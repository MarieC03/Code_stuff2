#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import csv
import os
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple, Optional

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import AutoMinorLocator
from scipy.signal import savgol_filter, get_window, detrend, find_peaks
from kuibit.simdir import SimDir


rcParams['text.usetex'] = True
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern']
rcParams['axes.labelsize'] = 12
rcParams['axes.linewidth'] = 1.0
rcParams['lines.markersize'] = 6
rcParams['xtick.labelsize'] = 14
rcParams['ytick.labelsize'] = 14

plt.rcParams.update({
    "font.family": "serif",
    "text.usetex": True,
    "pgf.rcfonts": False,
    "pgf.preamble": "\n".join([
        r"\usepackage{unicode-math}",
        r"\setmainfont{Computer Modern}",
    ])
})

MSUN_TIME_S = 4.925490947e-6
MSUN_TIME_MS = MSUN_TIME_S * 1e3
MPC_TO_LMSOL = 2.089152052322275e19

AXIS_FONTSIZE = 18
LEGEND_FONTSIZE = 16
TICK_FONTSIZE = 16

plt.rcParams.update({
    "font.size": AXIS_FONTSIZE,
    "axes.labelsize": AXIS_FONTSIZE,
    "xtick.labelsize": TICK_FONTSIZE,
    "ytick.labelsize": TICK_FONTSIZE,
    "legend.fontsize": LEGEND_FONTSIZE,
})


@dataclass
class ComplexSeries:
    t_code: np.ndarray
    y: np.ndarray

    def sorted_unique(self) -> "ComplexSeries":
        t = np.asarray(self.t_code, dtype=float)
        y = np.asarray(self.y, dtype=complex)
        order = np.argsort(t)
        t = t[order]
        y = y[order]
        keep = np.ones(len(t), dtype=bool)
        if len(t) > 1:
            keep[1:] = np.diff(t) > 0.0
        return ComplexSeries(t[keep], y[keep])

    @property
    def t_ms(self) -> np.ndarray:
        return self.t_code * MSUN_TIME_MS

    @property
    def dt_code(self) -> float:
        if len(self.t_code) < 2:
            raise ValueError("Need at least two time samples.")
        return float(np.median(np.diff(self.t_code)))

    @property
    def dt_s(self) -> float:
        return self.dt_code * MSUN_TIME_S


@dataclass
class SpectrumResult:
    f_khz: np.ndarray
    psd: np.ndarray
    char_strain: np.ndarray
    eff_strain: np.ndarray
    dominant_peak_khz: float
    peak_table: List[Tuple[float, float]]


@dataclass
class GWFluxHistory:
    t_ms_rel: np.ndarray
    e_gw: np.ndarray
    j_gw_z: np.ndarray


@dataclass
class GWProducts:
    psi4_22: ComplexSeries
    strain_22: ComplexSeries
    hp: np.ndarray
    hc: np.ndarray
    h_amp: np.ndarray
    merge_time_code: float
    merge_time_ms: float
    postmerger_spectrum: SpectrumResult
    inst_freq_khz: np.ndarray
    inst_freq_time_ms: np.ndarray
    inst_freq_valid: np.ndarray
    energy_mode_only: float
    angular_momentum_mode_only: float
    flux_history: GWFluxHistory
    energy_modes_used: List[Tuple[int, int]]


def load_psi4_mode(sim_path: str, det_radius: float, l: int, m: int) -> ComplexSeries:
    sim_dir = SimDir(sim_path)
    psi4 = sim_dir.gws[det_radius].get_psi4_lm(l, m)
    return ComplexSeries(
        np.asarray(psi4.t, dtype=float),
        np.asarray(psi4.y, dtype=complex),
    ).sorted_unique()


def merge_early_signal_from_reference(ref: ComplexSeries, restarted: ComplexSeries) -> ComplexSeries:
    if len(restarted.t_code) == 0:
        raise ValueError("Restarted signal is empty.")
    if len(ref.t_code) == 0:
        raise ValueError("Reference signal is empty.")
    t0 = restarted.t_code[0]
    early = ref.t_code < t0
    t = np.concatenate([ref.t_code[early], restarted.t_code])
    y = np.concatenate([ref.y[early], restarted.y])
    return ComplexSeries(t, y).sorted_unique()


def regularize_complex_series(series: ComplexSeries) -> ComplexSeries:
    s = series.sorted_unique()
    t = s.t_code
    y = s.y
    if len(t) < 4:
        return s
    dt = np.diff(t)
    dt_med = float(np.median(dt))
    if np.allclose(dt, dt_med, rtol=1e-5, atol=1e-12):
        return s
    t_uniform = np.arange(t[0], t[-1] + 0.5 * dt_med, dt_med)
    y_real = np.interp(t_uniform, t, y.real)
    y_imag = np.interp(t_uniform, t, y.imag)
    return ComplexSeries(t_uniform, y_real + 1j * y_imag)


def _window_array(kind: str, n: int) -> np.ndarray:
    name = (kind or "tukey").strip().lower()
    if name in ("none", "rect", "rectangular", "boxcar"):
        return np.ones(n, dtype=float)
    if name in ("hann", "hanning"):
        return get_window("hann", n)
    if name == "blackman":
        return get_window("blackman", n)
    return get_window(("tukey", 0.2), n)


def trim_to_valid_support(series: ComplexSeries, frac_of_peak: float = 1e-3, pad_points: int = 0) -> ComplexSeries:
    amp = np.abs(series.y)
    if len(amp) == 0:
        return series
    thr = frac_of_peak * float(np.max(amp))
    idx = np.where(amp >= thr)[0]
    if len(idx) == 0:
        return series
    i0 = max(0, int(idx[0]) - pad_points)
    i1 = min(len(amp), int(idx[-1]) + 1 + pad_points)
    return ComplexSeries(series.t_code[i0:i1], series.y[i0:i1])


def remove_edge_baseline(series: ComplexSeries, edge_points: int) -> ComplexSeries:
    y = np.asarray(series.y, dtype=complex).copy()
    n = len(y)
    if n < 4:
        return series
    m = max(2, min(edge_points, n // 4))
    x = np.linspace(0.0, 1.0, n)
    r0 = np.mean(y.real[:m])
    r1 = np.mean(y.real[-m:])
    i0 = np.mean(y.imag[:m])
    i1 = np.mean(y.imag[-m:])
    baseline = (r0 + (r1 - r0) * x) + 1j * (i0 + (i1 - i0) * x)
    y -= baseline
    return ComplexSeries(series.t_code, y)


def remove_slow_baseline(series: ComplexSeries, poly_order: int = 3) -> ComplexSeries:
    n = len(series.t_code)
    if n < max(20, poly_order + 5):
        return series
    x = np.linspace(-1.0, 1.0, n)
    cr = np.polyfit(x, series.y.real, deg=poly_order)
    ci = np.polyfit(x, series.y.imag, deg=poly_order)
    trend = np.polyval(cr, x) + 1j * np.polyval(ci, x)
    return ComplexSeries(series.t_code, series.y - trend)


def ffi_strain_from_psi4(
    psi4: ComplexSeries,
    cutoff_period_code: float,
    cutoff_fraction: float,
    window_function: str,
    support_threshold_frac: float = 1e-3,
) -> ComplexSeries:
    s = regularize_complex_series(psi4).sorted_unique()
    if len(s.t_code) < 8:
        return ComplexSeries(s.t_code, np.zeros_like(s.y))

    pad_pts = max(64, int(2.0 / max(s.dt_code, 1e-12)))
    s = trim_to_valid_support(s, frac_of_peak=support_threshold_frac, pad_points=pad_pts)

    t = s.t_code
    y = np.asarray(s.y, dtype=complex).copy()
    y = detrend(y.real, type="linear") + 1j * detrend(y.imag, type="linear")

    n = len(y)
    dt = float(np.median(np.diff(t)))
    y_in = y * _window_array(window_function, n)

    freqs_code = np.fft.fftfreq(n, d=dt)
    omega = 2.0 * np.pi * freqs_code
    omega0 = max(1e-12, 2.0 * np.pi * (cutoff_fraction / max(cutoff_period_code, 1e-300)))
    omega_eff = np.sign(omega) * np.maximum(np.abs(omega), omega0)
    omega_eff[0] = omega0

    psi4_fft = np.fft.fft(y_in)
    h_fft = -psi4_fft / (omega_eff ** 2)
    h = np.fft.ifft(h_fft)

    out = ComplexSeries(t, h)
    out = remove_edge_baseline(out, max(8, n // 20))
    out = remove_slow_baseline(out, 3)
    return out


def scale_strain_to_distance(strain_rh: ComplexSeries, distance_mpc: float) -> ComplexSeries:
    factor = 1.0 / (distance_mpc * MPC_TO_LMSOL)
    return ComplexSeries(strain_rh.t_code.copy(), strain_rh.y * factor)


def polarizations_from_scaled_strain(strain_h: ComplexSeries) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    hp = strain_h.y.real
    hc = -strain_h.y.imag
    h_amp = np.sqrt(hp**2 + hc**2)
    return hp, hc, h_amp


def crop_time_window(series: ComplexSeries, tmin_code: float, tmax_code: Optional[float]) -> ComplexSeries:
    mask = series.t_code >= tmin_code
    if tmax_code is not None:
        mask &= series.t_code <= tmax_code
    return ComplexSeries(series.t_code[mask], series.y[mask])


def smooth_curve(y: np.ndarray, window: int, poly: int, log_space: bool = False) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n < 7:
        return y.copy()
    w = int(window)
    if w % 2 == 0:
        w += 1
    w = min(w, n if n % 2 == 1 else n - 1)
    w = max(5, w)
    if w >= n:
        return y.copy()
    p = min(poly, w - 2)
    if log_space:
        z = np.log10(np.maximum(y, 1e-300))
        return 10.0 ** savgol_filter(z, w, p)
    out = savgol_filter(y, w, p)
    return np.clip(out, 0.0, None)


def phase_frequency_from_series(series: ComplexSeries, amp_threshold_frac: float, smooth_window: int, smooth_poly: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    s = regularize_complex_series(series).sorted_unique()
    if len(s.t_code) < 8:
        return s.t_code, np.full(len(s.t_code), np.nan), np.zeros(len(s.t_code), dtype=bool)

    amp = np.abs(s.y)
    thr = amp_threshold_frac * float(np.max(amp))
    valid = amp >= thr

    phase = np.unwrap(np.angle(s.y))
    w = int(smooth_window)
    if w % 2 == 0:
        w += 1
    if w >= len(phase):
        w = len(phase) - 1 if len(phase) % 2 == 0 else len(phase)
    if w >= 5:
        phase = savgol_filter(phase, w, min(smooth_poly, w - 2))

    omega = np.gradient(phase, s.t_code * MSUN_TIME_S)
    med_omega = np.nanmedian(omega[valid]) if np.any(valid) else np.nan
    if np.isfinite(med_omega) and med_omega < 0.0:
        omega = -omega

    f_hz = np.abs(omega) / (2.0 * np.pi)
    f_khz = f_hz / 1e3
    valid &= np.isfinite(f_khz) & (f_khz > 0.5) & (f_khz < 6.5)
    return s.t_code, f_khz, valid


def merger_time_from_strain_peak(strain_rh: ComplexSeries) -> Tuple[float, float]:
    s = regularize_complex_series(strain_rh).sorted_unique()
    amp = np.abs(s.y)
    if len(amp) == 0:
        raise ValueError("Empty strain series.")
    thr = 0.15 * float(np.max(amp))
    idx = np.where(amp >= thr)[0]
    if len(idx) == 0:
        i_peak = int(np.argmax(amp))
    else:
        i0 = int(idx[0])
        i_peak = i0 + int(np.argmax(amp[i0:]))
    t_merge = float(s.t_code[i_peak])
    return t_merge, t_merge * MSUN_TIME_MS


def instantaneous_frequency_from_phase(psi4_22: ComplexSeries, tmerge_code: float, pre_ms: float, post_ms: float, amp_threshold_frac: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    tmin = tmerge_code - pre_ms / MSUN_TIME_MS
    tmax = tmerge_code + post_ms / MSUN_TIME_MS
    sub = crop_time_window(psi4_22, tmin, tmax)
    t_code, f_khz, valid = phase_frequency_from_series(
        sub,
        amp_threshold_frac=amp_threshold_frac,
        smooth_window=11,
        smooth_poly=2,
    )
    t_ms = t_code * MSUN_TIME_MS - tmerge_code * MSUN_TIME_MS

    # Late-time spikes come from differentiating the phase where the signal amplitude
    # is too small. We therefore only interpolate inside the trusted interval and
    # never extrapolate beyond the last trusted point.
    amp = np.abs(sub.y)
    amp_thr = max(amp_threshold_frac * float(np.max(amp)), 0.0)
    trusted = np.isfinite(f_khz) & (f_khz > 0.5) & (f_khz < 6.5) & (amp >= amp_thr)

    f_plot = np.full_like(f_khz, np.nan, dtype=float)
    if np.count_nonzero(trusted) >= 4:
        idx = np.arange(len(f_khz), dtype=float)
        i_first = int(np.where(trusted)[0][0])
        i_last = int(np.where(trusted)[0][-1])
        support = np.arange(i_first, i_last + 1)
        f_support = np.interp(support, idx[trusted], f_khz[trusted])
        if len(f_support) >= 7:
            f_support = smooth_curve(f_support, window=9, poly=2, log_space=False)
        f_plot[support] = f_support

    return t_ms, f_plot, trusted


def tapered_periodogram(x: np.ndarray, dt_s: float, window_name: str = "hann", pad_factor: int = 32) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = len(x)
    nfft = int(max(2, pad_factor) * n)
    win = get_window(window_name, n)
    xw = x * win
    fft = np.fft.rfft(xw, n=nfft)
    f_hz = np.fft.rfftfreq(nfft, d=dt_s)
    fs = 1.0 / dt_s
    psd = (np.abs(fft) ** 2) / (fs * np.sum(win ** 2))
    if len(psd) > 2:
        psd[1:-1] *= 2.0
    htilde = dt_s * fft
    return f_hz, psd, htilde


def compute_spectrum(
    strain_h: ComplexSeries,
    tmerge_code: float,
    tstart_offset_ms: float,
    tend_offset_ms: Optional[float],
    fmin_khz: float,
    fmax_khz: float,
    peak_search_min_khz: float,
    spectrum_smooth_window: int,
    spectrum_smooth_poly: int,
    pad_factor: int,
) -> SpectrumResult:
    tmin = tmerge_code + tstart_offset_ms / MSUN_TIME_MS
    tmax = None if tend_offset_ms is None else tmerge_code + tend_offset_ms / MSUN_TIME_MS
    post = regularize_complex_series(crop_time_window(strain_h, tmin, tmax))
    if len(post.t_code) < 64:
        raise ValueError("Post-merger segment is too short for spectral analysis.")

    x = detrend(post.y.real.copy(), type="linear")
    f_hz, psd_raw, htilde = tapered_periodogram(x, post.dt_s, window_name="hann", pad_factor=pad_factor)
    f_khz = f_hz / 1e3

    band = (f_khz >= fmin_khz) & (f_khz <= fmax_khz)
    f_band = f_khz[band]
    psd_band = psd_raw[band]
    psd_plot = smooth_curve(psd_band, spectrum_smooth_window, spectrum_smooth_poly, log_space=True)

    hchar = 2.0 * f_hz * np.abs(htilde)
    heff = 2.0 * np.sqrt(np.maximum(f_hz, 0.0)) * np.abs(htilde)

    hchar_band = hchar[band]
    heff_band = heff[band]

    hchar_plot = smooth_curve(hchar_band, max(5, spectrum_smooth_window), spectrum_smooth_poly, log_space=True)
    heff_plot = smooth_curve(heff_band, max(5, spectrum_smooth_window), spectrum_smooth_poly, log_space=True)

    search = hchar_band.copy()
    search[f_band < max(peak_search_min_khz, fmin_khz)] = 0.0
    prominence = 0.03 * np.max(search) if np.max(search) > 0 else 0.0
    pidx, _ = find_peaks(search, prominence=prominence, distance=max(3, len(search) // 60))
    if len(pidx) == 0:
        dominant_idx = int(np.argmax(search))
        pidx = np.array([dominant_idx], dtype=int)
    else:
        dominant_idx = int(pidx[np.argmax(search[pidx])])

    peak_table = sorted(
        [(float(f_band[i]), float(hchar_band[i])) for i in pidx],
        key=lambda item: item[1],
        reverse=True,
    )

    return SpectrumResult(
        f_khz=f_band,
        psd=psd_plot,
        char_strain=hchar_plot,
        eff_strain=heff_plot,
        dominant_peak_khz=float(f_band[dominant_idx]),
        peak_table=peak_table,
    )


def cumtrapz_np(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    out = np.zeros_like(y, dtype=float)
    if len(y) < 2:
        return out
    dx = np.diff(x)
    out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * dx)
    return out


def compute_energy_angular_momentum_history(
    mode_series: Dict[Tuple[int, int], ComplexSeries],
    cutoff_period_code: float,
    cutoff_fraction: float,
    window_function: str,
    tmerge_code: float,
) -> Tuple[GWFluxHistory, float, float]:
    regularized_h: Dict[Tuple[int, int], ComplexSeries] = {}
    for mode, psi4 in mode_series.items():
        h_r = regularize_complex_series(
            ffi_strain_from_psi4(psi4, cutoff_period_code, cutoff_fraction, window_function)
        )
        regularized_h[mode] = h_r

    ref_mode = list(regularized_h.keys())[0]
    t_ref = regularized_h[ref_mode].t_code
    eflux_total = np.zeros_like(t_ref, dtype=float)
    jflux_total = np.zeros_like(t_ref, dtype=float)

    for (l, m), h_r in regularized_h.items():
        if not np.allclose(h_r.t_code, t_ref):
            hr_real = np.interp(t_ref, h_r.t_code, h_r.y.real)
            hr_imag = np.interp(t_ref, h_r.t_code, h_r.y.imag)
            h_y = hr_real + 1j * hr_imag
        else:
            h_y = h_r.y

        news = np.gradient(h_y, t_ref)
        eflux = np.abs(news) ** 2 / (16.0 * np.pi)
        jflux = (m / (16.0 * np.pi)) * np.imag(h_y * np.conjugate(news))

        eflux_total += eflux
        jflux_total += jflux

    e_cum = cumtrapz_np(eflux_total, t_ref)
    j_cum = cumtrapz_np(jflux_total, t_ref)
    t_rel_ms = t_ref * MSUN_TIME_MS - tmerge_code * MSUN_TIME_MS

    history = GWFluxHistory(
        t_ms_rel=t_rel_ms,
        e_gw=e_cum,
        j_gw_z=j_cum,
    )
    return history, float(e_cum[-1]), float(j_cum[-1])


def analyse_simulation_from_psi4(
    psi4_22: ComplexSeries,
    mode_series: Dict[Tuple[int, int], ComplexSeries],
    distance_mpc: float,
    gomega: float,
    cutoff_fraction: float,
    window_function: str,
    tstart_offset_ms: float,
    tend_offset_ms: Optional[float],
    fmin_khz: float,
    fmax_khz: float,
    peak_search_min_khz: float,
    spectrum_smooth_window: int,
    spectrum_smooth_poly: int,
    pad_factor: int,
    instfreq_pre_ms: float,
    instfreq_post_ms: float,
    instfreq_amp_threshold_frac: float,
    trim_strain_plot_frac: float,
) -> GWProducts:
    cutoff_period_code = 2.0 * np.pi / gomega

    strain_rh = ffi_strain_from_psi4(psi4_22, cutoff_period_code, cutoff_fraction, window_function)
    tmerge_code, tmerge_ms = merger_time_from_strain_peak(strain_rh)

    strain_h_full = scale_strain_to_distance(strain_rh, distance_mpc)

    pad_points = max(64, int(2.0 / max(strain_h_full.dt_code, 1e-12)))
    psi4_trim = trim_to_valid_support(regularize_complex_series(psi4_22), frac_of_peak=1e-3, pad_points=pad_points)
    mask = (strain_h_full.t_code >= psi4_trim.t_code[0]) & (strain_h_full.t_code <= psi4_trim.t_code[-1])
    strain_plot = ComplexSeries(strain_h_full.t_code[mask], strain_h_full.y[mask])
    strain_plot = trim_to_valid_support(strain_plot, frac_of_peak=trim_strain_plot_frac, pad_points=pad_points)

    hp, hc, h_amp = polarizations_from_scaled_strain(strain_plot)
    spec = compute_spectrum(
        strain_h_full, tmerge_code, tstart_offset_ms, tend_offset_ms, fmin_khz, fmax_khz,
        peak_search_min_khz, spectrum_smooth_window, spectrum_smooth_poly, pad_factor
    )
    tf_ms, f_inst_khz, inst_valid = instantaneous_frequency_from_phase(
        psi4_22, tmerge_code, instfreq_pre_ms, instfreq_post_ms, instfreq_amp_threshold_frac
    )
    flux_history, energy, angmom = compute_energy_angular_momentum_history(
        mode_series, cutoff_period_code, cutoff_fraction, window_function, tmerge_code
    )

    return GWProducts(
        psi4_22=psi4_22,
        strain_22=strain_plot,
        hp=hp,
        hc=hc,
        h_amp=h_amp,
        merge_time_code=tmerge_code,
        merge_time_ms=tmerge_ms,
        postmerger_spectrum=spec,
        inst_freq_khz=f_inst_khz,
        inst_freq_time_ms=tf_ms,
        inst_freq_valid=inst_valid,
        energy_mode_only=energy,
        angular_momentum_mode_only=angmom,
        flux_history=flux_history,
        energy_modes_used=list(mode_series.keys()),
    )


def enforce_identical_premerger_strain(m1: GWProducts, nom1: GWProducts) -> None:
    """
    Force the pre-merger strain to be identical between M1 and no-M1 when the
    no-M1 signal has been prepended with the early inspiral taken from the M1 run.
    """
    tm = m1.merge_time_ms
    t_m1 = m1.strain_22.t_ms - tm
    t_nu = nom1.strain_22.t_ms - tm

    pre_mask = t_nu < 0.0
    if not np.any(pre_mask):
        return

    nom1.hp[pre_mask] = np.interp(t_nu[pre_mask], t_m1, m1.hp)
    nom1.hc[pre_mask] = np.interp(t_nu[pre_mask], t_m1, m1.hc)
    nom1.h_amp = np.sqrt(nom1.hp**2 + nom1.hc**2)
    nom1.strain_22 = ComplexSeries(nom1.strain_22.t_code, nom1.hp - 1j * nom1.hc)


def _enforce_identical_premerger_cumulative_series(
    t_ref: np.ndarray,
    y_ref: np.ndarray,
    t_target: np.ndarray,
    y_target: np.ndarray,
) -> np.ndarray:
    """
    Replace the pre-merger cumulative history in ``y_target`` by the reference
    history sampled from ``y_ref``. For post-merger times, remove the constant
    offset at merger so the curve remains continuous.
    """
    y_out = np.asarray(y_target, dtype=float).copy()
    pre_mask = t_target < 0.0
    if not np.any(pre_mask):
        return y_out

    y_out[pre_mask] = np.interp(t_target[pre_mask], t_ref, y_ref)

    post_mask = ~pre_mask
    if np.any(post_mask):
        t0 = float(t_target[post_mask][0])
        ref_at_t0 = float(np.interp(t0, t_ref, y_ref))
        offset = y_out[post_mask][0] - ref_at_t0
        y_out[post_mask] -= offset

    return y_out


def enforce_identical_premerger_flux_history(m1: GWProducts, nom1: GWProducts) -> None:
    """
    When the no-M1 run uses the M1 inspiral before merger, the cumulative GW
    energy and angular-momentum histories must also coincide before merger.
    Small mismatches otherwise arise because the FFI/baseline-removal steps use
    the full time series and can imprint tiny post-merger differences back onto
    the inspiral part.
    """
    t_m1 = np.asarray(m1.flux_history.t_ms_rel, dtype=float)
    t_nu = np.asarray(nom1.flux_history.t_ms_rel, dtype=float)

    nom1_e = _enforce_identical_premerger_cumulative_series(
        t_m1, m1.flux_history.e_gw, t_nu, nom1.flux_history.e_gw
    )
    nom1_j = _enforce_identical_premerger_cumulative_series(
        t_m1, m1.flux_history.j_gw_z, t_nu, nom1.flux_history.j_gw_z
    )

    nom1.flux_history = GWFluxHistory(
        t_ms_rel=nom1.flux_history.t_ms_rel,
        e_gw=nom1_e,
        j_gw_z=nom1_j,
    )
    nom1.energy_mode_only = float(nom1_e[-1])
    nom1.angular_momentum_mode_only = float(nom1_j[-1])


def style_axes(ax, x_minor: bool = True, y_minor: bool = True) -> None:
    ax.tick_params(axis="both", which="major", direction="in", top=True, right=True,
                   length=7, width=1.0, labelsize=TICK_FONTSIZE)
    ax.tick_params(axis="both", which="minor", direction="in", top=True, right=True,
                   length=4, width=0.8)
    if x_minor:
        ax.xaxis.set_minor_locator(AutoMinorLocator())
    if y_minor:
        ax.yaxis.set_minor_locator(AutoMinorLocator())


def make_square_figure():
    return plt.subplots(figsize=(7.2, 7.2), dpi=160)


def plot_psi4_norm(products: Dict[str, GWProducts], outpath: str, colours: Dict[str, str],
                  xmin_ms: float = -10.0, xmax_ms: Optional[float] = 36.0) -> None:
    fig, ax = make_square_figure()
    for label, prod in products.items():
        t_rel = prod.psi4_22.t_ms - prod.merge_time_ms
        ax.plot(t_rel, np.abs(prod.psi4_22.y), label=label, linewidth=1.7, color=colours[label])
    ax.set_xlim(xmin_ms, xmax_ms)
    ax.set_xlabel(r"$t-t_{\rm mer}\,[\mathrm{ms}]$", fontsize=AXIS_FONTSIZE)
    ax.set_ylabel(r"$|\Psi_4^{22}|$", fontsize=AXIS_FONTSIZE)
    style_axes(ax)
    ax.legend(frameon=False, fontsize=LEGEND_FONTSIZE)
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def plot_strain(products: Dict[str, GWProducts], outpath: str, colours: Dict[str, str],
                xmin_ms: float = -10.0, xmax_ms: Optional[float] = 36.0) -> None:
    fig, ax = make_square_figure()
    for label, prod in products.items():
        t = prod.strain_22.t_ms - prod.merge_time_ms
        ax.plot(t, prod.hp * 1e21, label=label, linewidth=1.7, color=colours[label])
    ax.set_xlim(xmin_ms, xmax_ms)
    ax.set_xlabel(r"$t-t_{\rm mer}\,[\mathrm{ms}]$", fontsize=AXIS_FONTSIZE)
    ax.set_ylabel(r"$h_+^{22}\times 10^{21}$", fontsize=AXIS_FONTSIZE)
    style_axes(ax)
    ax.legend(frameon=False, fontsize=LEGEND_FONTSIZE)
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def plot_postmerger_psd(products: Dict[str, GWProducts], outpath: str, colours: Dict[str, str]) -> None:
    fig, ax = make_square_figure()
    for label, prod in products.items():
        c = colours[label]
        ax.semilogy(prod.postmerger_spectrum.f_khz, prod.postmerger_spectrum.psd,
                    label=label, linewidth=1.8, color=c)
        ax.axvline(prod.postmerger_spectrum.dominant_peak_khz, linestyle="--", linewidth=1.2, color=c)
    ax.set_xlabel(r"$f_{\rm GW}\,[\mathrm{kHz}]$", fontsize=AXIS_FONTSIZE)
    ax.set_ylabel(r"PSD $[1/\mathrm{Hz}]$", fontsize=AXIS_FONTSIZE)
    style_axes(ax, y_minor=False)
    ax.legend(frameon=False, fontsize=LEGEND_FONTSIZE)
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def plot_characteristic_spectrum(products: Dict[str, GWProducts], outpath: str, colours: Dict[str, str]) -> None:
    fig, ax = make_square_figure()
    for label, prod in products.items():
        c = colours[label]
        ax.semilogy(prod.postmerger_spectrum.f_khz, prod.postmerger_spectrum.char_strain,
                    label=label, linewidth=1.8, color=c)
        ax.axvline(prod.postmerger_spectrum.dominant_peak_khz, linestyle="--", linewidth=1.2, color=c)
    ax.set_xlabel(r"$f_{\rm GW}\,[\mathrm{kHz}]$", fontsize=AXIS_FONTSIZE)
    ax.set_ylabel(r"$2 f |\tilde h_+^{22}(f)|$", fontsize=AXIS_FONTSIZE)
    style_axes(ax, y_minor=False)
    ax.legend(frameon=False, fontsize=LEGEND_FONTSIZE)
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def plot_effective_spectrum(products: Dict[str, GWProducts], outpath: str, colours: Dict[str, str]) -> None:
    fig, ax = make_square_figure()
    for label, prod in products.items():
        c = colours[label]
        ax.semilogy(prod.postmerger_spectrum.f_khz, prod.postmerger_spectrum.eff_strain,
                    label=label, linewidth=1.8, color=c)
        ax.axvline(prod.postmerger_spectrum.dominant_peak_khz, linestyle="--", linewidth=1.2, color=c)
    ax.set_xlabel(r"$f\,[\mathrm{kHz}]$", fontsize=AXIS_FONTSIZE)
    ax.set_ylabel(r"$2 \sqrt{f}\, |\tilde h_+^{22}(f)|\,[\mathrm{Hz}^{-1/2}]$", fontsize=AXIS_FONTSIZE)
    style_axes(ax, y_minor=False)
    ax.legend(frameon=False, fontsize=LEGEND_FONTSIZE)
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def plot_instantaneous_frequency(products: Dict[str, GWProducts], outpath: str, colours: Dict[str, str],
                               xmin_ms: float = -10.0, xmax_ms: Optional[float] = 15.0,
                               ymax_khz: Optional[float] = 3.8) -> None:
    fig, ax = make_square_figure()
    for label, prod in products.items():
        ax.plot(prod.inst_freq_time_ms, prod.inst_freq_khz, label=label,
                linewidth=1.7, color=colours[label])
    ax.set_xlim(xmin_ms, xmax_ms)
    if ymax_khz is not None:
        ax.set_ylim(top=ymax_khz)
    ax.set_xlabel(r"$t-t_{\rm mer}\,[\mathrm{ms}]$", fontsize=AXIS_FONTSIZE)
    ax.set_ylabel(r"$f_{\rm GW}(t)\,[\mathrm{kHz}]$", fontsize=AXIS_FONTSIZE)
    style_axes(ax)
    ax.legend(frameon=False, fontsize=LEGEND_FONTSIZE)
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def plot_instantaneous_frequency_extended(products: Dict[str, GWProducts], outpath: str, colours: Dict[str, str],
                                        xmin_ms: float = -10.0, xmax_ms: Optional[float] = 36.0) -> None:
    fig, ax = make_square_figure()
    for label, prod in products.items():
        ax.plot(prod.inst_freq_time_ms, prod.inst_freq_khz, label=label,
                linewidth=1.7, color=colours[label])
    ax.set_xlim(xmin_ms, xmax_ms)
    ax.set_xlabel(r"$t-t_{\rm mer}\,[\mathrm{ms}]$", fontsize=AXIS_FONTSIZE)
    ax.set_ylabel(r"$f_{\rm GW}(t)\,[\mathrm{kHz}]$", fontsize=AXIS_FONTSIZE)
    style_axes(ax)
    ax.legend(frameon=False, fontsize=LEGEND_FONTSIZE)
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)

def plot_energy_vs_time(products: Dict[str, GWProducts], outpath: str, colours: Dict[str, str],
                        xmin_ms: float = -10.0, xmax_ms: Optional[float] = 36.0) -> None:
    fig, ax = make_square_figure()
    for label, prod in products.items():
        ax.plot(prod.flux_history.t_ms_rel, prod.flux_history.e_gw,
                label=label, linewidth=1.7, color=colours[label])
    ax.set_xlim(xmin_ms, xmax_ms)
    ax.set_xlabel(r"$t-t_{\rm mer}\,[\mathrm{ms}]$", fontsize=AXIS_FONTSIZE)
    ax.set_ylabel(r"$E_{\rm GW}(t)\,[M_\odot]$", fontsize=AXIS_FONTSIZE)
    style_axes(ax)
    ax.legend(frameon=False, fontsize=LEGEND_FONTSIZE)
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def plot_angmom_vs_time(products: Dict[str, GWProducts], outpath: str, colours: Dict[str, str],
                        xmin_ms: float = -10.0, xmax_ms: Optional[float] = 36.0) -> None:
    fig, ax = make_square_figure()
    for label, prod in products.items():
        ax.plot(prod.flux_history.t_ms_rel, prod.flux_history.j_gw_z,
                label=label, linewidth=1.7, color=colours[label])
    ax.set_xlim(xmin_ms, xmax_ms)
    ax.set_xlabel(r"$t-t_{\rm mer}\,[\mathrm{ms}]$", fontsize=AXIS_FONTSIZE)
    ax.set_ylabel(r"$J_{{\rm GW},z}(t)\,[M_\odot^2]$", fontsize=AXIS_FONTSIZE)
    style_axes(ax)
    ax.legend(frameon=False, fontsize=LEGEND_FONTSIZE)
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def write_summary_csv(products: Dict[str, GWProducts], out_csv: str) -> None:
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "simulation", "t_merge_ms", "dominant_postmerger_peak_kHz",
            "energy_radiated_Msun", "angular_momentum_radiated_Msun2",
            "energy_modes_used", "top_peak_1_kHz", "top_peak_2_kHz", "top_peak_3_kHz"
        ])
        for label, prod in products.items():
            peaks = [f"{pk[0]:.6f}" for pk in prod.postmerger_spectrum.peak_table[:3]]
            while len(peaks) < 3:
                peaks.append("")
            writer.writerow([
                label,
                f"{prod.merge_time_ms:.6f}",
                f"{prod.postmerger_spectrum.dominant_peak_khz:.6f}",
                f"{prod.energy_mode_only:.10e}",
                f"{prod.angular_momentum_mode_only:.10e}",
                " ".join([f"({l},{m})" for l, m in prod.energy_modes_used]),
                *peaks
            ])


def save_dat(path: str, header: str, columns: List[np.ndarray]) -> None:
    arr = np.column_stack(columns)
    np.savetxt(path, arr, header=header)


def main() -> None:
    parser = argparse.ArgumentParser(description="GW diagnostics for FIL M1/no-M1 from Psi4.")
    parser.add_argument("--root", required=True)
    parser.add_argument("--sim-m1", required=True)
    parser.add_argument("--sim-nom1", required=True)
    parser.add_argument("--dist", type=float, default=500.0)
    parser.add_argument("--distance-mpc", type=float, default=40.0)
    parser.add_argument("--gomega", type=float, default=0.009)
    parser.add_argument("--cutoff-fraction", type=float, default=0.001)
    parser.add_argument("--window-function", default="tukey", choices=["tukey", "hann", "hanning", "blackman", "none"])
    parser.add_argument("--prepend-early-m1-to-nom1", action="store_true")
    parser.add_argument("--postmerger-start-ms", type=float, default=0.0)
    parser.add_argument("--postmerger-end-ms", type=float, default=15.0)
    parser.add_argument("--fmin-khz", type=float, default=1.0)
    parser.add_argument("--fmax-khz", type=float, default=5.0)
    parser.add_argument("--peak-search-min-khz", type=float, default=1.2)
    parser.add_argument("--spectrum-smooth-window", type=int, default=5)
    parser.add_argument("--spectrum-smooth-poly", type=int, default=2)
    parser.add_argument("--pad-factor", type=int, default=32)
    parser.add_argument("--instfreq-pre-ms", type=float, default=3.0)
    parser.add_argument("--instfreq-post-ms", type=float, default=15.0)
    parser.add_argument("--time-plot-min-ms", type=float, default=-10.0,
                        help="Left x-limit for time-domain plots relative to merger.")
    parser.add_argument("--time-plot-max-ms", type=float, default=36.0,
                        help="Right x-limit for extended time-domain plots relative to merger.")
    parser.add_argument("--instfreq-plot-max-ms", type=float, default=36.0,
                        help="Right x-limit for the extended instantaneous-frequency plot.")
    parser.add_argument("--instfreq-ymax-khz", type=float, default=3.8,
                        help="Upper y-limit for the non-extended instantaneous-frequency plot.")
    parser.add_argument("--flux-smooth-frac", type=float, default=0.03,
                        help="Backward-compatible option kept for old command lines; not used in this version.")
    parser.add_argument("--flux-smooth-min", type=int, default=31,
                        help="Backward-compatible option kept for old command lines; not used in this version.")
    parser.add_argument("--instfreq-amp-threshold-frac", type=float, default=0.05)
    parser.add_argument("--trim-strain-plot-frac", type=float, default=0.01)
    parser.add_argument("--colour-M1", default="black")
    parser.add_argument("--colour-noM1", default="#D93C1E")
    parser.add_argument("--save-data", action="store_true")
    parser.add_argument("--energy-modes", nargs="*", default=["2,2"])
    parser.add_argument("--outdir", default="gw_analysis_output")
    parser.add_argument("--merger-time-code", type=float, default=None,
                        help="Optional manual override for merger time in code units (Msun).")
    args = parser.parse_args()

    sim_m1_path = args.sim_m1 if os.path.isabs(args.sim_m1) else os.path.join(args.root, args.sim_m1)
    sim_nom1_path = args.sim_nom1 if os.path.isabs(args.sim_nom1) else os.path.join(args.root, args.sim_nom1)
    os.makedirs(args.outdir, exist_ok=True)

    colours = {"M1": args.colour_M1, "no-M1": args.colour_noM1}

    requested_modes = []
    for item in args.energy_modes:
        l, m = item.split(",")
        requested_modes.append((int(l), int(m)))
    if (2, 2) not in requested_modes:
        requested_modes = [(2, 2)] + requested_modes

    m1_modes = {mode: load_psi4_mode(sim_m1_path, args.dist, mode[0], mode[1]) for mode in requested_modes}
    nom1_modes = {mode: load_psi4_mode(sim_nom1_path, args.dist, mode[0], mode[1]) for mode in requested_modes}

    if args.prepend_early_m1_to_nom1:
        nom1_modes = {mode: merge_early_signal_from_reference(m1_modes[mode], nom1_modes[mode]) for mode in requested_modes}

    pm_end = None if args.postmerger_end_ms <= 0.0 else args.postmerger_end_ms
    instfreq_pre_ms_eff = max(args.instfreq_pre_ms, max(0.0, -args.time_plot_min_ms))
    instfreq_post_ms_eff = max(args.instfreq_post_ms, args.instfreq_plot_max_ms, args.time_plot_max_ms)

    common = dict(
        distance_mpc=args.distance_mpc,
        gomega=args.gomega,
        cutoff_fraction=args.cutoff_fraction,
        window_function=args.window_function,
        tstart_offset_ms=args.postmerger_start_ms,
        tend_offset_ms=pm_end,
        fmin_khz=args.fmin_khz,
        fmax_khz=args.fmax_khz,
        peak_search_min_khz=args.peak_search_min_khz,
        spectrum_smooth_window=args.spectrum_smooth_window,
        spectrum_smooth_poly=args.spectrum_smooth_poly,
        pad_factor=args.pad_factor,
        instfreq_pre_ms=instfreq_pre_ms_eff,
        instfreq_post_ms=instfreq_post_ms_eff,
        instfreq_amp_threshold_frac=min(args.instfreq_amp_threshold_frac, 5e-3),
        trim_strain_plot_frac=args.trim_strain_plot_frac,
    )

    products = {
        "M1": analyse_simulation_from_psi4(
            psi4_22=m1_modes[(2, 2)],
            mode_series={mode: m1_modes[mode] for mode in requested_modes},
            **common,
        ),
        "no-M1": analyse_simulation_from_psi4(
            psi4_22=nom1_modes[(2, 2)],
            mode_series={mode: nom1_modes[mode] for mode in requested_modes},
            **common,
        ),
    }

    if args.merger_time_code is not None:
        tmerge_code = float(args.merger_time_code)
        tmerge_ms = tmerge_code * MSUN_TIME_MS
        for key in products:
            old_ms = products[key].merge_time_ms
            products[key].merge_time_code = tmerge_code
            products[key].merge_time_ms = tmerge_ms
            products[key].flux_history = GWFluxHistory(
                t_ms_rel=products[key].flux_history.t_ms_rel + (old_ms - tmerge_ms),
                e_gw=products[key].flux_history.e_gw,
                j_gw_z=products[key].flux_history.j_gw_z,
            )

    if args.prepend_early_m1_to_nom1:
        products["no-M1"].merge_time_code = products["M1"].merge_time_code
        products["no-M1"].merge_time_ms = products["M1"].merge_time_ms
        enforce_identical_premerger_strain(products["M1"], products["no-M1"])
        enforce_identical_premerger_flux_history(products["M1"], products["no-M1"])

    plot_psi4_norm(products, os.path.join(args.outdir, "psi4_22_norm_comparison.png"), colours,
                  xmin_ms=args.time_plot_min_ms, xmax_ms=args.time_plot_max_ms)
    plot_strain(products, os.path.join(args.outdir, "strain_22_comparison.png"), colours,
                xmin_ms=args.time_plot_min_ms, xmax_ms=args.time_plot_max_ms)
    plot_postmerger_psd(products, os.path.join(args.outdir, "postmerger_psd_comparison.png"), colours)
    plot_characteristic_spectrum(products, os.path.join(args.outdir, "postmerger_amplitude_spectrum_comparison.png"), colours)
    plot_effective_spectrum(products, os.path.join(args.outdir, "postmerger_effective_spectrum_comparison.png"), colours)
    plot_instantaneous_frequency(products, os.path.join(args.outdir, "instantaneous_gw_frequency_comparison.png"), colours,
                               xmin_ms=args.time_plot_min_ms, xmax_ms=args.instfreq_plot_max_ms, ymax_khz=args.instfreq_ymax_khz)
    plot_instantaneous_frequency_extended(products, os.path.join(args.outdir, "instantaneous_gw_frequency_comparison_extended.png"), colours,
                                        xmin_ms=args.time_plot_min_ms, xmax_ms=args.instfreq_plot_max_ms)
    plot_energy_vs_time(products, os.path.join(args.outdir, "gw_energy_vs_time.png"), colours,
                       xmin_ms=args.time_plot_min_ms, xmax_ms=args.time_plot_max_ms)
    plot_angmom_vs_time(products, os.path.join(args.outdir, "gw_angular_momentum_vs_time.png"), colours,
                       xmin_ms=args.time_plot_min_ms, xmax_ms=args.time_plot_max_ms)
    write_summary_csv(products, os.path.join(args.outdir, "gw_summary.csv"))

    if args.save_data:
        for label, prod in products.items():
            safe = label.replace("-", "").replace(" ", "_")
            save_dat(os.path.join(args.outdir, f"psi4_22_norm_{safe}.dat"),
                     "t_minus_tmerger_ms abs_psi4_22",
                     [prod.psi4_22.t_ms - prod.merge_time_ms, np.abs(prod.psi4_22.y)])
            save_dat(os.path.join(args.outdir, f"strain_22_{safe}.dat"),
                     "t_minus_tmerger_ms hp_22 hc_22 h_amp_22",
                     [prod.strain_22.t_ms - prod.merge_time_ms, prod.hp, prod.hc, prod.h_amp])
            save_dat(os.path.join(args.outdir, f"f_GW_time_{safe}.dat"),
                     "t_minus_tmerger_ms f_GW_kHz",
                     [prod.inst_freq_time_ms, prod.inst_freq_khz])
            save_dat(os.path.join(args.outdir, f"PSD_{safe}.dat"),
                     "f_kHz PSD_1_per_Hz",
                     [prod.postmerger_spectrum.f_khz, prod.postmerger_spectrum.psd])
            save_dat(os.path.join(args.outdir, f"char_spectrum_{safe}.dat"),
                     "f_kHz two_f_abs_h_tilde",
                     [prod.postmerger_spectrum.f_khz, prod.postmerger_spectrum.char_strain])
            save_dat(os.path.join(args.outdir, f"effective_spectrum_{safe}.dat"),
                     "f_kHz two_sqrtf_abs_h_tilde",
                     [prod.postmerger_spectrum.f_khz, prod.postmerger_spectrum.eff_strain])
            save_dat(os.path.join(args.outdir, f"gw_energy_vs_time_{safe}.dat"),
                     "t_minus_tmerger_ms E_GW_Msun",
                     [prod.flux_history.t_ms_rel, prod.flux_history.e_gw])
            save_dat(os.path.join(args.outdir, f"gw_angular_momentum_vs_time_{safe}.dat"),
                     "t_minus_tmerger_ms J_GW_z_Msun2",
                     [prod.flux_history.t_ms_rel, prod.flux_history.j_gw_z])

    for label, prod in products.items():
        print(f"\n=== {label} ===")
        print(f"merge time [code]               : {prod.merge_time_code:.6f}")
        print(f"merge time [ms]                 : {prod.merge_time_ms:.6f}")
        print(f"dominant post-merger peak [kHz] : {prod.postmerger_spectrum.dominant_peak_khz:.6f}")
        print(f"radiated GW energy [Msun]       : {prod.energy_mode_only:.10e}")
        print(f"radiated GW ang. mom. [Msun^2]  : {prod.angular_momentum_mode_only:.10e}")
        for freq, val in prod.postmerger_spectrum.peak_table[:5]:
            print(f"  {freq:.6f} kHz   char={val:.6e}")

    print(f"\nSaved outputs in: {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()
