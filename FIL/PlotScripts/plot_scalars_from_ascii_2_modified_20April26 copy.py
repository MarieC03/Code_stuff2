#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import re
import argparse
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import transforms

from UtilityUnits import *

FIGSIZE_SQUARE = (6.8, 6.4)

plt.rcParams.update({
    "text.usetex": True,
    "font.size": 12,
    "font.family": "Serif",
    "font.serif": "Computer Modern",
})

# ============================================================
# Utilities
# ============================================================

def natural_key(path_obj):
    s = str(path_obj)
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)]


def sanitize_stem(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_\-]+", "_", str(name)).strip("_")


def ensure_time_ms(df: pd.DataFrame) -> pd.DataFrame:
    if df is None:
        return pd.DataFrame()
    if df.empty:
        if "time_ms" not in df.columns:
            df = df.copy()
            df["time_ms"] = pd.Series(dtype=float)
        return df
    if "time_ms" not in df.columns:
        df = df.copy()
        if "time_code" in df.columns:
            df["time_ms"] = df["time_code"] * CU_to_ms
        else:
            raise KeyError("DataFrame has neither 'time_ms' nor 'time_code'.")
    return df




def remove_isolated_log_spikes(x: np.ndarray, y: np.ndarray, *, log10_jump: float = 0.75,
                              min_factor_vs_neighbors: float = 5.0,
                              neighbor_match_log10: float = 0.35) -> np.ndarray:
    """Return mask that removes isolated one-point glitches in a positive series.

    A point is rejected if it is finite/positive but differs strongly from both
    neighboring points in log10-space while the two neighbors are mutually
    consistent. This removes both isolated upward overshoots and isolated
    downward drops that can appear near the end of corrupted runs.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y) & (y > 0.0)
    idx = np.where(keep)[0]
    if idx.size < 3:
        return keep
    for j in range(1, idx.size - 1):
        i = idx[j]
        il = idx[j - 1]
        ir = idx[j + 1]
        yi = y[i]
        yl = y[il]
        yr = y[ir]
        if yl <= 0.0 or yr <= 0.0:
            continue

        logy_i = np.log10(yi)
        logy_l = np.log10(yl)
        logy_r = np.log10(yr)
        neighbor_match = abs(logy_l - logy_r) <= neighbor_match_log10
        if not neighbor_match:
            continue

        dev_left = abs(logy_i - logy_l)
        dev_right = abs(logy_i - logy_r)
        if dev_left <= log10_jump or dev_right <= log10_jump:
            continue

        # keep broad physical peaks/valleys if the central point is not really isolated
        large_neighbor = max(yl, yr)
        small_neighbor = min(yl, yr)
        broad_peak = yi <= min_factor_vs_neighbors * large_neighbor
        broad_dip = yi >= small_neighbor / min_factor_vs_neighbors

        if yi > large_neighbor:
            if not broad_peak:
                keep[i] = False
        elif yi < small_neighbor:
            if not broad_dip:
                keep[i] = False
    return keep

def style_axes(axes):
    for ax in axes:
        ax.tick_params(which="both", length=8, direction="in", right=True, top=True, labelsize=14)
        ax.minorticks_on()


def set_xlimits(axes, tmin, *x_arrays):
    valid_arrays = []
    for x in x_arrays:
        arr = np.asarray(x)
        if arr.size == 0:
            continue
        arr = arr[np.isfinite(arr)]
        if arr.size == 0:
            continue
        valid_arrays.append(arr)
    if not valid_arrays:
        return
    xmin = float(tmin) if tmin is not None else min(np.min(x) for x in valid_arrays)
    xmax = max(np.max(x) for x in valid_arrays)
    if xmax <= xmin:
        xmax = xmin + 1.0
    for ax in axes:
        ax.set_xlim(xmin, xmax)


def add_collapse_marker(ax, tcoll_plot, label=r"${\rm collapse}$",
                        color="gray", alpha=0.65, linewidth=1.2,
                        y_axes: float = 0.98, va: str = "top"):
    ax.axvline(
        tcoll_plot,
        color=color,
        linestyle="-",
        linewidth=linewidth,
        alpha=alpha,
        zorder=5,
    )
    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    ax.text(
        tcoll_plot,
        y_axes,
        label,
        transform=trans,
        color=color,
        alpha=min(alpha + 0.10, 0.95),
        rotation=90,
        ha="left",
        va=va,
        fontsize=13,
        clip_on=False,
    )


def maybe_apply_collapse(ax, tcoll_plot, label=r"${\rm collapse}$",
                         color="gray", alpha=0.65, linewidth=1.2,
                         y_axes: float = 0.98, va: str = "top"):
    if tcoll_plot is not None:
        add_collapse_marker(
            ax,
            tcoll_plot,
            label=label,
            color=color,
            alpha=alpha,
            linewidth=linewidth,
            y_axes=y_axes,
            va=va,
        )


def compute_tcoll_plot(tcoll_ms_arg: Optional[float], tcoll_code_arg: Optional[float],
                       tmerg_ms: Optional[float]) -> Optional[float]:
    if tcoll_ms_arg is not None:
        return tcoll_ms_arg
    if tcoll_code_arg is not None:
        tcoll_ms = tcoll_code_arg * CU_to_ms
        return tcoll_ms if tmerg_ms is None else tcoll_ms - tmerg_ms
    return None

# ============================================================
# Reading ASCII files
# ============================================================

EMPTY_COLUMNS = ["iteration", "time_code", "time_ms", "value", "source_file", "output_dir"]


def read_carpet_ascii_file(file_path: Path) -> pd.DataFrame:
    rows = []
    try:
        if (not file_path.exists()) or file_path.stat().st_size == 0:
            return pd.DataFrame(columns=EMPTY_COLUMNS)
    except OSError:
        return pd.DataFrame(columns=EMPTY_COLUMNS)

    is_scalar_file = "data_Scalar" in str(file_path)
    is_0d_file = "data_asc_0D" in str(file_path)

    with open(file_path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            try:
                vals = [float(x) for x in parts]
            except ValueError:
                continue
            if len(vals) < 3:
                continue
            try:
                iteration = int(vals[0])
            except (ValueError, OverflowError):
                continue

            if is_scalar_file:
                time_code = vals[1]
                value = vals[-1]
            elif is_0d_file and len(vals) >= 9:
                time_code = vals[8]
                value = vals[-1]
            elif len(vals) >= 9:
                time_code = vals[8]
                value = vals[-1]
            else:
                time_code = vals[1]
                value = vals[-1]

            if not np.isfinite(time_code) or not np.isfinite(value):
                continue
            rows.append({
                "iteration": iteration,
                "time_code": time_code,
                "value": value,
                "source_file": str(file_path),
                "output_dir": file_path.parent.parent.name,
            })

    if not rows:
        return pd.DataFrame(columns=EMPTY_COLUMNS)

    df = pd.DataFrame(rows)
    df["time_ms"] = df["time_code"] * CU_to_ms
    return df


def read_em_luminosity_file(file_path: Path) -> pd.DataFrame:
    arr = np.loadtxt(file_path, comments="#")
    if arr.ndim == 1:
        arr = arr[None, :]
    df = pd.DataFrame({
        "iteration": arr[:, 0].astype(int),
        "time_code": arr[:, 1],
        "time_ms": arr[:, 2],
        "radius_code": arr[:, 3],
        "mean_Sr_code": arr[:, 4],
        "L_EM_code": arr[:, 5],
        "mean_Sr_cgs": arr[:, 6],
        "L_EM_cgs": arr[:, 7],
    })
    df = df.sort_values(["iteration", "time_code"]).reset_index(drop=True)
    df = df.drop_duplicates(subset=["iteration"], keep="last").reset_index(drop=True)
    return df


def candidate_output_dirs(root: Path) -> List[Path]:
    return sorted([p for p in root.glob("output*") if p.is_dir()], key=natural_key)


def combine_series_across_outputs(root: Path, relfile_or_aliases: Sequence[str] | str) -> pd.DataFrame:
    aliases = [relfile_or_aliases] if isinstance(relfile_or_aliases, str) else list(relfile_or_aliases)
    files = []
    for outdir in candidate_output_dirs(root):
        picked = None
        for relfile in aliases:
            candidate = outdir / relfile
            if candidate.exists():
                picked = candidate
                break
        if picked is not None:
            files.append(picked)

    if not files:
        raise FileNotFoundError(f"No files found for any of: {aliases}")

    frames = []
    for f in files:
        df = read_carpet_ascii_file(f)
        if not df.empty:
            frames.append(df)

    if not frames:
        print(f"[WARNING] All matching files were empty/unreadable: {aliases} -> skipping.")
        return pd.DataFrame(columns=EMPTY_COLUMNS)

    df = pd.concat(frames, ignore_index=True)
    df = ensure_time_ms(df)
    df = df.sort_values(["iteration", "time_code", "source_file"]).reset_index(drop=True)
    df = df.drop_duplicates(subset=["iteration"], keep="last").reset_index(drop=True)
    return df


def load_optional_series(root: Path, aliases: Sequence[str] | str) -> pd.DataFrame:
    try:
        return combine_series_across_outputs(root, aliases)
    except FileNotFoundError:
        return pd.DataFrame(columns=EMPTY_COLUMNS)

# ============================================================
# Scalar definitions
# ============================================================

STAT_ORDER = ["minimum", "average", "maximum"]
STAT_SUFFIX = {"minimum": "min", "average": "ave", "maximum": "max"}


def stat_text(base_tex: str, stat: str) -> str:
    base_grouped = r"{" + base_tex + r"}"
    if stat == "minimum":
        return base_grouped + r"_{\min}"
    if stat == "maximum":
        return base_grouped + r"_{\max}"
    return r"\langle " + base_tex + r" \rangle"


def build_scalar_specs(include_m1: bool) -> List[Dict[str, object]]:
    quantities = [
        {"name": "rho_b", "display": r"\rho_b", "aliases": {s: ["rho_b"] for s in STAT_ORDER}, "factor": CU_to_densCGS, "yscale": "log", "ylim": None},
        {"name": "temp", "display": r"T", "aliases": {s: ["temp", "temp_b"] for s in STAT_ORDER}, "factor": 1.0, "yscale": "log", "ylim": None},
        {"name": "ye", "display": r"Y_e", "aliases": {s: ["ye"] for s in STAT_ORDER}, "factor": 1.0, "yscale": "linear", "ylim": (0.0, 1.0)},
        {"name": "alp", "display": r"\alpha", "aliases": {s: ["alp"] for s in STAT_ORDER}, "factor": 1.0, "yscale": "linear", "ylim": (0.0, 1.05)},
        {"name": "smallb2", "display": r"b^2", "aliases": {s: ["smallb2"] for s in STAT_ORDER}, "factor": 1.0, "yscale": "log", "ylim": None},
        {"name": "Bvec0", "display": r"B_x", "aliases": {s: ["Bvec[0]", "Bx"] for s in STAT_ORDER}, "factor": 1.0, "yscale": "linear", "ylim": None},
        {"name": "Bvec1", "display": r"B_y", "aliases": {s: ["Bvec[1]", "By"] for s in STAT_ORDER}, "factor": 1.0, "yscale": "linear", "ylim": None},
        {"name": "Bvec2", "display": r"B_z", "aliases": {s: ["Bvec[2]", "Bz"] for s in STAT_ORDER}, "factor": 1.0, "yscale": "linear", "ylim": None},
        {"name": "P", "display": r"P", "aliases": {s: ["P", "press"] for s in STAT_ORDER}, "factor": CU_to_pressCGS, "yscale": "log", "ylim": None},
    ]
    if include_m1:
        quantities.extend([
            {"name": "Enue", "display": r"E_{\nu_e}", "aliases": {s: ["Enue"] for s in STAT_ORDER}, "factor": 1.0, "yscale": "log", "ylim": None},
            {"name": "Enue_bar", "display": r"E_{\bar{\nu}_e}", "aliases": {s: ["Enue_bar"] for s in STAT_ORDER}, "factor": 1.0, "yscale": "log", "ylim": None},
            {"name": "Enux", "display": r"E_{\nu_x}", "aliases": {s: ["Enux"] for s in STAT_ORDER}, "factor": 1.0, "yscale": "log", "ylim": None},
        ])
    specs = []
    for q in quantities:
        for stat in STAT_ORDER:
            aliases = [f"data_Scalar/{root}.{stat}.asc" for root in q["aliases"][stat]]
            specs.append({
                "key": f"{q['name']}.{stat}",
                "outstem": f"{q['name']}_{STAT_SUFFIX[stat]}",
                "aliases": aliases,
                "ylabel": r"$" + stat_text(q["display"], stat) + r"$",
                "factor": q["factor"],
                "yscale": q["yscale"],
                "ylim": q["ylim"],
                "quantity_name": q["name"],
                "stat": stat,
            })
    return specs


def load_scalar_collection(root: Path, include_m1: bool) -> Dict[str, Dict[str, object]]:
    out: Dict[str, Dict[str, object]] = {}
    for spec in build_scalar_specs(include_m1):
        df = load_optional_series(root, spec["aliases"])
        if df.empty:
            continue
        out[spec["key"]] = {**spec, "df": df}
    return out

# ============================================================
# Curve prep and plotting
# ============================================================

def get_x(df: pd.DataFrame, tmerg_ms: Optional[float]) -> np.ndarray:
    df = ensure_time_ms(df)
    if df.empty:
        return np.array([])
    if tmerg_ms is None:
        return df["time_ms"].to_numpy()
    return (df["time_ms"] - tmerg_ms).to_numpy()


def prepare_scalar_curve(entry: Dict[str, object], tmerg_ms: Optional[float]) -> Tuple[np.ndarray, np.ndarray]:
    df = entry["df"]
    if df is None or df.empty or "value" not in df.columns:
        return np.array([]), np.array([])
    x = get_x(df, tmerg_ms)
    y = df["value"].to_numpy() * float(entry["factor"])
    return x, y


def finalize_axes(ax, xlabel: str, ylabel: str, x_arrays: Sequence[np.ndarray],
                  tmin: Optional[float], ylim: Optional[Tuple[float, float]],
                  tcoll_m1_plot: Optional[float] = None,
                  tcoll_nom1_plot: Optional[float] = None,
                  collapse_y_axes: float = 0.98,
                  collapse_va: str = "top"):
    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    if ylim is not None:
        ax.set_ylim(ylim)
    style_axes([ax])
    set_xlimits([ax], tmin, *x_arrays)
    maybe_apply_collapse(ax, tcoll_m1_plot, label=r"collapse-M1", color="dimgray", alpha=0.75, linewidth=1.2, y_axes=collapse_y_axes, va=collapse_va)
    maybe_apply_collapse(ax, tcoll_nom1_plot, label=r"collapse-no-$\nu$'s", color="gray", alpha=0.60, linewidth=1.0, y_axes=collapse_y_axes, va=collapse_va)


def plot_scalar_single(entry: Dict[str, object], tmerg_ms: Optional[float], outdir: Path,
                       prefix: str, xlabel: str, tmin: Optional[float], tcoll_plot: Optional[float],
                       colour_m1: str = "black") -> bool:
    x, y = prepare_scalar_curve(entry, tmerg_ms)
    if x.size == 0 or y.size == 0:
        print(f"[WARNING] No usable data for {entry['key']} -> skipping plot.")
        return False
    fig = plt.figure(figsize=FIGSIZE_SQUARE)
    ax = fig.add_subplot(111)
    if entry["yscale"] == "log":
        mask = np.isfinite(y) & (y > 0.0)
        ax.semilogy(x[mask], y[mask], color=colour_m1, linewidth=2)
    else:
        mask = np.isfinite(y)
        ax.plot(x[mask], y[mask], color=colour_m1, linewidth=2)
    if not np.any(mask):
        plt.close(fig)
        return False
    finalize_axes(ax, xlabel, entry["ylabel"], [x[mask]], tmin, entry["ylim"], tcoll_plot, None)
    fig.tight_layout()
    fname = outdir / f"{prefix}_{entry['outstem']}.pdf"
    fig.savefig(fname, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {fname}")
    return True


def plot_scalar_compare(entry_m1: Dict[str, object], entry_nom1: Dict[str, object],
                        tmerg_m1_ms: Optional[float], tmerg_nom1_ms: Optional[float], outdir: Path,
                        prefix: str, xlabel: str, tmin: Optional[float], tcoll_m1_plot: Optional[float],
                        tcoll_nom1_plot: Optional[float], colour_m1: str = "black",
                        colour_nom1: str = "gray", collapse_y_axes: float = 0.98,
                        collapse_va: str = "top") -> bool:
    x1, y1 = prepare_scalar_curve(entry_m1, tmerg_m1_ms)
    x0, y0 = prepare_scalar_curve(entry_nom1, tmerg_nom1_ms)
    if x1.size == 0 or y1.size == 0 or x0.size == 0 or y0.size == 0:
        return False
    fig = plt.figure(figsize=FIGSIZE_SQUARE)
    ax = fig.add_subplot(111)
    use_log = entry_m1["yscale"] == "log"
    if use_log:
        m1 = remove_isolated_log_spikes(x1, y1)
        m0 = remove_isolated_log_spikes(x0, y0)
        ax.semilogy(x1[m1], y1[m1], color=colour_m1, linewidth=2.0, label=r"${\rm M1}$")
        ax.semilogy(x0[m0], y0[m0], color=colour_nom1, linewidth=1.6, alpha=0.85, label=r"${\rm no~M1}$")
    else:
        m1 = np.isfinite(y1)
        m0 = np.isfinite(y0)
        ax.plot(x1[m1], y1[m1], color=colour_m1, linewidth=2.0, label=r"${\rm M1}$")
        ax.plot(x0[m0], y0[m0], color=colour_nom1, linewidth=1.6, alpha=0.85, label=r"${\rm no~M1}$")
    if not np.any(m1) and not np.any(m0):
        plt.close(fig)
        return False
    ax.legend(fontsize=13)
    finalize_axes(ax, xlabel, entry_m1["ylabel"], [x1[m1] if np.any(m1) else np.array([]), x0[m0] if np.any(m0) else np.array([])], tmin, entry_m1["ylim"], tcoll_m1_plot, tcoll_nom1_plot, collapse_y_axes=collapse_y_axes, collapse_va=collapse_va)
    fig.tight_layout()
    fname = outdir / f"{prefix}_{entry_m1['outstem']}_compare_M1_noM1.pdf"
    fig.savefig(fname, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {fname}")
    return True


def plot_series_single(df: pd.DataFrame, ycol: str, ylabel: str,
                       tmerg_ms: Optional[float], outdir: Path, prefix: str, suffix: str,
                       xlabel: str, tmin: Optional[float], tcoll_plot: Optional[float],
                       factor: float = 1.0, logscale: bool = True,
                       ylim: Optional[Tuple[float, float]] = None, colour_m1: str = "black",
                       collapse_y_axes: float = 0.98, collapse_va: str = "top") -> bool:
    if df.empty or ycol not in df.columns:
        return False
    x = get_x(df, tmerg_ms)
    y = df[ycol].to_numpy() * factor
    if logscale:
        m = remove_isolated_log_spikes(x, y)
    else:
        m = np.isfinite(x) & np.isfinite(y)
    if not np.any(m):
        return False
    fig = plt.figure(figsize=FIGSIZE_SQUARE)
    ax = fig.add_subplot(111)
    if logscale:
        ax.semilogy(x[m], y[m], color=colour_m1, linewidth=2.0)
    else:
        ax.plot(x[m], y[m], color=colour_m1, linewidth=2.0)
    finalize_axes(ax, xlabel, ylabel, [x[m]], tmin, ylim, tcoll_plot, None, collapse_y_axes=collapse_y_axes, collapse_va=collapse_va)
    fig.tight_layout()
    fname = outdir / f"{prefix}_{suffix}.pdf"
    fig.savefig(fname, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {fname}")
    return True


def plot_series_compare(df_m1: pd.DataFrame, df_nom1: pd.DataFrame,
                        ycol_m1: str, ycol_nom1: str, ylabel: str,
                        tmerg_m1_ms: Optional[float], tmerg_nom1_ms: Optional[float],
                        outdir: Path, prefix: str, suffix: str, xlabel: str,
                        tmin: Optional[float], tcoll_m1_plot: Optional[float],
                        tcoll_nom1_plot: Optional[float], factor_m1: float = 1.0,
                        factor_nom1: float = 1.0, logscale: bool = True,
                        ylim: Optional[Tuple[float, float]] = None,
                        label_m1: str = r"${\rm M1}$",
                        label_nom1: str = r"${\rm no~M1}$",
                        colour_m1: str = "black", colour_nom1: str = "gray",
                        collapse_y_axes: float = 0.98, collapse_va: str = "top") -> bool:
    x1 = get_x(df_m1, tmerg_m1_ms) if (not df_m1.empty and ycol_m1 in df_m1.columns) else np.array([])
    x0 = get_x(df_nom1, tmerg_nom1_ms) if (not df_nom1.empty and ycol_nom1 in df_nom1.columns) else np.array([])
    y1 = df_m1[ycol_m1].to_numpy() * factor_m1 if (not df_m1.empty and ycol_m1 in df_m1.columns) else np.array([])
    y0 = df_nom1[ycol_nom1].to_numpy() * factor_nom1 if (not df_nom1.empty and ycol_nom1 in df_nom1.columns) else np.array([])
    if logscale:
        m1 = remove_isolated_log_spikes(x1, y1)
        m0 = remove_isolated_log_spikes(x0, y0)
    else:
        m1 = np.isfinite(x1) & np.isfinite(y1)
        m0 = np.isfinite(x0) & np.isfinite(y0)
    if not np.any(m1) and not np.any(m0):
        return False
    fig = plt.figure(figsize=FIGSIZE_SQUARE)
    ax = fig.add_subplot(111)
    if np.any(m1):
        if logscale:
            ax.semilogy(x1[m1], y1[m1], color=colour_m1, linewidth=2.0, label=label_m1)
        else:
            ax.plot(x1[m1], y1[m1], color=colour_m1, linewidth=2.0, label=label_m1)
    if np.any(m0):
        if logscale:
            ax.semilogy(x0[m0], y0[m0], color=colour_nom1, linewidth=1.6, alpha=0.85, label=label_nom1)
        else:
            ax.plot(x0[m0], y0[m0], color=colour_nom1, linewidth=1.6, alpha=0.85, label=label_nom1)
    ax.legend(fontsize=13)
    finalize_axes(ax, xlabel, ylabel, [x1[m1] if np.any(m1) else np.array([]), x0[m0] if np.any(m0) else np.array([])], tmin, ylim, tcoll_m1_plot, tcoll_nom1_plot, collapse_y_axes=collapse_y_axes, collapse_va=collapse_va)
    fig.tight_layout()
    fname = outdir / f"{prefix}_{suffix}.pdf"
    fig.savefig(fname, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {fname}")
    return True


def plot_neutrino_luminosities_m1(df_lnue: pd.DataFrame, df_lnua: pd.DataFrame, df_lnux: pd.DataFrame,
                                  tmerg_m1_ms: Optional[float], outdir: Path, prefix: str,
                                  xlabel: str, tmin: Optional[float], tcoll_m1_plot: Optional[float],
                                  colour_m1: str = "black", collapse_y_axes: float = 0.15,
                                  collapse_va: str = "bottom") -> bool:
    if df_lnue.empty or df_lnua.empty or df_lnux.empty:
        return False
    x_nue = get_x(df_lnue, tmerg_m1_ms)
    x_nua = get_x(df_lnua, tmerg_m1_ms)
    x_nux = get_x(df_lnux, tmerg_m1_ms)
    conv = CU_to_energyCGS / CU_to_s
    with np.errstate(over='ignore', invalid='ignore'):
        Lnue = df_lnue['value'].to_numpy(dtype=float) * conv
        Lnua = df_lnua['value'].to_numpy(dtype=float) * conv
        Lnux = 4.0 * df_lnux['value'].to_numpy(dtype=float) * conv
        Lnu = Lnue + Lnua + Lnux
    m_tot = remove_isolated_log_spikes(x_nue, Lnu)
    m_nue = remove_isolated_log_spikes(x_nue, Lnue)
    m_nua = remove_isolated_log_spikes(x_nua, Lnua)
    m_nux = remove_isolated_log_spikes(x_nux, Lnux)
    if not (np.any(m_tot) or np.any(m_nue) or np.any(m_nua) or np.any(m_nux)):
        return False
    fig = plt.figure(figsize=FIGSIZE_SQUARE)
    ax = fig.add_subplot(111)
    if np.any(m_tot):
        ax.semilogy(x_nue[m_tot], Lnu[m_tot], color=colour_m1, linewidth=2.0, label=r'${\rm total}$')
    if np.any(m_nue):
        ax.semilogy(x_nue[m_nue], Lnue[m_nue], color='red', linestyle='--', linewidth=1.8, label=r'$\nu_e$')
    if np.any(m_nua):
        ax.semilogy(x_nua[m_nua], Lnua[m_nua], color='blue', linestyle='--', linewidth=1.8, label=r'$\bar{\nu}_e$')
    if np.any(m_nux):
        ax.semilogy(x_nux[m_nux], Lnux[m_nux], color='green', linestyle=':', linewidth=1.8, label=r'$\nu_x$')
    ax.legend(fontsize=13, loc='lower right')
    finalize_axes(ax, xlabel, r'$L_\nu~[{\rm erg/s}]$', [x_nue[m_tot] if np.any(m_tot) else np.array([])], tmin, (1e51, 2e54), tcoll_m1_plot, None, collapse_y_axes=collapse_y_axes, collapse_va=collapse_va)
    fig.tight_layout()
    fname = outdir / f"{prefix}_M1only_L_nu.pdf"
    fig.savefig(fname, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {fname}")
    return True

# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=str, default=".")
    parser.add_argument("--root-M1", type=str, default=None)
    parser.add_argument("--root-noM1", type=str, default=None)
    parser.add_argument("--compare-M1-noM1", action="store_true")
    parser.add_argument("--tmerg-code", type=float, default=None)
    parser.add_argument("--tmerg-code-M1", type=float, default=None)
    parser.add_argument("--tmerg-code-noM1", type=float, default=None)
    parser.add_argument("--tbreakout-code", type=float, default=None)
    parser.add_argument("--trec-code", type=float, default=None)
    parser.add_argument("--tmin", type=float, default=0.0)
    parser.add_argument("--tcoll", type=float, default=None)
    parser.add_argument("--tcoll-code", type=float, default=None)
    parser.add_argument("--tcoll-M1", type=float, default=None)
    parser.add_argument("--tcoll-code-M1", type=float, default=None)
    parser.add_argument("--tcoll-noM1", type=float, default=None)
    parser.add_argument("--tcoll-code-noM1", type=float, default=None)
    parser.add_argument("--em-file", type=str, default="em_luminosity_det3.dat")
    parser.add_argument("--em-file-M1", type=str, default=None)
    parser.add_argument("--em-file-noM1", type=str, default=None)
    parser.add_argument("--out", type=str, default="outflows_joined_from_ascii.pdf")
    parser.add_argument("--colour-M1", type=str, default="black")
    parser.add_argument("--colour-noM1", type=str, default="gray")
    parser.add_argument("--export-csv", action="store_true")
    parser.add_argument("--noM1", action="store_true")
    args = parser.parse_args()

    outpath = Path(args.out)
    outbase = outpath.with_suffix("")
    outdir = outbase.parent.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = outbase.name

    if not args.compare_M1_noM1:
        root = Path(args.root).resolve()
        em_file = Path(args.em_file).resolve()
        tmerg_ms = args.tmerg_code * CU_to_ms if args.tmerg_code is not None else None
        xlabel = r"$t~[{\rm ms}]$" if tmerg_ms is None else r"$t-t_{\rm mer}~[{\rm ms}]$"
        tcoll_plot = compute_tcoll_plot(args.tcoll, args.tcoll_code, tmerg_ms)

        df_bern = load_optional_series(root, "data_asc_0D/bernoulli_mass_ubound[2]..asc")
        df_emag = load_optional_series(root, "data_asc_0D/magnetic_Emag..asc")
        df_epol = load_optional_series(root, "data_asc_0D/magnetic_EmagPOL..asc")
        df_etor = load_optional_series(root, "data_asc_0D/magnetic_EmagTOR..asc")
        df_lem = read_em_luminosity_file(em_file) if em_file.exists() else pd.DataFrame()
        plot_series_single(df_emag, "value", r"$E_{\rm mag}~[{\rm erg}]$", tmerg_ms, outdir, prefix, "E_mag", xlabel, args.tmin, tcoll_plot, CU_to_energyCGS, True, None, colour_m1=args.colour_M1)
        plot_series_single(df_epol, "value", r"$E_{\rm mag,pol}~[{\rm erg}]$", tmerg_ms, outdir, prefix, "E_magPOL", xlabel, args.tmin, tcoll_plot, CU_to_energyCGS, True, None, colour_m1=args.colour_M1)
        plot_series_single(df_etor, "value", r"$E_{\rm mag,tor}~[{\rm erg}]$", tmerg_ms, outdir, prefix, "E_magTOR", xlabel, args.tmin, tcoll_plot, CU_to_energyCGS, True, None, colour_m1=args.colour_M1)
        plot_series_single(df_lem, "L_EM_cgs", r"$L_{\rm EM}~[{\rm erg/s}]$", tmerg_ms, outdir, prefix, "L_EM", xlabel, args.tmin, tcoll_plot, 1.0, True, (1e30, 1e48), colour_m1=args.colour_M1, collapse_y_axes=0.1, collapse_va="bottom")
        plot_series_single(df_bern, "value", r"$\dot{M}_{\rm ej}~[M_{\odot}/{\rm s}]$", tmerg_ms, outdir, prefix, "M_ej", xlabel, args.tmin, tcoll_plot, 1.0 / CU_to_s, True, (1e-4, 1e0), colour_m1=args.colour_M1, collapse_y_axes=0.1, collapse_va="bottom")
        if not args.noM1:
            df_lnue = load_optional_series(root, "data_asc_0D/L_nue[2]..asc")
            df_lnua = load_optional_series(root, "data_asc_0D/L_nue_bar[2]..asc")
            df_lnux = load_optional_series(root, "data_asc_0D/L_nux[2]..asc")
            plot_neutrino_luminosities_m1(df_lnue, df_lnua, df_lnux, tmerg_ms, outdir, prefix, xlabel, args.tmin, tcoll_plot, colour_m1=args.colour_M1, collapse_y_axes=0.15, collapse_va="bottom")
        scalar_data = load_scalar_collection(root, include_m1=not args.noM1)
        for key in sorted(scalar_data.keys()):
            plot_scalar_single(scalar_data[key], tmerg_ms, outdir, prefix, xlabel, args.tmin, tcoll_plot, colour_m1=args.colour_M1)
        if args.export_csv:
            exports = {
                "magnetic_Emag": df_emag,
                "magnetic_EmagPOL": df_epol,
                "magnetic_EmagTOR": df_etor,
                "L_EM": df_lem,
                "M_ej": df_bern,
            }
            if not args.noM1:
                exports.update({"L_nue": df_lnue, "L_nue_bar": df_lnua, "L_nux": df_lnux})
            for key, entry in scalar_data.items():
                exports[key] = entry["df"]
            for name, df in exports.items():
                if df is None or df.empty:
                    continue
                csvfile = outdir / f"{sanitize_stem(name)}.csv"
                df.to_csv(csvfile, index=False)
                print(f"Saved {csvfile}")
        return

    # compare mode
    root_m1 = Path(args.root_M1).resolve()
    root_nom1 = Path(args.root_noM1).resolve()
    tmerg_m1_ms = (args.tmerg_code_M1 if args.tmerg_code_M1 is not None else args.tmerg_code)
    tmerg_nom1_ms = (args.tmerg_code_noM1 if args.tmerg_code_noM1 is not None else args.tmerg_code)
    tmerg_m1_ms = None if tmerg_m1_ms is None else tmerg_m1_ms * CU_to_ms
    tmerg_nom1_ms = None if tmerg_nom1_ms is None else tmerg_nom1_ms * CU_to_ms
    xlabel = r"$t~[{\rm ms}]$" if (tmerg_m1_ms is None and tmerg_nom1_ms is None) else r"$t-t_{\rm mer}~[{\rm ms}]$"
    tcoll_m1_plot = compute_tcoll_plot(args.tcoll_M1, args.tcoll_code_M1, tmerg_m1_ms)
    tcoll_nom1_plot = compute_tcoll_plot(args.tcoll_noM1, args.tcoll_code_noM1, tmerg_nom1_ms)

    created_any = False

    df_bern_m1 = load_optional_series(root_m1, "data_asc_0D/bernoulli_mass_ubound[2]..asc")
    df_bern_nom1 = load_optional_series(root_nom1, "data_asc_0D/bernoulli_mass_ubound[2]..asc")
    created_any |= plot_series_compare(df_bern_m1, df_bern_nom1, "value", "value", r"$\dot{M}_{\rm ej}~[M_{\odot}/{\rm s}]$", tmerg_m1_ms, tmerg_nom1_ms, outdir, prefix, "M_ej_compare_M1_noM1", xlabel, args.tmin, tcoll_m1_plot, tcoll_nom1_plot, 1.0/CU_to_s, 1.0/CU_to_s, True, (1e-4,1e0), colour_m1=args.colour_M1, colour_nom1=args.colour_noM1, collapse_y_axes=0.1, collapse_va="bottom")

    df_emag_m1 = load_optional_series(root_m1, "data_asc_0D/magnetic_Emag..asc")
    df_emag_nom1 = load_optional_series(root_nom1, "data_asc_0D/magnetic_Emag..asc")
    created_any |= plot_series_compare(df_emag_m1, df_emag_nom1, "value", "value", r"$E_{\rm mag}~[{\rm erg}]$", tmerg_m1_ms, tmerg_nom1_ms, outdir, prefix, "E_mag_compare_M1_noM1", xlabel, args.tmin, tcoll_m1_plot, tcoll_nom1_plot, CU_to_energyCGS, CU_to_energyCGS, True, None, colour_m1=args.colour_M1, colour_nom1=args.colour_noM1)

    df_epol_m1 = load_optional_series(root_m1, "data_asc_0D/magnetic_EmagPOL..asc")
    df_epol_nom1 = load_optional_series(root_nom1, "data_asc_0D/magnetic_EmagPOL..asc")
    created_any |= plot_series_compare(df_epol_m1, df_epol_nom1, "value", "value", r"$E_{\rm mag,pol}~[{\rm erg}]$", tmerg_m1_ms, tmerg_nom1_ms, outdir, prefix, "E_magPOL_compare_M1_noM1", xlabel, args.tmin, tcoll_m1_plot, tcoll_nom1_plot, CU_to_energyCGS, CU_to_energyCGS, True, None, colour_m1=args.colour_M1, colour_nom1=args.colour_noM1)

    df_etor_m1 = load_optional_series(root_m1, "data_asc_0D/magnetic_EmagTOR..asc")
    df_etor_nom1 = load_optional_series(root_nom1, "data_asc_0D/magnetic_EmagTOR..asc")
    created_any |= plot_series_compare(df_etor_m1, df_etor_nom1, "value", "value", r"$E_{\rm mag,tor}~[{\rm erg}]$", tmerg_m1_ms, tmerg_nom1_ms, outdir, prefix, "E_magTOR_compare_M1_noM1", xlabel, args.tmin, tcoll_m1_plot, tcoll_nom1_plot, CU_to_energyCGS, CU_to_energyCGS, True, None, colour_m1=args.colour_M1, colour_nom1=args.colour_noM1)

    df_lem_m1 = read_em_luminosity_file(Path(args.em_file_M1).resolve()) if args.em_file_M1 else pd.DataFrame()
    df_lem_nom1 = read_em_luminosity_file(Path(args.em_file_noM1).resolve()) if args.em_file_noM1 else pd.DataFrame()
    created_any |= plot_series_compare(df_lem_m1, df_lem_nom1, "L_EM_cgs", "L_EM_cgs", r"$L_{\rm EM}~[{\rm erg/s}]$", tmerg_m1_ms, tmerg_nom1_ms, outdir, prefix, "L_EM_compare_M1_noM1", xlabel, args.tmin, tcoll_m1_plot, tcoll_nom1_plot, 1.0, 1.0, True, (1e30, 1e48), "M1", "no M1", colour_m1=args.colour_M1, colour_nom1=args.colour_noM1, collapse_y_axes=0.1, collapse_va="bottom")

    df_lnue_m1 = load_optional_series(root_m1, "data_asc_0D/L_nue[2]..asc")
    df_lnua_m1 = load_optional_series(root_m1, "data_asc_0D/L_nue_bar[2]..asc")
    df_lnux_m1 = load_optional_series(root_m1, "data_asc_0D/L_nux[2]..asc")
    created_any |= plot_neutrino_luminosities_m1(df_lnue_m1, df_lnua_m1, df_lnux_m1, tmerg_m1_ms, outdir, prefix, xlabel, args.tmin, tcoll_m1_plot, colour_m1=args.colour_M1, collapse_y_axes=0.15, collapse_va="bottom")

    scalar_m1 = load_scalar_collection(root_m1, include_m1=True)
    scalar_nom1 = load_scalar_collection(root_nom1, include_m1=False)
    common_keys = sorted(set(scalar_m1.keys()) & set(scalar_nom1.keys()))
    if not common_keys:
        print("[WARNING] No common scalar files were found between the M1 and no-M1 runs.")
    for key in common_keys:
        created_any |= plot_scalar_compare(scalar_m1[key], scalar_nom1[key], tmerg_m1_ms, tmerg_nom1_ms, outdir, prefix, xlabel, args.tmin, tcoll_m1_plot, tcoll_nom1_plot, colour_m1=args.colour_M1, colour_nom1=args.colour_noM1)
    m1_only_keys = sorted(set(scalar_m1.keys()) - set(scalar_nom1.keys()))
    for key in m1_only_keys:
        if scalar_m1[key]["quantity_name"].startswith("Enu"):
            created_any |= plot_scalar_single(scalar_m1[key], tmerg_m1_ms, outdir, f"{prefix}_M1only", xlabel, args.tmin, tcoll_m1_plot, colour_m1=args.colour_M1)

    if args.export_csv:
        compare_exports = {
            "M1_M_ej": df_bern_m1,
            "noM1_M_ej": df_bern_nom1,
            "M1_E_mag": df_emag_m1,
            "noM1_E_mag": df_emag_nom1,
            "M1_E_magPOL": df_epol_m1,
            "noM1_E_magPOL": df_epol_nom1,
            "M1_E_magTOR": df_etor_m1,
            "noM1_E_magTOR": df_etor_nom1,
            "M1_L_EM": df_lem_m1,
            "noM1_L_EM": df_lem_nom1,
            "M1_L_nue": df_lnue_m1,
            "M1_L_nue_bar": df_lnua_m1,
            "M1_L_nux": df_lnux_m1,
        }
        for _, entry in scalar_m1.items():
            compare_exports[f"M1_{entry['outstem']}"] = entry["df"]
        for _, entry in scalar_nom1.items():
            compare_exports[f"noM1_{entry['outstem']}"] = entry["df"]
        for name, df in compare_exports.items():
            if df is None or df.empty:
                continue
            csvfile = outdir / f"{sanitize_stem(name)}.csv"
            df.to_csv(csvfile, index=False)
            print(f"Saved {csvfile}")

    if not created_any:
        raise RuntimeError("No scalar or diagnostic plots could be created.")


if __name__ == "__main__":
    main()
