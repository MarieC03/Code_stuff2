#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import glob
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import matplotlib.pyplot as plt


def _require_healpy():
    try:
        import healpy as hp
    except Exception as exc:
        raise SystemExit(
            "This script needs healpy to reconstruct the HEALPix angular directions. "
            "Install it, e.g. with 'pip install healpy', or run in an environment where healpy is available."
        ) from exc
    return hp


def find_dataset(group: h5py.Group, candidates: list[str]) -> np.ndarray | None:
    for name in candidates:
        if name in group:
            return np.asarray(group[name][...])
    return None


def build_unit_vectors(npix: int, nest: bool = False):
    hp = _require_healpy()
    nside = hp.npix2nside(npix)
    pix = np.arange(npix)
    nx, ny, nz = hp.pix2vec(nside, pix, nest=nest)
    return np.asarray(nx), np.asarray(ny), np.asarray(nz), nside


def safe_flat(a: np.ndarray) -> np.ndarray:
    return np.asarray(a, dtype=float).reshape(-1)


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    wsum = np.sum(weights)
    if wsum <= 0.0:
        return np.nan
    return float(np.sum(weights * values) / wsum)


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    if values.size == 0:
        return np.nan
    idx = np.argsort(values)
    v = values[idx]
    w = weights[idx]
    cdf = np.cumsum(w)
    if cdf[-1] <= 0.0:
        return np.nan
    cdf /= cdf[-1]
    return float(np.interp(q, cdf, v))


def find_restart_number(path: str) -> int:
    m = re.search(r"output[-_]?(\d+)", path)
    return int(m.group(1)) if m else -1


def discover_files(root: str, detector: int | None, pattern: str | None) -> list[str]:
    if pattern is None:
        if detector is None:
            raise SystemExit("Please provide either --h5, or (--root and --detector), or --pattern.")
        patterns = [
            str(Path(root) / "output*" / "data" / f"healpix_det_{detector}_surf.h5"),
            str(Path(root) / "output*" / f"healpix_det_{detector}_surf.h5"),
            str(Path(root) / "**" / f"healpix_det_{detector}_surf.h5"),
        ]
    else:
        patterns = [pattern]

    found: list[str] = []
    for pat in patterns:
        recursive = "**" in pat
        found.extend(glob.glob(pat, recursive=recursive))

    files = sorted(set(found), key=lambda x: (find_restart_number(x), x))
    if not files:
        pats = "\n  ".join(patterns)
        raise FileNotFoundError(f"No detector files found. Looked for:\n  {pats}")
    return files


def parse_args():
    ap = argparse.ArgumentParser(
        description=(
            "Compute radial-velocity time series from FIL/outflow_plusplus "
            "HEALPix detector files such as healpix_det_3_surf.h5."
        )
    )
    src = ap.add_mutually_exclusive_group(required=False)
    src.add_argument("--h5", help="Single explicit path to healpix_det_*_surf.h5")
    src.add_argument(
        "--pattern",
        help=(
            "Explicit glob pattern for detector files. Example: "
            "'/path/to/sim/output*/data/healpix_det_3_surf.h5'"
        ),
    )
    ap.add_argument(
        "--root",
        default=".",
        help=(
            "Root directory under which detector files are searched. "
            "Default patterns are <root>/output*/data/healpix_det_<det>_surf.h5, "
            "<root>/output*/healpix_det_<det>_surf.h5, and recursive fallback."
        ),
    )
    ap.add_argument(
        "--detector",
        type=int,
        default=None,
        help="Detector index N in healpix_det_N_surf.h5. Needed with --root unless --h5 or --pattern is used.",
    )
    ap.add_argument(
        "--flux",
        default="mass",
        choices=["mass", "bernoulli", "geodesic", "none"],
        help=(
            "Weight / selection field. 'mass' uses outflow_plusplus::mass_flux_total_gf, "
            "'bernoulli' uses Bernoulli-unbound flux, 'geodesic' uses geodesic-unbound flux, "
            "'none' falls back to rho-weighting."
        ),
    )
    ap.add_argument(
        "--min-rho",
        type=float,
        default=0.0,
        help="Optional density floor for selecting pixels (code units).",
    )
    ap.add_argument(
        "--max-rho",
        type=float,
        default=np.inf,
        help="Optional density ceiling for selecting pixels (code units).",
    )
    ap.add_argument(
        "--vr-min",
        type=float,
        default=0.0,
        help="Only keep pixels with v_r > vr-min. Default: outward-moving only.",
    )
    ap.add_argument(
        "--nest",
        action="store_true",
        help="Use NESTED HEALPix ordering. Default assumes RING ordering.",
    )
    ap.add_argument(
        "--cu-to-ms",
        type=float,
        default=4.925490947e-3,
        help="Conversion factor from FIL code time units (M_sun) to ms.",
    )
    ap.add_argument(
        "--tmerger",
        type=float,
        default=None,
        help="Merger time in code units. If given, plot time as (t - t_merger) in ms.",
    )
    ap.add_argument(
        "--out-prefix",
        type=str,
        default="radial_velocity",
        help="Prefix for output CSV/PDF/PNG files.",
    )
    return ap.parse_args()


def process_one_file(
    h5_path: Path,
    args,
    times_code_all: dict[float, tuple[float, float, float, float, int]],
) -> None:
    flux_candidates = {
        "mass": ["outflow_plusplus::mass_flux_total_gf"],
        "bernoulli": [
            "outflow_plusplus::bernoulli_mass_ubound_gf",
            "outflow_plusplus::bernoulli_flux_ubound_gf",
        ],
        "geodesic": [
            "outflow_plusplus::geodesic_mass_ubound_gf",
            "outflow_plusplus::geodesic_flux_ubound_gf",
        ],
        "none": [],
    }

    velx_candidates = ["hydrobase::vel[0]", "vel[0]", "velx"]
    vely_candidates = ["hydrobase::vel[1]", "vel[1]", "vely"]
    velz_candidates = ["hydrobase::vel[2]", "vel[2]", "velz"]
    rho_candidates = ["hydrobase::rho", "rho", "hydrobase::rho_b", "rho_b"]
    time_candidates = ["time", "cctk_time"]

    with h5py.File(h5_path, "r") as f:
        if "data" not in f:
            print(f"[skip] {h5_path}: missing /data group", file=sys.stderr)
            return

        data_group = f["data"]
        keys = sorted(data_group.keys(), key=lambda s: int(s) if str(s).isdigit() else str(s))
        nx = ny = nz = None
        inferred_nside = None

        for key in keys:
            g = data_group[key]

            vx = find_dataset(g, velx_candidates)
            vy = find_dataset(g, vely_candidates)
            vz = find_dataset(g, velz_candidates)
            rho = find_dataset(g, rho_candidates)
            tarr = find_dataset(g, time_candidates)

            if vx is None or vy is None or vz is None:
                print(f"[skip] {h5_path} iteration {key}: missing velocity datasets", file=sys.stderr)
                continue
            if rho is None:
                print(f"[skip] {h5_path} iteration {key}: missing density dataset", file=sys.stderr)
                continue

            vx = safe_flat(vx)
            vy = safe_flat(vy)
            vz = safe_flat(vz)
            rho = safe_flat(rho)

            if not (vx.size == vy.size == vz.size == rho.size):
                print(f"[skip] {h5_path} iteration {key}: inconsistent dataset sizes", file=sys.stderr)
                continue

            if nx is None:
                nx, ny, nz, inferred_nside = build_unit_vectors(vx.size, nest=args.nest)
                print(f"[info] {h5_path.name}: detected HEALPix nside={inferred_nside}, npix={vx.size}")

            vr = vx * nx + vy * ny + vz * nz

            flux = None
            if args.flux != "none":
                flux = find_dataset(g, flux_candidates[args.flux])
                if flux is not None:
                    flux = safe_flat(flux)
                    if flux.size != vr.size:
                        print(f"[warn] {h5_path} iteration {key}: flux size mismatch, ignoring flux", file=sys.stderr)
                        flux = None

            mask = np.isfinite(vr) & np.isfinite(rho)
            mask &= (rho >= args.min_rho) & (rho <= args.max_rho)
            mask &= (vr > args.vr_min)

            if flux is not None:
                mask &= np.isfinite(flux)
                mask &= (flux > 0.0)
                weights = flux
            else:
                weights = rho

            if not np.any(mask):
                continue

            vsel = vr[mask]
            wsel = weights[mask]

            if tarr is None:
                tcode = float(key)
            else:
                tflat = safe_flat(tarr)
                tcode = float(tflat[0]) if tflat.size > 0 else float(key)

            # Later files overwrite earlier ones at identical times, matching restart logic.
            times_code_all[tcode] = (
                weighted_mean(vsel, wsel),
                weighted_quantile(vsel, wsel, 0.5),
                weighted_quantile(vsel, wsel, 0.9),
                float(np.nanmax(vsel)),
                int(vsel.size),
            )


def main():
    args = parse_args()

    if args.h5 is not None:
        files = [str(Path(args.h5))]
    else:
        files = discover_files(args.root, args.detector, args.pattern)

    print("[info] Using detector file(s):")
    for f in files:
        print(f"  {f}")

    series_by_time: dict[float, tuple[float, float, float, float, int]] = {}
    for fpath in files:
        process_one_file(Path(fpath), args, series_by_time)

    if len(series_by_time) == 0:
        raise SystemExit("No valid time samples were found.")

    times_code = np.array(sorted(series_by_time.keys()), dtype=float)
    vr_mean = np.array([series_by_time[t][0] for t in times_code])
    vr_median = np.array([series_by_time[t][1] for t in times_code])
    vr_p90 = np.array([series_by_time[t][2] for t in times_code])
    vr_max = np.array([series_by_time[t][3] for t in times_code])
    npixels_used = np.array([series_by_time[t][4] for t in times_code], dtype=int)

    if args.tmerger is None:
        times_ms = times_code * args.cu_to_ms
        xlabel = r"$t\ [{\rm ms}]$"
    else:
        times_ms = (times_code - args.tmerger) * args.cu_to_ms
        xlabel = r"$t - t_{\rm mer}\ [{\rm ms}]$"

    out_prefix = Path(args.out_prefix)
    csv_path = out_prefix.with_suffix(".csv")
    png_path = out_prefix.with_suffix(".png")
    pdf_path = out_prefix.with_suffix(".pdf")

    header = "time_code,time_ms,vr_mean,vr_median,vr_p90,vr_max,npixels_used"
    data = np.column_stack([times_code, times_ms, vr_mean, vr_median, vr_p90, vr_max, npixels_used])
    np.savetxt(csv_path, data, delimiter=",", header=header, comments="")

    plt.rcParams.update({
        "text.usetex": True,
        "font.size": 14,
        "font.family": "Serif",
        "font.serif": ["Computer Modern Roman", "CMU Serif", "DejaVu Serif"],
    })

    fig, ax = plt.subplots(figsize=(7.0, 6.0))
    ax.plot(times_ms, vr_mean, lw=2.2, label=rf"{args.flux}-weighted $\langle v_r \rangle$")
    ax.plot(times_ms, vr_median, lw=1.8, label=r"weighted median $v_r$")
    ax.plot(times_ms, vr_p90, lw=1.6, label=r"weighted 90th pct.")
    ax.set_xlabel(xlabel, fontsize=20)
    ax.set_ylabel(r"$v_r$", fontsize=20)
    ax.minorticks_on()
    ax.tick_params(direction="in", which="both", top=True, right=True, length=8)
    ax.legend(frameon=True, loc="best", fontsize=14)
    fig.tight_layout()
    fig.savefig(png_path, dpi=200)
    fig.savefig(pdf_path)

    print(f"Wrote {csv_path}")
    print(f"Wrote {png_path}")
    print(f"Wrote {pdf_path}")


if __name__ == "__main__":
    main()
