#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_bfield_from_npz.py
=======================
Reads pre-saved bundle_it*.npz files (produced by fil_2d_bfield_plot.py) and
reproduces the identical 5-panel b-field figure for every bundle found.

No HDF5 data, no physics computation, no external modules needed beyond
numpy and matplotlib.

Bundle keys expected in each .npz file
---------------------------------------
  X, Y               -- coordinate meshgrids [km]
  plasma_beta_inv    -- beta^{-1}  (= b^2 / 2p)
  b_pol              -- poloidal b-field [code units]
  b_tor              -- toroidal b-field [code units]
  W                  -- Lorentz factor
  Ye                 -- electron fraction
  rho_log_cgs        -- log10(rho  [g/cm^3])
  time_ms            -- simulation time [ms]
  time_plot_ms       -- plotted time (= time_ms - t_merg if merger time was set)
  plane              -- plane string, e.g. b'xz'
  axis0_name         -- e.g. b'x'
  axis1_name         -- e.g. b'z'

Usage
-----
  # Plot every bundle under the default search path:
  python plot_bfield_from_npz.py --npzdir derived_bfield_slices/xz

  # Plot a single file:
  python plot_bfield_from_npz.py --npzdir . --pattern bundle_it00012800.npz

  # Override conversion factor, disable TeX, custom output dir:
  python plot_bfield_from_npz.py --npzdir derived_bfield_slices/xz \\
      --cu-to-gauss 8.35195e19 --no-tex --outdir my_plots
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter

try:
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    _HAS_INSET = True
except ImportError:
    _HAS_INSET = False

try:
    import seaborn as sns

    def _rocket_cmap():
        return sns.color_palette("rocket", as_cmap=True)

except ImportError:

    def _rocket_cmap():
        return plt.get_cmap("magma")


# ============================================================================
# Physical constants (same as fil_2d_bfield_plot.py)
# ============================================================================
_DEFAULT_CU_TO_GAUSS = 8.35195e19 / np.sqrt(4.0 * np.pi)  # ~2.3548e19 G


# ============================================================================
# Helpers
# ============================================================================


def _natural_key(p):
    s = str(p)
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)]


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def safe_log10(arr, floor=1e-30):
    a = np.where(np.isfinite(arr), np.abs(arr), np.nan)
    a = np.where(a > floor, a, np.nan)
    return np.log10(a)


def _colorbar(ax, cf, label, ticks=None, fmt="%.2f"):
    if _HAS_INSET:
        axins = inset_axes(
            ax,
            width="80%",
            height="3.5%",
            loc="lower left",
            bbox_to_anchor=(0.125, 1.01, 1, 1),
            bbox_transform=ax.transAxes,
            borderpad=0,
        )
        cb = plt.colorbar(cf, cax=axins, orientation="horizontal")
    else:
        cb = plt.colorbar(
            cf, ax=ax, orientation="horizontal", pad=0.02, fraction=0.04, aspect=40
        )
    cb.set_label(label, fontsize=15)
    cb.ax.xaxis.set_ticks_position("top")
    cb.ax.xaxis.set_label_position("top")
    cb.ax.tick_params(rotation=45, labelsize=11)
    if ticks is not None:
        cb.set_ticks(ticks)
    cb.ax.xaxis.set_major_formatter(FormatStrFormatter(fmt))
    return cb


# ============================================================================
# Bundle discovery
# ============================================================================


def find_bundles(npzdir: Path, pattern: Optional[str]) -> List[Path]:
    """Return sorted list of .npz bundle files."""
    if pattern:
        files = sorted(npzdir.glob(pattern), key=_natural_key)
    else:
        files = sorted(npzdir.glob("bundle_it*.npz"), key=_natural_key)
    if not files:
        raise SystemExit(
            f"No bundle files found in {npzdir} "
            f"(pattern={'bundle_it*.npz' if not pattern else pattern})."
        )
    return files


# ============================================================================
# Core plot function  (identical logic to make_five_panel_plot in original)
# ============================================================================


def make_five_panel_plot(
    bundle: Dict,
    outpath: Path,
    title: str,
    CU_TO_GAUSS: float,
    use_tex: bool,
) -> None:
    if use_tex:
        plt.rcParams.update(
            {
                "text.usetex": True,
                "font.family": "Serif",
                "font.serif": "Computer Modern",
            }
        )
    else:
        plt.rcParams.update({"text.usetex": False, "font.family": "DejaVu Sans"})
    plt.rcParams["font.size"] = 12

    X = bundle["X"].astype(float)
    Y = bundle["Y"].astype(float)
    ax0n = bundle["ax0"]
    ax1n = bundle["ax1"]

    log_bi   = safe_log10(bundle["beta_inv"])
    log_bpol = safe_log10(bundle["bpol"] * CU_TO_GAUSS)
    log_btor = safe_log10(bundle["btor"] * CU_TO_GAUSS)
    W        = bundle["W"].astype(float)
    Ye       = bundle["Ye"].astype(float)
    rho_cgs  = bundle["rho_log_cgs"]

    # Fixed colormap limits (identical to original)
    bi_lo, bi_hi   = -2.0, 2.0    # log10(beta^-1)
    bp_lo, bp_hi   = 12.0, 16.0   # log10(b_pol [G])
    bt_lo, bt_hi   = 12.0, 16.0   # log10(b_tor [G])
    wmax           = 1.1           # Lorentz factor W upper limit
    ye_lo, ye_hi   = 0.01, 0.6    # electron fraction Y_e

    fig  = plt.figure(figsize=(30, 12))
    ax1p = fig.add_subplot(151)
    ax2p = fig.add_subplot(152, sharey=ax1p)
    ax3p = fig.add_subplot(153, sharey=ax2p)
    ax4p = fig.add_subplot(154, sharey=ax3p)
    ax5p = fig.add_subplot(155, sharey=ax4p)

    N = 100
    cf1 = ax1p.contourf(
        X, Y, log_bi,
        levels=np.linspace(bi_lo, bi_hi, N),
        cmap="seismic",
        extend="both",
    )
    cf2 = ax2p.contourf(
        X, Y, log_bpol,
        levels=np.linspace(bp_lo, bp_hi, N),
        cmap="plasma",
        extend="both",
    )
    cf3 = ax3p.contourf(
        X, Y, log_btor,
        levels=np.linspace(bt_lo, bt_hi, N),
        cmap=_rocket_cmap(),
        extend="both",
    )
    cf4 = ax4p.contourf(
        X, Y, W,
        levels=np.linspace(1.0, wmax, N),
        cmap="jet",
        extend="max",
    )
    cf5 = ax5p.contourf(
        X, Y, Ye,
        levels=np.linspace(ye_lo, ye_hi, N),
        cmap="PRGn",
        extend="both",
    )

    def L(s):
        return rf"${s}$" if use_tex else s

    _colorbar(ax1p, cf1, L(r"\log_{10}(\beta^{-1})"),        ticks=np.linspace(bi_lo, bi_hi, 9))
    _colorbar(ax2p, cf2, L(r"\log_{10}(b_{\rm pol}\ [G])"),  ticks=np.linspace(bp_lo, bp_hi, 9))
    _colorbar(ax3p, cf3, L(r"\log_{10}(b_{\rm tor}\ [G])"),  ticks=np.linspace(bt_lo, bt_hi, 9))
    _colorbar(ax4p, cf4, L(r"W"),                             ticks=np.linspace(1.0, wmax, 7))
    _colorbar(ax5p, cf5, L(r"Y_e"),                           ticks=np.linspace(ye_lo, ye_hi, 7))

    # Density contours overlaid on every panel
    for ax in [ax1p, ax2p, ax3p, ax4p, ax5p]:
        try:
            ax.contour(
                X, Y, rho_cgs,
                levels=[np.log10(1e10), 12, 13, 14],
                colors=["deepskyblue"],
                linewidths=1.2,
                linestyles=[":", "-.", "--", "-"],
            )
        except Exception:
            pass

    xlabel = L(rf"{ax0n}\ [\rm km]")
    ylabel = L(rf"{ax1n}\ [\rm km]")
    xmin, xmax = float(X[0, 0]),  float(X[0, -1])
    ymin, ymax = float(Y[0, 0]),  float(Y[-1, 0])

    ax1p.set_ylabel(ylabel, fontsize=22)
    for ax in [ax1p, ax2p, ax3p, ax4p, ax5p]:
        ax.set_xlabel(xlabel, fontsize=22)
        ax.tick_params(
            which="both",
            top=True, right=True,
            direction="in",
            length=6,
            color="white",
            labelsize=16,
        )
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
    for ax in [ax2p, ax3p, ax4p]:
        ax.tick_params(labelleft=False)
    ax5p.tick_params(labelleft=False, labelright=True)

    ax1p.text(
        xmin + 0.03 * (xmax - xmin),
        ymax - 0.08 * (ymax - ymin),
        title,
        color="white",
        fontsize=17,
        bbox=dict(boxstyle="round,pad=0.2", fc="black", alpha=0.45),
    )

    plt.subplots_adjust(wspace=0.02)
    ensure_dir(outpath.parent)
    fig.savefig(outpath, bbox_inches="tight", dpi=180)
    plt.close(fig)
    print(f"  -> {outpath}")


# ============================================================================
# Load a single .npz bundle and return the dict expected by make_five_panel_plot
# ============================================================================


def load_bundle(npz_path: Path) -> Tuple[Dict, str, int, float, float]:
    """
    Returns
    -------
    bundle : dict  with keys expected by make_five_panel_plot
    plane  : str
    iteration : int
    time_ms   : float
    time_plot_ms : float
    """
    data = np.load(npz_path, allow_pickle=True)

    def _scalar(key, default):
        v = data.get(key, default)
        if v is None:
            return default
        v = np.asarray(v)
        return v.item() if v.ndim == 0 else default

    iteration    = int(_scalar("iteration", -1))
    time_ms      = float(_scalar("time_ms", float("nan")))
    time_plot_ms = float(_scalar("time_plot_ms", time_ms))

    # Decode plane / axis names (stored as 0-d object arrays of bytes/str)
    def _str(key, default="?"):
        v = data.get(key, None)
        if v is None:
            return default
        v = np.asarray(v).item()
        return v.decode() if isinstance(v, bytes) else str(v)

    plane = _str("plane", "xz")
    ax0n  = _str("axis0_name", plane[0])
    ax1n  = _str("axis1_name", plane[1])

    bundle = {
        "X":           data["X"],
        "Y":           data["Y"],
        "beta_inv":    data["plasma_beta_inv"],
        "bpol":        data["b_pol"],
        "btor":        data["b_tor"],
        "W":           data["W"],
        "Ye":          data["Ye"],
        "rho_log_cgs": data["rho_log_cgs"],
        "ax0":         ax0n,
        "ax1":         ax1n,
    }
    return bundle, plane, iteration, time_ms, time_plot_ms


# ============================================================================
# Main
# ============================================================================


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--npzdir",
        required=True,
        help="Directory containing bundle_it*.npz files (e.g. derived_bfield_slices/xz).",
    )
    ap.add_argument(
        "--pattern",
        default=None,
        help="Glob pattern to filter bundles (default: 'bundle_it*.npz').",
    )
    ap.add_argument(
        "--outdir",
        default=None,
        help=(
            "Output directory for PNG files. "
            "Default: <npzdir>/plots"
        ),
    )
    ap.add_argument(
        "--tmerg-ms",
        type=float,
        default=None,
        help=(
            "Merger time in ms. If set the title shows t - t_merg; "
            "otherwise uses the time_plot_ms stored in each bundle."
        ),
    )
    ap.add_argument("--no-tex", action="store_true", help="Disable LaTeX rendering.")
    ap.add_argument(
        "--cu-to-gauss",
        type=float,
        default=_DEFAULT_CU_TO_GAUSS,
        help=(
            f"B-field code-unit → Gauss conversion factor "
            f"(default {_DEFAULT_CU_TO_GAUSS:.4e})."
        ),
    )
    args = ap.parse_args()

    use_tex      = not args.no_tex
    CU_TO_GAUSS  = args.cu_to_gauss
    npzdir       = Path(args.npzdir).resolve()
    outdir       = Path(args.outdir).resolve() if args.outdir else npzdir / "plots"
    ensure_dir(outdir)

    bundles = find_bundles(npzdir, args.pattern)
    print(f"Found {len(bundles)} bundle(s) in {npzdir}")

    for npz_path in bundles:
        print(f"\n[{npz_path.name}]")
        try:
            bundle, plane, iteration, time_ms, time_plot_ms = load_bundle(npz_path)
        except Exception as e:
            print(f"  [skip] failed to load: {e}", file=sys.stderr)
            continue

        # Build title
        if args.tmerg_ms is not None:
            t_plot = time_ms - args.tmerg_ms
            title = (
                f"t - t_merg = {t_plot:.2f} ms"
                if np.isfinite(t_plot)
                else (f"it={iteration}" if iteration >= 0 else npz_path.stem)
            )
        else:
            t_use = time_plot_ms if np.isfinite(time_plot_ms) else time_ms
            title = (
                f"t = {t_use:.2f} ms"
                if np.isfinite(t_use)
                else (f"it={iteration}" if iteration >= 0 else npz_path.stem)
            )

        it_str  = f"{iteration:08d}" if iteration >= 0 else npz_path.stem
        fig_path = outdir / f"five_panel_{plane}_it{it_str}.png"

        try:
            make_five_panel_plot(bundle, fig_path, title, CU_TO_GAUSS, use_tex)
        except Exception as e:
            print(f"  [skip] plotting failed: {e}", file=sys.stderr)
            continue

        print("  [ok]")

    print(f"\nDone. Figures written to {outdir}")


if __name__ == "__main__":
    main()
