# -*- coding: utf-8 -*-
"""
Plot FIL / kuibit grid variables on a chosen plane.

New features added:
  --var <name or pattern[,pattern2,...]>
      Select only chosen variable(s) for 2D plotting. Supports exact names,
      shell-style wildcards, and comma-separated lists.

  --simdir <path>
      Path to the FIL simulation directory. The script inspects available
      output-* directories, reports present/non-empty/empty files, and then
      loads the simulation from this path.

  --tmerger <time in code units = M_sun>
      Merger time in code units. Plot labels are shifted and shown as
      "t-t_mer = ... ms" in a white box at the top left of each 2D plot.

Existing features preserved:
  --Onedir <dirname>
      Save all plots of all variables into one common directory instead of
      creating one directory per variable.

  --log
      Plot log10(variable) instead of variable.
      Values <= 0 are masked as NaN.

Coordinate scaling:
  The simulation data are sampled on a UniformGrid in the ORIGINAL code units.
  After sampling, the grid-point coordinates are interpreted in physical units
  by multiplying all coordinate axes with

      Rscale = 1.47662504 km per code-length unit

So:
  - interpolation/resampling stays in code units
  - displayed coordinates are rescaled to km
"""

from __future__ import annotations

import os
import fnmatch
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable

from kuibit import grid_data as gd
from kuibit.grid_data import UniformGrid
from kuibit.simdir import SimDir


# ============================================================
# Defaults / constants
# ============================================================

DEFAULT_EOS_NAME = "_DD2_"
DEFAULT_NGRID = 400
DEFAULT_COLORMAP = "RdYlBu_r"
MSUN_TO_MS = 0.0049267398258

# km per code-length unit (GM_sun/c^2)
Rscale = 1.47662504


# ============================================================
# Plot style
# ============================================================

def setup_matplotlib() -> None:
    mpl.rcParams["axes.linewidth"] = 1.0
    mpl.rcParams["xtick.major.size"] = 6
    mpl.rcParams["xtick.major.width"] = 0.8
    mpl.rcParams["xtick.minor.size"] = 4.5
    mpl.rcParams["xtick.minor.width"] = 0.6
    mpl.rcParams["ytick.major.size"] = 6
    mpl.rcParams["ytick.major.width"] = 0.8
    mpl.rcParams["ytick.minor.size"] = 4.5
    mpl.rcParams["ytick.minor.width"] = 0.6

    matplotlib.rcParams.update({
        "figure.figsize": [8, 6],
        "legend.fontsize": 14,
        "text.usetex": True,
        "axes.titlesize": 18,
        "axes.labelsize": 18,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
    })


# ============================================================
# Helpers
# ============================================================

def safe_mkdir(path: str | Path) -> None:
    os.makedirs(path, exist_ok=True)


def sanitize_name(name: str) -> str:
    """Make variable names safe for file/folder names."""
    return (str(name).replace("/", "_").replace(" ", "_").replace("[", "_").replace("]", "_"))


def code_to_km(x: float | np.ndarray) -> float | np.ndarray:
    return x * Rscale


def choose_dimension(gf, plane: str):
    """Return (vars_dir, dim_name, is_3d_source)."""
    if plane == "xy":
        return gf.xy, "xy", False
    if plane == "xz":
        return gf.xz, "xz", False
    if plane == "xyz":
        return gf.xyz, "xyz", True
    raise ValueError(f"Unsupported plane '{plane}'")


def get_uniform_grids(args):
    """
    Build UniformGrid in ORIGINAL CODE UNITS.
    This is the correct grid for kuibit resampling/interpolation.
    """
    grid_2d = UniformGrid(
        [args.ngrid, args.ngrid],
        x0=[args.xmin, args.ymin],
        x1=[args.xmax, args.ymax],
    )
    grid_3d = UniformGrid(
        [args.ngrid, args.ngrid, args.ngrid],
        x0=[args.xmin, args.ymin, args.zmin],
        x1=[args.xmax, args.ymax, args.zmax],
    )
    return grid_2d, grid_3d


def get_plot_extent_km(args):
    """Plot/image extent in km."""
    return (
        code_to_km(args.plot_xmin),
        code_to_km(args.plot_xmax),
        code_to_km(args.plot_ymin),
        code_to_km(args.plot_ymax),
    )


def get_axis_labels(dim_name: str, is_3d_source: bool):
    if dim_name == "xy" or is_3d_source:
        return r"$x \, [\rm km]$", r"$y \, [\rm km]$", "_xy_it_"
    if dim_name == "xz":
        return r"$x \, [\rm km]$", r"$z \, [\rm km]$", "_xz_it_"
    raise ValueError(f"Unknown dim_name={dim_name}")


def get_iterations(one_var) -> np.ndarray:
    iterations = np.array(sorted(one_var.available_iterations), dtype=int)
    if len(iterations) == 0:
        raise ValueError("No iterations available.")
    return iterations


def get_merger_iteration(one_var, iterations: np.ndarray, it_merger: int) -> int:
    if it_merger == 0:
        return int(iterations[0])

    if it_merger in iterations:
        return int(it_merger)

    print(
        f"WARNING: merger iteration {it_merger} not found. "
        f"Using first available iteration {iterations[0]}."
    )
    return int(iterations[0])


def get_iteration_range(iterations: np.ndarray, i_start: int | None, i_end: int | None):
    if i_start is None:
        i0 = 0
    else:
        i0 = max(0, int(i_start))

    if i_end is None:
        i1 = len(iterations) - 1
    else:
        i1 = min(len(iterations) - 1, int(i_end))

    if i1 < i0:
        raise ValueError(f"Invalid iteration index range: [{i0}, {i1}]")

    return i0, i1


def time_ms_relative(one_var, iteration: int, merger_iteration: int) -> float:
    t_merger_ms = one_var.time_at_iteration(merger_iteration) * MSUN_TO_MS
    t_curr_ms = one_var.time_at_iteration(iteration) * MSUN_TO_MS - t_merger_ms
    return round(float(t_curr_ms), 4)


def time_ms_relative_from_tmerger(one_var, iteration: int, tmerger_code: float) -> float:
    t_curr_ms = one_var.time_at_iteration(iteration) * MSUN_TO_MS - tmerger_code * MSUN_TO_MS
    return round(float(t_curr_ms), 4)


def read_variable_as_2d(one_var, iteration: int, grid_2d, grid_3d, is_3d_source: bool):
    """
    Return (var_2d, data_2d).

    For xyz source data, read on 3D uniform grid (in code units) and take
    the z-slice at index 0, matching the original script's behavior.
    """
    if is_3d_source:
        var_3d = one_var.read_on_grid(iteration, grid_3d, resample=True)
        data_3d = np.array(var_3d.data_xyz)
        data_2d = data_3d[:, :, 0]
        var_2d = gd.UniformGridData(grid_2d, data_2d)
        return var_2d, data_2d

    var_2d = one_var.read_on_grid(iteration, grid_2d, resample=True)
    data_2d = np.array(var_2d.data_xyz)
    return var_2d, data_2d


def transform_for_plot(data_2d: np.ndarray, use_log: bool) -> np.ndarray:
    """
    Return data prepared for plotting/saving.

    If use_log:
      - values <= 0 are set to NaN
      - log10 is applied to positive values

    For linear plots, keep the raw values unchanged so quantities such as
    smallb2 and Y_e are shown in their native FIL units.
    """
    data = np.array(data_2d, dtype=float, copy=True)

    if not use_log:
        return data

    mask = data > 0.0
    out = np.full_like(data, np.nan, dtype=float)
    out[mask] = np.log10(data[mask])
    return out


def finite_minmax(data: np.ndarray):
    finite = np.asarray(data)[np.isfinite(data)]
    if finite.size == 0:
        return np.nan, np.nan
    return float(np.nanmin(finite)), float(np.nanmax(finite))


def finite_percentiles(data: np.ndarray, qlow: float = 1.0, qhigh: float = 99.0):
    finite = np.asarray(data)[np.isfinite(data)]
    if finite.size == 0:
        return np.nan, np.nan
    return float(np.nanpercentile(finite, qlow)), float(np.nanpercentile(finite, qhigh))


def display_label(quantity_name: str, use_log: bool) -> str:
    if use_log:
        return rf"$\log_{{10}}({quantity_name})$"
    return quantity_name


def save_flattened_data(quantity_name: str, data_2d: np.ndarray, tms_name: str,
                        data_dir: str, eos_name: str, use_log: bool) -> None:
    arr = data_2d.reshape(-1, 1)
    df = pd.DataFrame(arr)

    suffix = "_log10" if use_log else ""
    filename = os.path.join(
        data_dir,
        f"{sanitize_name(quantity_name)}{suffix}_2D_np_saved{eos_name}{tms_name}.txt"
    )
    df.to_csv(filename, header=False, index=False)


def get_output_dir(quantity_name: str, dim_name: str, onedir: str | None) -> str:
    if onedir is not None:
        return onedir
    return f"Plots_{sanitize_name(quantity_name)}_{dim_name}"


def format_time_label(t_ms: float) -> str:
    return rf"$t-t_{{\rm mer}} = {t_ms:.3f}\;{{\rm ms}}$"


def add_time_box(ax, t_ms: float) -> None:
    ax.text(
        0.03,
        0.97,
        format_time_label(t_ms),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=16,
        color="black",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.9, pad=4.0),
        zorder=50,
    )


def make_2d_plot(quantity_name: str, data_2d: np.ndarray, iteration: int, dim_name: str,
                 out_dir: str, args, is_3d_source: bool, t_label_ms: float | None = None) -> None:
    xlabel, ylabel, gridname = get_axis_labels(dim_name, is_3d_source)
    cbar_label = display_label(quantity_name, args.log)
    extent_km = get_plot_extent_km(args)

    fig = plt.figure(figsize=(6, 6), dpi=128)
    ax = fig.add_subplot()

    if args.fixed_cmap:
        im = ax.imshow(
            data_2d,
            origin="lower",
            cmap=args.colormap,
            extent=extent_km,
            vmin=args.vmin,
            vmax=args.vmax,
        )
    else:
        vmin_auto, vmax_auto = finite_percentiles(data_2d, args.auto_vmin_percentile, args.auto_vmax_percentile)
        if np.isfinite(vmin_auto) and np.isfinite(vmax_auto) and vmax_auto > vmin_auto:
            im = ax.imshow(
                data_2d,
                origin="lower",
                cmap=args.colormap,
                extent=extent_km,
                vmin=vmin_auto,
                vmax=vmax_auto,
            )
        else:
            im = ax.imshow(
                data_2d,
                origin="lower",
                cmap=args.colormap,
                extent=extent_km,
            )

    ax.get_yaxis().set_tick_params(which="both", direction="in")
    ax.get_xaxis().set_tick_params(which="both", direction="in")
    ax.xaxis.set_ticks_position("both")
    ax.yaxis.set_ticks_position("both")

    plt.setp(ax.get_xticklabels(), rotation=0, fontsize=15)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=15)

    ax.set_xlabel(xlabel, fontsize=22, style="italic")
    ax.set_ylabel(ylabel, fontsize=22, style="italic")

    if t_label_ms is not None:
        add_time_box(ax, t_label_ms)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size="3%", pad=0.3)
    cbar = fig.colorbar(im, orientation="horizontal", cax=cax)
    cbar.ax.tick_params(labelsize=12)
    cbar.formatter.set_powerlimits((-2, 2))
    cbar.update_ticks()
    cbar.ax.xaxis.set_label_position("top")
    cbar.set_label(
        cbar_label,
        fontsize=18,
        style="italic",
        fontweight="bold",
        color="black",
    )

    plt.tight_layout()

    if args.save_plots:
        suffix = "_log10" if args.log else ""
        filename = os.path.join(
            out_dir,
            f"{sanitize_name(quantity_name)}{suffix}{gridname}{iteration}.jpg"
        )
        plt.savefig(filename, dpi=400, pad_inches=0, bbox_inches="tight", format="jpg")

    plt.close(fig)


def make_1d_plots(quantity_name: str, var_2d, iteration: int, dim_name: str, out_dir: str, args) -> None:
    if not args.make_1d:
        return

    ylabel_txt = display_label(quantity_name, args.log)

    try:
        fig = plt.figure(figsize=(6, 6), dpi=128)
        ax = fig.add_subplot()

        yslice_code = 0.0
        qx = var_2d.sliced([None, yslice_code])

        xcoord_code = np.linspace(args.slice_min, args.slice_max, args.ngrid)
        xcoord_km = code_to_km(xcoord_code)

        qx_np = np.array(qx.data_xyz, dtype=float)

        if args.log:
            mask = qx_np > 0.0
            qx_plot = np.full_like(qx_np, np.nan, dtype=float)
            qx_plot[mask] = np.log10(qx_np[mask])
        else:
            qx_plot = qx_np

        ax.plot(xcoord_km, qx_plot, color="blue")
        ax.set_xlim(code_to_km(args.slice_plot_min), code_to_km(args.slice_plot_max))
        ax.set_xlabel(r"$x \, [\rm km]$", fontsize=22, style="italic")
        ax.set_ylabel(ylabel_txt, fontsize=22, style="italic")

        suffix = "_log10" if args.log else ""
        filename = os.path.join(out_dir, f"{sanitize_name(quantity_name)}{suffix}_x_it_{iteration}.jpg")
        if args.save_plots:
            plt.savefig(filename, dpi=400, pad_inches=0, bbox_inches="tight", format="jpg")
        plt.close(fig)
    except Exception as exc:
        print(f"  Could not make x-slice for {quantity_name} at it={iteration}: {exc}")

    try:
        fig = plt.figure(figsize=(6, 6), dpi=128)
        ax = fig.add_subplot()

        xslice_code = 0.0
        qy = var_2d.sliced([xslice_code, None])

        ycoord_code = np.linspace(args.slice_min, args.slice_max, args.ngrid)
        ycoord_km = code_to_km(ycoord_code)

        qy_np = np.array(qy.data_xyz, dtype=float)

        if args.log:
            mask = qy_np > 0.0
            qy_plot = np.full_like(qy_np, np.nan, dtype=float)
            qy_plot[mask] = np.log10(qy_np[mask])
        else:
            qy_plot = qy_np

        ax.plot(ycoord_km, qy_plot, color="blue")
        ax.set_xlim(code_to_km(args.slice_plot_min), code_to_km(args.slice_plot_max))

        if dim_name == "xz":
            gridname = "_z_it_"
            ax.set_xlabel(r"$z \, [\rm km]$", fontsize=22, style="italic")
        else:
            gridname = "_y_it_"
            ax.set_xlabel(r"$y \, [\rm km]$", fontsize=22, style="italic")

        ax.set_ylabel(ylabel_txt, fontsize=22, style="italic")

        suffix = "_log10" if args.log else ""
        filename = os.path.join(out_dir, f"{sanitize_name(quantity_name)}{suffix}{gridname}{iteration}.jpg")
        if args.save_plots:
            plt.savefig(filename, dpi=400, pad_inches=0, bbox_inches="tight", format="jpg")
        plt.close(fig)
    except Exception as exc:
        print(f"  Could not make second slice for {quantity_name} at it={iteration}: {exc}")


def split_var_patterns(var_argument: str | None) -> list[str]:
    if var_argument is None:
        return []
    return [item.strip() for item in var_argument.split(",") if item.strip()]


def _pattern_has_wildcards(pattern: str) -> bool:
    """Treat [] in FIL variable names literally; only * and ? trigger wildcard matching."""
    return ("*" in pattern) or ("?" in pattern)


def filter_variable_names(variable_names: list[str], patterns: list[str]) -> list[str]:
    if not patterns:
        return variable_names

    selected = []
    for name in variable_names:
        matched = False
        for pat in patterns:
            if _pattern_has_wildcards(pat):
                if fnmatch.fnmatchcase(name, pat):
                    matched = True
                    break
            else:
                if name == pat:
                    matched = True
                    break
        if matched:
            selected.append(name)

    return selected


def summarize_simdir(simdir: Path) -> None:
    print("Simulation directory scan:")
    print(f"  root = {simdir}")

    output_dirs = sorted([p for p in simdir.glob("output*") if p.is_dir()])
    if not output_dirs:
        print("  No output* directories found. kuibit may still load data from the root path.")
        print("------------------------------------------------")
        return

    print(f"  Found {len(output_dirs)} output directories:")
    total_files = 0
    total_empty = 0
    for odir in output_dirs:
        files = [p for p in odir.rglob("*") if p.is_file()]
        n_total = len(files)
        n_empty = sum(1 for p in files if p.stat().st_size == 0)
        n_nonempty = n_total - n_empty
        total_files += n_total
        total_empty += n_empty
        print(f"    {odir.name}: files={n_total}, non-empty={n_nonempty}, empty={n_empty}")

    print(f"  Total files     : {total_files}")
    print(f"  Total empty     : {total_empty}")
    print(f"  Total non-empty : {total_files - total_empty}")
    print("  Empty files are reported here and unreadable/empty data are skipped during plotting.")
    print("------------------------------------------------")


# ============================================================
# Argument parsing
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Loop over kuibit variables and plot FIL grid data."
    )

    parser.add_argument("--simdir", type=str, default=".", help="Path to FIL simulation directory.")
    parser.add_argument("--var", type=str, default=None,
                        help="Variable name or comma-separated list/patterns to plot, e.g. rho or 'rho,press,*Ye*'.")

    parser.add_argument(
        "--plane",
        choices=["xy", "xz", "xyz"],
        default="xy",
        help="Which grid-function set to use.",
    )

    parser.add_argument(
        "--ngrid",
        type=int,
        default=DEFAULT_NGRID,
        help="Uniform-grid resolution in each dimension.",
    )

    parser.add_argument("--eos-name", default=DEFAULT_EOS_NAME, help="String added to saved data filenames.")
    parser.add_argument("--colormap", default=DEFAULT_COLORMAP, help="Matplotlib colormap name.")

    parser.add_argument("--xmin", type=float, default=-30.0)
    parser.add_argument("--xmax", type=float, default=30.0)
    parser.add_argument("--ymin", type=float, default=-30.0)
    parser.add_argument("--ymax", type=float, default=30.0)
    parser.add_argument("--zmin", type=float, default=0.0)
    parser.add_argument("--zmax", type=float, default=30.0)

    parser.add_argument("--plot-xmin", type=float, default=None, help="Image extent xmin in code units. Defaults to --xmin.")
    parser.add_argument("--plot-xmax", type=float, default=None, help="Image extent xmax in code units. Defaults to --xmax.")
    parser.add_argument("--plot-ymin", type=float, default=None, help="Image extent ymin in code units. Defaults to --ymin.")
    parser.add_argument("--plot-ymax", type=float, default=None, help="Image extent ymax in code units. Defaults to --ymax.")

    parser.add_argument("--fixed-cmap", action="store_true", help="Use fixed colorbar range.")
    parser.add_argument("--vmin", type=float, default=0.0, help="Minimum colorbar value.")
    parser.add_argument("--vmax", type=float, default=0.0014, help="Maximum colorbar value.")
    parser.add_argument("--auto-vmin-percentile", type=float, default=1.0,
                        help="Lower percentile for automatic color scaling when --fixed-cmap is not used.")
    parser.add_argument("--auto-vmax-percentile", type=float, default=99.0,
                        help="Upper percentile for automatic color scaling when --fixed-cmap is not used.")

    parser.add_argument("--it-merger", type=int, default=0,
                        help="Merger iteration. Used only if --tmerger is not given. If 0, use first available iteration.")
    parser.add_argument("--tmerger", type=float, default=None,
                        help="Merger time in code units / solar masses. If given, times are shifted by this value.")
    parser.add_argument("--it-start", type=int, default=None, help="Start index in available-iteration list.")
    parser.add_argument("--it-end", type=int, default=None, help="End index in available-iteration list.")

    parser.add_argument("--save-plots", action="store_true", default=True, help="Save plots.")
    parser.add_argument("--no-save-plots", dest="save_plots", action="store_false", help="Do not save plots.")

    parser.add_argument("--save-data", action="store_true", default=True, help="Save flattened 2D data.")
    parser.add_argument("--no-save-data", dest="save_data", action="store_false", help="Do not save flattened 2D data.")

    parser.add_argument("--make-1d", action="store_true", help="Also create 1D slice plots.")

    parser.add_argument("--slice-min", type=float, default=0.0)
    parser.add_argument("--slice-max", type=float, default=50.0)
    parser.add_argument("--slice_plot_min", type=float, default=0.0)
    parser.add_argument("--slice_plot_max", type=float, default=50.0)

    parser.add_argument(
        "--Onedir",
        type=str,
        default=None,
        help="If given, save all plots of all variables into this one directory.",
    )

    parser.add_argument(
        "--log",
        action="store_true",
        help="Plot log10(variable). Values <= 0 are masked as NaN.",
    )

    args = parser.parse_args()

    if args.plot_xmin is None:
        args.plot_xmin = args.xmin
    if args.plot_xmax is None:
        args.plot_xmax = args.xmax
    if args.plot_ymin is None:
        args.plot_ymin = args.ymin
    if args.plot_ymax is None:
        args.plot_ymax = args.ymax

    return args


# ============================================================
# Main
# ============================================================

def main() -> None:
    args = parse_args()
    setup_matplotlib()

    simdir = Path(args.simdir).expanduser().resolve()
    if not simdir.exists() or not simdir.is_dir():
        raise FileNotFoundError(f"Simulation directory not found: {simdir}")

    safe_mkdir("Data_np")
    if args.Onedir is not None:
        safe_mkdir(args.Onedir)

    print("Using coordinate rescaling for display only:")
    print(f"  x_km = x_code * {Rscale}")
    print(f"  y_km = y_code * {Rscale}")
    print(f"  z_km = z_code * {Rscale}")
    print("Resampling/interpolation grid remains in code units.")
    print("------------------------------------------------")

    summarize_simdir(simdir)

    sim = SimDir(str(simdir))
    gf = sim.gridfunctions

    print("Available grid functions:")
    print(gf)
    print("------------------------------------------------")

    vars_dir, dim_name, is_3d_source = choose_dimension(gf, args.plane)

    all_variable_names = sorted(list(vars_dir.keys()))
    selected_patterns = split_var_patterns(args.var)
    variable_names = filter_variable_names(all_variable_names, selected_patterns)

    print(f"Selected dimension: {dim_name}")
    print(f"Total variables found    : {len(all_variable_names)}")
    if selected_patterns:
        print(f"Variable filter patterns : {selected_patterns}")
    print(f"Variables to be plotted  : {len(variable_names)}")
    print(variable_names)
    print("------------------------------------------------")

    if len(variable_names) == 0:
        raise RuntimeError(
            f"No variables selected for dimension '{dim_name}'. "
            f"Available examples: {all_variable_names[:20]}"
        )

    grid_2d, grid_3d = get_uniform_grids(args)

    for quantity_name in variable_names:
        print()
        print("=" * 80)
        print(f"Processing variable: {quantity_name}")
        print("=" * 80)

        out_dir = get_output_dir(quantity_name, dim_name, args.Onedir)
        safe_mkdir(out_dir)

        try:
            one_var = vars_dir[quantity_name]
        except Exception as exc:
            print(f"  Could not access variable '{quantity_name}': {exc}")
            continue

        try:
            iterations = get_iterations(one_var)
        except Exception as exc:
            print(f"  Could not get iterations for '{quantity_name}': {exc}")
            continue

        i0, i1 = get_iteration_range(iterations, args.it_start, args.it_end)

        if args.tmerger is None:
            merger_iteration = get_merger_iteration(one_var, iterations, args.it_merger)
            print(f"  Number of iterations: {len(iterations)}")
            print(f"  First iteration     : {iterations[0]}")
            print(f"  Last iteration      : {iterations[-1]}")
            print(f"  Merger iteration    : {merger_iteration}")
            print(f"  Looping index range : [{i0}, {i1}]")
        else:
            print(f"  Number of iterations: {len(iterations)}")
            print(f"  First iteration     : {iterations[0]}")
            print(f"  Last iteration      : {iterations[-1]}")
            print(f"  Merger time         : {args.tmerger} M_sun (code units)")
            print(f"  Looping index range : [{i0}, {i1}]")

        for idx in range(i0, i1 + 1):
            iteration = int(iterations[idx])

            try:
                if args.tmerger is None:
                    t_curr = time_ms_relative(one_var, iteration, merger_iteration)
                else:
                    t_curr = time_ms_relative_from_tmerger(one_var, iteration, args.tmerger)
                tms_name = str(t_curr)

                print(f"  -> iteration {iteration}, t-t_mer = {t_curr} ms")

                var_2d, data_2d_raw = read_variable_as_2d(
                    one_var=one_var,
                    iteration=iteration,
                    grid_2d=grid_2d,
                    grid_3d=grid_3d,
                    is_3d_source=is_3d_source,
                )

                data_2d_plot = transform_for_plot(data_2d_raw, args.log)

                raw_min, raw_max = finite_minmax(data_2d_raw)
                plot_min, plot_max = finite_minmax(data_2d_plot)
                p01, p99 = finite_percentiles(data_2d_plot, args.auto_vmin_percentile, args.auto_vmax_percentile)
                print(
                    f"     stats raw=[{raw_min:.6e}, {raw_max:.6e}] "
                    f"plot=[{plot_min:.6e}, {plot_max:.6e}] "
                    f"auto-clim=[{p01:.6e}, {p99:.6e}]"
                )

                if data_2d_plot.size == 0:
                    print(f"    Skipping iteration {iteration}: empty array returned.")
                    continue

                if not np.any(np.isfinite(data_2d_plot)):
                    print(f"    Skipping iteration {iteration}: array is empty/all-NaN after processing.")
                    continue

                if args.save_data:
                    save_flattened_data(
                        quantity_name=quantity_name,
                        data_2d=data_2d_plot,
                        tms_name=tms_name,
                        data_dir="Data_np",
                        eos_name=args.eos_name,
                        use_log=args.log,
                    )

                make_2d_plot(
                    quantity_name=quantity_name,
                    data_2d=data_2d_plot,
                    iteration=iteration,
                    dim_name=dim_name,
                    out_dir=out_dir,
                    args=args,
                    is_3d_source=is_3d_source,
                    t_label_ms=t_curr,
                )

                make_1d_plots(
                    quantity_name=quantity_name,
                    var_2d=var_2d,
                    iteration=iteration,
                    dim_name=dim_name,
                    out_dir=out_dir,
                    args=args,
                )

            except Exception as exc:
                print(f"    Skipping iteration {iteration} for '{quantity_name}': {exc}")
                continue

    print("\nDone.")


if __name__ == "__main__":
    main()
