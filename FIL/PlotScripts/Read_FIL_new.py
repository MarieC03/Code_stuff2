# -*- coding: utf-8 -*-
"""
Plot all available kuibit grid variables on a chosen plane.

Features:
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

  This means:
    if a structure is located at r = 15 in code units,
    it will be displayed at r = 15 * Rscale in km.

So:
  - interpolation/resampling stays in code units
  - displayed coordinates are rescaled to km

Examples:
  python3 plot_all_vars.py
  python3 plot_all_vars.py --plane xy
  python3 plot_all_vars.py --plane xz --Onedir AllPlots
  python3 plot_all_vars.py --plane xy --Onedir AllPlots --log
"""

import os
import argparse
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

def setup_matplotlib():
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

def safe_mkdir(path):
    os.makedirs(path, exist_ok=True)


def sanitize_name(name):
    """Make variable names safe for file/folder names."""
    return str(name).replace("/", "_").replace(" ", "_")


def code_to_km(x):
    return x * Rscale


def choose_dimension(gf, plane):
    """
    Return (vars_dir, dim_name, is_3d_source)
    """
    if plane == "xy":
        return gf.xy, "xy", False
    if plane == "xz":
        return gf.xz, "xz", False
    if plane == "xyz":
        return gf.xyz, "xyz", True
    raise ValueError(f"Unsupported plane '{plane}'")


def get_uniform_grids(args):
    """
    IMPORTANT:
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
    """
    Plot/image extent in km.
    The sampled data live on code-unit coordinates, but are DISPLAYED in km.
    """
    return (
        code_to_km(args.plot_xmin),
        code_to_km(args.plot_xmax),
        code_to_km(args.plot_ymin),
        code_to_km(args.plot_ymax),
    )


def get_axis_labels(dim_name, is_3d_source):
    if dim_name == "xy" or is_3d_source:
        return r"$x \,[\rm km]$", r"$y \,[\rm km]$", "_xy_it_"
    if dim_name == "xz":
        return r"$x \,[\rm km]$", r"$z \,[\rm km]$", "_xz_it_"
    raise ValueError(f"Unknown dim_name={dim_name}")


def get_iterations(one_var):
    iterations = np.array(sorted(one_var.available_iterations), dtype=int)
    if len(iterations) == 0:
        raise ValueError("No iterations available.")
    return iterations


def get_merger_iteration(one_var, iterations, it_merger):
    if it_merger == 0:
        return int(iterations[0])

    if it_merger in iterations:
        return int(it_merger)

    print(
        f"WARNING: merger iteration {it_merger} not found. "
        f"Using first available iteration {iterations[0]}."
    )
    return int(iterations[0])


def get_iteration_range(iterations, i_start, i_end):
    if i_start is None:
        i0 = 0
    else:
        i0 = max(0, int(i_start))

    if i_end is None:
        i1 = len(iterations) - 1
    else:
        i1 = min(len(iterations) - 1, int(i_end))

    return i0, i1


def time_ms_relative(one_var, iteration, merger_iteration):
    t_merger = one_var.time_at_iteration(merger_iteration) * MSUN_TO_MS
    t_curr = one_var.time_at_iteration(iteration) * MSUN_TO_MS - t_merger
    return round(t_curr, 4)


def read_variable_as_2d(one_var, iteration, grid_2d, grid_3d, is_3d_source):
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


def transform_for_plot(data_2d, use_log):
    """
    Return data prepared for plotting/saving.

    If use_log:
      - values <= 0 are set to NaN
      - log10 is applied to positive values
    """
    data = np.array(data_2d, dtype=float, copy=True)

    if not use_log:
        return data

    mask = data > 0.0
    out = np.full_like(data, np.nan, dtype=float)
    out[mask] = np.log10(data[mask])
    return out


def display_label(quantity_name, use_log):
    if use_log:
        return rf"$\log_{{10}}({quantity_name})$"
    return quantity_name


def save_flattened_data(quantity_name, data_2d, tms_name, data_dir, eos_name, use_log):
    arr = data_2d.reshape(-1, 1)
    df = pd.DataFrame(arr)

    suffix = "_log10" if use_log else ""
    filename = os.path.join(
        data_dir,
        f"{sanitize_name(quantity_name)}{suffix}_2D_np_saved{eos_name}{tms_name}.txt"
    )
    df.to_csv(filename, header=False, index=False)


def get_output_dir(quantity_name, dim_name, onedir):
    """
    Decide where to save plots:
    - if onedir is given: all variables go there
    - else: one folder per variable
    """
    if onedir is not None:
        return onedir
    return f"Plots_{sanitize_name(quantity_name)}_{dim_name}"


def make_2d_plot(quantity_name, data_2d, iteration, dim_name, out_dir, args, is_3d_source):
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

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size="3%", pad=0.3)
    cbar = fig.colorbar(im, orientation="horizontal", cax=cax)
    cbar.ax.tick_params(labelsize=12)
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


def make_1d_plots(quantity_name, var_2d, iteration, dim_name, out_dir, args):
    if not args.make_1d:
        return

    ylabel_txt = display_label(quantity_name, args.log)

    # x-slice at fixed second coordinate = 0 (IN CODE UNITS)
    try:
        fig = plt.figure(figsize=(6, 6), dpi=128)
        ax = fig.add_subplot()

        yslice_code = 0.0
        qx = var_2d.sliced([None, yslice_code])

        # coordinate array built in code units, then displayed in km
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
        ax.set_xlabel(r"$x \,[\rm km]$", fontsize=22, style="italic")
        ax.set_ylabel(ylabel_txt, fontsize=22, style="italic")

        suffix = "_log10" if args.log else ""
        filename = os.path.join(
            out_dir,
            f"{sanitize_name(quantity_name)}{suffix}_x_it_{iteration}.jpg"
        )
        if args.save_plots:
            plt.savefig(filename, dpi=400, pad_inches=0, bbox_inches="tight", format="jpg")
        plt.close(fig)
    except Exception as exc:
        print(f"  Could not make x-slice for {quantity_name} at it={iteration}: {exc}")

    # second-direction slice at fixed x=0 (IN CODE UNITS)
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
            ax.set_xlabel(r"$z \,[\rm km]$", fontsize=22, style="italic")
        else:
            gridname = "_y_it_"
            ax.set_xlabel(r"$y \,[\rm km]$", fontsize=22, style="italic")

        ax.set_ylabel(ylabel_txt, fontsize=22, style="italic")

        suffix = "_log10" if args.log else ""
        filename = os.path.join(
            out_dir,
            f"{sanitize_name(quantity_name)}{suffix}{gridname}{iteration}.jpg"
        )
        if args.save_plots:
            plt.savefig(filename, dpi=400, pad_inches=0, bbox_inches="tight", format="jpg")
        plt.close(fig)
    except Exception as exc:
        print(f"  Could not make second slice for {quantity_name} at it={iteration}: {exc}")


# ============================================================
# Argument parsing
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Loop over all available kuibit variables and plot them."
    )

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

    # Bounds stay in ORIGINAL CODE UNITS
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

    parser.add_argument("--it-merger", type=int, default=0, help="Merger iteration. If 0, use first available.")
    parser.add_argument("--it-start", type=int, default=None, help="Start index in available-iteration list.")
    parser.add_argument("--it-end", type=int, default=None, help="End index in available-iteration list.")

    parser.add_argument("--save-plots", action="store_true", default=True, help="Save plots.")
    parser.add_argument("--no-save-plots", dest="save_plots", action="store_false", help="Do not save plots.")

    parser.add_argument("--save-data", action="store_true", default=True, help="Save flattened 2D data.")
    parser.add_argument("--no-save-data", dest="save_data", action="store_false", help="Do not save flattened 2D data.")

    parser.add_argument("--make-1d", action="store_true", help="Also create 1D slice plots.")

    # Slice ranges stay in ORIGINAL CODE UNITS
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

def main():
    args = parse_args()
    setup_matplotlib()

    safe_mkdir("Data_np")
    if args.Onedir is not None:
        safe_mkdir(args.Onedir)

    print(f"Using coordinate rescaling for display only:")
    print(f"  x_km = x_code * {Rscale}")
    print(f"  y_km = y_code * {Rscale}")
    print(f"  z_km = z_code * {Rscale}")
    print("Resampling/interpolation grid remains in code units.")
    print("------------------------------------------------")

    sim = SimDir(".")
    gf = sim.gridfunctions

    print("Available grid functions:")
    print(gf)
    print("------------------------------------------------")

    vars_dir, dim_name, is_3d_source = choose_dimension(gf, args.plane)

    variable_names = sorted(list(vars_dir.keys()))
    print(f"Selected dimension: {dim_name}")
    print(f"Found {len(variable_names)} variables:")
    print(variable_names)
    print("------------------------------------------------")

    if len(variable_names) == 0:
        raise RuntimeError(f"No variables found for dimension '{dim_name}'.")

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

        merger_iteration = get_merger_iteration(one_var, iterations, args.it_merger)
        i0, i1 = get_iteration_range(iterations, args.it_start, args.it_end)

        print(f"  Number of iterations: {len(iterations)}")
        print(f"  First iteration     : {iterations[0]}")
        print(f"  Last iteration      : {iterations[-1]}")
        print(f"  Merger iteration    : {merger_iteration}")
        print(f"  Looping index range : [{i0}, {i1}]")

        for idx in range(i0, i1 + 1):
            iteration = int(iterations[idx])

            try:
                t_curr = time_ms_relative(one_var, iteration, merger_iteration)
                tms_name = str(t_curr)

                print(f"  -> iteration {iteration}, t_postmerger = {t_curr} ms")

                var_2d, data_2d_raw = read_variable_as_2d(
                    one_var=one_var,
                    iteration=iteration,
                    grid_2d=grid_2d,
                    grid_3d=grid_3d,
                    is_3d_source=is_3d_source,
                )

                data_2d_plot = transform_for_plot(data_2d_raw, args.log)

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