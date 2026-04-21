#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from kuibit.simdir import SimDir

# ============================================================
# Defaults
# ============================================================

dir_name = "/mnt/raarchive/cassing/FIL_runs/"
gomega = 0.009                  # FUKA orbital omega from .info file
dist = 500                      # GW detector radius in Msun
tMsolTomsec = 0.00493           # 1 Msun in ms

# orbital period, not directly used but kept from original script
T = 2.0 * np.pi / gomega


# PLOT
lw = 1.0
colours = ["black", "red"]

# ============================================================
# Helpers
# ============================================================

def get_psi4_22(sim_path, det_radius):
    """
    Load the l=m=2 Psi4 mode from a simulation directory.

    Returns
    -------
    t : np.ndarray
        Time array in code units (Msun).
    y : np.ndarray
        Complex Psi4(2,2) mode values.
    """
    sim_dir = SimDir(sim_path)
    psi4 = sim_dir.gws[det_radius].get_psi4_lm(2, 2)
    return np.asarray(psi4.t), np.asarray(psi4.y)


def merge_early_signal_from_reference(t_ref, y_ref, t_restart, y_restart):
    """
    Prepend the early part of a reference signal to a restarted signal.

    The restarted signal is assumed to begin later in time.
    We take all reference points with t < t_restart[0], then append
    the restarted signal.

    Parameters
    ----------
    t_ref, y_ref : arrays
        Reference signal (here: M1).
    t_restart, y_restart : arrays
        Restarted signal (here: no-M1).

    Returns
    -------
    t_merged, y_merged : arrays
        Combined time series.
    """
    if len(t_restart) == 0:
        raise ValueError("Restart signal is empty.")
    if len(t_ref) == 0:
        raise ValueError("Reference signal is empty.")

    t0_restart = t_restart[0]

    early_mask = t_ref < t0_restart

    t_early = t_ref[early_mask]
    y_early = y_ref[early_mask]

    t_merged = np.concatenate([t_early, t_restart])
    y_merged = np.concatenate([y_early, y_restart])

    return t_merged, y_merged


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Plot Psi4(2,2) for M1 and no-M1, and prepend early M1 inspiral to restarted no-M1."
    )
    parser.add_argument(
        "--root",
        default=dir_name,
        help="Root directory containing the simulations.",
    )
    parser.add_argument(
        "--sim-m1",
        default="M139_Bseed_HIGHRES_BNS_BHBL/Supermuc",
        help="Relative path to the M1 simulation.",
    )
    parser.add_argument(
        "--sim-nom1",
        default="noM1_M139_Bseed_HIGHRES_BNS_BHBL/Supermuc",
        help="Relative path to the no-M1 simulation.",
    )
    parser.add_argument(
        "--dist",
        type=float,
        default=dist,
        help="GW detector extraction radius.",
    )
    parser.add_argument(
        "--outdir",
        default=None,
        help="Output directory. If not given, uses <root>/<sim-m1>/Plot_Scripts/Psi4/",
    )
    parser.add_argument(
        "--outname",
        default="Psi4",
        help="Base name for output files.",
    )
    parser.add_argument(
        "--xmin",
        type=float,
        default=0.0,
        help="Minimum x in ms.",
    )
    parser.add_argument(
        "--xmax",
        type=float,
        default=70.0,
        help="Maximum x in ms.",
    )
    parser.add_argument(
        "--ymin",
        type=float,
        default=-0.75e-5,
        help="Minimum y.",
    )
    parser.add_argument(
        "--ymax",
        type=float,
        default=0.75e-5,
        help="Maximum y.",
    )
    args = parser.parse_args()

    sim_m1_path = os.path.join(args.root, args.sim_m1)
    sim_nom1_path = os.path.join(args.root, args.sim_nom1)

    print(f"Loading M1 simulation:    {sim_m1_path}")
    t_m1, y_m1 = get_psi4_22(sim_m1_path, args.dist)

    print(f"Loading no-M1 simulation: {sim_nom1_path}")
    t_nom1, y_nom1 = get_psi4_22(sim_nom1_path, args.dist)

    # Build combined no-M1 curve:
    # early inspiral from M1 + restarted no-M1 from its own start onward
    t_nom1_combined, y_nom1_combined = merge_early_signal_from_reference(
        t_m1, y_m1, t_nom1, y_nom1
    )

    # Plot
    fig, ax = plt.subplots(figsize=(10, 5), dpi=128)

    # M1 original
    ax.plot(
        t_m1 * tMsolTomsec,
        y_m1.real,
        label="M1",
        linewidth=lw,
        color=colours[0],
    )

    # no-M1 with prepended early inspiral from M1
    ax.plot(
        t_nom1_combined * tMsolTomsec,
        y_nom1_combined.real,
        label="no-M1",
        linewidth=lw,
        color=colours[1],
    )

    ax.set_ylabel(r"$\Psi_{4}^{22}$", fontsize=15)
    ax.set_xlabel(r"$t~[\rm{ms}]$", fontsize=15)
    ax.set_xlim(args.xmin, args.xmax)
    ax.set_ylim(args.ymin, args.ymax)

    ax.legend(fontsize=12, frameon=False)

    # Output path
    if args.outdir is None:
        outdir = os.path.join(args.root, args.sim_m1, "Plot_Scripts", "Psi4")
    else:
        outdir = args.outdir

    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, args.outname)

    fig.savefig(outpath + ".pdf", bbox_inches="tight", dpi=180)
    fig.savefig(outpath + ".png", bbox_inches="tight", dpi=180)

    print(f"Saved: {outpath}.pdf")
    print(f"Saved: {outpath}.png")


if __name__ == "__main__":
    main()
