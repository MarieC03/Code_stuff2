#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    "font.size": 12,
    "font.family": "Serif",
    "font.serif": "Computer Modern",
})


def resolve_var_to_csv(varname: str) -> str:
    mapping = {
        "rhob-max": "rho_b_maximum.csv",
        "rhob-min": "rho_b_minimum.csv",
        "rhob-ave": "rho_b_average.csv",
        "temp-max": "temp_maximum.csv",
        "temp-min": "temp_minimum.csv",
        "temp-ave": "temp_average.csv",
        "ye-max": "ye_maximum.csv",
        "ye-min": "ye_minimum.csv",
        "ye-ave": "ye_average.csv",
        "alp-max": "alp_maximum.csv",
        "alp-min": "alp_minimum.csv",
        "alp-ave": "alp_average.csv",
        "p-max": "P_maximum.csv",
        "p-min": "P_minimum.csv",
        "p-ave": "P_average.csv",
        "smallb2-max": "smallb2_maximum.csv",
        "smallb2-min": "smallb2_minimum.csv",
        "smallb2-ave": "smallb2_average.csv",
        "bx-max": "Bvec0_maximum.csv",
        "bx-min": "Bvec0_minimum.csv",
        "bx-ave": "Bvec0_average.csv",
        "by-max": "Bvec1_maximum.csv",
        "by-min": "Bvec1_minimum.csv",
        "by-ave": "Bvec1_average.csv",
        "bz-max": "Bvec2_maximum.csv",
        "bz-min": "Bvec2_minimum.csv",
        "bz-ave": "Bvec2_average.csv",
        "emag": "magnetic_Emag.csv",
        "emagpol": "magnetic_EmagPOL.csv",
        "emagtor": "magnetic_EmagTOR.csv",
        "mej": "M_ej.csv",
        "lem": "L_EM.csv",
    }

    key = varname.strip().lower()
    if key not in mapping:
        known = ", ".join(sorted(mapping.keys()))
        raise ValueError(f"Unknown --var '{varname}'. Known options: {known}")
    return mapping[key]


def make_ylabel(varname: str) -> str:
    labels = {
        "rhob-max": r"$\rho_b^{\max}$",
        "rhob-min": r"$\rho_b^{\min}$",
        "rhob-ave": r"$\langle \rho_b \rangle$",
        "temp-max": r"$T^{\max}$",
        "temp-min": r"$T^{\min}$",
        "temp-ave": r"$\langle T \rangle$",
        "ye-max": r"$Y_e^{\max}$",
        "ye-min": r"$Y_e^{\min}$",
        "ye-ave": r"$\langle Y_e \rangle$",
        "alp-max": r"$\alpha^{\max}$",
        "alp-min": r"$\alpha^{\min}$",
        "alp-ave": r"$\langle \alpha \rangle$",
        "p-max": r"$P^{\max}$",
        "p-min": r"$P^{\min}$",
        "p-ave": r"$\langle P \rangle$",
        "smallb2-max": r"$\left(b^2\right)^{\max}$",
        "smallb2-min": r"$\left(b^2\right)^{\min}$",
        "smallb2-ave": r"$\langle b^2 \rangle$",
        "bx-max": r"$B_x^{\max}$",
        "bx-min": r"$B_x^{\min}$",
        "bx-ave": r"$\langle B_x \rangle$",
        "by-max": r"$B_y^{\max}$",
        "by-min": r"$B_y^{\min}$",
        "by-ave": r"$\langle B_y \rangle$",
        "bz-max": r"$B_z^{\max}$",
        "bz-min": r"$B_z^{\min}$",
        "bz-ave": r"$\langle B_z \rangle$",
        "emag": r"$E_{\rm mag}$",
        "emagpol": r"$E_{\rm mag,pol}$",
        "emagtor": r"$E_{\rm mag,tor}$",
        "mej": r"$\dot{M}_{\rm ej}$",
        "lem": r"$L_{\rm EM}$",
    }
    return labels.get(varname.strip().lower(), "value")


def load_csv(filepath: Path) -> pd.DataFrame:
    df = pd.read_csv(filepath)

    required = {"value"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{filepath} is missing required columns: {sorted(missing)}")

    if "time_ms" not in df.columns:
        raise ValueError(f"{filepath} does not contain a 'time_ms' column.")

    df = df.copy()
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=["time_ms", "value"])
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--var", type=str, help="Shortcut variable name, e.g. rhob-max")
    parser.add_argument("--file", type=str, help="Direct CSV file path")
    parser.add_argument("--root", type=str, default=".", help="Directory containing the combined CSV files")
    parser.add_argument("--xmin", type=float, default=None)
    parser.add_argument("--xmax", type=float, default=None)
    parser.add_argument("--ymin", type=float, default=None)
    parser.add_argument("--ymax", type=float, default=None)
    parser.add_argument("--logy", action="store_true")
    parser.add_argument("--save", type=str, default=None, help="Optional output filename, e.g. rhob_max.pdf")
    args = parser.parse_args()

    if args.file:
        csvfile = Path(args.file).resolve()
        ylabel = "value"
        plot_name = csvfile.stem
    elif args.var:
        filename = resolve_var_to_csv(args.var)
        csvfile = (Path(args.root) / filename).resolve()
        ylabel = make_ylabel(args.var)
        plot_name = args.var.replace("-", "_")
    else:
        raise ValueError("Provide either --var or --file")

    if not csvfile.exists():
        raise FileNotFoundError(f"{csvfile} not found")

    print(f"Reading: {csvfile}")
    df = load_csv(csvfile)

    x = df["time_ms"].to_numpy(dtype=float)
    y = df["value"].to_numpy(dtype=float)

    fig = plt.figure(figsize=(6.8, 6.0))
    ax = fig.add_subplot(111)

    if args.logy:
        mask = np.isfinite(x) & np.isfinite(y) & (y > 0.0)
        ax.semilogy(x[mask], y[mask], linewidth=2.0)
    else:
        mask = np.isfinite(x) & np.isfinite(y)
        ax.plot(x[mask], y[mask], linewidth=2.0)

    ax.set_xlabel(r"$t\,[{\rm ms}]$", fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)

    ax.tick_params(which="both", direction="in", top=True, right=True, length=7)
    ax.minorticks_on()
    ax.grid(alpha=0.25)

    if args.xmin is not None or args.xmax is not None:
        ax.set_xlim(args.xmin, args.xmax)
    if args.ymin is not None or args.ymax is not None:
        ax.set_ylim(args.ymin, args.ymax)

    fig.tight_layout()

    if args.save:
        out = Path(args.save)
    else:
        out = Path(f"{plot_name}.pdf")

    fig.savefig(out, bbox_inches="tight")
    print(f"Saved {out}")

    plt.show()


if __name__ == "__main__":
    main()