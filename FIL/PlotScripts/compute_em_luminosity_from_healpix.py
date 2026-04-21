#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
compute_em_luminosity_from_healpix.py

Reconstruct electromagnetic Poynting flux and luminosity from FIL healpix
surface outputs of the form

    outputXXXX/data/healpix_det_<det>_surf.h5

using ideal-MHD fields stored on the detector sphere.

What it computes
----------------
At each healpix pixel:
    B2   = gamma_ij B^i B^j
    vdotB = gamma_ij v^i B^j

    S^i = alpha * ( B2 v^i - (v.B) B^i ) / (4*pi)

This is then projected onto the outward radial unit normal and integrated
over the sphere:
    L_EM = ∮ S_r dA

The script uses the spatial metric gamma_ij for contractions and for the
proper area element correction on the sphere.

Outputs
-------
ASCII file with columns:
  1  iteration
  2  time_code
  3  time_ms
  4  radius_code
  5  <S_r>_code           (proper-area-weighted average radial Poynting flux)
  6  L_EM_code            (code units)
  7  <S_r>_cgs            (erg s^-1 cm^-2)
  8  L_EM_cgs             (erg s^-1)

Notes
-----
- This is a postprocessed approximation from the stored fields.
- It is not guaranteed to be exactly identical to FIL's internal detector
  Poynting diagnostic, if such a thorn-level diagnostic was used.
- Default healpix ordering is RING. Change with --ordering nest if needed.
"""

from __future__ import annotations

import argparse
import glob
import os
import re
from pathlib import Path

import h5py
import numpy as np

try:
    import healpy as hp
except ImportError as e:
    raise SystemExit(
        "This script requires healpy.\n"
        "Install it, e.g. with:\n"
        "  pip install healpy\n"
        "or in conda:\n"
        "  conda install -c conda-forge healpy"
    ) from e


# ----------------------------------------------------------------------
# Unit conversions from UtilityUnits.py
# ----------------------------------------------------------------------
G     = 6.6738e-11
c     = 299792458.0
M_sun = 1.9885e30

CU_to_s         = (M_sun * G / (c**3))
CU_to_m         = M_sun * G / (c**2)
CU_to_cm        = CU_to_m * 100.0
CU_to_energyCGS = M_sun * c**2 * 1000.0 * 100.0**2   # erg

FOUR_PI = 4.0 * np.pi


def parse_args():
    p = argparse.ArgumentParser(
        description="Compute EM luminosity from FIL healpix surface outputs."
    )
    p.add_argument(
        "--root",
        default=".",
        help="Root directory under which output*/data/healpix_det_<det>_surf.h5 are searched."
    )
    p.add_argument(
        "--detector",
        type=int,
        required=True,
        help="Detector index N in healpix_det_N_surf.h5"
    )
    p.add_argument(
        "--pattern",
        default=None,
        help=(
            "Optional explicit glob pattern. "
            "Default: <root>/output*/data/healpix_det_<detector>_surf.h5"
        )
    )
    p.add_argument(
        "--ordering",
        choices=["ring", "nest"],
        default="ring",
        help="Healpix pixel ordering"
    )
    p.add_argument(
        "--no-alpha-factor",
        action="store_true",
        help="Do not multiply Poynting vector by lapse alpha"
    )
    p.add_argument(
        "--flat-area",
        action="store_true",
        help="Use flat-space pixel area r^2 dOmega instead of metric-corrected proper area"
    )
    p.add_argument(
        "--flat-normal",
        action="store_true",
        help="Use simple Euclidean radial projection instead of metric unit normal"
    )
    p.add_argument(
        "--output",
        default=None,
        help="Output ASCII filename"
    )
    return p.parse_args()


def find_restart_number(path: str) -> int:
    m = re.search(r"output(\d+)", path)
    return int(m.group(1)) if m else -1


def discover_files(root: str, detector: int, pattern: str | None):
    if pattern is None:
        pattern = os.path.join(root, "output*", "data", f"healpix_det_{detector}_surf.h5")
    files = glob.glob(pattern)
    files = sorted(files, key=lambda x: (find_restart_number(x), x))
    if not files:
        raise FileNotFoundError(f"No files found for pattern: {pattern}")
    return files


def metric_dot(gxx, gxy, gxz, gyy, gyz, gzz, ax, ay, az, bx, by, bz):
    return (
        gxx * ax * bx
        + gxy * (ax * by + ay * bx)
        + gxz * (ax * bz + az * bx)
        + gyy * ay * by
        + gyz * (ay * bz + az * by)
        + gzz * az * bz
    )


def build_healpix_geometry(npix: int, ordering: str):
    nside = hp.npix2nside(npix)
    ipix = np.arange(npix, dtype=np.int64)
    x, y, z = hp.pix2vec(nside, ipix, nest=(ordering == "nest"))
    theta, phi = hp.pix2ang(nside, ipix, nest=(ordering == "nest"))

    # Unit radial vector
    erx = x
    ery = y
    erz = z

    # Tangent basis on unit sphere:
    # u_theta = d n / d theta
    uthx = np.cos(theta) * np.cos(phi)
    uthy = np.cos(theta) * np.sin(phi)
    uthz = -np.sin(theta)

    # u_phi = (1/sin(theta)) d n / d phi
    uphx = -np.sin(phi)
    uphy =  np.cos(phi)
    uphz = np.zeros_like(theta)

    dOmega = FOUR_PI / npix
    return nside, dOmega, erx, ery, erz, uthx, uthy, uthz, uphx, uphy, uphz


def compute_one_iteration(
    grp,
    radius,
    dOmega,
    erx, ery, erz,
    uthx, uthy, uthz,
    uphx, uphy, uphz,
    use_alpha_factor=True,
    use_metric_area=True,
    use_metric_normal=True,
):
    alpha = grp["admbase::alp"][...]

    betax = grp["admbase::betax"][...]
    betay = grp["admbase::betay"][...]
    betaz = grp["admbase::betaz"][...]

    gxx = grp["admbase::gxx"][...]
    gxy = grp["admbase::gxy"][...]
    gxz = grp["admbase::gxz"][...]
    gyy = grp["admbase::gyy"][...]
    gyz = grp["admbase::gyz"][...]
    gzz = grp["admbase::gzz"][...]

    bx = grp["hydrobase::bvec[0]"][...]
    by = grp["hydrobase::bvec[1]"][...]
    bz = grp["hydrobase::bvec[2]"][...]

    vx = grp["hydrobase::vel[0]"][...]
    vy = grp["hydrobase::vel[1]"][...]
    vz = grp["hydrobase::vel[2]"][...]

    time_code = float(grp["time"][0])

    # Metric contractions
    B2 = metric_dot(gxx, gxy, gxz, gyy, gyz, gzz, bx, by, bz, bx, by, bz)
    vdotB = metric_dot(gxx, gxy, gxz, gyy, gyz, gzz, vx, vy, vz, bx, by, bz)

    # Ideal-MHD Poynting vector (code units)
    Sx = (B2 * vx - vdotB * bx) / FOUR_PI
    Sy = (B2 * vy - vdotB * by) / FOUR_PI
    Sz = (B2 * vz - vdotB * bz) / FOUR_PI

    if use_alpha_factor:
        Sx *= alpha
        Sy *= alpha
        Sz *= alpha

    # Radial projection
    if use_metric_normal:
        # covariant normal n_i = gamma_ij e_r^j, normalized
        nr_cov_x = gxx * erx + gxy * ery + gxz * erz
        nr_cov_y = gxy * erx + gyy * ery + gyz * erz
        nr_cov_z = gxz * erx + gyz * ery + gzz * erz

        nr_norm = np.sqrt(
            erx * nr_cov_x + ery * nr_cov_y + erz * nr_cov_z
        )
        nr_norm = np.where(nr_norm > 0.0, nr_norm, 1.0)

        Sr = (Sx * nr_cov_x + Sy * nr_cov_y + Sz * nr_cov_z) / nr_norm
    else:
        Sr = Sx * erx + Sy * ery + Sz * erz

    # Area element per solid angle
    if use_metric_area:
        # A_Omega = r^2 sqrt( (u_theta.u_theta)(u_phi.u_phi) - (u_theta.u_phi)^2 )
        uth2 = metric_dot(gxx, gxy, gxz, gyy, gyz, gzz,
                          uthx, uthy, uthz, uthx, uthy, uthz)
        uph2 = metric_dot(gxx, gxy, gxz, gyy, gyz, gzz,
                          uphx, uphy, uphz, uphx, uphy, uphz)
        uthuph = metric_dot(gxx, gxy, gxz, gyy, gyz, gzz,
                            uthx, uthy, uthz, uphx, uphy, uphz)

        geom = uth2 * uph2 - uthuph**2
        geom = np.where(geom > 0.0, geom, 0.0)
        area_per_solid_angle = radius**2 * np.sqrt(geom)
    else:
        area_per_solid_angle = np.full_like(Sr, radius**2)

    dA = area_per_solid_angle * dOmega

    # Proper-area weighted mean radial flux
    area_total = np.sum(dA)
    mean_Sr_code = np.sum(Sr * dA) / area_total if area_total > 0.0 else np.nan

    # Luminosity in code units
    L_code = np.sum(Sr * dA)

    # Convert to cgs
    # flux: energy / time / area
    flux_conv = CU_to_energyCGS / CU_to_s / (CU_to_cm**2)
    lum_conv  = CU_to_energyCGS / CU_to_s

    mean_Sr_cgs = mean_Sr_code * flux_conv
    L_cgs = L_code * lum_conv

    return time_code, mean_Sr_code, L_code, mean_Sr_cgs, L_cgs


def main():
    args = parse_args()

    files = discover_files(args.root, args.detector, args.pattern)

    # Keep latest restart content for duplicate iterations
    data_by_iteration = {}

    radius_ref = None
    geometry = None

    for fpath in files:
        restart_no = find_restart_number(fpath)
        print(f"[INFO] Reading restart output{restart_no:04d}: {fpath}")

        with h5py.File(fpath, "r") as f:
            radius = float(f["/radius"][0])

            if radius_ref is None:
                radius_ref = radius
            elif not np.isclose(radius, radius_ref):
                raise ValueError(
                    f"Radius mismatch: got {radius} in {fpath}, expected {radius_ref}"
                )

            # Find first iteration group to infer npix / geometry
            it_keys = sorted(f["/data"].keys(), key=lambda x: int(x))
            if not it_keys:
                continue

            if geometry is None:
                sample = f[f"/data/{it_keys[0]}/hydrobase::bvec[0]"]
                npix = sample.shape[0]
                geometry = build_healpix_geometry(npix, args.ordering)
                nside = geometry[0]
                print(f"[INFO] Healpix npix={npix}, inferred nside={nside}, ordering={args.ordering}")
                print(f"[INFO] Detector radius = {radius_ref}")

            _, dOmega, erx, ery, erz, uthx, uthy, uthz, uphx, uphy, uphz = geometry

            for it_str in it_keys:
                it = int(it_str)
                grp = f[f"/data/{it_str}"]

                result = compute_one_iteration(
                    grp=grp,
                    radius=radius_ref,
                    dOmega=dOmega,
                    erx=erx, ery=ery, erz=erz,
                    uthx=uthx, uthy=uthy, uthz=uthz,
                    uphx=uphx, uphy=uphy, uphz=uphz,
                    use_alpha_factor=(not args.no_alpha_factor),
                    use_metric_area=(not args.flat_area),
                    use_metric_normal=(not args.flat_normal),
                )

                # Later restarts overwrite earlier ones for identical iterations
                data_by_iteration[it] = result

    if not data_by_iteration:
        raise RuntimeError("No iterations found in the discovered files.")

    its = np.array(sorted(data_by_iteration.keys()), dtype=np.int64)

    time_code = np.array([data_by_iteration[it][0] for it in its])
    mean_Sr_code = np.array([data_by_iteration[it][1] for it in its])
    L_code = np.array([data_by_iteration[it][2] for it in its])
    mean_Sr_cgs = np.array([data_by_iteration[it][3] for it in its])
    L_cgs = np.array([data_by_iteration[it][4] for it in its])

    time_ms = time_code * 1000.0 * CU_to_s

    outname = args.output
    if outname is None:
        outname = f"em_luminosity_det{args.detector}.dat"

    header = (
        "1:iteration  "
        "2:time_code  "
        "3:time_ms  "
        "4:radius_code  "
        "5:mean_Sr_code  "
        "6:L_EM_code  "
        "7:mean_Sr_cgs[erg/s/cm^2]  "
        "8:L_EM_cgs[erg/s]"
    )

    arr = np.column_stack([
        its,
        time_code,
        time_ms,
        np.full_like(time_code, radius_ref, dtype=float),
        mean_Sr_code,
        L_code,
        mean_Sr_cgs,
        L_cgs,
    ])

    np.savetxt(outname, arr, header=header)
    print(f"[INFO] Wrote {outname}")
    print(f"[INFO] Number of unique iterations: {len(its)}")
    print(f"[INFO] Time range [code]: {time_code.min()} .. {time_code.max()}")
    print(f"[INFO] L_EM_cgs range: {np.nanmin(L_cgs):.6e} .. {np.nanmax(L_cgs):.6e} erg/s")


if __name__ == "__main__":
    main()
