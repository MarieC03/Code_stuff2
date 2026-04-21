#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fil_2d_bfield_plot.py
==================
Reads ALL required fields from the pre-sliced 2D Carpet HDF5 plane files
(data_hdf5_2D) and produces the 5-panel b-field figure.

No 3D data needed. No external UtilityUnits module needed.

Memory note: only tiles overlapping the requested --xmin/xmax/ymin/ymax window
are read into memory and the canvas is restricted to that window.  This avoids
OOM-kills on large high-resolution simulations.

Variables read from data_hdf5_2D/<var>.<plane>.h5:
  Bvec[0], Bvec[1], Bvec[2]   (magnetic field)
  vel[0],  vel[1],  vel[2]     (velocity)
  alp, betax, betay, betaz     (lapse, shift)
  gxx, gxy, gxz, gyy, gyz, gzz  (spatial metric)
  w_lorentz                    (Lorentz factor W)
  Y_e                          (electron fraction)
  rho                          (rest-mass density)
  smallb2                      (comoving b^2)
  press                        (pressure)

Physics computed:
  b^mu  -- comoving magnetic field (Valencia/FIL convention)
  b_pol = sqrt(b_R^2 + b_z^2)  using orthonormal cylindrical basis from metric
  b_tor = |b_phi|
  beta^{-1} = smallb2 / (2 * press)

Usage
-----
  python fil_2d_bfield_plot.py --root /path/to/output0020 --plane xz
  python fil_2d_bfield_plot.py --root /path/to/output0020 --plane xz \\
      --xmin -120 --xmax 120 --ymin 0 --ymax 280 \\
      --outdir bfield_out --no-tex
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
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
# Physical unit conversions  (ET geometric units, M_sun = 1)
# ============================================================================
CU_TO_MS = 4.925490947e-3
CU_TO_DENS_CGS = 6.17714470405638e17
_DEFAULT_CU_TO_GAUSS = 8.35195e19 / np.sqrt(4.0 * np.pi)  # ~2.3548e19 G
Rscale = 1.47662504  # 1 code length = Rscale km

# ============================================================================
# Carpet dataset-name regex
# ============================================================================
DSET_RE = re.compile(r"^(\S+)::(\S+) it=(\d+) tl=(\d+) rl=(\d+)(?: c=(\d+))?$")

# ============================================================================
# Variable aliases
# ============================================================================
ALIASES: Dict[str, List[str]] = {
    "Bx": ["Bvec[0]", "Bx", "bvecx"],
    "By": ["Bvec[1]", "By", "bvecy"],
    "Bz": ["Bvec[2]", "Bz", "bvecz"],
    "vx": ["vel[0]", "vx"],
    "vy": ["vel[1]", "vy"],
    "vz": ["vel[2]", "vz"],
    "alp": ["alp", "alpha"],
    "betax": ["betax"],
    "betay": ["betay"],
    "betaz": ["betaz"],
    "gxx": ["gxx"],
    "gxy": ["gxy"],
    "gxz": ["gxz"],
    "gyy": ["gyy"],
    "gyz": ["gyz"],
    "gzz": ["gzz"],
    "W": ["w_lorentz", "W"],
    "Ye": ["Y_e", "Ye", "ye"],
    "rho": ["rho", "rho_b"],
    "smallb2": ["smallb2", "b2", "bsq"],
    "press": ["press", "P"],
}

ALL_VARS = [
    "Bx",
    "By",
    "Bz",
    "vx",
    "vy",
    "vz",
    "alp",
    "betax",
    "betay",
    "betaz",
    "gxx",
    "gxy",
    "gxz",
    "gyy",
    "gyz",
    "gzz",
    "W",
    "Ye",
    "rho",
    "smallb2",
    "press",
]


# ============================================================================
# Data class
# ============================================================================
@dataclass
class Tile2D:
    file: str
    dset: str
    iteration: int
    tl: int
    rl: int
    c: int
    time: float
    origin: np.ndarray
    delta: np.ndarray
    shape_xy: np.ndarray
    corner: np.ndarray


# ============================================================================
# File-system helpers
# ============================================================================


def _natural_key(p):
    s = str(p)
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)]


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def find_data2d_dirs(root: Path) -> List[Path]:
    root = root.resolve()
    subnames = ["data_hdf5_2D", "data_hdf5", "data_h5", "data"]
    if root.is_dir() and root.name in subnames:
        return [root]
    dirs: List[Path] = []
    seen: set = set()
    if root.name.startswith("output"):
        output_dirs = [root]
    else:
        output_dirs = sorted(
            [p for p in root.iterdir() if p.is_dir() and p.name.startswith("output")],
            key=_natural_key,
        )
    for out in output_dirs:
        for sub in subnames:
            d = out / sub
            if d.is_dir() and d not in seen:
                dirs.append(d)
                seen.add(d)
    return dirs


def find_var_file(root: Path, plane: str, aliases: List[str]) -> List[Path]:
    found: List[Path] = []
    seen: set = set()
    for d in find_data2d_dirs(root):
        try:
            names = os.listdir(d)
        except OSError:
            continue
        for alias in aliases:
            fname = f"{alias}.{plane}.h5"
            if fname in names:
                p = d / fname
                if p not in seen:
                    found.append(p)
                    seen.add(p)
    return sorted(found, key=_natural_key)


# ============================================================================
# Variable store
# ============================================================================


class VariableStore2D:
    def __init__(self, root: Path, plane: str):
        self.root = root
        self.plane = plane
        self.index: Dict[str, Dict[int, Dict[int, Dict[int, Tile2D]]]] = {}

    @staticmethod
    def _read_tiles(fp: Path) -> List[Tile2D]:
        tiles: List[Tile2D] = []
        try:
            with h5py.File(fp, "r") as f:
                for key in f.keys():
                    m = DSET_RE.match(key)
                    if not m:
                        continue
                    _grp, _var, it, tl, rl, c = m.groups()
                    ds = f[key]
                    if ds.ndim != 2:
                        continue
                    origin = np.asarray(ds.attrs["origin"], dtype=float)
                    delta = np.asarray(ds.attrs["delta"], dtype=float)
                    shape_xy = np.asarray(ds.shape[::-1], dtype=int)
                    corner = origin + delta * (shape_xy - 1)
                    time = float(np.real(ds.attrs.get("time", np.nan)))
                    tiles.append(
                        Tile2D(
                            file=fp.as_posix(),
                            dset=key,
                            iteration=int(it),
                            tl=int(tl),
                            rl=int(rl),
                            c=int(c) if c is not None else -1,
                            time=time,
                            origin=origin,
                            delta=delta,
                            shape_xy=shape_xy,
                            corner=corner,
                        )
                    )
        except Exception as e:
            print(f"  [warn] could not read {fp}: {e}", file=sys.stderr)
        return tiles

    def build(self, names: Sequence[str]) -> None:
        for name in names:
            files = find_var_file(self.root, self.plane, ALIASES[name])
            if not files:
                continue
            idx: Dict[int, Dict[int, Dict[int, Tile2D]]] = defaultdict(
                lambda: defaultdict(dict)
            )
            for fp in files:
                for tile in self._read_tiles(fp):
                    if tile.tl == 0:
                        idx[tile.iteration][tile.rl][tile.c] = tile
            if idx:
                self.index[name] = idx

    def missing(self, names: Sequence[str]) -> List[str]:
        return [n for n in names if n not in self.index]

    def common_iterations(self, names: Sequence[str]) -> List[int]:
        sets = [set(self.index[n].keys()) for n in names if n in self.index]
        return sorted(set.intersection(*sets)) if sets else []


# ============================================================================
# Canvas helpers
# ============================================================================


def _centers_to_edges(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    if v.size == 1:
        return np.array([v[0] - 0.5, v[0] + 0.5])
    mids = 0.5 * (v[:-1] + v[1:])
    left = v[0] - 0.5 * (v[1] - v[0])
    right = v[-1] + 0.5 * (v[-1] - v[-2])
    return np.concatenate(([left], mids, [right]))


def _paste_tile(
    canvas: np.ndarray,
    x_global: np.ndarray,
    y_global: np.ndarray,
    x_tile: np.ndarray,
    y_tile: np.ndarray,
    vals: np.ndarray,
) -> None:
    """
    Vectorised paste.  For the common case where tile and canvas share the
    same resolution (fine tiles) this is a single slice assignment.
    For coarse tiles spanning multiple canvas cells it falls back to a
    row-wise loop with vectorised x handling.  ~100x faster than a pure
    Python double loop.
    """
    xce = _centers_to_edges(x_global)
    yce = _centers_to_edges(y_global)
    xte = _centers_to_edges(x_tile)
    yte = _centers_to_edges(y_tile)
    ny_t, nx_t = vals.shape

    x_lo = np.minimum(xte[:-1], xte[1:])
    x_hi = np.maximum(xte[:-1], xte[1:])
    y_lo = np.minimum(yte[:-1], yte[1:])
    y_hi = np.maximum(yte[:-1], yte[1:])

    gx0_all = np.clip(np.searchsorted(xce, x_lo, side="right") - 1, 0, canvas.shape[1])
    gx1_all = np.clip(np.searchsorted(xce, x_hi, side="left"), 0, canvas.shape[1])
    gy0_all = np.clip(np.searchsorted(yce, y_lo, side="right") - 1, 0, canvas.shape[0])
    gy1_all = np.clip(np.searchsorted(yce, y_hi, side="left"), 0, canvas.shape[0])

    # Fast path: every tile cell maps 1-to-1 to one canvas cell
    if np.all(gx1_all - gx0_all == 1) and np.all(gy1_all - gy0_all == 1):
        cx0 = int(gx0_all[0])
        cx1 = int(gx0_all[-1]) + 1
        cy0 = int(gy0_all[0])
        cy1 = int(gy0_all[-1]) + 1
        if cx1 > cx0 and cy1 > cy0:
            src = vals[: cy1 - cy0, : cx1 - cx0]
            mask = np.isfinite(src)
            dst = canvas[cy0:cy1, cx0:cx1]
            dst[mask] = src[mask]
            canvas[cy0:cy1, cx0:cx1] = dst
        return

    # General path: coarse tile, each cell may cover multiple canvas cells
    for j in range(ny_t):
        gy0 = int(gy0_all[j])
        gy1 = int(gy1_all[j])
        if gy1 <= gy0:
            continue
        row = vals[j, :]
        for i in np.where(np.isfinite(row))[0]:
            gx0 = int(gx0_all[i])
            gx1 = int(gx1_all[i])
            if gx1 > gx0:
                canvas[gy0:gy1, gx0:gx1] = row[i]


def join_var(
    store: VariableStore2D,
    name: str,
    iteration: int,
    requested_rl: Optional[int] = None,
    xmin_km: float = -1e30,
    xmax_km: float = 1e30,
    ymin_km: float = -1e30,
    ymax_km: float = 1e30,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """
    Join AMR tiles for (name, iteration) into one canvas.
    Only tiles overlapping the km window are read; canvas is clamped to window.
    Returns: X [km meshgrid], Y [km meshgrid], canvas [float32], time_code.
    """
    tiles: List[Tile2D] = []
    for rl in sorted(store.index[name][iteration].keys()):
        if requested_rl is not None and rl != requested_rl:
            continue
        for c in sorted(store.index[name][iteration][rl].keys()):
            tiles.append(store.index[name][iteration][rl][c])
    if not tiles:
        raise ValueError(f"No tiles for '{name}' at iteration {iteration}")

    tile_info = []
    for t in tiles:
        x_t = t.origin[0] * Rscale + t.delta[0] * Rscale * np.arange(int(t.shape_xy[0]))
        y_t = t.origin[1] * Rscale + t.delta[1] * Rscale * np.arange(int(t.shape_xy[1]))
        dx = abs(t.delta[0]) * Rscale
        dy = abs(t.delta[1]) * Rscale
        # Skip tiles that don't overlap the window
        if x_t[-1] < xmin_km or x_t[0] > xmax_km:
            continue
        if y_t[-1] < ymin_km or y_t[0] > ymax_km:
            continue
        tile_info.append((x_t, y_t, dx, dy, t))

    if not tile_info:
        raise ValueError(
            f"No tiles for '{name}' at it={iteration} overlap "
            f"x=[{xmin_km},{xmax_km}] y=[{ymin_km},{ymax_km}]"
        )

    # Canvas extents at finest resolution, clamped to window
    cxmin = max(min(x[0] for x, _, _, _, _ in tile_info), xmin_km)
    cxmax = min(max(x[-1] for x, _, _, _, _ in tile_info), xmax_km)
    cymin = max(min(y[0] for _, y, _, _, _ in tile_info), ymin_km)
    cymax = min(max(y[-1] for _, y, _, _, _ in tile_info), ymax_km)
    dx = min(d for _, _, d, _, _ in tile_info)
    dy = min(d for _, _, _, d, _ in tile_info)

    nx = int(round((cxmax - cxmin) / dx)) + 1
    ny = int(round((cymax - cymin) / dy)) + 1
    x_g = cxmin + dx * np.arange(nx)
    y_g = cymin + dy * np.arange(ny)

    canvas = np.full((ny, nx), np.nan, dtype=np.float32)
    for x_t, y_t, _dx, _dy, t in sorted(tile_info, key=lambda q: (q[4].rl, q[4].c)):
        with h5py.File(t.file, "r") as f:
            vals = np.asarray(f[t.dset], dtype=np.float32)
        _paste_tile(canvas, x_g, y_g, x_t, y_t, vals)

    X, Y = np.meshgrid(x_g, y_g)
    return X, Y, canvas, float(tiles[0].time)


def read_all_vars(
    store: VariableStore2D,
    iteration: int,
    requested_rl: Optional[int],
    xmin_km: float,
    xmax_km: float,
    ymin_km: float,
    ymax_km: float,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray], float]:
    """Read every required variable restricted to the km window."""
    X = Y = None
    time_ms = float("nan")
    fields: Dict[str, np.ndarray] = {}
    for name in ALL_VARS:
        Xv, Yv, arr, tc = join_var(
            store,
            name,
            iteration,
            requested_rl,
            xmin_km=xmin_km,
            xmax_km=xmax_km,
            ymin_km=ymin_km,
            ymax_km=ymax_km,
        )
        if X is None:
            X, Y = Xv, Yv
            time_ms = tc * CU_TO_MS
        fields[name] = arr
    fields["X"] = X
    fields["Y"] = Y
    return X, Y, fields, time_ms


# ============================================================================
# GRMHD: comoving b-field + cylindrical decomposition
# ============================================================================


def _sdot(ax, ay, az, bx, by, bz, gxx, gxy, gxz, gyy, gyz, gzz):
    return (
        gxx * ax * bx
        + gyy * ay * by
        + gzz * az * bz
        + gxy * (ax * by + ay * bx)
        + gxz * (ax * bz + az * bx)
        + gyz * (ay * bz + az * by)
    )


def _norm(rx, ry, rz, gxx, gxy, gxz, gyy, gyz, gzz):
    n2 = _sdot(rx, ry, rz, rx, ry, rz, gxx, gxy, gxz, gyy, gyz, gzz)
    n = np.where(n2 > 1e-28, np.sqrt(n2), np.nan)
    return rx / n, ry / n, rz / n


def compute_bpol_btor(f: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    Bx = f["Bx"].astype(float)
    By = f["By"].astype(float)
    Bz = f["Bz"].astype(float)
    vx = f["vx"].astype(float)
    vy = f["vy"].astype(float)
    vz = f["vz"].astype(float)
    alp = f["alp"].astype(float)
    bex = f["betax"].astype(float)
    bey = f["betay"].astype(float)
    bez = f["betaz"].astype(float)
    gxx = f["gxx"].astype(float)
    gxy = f["gxy"].astype(float)
    gxz = f["gxz"].astype(float)
    gyy = f["gyy"].astype(float)
    gyz = f["gyz"].astype(float)
    gzz = f["gzz"].astype(float)
    W = f["W"].astype(float)
    X_km = f["X"].astype(float)
    Y_km = f["Y"].astype(float)

    # Lorentz factor: prefer stored W, fall back to computing from v^i
    vsq = _sdot(vx, vy, vz, vx, vy, vz, gxx, gxy, gxz, gyy, gyz, gzz)
    vsq = np.clip(vsq, 0.0, 1.0 - 1e-12)
    W_calc = 1.0 / np.sqrt(1.0 - vsq)
    W = np.where(np.isfinite(W) & (W > 0.5), W, W_calc)

    # Comoving b^mu  (Valencia convention)
    alp_s = np.where(np.abs(alp) > 1e-14, alp, np.nan)
    vdotB = _sdot(Bx, By, Bz, vx, vy, vz, gxx, gxy, gxz, gyy, gyz, gzz)
    b0 = W * vdotB / alp_s
    ux = W * (vx - bex / alp_s)
    uy = W * (vy - bey / alp_s)
    uz = W * (vz - bez / alp_s)
    bx = (Bx + alp_s * b0 * ux) / W
    by = (By + alp_s * b0 * uy) / W
    bz = (Bz + alp_s * b0 * uz) / W
    b2 = _sdot(bx, by, bz, bx, by, bz, gxx, gxy, gxz, gyy, gyz, gzz) - alp_s**2 * b0**2

    # Orthonormal cylindrical basis from spatial metric
    Xc = X_km / Rscale
    Yc = Y_km / Rscale
    z0 = np.zeros_like(Xc)
    z1 = np.ones_like(Xc)

    eRx, eRy, eRz = _norm(Xc, Yc, z0, gxx, gxy, gxz, gyy, gyz, gzz)
    epx, epy, epz = _norm(-Yc, Xc, z0, gxx, gxy, gxz, gyy, gyz, gzz)

    zeR = _sdot(z0, z0, z1, eRx, eRy, eRz, gxx, gxy, gxz, gyy, gyz, gzz)
    zephi = _sdot(z0, z0, z1, epx, epy, epz, gxx, gxy, gxz, gyy, gyz, gzz)
    ezx, ezy, ezz = _norm(
        z0 - zeR * eRx - zephi * epx,
        z0 - zeR * eRy - zephi * epy,
        z1 - zeR * eRz - zephi * epz,
        gxx,
        gxy,
        gxz,
        gyy,
        gyz,
        gzz,
    )

    bR = _sdot(bx, by, bz, eRx, eRy, eRz, gxx, gxy, gxz, gyy, gyz, gzz)
    bphi = _sdot(bx, by, bz, epx, epy, epz, gxx, gxy, gxz, gyy, gyz, gzz)
    bzp = _sdot(bx, by, bz, ezx, ezy, ezz, gxx, gxy, gxz, gyy, gyz, gzz)

    bpol = np.sqrt(np.maximum(bR**2 + bzp**2, 0.0))
    btor = np.abs(bphi)

    return {
        "b2": b2,
        "bpol": bpol,
        "btor": btor,
        "bx_comov": bx,
        "by_comov": by,
        "bz_comov": bz,
    }


# ============================================================================
# Iteration selection
# ============================================================================


def _time_map(store: VariableStore2D) -> Dict[int, float]:
    ref = next((v for v in ["rho", "W", "Ye"] if v in store.index), None)
    if ref is None:
        return {}
    out: Dict[int, float] = {}
    for it, rl_map in store.index[ref].items():
        for rl in sorted(rl_map):
            for c in sorted(rl_map[rl]):
                t = rl_map[rl][c].time
                if np.isfinite(t):
                    out[it] = t * CU_TO_MS
                break
            if it in out:
                break
    return out


def select_iterations(
    store: VariableStore2D, args, tmerg_ms: Optional[float]
) -> List[int]:
    common = store.common_iterations(ALL_VARS)
    if not common:
        raise SystemExit("No iterations found with all required variables present.")
    tmap = _time_map(store)
    result = []
    for it in common:
        if args.itmin is not None and it < args.itmin:
            continue
        if args.itmax is not None and it > args.itmax:
            continue
        t_ms = tmap.get(it, float("nan"))
        if args.tmin is not None or args.tmax is not None:
            if not np.isfinite(t_ms):
                continue
            t_plot = t_ms if tmerg_ms is None else t_ms - tmerg_ms
            if args.tmin is not None and t_plot < args.tmin:
                continue
            if args.tmax is not None and t_plot > args.tmax:
                continue
        result.append(it)
    if not result:
        raise SystemExit("No iterations remain after applying filters.")
    return result


# ============================================================================
# Plotting
# ============================================================================


def safe_log10(arr, floor=1e-30):
    a = np.where(np.isfinite(arr), np.abs(arr), np.nan)
    a = np.where(a > floor, a, np.nan)
    return np.log10(a)


def _robust(arr, qlo=1.0, qhi=99.0, dlo=0.0, dhi=1.0):
    a = arr[np.isfinite(arr)]
    if a.size == 0:
        return dlo, dhi
    lo, hi = float(np.percentile(a, qlo)), float(np.percentile(a, qhi))
    if not (np.isfinite(lo) and np.isfinite(hi) and hi > lo):
        return dlo, dhi
    return lo, hi


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


def _axis_names(plane: str) -> Tuple[str, str]:
    return {"xy": ("x", "y"), "xz": ("x", "z"), "yz": ("y", "z")}[plane]


def make_five_panel_plot(
    bundle: Dict, outpath: Path, title: str, CU_TO_GAUSS: float, use_tex: bool
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

    log_bi = safe_log10(bundle["beta_inv"])
    log_bpol = safe_log10(bundle["bpol"] * CU_TO_GAUSS)
    log_btor = safe_log10(bundle["btor"] * CU_TO_GAUSS)
    W = bundle["W"].astype(float)
    Ye = bundle["Ye"].astype(float)
    rho_cgs = bundle["rho_log_cgs"]

    # Fixed colormap limits (applied uniformly across all iterations)
    bi_lo, bi_hi = -2.0, 2.0        # log10(beta^-1)
    bp_lo, bp_hi = 12.0, 16.0       # log10(b_pol [G])
    bt_lo, bt_hi = 12.0, 16.0       # log10(b_tor [G])
    wmax        = 1.1               # Lorentz factor W upper limit (lower=1.0)
    ye_lo, ye_hi = 0.01, 0.6        # electron fraction Y_e

    fig = plt.figure(figsize=(30, 12))
    ax1p = fig.add_subplot(151)
    ax2p = fig.add_subplot(152, sharey=ax1p)
    ax3p = fig.add_subplot(153, sharey=ax2p)
    ax4p = fig.add_subplot(154, sharey=ax3p)
    ax5p = fig.add_subplot(155, sharey=ax4p)

    N = 100
    cf1 = ax1p.contourf(
        X, Y, log_bi, levels=np.linspace(bi_lo, bi_hi, N), cmap="seismic", extend="both"
    )
    cf2 = ax2p.contourf(
        X,
        Y,
        log_bpol,
        levels=np.linspace(bp_lo, bp_hi, N),
        cmap="plasma",
        extend="both",
    )
    cf3 = ax3p.contourf(
        X,
        Y,
        log_btor,
        levels=np.linspace(bt_lo, bt_hi, N),
        cmap=_rocket_cmap(),
        extend="both",
    )
    cf4 = ax4p.contourf(
        X, Y, W, levels=np.linspace(1.0, wmax, N), cmap="jet", extend="max"
    )
    cf5 = ax5p.contourf(
        X, Y, Ye, levels=np.linspace(ye_lo, ye_hi, N), cmap="PRGn", extend="both"
    )

    def L(s):
        return rf"${s}$" if use_tex else s

    _colorbar(
        ax1p, cf1, L(r"\log_{10}(\beta^{-1})"), ticks=np.linspace(bi_lo, bi_hi, 9)
    )
    _colorbar(
        ax2p, cf2, L(r"\log_{10}(b_{\rm pol}\ [G])"), ticks=np.linspace(bp_lo, bp_hi, 9)
    )
    _colorbar(
        ax3p, cf3, L(r"\log_{10}(b_{\rm tor}\ [G])"), ticks=np.linspace(bt_lo, bt_hi, 9)
    )
    _colorbar(ax4p, cf4, L(r"W"), ticks=np.linspace(1.0, wmax, 7))
    _colorbar(ax5p, cf5, L(r"Y_e"), ticks=np.linspace(ye_lo, ye_hi, 7))

    for ax in [ax1p, ax2p, ax3p, ax4p, ax5p]:
        try:
            ax.contour(
                X,
                Y,
                rho_cgs,
                levels=[np.log10(1e10), 12, 13, 14],
                colors=["deepskyblue"],
                linewidths=1.2,
                linestyles=[":", "-.", "--", "-"],
            )
        except Exception:
            pass

    xlabel = L(rf"{ax0n}\ [\rm km]")
    ylabel = L(rf"{ax1n}\ [\rm km]")
    xmin, xmax = float(X[0, 0]), float(X[0, -1])
    ymin, ymax = float(Y[0, 0]), float(Y[-1, 0])

    ax1p.set_ylabel(ylabel, fontsize=22)
    for ax in [ax1p, ax2p, ax3p, ax4p, ax5p]:
        ax.set_xlabel(xlabel, fontsize=22)
        ax.tick_params(
            which="both",
            top=True,
            right=True,
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
# Main
# ============================================================================


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--root", required=True)
    ap.add_argument("--plane", default="xz", choices=["xy", "xz", "yz"])
    ap.add_argument(
        "--rl",
        type=int,
        default=None,
        help="Force single refinement level (default: all, fine overwrites coarse)",
    )
    ap.add_argument("--tmerg-code", type=float, default=None)
    ap.add_argument("--xmin", type=float, default=-120.0)
    ap.add_argument("--xmax", type=float, default=120.0)
    ap.add_argument("--ymin", type=float, default=0.0)
    ap.add_argument("--ymax", type=float, default=280.0)
    ap.add_argument("--itmin", type=int, default=None)
    ap.add_argument("--itmax", type=int, default=None)
    ap.add_argument("--tmin", type=float, default=None)
    ap.add_argument("--tmax", type=float, default=None)
    ap.add_argument("--outdir", default="derived_bfield_slices")
    ap.add_argument("--no-tex", action="store_true")
    ap.add_argument(
        "--cu-to-gauss",
        type=float,
        default=_DEFAULT_CU_TO_GAUSS,
        help=(
            f"B-field code-unit -> Gauss factor "
            f"(default {_DEFAULT_CU_TO_GAUSS:.4e}). "
            f"Use 8.35195e19 for the sqrt(4pi)-rescaled convention."
        ),
    )
    args = ap.parse_args()

    use_tex = not args.no_tex
    CU_TO_GAUSS = args.cu_to_gauss
    tmerg_ms = args.tmerg_code * CU_TO_MS if args.tmerg_code is not None else None
    root = Path(args.root).resolve()
    outdir = Path(args.outdir).resolve()
    ensure_dir(outdir)
    ax0n, ax1n = _axis_names(args.plane)

    print(f"Building 2D variable index for plane={args.plane} …")
    store = VariableStore2D(root, args.plane)
    store.build(ALL_VARS)

    missing = store.missing(ALL_VARS)
    if missing:
        print("\nData directories searched:")
        for d in find_data2d_dirs(root):
            print(f"  {d}")
        raise SystemExit(f"\nMissing variables (no matching file found): {missing}")

    common = store.common_iterations(ALL_VARS)
    print(f"  Found {len(common)} iteration(s) with all variables present.")

    iterations = select_iterations(store, args, tmerg_ms)
    print(f"  {len(iterations)} iteration(s) after filters.")

    plane_dir = outdir / args.plane
    plots_dir = plane_dir / "plots"
    ensure_dir(plots_dir)
    tmap = _time_map(store)

    for it in iterations:
        t_ms = tmap.get(it, float("nan"))
        t_plot = t_ms if tmerg_ms is None else t_ms - tmerg_ms

        if tmerg_ms is None:
            title = f"t = {t_ms:.2f} ms" if np.isfinite(t_ms) else f"it={it}"
        else:
            title = (
                f"t - t_merg = {t_plot:.2f} ms" if np.isfinite(t_plot) else f"it={it}"
            )

        print(f"\n[it={it}]  {title}")
        print("  Reading fields (window-restricted) …")

        try:
            X, Y, fields, _ = read_all_vars(
                store,
                it,
                args.rl,
                xmin_km=args.xmin,
                xmax_km=args.xmax,
                ymin_km=args.ymin,
                ymax_km=args.ymax,
            )
        except Exception as e:
            print(f"  [skip] read failed: {e}", file=sys.stderr)
            continue

        print(f"  Canvas size: {X.shape[1]} x {X.shape[0]} pixels")

        print("  Computing comoving b-field …")
        try:
            derived = compute_bpol_btor(fields)
        except Exception as e:
            print(f"  [skip] physics failed: {e}", file=sys.stderr)
            continue

        press_safe = np.where(
            np.abs(fields["press"]) > 0.0, fields["press"].astype(float), np.nan
        )
        beta_inv = fields["smallb2"].astype(float) / (2.0 * press_safe)
        rho_log = safe_log10(fields["rho"].astype(float) * CU_TO_DENS_CGS)

        bundle_path = plane_dir / f"bundle_it{it:08d}.npz"
        np.savez_compressed(
            bundle_path,
            iteration=it,
            time_ms=t_ms,
            time_plot_ms=t_plot,
            X=X.astype(np.float32),
            Y=Y.astype(np.float32),
            b_pol=derived["bpol"].astype(np.float32),
            b_tor=derived["btor"].astype(np.float32),
            b2_comov=derived["b2"].astype(np.float32),
            bx_comov=derived["bx_comov"].astype(np.float32),
            by_comov=derived["by_comov"].astype(np.float32),
            bz_comov=derived["bz_comov"].astype(np.float32),
            plasma_beta_inv=beta_inv.astype(np.float32),
            W=fields["W"].astype(np.float32),
            Ye=fields["Ye"].astype(np.float32),
            rho_log_cgs=rho_log.astype(np.float32),
            smallb2=fields["smallb2"].astype(np.float32),
            press=fields["press"].astype(np.float32),
            plane=np.array(args.plane),
            axis0_name=np.array(ax0n),
            axis1_name=np.array(ax1n),
        )
        print(f"  bundle -> {bundle_path.name}")

        fig_path = plots_dir / f"five_panel_{args.plane}_it{it:08d}.png"
        make_five_panel_plot(
            {
                "X": X,
                "Y": Y,
                "ax0": ax0n,
                "ax1": ax1n,
                "beta_inv": beta_inv,
                "bpol": derived["bpol"],
                "btor": derived["btor"],
                "W": fields["W"],
                "Ye": fields["Ye"],
                "rho_log_cgs": rho_log,
            },
            fig_path,
            title,
            CU_TO_GAUSS,
            use_tex,
        )
        print(f"  [ok]")

        # Explicitly free memory before next iteration
        del fields, derived, beta_inv, rho_log, X, Y

    print(f"\nDone. Outputs in {plane_dir}")


if __name__ == "__main__":
    main()
