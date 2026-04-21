#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

import argparse
import os
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

try:
    import seaborn as sns
except Exception:
    sns = None

from UtilityUnits import CU_to_Gauss, CU_to_densCGS, CU_to_ms

plt.rcParams.update({
    "text.usetex": True,
    "font.size": 12,
    "font.family": "Serif",
    "font.serif": "Computer Modern",
})

Rscale = 1.47662504  # km
DSET_RE = re.compile(r"(\S*)::(\S*) it=(\d+) tl=(\d+) rl=(\d+)(?: c=(\d+))?$")

ALIASES_3D: Dict[str, List[str]] = {
    "rho": ["rho", "rho_b"],
    "W": ["w_lorentz", "W"],
    "Ye": ["Y_e", "Ye", "ye"],
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
}

ALIASES_2D: Dict[str, List[str]] = {
    "W": ["w_lorentz", "W"],
    "Ye": ["Y_e", "Ye", "ye"],
    "rho": ["rho", "rho_b"],
    "smallb2": ["smallb2", "b2", "bsq"],
    "press": ["press", "P"],
}

REQUIRED_3D = ["rho", "Bx", "By", "Bz", "vx", "vy", "vz", "alp", "betax", "betay", "betaz", "gxx", "gxy", "gxz", "gyy", "gyz", "gzz"]
OPTIONAL_3D = ["W"]
REQUIRED_2D = ["W", "Ye", "rho", "smallb2", "press"]


@dataclass
class Tile2D:
    file: str
    dset: str
    group: str
    var: str
    iteration: int
    tl: int
    rl: int
    c: int
    time: float
    origin: np.ndarray
    delta: np.ndarray
    shape_xy: np.ndarray
    corner: np.ndarray


@dataclass
class Tile3D:
    file: str
    dset: str
    group: str
    var: str
    iteration: int
    tl: int
    rl: int
    c: int
    time: float
    origin: np.ndarray
    delta: np.ndarray
    shape_xyz: np.ndarray
    corner: np.ndarray


def natural_key(path_obj: Path | str):
    s = str(path_obj)
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)]


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def safe_log10(arr: np.ndarray, floor: float = 1e-300) -> np.ndarray:
    arr = np.asarray(arr, dtype=float)
    arr = np.where(np.isfinite(arr), arr, np.nan)
    arr = np.where(np.abs(arr) > floor, np.abs(arr), np.nan)
    return np.log10(arr)


def get_rocket_cmap():
    if sns is not None:
        return sns.color_palette("rocket", as_cmap=True)
    return plt.get_cmap("magma")


def choose_plane_axes(plane: str) -> Tuple[int, Tuple[int, int], str, str]:
    if plane == "xy":
        return 2, (0, 1), "x", "y"
    if plane == "xz":
        return 1, (0, 2), "x", "z"
    if plane == "yz":
        return 0, (1, 2), "y", "z"
    raise ValueError(f"Unsupported plane {plane}")


def candidate_output_dirs(root: Path) -> List[Path]:
    root = root.resolve()
    if root.is_dir() and root.name.startswith("output"):
        return [root]
    if root.is_dir() and root.name.startswith("data_hdf5_"):
        parent = root.parent
        if parent.is_dir() and parent.name.startswith("output"):
            return [parent]
    return [p for p in sorted(root.glob("output*"), key=natural_key) if p.is_dir() and p.name.startswith("output")]


def candidate_data_dirs(root: Path, ndim: int) -> List[Path]:
    subnames = [f"data_hdf5_{ndim}D", "data_hdf5", "data_h5", "data"]
    root = root.resolve()
    if root.is_dir() and root.name in subnames:
        return [root]
    if root.is_dir() and root.name.startswith("output"):
        return [d for sub in subnames if (d := root / sub).is_dir()]
    dirs: List[Path] = []
    for out in candidate_output_dirs(root):
        for sub in subnames:
            d = out / sub
            if d.is_dir() and d not in dirs:
                dirs.append(d)
    return dirs


def filename_matches_alias(filename: str, alias: str, plane: Optional[str] = None) -> bool:
    if not filename.endswith(".h5"):
        return False
    if plane is not None:
        return filename.startswith(f"{alias}.{plane}")
    return filename.startswith(f"{alias}.xyz") or filename == f"{alias}.h5"


class VariableStore2D:
    def __init__(self, root: Path, plane: str):
        self.root = root
        self.plane = plane
        self.index: Dict[str, Dict[int, Dict[int, Dict[int, Tile2D]]]] = {}

    def find_variable_files(self, aliases: Sequence[str]) -> List[Path]:
        found: List[Path] = []
        seen = set()
        for d in candidate_data_dirs(self.root, 2):
            try:
                names = os.listdir(d)
            except Exception:
                continue
            for alias in aliases:
                for name in names:
                    if filename_matches_alias(name, alias, plane=self.plane):
                        p = d / name
                        if p not in seen:
                            found.append(p)
                            seen.add(p)
        return sorted(found, key=natural_key)

    @staticmethod
    def read_file_tiles(file_path: Path) -> List[Tile2D]:
        tiles: List[Tile2D] = []
        with h5py.File(file_path, "r") as f:
            for key in f.keys():
                m = DSET_RE.match(key)
                if not m:
                    continue
                group, var, it, tl, rl, c = m.groups()
                ds = f[key]
                if ds.ndim != 2:
                    continue
                origin = np.asarray(ds.attrs["origin"], dtype=float)
                delta = np.asarray(ds.attrs["delta"], dtype=float)
                shape_xy = np.asarray(ds.shape[::-1], dtype=int)
                corner = origin + delta * (shape_xy - 1)
                time = float(np.real(ds.attrs.get("time", np.nan)))
                tiles.append(Tile2D(file=file_path.as_posix(), dset=key, group=group, var=var,
                                    iteration=int(it), tl=int(tl), rl=int(rl), c=int(c) if c is not None else -1,
                                    time=time, origin=origin, delta=delta, shape_xy=shape_xy, corner=corner))
        return tiles

    def build(self, names: Sequence[str]) -> None:
        for name in names:
            files = self.find_variable_files(ALIASES_2D[name])
            if not files:
                continue
            idx: Dict[int, Dict[int, Dict[int, Tile2D]]] = defaultdict(lambda: defaultdict(dict))
            for fp in files:
                for tile in self.read_file_tiles(fp):
                    if tile.tl == 0:
                        idx[tile.iteration][tile.rl][tile.c] = tile
            if idx:
                self.index[name] = idx

    def missing(self, names: Sequence[str]) -> List[str]:
        return [n for n in names if n not in self.index]

    def common_iterations(self, names: Sequence[str]) -> List[int]:
        sets = [set(self.index[n].keys()) for n in names if n in self.index]
        return sorted(set.intersection(*sets)) if sets else []


class VariableStore3D:
    def __init__(self, root: Path):
        self.root = root
        self.index: Dict[str, Dict[int, Dict[int, Dict[int, Tile3D]]]] = {}
        self.files_used: Dict[str, List[str]] = {}

    def find_variable_files(self, aliases: Sequence[str]) -> List[Path]:
        found: List[Path] = []
        seen = set()
        for d in candidate_data_dirs(self.root, 3):
            try:
                names = os.listdir(d)
            except Exception:
                continue
            for alias in aliases:
                for name in names:
                    if filename_matches_alias(name, alias):
                        p = d / name
                        if p not in seen:
                            found.append(p)
                            seen.add(p)
        return sorted(found, key=natural_key)

    @staticmethod
    def read_file_tiles(file_path: Path) -> List[Tile3D]:
        tiles: List[Tile3D] = []
        with h5py.File(file_path, "r") as f:
            for key in f.keys():
                m = DSET_RE.match(key)
                if not m:
                    continue
                group, var, it, tl, rl, c = m.groups()
                ds = f[key]
                if ds.ndim != 3:
                    continue
                origin = np.asarray(ds.attrs["origin"], dtype=float)
                delta = np.asarray(ds.attrs["delta"], dtype=float)
                shape_xyz = np.asarray(ds.shape[::-1], dtype=int)
                corner = origin + delta * (shape_xyz - 1)
                time = float(np.real(ds.attrs.get("time", np.nan)))
                tiles.append(Tile3D(file=file_path.as_posix(), dset=key, group=group, var=var,
                                    iteration=int(it), tl=int(tl), rl=int(rl), c=int(c) if c is not None else -1,
                                    time=time, origin=origin, delta=delta, shape_xyz=shape_xyz, corner=corner))
        return tiles

    def build(self, names: Sequence[str]) -> None:
        for name in names:
            files = self.find_variable_files(ALIASES_3D[name])
            if not files:
                continue
            idx: Dict[int, Dict[int, Dict[int, Tile3D]]] = defaultdict(lambda: defaultdict(dict))
            used: List[str] = []
            for fp in files:
                tiles = self.read_file_tiles(fp)
                if not tiles:
                    continue
                used.append(fp.as_posix())
                for tile in tiles:
                    if tile.tl == 0:
                        idx[tile.iteration][tile.rl][tile.c] = tile
            if idx:
                self.index[name] = idx
                self.files_used[name] = used

    def missing(self, names: Sequence[str]) -> List[str]:
        return [n for n in names if n not in self.index]

    def common_iterations(self, names: Sequence[str]) -> List[int]:
        sets = [set(self.index[n].keys()) for n in names if n in self.index]
        return sorted(set.intersection(*sets)) if sets else []

    def common_rls(self, names: Sequence[str], iteration: int) -> List[int]:
        sets = []
        for name in names:
            if name not in self.index or iteration not in self.index[name]:
                return []
            sets.append(set(self.index[name][iteration].keys()))
        return sorted(set.intersection(*sets)) if sets else []


def collect_tiles_2d(store: VariableStore2D, name: str, iteration: int, requested_rl: Optional[int]) -> List[Tile2D]:
    return [store.index[name][iteration][rl][c]
            for rl in sorted(store.index[name][iteration].keys())
            if requested_rl is None or rl == requested_rl
            for c in sorted(store.index[name][iteration][rl].keys())]


def collect_tiles_3d(store: VariableStore3D, name: str, iteration: int, requested_rl: Optional[int]) -> List[Tile3D]:
    return [store.index[name][iteration][rl][c]
            for rl in sorted(store.index[name][iteration].keys())
            if requested_rl is None or rl == requested_rl
            for c in sorted(store.index[name][iteration][rl].keys())]


def tile_bounds_km_2d(tile: Tile2D) -> Tuple[float, float, float, float]:
    return (float(min(tile.origin[0], tile.corner[0]) * Rscale),
            float(max(tile.origin[0], tile.corner[0]) * Rscale),
            float(min(tile.origin[1], tile.corner[1]) * Rscale),
            float(max(tile.origin[1], tile.corner[1]) * Rscale))


def tile_bounds_km_3d(tile: Tile3D, plane: str) -> Tuple[float, float, float, float]:
    _, axes_keep, _, _ = choose_plane_axes(plane)
    ia, ib = axes_keep
    return (float(min(tile.origin[ia], tile.corner[ia]) * Rscale),
            float(max(tile.origin[ia], tile.corner[ia]) * Rscale),
            float(min(tile.origin[ib], tile.corner[ib]) * Rscale),
            float(max(tile.origin[ib], tile.corner[ib]) * Rscale))


def build_uniform_mesh(xmin: float, xmax: float, ymin: float, ymax: float, dx: float, dy: float):
    nx = int(np.floor((xmax - xmin) / dx + 0.5)) + 1
    ny = int(np.floor((ymax - ymin) / dy + 0.5)) + 1
    x = xmin + dx * np.arange(nx, dtype=np.float32)
    y = ymin + dy * np.arange(ny, dtype=np.float32)
    return np.meshgrid(x, y)


def compute_target_spacing_km(tile_deltas_code: List[Tuple[float, float]], requested_dx: float, requested_dy: Optional[float],
                              max_points: int, xmin: float, xmax: float, ymin: float, ymax: float) -> Tuple[float, float]:
    min_dx = min(abs(dx_) for dx_, _ in tile_deltas_code) * Rscale
    min_dy = min(abs(dy_) for _, dy_ in tile_deltas_code) * Rscale
    dx = max(requested_dx, min_dx)
    dy = max(requested_dy if requested_dy is not None else requested_dx, min_dy)
    nx = int(np.floor((xmax - xmin) / dx + 0.5)) + 1
    ny = int(np.floor((ymax - ymin) / dy + 0.5)) + 1
    while nx > max_points or ny > max_points:
        dx *= 1.25
        dy *= 1.25
        nx = int(np.floor((xmax - xmin) / dx + 0.5)) + 1
        ny = int(np.floor((ymax - ymin) / dy + 0.5)) + 1
    return float(dx), float(dy)


def compute_global_mesh(store3d: VariableStore3D, store2d: VariableStore2D,
                        it3d: int, it2d: int, plane: str,
                        requested_rl3d: Optional[int], requested_rl2d: Optional[int],
                        save_dx_km: float, save_dy_km: Optional[float], max_points: int):
    bounds: List[Tuple[float, float, float, float]] = []
    deltas: List[Tuple[float, float]] = []
    _, axes_keep, _, _ = choose_plane_axes(plane)
    for name in REQUIRED_3D + [n for n in OPTIONAL_3D if n in store3d.index]:
        if name not in store3d.index or it3d not in store3d.index[name]:
            continue
        for t in collect_tiles_3d(store3d, name, it3d, requested_rl3d):
            bounds.append(tile_bounds_km_3d(t, plane))
            deltas.append((float(t.delta[axes_keep[0]]), float(t.delta[axes_keep[1]])))
    for name in REQUIRED_2D:
        if name not in store2d.index or it2d not in store2d.index[name]:
            continue
        for t in collect_tiles_2d(store2d, name, it2d, requested_rl2d):
            bounds.append(tile_bounds_km_2d(t))
            deltas.append((float(t.delta[0]), float(t.delta[1])))
    if not bounds:
        raise ValueError("No bounds available for common mesh.")
    xmin = min(b[0] for b in bounds); xmax = max(b[1] for b in bounds)
    ymin = min(b[2] for b in bounds); ymax = max(b[3] for b in bounds)
    dx, dy = compute_target_spacing_km(deltas, save_dx_km, save_dy_km, max_points, xmin, xmax, ymin, ymax)
    X, Y = build_uniform_mesh(xmin, xmax, ymin, ymax, dx, dy)
    return X, Y, dx, dy


def centers_to_edges(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    if v.size == 1:
        return np.array([v[0] - 0.5, v[0] + 0.5], dtype=float)
    mids = 0.5 * (v[:-1] + v[1:])
    left = v[0] - 0.5 * (v[1] - v[0])
    right = v[-1] + 0.5 * (v[-1] - v[-2])
    return np.concatenate(([left], mids, [right]))


def paste_regular_tile(canvas: np.ndarray, x_centers: np.ndarray, y_centers: np.ndarray,
                       x_tile: np.ndarray, y_tile: np.ndarray, vals: np.ndarray) -> None:
    x_canvas_edges = centers_to_edges(x_centers)
    y_canvas_edges = centers_to_edges(y_centers)
    x_tile_edges = centers_to_edges(x_tile)
    y_tile_edges = centers_to_edges(y_tile)
    ny_loc, nx_loc = vals.shape
    for j in range(ny_loc):
        ylo = min(y_tile_edges[j], y_tile_edges[j + 1])
        yhi = max(y_tile_edges[j], y_tile_edges[j + 1])
        gy0 = max(0, np.searchsorted(y_canvas_edges, ylo, side="right") - 1)
        gy1 = min(canvas.shape[0], np.searchsorted(y_canvas_edges, yhi, side="left"))
        if gy1 <= gy0:
            continue
        for i in range(nx_loc):
            v = vals[j, i]
            if not np.isfinite(v):
                continue
            xlo = min(x_tile_edges[i], x_tile_edges[i + 1])
            xhi = max(x_tile_edges[i], x_tile_edges[i + 1])
            gx0 = max(0, np.searchsorted(x_canvas_edges, xlo, side="right") - 1)
            gx1 = min(canvas.shape[1], np.searchsorted(x_canvas_edges, xhi, side="left"))
            if gx1 <= gx0:
                continue
            canvas[gy0:gy1, gx0:gx1] = v


def read_dataset_values(tile) -> np.ndarray:
    with h5py.File(tile.file, "r") as f:
        return np.asarray(f[tile.dset], dtype=np.float32)


ODD_Z_VARS = {"vel[2]", "vz", "betaz", "gxz", "gyz", "Bvec[0]", "Bvec[1]", "Bx", "By"}


def extract_2d_values_from_3d(tile: Tile3D, values3d: np.ndarray, plane: str, coord0: float, mirror_z: bool):
    x = tile.origin[0] + tile.delta[0] * np.arange(values3d.shape[2], dtype=np.float32)
    y = tile.origin[1] + tile.delta[1] * np.arange(values3d.shape[1], dtype=np.float32)
    z = tile.origin[2] + tile.delta[2] * np.arange(values3d.shape[0], dtype=np.float32)

    if mirror_z and plane in ("xz", "yz"):
        zmin, zmax = np.nanmin(z), np.nanmax(z)
        if zmin >= -1e-12 and zmax > 0 and values3d.shape[0] > 1:
            vals_neg = values3d[1:][::-1].copy()
            z_neg = -z[1:][::-1]
            if tile.var in ODD_Z_VARS:
                vals_neg *= -1.0
            values3d = np.concatenate([vals_neg, values3d], axis=0)
            z = np.concatenate([z_neg, z], axis=0)

    dir_idx, _, _, _ = choose_plane_axes(plane)
    coords = [x, y, z]
    idx = int(np.argmin(np.abs(coords[dir_idx] - coord0)))
    if plane == "xy":
        vals2d = values3d[idx, :, :]
        Xt, Yt = np.meshgrid(x * Rscale, y * Rscale)
    elif plane == "xz":
        vals2d = values3d[:, idx, :]
        Xt, Yt = np.meshgrid(x * Rscale, z * Rscale)
    else:
        vals2d = values3d[:, :, idx]
        Xt, Yt = np.meshgrid(y * Rscale, z * Rscale)
    return Xt.astype(np.float32), Yt.astype(np.float32), np.asarray(vals2d, dtype=np.float32)


def join_iteration_2d_to_mesh(store: VariableStore2D, name: str, iteration: int, X: np.ndarray, Y: np.ndarray,
                              requested_rl: Optional[int]) -> Tuple[np.ndarray, float]:
    tiles = collect_tiles_2d(store, name, iteration, requested_rl)
    if not tiles:
        raise ValueError(f"No 2D tiles found for {name} at iteration {iteration}")
    canvas = np.full(X.shape, np.nan, dtype=np.float32)
    x_centers = X[0, :]
    y_centers = Y[:, 0]
    for tile in sorted(tiles, key=lambda t: (t.rl, t.c, t.file)):
        vals = read_dataset_values(tile)
        x0 = tile.origin[0] * Rscale
        y0 = tile.origin[1] * Rscale
        x_tile = x0 + tile.delta[0] * Rscale * np.arange(vals.shape[1], dtype=np.float32)
        y_tile = y0 + tile.delta[1] * Rscale * np.arange(vals.shape[0], dtype=np.float32)
        paste_regular_tile(canvas, x_centers, y_centers, x_tile, y_tile, vals)
    return canvas, float(tiles[0].time)


def join_iteration_3d_to_mesh(store: VariableStore3D, name: str, iteration: int, plane: str, coord0: float,
                              X: np.ndarray, Y: np.ndarray, requested_rl: Optional[int], mirror_z: bool) -> Tuple[np.ndarray, float]:
    tiles = collect_tiles_3d(store, name, iteration, requested_rl)
    if not tiles:
        raise ValueError(f"No 3D tiles found for {name} at iteration {iteration}")
    canvas = np.full(X.shape, np.nan, dtype=np.float32)
    x_centers = X[0, :]
    y_centers = Y[:, 0]
    for tile in sorted(tiles, key=lambda t: (t.rl, t.c, t.file)):
        vals3d = read_dataset_values(tile)
        Xt, Yt, vals2d = extract_2d_values_from_3d(tile, vals3d, plane, coord0, mirror_z)
        paste_regular_tile(canvas, x_centers, y_centers, Xt[0, :], Yt[:, 0], vals2d)
    return canvas, float(tiles[0].time)


def spatial_dot(ax, ay, az, bx, by, bz, gxx, gxy, gxz, gyy, gyz, gzz):
    return (gxx * ax * bx + gyy * ay * by + gzz * az * bz
            + gxy * (ax * by + ay * bx) + gxz * (ax * bz + az * bx) + gyz * (ay * bz + az * by))


def maybe_compute_W(data: Dict[str, np.ndarray]) -> np.ndarray:
    if "W" in data:
        W = np.asarray(data["W"], dtype=np.float32)
        if np.any(np.isfinite(W)):
            return W
    vsq = spatial_dot(data["vx"], data["vy"], data["vz"], data["vx"], data["vy"], data["vz"],
                      data["gxx"], data["gxy"], data["gxz"], data["gyy"], data["gyz"], data["gzz"])
    vsq = np.clip(vsq, 0.0, 1.0 - 1e-12)
    return (1.0 / np.sqrt(1.0 - vsq)).astype(np.float32)


def compute_comoving_b(data: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    Bx, By, Bz = data["Bx"], data["By"], data["Bz"]
    vx, vy, vz = data["vx"], data["vy"], data["vz"]
    alp = data["alp"]
    betax, betay, betaz = data["betax"], data["betay"], data["betaz"]
    W = maybe_compute_W(data)
    alp_safe = np.where(np.abs(alp) > 1e-14, alp, np.nan)
    vdotB = spatial_dot(Bx, By, Bz, vx, vy, vz, data["gxx"], data["gxy"], data["gxz"], data["gyy"], data["gyz"], data["gzz"])
    b0 = W * vdotB / alp_safe
    ux = W * (vx - betax / alp_safe)
    uy = W * (vy - betay / alp_safe)
    uz = W * (vz - betaz / alp_safe)
    bx = (Bx + alp * b0 * ux) / W
    by = (By + alp * b0 * uy) / W
    bz = (Bz + alp * b0 * uz) / W
    b2 = spatial_dot(bx, by, bz, bx, by, bz, data["gxx"], data["gxy"], data["gxz"], data["gyy"], data["gyz"], data["gzz"]) - (alp ** 2) * (b0 ** 2)
    return {"b0": b0.astype(np.float32), "bx": bx.astype(np.float32), "by": by.astype(np.float32), "bz": bz.astype(np.float32),
            "b2": b2.astype(np.float32), "W": W.astype(np.float32)}


def normalize_basis(rx, ry, rz, gxx, gxy, gxz, gyy, gyz, gzz):
    norm = np.sqrt(np.maximum(spatial_dot(rx, ry, rz, rx, ry, rz, gxx, gxy, gxz, gyy, gyz, gzz), 0.0))
    norm = np.where(norm > 1e-14, norm, np.nan)
    return rx / norm, ry / norm, rz / norm


def panel_coords_to_cartesian(X, Y, plane: str, coord0: float):
    Xc = np.asarray(X, dtype=np.float32) / Rscale
    Yc = np.asarray(Y, dtype=np.float32) / Rscale
    c0 = np.float32(coord0)
    if plane == "xy":
        return Xc, Yc, np.full_like(Xc, c0)
    if plane == "xz":
        return Xc, np.full_like(Xc, c0), Yc
    return np.full_like(Xc, c0), Xc, Yc


def orthonormal_cylindrical_basis(x, y, z, gxx, gxy, gxz, gyy, gyz, gzz):
    raw_Rx = np.asarray(x, dtype=np.float32)
    raw_Ry = np.asarray(y, dtype=np.float32)
    raw_Rz = np.zeros_like(raw_Rx)
    r_cyl = np.sqrt(raw_Rx ** 2 + raw_Ry ** 2)
    mask_axis = r_cyl <= 1e-14
    eRx, eRy, eRz = normalize_basis(raw_Rx, raw_Ry, raw_Rz, gxx, gxy, gxz, gyy, gyz, gzz)
    raw_phix = -raw_Ry
    raw_phiy = raw_Rx
    raw_phiz = np.zeros_like(raw_phix)
    ephix, ephiy, ephiz = normalize_basis(raw_phix, raw_phiy, raw_phiz, gxx, gxy, gxz, gyy, gyz, gzz)
    zzx = np.zeros_like(raw_Rx)
    zzy = np.zeros_like(raw_Rx)
    zzz = np.ones_like(raw_Rx)
    z_dot_eR = spatial_dot(zzx, zzy, zzz, eRx, eRy, eRz, gxx, gxy, gxz, gyy, gyz, gzz)
    z_dot_ephi = spatial_dot(zzx, zzy, zzz, ephix, ephiy, ephiz, gxx, gxy, gxz, gyy, gyz, gzz)
    raw_ezx = zzx - z_dot_eR * eRx - z_dot_ephi * ephix
    raw_ezy = zzy - z_dot_eR * eRy - z_dot_ephi * ephiy
    raw_ezz = zzz - z_dot_eR * eRz - z_dot_ephi * ephiz
    ezx, ezy, ezz = normalize_basis(raw_ezx, raw_ezy, raw_ezz, gxx, gxy, gxz, gyy, gyz, gzz)
    for arr in (eRx, eRy, eRz, ephix, ephiy, ephiz):
        arr[mask_axis] = np.nan
    return {"eR": (eRx, eRy, eRz), "ephi": (ephix, ephiy, ephiz), "ez": (ezx, ezy, ezz)}


def compute_derived_slice_fields(slice_data: Dict[str, np.ndarray], plane: str, coord0: float) -> Dict[str, np.ndarray]:
    comov = compute_comoving_b(slice_data)
    x, y, z = panel_coords_to_cartesian(slice_data["X"], slice_data["Y"], plane, coord0)
    basis = orthonormal_cylindrical_basis(x, y, z, slice_data["gxx"], slice_data["gxy"], slice_data["gxz"],
                                          slice_data["gyy"], slice_data["gyz"], slice_data["gzz"])
    bR = spatial_dot(comov["bx"], comov["by"], comov["bz"], *basis["eR"],
                     slice_data["gxx"], slice_data["gxy"], slice_data["gxz"], slice_data["gyy"], slice_data["gyz"], slice_data["gzz"])
    bphi = spatial_dot(comov["bx"], comov["by"], comov["bz"], *basis["ephi"],
                       slice_data["gxx"], slice_data["gxy"], slice_data["gxz"], slice_data["gyy"], slice_data["gyz"], slice_data["gzz"])
    bzproj = spatial_dot(comov["bx"], comov["by"], comov["bz"], *basis["ez"],
                         slice_data["gxx"], slice_data["gxy"], slice_data["gxz"], slice_data["gyy"], slice_data["gyz"], slice_data["gzz"])
    bpol = np.sqrt(np.maximum(bR ** 2 + bzproj ** 2, 0.0)).astype(np.float32)
    return {"bx_comov": comov["bx"], "by_comov": comov["by"], "bz_comov": comov["bz"],
            "b2": comov["b2"], "btor": np.abs(bphi).astype(np.float32), "bpol": bpol, "W3D": comov["W"]}


def top_cbar(ax, cf, label, ticks=None, fmt=r"$%.2f$"):
    axins = inset_axes(ax, width="80%", height="3.5%", loc="lower left",
                       bbox_to_anchor=(0.125, 1.01, 1, 1), bbox_transform=ax.transAxes, borderpad=0)
    cb = plt.colorbar(cf, cax=axins, orientation="horizontal")
    cb.set_label(label, fontsize=18)
    cb.ax.xaxis.set_ticks_position("top")
    cb.ax.xaxis.set_label_position("top")
    cb.ax.tick_params(rotation=45)
    if ticks is not None:
        cb.set_ticks(ticks)
    cb.ax.xaxis.set_major_formatter(FormatStrFormatter(fmt))
    return cb


def overlay_rho_contours(axs, X, Y, rho_log_cgs):
    levels = [np.log10(1e10), 12, 13, 14]
    for ax in axs:
        ax.contour(X, Y, rho_log_cgs, levels=levels, colors=["deepskyblue"], linewidths=[1.2], linestyles=[":", "-.", "--", "-"])


def choose_ranges(fields: Dict[str, np.ndarray]) -> Dict[str, Tuple[float, float]]:
    def robust(arr, default_min, default_max, qlo=1.0, qhi=99.0):
        arr = np.asarray(arr, dtype=float)
        arr = arr[np.isfinite(arr)]
        if arr.size == 0:
            return default_min, default_max
        lo = float(np.nanpercentile(arr, qlo))
        hi = float(np.nanpercentile(arr, qhi))
        if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
            return default_min, default_max
        return lo, hi
    beta_min, beta_max = robust(fields["log10_beta_inv"], -2.0, 2.0)
    btor_min, btor_max = robust(fields["log10_btor"], 12.0, 16.0)
    bpol_min, bpol_max = robust(fields["log10_bpol"], 12.0, 16.0)
    beta_min = min(-2.0, np.floor(beta_min)); beta_max = max(2.0, np.ceil(beta_max))
    btor_min = min(12.0, np.floor(btor_min)); btor_max = max(16.0, np.ceil(btor_max))
    bpol_min = min(12.0, np.floor(bpol_min)); bpol_max = max(16.0, np.ceil(bpol_max))
    W = np.asarray(fields["W"]); Ye = np.asarray(fields["Ye"])
    wmax = max(1.1, float(np.nanpercentile(W[np.isfinite(W)], 99.5)) if np.any(np.isfinite(W)) else 1.1)
    yemax = min(0.6, max(0.5, float(np.nanpercentile(Ye[np.isfinite(Ye)], 99.5)) if np.any(np.isfinite(Ye)) else 0.5))
    return {"beta_inv": (beta_min, beta_max), "btor": (btor_min, btor_max), "bpol": (bpol_min, bpol_max), "W": (1.0, wmax), "Ye": (0.0, yemax)}


def make_five_panel_plot(slice_all: Dict[str, np.ndarray], outpath: Path, xmin: float, xmax: float, ymin: float, ymax: float):
    X = slice_all["X"]; Y = slice_all["Y"]
    fields = {
        "log10_beta_inv": safe_log10(slice_all["beta_inv"]),
        "log10_btor": safe_log10(slice_all["btor"] * CU_to_Gauss),
        "log10_bpol": safe_log10(slice_all["bpol"] * CU_to_Gauss),
        "W": slice_all["W"], "Ye": slice_all["Ye"],
    }
    ranges = choose_ranges(fields)
    fig = plt.figure(figsize=(30, 12))
    ax1 = fig.add_subplot(151); ax2 = fig.add_subplot(152, sharey=ax1); ax3 = fig.add_subplot(153, sharey=ax2)
    ax4 = fig.add_subplot(154, sharey=ax3); ax5 = fig.add_subplot(155, sharey=ax4)
    cf1 = ax1.contourf(X, Y, fields["log10_beta_inv"], levels=np.linspace(*ranges["beta_inv"], 100), cmap="seismic", extend="both")
    cf2 = ax2.contourf(X, Y, fields["log10_bpol"], levels=np.linspace(*ranges["bpol"], 100), cmap="plasma", extend="both")
    cf3 = ax3.contourf(X, Y, fields["log10_btor"], levels=np.linspace(*ranges["btor"], 100), cmap=get_rocket_cmap(), extend="both")
    cf4 = ax4.contourf(X, Y, fields["W"], levels=np.linspace(*ranges["W"], 100), cmap="jet", extend="max")
    cf5 = ax5.contourf(X, Y, fields["Ye"], levels=np.linspace(*ranges["Ye"], 100), cmap="PRGn", extend="both")
    top_cbar(ax1, cf1, r"$\log_{10}(\beta^{-1})$", ticks=np.linspace(*ranges["beta_inv"], 9))
    top_cbar(ax2, cf2, r"$\log_{10}(b_{\rm pol}[{\rm G}])$", ticks=np.linspace(*ranges["bpol"], 9))
    top_cbar(ax3, cf3, r"$\log_{10}(b_{\rm tor}[{\rm G}])$", ticks=np.linspace(*ranges["btor"], 9))
    top_cbar(ax4, cf4, r"$W$", ticks=np.linspace(*ranges["W"], 7))
    top_cbar(ax5, cf5, r"$Y_e$", ticks=np.linspace(*ranges["Ye"], 7))
    overlay_rho_contours([ax1, ax2, ax3, ax4, ax5], X, Y, slice_all["rho_log_cgs"])
    y_label = str(slice_all["axis1_name"]); x_label = str(slice_all["axis0_name"])
    ax1.set_ylabel(rf"${y_label}~[{{\rm km}}]$", fontsize=24)
    for ax in [ax1, ax2, ax3, ax4, ax5]:
        ax.set_xlabel(rf"${x_label}~[{{\rm km}}]$", fontsize=24)
        ax.tick_params(which="both", top=True, right=True, direction="in", length=6, color="white", labelsize=18)
        ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
    for ax in [ax2, ax3, ax4]:
        ax.tick_params(labelleft=False)
    ax5.tick_params(labelleft=False, labelright=True)
    plt.subplots_adjust(wspace=0.02)
    fig.savefig(outpath, bbox_inches="tight", dpi=180)
    plt.close(fig)


def crop_mesh_and_fields(X: np.ndarray, Y: np.ndarray, fields: Dict[str, np.ndarray], xmin: float, xmax: float, ymin: float, ymax: float):
    x = X[0, :]
    y = Y[:, 0]
    ix = np.where((x >= xmin) & (x <= xmax))[0]
    iy = np.where((y >= ymin) & (y <= ymax))[0]
    if ix.size == 0 or iy.size == 0:
        raise ValueError("Requested plotting range does not overlap available data.")
    Xc = X[np.ix_(iy, ix)]
    Yc = Y[np.ix_(iy, ix)]
    out = {k: np.asarray(v)[np.ix_(iy, ix)] for k, v in fields.items()}
    return Xc, Yc, out


def write_quantity_npz(outdir: Path, quantity: str, iteration: int, time_ms: float, time_plot_ms: float,
                       X: np.ndarray, Y: np.ndarray, values: np.ndarray, axis0_name: str, axis1_name: str) -> None:
    ensure_dir(outdir)
    np.savez_compressed(outdir / f"{quantity}_it{iteration:08d}.npz",
                        iteration=int(iteration), time_ms=float(time_ms), time_plot_ms=float(time_plot_ms),
                        X=X.astype(np.float32), Y=Y.astype(np.float32), values=np.asarray(values, dtype=np.float32),
                        quantity=np.array(quantity), axis0_name=np.array(axis0_name), axis1_name=np.array(axis1_name))


def _time_map_3d_ms(store3d: VariableStore3D, names: Sequence[str]) -> Dict[int, float]:
    out: Dict[int, float] = {}
    for name in names:
        if name not in store3d.index:
            continue
        for it, rl_map in store3d.index[name].items():
            if it in out:
                continue
            found = None
            for rl in sorted(rl_map.keys()):
                for c in sorted(rl_map[rl].keys()):
                    found = rl_map[rl][c]; break
                if found is not None:
                    break
            if found is not None and np.isfinite(found.time):
                out[int(it)] = float(found.time) * CU_to_ms
    return out


def _time_map_2d_ms(store2d: VariableStore2D) -> Dict[int, float]:
    out: Dict[int, float] = {}
    if "rho" not in store2d.index:
        return out
    for it, rl_map in store2d.index["rho"].items():
        found = None
        for rl in sorted(rl_map.keys()):
            for c in sorted(rl_map[rl].keys()):
                found = rl_map[rl][c]; break
            if found is not None:
                break
        if found is not None and np.isfinite(found.time):
            out[int(it)] = float(found.time) * CU_to_ms
    return out


def determine_iteration_pairs(store3d: VariableStore3D, store2d: VariableStore2D, args):
    common_3d = store3d.common_iterations(REQUIRED_3D)
    if not common_3d:
        raise SystemExit("No common 3D iterations found across required variables.")
    missing_2d = store2d.missing(REQUIRED_2D)
    if missing_2d:
        raise SystemExit(f"Missing required 2D variables: {missing_2d}")
    common_2d = store2d.common_iterations(REQUIRED_2D)
    if not common_2d:
        raise SystemExit("No common 2D iterations found across required variables.")
    tmerg_ms = args.tmerg_code * CU_to_ms if args.tmerg_code is not None else None
    time3d_ms = _time_map_3d_ms(store3d, REQUIRED_3D)
    time2d_ms = _time_map_2d_ms(store2d)

    def pass_filters(it: int, time_map: Dict[int, float]) -> bool:
        if args.itmin is not None and it < args.itmin:
            return False
        if args.itmax is not None and it > args.itmax:
            return False
        t_ms = time_map.get(int(it))
        if args.tmin is not None or args.tmax is not None:
            if t_ms is None:
                return False
            t_plot = t_ms if tmerg_ms is None else (t_ms - tmerg_ms)
            if args.tmin is not None and t_plot < args.tmin:
                return False
            if args.tmax is not None and t_plot > args.tmax:
                return False
        return True

    eligible_3d = [int(it) for it in common_3d if pass_filters(int(it), time3d_ms)]
    eligible_2d = [int(it) for it in common_2d if pass_filters(int(it), time2d_ms)]
    if not eligible_3d:
        raise SystemExit("No 3D iterations left after filters.")
    if not eligible_2d:
        raise SystemExit("No 2D iterations left after filters.")

    exact = sorted(set(eligible_3d).intersection(eligible_2d))
    if exact:
        return [(it, it) for it in exact]

    pairs = []
    for it3d in eligible_3d:
        if it3d in time3d_ms and time2d_ms:
            t3 = time3d_ms[it3d]
            best = min(eligible_2d, key=lambda it2d: (abs(time2d_ms.get(it2d, np.inf) - t3), abs(it2d - it3d), it2d))
        else:
            best = min(eligible_2d, key=lambda it2d: (abs(it2d - it3d), it2d))
        pairs.append((it3d, best))
    return pairs


def choose_rl_for_iteration(store3d: VariableStore3D, iteration: int, requested_rl: Optional[int]) -> Optional[int]:
    common_rls = store3d.common_rls(REQUIRED_3D, iteration)
    if not common_rls:
        raise SystemExit(f"No common 3D refinement levels for iteration {iteration}.")
    if requested_rl is not None and requested_rl not in common_rls:
        raise SystemExit(f"Requested 3D rl={requested_rl} not available for iteration {iteration}. Common RLs: {common_rls}")
    return requested_rl


def main():
    ap = argparse.ArgumentParser(description="Compute full-plane b_pol and b_tor on one common mesh, then crop only for plotting.")
    ap.add_argument("--root", required=True)
    ap.add_argument("--plane", default="xz", choices=["xy", "xz", "yz"])
    ap.add_argument("--coord0", type=float, default=0.0)
    ap.add_argument("--rl3d", type=int, default=None)
    ap.add_argument("--rl2d", type=int, default=None)
    ap.add_argument("--mirror-z", action="store_true")
    ap.add_argument("--tmerg-code", type=float, default=None)
    ap.add_argument("--xmin", type=float, default=-120.0)
    ap.add_argument("--xmax", type=float, default=120.0)
    ap.add_argument("--ymin", type=float, default=0.0)
    ap.add_argument("--ymax", type=float, default=300.0)
    ap.add_argument("--itmin", type=int, default=None)
    ap.add_argument("--itmax", type=int, default=None)
    ap.add_argument("--tmin", type=float, default=None)
    ap.add_argument("--tmax", type=float, default=None)
    ap.add_argument("--outdir", default="derived_bfield_panels")
    ap.add_argument("--save-dx-km", type=float, default=1.0)
    ap.add_argument("--save-dy-km", type=float, default=None)
    ap.add_argument("--max-save-points", type=int, default=4096)
    args = ap.parse_args()

    root = Path(args.root).resolve()
    outdir = Path(args.outdir).resolve()
    ensure_dir(outdir)

    store3d = VariableStore3D(root)
    store3d.build(REQUIRED_3D + OPTIONAL_3D)
    missing_3d = store3d.missing(REQUIRED_3D)
    if missing_3d:
        raise SystemExit(f"Missing required 3D variables: {missing_3d}")
    store2d = VariableStore2D(root, args.plane)
    store2d.build(REQUIRED_2D)
    pairs = determine_iteration_pairs(store3d, store2d, args)

    plane_dir = outdir / args.plane
    plots_dir = plane_dir / "plots"
    ensure_dir(plots_dir)
    _, _, axis0_name, axis1_name = choose_plane_axes(args.plane)
    tmerg_ms = args.tmerg_code * CU_to_ms if args.tmerg_code is not None else None

    for it3d, it2d in pairs:
        rl3d = choose_rl_for_iteration(store3d, int(it3d), args.rl3d)
        X, Y, dx_used, dy_used = compute_global_mesh(store3d, store2d, int(it3d), int(it2d), args.plane,
                                                     rl3d, args.rl2d, args.save_dx_km, args.save_dy_km, args.max_save_points)

        slice3d: Dict[str, np.ndarray] = {"X": X, "Y": Y}
        time3d_code = np.nan
        for name in REQUIRED_3D + [n for n in OPTIONAL_3D if n in store3d.index]:
            arr, tcode = join_iteration_3d_to_mesh(store3d, name, int(it3d), args.plane, args.coord0, X, Y, rl3d, args.mirror_z)
            slice3d[name] = arr
            if not np.isfinite(time3d_code):
                time3d_code = tcode

        derived = compute_derived_slice_fields(slice3d, args.plane, args.coord0)

        W2d, t2d_code = join_iteration_2d_to_mesh(store2d, "W", int(it2d), X, Y, args.rl2d)
        Ye2d, _ = join_iteration_2d_to_mesh(store2d, "Ye", int(it2d), X, Y, args.rl2d)
        rho2d, _ = join_iteration_2d_to_mesh(store2d, "rho", int(it2d), X, Y, args.rl2d)
        smallb2_2d, _ = join_iteration_2d_to_mesh(store2d, "smallb2", int(it2d), X, Y, args.rl2d)
        press2d, _ = join_iteration_2d_to_mesh(store2d, "press", int(it2d), X, Y, args.rl2d)

        beta_inv = smallb2_2d / (2.0 * np.where(np.abs(press2d) > 0.0, press2d, np.nan))
        rho_log_cgs = safe_log10(rho2d * CU_to_densCGS).astype(np.float32)
        time_ms = float(t2d_code) * CU_to_ms
        time_plot_ms = time_ms if tmerg_ms is None else (time_ms - tmerg_ms)

        for qname, qval in {"b_pol": derived["bpol"], "b_tor": derived["btor"], "plasma_beta_inv": beta_inv, "W": W2d, "Ye": Ye2d}.items():
            write_quantity_npz(plane_dir / qname, qname, int(it3d), time_ms, time_plot_ms, X, Y, qval, axis0_name, axis1_name)

        np.savez_compressed(
            plane_dir / f"bundle_it3d{int(it3d):08d}_it2d{int(it2d):08d}.npz",
            iteration_3d=int(it3d), iteration_2d=int(it2d),
            time_ms=float(time_ms), time_plot_ms=float(time_plot_ms),
            X=X.astype(np.float32), Y=Y.astype(np.float32),
            b_pol=derived["bpol"], b_tor=derived["btor"], plasma_beta_inv=beta_inv.astype(np.float32),
            W=W2d.astype(np.float32), Ye=Ye2d.astype(np.float32), rho_log_cgs=rho_log_cgs,
            bx_comov=derived["bx_comov"], by_comov=derived["by_comov"], bz_comov=derived["bz_comov"],
            b2_comov_3d=derived["b2"], smallb2_2d=smallb2_2d, press_2d=press2d,
            plane=np.array(args.plane), axis0_name=np.array(axis0_name), axis1_name=np.array(axis1_name),
            rl3d=np.array(-1 if rl3d is None else int(rl3d)), rl2d=np.array(-1 if args.rl2d is None else int(args.rl2d)),
            dx_km=np.array(dx_used), dy_km=np.array(dy_used),
        )

        Xp, Yp, cropped = crop_mesh_and_fields(X, Y,
            {"W": W2d, "Ye": Ye2d, "rho_log_cgs": rho_log_cgs, "bpol": derived["bpol"], "btor": derived["btor"], "beta_inv": beta_inv},
            args.xmin, args.xmax, args.ymin, args.ymax)
        fig_name = plots_dir / f"five_panel_{args.plane}_it3d{int(it3d):08d}_it2d{int(it2d):08d}.png"
        make_five_panel_plot({"X": Xp, "Y": Yp, "axis0_name": np.array(axis0_name), "axis1_name": np.array(axis1_name), **cropped},
                             fig_name, args.xmin, args.xmax, args.ymin, args.ymax)
        print(f"[ok] it3d={it3d} it2d={it2d} rl3d={'all' if rl3d is None else rl3d} dx={dx_used:.3f} km -> {fig_name}")

    print(f"Done. Outputs written to {plane_dir}")


if __name__ == "__main__":
    main()
