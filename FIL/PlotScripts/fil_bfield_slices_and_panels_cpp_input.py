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
BFIELD_FILE_RE = re.compile(r"bfield_slice_(xy|xz|yz)_it(\d{8})\.h5$")

ALIASES_2D: Dict[str, List[str]] = {
    "W": ["w_lorentz", "W"],
    "Ye": ["Y_e", "Ye", "ye"],
    "rho": ["rho", "rho_b"],
    "smallb2": ["smallb2", "b2", "bsq"],
    "press": ["press", "P"],
}
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
class BFieldSlice:
    file: str
    iteration: int
    time_ms: float
    plane: str
    axis0_name: str
    axis1_name: str
    X: np.ndarray
    Y: np.ndarray
    bpol: np.ndarray
    btor: np.ndarray
    W: Optional[np.ndarray]
    bx_comov: Optional[np.ndarray]
    by_comov: Optional[np.ndarray]
    bz_comov: Optional[np.ndarray]
    b2_comov: Optional[np.ndarray]


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
                tiles.append(Tile2D(
                    file=file_path.as_posix(), dset=key, group=group, var=var,
                    iteration=int(it), tl=int(tl), rl=int(rl), c=int(c) if c is not None else -1,
                    time=time, origin=origin, delta=delta, shape_xy=shape_xy, corner=corner,
                ))
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


def read_dataset_values(tile) -> np.ndarray:
    with h5py.File(tile.file, "r") as f:
        return np.asarray(f[tile.dset], dtype=np.float32)


def collect_tiles_2d(store: VariableStore2D, name: str, iteration: int, requested_rl: Optional[int]) -> List[Tile2D]:
    return [
        store.index[name][iteration][rl][c]
        for rl in sorted(store.index[name][iteration].keys())
        if requested_rl is None or rl == requested_rl
        for c in sorted(store.index[name][iteration][rl].keys())
    ]


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


class BFieldSliceStore:
    def __init__(self, root: Path, plane: str):
        self.root = root.resolve()
        self.plane = plane
        self.index: Dict[int, BFieldSlice] = {}

    def candidate_dirs(self) -> List[Path]:
        dirs: List[Path] = []
        # explicit derived output dir
        for p in [self.root, self.root / self.plane, self.root / "derived_bfield_slices_cpp", self.root / "derived_bfield_slices_cpp" / self.plane]:
            if p.is_dir() and p not in dirs:
                dirs.append(p)
        # look recursively but keep it modest
        for pat in [f"**/bfield_slice_{self.plane}_it*.h5", f"**/{self.plane}/bfield_slice_{self.plane}_it*.h5"]:
            for fp in self.root.glob(pat):
                parent = fp.parent
                if parent.is_dir() and parent not in dirs:
                    dirs.append(parent)
        return sorted(dirs, key=natural_key)

    def build(self) -> None:
        seen: set[Path] = set()
        for d in self.candidate_dirs():
            for fp in sorted(d.glob(f"bfield_slice_{self.plane}_it*.h5"), key=natural_key):
                if fp in seen:
                    continue
                seen.add(fp)
                sl = self.read_one(fp)
                if sl is not None:
                    self.index[sl.iteration] = sl

    @staticmethod
    def _read_optional_dataset(f: h5py.File, name: str) -> Optional[np.ndarray]:
        if name in f:
            return np.asarray(f[name], dtype=np.float32)
        return None

    def read_one(self, file_path: Path) -> Optional[BFieldSlice]:
        m = BFIELD_FILE_RE.match(file_path.name)
        if not m:
            return None
        with h5py.File(file_path, "r") as f:
            if "axis0_km" not in f or "axis1_km" not in f or "bpol" not in f or "btor" not in f:
                return None
            x = np.asarray(f["axis0_km"], dtype=np.float32)
            y = np.asarray(f["axis1_km"], dtype=np.float32)
            X, Y = np.meshgrid(x, y)
            it = int(f.attrs.get("iteration", int(m.group(2))))
            time_ms = float(f.attrs.get("time_ms", np.nan))
            plane = str(f.attrs.get("plane", m.group(1)))
            axis0_name = f.attrs.get("axis0_name", "x")
            axis1_name = f.attrs.get("axis1_name", "z")
            if isinstance(axis0_name, bytes):
                axis0_name = axis0_name.decode()
            if isinstance(axis1_name, bytes):
                axis1_name = axis1_name.decode()
            return BFieldSlice(
                file=file_path.as_posix(), iteration=it, time_ms=time_ms, plane=plane,
                axis0_name=str(axis0_name), axis1_name=str(axis1_name),
                X=X.astype(np.float32), Y=Y.astype(np.float32),
                bpol=np.asarray(f["bpol"], dtype=np.float32),
                btor=np.asarray(f["btor"], dtype=np.float32),
                W=self._read_optional_dataset(f, "W"),
                bx_comov=self._read_optional_dataset(f, "bx_comov"),
                by_comov=self._read_optional_dataset(f, "by_comov"),
                bz_comov=self._read_optional_dataset(f, "bz_comov"),
                b2_comov=self._read_optional_dataset(f, "b2_comov"),
            )

    def common_iterations(self) -> List[int]:
        return sorted(self.index.keys())



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

    beta_min = min(-2.0, np.floor(beta_min))
    beta_max = max(2.0, np.ceil(beta_max))
    btor_min = min(12.0, np.floor(btor_min))
    btor_max = max(16.0, np.ceil(btor_max))
    bpol_min = min(12.0, np.floor(bpol_min))
    bpol_max = max(16.0, np.ceil(bpol_max))

    W = np.asarray(fields["W"], dtype=float)
    Wf = W[np.isfinite(W)]
    wmax = max(1.1, float(np.nanpercentile(Wf, 99.5)) if Wf.size else 1.1)

    Ye = np.asarray(fields["Ye"], dtype=float)
    Yef = Ye[np.isfinite(Ye)]
    yemax = min(0.6, max(0.06, float(np.nanpercentile(Yef, 99.5)) if Yef.size else 0.06))
    yemin = max(0.0, min(0.01, float(np.nanpercentile(Yef, 0.5)) if Yef.size else 0.0))

    return {
        "beta_inv": (float(beta_min), float(beta_max)),
        "btor": (float(btor_min), float(btor_max)),
        "bpol": (float(bpol_min), float(bpol_max)),
        "W": (1.0, float(wmax)),
        "Ye": (float(yemin), float(yemax)),
    }


def make_five_panel_plot(slice_all: Dict[str, np.ndarray], outpath: Path, xmin: float, xmax: float, ymin: float, ymax: float):
    X = slice_all["X"]
    Y = slice_all["Y"]
    fields = {
        "log10_beta_inv": safe_log10(slice_all["beta_inv"]),
        "log10_btor": safe_log10(slice_all["btor"] * CU_to_Gauss),
        "log10_bpol": safe_log10(slice_all["bpol"] * CU_to_Gauss),
        "W": slice_all["W"],
        "Ye": slice_all["Ye"],
    }
    ranges = choose_ranges(fields)

    fig = plt.figure(figsize=(30, 12))
    ax1 = fig.add_subplot(151)
    ax2 = fig.add_subplot(152, sharey=ax1)
    ax3 = fig.add_subplot(153, sharey=ax2)
    ax4 = fig.add_subplot(154, sharey=ax3)
    ax5 = fig.add_subplot(155, sharey=ax4)

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

    y_label = str(slice_all["axis1_name"])
    x_label = str(slice_all["axis0_name"])
    ax1.set_ylabel(rf"${y_label}~[{{\rm km}}]$", fontsize=24)
    for ax in [ax1, ax2, ax3, ax4, ax5]:
        ax.set_xlabel(rf"${x_label}~[{{\rm km}}]$", fontsize=24)
        ax.tick_params(which="both", top=True, right=True, direction="in", length=6, color="white", labelsize=18)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
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
    np.savez_compressed(
        outdir / f"{quantity}_it{iteration:08d}.npz",
        iteration=int(iteration), time_ms=float(time_ms), time_plot_ms=float(time_plot_ms),
        X=X.astype(np.float32), Y=Y.astype(np.float32), values=np.asarray(values, dtype=np.float32),
        quantity=np.array(quantity), axis0_name=np.array(axis0_name), axis1_name=np.array(axis1_name),
    )


def _time_map_2d_ms(store2d: VariableStore2D) -> Dict[int, float]:
    out: Dict[int, float] = {}
    if "rho" not in store2d.index:
        return out
    for it, rl_map in store2d.index["rho"].items():
        found = None
        for rl in sorted(rl_map.keys()):
            for c in sorted(rl_map[rl].keys()):
                found = rl_map[rl][c]
                break
            if found is not None:
                break
        if found is not None and np.isfinite(found.time):
            out[int(it)] = float(found.time) * CU_to_ms
    return out


def determine_iteration_pairs(store_bf: BFieldSliceStore, store2d: VariableStore2D, args):
    common_bf = store_bf.common_iterations()
    if not common_bf:
        raise SystemExit("No precomputed C++ bfield slices found.")

    missing_2d = store2d.missing(REQUIRED_2D)
    if missing_2d:
        raise SystemExit(f"Missing required 2D variables: {missing_2d}")
    common_2d = store2d.common_iterations(REQUIRED_2D)
    if not common_2d:
        raise SystemExit("No common 2D iterations found across required variables.")

    tmerg_ms = args.tmerg_code * CU_to_ms if args.tmerg_code is not None else None
    timebf_ms = {it: store_bf.index[it].time_ms for it in common_bf}
    time2d_ms = _time_map_2d_ms(store2d)

    def pass_filters(it: int, time_map: Dict[int, float]) -> bool:
        if args.itmin is not None and it < args.itmin:
            return False
        if args.itmax is not None and it > args.itmax:
            return False
        t_ms = time_map.get(int(it))
        if args.tmin is not None or args.tmax is not None:
            if t_ms is None or not np.isfinite(t_ms):
                return False
            t_plot = t_ms if tmerg_ms is None else (t_ms - tmerg_ms)
            if args.tmin is not None and t_plot < args.tmin:
                return False
            if args.tmax is not None and t_plot > args.tmax:
                return False
        return True

    eligible_bf = [int(it) for it in common_bf if pass_filters(int(it), timebf_ms)]
    eligible_2d = [int(it) for it in common_2d if pass_filters(int(it), time2d_ms)]
    if not eligible_bf:
        raise SystemExit("No C++ bfield iterations left after filters.")
    if not eligible_2d:
        raise SystemExit("No 2D iterations left after filters.")

    exact = sorted(set(eligible_bf).intersection(eligible_2d))
    if exact:
        return [(it, it) for it in exact]

    pairs = []
    for itbf in eligible_bf:
        tbf = timebf_ms.get(itbf, np.nan)
        if np.isfinite(tbf) and time2d_ms:
            best = min(eligible_2d, key=lambda it2d: (abs(time2d_ms.get(it2d, np.inf) - tbf), abs(it2d - itbf), it2d))
        else:
            best = min(eligible_2d, key=lambda it2d: (abs(it2d - itbf), it2d))
        pairs.append((itbf, best))
    return pairs


def main():
    ap = argparse.ArgumentParser(description="Read precomputed C++ bpol/btor slices and make the five-panel figure together with 2D FIL fields.")
    ap.add_argument("--root", required=True, help="Simulation root containing output*/data_hdf5_2D.")
    ap.add_argument("--bfield-root", default=None, help="Directory containing C++ bfield_slice_*.h5 outputs. Default: search under --root.")
    ap.add_argument("--plane", default="xz", choices=["xy", "xz", "yz"])
    ap.add_argument("--rl2d", type=int, default=None)
    ap.add_argument("--tmerg-code", type=float, default=None)
    ap.add_argument("--xmin", type=float, default=-120.0)
    ap.add_argument("--xmax", type=float, default=120.0)
    ap.add_argument("--ymin", type=float, default=0.0)
    ap.add_argument("--ymax", type=float, default=300.0)
    ap.add_argument("--itmin", type=int, default=None)
    ap.add_argument("--itmax", type=int, default=None)
    ap.add_argument("--tmin", type=float, default=None)
    ap.add_argument("--tmax", type=float, default=None)
    ap.add_argument("--outdir", default="derived_bfield_panels_from_cpp")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    bfield_root = Path(args.bfield_root).resolve() if args.bfield_root is not None else root
    outdir = Path(args.outdir).resolve()
    ensure_dir(outdir)

    store_bf = BFieldSliceStore(bfield_root, args.plane)
    store_bf.build()

    store2d = VariableStore2D(root, args.plane)
    store2d.build(REQUIRED_2D)

    pairs = determine_iteration_pairs(store_bf, store2d, args)
    plane_dir = outdir / args.plane
    plots_dir = plane_dir / "plots"
    ensure_dir(plots_dir)
    tmerg_ms = args.tmerg_code * CU_to_ms if args.tmerg_code is not None else None

    for itbf, it2d in pairs:
        bf = store_bf.index[int(itbf)]
        X = bf.X.astype(np.float32)
        Y = bf.Y.astype(np.float32)

        W2d, t2d_code = join_iteration_2d_to_mesh(store2d, "W", int(it2d), X, Y, args.rl2d)
        Ye2d, _ = join_iteration_2d_to_mesh(store2d, "Ye", int(it2d), X, Y, args.rl2d)
        rho2d, _ = join_iteration_2d_to_mesh(store2d, "rho", int(it2d), X, Y, args.rl2d)
        smallb2_2d, _ = join_iteration_2d_to_mesh(store2d, "smallb2", int(it2d), X, Y, args.rl2d)
        press2d, _ = join_iteration_2d_to_mesh(store2d, "press", int(it2d), X, Y, args.rl2d)

        beta_inv = smallb2_2d / (2.0 * np.where(np.abs(press2d) > 0.0, press2d, np.nan))
        rho_log_cgs = safe_log10(rho2d * CU_to_densCGS).astype(np.float32)

        time_ms = bf.time_ms if np.isfinite(bf.time_ms) else float(t2d_code) * CU_to_ms
        time_plot_ms = time_ms if tmerg_ms is None else (time_ms - tmerg_ms)

        Wplot = W2d
        if bf.W is not None:
            # keep 2D W as the plotted quantity for consistency with the original panel,
            # but store the C++ W alongside the bundle for debugging/comparison.
            Wcpp = bf.W.astype(np.float32)
        else:
            Wcpp = None

        for qname, qval in {
            "b_pol": bf.bpol,
            "b_tor": bf.btor,
            "plasma_beta_inv": beta_inv,
            "W": Wplot,
            "Ye": Ye2d,
        }.items():
            write_quantity_npz(plane_dir / qname, qname, int(itbf), time_ms, time_plot_ms, X, Y, qval, bf.axis0_name, bf.axis1_name)

        bundle = {
            "iteration_bfield": np.array(int(itbf)),
            "iteration_2d": np.array(int(it2d)),
            "time_ms": np.array(float(time_ms)),
            "time_plot_ms": np.array(float(time_plot_ms)),
            "X": X.astype(np.float32),
            "Y": Y.astype(np.float32),
            "b_pol": bf.bpol.astype(np.float32),
            "b_tor": bf.btor.astype(np.float32),
            "plasma_beta_inv": beta_inv.astype(np.float32),
            "W": Wplot.astype(np.float32),
            "Ye": Ye2d.astype(np.float32),
            "rho_log_cgs": rho_log_cgs,
            "smallb2_2d": smallb2_2d.astype(np.float32),
            "press_2d": press2d.astype(np.float32),
            "plane": np.array(args.plane),
            "axis0_name": np.array(bf.axis0_name),
            "axis1_name": np.array(bf.axis1_name),
            "rl2d": np.array(-1 if args.rl2d is None else int(args.rl2d)),
        }
        if bf.b2_comov is not None:
            bundle["b2_comov_cpp"] = bf.b2_comov.astype(np.float32)
        if Wcpp is not None:
            bundle["W_cpp"] = Wcpp
        if bf.bx_comov is not None:
            bundle["bx_comov"] = bf.bx_comov.astype(np.float32)
        if bf.by_comov is not None:
            bundle["by_comov"] = bf.by_comov.astype(np.float32)
        if bf.bz_comov is not None:
            bundle["bz_comov"] = bf.bz_comov.astype(np.float32)

        np.savez_compressed(plane_dir / f"bundle_itbf{int(itbf):08d}_it2d{int(it2d):08d}.npz", **bundle)

        Xp, Yp, cropped = crop_mesh_and_fields(
            X, Y,
            {
                "W": Wplot,
                "Ye": Ye2d,
                "rho_log_cgs": rho_log_cgs,
                "bpol": bf.bpol,
                "btor": bf.btor,
                "beta_inv": beta_inv,
            },
            args.xmin, args.xmax, args.ymin, args.ymax,
        )
        fig_name = plots_dir / f"five_panel_{args.plane}_itbf{int(itbf):08d}_it2d{int(it2d):08d}.png"
        make_five_panel_plot(
            {
                "X": Xp,
                "Y": Yp,
                "axis0_name": np.array(bf.axis0_name),
                "axis1_name": np.array(bf.axis1_name),
                **cropped,
            },
            fig_name, args.xmin, args.xmax, args.ymin, args.ymax,
        )
        print(f"[ok] itbf={itbf} it2d={it2d} -> {fig_name}")

    print(f"Done. Outputs written to {plane_dir}")


if __name__ == "__main__":
    main()
