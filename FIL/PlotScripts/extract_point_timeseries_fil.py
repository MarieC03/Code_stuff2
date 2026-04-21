#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import fnmatch
import os
import re
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import h5py
import numpy as np

# FIL / Cactus code-unit time conversion: 1 M_sun = 4.925490947e-3 ms
CU_TO_MS = 4.925490947e-3

OUTPUT_RE = re.compile(r"output-?\d{4}")
DSET_RE = re.compile(
    r"^(?P<group>.+?)::(?P<var>.+?)\s+it=(?P<it>\d+)\s+tl=(?P<tl>\d+)\s+rl=(?P<rl>\d+)(?:\s+c=(?P<c>-?\d+))?$"
)
SLICE_FILE_RE = re.compile(r"\.(?P<kind>xy|xz|yz|x|y|z|d)\.(?:file_\d+\.)?(?:h5|hdf5)$")

THREE_D_KIND = "xyz"
TWO_D_KINDS = {"xy", "xz", "yz"}
ONE_D_KINDS = {"x", "y", "z", "d"}
ALL_KINDS = {THREE_D_KIND, *TWO_D_KINDS, *ONE_D_KINDS}


@dataclass(frozen=True)
class TileRecord:
    file_path: str
    dset_name: str
    iteration: int
    time_msun: float
    rl: int
    tl: int
    component: int
    origin: np.ndarray
    delta: np.ndarray
    shape_xyz: np.ndarray
    corner: np.ndarray
    ndim: int
    slice_kind: str

    def contains_point(self, point: np.ndarray, tol: float = 1.0e-12) -> bool:
        lower = np.minimum(self.origin, self.corner) - tol
        upper = np.maximum(self.origin, self.corner) + tol

        if self.slice_kind == THREE_D_KIND:
            return bool(np.all(point >= lower) and np.all(point <= upper))

        if self.slice_kind == "xy":
            if abs(point[2] - self.origin[2]) > tol:
                return False
            return bool(np.all(point[:2] >= lower[:2]) and np.all(point[:2] <= upper[:2]))

        if self.slice_kind == "xz":
            if abs(point[1] - self.origin[1]) > tol:
                return False
            return bool(
                lower[0] <= point[0] <= upper[0] and lower[2] <= point[2] <= upper[2]
            )

        if self.slice_kind == "yz":
            if abs(point[0] - self.origin[0]) > tol:
                return False
            return bool(
                lower[1] <= point[1] <= upper[1] and lower[2] <= point[2] <= upper[2]
            )

        if self.slice_kind == "x":
            if abs(point[1] - self.origin[1]) > tol or abs(point[2] - self.origin[2]) > tol:
                return False
            return bool(lower[0] <= point[0] <= upper[0])

        if self.slice_kind == "y":
            if abs(point[0] - self.origin[0]) > tol or abs(point[2] - self.origin[2]) > tol:
                return False
            return bool(lower[1] <= point[1] <= upper[1])

        if self.slice_kind == "z":
            if abs(point[0] - self.origin[0]) > tol or abs(point[1] - self.origin[1]) > tol:
                return False
            return bool(lower[2] <= point[2] <= upper[2])

        if self.slice_kind == "d":
            return self._diagonal_index(point, tol=tol) is not None

        return False

    def _diagonal_index(self, point: np.ndarray, tol: float = 1.0e-10) -> Optional[int]:
        if self.slice_kind != "d":
            return None

        rels: List[float] = []
        for dim in range(3):
            if abs(self.delta[dim]) <= tol:
                if abs(point[dim] - self.origin[dim]) > tol:
                    return None
                continue
            rels.append((point[dim] - self.origin[dim]) / self.delta[dim])

        if not rels:
            return 0

        idxf = rels[0]
        if not all(abs(r - idxf) <= tol for r in rels[1:]):
            return None

        idx = int(np.rint(idxf))
        if idx < 0 or idx >= int(self.shape_xyz[0]):
            return None

        reconstructed = self.origin + idx * self.delta
        if not np.allclose(reconstructed, point, atol=tol, rtol=0.0):
            return None
        return idx

    def grid_indices(self, point: np.ndarray, tol: float = 1.0e-10) -> Optional[Tuple[int, ...]]:
        if self.slice_kind == THREE_D_KIND:
            rel = (point - self.origin) / self.delta
            idx = np.rint(rel).astype(int)
            if np.any(idx < 0) or np.any(idx >= self.shape_xyz):
                return None
            reconstructed = self.origin + idx * self.delta
            if not np.allclose(reconstructed, point, atol=tol, rtol=0.0):
                return None
            return tuple(int(i) for i in idx)

        if self.slice_kind == "xy":
            ix = _grid_index_1d(point[0], self.origin[0], self.delta[0], int(self.shape_xyz[0]), tol)
            iy = _grid_index_1d(point[1], self.origin[1], self.delta[1], int(self.shape_xyz[1]), tol)
            if ix is None or iy is None or abs(point[2] - self.origin[2]) > tol:
                return None
            return (ix, iy, 0)

        if self.slice_kind == "xz":
            ix = _grid_index_1d(point[0], self.origin[0], self.delta[0], int(self.shape_xyz[0]), tol)
            iz = _grid_index_1d(point[2], self.origin[2], self.delta[2], int(self.shape_xyz[2]), tol)
            if ix is None or iz is None or abs(point[1] - self.origin[1]) > tol:
                return None
            return (ix, 0, iz)

        if self.slice_kind == "yz":
            iy = _grid_index_1d(point[1], self.origin[1], self.delta[1], int(self.shape_xyz[1]), tol)
            iz = _grid_index_1d(point[2], self.origin[2], self.delta[2], int(self.shape_xyz[2]), tol)
            if iy is None or iz is None or abs(point[0] - self.origin[0]) > tol:
                return None
            return (0, iy, iz)

        if self.slice_kind == "x":
            ix = _grid_index_1d(point[0], self.origin[0], self.delta[0], int(self.shape_xyz[0]), tol)
            if ix is None or abs(point[1] - self.origin[1]) > tol or abs(point[2] - self.origin[2]) > tol:
                return None
            return (ix, 0, 0)

        if self.slice_kind == "y":
            iy = _grid_index_1d(point[1], self.origin[1], self.delta[1], int(self.shape_xyz[1]), tol)
            if iy is None or abs(point[0] - self.origin[0]) > tol or abs(point[2] - self.origin[2]) > tol:
                return None
            return (0, iy, 0)

        if self.slice_kind == "z":
            iz = _grid_index_1d(point[2], self.origin[2], self.delta[2], int(self.shape_xyz[2]), tol)
            if iz is None or abs(point[0] - self.origin[0]) > tol or abs(point[1] - self.origin[1]) > tol:
                return None
            return (0, 0, iz)

        if self.slice_kind == "d":
            idx = self._diagonal_index(point, tol=tol)
            if idx is None:
                return None
            return (idx, idx, idx)

        return None


@dataclass
class IterationValue:
    iteration: int
    time_msun: float
    rl: int
    value: float


class PointNotOnGridError(RuntimeError):
    pass


class VariableNotFoundError(RuntimeError):
    pass


def _grid_index_1d(coord: float, origin: float, delta: float, n: int, tol: float) -> Optional[int]:
    if abs(delta) <= tol:
        return 0 if abs(coord - origin) <= tol else None
    rel = (coord - origin) / delta
    idx = int(np.rint(rel))
    if idx < 0 or idx >= n:
        return None
    if abs((origin + idx * delta) - coord) > tol:
        return None
    return idx


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Extract one-point FIL time series from Carpet HDF5 output files and "
            "write ASCII files with columns: iteration, time(Msun), time(ms), value."
        )
    )
    p.add_argument("--simdir", required=True, help="Simulation root directory.")
    p.add_argument("--x", type=float, required=True, help="Target x coordinate in code units (Msun).")
    p.add_argument("--y", type=float, required=True, help="Target y coordinate in code units (Msun).")
    p.add_argument("--z", type=float, required=True, help="Target z coordinate in code units (Msun).")
    p.add_argument("--vars", nargs="+", required=True, help="Variable names, e.g. rho temp ye alp eps")
    p.add_argument("--outdir", default="point_timeseries", help="Output directory.")
    p.add_argument(
        "--rl-mode",
        choices=["finest", "coarsest"],
        default="finest",
        help="How to choose between multiple refinement levels in the same iteration.",
    )
    p.add_argument(
        "--tol",
        type=float,
        default=1.0e-10,
        help="Absolute tolerance for checking whether the point lies on a grid point.",
    )
    p.add_argument(
        "--plane",
        choices=["xyz", "xy", "xz", "yz", "x", "y", "z", "d"],
        default="xyz",
        help=(
            "Restrict search to full 3D files ('xyz'), a 2D plane ('xy','xz','yz'), "
            "or a 1D line file ('x','y','z','d'). Here 'd' denotes the diagonal 1D output."
        ),
    )
    p.add_argument(
        "--recursive-glob",
        default="*.h5",
        help="Filename pattern searched recursively below output directories. Default: %(default)s",
    )
    p.add_argument(
        "--verbose-skip",
        action="store_true",
        help="Print skipped/corrupt files while scanning.",
    )
    return p.parse_args()


def path_is_under_output_tree(path: str) -> bool:
    parts = os.path.normpath(path).split(os.sep)
    return any(OUTPUT_RE.fullmatch(part) for part in parts)


def guess_slice_kind_from_path(path: str) -> str:
    m = SLICE_FILE_RE.search(path)
    if m:
        return m.group("kind")
    return THREE_D_KIND


def slice_file_is_allowed(path: str, requested_kind: str) -> bool:
    file_kind = guess_slice_kind_from_path(path)
    if requested_kind == THREE_D_KIND:
        return file_kind == THREE_D_KIND
    return file_kind == requested_kind


def discover_candidate_files(simdir: str, file_glob: str, requested_kind: str) -> List[str]:
    candidates: List[str] = []
    for root, _, files in os.walk(simdir):
        if not path_is_under_output_tree(root):
            continue
        for fname in files:
            if not fnmatch.fnmatch(fname, file_glob):
                continue
            path = os.path.join(root, fname)
            if not slice_file_is_allowed(path, requested_kind):
                continue
            candidates.append(path)
    candidates.sort()
    return candidates


def safe_open_h5(path: str) -> Optional[h5py.File]:
    try:
        if os.path.getsize(path) == 0:
            return None
    except OSError:
        return None
    try:
        return h5py.File(path, "r")
    except Exception:
        return None


def get_var_from_dataset_name(dset_name: str) -> Optional[str]:
    m = DSET_RE.match(dset_name)
    return None if m is None else m.group("var")


def iter_matching_dataset_names(h5f: h5py.File, wanted_vars: Optional[set[str]] = None) -> Iterable[str]:
    """Robustly iterate dataset names and survive partially corrupt HDF5 files."""
    try:
        root_keys = list(h5f.keys())
    except Exception:
        return

    for name in root_keys:
        try:
            obj = h5f.get(name, default=None)
        except Exception:
            continue
        if obj is None:
            continue
        if not isinstance(obj, h5py.Dataset):
            continue
        var = get_var_from_dataset_name(name)
        if var is None:
            continue
        if wanted_vars is not None and var not in wanted_vars:
            continue
        yield name


def build_var_file_map(
    files: Sequence[str], wanted_vars: Sequence[str], verbose_skip: bool = False
) -> Tuple[Dict[str, List[str]], Dict[str, int]]:
    wanted = set(wanted_vars)
    found: Dict[str, List[str]] = {v: [] for v in wanted_vars}
    stats = {
        "empty_or_unopenable": 0,
        "corrupt_root": 0,
        "no_matching_datasets": 0,
        "used": 0,
    }

    for path in files:
        h5f = safe_open_h5(path)
        if h5f is None:
            stats["empty_or_unopenable"] += 1
            if verbose_skip:
                print(f"[SKIP] unreadable or empty HDF5 file: {path}")
            continue
        try:
            hit_vars = set(iter_matching_dataset_names(h5f, wanted))
            if not hit_vars:
                stats["no_matching_datasets"] += 1
                continue
            resolved_vars = {get_var_from_dataset_name(name) for name in hit_vars}
            for var in sorted(v for v in resolved_vars if v is not None):
                found[var].append(path)
            stats["used"] += 1
        except Exception:
            stats["corrupt_root"] += 1
            if verbose_skip:
                print(f"[SKIP] corrupt HDF5 symbol table / root group: {path}")
        finally:
            try:
                h5f.close()
            except Exception:
                pass
    return found, stats


def parse_tile_record(file_path: str, dset_name: str, dset: h5py.Dataset) -> Optional[TileRecord]:
    m = DSET_RE.match(dset_name)
    if m is None:
        return None
    try:
        origin = np.asarray(dset.attrs["origin"], dtype=float)
        delta = np.asarray(dset.attrs["delta"], dtype=float)
        time_msun = float(np.asarray(dset.attrs["time"]))
    except Exception:
        return None

    ndim = dset.ndim
    slice_kind = guess_slice_kind_from_path(file_path)

    if origin.shape[0] != 3 or delta.shape[0] != 3:
        return None

    if ndim == 3:
        if slice_kind != THREE_D_KIND:
            slice_kind = THREE_D_KIND
        shape_xyz = np.asarray(dset.shape[::-1], dtype=int)
    elif ndim == 2:
        if slice_kind == "xy":
            shape_xyz = np.array([dset.shape[1], dset.shape[0], 1], dtype=int)
        elif slice_kind == "xz":
            shape_xyz = np.array([dset.shape[1], 1, dset.shape[0]], dtype=int)
        elif slice_kind == "yz":
            shape_xyz = np.array([1, dset.shape[1], dset.shape[0]], dtype=int)
        else:
            return None
    elif ndim == 1:
        if slice_kind == "x":
            shape_xyz = np.array([dset.shape[0], 1, 1], dtype=int)
        elif slice_kind == "y":
            shape_xyz = np.array([1, dset.shape[0], 1], dtype=int)
        elif slice_kind == "z":
            shape_xyz = np.array([1, 1, dset.shape[0]], dtype=int)
        elif slice_kind == "d":
            shape_xyz = np.array([dset.shape[0], dset.shape[0], dset.shape[0]], dtype=int)
        else:
            # Fallback: a plain 1D dataset without suffix is ambiguous, so ignore it.
            return None
    else:
        return None

    if ndim == 1 and slice_kind in ONE_D_KINDS:
        if slice_kind == "x":
            corner = origin + np.array([delta[0] * (shape_xyz[0] - 1), 0.0, 0.0])
        elif slice_kind == "y":
            corner = origin + np.array([0.0, delta[1] * (shape_xyz[1] - 1), 0.0])
        elif slice_kind == "z":
            corner = origin + np.array([0.0, 0.0, delta[2] * (shape_xyz[2] - 1)])
        else:  # diagonal
            n = int(dset.shape[0])
            corner = origin + delta * (n - 1)
    else:
        corner = origin + delta * (shape_xyz - 1)

    return TileRecord(
        file_path=file_path,
        dset_name=dset_name,
        iteration=int(m.group("it")),
        time_msun=time_msun,
        rl=int(m.group("rl")),
        tl=int(m.group("tl")),
        component=int(m.group("c")) if m.group("c") is not None else -1,
        origin=origin,
        delta=delta,
        shape_xyz=shape_xyz,
        corner=corner,
        ndim=ndim,
        slice_kind=slice_kind,
    )


def collect_tiles_for_var(
    var: str, files: Sequence[str], point: np.ndarray, verbose_skip: bool = False
) -> Dict[int, List[TileRecord]]:
    by_iteration: Dict[int, List[TileRecord]] = {}
    for path in files:
        h5f = safe_open_h5(path)
        if h5f is None:
            continue
        try:
            for dset_name in iter_matching_dataset_names(h5f, {var}):
                try:
                    dset = h5f.get(dset_name, default=None)
                except Exception:
                    continue
                if dset is None or not isinstance(dset, h5py.Dataset):
                    continue
                rec = parse_tile_record(path, dset_name, dset)
                if rec is None:
                    continue
                if rec.contains_point(point):
                    by_iteration.setdefault(rec.iteration, []).append(rec)
        except Exception:
            if verbose_skip:
                print(f"[SKIP] failed while scanning variable '{var}' in file: {path}")
        finally:
            try:
                h5f.close()
            except Exception:
                pass
    return by_iteration


def choose_record(records: Sequence[TileRecord], rl_mode: str) -> TileRecord:
    if rl_mode == "finest":
        return max(records, key=lambda r: (r.rl, r.tl, r.component))
    return min(records, key=lambda r: (r.rl, r.tl, r.component))


def read_value_from_record(rec: TileRecord, point: np.ndarray, tol: float) -> float:
    idx_xyz = rec.grid_indices(point, tol=tol)
    if idx_xyz is None:
        raise PointNotOnGridError(
            f"Point {tuple(point)} is inside tile bounds but not on a grid point for {rec.dset_name}"
        )
    ix, iy, iz = idx_xyz
    with h5py.File(rec.file_path, "r") as h5f:
        arr = h5f[rec.dset_name]
        if arr.ndim == 3:
            value = arr[iz, iy, ix]
        elif arr.ndim == 2:
            if rec.slice_kind == "xy":
                value = arr[iy, ix]
            elif rec.slice_kind == "xz":
                value = arr[iz, ix]
            elif rec.slice_kind == "yz":
                value = arr[iz, iy]
            else:
                raise PointNotOnGridError(f"2D dataset without identifiable plane: {rec.dset_name}")
        elif arr.ndim == 1:
            if rec.slice_kind == "x":
                value = arr[ix]
            elif rec.slice_kind == "y":
                value = arr[iy]
            elif rec.slice_kind == "z":
                value = arr[iz]
            elif rec.slice_kind == "d":
                value = arr[ix]
            else:
                raise PointNotOnGridError(f"1D dataset without identifiable line kind: {rec.dset_name}")
        else:
            raise RuntimeError(f"Unsupported dataset dimension {arr.ndim} in {rec.dset_name}")
    return float(np.asarray(value))


def extract_var_series(
    var: str,
    files: Sequence[str],
    point: np.ndarray,
    rl_mode: str,
    tol: float,
    verbose_skip: bool = False,
) -> List[IterationValue]:
    by_iteration = collect_tiles_for_var(var, files, point, verbose_skip=verbose_skip)
    if not by_iteration:
        raise VariableNotFoundError(f"No tiles for variable '{var}' contain the point {tuple(point)}")

    series: List[IterationValue] = []
    for it in sorted(by_iteration):
        rec = choose_record(by_iteration[it], rl_mode=rl_mode)
        try:
            value = read_value_from_record(rec, point, tol=tol)
        except PointNotOnGridError:
            continue
        series.append(IterationValue(iteration=rec.iteration, time_msun=rec.time_msun, rl=rec.rl, value=value))
    return series


def sanitize_varname(var: str) -> str:
    out = var.replace("/", "_")
    out = out.replace("[", "_").replace("]", "")
    out = out.replace(" ", "_")
    return out


def write_ascii(var: str, series: Sequence[IterationValue], outdir: str, point: np.ndarray, plane: str) -> str:
    os.makedirs(outdir, exist_ok=True)
    suffix = f"_{plane}" if plane != THREE_D_KIND else ""
    outname = f"{sanitize_varname(var)}_x{point[0]:g}_y{point[1]:g}_z{point[2]:g}{suffix}.asc"
    outpath = os.path.join(outdir, outname)
    with open(outpath, "w", encoding="utf-8") as f:
        f.write(f"# variable = {var}\n")
        f.write(f"# plane = {plane}\n")
        f.write("# columns: iteration time_Msun time_ms value\n")
        for row in series:
            t_ms = row.time_msun * CU_TO_MS
            f.write(f"{row.iteration:12d} {row.time_msun:24.16e} {t_ms:24.16e} {row.value:24.16e}\n")
    return outpath


if __name__ == "__main__":
    args = parse_args()
    point = np.array([args.x, args.y, args.z], dtype=float)

    if args.plane == "xy" and abs(args.z) > args.tol:
        raise SystemExit("For --plane xy you must set --z 0 (within --tol).")
    if args.plane == "xz" and abs(args.y) > args.tol:
        raise SystemExit("For --plane xz you must set --y 0 (within --tol).")
    if args.plane == "yz" and abs(args.x) > args.tol:
        raise SystemExit("For --plane yz you must set --x 0 (within --tol).")

    all_files = discover_candidate_files(args.simdir, args.recursive_glob, args.plane)
    if not all_files:
        raise SystemExit(
            f"No readable HDF5 files matching '{args.recursive_glob}' for selection '{args.plane}' found below {args.simdir}."
        )

    var_to_files, scan_stats = build_var_file_map(all_files, args.vars, verbose_skip=args.verbose_skip)

    print(f"Found {len(all_files)} candidate HDF5 files below: {args.simdir}")
    print(f"Target point in code units: x={args.x}, y={args.y}, z={args.z}")
    print(f"Slice selection: {args.plane}")
    print(f"Refinement-level selection mode: {args.rl_mode}")
    print(
        "Scan summary: "
        f"used={scan_stats['used']}, "
        f"no_matching_datasets={scan_stats['no_matching_datasets']}, "
        f"empty_or_unopenable={scan_stats['empty_or_unopenable']}, "
        f"corrupt_root={scan_stats['corrupt_root']}\n"
    )

    failures: List[str] = []
    for var in args.vars:
        files = var_to_files.get(var, [])
        if not files:
            failures.append(var)
            print(f"[WARNING] Variable '{var}' was not found in any readable/scannable HDF5 file. Skipping.")
            continue
        try:
            series = extract_var_series(
                var, files, point, args.rl_mode, args.tol, verbose_skip=args.verbose_skip
            )
        except VariableNotFoundError as exc:
            failures.append(var)
            print(f"[WARNING] {exc}")
            continue
        if not series:
            failures.append(var)
            print(f"[WARNING] Variable '{var}' had no valid samples at the requested point. Skipping.")
            continue
        outpath = write_ascii(var, series, args.outdir, point, args.plane)
        print(f"[OK] {var:>16s} -> {outpath}   ({len(series)} rows)")

    if failures:
        print("\nFinished with warnings for:", ", ".join(failures))
