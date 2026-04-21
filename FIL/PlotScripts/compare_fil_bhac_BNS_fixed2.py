#!/usr/bin/env python3
"""
Compare FIL HDF5 2D slices with BHAC VTU data along a 1D line.

This version fixes:
- FIL aliases: rho -> rho_b, temperature -> temp, Y_e -> ye
- correct coordinate handling for xy/xz/yz in FIL
- auto-detection of FIL refinement levels
- BHAC plane handling:
    * understands the user's BHAC naming convention:
        d1 -> yz slice
        d2 -> xz slice
        d3 -> xy slice
      even though filenames contain x+0.00...
    * falls back to full 3D testNNNN.vtu if no matching d# slice exists
- robust BHAC variable lookup, including Ye = dye/d when needed
"""

import argparse
import glob
import os
import re
import xml.etree.ElementTree as ET

import h5py
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import rcParams

# ---------------------------------------------------------------------
# Matplotlib style
# ---------------------------------------------------------------------
rcParams["text.usetex"] = True
rcParams["font.family"] = "serif"
rcParams["font.serif"] = ["Computer Modern"]
rcParams["axes.labelsize"] = 18
rcParams["axes.linewidth"] = 1.5
rcParams["lines.markersize"] = 6
rcParams["xtick.labelsize"] = 18
rcParams["ytick.labelsize"] = 18
rcParams["font.size"] = 18


# ---------------------------------------------------------------------
# Canonical variable aliases
# ---------------------------------------------------------------------
FIL_ALIASES = {
    "rho": "rho_b",
    "rho_b": "rho_b",
    "density": "rho_b",
    "temperature": "temp",
    "temp": "temp",
    "t_eps": "temp",
    "t": "temp",
    "y_e": "ye",
    "ye": "ye",
    "enue": "Enue",
    "enue_bar": "Enue_bar",
    "enux": "Enux",
    "nnue": "Nnue",
    "nnue_bar": "Nnue_bar",
    "nnux": "Nnux",
    "qnue": "Qnue",
    "qnue_bar": "Qnue_bar",
    "qnux": "Qnux",
    "rnue": "Rnue",
    "rnue_bar": "Rnue_bar",
    "rnux": "Rnux",
    "kappa_nue_a": "kappa_nue_a",
    "kappa_nue_bar_a": "kappa_nue_bar_a",
    "kappa_nux_a": "kappa_nux_a",
    "kappa_nue_s": "kappa_nue_s",
    "kappa_nue_bar_s": "kappa_nue_bar_s",
    "kappa_nux_s": "kappa_nux_s",
    "kappa_nue_n": "kappa_nue_n",
    "kappa_nue_bar_n": "kappa_nue_bar_n",
    "kappa_nux_n": "kappa_nux_n",
}

FIL_TO_BHAC = {
    "rho_b": "rho",
    "ye": "ye",
    "temp": "temperature",
    "Nnue": "nrad1_prim",
    "Nnue_bar": "nrad2_prim",
    "Nnux": "nrad3_prim",
    "Enue": "erad1_prim",
    "Enue_bar": "erad2_prim",
    "Enux": "erad3_prim",
    "Fnue_x": "frad11",
    "Fnue_bar_x": "frad21",
    "Fnux_x": "frad31",
    "Fnue_y": "frad12",
    "Fnue_bar_y": "frad22",
    "Fnux_y": "frad32",
    "Fnue_z": "frad13",
    "Fnue_bar_z": "frad23",
    "Fnux_z": "frad33",
    "Qnue": "Q_er1",
    "Qnue_bar": "Q_er2",
    "Qnux": "Q_er3",
    "Rnue": "Q_nr1",
    "Rnue_bar": "Q_nr2",
    "Rnux": "Q_nr3",
    "kappa_nue_a": "kappa_a1",
    "kappa_nue_bar_a": "kappa_a2",
    "kappa_nux_a": "kappa_a3",
    "kappa_nue_s": "kappa_s1",
    "kappa_nue_bar_s": "kappa_s2",
    "kappa_nux_s": "kappa_s3",
    "kappa_nue_n": "kappa_n1",
    "kappa_nue_bar_n": "kappa_n2",
    "kappa_nux_n": "kappa_n3",
}

# ---------------------------------------------------------------------
# LaTeX y-axis labels for physical variables
# ---------------------------------------------------------------------
YLABEL_MAP = {
    "rho_b":           r"$\rho_{\rm b}$",
    "temp":            r"$T~[\rm MeV]$",
    "ye":              r"$Y_e$",
    "Enue":            r"$E_{\nu_e}$",
    "Enue_bar":        r"$E_{\bar{\nu}_e}$",
    "Enux":            r"$E_{\nu_x}$",
    "Nnue":            r"$N_{\nu_e}$",
    "Nnue_bar":        r"$N_{\bar{\nu}_e}$",
    "Nnux":            r"$N_{\nu_x}$",
    "Qnue":            r"$Q_{\nu_e}$",
    "Qnue_bar":        r"$Q_{\bar{\nu}_e}$",
    "Qnux":            r"$Q_{\nu_x}$",
    "Rnue":            r"$R_{\nu_e}$",
    "Rnue_bar":        r"$R_{\bar{\nu}_e}$",
    "Rnux":            r"$R_{\nu_x}$",
    "kappa_nue_a":     r"$\kappa^{a}_{\nu_e}$",
    "kappa_nue_bar_a": r"$\kappa^{a}_{\bar{\nu}_e}$",
    "kappa_nux_a":     r"$\kappa^{a}_{\nu_x}$",
    "kappa_nue_s":     r"$\kappa^{s}_{\nu_e}$",
    "kappa_nue_bar_s": r"$\kappa^{s}_{\bar{\nu}_e}$",
    "kappa_nux_s":     r"$\kappa^{s}_{\nu_x}$",
    "kappa_nue_n":     r"$\kappa^{n}_{\nu_e}$",
    "kappa_nue_bar_n": r"$\kappa^{n}_{\bar{\nu}_e}$",
    "kappa_nux_n":     r"$\kappa^{n}_{\nu_x}$",
    "Fnue_x":          r"$F^{x}_{\nu_e}$",
    "Fnue_bar_x":      r"$F^{x}_{\bar{\nu}_e}$",
    "Fnux_x":          r"$F^{x}_{\nu_x}$",
    "Fnue_y":          r"$F^{y}_{\nu_e}$",
    "Fnue_bar_y":      r"$F^{y}_{\bar{\nu}_e}$",
    "Fnux_y":          r"$F^{y}_{\nu_x}$",
    "Fnue_z":          r"$F^{z}_{\nu_e}$",
    "Fnue_bar_z":      r"$F^{z}_{\bar{\nu}_e}$",
    "Fnux_z":          r"$F^{z}_{\nu_x}$",
}

BHAC_ALIASES = {
    "rho": ["rho", "rhoprim", "rhoc", "d"],
    "ye": ["ye", "Y_e", "Ye", "dye_over_d"],
    "temperature": ["temperature", "temp", "T", "T_eps", "eps"],
    "nrad1_prim": ["nrad1_prim", "nrad1"],
    "nrad2_prim": ["nrad2_prim", "nrad2"],
    "nrad3_prim": ["nrad3_prim", "nrad3"],
    "erad1_prim": ["erad1_prim", "erad1"],
    "erad2_prim": ["erad2_prim", "erad2"],
    "erad3_prim": ["erad3_prim", "erad3"],
    "frad11": ["frad11"],
    "frad21": ["frad21"],
    "frad31": ["frad31"],
    "frad12": ["frad12"],
    "frad22": ["frad22"],
    "frad32": ["frad32"],
    "frad13": ["frad13"],
    "frad23": ["frad23"],
    "frad33": ["frad33"],
    "Q_er1": ["Q_er1"],
    "Q_er2": ["Q_er2"],
    "Q_er3": ["Q_er3"],
    "Q_nr1": ["Q_nr1"],
    "Q_nr2": ["Q_nr2"],
    "Q_nr3": ["Q_nr3"],
    "kappa_a1": ["kappa_a1"],
    "kappa_a2": ["kappa_a2"],
    "kappa_a3": ["kappa_a3"],
    "kappa_s1": ["kappa_s1"],
    "kappa_s2": ["kappa_s2"],
    "kappa_s3": ["kappa_s3"],
    "kappa_n1": ["kappa_n1"],
    "kappa_n2": ["kappa_n2"],
    "kappa_n3": ["kappa_n3"],
}


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
def normalize_name(name: str) -> str:
    return name.strip().replace("[", "_").replace("]", "").replace("-", "_").lower()


def canonical_fil_name(user_name: str) -> str:
    return FIL_ALIASES.get(normalize_name(user_name), user_name)


def plane_axes(plane: str):
    plane = plane.lower()
    if plane == "xy":
        return 0, 1, 2, "x", "y", "z"
    if plane == "xz":
        return 0, 2, 1, "x", "z", "y"
    if plane == "yz":
        return 1, 2, 0, "y", "z", "x"
    raise ValueError(f"Unsupported plane '{plane}'")


# User's BHAC convention
DSET_TO_PLANE = {
    1: "yz",
    2: "xz",
    3: "xy",
}
PLANE_TO_DSET = {v: k for k, v in DSET_TO_PLANE.items()}


def parse_bhac_slice_filename(filename: str):
    """
    Parse filenames like:
      test_d1_x+0.00D+00_n0000.vtu
      test_d2_x+0.00D+00_n0000.vtu
      test_d3_x+0.00D+00_n0000.vtu

    In this user's setup:
      d1 -> yz
      d2 -> xz
      d3 -> xy
    so the d# index determines the geometric plane.
    """
    base = os.path.basename(filename)
    m = re.search(r"_d(\d+)_([xyz])([+-]\d+\.\d+D[+-]\d+)_n(\d+)\.vtu$", base)
    if not m:
        return None
    dset = int(m.group(1))
    plane = DSET_TO_PLANE.get(dset)
    return {
        "kind": "slice",
        "dataset_index": dset,
        "plane": plane,
        "axis_tag": m.group(2),
        "fixed_value": float(m.group(3).replace("D", "E")),
        "iteration": int(m.group(4)),
        "filename": filename,
    }


def parse_bhac_volume_filename(filename: str):
    base = os.path.basename(filename)
    m = re.fullmatch(r"test(\d{4})\.vtu", base)
    if not m:
        return None
    return {
        "kind": "volume",
        "iteration": int(m.group(1)),
        "filename": filename,
    }


# ---------------------------------------------------------------------
# FIL loading
# ---------------------------------------------------------------------
def fil_file_map(data_dir: str, plane: str):
    mapping = {}
    pattern = os.path.join(data_dir, f"*.{plane}.h5")
    for path in glob.glob(pattern):
        base = os.path.basename(path)[: -len(f".{plane}.h5")]
        mapping[normalize_name(base)] = path
    return mapping


def fil_resolve_filename(data_dir: str, plane: str, user_var: str):
    canonical = canonical_fil_name(user_var)
    fmap = fil_file_map(data_dir, plane)
    key = normalize_name(canonical)
    if key in fmap:
        return canonical, fmap[key]
    key2 = normalize_name(user_var)
    if key2 in fmap:
        return user_var, fmap[key2]
    return canonical, None


def fil_available_rls(filename: str, iteration: int):
    rls = set()
    with h5py.File(filename, "r") as f:
        for k in f.keys():
            if f"it={iteration}" not in k:
                continue
            m = re.search(r"rl=(\d+)", k)
            if m:
                rls.add(int(m.group(1)))
    return sorted(rls)


def fil_extract_time(filename: str, iteration: int, rl: int):
    with h5py.File(filename, "r") as f:
        for k in f.keys():
            if f"it={iteration}" in k and f"rl={rl}" in k:
                dset = f[k]
                if "time" in dset.attrs:
                    return float(dset.attrs["time"])
                break
    return None


def fil_load_patches(filename: str, iteration: int, rl: int, plane: str):
    ax0, ax1, _, _, _, _ = plane_axes(plane)
    patches = []
    with h5py.File(filename, "r") as f:
        keys = list(f.keys())
        matched = [k for k in keys if f"it={iteration}" in k and f"rl={rl}" in k]
        if not matched:
            return patches

        full_var_name = matched[0].split(" it=")[0]
        prefix = f"{full_var_name} it={iteration} tl=0 rl={rl}"
        components = [k for k in keys if k.startswith(prefix)]

        for comp in components:
            dset = f[comp]
            data = np.asarray(dset[:], dtype=float)
            attrs = dict(dset.attrs)

            origin = np.asarray(attrs.get("origin", [0, 0, 0]), dtype=float)
            delta = np.asarray(attrs.get("delta", [1, 1, 1]), dtype=float)

            n2, n1 = data.shape
            dc1 = float(delta[ax0])
            dc2 = float(delta[ax1])
            c1 = origin[ax0] + np.arange(n1) * dc1
            c2 = origin[ax1] + np.arange(n2) * dc2

            patches.append(
                {
                    "data": data,
                    "c1": c1,
                    "c2": c2,
                    "dc1": dc1,
                    "dc2": dc2,
                    "c1_min": c1[0] - 0.5 * dc1,
                    "c1_max": c1[-1] + 0.5 * dc1,
                    "c2_min": c2[0] - 0.5 * dc2,
                    "c2_max": c2[-1] + 0.5 * dc2,
                }
            )
    return patches


def fil_find_patch(a, b, patches):
    candidates = []
    for p in patches:
        if p["c1_min"] <= a <= p["c1_max"] and p["c2_min"] <= b <= p["c2_max"]:
            candidates.append(p)
    if not candidates:
        return None
    return min(candidates, key=lambda p: p["dc1"] * p["dc2"])


def fil_bilinear_interpolate(a, b, patch):
    c1 = patch["c1"]
    c2 = patch["c2"]
    data = patch["data"]

    i = np.searchsorted(c1, a) - 1
    j = np.searchsorted(c2, b) - 1
    i = max(0, min(i, len(c1) - 2))
    j = max(0, min(j, len(c2) - 2))

    a0, a1 = c1[i], c1[i + 1]
    b0, b1 = c2[j], c2[j + 1]
    q00 = data[j, i]
    q10 = data[j, i + 1]
    q01 = data[j + 1, i]
    q11 = data[j + 1, i + 1]

    t = 0.0 if a1 == a0 else (a - a0) / (a1 - a0)
    s = 0.0 if b1 == b0 else (b - b0) / (b1 - b0)
    t = np.clip(t, 0.0, 1.0)
    s = np.clip(s, 0.0, 1.0)

    return (1 - t) * (1 - s) * q00 + t * (1 - s) * q10 + (1 - t) * s * q01 + t * s * q11


def fil_interpolate_line(patches, a_line, b_line):
    out = np.full(len(a_line), np.nan)
    for i, (a, b) in enumerate(zip(a_line, b_line)):
        patch = fil_find_patch(a, b, patches)
        if patch is not None:
            out[i] = fil_bilinear_interpolate(a, b, patch)
    return out


# ---------------------------------------------------------------------
# BHAC VTU parsing
# ---------------------------------------------------------------------
def _parse_piece_points(piece):
    n_points = int(piece.get("NumberOfPoints"))
    points_text = piece.find("Points").find("DataArray").text.strip().split()
    return np.asarray(points_text, dtype=float).reshape(n_points, -1)


def _parse_piece_cells(piece):
    cells_elem = piece.find("Cells")
    connectivity = np.asarray(
        cells_elem.find("DataArray[@Name='connectivity']").text.strip().split(),
        dtype=int
    )
    offsets = np.asarray(
        cells_elem.find("DataArray[@Name='offsets']").text.strip().split(),
        dtype=int
    )
    return connectivity, offsets


def _parse_piece_cell_data(piece):
    out = {}
    cell_data = piece.find("CellData")
    if cell_data is None:
        return out
    for da in cell_data.findall("DataArray"):
        name = da.get("Name")
        text = da.text.strip() if da.text else ""
        out[name] = np.asarray(text.split(), dtype=float) if text else np.array([], dtype=float)
    return out


def bhac_parse_vtu_full_geometry(filename: str):
    tree = ET.parse(filename)
    root = tree.getroot()

    time = None
    field_data = root.find(".//FieldData")
    if field_data is not None:
        for arr in field_data.findall("DataArray"):
            if arr.get("Name") == "TIME" and arr.text:
                time = float(arr.text.strip())
                break

    all_centers = []
    all_sizes = []
    all_data_lists = {}

    for piece in root.findall(".//Piece"):
        points = _parse_piece_points(piece)
        n_cells = int(piece.get("NumberOfCells"))
        connectivity, offsets = _parse_piece_cells(piece)

        centers = np.zeros((n_cells, 3), dtype=float)
        sizes = np.zeros((n_cells, 3), dtype=float)

        start = 0
        for icell, end in enumerate(offsets):
            vert_ids = connectivity[start:end]
            verts = points[vert_ids]
            centers[icell, :] = np.mean(verts[:, :3], axis=0)
            sizes[icell, 0] = verts[:, 0].max() - verts[:, 0].min()
            sizes[icell, 1] = verts[:, 1].max() - verts[:, 1].min()
            sizes[icell, 2] = verts[:, 2].max() - verts[:, 2].min()
            start = end

        all_centers.append(centers)
        all_sizes.append(sizes)

        piece_data = _parse_piece_cell_data(piece)
        for name, vals in piece_data.items():
            all_data_lists.setdefault(name, []).append(vals)

    centers_xyz = np.vstack(all_centers)
    sizes_xyz = np.vstack(all_sizes)
    data = {k: np.concatenate(v) for k, v in all_data_lists.items()}
    return centers_xyz, sizes_xyz, data, time


def bhac_available_slice_files(data_dir: str):
    files = sorted(glob.glob(os.path.join(data_dir, "*.vtu")))
    out = []
    for f in files:
        p = parse_bhac_slice_filename(f)
        if p is not None:
            out.append(p)
    return out


def bhac_select_best_source(data_dir: str, plane: str, iteration: int, fixed_value: float):
    """
    Prefer the user's dedicated d# slice files using:
      d1 -> yz
      d2 -> xz
      d3 -> xy
    If absent, fall back to full 3D testNNNN.vtu.
    """
    expected_dset = PLANE_TO_DSET[plane]
    slice_candidates = [
        p for p in bhac_available_slice_files(data_dir)
        if p["iteration"] == iteration and p["dataset_index"] == expected_dset
    ]
    if slice_candidates:
        slice_candidates.sort(key=lambda p: abs(p["fixed_value"] - fixed_value))
        best = slice_candidates[0]
        return {
            "mode": "slice_vtu",
            "filename": best["filename"],
            "plane": best["plane"],
            "slice_value": best["fixed_value"],
            "dataset_index": best["dataset_index"],
        }

    vol_name = os.path.join(data_dir, f"test{iteration:04d}.vtu")
    if os.path.exists(vol_name):
        return {
            "mode": "volume_vtu",
            "filename": vol_name,
            "plane": plane,
            "slice_value": fixed_value,
        }

    all_vtu = sorted(glob.glob(os.path.join(data_dir, "*.vtu")))
    parsed_vols = [parse_bhac_volume_filename(f) for f in all_vtu]
    parsed_vols = [p for p in parsed_vols if p is not None and p["iteration"] == iteration]
    if parsed_vols:
        return {
            "mode": "volume_vtu",
            "filename": parsed_vols[0]["filename"],
            "plane": plane,
            "slice_value": fixed_value,
        }

    return None


def bhac_extract_plane_geometry_from_source(source_info: dict, plane: str):
    centers_xyz, sizes_xyz, data, time = bhac_parse_vtu_full_geometry(source_info["filename"])
    ax0, ax1, axfix, _, _, _ = plane_axes(plane)

    if source_info["mode"] == "volume_vtu":
        slice_value = source_info["slice_value"]
        tol = np.maximum(0.51 * sizes_xyz[:, axfix], 1.0e-12)
        mask = np.abs(centers_xyz[:, axfix] - slice_value) <= tol
        if not np.any(mask):
            tol = np.maximum(1.01 * sizes_xyz[:, axfix], 1.0e-10)
            mask = np.abs(centers_xyz[:, axfix] - slice_value) <= tol

        centers2d = centers_xyz[mask][:, [ax0, ax1]]
        sizes2d = np.maximum(sizes_xyz[mask][:, ax0], sizes_xyz[mask][:, ax1])
        data2d = {k: v[mask] for k, v in data.items()}
    else:
        # Dedicated slice already matches requested plane by d# convention
        centers2d = centers_xyz[:, [ax0, ax1]]
        sizes2d = np.maximum(sizes_xyz[:, ax0], sizes_xyz[:, ax1])
        data2d = data

    return centers2d, sizes2d, data2d, time


def bhac_get_values(bhac_data: dict, requested_name: str):
    if requested_name in bhac_data:
        return requested_name, bhac_data[requested_name]

    aliases = BHAC_ALIASES.get(requested_name, [requested_name])
    for alias in aliases:
        if alias == "dye_over_d":
            if "dye" in bhac_data and "d" in bhac_data:
                d = bhac_data["d"]
                with np.errstate(divide="ignore", invalid="ignore"):
                    ye = np.where(np.abs(d) > 0.0, bhac_data["dye"] / d, np.nan)
                return "dye/d", ye
        elif alias in bhac_data:
            return alias, bhac_data[alias]
    return requested_name, None


def bhac_interpolate_line(centers2d, sizes2d, values, a_line, b_line):
    out = np.full(len(a_line), np.nan)
    c1 = centers2d[:, 0]
    c2 = centers2d[:, 1]

    for i, (a, b) in enumerate(zip(a_line, b_line)):
        dist2 = (c1 - a) ** 2 + (c2 - b) ** 2
        idx_near = np.argsort(dist2)[:64]
        if idx_near.size == 0:
            continue

        local_sizes = sizes2d[idx_near]
        finest = np.min(local_sizes)
        mask = (sizes2d <= 1.5 * finest) & (dist2 <= max((3.0 * finest) ** 2, 1.0e-14))
        use = np.where(mask)[0]
        if use.size < 4:
            use = idx_near[: min(16, idx_near.size)]

        uu1 = c1[use]
        uu2 = c2[use]
        uval = values[use]

        left = uu1 <= a
        right = uu1 >= a
        below = uu2 <= b
        above = uu2 >= b

        def pick(maskq):
            ids = np.where(maskq)[0]
            if ids.size == 0:
                return None
            j = ids[np.argmin((uu1[ids] - a) ** 2 + (uu2[ids] - b) ** 2)]
            return uu1[j], uu2[j], uval[j]

        q00 = pick(left & below)
        q10 = pick(right & below)
        q01 = pick(left & above)
        q11 = pick(right & above)

        if None in (q00, q10, q01, q11):
            j = np.argmin((uu1 - a) ** 2 + (uu2 - b) ** 2)
            out[i] = uval[j]
            continue

        a0, b0, v00 = q00
        a1, _, v10 = q10
        _, b1, v01 = q01
        _, _, v11 = q11

        t = 0.5 if a1 == a0 else (a - a0) / (a1 - a0)
        s = 0.5 if b1 == b0 else (b - b0) / (b1 - b0)
        t = np.clip(t, 0.0, 1.0)
        s = np.clip(s, 0.0, 1.0)
        out[i] = (1 - t) * (1 - s) * v00 + t * (1 - s) * v10 + (1 - t) * s * v01 + t * s * v11

    return out


# ---------------------------------------------------------------------
# Line definition
# ---------------------------------------------------------------------
def get_line_coordinates(line_type, n_points, line_range, offset=0.0, plane="xy"):
    amin, amax = line_range
    _, _, _, label_a, label_b, _ = plane_axes(plane)

    if line_type == "axis0":
        a_line = np.linspace(amin, amax, n_points)
        b_line = np.full(n_points, offset)
        coord_values = a_line
        coord_name = label_a
    elif line_type == "axis1":
        b_line = np.linspace(amin, amax, n_points)
        a_line = np.full(n_points, offset)
        coord_values = b_line
        coord_name = label_b
    elif line_type == "diagonal":
        a_line = np.linspace(amin, amax, n_points)
        b_line = np.linspace(amin, amax, n_points)
        coord_values = np.sqrt(a_line**2 + b_line**2)
        coord_name = "distance"
    else:
        raise ValueError(f"Unknown line type '{line_type}'")

    return a_line, b_line, coord_values, coord_name


# ---------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------
def save_dat(coord_values, fil_values, bhac_values, fil_name, coord_name, fil_iteration, bhac_iteration, output_dir):
    """Save comparison data to a .dat file with columns: coord, BHAC, FIL."""
    os.makedirs(output_dir, exist_ok=True)
    var_tag = normalize_name(fil_name)
    filename = os.path.join(output_dir, f"{var_tag}_itbhac{bhac_iteration:04d}_itfil{fil_iteration:04d}.dat")

    n = len(coord_values)
    bhac_col = bhac_values if bhac_values is not None else np.full(n, np.nan)
    fil_col  = fil_values  if fil_values  is not None else np.full(n, np.nan)

    header = f"{coord_name:<24s}  {'BHAC':<24s}  {'FIL':<24s}"
    data   = np.column_stack([coord_values, bhac_col, fil_col])
    np.savetxt(filename, data, header=header, fmt="%24.16e", comments="# ")
    print(f"  Saved dat: {filename}")


def plot_comparison(coord_values, fil_values, bhac_values, fil_name, bhac_name, fil_time, bhac_time, coord_name, output_file):
    fig, ax = plt.subplots(figsize=(7, 7))

    if fil_values is not None:
        m = np.isfinite(fil_values)
        ax.plot(coord_values[m], fil_values[m], lw=2, label="FIL")

    if bhac_values is not None:
        m = np.isfinite(bhac_values)
        ax.plot(coord_values[m], bhac_values[m], "--", lw=2, label="BHAC")

    # LaTeX y-axis label: look up by FIL canonical name, fall back to plain name
    ylabel = YLABEL_MAP.get(fil_name, YLABEL_MAP.get(bhac_name, fil_name))
    ax.set_xlabel(f"${coord_name}$")
    ax.set_ylabel(ylabel)
    ax.set_title(f"{fil_name} (FIL) vs {bhac_name} (BHAC)")
    ax.grid(True, alpha=0.3, linestyle="--")
    ax.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(output_file, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_file}")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Compare FIL HDF5 and BHAC VTU data along a 1D line")
    parser.add_argument("--fil-dir", required=True)
    parser.add_argument("--bhac-dir", required=True)
    parser.add_argument("--fil-iteration", type=int, required=True)
    parser.add_argument("--bhac-iteration", type=int, required=True)
    parser.add_argument("--variables", nargs="+", required=True)
    parser.add_argument("--output-dir", default="plots_comparison")
    parser.add_argument("--plane", choices=["xy", "xz", "yz"], default="xy")
    parser.add_argument("--line-type", choices=["axis0", "axis1", "diagonal"], default="axis0")
    parser.add_argument("--line-range", nargs=2, type=float, default=[-11.0, 11.0])
    parser.add_argument("--offset", type=float, default=0.0,
                        help="Fixed in-plane coordinate for axis0/axis1 lines")
    parser.add_argument("--plane-value", type=float, default=0.0,
                        help="Fixed out-of-plane coordinate value")
    parser.add_argument("--n-points", type=int, default=2000)
    parser.add_argument("--fil-rl", type=int, default=None)

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    data_dir = os.path.join(args.output_dir, "data")
    os.makedirs(data_dir, exist_ok=True)

    _, _, _, label_a, label_b, label_fix = plane_axes(args.plane)

    print("=" * 72)
    print("FIL vs BHAC comparison")
    print("=" * 72)
    print(f"FIL dir          : {args.fil_dir}")
    print(f"BHAC dir         : {args.bhac_dir}")
    print(f"FIL iteration    : {args.fil_iteration}")
    print(f"BHAC iteration   : {args.bhac_iteration}")
    print(f"Plane            : {args.plane} ({label_fix} = {args.plane_value})")
    print(f"Line type        : {args.line_type}")
    print(f"Line range       : {args.line_range}")
    print(f"Offset           : {args.offset}")
    print(f"Variables        : {args.variables}")
    print("=" * 72)

    a_line, b_line, coord_values, coord_name = get_line_coordinates(
        args.line_type, args.n_points, tuple(args.line_range), offset=args.offset, plane=args.plane
    )

    print("\nLoading BHAC data...")
    source_info = bhac_select_best_source(
        args.bhac_dir, args.plane, args.bhac_iteration, fixed_value=args.plane_value
    )
    if source_info is None:
        print("ERROR: no compatible BHAC data source found.")
        print(f"Tried dedicated d# slice files for plane {args.plane} and full 3D test{args.bhac_iteration:04d}.vtu.")
        return

    print(f"  BHAC source mode : {source_info['mode']}")
    print(f"  BHAC source file : {source_info['filename']}")
    if source_info["mode"] == "slice_vtu":
        print(f"  BHAC plane       : {source_info['plane']} via d{source_info['dataset_index']}")
        print(f"  Slice value      : {source_info['slice_value']}")
    else:
        print(f"  Extracted plane  : {args.plane} from full 3D VTU at {label_fix}={source_info['slice_value']}")

    bhac_centers2d, bhac_sizes2d, bhac_data, bhac_time = bhac_extract_plane_geometry_from_source(source_info, args.plane)
    print(f"  BHAC cells in plane : {len(bhac_centers2d)}")
    print(f"  BHAC time           : {bhac_time}")
    print(f"  BHAC vars available : {len(bhac_data)}")

    for user_var in args.variables:
        print(f"\n{user_var}:")
        fil_var, fil_filename = fil_resolve_filename(args.fil_dir, args.plane, user_var)
        bhac_requested = FIL_TO_BHAC.get(fil_var, user_var)

        print(f"  FIL canonical  : {fil_var}")
        print(f"  BHAC requested : {bhac_requested}")

        fil_values = None
        fil_time = None

        if fil_filename is None:
            print(f"  FIL: file for variable '{user_var}' / '{fil_var}' not found in plane {args.plane}")
        else:
            rls = fil_available_rls(fil_filename, args.fil_iteration)
            if not rls:
                print(f"  FIL: iteration {args.fil_iteration} not found in {os.path.basename(fil_filename)}")
            else:
                rl = args.fil_rl if args.fil_rl is not None else max(rls)
                if rl not in rls:
                    print(f"  FIL: requested rl={rl} not present. Available rls: {rls}")
                else:
                    print(f"  FIL file       : {os.path.basename(fil_filename)}")
                    print(f"  FIL RL         : {rl} (available: {rls})")
                    patches = fil_load_patches(fil_filename, args.fil_iteration, rl, args.plane)
                    fil_time = fil_extract_time(fil_filename, args.fil_iteration, rl)
                    if len(patches) == 0:
                        print("  FIL: no patches loaded")
                    else:
                        print(f"  FIL patches    : {len(patches)}")
                        fil_values = fil_interpolate_line(patches, a_line, b_line)

        bhac_resolved, bhac_raw = bhac_get_values(bhac_data, bhac_requested)
        bhac_values = None
        if bhac_raw is None:
            print(f"  BHAC: variable '{bhac_requested}' not found")
            sample = sorted(list(bhac_data.keys()))[:25]
            print(f"  BHAC available sample: {sample}")
        else:
            print(f"  BHAC resolved  : {bhac_resolved}")
            bhac_values = bhac_interpolate_line(bhac_centers2d, bhac_sizes2d, bhac_raw, a_line, b_line)

        if fil_values is None and bhac_values is None:
            print("  Skipping plot: neither side available.")
            continue

        save_dat(coord_values, fil_values, bhac_values, fil_var, coord_name,
                 args.fil_iteration, args.bhac_iteration, data_dir)

        out = os.path.join(
            args.output_dir,
            f"compare_{normalize_name(user_var)}_{args.plane}_fil{args.fil_iteration}_bhac{args.bhac_iteration}.png",
        )
        plot_comparison(coord_values, fil_values, bhac_values, fil_var, bhac_resolved, fil_time, bhac_time, coord_name, out)

    print("\nDone.")
    print(f"Plots written to: {args.output_dir}")


if __name__ == "__main__":
    main()


