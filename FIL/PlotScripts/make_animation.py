#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, Optional

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.image as mpimg


# 1 code time unit (1 Msun in geometric units) in ms
CU_TO_MS = 4.925490947e-3


def extract_iteration(path: Path) -> int:
    """
    Extract iteration number from filenames like:
    five_panel_xz_it00128000.png
    """
    m = re.search(r"_it(\d+)\.png$", path.name)
    if m is None:
        raise ValueError(f"Could not extract iteration from filename: {path.name}")
    return int(m.group(1))


def natural_sort_images(files):
    return sorted(files, key=extract_iteration)


def detect_run_kind(files, indir: Path) -> str:
    """
    Detect whether this is an M1 or noM1 case.
    We must check 'noM1' before 'M1' because 'M1' is contained in 'noM1'.
    """
    texts = [str(indir)] + [f.name for f in files]

    for text in texts:
        if "noM1" in text:
            return "noM1"

    for text in texts:
        if "M1" in text:
            return "M1"

    # fallback
    return "M1"


def try_read_table(path: Path) -> Optional[pd.DataFrame]:
    """
    Try to read a time table from CSV or whitespace-separated ASCII.
    Expected useful columns include:
      iteration
      time_code
      time_ms
    """
    suffix = path.suffix.lower()

    try:
        if suffix == ".csv":
            df = pd.read_csv(path)
        else:
            # fallback for ASCII-like files
            df = pd.read_csv(path, comment="#", delim_whitespace=True)
    except Exception:
        return None

    # normalize column names
    df.columns = [str(c).strip() for c in df.columns]

    if "iteration" not in df.columns:
        return None

    if ("time_code" not in df.columns) and ("time_ms" not in df.columns):
        return None

    return df


def find_time_file(indir: Path, run_kind: str) -> Path:
    """
    Find a CSV/ASCII file in indir that contains iteration and time information.
    Preference:
      1) matching run kind in filename
      2) files containing rho and/or max
      3) any file with the required columns
    """
    candidates = []

    for ext in ("*.csv", "*.asc", "*.txt", "*.dat"):
        candidates.extend(indir.glob(ext))

    if not candidates:
        raise FileNotFoundError(f"No CSV/ASCII files found in {indir.resolve()}")

    good = []
    for path in candidates:
        name = path.name

        # strict run-kind filtering
        if run_kind == "noM1":
            if "noM1" not in name:
                continue
        elif run_kind == "M1":
            # exclude noM1 here
            if ("M1" not in name) or ("noM1" in name):
                continue

        df = try_read_table(path)
        if df is not None:
            score = 0
            lname = name.lower()
            if "rho" in lname:
                score += 3
            if "max" in lname or "maximum" in lname:
                score += 2
            if "time" in lname:
                score += 1
            good.append((score, path))

    if good:
        good.sort(key=lambda x: (-x[0], x[1].name))
        return good[0][1]

    # fallback: try any readable table
    for path in candidates:
        df = try_read_table(path)
        if df is not None:
            return path

    raise FileNotFoundError(
        f"Could not find a readable time table in {indir.resolve()} for run kind '{run_kind}'"
    )


def build_time_lookup(time_file: Path, t_merger_code: float) -> Dict[int, float]:
    """
    Build mapping:
      iteration -> (time_code - t_merger_code) * CU_TO_MS  [ms]
    """
    df = try_read_table(time_file)
    if df is None:
        raise ValueError(f"Could not read time table from {time_file}")

    # clean duplicates if any
    df = df.dropna(subset=["iteration"]).copy()
    df["iteration"] = df["iteration"].astype(int)

    if "time_code" in df.columns:
        df = df.dropna(subset=["time_code"]).copy()
        df["time_rel_ms"] = (df["time_code"].astype(float) - t_merger_code) * CU_TO_MS
    elif "time_ms" in df.columns:
        # if only time_ms exists, convert merger time to ms and subtract
        t_merger_ms = t_merger_code * CU_TO_MS
        df = df.dropna(subset=["time_ms"]).copy()
        df["time_rel_ms"] = df["time_ms"].astype(float) - t_merger_ms
    else:
        raise ValueError(
            f"Time table {time_file} has 'iteration' but neither 'time_code' nor 'time_ms'"
        )

    # keep last occurrence if duplicate iterations exist
    df = df.drop_duplicates(subset=["iteration"], keep="last")
    return dict(zip(df["iteration"], df["time_rel_ms"]))


def format_title(t_rel_ms: float) -> str:
    return rf"$t - t_{{\rm mer}} = {t_rel_ms:.2f}\,\mathrm{{ms}}$"


def main():
    parser = argparse.ArgumentParser(
        description="Create an animation from plot images, labeling frames with t-t_mer in ms."
    )
    parser.add_argument(
        "--indir",
        type=str,
        default=".",
        help="Directory containing the images and time table"
    )
    parser.add_argument(
        "--pattern",
        type=str,
        default="five_panel_xz_it*.png",
        help="Glob pattern for input images"
    )
    parser.add_argument(
        "--out",
        type=str,
        default="five_panel_xz_animation.mp4",
        help="Output animation filename (.mp4 or .gif)"
    )
    parser.add_argument(
        "--fps",
        type=int,
        default=10,
        help="Frames per second"
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="Output DPI"
    )
    parser.add_argument(
        "--t-merger",
        type=float,
        required=True,
        help="Merger time in code units (Msun units)"
    )
    parser.add_argument(
        "--time-file",
        type=str,
        default=None,
        help="Optional explicit CSV/ASCII file with iteration/time mapping"
    )
    parser.add_argument(
        "--title-fontsize",
        type=int,
        default=18,
        help="Font size of the time label"
    )

    args = parser.parse_args()

    indir = Path(args.indir)
    files = list(indir.glob(args.pattern))

    if not files:
        raise FileNotFoundError(
            f"No files found in {indir.resolve()} matching pattern: {args.pattern}"
        )

    files = natural_sort_images(files)

    print(f"Found {len(files)} frames.")
    print("First file:", files[0].name)
    print("Last file: ", files[-1].name)

    run_kind = detect_run_kind(files, indir)
    print(f"Detected run kind: {run_kind}")

    if args.time_file is not None:
        time_file = Path(args.time_file)
        if not time_file.is_absolute():
            time_file = indir / time_file
        if not time_file.exists():
            raise FileNotFoundError(f"Specified --time-file does not exist: {time_file}")
    else:
        time_file = find_time_file(indir, run_kind)

    print(f"Using time table: {time_file}")

    time_lookup = build_time_lookup(time_file, args.t_merger)

    # Read first image to set up figure size
    first_img = mpimg.imread(files[0])
    ny, nx = first_img.shape[:2]

    fig_w = nx / 100.0
    fig_h = ny / 100.0

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off")

    im = ax.imshow(first_img)

    # initial title
    first_it = extract_iteration(files[0])
    if first_it not in time_lookup:
        raise KeyError(
            f"Iteration {first_it} from image {files[0].name} not found in time table {time_file}"
        )
    title = ax.set_title(
        format_title(time_lookup[first_it]),
        fontsize=args.title_fontsize,
        pad=12
    )

    def update(frame_idx):
        path = files[frame_idx]
        iteration = extract_iteration(path)

        if iteration not in time_lookup:
            raise KeyError(
                f"Iteration {iteration} from image {path.name} not found in time table {time_file}"
            )

        img = mpimg.imread(path)
        im.set_data(img)
        title.set_text(format_title(time_lookup[iteration]))
        return [im, title]

    ani = animation.FuncAnimation(
        fig,
        update,
        frames=len(files),
        interval=1000 / args.fps,
        blit=True
    )

    out_path = Path(args.out)
    suffix = out_path.suffix.lower()

    try:
        if suffix == ".gif":
            ani.save(out_path, writer="pillow", fps=args.fps, dpi=args.dpi)
        else:
            ani.save(out_path, writer="ffmpeg", fps=args.fps, dpi=args.dpi)
        print(f"Saved animation to: {out_path}")
    except Exception as e:
        print(f"Primary save failed: {e}")
        fallback = out_path.with_suffix(".gif")
        print(f"Trying GIF fallback: {fallback}")
        ani.save(fallback, writer="pillow", fps=args.fps, dpi=args.dpi)
        print(f"Saved animation to: {fallback}")

    plt.close(fig)


if __name__ == "__main__":
    main()