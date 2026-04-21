#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import shutil
import sys
from pathlib import Path
import re

import h5py


ITER_RE = re.compile(r"\bit=(\d+)\b")


def parse_args():
    p = argparse.ArgumentParser(
        description="Extract one iteration from FIL/Carpet split 3D HDF5 files."
    )
    p.add_argument("--input-dir", required=True,
                   help="Directory with original .h5 files")
    p.add_argument("--output-dir", required=True,
                   help="Directory for reduced .h5 files")
    p.add_argument("--iteration", type=int,
                   help="Iteration to keep, e.g. 196608")
    p.add_argument("--pattern", default="*.h5",
                   help="Glob pattern for files (default: *.h5)")
    p.add_argument("--overwrite", action="store_true",
                   help="Overwrite output files")
    p.add_argument("--dry-run", action="store_true",
                   help="Do not write, only print actions")
    p.add_argument("--list-iterations", action="store_true",
                   help="List iterations found and exit")
    p.add_argument("--max-files-list", type=int, default=10,
                   help="How many files to inspect for --list-iterations")
    p.add_argument("--verbose", action="store_true",
                   help="Verbose output")
    return p.parse_args()


def extract_iteration_from_name(name):
    m = ITER_RE.search(name)
    if m:
        return int(m.group(1))
    return None


def should_keep_object(name, obj, selected_it):
    """
    Keep:
      - datasets whose name contains it=<selected_it>
      - the group 'Parameters and Global Attributes'
    """
    if isinstance(obj, h5py.Group):
        if name == "Parameters and Global Attributes":
            return True
        return False

    if isinstance(obj, h5py.Dataset):
        it = extract_iteration_from_name(name)
        return it == selected_it

    return False


def copy_attrs(src_obj, dst_obj):
    for k, v in src_obj.attrs.items():
        try:
            dst_obj.attrs[k] = v
        except Exception:
            pass


def list_iterations_in_file(h5file):
    found = set()
    with h5py.File(h5file, "r") as f:
        for name in f.keys():
            it = extract_iteration_from_name(name)
            if it is not None:
                found.add(it)
    return found


def copy_selected_content(src_path, dst_path, selected_it, dry_run=False, verbose=False):
    n_copied = 0

    with h5py.File(src_path, "r") as src:
        if dry_run:
            print(f"[DRY] {src_path.name}")
            for name in src.keys():
                obj = src[name]
                if should_keep_object(name, obj, selected_it):
                    print(f"  keep: {name}")
            return 0

        with h5py.File(dst_path, "w") as dst:
            copy_attrs(src, dst)

            for name in src.keys():
                obj = src[name]
                if should_keep_object(name, obj, selected_it):
                    src.copy(name, dst)
                    n_copied += 1
                    if verbose:
                        print(f"    copied: {name}")

            # Also copy attrs of copied top-level objects
            for name in dst.keys():
                copy_attrs(src[name], dst[name])

    return n_copied


def main():
    args = parse_args()

    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()

    if not input_dir.is_dir():
        print(f"ERROR: input directory does not exist: {input_dir}", file=sys.stderr)
        sys.exit(1)

    files = sorted(input_dir.glob(args.pattern))
    if not files:
        print(f"ERROR: no files found in {input_dir} matching {args.pattern}", file=sys.stderr)
        sys.exit(1)

    if args.list_iterations:
        all_iters = set()
        n = min(args.max_files_list, len(files))
        print(f"Scanning up to {n} files for iterations...")
        for f in files[:n]:
            try:
                its = list_iterations_in_file(f)
                print(f"  {f.name}: {sorted(its)}")
                all_iters.update(its)
            except Exception as e:
                print(f"  WARNING: failed on {f.name}: {e}")
        print("Discovered iterations:")
        print(" ".join(str(it) for it in sorted(all_iters)))
        return

    if args.iteration is None:
        print("ERROR: please provide --iteration", file=sys.stderr)
        sys.exit(1)

    if output_dir.exists():
        if not output_dir.is_dir():
            print(f"ERROR: output path exists and is not a directory: {output_dir}", file=sys.stderr)
            sys.exit(1)
        if any(output_dir.iterdir()) and not args.overwrite:
            print(f"ERROR: output directory already contains files: {output_dir}", file=sys.stderr)
            print("Use --overwrite, or choose another --output-dir", file=sys.stderr)
            sys.exit(1)
        if args.overwrite:
            for p in output_dir.iterdir():
                if p.is_file() or p.is_symlink():
                    p.unlink()
                elif p.is_dir():
                    shutil.rmtree(p)
    else:
        output_dir.mkdir(parents=True)

    total = len(files)
    copied_files = 0
    empty_files = 0
    failed_files = 0

    print(f"Input directory : {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Iteration       : {args.iteration}")
    print(f"Files           : {total}")
    print()

    for i, f in enumerate(files, start=1):
        out = output_dir / f.name
        print(f"[{i}/{total}] {f.name}")
        try:
            n = copy_selected_content(
                f, out, args.iteration,
                dry_run=args.dry_run,
                verbose=args.verbose
            )
            if not args.dry_run:
                if n > 0:
                    copied_files += 1
                    print(f"  -> copied {n} objects")
                else:
                    empty_files += 1
                    print("  -> no matching datasets")
        except Exception as e:
            failed_files += 1
            print(f"  -> ERROR: {e}")

    print()
    print("Done.")
    print(f"  copied files : {copied_files}")
    print(f"  empty files  : {empty_files}")
    print(f"  failed files : {failed_files}")


if __name__ == "__main__":
    main()