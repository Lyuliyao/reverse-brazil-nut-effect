#!/usr/bin/env python
"""Consolidate per-case analysis/z_profiles/*/dump*.z_profile.tsv into a single HDF5.

For each case directory passed on the command line, reads every
``z_profiles/<segment>/dump<step>.z_profile.tsv`` and writes a single
``z_profiles.h5`` next to the ``z_profiles/`` folder.

HDF5 layout::

    /profile
        attrs["columns"] = list of 7 numeric column names
        attrs["n_rows"]  = total number of rows
        /segment   (N,) variable-length strings
        /dump_step (N,) int64
        /values    (N, 7) float64

Does NOT delete the originals. Verifies row counts against the source files
and aborts (per case) if anything looks off so that nothing gets quietly lost.
"""

from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

NUMERIC_COLS = [
    "bin_index",
    "z_center",
    "phi_small",
    "phi_large",
    "phi_total",
    "count_small",
    "count_large",
]

DUMP_RE = re.compile(r"^dump(\d+)\.z_profile\.tsv$")


def pack_case(case_dir: Path, out_path: Path | None = None, overwrite: bool = False) -> dict:
    zp_dir = case_dir / "analysis" / "z_profiles"
    if not zp_dir.is_dir():
        raise FileNotFoundError(f"no z_profiles dir under {case_dir}")

    if out_path is None:
        out_path = case_dir / "analysis" / "z_profiles.h5"
    if out_path.exists() and not overwrite:
        raise FileExistsError(f"{out_path} already exists; pass --overwrite to replace")

    segments = sorted(p for p in zp_dir.iterdir() if p.is_dir())
    if not segments:
        raise RuntimeError(f"no segment dirs inside {zp_dir}")

    seg_list: list[str] = []
    step_list: list[int] = []
    value_blocks: list[np.ndarray] = []
    file_count = 0

    expected_cols = NUMERIC_COLS  # exact column order to enforce

    for seg_dir in segments:
        files = sorted(seg_dir.glob("dump*.z_profile.tsv"))
        for f in files:
            m = DUMP_RE.match(f.name)
            if not m:
                raise RuntimeError(f"unexpected filename: {f}")
            step = int(m.group(1))
            df = pd.read_csv(f, sep="\t")
            if list(df.columns) != expected_cols:
                raise RuntimeError(
                    f"column mismatch in {f}: got {list(df.columns)}, "
                    f"expected {expected_cols}"
                )
            value_blocks.append(df.to_numpy(dtype=np.float64))
            n = len(df)
            seg_list.extend([seg_dir.name] * n)
            step_list.extend([step] * n)
            file_count += 1

    if not value_blocks:
        raise RuntimeError(f"no dump*.z_profile.tsv files found under {zp_dir}")

    values = np.vstack(value_blocks)
    segs = np.array(seg_list, dtype=object)
    steps = np.array(step_list, dtype=np.int64)
    assert values.shape == (len(segs), len(NUMERIC_COLS))
    assert len(segs) == len(steps)

    tmp = out_path.with_suffix(out_path.suffix + ".tmp")
    if tmp.exists():
        tmp.unlink()
    with h5py.File(tmp, "w") as h:
        g = h.create_group("profile")
        g.attrs["columns"] = NUMERIC_COLS
        g.attrs["n_rows"] = int(values.shape[0])
        g.attrs["n_source_files"] = int(file_count)
        g.attrs["source_dir"] = str(zp_dir)
        g.create_dataset(
            "values",
            data=values,
            compression="gzip",
            compression_opts=4,
            shuffle=True,
        )
        g.create_dataset(
            "dump_step",
            data=steps,
            compression="gzip",
            compression_opts=4,
        )
        g.create_dataset(
            "segment",
            data=segs.astype(h5py.string_dtype()),
            compression="gzip",
            compression_opts=4,
        )

    # round-trip read-back check
    with h5py.File(tmp, "r") as h:
        rb_values = h["profile/values"][...]
        rb_steps = h["profile/dump_step"][...]
        rb_segs = h["profile/segment"].asstr()[...]
        if rb_values.shape != values.shape:
            raise RuntimeError("readback shape mismatch")
        if not np.array_equal(rb_steps, steps):
            raise RuntimeError("readback dump_step mismatch")
        if not np.array_equal(rb_segs, segs.astype(str)):
            raise RuntimeError("readback segment mismatch")
        # check a few random rows match exactly
        rng = np.random.default_rng(0)
        idx = rng.choice(values.shape[0], size=min(64, values.shape[0]), replace=False)
        if not np.allclose(rb_values[idx], values[idx], equal_nan=True):
            raise RuntimeError("readback value mismatch")

    tmp.replace(out_path)

    return {
        "case": case_dir.name,
        "out": str(out_path),
        "n_files": file_count,
        "n_rows": int(values.shape[0]),
        "n_segments": len(segments),
        "size_bytes": out_path.stat().st_size,
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "cases",
        nargs="+",
        type=Path,
        help="case directories (each must contain analysis/z_profiles/)",
    )
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="overwrite existing z_profiles.h5",
    )
    p.add_argument(
        "--quiet",
        action="store_true",
        help="only print errors and final summary",
    )
    args = p.parse_args()

    total_files = 0
    total_rows = 0
    total_size = 0
    failed: list[tuple[str, str]] = []
    for case in args.cases:
        t0 = time.time()
        try:
            info = pack_case(case, overwrite=args.overwrite)
        except Exception as e:  # noqa: BLE001
            failed.append((str(case), repr(e)))
            print(f"[FAIL] {case}: {e}", file=sys.stderr)
            continue
        dt = time.time() - t0
        total_files += info["n_files"]
        total_rows += info["n_rows"]
        total_size += info["size_bytes"]
        if not args.quiet:
            mb = info["size_bytes"] / (1024 * 1024)
            print(
                f"[ok] {info['case']}: {info['n_files']} files, "
                f"{info['n_rows']} rows, {info['n_segments']} segs, "
                f"{mb:.1f} MiB, {dt:.1f}s -> {info['out']}"
            )

    print(
        f"\nDONE: packed {total_files} files / {total_rows} rows / "
        f"{total_size / (1024 * 1024):.1f} MiB across {len(args.cases) - len(failed)} cases"
    )
    if failed:
        print(f"FAILED cases: {len(failed)}", file=sys.stderr)
        for c, e in failed:
            print(f"  {c}: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
