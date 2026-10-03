#!/usr/bin/env python
"""Sample-verify packed z_profiles.h5 against original *.z_profile.tsv files.

For each case:
  * count tsv files under analysis/z_profiles/ and rows in z_profiles.h5;
    the row total must equal sum(rows in every tsv).
  * pick N random tsv files, locate the rows in the h5 by
    (segment, dump_step), compare values bit-for-bit.
Exits non-zero if any case fails.
"""

from __future__ import annotations

import argparse
import re
import sys
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


def verify_case(case_dir: Path, n_samples: int, rng: np.random.Generator) -> tuple[bool, str]:
    zp = case_dir / "analysis" / "z_profiles"
    h5p = case_dir / "analysis" / "z_profiles.h5"
    if not zp.is_dir():
        return False, f"no z_profiles dir: {zp}"
    if not h5p.exists():
        return False, f"no h5 file: {h5p}"

    tsvs: list[tuple[str, int, Path]] = []
    for seg in sorted(p for p in zp.iterdir() if p.is_dir()):
        for f in seg.glob("dump*.z_profile.tsv"):
            m = DUMP_RE.match(f.name)
            if not m:
                return False, f"bad name: {f}"
            tsvs.append((seg.name, int(m.group(1)), f))
    if not tsvs:
        return False, f"no tsv files under {zp}"

    with h5py.File(h5p, "r") as h:
        g = h["profile"]
        cols = [c.decode() if isinstance(c, bytes) else c for c in list(g.attrs["columns"])]
        if cols != NUMERIC_COLS:
            return False, f"h5 column mismatch: {cols}"
        values = g["values"][...]
        steps = g["dump_step"][...]
        segs = g["segment"].asstr()[...]

    # row-count check: total tsv lines (minus 1 header each) == n_rows in h5
    total_tsv_rows = 0
    # only sum line counts from sampled set + a fast estimate via first file's bin count
    # to keep this cheap, assume all files share same bin count -> verify against first
    ref_rows = sum(1 for _ in open(tsvs[0][2])) - 1
    # but bins per file may differ between cases, so count for THIS case using first file's count assumption,
    # then we still verify per-sample rows below
    expected_total = ref_rows * len(tsvs)
    if expected_total != values.shape[0]:
        # fall back to exact count
        actual = 0
        for _, _, f in tsvs:
            actual += sum(1 for _ in open(f)) - 1
        if actual != values.shape[0]:
            return False, f"row total mismatch: tsv={actual} h5={values.shape[0]}"
        total_tsv_rows = actual
    else:
        total_tsv_rows = expected_total

    # spot-check N random tsv files
    n = min(n_samples, len(tsvs))
    idx = rng.choice(len(tsvs), size=n, replace=False)
    for i in idx:
        seg, step, f = tsvs[i]
        src = pd.read_csv(f, sep="\t")
        if list(src.columns) != NUMERIC_COLS:
            return False, f"column mismatch in {f}"
        mask = (segs == seg) & (steps == step)
        m_idx = np.where(mask)[0]
        if m_idx.size == 0:
            return False, f"missing in h5: seg={seg} step={step}"
        if m_idx.size != len(src):
            return False, f"row count mismatch for seg={seg} step={step}: tsv={len(src)} h5={m_idx.size}"
        # h5 rows are stored in the order they were read; bin_index column is monotonic
        h5_block = values[m_idx]
        src_arr = src[NUMERIC_COLS].to_numpy(dtype=np.float64)
        if not np.array_equal(h5_block, src_arr):
            # find first mismatch for diagnostics
            diff = np.where(~np.isclose(h5_block, src_arr, atol=0, rtol=0))
            r, c = int(diff[0][0]), int(diff[1][0])
            return False, (
                f"value mismatch seg={seg} step={step} at row={r} col={NUMERIC_COLS[c]}: "
                f"tsv={src_arr[r, c]} h5={h5_block[r, c]}"
            )

    return True, f"{len(tsvs)} tsv / {total_tsv_rows} rows, sampled {n}, all match"


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cases", nargs="+", type=Path)
    p.add_argument("--samples", type=int, default=8, help="random tsv files to verify per case")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--quiet", action="store_true")
    args = p.parse_args()
    rng = np.random.default_rng(args.seed)

    fails: list[tuple[str, str]] = []
    for c in args.cases:
        ok, msg = verify_case(c, args.samples, rng)
        tag = "ok" if ok else "FAIL"
        if not ok or not args.quiet:
            print(f"[{tag}] {c.name}: {msg}")
        if not ok:
            fails.append((str(c), msg))
    print(f"\nverified {len(args.cases) - len(fails)}/{len(args.cases)} cases ok")
    if fails:
        print("FAILURES:", file=sys.stderr)
        for c, m in fails:
            print(f"  {c}: {m}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
