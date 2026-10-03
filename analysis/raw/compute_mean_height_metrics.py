#!/usr/bin/env python3
"""Compute species mean-height time series and RBNE/BNE metrics for one case.

Sign convention (project-wide): z is positive upward (wall-normal to the
incline); dz = zbar_L - zbar_S.  dz < 0 means large particles are below small
particles, i.e. transient RBNE.  dz > 0 is conventional BNE.

Outputs
-------
* per-frame CSV (default results/<case>_meanz.csv) with columns:
    time, zbar_small, zbar_large, dz, vz_small, vz_large,
    cc_small, cc_large, n_small, n_large
  (vz_* are NaN when the dump lacks velocities; cc_* are NaN when c_cc absent.)
* metrics JSON (results/<case>_meanz_metrics.json) with:
    dz0, early_slope, A_RBNE, t_min, t_cross, dz_final, R_BNE, plus bookkeeping.

t_cross is, by definition here, the FIRST frame time with dz >= 0 that occurs
strictly AFTER t_min.  We never read any precomputed cross field from elsewhere.

This script does NOT require velocity or contact data: vz columns are reported
when present but are not needed for the height metrics.
"""

import argparse
import csv
import json
from pathlib import Path

import numpy as np

import io_common


def _safe_mean(values, mask):
    sel = values[mask]
    if sel.size == 0:
        return float("nan")
    return float(np.mean(sel))


def compute_series(case_dir):
    """Return a list of per-frame dicts (ascending time)."""
    rows = []
    for time, cols, _box in io_common.iter_frames(case_dir):
        mask_s, mask_l, _rs, _rl = io_common.species_masks(cols)
        z = cols.get("z")
        if z is None:
            raise ValueError("Frame has no 'z' column; cannot compute mean heights")

        vz = cols.get("vz")
        cc = cols.get("c_cc")

        rows.append({
            "time": float(time),
            "zbar_small": _safe_mean(z, mask_s),
            "zbar_large": _safe_mean(z, mask_l),
            "dz": _safe_mean(z, mask_l) - _safe_mean(z, mask_s),
            "vz_small": _safe_mean(vz, mask_s) if vz is not None else float("nan"),
            "vz_large": _safe_mean(vz, mask_l) if vz is not None else float("nan"),
            "cc_small": _safe_mean(cc, mask_s) if cc is not None else float("nan"),
            "cc_large": _safe_mean(cc, mask_l) if cc is not None else float("nan"),
            "n_small": int(np.count_nonzero(mask_s)),
            "n_large": int(np.count_nonzero(mask_l)),
        })
    return rows


def _linfit_slope(t, y):
    """Slope of a least-squares line y ~ a + b t, ignoring NaNs. NaN if <2 pts."""
    t = np.asarray(t, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    good = np.isfinite(t) & np.isfinite(y)
    if np.count_nonzero(good) < 2:
        return float("nan")
    tt = t[good]
    yy = y[good]
    if np.allclose(tt, tt[0]):
        return float("nan")
    b = np.polyfit(tt, yy, 1)[0]
    return float(b)


def compute_metrics(rows):
    """Compute RBNE/BNE metrics from a per-frame series.

    early_slope: linear fit of dz over the first ~10% of the time span OR the
    first 2 s, whichever covers more frames (with at least 2 points).
    """
    if not rows:
        return {
            "dz0": None, "early_slope": None, "A_RBNE": None,
            "t_min": None, "t_cross": None, "dz_final": None,
            "R_BNE": None, "n_frames": 0,
        }

    times = np.array([r["time"] for r in rows], dtype=np.float64)
    dz = np.array([r["dz"] for r in rows], dtype=np.float64)

    t0 = float(times[0])
    T = float(times[-1])
    span = T - t0

    dz0 = float(dz[0])
    dz_final = float(dz[-1])

    # Early window: first max(10% span, 2 s).
    early_cut = t0 + max(0.1 * span, 2.0)
    early_mask = times <= early_cut
    if np.count_nonzero(early_mask) < 2:
        # Fall back to the first two frames.
        early_mask = np.zeros_like(times, dtype=bool)
        early_mask[:2] = True
    early_slope = _linfit_slope(times[early_mask], dz[early_mask])

    # A_RBNE and t_min from the global minimum of dz.
    finite = np.isfinite(dz)
    if not np.any(finite):
        return {
            "dz0": dz0, "early_slope": early_slope, "A_RBNE": None,
            "t_min": None, "t_cross": None, "dz_final": dz_final,
            "R_BNE": None, "n_frames": len(rows),
        }
    finite_idx = np.where(finite)[0]
    imin_local = int(np.argmin(dz[finite_idx]))
    imin = int(finite_idx[imin_local])
    dz_min = float(dz[imin])
    t_min = float(times[imin])
    A_RBNE = max(0.0, -dz_min)

    # t_cross: FIRST frame time with dz >= 0 strictly AFTER t_min.
    t_cross = None
    for k in range(imin + 1, len(rows)):
        if np.isfinite(dz[k]) and dz[k] >= 0.0:
            t_cross = float(times[k])
            break

    # R_BNE: slope of dz over [t_min, min(t_min + 10 s, T)].
    rec_hi = min(t_min + 10.0, T)
    rec_mask = (times >= t_min) & (times <= rec_hi)
    R_BNE = _linfit_slope(times[rec_mask], dz[rec_mask])

    return {
        "dz0": dz0,
        "early_slope": early_slope,
        "A_RBNE": float(A_RBNE),
        "t_min": t_min,
        "t_cross": t_cross,
        "dz_final": dz_final,
        "R_BNE": R_BNE,
        "n_frames": len(rows),
        "t_start": t0,
        "t_end": T,
        "dz_min": dz_min,
    }


CSV_FIELDS = [
    "time", "zbar_small", "zbar_large", "dz",
    "vz_small", "vz_large", "cc_small", "cc_large",
    "n_small", "n_large",
]


def write_csv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row[k] for k in CSV_FIELDS})


def main():
    parser = argparse.ArgumentParser(
        description="Mean-height time series and RBNE/BNE metrics for one case."
    )
    parser.add_argument("--case", required=True, help="Path to a case directory")
    parser.add_argument(
        "--out", default=None,
        help="Per-frame CSV path (default results/<case>_meanz.csv)",
    )
    args = parser.parse_args()

    case_dir = Path(args.case)
    if not case_dir.is_dir():
        raise SystemExit(f"error: case dir not found: {case_dir}")

    metadata = io_common.read_metadata(case_dir)
    name = metadata.get("case_name", case_dir.name)

    analysis_root = Path(__file__).resolve().parent
    out_csv = Path(args.out) if args.out else analysis_root / "results" / f"{name}_meanz.csv"
    out_metrics = analysis_root / "results" / f"{name}_meanz_metrics.json"

    rows = compute_series(case_dir)
    if not rows:
        raise SystemExit(f"error: no particle frames found under {case_dir}")

    write_csv(out_csv, rows)
    metrics = compute_metrics(rows)
    metrics.update({
        "case_name": name,
        "case_path": str(case_dir.resolve()),
        "control_type": metadata.get("control_type"),
        "csv_path": str(out_csv.resolve()),
    })

    out_metrics.parent.mkdir(parents=True, exist_ok=True)
    out_metrics.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"case: {name}")
    print(f"frames: {metrics['n_frames']}")
    print(f"dz0={metrics['dz0']}  A_RBNE={metrics['A_RBNE']}  "
          f"t_min={metrics['t_min']}  t_cross={metrics['t_cross']}  "
          f"dz_final={metrics['dz_final']}  R_BNE={metrics['R_BNE']}")
    print(f"per-frame CSV: {out_csv}")
    print(f"metrics JSON:  {out_metrics}")


if __name__ == "__main__":
    main()
