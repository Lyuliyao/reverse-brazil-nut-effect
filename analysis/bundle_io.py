"""Portable readers for the paper's compact data bundle (all stored units SI)."""
import os
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
CASES = REPO / "data" / "cases"
OUTPUT = Path(os.environ.get("RBNE_OUTPUT", str(REPO / "outputs" / "figures"))).resolve()
OUTPUT.mkdir(parents=True, exist_ok=True)

def load_snapshot(case, target):
    paths = list((Path(case) / "snapshots").glob("t*.npz"))
    if not paths:
        raise FileNotFoundError(f"No bundled snapshots in {case}")
    selected = min(paths, key=lambda p: abs(float(p.stem[1:]) - target))
    with np.load(selected, allow_pickle=False) as f:
        return f["xyz"], f["types"], f["radii"], float(f["time"]), f["box"]

def read_height_series(case):
    df = pd.read_csv(Path(case) / "time_series.tsv", sep="\t").sort_values("time")
    return df["time"].to_numpy(), df["z_gap_large_minus_small"].to_numpy() * 1e3

def read_fluid_profile(case, target):
    paths = list(Path(case).glob("fluid_profile_*.csv"))
    if not paths:
        raise FileNotFoundError(f"No bundled fluid profiles in {case}")
    selected = min(paths, key=lambda p: abs(float(p.stem.split("_")[-1]) - target))
    df = pd.read_csv(selected)
    return (float(selected.stem.split("_")[-1]), df["z_m"].to_numpy() * 1e3,
            df["ux_m_s"].to_numpy() * 1e3, df["uz_m_s"].to_numpy() * 1e3)
