#!/usr/bin/env python3
"""Prepare a fresh run from an archived paper case, without altering the archive."""
import argparse
import json
from pathlib import Path
import re
import shutil

ROOT = Path(__file__).resolve().parents[1]

def replace_entry(text, key, value):
    text, count = re.subn(rf"(?m)^(\s*){key}\s+[^;]+;", rf"\g<1>{key} {value};", text)
    if count != 1:
        raise ValueError(f"Expected one {key} entry; found {count}")
    return text

def prepare_case(name, dest, end_time=40., nprocs=8):
    source = ROOT / "cases" / name
    if not source.is_dir() or source.parent != ROOT / "cases":
        raise ValueError(f"Unknown case: {name}")
    dest = Path(dest).resolve()
    if dest.exists():
        raise FileExistsError(f"Run directory already exists: {dest}")
    if end_time <= 0 or nprocs < 1:
        raise ValueError("end_time and nprocs must be positive")
    for field in ("U", "Us", "rho", "p", "Ksl", "voidfraction"):
        if not (source / "CFD" / "0" / field).is_file():
            raise FileNotFoundError(f"Missing initial field: {field}")
    shutil.copytree(source, dest)
    for relative in ("DEM/post/restart", "DEM/post/reports",
                     "DEM/post/forcebalance_segments/t0_to_40",
                     "DEM/post/contact_segments/t0_to_40", "CFD/couplingFiles"):
        (dest / relative).mkdir(parents=True, exist_ok=True)
    control = dest / "CFD/system/controlDict"
    text = control.read_text()
    for key, value in (("startFrom", "startTime"), ("startTime", "0"),
                       ("endTime", f"{end_time:g}"), ("writeInterval", "0.5"),
                       ("purgeWrite", "0")):
        text = replace_entry(text, key, value)
    control.write_text(text)
    decompose = dest / "CFD/system/decomposeParDict"
    decompose.write_text(replace_entry(decompose.read_text(), "numberOfSubdomains", str(nprocs)))
    resume = dest / "DEM/in.liggghts_resume"
    text = resume.read_text()
    text, count = re.subn(r"(?m)^read_restart\s+.*$",
                          "read_restart ../DEM/post/restart/liggghts.restart", text)
    if count != 1:
        raise ValueError("Expected one DEM restart input")
    text = re.sub(r"forcebalance_segments/t[^/\s]+/", "forcebalance_segments/t0_to_40/", text)
    text = re.sub(r"contact_segments/t[^/\s]+/", "contact_segments/t0_to_40/", text)
    resume.write_text(text)
    parameters = json.loads((dest / "case_params.json").read_text())
    parameters.update(end_time=end_time, nprocs=nprocs, purge_write=0,
                      forcebalance_segment_label="t0_to_40",
                      forcebalance_segment_dir=str(dest / "DEM/post/forcebalance_segments/t0_to_40"))
    (dest / "case_params.json").write_text(json.dumps(parameters, indent=2) + "\n")
    record = {"source_case": name, "end_time_s": end_time, "nprocs": nprocs,
              "mode": "fresh initialization from archived DEM insertion seeds; no restart binary copied"}
    (dest / "prepared_case.json").write_text(json.dumps(record, indent=2) + "\n")
    return dest

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", required=True, choices=sorted(p.name for p in (ROOT / "cases").iterdir()))
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--end-time", type=float, default=40.)
    parser.add_argument("--nprocs", type=int, default=8)
    args = parser.parse_args()
    print(prepare_case(args.case, args.out, args.end_time, args.nprocs))

if __name__ == "__main__":
    main()
