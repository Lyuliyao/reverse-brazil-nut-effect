#!/usr/bin/env python3
"""Regenerate provenance/checksums.json for the research artifacts.

Only the data, simulation inputs and solver patches are
checksummed. Documentation and code are versioned by git and may be edited
without touching the manifest. Run this after intentionally changing a file
under one of the CHECKED prefixes, then commit the updated manifest.
"""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHECKED = ("data/", "cases/", "simulation/patches/")
MANIFEST = ROOT / "provenance" / "checksums.json"


def tracked_files():
    out = subprocess.run(["git", "ls-files", "-z"], cwd=ROOT, check=True,
                         capture_output=True).stdout.decode("utf-8")
    return sorted(f for f in out.split("\0") if f.startswith(CHECKED))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true",
                        help="only report differences; do not write the manifest")
    args = parser.parse_args()
    new = {f: hashlib.sha256((ROOT / f).read_bytes()).hexdigest() for f in tracked_files()}
    old = json.loads(MANIFEST.read_text()) if MANIFEST.exists() else {}
    changed = sorted(k for k in new if old.get(k) != new[k])
    removed = sorted(set(old) - set(new))
    print(f"{len(new)} files checksummed; {len(changed)} new/changed, {len(removed)} dropped")
    if args.check:
        raise SystemExit(1 if changed or removed else 0)
    MANIFEST.write_text(json.dumps(new, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
