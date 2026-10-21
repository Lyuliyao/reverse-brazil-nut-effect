#!/usr/bin/env python3
"""Reproduce selected numerical figures using repository-relative data paths."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figures", nargs="+", type=int, choices=range(6, 11), default=list(range(6, 11)))
    parser.add_argument("--out", type=Path, default=root / "outputs" / "figures")
    args = parser.parse_args()
    pythonpath = os.pathsep.join(filter(None, [str(root / "analysis"), os.environ.get("PYTHONPATH")]))
    env = dict(os.environ, MPLBACKEND="Agg", RBNE_OUTPUT=str(args.out.resolve()), PYTHONPATH=pythonpath)
    for number in args.figures:
        subprocess.run([sys.executable, str(root / "analysis" / f"figure{number:02}.py")],
                       env=env, check=True, cwd=root)

if __name__ == "__main__":
    main()
