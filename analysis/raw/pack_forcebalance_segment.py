#!/usr/bin/env python3
import argparse
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import h5py
import numpy as np


DUMP_RE = re.compile(r"^dump([0-9]+)\.liggghts_forcebalance$")
SEGMENT_RE = re.compile(r"^t([0-9]+)_to_([0-9]+)$")
DEFAULT_H5_NAME = "forcebalance.h5"


def natural_dump_key(path):
    match = DUMP_RE.match(path.name)
    if match:
        return int(match.group(1)), path.name
    return sys.maxsize, path.name


def segment_key(path):
    match = SEGMENT_RE.match(path.name)
    if match:
        return int(match.group(1)), int(match.group(2)), path.name
    return -1, -1, path.name


def is_latest_segment(segment_dir):
    siblings = [path for path in segment_dir.parent.iterdir() if path.is_dir() and SEGMENT_RE.match(path.name)]
    if not siblings:
        return False
    return segment_dir == max(siblings, key=segment_key)


def decode_attr(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def read_dump(path):
    timestep = None
    n_atoms = None
    box_bounds = []
    columns = None
    rows = []

    with path.open("r", encoding="utf-8") as handle:
        line_iter = iter(handle)
        for raw_line in line_iter:
            line = raw_line.strip()
            if line == "ITEM: TIMESTEP":
                try:
                    timestep = int(next(line_iter).strip())
                except StopIteration as exc:
                    raise ValueError(f"Unexpected EOF after TIMESTEP in {path}") from exc
                continue
            if line == "ITEM: NUMBER OF ATOMS":
                try:
                    n_atoms = int(next(line_iter).strip())
                except StopIteration as exc:
                    raise ValueError(f"Unexpected EOF after NUMBER OF ATOMS in {path}") from exc
                continue
            if line.startswith("ITEM: BOX BOUNDS"):
                box_bounds = []
                for _ in range(3):
                    try:
                        bounds = np.fromstring(next(line_iter).strip(), sep=" ", dtype=np.float64)
                    except StopIteration as exc:
                        raise ValueError(f"Unexpected EOF in BOX BOUNDS in {path}") from exc
                    if bounds.size < 2:
                        raise ValueError(f"Malformed BOX BOUNDS row in {path}")
                    box_bounds.append(bounds[:2])
                continue
            if line.startswith("ITEM: ATOMS "):
                columns = line.split()[2:]
                if not columns:
                    raise ValueError(f"Missing ATOMS columns in {path}")
                if n_atoms is None:
                    raise ValueError(f"ATOMS section appeared before atom count in {path}")
                rows = np.empty((n_atoms, len(columns)), dtype=np.float64)
                for index in range(n_atoms):
                    try:
                        atom_line = next(line_iter)
                    except StopIteration as exc:
                        raise ValueError(f"Unexpected EOF in ATOMS data in {path}") from exc
                    values = np.fromstring(atom_line, sep=" ", dtype=np.float64)
                    if values.size != len(columns):
                        raise ValueError(
                            f"Expected {len(columns)} ATOMS values, got {values.size} in {path} line {index + 1}"
                        )
                    rows[index, :] = values
                break

    if timestep is None or n_atoms is None or columns is None:
        raise ValueError(f"Failed to parse required dump sections in {path}")
    if len(box_bounds) != 3:
        raise ValueError(f"Failed to parse 3 BOX BOUNDS rows in {path}")
    return timestep, np.asarray(box_bounds, dtype=np.float64), columns, rows


def h5_columns(h5_path):
    with h5py.File(h5_path, "r") as h5_file:
        if "particles" not in h5_file or "frames" not in h5_file:
            raise ValueError(f"{h5_path} is missing required datasets")
        raw_columns = h5_file["particles"].attrs.get("columns")
        if raw_columns is None:
            raise ValueError(f"{h5_path} is missing particles columns metadata")
        return [decode_attr(value) for value in raw_columns]


def validate_h5(h5_path, expected_frames=None, expected_rows=None, expected_sources=None):
    with h5py.File(h5_path, "r") as h5_file:
        if "particles" not in h5_file or "frames" not in h5_file or "box_bounds" not in h5_file:
            raise ValueError(f"{h5_path} is missing one of particles, frames, box_bounds")
        particles = h5_file["particles"]
        frames = h5_file["frames"]
        box_bounds = h5_file["box_bounds"]
        if particles.ndim != 2:
            raise ValueError(f"{h5_path}: particles dataset must be 2-D")
        if len(frames) != len(box_bounds):
            raise ValueError(f"{h5_path}: frames and box_bounds lengths differ")
        if expected_frames is not None and len(frames) != expected_frames:
            raise ValueError(f"{h5_path}: expected {expected_frames} frames, found {len(frames)}")
        if expected_sources is not None:
            actual_sources = [decode_attr(frame["source_name"]) for frame in frames]
            if actual_sources != list(expected_sources):
                raise ValueError(f"{h5_path}: frame source names do not match current raw dumps")
        frame_rows = int(np.sum(frames["row_count"])) if len(frames) else 0
        if particles.shape[0] != frame_rows:
            raise ValueError(f"{h5_path}: particles rows {particles.shape[0]} != frame row total {frame_rows}")
        if expected_rows is not None and particles.shape[0] != expected_rows:
            raise ValueError(f"{h5_path}: expected {expected_rows} particle rows, found {particles.shape[0]}")
        if particles.shape[1] != len(particles.attrs.get("columns", [])):
            raise ValueError(f"{h5_path}: particle column count metadata mismatch")
    return True


def write_h5(segment_dir, dump_paths, h5_path, force):
    tmp_path = h5_path.with_name(f".{h5_path.name}.tmp.{os.getpid()}")
    if tmp_path.exists():
        tmp_path.unlink()

    columns = None
    frame_rows = []
    total_rows = 0
    string_dtype = h5py.string_dtype(encoding="utf-8", length=128)
    frame_dtype = np.dtype(
        [
            ("timestep", "<i8"),
            ("row_start", "<i8"),
            ("row_count", "<i8"),
            ("source_name", string_dtype),
        ]
    )

    try:
        with h5py.File(tmp_path, "w") as h5_file:
            h5_file.attrs["format"] = "bne_forcebalance_segment_v1"
            h5_file.attrs["segment_label"] = segment_dir.name
            h5_file.attrs["source_dir"] = str(segment_dir)
            h5_file.attrs["created_unix_time"] = int(time.time())
            h5_file.attrs["created_by"] = Path(__file__).name
            h5_file.attrs["source_dump_count"] = len(dump_paths)

            particles = None
            box_bounds_dataset = h5_file.create_dataset(
                "box_bounds",
                shape=(0, 3, 2),
                maxshape=(None, 3, 2),
                chunks=(1, 3, 2),
                dtype=np.float64,
            )

            for dump_index, dump_path in enumerate(dump_paths):
                timestep, box_bounds, dump_columns, rows = read_dump(dump_path)
                if columns is None:
                    columns = dump_columns
                    row_chunk = max(1024, min(65536, max(rows.shape[0], 1)))
                    particles = h5_file.create_dataset(
                        "particles",
                        shape=(0, len(columns)),
                        maxshape=(None, len(columns)),
                        chunks=(row_chunk, len(columns)),
                        dtype=np.float64,
                        compression="lzf",
                        shuffle=True,
                    )
                    particles.attrs["columns"] = np.asarray(columns, dtype=h5py.string_dtype(encoding="utf-8"))
                elif dump_columns != columns:
                    raise ValueError(f"{dump_path} columns differ from first dump")

                row_start = total_rows
                row_count = rows.shape[0]
                particles.resize((total_rows + row_count, len(columns)))
                particles[row_start : row_start + row_count, :] = rows

                box_bounds_dataset.resize((dump_index + 1, 3, 2))
                box_bounds_dataset[dump_index, :, :] = box_bounds
                frame_rows.append((timestep, row_start, row_count, dump_path.name))
                total_rows += row_count

            if particles is None:
                raise ValueError(f"No dump rows found in {segment_dir}")

            frames = h5_file.create_dataset("frames", data=np.asarray(frame_rows, dtype=frame_dtype))
            frames.attrs["columns"] = np.asarray(
                ["timestep", "row_start", "row_count", "source_name"],
                dtype=h5py.string_dtype(encoding="utf-8"),
            )

        validate_h5(tmp_path, expected_frames=len(dump_paths), expected_rows=total_rows)
        if h5_path.exists() and force:
            h5_path.unlink()
        tmp_path.replace(h5_path)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()

    return total_rows


def archive_segment(segment_dir, archive_path, force, h5_name, expected_dump_count, expected_dump_names):
    tmp_archive = archive_path.with_name(f".{archive_path.name}.tmp.{os.getpid()}")
    if tmp_archive.exists():
        tmp_archive.unlink()
    if archive_path.exists():
        if not force:
            validate_archive(
                archive_path,
                segment_name=segment_dir.name,
                h5_name=h5_name,
                expected_dump_count=expected_dump_count,
                expected_dump_names=expected_dump_names,
            )
            return
        archive_path.unlink()

    threads = os.environ.get("BNE_ZSTD_THREADS", "1")
    level = os.environ.get("BNE_ZSTD_LEVEL", "3")
    compressor = f"zstd -{level} -T{threads}"
    cmd = ["tar", "-I", compressor, "-cf", str(tmp_archive.name), segment_dir.name]
    subprocess.run(cmd, cwd=segment_dir.parent, check=True)
    validate_archive(
        tmp_archive,
        segment_name=segment_dir.name,
        h5_name=h5_name,
        expected_dump_count=expected_dump_count,
        expected_dump_names=expected_dump_names,
    )
    tmp_archive.replace(archive_path)


def validate_archive(
    archive_path,
    segment_name=None,
    h5_name=None,
    expected_dump_count=None,
    expected_dump_names=None,
):
    result = subprocess.run(
        ["tar", "-I", "zstd", "-tf", str(archive_path)],
        stdout=subprocess.PIPE,
        text=True,
        check=True,
    )
    entries = set(result.stdout.splitlines())
    if segment_name is not None and h5_name is not None:
        h5_entry = f"{segment_name}/{h5_name}"
        if h5_entry not in entries:
            raise ValueError(f"{archive_path} does not contain {h5_entry}")
    if segment_name is not None and expected_dump_count is not None:
        dump_prefix = f"{segment_name}/dump"
        dump_names = sorted(
            Path(entry).name
            for entry in entries
            if entry.startswith(dump_prefix) and entry.endswith(".liggghts_forcebalance")
        )
        dump_count = len(dump_names)
        if dump_count != expected_dump_count:
            raise ValueError(
                f"{archive_path} contains {dump_count} raw dumps, expected {expected_dump_count}"
            )
        if expected_dump_names is not None and dump_names != sorted(expected_dump_names):
            raise ValueError(f"{archive_path} raw dump names do not match current raw dumps")


def delete_raw_dumps(dump_paths):
    for dump_path in dump_paths:
        dump_path.unlink()


def pack_segment(segment_dir, h5_name=DEFAULT_H5_NAME, force=False, dry_run=False, keep_dumps=False):
    segment_dir = Path(segment_dir).resolve()
    if not segment_dir.is_dir():
        raise ValueError(f"Segment directory does not exist: {segment_dir}")

    dump_paths = sorted(segment_dir.glob("dump*.liggghts_forcebalance"), key=natural_dump_key)
    h5_path = segment_dir / h5_name
    archive_path = segment_dir.parent / f"{segment_dir.name}.tar.zst"

    summary = {
        "segment_dir": str(segment_dir),
        "dump_count": len(dump_paths),
        "h5_path": str(h5_path),
        "archive_path": str(archive_path),
        "dry_run": dry_run,
        "keep_dumps": keep_dumps,
    }

    if dry_run:
        print(json.dumps(summary, indent=2, sort_keys=True))
        return summary

    if not dump_paths:
        if h5_path.exists():
            validate_h5(h5_path)
            print(f"skip: no raw dumps remain; validated {h5_path}")
            return summary
        raise ValueError(f"No raw dump files found in {segment_dir}")

    if os.environ.get("BNE_PACK_ALLOW_LATEST_SEGMENT", "0") != "1" and is_latest_segment(segment_dir):
        print(f"skip: latest segment may still be active; left {len(dump_paths)} raw dump files in {segment_dir}")
        return summary

    expected_sources = [path.name for path in dump_paths]
    expected_rows = None
    if h5_path.exists() and not force:
        validate_h5(h5_path, expected_frames=len(dump_paths), expected_sources=expected_sources)
        print(f"reuse: validated existing {h5_path}")
    else:
        expected_rows = write_h5(segment_dir, dump_paths, h5_path, force=force)
        print(f"wrote: {h5_path} ({len(dump_paths)} frames, {expected_rows} particle rows)")

    archive_segment(
        segment_dir,
        archive_path,
        force=force,
        h5_name=h5_path.name,
        expected_dump_count=len(dump_paths),
        expected_dump_names=expected_sources,
    )
    print(f"archived: {archive_path}")

    validate_h5(
        h5_path,
        expected_frames=len(dump_paths),
        expected_rows=expected_rows,
        expected_sources=expected_sources,
    )
    validate_archive(
        archive_path,
        segment_name=segment_dir.name,
        h5_name=h5_path.name,
        expected_dump_count=len(dump_paths),
        expected_dump_names=expected_sources,
    )

    if not keep_dumps:
        delete_raw_dumps(dump_paths)
        print(f"deleted: {len(dump_paths)} raw dump files from {segment_dir}")
    else:
        print(f"kept: {len(dump_paths)} raw dump files in {segment_dir}")

    return summary


def main():
    parser = argparse.ArgumentParser(description="Pack one forcebalance segment into HDF5 and tar.zst.")
    parser.add_argument("segment_dir", help="Path to DEM/post/forcebalance_segments/t*_to_*")
    parser.add_argument("--h5-name", default=DEFAULT_H5_NAME, help="HDF5 filename to keep inside the segment")
    parser.add_argument("--force", action="store_true", help="Rebuild existing HDF5/archive if present")
    parser.add_argument("--dry-run", action="store_true", help="Report intended work without writing or deleting")
    parser.add_argument("--keep-dumps", action="store_true", help="Keep raw dump files after HDF5/archive validation")
    args = parser.parse_args()

    try:
        pack_segment(
            args.segment_dir,
            h5_name=args.h5_name,
            force=args.force,
            dry_run=args.dry_run,
            keep_dumps=args.keep_dumps,
        )
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
