#!/usr/bin/env python3
"""Shared I/O helpers for the CFD-DEM mechanism-isolation campaign analysis.

This module centralises every assumption about on-disk data formats so that the
individual analysis scripts (compute_mean_height_metrics.py, compute_vertical_flux.py,
etc.) only deal with numpy arrays keyed by column name.

Authoritative data-format notes (see campaign CLAUDE.md):

* metadata.json (per case dir) carries control_type, small_radius, large_radius,
  small_diameter_um, large_diameter_um, volfrac, angle_deg, rhop, rhof(=1000),
  g_mag, gx, gy, gz, g_eff, g_par, g_perp, dem_timestep(=1e-7),
  particle_dump_every, contact_dump_every, domain{xlo,xhi,ylo,yhi,zlo,zhi},
  seeds, seed_index.  Some older cases use case_params.json instead.

* PARTICLE DATA lives under DEM/post/forcebalance_segments/<label>/ as either
  raw `dump*.liggghts_forcebalance` text files OR a packed `forcebalance.h5`.
  Segment labels encode a time window: `t0_to_1`, `t12_to_13`; encoding maps
  'm'->'-' and 'p'->'.'.  Columns are read BY NAME from the `ITEM: ATOMS` header
  (the EXTENDED 32-column layout has vx vy vz; older runs have 29 columns
  without velocities -- this reader handles both).

* CONTACT DATA (this campaign only) lives under DEM/post/contact_segments/<label>/
  as `contact*.liggghts_local` (LIGGGHTS `dump local`).  The 13 columns, in
  order, are: id1 id2 periodicFlag fx fy fz fnx fny fnz delta cpx cpy cpz.
  fx,fy,fz is the pairwise contact force ON particle id1 FROM id2 (Newton 3rd:
  on id2 it is the negative).  A contact frame at a given timestep coincides
  with a particle dump at the SAME timestep (used for id->type mapping).

The physical-time convention (parse_segment_time_label /
physical_time_from_segment_timestep) is copied verbatim from
brazil_nut/tools/analyze_bne_case.py so that frame times are consistent across
the whole project.

Only numpy + h5py + matplotlib + the standard library are used.
"""

import json
import re
from pathlib import Path

import numpy as np

try:
    import h5py
except ImportError:  # h5py is optional for raw-only cases.
    h5py = None


# Canonical column order of the contact `dump local` files for this campaign.
CONTACT_COLUMNS = [
    "id1", "id2", "periodicFlag",
    "fx", "fy", "fz",
    "fnx", "fny", "fnz",
    "delta", "cpx", "cpy", "cpz",
]


# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------
def read_metadata(case_dir):
    """Return the metadata dict for a case.

    Prefers metadata.json; falls back to case_params.json (older convention).
    Returns {} if neither file exists (callers should be defensive).
    """
    case_dir = Path(case_dir)
    for name in ("metadata.json", "case_params.json"):
        path = case_dir / name
        if path.exists():
            try:
                return json.loads(path.read_text(encoding="utf-8"))
            except (ValueError, OSError):
                continue
    return {}


def case_box(metadata):
    """Return a box dict {xlo,xhi,ylo,yhi,zlo,zhi} from metadata['domain'].

    Falls back to the campaign default domain (0..0.01, 0..0.01, 0..0.005) if
    missing, and warns via the returned 'box_is_default' flag is NOT added --
    callers can detect missing keys themselves.  Kept simple and defensive.
    """
    domain = {}
    if isinstance(metadata, dict):
        domain = metadata.get("domain") or {}
    defaults = {"xlo": 0.0, "xhi": 0.01, "ylo": 0.0, "yhi": 0.01, "zlo": 0.0, "zhi": 0.005}
    box = {}
    for key, default in defaults.items():
        try:
            box[key] = float(domain.get(key, default))
        except (TypeError, ValueError):
            box[key] = default
    return box


# ---------------------------------------------------------------------------
# Segment-time convention (copied verbatim from analyze_bne_case.py)
# ---------------------------------------------------------------------------
def parse_segment_time_label(label):
    """Parse a segment label like 't0_to_1' or 't12_to_13'.

    Returns (start, end) floats, or None if the label is not parseable
    (e.g. 't0_to_initial').  Encoding: 'm'->'-', 'p'->'.'.
    """
    match = re.match(r"^t(.+)_to_(.+)$", label)
    if not match:
        return None

    def decode(value):
        return float(value.replace("m", "-").replace("p", "."))

    try:
        return decode(match.group(1)), decode(match.group(2))
    except ValueError:
        return None


def physical_time_from_segment_timestep(segment_label, timestep, dem_dt):
    """Map a DEM timestep index to physical time using the segment label.

    Convention copied verbatim from analyze_bne_case.py:
      raw_time = timestep * dem_dt
      if start <= raw_time <= end (with 1e-9 tolerance): use raw_time
      else: start + raw_time
    """
    parsed = parse_segment_time_label(segment_label)
    segment_start = parsed[0] if parsed is not None else 0.0
    segment_end = parsed[1] if parsed is not None else None
    raw_time = timestep * dem_dt
    if (
        segment_end is not None
        and segment_start - 1.0e-9 <= raw_time <= segment_end + 1.0e-9
    ):
        return raw_time
    return segment_start + raw_time


def _decode_attr(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


# ---------------------------------------------------------------------------
# Raw dump parsing
# ---------------------------------------------------------------------------
def _read_raw_dump(path):
    """Parse a raw `dump*.liggghts_forcebalance` text file.

    Returns (timestep:int, box_bounds:np.ndarray(3,2), columns:list[str],
    data:np.ndarray(n_atoms, n_cols) float64).  Robust to a variable number of
    columns -- the caller maps by name.
    """
    path = Path(path)
    timestep = None
    n_atoms = None
    box_bounds = []
    columns = None
    data = None

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
                data = np.empty((n_atoms, len(columns)), dtype=np.float64)
                for index in range(n_atoms):
                    try:
                        atom_line = next(line_iter)
                    except StopIteration as exc:
                        raise ValueError(f"Unexpected EOF in ATOMS data in {path}") from exc
                    values = np.fromstring(atom_line, sep=" ", dtype=np.float64)
                    if values.size != len(columns):
                        raise ValueError(
                            f"Expected {len(columns)} ATOMS values, got {values.size} "
                            f"in {path} line {index + 1}"
                        )
                    data[index, :] = values
                break

    if timestep is None or columns is None or data is None:
        raise ValueError(f"Failed to parse required dump sections in {path}")
    if len(box_bounds) != 3:
        # Box bounds are not always strictly needed by callers; supply NaNs.
        box_bounds = [[np.nan, np.nan]] * 3
    return timestep, np.asarray(box_bounds, dtype=np.float64), columns, data


def _columns_to_dict(columns, data):
    """Map a 2-D float array to a dict {column_name: 1-D array}."""
    out = {}
    for index, name in enumerate(columns):
        out[name] = data[:, index]
    return out


def _box_dict_from_bounds(box_bounds, metadata):
    """Build a box dict from parsed BOX BOUNDS, falling back to metadata domain."""
    box = case_box(metadata)
    if box_bounds is not None and np.all(np.isfinite(box_bounds)):
        box = {
            "xlo": float(box_bounds[0, 0]), "xhi": float(box_bounds[0, 1]),
            "ylo": float(box_bounds[1, 0]), "yhi": float(box_bounds[1, 1]),
            "zlo": float(box_bounds[2, 0]), "zhi": float(box_bounds[2, 1]),
        }
    return box


# ---------------------------------------------------------------------------
# Frame collection (h5 + raw, deduped, time-sorted)
# ---------------------------------------------------------------------------
def _collect_particle_records(case_dir, dem_dt):
    """Collect particle-frame records across all segments.

    Returns a list of record dicts sorted by (physical_time, timestep, segment,
    name).  Each record carries enough info to lazily load its rows.  Raw dumps
    take priority over packed h5 for an identical (segment, name) key so that
    re-packing in progress does not duplicate frames.
    """
    segment_root = Path(case_dir) / "DEM" / "post" / "forcebalance_segments"
    records_by_key = {}

    def add(record):
        key = (record["segment_label"], record["dump_name"])
        existing = records_by_key.get(key)
        if existing is None or (existing["kind"] == "hdf5" and record["kind"] == "raw"):
            records_by_key[key] = record

    if segment_root.exists():
        for segment_dir in sorted(p for p in segment_root.iterdir() if p.is_dir()):
            label = segment_dir.name
            h5_path = segment_dir / "forcebalance.h5"
            if h5_path.exists():
                if h5py is None:
                    raise ImportError(
                        f"h5py is required to read packed forcebalance data: {h5_path}"
                    )
                with h5py.File(h5_path, "r") as h5_file:
                    if "frames" not in h5_file:
                        continue
                    frames = h5_file["frames"]
                    for frame_index, frame in enumerate(frames):
                        timestep = int(frame["timestep"])
                        add({
                            "kind": "hdf5",
                            "physical_time": physical_time_from_segment_timestep(label, timestep, dem_dt),
                            "timestep": timestep,
                            "h5_path": h5_path,
                            "frame_index": frame_index,
                            "row_start": int(frame["row_start"]),
                            "row_count": int(frame["row_count"]),
                            "dump_name": _decode_attr(frame["source_name"]),
                            "segment_label": label,
                        })
            for dump_path in segment_dir.glob("dump*.liggghts_forcebalance"):
                ts, _, _, _ = _read_raw_dump(dump_path)
                add({
                    "kind": "raw",
                    "physical_time": physical_time_from_segment_timestep(label, ts, dem_dt),
                    "timestep": ts,
                    "path": dump_path,
                    "dump_name": dump_path.name,
                    "segment_label": label,
                })

    records = list(records_by_key.values())
    # Legacy fallback: dumps directly under DEM/post.
    if not records:
        legacy_dir = Path(case_dir) / "DEM" / "post"
        for dump_path in legacy_dir.glob("dump*.liggghts_forcebalance"):
            ts, _, _, _ = _read_raw_dump(dump_path)
            records.append({
                "kind": "raw",
                "physical_time": ts * dem_dt,
                "timestep": ts,
                "path": dump_path,
                "dump_name": dump_path.name,
                "segment_label": "legacy",
            })

    records.sort(
        key=lambda r: (r["physical_time"], r["timestep"], r["segment_label"], r["dump_name"])
    )
    return records


def iter_frames(case_dir):
    """Yield (time, columns_dict, box) for every particle frame in a case.

    * time: float physical time (s).
    * columns_dict: {column_name: 1-D np.ndarray}, including whatever columns are
      present in the dump header.  vx/vy/vz may be absent for older runs --
      callers must check with has_velocity().
    * box: dict {xlo,xhi,ylo,yhi,zlo,zhi} (from BOX BOUNDS when available, else
      metadata domain).

    Frames across all segments are merged, deduped at segment boundaries
    (identical (segment, name) keys), and yielded in ascending time order.
    Additionally, frames with an identical timestep at segment boundaries are
    de-duplicated here so each physical timestep is emitted once.
    """
    case_dir = Path(case_dir)
    metadata = read_metadata(case_dir)
    dem_dt = float(metadata.get("dem_timestep", 1.0e-7))
    records = _collect_particle_records(case_dir, dem_dt)

    seen_timesteps = set()
    open_h5 = None
    open_h5_path = None
    open_columns = None
    open_particles = None
    open_box_bounds = None
    try:
        for record in records:
            ts = record["timestep"]
            if ts in seen_timesteps:
                continue
            seen_timesteps.add(ts)

            if record["kind"] == "hdf5":
                if h5py is None:
                    raise ImportError("h5py is required to read packed forcebalance data")
                if open_h5_path != record["h5_path"]:
                    if open_h5 is not None:
                        open_h5.close()
                    open_h5 = h5py.File(record["h5_path"], "r")
                    open_h5_path = record["h5_path"]
                    open_particles = open_h5["particles"]
                    open_columns = [_decode_attr(v) for v in open_particles.attrs["columns"]]
                    open_box_bounds = open_h5["box_bounds"] if "box_bounds" in open_h5 else None
                rs = record["row_start"]
                re_ = rs + record["row_count"]
                data = np.asarray(open_particles[rs:re_, :], dtype=np.float64)
                columns = open_columns
                bounds = None
                if open_box_bounds is not None:
                    try:
                        bounds = np.asarray(open_box_bounds[record["frame_index"]], dtype=np.float64)
                    except (IndexError, ValueError):
                        bounds = None
            else:
                _, bounds, columns, data = _read_raw_dump(record["path"])

            box = _box_dict_from_bounds(bounds, metadata)
            yield record["physical_time"], _columns_to_dict(columns, data), box
    finally:
        if open_h5 is not None:
            open_h5.close()


def has_velocity(columns_dict):
    """True if the frame carries vx/vy/vz columns."""
    return all(key in columns_dict for key in ("vx", "vy", "vz"))


# ---------------------------------------------------------------------------
# Contact frames
# ---------------------------------------------------------------------------
def _collect_contact_records(case_dir, dem_dt):
    """Collect contact-frame records across all contact segments.

    Each `contact*.liggghts_local` file may itself contain multiple snapshots;
    we record (path, label) and parse on demand.  Returns a list sorted by
    (segment_label, name).
    """
    segment_root = Path(case_dir) / "DEM" / "post" / "contact_segments"
    records = []
    if segment_root.exists():
        for segment_dir in sorted(p for p in segment_root.iterdir() if p.is_dir()):
            for path in sorted(segment_dir.glob("contact*.liggghts_local")):
                records.append({"path": path, "segment_label": segment_dir.name})
    return records


def _iter_local_snapshots(path):
    """Yield (timestep:int, data:np.ndarray(n_entries, 13)) per snapshot.

    Parses a LIGGGHTS `dump local` file.  Header per snapshot:
      ITEM: TIMESTEP / <step>
      ITEM: NUMBER OF ENTRIES / <n>
      ITEM: BOX BOUNDS / 3 lines
      ITEM: ENTRIES c_cpgl[1] ... c_cpgl[13] / <n rows>
    A file may contain more than one snapshot; all are yielded.
    """
    # NOTE: robust to a truncated TRAILING frame (the file may be appended to by a
    # live LIGGGHTS run while we read). Any EOF / short row mid-frame cleanly STOPS
    # iteration, silently dropping only the incomplete final frame. All COMPLETE
    # frames already yielded are returned.
    path = Path(path)
    with path.open("r", encoding="utf-8") as handle:
        line_iter = iter(handle)
        timestep = None
        n_entries = None
        for raw_line in line_iter:
            line = raw_line.strip()
            if line == "ITEM: TIMESTEP":
                try:
                    timestep = int(next(line_iter).strip())
                except StopIteration:
                    return
                continue
            if line == "ITEM: NUMBER OF ENTRIES":
                try:
                    n_entries = int(next(line_iter).strip())
                except StopIteration:
                    return
                continue
            if line.startswith("ITEM: BOX BOUNDS"):
                truncated = False
                for _ in range(3):
                    try:
                        next(line_iter)
                    except StopIteration:
                        truncated = True
                        break
                if truncated:
                    return
                continue
            if line.startswith("ITEM: ENTRIES"):
                if n_entries is None:
                    raise ValueError(f"ENTRIES section appeared before count in {path}")
                ncol = len(CONTACT_COLUMNS)
                data = np.empty((n_entries, ncol), dtype=np.float64)
                truncated = False
                for index in range(n_entries):
                    try:
                        entry_line = next(line_iter)
                    except StopIteration:
                        truncated = True
                        break
                    values = np.fromstring(entry_line, sep=" ", dtype=np.float64)
                    if values.size < ncol:
                        # incomplete trailing row (live append) -> drop this frame
                        truncated = True
                        break
                    data[index, :] = values[:ncol]
                if truncated or timestep is None:
                    return
                yield timestep, data
                timestep = None
                n_entries = None


def read_contact_frames(case_dir):
    """Yield (time, contact_array, timestep) for every contact snapshot.

    * time: physical time (s), via physical_time_from_segment_timestep using the
      segment label and dem_timestep from metadata.
    * contact_array: np.ndarray(n_entries, 13) in CONTACT_COLUMNS order.
    * timestep: int DEM timestep (so callers can match the coincident particle
      dump at the SAME timestep for id->type mapping).

    Snapshots are deduped by timestep and yielded in ascending time order.
    """
    case_dir = Path(case_dir)
    metadata = read_metadata(case_dir)
    dem_dt = float(metadata.get("dem_timestep", 1.0e-7))
    records = _collect_contact_records(case_dir, dem_dt)

    collected = []  # (time, timestep, path, snapshot_index)
    for record in records:
        label = record["segment_label"]
        for snap_index, (timestep, data) in enumerate(_iter_local_snapshots(record["path"])):
            time = physical_time_from_segment_timestep(label, timestep, dem_dt)
            collected.append((time, timestep, record["path"], snap_index, data))

    collected.sort(key=lambda item: (item[0], item[1], str(item[2]), item[3]))
    seen = set()
    for time, timestep, _path, _snap, data in collected:
        if timestep in seen:
            continue
        seen.add(timestep)
        yield time, data, timestep


def contact_columns_index():
    """Return {column_name: index} for the 13 contact columns."""
    return {name: idx for idx, name in enumerate(CONTACT_COLUMNS)}


# ---------------------------------------------------------------------------
# Species masks
# ---------------------------------------------------------------------------
def species_masks(columns_dict):
    """Return (mask_small, mask_large, small_radius, large_radius) for a frame.

    Prefers the `type` column (type 1 = small, type 2 = large).  Falls back to
    the `radius` column (small = min radius, large = max radius) when type is
    absent or degenerate.  Defensive: if only one species is present, the other
    mask is all-False.
    """
    n = None
    for key in ("id", "type", "radius", "x", "z"):
        if key in columns_dict:
            n = columns_dict[key].shape[0]
            break
    if n is None:
        raise ValueError("Cannot determine frame size: no recognised columns present")

    small_radius = None
    large_radius = None

    radius = columns_dict.get("radius")
    if radius is not None:
        unique = np.unique(radius)
        if unique.size >= 2:
            small_radius = float(unique.min())
            large_radius = float(unique.max())

    types = columns_dict.get("type")
    if types is not None and np.unique(types).size >= 2:
        mask_small = types == 1
        mask_large = types == 2
        # Reconcile radii with type assignment when both are available.
        if radius is not None:
            if np.any(mask_small):
                small_radius = float(np.mean(radius[mask_small]))
            if np.any(mask_large):
                large_radius = float(np.mean(radius[mask_large]))
        return mask_small, mask_large, small_radius, large_radius

    if radius is not None and small_radius is not None:
        mask_small = radius == small_radius
        mask_large = radius == large_radius
        return mask_small, mask_large, small_radius, large_radius

    # Only one species detectable.
    mask_small = np.zeros(n, dtype=bool)
    mask_large = np.zeros(n, dtype=bool)
    if radius is not None:
        small_radius = large_radius = float(np.unique(radius)[0])
        mask_small = np.ones(n, dtype=bool)
    return mask_small, mask_large, small_radius, large_radius


def particle_volume(radius):
    """Sphere volume V = (4/3) pi r^3 (works on scalars or arrays)."""
    return (4.0 / 3.0) * np.pi * np.asarray(radius, dtype=np.float64) ** 3


def normal_gravity(metadata):
    """Return g_perp for the case, falling back to g_eff/g_mag if g_perp==0.

    Returns (g_value, used_fallback_bool).  Used by contact-force normalisation.
    """
    g_perp = float(metadata.get("g_perp", 0.0) or 0.0)
    if g_perp > 0.0:
        return g_perp, False
    for key in ("g_eff", "g_mag"):
        val = metadata.get(key)
        if val is not None and float(val) > 0.0:
            return float(val), True
    return 9.81, True


if __name__ == "__main__":  # pragma: no cover - smoke usage hint
    import argparse

    parser = argparse.ArgumentParser(
        description="Quick inventory of frames/contacts for a case (read-only)."
    )
    parser.add_argument("--case", required=True, help="Path to a case directory")
    args = parser.parse_args()

    meta = read_metadata(args.case)
    print(f"case: {args.case}")
    print(f"control_type: {meta.get('control_type')}")
    n_frames = 0
    has_vel = None
    for time, cols, box in iter_frames(args.case):
        n_frames += 1
        if has_vel is None:
            has_vel = has_velocity(cols)
    print(f"particle frames: {n_frames}; has velocity columns: {has_vel}")
    n_contact = sum(1 for _ in read_contact_frames(args.case))
    print(f"contact snapshots: {n_contact}")
