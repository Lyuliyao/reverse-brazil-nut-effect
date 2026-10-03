#!/usr/bin/env python3
import argparse
import csv
import json
import math
import re
from pathlib import Path

try:
    import h5py
except ImportError:
    h5py = None


def parse_case_params(case_path):
    params_path = case_path / "case_params.json"
    if params_path.exists():
        return json.loads(params_path.read_text(encoding="utf-8"))
    return {}


def parse_liggghts_timestep(case_path):
    text = (case_path / "DEM" / "in.liggghts_setting").read_text(encoding="utf-8")
    match = re.search(r"^\s*timestep\s+([0-9.eE+-]+)\s*$", text, flags=re.MULTILINE)
    if not match:
        raise ValueError(f"Could not find timestep in {case_path / 'DEM' / 'in.liggghts_setting'}")
    return float(match.group(1))


def parse_box(case_path):
    values = {}
    pattern = re.compile(r"^\s*variable\s+(\w+)\s+equal\s+([-+0-9.eE]+)\s*$")
    for line in (case_path / "DEM" / "in.liggghts_variable").read_text(encoding="utf-8").splitlines():
        match = pattern.match(line)
        if match:
            values[match.group(1)] = float(match.group(2))
    return values


def read_dump_timestep(path):
    with path.open(encoding="utf-8") as handle:
        line_iter = iter(handle)
        for raw_line in line_iter:
            if raw_line.strip().startswith("ITEM: TIMESTEP"):
                try:
                    return int(next(line_iter).strip())
                except StopIteration as exc:
                    raise ValueError(f"Unexpected EOF after TIMESTEP in {path}") from exc
    raise ValueError(f"Failed to parse timestep from dump: {path}")


def parse_segment_time_label(label):
    match = re.match(r"^t(.+)_to_(.+)$", label)
    if not match:
        return None

    def decode(value):
        return float(value.replace("m", "-").replace("p", "."))

    try:
        return decode(match.group(1)), decode(match.group(2))
    except ValueError:
        return None


def decode_h5_string(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def physical_time_from_segment_timestep(segment_label, timestep, dem_dt):
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


def collect_frame_records(case_path, dem_dt):
    segment_root = case_path / "DEM" / "post" / "forcebalance_segments"
    records_by_key = {}

    def add_record(record):
        key = (record["segment_label"], record["dump_name"])
        existing = records_by_key.get(key)
        if existing is None or (existing["kind"] == "hdf5" and record["kind"] == "raw"):
            records_by_key[key] = record

    if segment_root.exists():
        for segment_dir in sorted(path for path in segment_root.iterdir() if path.is_dir()):
            h5_path = segment_dir / "forcebalance.h5"
            if h5_path.exists():
                if h5py is None:
                    raise ImportError(
                        f"h5py is required to read packed forcebalance data: {h5_path}"
                    )
                with h5py.File(h5_path, "r") as h5_file:
                    frames = h5_file["frames"]
                    for frame_index, frame in enumerate(frames):
                        timestep = int(frame["timestep"])
                        add_record(
                            {
                                "kind": "hdf5",
                                "physical_time": physical_time_from_segment_timestep(segment_dir.name, timestep, dem_dt),
                                "timestep": timestep,
                                "h5_path": h5_path,
                                "frame_index": frame_index,
                                "row_start": int(frame["row_start"]),
                                "row_count": int(frame["row_count"]),
                                "dump_name": decode_h5_string(frame["source_name"]),
                                "segment_label": segment_dir.name,
                            }
                        )

            for dump_path in segment_dir.glob("dump*.liggghts_forcebalance"):
                timestep = read_dump_timestep(dump_path)
                add_record(
                    {
                        "kind": "raw",
                        "physical_time": physical_time_from_segment_timestep(segment_dir.name, timestep, dem_dt),
                        "timestep": timestep,
                        "path": dump_path,
                        "dump_name": dump_path.name,
                        "segment_label": segment_dir.name,
                    }
                )

    records = list(records_by_key.values())

    # Legacy location used before segment-specific dump directories were added.
    if not records:
        for dump_path in (case_path / "DEM" / "post").glob("dump*.liggghts_forcebalance"):
            timestep = read_dump_timestep(dump_path)
            records.append(
                {
                    "kind": "raw",
                    "physical_time": timestep * dem_dt,
                    "timestep": timestep,
                    "path": dump_path,
                    "dump_name": dump_path.name,
                    "segment_label": "legacy",
                }
            )

    records.sort(key=lambda item: (item["physical_time"], item["timestep"], item["segment_label"], item["dump_name"]))
    return records


def read_dump(path):
    timestep = None
    header = None
    rows = []
    in_atoms_section = False

    # Stream each dump instead of loading the full file into memory first.
    with path.open(encoding="utf-8") as handle:
        line_iter = iter(handle)
        for raw_line in line_iter:
            line = raw_line.strip()
            if line.startswith("ITEM: TIMESTEP"):
                try:
                    timestep = int(next(line_iter).strip())
                except StopIteration as exc:
                    raise ValueError(f"Unexpected EOF after TIMESTEP in {path}") from exc
                in_atoms_section = False
                continue
            if line.startswith("ITEM: ATOMS "):
                header = line.split()[2:]
                in_atoms_section = True
                continue
            if line.startswith("ITEM: "):
                in_atoms_section = False
                continue
            if not in_atoms_section or not line:
                continue

            values = line.split()
            entry = {}
            for key, value in zip(header, values):
                if key in {"id", "type"}:
                    entry[key] = int(value)
                else:
                    entry[key] = float(value)
            rows.append(entry)

    if header is None or timestep is None:
        raise ValueError(f"Failed to parse dump: {path}")
    return timestep, rows


def rows_from_h5_array(data, columns):
    rows = []
    for values in data:
        rows.append({column: float(values[index]) for index, column in enumerate(columns)})
    return rows

def mean(values):
    return sum(values) / len(values) if values else 0.0


def percentile(values, frac):
    if not values:
        return 0.0
    sorted_values = sorted(values)
    idx = frac * (len(sorted_values) - 1)
    lo = int(math.floor(idx))
    hi = int(math.ceil(idx))
    if lo == hi:
        return sorted_values[lo]
    weight = idx - lo
    return sorted_values[lo] * (1.0 - weight) + sorted_values[hi] * weight


def particle_volume(radius):
    return 4.0 * math.pi * radius ** 3 / 3.0


def species_rows(rows):
    radii = sorted({row["radius"] for row in rows})
    if len(radii) != 2:
        raise ValueError(f"Expected exactly 2 radii, found {radii}")
    small_radius, large_radius = radii
    small = [row for row in rows if row["radius"] == small_radius]
    large = [row for row in rows if row["radius"] == large_radius]
    return small_radius, large_radius, small, large


def vec_component(row, prefix):
    return (
        row[f"{prefix}[1]"],
        row[f"{prefix}[2]"],
        row[f"{prefix}[3]"],
    )


def z_profile(rows, small_radius, large_radius, box, n_bins):
    area = (box["xhi"] - box["xlo"]) * (box["yhi"] - box["ylo"])
    dz = (box["zhi"] - box["zlo"]) / n_bins
    bins = []
    for ibin in range(n_bins):
        z0 = box["zlo"] + ibin * dz
        z1 = z0 + dz
        bins.append(
            {
                "bin_index": ibin,
                "z_center": z0 + 0.5 * dz,
                "phi_small": 0.0,
                "phi_large": 0.0,
                "phi_total": 0.0,
                "count_small": 0,
                "count_large": 0,
            }
        )
    for row in rows:
        ibin = int((row["z"] - box["zlo"]) / dz)
        ibin = max(0, min(n_bins - 1, ibin))
        volume = particle_volume(row["radius"])
        phi = volume / (area * dz)
        if row["radius"] == small_radius:
            bins[ibin]["phi_small"] += phi
            bins[ibin]["count_small"] += 1
        else:
            bins[ibin]["phi_large"] += phi
            bins[ibin]["count_large"] += 1
        bins[ibin]["phi_total"] += phi
    return bins


def local_phi_near_large(large_rows, all_rows, box):
    area = (box["xhi"] - box["xlo"]) * (box["yhi"] - box["ylo"])
    if not large_rows:
        return {
            "phi_total_below_large": 0.0,
            "phi_small_below_large": 0.0,
            "phi_total_above_large": 0.0,
            "phi_small_above_large": 0.0,
        }
    large_radius = max(row["radius"] for row in large_rows)
    small_rows = [row for row in all_rows if row["radius"] < large_radius]

    phi_total_below = []
    phi_small_below = []
    phi_total_above = []
    phi_small_above = []

    for large in large_rows:
        slab = 2.0 * large["radius"]
        below_lo = max(box["zlo"], large["z"] - slab)
        below_hi = large["z"]
        above_lo = large["z"]
        above_hi = min(box["zhi"], large["z"] + slab)
        below_width = max(1.0e-16, below_hi - below_lo)
        above_width = max(1.0e-16, above_hi - above_lo)

        total_below = sum(
            particle_volume(row["radius"])
            for row in all_rows
            if below_lo <= row["z"] < below_hi
        )
        small_below = sum(
            particle_volume(row["radius"])
            for row in small_rows
            if below_lo <= row["z"] < below_hi
        )
        total_above = sum(
            particle_volume(row["radius"])
            for row in all_rows
            if above_lo <= row["z"] < above_hi
        )
        small_above = sum(
            particle_volume(row["radius"])
            for row in small_rows
            if above_lo <= row["z"] < above_hi
        )

        phi_total_below.append(total_below / (area * below_width))
        phi_small_below.append(small_below / (area * below_width))
        phi_total_above.append(total_above / (area * above_width))
        phi_small_above.append(small_above / (area * above_width))

    return {
        "phi_total_below_large": mean(phi_total_below),
        "phi_small_below_large": mean(phi_small_below),
        "phi_total_above_large": mean(phi_total_above),
        "phi_small_above_large": mean(phi_small_above),
    }


def summarize_species(rows, gravity_vector):
    z_values = [row["z"] for row in rows]
    f_total_z = [row["fz"] for row in rows]
    f_contact_z = [row["f_contactforce[3]"] for row in rows]
    f_wall_z = [row["f_wallforce[3]"] for row in rows]
    f_fluid_z = [row["f_fluidforce[3]"] for row in rows]
    coordination = [row["c_cc"] for row in rows]
    masses = [row["mass"] for row in rows]
    return {
        "count": len(rows),
        "z_mean": mean(z_values),
        "z_median": percentile(z_values, 0.5),
        "z_p10": percentile(z_values, 0.1),
        "z_p90": percentile(z_values, 0.9),
        "z_min": min(z_values) if z_values else 0.0,
        "z_max": max(z_values) if z_values else 0.0,
        "fz_total_mean": mean(f_total_z),
        "fz_contact_mean": mean(f_contact_z),
        "fz_wall_mean": mean(f_wall_z),
        "fz_fluid_mean": mean(f_fluid_z),
        "fz_gravity_mean": mean(masses) * gravity_vector[2],
        "coordination_mean": mean(coordination),
    }


SERIES_FIELDNAMES = [
    "dump_name",
    "segment",
    "dem_timestep",
    "time",
    "count_small",
    "count_large",
    "z_mean_small",
    "z_mean_large",
    "z_median_small",
    "z_median_large",
    "z_gap_large_minus_small",
    "large_above_small_median_fraction",
    "fz_total_small",
    "fz_contact_small",
    "fz_wall_small",
    "fz_fluid_small",
    "fz_gravity_small",
    "fz_total_large",
    "fz_contact_large",
    "fz_wall_large",
    "fz_fluid_large",
    "fz_gravity_large",
    "fz_fluid_to_abs_gravity_ratio_large",
    "coordination_small",
    "coordination_large",
    "phi_total_below_large",
    "phi_small_below_large",
    "phi_total_above_large",
    "phi_small_above_large",
]

INT_SERIES_FIELDS = {"dem_timestep", "count_small", "count_large"}
FLOAT_SERIES_FIELDS = set(SERIES_FIELDNAMES) - {"dump_name", "segment"} - INT_SERIES_FIELDS


def record_key(record):
    return record["segment_label"], record["dump_name"]


def series_key(row):
    return row["segment"], row["dump_name"]


def coerce_series_row(row):
    coerced = {}
    for field in SERIES_FIELDNAMES:
        value = row.get(field, "")
        if field in INT_SERIES_FIELDS:
            coerced[field] = int(float(value)) if value not in {"", None} else 0
        elif field in FLOAT_SERIES_FIELDS:
            coerced[field] = float(value) if value not in {"", None} else 0.0
        else:
            coerced[field] = value
    return coerced


def read_existing_series(series_path):
    if not series_path.exists():
        return {}
    with series_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not set(SERIES_FIELDNAMES).issubset(reader.fieldnames):
            return {}
        rows = {}
        for row in reader:
            coerced = coerce_series_row(row)
            rows[series_key(coerced)] = coerced
        return rows


def sorted_series_rows(rows):
    return sorted(
        rows,
        key=lambda item: (item["time"], item["dem_timestep"], item["segment"], item["dump_name"]),
    )


def profile_path_for_record(profiles_dir, record):
    return (
        profiles_dir
        / record["segment_label"]
        / f"{Path(record['dump_name']).stem}.z_profile.tsv"
    )


def write_profile(profile_path, profile_rows):
    profile_path.parent.mkdir(parents=True, exist_ok=True)
    with profile_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "bin_index",
                "z_center",
                "phi_small",
                "phi_large",
                "phi_total",
                "count_small",
                "count_large",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in profile_rows:
            writer.writerow(row)


def write_time_series(series_path, series_rows):
    with series_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SERIES_FIELDNAMES, delimiter="\t")
        writer.writeheader()
        for row in series_rows:
            writer.writerow({field: row[field] for field in SERIES_FIELDNAMES})


def summarize_frame(record, rows, gravity, box):
    small_radius, large_radius, small_rows, large_rows = species_rows(rows)
    small_summary = summarize_species(small_rows, gravity)
    large_summary = summarize_species(large_rows, gravity)
    local_phi = local_phi_near_large(large_rows, rows, box)

    small_median = small_summary["z_median"]
    large_above_small_median_fraction = mean([1.0 if row["z"] > small_median else 0.0 for row in large_rows])
    z_gap = large_summary["z_mean"] - small_summary["z_mean"]

    fluid_to_gravity_ratio = 0.0
    if abs(large_summary["fz_gravity_mean"]) > 0.0:
        fluid_to_gravity_ratio = large_summary["fz_fluid_mean"] / abs(large_summary["fz_gravity_mean"])

    return (
        {
            "dump_name": record["dump_name"],
            "segment": record["segment_label"],
            "dem_timestep": record["timestep"],
            "time": record["physical_time"],
            "count_small": small_summary["count"],
            "count_large": large_summary["count"],
            "z_mean_small": small_summary["z_mean"],
            "z_mean_large": large_summary["z_mean"],
            "z_median_small": small_summary["z_median"],
            "z_median_large": large_summary["z_median"],
            "z_gap_large_minus_small": z_gap,
            "large_above_small_median_fraction": large_above_small_median_fraction,
            "fz_total_small": small_summary["fz_total_mean"],
            "fz_contact_small": small_summary["fz_contact_mean"],
            "fz_wall_small": small_summary["fz_wall_mean"],
            "fz_fluid_small": small_summary["fz_fluid_mean"],
            "fz_gravity_small": small_summary["fz_gravity_mean"],
            "fz_total_large": large_summary["fz_total_mean"],
            "fz_contact_large": large_summary["fz_contact_mean"],
            "fz_wall_large": large_summary["fz_wall_mean"],
            "fz_fluid_large": large_summary["fz_fluid_mean"],
            "fz_gravity_large": large_summary["fz_gravity_mean"],
            "fz_fluid_to_abs_gravity_ratio_large": fluid_to_gravity_ratio,
            "coordination_small": small_summary["coordination_mean"],
            "coordination_large": large_summary["coordination_mean"],
            **local_phi,
        },
        small_radius,
        large_radius,
    )


def analyze_case(case_path, n_bins, full_rebuild=False):
    case_path = Path(case_path).resolve()
    params = parse_case_params(case_path)
    box = parse_box(case_path)
    dem_dt = parse_liggghts_timestep(case_path)
    gravity = (box["gx"], box["gy"], box["gz"])

    analysis_dir = case_path / "analysis"
    profiles_dir = analysis_dir / "z_profiles"
    series_path = analysis_dir / "time_series.tsv"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    profiles_dir.mkdir(parents=True, exist_ok=True)

    dump_records = collect_frame_records(case_path, dem_dt)

    if not dump_records:
        raise ValueError(f"No forcebalance dump frames found under {case_path / 'DEM' / 'post'}")

    record_keys = {record_key(record) for record in dump_records}
    if full_rebuild:
        series_by_key = {}
    else:
        series_by_key = {
            key: row
            for key, row in read_existing_series(series_path).items()
            if key in record_keys
        }

    records_to_process = [
        record
        for record in dump_records
        if record_key(record) not in series_by_key
        or not profile_path_for_record(profiles_dir, record).exists()
    ]

    open_h5 = None
    open_h5_path = None
    open_h5_columns = None
    open_h5_particles = None

    try:
        for record in records_to_process:
            if record["kind"] == "hdf5":
                if h5py is None:
                    raise ImportError("h5py is required to read packed forcebalance data")
                if open_h5_path != record["h5_path"]:
                    if open_h5 is not None:
                        open_h5.close()
                    open_h5 = h5py.File(record["h5_path"], "r")
                    open_h5_path = record["h5_path"]
                    open_h5_particles = open_h5["particles"]
                    open_h5_columns = [
                        decode_h5_string(value)
                        for value in open_h5_particles.attrs["columns"]
                    ]
                row_start = record["row_start"]
                row_end = row_start + record["row_count"]
                rows = rows_from_h5_array(open_h5_particles[row_start:row_end, :], open_h5_columns)
            else:
                _, rows = read_dump(record["path"])

            series_row, small_radius, large_radius = summarize_frame(record, rows, gravity, box)
            series_by_key[record_key(record)] = series_row
            profile_rows = z_profile(rows, small_radius, large_radius, box, n_bins)
            write_profile(profile_path_for_record(profiles_dir, record), profile_rows)
    finally:
        if open_h5 is not None:
            open_h5.close()

    series_rows = sorted_series_rows(series_by_key.values())
    write_time_series(series_path, series_rows)

    initial = series_rows[0]
    final = series_rows[-1]
    cross_time = next(
        (row["time"] for row in series_rows if row["z_gap_large_minus_small"] > 0.0),
        None,
    )
    summary = {
        "case_name": params.get("case_name", case_path.name),
        "case_path": str(case_path),
        "n_dumps": len(series_rows),
        "initial_time": initial["time"],
        "final_time": final["time"],
        "initial_z_gap_large_minus_small": initial["z_gap_large_minus_small"],
        "final_z_gap_large_minus_small": final["z_gap_large_minus_small"],
        "delta_z_gap_large_minus_small": final["z_gap_large_minus_small"] - initial["z_gap_large_minus_small"],
        "initial_large_above_small_median_fraction": initial["large_above_small_median_fraction"],
        "final_large_above_small_median_fraction": final["large_above_small_median_fraction"],
        "cross_time_large_mean_above_small_mean": cross_time,
        "final_fz_fluid_large": final["fz_fluid_large"],
        "final_fz_gravity_large": final["fz_gravity_large"],
        "final_fz_contact_large": final["fz_contact_large"],
        "final_fz_wall_large": final["fz_wall_large"],
        "final_fz_fluid_to_abs_gravity_ratio_large": final["fz_fluid_to_abs_gravity_ratio_large"],
        "final_phi_total_below_large": final["phi_total_below_large"],
        "final_phi_total_above_large": final["phi_total_above_large"],
        "final_phi_small_below_large": final["phi_small_below_large"],
        "final_phi_small_above_large": final["phi_small_above_large"],
        "final_coordination_large": final["coordination_large"],
        "final_coordination_small": final["coordination_small"],
    }
    summary.update({k: v for k, v in params.items() if k not in summary})

    summary_json_path = analysis_dir / "summary.json"
    summary_txt_path = analysis_dir / "summary.txt"
    summary_json_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    summary_txt_path.write_text(
        "\n".join(
            [
                f"case_name: {summary['case_name']}",
                f"case_path: {summary['case_path']}",
                f"n_dumps: {summary['n_dumps']}",
                f"initial_z_gap_large_minus_small: {summary['initial_z_gap_large_minus_small']:.12e}",
                f"final_z_gap_large_minus_small: {summary['final_z_gap_large_minus_small']:.12e}",
                f"delta_z_gap_large_minus_small: {summary['delta_z_gap_large_minus_small']:.12e}",
                f"initial_large_above_small_median_fraction: {summary['initial_large_above_small_median_fraction']:.12e}",
                f"final_large_above_small_median_fraction: {summary['final_large_above_small_median_fraction']:.12e}",
                f"cross_time_large_mean_above_small_mean: {summary['cross_time_large_mean_above_small_mean']}",
                f"final_fz_fluid_large: {summary['final_fz_fluid_large']:.12e}",
                f"final_fz_gravity_large: {summary['final_fz_gravity_large']:.12e}",
                f"final_fz_contact_large: {summary['final_fz_contact_large']:.12e}",
                f"final_fz_wall_large: {summary['final_fz_wall_large']:.12e}",
                f"final_fz_fluid_to_abs_gravity_ratio_large: {summary['final_fz_fluid_to_abs_gravity_ratio_large']:.12e}",
                f"final_phi_total_below_large: {summary['final_phi_total_below_large']:.12e}",
                f"final_phi_total_above_large: {summary['final_phi_total_above_large']:.12e}",
                f"final_phi_small_below_large: {summary['final_phi_small_below_large']:.12e}",
                f"final_phi_small_above_large: {summary['final_phi_small_above_large']:.12e}",
                f"final_coordination_large: {summary['final_coordination_large']:.12e}",
                f"final_coordination_small: {summary['final_coordination_small']:.12e}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"Analysis written under {analysis_dir}")
    print(f"Incremental frames analyzed: {len(records_to_process)}")
    print(f"Existing frames reused:      {len(series_rows) - len(records_to_process)}")
    print(f"Time series: {series_path}")
    print(f"Summary:     {summary_txt_path}")


def main():
    parser = argparse.ArgumentParser(description="Analyze one Brazil-nut CFDEM case.")
    parser.add_argument("case_path", help="Absolute path to the case")
    parser.add_argument("--n-bins", type=int, default=32, help="Number of z bins for volume-fraction profiles")
    parser.add_argument("--full-rebuild", action="store_true", help="Ignore existing analysis outputs and recompute all frames")
    args = parser.parse_args()
    analyze_case(Path(args.case_path), args.n_bins, full_rebuild=args.full_rebuild)


if __name__ == "__main__":
    main()
