#!/usr/bin/env python3
import argparse
import csv
import glob
import os
import re
import sys


def parse_variables(path):
    values = {}
    pattern = re.compile(r"^\s*variable\s+(\w+)\s+equal\s+([-+0-9.eE]+)\s*$")
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            match = pattern.match(line)
            if match:
                values[match.group(1)] = float(match.group(2))
    missing = [name for name in ("gx", "gy", "gz", "zlo", "zhi") if name not in values]
    if missing:
        raise ValueError(f"Missing variables in {path}: {', '.join(missing)}")
    return values


def vec_add(a, b):
    return [a[i] + b[i] for i in range(3)]


def vec_sub(a, b):
    return [a[i] - b[i] for i in range(3)]


def vec_inf_norm(a):
    return max(abs(a[0]), abs(a[1]), abs(a[2]))


REQUIRED_COLUMNS = [
    "id", "type", "mass", "x", "y", "z", "fx", "fy", "fz",
    "f_contactforce[1]", "f_contactforce[2]", "f_contactforce[3]",
    "f_gradPForce[1]", "f_gradPForce[2]", "f_gradPForce[3]",
    "f_viscForce[1]", "f_viscForce[2]", "f_viscForce[3]",
    "f_fluidforce[1]", "f_fluidforce[2]", "f_fluidforce[3]",
    "f_dragforce[1]", "f_dragforce[2]", "f_dragforce[3]",
    "radius",
]

OPTIONAL_WALLFORCE_COLUMNS = [
    "f_wallforce[1]", "f_wallforce[2]", "f_wallforce[3]",
]


def read_dump(path):
    with open(path, "r", encoding="utf-8") as handle:
        lines = [line.rstrip("\n") for line in handle]

    timestep = None
    header = None
    rows = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("ITEM: TIMESTEP") and i + 1 < len(lines):
            timestep = lines[i + 1].strip()
            i += 2
            continue
        if line.startswith("ITEM: ATOMS "):
            header = line.split()[2:]
            i += 1
            while i < len(lines) and not lines[i].startswith("ITEM:"):
                if lines[i].strip():
                    rows.append(lines[i].split())
                i += 1
            continue
        i += 1

    if header is None:
        raise ValueError(f"Could not find ATOMS section in {path}")

    index = {name: idx for idx, name in enumerate(header)}
    missing = [name for name in REQUIRED_COLUMNS if name not in index]
    if missing:
        raise ValueError(f"Missing required dump columns in {path}: {', '.join(missing)}")

    has_explicit_wallforce = all(name in index for name in OPTIONAL_WALLFORCE_COLUMNS)

    parsed = []
    for row in rows:
        entry = {}
        for key, idx in index.items():
            value = row[idx]
            if key in {"id", "type"}:
                entry[key] = value
            else:
                entry[key] = float(value)
        parsed.append(entry)

    return timestep, parsed, has_explicit_wallforce


def analyze_dump(dump_path, values, output_dir, balance_tol, wall_force_tol, wall_gap_tol):
    timestep, rows, has_explicit_wallforce = read_dump(dump_path)
    basename = os.path.basename(dump_path)
    stem = basename.rsplit(".", 1)[0]
    tsv_path = os.path.join(output_dir, f"{stem}.force_balance.tsv")
    report_path = os.path.join(output_dir, f"{stem}.force_balance_report.txt")

    gx, gy, gz = values["gx"], values["gy"], values["gz"]
    zlo, zhi = values["zlo"], values["zhi"]

    max_primary_balance = 0.0
    max_full_balance = 0.0
    max_wall_force = 0.0
    max_wall_balance_residual = 0.0
    max_wall_balance_error = 0.0
    max_fluid_residual = 0.0
    wall_force_off_wall = []
    significant_wall_force = 0

    fieldnames = [
        "id", "type", "mass", "x", "y", "z", "radius",
        "f_total_x", "f_total_y", "f_total_z",
        "f_contact_pp_x", "f_contact_pp_y", "f_contact_pp_z",
        "f_wall_x", "f_wall_y", "f_wall_z",
        "f_wall_balance_residual_x", "f_wall_balance_residual_y", "f_wall_balance_residual_z",
        "f_wall_balance_error_x", "f_wall_balance_error_y", "f_wall_balance_error_z",
        "f_gravity_x", "f_gravity_y", "f_gravity_z",
        "f_fluid_total_x", "f_fluid_total_y", "f_fluid_total_z",
        "f_drag_x", "f_drag_y", "f_drag_z",
        "f_gradP_x", "f_gradP_y", "f_gradP_z",
        "f_visc_x", "f_visc_y", "f_visc_z",
        "f_fluid_residual_x", "f_fluid_residual_y", "f_fluid_residual_z",
        "primary_balance_x", "primary_balance_y", "primary_balance_z",
        "full_balance_x", "full_balance_y", "full_balance_z",
        "touch_lower_wall", "touch_upper_wall",
    ]

    with open(tsv_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for row in rows:
            total = [row["fx"], row["fy"], row["fz"]]
            contact_pp = [row["f_contactforce[1]"], row["f_contactforce[2]"], row["f_contactforce[3]"]]
            if has_explicit_wallforce:
                wall_force = [row["f_wallforce[1]"], row["f_wallforce[2]"], row["f_wallforce[3]"]]
            else:
                wall_force = [0.0, 0.0, 0.0]
            grad_p = [row["f_gradPForce[1]"], row["f_gradPForce[2]"], row["f_gradPForce[3]"]]
            visc = [row["f_viscForce[1]"], row["f_viscForce[2]"], row["f_viscForce[3]"]]
            fluid_total = [row["f_fluidforce[1]"], row["f_fluidforce[2]"], row["f_fluidforce[3]"]]
            drag = [row["f_dragforce[1]"], row["f_dragforce[2]"], row["f_dragforce[3]"]]
            gravity = [row["mass"] * gx, row["mass"] * gy, row["mass"] * gz]
            fluid_residual = vec_sub(vec_sub(vec_sub(fluid_total, drag), grad_p), visc)
            wall_balance_residual = vec_sub(vec_sub(vec_sub(total, contact_pp), gravity), fluid_total)
            if not has_explicit_wallforce:
                wall_force = wall_balance_residual
            wall_balance_error = vec_sub(wall_force, wall_balance_residual)
            primary_balance = vec_sub(total, vec_add(vec_add(vec_add(contact_pp, wall_force), gravity), fluid_total))
            full_balance = vec_sub(
                total,
                vec_add(
                    vec_add(vec_add(contact_pp, wall_force), gravity),
                    vec_add(vec_add(drag, grad_p), vec_add(visc, fluid_residual)),
                ),
            )

            touch_lower_wall = (row["z"] - row["radius"]) <= (zlo + wall_gap_tol)
            touch_upper_wall = (row["z"] + row["radius"]) >= (zhi - wall_gap_tol)
            wall_mag = vec_inf_norm(wall_force)
            wall_balance_residual_mag = vec_inf_norm(wall_balance_residual)
            wall_balance_error_mag = vec_inf_norm(wall_balance_error)
            fluid_res_mag = vec_inf_norm(fluid_residual)
            primary_balance_mag = vec_inf_norm(primary_balance)
            full_balance_mag = vec_inf_norm(full_balance)

            max_wall_force = max(max_wall_force, wall_mag)
            max_wall_balance_residual = max(max_wall_balance_residual, wall_balance_residual_mag)
            max_wall_balance_error = max(max_wall_balance_error, wall_balance_error_mag)
            max_fluid_residual = max(max_fluid_residual, fluid_res_mag)
            max_primary_balance = max(max_primary_balance, primary_balance_mag)
            max_full_balance = max(max_full_balance, full_balance_mag)

            if wall_mag > wall_force_tol:
                significant_wall_force += 1
                if not (touch_lower_wall or touch_upper_wall):
                    wall_force_off_wall.append(
                        {
                            "id": row["id"],
                            "x": row["x"],
                            "y": row["y"],
                            "z": row["z"],
                            "radius": row["radius"],
                            "wall_force": wall_force,
                        }
                    )

            writer.writerow(
                {
                    "id": row["id"],
                    "type": row["type"],
                    "mass": row["mass"],
                    "x": row["x"],
                    "y": row["y"],
                    "z": row["z"],
                    "radius": row["radius"],
                    "f_total_x": total[0],
                    "f_total_y": total[1],
                    "f_total_z": total[2],
                    "f_contact_pp_x": contact_pp[0],
                    "f_contact_pp_y": contact_pp[1],
                    "f_contact_pp_z": contact_pp[2],
                    "f_wall_x": wall_force[0],
                    "f_wall_y": wall_force[1],
                    "f_wall_z": wall_force[2],
                    "f_wall_balance_residual_x": wall_balance_residual[0],
                    "f_wall_balance_residual_y": wall_balance_residual[1],
                    "f_wall_balance_residual_z": wall_balance_residual[2],
                    "f_wall_balance_error_x": wall_balance_error[0],
                    "f_wall_balance_error_y": wall_balance_error[1],
                    "f_wall_balance_error_z": wall_balance_error[2],
                    "f_gravity_x": gravity[0],
                    "f_gravity_y": gravity[1],
                    "f_gravity_z": gravity[2],
                    "f_fluid_total_x": fluid_total[0],
                    "f_fluid_total_y": fluid_total[1],
                    "f_fluid_total_z": fluid_total[2],
                    "f_drag_x": drag[0],
                    "f_drag_y": drag[1],
                    "f_drag_z": drag[2],
                    "f_gradP_x": grad_p[0],
                    "f_gradP_y": grad_p[1],
                    "f_gradP_z": grad_p[2],
                    "f_visc_x": visc[0],
                    "f_visc_y": visc[1],
                    "f_visc_z": visc[2],
                    "f_fluid_residual_x": fluid_residual[0],
                    "f_fluid_residual_y": fluid_residual[1],
                    "f_fluid_residual_z": fluid_residual[2],
                    "primary_balance_x": primary_balance[0],
                    "primary_balance_y": primary_balance[1],
                    "primary_balance_z": primary_balance[2],
                    "full_balance_x": full_balance[0],
                    "full_balance_y": full_balance[1],
                    "full_balance_z": full_balance[2],
                    "touch_lower_wall": int(touch_lower_wall),
                    "touch_upper_wall": int(touch_upper_wall),
                }
            )

    with open(report_path, "w", encoding="utf-8") as handle:
        handle.write(f"case_dump: {dump_path}\n")
        handle.write(f"timestep: {timestep}\n")
        handle.write(f"num_particles: {len(rows)}\n")
        handle.write(f"gravity_vector: ({gx}, {gy}, {gz})\n")
        handle.write(f"explicit_wallforce_columns: {'yes' if has_explicit_wallforce else 'no'}\n")
        if has_explicit_wallforce:
            handle.write("assumption: total force is decomposed as particle-particle contact + explicit wallforce + gravity + fluid total.\n")
            handle.write("check: explicit wallforce is also compared against the balance-derived wall residual.\n")
        else:
            handle.write("assumption: for this case, total force is decomposed as particle-particle contact + wall residual + gravity + fluid total.\n")
            handle.write("assumption: wall residual is interpreted as wall force because the active DEM-side force sources are pair/gran, wall/gran, gravity, and couple/cfd/force.\n")
        handle.write(f"max_primary_balance_inf_norm: {max_primary_balance:.12e}\n")
        handle.write(f"max_full_balance_inf_norm: {max_full_balance:.12e}\n")
        handle.write(f"max_wall_force_inf_norm: {max_wall_force:.12e}\n")
        handle.write(f"max_wall_balance_residual_inf_norm: {max_wall_balance_residual:.12e}\n")
        handle.write(f"max_wall_balance_error_inf_norm: {max_wall_balance_error:.12e}\n")
        handle.write(f"max_fluid_residual_inf_norm: {max_fluid_residual:.12e}\n")
        handle.write(f"significant_wall_force_particles: {significant_wall_force}\n")
        handle.write(f"wall_force_off_wall_particles: {len(wall_force_off_wall)}\n")
        if wall_force_off_wall:
            handle.write("off_wall_examples:\n")
            for item in wall_force_off_wall[:10]:
                wx, wy, wz = item["wall_force"]
                handle.write(
                    f"  id={item['id']} x={item['x']:.12e} y={item['y']:.12e} z={item['z']:.12e} radius={item['radius']:.12e} "
                    f"wallForce=({wx:.12e},{wy:.12e},{wz:.12e})\n"
                )

    balance_ok = (
        max_primary_balance <= balance_tol
        and max_full_balance <= balance_tol
        and max_wall_balance_error <= balance_tol
    )
    off_wall_ok = not wall_force_off_wall
    ok = balance_ok and off_wall_ok
    return {
        "dump": dump_path,
        "report": report_path,
        "table": tsv_path,
        "ok": ok,
        "balance_ok": balance_ok,
        "off_wall_ok": off_wall_ok,
        "has_explicit_wallforce": has_explicit_wallforce,
        "max_primary_balance": max_primary_balance,
        "max_full_balance": max_full_balance,
        "max_wall_balance_error": max_wall_balance_error,
        "off_wall_count": len(wall_force_off_wall),
    }


def decode_segment_time(text):
    try:
        return float(text.replace("m", "-").replace("p", "."))
    except ValueError:
        return None


def segment_times_from_path(path):
    parts = os.path.normpath(path).split(os.sep)
    try:
        idx = parts.index("forcebalance_segments")
    except ValueError:
        return False, None, None
    if idx + 1 >= len(parts):
        return True, None, None
    match = re.match(r"^t(.+)_to_(.+)$", parts[idx + 1])
    if not match:
        return True, None, None
    return True, decode_segment_time(match.group(1)), decode_segment_time(match.group(2))


def dump_timestep_from_path(path):
    match = re.search(r"dump([0-9]+)\.liggghts_forcebalance$", os.path.basename(path))
    if not match:
        return -1
    return int(match.group(1))


def dump_sort_key(path):
    in_segment_dir, segment_start, segment_end = segment_times_from_path(path)
    timestep = dump_timestep_from_path(path)
    try:
        mtime = os.path.getmtime(path)
    except OSError:
        mtime = 0.0
    if segment_end is not None:
        return (2, segment_end, segment_start if segment_start is not None else 0.0, timestep, path)
    if in_segment_dir:
        return (1, mtime, 0.0, timestep, path)
    return (0, mtime, 0.0, timestep, path)


def collect_dump_files(case_path):
    patterns = [
        os.path.join(case_path, "DEM", "post", "forcebalance_segments", "*", "dump*.liggghts_forcebalance"),
        os.path.join(case_path, "DEM", "post", "dump*.liggghts_forcebalance"),
    ]
    dump_files = []
    seen = set()
    for pattern in patterns:
        for dump_path in glob.glob(pattern):
            if dump_path not in seen:
                dump_files.append(dump_path)
                seen.add(dump_path)
    return sorted(dump_files, key=dump_sort_key), patterns


def main():
    parser = argparse.ArgumentParser(description="Check per-particle force balance for a CFDEM case dump.")
    parser.add_argument("case_path", help="Absolute path to the case root")
    parser.add_argument("--latest", action="store_true", help="Analyze only the latest dump file")
    parser.add_argument("--all", action="store_true", help="Analyze all matching dump files")
    parser.add_argument("--balance-tol", type=float, default=1e-8)
    parser.add_argument("--wall-force-tol", type=float, default=1e-10)
    parser.add_argument("--wall-gap-tol", type=float, default=1e-10)
    parser.add_argument(
        "--warn-only",
        action="store_true",
        help="Do not fail the command on balance residuals above --balance-tol; still fail on off-wall wallforce.",
    )
    args = parser.parse_args()

    case_path = os.path.abspath(args.case_path)
    variables = parse_variables(os.path.join(case_path, "DEM", "in.liggghts_variable"))
    dump_files, dump_patterns = collect_dump_files(case_path)
    if not dump_files:
        print("No dump files found matching:", file=sys.stderr)
        for pattern in dump_patterns:
            print(f"  {pattern}", file=sys.stderr)
        return 1

    if args.latest and args.all:
        print("Use either --latest or --all, not both", file=sys.stderr)
        return 1
    if args.latest or not args.all:
        dump_files = [dump_files[-1]]

    output_dir = os.path.join(case_path, "DEM", "post", "reports")
    os.makedirs(output_dir, exist_ok=True)

    overall_ok = True
    overall_balance_ok = True
    overall_off_wall_ok = True
    for dump_path in dump_files:
        result = analyze_dump(
            dump_path,
            variables,
            output_dir,
            args.balance_tol,
            args.wall_force_tol,
            args.wall_gap_tol,
        )
        print(f"force-balance report: {result['report']}")
        print(f"force-balance table:  {result['table']}")
        print(
            "summary: "
            f"explicit_wallforce={'yes' if result['has_explicit_wallforce'] else 'no'} "
            f"max_primary_balance={result['max_primary_balance']:.12e} "
            f"max_full_balance={result['max_full_balance']:.12e} "
            f"max_wall_balance_error={result['max_wall_balance_error']:.12e} "
            f"off_wall_count={result['off_wall_count']}"
        )
        if not result["balance_ok"]:
            if args.warn_only:
                print(
                    "force-balance warning: tolerance exceeded but continuing "
                    f"(warn-only mode, balance_tol={args.balance_tol:.12e}) "
                    f"max_primary_balance={result['max_primary_balance']:.12e} "
                    f"max_full_balance={result['max_full_balance']:.12e} "
                    f"max_wall_balance_error={result['max_wall_balance_error']:.12e}"
                )
            else:
                print(
                    "force-balance failure: tolerance exceeded "
                    f"(balance_tol={args.balance_tol:.12e}) "
                    f"max_primary_balance={result['max_primary_balance']:.12e} "
                    f"max_full_balance={result['max_full_balance']:.12e} "
                    f"max_wall_balance_error={result['max_wall_balance_error']:.12e}",
                    file=sys.stderr,
                )
        if not result["off_wall_ok"]:
            print(
                "force-balance failure: explicit wallforce detected away from walls "
                f"(off_wall_count={result['off_wall_count']})",
                file=sys.stderr,
            )
        overall_ok = overall_ok and result["ok"]
        overall_balance_ok = overall_balance_ok and result["balance_ok"]
        overall_off_wall_ok = overall_off_wall_ok and result["off_wall_ok"]

    if overall_ok:
        return 0
    if args.warn_only and overall_off_wall_ok:
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(main())
