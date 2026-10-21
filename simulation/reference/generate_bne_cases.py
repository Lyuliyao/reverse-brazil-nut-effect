#!/usr/bin/env python3
"""Generate the parameter-sweep CFD-DEM cases (Figures 6-9).

Reference implementation used to create the archived cases in cases/. Paths are
taken from the environment: RBNE_CAMPAIGN_ROOT (output root; default runs/sweep)
and RBNE_BASE_CASE (a fully populated template case to copy CFD/DEM inputs from;
default: the bundled canonical case).
"""
import csv
import json
import math
import os
import re
import shutil
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN_ROOT = Path(os.environ.get("RBNE_CAMPAIGN_ROOT", REPO_ROOT / "runs" / "sweep"))
CASES_ROOT = CAMPAIGN_ROOT / "cases"
BASE_CASE = Path(os.environ.get(
    "RBNE_BASE_CASE",
    REPO_ROOT / "cases" / "bne_vf0p1_rs100um_rl300um_rho2470_cor0p1_mu0p5_ang30",
))
MANIFEST_PATH = CAMPAIGN_ROOT / "case_manifest.csv"
MANIFEST_FIELDNAMES = [
    "case_name",
    "case_path",
    "volfrac",
    "initial_state",
    "small_radius",
    "large_radius",
    "rhop",
    "restitution",
    "friction",
    "angle_deg",
    "gx",
    "gy",
    "gz",
    "end_time",
    "write_interval",
    "dump_every",
    "dump_phys_interval",
    "nprocs",
    "g_mag",
]

G_MAG = 9.796
DEFAULT_END_TIME = "5e-4"
DEFAULT_WRITE_INTERVAL = "5e-4"
DEM_TIMESTEP = 1.0e-7
LOW_VOLFRAC_DUMP_PHYS_INTERVAL = 0.01
HIGH_VOLFRAC_DUMP_PHYS_INTERVAL = 0.005
DEFAULT_NPROCS = 8
DEFAULT_BLOCK_CELLS = (32, 32, 16)
DUMP_COLUMNS = (
    "id type mass x y z fx fy fz "
    "f_contactforce[1] f_contactforce[2] f_contactforce[3] "
    "f_wallforce[1] f_wallforce[2] f_wallforce[3] "
    "f_gradPForce[1] f_gradPForce[2] f_gradPForce[3] "
    "f_viscForce[1] f_viscForce[2] f_viscForce[3] "
    "f_fluidforce[1] f_fluidforce[2] f_fluidforce[3] "
    "f_dragforce[1] f_dragforce[2] f_dragforce[3] radius c_cc"
)

VOLFRACS = [0.05, 0.2, 0.4]
RADIUS_PAIRS = [(1.0e-4, 3.0e-4), (1.0e-4, 6.0e-4)]
RHOPS = [2470, 4960]
RESTITUTIONS = [0.9, 0.1]
FRICTIONS = [0.1, 0.5]
ANGLES_DEG = [10, 20, 30]

ROOT_FILES = [
    "Allrun.sh",
    "parDEMrun.sh",
    "parCFDDEMrun.sh",
    "check_force_balance.py",
]

DEM_FILES = [
    "in.liggghts_resume",
    "in.liggghts_setting",
    "in.liggghts_variable",
]

CFD_ZERO_FILES = [
    "Us",
    "U",
    "rho",
    "Ksl",
    "p",
    "voidfraction",
]

CFD_CONSTANT_FILES = [
    "RASProperties",
    "couplingProperties",
    "liggghtsCommands",
    "turbulenceProperties",
    "g",
    "transportProperties",
]

CFD_SYSTEM_FILES = [
    "decomposeParDict",
    "fvSchemes",
    "fvOptions",
    "blockMeshDict",
    "controlDict",
    "fvSolution",
]


def format_name_value(value):
    return str(value).replace(".", "p").replace("-", "m")


def format_float(value):
    text = f"{value:.15g}"
    if "e" not in text and "." not in text:
        text += ".0"
    return text


def case_name(volfrac, small_radius, large_radius, rhop, restitution, friction, angle_deg):
    return (
        f"bne_vf{format_name_value(volfrac)}"
        f"_rs{int(round(small_radius * 1e6))}um"
        f"_rl{int(round(large_radius * 1e6))}um"
        f"_rho{rhop}"
        f"_cor{format_name_value(restitution)}"
        f"_mu{format_name_value(friction)}"
        f"_ang{angle_deg}"
    )


def dump_phys_interval_for_volfrac(volfrac):
    if float(volfrac) <= 0.05 + 1.0e-12:
        return LOW_VOLFRAC_DUMP_PHYS_INTERVAL
    return HIGH_VOLFRAC_DUMP_PHYS_INTERVAL


def dump_every_for_volfrac(volfrac):
    return max(1, int(round(dump_phys_interval_for_volfrac(volfrac) / DEM_TIMESTEP)))


def build_case_metadata(volfrac, small_radius, large_radius, rhop, restitution, friction, angle_deg):
    name = case_name(volfrac, small_radius, large_radius, rhop, restitution, friction, angle_deg)
    case_dir = CASES_ROOT / name
    return {
        "case_name": name,
        "case_path": str(case_dir),
        "volfrac": volfrac,
        "initial_state": "well_mixed",
        "small_radius": small_radius,
        "large_radius": large_radius,
        "rhop": rhop,
        "restitution": restitution,
        "friction": friction,
        "angle_deg": angle_deg,
        "end_time": DEFAULT_END_TIME,
        "write_interval": DEFAULT_WRITE_INTERVAL,
        "dump_every": dump_every_for_volfrac(volfrac),
        "dump_phys_interval": dump_phys_interval_for_volfrac(volfrac),
        "nprocs": DEFAULT_NPROCS,
        "g_mag": G_MAG,
    }


def iter_case_metadata():
    for volfrac in VOLFRACS:
        for small_radius, large_radius in RADIUS_PAIRS:
            for rhop in RHOPS:
                for restitution in RESTITUTIONS:
                    for friction in FRICTIONS:
                        for angle_deg in ANGLES_DEG:
                            yield build_case_metadata(
                                volfrac,
                                small_radius,
                                large_radius,
                                rhop,
                                restitution,
                                friction,
                                angle_deg,
                            )


def metadata_by_case_name():
    return {row["case_name"]: row for row in iter_case_metadata()}


def read_text(path):
    return path.read_text(encoding="utf-8")


def write_text(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def replace_pattern(text, pattern, replacement):
    new_text, count = re.subn(pattern, replacement, text, flags=re.MULTILINE)
    if count == 0:
        raise ValueError(f"Pattern not found: {pattern}")
    return new_text


def build_liggghts_variable(rhop, gx, gz):
    return (
        "variable           xlo equal 0.0\n"
        "variable           xhi equal 1e-2\n"
        "variable           ylo equal 0.0\n"
        "variable           yhi equal 1e-2\n"
        "variable           zlo equal 0.0\n"
        "variable           zhi equal 5e-3\n"
        f"variable           rhop equal {rhop}\n"
        f"variable           gx equal {format_float(gx)}\n"
        "variable           gy equal 0.0\n"
        f"variable           gz equal {format_float(gz)}\n"
    )


def build_liggghts_setting(restitution, friction):
    return (
        "fix  matE   all property/global youngsModulus  peratomtype  5e8     5e8\n"
        "fix  matNu  all property/global poissonsRatio  peratomtype  0.45    0.45\n"
        "\n"
        "fix cor all property/global coefficientRestitution peratomtypepair 2  "
        f"{restitution}  {restitution}  {restitution} {restitution}\n"
        "\n"
        "fix mu  all property/global coefficientFriction    peratomtypepair 2  "
        f"{friction}  {friction}  {friction}  {friction}\n"
        "\n"
        "fix mur all property/global coefficientRollingFriction peratomtypepair 2  0.0  0.0  0.0  0.0\n"
        "\n"
        "pair_style  gran model hertz tangential history rolling_friction cdt\n"
        "pair_coeff  * *\n"
        "\n"
        "timestep            1e-7\n"
        "fix                 grav all gravity 9.796 vector ${gx} ${gy} ${gz}\n"
        "fix                 zwalls1 all wall/gran model hertz tangential history rolling_friction cdt primitive type 1  zplane 0.0\n"
        "fix                 zwalls2 all wall/gran model hertz tangential history rolling_friction cdt primitive type 1 zplane 5e-3\n"
        "\n"
        "\n"
    )


def build_liggghts_init(volfrac, small_radius, large_radius):
    pad = 1.0e-5
    z_insert_lo = 0.0 + large_radius + pad
    z_insert_hi = 5.0e-3 - large_radius - pad
    return (
        "echo            both\n"
        "log             ../DEM/log.liggghts\n"
        "\n"
        "atom_style      granular\n"
        "atom_modify     map array\n"
        "communicate     single vel yes\n"
        "\n"
        "boundary        p p f\n"
        "newton          off\n"
        "\n"
        "units           si\n"
        "\n"
        "include          ../DEM/in.liggghts_variable\n"
        "\n"
        "region             regDomain block ${xlo} ${xhi} ${ylo} ${yhi} ${zlo} ${zhi} units box\n"
        "create_box         2 regDomain\n"
        "\n"
        "include ../DEM/in.liggghts_setting\n"
        "\n"
        f"variable           rSmall equal {format_float(small_radius)}\n"
        f"variable           rLarge equal {format_float(large_radius)}\n"
        f"region             regInsertMix block ${{xlo}} ${{xhi}} ${{ylo}} ${{yhi}} {format_float(z_insert_lo)} {format_float(z_insert_hi)} units box\n"
        "\n"
        "fix                 ptsSmall all particletemplate/sphere 15485867 atom_type 1 density constant ${rhop} radius constant ${rSmall}\n"
        "fix                 ptsLarge all particletemplate/sphere 32452867 atom_type 2 density constant ${rhop} radius constant ${rLarge}\n"
        "fix                 pddMix all particledistribution/discrete 86028121 2 ptsSmall 0.5 ptsLarge 0.5\n"
        f"fix                 insMix all insert/pack seed 49979687 distributiontemplate pddMix vel constant 0. 0. 0. insert_every once region regInsertMix volumefraction_region {format_float(volfrac)} ntry_mc 100000\n"
        "\n"
        "fix                 integr all nve/sphere\n"
        "\n"
        "compute         1 all erotate/sphere\n"
        "thermo_style    custom step atoms ke c_1 vol\n"
        "thermo          1000\n"
        "thermo_modify   lost ignore norm no\n"
        "\n"
        "run     10 upto\n"
        "write_restart post/restart/liggghts.restart\n"
        "\n"
        "\n"
        "\n"
        "\n"
        "\n"
    )


def update_resume_dump(text, dump_every):
    replacement = (
        "dump        dmp all custom "
        f"{dump_every} "
        "../DEM/post/forcebalance_segments/t0_to_initial/dump*.liggghts_forcebalance "
        f"{DUMP_COLUMNS}"
    )
    return replace_pattern(text, r"^dump\s+dmp all custom .*$", replacement)


def update_par_cfd_dem_run(text):
    warn_only_line = 'python3 $casePath/check_force_balance.py "$casePath" --latest --balance-tol 1e-8 --warn-only'

    if 'nrProcs="${CFDEM_NPROCS:-8}"' not in text:
        updated = re.sub(r'^nrProcs=".*"$', 'nrProcs="${CFDEM_NPROCS:-8}"', text, count=1, flags=re.MULTILINE)
        if updated == text:
            raise ValueError('Could not normalize nrProcs line in parCFDDEMrun.sh')
        text = updated
    if 'separateDEM=' not in text:
        text = replace_pattern(
            text,
            r'^(debugMode="off".*)$',
            '\\1\nseparateDEM="false"\nreconstuctCase="true"',
        )
    if '$separateDEM $reconstuctCase' not in text:
        text = replace_pattern(
            text,
            r'^parCFDDEMrun \$logpath \$logfileName \$casePath \$headerText \$solverName \$nrProcs \$machineFileName \$debugMode$',
            'parCFDDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode $separateDEM $reconstuctCase',
        )
    if warn_only_line not in text:
        text = replace_pattern(
            text,
            r'^python3 \$casePath/check_force_balance\.py "\$casePath" --latest --balance-tol 1e-8$',
            warn_only_line,
        )
    return text

def update_control_dict(text, end_time, write_interval):
    text = replace_pattern(text, r"^endTime\s+.*;$", f"endTime         {end_time};")
    text = replace_pattern(text, r"^writeInterval\s+.*;$", f"writeInterval   {write_interval};")
    return text


def update_decompose_dict(text, nprocs):
    return replace_pattern(text, r"^numberOfSubdomains\s+\d+;$", f"numberOfSubdomains {nprocs};")


def update_block_mesh(text, block_cells):
    nx, ny, nz = block_cells
    replacement = f"    hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1)"
    return replace_pattern(
        text,
        r"^\s*hex\s+\(0 1 2 3 4 5 6 7\)\s+\(\d+\s+\d+\s+\d+\)\s+simpleGrading\s+\(1 1 1\)\s*$",
        replacement,
    )


def update_g_file(text, gx, gz, angle_deg):
    comment = f"// {angle_deg} deg: (g sin(theta), 0, -g cos(theta))"
    replacement = f"value           ({format_float(gx)} 0 {format_float(gz)});   {comment}"
    return replace_pattern(text, r"^value\s+\(.*\);\s*.*$", replacement)


def update_fv_options(text, gx, gz):
    replacement = f"            U   ( ({format_float(gx)} 0 {format_float(gz)}) 0 ); "
    return replace_pattern(text, r"^\s*U\s+\(\s*\(.*\)\s*0\s*\)\s*;\s*$", replacement)


def build_liggghts_commands():
    return (
        "/*---------------------------------------------------------------------------*\\\\\n"
        "| =========                 |                                                 |\n"
        "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
        "|  \\\\    /   O peration     | Version:  1.4                                   |\n"
        "|   \\\\  /    A nd           | Web:      http://www.openfoam.org               |\n"
        "|    \\\\/     M anipulation  |                                                 |\n"
        "\\*---------------------------------------------------------------------------*/\n"
        "\n"
        "\n"
        "FoamFile\n"
        "{\n"
        "    version         2.0;\n"
        "    format          ascii;\n"
        "\n"
        "    root            \"\";\n"
        "    case            \"\";\n"
        "    instance        \"\";\n"
        "    local           \"\";\n"
        "\n"
        "    class           dictionary;\n"
        "    object          liggghtsCommands;\n"
        "}\n"
        "\n"
        "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n"
        "\n"
        "liggghtsCommandModels\n"
        "(\n"
        "    runLiggghts\n"
        "    writeLiggghts\n"
        ");\n"
        "\n"
        "runLiggghtsProps\n"
        "{\n"
        "    preNo false;\n"
        "}\n"
        "\n"
        "writeLiggghtsProps\n"
        "{\n"
        "    path \"../DEM/post/restart\";\n"
        "    writeName \"liggghts.restart\";\n"
        "    overwrite off;\n"
        "    writeLastOnly on;\n"
        "}\n"
        "\n"
        "// ************************************************************************* //\n"
    )


def copy_case_skeleton(case_dir):
    if case_dir.exists():
        shutil.rmtree(case_dir)

    (case_dir / "DEM" / "post" / "restart").mkdir(parents=True, exist_ok=True)
    (case_dir / "DEM" / "post" / "reports").mkdir(parents=True, exist_ok=True)
    (case_dir / "DEM" / "post" / "forcebalance_segments" / "t0_to_initial").mkdir(parents=True, exist_ok=True)
    (case_dir / "CFD" / "0").mkdir(parents=True, exist_ok=True)
    (case_dir / "CFD" / "constant").mkdir(parents=True, exist_ok=True)
    (case_dir / "CFD" / "system").mkdir(parents=True, exist_ok=True)

    for filename in ROOT_FILES:
        shutil.copy2(BASE_CASE / filename, case_dir / filename)
    for filename in DEM_FILES:
        shutil.copy2(BASE_CASE / "DEM" / filename, case_dir / "DEM" / filename)
    for filename in CFD_ZERO_FILES:
        shutil.copy2(BASE_CASE / "CFD" / "0" / filename, case_dir / "CFD" / "0" / filename)
    for filename in CFD_CONSTANT_FILES:
        shutil.copy2(BASE_CASE / "CFD" / "constant" / filename, case_dir / "CFD" / "constant" / filename)
    for filename in CFD_SYSTEM_FILES:
        shutil.copy2(BASE_CASE / "CFD" / "system" / filename, case_dir / "CFD" / "system" / filename)


def write_case(case_dir, metadata):
    gx = G_MAG * math.sin(math.radians(metadata["angle_deg"]))
    gz = -G_MAG * math.cos(math.radians(metadata["angle_deg"]))
    metadata["gx"] = gx
    metadata["gy"] = 0.0
    metadata["gz"] = gz

    copy_case_skeleton(case_dir)

    write_text(case_dir / "DEM" / "in.liggghts_variable", build_liggghts_variable(metadata["rhop"], gx, gz))
    write_text(case_dir / "DEM" / "in.liggghts_setting", build_liggghts_setting(metadata["restitution"], metadata["friction"]))
    write_text(
        case_dir / "DEM" / "in.liggghts_init",
        build_liggghts_init(metadata["volfrac"], metadata["small_radius"], metadata["large_radius"]),
    )

    resume_text = read_text(case_dir / "DEM" / "in.liggghts_resume")
    write_text(case_dir / "DEM" / "in.liggghts_resume", update_resume_dump(resume_text, metadata["dump_every"]))

    par_cfd_dem_run = read_text(case_dir / "parCFDDEMrun.sh")
    write_text(case_dir / "parCFDDEMrun.sh", update_par_cfd_dem_run(par_cfd_dem_run))

    control_dict = read_text(case_dir / "CFD" / "system" / "controlDict")
    write_text(
        case_dir / "CFD" / "system" / "controlDict",
        update_control_dict(control_dict, metadata["end_time"], metadata["write_interval"]),
    )

    decompose_dict = read_text(case_dir / "CFD" / "system" / "decomposeParDict")
    write_text(case_dir / "CFD" / "system" / "decomposeParDict", update_decompose_dict(decompose_dict, metadata["nprocs"]))

    block_mesh = read_text(case_dir / "CFD" / "system" / "blockMeshDict")
    write_text(
        case_dir / "CFD" / "system" / "blockMeshDict",
        update_block_mesh(block_mesh, DEFAULT_BLOCK_CELLS),
    )

    g_text = read_text(case_dir / "CFD" / "constant" / "g")
    write_text(case_dir / "CFD" / "constant" / "g", update_g_file(g_text, gx, gz, metadata["angle_deg"]))

    write_text(case_dir / "CFD" / "constant" / "liggghtsCommands", build_liggghts_commands())

    fv_options = read_text(case_dir / "CFD" / "system" / "fvOptions")
    write_text(case_dir / "CFD" / "system" / "fvOptions", update_fv_options(fv_options, gx, gz))

    write_text(case_dir / "case_params.json", json.dumps(metadata, indent=2, sort_keys=True) + "\n")


def write_manifest(rows):
    with MANIFEST_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDNAMES)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def generate_cases(case_names=None):
    CASES_ROOT.mkdir(parents=True, exist_ok=True)
    selected_names = None if case_names is None else set(case_names)
    rows = []

    for metadata in iter_case_metadata():
        if selected_names is not None and metadata["case_name"] not in selected_names:
            continue
        write_case(CASES_ROOT / metadata["case_name"], metadata)
        rows.append(metadata)

    if selected_names is None:
        write_manifest(rows)

    return rows


def main():
    rows = generate_cases()

    print(f"Generated {len(rows)} cases under {CASES_ROOT}")
    print(f"Manifest: {MANIFEST_PATH}")


if __name__ == "__main__":
    main()
