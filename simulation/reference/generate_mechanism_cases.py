#!/usr/bin/env python3
"""Generate the mechanism-isolation CFD-DEM cases (Figure 10 controls).

These cases are generated separately from the parameter-sweep cases. It never
overwrites existing cases (refuses to write into an existing case dir unless
--force) and never rmtree's by default.

Differences from the parameter-sweep generator (generate_bne_cases.py):
  * Direct gravity components (gx, gy, gz) with the DEM `fix gravity` MAGNITUDE set to
    g_eff = sqrt(gx^2+gy^2+gz^2)  (LIGGGHTS normalizes the vector then scales by the
    magnitude; fix_gravity.cpp:379-397). OpenFOAM constant/g and fvOptions get the RAW
    (gx,0,gz) acceleration vector.
  * Extended per-particle dump: adds vx vy vz (for species vertical flux and contact work).
  * Per-contact pair dump via `compute pair/gran/local` + `dump local`
    (for contact-type decomposition F^z_{L<-S} vs F^z_{L<-L} and contact work/impulse).
  * Optional per-species vertical body force via `fix addforce` (equal-settling probe).
  * Parameterized insertion seeds (multi-seed; paired ICs reuse one restart per seed).
  * Rich metadata.json (control type, g components, g_eff, seeds, ...).

Notation: case names/metadata record RADII (rs/rl, micron) AND diameters. The
manuscript canonical d_L=600um corresponds to large RADIUS 300um (rl300).

Usage examples (this script only WRITES case dirs; submission is separate):
  python3 generate_mechanism_cases.py --control baseline       --seed-index 1
  python3 generate_mechanism_cases.py --control settling_only  --seed-index 1
  python3 generate_mechanism_cases.py --control shear_only      --seed-index 1
  python3 generate_mechanism_cases.py --control equal_settling  --seed-index 1 \
        --addforce-large <Fz_N>   # upward (+) per-large-particle correction
  python3 generate_mechanism_cases.py --control gravity_split   --seed-index 1 \
        --gpar 0.5 --gperp 0.866
  python3 generate_mechanism_cases.py --control const_ratio --dl-um 900 --ds-um 300 --seed-index 1
  python3 generate_mechanism_cases.py --control settling_velocity_test --species large --seed-index 1
"""
import argparse
import json
import math
import re
import shutil
import subprocess
import os
from pathlib import Path

# --------------------------------------------------------------------------- #
REPO_ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN_ROOT = Path(os.environ.get("RBNE_CAMPAIGN_ROOT", REPO_ROOT / "runs" / "mechanism"))
CASES_ROOT = CAMPAIGN_ROOT / "cases"
# Template = a fully populated canonical case (base material, alpha=30, phi=0.10).
TEMPLATE_CASE = Path(os.environ.get(
    "RBNE_TEMPLATE_CASE",
    REPO_ROOT / "cases" / "bne_vf0p1_rs100um_rl300um_rho2470_cor0p1_mu0p5_ang30",
))

G_MAG = 9.796
DEM_TIMESTEP = 1.0e-7
G_NODE = "9.81"  # informational only

# Canonical defaults (base material). RADII in metres.
DEF_RS = 1.0e-4          # small radius -> d_S = 200 um
DEF_RL = 3.0e-4          # large radius -> d_L = 600 um (manuscript canonical)
DEF_PHI = 0.10
DEF_ANGLE = 30
DEF_RHOP = 2470
DEF_REST = 0.1
DEF_FRIC = 0.5
DEF_RHOF = 1000.0        # fluid density (water); transportProperties uses this for buoyancy

# Domain (matches canonical template)
XLO, XHI, YLO, YHI, ZLO, ZHI = 0.0, 1e-2, 0.0, 1e-2, 0.0, 5e-3

# Dump intervals
DUMP_PHYS_INTERVAL = 0.005                                   # per-particle dump cadence (s)
PARTICLE_DUMP_EVERY = max(1, round(DUMP_PHYS_INTERVAL / DEM_TIMESTEP))   # 50000
CONTACT_DUMP_MULTIPLE = 10                                   # pair/local every k * particle dump (0.05 s)
CONTACT_DUMP_EVERY = PARTICLE_DUMP_EVERY * CONTACT_DUMP_MULTIPLE

# Extended per-particle dump columns: adds vx vy vz right after position.
PARTICLE_DUMP_COLUMNS = (
    "id type mass x y z vx vy vz fx fy fz "
    "f_contactforce[1] f_contactforce[2] f_contactforce[3] "
    "f_wallforce[1] f_wallforce[2] f_wallforce[3] "
    "f_gradPForce[1] f_gradPForce[2] f_gradPForce[3] "
    "f_viscForce[1] f_viscForce[2] f_viscForce[3] "
    "f_fluidforce[1] f_fluidforce[2] f_fluidforce[3] "
    "f_dragforce[1] f_dragforce[2] f_dragforce[3] radius c_cc"
)
# Per-contact pair compute keywords (order in output is fixed by LIGGGHTS, not by us):
# id(3) force(3) force_normal(3) delta(1) contactPoint(3) = 13 cols -> c_cpgl[1..13]
CONTACT_COMPUTE_KEYWORDS = "id force force_normal delta contactPoint"
CONTACT_N_COLS = 13


def fmt(value):
    text = f"{value:.15g}"
    if "e" not in text and "E" not in text and "." not in text:
        text += ".0"
    return text


def name_val(value):
    return str(value).replace(".", "p").replace("-", "m")


# --------------------------------------------------------------------------- #
def gravity_components(control, angle_deg, gpar=None, gperp=None):
    """Return (gx, gy, gz, g_eff) for the requested control.

    Convention: x = downslope, z = wall-normal (positive up). g_par = gx >= 0 is the
    downslope drive; g_perp = -gz > 0 is the wall-normal settling drive.
    """
    a = math.radians(angle_deg)
    if control in ("baseline", "equal_settling", "const_ratio", "same_size_density"):
        gx = G_MAG * math.sin(a)
        gz = -G_MAG * math.cos(a)
    elif control == "settling_only":
        # no downslope drive; keep the physical normal component g*cos(alpha)
        gx = 0.0
        gz = -G_MAG * math.cos(a)
    elif control in ("settling_velocity_test",):
        # pure normal gravity for terminal-velocity measurement
        gx = 0.0
        gz = -G_MAG * math.cos(a)
    elif control == "shear_only":
        # no wall-normal settling drive; keep the physical downslope component g*sin(alpha)
        gx = G_MAG * math.sin(a)
        gz = 0.0
    elif control == "gravity_split":
        if gpar is None or gperp is None:
            raise SystemExit("gravity_split requires --gpar and --gperp (in units of g)")
        gx = gpar * G_MAG
        gz = -gperp * G_MAG
    else:
        raise SystemExit(f"unknown control: {control}")
    g_eff = math.sqrt(gx * gx + gz * gz)
    return gx, 0.0, gz, g_eff


def _is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _next_prime(n):
    n = max(int(n), 10007)
    if n % 2 == 0:
        n += 1
    while not _is_prime(n):
        n += 2
    return n


def seed_set(seed_index):
    """Four DISTINCT primes > 1e7, deterministic per seed_index (paired ICs reuse the set).

    LIGGGHTS requires random seeds that are prime and > 10000 (random.cpp:93).
    """
    cur = 10_000_019 + seed_index * 1_000_003
    seeds = []
    while len(seeds) < 4:
        p = _next_prime(cur)
        seeds.append(p)
        cur = p + 2
    return seeds


# --------------------------------------------------------------------------- #
def build_variable(rhop, gx, gy, gz, g_eff, rhop_small=None, rhop_large=None):
    rs = rhop if rhop_small is None else rhop_small
    rl = rhop if rhop_large is None else rhop_large
    return (
        f"variable           xlo equal {fmt(XLO)}\n"
        f"variable           xhi equal {fmt(XHI)}\n"
        f"variable           ylo equal {fmt(YLO)}\n"
        f"variable           yhi equal {fmt(YHI)}\n"
        f"variable           zlo equal {fmt(ZLO)}\n"
        f"variable           zhi equal {fmt(ZHI)}\n"
        f"variable           rhop equal {rhop}\n"
        f"variable           rhopSmall equal {rs}\n"
        f"variable           rhopLarge equal {rl}\n"
        f"variable           gx equal {fmt(gx)}\n"
        f"variable           gy equal {fmt(gy)}\n"
        f"variable           gz equal {fmt(gz)}\n"
        f"variable           geff equal {fmt(g_eff)}\n"
    )


def build_setting(restitution, friction, addforce_small=0.0, addforce_large=0.0):
    """Contacts + gravity (magnitude g_eff!) + walls + species groups + optional addforce."""
    txt = (
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
        # IMPORTANT: magnitude = g_eff = |(gx,gy,gz)|, because fix gravity normalizes the vector.
        "fix                 grav all gravity ${geff} vector ${gx} ${gy} ${gz}\n"
        "fix                 zwalls1 all wall/gran model hertz tangential history rolling_friction cdt primitive type 1  zplane 0.0\n"
        "fix                 zwalls2 all wall/gran model hertz tangential history rolling_friction cdt primitive type 1 zplane 5e-3\n"
        "\n"
        "# species groups (type 1 = small, type 2 = large)\n"
        "group               gsmall type 1\n"
        "group               glarge type 2\n"
    )
    if abs(addforce_small) > 0 or abs(addforce_large) > 0:
        txt += (
            "\n# --- equal-settling probe: per-species vertical body force (ARTIFICIAL) ---\n"
            f"fix                 addsmall gsmall addforce 0.0 0.0 {fmt(addforce_small)}\n"
            f"fix                 addlarge glarge addforce 0.0 0.0 {fmt(addforce_large)}\n"
        )
    txt += "\n\n"
    return txt


def build_init(phi, small_radius, large_radius, seeds, single_species=None, chi=0.5):
    """Insertion script. single_species in {None,'small','large'} for settling tests.

    chi is the LARGE-species VOLUME (=mass, equal density) fraction of the solids,
    chi = phi_L/phi.  The LIGGGHTS particledistribution/discrete fractions are
    volume/mass-based here (empirically: 0.5/0.5 -> N_S/N_L = (d_L/d_S)^3 -> chi=0.5),
    so the large-template weight is chi and the small-template weight is 1-chi.
    chi=0.5 reproduces the canonical 0.5/0.5 line exactly (baseline unchanged)."""
    pad = 1.0e-5
    z_lo = ZLO + large_radius + pad
    z_hi = ZHI - large_radius - pad
    s_pts, l_pts, s_pdd, s_ins = seeds
    head = (
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
        f"variable           rSmall equal {fmt(small_radius)}\n"
        f"variable           rLarge equal {fmt(large_radius)}\n"
        f"region             regInsertMix block ${{xlo}} ${{xhi}} ${{ylo}} ${{yhi}} {fmt(z_lo)} {fmt(z_hi)} units box\n"
        "\n"
        f"fix                 ptsSmall all particletemplate/sphere {s_pts} atom_type 1 density constant ${{rhopSmall}} radius constant ${{rSmall}}\n"
        f"fix                 ptsLarge all particletemplate/sphere {l_pts} atom_type 2 density constant ${{rhopLarge}} radius constant ${{rLarge}}\n"
    )
    if single_species == "small":
        dist = f"fix                 pddMix all particledistribution/discrete {s_pdd} 1 ptsSmall 1.0\n"
    elif single_species == "large":
        dist = f"fix                 pddMix all particledistribution/discrete {s_pdd} 1 ptsLarge 1.0\n"
    else:
        dist = (f"fix                 pddMix all particledistribution/discrete {s_pdd} 2 "
                f"ptsSmall {fmt(1.0 - chi)} ptsLarge {fmt(chi)}\n")
    tail = (
        dist
        + f"fix                 insMix all insert/pack seed {s_ins} distributiontemplate pddMix vel constant 0. 0. 0. insert_every once region regInsertMix volumefraction_region {fmt(phi)} ntry_mc 100000\n"
        "\n"
        "fix                 integr all nve/sphere\n"
        "\n"
        "compute         1 all erotate/sphere\n"
        "thermo_style    custom step atoms ke c_1 vol\n"
        "thermo          1000\n"
        "thermo_modify   lost ignore norm no\n"
        "\n"
        "run     10 upto\n"
        "# IC snapshot (id+type+pos) for composition (phi,chi) verification and reproducibility\n"
        "write_dump all custom post/ic_snapshot.liggghts id type x y z\n"
        "write_restart post/restart/liggghts.restart\n"
    )
    return head + tail


def build_resume(particle_dump_dir, contact_dump_dir, with_contacts=True):
    """Coupled-run script. Extended particle dump + (optionally) pair/local contact dump."""
    lines = [
        "echo            both",
        "log             ../DEM/log.liggghts",
        "",
        "atom_style\t    granular",
        "atom_modify\t    map array",
        "communicate\t    single vel yes",
        "",
        "boundary        p p f",
        "newton\t\t    off",
        "units\t\t    si",
        "",
        "include         ../DEM/in.liggghts_variable",
        "read_restart \t../DEM/post/restart/liggghts.restart",
        "",
        "neighbor        5.0e-4 bin",
        "neigh_modify    every 1 delay 0 check yes",
        "communicate     single vel yes",
        "newton          off",
        "",
        "include         ../DEM/in.liggghts_setting",
        "",
        "fix\t\tcfd all couple/cfd couple_every 100 mpi",
        "fix\t\tcfd2 all couple/cfd/force/implicit",
        "fix\t\tintegr all nve/sphere",
        "",
        "compute         centerOfMass all com",
        "compute         cc all contact/atom",
        "variable    dragtotalX equal f_cfd2[1]",
        "variable    dragtotalY equal f_cfd2[2]",
        "variable    dragtotalZ equal f_cfd2[3]",
        "",
        "compute\t\t1 all erotate/sphere",
        "thermo_style\tcustom step atoms ke c_1 vol c_centerOfMass[3] v_dragtotalX v_dragtotalY v_dragtotalZ",
        "thermo\t\t10000",
        "thermo_modify\tlost ignore norm no",
        "",
        f"dump        dmp all custom {PARTICLE_DUMP_EVERY} {particle_dump_dir}/dump*.liggghts_forcebalance {PARTICLE_DUMP_COLUMNS}",
    ]
    if with_contacts:
        lines += [
            f"compute     cpgl all pair/gran/local {CONTACT_COMPUTE_KEYWORDS}",
            "dump        dlocal all local " + str(CONTACT_DUMP_EVERY) + f" {contact_dump_dir}/contact.liggghts_local "
            + " ".join(f"c_cpgl[{i+1}]" for i in range(CONTACT_N_COLS)),
        ]
    lines += ["run\t\t    0 "]
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
def patch_openfoam_g(text, gx, gz, label):
    repl = f"value           ({fmt(gx)} 0 {fmt(gz)});   // mechanism probe {label}"
    return re.sub(r"^value\s+\(.*\);\s*.*$", repl, text, flags=re.MULTILINE)


def patch_fvoptions(text, gx, gz):
    repl = f"            U   ( ({fmt(gx)} 0 {fmt(gz)}) 0 ); "
    return re.sub(r"^\s*U\s+\(\s*\(.*\)\s*0\s*\)\s*;\s*$", repl, text, flags=re.MULTILINE)


def patch_controldict(text, end_time, write_interval, purge_write):
    text = re.sub(r"^endTime\s+.*;$", f"endTime         {end_time};", text, flags=re.MULTILINE)
    text = re.sub(r"^writeInterval\s+.*;$", f"writeInterval   {write_interval};", text, flags=re.MULTILINE)
    text = re.sub(r"^purgeWrite\s+.*;$", f"purgeWrite      {purge_write};", text, flags=re.MULTILINE)
    return text


def patch_coupling_ddt(text, mode, steps=None):
    """Insert `useDDTvoidfraction <mode>;` (a|b|c|off) after the modelType line. Off/None = solver
    default (keyword absent -> cfdemCloud treats it as 'off'). Mode c (requires a locally patched solver build)
    additionally writes `ddtVoidfractionSteps <K>;` — the production solver FatalErrors
    on mode c, which is the intended fail-closed behaviour."""
    if mode is None:
        return text
    ins = f"useDDTvoidfraction {mode};   // control: d(eps)/dt source in continuity"
    if mode == "c":
        if not steps or steps < 2:
            raise SystemExit("--ddt-voidfraction c requires --ddt-voidfraction-steps >= 2")
        ins += f"\nddtVoidfractionSteps {steps};   // mode c: running-mean window (CFD steps)"
    if re.search(r"^\s*useDDTvoidfraction\b", text, flags=re.MULTILINE):
        text = re.sub(r"^\s*ddtVoidfractionSteps\s+\S+;.*$\n?", "", text, flags=re.MULTILINE)
        return re.sub(r"^\s*useDDTvoidfraction\s+\S+;.*$", ins, text, flags=re.MULTILINE)
    return re.sub(r'^(modelType\s+"[AB]"\s*;.*)$',
                  lambda m: m.group(1) + "\n" + ins,
                  text, count=1, flags=re.MULTILINE)


def patch_blockmesh(text, mesh):
    """Replace the (nx ny nz) cell counts of the single hex block; None = template mesh unchanged."""
    if mesh is None:
        return text
    nx, ny, nz = mesh
    new, n = re.subn(r"(hex\s*\([^)]*\)\s*)\(\s*\d+\s+\d+\s+\d+\s*\)",
                     lambda m: m.group(1) + f"({nx} {ny} {nz})", text, count=1)
    if n != 1:
        raise SystemExit("blockMeshDict hex line not found for mesh patch")
    return new


def patch_coupling_debug(text, on):
    """Insert `debug on;` after the modelType line (cfdemCloud.C:306 reads top-level `debug`), which makes
    the solver write ddtVoidfraction (cfdemCloud.C:309) and submodel debug fields; None = unchanged."""
    if not on:
        return text
    return re.sub(r'^(modelType\s+"[AB]"\s*;.*)$',
                  lambda m: m.group(1) + "\ndebug on;   // diagnostics: write ddtVoidfraction etc.",
                  text, count=1, flags=re.MULTILINE)


def patch_write_precision(text, prec):
    """controlDict writePrecision; None = template default (6)."""
    if prec is None:
        return text
    return re.sub(r"^writePrecision\s+.*;$", f"writePrecision  {prec};", text, flags=re.MULTILINE)


def patch_io_model(text, model):
    """couplingProperties IOModel (template: basicIO); None = unchanged. noIO removes CFDEM's global-path
    lagrangian output (IOModel.C:65 mkDir on the GLOBAL case time dir from every rank -> mkdir race with
    per-step writes)."""
    if model is None:
        return text
    return re.sub(r"^IOModel\s+\S+;.*$", f'IOModel "{model}";   // avoid global-path mkdir race', text,
                  count=1, flags=re.MULTILINE)


def patch_time_format(text, tfmt, tprec):
    if tfmt is not None:
        text = re.sub(r"^timeFormat\s+.*;$", f"timeFormat      {tfmt};", text, flags=re.MULTILINE)
    if tprec is not None:
        text = re.sub(r"^timePrecision\s+.*;$", f"timePrecision   {tprec};", text, flags=re.MULTILINE)
    return text


SERIAL_PARCFDDEM = """#!/bin/bash
# TWO-RANK coupled run (diagnostics). Constraints: OpenFOAM-6 aborts -parallel on
# 1 rank (UPstream.C:110); CFDEM twoWayMPI has no MPI_Initialized guard (needs mpirun+-parallel); the
# chain runner rewrites numberOfSubdomains from CFDEM_NPROCS. So: CFDEM_NPROCS=2 -> decomposePar (2) ->
# mpirun -np 2 -parallel -> rc-checked reconstructPar -> only then rm processor*.
# Written by generate_mechanism_cases.py --serial-cfddem; replaces the template parCFDDEMrun.sh.
source "${CFDEM_ENV_SCRIPT:?set CFDEM_ENV_SCRIPT to your CFDEM environment script}" || { echo "env source FAILED"; exit 90; }
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
logfile="$casePath/log_run_parallel_cfdemSolverPiso_periodicChannel_CFDDEM"
nprocs="${CFDEM_NPROCS:-2}"
# ---- fail-closed binary/lib/rank assertions ----
if [ "$nprocs" != "2" ]; then echo "ABORT: CFDEM_NPROCS=$nprocs != 2" | tee -a "$logfile"; exit 91; fi
SOLVER_BIN="$(command -v cfdemSolverPiso)" || { echo "ABORT: solver not on PATH" | tee -a "$logfile"; exit 92; }
SOLVER_SHA=$(sha256sum "$SOLVER_BIN" | cut -d" " -f1)
LIB_PATH=$(ldd "$SOLVER_BIN" | awk "/lagrangianCFDEM/ {print \$3}")
LIB_SHA=$(sha256sum "$LIB_PATH" | cut -d" " -f1)
echo "resolved solver=$SOLVER_BIN sha256=$SOLVER_SHA" | tee -a "$logfile"
echo "resolved lib=$LIB_PATH sha256=$LIB_SHA" | tee -a "$logfile"
EXPECT_SOLVER_SHA="__EXPECT_SOLVER_SHA__"
EXPECT_LIB_SHA="__EXPECT_LIB_SHA__"
if [ -n "$EXPECT_SOLVER_SHA" ] && [ "$SOLVER_SHA" != "$EXPECT_SOLVER_SHA" ]; then
    echo "ABORT: solver sha256 mismatch (expected $EXPECT_SOLVER_SHA)" | tee -a "$logfile"; exit 93; fi
if [ -n "$EXPECT_LIB_SHA" ] && [ "$LIB_SHA" != "$EXPECT_LIB_SHA" ]; then
    echo "ABORT: lib sha256 mismatch (expected $EXPECT_LIB_SHA)" | tee -a "$logfile"; exit 94; fi
cd "$casePath/CFD"
echo "TWO-RANK cfdemSolverPiso run (np=$nprocs)" | tee -a "$logfile"
decomposePar -force 2>&1 | tee -a "$logfile"
if [ ${PIPESTATUS[0]} -ne 0 ]; then echo "decomposePar failed" | tee -a "$logfile"; exit 1; fi
mpirun -np "$nprocs" cfdemSolverPiso -parallel 2>&1 | tee -a "$logfile"
rc=${PIPESTATUS[0]}
if [ $rc -ne 0 ]; then echo "solver failed rc=$rc" | tee -a "$logfile"; exit $rc; fi
reconstructPar -noLagrangian 2>&1 | tee -a "$logfile"
rrc=${PIPESTATUS[0]}
nglob=$(ls -d [0-9]*[0-9.] 2>/dev/null | wc -l)
if [ $rrc -eq 0 ] && [ "$nglob" -gt 2 ]; then
    echo "reconstructPar OK ($nglob global time dirs): removing processor dirs" | tee -a "$logfile"
    rm -rf processor*
else
    echo "reconstructPar rc=$rrc nglob=$nglob: KEEPING processor dirs" | tee -a "$logfile"
    exit 1
fi
python3 "$casePath/check_force_balance.py" "$casePath" --latest --balance-tol 1e-8 --warn-only
exit $?
"""


SMOOTH_FIELD_S = """FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      sSmoothField;
}
dimensions      [0 0 0 0 0 0 0];
internalField   uniform 0;
boundaryField
{
    "(xMin|xMax|yMin|yMax)"
    {
        type            cyclic;
    }
    "(top|bottom)"
    {
        type            zeroGradient;
    }
}
"""

SMOOTH_FIELD_V = SMOOTH_FIELD_S.replace("volScalarField", "volVectorField").replace(
    "sSmoothField", "vSmoothField").replace("uniform 0;", "uniform (0 0 0);")


def patch_coupling_smoothing(text, length):
    """Optional voidfraction smoothing: smoothingModel constDiffSmoothing with
    smoothingLength <length> (m), lowerLimit 0.0, upperLimit 1.0. None = template unchanged.
    NOTE: this smooths voidFractionNext BEFORE calcDdtVoidfraction, so every model sees the
    smoothed eps."""
    if length is None:
        return text
    new, n = re.subn(r"^smoothingModel\s+[^;]+;.*$",
                     "smoothingModel constDiffSmoothing;   // voidfraction smoothing", text,
                     count=1, flags=re.MULTILINE)
    if n != 1:
        raise SystemExit("smoothingModel line not found for smoothing patch")
    props = ("\nconstDiffSmoothingProps\n{\n"
             "    lowerLimit 0.0;\n    upperLimit 1.0;\n"
             f"    smoothingLength {fmt(length)};\n}}\n")
    return new + props


def patch_divided_alphamin(text, alpha_min):
    """Set dividedProps { alphaMin <v>; }; None = model default (0.1, dividedVoidFraction.C:67)."""
    if alpha_min is None:
        return text
    new, n = re.subn(r"dividedProps\s*\{[^}]*\}",
                     f"dividedProps\n{{\n    alphaMin {fmt(alpha_min)};   // control: void-fraction floor\n}}",
                     text, count=1)
    if n != 1:
        raise SystemExit("dividedProps block not found for alphaMin patch")
    return new


def template_mesh():
    txt = (TEMPLATE_CASE / "CFD" / "system" / "blockMeshDict").read_text(encoding="utf-8")
    m = re.search(r"hex\s*\([^)]*\)\s*\(\s*(\d+)\s+(\d+)\s+(\d+)\s*\)", txt)
    return [int(m.group(1)), int(m.group(2)), int(m.group(3))] if m else None


HEAVY = ("CFD/processor", "CFD/0.", "CFD/1", "CFD/2", "CFD/3", "CFD/4", "CFD/5", "CFD/6",
         "CFD/7", "CFD/8", "CFD/9", "DEM/post", "slurm", "analysis")


def clone_template(case_dir, force):
    if case_dir.exists():
        if not force:
            raise SystemExit(f"REFUSING to overwrite existing case dir: {case_dir} (use --force)")
        shutil.rmtree(case_dir)
    # Copy ONLY the skeleton (inputs), not the heavy CFD time dirs / processors / post.
    # This is fast and avoids transient disk churn.
    (case_dir / "DEM" / "post" / "restart").mkdir(parents=True, exist_ok=True)
    (case_dir / "DEM" / "post" / "reports").mkdir(parents=True, exist_ok=True)
    (case_dir / "CFD" / "0").mkdir(parents=True, exist_ok=True)
    (case_dir / "CFD" / "constant").mkdir(parents=True, exist_ok=True)
    (case_dir / "CFD" / "system").mkdir(parents=True, exist_ok=True)
    (case_dir / "slurm").mkdir(parents=True, exist_ok=True)
    # root run scripts + force-balance checker
    for fn in ("Allrun.sh", "parDEMrun.sh", "parCFDDEMrun.sh", "check_force_balance.py"):
        src = TEMPLATE_CASE / fn
        if src.exists():
            shutil.copy2(src, case_dir / fn)
    # DEM inputs (will be overwritten by the generator, but copy for completeness)
    for fn in ("in.liggghts_resume", "in.liggghts_setting", "in.liggghts_variable", "in.liggghts_init"):
        src = TEMPLATE_CASE / "DEM" / fn
        if src.exists():
            shutil.copy2(src, case_dir / "DEM" / fn)
    # CFD/0 fields, constant (incl. polyMesh skipped — blockMesh regenerates), system dicts
    for fn in (TEMPLATE_CASE / "CFD" / "0").iterdir():
        if fn.is_file():
            shutil.copy2(fn, case_dir / "CFD" / "0" / fn.name)
    for fn in (TEMPLATE_CASE / "CFD" / "constant").iterdir():
        if fn.is_file():
            shutil.copy2(fn, case_dir / "CFD" / "constant" / fn.name)
    for fn in (TEMPLATE_CASE / "CFD" / "system").iterdir():
        if fn.is_file():
            shutil.copy2(fn, case_dir / "CFD" / "system" / fn.name)


def git_status():
    return "git: not a repo"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--control", required=True,
                    choices=["baseline", "settling_only", "shear_only", "gravity_split",
                             "equal_settling", "const_ratio", "settling_velocity_test",
                             "same_size_density"])
    ap.add_argument("--seed-index", type=int, default=1)
    ap.add_argument("--angle", type=int, default=DEF_ANGLE)
    ap.add_argument("--phi", type=float, default=DEF_PHI)
    ap.add_argument("--chi", type=float, default=0.5,
                    help="large-species VOLUME fraction of solids chi=phi_L/phi (default 0.5 = canonical)")
    ap.add_argument("--ds-um", type=float, default=None, help="small DIAMETER in um (default 200)")
    ap.add_argument("--dl-um", type=float, default=None, help="large DIAMETER in um (default 600)")
    ap.add_argument("--rhop", type=float, default=DEF_RHOP)
    ap.add_argument("--rhop-small", type=float, default=None,
                    help="density of species 1 (type 1). For same_size_density (§5.7).")
    ap.add_argument("--rhop-large", type=float, default=None,
                    help="density of species 2 (type 2). For same_size_density (§5.7).")
    ap.add_argument("--restitution", type=float, default=DEF_REST)
    ap.add_argument("--friction", type=float, default=DEF_FRIC)
    ap.add_argument("--gpar", type=float, default=None, help="g_par/g for gravity_split")
    ap.add_argument("--gperp", type=float, default=None, help="g_perp/g for gravity_split")
    ap.add_argument("--addforce-small", type=float, default=0.0, help="per-small-particle Fz (N)")
    ap.add_argument("--addforce-large", type=float, default=0.0, help="per-large-particle Fz (N)")
    ap.add_argument("--species", choices=["small", "large"], default=None,
                    help="for settling_velocity_test: which single species")
    ap.add_argument("--phi-test", type=float, default=0.01, help="dilute phi for settling test")
    ap.add_argument("--end-time", type=float, default=None, help="CFD endTime (s)")
    ap.add_argument("--write-interval", type=float, default=0.5)
    ap.add_argument("--purge-write", type=int, default=0, help="0 keeps all CFD writes (for J_rel)")
    ap.add_argument("--name", default=None, help="explicit case name override")
    ap.add_argument("--force", action="store_true")
    # --- CFD-coupling controls (defaults reproduce the archived cases byte-for-byte) ---
    ap.add_argument("--ddt-voidfraction", choices=["off", "a", "b", "c"], default=None,
                    help="CFDEM useDDTvoidfraction: a = div(Us(1-eps)), b = ddt(eps), c = running mean (patched build only), off; default = keyword absent")
    ap.add_argument("--ddt-voidfraction-steps", type=int, default=None,
                    help="mode c: running-mean window K in CFD steps (required with c)")
    ap.add_argument("--smoothing-length", type=float, default=None,
                    help="voidfraction smoothing: constDiffSmoothing smoothingLength (m); default = template smoothingModel off")
    ap.add_argument("--env-script", default=None,
                    help="bake this CFDEM env script path as the parCFDDEMrun.sh default")
    ap.add_argument("--expect-solver-sha", default=None,
                    help="fail-closed: abort the node-side run unless the resolved cfdemSolverPiso sha256 matches")
    ap.add_argument("--expect-lib-sha", default=None,
                    help="fail-closed: abort unless the ldd-resolved liblagrangianCFDEM sha256 matches")
    ap.add_argument("--cfd-mesh", type=int, nargs=3, metavar=("NX", "NY", "NZ"), default=None,
                    help="blockMesh cell counts (default: template 32 32 16)")
    ap.add_argument("--alpha-min", type=float, default=None,
                    help="divided void-fraction model floor alphaMin (default: model default 0.1)")
    ap.add_argument("--cases-root", default=None, help="override the cases root (tests / scratch)")
    ap.add_argument("--coupling-debug", action="store_true", help="couplingProperties `debug on;` (writes ddtVoidfraction)")
    ap.add_argument("--write-precision", type=int, default=None, help="controlDict writePrecision (default template 6)")
    ap.add_argument("--io-model", default=None, help="couplingProperties IOModel; NOTE this build registers the no-op model as \"off\" (source dir noIO); default template basicIO")
    ap.add_argument("--time-format", default=None, choices=["general", "fixed", "scientific"])
    ap.add_argument("--time-precision", type=int, default=None)
    ap.add_argument("--serial-cfddem", action="store_true",
                    help="replace parCFDDEMrun.sh with a SERIAL coupled run (global writes only; diagnostics)")
    args = ap.parse_args()

    rs = (args.ds_um * 1e-6 / 2.0) if args.ds_um else DEF_RS
    rl = (args.dl_um * 1e-6 / 2.0) if args.dl_um else DEF_RL
    dS_um = rs * 2e6
    dL_um = rl * 2e6
    phi = args.phi_test if args.control == "settling_velocity_test" else args.phi
    end_time = args.end_time if args.end_time is not None else (
        2.0 if args.control == "settling_velocity_test" else 40.0)

    gx, gy, gz, g_eff = gravity_components(args.control, args.angle, args.gpar, args.gperp)

    # Case name
    if args.name:
        name = args.name
    else:
        tag = args.control
        extra = ""
        if args.control == "gravity_split":
            extra = f"_gpar{name_val(round(args.gpar,3))}_gperp{name_val(round(args.gperp,3))}"
        if args.control == "settling_velocity_test":
            extra = f"_{args.species}"
        name = (f"mech_{tag}{extra}_dS{int(round(dS_um))}_dL{int(round(dL_um))}"
                f"_phi{name_val(phi)}_ang{args.angle}_seed{args.seed_index}")

    cases_root = Path(args.cases_root).resolve() if args.cases_root else CASES_ROOT
    case_dir = cases_root / name
    clone_template(case_dir, args.force)

    seeds = seed_set(args.seed_index)
    single = args.species if args.control == "settling_velocity_test" else None
    with_contacts = args.control != "settling_velocity_test"

    # --- DEM inputs ---
    (case_dir / "DEM" / "in.liggghts_variable").write_text(
        build_variable(args.rhop, gx, gy, gz, g_eff,
                       rhop_small=args.rhop_small, rhop_large=args.rhop_large), encoding="utf-8")
    (case_dir / "DEM" / "in.liggghts_setting").write_text(
        build_setting(args.restitution, args.friction, args.addforce_small, args.addforce_large),
        encoding="utf-8")
    (case_dir / "DEM" / "in.liggghts_init").write_text(
        build_init(phi, rs, rl, seeds, single_species=single, chi=args.chi), encoding="utf-8")
    pdir = "../DEM/post/forcebalance_segments/t0_to_initial"
    cdir = "../DEM/post/contact_segments/t0_to_initial"
    (case_dir / "DEM" / "post" / "forcebalance_segments" / "t0_to_initial").mkdir(parents=True, exist_ok=True)
    if with_contacts:
        (case_dir / "DEM" / "post" / "contact_segments" / "t0_to_initial").mkdir(parents=True, exist_ok=True)
    (case_dir / "DEM" / "in.liggghts_resume").write_text(
        build_resume(pdir, cdir, with_contacts=with_contacts), encoding="utf-8")

    # --- CFD inputs ---
    gpath = case_dir / "CFD" / "constant" / "g"
    gpath.write_text(patch_openfoam_g(gpath.read_text(encoding="utf-8"), gx, gz, name), encoding="utf-8")
    fvo = case_dir / "CFD" / "system" / "fvOptions"
    fvo.write_text(patch_fvoptions(fvo.read_text(encoding="utf-8"), gx, gz), encoding="utf-8")
    cd = case_dir / "CFD" / "system" / "controlDict"
    cd_txt = patch_controldict(cd.read_text(encoding="utf-8"), fmt(end_time),
                               fmt(args.write_interval), args.purge_write)
    cd_txt = patch_write_precision(cd_txt, args.write_precision)
    cd.write_text(patch_time_format(cd_txt, args.time_format, args.time_precision), encoding="utf-8")
    cp = case_dir / "CFD" / "constant" / "couplingProperties"
    cp_txt = patch_coupling_ddt(cp.read_text(encoding="utf-8"), args.ddt_voidfraction,
                                args.ddt_voidfraction_steps)
    cp_txt = patch_coupling_smoothing(cp_txt, args.smoothing_length)
    cp_txt = patch_coupling_debug(cp_txt, args.coupling_debug)
    if args.io_model == "noIO":      # normalise the alias: the no-op IOModel registers as "off"
        args.io_model = "off"
    cp_txt = patch_io_model(cp_txt, args.io_model)
    cp_txt = patch_divided_alphamin(cp_txt, args.alpha_min)
    if cp_txt != cp.read_text(encoding="utf-8"):
        cp.write_text(cp_txt, encoding="utf-8")
    bm = case_dir / "CFD" / "system" / "blockMeshDict"
    bm_txt = patch_blockmesh(bm.read_text(encoding="utf-8"), args.cfd_mesh)
    if bm_txt != bm.read_text(encoding="utf-8"):
        bm.write_text(bm_txt, encoding="utf-8")

    if args.serial_cfddem:
        run_txt = SERIAL_PARCFDDEM
        if args.env_script:
            run_txt = run_txt.replace(
                "${CFDEM_ENV_SCRIPT:?set CFDEM_ENV_SCRIPT to your CFDEM environment script}",
                "${CFDEM_ENV_SCRIPT:-" + args.env_script + "}")
        run_txt = run_txt.replace("__EXPECT_SOLVER_SHA__", args.expect_solver_sha or "")
        run_txt = run_txt.replace("__EXPECT_LIB_SHA__", args.expect_lib_sha or "")
        (case_dir / "parCFDDEMrun.sh").write_text(run_txt, encoding="utf-8")
        (case_dir / "parCFDDEMrun.sh").chmod(0o755)
    if args.smoothing_length is not None:
        (case_dir / "CFD" / "0" / "sSmoothField").write_text(SMOOTH_FIELD_S, encoding="utf-8")
        (case_dir / "CFD" / "0" / "vSmoothField").write_text(SMOOTH_FIELD_V, encoding="utf-8")
        # constDiffSmoothing solves fvm equations for the smooth fields -> solver entries required
        fvs = case_dir / "CFD" / "system" / "fvSolution"
        txt = fvs.read_text(encoding="utf-8")
        if "sSmoothField" not in txt:
            import re as _re
            txt2, n = _re.subn(r"(solvers\s*\{)",
                r"""\1
    "(sSmoothField|vSmoothField|fSmooth)"
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-09;
        relTol          1e-06;
    }
""", txt, count=1)
            if n != 1:
                raise SystemExit("fvSolution solvers block not found for smoothing patch")
            fvs.write_text(txt2, encoding="utf-8")

    # --- metadata ---
    meta = {
        "case_name": name,
        "campaign": "mechanism_campaign_20260611",
        "control_type": args.control,
        "ddt_voidfraction": args.ddt_voidfraction if args.ddt_voidfraction else "off (solver default; keyword absent)",
        "ddt_voidfraction_steps": args.ddt_voidfraction_steps,
        "smoothing_length": args.smoothing_length,
        "env_script": args.env_script or "default ($CFDEM_ENV_SCRIPT)",
        "expect_solver_sha": args.expect_solver_sha,
        "expect_lib_sha": args.expect_lib_sha,
        "cfd_mesh": list(args.cfd_mesh) if args.cfd_mesh else template_mesh(),
        "alpha_min": args.alpha_min if args.alpha_min is not None else 0.1,
        "alpha_min_source": "dividedProps" if args.alpha_min is not None else "model default (dividedVoidFraction.C:67)",
        "coupling_debug": bool(args.coupling_debug),
        "write_precision": args.write_precision if args.write_precision is not None else 6,
        "io_model": args.io_model if args.io_model else "basicIO (template)",
        "serial_cfddem": bool(args.serial_cfddem),
        "time_format": args.time_format or "general (template)", "time_precision": args.time_precision if args.time_precision is not None else 6,
        "template_case": str(TEMPLATE_CASE),
        "git": git_status(),
        "seed_index": args.seed_index,
        "seeds": seeds,
        "small_radius": rs, "large_radius": rl,
        "small_diameter_um": dS_um, "large_diameter_um": dL_um,
        "size_ratio_dL_over_dS": dL_um / dS_um,
        "volfrac": phi,
        "chi_large_volfrac": args.chi,
        "insert_frac_small": 1.0 - args.chi,
        "insert_frac_large": args.chi,
        "angle_deg": args.angle,
        "rhop": args.rhop, "rhof": DEF_RHOF,
        "rhop_small": args.rhop_small if args.rhop_small is not None else args.rhop,
        "rhop_large": args.rhop_large if args.rhop_large is not None else args.rhop,
        "restitution": args.restitution, "friction": args.friction,
        "g_mag": G_MAG, "gx": gx, "gy": gy, "gz": gz, "g_eff": g_eff,
        "g_par": gx, "g_perp": -gz,
        "addforce_small_z": args.addforce_small, "addforce_large_z": args.addforce_large,
        "single_species": single,
        "dem_timestep": DEM_TIMESTEP,
        "particle_dump_every": PARTICLE_DUMP_EVERY,
        "contact_dump_every": CONTACT_DUMP_EVERY if with_contacts else None,
        "dump_phys_interval": DUMP_PHYS_INTERVAL,
        "particle_dump_columns": PARTICLE_DUMP_COLUMNS,
        "contact_compute_keywords": CONTACT_COMPUTE_KEYWORDS if with_contacts else None,
        "contact_n_cols": CONTACT_N_COLS if with_contacts else None,
        "end_time": end_time, "write_interval": args.write_interval, "purge_write": args.purge_write,
        "domain": {"xlo": XLO, "xhi": XHI, "ylo": YLO, "yhi": YHI, "zlo": ZLO, "zhi": ZHI},
    }
    (case_dir / "metadata.json").write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    # also write case_params.json for compatibility with existing analysis tools
    (case_dir / "case_params.json").write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"WROTE {case_dir}")
    print(f"  control={args.control} gx={gx:.6g} gz={gz:.6g} g_eff={g_eff:.6g} "
          f"dS={dS_um:.0f}um dL={dL_um:.0f}um phi={phi} ang={args.angle} seed_index={args.seed_index}")


if __name__ == "__main__":
    main()
