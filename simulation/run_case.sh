#!/usr/bin/env bash
# Run only inside a fresh directory created by prepare_case.py.
set -eo pipefail
if [ "$#" -ne 2 ]; then
    echo "Usage: bash simulation/run_case.sh RUN_DIRECTORY CFDEM_ENV_SCRIPT" >&2
    exit 2
fi
rbne_run_dir="$(cd "$1" && pwd -P)"
rbne_env_script="$(realpath "$2")"
test -f "$rbne_run_dir/prepared_case.json"
if [ -e "$rbne_run_dir/.run-started" ]; then
    echo "This run has already been started; its results are preserved." >&2
    exit 1
fi
source "$rbne_env_script"
command -v blockMesh >/dev/null
command -v decomposePar >/dev/null
command -v cfdemSolverPiso >/dev/null
command -v mpirun >/dev/null
test -x "${CFDEM_LIGGGHTS_EXEC:?CFDEM environment must define CFDEM_LIGGGHTS_EXEC}"
rbne_nprocs="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["nprocs"])' "$rbne_run_dir/prepared_case.json")"
touch "$rbne_run_dir/.run-started"
cd "$rbne_run_dir/CFD"
blockMesh 2>&1 | tee "$rbne_run_dir/log.blockMesh"
cd "$rbne_run_dir/DEM"
mpirun -np 1 "$CFDEM_LIGGGHTS_EXEC" -in in.liggghts_init 2>&1 | tee "$rbne_run_dir/log.initialization"
cd "$rbne_run_dir/CFD"
decomposePar -force 2>&1 | tee "$rbne_run_dir/log.decomposition"
mpirun -np "$rbne_nprocs" cfdemSolverPiso -parallel 2>&1 | tee "$rbne_run_dir/log.coupled"
reconstructPar -noLagrangian 2>&1 | tee "$rbne_run_dir/log.reconstruction"
