# CFD–DEM simulations

The compact data bundle is sufficient to replot the paper. New trajectories
require the coupled solver stack below and substantial computation.

## Solver versions and local modifications

| Component | Upstream repository | Recorded commit |
| --- | --- | --- |
| OpenFOAM 6 | [OpenFOAM/OpenFOAM-6](https://github.com/OpenFOAM/OpenFOAM-6) | `af7d7f427be78e9b9beb6aceca8fe7d5d4636876` |
| CFDEMcoupling PUBLIC | [CFDEMproject/CFDEMcoupling-PUBLIC](https://github.com/CFDEMproject/CFDEMcoupling-PUBLIC) | `abedb0712debed5c055be5aa7262fef62fcf3a62` |
| LIGGGHTS PUBLIC | [CFDEMproject/LIGGGHTS-PUBLIC](https://github.com/CFDEMproject/LIGGGHTS-PUBLIC) | `3d5c00f20519e6bb6eb6756f51f1ad36564e649d` |

Use consistent compiler and MPI libraries across all three components.
The source installation used GCC 13.2.0 and OpenMPI 4.1.6. Configure each
upstream installation for your Linux / cluster environment, following its
installation documentation; solver binaries are not bundled.

Before building, check out the recorded commits and apply the patches from
this repository. With `REPO` set to the absolute path of this checkout, run
these commands from the respective **upstream** checkouts:

```bash
# OpenFOAM-6 checkout: compatibility for the recorded build environment
 git apply --check "$REPO/simulation/patches/openfoam6-build-compat.patch"
 git apply "$REPO/simulation/patches/openfoam6-build-compat.patch"
# CFDEMcoupling-PUBLIC checkout: diagnostic force exchange
 git apply --check "$REPO/simulation/patches/cfdem-force-diagnostics.patch"
 git apply "$REPO/simulation/patches/cfdem-force-diagnostics.patch"
# LIGGGHTS-PUBLIC checkout: per-particle contact / wall / fluid force output
 git apply --check "$REPO/simulation/patches/liggghts-force-diagnostics.patch"
 git apply "$REPO/simulation/patches/liggghts-force-diagnostics.patch"
```

`cfdem-run-compat.patch` records the source installation's MPI / Slurm helper
changes. It is optional for the direct `mpirun` runner supplied here. The
patch manifest records each modified source path and its base commit.
The patches were checked against the recorded upstream source files during
repository preparation. A new solver build and a full CFD–DEM run have not
been performed as part of preparing this repository.

Create an environment script that loads your compiler / MPI, sources
OpenFOAM and CFDEMcoupling, and sets the matching LIGGGHTS executable and
library paths. It must expose `blockMesh`, `decomposePar`, `cfdemSolverPiso`,
`reconstructPar`, `mpirun`, and `CFDEM_LIGGGHTS_EXEC`.

## Prepare and run a fresh case

Original archived inputs are in `cases/`. Some archive control dictionaries
retain their last segment settings. **Prepare a fresh copy before running**:

```bash
python simulation/prepare_case.py \
  --case bne_vf0p1_rs100um_rl300um_rho2470_cor0p1_mu0p5_ang30 \
  --out runs/canonical --end-time 40 --nprocs 8
bash simulation/run_case.sh runs/canonical /absolute/path/to/your/cfdem_env.sh
```

`prepare_case.py` preserves the recorded particle sizes, material parameters,
gravity, contact law, force models, and insertion seeds. It resets only the
execution settings for a fresh run: CFD start time, target end time, write /
purge settings, process count, restart location, and dump output paths.
The runner generates the initial DEM restart using one MPI rank, then runs
the coupled CFD–DEM solver with the requested ranks. It preserves its output
and refuses to start a directory that has already been started. The script
does not submit or chain scheduler jobs; execute it inside a suitable
allocation on a cluster.

For Figure 10, choose one of the `mech_*_seed1` cases. The supplied
`mech_equal_settling_stokes_*` input already contains its analytic body force.
`tools/compute_equal_settling_addforce.py --help` explains the calculation.
`reference/` retains the original campaign-generation implementation for
provenance; these historical scripts include campaign-specific paths and
options beyond the paper's selected cases.

Fresh runs can vary with the MPI partition, solver / compiler build, and
insertion implementation. Use the bundled recorded data for direct
comparison with the manuscript figures.

## Analyze new outputs

After a run has finished, preserve the raw dumps while packing and checking:

```bash
BNE_PACK_ALLOW_LATEST_SEGMENT=1 python analysis/raw/pack_forcebalance_segment.py \
  runs/canonical/DEM/post/forcebalance_segments/t0_to_40 --keep-dumps
python analysis/raw/analyze_bne_case.py runs/canonical --full-rebuild
python analysis/raw/compute_mean_height_metrics.py --case runs/canonical
```

The first two tools implement the original force / height analysis. The last
one implements the mechanism-campaign mean-height metrics. The paper figure
scripts read the curated `data/` bundle; new run outputs are written under
`runs/` and do not replace those reference data automatically.
