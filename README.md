# Transient Reverse Brazil Nut Effect in Particle-Laden Inclined Flows

Code and essential data accompanying the paper by **Liyao Lyu, Lingyun Ding,
Yifan Gu, Yanming Ba, and Andrea Bertozzi**.

The author-supplied final manuscript is [paper/reverse_Brazil_nut_effect.pdf](paper/reverse_Brazil_nut_effect.pdf).
This repository contains a compact selection of the original CFD–DEM research
project, organized around the figures in that manuscript.

## Reproduce the numerical figures

Python 3.9 or later is required. From a fresh checkout:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/validate_bundle.py
python scripts/reproduce_figures.py
```

The last command writes Figures 6–10 as PDFs and PNGs to `outputs/figures/`.
To select individual figures or change the output directory:

```bash
python scripts/reproduce_figures.py --figures 6 9 10 --out outputs/selected
```

The figure scripts retain the original plotting conventions, with portable
readers for the bundled data. Exact figure images embedded in the final paper
are retained separately in `paper/figure_assets/` for reference.

## Figure-to-data map

| Paper figure | Contents | Code / data |
| --- | --- | --- |
| 1 | Experimental schematic | `paper/figure_assets/figure01.pdf` |
| 2–5 | Experimental image sequences | `data/experiments/frames/`; panel provenance in `data/experiments/manifest.json` |
| 6 | Canonical transient segregation and fluid profiles | `analysis/figure06.py`; canonical case in `data/cases/` |
| 7 | Volume fraction and incline angle | `analysis/figure07.py`; five volume fractions, angles 20° and 30°, large diameter 600 µm |
| 8 | Particle size ratio | `analysis/figure08.py`; five volume fractions, large diameters 600 and 1200 µm, angle 30° |
| 9 | Species mean heights and vertical forces | `analysis/figure09.py`; canonical case time series and parameters |
| 10 | Baseline and mechanism controls | `analysis/figure10.py`; four `data/mechanisms/*_meanz.csv` files |

The canonical numerical case uses small diameter **200 µm**, large diameter
**600 µm**, volume fraction **0.10**, and incline angle **30°**. Archive names
encode **radii**: `rs100um_rl300um` corresponds to these diameters.

The parameter sweep includes volume fractions 0.05, 0.10, 0.20, 0.30, and 0.40.
The mechanism comparison uses its own recorded initialization; its baseline
curve is therefore a separate realization from the canonical case in Figure 6.
The equal-settling control in Figure 10 is the analytic Stokes control,
with upward force `F_add = (8/9) W'_L` on each large particle.

## What is included

- Final manuscript, reference figure assets, and the 18 experimental panels
  embedded in the final manuscript.
- Per-frame species mean-height curves, the force quantities used in Figure 9,
  horizontally averaged fluid profiles, and the particle snapshots used in
  Figures 6–8.
- Original inputs for the 15 sweep cases and four mechanism controls, including
  the required OpenFOAM `CFD/0/` initial fields.
- Figure scripts, original raw-output analysis tools, simulation preparation
  utilities, and the local solver modifications needed for diagnostic output.
- Source locations, upstream version pins, extraction metadata, and SHA-256
  checksums for the curated bundle.

The multi-terabyte CFD fields, full particle trajectories, restart binaries,
compiled software, temporary jobs, and unrelated parameter campaigns remain
in the source archive. Full laboratory movies and particle-tracking data are
not included in this bundle. The experimental panels are manuscript-resolution
image data; they do not constitute raw video measurements.

## Run new CFD–DEM simulations

See [simulation/README.md](simulation/README.md) for the solver stack, patches,
case preparation, and execution instructions. Replotting the included data
only requires Python; running new simulations requires Linux, MPI, OpenFOAM 6,
CFDEMcoupling, and LIGGGHTS.

## Data conventions and provenance

`z` points upward, normal to the incline. `Δz = mean(z_large) - mean(z_small)`:
negative values indicate RBNE and positive values indicate BNE. Stored lengths,
times, velocities, and forces use SI units; plots convert lengths to mm and
velocities to mm/s. Particle type 1 is small and type 2 is large.

See [docs/data_dictionary.md](docs/data_dictionary.md) and
[docs/provenance.md](docs/provenance.md) for file formats and source mapping.
The manuscript is included unchanged. No new CFD–DEM trajectories were
computed during repository preparation.

Please cite the paper when using the material; citation metadata are in
[CITATION.cff](CITATION.cff). For rights and upstream licensing, see
[RIGHTS.md](RIGHTS.md).
