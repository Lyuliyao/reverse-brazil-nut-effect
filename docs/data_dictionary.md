# Data dictionary

All numerical data are stored in SI units unless the column name says otherwise.

## Sweep / canonical cases

`data/cases/<case>/time_series.tsv` is a tab-separated per-frame table.

| Column | Meaning | Unit |
| --- | --- | --- |
| `time` | Physical simulation time | s |
| `z_mean_small`, `z_mean_large` | Unweighted arithmetic mean particle-center height per species | m |
| `z_gap_large_minus_small` | Large minus small species mean height | m |
| `count_small`, `count_large` | Particle counts per species | dimensionless |
| `fz_contact_small`, `fz_contact_large` | Mean vertical particle-particle contact force per particle | N |
| `fz_wall_small`, `fz_wall_large` | Mean vertical wall-contact force per particle | N |
| `fz_fluid_small`, `fz_fluid_large` | Mean vertical fluid force per particle | N |

Force columns are needed only for the canonical case / Figure 9.
The fluid-force trace includes its hydrostatic pressure-gradient contribution.
Figure 9 subtracts the buoyancy contribution and normalizes by the vertical
submerged particle weight. A centered 1 s rolling average is used for its
smoothed curves, following the original figure script.

`fluid_profile_1.5.csv` and `fluid_profile_28.5.csv` contain `z_m`, `ux_m_s`,
`uy_m_s`, and `uz_m_s`: cell-center height and horizontal arithmetic mean
velocity over a 32 × 32 grid in each of the 16 vertical layers.

`snapshots/t<time>.npz` contains:

| Array | Shape / meaning |
| --- | --- |
| `xyz` | N × 3 particle-center coordinates in m |
| `types` | N particle types: 1 = small; 2 = large |
| `radii` | N radii in m |
| `box` | 3 × 2 lower / upper box coordinates in m |
| `time` | Scalar physical time in s |

Each snapshot is an extracted recorded frame, not an interpolated particle
configuration. Its source HDF5 member and frame index are recorded in
`provenance/archive_export.json`.

## Mechanism controls

`data/mechanisms/*_meanz.csv` contains the original per-frame records:
`time`, `zbar_small`, `zbar_large`, `dz`, `vz_small`, `vz_large`, `cc_small`,
`cc_large`, `n_small`, and `n_large`. Heights are in m; velocities in m/s;
`cc_*` is mean contact count. Missing quantities are written as `nan`.
`*_meanz_metrics.json` retains the original summary metrics and their definitions.
Only the mean-height columns are required for Figure 10.

`data/canonical_mean_heights.csv` is an independent archived recalculation
of the canonical mean heights from the packed particle records, used to check
the bundled canonical time series.

## Experimental panels and reference assets

`data/experiments/manifest.json` maps each image panel to its figure number,
panel letter, one-based PDF page, image XObject, and SHA-256 checksum.
`paper/figure_assets/manifest.json` identifies the numerical figure images
embedded in the final manuscript. Their pixels preserve the manuscript's
reference presentation; generated plots may differ slightly in font rendering
between Python / Matplotlib versions.
