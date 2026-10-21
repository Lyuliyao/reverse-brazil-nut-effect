# Source provenance and scope

All data in this repository were exported from the authors' original CFD–DEM
research archive. No simulations
were re-run to produce this repository; absolute paths that appear in the JSON
records refer to locations in that archive.

## Included sources

- `paper/reverse_Brazil_nut_effect.pdf`: the manuscript PDF.
- Experimental panels and numerical reference images: extracted from that
  PDF's image objects; manifests identify figure, panel, PDF page, and checksum.
- Figure 1 schematic: `paper/pic/diagram.pdf` in the source project.
- Figures 6–9 code: `CFDEM_PROJECT/brazil_nut/results/paper_figs/`.
- Figure 10 code: adapted from the authors' development repository
  (commit `408f94f91d1c0bc8b8f8fa267b9d3b6507cc977a`) to bundled paths.
- Figure 10 CSVs: `CFDEM_PROJECT/mechanism_campaign_20260611/analysis/results/`,
  including the analytic `equal_settling_stokes` control.
- Sweep / canonical inputs, time series, and selected frames: the original
  `.tar` case archives under `CFDEM_PROJECT/` or
  `CFDEM_PROJECT/brazil_nut/cases/`.
- Mechanism inputs: original case archives under
  `CFDEM_PROJECT/mechanism_campaign_20260611/cases/`.
- Solver modifications: source checkouts under `CFDEM_INSTALL/`, represented
  as patches against recorded upstream commits.

`provenance/source_files.json` records source script / data hashes.
`provenance/archive_export.json` records archive names and sizes, extracted
members, selected frame indices, and curve derivations. Fields not needed
by the paper are omitted from the noncanonical sweep tables, preserving the
numeric strings in the retained columns. Tables whose archived analysis had
stopped before 40 s are recalculated from every recorded particle frame,
using the original figure reader's arithmetic species means. No heights or
particle configurations are reconstructed from the figure images.

Two archived cases had empty CFD input folders. Their common `CFD/0/`,
`constant/`, and `system/` templates were restored from the canonical case,
following the original case generator. All six initial fields in the other
17 cases are byte-identical; the physical CFD dictionaries match the other
sweep cases at the same angle. Fresh-run decomposition settings are selected
by `prepare_case.py`.
`provenance/initial_field_restoration.json` records this restoration and hashes.
All remaining case inputs retain their archived values; use the preparation
utility to make a fresh runnable copy.

## Validation and retained archive

`provenance/checksums.json` covers the research artifacts under `data/`,
`cases/`, `paper/` and `simulation/patches/`. Documentation and code are
versioned by git and are not checksummed; after an intentional change to a
checksummed file, run `python scripts/update_checksums.py` (CI runs it with
`--check`).

The bundle validator checks checksums, case / species identities, time coverage,
height definitions, snapshot means, canonical cross-check data, fluid-profile
maxima, and the Stokes-control force. Figure scripts generate PDFs and PNGs;
reference assets retain the final manuscript's presentation.

The raw multi-terabyte archive is not a Git dependency. It retains full CFD
fields, complete particle records and contact dumps, restart binaries, compiled
software, and campaigns outside the paper's scope. Extraction tools are kept
under `scripts/` for use where those source archives are available.

Full experimental videos are outside this curated bundle. The included
experimental images are the precise panels embedded in the final manuscript,
and their panel mapping is the record of what is shared here.
