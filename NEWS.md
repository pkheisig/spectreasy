# spectreasy 0.99.0

## New

- Added `unmix_controls()` as the main SCC workflow entry point.
- Added `unmix_samples()` support for loading a saved unmixing matrix file.
- Added structured control-file preflight validation and auto-generation support.

## Changes

- WLS now uses SCC-derived detector variances only; it no longer guesses weights from the reference matrix or sample brightness.
- `unmix_samples(method = "WLS")` can recompute missing WLS variances from SCC files and save a regenerated `scc_variances.csv`.
- Multi-AF WLS now uses SCC-derived variances for AF-band selection and coefficient fitting.
- `control.type` in the control mapping now controls bead/cell FSC/SSC gating, with filename guessing used only when the column is empty.
- `universal.negative` can now point to a specific negative FCS file for SCC subtraction.
- Removed the unused `background_noise` argument from spectral spread matrix calculation.
- Added AF basis-band extraction via `af_n_bands`/`af_max_cells` in `build_reference_matrix()`.
- Added optional deterministic `seed` support in SCC report/matrix/control workflows.
- Updated static matrix export in `unmix_controls()` to follow selected method (`OLS`, `WLS`, `NNLS` proxy).
- Updated SCC/Sample QC reports to exclude AF bands from spectra overlays, SSM, and NPS pages.
- Removed the conclusions/recommendations page from sample QC report.
- Improved detector fallback ordering in spectra overlays to laser-first order (`UV`, `V`, `B`, `YG`, `R`).
- Split control-stage and sample-stage unmixing workflows.
- Improved SCC histogram gating and unmixing scatter matrix plotting.
- Switched CSV I/O to base R readers/writers where applicable.
- Moved frontend assets to `inst/gui` for package build compliance.

## Maintenance

- Removed deprecated wrapper scripts and duplicate backend sources.
- Updated package metadata for Bioconductor submission preparation.
