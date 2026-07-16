# spectreasy 1.0.0

- Renamed the canonical `unmix_controls()` arguments `save_qc_plots` to
  `save_qc_png` and `refine` to `autospectral_refine`, and removed the former
  `unmix_threads` argument in favor of `n_threads`. The old names are no longer
  accepted by `unmix_controls()`.

- `n_threads` now covers event-wise AutoSpectral/Spectreasy AF assignment,
  NNLS, WLS, and RWLS fitting, including multi-AF OLS/NNLS candidate fitting.
  Dense fixed-reference OLS remains a vectorized matrix operation.

- Control HTML reports now save their complete NxN scatter matrix as a linked,
  high-resolution PNG instead of creating a standalone NxN HTML viewer. Points
  are larger and more opaque so sparse events remain visible when zooming.

- `unmix_controls()` now uses `gating_mode = "interactive"`, `"reuse"`, or
  `"automatic"` in place of `manual_gating`. Interactive mode preloads an
  existing gate CSV, reuse mode requires it, and automatic mode uses GMM gates.

- Report APIs now use `report_format` instead of `output_format`, with HTML as
  the default. Enumerated string arguments are matched case-insensitively
  across the public workflow, reporting, gating, and plotting interfaces.

## New

- Added the hosted Spectreasy cockpit with native project selection, schema-aligned control mapping, embedded analysis applets, and automatically persisted GUI settings.
- Added `unmix_controls()` as the main SCC workflow entry point.
- Added `unmix_samples()` support for loading a saved unmixing matrix file.
- Added structured control-file preflight validation and auto-generation support.

## Changes

- Changed the default `spectreasy_weight_quantile` from `0.9` to `0.65`.
- WLS now uses event-wise detector-error weights with detector-specific noise floors estimated from SCC low-signal tails by default.
- Increased the default WLS maximum detector weight ratio from 100 to 1600 based on cross-dataset benchmark results.
- `unmix_controls()` writes `scc_detector_noise.csv`; `unmix_samples(method = "WLS")` loads it beside `scc_reference_matrix.csv`, estimates from `scc_dir` when available, or falls back to a scalar noise floor of 125.
- Multi-AF WLS now uses the same event-wise detector-error weights for AF-band selection and coefficient fitting.
- `control.type` in the control mapping now controls bead/cell FSC/SSC gating, with filename guessing used only when the column is empty.
- `universal.negative` can now point to a specific negative FCS file for SCC subtraction.
- Removed legacy spectral spread matrix helpers and unused spread-score QC outputs.
- Added AF basis-band extraction via `af_n_bands`/`af_max_cells` in `build_reference_matrix()`.
- `unmix_controls()`, `unmix_samples()`, and `calc_residuals()` now default to `unmixing_method = "Spectreasy"` / `method = "Spectreasy"`.
- `af_n_bands` now defaults to a broad fixed bank of 100 AF signatures.
- `af_n_bands` is now an exact contract: requested bands are never silently
  discarded, and impossible requests fail rather than returning a smaller bank.
- Added optional deterministic `seed` support in SCC report/matrix/control workflows.
- Updated static matrix export in `unmix_controls()` to follow selected method (`OLS`, `WLS`, `NNLS` proxy).
- Updated SCC/Sample QC reports to exclude AF bands from spectra overlays and NPS pages.
- Removed the conclusions/recommendations page from sample QC report.
- Improved detector fallback ordering in spectra overlays to laser-first order (`UV`, `V`, `B`, `YG`, `R`).
- Split control-stage and sample-stage unmixing workflows.
- Improved SCC histogram gating and unmixing scatter matrix plotting.
- Switched CSV I/O to base R readers/writers where applicable.
- Moved frontend assets to `inst/gui` for package build compliance.

## Maintenance

- Removed deprecated wrapper scripts and duplicate backend sources.
- Updated package metadata for Bioconductor submission preparation.
