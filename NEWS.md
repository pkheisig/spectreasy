# spectreasy 0.99.0

## New

- Added `unmix_controls()` as the main SCC workflow entry point.
- Added `unmix_samples()` support for loading a saved unmixing matrix file.
- Added structured control-file preflight validation and auto-generation support.
- Added `build_spectral_panel()` as a browser-based spectral panel builder for
  Aurora, Discover, ID7000, and Xenith theoretical spectra.
- Renamed the matrix-adjustment GUI entry point to `adjust_matrix()`.

## Changes

- WLS now uses event-wise detector-error weights with detector-specific noise floors estimated from SCC low-signal tails by default. SCC variances are retained as reference QC metadata.
- Increased the default WLS maximum detector weight ratio from 100 to 1600 based on cross-dataset benchmark results.
- `unmix_controls()` writes `scc_detector_noise.csv`; `unmix_samples(method = "WLS")` loads it beside `scc_reference_matrix.csv`, estimates from `scc_dir` when available, or falls back to a scalar noise floor of 125.
- Multi-AF event assignment now uses a joint covariance + residual score to choose the AF profile per event. WLS/RWLS are used for the final coefficient fit after AF assignment.
- Multi-AF event-wise WLS/RWLS can now run through `RcppParallel` when `multithreading = TRUE`. `n_threads = "auto"` uses `RcppParallel::defaultNumThreads()`, while oversized integer requests are clipped to the available thread count.
- `unmix_samples(verbose = TRUE)` now shows an interactive console progress bar across samples when available.
- `control.type` in the control mapping now controls bead/cell FSC/SSC gating, with filename guessing used only when the column is empty.
- Removed the unused `background_noise` argument from spectral spread matrix calculation.
- Added AF basis-band extraction via `af_n_bands`/`af_max_cells` in `build_reference_matrix()`.
- Added FlowSOM-based multi-AF bank extraction via `af_n_bands`/`af_auto_max_bands`; the default auto bank uses up to 100 SOM nodes plus a prepended mean AF row.
- `af_n_bands = "auto"` is now the default in reference-matrix workflows, and auto-built AF banks keep the SOM centers for covariance-based event assignment.
- Added optional `af_refine = TRUE` second-pass AF-library refinement, which appends modulated AF spectra from high-error unstained cells for benchmarking difficult AF cases.
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
