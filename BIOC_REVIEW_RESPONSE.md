# Bioconductor review response checklist for `spectreasy`

This file maps each major reviewer comment to the implemented fix in the package source.

## DESCRIPTION
- ORCID in `Authors@R`
  - Implemented in `DESCRIPTION`.
- `nnls` import / optional dependency mismatch
  - `nnls` checks were removed.
  - Dynamic NNLS and WLS are implemented via internal Rcpp kernels in `src/unmix_kernels.cpp`.
  - Relevant files: `R/calc_residuals.R`, `R/unmixing_matrix.R`, `src/unmix_kernels.cpp`.
- `plumber` in Suggests but used unguarded
  - Added explicit `requireNamespace("plumber", quietly = TRUE)` check in `R/launch_gui.R`.
- Funding role
  - Not addressed here.

## Vignette
- Vignette should run real package code
  - `vignettes/getting_started.Rmd` now demonstrates actual package workflows.
  - It covers matrix derivation, `flowSet` workflows, `SingleCellExperiment` workflows, plotting, and QC report generation.

## R code / maintainability
- Very long functions
  - Refactored into helper-driven orchestration:
    - `create_control_file()` → `R/control_interface.R`
    - `build_reference_matrix()` → `R/build_reference_matrix.R`
    - `autounmix_controls()` → `R/autounmix_controls.R`
    - `validate_control_file_mapping()` → `R/control_validation.R`
    - `plot_unmixing_scatter_matrix()` → `R/plot_functions.R`
    - `generate_scc_report()` → `R/scc_check.R`
    - `unmix_samples()` → `R/unmixing.R`
    - `plot_detector_residuals()` → `R/diagnostic_plots.R`
- Cyclomatic complexity hotspots
  - Substantially reduced by extracting internal helpers in the files above.
- Use of `<<-`
  - Removed; no remaining `<<-` in `R/`.
- `launch_gui()` should use installed-package paths
  - Uses `system.file()` in `R/launch_gui.R`.
- Check for npm before dev-mode GUI launch
  - Implemented in `R/launch_gui.R`.
- Avoid saving files without consent
  - Plot/report helpers now save only when explicit output paths are supplied.
  - `generate_scc_report()` no longer retains intermediate QC PNGs by default.
- GUI host security
  - `launch_gui()` runs on `127.0.0.1`.
- `unmix_samples()` overwriting FCS outputs
  - Default is `write_fcs = FALSE`.
  - Safe collision-avoiding filenames are used when writing output.
- Heavy R loops in NNLS/WLS
  - Moved to Rcpp in `src/unmix_kernels.cpp`.
- Bioconductor-native classes
  - `unmix_samples()` accepts `flowSet` and `SingleCellExperiment`.
  - It can return `flowSet` or `SingleCellExperiment`.
  - Detector residuals are attached in `attr(x, "spectreasy_residuals")` for `flowSet`
    and in `altExp(x, "detector_residuals")` for `SingleCellExperiment`.

## Tests / coverage
- Coverage concerns
  - Added tests for:
    - SCC workflows
    - GUI helpers
    - AF exclusion
    - reports
    - `flowSet` return paths
    - `SingleCellExperiment` input/output
  - Current test suite passes with `testthat::test_local()`.

## Packaging / release hygiene
- GUI payload too large
  - Removed `inst/gui/node_modules` from the packaged source.
  - Bundled production assets remain in `inst/gui/dist`.
- Local build/check artifacts included
  - `.Rbuildignore` updated to exclude `.Rcheck`, tarballs, generated vignette outputs,
    local session artifacts, and compiled objects.
- Vignette packaging
  - Removed generated `vignettes/getting_started.html` from source tree.
- GUI OpenAPI warning
  - Removed `NULL` default in the Plumber endpoint in `inst/api/gui_api.R`.

## Notes
- The source tarball builds cleanly enough for code/testing checks when built from a clean source tree.
- Full vignette rendering was not run here.
