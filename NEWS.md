# spectreasy 0.99.0

## New

- Added `autounmix_controls()` as the main SCC workflow entry point.
- Added `unmix_samples()` support for loading a saved unmixing matrix file.
- Added structured control-file preflight validation and auto-generation support.

## Changes

- Split control-stage and sample-stage unmixing workflows.
- Improved SCC histogram gating and unmixing scatter matrix plotting.
- Switched CSV I/O to base R readers/writers where applicable.
- Moved frontend assets to `inst/gui` for package build compliance.

## Maintenance

- Removed deprecated wrapper scripts and duplicate backend sources.
- Updated package metadata for Bioconductor submission preparation.
