# spectreasy: Full Spectrum Flow Cytometry Quality Control

`spectreasy` is an R package for reviewing single-color controls, building spectral reference matrices, and unmixing experimental samples.

## Key Features

- **Automated Gating**: Isolate positive populations from beads or cells using Gaussian Mixture Models
- **Background Subtraction**: Automatically subtract internal negative populations to isolate pure fluorophore signatures
- **Pre-Unmix SCC Review**: Generate a PDF with per-control gating, histogram, and spectrum diagnostics before unmixing
- **SCC-Variance WLS Unmixing**: Weighted unmixing using detector noise measured from the single-color controls
- **SCC Diagnostics & Visualization**: Spectra, gating plots, and SCC unmixing scatter outputs for control-stage QC
- **Interactive GUI**: Web-based interface for manual matrix adjustment
- **Bioconductor-Native In-Memory Workflows**: `unmix_samples()` accepts `flowSet` and `SingleCellExperiment`, and can return either container

---

## Installation

For the released Bioconductor version, install with:

```r
BiocManager::install("spectreasy")
```

If you need the development version from GitHub, install with:

```r
remotes::install_github("pkheisig/spectreasy")
```

`remotes` is used here rather than `devtools` because only GitHub installation is needed.

---

# Example workflow

This walkthrough demonstrates the primary `spectreasy` workflow on the release-hosted example dataset. The example project contains:

- single-color controls in `scc/`
- one experimental sample in `sample/sample.fcs`

The user-facing workflow is:

1. download the example data into a project directory
2. run `unmix_controls()`
3. review and supplement the generated `fcs_mapping.csv`
4. confirm the control file in the console so `unmix_controls()` can finish
5. run `unmix_samples()`
6. generate QC reports for the controls and unmixed samples (optional)

## 1. Download the example data

`spectreasy_example_data()` downloads the example archive once, caches it under the R user cache, and can copy the extracted files into a project directory for local work.

```r
library(spectreasy)

project_dir <- file.path(tempdir(), "spectreasy_vignette_project")
if (dir.exists(project_dir)) {
  unlink(project_dir, recursive = TRUE, force = TRUE)
}

example_paths <- spectreasy_example_data(dest_dir = project_dir)

list.files(project_dir, recursive = TRUE)
#> [1] "sample/sample.fcs"               "scc/Alexa Fluor 700 (Beads).fcs"
#> [3] "scc/BUV395 (Beads).fcs"          "scc/BV510 (Beads).fcs"          
#> [5] "scc/FITC (Beads).fcs"            "scc/LIVE DEAD NIR (Cells).fcs"  
#> [7] "scc/PerCP-Cy5.5 (Beads).fcs"     "scc/Unstained (Cells).fcs"
```

For the remainder of this walkthrough, the commands are shown as they would be run from the project directory created above.

## 2. Start the control-stage workflow

Run `unmix_controls()` first. If `fcs_mapping.csv` is missing, `auto_create_control = TRUE` creates it automatically and then pauses for review.

```r
setwd(project_dir)

unmix_controls(
  scc_dir = "scc",
  auto_create_control = TRUE,
  cytometer = "Aurora",
  auto_unknown_fluor_policy = "by_channel",
  unmix_method = "WLS",
  unmix_scatter_panel_size_mm = 30,
  save_qc_plots = TRUE
)
```

After the control file is created, `unmix_controls()` prints a confirmation prompt and waits:

```text
Proceed with unmix_controls now? [y/n]:
```

## 3. Review and supplement `fcs_mapping.csv`

Open the generated `fcs_mapping.csv` in the project directory and complete the panel annotation before continuing. At minimum, review these columns:

- `fluorophore`
- `marker`
- `channel`
- `control.type`
- `is.viability`

For the example dataset, the reviewed control file looks like this:

|filename                    |fluorophore     |marker           |channel |control.type |universal.negative |large.gate |is.viability |
|:---------------------------|:---------------|:----------------|:-------|:------------|:------------------|:----------|:------------|
|Alexa Fluor 700 (Beads).fcs |Alexa Fluor 700 |CD3              |R4-A    |beads        |                   |           |             |
|BUV395 (Beads).fcs          |BUV395          |CD45RA           |UV2-A   |beads        |                   |           |             |
|BV510 (Beads).fcs           |BV510           |CD27             |V7-A    |beads        |                   |           |             |
|FITC (Beads).fcs            |FITC            |CD8              |B2-A    |beads        |                   |           |             |
|LIVE DEAD NIR (Cells).fcs   |LIVE DEAD NIR   |Live             |R7-A    |cells        |                   |           |TRUE         |
|PerCP-Cy5.5 (Beads).fcs     |PerCP-Cy5.5     |CCR7             |B9-A    |beads        |                   |           |             |
|Unstained (Cells).fcs       |AF              |Autofluorescence |UV7-A   |cells        |                   |           |             |

## 4. Return to the console and confirm with `y`

Once `fcs_mapping.csv` has been reviewed and saved, return to the console where `unmix_controls()` is waiting and enter:

```text
y
```

The same `unmix_controls()` call then continues and writes the control-stage outputs to `spectreasy_outputs/unmix_controls/`.

```
#>  [1] "fsc_ssc/Alexa Fluor 700 (Beads)_fsc_ssc.png"
#>  [2] "fsc_ssc/BUV395 (Beads)_fsc_ssc.png"
#>  [3] "fsc_ssc/BV510 (Beads)_fsc_ssc.png"
#>  [4] "fsc_ssc/FITC (Beads)_fsc_ssc.png"
#>  [5] "fsc_ssc/LIVE DEAD NIR (Cells)_fsc_ssc.png"
#>  [6] "fsc_ssc/PerCP-Cy5.5 (Beads)_fsc_ssc.png"
#>  [7] "intensity_scatter/Alexa Fluor 700 (Beads)_intensity_scatter.png"
#>  [8] "intensity_scatter/BUV395 (Beads)_intensity_scatter.png"
#>  [9] "intensity_scatter/BV510 (Beads)_intensity_scatter.png"
#> [10] "intensity_scatter/FITC (Beads)_intensity_scatter.png"
#> [11] "intensity_scatter/LIVE DEAD NIR (Cells)_intensity_scatter.png"
#> [12] "intensity_scatter/PerCP-Cy5.5 (Beads)_intensity_scatter.png"
#> [13] "scc_af_spectra.png"
#> [14] "scc_detector_noise.csv"
#> [15] "scc_reference_matrix.csv"
#> [16] "scc_spectra.png"
#> [17] "scc_unmixing_matrix.csv"
#> [18] "scc_unmixing_scatter_matrix.png"
#> [19] "scc_variances.csv"
#> [20] "spectrum/Alexa Fluor 700 (Beads)_spectrum.png"
#> [21] "spectrum/BUV395 (Beads)_spectrum.png"
#> [22] "spectrum/BV510 (Beads)_spectrum.png"
#> [23] "spectrum/FITC (Beads)_spectrum.png"
#> [24] "spectrum/LIVE DEAD NIR (Cells)_spectrum.png"
#> [25] "spectrum/PerCP-Cy5.5 (Beads)_spectrum.png"
#> [26] "unmixed_fcs/Alexa Fluor 700 (Beads)_unmixed.fcs"
#> [27] "unmixed_fcs/BUV395 (Beads)_unmixed.fcs"
#> [28] "unmixed_fcs/BV510 (Beads)_unmixed.fcs"
#> [29] "unmixed_fcs/FITC (Beads)_unmixed.fcs"
#> [30] "unmixed_fcs/LIVE DEAD NIR (Cells)_unmixed.fcs"
#> [31] "unmixed_fcs/PerCP-Cy5.5 (Beads)_unmixed.fcs"
#> [32] "unmixed_fcs/Unstained (Cells)_unmixed.fcs"
```

Key outputs from this step include:

- `fcs_mapping.csv`
- `spectreasy_outputs/unmix_controls/scc_detector_noise.csv`
- `spectreasy_outputs/unmix_controls/scc_reference_matrix.csv`
- `spectreasy_outputs/unmix_controls/scc_variances.csv`
- `spectreasy_outputs/unmix_controls/scc_spectra.png`
- `spectreasy_outputs/unmix_controls/scc_unmixing_matrix.csv`
- `spectreasy_outputs/unmix_controls/scc_unmixing_scatter_matrix.png`
- `spectreasy_outputs/unmix_controls/fsc_ssc/*.png`
- `spectreasy_outputs/unmix_controls/intensity_scatter/*.png`
- `spectreasy_outputs/unmix_controls/spectrum/*.png`
- `spectreasy_outputs/unmix_controls/unmixed_fcs/*.fcs`

The control-stage run also writes visual checks for each single-color control. For one color, the three plots below show the FSC/SSC gate, the intensity-vs-FSC scatter gate, and the detector spectrum used to build the reference matrix:

<p align="center">
  <img src="man/figures/vignette_fsc_ssc.png" width="48%" />
  <img src="man/figures/vignette_intensity_scatter.png" width="48%" />
</p>

<p align="center">
  <img src="man/figures/vignette_spectrum.png" width="100%" />
</p>

The same run creates the NxN scatter matrix for the single-color controls. Each row is one control, and each column checks how much signal appears in the other unmixed markers.

<p align="center">
  <img src="man/figures/vignette_scatter_matrix.png" width="100%" />
</p>

### What WLS Uses

When `method = "WLS"`, `spectreasy` uses an event-wise detector-error model: detectors with higher non-negative signal in an event get lower weight for that event. The detector noise floor is estimated from the low-signal tail of the SCC files and written to `scc_detector_noise.csv`; if no estimate is available, `spectreasy` falls back to a scalar floor of 125. The SCC population variances in `scc_variances.csv` are still written as reference QC metadata, but they are not used as default WLS detector weights.

The `control.type` column in `fcs_mapping.csv` also matters for this step. It tells `spectreasy` whether each control should be gated as `beads` or `cells`. If `control.type` is empty, `spectreasy` falls back to filename-based guessing.

## 5. Unmix the experimental sample

After the control-stage workflow has completed, unmix the experimental files with `unmix_samples()`. The reference matrix written by `unmix_controls()` is loaded by default.

```r
unmixed <- unmix_samples(
  sample_dir = "sample",
  output_dir = "spectreasy_outputs/unmix_samples"
)
```

For the example dataset, this writes:

- `spectreasy_outputs/unmix_samples/sample_unmixed.fcs`

and returns a named list with one element per sample.

| Alexa Fluor 700|      BUV395|     BV510|        FITC| LIVE DEAD NIR| PerCP-Cy5.5|         AF|File   |
|---------------:|-----------:|---------:|-----------:|-------------:|-----------:|----------:|:------|
|        69.30757|    92.68504|  815.3799|   363.95801|     146.16334|   -19.78563|     0.0000|sample |
|        30.20207|   669.77978|  -69.6220| -1712.11376|     156.65262|   695.78227| 11494.0203|sample |
|       -19.07283|  -261.43980|  616.0712|   182.37640|     -56.23331|    93.10464|     0.0000|sample |
|      -103.82067|    78.38529|  496.1857|   146.58977|      67.91341|   -61.55665|  -442.6262|sample |
|      8233.14537| 16461.14920| 1683.2967|  1009.26887|      32.07973|  -257.83746|  3400.9268|sample |
|       -42.85445|  -235.11592| -247.3907|   -55.22978|     218.27591|  -157.39362|   263.2533|sample |

## 6. Generate quality control reports (optional)

After unmixing, you can generate comprehensive PDF reports to inspect the quality of both the single-color controls and the unmixed experimental samples.

### Single-Color Control (SCC) Report

The SCC report reviews gating, peak channels, and signal distributions for each control file. By default, it writes the report to `"spectreasy_outputs/unmix_samples/qc_controls_report.pdf"`.

```r
qc_controls(
  scc_dir = "scc",
  cytometer = "Aurora",
  seed = 1
)
```

### Samples Report

The overall sample report visualizes unmixing quality across samples, including spectra overlays, detector residuals, spread matrices, and marker scatter plots. By default, it writes the report to `"spectreasy_outputs/unmix_samples/qc_samples_report.pdf"`.

```r
qc_samples(
  results = unmixed,
  M = ctrl$M
)
```

# Optional steps

The sections below are useful extensions, but they are not required for the core `unmix_controls()` -> `unmix_samples()` workflow shown above.

## Per-cell Autofluorescence (AF) Extraction

By default, `unmix_controls()` uses `af_n_bands = "auto"` to build a FlowSOM autofluorescence bank from pooled unstained/AF control events. If your cells have different AF shapes from cell to cell, the SOM bank represents those shapes as multiple candidate AF signatures.

Use the two multi-AF settings in the control-stage call. `af_n_bands` controls the size of the shared pooled AF bank, while `include_multi_af` tells `spectreasy` to include additional AF controls from the `af/` directory when those files are available.

```r
ctrl_multi_af <- unmix_controls(
  scc_dir = "scc",
  control_file = "fcs_mapping.csv",
  cytometer = "Aurora",
  output_dir = "spectreasy_outputs/unmix_controls_multi_af",
  unmix_method = "WLS",
  include_multi_af = TRUE,
  af_n_bands = "auto",
  af_auto_max_bands = 100,
  seed = 1
)
```

Then pass the saved control-stage matrix to `unmix_samples()`. `unmix_samples()` does not rebuild missing matrices from SCC files; if the matrix is absent, run `unmix_controls()` first.

```r
unmixed_multi_af <- unmix_samples(
  sample_dir = "sample",
  unmixing_matrix_file = ctrl_multi_af$reference_matrix_file,
  method = "WLS",
  output_dir = "spectreasy_outputs/unmix_samples_multi_af"
)
```

`af_n_bands` is like choosing how many AF "flavors" to model. With `af_n_bands = "auto"`, spectreasy builds a SOM bank with up to `af_auto_max_bands = 100` SOM nodes by default, then prepends a mean AF row. That gives 101 AF rows before any contaminant QC removals. During unmixing, each event chooses one AF profile with a joint covariance + residual score, so the larger AF bank can describe varied autofluorescence without using every AF profile in the final fit.

For direct control over the multi-AF bank size, set `af_auto_max_bands` or pass an explicit integer to `af_n_bands`. For difficult samples with structured AF left after the first pass, set `af_refine = TRUE` to append second-pass modulated AF spectra from high-error unstained cells. Keep it off unless benchmark metrics show it helps your panel.

## Use a reviewed control CSV in non-interactive workflows

For scripts, reports, or CI jobs, you can supply a pre-existing, reviewed control CSV file via `control_file` to `unmix_controls()` to skip the confirmation prompt.

```r
ctrl_noninteractive <- unmix_controls(
  scc_dir = file.path(project_dir, "scc"),
  control_file = control_file,
  cytometer = "Aurora",
  output_dir = file.path(project_dir, "spectreasy_outputs", "unmix_controls_noninteractive"),
  unmix_method = "WLS",
  seed = 1
)

dim(ctrl_noninteractive$M)
#> [1] 9 64
```

## Pass the in-memory reference matrix directly

You can pass the in-memory reference matrix returned by `unmix_controls()` directly to `unmix_samples()` instead of loading it from the saved CSV file.

For WLS, `unmix_samples()` will also load `scc_detector_noise.csv` beside the saved reference matrix when it is available. The optional `scc_variances.csv` remains useful as control-spread QC metadata.

```r
fluor_reference_matrix <- ctrl$M
marker_map <- stats::setNames(control_df$marker, control_df$fluorophore)
reference_matrix <- fluor_reference_matrix
mapped_names <- marker_map[rownames(reference_matrix)]
na_idx <- is.na(mapped_names)
mapped_names[na_idx] <- rownames(reference_matrix)[na_idx]
rownames(reference_matrix) <- unname(mapped_names)

unmixed_direct <- unmix_samples(
  sample_dir = file.path(project_dir, "sample"),
  M = reference_matrix,
  method = "WLS",
  output_dir = file.path(project_dir, "spectreasy_outputs", "unmix_samples_direct")
)

names(unmixed_direct)
#> [1] "sample"
```

## Inspect quick QC plots

Reference spectra and spectral spread plots should be interpreted in fluorophore space, so the original fluorophore-labeled control matrix is used here.

```r
reference_matrix_no_af <- fluor_reference_matrix[!grepl("^AF($|_)", rownames(fluor_reference_matrix), ignore.case = TRUE), , drop = FALSE]

plot_spectra(reference_matrix_no_af, output_file = NULL)
```

<p align="center">
  <img src="man/figures/vignette_spectra.png" width="80%" />
</p>

```r
plot_ssm(calculate_ssm(reference_matrix_no_af), output_file = NULL)
```

<p align="center">
  <img src="man/figures/vignette_ssm.png" width="80%" />
</p>

```r
plot_nps(calculate_nps(sample_results), output_file = NULL)
```

<p align="center">
  <img src="man/figures/vignette_nps.png" width="80%" />
</p>

---

**Author**: Paul Heisig  
**Email**: p.k.s.heisig@amsterdamumc.nl
