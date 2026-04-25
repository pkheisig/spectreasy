# spectreasy: Full Spectrum Flow Cytometry Quality Control

`spectreasy` is an R package for reviewing single-color controls, building spectral reference matrices, and unmixing experimental samples.

## Key Features

- **Automated Gating**: Isolate positive populations from beads or cells using Gaussian Mixture Models
- **Background Subtraction**: Automatically subtract internal negative populations to isolate pure fluorophore signatures
- **Pre-Unmix SCC Review**: Generate a PDF with per-control gating, histogram, and spectrum diagnostics before unmixing
- **Per-Cell WLS Unmixing**: High-accuracy unmixing using photon-counting variance weighting
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

## Example data

`spectreasy_example_data()` downloads the example FCS dataset from the GitHub
release asset matching the installed package version, caches it locally, and
unzips it automatically.

```r
paths <- spectreasy_example_data()
paths$scc_dir
paths$sample_dir

local_paths <- spectreasy_example_data(dest_dir = getwd())
local_paths$sample_dir
```

---

## Project Setup

Organize your data with this folder structure:

```
my_project/
├── scc/                          # Single-color control FCS files
│   ├── FITC_beads.fcs
│   ├── PE_beads.fcs
│   └── Unstained.fcs             # Autofluorescence control
├── samples/                      # Experimental FCS files
│   └── Sample1.fcs
└── fcs_mapping.csv          # Control file mapping (optional)
```

### Control File

<img width="1102" height="201" alt="Screenshot 2026-02-20 at 11 00 23" src="https://github.com/user-attachments/assets/a8c5e253-e1c8-4592-8ea2-e0f404932a54" />


The control file maps FCS filenames to fluorophores and markers. Generate it with `spectreasy`:

```r
create_control_file(
  input_folder = "scc",
  cytometer = "Aurora",
  output_file = "fcs_mapping.csv"
)
```

Or create manually with columns: `filename`, `fluorophore`, `marker`, `channel` (`universal.negative` is optional).

For auto-generation, `spectreasy` uses shipped dictionaries:
- `inst/extdata/fluorophore_dictionary.csv`
- `inst/extdata/marker_dictionary.csv`

Logic:
- detect `fluorophore` and `marker` from filename aliases first
- if no fluorophore is detected, fallback to peak-channel mapping (for example `YG1 -> PE`, `UV2 -> BUV395`)
- if a file is unlabeled (no marker/fluor match) but matches cytometer AF-channel behavior, it is auto-tagged as `AF`
- if filename contains viability cues (`live`, `dead`, `viability`), `is.viability` is auto-set to `TRUE`
- if no marker is detected, marker is left empty

---

## Workflow

### Path A: Recommended Quick Workflow

#### Step 1: Generate the SCC review report with `generate_scc_report()`:

```r
library(spectreasy)

generate_scc_report(
  scc_dir = "scc",
  control_file = "fcs_mapping.csv",
  cytometer = "Aurora",
  output_file = file.path("spectreasy_outputs", "SCC_QC_Report.pdf"),
  seed = 1
)
```

Use this report to inspect each SCC before unmixing. The PDF includes:
- FSC/SSC auto-gating per control
- peak-channel histogram gating per control
- per-control spectrum distributions
- global reference spectra overlay
- spectral spread matrix

AF handling in reports: AF basis rows (`AF`, `AF_2`, ...) are excluded from spectra overlays, spectral spread matrix pages, and NPS pages to keep summaries focused on non-AF biological channels.

By default, the report does not retain intermediate PNG QC assets. If you want to keep them, set `save_qc_pngs = TRUE` and choose `qc_plot_dir`.

#### Step 2: Run controls with `autounmix_controls()`:

```r
ctrl <- autounmix_controls(
  scc_dir = "scc",
  control_file = "fcs_mapping.csv",
  auto_create_control = TRUE,
  cytometer = "Aurora",
  auto_unknown_fluor_policy = "by_channel",
  output_dir = "spectreasy_outputs/autounmix_controls",
  exclude_af = FALSE,                 # set TRUE to ignore unstained/AF controls
  unmix_method = "WLS",
  af_n_bands = 5,                     # optional: derive multiple AF basis signatures
  build_qc_plots = TRUE,
  unmix_scatter_panel_size_mm = 30,
  seed = 1
)
```

If `fcs_mapping.csv` is missing and `auto_create_control = TRUE`, `autounmix_controls()` auto-generates a control file (filename, marker, fluorophore, and detected peak channel), then asks for confirmation before continuing.
`cytometer` is used for channel-aware fluorophore inference via spectreasy's shipped dictionaries and detector metadata.

For newly auto-created files:
- `control.type` is auto-detected from filename tokens (`beads`/`cells`) when possible; AF rows are always `cells`
- `is.viability` is auto-detected from filename cues (`live`/`dead`/`viability`) when possible
- `universal.negative` is left empty by default for all rows
- `autounmix_controls()` pauses and asks for `y/n` confirmation so you can review/edit the file first

`autounmix_controls()` writes:
- `scc_reference_matrix.csv`
- `scc_spectra.png` (reference spectra overlay)
- `scc_unmixing_matrix.png` and `scc_unmixing_matrix.csv`
- `scc_unmixing_scatter_matrix.png` (lower-triangle scatter matrix, one single-stain file per row, with x=0/y=0 guides)

Note: `scc_unmixing_matrix.csv` is now exported using the selected `unmix_method` (`OLS`, `WLS`, or `NNLS`).  
For `NNLS`, the exported static matrix is a deterministic linear proxy; exact NNLS behavior remains available via dynamic unmixing (`unmix_samples(..., M = ..., method = "NNLS")`).

Set `unmix_scatter_panel_size_mm` higher (for example `40`) if you want larger per-panel scatter plots.

`autounmix_controls()` also runs a strict preflight check before processing:
- every SCC file must be mapped in `fcs_mapping.csv`
- non-AF rows must define a valid `channel`
- if `universal.negative` is present, values for active SCC rows must be empty/keyword or reference a file present in your selected SCC/AF directories

If you want to keep an unstained/AF file in the SCC folder but leave it out of the control-stage matrix, set `exclude_af = TRUE`. This safely skips AF/unstained rows even when no AF row exists in `fcs_mapping.csv`.

---

#### Step 3: Launch the GUI if you need manual matrix adjustment.

```r
launch_gui(
  matrix_dir = "spectreasy_outputs/autounmix_controls",
  samples_dir = "samples"
)
```

This starts both the backend API and bundled frontend on one port (default `http://localhost:8000`) and opens it in your browser.

#### Step 4: Unmix samples using the unmixing matrix generated in the autounmix_controls step.

```r
# Uses saved unmixing matrix by filepath (default points to autounmix_controls output)
unmixed <- unmix_samples(
  sample_dir = "samples",
  unmixing_matrix_file = "spectreasy_outputs/autounmix_controls/scc_unmixing_matrix.csv",
  output_dir = "spectreasy_outputs/unmix_samples",
  write_fcs = TRUE
)
```

When writing unmixed FCS files, primary feature names come from matrix row names (fluorophores).  
Secondary feature names are taken from `fcs_mapping.csv` (`marker` column) when available.

For in-memory Bioconductor workflows, `unmix_samples()` also accepts a `flowSet` or
`SingleCellExperiment`, and can return `return_type = "flowSet"` or
`return_type = "SingleCellExperiment"`.

### Path B: Step-wise manual workflow

#### Step 1: Build Reference Matrix

Extract spectral signatures from single-color controls:

```r
library(spectreasy)

# Load control file (optional but recommended)
control_df <- read.csv("fcs_mapping.csv", stringsAsFactors = FALSE, check.names = FALSE)

# Build reference matrix from SCC files
M <- build_reference_matrix(
  input_folder = "scc",
  output_folder = "gating_plots",
  save_qc_plots = TRUE,
  control_df = control_df,
  default_sample_type = "beads",
  cytometer = "Aurora",
  seed = 1
)
```

This saves gating/spectrum plots to `gating_plots/` and returns the reference matrix in memory.

For per-cell AF extraction with multiple AF basis signatures, increase `af_n_bands` (for example `af_n_bands = 10`).

For reproducible SCC gating/subsampling and AF-band extraction, set `seed` in `generate_scc_report()`, `build_reference_matrix()`, or `autounmix_controls()`.

---

Play around with gating parameters if auto-gating fails:

```r
M <- build_reference_matrix(
  input_folder = "scc",
  output_folder = "gating_plots",
  save_qc_plots = TRUE,
  histogram_pct_beads = 0.98,         # Width of positive gate for beads
  histogram_pct_cells = 0.35,         # Width of positive gate for cells
  max_clusters = 6,                   # Maximum GMM clusters
  gate_contour_beads = 0.999999999,   # Contour level for bead gates
  gate_contour_cells = 0.95           # Contour level for cell gates
)
```

---

#### Step 2: Unmix Experimental Samples

Apply the reference matrix to your samples:

```r
# Option 1: dynamic unmixing directly from reference matrix (M)
unmixed <- unmix_samples(
  sample_dir = "samples",
  M = M,
  method = "WLS",                     # "OLS", "WLS", or "NNLS"
  cytometer = "Aurora",
  output_dir = "spectreasy_outputs/unmix_samples",
  write_fcs = TRUE
)

# Option 2: static unmixing from saved unmixing matrix (W)
unmixed_w <- unmix_samples(
  sample_dir = "samples",
  unmixing_matrix_file = "spectreasy_outputs/autounmix_controls/scc_unmixing_matrix.csv",
  output_dir = "spectreasy_outputs/unmix_samples_w",
  write_fcs = TRUE
)
```

**Methods:**
- **OLS**: Ordinary least squares — fast, suitable for most panels
- **WLS**: Weighted least squares — accounts for photon-counting noise, best accuracy
- **NNLS**: Non-negative least squares — forces positive abundances (dynamic mode). Static NNLS matrix export is a linear proxy.

---

### Optional: Interactive Matrix Adjustment (before sample unmixing)

For manual fine-tuning, use the web interface.

No terminal setup is required for end users. The bundled GUI is served directly by `launch_gui()`.

#### Launch the GUI:

```r
launch_gui(
  matrix_dir = getwd(),
  samples_dir = "samples"
)
```

Developer mode (optional, for GUI hacking):

```r
launch_gui(
  matrix_dir = getwd(),
  samples_dir = "samples",
  dev_mode = TRUE
)
```

#### What To Do After GUI

1. Save your adjusted matrix CSV (typically `scc_unmixing_matrix.csv` or `scc_reference_matrix.csv`).
2. Run `unmix_samples(...)` on your experimental samples:
   - if you edited `scc_unmixing_matrix.csv`, pass `unmixing_matrix_file = "..."`
   - if you edited a reference matrix, load it as `M` and pass `M = ...`

### Output Directories

- `spectreasy_outputs/autounmix_controls/scc_reference_matrix.csv`: Reference matrix written by `autounmix_controls(...)`
- `spectreasy_outputs/autounmix_controls/`: SCC control-stage outputs (`scc_reference_matrix.csv`, `scc_unmixing_matrix.csv/.png`, `scc_spectra.png`, `scc_unmixing_scatter_matrix.png`)
- `spectreasy_outputs/autounmix_controls/scc_unmixed/`: Unmixed SCC control files (FCS format)
- `spectreasy_outputs/unmix_samples/`: Unmixed experimental data (FCS format)

If you run the manual path with `build_reference_matrix(...)`, `output_folder` (for example `gating_plots/`) is used for build-stage QC plots.

---

## Report APIs

```r
# SCC review report (recommended before autounmix_controls)
generate_scc_report(
  scc_dir = "scc",
  control_file = "fcs_mapping.csv",
  output_file = file.path("spectreasy_outputs", "SCC_QC_Report.pdf")
)

# Full sample-level report
generate_qc_report(
      results_df = do.call(rbind, lapply(unmixed, `[[`, "data")),
      M = ctrl$M,  # matrix used for unmixing context
      output_file = file.path("spectreasy_outputs", "Sample_QC_Report.pdf"),
      sample_nxn_rows_per_page = 10,
      nxn_all_samples = FALSE
)
```

`generate_qc_report()` now adds per-sample NxN marker scatter pages after the
summary diagnostics. By default it shows 10 marker rows/columns per page block,
which keeps the panels square and standardized across pages. By default only the
first sample gets NxN pages; set `nxn_all_samples = TRUE` to include all samples.

---

## QC Plot Interpretation Guide

Use the same rule everywhere: first look for consistency across files/markers, then look for outlier structures that repeat in specific channels or populations.

- `Reference Spectra Overlay`: Good = each fluorophore shows a clear dominant detector profile with smooth shape. Bad = noisy, flattened, or unexpectedly broad/overlapping profiles. Action = recheck SCC gating, fluorophore labels, and control quality.
- `FSC/SSC Auto-Gate`: Good = the gate captures the main bead/cell cloud and excludes obvious debris/outliers. Bad = clipped main populations or large debris regions retained. Action = review sample type, scatter distributions, and gating parameters.
- `Peak-Channel Histogram Gate`: Good = the kept interval isolates the bright positive population. Bad = the interval covers mostly negative/background events or misses the dominant bright peak. Action = verify the control-file channel and the single-stain identity.
- `Per-Event Spectrum Distribution`: Good = the spectral shape is smooth and concentrated in the expected detector region. Bad = broad, noisy, or multimodal distributions that do not match the assigned fluorophore. Action = inspect labeling, staining quality, and instrument setup.
- `Unmixing Matrix Coefficients`: Good = strongest coefficients align with expected detector-marker relationships and off-target coefficients are moderate. Bad = many large off-target coefficients across rows. Action = inspect collinear markers, low-quality controls, or unstable matrix inversion.
- `Unmixing Scatter Matrix` (SCC controls): Good = row-stain events are high on Y while X (other markers) stays near zero. Bad = large off-axis/off-target clouds in lower-triangle panels. Action = inspect spillover-heavy pairs and verify single-stain identity.
- `Sample NxN Scatter Pages` (experimental samples): Good = each sample shows compact, interpretable pairwise structure without widespread off-axis spread. Bad = broad diagonal/off-axis clouds across many panels suggest unstable unmixing, high spread, or problematic marker combinations.
- `Spectral Spread Matrix`: Good = mostly low off-diagonal spread values. Bad = bright/high off-diagonal cells for specific marker pairs. Action = avoid those pairs for dim co-expression readouts or adjust panel design.
- `Residual Contributions`: Good = median residual per detector stays near zero. Bad = systematic positive/negative shifts in specific detectors. Action = investigate missing fluorophores, detector drift, or matrix mismatch in affected detector groups.
- `Negative Population Spread (NPS)`: Good = low and comparable MAD across files and markers. Bad = isolated high bars in selected markers/files. Action = identify bright spreaders into those channels and rebalance panel usage.

---

**Author**: Paul Heisig  
**Email**: p.k.s.heisig@amsterdamumc.nl
