#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(flowCore)
  library(pkgload)
})

repo_dir <- "/Users/pkheisig/Building/spectreasy"
data_dir <- "/Users/pkheisig/Library/CloudStorage/GoogleDrive-pkheisig@gmail.com/My Drive/PhD/Project_Erlangen/Erlangen_Data/PBMC_test_culture_titrations"
out_root <- file.path(data_dir, "spectreasy_outputs", "wls_comparison_20260605_111902")
subsample_n <- 5000L

pkgload::load_all(repo_dir, quiet = TRUE)

read_reference_matrix <- function(path) {
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  marker <- df[[1]]
  mat <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- marker
  mat
}

auto_fixed_matrix_for_sample <- function(M, weights) {
  noise_floor <- 1 / weights
  noise_floor[!is.finite(noise_floor) | noise_floor <= 0] <- 1
  M_auto <- M
  attr(M_auto, "detector_noise") <- data.frame(
    detector = colnames(M),
    noise_floor = as.numeric(noise_floor),
    signal_scale = rep(0, ncol(M)),
    stringsAsFactors = FALSE
  )
  M_auto
}

reference_file <- file.path(out_root, "reference", "scc_reference_matrix.csv")
detector_noise_file <- file.path(out_root, "reference", "scc_detector_noise.csv")
weights_file <- file.path(out_root, "tables", "autospectral_fixed_weights_by_sample.csv")

M_base <- read_reference_matrix(reference_file)
detector_noise <- utils::read.csv(detector_noise_file, stringsAsFactors = FALSE, check.names = FALSE)
M_eventwise <- spectreasy:::.attach_detector_noise(M_base, detector_noise, source = detector_noise_file)
weights <- utils::read.csv(weights_file, stringsAsFactors = FALSE, check.names = FALSE)

sample_files <- sort(list.files(file.path(data_dir, "samples"), pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE))
results_eventwise <- list()
results_auto_fixed <- list()

set.seed(20260605)
for (sample_file in sample_files) {
  sample_name <- tools::file_path_sans_ext(basename(sample_file))
  message("Subsampled report residuals: ", sample_name)
  ff <- flowCore::read.FCS(sample_file, transformation = FALSE, truncate_max_range = FALSE)
  expr <- flowCore::exprs(ff)
  keep <- if (nrow(expr) > subsample_n) sort(sample(seq_len(nrow(expr)), subsample_n)) else seq_len(nrow(expr))
  ff_sub <- flowCore::flowFrame(expr[keep, , drop = FALSE])

  results_eventwise[[sample_name]] <- spectreasy::calc_residuals(
    ff_sub,
    M_eventwise,
    method = "WLS",
    file_name = sample_name,
    return_residuals = TRUE
  )

  sample_weights <- weights[weights$sample == sample_name, , drop = FALSE]
  sample_weights <- sample_weights[match(colnames(M_base), sample_weights$detector), , drop = FALSE]
  if (any(is.na(sample_weights$autospectral_fixed_weight))) {
    stop("Missing AutoSpectral-style weights for sample: ", sample_name)
  }
  M_auto <- auto_fixed_matrix_for_sample(M_base, sample_weights$autospectral_fixed_weight)
  results_auto_fixed[[sample_name]] <- spectreasy::calc_residuals(
    ff_sub,
    M_auto,
    method = "WLS",
    file_name = sample_name,
    return_residuals = TRUE,
    wls_max_weight_ratio = 1e12
  )
}

attr(results_eventwise, "method") <- "WLS"
attr(results_auto_fixed, "method") <- "WLS"
class(results_eventwise) <- c("spectreasy_unmixed_results", "list")
class(results_auto_fixed) <- c("spectreasy_unmixed_results", "list")

spectreasy::qc_samples(
  results = results_eventwise,
  M = M_eventwise,
  output_file = file.path(out_root, "spectreasy_eventwise_wls", "qc_samples_report_eventwise_wls_subsampled.pdf"),
  method = "WLS",
  sample_nxn_max_points = 2000,
  sample_nxn_rows_per_page = 6,
  nxn_all_samples = FALSE
)

spectreasy::qc_samples(
  results = results_auto_fixed,
  M = M_base,
  output_file = file.path(out_root, "autospectral_fixed_mean_wls", "qc_samples_report_autospectral_fixed_wls_subsampled.pdf"),
  method = "WLS",
  sample_nxn_max_points = 2000,
  sample_nxn_rows_per_page = 6,
  nxn_all_samples = FALSE
)

paired <- utils::read.csv(file.path(out_root, "tables", "paired_sample_metric_deltas.csv"), stringsAsFactors = FALSE)
summary_lines <- c(
  "# WLS comparison summary",
  "",
  "Methods compared:",
  "- `spectreasy_eventwise_wls`: SCC-derived detector noise floor plus event signal; weights recalculated per event.",
  "- `autospectral_fixed_mean_wls`: AutoSpectral 0.8.6-style fixed detector weights; weights = `1 / (sample channel mean + 1e-6)`.",
  "",
  "Important comparison detail:",
  "Both methods used the same spectreasy reference matrix and the same multi-AF candidate selection. This isolates the WLS weighting model difference.",
  "",
  "Outputs:",
  "- Full unmixed FCS files were written for both methods.",
  "- Full-dataset residual metrics are in `tables/`.",
  "- QC PDFs are subsampled to 5,000 events per sample to avoid R's 16 GB vector limit during PDF generation.",
  "",
  "Paired RMSE/RSS deltas:",
  "Positive values mean spectreasy event-wise WLS had lower residuals than AutoSpectral-style fixed WLS.",
  "",
  paste(capture.output(print(paired[, c("sample", "rmse_pct_better_spectreasy", "rss_pct_better_spectreasy")], row.names = FALSE)), collapse = "\n"),
  "",
  paste0("Mean RMSE improvement: ", round(mean(paired$rmse_pct_better_spectreasy), 3), "%"),
  paste0("Median RMSE improvement: ", round(stats::median(paired$rmse_pct_better_spectreasy), 3), "%")
)
writeLines(summary_lines, file.path(out_root, "README_wls_comparison.md"))
writeLines(
  c(
    "The initial all-event QC report attempt hit R's 16 GB vector memory limit.",
    paste0("These replacement QC reports use ", subsample_n, " events per sample."),
    "Full-dataset metrics and full unmixed FCS outputs were already generated before the report memory limit."
  ),
  file.path(out_root, "report_generation_notes.txt")
)

message("Done: ", out_root)
