#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(flowCore)
  library(ggplot2)
  library(gridExtra)
  library(pkgload)
})

repo_dir <- "/Users/pkheisig/Building/spectreasy"
data_dir <- "/Users/pkheisig/Library/CloudStorage/GoogleDrive-pkheisig@gmail.com/My Drive/PhD/Project_Erlangen/Erlangen_Data/PBMC_test_culture_titrations"

pkgload::load_all(repo_dir, quiet = TRUE)

run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_root <- file.path(data_dir, "spectreasy_outputs", paste0("wls_comparison_", run_id))
dirs <- list(
  root = out_root,
  reference = file.path(out_root, "reference"),
  spectreasy = file.path(out_root, "spectreasy_eventwise_wls"),
  autospectral = file.path(out_root, "autospectral_fixed_mean_wls"),
  plots = file.path(out_root, "comparison_plots"),
  tables = file.path(out_root, "tables")
)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(dirs$spectreasy, "unmixed_fcs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(dirs$autospectral, "unmixed_fcs"), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_root, "run_log.txt")
log_msg <- function(...) {
  msg <- paste0(...)
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg, "\n")
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg, "\n", file = log_file, append = TRUE)
}

read_reference_matrix <- function(path) {
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  marker <- df[[1]]
  mat <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- marker
  mat
}

write_reference_matrix <- function(M, path) {
  df <- as.data.frame(M, check.names = FALSE)
  df <- cbind(Marker = rownames(M), df)
  utils::write.csv(df, path, row.names = FALSE, quote = TRUE)
}

copy_if_exists <- function(src, dest) {
  if (file.exists(src)) file.copy(src, dest, overwrite = FALSE)
}

reference_file <- file.path(data_dir, "spectreasy_outputs", "unmix_controls", "scc_reference_matrix.csv")
detector_noise_file <- file.path(data_dir, "spectreasy_outputs", "unmix_controls", "scc_detector_noise.csv")

if (!file.exists(reference_file)) stop("Reference matrix not found: ", reference_file)
if (!file.exists(detector_noise_file)) stop("Detector noise file not found: ", detector_noise_file)

copy_if_exists(reference_file, file.path(dirs$reference, "scc_reference_matrix.csv"))
copy_if_exists(detector_noise_file, file.path(dirs$reference, "scc_detector_noise.csv"))

M_base <- read_reference_matrix(reference_file)
detector_noise <- utils::read.csv(detector_noise_file, stringsAsFactors = FALSE, check.names = FALSE)
M_eventwise <- spectreasy:::.attach_detector_noise(M_base, detector_noise, source = detector_noise_file)

sample_files <- list.files(file.path(data_dir, "samples"), pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
if (length(sample_files) == 0) stop("No sample FCS files found.")
sample_files <- sort(sample_files)

log_msg("Output root: ", out_root)
log_msg("Samples: ", length(sample_files))
log_msg("Reference rows: ", nrow(M_base), "; detectors: ", ncol(M_base))
log_msg("AutoSpectral installed: ", requireNamespace("AutoSpectral", quietly = TRUE),
        if (requireNamespace("AutoSpectral", quietly = TRUE)) paste0(" v", utils::packageVersion("AutoSpectral")) else "")

sample_metrics <- list()
detector_metrics <- list()
weight_rows <- list()
event_metric_rows <- list()
results_eventwise <- list()
results_auto_fixed <- list()

marker_columns <- rownames(M_base)
metadata_columns <- function(x) {
  intersect(colnames(x), c("Time", grep("^FSC|^SSC", colnames(x), value = TRUE), "File"))
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

summarise_result <- function(sample_name, method_label, Y, residuals) {
  event_rss <- rowSums(residuals^2)
  event_rmse <- sqrt(rowMeans(residuals^2))
  raw_event_rms <- sqrt(rowMeans(Y^2))
  rel_event_rmse <- event_rmse / pmax(raw_event_rms, 1e-6)
  data.frame(
    sample = sample_name,
    method = method_label,
    events = nrow(Y),
    detectors = ncol(Y),
    rss = sum(event_rss),
    mse = mean(residuals^2),
    rmse = sqrt(mean(residuals^2)),
    median_event_rmse = stats::median(event_rmse),
    q95_event_rmse = as.numeric(stats::quantile(event_rmse, 0.95, names = FALSE)),
    median_relative_event_rmse = stats::median(rel_event_rmse),
    q95_relative_event_rmse = as.numeric(stats::quantile(rel_event_rmse, 0.95, names = FALSE)),
    stringsAsFactors = FALSE
  )
}

summarise_detectors <- function(sample_name, method_label, residuals) {
  data.frame(
    sample = sample_name,
    method = method_label,
    detector = colnames(residuals),
    detector_rmse = sqrt(colMeans(residuals^2)),
    detector_median_abs_residual = apply(abs(residuals), 2, stats::median),
    stringsAsFactors = FALSE
  )
}

write_unmixed_fcs <- function(res_obj, output_dir, sample_name) {
  data <- res_obj$data
  cols_to_write <- setdiff(colnames(data), "File")
  exprs_mat <- as.matrix(data[, cols_to_write, drop = FALSE])
  storage.mode(exprs_mat) <- "numeric"
  ff <- flowCore::flowFrame(exprs_mat)
  out_file <- file.path(output_dir, paste0(sample_name, "_unmixed.fcs"))
  flowCore::write.FCS(ff, filename = out_file)
  out_file
}

for (sample_file in sample_files) {
  sample_name <- tools::file_path_sans_ext(basename(sample_file))
  log_msg("Processing ", sample_name)
  ff <- flowCore::read.FCS(sample_file, transformation = FALSE, truncate_max_range = FALSE)
  expr <- flowCore::exprs(ff)
  missing_detectors <- setdiff(colnames(M_base), colnames(expr))
  if (length(missing_detectors) > 0) {
    stop("Sample ", sample_name, " is missing detectors: ", paste(missing_detectors, collapse = ", "))
  }
  Y <- expr[, colnames(M_base), drop = FALSE]

  res_eventwise <- spectreasy::calc_residuals(
    ff,
    M_eventwise,
    method = "WLS",
    file_name = sample_name,
    return_residuals = TRUE
  )

  auto_weights <- 1 / (colMeans(Y) + 1e-6)
  auto_weights[!is.finite(auto_weights) | auto_weights <= 0] <- NA_real_
  if (anyNA(auto_weights)) {
    fallback <- stats::median(auto_weights[is.finite(auto_weights) & auto_weights > 0], na.rm = TRUE)
    if (!is.finite(fallback) || fallback <= 0) fallback <- 1
    auto_weights[!is.finite(auto_weights) | auto_weights <= 0] <- fallback
  }
  M_auto_fixed <- auto_fixed_matrix_for_sample(M_base, auto_weights)
  res_auto_fixed <- spectreasy::calc_residuals(
    ff,
    M_auto_fixed,
    method = "WLS",
    file_name = sample_name,
    return_residuals = TRUE,
    wls_max_weight_ratio = 1e12
  )

  results_eventwise[[sample_name]] <- res_eventwise
  results_auto_fixed[[sample_name]] <- res_auto_fixed

  write_unmixed_fcs(res_eventwise, file.path(dirs$spectreasy, "unmixed_fcs"), sample_name)
  write_unmixed_fcs(res_auto_fixed, file.path(dirs$autospectral, "unmixed_fcs"), sample_name)

  sample_metrics[[paste0(sample_name, "_eventwise")]] <- summarise_result(
    sample_name, "spectreasy_eventwise_wls", Y, res_eventwise$residuals
  )
  sample_metrics[[paste0(sample_name, "_auto")]] <- summarise_result(
    sample_name, "autospectral_fixed_mean_wls", Y, res_auto_fixed$residuals
  )
  detector_metrics[[paste0(sample_name, "_eventwise")]] <- summarise_detectors(
    sample_name, "spectreasy_eventwise_wls", res_eventwise$residuals
  )
  detector_metrics[[paste0(sample_name, "_auto")]] <- summarise_detectors(
    sample_name, "autospectral_fixed_mean_wls", res_auto_fixed$residuals
  )
  weight_rows[[sample_name]] <- data.frame(
    sample = sample_name,
    detector = names(auto_weights),
    autospectral_fixed_weight = as.numeric(auto_weights),
    stringsAsFactors = FALSE
  )

  set.seed(1001)
  keep <- if (nrow(Y) > 5000) sample(seq_len(nrow(Y)), 5000) else seq_len(nrow(Y))
  event_metric_rows[[paste0(sample_name, "_eventwise")]] <- data.frame(
    sample = sample_name,
    method = "spectreasy_eventwise_wls",
    event_rmse = sqrt(rowMeans(res_eventwise$residuals[keep, , drop = FALSE]^2)),
    stringsAsFactors = FALSE
  )
  event_metric_rows[[paste0(sample_name, "_auto")]] <- data.frame(
    sample = sample_name,
    method = "autospectral_fixed_mean_wls",
    event_rmse = sqrt(rowMeans(res_auto_fixed$residuals[keep, , drop = FALSE]^2)),
    stringsAsFactors = FALSE
  )

  rm(ff, expr, Y, res_eventwise, res_auto_fixed, M_auto_fixed)
  gc()
}

attr(results_eventwise, "method") <- "WLS"
attr(results_auto_fixed, "method") <- "WLS"
class(results_eventwise) <- c("spectreasy_unmixed_results", "list")
class(results_auto_fixed) <- c("spectreasy_unmixed_results", "list")

metrics <- do.call(rbind, sample_metrics)
detectors <- do.call(rbind, detector_metrics)
weights <- do.call(rbind, weight_rows)
event_metrics <- do.call(rbind, event_metric_rows)

wide_metrics <- reshape(
  metrics[, c("sample", "method", "rss", "rmse", "median_event_rmse", "q95_event_rmse", "median_relative_event_rmse")],
  idvar = "sample",
  timevar = "method",
  direction = "wide"
)
wide_metrics$rmse_delta_auto_minus_spectreasy <- wide_metrics$rmse.autospectral_fixed_mean_wls - wide_metrics$rmse.spectreasy_eventwise_wls
wide_metrics$rmse_pct_better_spectreasy <- 100 * wide_metrics$rmse_delta_auto_minus_spectreasy / wide_metrics$rmse.autospectral_fixed_mean_wls
wide_metrics$rss_delta_auto_minus_spectreasy <- wide_metrics$rss.autospectral_fixed_mean_wls - wide_metrics$rss.spectreasy_eventwise_wls
wide_metrics$rss_pct_better_spectreasy <- 100 * wide_metrics$rss_delta_auto_minus_spectreasy / wide_metrics$rss.autospectral_fixed_mean_wls

utils::write.csv(metrics, file.path(dirs$tables, "sample_residual_metrics.csv"), row.names = FALSE)
utils::write.csv(wide_metrics, file.path(dirs$tables, "paired_sample_metric_deltas.csv"), row.names = FALSE)
utils::write.csv(detectors, file.path(dirs$tables, "detector_residual_metrics.csv"), row.names = FALSE)
utils::write.csv(weights, file.path(dirs$tables, "autospectral_fixed_weights_by_sample.csv"), row.names = FALSE)
write_reference_matrix(M_base, file.path(dirs$reference, "reference_matrix_used.csv"))

method_labels <- c(
  spectreasy_eventwise_wls = "spectreasy event-wise WLS",
  autospectral_fixed_mean_wls = "AutoSpectral fixed mean WLS"
)

p_rmse <- ggplot(metrics, aes(x = reorder(sample, rmse), y = rmse, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = c(spectreasy_eventwise_wls = "#2563eb", autospectral_fixed_mean_wls = "#dc2626"), labels = method_labels) +
  coord_flip() +
  labs(title = "Sample RMSE: lower is better", x = NULL, y = "RMSE across detector residuals", fill = NULL) +
  theme_minimal(base_size = 11)

p_delta <- ggplot(wide_metrics, aes(x = reorder(sample, rmse_pct_better_spectreasy), y = rmse_pct_better_spectreasy)) +
  geom_hline(yintercept = 0, color = "grey45") +
  geom_col(fill = "#2563eb", width = 0.75) +
  coord_flip() +
  labs(
    title = "RMSE improvement of spectreasy event-wise WLS over AutoSpectral-style WLS",
    subtitle = "Positive means spectreasy has lower residual RMSE",
    x = NULL,
    y = "% lower RMSE"
  ) +
  theme_minimal(base_size = 11)

p_event <- ggplot(event_metrics, aes(x = method, y = event_rmse, fill = method)) +
  geom_boxplot(outlier.alpha = 0.08, width = 0.7) +
  scale_y_continuous(trans = "log10") +
  scale_fill_manual(values = c(spectreasy_eventwise_wls = "#2563eb", autospectral_fixed_mean_wls = "#dc2626"), labels = method_labels) +
  facet_wrap(~ sample, scales = "free_y") +
  labs(title = "Per-event residual RMSE distribution", x = NULL, y = "Event RMSE, log10 scale", fill = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

top_detectors <- detectors
top_detectors$sample_method <- paste(top_detectors$sample, top_detectors$method)
p_detector <- ggplot(top_detectors, aes(x = detector, y = sample, fill = detector_rmse)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "sqrt") +
  facet_wrap(~ method, ncol = 1) +
  labs(title = "Detector residual RMSE heatmap", x = "Detector", y = NULL, fill = "RMSE") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(file.path(dirs$plots, "sample_rmse_barplot.png"), p_rmse, width = 10, height = 7, dpi = 220)
ggsave(file.path(dirs$plots, "spectreasy_rmse_improvement.png"), p_delta, width = 10, height = 7, dpi = 220)
ggsave(file.path(dirs$plots, "event_rmse_distribution.png"), p_event, width = 12, height = 9, dpi = 220)
ggsave(file.path(dirs$plots, "detector_rmse_heatmap.png"), p_detector, width = 13, height = 9, dpi = 220)

pdf(file.path(dirs$root, "wls_comparison_residual_overview.pdf"), width = 11, height = 8.5)
gridExtra::grid.arrange(p_rmse, p_delta, nrow = 1)
gridExtra::grid.arrange(p_event, nrow = 1)
gridExtra::grid.arrange(p_detector, nrow = 1)
dev.off()

log_msg("Generating spectreasy event-wise WLS QC report")
spectreasy::qc_samples(
  results = results_eventwise,
  M = M_eventwise,
  output_file = file.path(dirs$spectreasy, "qc_samples_report_eventwise_wls.pdf"),
  method = "WLS",
  nxn_all_samples = FALSE
)

log_msg("Generating AutoSpectral-style fixed WLS QC report")
spectreasy::qc_samples(
  results = results_auto_fixed,
  M = M_base,
  output_file = file.path(dirs$autospectral, "qc_samples_report_autospectral_fixed_wls.pdf"),
  method = "WLS",
  nxn_all_samples = FALSE
)

summary_lines <- c(
  "# WLS comparison summary",
  "",
  paste0("Run ID: ", run_id),
  paste0("Data directory: ", data_dir),
  paste0("Output directory: ", out_root),
  "",
  "Methods compared:",
  "- spectreasy_eventwise_wls: SCC-derived detector noise floor plus event signal, weights recalculated per event.",
  "- autospectral_fixed_mean_wls: AutoSpectral 0.8.6-style fixed detector weights, weights = 1 / (sample channel mean + 1e-6).",
  "",
  "Important note:",
  "Both methods used the same spectreasy reference matrix and the same multi-AF candidate selection. This isolates the weighting model difference.",
  "",
  "Top-level paired deltas:",
  paste(capture.output(print(wide_metrics[, c("sample", "rmse_pct_better_spectreasy", "rss_pct_better_spectreasy")], row.names = FALSE)), collapse = "\n")
)
writeLines(summary_lines, file.path(dirs$root, "README_wls_comparison.md"))

log_msg("Done.")
cat(out_root, "\n")
