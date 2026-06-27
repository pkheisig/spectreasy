#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(flowCore)
  library(spectreasy)
})

benchmark_dir <- "/Users/pkheisig/Building/Project_Spectreasy/benchmarking"

source(file.path(benchmark_dir, "run_wls_benchmark.R"), chdir = TRUE)

out_dir <- file.path(benchmark_dir, "runs", "BLIND-AF-SPIKE")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

max_train_per_sample <- as.integer(Sys.getenv("BLIND_AF_TRAIN_EVENTS", "2500"))
max_eval_per_sample <- as.integer(Sys.getenv("BLIND_AF_EVAL_EVENTS", "2500"))
max_samples <- as.integer(Sys.getenv("BLIND_AF_MAX_SAMPLES", "6"))
parse_quantiles <- function(x) {
  vals <- suppressWarnings(as.numeric(strsplit(x, ",", fixed = TRUE)[[1]]))
  vals <- vals[is.finite(vals) & vals > 0 & vals < 1]
  if (!length(vals)) vals <- c(0.80, 0.90, 0.95)
  unique(vals)
}
candidate_quantiles <- parse_quantiles(Sys.getenv("BLIND_AF_CANDIDATE_QUANTILES", "0.80,0.90,0.95"))

read_reference <- function(path) {
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  marker <- df[[1]]
  M <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(M) <- "numeric"
  rownames(M) <- marker
  M
}

af_rows <- function(M) grepl("^AF($|_)", rownames(M), ignore.case = TRUE)

drop_af <- function(M) M[!af_rows(M), , drop = FALSE]

row_normalize_max <- function(x) {
  x <- as.matrix(x)
  x[!is.finite(x)] <- 0
  x <- pmax(x, 0)
  mx <- apply(x, 1, max, na.rm = TRUE)
  keep <- is.finite(mx) & mx > 0
  x <- x[keep, , drop = FALSE]
  mx <- mx[keep]
  sweep(x, 1, mx, `/`)
}

cosine_matrix <- function(A, B) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  An <- sqrt(rowSums(A^2))
  Bn <- sqrt(rowSums(B^2))
  denom <- outer(An, Bn)
  out <- A %*% t(B) / denom
  out[!is.finite(out)] <- NA_real_
  out
}

prune_centers <- function(centers, marker_M, max_marker_cosine = 0.995, max_center_cosine = 0.985) {
  centers <- row_normalize_max(centers)
  if (!nrow(centers)) return(centers)

  marker_cos <- cosine_matrix(centers, marker_M)
  marker_max <- apply(marker_cos, 1, max, na.rm = TRUE)
  centers <- centers[!is.finite(marker_max) | marker_max < max_marker_cosine, , drop = FALSE]
  if (nrow(centers) <= 1) return(centers)

  keep <- integer()
  for (i in seq_len(nrow(centers))) {
    if (!length(keep)) {
      keep <- i
      next
    }
    center_cos <- cosine_matrix(centers[i, , drop = FALSE], centers[keep, , drop = FALSE])
    if (!any(is.finite(center_cos) & center_cos >= max_center_cosine)) {
      keep <- c(keep, i)
    }
  }
  centers[keep, , drop = FALSE]
}

subset_indices <- function(n, train_n, eval_n, seed) {
  set.seed(seed)
  total <- min(n, train_n + eval_n)
  idx <- sample.int(n, total)
  list(
    train = idx[seq_len(min(train_n, total))],
    eval = idx[seq.int(min(train_n, total) + 1L, total)]
  )
}

read_split_events <- function(path, M, train_n, eval_n, seed) {
  ff <- flowCore::read.FCS(path, transformation = FALSE, truncate_max_range = FALSE)
  detectors <- detector_columns(ff, M)
  n <- nrow(flowCore::exprs(ff))
  idx <- subset_indices(n, train_n, eval_n, seed)
  list(
    path = path,
    train = ff[sort(idx$train), ],
    eval = ff[sort(idx$eval), ],
    detectors = detectors
  )
}

learn_blind_af <- function(splits,
                           marker_M,
                           detector_noise,
                           k = 1L,
                           learn_method = "NNLS",
                           candidate_q = 0.80,
                           seed = 1L) {
  residual_shapes <- list()
  signal_values <- numeric()

  M_learn <- attach_detector_noise(marker_M, detector_noise)
  for (i in seq_along(splits)) {
    ff <- splits[[i]]$train
    if (nrow(flowCore::exprs(ff)) < 50) next
    res <- spectreasy::calc_residuals(
      ff,
      M_learn,
      method = learn_method,
      return_residuals = TRUE,
      rwls_max_iter = 3L
    )
    R <- as.matrix(res$residuals)
    Rpos <- pmax(R, 0)
    signal <- rowSums(Rpos)
    positive_fraction <- signal / (rowSums(abs(R)) + 1e-9)
    cutoff <- stats::quantile(signal[is.finite(signal)], candidate_q, na.rm = TRUE, names = FALSE)
    take <- is.finite(signal) & signal >= cutoff & positive_fraction >= 0.55
    if (sum(take) < 25) {
      take <- is.finite(signal) & signal >= cutoff
    }
    shapes <- row_normalize_max(Rpos[take, , drop = FALSE])
    if (nrow(shapes)) {
      residual_shapes[[length(residual_shapes) + 1L]] <- shapes
      signal_values <- c(signal_values, signal[take][seq_len(nrow(shapes))])
    }
  }

  shape_mat <- do.call(rbind, residual_shapes)
  if (is.null(shape_mat) || nrow(shape_mat) < 25) {
    return(NULL)
  }

  set.seed(seed)
  k <- min(as.integer(k), nrow(shape_mat))
  if (k <= 1L) {
    centers <- matrix(colMeans(shape_mat, na.rm = TRUE), nrow = 1)
  } else {
    km <- stats::kmeans(shape_mat, centers = k, nstart = 8, iter.max = 100)
    centers <- km$centers
  }
  colnames(centers) <- colnames(marker_M)
  centers <- prune_centers(centers, marker_M)
  if (!nrow(centers)) return(NULL)
  rownames(centers) <- paste0("AF_blind_", seq_len(nrow(centers)))
  rbind(marker_M, centers)
}

calc_metrics <- function(residuals, Y, detector_noise) {
  common <- calc_common_metrics(residuals, Y, detector_noise = detector_noise)
  common$Median_Event_RMSE <- stats::median(sqrt(rowMeans(residuals^2, na.rm = TRUE)), na.rm = TRUE)
  common$Q95_Event_RMSE <- stats::quantile(sqrt(rowMeans(residuals^2, na.rm = TRUE)), 0.95, na.rm = TRUE, names = FALSE)
  common
}

score_model <- function(dataset, splits, M, detector_noise, model_label, method_label = benchmark_method_rwls) {
  rows <- list()
  M_use <- attach_detector_noise(M, detector_noise)
  for (i in seq_along(splits)) {
    ff <- splits[[i]]$eval
    if (nrow(flowCore::exprs(ff)) < 50) next
    detectors <- detector_columns(ff, M)
    Y <- flowCore::exprs(ff)[, detectors, drop = FALSE]
    elapsed <- system.time({
      res <- unmix_benchmark_method(ff, M_use, detector_noise, method_label)
    })
    metrics <- calc_metrics(res$residuals, Y, detector_noise)
    rows[[length(rows) + 1L]] <- cbind(
      data.frame(
        Dataset = dataset,
        Sample = basename(splits[[i]]$path),
        Model = model_label,
        Method = method_label,
        Events = nrow(Y),
        Detectors = ncol(Y),
        Reference_Rows = nrow(M),
        AF_Rows = sum(af_rows(M)),
        Runtime_Sec = round(as.numeric(elapsed["elapsed"]), 3),
        stringsAsFactors = FALSE
      ),
      metrics
    )
  }
  do.call(rbind, rows)
}

summarize_results <- function(results) {
  metric_cols <- c("Raw_RMSE", "Raw_MAE", "Noise_Scaled_RMSE", "Signal_Scaled_RMSE", "Median_Event_RMSE", "Q95_Event_RMSE")
  aggregate(results[, metric_cols], results[c("Dataset", "Model", "AF_Rows")], mean, na.rm = TRUE)
}

compare_to_conventional <- function(summary) {
  out <- summary
  out$Raw_RMSE_vs_Conventional_Pct <- NA_real_
  out$Signal_Scaled_vs_Conventional_Pct <- NA_real_
  for (dataset in unique(out$Dataset)) {
    idx <- out$Dataset == dataset
    conv <- out[idx & out$Model == "conventional_multi_af", , drop = FALSE]
    if (!nrow(conv)) conv <- out[idx & out$Model == "conventional_single_af", , drop = FALSE]
    if (!nrow(conv)) next
    out$Raw_RMSE_vs_Conventional_Pct[idx] <- 100 * (out$Raw_RMSE[idx] - conv$Raw_RMSE[1]) / conv$Raw_RMSE[1]
    out$Signal_Scaled_vs_Conventional_Pct[idx] <- 100 * (out$Signal_Scaled_RMSE[idx] - conv$Signal_Scaled_RMSE[1]) / conv$Signal_Scaled_RMSE[1]
  }
  out
}

af_similarity <- function(dataset, conventional_M, blind_models) {
  conventional_af <- conventional_M[af_rows(conventional_M), , drop = FALSE]
  rows <- list()
  for (name in names(blind_models)) {
    M <- blind_models[[name]]
    if (is.null(M)) next
    blind_af <- M[grepl("^AF_blind_", rownames(M), ignore.case = TRUE), , drop = FALSE]
    if (!nrow(blind_af) || !nrow(conventional_af)) next
    cos <- cosine_matrix(blind_af, conventional_af)
    rows[[length(rows) + 1L]] <- data.frame(
      Dataset = dataset,
      Model = name,
      Blind_AF = rownames(blind_af),
      Best_Conventional_AF = colnames(t(conventional_af))[max.col(cos, ties.method = "first")],
      Best_Cosine = apply(cos, 1, max, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

run_dataset <- function(config) {
  message("\n=== ", config$name, " ===")
  conventional_multi <- read_reference(config$multi_ref)
  conventional_single <- read_reference(config$single_ref)
  detector_noise <- read_detector_noise(config$noise)
  marker_M <- drop_af(conventional_multi)

  sample_paths <- sort(config$samples)
  if (length(sample_paths) > max_samples) {
    sample_paths <- sample_paths[round(seq(1, length(sample_paths), length.out = max_samples))]
  }

  splits <- lapply(seq_along(sample_paths), function(i) {
    message("Loading sampled train/eval events: ", basename(sample_paths[i]))
    read_split_events(sample_paths[i], marker_M, max_train_per_sample, max_eval_per_sample, seed = 900 + i)
  })

  blind_models <- list()
  for (candidate_q in candidate_quantiles) {
    q_label <- paste0("q", sprintf("%02d", as.integer(round(candidate_q * 100))))
    for (k in c(1L, 3L, 6L, 10L)) {
      model_name <- paste0("blind_residual_", q_label, "_k", k)
      message("Learning ", model_name)
      blind_models[[model_name]] <- learn_blind_af(
        splits,
        marker_M,
        detector_noise,
        k = k,
        candidate_q = candidate_q,
        seed = 1000L + as.integer(round(candidate_q * 100)) + k
      )
    }
  }

  models <- c(
    list(
      marker_only_no_af = marker_M,
      conventional_single_af = conventional_single,
      conventional_multi_af = conventional_multi
    ),
    blind_models
  )
  models <- models[!vapply(models, is.null, logical(1))]

  result_rows <- list()
  for (model_name in names(models)) {
    message("Scoring ", model_name)
    result_rows[[length(result_rows) + 1L]] <- score_model(
      dataset = config$name,
      splits = splits,
      M = models[[model_name]],
      detector_noise = detector_noise,
      model_label = model_name
    )
  }

  list(
    results = do.call(rbind, result_rows),
    similarity = af_similarity(config$name, conventional_multi, blind_models),
    blind_models = blind_models
  )
}

pbmc_source <- "/Users/pkheisig/Library/CloudStorage/GoogleDrive-pkheisig@gmail.com/My Drive/PhD/Project_Erlangen/Erlangen_Data/PBMC_test_culture_titrations_second"
configs <- list(
  list(
    name = "PBMC-CULTURE-TITRATION",
    samples = list.files(file.path(pbmc_source, "samples"), pattern = "[.]fcs$", full.names = TRUE, ignore.case = TRUE),
    single_ref = file.path(benchmark_dir, "runs/PBMC-CULTURE-TITRATION/benchmark_outputs/references/single_af/scc_reference_matrix.csv"),
    multi_ref = file.path(benchmark_dir, "runs/PBMC-CULTURE-TITRATION/benchmark_outputs/references/multi_af/scc_reference_matrix.csv"),
    noise = file.path(benchmark_dir, "runs/PBMC-CULTURE-TITRATION/benchmark_outputs/references/multi_af/scc_detector_noise.csv")
  ),
  list(
    name = "MONOCYTE-XENITH-CELLS",
    samples = list.files(file.path(benchmark_dir, "runs/MONOCYTE-XENITH-CELLS/data/samples"), pattern = "[.]fcs$", full.names = TRUE, ignore.case = TRUE),
    single_ref = file.path(benchmark_dir, "runs/MONOCYTE-XENITH-CELLS/benchmark_outputs/references/single_af/scc_reference_matrix.csv"),
    multi_ref = file.path(benchmark_dir, "runs/MONOCYTE-XENITH-CELLS/benchmark_outputs/references/multi_af/scc_reference_matrix.csv"),
    noise = file.path(benchmark_dir, "runs/MONOCYTE-XENITH-CELLS/benchmark_outputs/references/multi_af/scc_detector_noise.csv")
  )
)

all_runs <- lapply(configs, run_dataset)
results <- do.call(rbind, lapply(all_runs, `[[`, "results"))
summary <- compare_to_conventional(summarize_results(results))
similarity <- do.call(rbind, lapply(all_runs, `[[`, "similarity"))

utils::write.csv(results, file.path(out_dir, "blind_af_event_metrics.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(out_dir, "blind_af_summary.csv"), row.names = FALSE)
utils::write.csv(similarity, file.path(out_dir, "blind_af_similarity_to_control_af.csv"), row.names = FALSE)

message("\nSummary:")
print(summary[order(summary$Dataset, summary$Signal_Scaled_RMSE), ], row.names = FALSE)
message("\nWrote blind-AF spike results to: ", out_dir)
