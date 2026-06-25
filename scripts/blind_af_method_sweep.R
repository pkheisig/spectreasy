#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(flowCore)
  library(spectreasy)
  library(nnls)
})

benchmark_dir <- "/Users/pkheisig/Building/Project_Spectreasy/benchmarking"
source(file.path(benchmark_dir, "run_wls_benchmark.R"), chdir = TRUE)

out_dir <- file.path(benchmark_dir, "runs", "BLIND-AF-METHOD-SWEEP")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

max_train_per_sample <- as.integer(Sys.getenv("BLIND_AF_TRAIN_EVENTS", "1800"))
max_eval_per_sample <- as.integer(Sys.getenv("BLIND_AF_EVAL_EVENTS", "1800"))
max_samples <- as.integer(Sys.getenv("BLIND_AF_MAX_SAMPLES", "6"))
nmf_max_train_events <- as.integer(Sys.getenv("BLIND_AF_NMF_EVENTS", "3000"))
final_method <- Sys.getenv("BLIND_AF_FINAL_METHOD", benchmark_method_event_wise_wls)

read_reference <- function(path) {
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  M <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(M) <- "numeric"
  rownames(M) <- df[[1]]
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
  if (!any(keep)) return(x[FALSE, , drop = FALSE])
  sweep(x[keep, , drop = FALSE], 1, mx[keep], `/`)
}

cosine_matrix <- function(A, B) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  out <- A %*% t(B) / outer(sqrt(rowSums(A^2)), sqrt(rowSums(B^2)))
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
    if (!any(is.finite(center_cos) & center_cos >= max_center_cosine)) keep <- c(keep, i)
  }
  centers[keep, , drop = FALSE]
}

cluster_shapes <- function(shape_mat, k, seed = 1L, min_cluster_events = 20L) {
  shape_mat <- row_normalize_max(shape_mat)
  if (is.null(shape_mat) || nrow(shape_mat) < min_cluster_events) return(NULL)
  distinct_n <- nrow(unique(as.data.frame(round(shape_mat, 6))))
  k <- min(as.integer(k), nrow(shape_mat), max(distinct_n, 1L))
  set.seed(seed)
  if (k <= 1L) {
    centers <- matrix(colMeans(shape_mat, na.rm = TRUE), nrow = 1)
  } else {
    km <- stats::kmeans(shape_mat, centers = k, nstart = 8, iter.max = 100)
    cluster_sizes <- as.numeric(table(factor(km$cluster, levels = seq_len(k))))
    keep <- which(cluster_sizes >= min_cluster_events)
    if (!length(keep)) keep <- which.max(cluster_sizes)
    centers <- km$centers[keep, , drop = FALSE]
    centers <- centers[order(cluster_sizes[keep], decreasing = TRUE), , drop = FALSE]
  }
  colnames(centers) <- colnames(shape_mat)
  row_normalize_max(centers)
}

append_af <- function(marker_M, centers) {
  centers <- row_normalize_max(centers)
  if (is.null(centers) || !nrow(centers)) return(NULL)
  rownames(centers) <- paste0("AF_blind_", seq_len(nrow(centers)))
  colnames(centers) <- colnames(marker_M)
  rbind(marker_M, centers)
}

subset_indices <- function(n, train_n, eval_n, seed) {
  set.seed(seed)
  total <- min(n, train_n + eval_n)
  idx <- sample.int(n, total)
  list(train = idx[seq_len(min(train_n, total))], eval = idx[seq.int(min(train_n, total) + 1L, total)])
}

read_split_events <- function(path, M, train_n, eval_n, seed) {
  ff <- flowCore::read.FCS(path, transformation = FALSE, truncate_max_range = FALSE)
  detector_columns(ff, M)
  n <- nrow(flowCore::exprs(ff))
  idx <- subset_indices(n, train_n, eval_n, seed)
  list(path = path, train = ff[sort(idx$train), ], eval = ff[sort(idx$eval), ])
}

candidate_shapes_from_residuals <- function(res, candidate_q = 0.95, strategy = "residual") {
  R <- as.matrix(res$residuals)
  Rpos <- pmax(R, 0)
  signal <- rowSums(Rpos)
  positive_fraction <- signal / (rowSums(abs(R)) + 1e-9)
  score <- signal
  if (identical(strategy, "low_marker")) {
    marker_cols <- setdiff(colnames(res$data), grep("^(Time|FSC|SSC|File)", colnames(res$data), value = TRUE, ignore.case = TRUE))
    burden <- if (length(marker_cols)) rowSums(pmax(as.matrix(res$data[, marker_cols, drop = FALSE]), 0), na.rm = TRUE) else 0
    score <- signal / sqrt(burden + 1)
  }
  score_finite <- score[is.finite(score)]
  if (length(score_finite) < 25) return(NULL)
  effective_q <- min(candidate_q, max(0, 1 - 25 / length(score_finite)))
  cutoff <- stats::quantile(score_finite, effective_q, na.rm = TRUE, names = FALSE)
  take <- is.finite(score) & score >= cutoff & positive_fraction >= 0.55
  if (sum(take) < 25) take <- is.finite(score) & score >= cutoff
  row_normalize_max(Rpos[take, , drop = FALSE])
}

learn_residual_global <- function(splits, marker_M, detector_noise, k, candidate_q, learn_method, strategy, seed) {
  M_learn <- attach_detector_noise(marker_M, detector_noise)
  shapes <- list()
  for (i in seq_along(splits)) {
    method <- learn_method
    res <- spectreasy::calc_residuals(splits[[i]]$train, M_learn, method = method, return_residuals = TRUE)
    s <- candidate_shapes_from_residuals(res, candidate_q = candidate_q, strategy = strategy)
    if (!is.null(s) && nrow(s)) shapes[[length(shapes) + 1L]] <- s
  }
  shape_mat <- do.call(rbind, shapes)
  centers <- cluster_shapes(shape_mat, k = k, seed = seed)
  centers <- prune_centers(centers, marker_M)
  append_af(marker_M, centers)
}

learn_residual_per_sample <- function(splits, marker_M, detector_noise, k_total, candidate_q, learn_method, seed) {
  M_learn <- attach_detector_noise(marker_M, detector_noise)
  centers <- list()
  k_each <- max(1L, ceiling(k_total / max(1L, length(splits))))
  for (i in seq_along(splits)) {
    res <- spectreasy::calc_residuals(splits[[i]]$train, M_learn, method = learn_method, return_residuals = TRUE)
    s <- candidate_shapes_from_residuals(res, candidate_q = candidate_q, strategy = "residual")
    c_i <- cluster_shapes(s, k = k_each, seed = seed + i)
    if (!is.null(c_i) && nrow(c_i)) centers[[length(centers) + 1L]] <- c_i
  }
  pooled <- do.call(rbind, centers)
  final <- cluster_shapes(pooled, k = k_total, seed = seed + 100L, min_cluster_events = 1L)
  final <- prune_centers(final, marker_M)
  append_af(marker_M, final)
}

scatter_groups <- function(ff, max_groups = 3L) {
  X <- flowCore::exprs(ff)
  scatter <- grep("^(FSC|SSC)", colnames(X), ignore.case = TRUE, value = TRUE)
  if (length(scatter) < 2 || nrow(X) < 100) return(rep(1L, nrow(X)))
  scatter <- scatter[seq_len(2)]
  Z <- scale(log1p(pmax(X[, scatter, drop = FALSE], 0)))
  Z[!is.finite(Z)] <- 0
  k <- min(max_groups, max(1L, nrow(unique(as.data.frame(round(Z, 3))))))
  if (k <= 1L) return(rep(1L, nrow(X)))
  stats::kmeans(Z, centers = k, nstart = 4, iter.max = 50)$cluster
}

learn_residual_scatter_population <- function(splits, marker_M, detector_noise, k_total, candidate_q, learn_method, seed) {
  M_learn <- attach_detector_noise(marker_M, detector_noise)
  centers <- list()
  for (i in seq_along(splits)) {
    ff <- splits[[i]]$train
    groups <- scatter_groups(ff)
    for (g in sort(unique(groups))) {
      idx <- which(groups == g)
      if (length(idx) < 80) next
      res <- spectreasy::calc_residuals(ff[idx, ], M_learn, method = learn_method, return_residuals = TRUE)
      s <- candidate_shapes_from_residuals(res, candidate_q = candidate_q, strategy = "residual")
      c_i <- cluster_shapes(s, k = 2L, seed = seed + i + g)
      if (!is.null(c_i) && nrow(c_i)) centers[[length(centers) + 1L]] <- c_i
    }
  }
  pooled <- do.call(rbind, centers)
  final <- cluster_shapes(pooled, k = k_total, seed = seed + 200L, min_cluster_events = 1L)
  final <- prune_centers(final, marker_M)
  append_af(marker_M, final)
}

sample_train_matrix <- function(splits, marker_M, max_events, seed = 1L) {
  Y_parts <- lapply(splits, function(s) flowCore::exprs(s$train)[, colnames(marker_M), drop = FALSE])
  Y <- do.call(rbind, Y_parts)
  if (nrow(Y) > max_events) {
    set.seed(seed)
    Y <- Y[sort(sample.int(nrow(Y), max_events)), , drop = FALSE]
  }
  Y
}

nnls_coefficients <- function(Y, M) {
  design <- t(M)
  A <- matrix(0, nrow = nrow(Y), ncol = nrow(M), dimnames = list(NULL, rownames(M)))
  for (i in seq_len(nrow(Y))) {
    A[i, ] <- nnls::nnls(design, pmax(as.numeric(Y[i, ]), 0))$x
  }
  A
}

learn_semisupervised_nnls <- function(splits, marker_M, detector_noise, k, candidate_q, seed) {
  init <- learn_residual_global(splits, marker_M, detector_noise, k = k, candidate_q = candidate_q, learn_method = "WLS", strategy = "residual", seed = seed)
  if (is.null(init)) return(NULL)
  B <- init[af_rows(init), , drop = FALSE]
  if (!nrow(B)) return(NULL)
  Y <- sample_train_matrix(splits, marker_M, max_events = nmf_max_train_events, seed = seed)

  for (iter in seq_len(4L)) {
    M_all <- rbind(marker_M, B)
    A <- nnls_coefficients(Y, M_all)
    A_marker <- A[, seq_len(nrow(marker_M)), drop = FALSE]
    A_af <- A[, seq.int(nrow(marker_M) + 1L, ncol(A)), drop = FALSE]
    residual <- pmax(Y - A_marker %*% marker_M, 0)
    B_new <- matrix(0, nrow = ncol(A_af), ncol = ncol(marker_M), dimnames = list(rownames(B), colnames(marker_M)))
    usable <- colSums(A_af) > 1e-8
    if (!any(usable)) break
    for (d in seq_len(ncol(residual))) {
      B_new[usable, d] <- nnls::nnls(A_af[, usable, drop = FALSE], residual[, d])$x
    }
    B <- row_normalize_max(B_new)
    if (!nrow(B)) break
  }

  B <- prune_centers(B, marker_M)
  append_af(marker_M, B)
}

calc_metrics2 <- function(residuals, Y, detector_noise) {
  common <- calc_common_metrics(residuals, Y, detector_noise = detector_noise)
  common$Median_Event_RMSE <- stats::median(sqrt(rowMeans(residuals^2, na.rm = TRUE)), na.rm = TRUE)
  common$Q95_Event_RMSE <- stats::quantile(sqrt(rowMeans(residuals^2, na.rm = TRUE)), 0.95, na.rm = TRUE, names = FALSE)
  common
}

score_model <- function(dataset, splits, M, detector_noise, model_label, method_label = final_method) {
  rows <- list()
  M_use <- attach_detector_noise(M, detector_noise)
  for (i in seq_along(splits)) {
    ff <- splits[[i]]$eval
    if (nrow(flowCore::exprs(ff)) < 50) next
    Y <- flowCore::exprs(ff)[, colnames(M), drop = FALSE]
    elapsed <- system.time(res <- unmix_benchmark_method(ff, M_use, detector_noise, method_label))
    metrics <- calc_metrics2(res$residuals, Y, detector_noise)
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
    conv <- out[idx & out$Model == "conventional_multi_af_10", , drop = FALSE]
    if (!nrow(conv)) next
    out$Raw_RMSE_vs_Conventional_Pct[idx] <- 100 * (out$Raw_RMSE[idx] - conv$Raw_RMSE[1]) / conv$Raw_RMSE[1]
    out$Signal_Scaled_vs_Conventional_Pct[idx] <- 100 * (out$Signal_Scaled_RMSE[idx] - conv$Signal_Scaled_RMSE[1]) / conv$Signal_Scaled_RMSE[1]
  }
  out
}

af_similarity <- function(dataset, conventional_M, models) {
  conventional_af <- conventional_M[af_rows(conventional_M), , drop = FALSE]
  rows <- list()
  for (name in names(models)) {
    M <- models[[name]]
    if (is.null(M) || name %in% c("marker_only_no_af", "conventional_single_af", "conventional_multi_af_10")) next
    blind_af <- M[af_rows(M), , drop = FALSE]
    if (!nrow(blind_af) || !nrow(conventional_af)) next
    cos <- cosine_matrix(blind_af, conventional_af)
    rows[[length(rows) + 1L]] <- data.frame(
      Dataset = dataset,
      Model = name,
      Blind_AF = rownames(blind_af),
      Best_Conventional_AF = rownames(conventional_af)[max.col(cos, ties.method = "first")],
      Best_Cosine = apply(cos, 1, max, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

build_models <- function(splits, marker_M, conventional_single, conventional_multi, detector_noise) {
  models <- list(
    marker_only_no_af = marker_M,
    conventional_single_af = conventional_single,
    conventional_multi_af_10 = conventional_multi
  )

  q_grid <- c(0.90, 0.95)
  k_grid <- c(6L, 10L)
  for (q in q_grid) {
    q_label <- paste0("q", as.integer(q * 100))
    for (k in k_grid) {
      message("Learning residual NNLS ", q_label, " k", k)
      models[[paste0("residual_nnls_", q_label, "_k", k)]] <- learn_residual_global(splits, marker_M, detector_noise, k, q, "NNLS", "residual", seed = 100 + k)
      message("Learning residual WLS ", q_label, " k", k)
      models[[paste0("residual_wls_", q_label, "_k", k)]] <- learn_residual_global(splits, marker_M, detector_noise, k, q, "WLS", "residual", seed = 200 + k)
      message("Learning residual WLS low-marker ", q_label, " k", k)
      models[[paste0("residual_wls_low_marker_", q_label, "_k", k)]] <- learn_residual_global(splits, marker_M, detector_noise, k, q, "WLS", "low_marker", seed = 300 + k)
    }
  }

  message("Learning per-sample WLS")
  models$residual_wls_per_sample_q95_k10 <- learn_residual_per_sample(splits, marker_M, detector_noise, 10L, 0.95, "WLS", seed = 400)

  message("Learning scatter-population WLS")
  models$residual_wls_scatter_q95_k10 <- learn_residual_scatter_population(splits, marker_M, detector_noise, 10L, 0.95, "WLS", seed = 500)

  message("Learning semi-supervised NNLS")
  models$semisupervised_nnls_q95_k10 <- learn_semisupervised_nnls(splits, marker_M, detector_noise, 10L, 0.95, seed = 600)

  models[!vapply(models, is.null, logical(1))]
}

run_dataset <- function(config) {
  message("\n=== ", config$name, " ===")
  conventional_multi <- read_reference(config$multi_ref)
  conventional_single <- read_reference(config$single_ref)
  detector_noise <- read_detector_noise(config$noise)
  marker_M <- drop_af(conventional_multi)

  sample_paths <- sort(config$samples)
  if (length(sample_paths) > max_samples) sample_paths <- sample_paths[round(seq(1, length(sample_paths), length.out = max_samples))]
  splits <- lapply(seq_along(sample_paths), function(i) {
    message("Loading sampled train/eval events: ", basename(sample_paths[i]))
    read_split_events(sample_paths[i], marker_M, max_train_per_sample, max_eval_per_sample, seed = 1900 + i)
  })

  models <- build_models(splits, marker_M, conventional_single, conventional_multi, detector_noise)
  results <- do.call(rbind, lapply(names(models), function(model_name) {
    message("Scoring ", model_name)
    score_model(config$name, splits, models[[model_name]], detector_noise, model_name)
  }))

  list(
    results = results,
    similarity = af_similarity(config$name, conventional_multi, models)
  )
}

pbmc_source <- "/Users/pkheisig/Library/CloudStorage/GoogleDrive-pkheisig@gmail.com/My Drive/PhD/Project_Erlangen/Erlangen_Data/PBMC_test_culture_titrations_second"
configs <- list(
  list(
    name = "PBMC-CULTURE-TITRATION",
    samples = list.files(file.path(pbmc_source, "samples"), pattern = "[.]fcs$", full.names = TRUE, ignore.case = TRUE),
    single_ref = file.path(benchmark_dir, "runs/PBMC-CULTURE-TITRATION/benchmark_outputs/references/single_af/scc_reference_matrix.csv"),
    multi_ref = file.path(benchmark_dir, "runs/PBMC-CULTURE-TITRATION/benchmark_outputs/references/multi_af_10/scc_reference_matrix.csv"),
    noise = file.path(benchmark_dir, "runs/PBMC-CULTURE-TITRATION/benchmark_outputs/references/multi_af_10/scc_detector_noise.csv")
  ),
  list(
    name = "MONOCYTE-XENITH-CELLS",
    samples = list.files(file.path(benchmark_dir, "runs/MONOCYTE-XENITH-CELLS/data/samples"), pattern = "[.]fcs$", full.names = TRUE, ignore.case = TRUE),
    single_ref = file.path(benchmark_dir, "runs/MONOCYTE-XENITH-CELLS/benchmark_outputs/references/single_af/scc_reference_matrix.csv"),
    multi_ref = file.path(benchmark_dir, "runs/MONOCYTE-XENITH-CELLS/benchmark_outputs/references/multi_af_10/scc_reference_matrix.csv"),
    noise = file.path(benchmark_dir, "runs/MONOCYTE-XENITH-CELLS/benchmark_outputs/references/multi_af_10/scc_detector_noise.csv")
  )
)

all_runs <- lapply(configs, run_dataset)
results <- do.call(rbind, lapply(all_runs, `[[`, "results"))
summary <- compare_to_conventional(summarize_results(results))
similarity <- do.call(rbind, lapply(all_runs, `[[`, "similarity"))

utils::write.csv(results, file.path(out_dir, "method_sweep_event_metrics.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(out_dir, "method_sweep_summary.csv"), row.names = FALSE)
utils::write.csv(similarity, file.path(out_dir, "method_sweep_similarity_to_control_af.csv"), row.names = FALSE)

message("\nTop models by Signal_Scaled_RMSE:")
for (dataset in unique(summary$Dataset)) {
  cat("\n", dataset, "\n", sep = "")
  print(
    head(summary[summary$Dataset == dataset, ][order(summary$Signal_Scaled_RMSE[summary$Dataset == dataset]), ], 10),
    row.names = FALSE
  )
}
message("\nWrote method sweep results to: ", out_dir)
