.blind_af_row_mask <- function(M) {
    grepl("^AF($|_)", rownames(M), ignore.case = TRUE)
}

.subset_reference_rows_preserve_attrs <- function(M, rows) {
    M <- .as_reference_matrix(M, "M")
    rows <- as.logical(rows)
    out <- M[rows, , drop = FALSE]

    variances <- attr(M, "variances")
    if (!is.null(variances)) {
        attr(out, "variances") <- variances[rownames(out), colnames(out), drop = FALSE]
    }

    detector_noise <- attr(M, "detector_noise")
    if (!is.null(detector_noise)) {
        attr(out, "detector_noise") <- detector_noise
    }

    out
}

.normalize_blind_af_shapes <- function(x) {
    x <- as.matrix(x)
    x[!is.finite(x)] <- 0
    x <- pmax(x, 0)
    row_scale <- apply(x, 1, max, na.rm = TRUE)
    keep <- is.finite(row_scale) & row_scale > 0
    if (!any(keep)) {
        return(x[FALSE, , drop = FALSE])
    }
    sweep(x[keep, , drop = FALSE], 1, row_scale[keep], `/`)
}

.blind_af_cosine_matrix <- function(A, B) {
    A <- as.matrix(A)
    B <- as.matrix(B)
    An <- sqrt(rowSums(A^2))
    Bn <- sqrt(rowSums(B^2))
    out <- A %*% t(B) / outer(An, Bn)
    out[!is.finite(out)] <- NA_real_
    out
}

.prune_blind_af_centers <- function(centers,
                                    marker_M,
                                    max_marker_cosine = 0.99,
                                    max_center_cosine = 0.985) {
    centers <- .normalize_blind_af_shapes(centers)
    if (nrow(centers) == 0) {
        return(centers)
    }

    marker_cos <- .blind_af_cosine_matrix(centers, marker_M)
    marker_max <- apply(marker_cos, 1, max, na.rm = TRUE)
    centers <- centers[!is.finite(marker_max) | marker_max < max_marker_cosine, , drop = FALSE]
    if (nrow(centers) <= 1) {
        return(centers)
    }

    keep <- integer()
    for (i in seq_len(nrow(centers))) {
        if (length(keep) == 0) {
            keep <- i
            next
        }
        center_cos <- .blind_af_cosine_matrix(centers[i, , drop = FALSE], centers[keep, , drop = FALSE])
        if (!any(is.finite(center_cos) & center_cos >= max_center_cosine)) {
            keep <- c(keep, i)
        }
    }

    centers[keep, , drop = FALSE]
}

.sample_blind_af_training_frame <- function(ff, max_events = 20000, seed = NULL) {
    n <- nrow(flowCore::exprs(ff))
    if (!is.finite(max_events) || is.na(max_events) || max_events <= 0 || n <= max_events) {
        return(ff)
    }
    if (!is.null(seed)) set.seed(seed)
    ff[sort(sample.int(n, max_events)), ]
}

.blind_af_training_frame_from_entry <- function(entry, max_events = 20000, seed = NULL) {
    ff <- if (inherits(entry$flow_frame, "flowFrame")) {
        entry$flow_frame
    } else {
        flowCore::read.FCS(entry$file_path, transformation = FALSE, truncate_max_range = FALSE)
    }
    .sample_blind_af_training_frame(ff, max_events = max_events, seed = seed)
}

.blind_af_split_frame <- function(ff,
                                  max_training_events = 20000L,
                                  max_evaluation_events = 5000L,
                                  train_fraction = 0.70,
                                  seed = NULL) {
    n <- nrow(flowCore::exprs(ff))
    if (n < 2L) {
        return(list(train = ff, eval = ff))
    }

    max_total <- max_training_events + max_evaluation_events
    total <- min(n, max_total)
    if (!is.null(seed)) set.seed(seed)
    idx <- sample.int(n, total)
    n_train <- min(max_training_events, max(1L, floor(total * train_fraction)))
    n_eval <- min(max_evaluation_events, total - n_train)
    if (n_eval < 1L) {
        n_eval <- min(1L, total)
        n_train <- max(1L, total - n_eval)
    }

    train_idx <- sort(idx[seq_len(n_train)])
    eval_idx <- sort(idx[seq.int(n_train + 1L, n_train + n_eval)])
    list(train = ff[train_idx, ], eval = ff[eval_idx, ])
}

.blind_af_split_entry <- function(entry,
                                  max_training_events = 20000L,
                                  max_evaluation_events = 5000L,
                                  seed = NULL) {
    ff <- if (inherits(entry$flow_frame, "flowFrame")) {
        entry$flow_frame
    } else {
        flowCore::read.FCS(entry$file_path, transformation = FALSE, truncate_max_range = FALSE)
    }
    .blind_af_split_frame(
        ff,
        max_training_events = max_training_events,
        max_evaluation_events = max_evaluation_events,
        seed = seed
    )
}

.cluster_blind_af_shapes <- function(shape_mat,
                                     n_bands = 10L,
                                     min_cluster_events = 20L,
                                     min_cluster_proportion = 0.005,
                                     seed = NULL) {
    shape_mat <- .normalize_blind_af_shapes(shape_mat)
    if (nrow(shape_mat) == 0) {
        return(shape_mat)
    }

    rounded <- as.data.frame(round(shape_mat, digits = 6))
    distinct_n <- nrow(unique(rounded))
    n_eff <- min(as.integer(n_bands), nrow(shape_mat), max(distinct_n, 1L))
    if (!is.finite(n_eff) || is.na(n_eff) || n_eff < 1L) {
        n_eff <- 1L
    }

    if (n_eff == 1L) {
        centers <- matrix(colMeans(shape_mat, na.rm = TRUE), nrow = 1)
        colnames(centers) <- colnames(shape_mat)
        return(.normalize_blind_af_shapes(centers))
    }

    if (!is.null(seed)) set.seed(seed)
    km <- stats::kmeans(shape_mat, centers = n_eff, nstart = 8, iter.max = 100)
    cluster_sizes <- as.numeric(table(factor(km$cluster, levels = seq_len(n_eff))))
    min_cluster_size <- .reference_min_af_cluster_size(
        n_events = nrow(shape_mat),
        min_cluster_events = min_cluster_events,
        min_cluster_proportion = min_cluster_proportion
    )
    keep_idx <- which(cluster_sizes >= min_cluster_size)
    if (length(keep_idx) == 0) {
        keep_idx <- which.max(cluster_sizes)
    }

    centers <- km$centers[keep_idx, , drop = FALSE]
    cluster_sizes <- cluster_sizes[keep_idx]
    centers <- centers[order(cluster_sizes, decreasing = TRUE), , drop = FALSE]
    colnames(centers) <- colnames(shape_mat)
    .normalize_blind_af_shapes(centers)
}

.append_blind_af_to_reference <- function(marker_M, af_centers) {
    marker_M <- .as_reference_matrix(marker_M, "marker_M")
    af_centers <- .normalize_blind_af_shapes(af_centers)
    if (nrow(af_centers) == 0) {
        return(marker_M)
    }

    rownames(af_centers) <- c("AF", if (nrow(af_centers) > 1L) paste0("AF_", seq.int(2L, nrow(af_centers))) else NULL)
    colnames(af_centers) <- colnames(marker_M)
    out <- rbind(marker_M, af_centers)

    variances <- attr(marker_M, "variances")
    if (!is.null(variances)) {
        af_variances <- matrix(
            0,
            nrow = nrow(af_centers),
            ncol = ncol(marker_M),
            dimnames = list(rownames(af_centers), colnames(marker_M))
        )
        attr(out, "variances") <- rbind(variances[rownames(marker_M), colnames(marker_M), drop = FALSE], af_variances)
    }

    detector_noise <- attr(marker_M, "detector_noise")
    if (!is.null(detector_noise)) {
        attr(out, "detector_noise") <- detector_noise
    }

    out
}

.blind_af_candidate_shapes <- function(res_obj,
                                       candidate_quantile = 0.90,
                                       strategy = c("residual", "low_marker"),
                                       min_events = 25L) {
    strategy <- match.arg(strategy)
    R <- as.matrix(res_obj$residuals)
    Rpos <- pmax(R, 0)
    signal <- rowSums(Rpos)
    total_residual <- rowSums(abs(R)) + 1e-9
    positive_fraction <- signal / total_residual

    score <- signal
    if (identical(strategy, "low_marker")) {
        marker_cols <- setdiff(
            colnames(res_obj$data),
            c("File", grep("^(Time|FSC|SSC)", colnames(res_obj$data), ignore.case = TRUE, value = TRUE))
        )
        if (length(marker_cols) > 0) {
            marker_burden <- rowSums(pmax(as.matrix(res_obj$data[, marker_cols, drop = FALSE]), 0), na.rm = TRUE)
            score <- signal / sqrt(marker_burden + 1)
        }
    }

    score_finite <- score[is.finite(score)]
    if (length(score_finite) < min_events) {
        return(NULL)
    }

    effective_quantile <- min(candidate_quantile, max(0, 1 - (min_events / length(score_finite))))
    cutoff <- stats::quantile(score_finite, effective_quantile, na.rm = TRUE, names = FALSE)
    take <- is.finite(score) & score >= cutoff & positive_fraction >= 0.55
    if (sum(take) < min_events) {
        take <- is.finite(score) & score >= cutoff
    }

    .normalize_blind_af_shapes(Rpos[take, , drop = FALSE])
}

.learn_blind_af_model <- function(splits,
                                  marker_M,
                                  model_id = "residual_wls_q90_k10",
                                  n_bands = 10L,
                                  candidate_quantile = 0.90,
                                  strategy = c("residual", "low_marker"),
                                  min_events = 25L,
                                  seed = NULL) {
    strategy <- match.arg(strategy)
    candidate_shapes <- list()
    candidate_counts <- integer(length(splits))

    for (i in seq_along(splits)) {
        ff <- splits[[i]]$train
        if (nrow(flowCore::exprs(ff)) < min_events) {
            next
        }

        res <- tryCatch(
            calc_residuals(ff, marker_M, method = "WLS", return_residuals = TRUE),
            error = function(e) NULL
        )
        if (is.null(res) || is.null(res$residuals)) {
            next
        }

        shapes <- .blind_af_candidate_shapes(
            res,
            candidate_quantile = candidate_quantile,
            strategy = strategy,
            min_events = min_events
        )
        if (!is.null(shapes) && nrow(shapes) > 0) {
            candidate_shapes[[length(candidate_shapes) + 1L]] <- shapes
            candidate_counts[[i]] <- nrow(shapes)
        }
    }

    shape_mat <- do.call(rbind, candidate_shapes)
    if (is.null(shape_mat) || nrow(shape_mat) < min_events) {
        return(NULL)
    }

    centers <- .cluster_blind_af_shapes(
        shape_mat,
        n_bands = n_bands,
        min_cluster_events = 20L,
        min_cluster_proportion = 0.005,
        seed = seed
    )
    centers <- .prune_blind_af_centers(centers, marker_M)
    if (nrow(centers) == 0) {
        return(NULL)
    }

    out <- .append_blind_af_to_reference(marker_M, centers)
    attr(out, "blind_af_info") <- list(
        model_id = model_id,
        requested_bands = as.integer(n_bands),
        derived_bands = sum(.blind_af_row_mask(out)),
        candidate_events = nrow(shape_mat),
        candidate_quantile = candidate_quantile,
        candidate_strategy = strategy,
        sample_candidate_counts = candidate_counts
    )
    out
}

.blind_af_signal_scaled_rmse <- function(residuals, Y, M) {
    noise <- .resolve_wls_noise_parameters(M)
    denom <- sweep(
        pmax(Y, 0),
        2,
        noise$signal_scale,
        FUN = `*`
    )
    denom <- sweep(
        denom,
        2,
        noise$noise_floor,
        FUN = `+`
    )
    denom[!is.finite(denom) | denom <= 0] <- .default_wls_background_noise()
    sqrt(mean((residuals^2) / denom, na.rm = TRUE))
}

.score_blind_af_model <- function(M, splits, af_assignment = "projection") {
    scores <- numeric()
    for (split in splits) {
        ff <- split$eval
        if (nrow(flowCore::exprs(ff)) < 1L) {
            next
        }
        res <- tryCatch(
            calc_residuals(ff, M, method = "WLS", return_residuals = TRUE, af_assignment = af_assignment),
            error = function(e) NULL
        )
        if (is.null(res) || is.null(res$residuals)) {
            next
        }
        Y <- flowCore::exprs(ff)[, colnames(M), drop = FALSE]
        scores <- c(scores, .blind_af_signal_scaled_rmse(res$residuals, Y, M))
    }
    if (length(scores) == 0) {
        return(NA_real_)
    }
    mean(scores, na.rm = TRUE)
}

.estimate_blind_af_reference <- function(M,
                                         sample_entries,
                                         n_bands = 10L,
                                         candidate_quantile = 0.90,
                                         max_training_events = 20000L,
                                         max_evaluation_events = 5000L,
                                         min_events = 25L,
                                         seed = NULL,
                                         af_assignment = "projection",
                                         verbose = TRUE) {
    M <- .as_reference_matrix(M, "M")
    af_assignment <- .normalize_af_assignment(af_assignment, choices = c("projection", "residual_alignment", "legacy_residual"))
    existing_af <- .blind_af_row_mask(M)
    if (all(existing_af)) {
        warning("estimate_af = TRUE requires at least one non-AF marker row; proceeding with the existing reference matrix.", call. = FALSE)
        return(M)
    }
    if (any(existing_af)) {
        warning("estimate_af = TRUE: replacing existing AF rows with AF estimated from stained samples.", call. = FALSE)
    }

    marker_M <- .subset_reference_rows_preserve_attrs(M, !existing_af)

    splits <- lapply(seq_along(sample_entries), function(i) {
        .blind_af_split_entry(
            sample_entries[[i]],
            max_training_events = max_training_events,
            max_evaluation_events = max_evaluation_events,
            seed = if (!is.null(seed)) seed + i else NULL
        )
    })

    model_specs <- list(
        list(id = "residual_wls_q90_k10", strategy = "residual"),
        list(id = "residual_wls_low_marker_q90_k10", strategy = "low_marker")
    )
    models <- list()
    score_rows <- list()
    for (spec in model_specs) {
        model <- .learn_blind_af_model(
            splits = splits,
            marker_M = marker_M,
            model_id = spec$id,
            n_bands = n_bands,
            candidate_quantile = candidate_quantile,
            strategy = spec$strategy,
            min_events = min_events,
            seed = seed
        )
        if (is.null(model)) {
            score_rows[[length(score_rows) + 1L]] <- data.frame(
                model_id = spec$id,
                heldout_signal_scaled_rmse = NA_real_,
                derived_bands = 0L,
                stringsAsFactors = FALSE
            )
            next
        }
        score <- .score_blind_af_model(
            model,
            splits,
            af_assignment = if (identical(af_assignment, "legacy_residual")) "legacy" else af_assignment
        )
        models[[spec$id]] <- model
        score_rows[[length(score_rows) + 1L]] <- data.frame(
            model_id = spec$id,
            heldout_signal_scaled_rmse = score,
            derived_bands = sum(.blind_af_row_mask(model)),
            stringsAsFactors = FALSE
        )
    }

    score_table <- do.call(rbind, score_rows)
    valid <- is.finite(score_table$heldout_signal_scaled_rmse)
    if (!any(valid)) {
        fallback <- models[["residual_wls_q90_k10"]]
        if (!is.null(fallback)) {
            out <- fallback
            info <- attr(out, "blind_af_info")
            info$selected_by <- "fallback"
            info$af_assignment <- af_assignment
            info$heldout_scores <- score_table
            attr(out, "blind_af_info") <- info
            return(out)
        }
        warning("estimate_af = TRUE could not find enough AF-like residual events; proceeding without estimated AF.", call. = FALSE)
        return(marker_M)
    }

    selected_id <- score_table$model_id[valid][which.min(score_table$heldout_signal_scaled_rmse[valid])]
    out <- models[[selected_id]]
    info <- attr(out, "blind_af_info")
    info$selected_by <- if (identical(af_assignment, "legacy_residual")) {
        "heldout_event_wise_wls"
    } else {
        paste0("heldout_", af_assignment, "_wls")
    }
    info$af_assignment <- af_assignment
    info$heldout_signal_scaled_rmse <- score_table$heldout_signal_scaled_rmse[match(selected_id, score_table$model_id)]
    info$heldout_scores <- score_table
    attr(out, "blind_af_info") <- info

    if (isTRUE(verbose)) {
        message(
            "  Estimated ",
            sum(.blind_af_row_mask(out)),
            " AF signature(s) from stained sample residuals using ",
            selected_id,
            "."
        )
    }

    out
}
