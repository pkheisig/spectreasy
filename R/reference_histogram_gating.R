.compute_reference_histogram_gate <- function(peak_vals,
                                              sample_type,
                                              histogram_pct_beads,
                                              histogram_direction_beads,
                                              histogram_pct_cells,
                                              histogram_direction_cells,
                                              is_viability = FALSE) {
    vals_log <- log10(pmax(peak_vals, 1))
    vals_log <- vals_log[is.finite(vals_log)]
    positive_gate_present <- !(sample_type %in% c("unstained"))
    vals_for_gate <- if (sample_type %in% c("unstained", "cells")) vals_log else .select_reference_hist_peak_population(vals_log)
    pct <- if (sample_type %in% c("unstained", "cells")) histogram_pct_cells else histogram_pct_beads
    dir <- if (sample_type %in% c("unstained", "cells")) histogram_direction_cells else histogram_direction_beads
    pct <- min(max(pct, 0.01), 0.999)

    neg_log_min <- min(vals_log, na.rm = TRUE)
    neg_log_max <- as.numeric(stats::quantile(vals_log, 0.15, na.rm = TRUE))
    neg_gate_method <- "density lower-tail fallback"
    negative_gate_present <- positive_gate_present
    negative_components <- NULL
    density_gate <- NULL

    if (positive_gate_present) {
        density <- .reference_histogram_density_gate(vals_log, sample_type, is_viability, dir)
        if (!is.null(density)) {
            density_gate <- density$gate
            neg_log_min <- density$neg_log_min
            neg_log_max <- density$neg_log_max
            neg_gate_method <- density$method
        }
    }

    negative <- list(
        min = neg_log_min,
        max = neg_log_max,
        method = neg_gate_method,
        components = negative_components
    )
    if (positive_gate_present && is.null(density_gate)) {
        negative <- .reference_histogram_negative_gmm(vals_log, negative)
        neg_log_min <- negative$min
        neg_log_max <- negative$max
        neg_gate_method <- negative$method
        negative_components <- negative$components
    }

    fit_vals <- vals_for_gate[is.finite(vals_for_gate)]
    model_info <- list(method = "quantile fallback", components = negative_components)
    gmm_upper <- NA_real_
    if (positive_gate_present) {
        bright <- .reference_histogram_bright_gmm(fit_vals, negative)
        model_info <- bright$info
        gmm_upper <- bright$upper
    }

    if (dir == "right") {
        lq <- 0.5
        uq <- min(1, 0.5 + pct)
    } else if (dir == "left") {
        lq <- max(0, 0.5 - pct)
        uq <- 0.5
    } else {
        lq <- 0.5 - pct / 2
        uq <- 0.5 + pct / 2
    }
    q_gate <- c(
        lower = as.numeric(stats::quantile(vals_for_gate, max(0, lq), na.rm = TRUE)),
        upper = as.numeric(stats::quantile(vals_for_gate, min(1, uq), na.rm = TRUE))
    )
    if (!positive_gate_present) {
        lower_log <- min(vals_log, na.rm = TRUE)
        upper_log <- max(vals_log, na.rm = TRUE)
        negative_gate_present <- FALSE
        model_info <- list(
            method = "AF/unstained control: scatter-gated only; no positive histogram gate",
            components = negative_components
        )
    } else if (!is.null(density_gate) && is.finite(density_gate$pos_lower) && is.finite(density_gate$pos_upper) && density_gate$pos_upper > density_gate$pos_lower) {
        lower_log <- density_gate$pos_lower
        upper_log <- density_gate$pos_upper
        component_df <- data.frame(
            component = seq_along(density_gate$component_means),
            n = NA_integer_,
            prop = NA_real_,
            mean = density_gate$component_means,
            sd = NA_real_,
            q005 = NA_real_,
            q995 = NA_real_
        )
        model_info <- list(
            method = paste0(
                "density/GMM mode gate: positive peak at ", round(density_gate$pos_peak, 2),
                "; negative: ", neg_gate_method
            ),
            components = component_df
        )
    } else if (dir == "right") {
        max_val_log <- max(vals_for_gate, na.rm = TRUE)
        med_val_log <- stats::median(vals_for_gate, na.rm = TRUE)
        is_cut_off <- (max_val_log - med_val_log) < 0.25
        if (is_cut_off) {
            lower_log <- as.numeric(stats::quantile(vals_for_gate, max(0, 0.5 - pct / 2), na.rm = TRUE))
            upper_log <- min(as.numeric(stats::quantile(vals_for_gate, min(1, 0.5 + pct / 2), na.rm = TRUE)), if (is.finite(gmm_upper)) gmm_upper else Inf, na.rm = TRUE)
        } else {
            lower_log <- med_val_log
            upper_log <- if (is.finite(gmm_upper)) gmm_upper else q_gate[["upper"]]
        }
    } else if (dir == "left") {
        lower_log <- q_gate[["lower"]]
        upper_log <- stats::median(vals_for_gate, na.rm = TRUE)
    } else {
        lower_log <- q_gate[["lower"]]
        upper_log <- min(q_gate[["upper"]], if (is.finite(gmm_upper)) gmm_upper else Inf, na.rm = TRUE)
    }

    retained_fraction <- mean(peak_vals >= 10^lower_log & peak_vals <= 10^upper_log, na.rm = TRUE)
    if (is.null(density_gate) && dir == "right" && (!is.finite(retained_fraction) || retained_fraction < 0.15)) {
        upper_log <- min(
            as.numeric(stats::quantile(fit_vals, 0.995, na.rm = TRUE)),
            q_gate[["upper"]],
            na.rm = TRUE
        )
    }
    if (!is.finite(lower_log) || !is.finite(upper_log) || upper_log <= lower_log) {
        lower_log <- q_gate[["lower"]]
        upper_log <- q_gate[["upper"]]
    }

    attr(vals_log, "neg_log_min") <- neg_log_min
    attr(vals_log, "neg_log_max") <- neg_log_max
    attr(vals_log, "negative_gate_present") <- negative_gate_present
    attr(vals_log, "positive_gate_present") <- positive_gate_present
    attr(vals_log, "gate_method") <- model_info$method
    attr(vals_log, "gmm_components") <- model_info$components
    attr(vals_log, "gate_type") <- "histogram"

    list(vals_log = vals_log, gate_min = 10^lower_log, gate_max = 10^upper_log)
}
.reference_histogram_posterior_boundary <- function(fit, left_k, right_k, left_mean, right_mean) {
    if (!is.finite(left_mean) || !is.finite(right_mean) || right_mean <= left_mean) {
        return(NA_real_)
    }
    grid <- seq(left_mean, right_mean, length.out = 512)
    pred <- tryCatch(stats::predict(fit, newdata = grid)$z, error = function(e) NULL)
    if (is.null(pred) || ncol(pred) < max(left_k, right_k)) {
        return(mean(c(left_mean, right_mean)))
    }
    cross <- which(pred[, left_k] - pred[, right_k] <= 0)
    if (length(cross) == 0) {
        return(mean(c(left_mean, right_mean)))
    }
    grid[min(cross)]
}

.reference_histogram_component_stats <- function(x, cls) {
    do.call(rbind, lapply(sort(unique(cls)), function(k) {
        xk <- x[cls == k]
        data.frame(
            component = k,
            n = length(xk),
            prop = length(xk) / length(x),
            mean = mean(xk, na.rm = TRUE),
            sd = stats::sd(xk, na.rm = TRUE),
            q005 = as.numeric(stats::quantile(xk, 0.005, na.rm = TRUE)),
            q995 = as.numeric(stats::quantile(xk, 0.995, na.rm = TRUE)),
            stringsAsFactors = FALSE
        )
    }))
}

.reference_histogram_nearest_trough <- function(d, trough_idx, peak, side) {
    candidates <- if (identical(side, "left")) trough_idx[trough_idx < peak] else trough_idx[trough_idx > peak]
    if (length(candidates) == 0) {
        return(if (identical(side, "left")) min(d$x, na.rm = TRUE) else max(d$x, na.rm = TRUE))
    }
    d$x[if (identical(side, "left")) max(candidates) else min(candidates)]
}

.reference_histogram_viability_density_gate <- function(d, peak_idx, trough_idx, vals_log) {
    sig_peaks <- peak_idx[
        d$y[peak_idx] >= max(d$y, na.rm = TRUE) * 0.03 &
            d$x[peak_idx] > 0.75
    ]
    sig_peaks <- sig_peaks[order(d$x[sig_peaks])]
    if (length(sig_peaks) < 2) return(NULL)

    neg_peak <- sig_peaks[length(sig_peaks) - 1L]
    pos_peak <- sig_peaks[length(sig_peaks)]
    boundary_candidates <- trough_idx[trough_idx > neg_peak & trough_idx < pos_peak]
    boundary <- if (length(boundary_candidates) > 0) {
        d$x[boundary_candidates[which.min(d$y[boundary_candidates])]]
    } else {
        mean(c(d$x[neg_peak], d$x[pos_peak]))
    }
    neg_left <- trough_idx[trough_idx < neg_peak]
    neg_log_min <- if (length(neg_left) > 0) {
        d$x[max(neg_left)]
    } else {
        as.numeric(stats::quantile(vals_log[vals_log <= boundary], 0.005, na.rm = TRUE))
    }
    pos_upper <- as.numeric(stats::quantile(vals_log[vals_log >= boundary], 0.999, na.rm = TRUE))
    valid <- is.finite(boundary) && is.finite(pos_upper) && pos_upper > boundary &&
        is.finite(neg_log_min) && boundary > neg_log_min
    if (!valid) return(NULL)

    list(
        gate = list(
            pos_lower = boundary,
            pos_upper = pos_upper,
            pos_peak = d$x[pos_peak],
            neg_peak = d$x[neg_peak],
            component_means = d$x[sig_peaks]
        ),
        neg_log_min = neg_log_min,
        neg_log_max = boundary,
        method = paste0(
            "viability gate: negative low mode at ", round(d$x[neg_peak], 2),
            "; positive high mode at ", round(d$x[pos_peak], 2)
        )
    )
}

.reference_histogram_negative_peak_group <- function(d, peak_idx, trough_idx, dominant_peak) {
    candidates <- peak_idx[
        peak_idx < dominant_peak &
            d$y[peak_idx] >= max(d$y, na.rm = TRUE) * 0.03 &
            d$x[peak_idx] > 0.75 &
            d$x[peak_idx] <= d$x[dominant_peak] - 0.75
    ]
    if (length(candidates) == 0) return(integer())

    group <- candidates[which.max(d$x[candidates])]
    previous <- rev(candidates[candidates < min(group)])
    for (peak in previous) {
        between <- trough_idx[trough_idx > peak & trough_idx < min(group)]
        if (length(between) == 0) break
        valley_y <- min(d$y[between], na.rm = TRUE)
        adjacent_y <- min(d$y[c(peak, min(group))], na.rm = TRUE)
        if (!is.finite(valley_y) || !is.finite(adjacent_y) || valley_y < adjacent_y * 0.60) break
        group <- c(peak, group)
    }
    sort(group)
}

.reference_histogram_bright_cell_density_gate <- function(d, peak_idx, trough_idx, vals_log, dominant_peak) {
    left <- trough_idx[trough_idx < dominant_peak]
    neg_group <- .reference_histogram_negative_peak_group(d, peak_idx, trough_idx, dominant_peak)
    neg_peak <- if (length(neg_group) > 0) max(neg_group) else NA_integer_
    neg_log_min <- min(vals_log, na.rm = TRUE)
    neg_log_max <- as.numeric(stats::quantile(vals_log, 0.15, na.rm = TRUE))
    method <- paste0("density dominant bright cell peak at ", round(d$x[dominant_peak], 2))

    if (!is.na(neg_peak)) {
        neg_left <- trough_idx[trough_idx < min(neg_group)]
        neg_right <- trough_idx[trough_idx > max(neg_group) & trough_idx < dominant_peak]
        neg_log_max <- if (length(neg_right) > 0) d$x[min(neg_right)] else d$x[max(left)]
        density_floor <- max(d$y[neg_group], na.rm = TRUE) * 0.10
        left_drop <- which(d$x < d$x[min(neg_group)] & d$y <= density_floor)
        neg_log_min <- if (length(left_drop) > 0) {
            d$x[max(left_drop)]
        } else if (length(neg_left) > 0) {
            d$x[max(neg_left)]
        } else {
            as.numeric(stats::quantile(vals_log[vals_log <= neg_log_max], 0.005, na.rm = TRUE))
        }
        method <- paste0("density middle negative cell mode at ", paste(round(d$x[neg_group], 2), collapse = "/"))
    }

    list(
        gate = list(
            pos_lower = if (length(left) > 0) d$x[max(left)] else as.numeric(stats::quantile(vals_log, 0.005, na.rm = TRUE)),
            pos_upper = as.numeric(stats::quantile(vals_log, 0.999, na.rm = TRUE)),
            pos_peak = d$x[dominant_peak],
            neg_peak = if (!is.na(neg_peak)) d$x[neg_peak] else NA_real_,
            component_means = d$x[peak_idx[d$y[peak_idx] >= max(d$y, na.rm = TRUE) * 0.005]]
        ),
        neg_log_min = neg_log_min,
        neg_log_max = neg_log_max,
        method = method
    )
}

.reference_histogram_cell_density_gate <- function(d, peak_idx, trough_idx, vals_log) {
    dominant_peak <- peak_idx[which.max(d$y[peak_idx])]
    left <- trough_idx[trough_idx < dominant_peak]
    right <- trough_idx[trough_idx > dominant_peak]
    bright_dominant <- d$x[dominant_peak] >= stats::median(vals_log, na.rm = TRUE) &&
        d$x[dominant_peak] >= max(vals_log, na.rm = TRUE) - 0.75
    if (bright_dominant) {
        return(.reference_histogram_bright_cell_density_gate(d, peak_idx, trough_idx, vals_log, dominant_peak))
    }

    neg_log_min <- if (length(left) > 0) d$x[max(left)] else as.numeric(stats::quantile(vals_log, 0.005, na.rm = TRUE))
    neg_log_max <- if (length(right) > 0) d$x[min(right)] else as.numeric(stats::quantile(vals_log, 0.85, na.rm = TRUE))
    pos_lower <- max(neg_log_max, as.numeric(stats::quantile(vals_log, 0.95, na.rm = TRUE)), na.rm = TRUE)
    pos_upper <- as.numeric(stats::quantile(vals_log, 0.999, na.rm = TRUE))
    if (!is.finite(pos_lower) || !is.finite(pos_upper) || pos_upper <= pos_lower) return(NULL)

    list(
        gate = list(
            pos_lower = pos_lower,
            pos_upper = pos_upper,
            pos_peak = d$x[peak_idx[which.max(d$x[peak_idx])]],
            neg_peak = d$x[dominant_peak],
            component_means = d$x[peak_idx[d$y[peak_idx] >= max(d$y, na.rm = TRUE) * 0.005]]
        ),
        neg_log_min = neg_log_min,
        neg_log_max = neg_log_max,
        method = paste0("density dominant negative peak at ", round(d$x[dominant_peak], 2))
    )
}

.reference_histogram_group_peaks <- function(d, sig_peaks) {
    if (length(sig_peaks) <= 1) return(sig_peaks)
    grouped <- list()
    current <- sig_peaks[1]
    for (idx in sig_peaks[-1]) {
        if ((d$x[idx] - d$x[current[length(current)]]) <= 0.35) {
            current <- c(current, idx)
        } else {
            grouped[[length(grouped) + 1L]] <- current
            current <- idx
        }
    }
    grouped[[length(grouped) + 1L]] <- current
    as.integer(vapply(grouped, function(group) group[which.max(d$y[group])], numeric(1)))
}

.reference_histogram_generic_density_gate <- function(d, peak_idx, trough_idx, vals_log, direction) {
    sig_peaks <- peak_idx[d$y[peak_idx] >= max(d$y, na.rm = TRUE) * 0.03 & d$x[peak_idx] > 0.75]
    sig_peaks <- .reference_histogram_group_peaks(d, sig_peaks)
    if (length(sig_peaks) == 0) return(NULL)

    sig_peaks <- sig_peaks[order(d$x[sig_peaks])]
    pos_peak <- sig_peaks[length(sig_peaks)]
    neg_peak <- if (length(sig_peaks) >= 2) sig_peaks[length(sig_peaks) - 1L] else NA_integer_
    left_trough <- .reference_histogram_nearest_trough(d, trough_idx, pos_peak, "left")
    right_trough <- .reference_histogram_nearest_trough(d, trough_idx, pos_peak, "right")
    pos_vals <- vals_log[vals_log >= left_trough]
    pos_upper <- as.numeric(stats::quantile(pos_vals[pos_vals >= d$x[pos_peak]], 0.999, na.rm = TRUE))
    if (is.finite(right_trough) && right_trough > d$x[pos_peak] && right_trough < max(vals_log, na.rm = TRUE)) {
        pos_upper <- min(pos_upper, right_trough, na.rm = TRUE)
    }
    cut_off <- max(vals_log, na.rm = TRUE) - d$x[pos_peak] < 0.25
    result <- list(
        gate = list(
            pos_lower = if (identical(direction, "right") && !cut_off) d$x[pos_peak] else left_trough,
            pos_upper = pos_upper,
            pos_peak = d$x[pos_peak],
            neg_peak = if (!is.na(neg_peak)) d$x[neg_peak] else NA_real_,
            component_means = d$x[sig_peaks]
        ),
        neg_log_min = min(vals_log, na.rm = TRUE),
        neg_log_max = as.numeric(stats::quantile(vals_log, 0.15, na.rm = TRUE)),
        method = "density lower-tail fallback"
    )
    if (!is.na(neg_peak)) {
        result$neg_log_min <- .reference_histogram_nearest_trough(d, trough_idx, neg_peak, "left")
        result$neg_log_max <- .reference_histogram_nearest_trough(d, trough_idx, neg_peak, "right")
        result$method <- paste0("density negative peak at ", round(d$x[neg_peak], 2))
    }
    result
}

.reference_histogram_density_gate <- function(vals_log, sample_type, is_viability, direction) {
    if (length(vals_log) < 80) return(NULL)
    d <- stats::density(vals_log, n = 2048)
    peak_idx <- which(diff(sign(diff(d$y))) == -2) + 1
    if (length(peak_idx) == 0) return(NULL)
    trough_idx <- which(diff(sign(diff(d$y))) == 2) + 1

    if (sample_type %in% "cells" && isTRUE(is_viability)) {
        result <- .reference_histogram_viability_density_gate(d, peak_idx, trough_idx, vals_log)
        if (!is.null(result)) return(result)
    }
    if (sample_type %in% "cells") {
        result <- .reference_histogram_cell_density_gate(d, peak_idx, trough_idx, vals_log)
        if (!is.null(result)) return(result)
    }
    .reference_histogram_generic_density_gate(d, peak_idx, trough_idx, vals_log, direction)
}

.reference_histogram_negative_gmm <- function(vals_log, defaults) {
    if (length(vals_log) < 80 || !requireNamespace("mclust", quietly = TRUE)) return(defaults)
    fit <- tryCatch(
        mclust::Mclust(vals_log, G = seq_len(min(6, max(1, floor(length(vals_log) / 40)))), verbose = FALSE),
        error = function(e) NULL
    )
    if (is.null(fit) || is.null(fit$classification)) return(defaults)

    component_stats <- .reference_histogram_component_stats(vals_log, fit$classification)
    defaults$components <- component_stats
    floor_component <- component_stats$mean <= 0.1 | component_stats$q995 <= 0.25
    high_artifact <- component_stats$prop < 0.02 & component_stats$mean > stats::median(vals_log, na.rm = TRUE)
    real_idx <- which(!floor_component & !high_artifact & component_stats$prop >= 0.02)
    if (length(real_idx) < 2) return(defaults)

    positive_row <- real_idx[which.max(component_stats$mean[real_idx])]
    candidates <- real_idx[component_stats$mean[real_idx] < component_stats$mean[positive_row]]
    if (length(candidates) == 0) return(defaults)
    negative_row <- candidates[which.max(component_stats$mean[candidates])]
    negative_k <- component_stats$component[negative_row]
    positive_k <- component_stats$component[positive_row]
    left_bound <- component_stats$q005[negative_row]
    left_rows <- which(component_stats$mean < component_stats$mean[negative_row])
    if (length(left_rows) > 0) {
        left_row <- left_rows[which.max(component_stats$mean[left_rows])]
        left_cross <- .reference_histogram_posterior_boundary(
            fit, component_stats$component[left_row], negative_k,
            component_stats$mean[left_row], component_stats$mean[negative_row]
        )
        if (is.finite(left_cross)) left_bound <- max(left_bound, left_cross)
    }
    right_bound <- component_stats$q995[negative_row]
    right_cross <- .reference_histogram_posterior_boundary(
        fit, negative_k, positive_k,
        component_stats$mean[negative_row], component_stats$mean[positive_row]
    )
    if (is.finite(right_cross)) right_bound <- min(right_bound, right_cross)
    if (!is.finite(left_bound) || !is.finite(right_bound) || right_bound <= left_bound) return(defaults)

    defaults$min <- left_bound
    defaults$max <- right_bound
    defaults$method <- paste0("GMM negative component ", negative_k, " left of positive component ", positive_k)
    defaults
}

.reference_histogram_bright_gmm <- function(fit_vals, negative) {
    result <- list(
        upper = NA_real_,
        info = list(method = "quantile fallback", components = negative$components)
    )
    if (length(fit_vals) < 80 || !requireNamespace("mclust", quietly = TRUE)) return(result)
    fit <- tryCatch(
        mclust::Mclust(fit_vals, G = seq_len(min(6, max(1, floor(length(fit_vals) / 40)))), verbose = FALSE),
        error = function(e) NULL
    )
    if (is.null(fit) || is.null(fit$classification)) return(result)

    component_stats <- .reference_histogram_component_stats(fit_vals, fit$classification)
    q995 <- as.numeric(stats::quantile(fit_vals, 0.995, na.rm = TRUE))
    q999 <- as.numeric(stats::quantile(fit_vals, 0.999, na.rm = TRUE))
    high_artifact <- (component_stats$prop < 0.02 & component_stats$mean > stats::median(fit_vals, na.rm = TRUE)) |
        (component_stats$prop < 0.05 & component_stats$mean >= q995) |
        (is.finite(component_stats$sd) & component_stats$sd < 0.005 & component_stats$mean >= q999)
    real_idx <- which(!high_artifact & component_stats$prop >= 0.02)
    if (length(real_idx) == 0) real_idx <- which(!high_artifact)
    if (length(real_idx) == 0) return(result)

    bright_row <- real_idx[which.max(component_stats$mean[real_idx])]
    artifact_rows <- which(component_stats$mean > component_stats$mean[bright_row] & high_artifact)
    if (length(artifact_rows) > 0) result$upper <- q995
    result$info <- list(
        method = paste0(
            "GMM bright component ", component_stats$component[bright_row],
            if (length(artifact_rows) > 0) "; high-end artifact components excluded" else "; no separated high-end artifact detected",
            "; negative: ", negative$method
        ),
        components = if (!is.null(negative$components)) negative$components else component_stats
    )
    result
}
