#' Calculate Negative Population Spread (NPS)
#' 
#' Quantifies the spreading of negative populations after unmixing. 
#' High spread indicates poor unmixing or high spillover noise.
#' 
#' @param data Unmixed data frame
#' @param markers Vector of marker names to analyze
#' @return A data frame with NPS values (MAD) per marker per file
#' @examples
#' demo_df <- data.frame(
#'   File = rep(c("A", "B"), each = 50),
#'   FITC = rnorm(100, sd = 0.3),
#'   PE = rnorm(100, sd = 0.4)
#' )
#' nps <- calculate_nps(demo_df)
#' head(nps)
#' @export
calculate_nps <- function(data, markers = NULL) {
    if (is.null(markers)) {
        exclude <- .get_result_metadata_columns(colnames(data))
        markers <- setdiff(colnames(data), exclude)
        markers <- markers[!grepl("^AF($|_)", markers, ignore.case = TRUE)]
    }
    
    # For each marker, isolate the negative population (intensity < threshold)
    # We use a robust estimate (Median Absolute Deviation) of the spread
    nps_results <- data |> 
        dplyr::group_by(File) |> 
        dplyr::summarize(dplyr::across(dplyr::all_of(markers), function(x) {
            # Heuristic: negative population is around 0. 
            # We take values between -2SD and +2SD to avoid true positives
            # But simpler: just take the MAD of all values < quantile(0.2)
            neg_subset <- x[x < quantile(x, 0.2)]
            stats::mad(neg_subset, na.rm = TRUE)
        })) |> 
        tidyr::pivot_longer(cols = -File, names_to = "Marker", values_to = "NPS")
    
    return(nps_results)
}

#' Plot Negative Population Spread
#' @param nps_results Output of [calculate_nps()].
#' @param output_file Path to save the plot (set `NULL` to skip saving).
#' @param width Plot width.
#' @return A `ggplot` object.
#' @examples
#' demo_df <- data.frame(
#'   File = rep(c("sample_a", "sample_b"), each = 60),
#'   FITC = c(rnorm(60, 0.2, 0.2), rnorm(60, 0.1, 0.15)),
#'   PE = c(rnorm(60, 0.3, 0.25), rnorm(60, 0.2, 0.2))
#' )
#' nps <- calculate_nps(demo_df)
#' p <- plot_nps(nps, output_file = NULL)
#' print(p)
#' @export
plot_nps <- function(nps_results, output_file = NULL, width = 200) {
    n_files <- length(unique(nps_results$File))
    
    p <- ggplot2::ggplot(nps_results, ggplot2::aes(Marker, NPS))
    
    if (n_files > 1) {
        p <- p + 
            ggplot2::geom_boxplot(outlier.shape = NA, fill = "grey95", color = "grey40", width = 0.5) +
            ggplot2::geom_jitter(ggplot2::aes(color = File), width = 0.15, height = 0, alpha = 0.8, size = 1.2)
    } else {
        p <- p + 
            ggplot2::geom_bar(stat = "identity", fill = "#E06666", width = 0.6)
    }
    
    p <- p +
        ggplot2::labs(
            title = "Negative Population Spread (Unmixing Noise Floor)",
            subtitle = "Lower is better. High spread can worsen separation between pos. and neg. populations.\nCell-based control colors will show higher MAD.\nUse colors with high spread on abundant markers whenever possible.",
            y = "Spread (MAD)",
            x = "Unmixed Marker"
        ) +
        ggplot2::theme_minimal(base_size = 13.75) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            plot.subtitle = ggplot2::element_text(size = 13.2, lineheight = 1.1)
        )
        
    if (n_files == 1) {
        p <- p + ggplot2::theme(legend.position = "none")
    }
    
    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, p, width = width, height = 120, units = "mm", dpi = 300)
    }
    return(p)
}

#' Calculate Cosine Similarity Matrix
#'
#' @param M Reference matrix (Markers x Detectors)
#' @return A square matrix representing pairwise cosine similarities.
#' @export
calculate_similarity_matrix <- function(M) {
    M <- .as_reference_matrix(M, "M")
    norms <- sqrt(rowSums(M^2, na.rm = TRUE))
    norms[norms == 0] <- 1e-6
    M_norm <- M / norms
    sim_mat <- M_norm %*% t(M_norm)
    sim_mat[sim_mat > 1] <- 1
    sim_mat[sim_mat < 0] <- 0
    return(sim_mat)
}

#' Plot Cosine Similarity Matrix
#'
#' @param similarity_matrix Matrix returned by calculate_similarity_matrix
#' @param output_file Optional path to save the plot. Set `NULL` to return the plot.
#' @param width Width of plot in mm
#' @param height Height of plot in mm
#' @return A `ggplot` object.
#' @export
plot_similarity_matrix <- function(similarity_matrix, output_file = NULL, width = 180, height = 160) {
    sim_tri <- similarity_matrix
    is_square_same_markers <- nrow(sim_tri) == ncol(sim_tri) &&
        identical(rownames(sim_tri), colnames(sim_tri))
    marker_order <- if (is_square_same_markers) rownames(sim_tri) else colnames(sim_tri)
    if (is_square_same_markers) {
        sim_tri <- sim_tri[marker_order, marker_order, drop = FALSE]
        sim_tri[upper.tri(sim_tri, diag = FALSE)] <- NA
    }
    
    long <- as.data.frame(sim_tri, check.names = FALSE)
    long$Marker1 <- rownames(sim_tri)
    long <- tidyr::pivot_longer(long, cols = -Marker1, names_to = "Marker2", values_to = "Similarity")
    long <- long[!is.na(long$Similarity), ]
    
    row_markers <- if (is_square_same_markers) marker_order else rownames(sim_tri)
    col_markers <- if (is_square_same_markers) marker_order else colnames(sim_tri)
    long$Marker1 <- factor(long$Marker1, levels = rev(row_markers))
    long$Marker2 <- factor(long$Marker2, levels = col_markers)
    long$is_diagonal <- as.character(long$Marker1) == as.character(long$Marker2)
    diag_long <- long[long$is_diagonal, , drop = FALSE]
    offdiag_long <- long[!long$is_diagonal, , drop = FALSE]
    offdiag_long$Similarity <- pmax(0, pmin(1, offdiag_long$Similarity))
    offdiag_long$SimilarityFill <- pmin(offdiag_long$Similarity, 0.99)
    
    n_markers <- max(length(row_markers), length(col_markers))
    text_size <- max(2.4, min(4.8, 36 / max(1, n_markers)))
    
    p <- ggplot2::ggplot() +
        ggplot2::geom_tile(
            data = offdiag_long,
            ggplot2::aes(Marker2, Marker1, fill = SimilarityFill),
            color = "white",
            linewidth = 0.1
        ) +
        ggplot2::geom_tile(
            data = diag_long,
            ggplot2::aes(Marker2, Marker1),
            fill = "#E6E8EB",
            color = "white",
            linewidth = 0.1
        ) +
        ggplot2::scale_fill_gradientn(
            colors = c("#FFFFFF", "#FEE5D9", "#FCAE91", "#FB6A4A", "#CB181D"),
            values = c(0, 0.5, 0.75, 0.9, 1.0),
            limits = c(0, 0.99),
            name = "Similarity"
        ) +
        ggplot2::geom_text(
            data = offdiag_long,
            ggplot2::aes(
                Marker2,
                Marker1,
                label = sprintf("%.2f", Similarity)
            ),
            size = text_size,
            color = ifelse(offdiag_long$Similarity > 0.8, "white", "black"),
            show.legend = FALSE
        ) +
        ggplot2::labs(
            title = "Fluorophore Spectral Similarity",
            subtitle = "Cosine similarity of reference signatures (0 = orthogonal, 1 = identical).\nHigh similarity can indicate potential spillover.",
            x = NULL, y = NULL
        ) +
        ggplot2::theme_minimal(base_size = 13.75) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            panel.grid = ggplot2::element_blank(),
            plot.subtitle = ggplot2::element_text(size = 13.2, lineheight = 1.1)
        )
        
    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, p, width = width, height = height, units = "mm", dpi = 300)
    }
    return(p)
}

.residual_detector_laser_group <- function(detectors) {
    key <- sub("-A$", "", as.character(detectors), ignore.case = TRUE)
    dplyr::case_when(
        grepl("^UV", key, ignore.case = TRUE) ~ "UV",
        grepl("^YG", key, ignore.case = TRUE) ~ "YG",
        grepl("^V", key, ignore.case = TRUE) ~ "Violet",
        grepl("^B", key, ignore.case = TRUE) ~ "Blue",
        grepl("^R", key, ignore.case = TRUE) ~ "Red",
        TRUE ~ "Other"
    )
}

.residual_detector_channel_order <- function(detectors) {
    detectors <- as.character(detectors)
    key <- toupper(trimws(sub("-A$", "", detectors, ignore.case = TRUE)))
    has_prefix <- grepl("^(UV|YG|V|B|Y|G|R)[0-9]+", key, perl = TRUE)
    prefix <- rep("", length(key))
    prefix[has_prefix] <- sub("^(UV|YG|V|B|Y|G|R)([0-9]+).*$", "\\1", key[has_prefix], perl = TRUE)

    channel <- rep(.Machine$integer.max, length(key))
    channel[has_prefix] <- suppressWarnings(as.integer(
        sub("^(UV|YG|V|B|Y|G|R)([0-9]+).*$", "\\2", key[has_prefix], perl = TRUE)
    ))
    missing_channel <- !is.finite(channel)
    if (any(missing_channel)) {
        channel[missing_channel] <- suppressWarnings(as.integer(
            sub("^[^0-9]*([0-9]+).*$", "\\1", key[missing_channel], perl = TRUE)
        ))
    }
    channel[!is.finite(channel)] <- .Machine$integer.max

    laser_rank <- dplyr::case_when(
        prefix == "UV" ~ 1L,
        prefix == "V" ~ 2L,
        prefix == "B" ~ 3L,
        prefix %in% c("YG", "Y", "G") ~ 4L,
        prefix == "R" ~ 5L,
        TRUE ~ 99L
    )
    order(laser_rank, channel, key, detectors, na.last = TRUE)
}

.collect_report_residual_matrix <- function(results, detector_names = NULL) {
    if (!is.list(results) || length(results) == 0) {
        return(NULL)
    }

    mats <- lapply(results, function(res_obj) {
        if (!is.list(res_obj) || is.null(res_obj$residuals)) return(NULL)
        R <- as.matrix(res_obj$residuals)
        if (nrow(R) == 0 || ncol(R) == 0) return(NULL)
        R
    })
    mats <- mats[!vapply(mats, is.null, logical(1))]
    if (length(mats) == 0) {
        return(NULL)
    }

    if (is.null(detector_names)) {
        detector_names <- unique(unlist(lapply(mats, colnames), use.names = FALSE))
    }
    detector_names <- detector_names[nzchar(detector_names)]
    if (length(detector_names) == 0) {
        return(NULL)
    }

    mats <- lapply(mats, function(R) {
        common <- intersect(detector_names, colnames(R))
        if (length(common) == 0) return(NULL)
        out <- matrix(NA_real_, nrow = nrow(R), ncol = length(detector_names), dimnames = list(NULL, detector_names))
        out[, common] <- R[, common, drop = FALSE]
        out
    })
    mats <- mats[!vapply(mats, is.null, logical(1))]
    if (length(mats) == 0) {
        return(NULL)
    }

    do.call(rbind, mats)
}

.resolve_residual_metric_method <- function(results, unmixing_method = NULL) {
    method <- unmixing_method
    if (is.null(method)) {
        method <- attr(results, "method")
    }
    if (is.null(method) && is.list(results) && length(results) > 0) {
        method <- attr(results[[1]], "method")
    }
    if (is.null(method)) {
        return("OLS")
    }
    .normalize_unmix_method(method, choices = c("AutoSpectral", "OLS", "NNLS", "WLS", "RWLS"))
}

.resolve_residual_metric_matrix <- function(results, M = NULL) {
    if (!is.null(M)) {
        return(.as_reference_matrix(M, "M"))
    }
    matrix_attr <- attr(results, "reference_matrix")
    if (!is.null(matrix_attr)) {
        return(.as_reference_matrix(matrix_attr, "M"))
    }
    if (is.list(results) && length(results) > 0) {
        matrix_attr <- attr(results[[1]], "reference_matrix")
        if (!is.null(matrix_attr)) {
            return(.as_reference_matrix(matrix_attr, "M"))
        }
    }
    NULL
}

.reconstruct_residual_raw_signal <- function(res_obj, M) {
    if (!is.list(res_obj) || is.null(res_obj$data) || is.null(res_obj$residuals) || is.null(M)) {
        return(NULL)
    }
    residuals <- as.matrix(res_obj$residuals)
    markers <- intersect(rownames(M), colnames(res_obj$data))
    detectors <- intersect(colnames(M), colnames(residuals))
    if (length(markers) == 0 || length(detectors) == 0) {
        return(NULL)
    }
    fitted <- as.matrix(res_obj$data[, markers, drop = FALSE]) %*% M[markers, detectors, drop = FALSE]
    residuals <- residuals[, detectors, drop = FALSE]
    fitted + residuals
}

.residual_metric_weights <- function(res_obj, M, unmixing_method = NULL) {
    method <- .resolve_residual_metric_method(list(res_obj), unmixing_method = unmixing_method)
    if (!(method %in% c("WLS", "RWLS")) || is.null(M)) {
        return(NULL)
    }
    Y <- .reconstruct_residual_raw_signal(res_obj, M)
    if (is.null(Y) || nrow(Y) == 0 || ncol(Y) == 0) {
        return(NULL)
    }
    params <- .resolve_wls_noise_parameters(M[, colnames(Y), drop = FALSE])
    weights <- t(apply(Y, 1, .wls_event_weights,
        noise_floor = params$noise_floor,
        signal_scale = params$signal_scale,
        max_weight_ratio = params$max_weight_ratio
    ))
    colnames(weights) <- colnames(Y)
    weights
}

.collect_report_residual_metric_matrix <- function(results, M = NULL, detector_names = NULL, unmixing_method = NULL) {
    if (!is.list(results) || length(results) == 0) {
        return(NULL)
    }
    M <- .resolve_residual_metric_matrix(results, M = M)
    method <- .resolve_residual_metric_method(results, unmixing_method = unmixing_method)
    mats <- lapply(results, function(res_obj) {
        if (!is.list(res_obj) || is.null(res_obj$residuals)) return(NULL)
        R <- as.matrix(res_obj$residuals)
        if (nrow(R) == 0 || ncol(R) == 0) return(NULL)
        if (method %in% c("WLS", "RWLS") && !is.null(M)) {
            weights <- .residual_metric_weights(res_obj, M = M, unmixing_method = method)
            if (!is.null(weights)) {
                common <- intersect(colnames(R), colnames(weights))
                R[, common] <- sqrt(R[, common, drop = FALSE]^2 * weights[, common, drop = FALSE])
            } else {
                R <- abs(R)
            }
        } else {
            R <- abs(R)
        }
        R
    })
    mats <- mats[!vapply(mats, is.null, logical(1))]
    if (length(mats) == 0) {
        return(NULL)
    }
    if (is.null(detector_names)) {
        detector_names <- unique(unlist(lapply(mats, colnames), use.names = FALSE))
    }
    detector_names <- detector_names[nzchar(detector_names)]
    if (length(detector_names) == 0) {
        return(NULL)
    }
    mats <- lapply(mats, function(R) {
        common <- intersect(detector_names, colnames(R))
        if (length(common) == 0) return(NULL)
        out <- matrix(NA_real_, nrow = nrow(R), ncol = length(detector_names), dimnames = list(NULL, detector_names))
        out[, common] <- R[, common, drop = FALSE]
        out
    })
    mats <- mats[!vapply(mats, is.null, logical(1))]
    if (length(mats) == 0) {
        return(NULL)
    }
    do.call(rbind, mats)
}

plot_detector_rms_residuals <- function(results, M = NULL, pd = NULL, output_file = NULL, width = 250, height = 120, unmixing_method = NULL) {
    detector_names <- if (!is.null(M)) colnames(.as_reference_matrix(M, "M")) else NULL
    method <- .resolve_residual_metric_method(results, unmixing_method = unmixing_method)
    residuals <- .collect_report_residual_metric_matrix(results, M = M, detector_names = detector_names, unmixing_method = unmixing_method)
    if (is.null(residuals) || nrow(residuals) == 0) {
        warning("No residuals available to plot detector RMS residuals.")
        return(NULL)
    }

    detector_names <- colnames(residuals)
    levels_sorted <- detector_names[.residual_detector_channel_order(detector_names)]
    if (length(levels_sorted) == 0) {
        levels_sorted <- detector_names
    }

    rms <- sqrt(colMeans(residuals[, levels_sorted, drop = FALSE]^2, na.rm = TRUE))
    labels <- sub("-A$", "", levels_sorted, ignore.case = TRUE)
    laser <- .residual_detector_laser_group(levels_sorted)
    plot_df <- data.frame(
        Detector = levels_sorted,
        DetectorIndex = seq_along(levels_sorted),
        Label = labels,
        Laser = factor(laser, levels = c("UV", "Violet", "Blue", "YG", "Red", "Other")),
        RMS = as.numeric(rms),
        stringsAsFactors = FALSE
    )
    plot_df <- plot_df[is.finite(plot_df$RMS), , drop = FALSE]
    if (nrow(plot_df) == 0) {
        warning("No finite detector RMS residuals available.")
        return(NULL)
    }

    y_label <- if (method %in% c("WLS", "RWLS")) "WLS-weighted RMS residual" else "RMS residual (raw detector units)"
    separators <- which(plot_df$Laser[-1] != plot_df$Laser[-nrow(plot_df)]) + 0.5
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(DetectorIndex, RMS)) +
        ggplot2::geom_col(width = 0.82, fill = "grey35", color = "grey25", linewidth = 0.15) +
        ggplot2::geom_vline(xintercept = separators, color = "grey45", linewidth = 0.25) +
        ggplot2::scale_x_continuous(
            breaks = plot_df$DetectorIndex,
            labels = plot_df$Label,
            expand = ggplot2::expansion(mult = c(0.005, 0.005))
        ) +
        ggplot2::scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
        ggplot2::labs(
            title = "Reconstruction Error per Detector",
            subtitle = "Root Mean Square (RMS) of reconstruction residuals per detector. Lower is better.",
            x = "Detector",
            y = y_label
        ) +
        ggplot2::theme_minimal(base_size = 13.75) +
        ggplot2::theme(
            legend.position = "none",
            axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.2),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.minor.x = ggplot2::element_blank(),
            plot.subtitle = ggplot2::element_text(size = 12.5, lineheight = 1.1)
        )

    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, p, width = width, height = height, units = "mm", dpi = 300)
    }
    p
}

#' Plot Overall Detector Reconstruction Error Across Samples
#'
#' @param results List of unmixed results from unmix_samples()
#' @param M Optional reference matrix to reconstruct raw intensities and compute relative error
#' @param output_file Optional path to save the plot. Set `NULL` to return the plot.
#' @param width Width of plot in mm
#' @param height Height of plot in mm
#' @param unmixing_method Optional unmixing method used to score residuals.
#'   When `"WLS"` or `"RWLS"`, residual RMS values use the WLS detector-noise
#'   weights; otherwise raw detector residuals are used.
#' @return A `ggplot` object.
#' @export
plot_sample_rms_residuals <- function(results, M = NULL, output_file = NULL, width = 225, height = 125, unmixing_method = NULL) {
    if (!is.list(results) || length(results) == 0) {
        stop("results must be a non-empty list of unmixed results.")
    }
    M_mat <- .resolve_residual_metric_matrix(results, M = M)
    method <- .resolve_residual_metric_method(results, unmixing_method = unmixing_method)
    y_label <- if (method %in% c("WLS", "RWLS")) "WLS-weighted RMS residual per cell" else "RMS residual per cell"
    
    sample_dfs <- list()
    for (sn in names(results)) {
        res_obj <- results[[sn]]
        if (is.list(res_obj) && !is.null(res_obj$residuals)) {
            residuals <- as.matrix(res_obj$residuals)
            weights <- if (method %in% c("WLS", "RWLS")) .residual_metric_weights(res_obj, M = M_mat, unmixing_method = method) else NULL
            if (!is.null(weights)) {
                common <- intersect(colnames(residuals), colnames(weights))
                if (length(common) > 0) {
                    metric_sq <- residuals[, common, drop = FALSE]^2 * weights[, common, drop = FALSE]
                    rms <- sqrt(rowMeans(metric_sq, na.rm = TRUE))
                } else {
                    rms <- sqrt(rowMeans(residuals^2, na.rm = TRUE))
                }
            } else {
                rms <- sqrt(rowMeans(residuals^2, na.rm = TRUE))
            }
            median_rms <- stats::median(rms, na.rm = TRUE)
            
            # Estimate peak raw signal for relative error context
            if (!is.null(M_mat)) {
                markers <- intersect(rownames(M_mat), colnames(res_obj$data))
                detectors <- intersect(colnames(M_mat), colnames(residuals))
                Fitted <- as.matrix(res_obj$data[, markers, drop = FALSE]) %*% M_mat[markers, detectors, drop = FALSE]
                Y <- Fitted + residuals[, detectors, drop = FALSE]
                peak_signal <- stats::quantile(Y, 0.995, na.rm = TRUE)
            } else {
                expr_cols <- setdiff(colnames(res_obj$data), c("File", "Time"))
                expr_cols <- expr_cols[!grepl("^FSC|^SSC", expr_cols)]
                expr_mat <- as.matrix(res_obj$data[, expr_cols, drop = FALSE])
                peak_signal <- stats::quantile(expr_mat, 0.995, na.rm = TRUE)
            }
            
            error_ratio <- (median_rms / max(peak_signal, 100)) * 100
            label_text <- sprintf("%s\n(Med: %.0f, Err: %.1f%%)", sn, median_rms, error_ratio)
            
            sample_dfs[[sn]] <- data.frame(
                Sample = label_text,
                RMS = rms,
                stringsAsFactors = FALSE
            )
        }
    }
    
    if (length(sample_dfs) == 0) {
        warning("No residuals available to plot RMS residuals.")
        return(NULL)
    }
    
    df <- do.call(rbind, sample_dfs)
    
    p <- ggplot2::ggplot(df, ggplot2::aes(Sample, RMS, fill = Sample)) +
        ggplot2::geom_violin(alpha = 0.5, color = "grey60", scale = "width") +
        ggplot2::geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, color = "grey30") +
        ggplot2::scale_y_continuous(
            trans = scales::pseudo_log_trans(base = 10, sigma = 100),
            breaks = c(0, 100, 500, 1000, 5000, 10000, 20000)
        ) +
        ggplot2::labs(
            title = "Reconstruction Error per Sample",
            subtitle = "Cell-level Root Mean Square (RMS) residuals norm. against the peak signal. Lower is better.",
            x = "Sample",
            y = y_label
        ) +
        ggplot2::theme_minimal(base_size = 13.75) +
        ggplot2::theme(
            legend.position = "none",
            axis.text.x = ggplot2::element_text(size = 9.375, angle = 45, hjust = 1),
            plot.subtitle = ggplot2::element_text(size = 13.2, lineheight = 1.1)
        )
        
    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, p, width = width, height = height, units = "mm", dpi = 300)
    }
    return(p)
}
