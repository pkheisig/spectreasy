# Internal helpers for QC report generation.
.draw_qc_report_plot_page <- function(p, square = FALSE, height_ratio = 1.0, width_ratio = 1.0, newpage = TRUE) {
    if (is.null(p)) return(invisible(NULL))
    if (isTRUE(newpage)) {
        grid::grid.newpage()
    }
    grob <- ggplot2::ggplotGrob(p)
    if (!isTRUE(square) && height_ratio == 1.0 && width_ratio == 1.0) {
        grid::grid.draw(grob)
        return(invisible(NULL))
    }

    vp_h <- if (isTRUE(square)) 0.99 else height_ratio
    vp_w <- if (isTRUE(square)) 0.78 else (0.95 * width_ratio)

    grid::grid.draw(
        grid::editGrob(
            grob,
            vp = grid::viewport(
                x = 0.5,
                y = 0.5,
                width = vp_w,
                height = vp_h
            )
        )
    )
    invisible(NULL)
}

.extract_qc_report_legend <- function(p) {
    grob <- ggplot2::ggplotGrob(p)
    grob_names <- vapply(grob$grobs, function(x) {
        if (is.null(x$name)) "" else x$name
    }, character(1))
    guide_idx <- which(grepl("^guide-box", grob_names))
    if (length(guide_idx) == 0) {
        return(grid::nullGrob())
    }
    grob$grobs[[guide_idx[1]]]
}

.draw_qc_report_spectra_page <- function(p, newpage = TRUE) {
    if (is.null(p)) return(invisible(NULL))

    n_labels <- 0L
    if (is.data.frame(p$data) && "Fluorophore" %in% colnames(p$data)) {
        n_labels <- length(unique(as.character(p$data$Fluorophore)))
    }

    legend_text_size <- if (n_labels > 34L) {
        5.8
    } else if (n_labels > 28L) {
        6.5
    } else if (n_labels > 22L) {
        7.5
    } else if (n_labels > 18L) {
        8.5
    } else {
        10
    }
    legend_title_size <- if (n_labels > 28L) 12 else if (n_labels > 22L) 14 else 16
    legend_key_size <- if (n_labels > 34L) {
        2.6
    } else if (n_labels > 28L) {
        3.0
    } else if (n_labels > 22L) {
        3.4
    } else if (n_labels > 18L) {
        3.8
    } else {
        4.3
    }
    legend_spacing <- if (n_labels > 22L) 0.1 else 0.35

    legend_plot <- p +
        ggplot2::guides(color = ggplot2::guide_legend(ncol = 1, byrow = FALSE, title.position = "top")) +
        ggplot2::theme(
            plot.margin = ggplot2::margin(0, 0, 0, 0),
            axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 3)),
            axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 3)),
            plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 4)),
            legend.position = "right",
            legend.title = ggplot2::element_text(size = legend_title_size),
            legend.text = ggplot2::element_text(size = legend_text_size),
            legend.key.size = grid::unit(legend_key_size, "mm"),
            legend.spacing.y = grid::unit(legend_spacing, "mm"),
            legend.margin = ggplot2::margin(0, 0, 0, 0),
            legend.box.margin = ggplot2::margin(0, 0, 0, 0)
        )
    plot_grob <- ggplot2::ggplotGrob(
        p +
            ggplot2::theme(
                legend.position = "none",
                plot.margin = ggplot2::margin(0, 0, 0, 0),
                axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 3)),
                axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 3)),
                plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 4))
            )
    )
    legend_grob <- .extract_qc_report_legend(legend_plot)

    if (isTRUE(newpage)) {
        grid::grid.newpage()
    }
    grid::grid.draw(
        grid::editGrob(
            plot_grob,
            vp = grid::viewport(
                x = 0.42,
                y = 0.49,
                width = 0.82,
                height = 0.66
            )
        )
    )
    grid::grid.draw(
        grid::editGrob(
            legend_grob,
            vp = grid::viewport(
                x = 0.915,
                y = 0.49,
                width = 0.16,
                height = 0.90
            )
        )
    )
    invisible(NULL)
}

.prepare_qc_report_png_dir <- function(qc_plot_dir = NULL, save_qc_pngs = FALSE, output_file = NULL) {
    if (!isTRUE(save_qc_pngs)) {
        return(NULL)
    }
    if (is.null(qc_plot_dir) || !nzchar(trimws(as.character(qc_plot_dir)[1]))) {
        report_dir <- if (!is.null(output_file)) dirname(output_file) else "."
        qc_plot_dir <- file.path(report_dir, "qc_samples_plots")
    }
    qc_plot_dir <- normalizePath(qc_plot_dir, mustWork = FALSE)
    dir.create(qc_plot_dir, showWarnings = FALSE, recursive = TRUE)
    qc_plot_dir
}

.save_qc_report_png <- function(p, qc_plot_dir, filename, width = 11, height = 8.5) {
    if (is.null(p) || is.null(qc_plot_dir)) {
        return(invisible(NULL))
    }
    .with_known_qc_plot_warnings_suppressed(
        ggplot2::ggsave(
            filename = file.path(qc_plot_dir, filename),
            plot = p,
            width = width,
            height = height,
            units = "in",
            dpi = 200
        )
    )
    invisible(file.path(qc_plot_dir, filename))
}

.is_known_qc_plot_warning <- function(message) {
    grepl("Binning grid too coarse for current \\(small\\) bandwidth", message) ||
        (
            grepl("Removed [0-9]+ rows? containing missing values or values outside the scale range", message) &&
                grepl("`geom_(tile|point)\\(\\)`", message)
        )
}

.with_known_qc_plot_warnings_suppressed <- function(expr) {
    withCallingHandlers(
        expr,
        warning = function(w) {
            if (.is_known_qc_plot_warning(conditionMessage(w))) {
                invokeRestart("muffleWarning")
            }
        }
    )
}

.write_qc_report_csv <- function(x, path, row_id = NULL) {
    if (is.null(path) || length(path) == 0 || is.na(path[1]) || !nzchar(trimws(as.character(path[1])))) {
        return(invisible(NULL))
    }
    out_dir <- dirname(path)
    if (!is.na(out_dir) && nzchar(out_dir) && out_dir != ".") {
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    }
    x <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
    if (!is.null(row_id)) {
        row_values <- rownames(x)
        if (is.null(row_values)) {
            row_values <- seq_len(nrow(x))
        }
        x[[row_id]] <- row_values
        x <- x[, c(row_id, setdiff(colnames(x), row_id)), drop = FALSE]
    }
    utils::write.csv(x, path, row.names = FALSE, quote = TRUE)
    invisible(path)
}

.write_qc_report_matrix_metric <- function(mat, path, row_id = "marker") {
    if (is.null(path) || length(path) == 0 || is.null(mat) || length(mat) == 0) {
        return(invisible(NULL))
    }
    .write_qc_report_csv(as.data.frame(mat, check.names = FALSE), path, row_id = row_id)
}

.prepare_qc_report_metrics_dir <- function(qc_metrics_dir = NULL) {
    if (!is.null(qc_metrics_dir) && length(qc_metrics_dir) > 0 &&
        !is.na(qc_metrics_dir[1]) && nzchar(trimws(as.character(qc_metrics_dir)[1]))) {
        qc_metrics_dir <- as.character(qc_metrics_dir)[1]
        dir.create(qc_metrics_dir, showWarnings = FALSE, recursive = TRUE)
        obsolete_files <- file.path(
            qc_metrics_dir,
            c(
                "sample_qc_summary.csv",
                "spectral_spread_matrix.csv",
                "directional_spread_score.csv",
                "overall_detector_reconstruction_error_per_sample.csv"
            )
        )
        unlink(obsolete_files[file.exists(obsolete_files)], force = TRUE)
        return(qc_metrics_dir)
    }
    NULL
}

.compute_qc_report_detector_rms <- function(res_list, M = NULL, pd = NULL, unmixing_method = NULL) {
    detector_names <- if (!is.null(M)) colnames(.as_reference_matrix(M, "M")) else NULL
    method <- .resolve_residual_metric_method(res_list, unmixing_method = unmixing_method)
    metric_name <- if (.uses_wls_residual_metric(method)) "wls_weighted_rms" else "raw_rms"
    residuals <- .collect_report_residual_metric_matrix(
        res_list,
        M = M,
        detector_names = detector_names,
        unmixing_method = unmixing_method
    )
    if (is.null(residuals) || nrow(residuals) == 0) {
        return(NULL)
    }
    detector_names <- colnames(residuals)
    levels_sorted <- detector_names[.residual_detector_channel_order(detector_names)]
    if (length(levels_sorted) == 0) {
        levels_sorted <- detector_names
    }
    rms <- sqrt(colMeans(residuals[, levels_sorted, drop = FALSE]^2, na.rm = TRUE))
    laser <- .residual_detector_laser_group(levels_sorted)
    out <- data.frame(
        detector = levels_sorted,
        label = sub("-A$", "", levels_sorted, ignore.case = TRUE),
        laser = laser,
        rms_residual = as.numeric(rms),
        residual_metric = metric_name,
        unmixing_method = method,
        stringsAsFactors = FALSE
    )
    out[is.finite(out$rms_residual), , drop = FALSE]
}

.compute_qc_report_sample_rms <- function(res_list, M = NULL, unmixing_method = NULL) {
    if (!is.list(res_list) || length(res_list) == 0) {
        return(NULL)
    }
    M_mat <- .resolve_residual_metric_matrix(res_list, M = M)
    method <- .resolve_residual_metric_method(res_list, unmixing_method = unmixing_method)
    metric_name <- if (.uses_wls_residual_metric(method)) "wls_weighted_rms" else "raw_rms"
    rows <- lapply(names(res_list), function(sn) {
        res_obj <- res_list[[sn]]
        if (!is.list(res_obj) || is.null(res_obj$residuals)) {
            return(NULL)
        }
        residuals <- as.matrix(res_obj$residuals)
        weights <- if (.uses_wls_residual_metric(method)) .residual_metric_weights(res_obj, M = M_mat, unmixing_method = method) else NULL
        if (!is.null(weights)) {
            common <- intersect(colnames(residuals), colnames(weights))
            if (length(common) > 0) {
                rms <- sqrt(rowMeans(residuals[, common, drop = FALSE]^2 * weights[, common, drop = FALSE], na.rm = TRUE))
            } else {
                rms <- sqrt(rowMeans(residuals^2, na.rm = TRUE))
            }
        } else {
            rms <- sqrt(rowMeans(residuals^2, na.rm = TRUE))
        }
        if (length(rms) == 0 || all(!is.finite(rms))) {
            return(NULL)
        }
        if (!is.null(M_mat)) {
            markers <- intersect(rownames(M_mat), colnames(res_obj$data))
            detectors <- intersect(colnames(M_mat), colnames(residuals))
            Fitted <- as.matrix(res_obj$data[, markers, drop = FALSE]) %*% M_mat[markers, detectors, drop = FALSE]
            Y <- Fitted + residuals[, detectors, drop = FALSE]
            peak_signal <- stats::quantile(Y, 0.995, na.rm = TRUE, names = FALSE)
        } else {
            expr_cols <- setdiff(colnames(res_obj$data), c("File", "Time"))
            expr_cols <- expr_cols[!grepl("^FSC|^SSC", expr_cols)]
            expr_mat <- as.matrix(res_obj$data[, expr_cols, drop = FALSE])
            peak_signal <- stats::quantile(expr_mat, 0.995, na.rm = TRUE, names = FALSE)
        }
        median_rms <- stats::median(rms, na.rm = TRUE)
        mean_rms <- mean(rms, na.rm = TRUE)
        error_ratio <- (median_rms / max(peak_signal, 100)) * 100
        data.frame(
            sample = sn,
            n_events = length(rms),
            median_rms_residual = median_rms,
            mean_rms_residual = mean_rms,
            q95_rms_residual = stats::quantile(rms, 0.95, na.rm = TRUE, names = FALSE),
            q995_peak_signal = peak_signal,
            median_rms_percent_of_peak = error_ratio,
            residual_metric = metric_name,
            unmixing_method = method,
            stringsAsFactors = FALSE
        )
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows) == 0) {
        return(NULL)
    }
    do.call(rbind, rows)
}

.split_qc_report_batches <- function(values, max_per_page = 15) {
    values <- as.character(values)
    values <- values[!is.na(values) & values != ""]
    if (length(values) == 0) {
        return(list())
    }

    max_per_page <- suppressWarnings(as.integer(max_per_page[1]))
    if (!is.finite(max_per_page) || max_per_page <= 0) {
        max_per_page <- 15L
    }

    split(values, ceiling(seq_along(values) / max_per_page))
}

.label_qc_report_batch_page <- function(p, page_idx, page_total, item_label = "Batch") {
    if (is.null(p) || page_total <= 1) {
        return(p)
    }

    p + ggplot2::labs(caption = paste0(item_label, " ", page_idx, " of ", page_total))
}

.split_qc_report_matrix_marker_batches <- function(marker_names, max_markers_per_page = 20) {
    marker_names <- as.character(marker_names)
    marker_names <- marker_names[!is.na(marker_names) & marker_names != ""]
    n_markers <- length(marker_names)
    if (n_markers == 0) {
        return(list())
    }

    max_markers_per_page <- suppressWarnings(as.integer(max_markers_per_page[1]))
    if (!is.finite(max_markers_per_page) || max_markers_per_page <= 0) {
        max_markers_per_page <- 20L
    }

    n_pages <- max(1L, ceiling(n_markers / max_markers_per_page))

    base_size <- n_markers %/% n_pages
    remainder <- n_markers %% n_pages
    sizes <- rep(base_size, n_pages)
    if (remainder > 0) {
        sizes[seq_len(remainder)] <- sizes[seq_len(remainder)] + 1L
    }

    ends <- cumsum(sizes)
    starts <- c(1L, head(ends, -1L) + 1L)
    Map(function(start, end) marker_names[seq.int(start, end)], starts, ends)
}

.build_qc_report_rms_pages <- function(res_list, M = NULL, max_files_per_page = 15, unmixing_method = NULL) {
    if (!is.list(res_list) || length(res_list) == 0) {
        return(list())
    }

    sample_names <- names(res_list)
    if (is.null(sample_names) || any(sample_names == "")) {
        p <- plot_sample_rms_residuals(res_list, M = M, output_file = NULL, unmixing_method = unmixing_method)
        return(if (is.null(p)) list() else list(p))
    }

    batches <- .split_qc_report_batches(sample_names, max_per_page = max_files_per_page)
    pages <- list()
    k <- 1L
    for (page_idx in seq_along(batches)) {
        p <- plot_sample_rms_residuals(res_list[batches[[page_idx]]], M = M, output_file = NULL, unmixing_method = unmixing_method)
        if (is.null(p)) next
        pages[[k]] <- .label_qc_report_batch_page(
            p,
            page_idx = page_idx,
            page_total = length(batches),
            item_label = "Files"
        )
        k <- k + 1L
    }

    pages
}

.build_qc_report_nps_pages <- function(nps_scores, max_files_per_page = 15) {
    if (!is.data.frame(nps_scores) || !("File" %in% colnames(nps_scores)) || nrow(nps_scores) == 0) {
        return(list())
    }
    if (!("NPS" %in% colnames(nps_scores))) {
        return(list())
    }
    nps_scores <- nps_scores[is.finite(nps_scores$NPS), , drop = FALSE]
    if (nrow(nps_scores) == 0L) {
        return(list())
    }

    batches <- .split_qc_report_batches(unique(as.character(nps_scores$File)), max_per_page = max_files_per_page)
    pages <- list()
    k <- 1L
    for (page_idx in seq_along(batches)) {
        page_scores <- nps_scores[nps_scores$File %in% batches[[page_idx]], , drop = FALSE]
        if (nrow(page_scores) == 0) next
        p <- plot_nps(page_scores, output_file = NULL)
        pages[[k]] <- .label_qc_report_batch_page(
            p,
            page_idx = page_idx,
            page_total = length(batches),
            item_label = "Files"
        )
        k <- k + 1L
    }

    pages
}

.build_qc_report_matrix_pages <- function(mat, plot_fun, max_markers_per_page = 20, item_label = "Markers") {
    if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
        return(list())
    }

    marker_names <- rownames(mat)
    if (is.null(marker_names) || any(marker_names == "")) {
        marker_names <- seq_len(nrow(mat))
    }

    shared_markers <- intersect(marker_names, colnames(mat))
    if (length(shared_markers) == 0) {
        p <- plot_fun(mat, output_file = NULL)
        return(if (is.null(p)) list() else list(p))
    }

    batches <- .split_qc_report_matrix_marker_batches(shared_markers, max_markers_per_page = max_markers_per_page)
    pages <- list()
    k <- 1L
    for (page_idx in seq_along(batches)) {
        markers <- batches[[page_idx]]
        p <- plot_fun(mat[markers, markers, drop = FALSE], output_file = NULL)
        if (is.null(p)) next
        pages[[k]] <- .label_qc_report_batch_page(
            p,
            page_idx = page_idx,
            page_total = length(batches),
            item_label = item_label
        )
        k <- k + 1L
    }

    pages
}

.normalize_qc_report_results_df <- function(results_df) {
    if (is.data.frame(results_df)) {
        return(as.data.frame(results_df, stringsAsFactors = FALSE, check.names = FALSE))
    }

    if (inherits(results_df, "spectreasy_unmixed_results")) {
        return(as.data.frame(results_df, stringsAsFactors = FALSE, check.names = FALSE))
    }

    if (is.list(results_df) && length(results_df) > 0) {
        looks_like_unmixed <- all(vapply(
            results_df,
            function(x) is.list(x) && is.data.frame(x$data),
            logical(1)
        ))
        if (looks_like_unmixed) {
            return(.as_unmixed_results_data_frame(results_df, arg_name = "results_df"))
        }
    }

    as.data.frame(results_df, stringsAsFactors = FALSE, check.names = FALSE)
}

.qc_report_cap_indices <- function(n, max_events_per_sample = NULL) {
    if (is.null(max_events_per_sample)) {
        return(seq_len(n))
    }
    cap <- suppressWarnings(as.integer(max_events_per_sample[1]))
    if (!is.finite(cap) || cap <= 0 || n <= cap) {
        return(seq_len(n))
    }
    unique(as.integer(round(seq(1, n, length.out = cap))))
}

.cap_qc_report_results_df <- function(results_df, max_events_per_sample = NULL) {
    if (is.null(max_events_per_sample) || !("File" %in% colnames(results_df))) {
        return(results_df)
    }

    groups <- split(seq_len(nrow(results_df)), results_df$File, drop = TRUE)
    keep <- unlist(lapply(groups, function(idx) {
        idx[.qc_report_cap_indices(length(idx), max_events_per_sample)]
    }), use.names = FALSE)
    results_df[sort(keep), , drop = FALSE]
}

.cap_qc_report_results_list <- function(results, max_events_per_sample = NULL) {
    if (is.null(max_events_per_sample) || !is.list(results) || is.data.frame(results)) {
        return(results)
    }

    lapply(results, function(res_obj) {
        if (!is.list(res_obj) || !is.data.frame(res_obj$data)) {
            return(res_obj)
        }
        idx <- .qc_report_cap_indices(nrow(res_obj$data), max_events_per_sample)
        res_obj$data <- res_obj$data[idx, , drop = FALSE]
        if (!is.null(res_obj$residuals) && nrow(res_obj$residuals) >= max(idx)) {
            res_obj$residuals <- res_obj$residuals[idx, , drop = FALSE]
        }
        res_obj
    })
}

.qc_report_file_counts <- function(results) {
    if (is.data.frame(results)) {
        if ("File" %in% colnames(results)) {
            return(table(results$File))
        }
        return(integer(0))
    }

    if (is.list(results) && length(results) > 0) {
        looks_like_unmixed <- all(vapply(
            results,
            function(x) is.list(x) && is.data.frame(x$data),
            logical(1)
        ))
        if (looks_like_unmixed) {
            counts <- vapply(results, function(x) nrow(x$data), integer(1))
            names(counts) <- names(results)
            return(counts)
        }
    }

    integer(0)
}

.get_qc_report_marker_labels <- function(markers) {
    stats::setNames(markers, markers)
}

.compute_qc_report_scatter_limits <- function(v) {
    v <- v[is.finite(v)]
    if (length(v) == 0) return(c(-1, 1))
    q <- stats::quantile(v, probs = c(0.005, 0.995), na.rm = TRUE, names = FALSE, type = 7)
    lo <- q[1]
    hi <- q[2]
    if (!is.finite(lo) || !is.finite(hi)) {
        lo <- min(v, na.rm = TRUE)
        hi <- max(v, na.rm = TRUE)
    }
    lo <- min(lo, 0)
    hi <- max(hi, 0)
    if (lo == hi) {
        pad <- max(1e-6, abs(lo) * 0.15)
        lo <- lo - pad
        hi <- hi + pad
    } else {
        pad <- (hi - lo) * 0.12
        lo <- lo - pad
        hi <- hi + pad
    }
    c(lo, hi)
}

.compute_qc_report_fixed_scatter_limits <- function(axis_limit,
                                                    transform = c("none", "asinh"),
                                                    asinh_cofactor = 150) {
    if (is.null(axis_limit)) {
        return(NULL)
    }
    transform <- .match_arg_ci(transform, c("none", "asinh"), "transform")
    axis_limit <- suppressWarnings(as.numeric(axis_limit)[1])
    if (!is.finite(axis_limit) || axis_limit <= 0) {
        return(NULL)
    }
    limits <- c(-axis_limit, axis_limit)
    if (transform == "asinh") {
        limits <- asinh(limits / asinh_cofactor)
    }
    limits
}

.add_qc_report_density_colors <- function(plot_df) {
    plot_df$color <- .scatter_density_base_color()
    groups <- split(seq_len(nrow(plot_df)), interaction(plot_df$panel_row, plot_df$panel_col, drop = TRUE))
    ramp <- .scatter_density_color_ramp()
    for (idx in groups) {
        if (length(idx) > 1 && stats::var(plot_df$x[idx]) > 0 && stats::var(plot_df$y[idx]) > 0) {
            plot_df$color[idx] <- .with_known_qc_plot_warnings_suppressed(
                grDevices::densCols(plot_df$x[idx], plot_df$y[idx], colramp = ramp)
            )
        }
    }
    plot_df
}

.pad_qc_report_page_levels <- function(real_names, block_size = 10, prefix = "pad") {
    real_names <- as.character(real_names)
    if (length(real_names) >= block_size) {
        return(real_names[seq_len(block_size)])
    }
    c(real_names, paste0("..", prefix, seq_len(block_size - length(real_names))))
}

.build_qc_report_sample_page_defs <- function(markers, block_size = 10) {
    n_markers <- length(markers)
    if (n_markers < 2) {
        return(list())
    }

    row_starts <- seq(1L, n_markers, by = block_size)
    page_defs <- list()
    k <- 1L

    for (row_start in row_starts) {
        row_end <- min(row_start + block_size - 1L, n_markers)
        col_starts <- seq(1L, row_end, by = block_size)
        for (col_start in col_starts) {
            col_end <- min(col_start + block_size - 1L, n_markers)
            row_idx <- seq.int(row_start, row_end)
            col_idx <- seq.int(col_start, col_end)
            keep_pairs <- outer(col_idx, row_idx, FUN = function(ci, ri) ci < ri)
            if (!any(keep_pairs)) next
            page_defs[[k]] <- list(row_idx = row_idx, col_idx = col_idx)
            k <- k + 1L
        }
    }

    page_defs
}
