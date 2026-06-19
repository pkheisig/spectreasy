# Internal helpers for QC report generation.
.draw_qc_report_plot_page <- function(p, square = FALSE, height_ratio = 1.0, width_ratio = 1.0) {
    if (is.null(p)) return(invisible(NULL))
    grid::grid.newpage()
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

.build_qc_report_rms_pages <- function(res_list, M = NULL, max_files_per_page = 15) {
    if (!is.list(res_list) || length(res_list) == 0) {
        return(list())
    }

    sample_names <- names(res_list)
    if (is.null(sample_names) || any(sample_names == "")) {
        p <- plot_sample_rms_residuals(res_list, M = M, output_file = NULL)
        return(if (is.null(p)) list() else list(p))
    }

    batches <- .split_qc_report_batches(sample_names, max_per_page = max_files_per_page)
    pages <- list()
    k <- 1L
    for (page_idx in seq_along(batches)) {
        p <- plot_sample_rms_residuals(res_list[batches[[page_idx]]], M = M, output_file = NULL)
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

.build_qc_report_matrix_pages <- function(mat, plot_fun, max_markers_per_page = 15, item_label = "Markers") {
    if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
        return(list())
    }

    marker_names <- rownames(mat)
    if (is.null(marker_names) || any(marker_names == "")) {
        marker_names <- seq_len(nrow(mat))
    }

    batches <- .split_qc_report_batches(marker_names, max_per_page = max_markers_per_page)
    pages <- list()
    k <- 1L
    for (page_idx in seq_along(batches)) {
        markers <- batches[[page_idx]]
        markers <- intersect(markers, intersect(rownames(mat), colnames(mat)))
        if (length(markers) == 0) next
        p <- plot_fun(mat[markers, markers, drop = FALSE], output_file = NULL)
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
    transform <- match.arg(transform)
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
    plot_df$color <- "#0000FF"
    groups <- split(seq_len(nrow(plot_df)), interaction(plot_df$panel_row, plot_df$panel_col, drop = TRUE))
    ramp <- grDevices::colorRampPalette(c("#0000FF", "#00FFFF", "#00FF00", "#FFFF00", "#FF0000"))
    for (idx in groups) {
        if (length(idx) > 1 && stats::var(plot_df$x[idx]) > 0 && stats::var(plot_df$y[idx]) > 0) {
            plot_df$color[idx] <- grDevices::densCols(plot_df$x[idx], plot_df$y[idx], colramp = ramp)
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

.build_sample_scatter_page_panels <- function(sample_df,
                                              markers,
                                              row_markers,
                                              col_markers,
                                              block_size = 10,
                                              transform = c("none", "asinh"),
                                              asinh_cofactor = 150,
                                              max_points = 3000,
                                              fixed_limits = NULL) {
    transform <- match.arg(transform)
    row_markers <- intersect(as.character(row_markers), markers)
    col_markers <- intersect(as.character(col_markers), markers)
    if (length(row_markers) == 0 || length(col_markers) == 0) {
        return(NULL)
    }

    panel_data <- list()
    panel_limits <- list()
    k <- 1L
    lim_k <- 1L

    for (row_stain in row_markers) {
        if (!(row_stain %in% colnames(sample_df))) next
        row_idx <- match(row_stain, markers)
        if (is.na(row_idx) || row_idx <= 1L) next

        valid_x <- col_markers[match(col_markers, markers) < row_idx]
        valid_x <- valid_x[valid_x %in% colnames(sample_df)]
        if (length(valid_x) == 0) next

        for (xm in valid_x) {
            d_pair <- sample_df[, c(xm, row_stain), drop = FALSE]
            d_pair <- d_pair[stats::complete.cases(d_pair), , drop = FALSE]
            if (nrow(d_pair) == 0) next
            if (nrow(d_pair) > max_points) {
                d_pair <- d_pair[sample.int(nrow(d_pair), max_points), , drop = FALSE]
            }

            x_vals <- d_pair[[xm]]
            y_vals <- d_pair[[row_stain]]
            if (transform == "asinh") {
                x_vals <- asinh(x_vals / asinh_cofactor)
                y_vals <- asinh(y_vals / asinh_cofactor)
            }

            if (is.null(fixed_limits)) {
                x_lim <- .compute_qc_report_scatter_limits(x_vals)
                y_lim <- .compute_qc_report_scatter_limits(y_vals)
            } else {
                x_lim <- fixed_limits
                y_lim <- fixed_limits
            }
            x_plot <- pmax(pmin(x_vals, x_lim[2]), x_lim[1])
            y_plot <- pmax(pmin(y_vals, y_lim[2]), y_lim[1])

            panel_data[[k]] <- data.frame(
                x = x_plot,
                y = y_plot,
                panel_col = xm,
                panel_row = row_stain,
                stringsAsFactors = FALSE
            )
            panel_limits[[lim_k]] <- data.frame(
                panel_col = xm,
                panel_row = row_stain,
                x_low = x_lim[1],
                x_high = x_lim[2],
                y_low = y_lim[1],
                y_high = y_lim[2],
                stringsAsFactors = FALSE
            )
            k <- k + 1L
            lim_k <- lim_k + 1L
        }
    }

    if (length(panel_data) == 0) {
        return(NULL)
    }

    plot_df <- do.call(rbind, panel_data)
    lim_df <- do.call(rbind, panel_limits)
    plot_df <- .add_qc_report_density_colors(plot_df)
    page_rows <- .pad_qc_report_page_levels(row_markers, block_size = block_size, prefix = "row")
    page_cols <- .pad_qc_report_page_levels(col_markers, block_size = block_size, prefix = "col")

    plot_df$panel_col <- factor(plot_df$panel_col, levels = page_cols)
    plot_df$panel_row <- factor(plot_df$panel_row, levels = page_rows)
    lim_df$panel_col <- factor(lim_df$panel_col, levels = page_cols)
    lim_df$panel_row <- factor(lim_df$panel_row, levels = page_rows)

    list(plot_df = plot_df, lim_df = lim_df, page_cols = page_cols, page_rows = page_rows)
}

.build_sample_scatter_page_plot <- function(panel_info,
                                            marker_labels,
                                            sample_name,
                                            page_idx,
                                            page_total,
                                            point_color = "#163B5C",
                                            axis_color = "#C75000") {
    panel_used <- unique(panel_info$lim_df[, c("panel_row", "panel_col"), drop = FALSE])
    panel_grid <- expand.grid(
        panel_row = factor(panel_info$page_rows, levels = panel_info$page_rows),
        panel_col = factor(panel_info$page_cols, levels = panel_info$page_cols),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )

    ggplot2::ggplot() +
        ggplot2::geom_blank(data = panel_grid, ggplot2::aes(x = 0, y = 0)) +
        ggplot2::geom_blank(data = panel_info$lim_df, ggplot2::aes(x = x_low, y = y_low)) +
        ggplot2::geom_blank(data = panel_info$lim_df, ggplot2::aes(x = x_high, y = y_high)) +
        ggplot2::geom_hline(
            data = panel_used,
            ggplot2::aes(yintercept = 0),
            inherit.aes = FALSE,
            color = axis_color,
            linewidth = 0.15
        ) +
        ggplot2::geom_vline(
            data = panel_used,
            ggplot2::aes(xintercept = 0),
            inherit.aes = FALSE,
            color = axis_color,
            linewidth = 0.15
        ) +
        ggplot2::geom_point(
            data = panel_info$plot_df,
            ggplot2::aes(x = x, y = y, color = color),
            alpha = 0.55,
            size = 0.2,
            stroke = 0
        ) +
        ggplot2::scale_color_identity() +
        ggplot2::facet_grid(
            panel_row ~ panel_col,
            drop = FALSE,
            switch = "y",
            scales = "free",
            labeller = ggplot2::labeller(
                panel_row = marker_labels,
                panel_col = marker_labels
            )
        ) +
        ggplot2::labs(
            title = paste0("Sample NxN Scatter Matrix: ", sample_name),
            subtitle = paste0("Page ", page_idx, " of ", page_total)
        ) +
        ggplot2::theme_bw(base_size = 6.5) +
        ggplot2::theme(
            legend.position = "none",
            panel.grid = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(color = "grey60", linewidth = 0.2),
            strip.background = ggplot2::element_rect(fill = "grey95", color = "grey80"),
            strip.text = ggplot2::element_text(size = 5.625, face = "bold"),
            axis.title = ggplot2::element_blank(),
            axis.text = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            axis.line = ggplot2::element_blank(),
            panel.spacing = grid::unit(0.25, "mm"),
            plot.title = ggplot2::element_text(size = 11.25, face = "bold", hjust = 0.5),
            plot.subtitle = ggplot2::element_text(size = 9.3, hjust = 0.5),
            plot.margin = grid::unit(c(1, 1, 1, 1), "mm")
        )
}

.build_qc_report_sample_scatter_pages <- function(results_df,
                                                  markers,
                                                  rows_per_page = 10,
                                                  max_points = 3000,
                                                  transform = c("none", "asinh"),
                                                  asinh_cofactor = 150,
                                                  axis_limit = NULL,
                                                  all_samples = FALSE) {
    transform <- match.arg(transform)
    if (!("File" %in% colnames(results_df))) {
        return(list())
    }
    markers <- intersect(markers, colnames(results_df))
    if (length(markers) < 2) {
        return(list())
    }

    marker_labels <- .get_qc_report_marker_labels(markers)
    sample_names <- unique(as.character(results_df$File))
    sample_names <- sample_names[!is.na(sample_names) & sample_names != ""]
    if (!isTRUE(all_samples) && length(sample_names) > 0) {
        sample_names <- sample_names[1]
    }

    pages <- list()
    k <- 1L
    fixed_limits <- .compute_qc_report_fixed_scatter_limits(
        axis_limit = axis_limit,
        transform = transform,
        asinh_cofactor = asinh_cofactor
    )
    for (sample_name in sample_names) {
        sample_df <- results_df[results_df$File == sample_name, markers, drop = FALSE]
        page_defs <- .build_qc_report_sample_page_defs(markers, block_size = rows_per_page)
        page_total <- length(page_defs)
        if (page_total == 0) next

        for (page_idx in seq_along(page_defs)) {
            def <- page_defs[[page_idx]]
            row_markers <- markers[def$row_idx]
            col_markers <- markers[def$col_idx]
            panel_info <- .build_sample_scatter_page_panels(
                sample_df = sample_df,
                markers = markers,
                row_markers = row_markers,
                col_markers = col_markers,
                block_size = rows_per_page,
                transform = transform,
                asinh_cofactor = asinh_cofactor,
                max_points = max_points,
                fixed_limits = fixed_limits
            )
            if (is.null(panel_info)) next

            page_labels <- marker_labels
            missing_rows <- setdiff(panel_info$page_rows, names(page_labels))
            missing_cols <- setdiff(panel_info$page_cols, names(page_labels))
            if (length(missing_rows) > 0) page_labels[missing_rows] <- ""
            if (length(missing_cols) > 0) page_labels[missing_cols] <- ""

            pages[[k]] <- .build_sample_scatter_page_plot(
                panel_info = panel_info,
                marker_labels = page_labels,
                sample_name = sample_name,
                page_idx = page_idx,
                page_total = page_total
            )
            k <- k + 1L
        }
    }

    pages
}

#' Generate a Full Sample PDF Report
#'
#' Creates a multi-page report summarizing unmixing quality, including spectra,
#' detector residuals, spread matrix, NPS, and per-sample NxN marker scatter pages.
#'
#' `qc_samples()` expects a combined data frame or the raw list returned
#' by [unmix_samples()]. In the usual workflow, pass the unmixed results object
#' directly.
#'
#' `M` should be the reference matrix used for the unmixing context, supplied as
#' a numeric matrix or detector-column data frame. In the usual
#' `unmix_controls()` workflow, pass `ctrl$M` or load
#' `scc_reference_matrix.csv`. Do not pass the path to
#' `scc_unmixing_matrix.csv` here.
#'
#' @param results Combined unmixed data frame, or the list returned by
#'   [unmix_samples()]. `qc_samples()` will automatically bind
#'   per-sample `$data` elements when needed.
#' @param M Reference matrix used for report context. Must be a numeric matrix or
#'   a data frame with detector columns, for example `ctrl$M` from
#'   `unmix_controls()` or `utils::read.csv("scc_reference_matrix.csv",
#'   check.names = FALSE)`.
#' @param unmixing_matrix_file Optional CSV path to a saved reference matrix.
#'   Used when `M` is not supplied. By default this points to the reference matrix
#'   produced by [unmix_controls()] (`"scc_reference_matrix.csv"`).
#' @param output_file Output PDF file path. Defaults to `"spectreasy_outputs/unmix_samples/qc_samples_report.pdf"`.
#' @param method Unmixing method used to create `results` (`"WLS"`, `"RWLS"`,
#'   `"OLS"`, or `"NNLS"`). When `"NNLS"`, the negative population spread page is skipped
#'   because constrained NNLS results are non-negative by construction.
#' @param res_list Optional residual object/list from `calc_residuals(..., return_residuals = TRUE)`.
#' @param png_dir Deprecated and ignored (kept for backward compatibility).
#' @param pd Optional detector metadata (`flowCore::pData(parameters(ff))`) for axis labels.
#'   If omitted, `attr(M, "detector_pd")` is used when available.
#' @param max_events_per_sample Maximum events per sample used for report-wide
#'   plots and diagnostics. Defaults to 1000 to keep large FCS reports
#'   responsive. Set `NULL` to use all events.
#' @param overview_files_per_page Maximum files shown on each RMS residual and
#'   NPS overview page.
#' @param matrix_markers_per_page Maximum markers shown on each similarity and
#'   spread matrix page.
#' @param sample_nxn_rows_per_page Number of marker rows and columns to show per per-sample NxN page block.
#'   Defaults to 10, which standardizes geometry across pages and samples.
#' @param sample_nxn_max_points Maximum cells sampled per sample for each NxN page.
#' @param sample_nxn_transform One of `"none"` or `"asinh"` for per-sample NxN pages.
#' @param sample_nxn_asinh_cofactor Cofactor used when `sample_nxn_transform = "asinh"`.
#' @param sample_nxn_axis_limit Optional fixed symmetric NxN scatter axis limit.
#'   The default `NULL` uses local per-panel ranges. Use `1e5` for
#'   `c(-1e5, 1e5)` on every NxN panel.
#' @param nxn_all_samples Logical; if `TRUE`, include per-sample NxN pages for all samples.
#'   If `FALSE` (default), only include NxN pages for the first sample in `results`.
#'
#' @return Invisibly returns `NULL`; writes report to disk.
#' @export
#' @examples
#' M_demo <- rbind(
#'   FITC = c(1.00, 0.20, 0.05),
#'   PE = c(0.10, 1.00, 0.20),
#'   APC = c(0.05, 0.15, 1.00)
#' )
#' colnames(M_demo) <- c("B2-A", "YG1-A", "R1-A")
#'
#' results <- data.frame(
#'   File = rep(c("sample_a", "sample_b"), each = 120),
#'   FITC = c(rnorm(120, 2, 0.4), rnorm(120, 0.1, 0.2)),
#'   PE = c(rnorm(120, 0.2, 0.2), rnorm(120, 2.5, 0.5)),
#'   APC = rnorm(240, 0.3, 0.3)
#' )
#'
#' # Typical workflow after unmix_samples():
#' # qc_samples(
#' #   results = unmixed,
#' #   M = ctrl$M
#' # )
#'
#' pdf_file <- tempfile(fileext = ".pdf")
#' qc_samples(results = results, M = M_demo, output_file = pdf_file)
#' file.exists(pdf_file)
qc_samples <- function(results,
                       M = NULL,
                       unmixing_matrix_file = file.path("spectreasy_outputs", "unmix_controls", "scc_reference_matrix.csv"),
                       output_file = "spectreasy_outputs/unmix_samples/qc_samples_report.pdf",
                       method = NULL,
                       res_list = NULL,
                       png_dir = NULL,
                       pd = NULL,
                       max_events_per_sample = 1000,
                       overview_files_per_page = 15,
                       matrix_markers_per_page = 15,
                       sample_nxn_rows_per_page = 10,
                       sample_nxn_max_points = max_events_per_sample,
                       sample_nxn_transform = c("none", "asinh"),
                       sample_nxn_asinh_cofactor = 150,
                       sample_nxn_axis_limit = NULL,
                       nxn_all_samples = FALSE) {
    if (is.null(output_file) || !nzchar(trimws(as.character(output_file)[1]))) {
        stop("Please supply output_file to save the QC PDF report.", call. = FALSE)
    }
    sample_nxn_transform <- match.arg(sample_nxn_transform)
    method_attr <- attr(results, "method")
    if (is.null(method)) {
        method <- if (!is.null(method_attr)) method_attr else "WLS"
    }
    method <- toupper(as.character(method)[1])
    if (!(method %in% c("WLS", "RWLS", "OLS", "NNLS"))) {
        stop("method must be one of: WLS, RWLS, OLS, NNLS", call. = FALSE)
    }

    message("Generating spectreasy Summary Report...")
    if (!is.null(png_dir)) {
        warning("png_dir is deprecated and ignored; report output is PDF-only.")
    }

    if (!is.null(M)) {
        M <- .as_reference_matrix(M, "M")
    } else if (!is.null(unmixing_matrix_file)) {
        if (!file.exists(unmixing_matrix_file)) {
            stop("unmixing_matrix_file not found: ", unmixing_matrix_file)
        }
        .stop_if_static_unmixing_matrix_path(unmixing_matrix_file, arg_name = "unmixing_matrix_file")
        M <- .read_unmixing_matrix_csv(unmixing_matrix_file)
        M <- .as_reference_matrix(M, "M")
    }

    if (is.null(M)) {
        stop("No reference matrix provided. Supply either M or a valid unmixing_matrix_file.")
    }

    if (is.null(pd)) {
        pd_attr <- attr(M, "detector_pd")
        if (is.data.frame(pd_attr)) {
            pd <- pd_attr
        }
    }

    report_results <- if (is.data.frame(results)) {
        .cap_qc_report_results_df(results, max_events_per_sample = max_events_per_sample)
    } else {
        .cap_qc_report_results_list(results, max_events_per_sample = max_events_per_sample)
    }

    if (is.null(res_list) && is.list(report_results) && !is.data.frame(report_results)) {
        has_residuals <- any(vapply(report_results, function(x) is.list(x) && !is.null(x$residuals), logical(1)))
        if (has_residuals) {
            res_list <- report_results
        }
    } else if (!is.null(res_list)) {
        res_list <- .cap_qc_report_results_list(res_list, max_events_per_sample = max_events_per_sample)
    }

    file_counts <- .qc_report_file_counts(results)
    results_df <- .normalize_qc_report_results_df(report_results)
    if (!("File" %in% colnames(results_df))) {
        stop("results must contain a 'File' column.")
    }
    out_dir <- dirname(output_file)
    if (!is.na(out_dir) && nzchar(out_dir) && out_dir != ".") {
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    }
    grDevices::pdf(output_file, width = 11, height = 8.5)
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

    grid::grid.newpage()
    if (!is.null(max_events_per_sample) && any(as.numeric(file_counts) > as.numeric(max_events_per_sample)[1], na.rm = TRUE)) {
        message(
            "  - Report diagnostics capped at ",
            as.integer(max_events_per_sample[1]),
            " events per sample. Set max_events_per_sample = NULL to use all events."
        )
    }

    keep_non_af <- !grepl("^AF($|_)", rownames(M), ignore.case = TRUE)
    M_no_af <- M[keep_non_af, , drop = FALSE]

    summary_txt <- paste0(
        "spectreasy: Spectral Unmixing Quality Control Report\n",
        "Generated on: ", Sys.time(), "\n\n",
        "Files Processed: ", length(unique(results_df$File)), "\n",
        "Total Markers (non-AF): ", nrow(M_no_af)
    )
    grid::grid.text(summary_txt, x = 0.5, y = 0.6, just = "center", gp = grid::gpar(fontsize = 15))

    if (!is.null(res_list)) {
        message("  - Adding overall RMS residuals page...")
        rms_pages <- .build_qc_report_rms_pages(
            res_list,
            M = M,
            max_files_per_page = overview_files_per_page
        )
        for (p in rms_pages) {
            .draw_qc_report_plot_page(
                p,
                height_ratio = 0.72,
                width_ratio = 0.72
            )
        }
    }

    message("  - Adding spectra overlay...")
    if (nrow(M_no_af) > 0) {
        .draw_qc_report_plot_page(plot_spectra(M_no_af, pd = pd, output_file = NULL), height_ratio = 0.6)
    }

    if (nrow(M_no_af) > 1) {
        message("  - Adding Fluorophore Similarity Matrix...")
        sim_mat <- calculate_similarity_matrix(M_no_af)
        sim_pages <- .build_qc_report_matrix_pages(
            sim_mat,
            plot_fun = plot_similarity_matrix,
            max_markers_per_page = matrix_markers_per_page,
            item_label = "Markers"
        )
        for (p in sim_pages) {
            .draw_qc_report_plot_page(p)
        }
    }

    message("  - Adding Spread Matrix...")
    if (nrow(M_no_af) > 1) {
        ssm_method <- if (method %in% c("NNLS", "RWLS")) {
            if (identical(method, "RWLS")) "WLS" else "OLS"
        } else {
            method
        }
        ssm <- calculate_ssm(M_no_af, method = ssm_method)
        ssm_pages <- .build_qc_report_matrix_pages(
            ssm,
            plot_fun = plot_ssm,
            max_markers_per_page = matrix_markers_per_page,
            item_label = "Markers"
        )
        for (p in ssm_pages) {
            .draw_qc_report_plot_page(p)
        }
    }

    if (identical(method, "NNLS")) {
        message("  - Skipping NPS diagnostics for NNLS...")
    } else {
        message("  - Adding NPS diagnostics...")
        nps_scores <- calculate_nps(results_df)
        nps_scores <- nps_scores[!grepl("^AF($|_)", nps_scores$Marker, ignore.case = TRUE), , drop = FALSE]
        if (nrow(nps_scores) > 0) {
            nps_pages <- .build_qc_report_nps_pages(
                nps_scores,
                max_files_per_page = overview_files_per_page
            )
            for (p in nps_pages) {
                .draw_qc_report_plot_page(p)
            }
        }
    }

    message("  - Adding per-sample NxN scatter pages...")
    sample_markers <- setdiff(colnames(results_df), .get_result_metadata_columns(colnames(results_df)))
    sample_markers <- sample_markers[!grepl("^AF($|_)", sample_markers, ignore.case = TRUE)]
    scatter_pages <- .build_qc_report_sample_scatter_pages(
        results_df = results_df,
        markers = sample_markers,
        rows_per_page = sample_nxn_rows_per_page,
        max_points = sample_nxn_max_points,
        transform = sample_nxn_transform,
        asinh_cofactor = sample_nxn_asinh_cofactor,
        axis_limit = sample_nxn_axis_limit,
        all_samples = nxn_all_samples
    )
    for (p in scatter_pages) {
        .draw_qc_report_plot_page(p, square = TRUE)
    }

    message("Report saved to: ", output_file)
}
