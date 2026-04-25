# Internal helpers for QC report generation.
.draw_qc_report_plot_page <- function(p, square = FALSE) {
    if (is.null(p)) return(invisible(NULL))
    grid::grid.newpage()
    grob <- ggplot2::ggplotGrob(p)
    if (!isTRUE(square)) {
        grid::grid.draw(grob)
        return(invisible(NULL))
    }

    grid::grid.draw(
        grid::editGrob(
            grob,
            vp = grid::viewport(
                x = 0.5,
                y = 0.5,
                width = 0.78,
                height = 0.99
            )
        )
    )
    invisible(NULL)
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
                                              max_points = 3000) {
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

            x_lim <- .compute_qc_report_scatter_limits(x_vals)
            y_lim <- .compute_qc_report_scatter_limits(y_vals)
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
            ggplot2::aes(x = x, y = y),
            color = point_color,
            alpha = 0.55,
            size = 0.2,
            stroke = 0
        ) +
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
        ggplot2::theme_bw(base_size = 5.2) +
        ggplot2::theme(
            legend.position = "none",
            panel.grid = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(color = "grey60", linewidth = 0.2),
            strip.background = ggplot2::element_rect(fill = "grey95", color = "grey80"),
            strip.text = ggplot2::element_text(size = 4.5, face = "bold"),
            axis.title = ggplot2::element_blank(),
            axis.text = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            axis.line = ggplot2::element_blank(),
            panel.spacing = grid::unit(0.25, "mm"),
            plot.title = ggplot2::element_text(size = 9, face = "bold", hjust = 0.5),
            plot.subtitle = ggplot2::element_text(size = 6.2, hjust = 0.5),
            plot.margin = grid::unit(c(1, 1, 1, 1), "mm")
        )
}

.build_qc_report_sample_scatter_pages <- function(results_df,
                                                  markers,
                                                  rows_per_page = 10,
                                                  max_points = 3000,
                                                  transform = c("none", "asinh"),
                                                  asinh_cofactor = 150,
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
                max_points = max_points
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

#' Generate a Full QC PDF Report
#'
#' Creates a multi-page report summarizing unmixing quality, including spectra,
#' detector residuals, spread matrix, NPS, and per-sample NxN marker scatter pages.
#'
#' @param results_df Combined unmixed data frame (typically `rbind` of sample results).
#' @param M Reference matrix used for unmixing.
#' @param output_file Output PDF file path. Must be supplied explicitly.
#' @param res_list Optional residual object/list from `calc_residuals(..., return_residuals = TRUE)`.
#' @param png_dir Deprecated and ignored (kept for backward compatibility).
#' @param pd Optional detector metadata (`flowCore::pData(parameters(ff))`) for axis labels.
#'   If omitted, `attr(M, "detector_pd")` is used when available.
#' @param sample_nxn_rows_per_page Number of marker rows and columns to show per per-sample NxN page block.
#'   Defaults to 10, which standardizes panel geometry across pages and samples.
#' @param sample_nxn_max_points Maximum cells sampled per sample for each NxN page.
#' @param sample_nxn_transform One of `"none"` or `"asinh"` for per-sample NxN pages.
#' @param sample_nxn_asinh_cofactor Cofactor used when `sample_nxn_transform = "asinh"`.
#' @param nxn_all_samples Logical; if `TRUE`, include per-sample NxN pages for all samples.
#'   If `FALSE` (default), only include NxN pages for the first sample in `results_df`.
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
#' results_df <- data.frame(
#'   File = rep(c("sample_a", "sample_b"), each = 120),
#'   FITC = c(rnorm(120, 2, 0.4), rnorm(120, 0.1, 0.2)),
#'   PE = c(rnorm(120, 0.2, 0.2), rnorm(120, 2.5, 0.5)),
#'   APC = rnorm(240, 0.3, 0.3)
#' )
#'
#' pdf_file <- tempfile(fileext = ".pdf")
#' generate_qc_report(results_df = results_df, M = M_demo, output_file = pdf_file)
#' file.exists(pdf_file)
generate_qc_report <- function(results_df,
                               M,
                               output_file = NULL,
                               res_list = NULL,
                               png_dir = NULL,
                               pd = NULL,
                               sample_nxn_rows_per_page = 10,
                               sample_nxn_max_points = 3000,
                               sample_nxn_transform = c("none", "asinh"),
                               sample_nxn_asinh_cofactor = 150,
                               nxn_all_samples = FALSE) {
    if (is.null(output_file) || !nzchar(trimws(as.character(output_file)[1]))) {
        stop("Please supply output_file to save the QC PDF report.", call. = FALSE)
    }
    sample_nxn_transform <- match.arg(sample_nxn_transform)

    message("Generating spectreasy Summary Report...")
    if (!is.null(png_dir)) {
        warning("png_dir is deprecated and ignored; report output is PDF-only.")
    }

    M <- .as_reference_matrix(M, "M")
    if (is.null(pd)) {
        pd_attr <- attr(M, "detector_pd")
        if (is.data.frame(pd_attr)) {
            pd <- pd_attr
        }
    }

    results_df <- as.data.frame(results_df, stringsAsFactors = FALSE)
    if (!("File" %in% colnames(results_df))) {
        stop("results_df must contain a 'File' column.")
    }
    out_dir <- dirname(output_file)
    if (!is.na(out_dir) && nzchar(out_dir) && out_dir != ".") {
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    }
    grDevices::pdf(output_file, width = 11, height = 8.5)
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

    grid::grid.newpage()
    file_counts <- table(results_df$File)
    median_events_txt <- if (length(file_counts) > 0) {
        as.character(stats::median(as.numeric(file_counts)))
    } else {
        "NA"
    }

    keep_non_af <- !grepl("^AF($|_)", rownames(M), ignore.case = TRUE)
    M_no_af <- M[keep_non_af, , drop = FALSE]

    summary_txt <- paste0(
        "spectreasy: Spectral Unmixing Quality Control Report\n",
        "Generated on: ", Sys.time(), "\n\n",
        "Files Processed: ", length(unique(results_df$File)), "\n",
        "Total Markers (non-AF): ", nrow(M_no_af), "\n",
        "Median Events per File: ", median_events_txt
    )
    grid::grid.text(summary_txt, x = 0.5, y = 0.6, just = "center", gp = grid::gpar(fontsize = 15))

    message("  - Adding spectra overlay...")
    if (nrow(M_no_af) > 0) {
        .draw_qc_report_plot_page(plot_spectra(M_no_af, pd = pd, output_file = NULL))
    }

    if (!is.null(res_list)) {
        message("  - Adding detector-level residual diagnostics...")
        rep_res <- if (!is.null(res_list$residuals)) res_list else res_list[[1]]
        .draw_qc_report_plot_page(
            plot_detector_residuals(rep_res, M = if (nrow(M_no_af) > 0) M_no_af else M, top_n = 50, output_file = NULL)
        )
    }

    message("  - Adding Spread Matrix...")
    if (nrow(M_no_af) > 1) {
        ssm <- calculate_ssm(M_no_af)
        .draw_qc_report_plot_page(plot_ssm(ssm, output_file = NULL))
    }

    message("  - Adding NPS diagnostics...")
    nps_scores <- calculate_nps(results_df)
    nps_scores <- nps_scores[!grepl("^AF($|_)", nps_scores$Marker, ignore.case = TRUE), , drop = FALSE]
    if (nrow(nps_scores) > 0) {
        .draw_qc_report_plot_page(plot_nps(nps_scores, output_file = NULL))
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
        all_samples = nxn_all_samples
    )
    for (p in scatter_pages) {
        .draw_qc_report_plot_page(p, square = TRUE)
    }

    message("Report saved to: ", output_file)
}
