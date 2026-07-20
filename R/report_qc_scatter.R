.build_sample_scatter_page_panels <- function(sample_df,
                                              markers,
                                              row_markers,
                                              col_markers,
                                              block_size = 10,
                                              transform = c("none", "asinh"),
                                              asinh_cofactor = 150,
                                              max_points = 3000,
                                              fixed_limits = NULL) {
    transform <- .match_arg_ci(transform, c("none", "asinh"), "transform")
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
    page_rows <- unique(as.character(row_markers))
    page_cols <- unique(as.character(col_markers))

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
                                            axis_color = "#C75000",
                                            point_size = 0.2,
                                            point_alpha = 0.55,
                                            title = NULL) {
    if (is.null(title)) {
        title <- paste0("Sample NxN Scatter Matrix: ", sample_name)
    }
    panel_used <- unique(panel_info$lim_df[, c("panel_row", "panel_col"), drop = FALSE])
    panel_grid <- expand.grid(
        panel_row = factor(panel_info$page_rows, levels = panel_info$page_rows),
        panel_col = factor(panel_info$page_cols, levels = panel_info$page_cols),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )

    p <- ggplot2::ggplot() +
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
            alpha = point_alpha,
            size = point_size,
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
            title = title,
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
    attr(p, "spectreasy_nxn_rows") <- length(panel_info$page_rows)
    attr(p, "spectreasy_nxn_cols") <- length(panel_info$page_cols)
    p
}

.build_qc_report_sample_scatter_pages <- function(results_df,
                                                  markers,
                                                  rows_per_page = 10,
                                                  max_points = 3000,
                                                  transform = c("none", "asinh"),
                                                  asinh_cofactor = 150,
                                                  axis_limit = NULL,
                                                  all_samples = FALSE) {
    transform <- .match_arg_ci(transform, c("none", "asinh"), "transform")
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

.build_qc_report_control_scatter_pages <- function(unmixed_list,
                                                   sample_to_marker = NULL,
                                                   markers = NULL,
                                                   marker_display = NULL,
                                                   rows_per_page = 10,
                                                   max_points_per_sample = 1000,
                                                   transform = c("none", "asinh"),
                                                   asinh_cofactor = 150,
                                                   axis_limit = NULL,
                                                   point_size = 0.2,
                                                   point_alpha = 0.55,
                                                   seed = NULL) {
    transform <- .match_arg_ci(transform, c("none", "asinh"), "transform")
    .with_optional_seed(seed)

    data_list <- .extract_unmix_scatter_data_list(unmixed_list)
    sample_stains <- .normalize_unmix_scatter_mapping(names(data_list), sample_to_marker = sample_to_marker)
    marker_info <- .resolve_unmix_scatter_markers(sample_stains, markers = markers, marker_display = marker_display)
    markers <- marker_info$markers
    marker_labels <- marker_info$marker_labels
    panel_info <- .build_unmix_scatter_panel_data(
        data_list = data_list,
        sample_stains = sample_stains,
        markers = markers,
        max_points_per_sample = max_points_per_sample,
        transform = transform,
        asinh_cofactor = asinh_cofactor,
        axis_limit = axis_limit
    )

    page_defs <- .build_qc_report_sample_page_defs(markers, block_size = rows_per_page)
    page_total <- length(page_defs)
    if (page_total == 0) {
        return(list())
    }

    pages <- list()
    k <- 1L
    for (page_idx in seq_along(page_defs)) {
        def <- page_defs[[page_idx]]
        row_markers <- markers[def$row_idx]
        col_markers <- markers[def$col_idx]
        keep_rows <- as.character(panel_info$plot_df$panel_row) %in% row_markers
        keep_cols <- as.character(panel_info$plot_df$panel_col) %in% col_markers
        keep_lim_rows <- as.character(panel_info$lim_df$panel_row) %in% row_markers
        keep_lim_cols <- as.character(panel_info$lim_df$panel_col) %in% col_markers
        plot_df <- panel_info$plot_df[keep_rows & keep_cols, , drop = FALSE]
        lim_df <- panel_info$lim_df[keep_lim_rows & keep_lim_cols, , drop = FALSE]
        if (nrow(plot_df) == 0 || nrow(lim_df) == 0) {
            next
        }

        page_rows <- unique(as.character(row_markers))
        page_cols <- unique(as.character(col_markers))
        plot_df$panel_row <- factor(as.character(plot_df$panel_row), levels = page_rows)
        plot_df$panel_col <- factor(as.character(plot_df$panel_col), levels = page_cols)
        lim_df$panel_row <- factor(as.character(lim_df$panel_row), levels = page_rows)
        lim_df$panel_col <- factor(as.character(lim_df$panel_col), levels = page_cols)

        page_labels <- marker_labels
        missing_rows <- setdiff(page_rows, names(page_labels))
        missing_cols <- setdiff(page_cols, names(page_labels))
        if (length(missing_rows) > 0) page_labels[missing_rows] <- ""
        if (length(missing_cols) > 0) page_labels[missing_cols] <- ""

        pages[[k]] <- .build_sample_scatter_page_plot(
            panel_info = list(plot_df = plot_df, lim_df = lim_df, page_cols = page_cols, page_rows = page_rows),
            marker_labels = page_labels,
            sample_name = "Controls",
            page_idx = page_idx,
            page_total = page_total,
            point_size = point_size,
            point_alpha = point_alpha,
            title = "Control NxN Scatter Matrix"
        )
        k <- k + 1L
    }

    pages
}

.qc_report_sample_names <- function(results) {
    results_df <- .normalize_qc_report_results_df(results)
    if (!("File" %in% colnames(results_df))) return(character())
    sample_names <- unique(as.character(results_df$File))
    sample_names[!is.na(sample_names) & nzchar(sample_names)]
}

.subset_qc_report_sample <- function(results, sample_name) {
    if (is.null(results)) return(NULL)
    original_attributes <- attributes(results)
    if (is.data.frame(results)) {
        if (!("File" %in% colnames(results))) return(results)
        out <- results[as.character(results$File) == sample_name, , drop = FALSE]
    } else if (is.list(results)) {
        result_names <- names(results)
        if (is.null(result_names)) return(results)
        out <- results[result_names == sample_name]
    } else {
        return(results)
    }
    preserved <- setdiff(names(original_attributes), c("names", "row.names", "class", "dim", "dimnames"))
    for (attribute_name in preserved) attr(out, attribute_name) <- original_attributes[[attribute_name]]
    out
}
