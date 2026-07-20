.reference_qc_plot_context <- function(raw_data,
                                       fsc,
                                       ssc,
                                       peak_channel,
                                       final_gate,
                                       manual_gate_info,
                                       source_filename,
                                       sn) {
    gate_settings <- if (!is.null(manual_gate_info) && !is.null(manual_gate_info$settings)) {
        manual_gate_info$settings
    } else {
        .reference_manual_gate_settings(NULL)
    }
    if (is.null(source_filename) || !nzchar(as.character(source_filename)[1])) {
        source_filename <- paste0(sn, ".fcs")
    }
    cell <- if (!is.null(manual_gate_info) && !is.null(manual_gate_info$cell)) manual_gate_info$cell else NULL
    singlet <- if (!is.null(manual_gate_info) && !is.null(manual_gate_info$singlet)) manual_gate_info$singlet else NULL
    cell_x <- if (!is.null(cell$x_channel) && cell$x_channel %in% colnames(raw_data)) cell$x_channel else fsc
    cell_y <- if (!is.null(cell$y_channel) && cell$y_channel %in% colnames(raw_data)) cell$y_channel else ssc
    singlet_x <- if (!is.null(singlet$x_channel) && singlet$x_channel %in% colnames(raw_data)) singlet$x_channel else NA_character_
    singlet_y <- if (!is.null(singlet$y_channel) && singlet$y_channel %in% colnames(raw_data)) singlet$y_channel else NA_character_
    list(
        settings = gate_settings,
        cell = cell,
        singlet = singlet,
        cell_x = cell_x,
        cell_y = cell_y,
        cell_vertices = if (!is.null(cell$vertices)) cell$vertices else final_gate,
        singlet_x = singlet_x,
        singlet_y = singlet_y,
        singlet_vertices = singlet$vertices,
        display_data = .reference_gui_display_data(
            raw_data = raw_data,
            filename = source_filename,
            channels = c(cell_x, cell_y, singlet_x, singlet_y, peak_channel),
            display_max_points = gate_settings$max_points,
            preload_max_points = 100000L
        )
    )
}

.reference_qc_cell_source <- function(context, raw_data, gated_data, pd, gui_scatter_domains, sample_type, sn, out_path) {
    keep_count <- if (!is.null(context$cell$keep_count) && is.finite(context$cell$keep_count)) {
        as.numeric(context$cell$keep_count)
    } else {
        nrow(gated_data)
    }
    plot <- .reference_qc_scatter_density_plot(
        data = context$display_data,
        x_channel = context$cell_x,
        y_channel = context$cell_y,
        pd = pd,
        gate_vertices = context$cell_vertices,
        title = paste0(sn, " - Cell Gate"),
        subtitle = paste0(round(100 * keep_count / nrow(raw_data), 1), "% gated"),
        max_points = nrow(context$display_data),
        point_size = context$settings$point_size,
        x_domain = .reference_gui_domain_for_channel(gui_scatter_domains, context$cell_x, raw_data[, context$cell_x]),
        y_domain = .reference_gui_domain_for_channel(gui_scatter_domains, context$cell_y, raw_data[, context$cell_y]),
        sample_type = sample_type
    )
    .save_reference_ggsave(
        file.path(out_path, "fsc_ssc", paste0(sn, "_fsc_ssc.png")),
        plot, sn, "FSC/SSC", width = 5.2, height = 4.2, dpi = 300
    )
    vertices <- context$cell_vertices
    if (is.null(vertices) || nrow(vertices) < 3 ||
        !context$cell_x %in% colnames(context$display_data) ||
        !context$cell_y %in% colnames(context$display_data)) {
        return(context$display_data)
    }
    keep <- sp::point.in.polygon(
        context$display_data[, context$cell_x], context$display_data[, context$cell_y],
        vertices$x, vertices$y
    ) > 0
    context$display_data[keep, , drop = FALSE]
}

.reference_qc_singlet_source <- function(context, histogram_source, raw_data, gated_data, pd, gui_scatter_domains, sample_type, sn, out_path) {
    vertices <- context$singlet_vertices
    if (is.null(vertices) || nrow(vertices) < 3 || is.na(context$singlet_x) || is.na(context$singlet_y)) {
        return(histogram_source)
    }
    keep_count <- if (!is.null(context$singlet$keep_count) && is.finite(context$singlet$keep_count)) {
        as.numeric(context$singlet$keep_count)
    } else {
        nrow(gated_data)
    }
    source_count <- if (!is.null(context$singlet$source_count) && is.finite(context$singlet$source_count)) {
        as.numeric(context$singlet$source_count)
    } else {
        nrow(histogram_source)
    }
    plot <- .reference_qc_scatter_density_plot(
        data = histogram_source,
        x_channel = context$singlet_x,
        y_channel = context$singlet_y,
        pd = pd,
        gate_vertices = vertices,
        title = paste0(sn, " - Singlet Gate"),
        subtitle = paste0(round(100 * keep_count / max(source_count, 1), 1), "% gated"),
        max_points = nrow(histogram_source),
        point_size = context$settings$point_size,
        x_domain = .reference_gui_domain_for_channel(gui_scatter_domains, context$singlet_x, raw_data[, context$singlet_x]),
        y_domain = .reference_gui_domain_for_channel(gui_scatter_domains, context$singlet_y, raw_data[, context$singlet_y]),
        sample_type = sample_type
    )
    .save_reference_ggsave(
        file.path(out_path, "singlet", paste0(sn, "_singlet.png")),
        plot, sn, "singlet", width = 5.2, height = 4.2, dpi = 300
    )
    if (!context$singlet_x %in% colnames(histogram_source) || !context$singlet_y %in% colnames(histogram_source)) {
        return(histogram_source[0, , drop = FALSE])
    }
    keep <- sp::point.in.polygon(
        histogram_source[, context$singlet_x], histogram_source[, context$singlet_y],
        vertices$x, vertices$y
    ) > 0
    histogram_source[keep, , drop = FALSE]
}

.reference_qc_histogram_stage <- function(histogram_source,
                                          final_gated_data,
                                          pd,
                                          peak_channel,
                                          vals_log,
                                          gate_min,
                                          gate_max,
                                          gate_settings,
                                          hist_info,
                                          histogram_domain,
                                          sn,
                                          out_path) {
    peak_vals <- if (peak_channel %in% colnames(histogram_source)) histogram_source[, peak_channel] else numeric()
    positive_valid <- isTRUE(attr(vals_log, "positive_gate_present")) &&
        is.finite(gate_min) && is.finite(gate_max) && gate_max > gate_min
    spectrum_data <- if (positive_valid && peak_channel %in% colnames(histogram_source)) {
        histogram_source[peak_vals >= gate_min & peak_vals <= gate_max, , drop = FALSE]
    } else {
        final_gated_data
    }
    plot <- .reference_qc_histogram_gui_plot(
        peak_vals = peak_vals,
        pd = pd,
        peak_channel = peak_channel,
        vals_log = vals_log,
        gate_min = gate_min,
        gate_max = gate_max,
        settings = gate_settings,
        hist_info = hist_info,
        x_domain = histogram_domain
    )
    .save_reference_ggsave(
        file.path(out_path, "histogram", paste0(sn, "_histogram.png")),
        plot, sn, "histogram", width = 6.5, height = 4, dpi = 300
    )
    spectrum_data
}

.reference_spectrum_histogram_data <- function(data, detector_names) {
    log_mat <- if (nrow(data) > 0L && all(detector_names %in% colnames(data))) {
        log10(pmax(data[, detector_names, drop = FALSE], 1e-3))
    } else {
        matrix(numeric(), nrow = 0L, ncol = length(detector_names), dimnames = list(NULL, detector_names))
    }
    finite <- log_mat[is.finite(log_mat)]
    if (length(finite) == 0L) finite <- c(0, 1)
    min_y <- floor(min(finite, na.rm = TRUE))
    max_y <- ceiling(max(finite, na.rm = TRUE))
    breaks <- seq(min_y, max_y, length.out = 151)
    bin_mid <- (breaks[-1] + breaks[-length(breaks)]) / 2
    counts <- vapply(
        seq_len(ncol(log_mat)),
        function(j) as.numeric(graphics::hist(log_mat[, j], breaks = breaks, plot = FALSE)$counts),
        numeric(length(bin_mid))
    )
    rownames(counts) <- as.character(seq_along(bin_mid))
    colnames(counts) <- as.character(seq_len(ncol(log_mat)))
    data <- data.table::as.data.table(as.table(counts))
    data.table::setnames(data, c("bin_idx", "ch_idx", "count"))
    data$bin_idx <- as.integer(as.character(data$bin_idx))
    data$ch_idx <- as.integer(as.character(data$ch_idx))
    data$y_orig <- bin_mid[data$bin_idx]
    data$fill <- log10(data$count + 1)
    data <- data[data$count >= 3, ]
    data$y <- data$y_orig^1.5
    list(data = data, min_y = min_y, max_y = max_y, bin_height = breaks[2] - breaks[1])
}

.save_reference_spectrum_qc <- function(data, detector_names, detector_labels, det_info, sn, out_path) {
    histogram <- .reference_spectrum_histogram_data(data, detector_names)
    plot_data <- histogram$data
    if (nrow(plot_data) == 0) {
        .spectreasy_console_step("QC plot", paste0(sn, " spectrum histogram is empty"))
    }
    fill_lo <- suppressWarnings(min(plot_data$fill, na.rm = TRUE))
    fill_hi <- suppressWarnings(stats::quantile(plot_data$fill, 0.96, na.rm = TRUE))
    if (!is.finite(fill_lo) || !is.finite(fill_hi) || fill_hi <= fill_lo) {
        fill_lo <- 0
        fill_hi <- 1
    }
    y_breaks <- 0:ceiling(histogram$max_y)
    plot <- ggplot2::ggplot(plot_data, ggplot2::aes(ch_idx, y, fill = fill)) +
        ggplot2::geom_tile(width = 0.7, height = histogram$bin_height * 3) +
        ggplot2::scale_fill_gradientn(
            colors = c("#0000FF", "#00FFFF", "#00FF00", "#FFFF00", "#FF0000"),
            limits = c(fill_lo, fill_hi), oob = scales::squish
        ) +
        ggplot2::scale_x_continuous(breaks = seq_along(detector_names), labels = detector_labels) +
        ggplot2::scale_y_continuous(
            limits = c(0, (histogram$max_y + 0.5)^1.5),
            breaks = y_breaks^1.5,
            labels = paste0("10^", y_breaks)
        ) +
        ggplot2::coord_cartesian(expand = FALSE) +
        ggplot2::labs(title = paste0(sn, " - Spectrum"), x = NULL, y = "Intensity") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
            legend.position = "none",
            panel.background = ggplot2::element_rect(fill = "white", color = NA)
        )
    .save_reference_ggsave(
        file.path(out_path, "spectrum", paste0(sn, "_spectrum.png")),
        plot, sn, "spectrum", width = 300, height = 120, units = "mm", dpi = 600
    )
}
