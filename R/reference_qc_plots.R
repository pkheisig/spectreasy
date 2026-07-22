.get_reference_axis_label <- function(ch_name, pd_tbl) {
    idx <- match(ch_name, as.character(pd_tbl$name))
    if (!is.na(idx) && "desc" %in% colnames(pd_tbl)) {
        desc_val <- trimws(as.character(pd_tbl$desc[idx]))
        if (!is.na(desc_val) && nzchar(desc_val)) {
            return(desc_val)
        }
    }
    ch_name
}

.warn_reference_qc_plot_failure <- function(sn, plot_type, condition) {
    warning(
        "SCC QC plot skipped for ", sn, " (", plot_type, "): ",
        conditionMessage(condition),
        ". Continuing unmixing.",
        call. = FALSE
    )
    invisible(FALSE)
}

.save_reference_ggsave <- function(filename, plot, sn, plot_type, ...) {
    tryCatch(
        {
            ggplot2::ggsave(filename, plot, ...)
            invisible(TRUE)
        },
        error = function(e) .warn_reference_qc_plot_failure(sn, plot_type, e)
    )
}

.save_reference_qc_plots_safely <- function(...) {
    args <- list(...)
    sn <- if (!is.null(args$sn)) args$sn else "control"
    tryCatch(
        do.call(.save_reference_qc_plots, args),
        error = function(e) .warn_reference_qc_plot_failure(sn, "plot bundle", e)
    )
}

.reference_even_indices <- function(n, max_points = 50000L) {
    n <- as.integer(n)
    max_points <- suppressWarnings(as.integer(max_points))
    if (!is.finite(max_points) || is.na(max_points) || max_points <= 0L) max_points <- 50000L
    if (n <= max_points) return(seq_len(n))
    pmax(1L, pmin(n, floor((seq_len(max_points) - 1L) * n / max_points) + 1L))
}

.reference_density_palette <- function(size = 64L) {
    vapply(seq_len(size), function(i) {
        v <- (i - 1) / max(size - 1, 1)
        r <- max(0, min(1, min(4 * v - 1.5, -4 * v + 4.5)))
        g <- max(0, min(1, min(4 * v - 0.5, -4 * v + 3.5)))
        b <- max(0, min(1, min(4 * v + 0.5, -4 * v + 2.5)))
        grDevices::rgb(r, g, b, alpha = 0.70)
    }, character(1))
}

.reference_point_density <- function(x, y) {
    n <- length(x)
    if (n == 0) return(numeric())
    min_x <- min(x, na.rm = TRUE)
    max_x <- max(x, na.rm = TRUE)
    min_y <- min(y, na.rm = TRUE)
    max_y <- max(y, na.rm = TRUE)
    rx <- max(max_x - min_x, 1)
    ry <- max(max_y - min_y, 1)
    num_bins <- 160L
    grid_side <- num_bins + 1L
    bx <- pmax(0L, pmin(num_bins, floor(((x - min_x) / rx) * num_bins)))
    by <- pmax(0L, pmin(num_bins, floor(((y - min_y) / ry) * num_bins)))
    grid <- tabulate(by * grid_side + bx + 1L, nbins = grid_side * grid_side)
    radius <- 5L
    sigma <- 2
    offsets <- expand.grid(dx = -radius:radius, dy = -radius:radius)
    offsets$weight <- exp(-(offsets$dx^2 + offsets$dy^2) / (2 * sigma^2))
    densities <- numeric(n)
    for (i in seq_len(n)) {
        nx <- bx[i] + offsets$dx
        ny <- by[i] + offsets$dy
        ok <- nx >= 0L & nx <= num_bins & ny >= 0L & ny <= num_bins
        grid_idx <- ny[ok] * grid_side + nx[ok] + 1L
        densities[i] <- sum(grid[grid_idx] * offsets$weight[ok])
    }
    densities
}

.reference_pretty_k_label <- function(x) {
    x <- as.numeric(x)
    out <- ifelse(is.finite(x), ifelse(abs(x) >= 1000, paste0(round(x / 1000), "K"), as.character(round(x))), "")
    million <- is.finite(x) & abs(x) >= 1000000
    out[million] <- paste0(sub("\\.0$", "", sprintf("%.1f", x[million] / 1000000)), "M")
    out
}

.reference_gui_plot_aspect <- function() {
    (520 - 54 - 18) / (420 - 18 - 46)
}

.reference_gui_coord_ratio <- function(x_lim, y_lim) {
    x_span <- diff(as.numeric(x_lim[1:2]))
    y_span <- diff(as.numeric(y_lim[1:2]))
    if (!is.finite(x_span) || !is.finite(y_span) || x_span <= 0 || y_span <= 0) return(1)
    x_span / (y_span * .reference_gui_plot_aspect())
}

.reference_qc_scatter_density_plot <- function(data,
                                               x_channel,
                                               y_channel,
                                               pd,
                                               gate_vertices = NULL,
                                               title,
                                               subtitle = NULL,
                                               max_points = 50000L,
                                               point_size = 1.5,
                                               x_domain = NULL,
                                               y_domain = NULL,
                                               sample_type = NULL) {
    x_desc <- .get_reference_axis_label(x_channel, pd)
    y_desc <- .get_reference_axis_label(y_channel, pd)
    x_vals <- data[, x_channel]
    y_vals <- data[, y_channel]
    keep <- is.finite(x_vals) & is.finite(y_vals)
    x_vals <- x_vals[keep]
    y_vals <- y_vals[keep]
    idx <- .reference_even_indices(length(x_vals), max_points = max_points)
    plot_df <- data.frame(x = x_vals[idx], y = y_vals[idx])
    density <- .reference_point_density(plot_df$x, plot_df$y)
    palette <- .reference_density_palette()
    density[!is.finite(density)] <- 0
    max_density <- if (length(density) > 0) max(density, na.rm = TRUE) else 0
    bucket <- if (length(density) > 0 && is.finite(max_density) && max_density > 0) {
        pmax(1L, pmin(length(palette), floor((density / max_density) * (length(palette) - 1L)) + 1L))
    } else {
        rep(1L, nrow(plot_df))
    }
    bucket[!is.finite(bucket) | is.na(bucket)] <- 1L
    plot_df$color <- palette[bucket]
    x_lim <- .reference_gui_extent(x_vals)
    y_lim <- .reference_gui_extent(y_vals)
    if (!is.null(x_domain) && length(x_domain) >= 2L && all(is.finite(x_domain[1:2])) && x_domain[2] > x_domain[1]) {
        x_lim <- as.numeric(x_domain[1:2])
    }
    if (!is.null(y_domain) && length(y_domain) >= 2L && all(is.finite(y_domain[1:2])) && y_domain[2] > y_domain[1]) {
        y_lim <- as.numeric(y_domain[1:2])
    }
    coord_ratio <- .reference_gui_coord_ratio(x_lim, y_lim)
    point_size <- max(0.3, min(1.8, as.numeric(point_size) * 0.55))
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x, y)) +
        ggplot2::geom_point(color = plot_df$color, size = point_size, alpha = 0.95, stroke = 0) +
        ggplot2::labs(title = title, subtitle = subtitle, x = x_desc, y = y_desc) +
        ggplot2::scale_x_continuous(labels = .reference_pretty_k_label) +
        ggplot2::scale_y_continuous(labels = .reference_pretty_k_label) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            legend.position = "none",
            panel.grid = ggplot2::element_blank(),
            panel.background = ggplot2::element_rect(fill = "#fbfaf6", color = "#ded9cf", linewidth = 0.35),
            axis.line = ggplot2::element_line(color = "#6c746f", linewidth = 0.35),
            axis.ticks = ggplot2::element_line(color = "#6c746f", linewidth = 0.3),
            axis.title = ggplot2::element_text(color = "#626a65", face = "bold"),
            axis.text = ggplot2::element_text(color = "#626a65", face = "bold"),
            plot.title = ggplot2::element_text(color = "#17201c", face = "plain", size = 14),
            plot.subtitle = ggplot2::element_text(color = "#69716b", size = 10.6)
        ) +
        ggplot2::coord_fixed(ratio = coord_ratio, xlim = x_lim, ylim = y_lim, expand = FALSE)
    if (!is.null(gate_vertices) && nrow(gate_vertices) >= 3) {
        is_bead_gate <- identical(tolower(trimws(as.character(sample_type)[1])), "beads")
        gate_color <- if (is_bead_gate) "#56b4e9" else "#d65238"
        gate_alpha <- if (is_bead_gate) 0.16 else 0.16
        gate_closed <- rbind(gate_vertices, gate_vertices[1, , drop = FALSE])
        p <- p +
            ggplot2::geom_polygon(data = gate_vertices, ggplot2::aes(x, y), inherit.aes = FALSE, fill = gate_color, alpha = gate_alpha) +
            ggplot2::geom_path(data = gate_closed, ggplot2::aes(x, y), inherit.aes = FALSE, color = gate_color, linewidth = 1.05) +
            ggplot2::geom_point(data = gate_vertices, ggplot2::aes(x, y), inherit.aes = FALSE, shape = 21, size = 2.4, stroke = 0.9, color = "#d65238", fill = "#fbfaf6")
    }
    p
}

.reference_histogram_auto_cofactor <- function(values) {
    finite <- sort(as.numeric(values[is.finite(values)]))
    if (length(finite) < 2L) return(1)
    q <- as.numeric(stats::quantile(finite, c(0.001, 0.05, 0.10, 0.999), names = FALSE, na.rm = TRUE))
    robust_span <- max(q[4] - q[1], .Machine$double.eps)
    lower_scale <- max(abs(q[2]), abs(q[3]), abs(q[3] - q[2]), robust_span / 1000)
    cofactor <- lower_scale / sinh(2.5)
    max(robust_span / 1e6, min(robust_span / 2, cofactor))
}

.reference_histogram_transform_values <- function(values, transform = "auto", cofactor = NULL) {
    transform <- tolower(trimws(as.character(transform)[1]))
    values <- as.numeric(values)
    if (identical(transform, "auto")) {
        if (is.null(cofactor) || !is.finite(cofactor) || cofactor <= 0) {
            cofactor <- .reference_histogram_auto_cofactor(values)
        }
        return(asinh(values / cofactor))
    }
    if (identical(transform, "log10")) {
        return(log10(pmax(values, 1)))
    }
    if (identical(transform, "asinh")) {
        return(asinh(values / 150))
    }
    if (identical(transform, "biexponential")) {
        return(sign(values) * log10(1 + abs(values) / 50))
    }
    values
}

.reference_histogram_inverse_values <- function(values, transform = "auto", cofactor = NULL) {
    transform <- tolower(trimws(as.character(transform)[1]))
    values <- as.numeric(values)
    if (identical(transform, "auto")) {
        if (is.null(cofactor) || !is.finite(cofactor) || cofactor <= 0) cofactor <- 1
        return(sinh(values) * cofactor)
    }
    if (identical(transform, "log10")) {
        return(10^values)
    }
    if (identical(transform, "asinh")) {
        return(sinh(values) * 150)
    }
    if (identical(transform, "biexponential")) {
        return(sign(values) * 50 * (10^abs(values) - 1))
    }
    values
}

.reference_gui_extent <- function(values) {
    finite <- as.numeric(values[is.finite(values)])
    if (length(finite) == 0L) return(c(0, 1))
    val_max <- as.numeric(stats::quantile(finite, 0.998, na.rm = TRUE, names = FALSE))
    if (!is.finite(val_max) || val_max <= 0) val_max <- max(finite, na.rm = TRUE)
    if (!is.finite(val_max) || val_max <= 0) val_max <- 1
    c(0, val_max + val_max * 0.04)
}

.reference_gui_ticks <- function(domain, count = 5L) {
    domain <- as.numeric(domain[1:2])
    if (!all(is.finite(domain))) return(c(0, 1))
    if (domain[1] == domain[2]) return(domain[1])
    seq(domain[1], domain[2], length.out = count)
}

.reference_gui_domain_for_channel <- function(domains, channel, fallback_values = NULL) {
    if (!is.null(domains) && !is.null(channel) && nzchar(channel) && !is.null(domains[[channel]])) {
        domain <- as.numeric(domains[[channel]][1:2])
        if (length(domain) >= 2L && all(is.finite(domain)) && domain[2] > domain[1]) return(domain)
    }
    .reference_gui_extent(fallback_values)
}

.reference_gui_channel_extent <- function(raw_data, channel, filename, max_points = 100000L) {
    if (is.null(channel) || !nzchar(channel) || !channel %in% colnames(raw_data)) return(c(0, 1))
    idx <- .reference_gui_payload_indices(filename, nrow(raw_data), max_points = max_points)
    .reference_gui_extent(as.numeric(raw_data[idx, channel]))
}

.reference_gui_payload_indices <- function(filename, n, max_points = 100000L) {
    if (!is.finite(n) || n <= 0L) return(integer())
    max_points <- suppressWarnings(as.integer(max_points[1]))
    if (!is.finite(max_points) || is.na(max_points) || max_points <= 0L || n <= max_points) return(seq_len(n))
    seed <- sum(utf8ToInt(basename(filename))) %% .Machine$integer.max
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
        if (is.null(old_seed)) {
            if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
        } else {
            assign(".Random.seed", old_seed, envir = .GlobalEnv)
        }
    }, add = TRUE)
    set.seed(seed)
    sort(sample.int(n, max_points))
}

.reference_gui_display_data <- function(raw_data,
                                        filename,
                                        channels,
                                        display_max_points,
                                        preload_max_points = 100000L) {
    if (is.null(raw_data) || nrow(raw_data) == 0L) return(raw_data[0, , drop = FALSE])
    idx <- .reference_gui_payload_indices(filename, nrow(raw_data), max_points = preload_max_points)
    payload <- raw_data[idx, , drop = FALSE]
    channels <- unique(stats::na.omit(as.character(channels)))
    channels <- intersect(channels[nzchar(channels)], colnames(payload))
    if (length(channels) > 0L) {
        payload <- payload[stats::complete.cases(payload[, channels, drop = FALSE]), , drop = FALSE]
    }
    display_idx <- .reference_even_indices(nrow(payload), max_points = display_max_points)
    payload[display_idx, , drop = FALSE]
}

.reference_scatter_channels <- function(col_names) {
    unique(grep("^(FSC|SSC).*-(A|H|W)$", col_names, value = TRUE, ignore.case = TRUE))
}

.reference_compute_gui_scatter_domains <- function(fcs_files, max_points = 100000L) {
    domains <- list()
    for (fcs_file in fcs_files) {
        ff <- tryCatch(.spectreasy_read_fcs(fcs_file, label = "SCC FCS file"), error = function(e) NULL)
        if (is.null(ff)) next
        expr <- flowCore::exprs(ff)
        idx <- .reference_gui_payload_indices(basename(fcs_file), nrow(expr), max_points = max_points)
        channels <- .reference_scatter_channels(colnames(expr))
        for (channel in channels) {
            if (!channel %in% colnames(expr)) next
            values <- as.numeric(expr[idx, channel])
            domain <- .reference_gui_extent(values[is.finite(values)])
            if (is.null(domains[[channel]])) {
                domains[[channel]] <- domain
            } else {
                domains[[channel]] <- c(
                    min(domains[[channel]][1], domain[1], na.rm = TRUE),
                    max(domains[[channel]][2], domain[2], na.rm = TRUE)
                )
            }
        }
    }
    domains
}

.reference_histogram_density_curve <- function(values, domain, bins = 100L) {
    bins <- suppressWarnings(as.integer(bins))
    if (!is.finite(bins) || is.na(bins)) bins <- 100L
    bins <- min(max(bins, 5L), 500L)
    values <- values[is.finite(values)]
    if (length(values) == 0 || !all(is.finite(domain)) || domain[2] <= domain[1]) {
        return(data.frame(x = numeric(), y = numeric()))
    }
    idx <- pmax(0L, pmin(bins, floor(((values - domain[1]) / (domain[2] - domain[1])) * bins)))
    counts <- tabulate(idx + 1L, nbins = bins + 1L)
    smooth_radius <- max(1L, min(3L, floor(bins / 40L)))
    y <- numeric(bins + 1L)
    for (i in 0:bins) {
        js <- (-smooth_radius):smooth_radius
        bin_idx <- i + js
        ok <- bin_idx >= 0L & bin_idx <= bins
        weights <- exp(-(js[ok] * js[ok]) / 4)
        y[i + 1L] <- sum(counts[bin_idx[ok] + 1L] * weights) / max(sum(weights), 1)
    }
    data.frame(x = seq(domain[1], domain[2], length.out = bins + 1L), y = y)
}

.reference_histogram_qc_bounds <- function(vals_log, hist_info, gate_min, gate_max) {
    scalar_finite <- function(x) length(x) == 1L && is.finite(x)
    positive_min <- if (!is.null(hist_info$positive_raw_min)) hist_info$positive_raw_min else attr(vals_log, "pos_raw_min")
    positive_max <- if (!is.null(hist_info$positive_raw_max)) hist_info$positive_raw_max else attr(vals_log, "pos_raw_max")
    if (!scalar_finite(positive_min)) positive_min <- gate_min
    if (!scalar_finite(positive_max)) positive_max <- gate_max
    negative_min <- if (!is.null(hist_info$negative_raw_min)) hist_info$negative_raw_min else attr(vals_log, "neg_raw_min")
    negative_max <- if (!is.null(hist_info$negative_raw_max)) hist_info$negative_raw_max else attr(vals_log, "neg_raw_max")
    if (!scalar_finite(negative_min) && scalar_finite(attr(vals_log, "neg_log_min"))) {
        negative_min <- 10^attr(vals_log, "neg_log_min")
    }
    if (!scalar_finite(negative_max) && scalar_finite(attr(vals_log, "neg_log_max"))) {
        negative_max <- 10^attr(vals_log, "neg_log_max")
    }
    list(
        negative = c(negative_min, negative_max),
        positive = c(positive_min, positive_max),
        negative_present = isTRUE(attr(vals_log, "negative_gate_present")),
        positive_present = isTRUE(attr(vals_log, "positive_gate_present"))
    )
}

.reference_add_histogram_qc_gate <- function(plot, bounds, present, transform, ymax, color, cofactor = NULL) {
    if (!isTRUE(present) || length(bounds) != 2L || !all(is.finite(bounds)) || bounds[2] <= bounds[1]) return(plot)
    transformed <- .reference_histogram_transform_values(bounds, transform = transform, cofactor = cofactor)
    handles <- data.frame(x = transformed, y = rep(ymax * 0.5, length(transformed)))
    plot +
        ggplot2::annotate("rect", xmin = min(transformed), xmax = max(transformed), ymin = -Inf, ymax = Inf, alpha = 0.13, fill = color) +
        ggplot2::geom_vline(xintercept = transformed, color = color, linewidth = 1) +
        ggplot2::geom_point(
            data = handles, ggplot2::aes(x, y), inherit.aes = FALSE,
            shape = 21, size = 2.4, stroke = 0.9, color = color, fill = "#fbfaf6"
        )
}

.reference_qc_histogram_gui_plot <- function(peak_vals,
                                             pd,
                                             peak_channel,
                                             vals_log,
                                             gate_min,
                                             gate_max,
                                             settings,
                                             hist_info = NULL,
                                             x_domain = NULL) {
    transform <- settings$histogram_transform
    bins <- settings$histogram_bins
    cofactor <- if (identical(transform, "auto")) .reference_histogram_auto_cofactor(peak_vals) else NULL
    values_t <- .reference_histogram_transform_values(peak_vals, transform = transform, cofactor = cofactor)
    bounds <- .reference_histogram_qc_bounds(vals_log, hist_info, gate_min, gate_max)
    raw_domain <- if (!is.null(x_domain) && length(x_domain) >= 2L && all(is.finite(x_domain[1:2])) && x_domain[2] > x_domain[1]) {
        as.numeric(x_domain[1:2])
    } else {
        .reference_gui_extent(peak_vals)
    }
    domain <- .reference_histogram_transform_values(raw_domain, transform = transform, cofactor = cofactor)
    if (!all(is.finite(domain)) || domain[2] <= domain[1]) domain <- c(0, 1)
    x_breaks <- .reference_gui_ticks(domain, 5L)
    x_labels <- .reference_pretty_k_label(.reference_histogram_inverse_values(x_breaks, transform = transform, cofactor = cofactor))
    curve <- .reference_histogram_density_curve(values_t, domain, bins = bins)
    ymax <- max(curve$y, 1, na.rm = TRUE)
    y_lim <- c(0, ymax * 1.12)
    coord_ratio <- .reference_gui_coord_ratio(domain, y_lim)
    x_desc <- .get_reference_axis_label(peak_channel, pd)
    p <- ggplot2::ggplot(curve, ggplot2::aes(x, y)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = y), fill = "#263f73", alpha = 0.28) +
        ggplot2::geom_line(color = "#263f73", linewidth = 0.85) +
        ggplot2::labs(title = "Histogram", x = x_desc, y = "Events per bin") +
        ggplot2::scale_x_continuous(breaks = x_breaks, labels = x_labels) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            legend.position = "none",
            panel.grid = ggplot2::element_blank(),
            panel.background = ggplot2::element_rect(fill = "#fbfaf6", color = "#ded9cf", linewidth = 0.35),
            axis.line = ggplot2::element_line(color = "#6c746f", linewidth = 0.35),
            axis.ticks = ggplot2::element_line(color = "#6c746f", linewidth = 0.3),
            axis.title = ggplot2::element_text(color = "#626a65", face = "bold"),
            axis.text = ggplot2::element_text(color = "#626a65", face = "bold"),
            plot.title = ggplot2::element_text(color = "#17201c", face = "plain", size = 14),
            plot.margin = ggplot2::margin(5.5, 24, 5.5, 5.5)
        ) +
        ggplot2::coord_fixed(ratio = coord_ratio, xlim = domain, ylim = y_lim, expand = FALSE, clip = "off")
    p <- .reference_add_histogram_qc_gate(p, bounds$negative, bounds$negative_present, transform, ymax, "#263f73", cofactor = cofactor)
    .reference_add_histogram_qc_gate(p, bounds$positive, bounds$positive_present, transform, ymax, "#d65238", cofactor = cofactor)
}

# Generates and saves PDF quality control (QC) plots for a processed sample.
# Produces three plots: a 2D density plot of FSC-SSC gating, a 1D density histogram of positive/negative
# peak gating, and a spectral profile plot of the normalized signature.
# Writes files to the output folder.
.save_reference_qc_plots <- function(sn,
                                     raw_data,
                                     gated_data,
                                     final_gated_data,
                                     pd,
                                     fsc,
                                     ssc,
                                     fsc_max,
                                     ssc_max,
                                     final_gate,
                                     vals_log,
                                     peak_vals,
                                     gate_min,
                                     gate_max,
                                     peak_channel,
                                     detector_names,
                                     detector_labels,
                                     det_info,
                                     out_path,
                                     manual_gate_info = NULL,
                                     hist_info = NULL,
                                     gui_scatter_domains = NULL,
                                     histogram_domain = NULL,
                                     source_filename = NULL,
                                     sample_type = NULL) {
    context <- .reference_qc_plot_context(
        raw_data = raw_data,
        fsc = fsc,
        ssc = ssc,
        peak_channel = peak_channel,
        final_gate = final_gate,
        manual_gate_info = manual_gate_info,
        source_filename = source_filename,
        sn = sn
    )
    histogram_source <- .reference_qc_cell_source(
        context = context,
        raw_data = raw_data,
        gated_data = gated_data,
        pd = pd,
        gui_scatter_domains = gui_scatter_domains,
        sample_type = sample_type,
        sn = sn,
        out_path = out_path
    )
    histogram_source <- .reference_qc_singlet_source(
        context = context,
        histogram_source = histogram_source,
        raw_data = raw_data,
        gated_data = gated_data,
        pd = pd,
        gui_scatter_domains = gui_scatter_domains,
        sample_type = sample_type,
        sn = sn,
        out_path = out_path
    )
    spectrum_data <- .reference_qc_histogram_stage(
        histogram_source = histogram_source,
        final_gated_data = final_gated_data,
        pd = pd,
        peak_channel = peak_channel,
        vals_log = vals_log,
        gate_min = gate_min,
        gate_max = gate_max,
        gate_settings = context$settings,
        hist_info = hist_info,
        histogram_domain = histogram_domain,
        sn = sn,
        out_path = out_path
    )
    .save_reference_spectrum_qc(
        data = spectrum_data,
        detector_names = detector_names,
        detector_labels = detector_labels,
        det_info = det_info,
        sn = sn,
        out_path = out_path
    )
    invisible(NULL)
}

# Processes a single FCS file through the entire reference construction pipeline.
# Reads the file, performs scatter gating, identifies the peak channel, runs histogram gating
# to separate positive and negative events, calculates the normalized spectrum, and saves QC plots.
# Returns a list containing the processed spectrum, QC summary row, and sample name metadata.
