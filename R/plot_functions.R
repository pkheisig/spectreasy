.natural_reference_row_order <- function(x) {
    x <- as.character(x)
    af_match <- grepl("^AF(?:_(\\d+))?$", x, ignore.case = TRUE, perl = TRUE)
    if (!any(af_match)) {
        return(seq_along(x))
    }
    af_names <- x[af_match]
    af_suffix <- suppressWarnings(as.integer(sub("^AF_(\\d+)$", "\\1", af_names, ignore.case = TRUE, perl = TRUE)))
    af_suffix[is.na(af_suffix)] <- 1L
    af_order <- order(af_suffix, af_names, na.last = TRUE)
    out <- seq_along(x)
    out[which(af_match)] <- which(af_match)[af_order]
    out
}

.natural_reference_row_levels <- function(x) {
    x <- as.character(x)
    x[.natural_reference_row_order(x)]
}

.abbreviate_reference_peak_label <- function(x) {
    x <- trimws(as.character(x))
    x[is.na(x)] <- ""
    out <- x

    out <- gsub("\\bAlexa\\s+Fluor\\s*([0-9]+)\\b", "Alexa\\1", out, ignore.case = TRUE, perl = TRUE)
    out <- gsub("\\bBrilliant\\s+Ultra\\s+Violet\\s*([0-9]+)\\b", "BUV\\1", out, ignore.case = TRUE, perl = TRUE)
    out <- gsub("\\bBrilliant\\s+Violet\\s*([0-9]+)\\b", "BV\\1", out, ignore.case = TRUE, perl = TRUE)
    out <- gsub("\\bBrilliant\\s+Blue\\s*([0-9]+)\\b", "BB\\1", out, ignore.case = TRUE, perl = TRUE)
    out <- gsub("\\beFluor\\s*([0-9]+)\\b", "eF\\1", out, ignore.case = TRUE, perl = TRUE)
    out <- gsub("\\bLIVE\\s*/?\\s*DEAD\\b", "LD", out, ignore.case = TRUE, perl = TRUE)
    out <- gsub("\\s+", "", out, perl = TRUE)
    out
}

.build_reference_peak_label_df <- function(ref_matrix,
                                           detectors,
                                           labels,
                                           fluor_levels,
                                           min_peak_intensity = 0.2) {
    if (nrow(ref_matrix) == 0L || ncol(ref_matrix) == 0L) {
        return(data.frame())
    }
    mat <- ref_matrix[fluor_levels, detectors, drop = FALSE]
    rows <- lapply(seq_len(nrow(mat)), function(i) {
        vals <- as.numeric(mat[i, ])
        vals[!is.finite(vals)] <- -Inf
        if (!any(is.finite(vals)) || max(vals, na.rm = TRUE) < min_peak_intensity) {
            return(NULL)
        }
        peak_idx <- which.max(vals)
        data.frame(
            Fluorophore = rownames(mat)[i],
            Detector = detectors[peak_idx],
            DetectorLabel = labels[peak_idx],
            Intensity = vals[peak_idx],
            Label = .abbreviate_reference_peak_label(rownames(mat)[i]),
            stringsAsFactors = FALSE
        )
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (!length(rows)) {
        return(data.frame())
    }
    out <- do.call(rbind, rows)
    out$Fluorophore <- factor(out$Fluorophore, levels = fluor_levels)
    out$Detector <- factor(out$Detector, levels = detectors, labels = labels)

    split_idx <- split(seq_len(nrow(out)), as.character(out$Detector))
    out$LabelY <- pmax(out$Intensity + 0.035, 1.035)
    for (idx in split_idx) {
        if (length(idx) > 1L) {
            idx <- idx[order(as.integer(out$Fluorophore[idx]))]
            out$LabelY[idx] <- 1.035 + (seq_along(idx) - 1L) * 0.055
        }
    }
    out
}

#' Plot Spectral Overlays
#' 
#' @param ref_matrix Reference matrix (Markers x Detectors)
#' @param pd Optional pData for descriptive labels
#' @param output_file Optional path to save the plot. Set `NULL` to return the plot without writing a file.
#' @param width Plot width.
#' @param height Plot height.
#' @param unit Plot size unit.
#' @param dpi Output resolution.
#' @param theme_custom Optional ggplot theme object (reserved).
#' @param annotate_peaks `"auto"` (default), `TRUE`, or `FALSE`. In auto mode,
#'   compact peak labels are added for regular marker spectra but skipped for
#'   all-AF-bank overlays.
#' @param peak_label_max Maximum number of spectra to label in auto mode.
#' @param peak_label_size Text size for peak labels.
#' @return A `ggplot` object.
#' @export
#' @examples
#' m <- diag(3)
#' rownames(m) <- c("FITC", "PE", "APC")
#' colnames(m) <- c("B2-A", "YG1-A", "R1-A")
#' p <- plot_spectra(ref_matrix = m, output_file = NULL)
#' print(p)
plot_spectra <- function(ref_matrix,
                         pd = NULL,
                         output_file = NULL,
                         width = 250,
                         height = 100,
                         unit = "mm",
                         dpi = 600,
                         theme_custom = NULL,
                         annotate_peaks = "auto",
                         peak_label_max = 40L,
                         peak_label_size = 2.6) {
    ref_matrix <- .as_reference_matrix(ref_matrix, "ref_matrix")
    detectors <- colnames(ref_matrix)
    
    # 1. Get Sorted Detectors and Labels
    if (!is.null(pd)) {
        det_info <- get_sorted_detectors(pd)
        # Filter ref_matrix to match
        common <- intersect(det_info$names, detectors)
        
        if (length(common) == 0) {
            warning("No matching detectors found between reference matrix and provided metadata. Ignoring metadata.")
            message("Ref Matrix cols (first 5): ", paste(utils::head(detectors, 5), collapse=", "))
            message("Metadata names (first 5): ", paste(utils::head(det_info$names, 5), collapse=", "))

            # Fallback to numerical sort
            nums <- as.numeric(gsub("[^0-9]", "", detectors))
            ord <- order(nums)
            ref_matrix <- ref_matrix[, ord, drop = FALSE]
            detectors <- colnames(ref_matrix)
            labels <- detectors
        } else {
            ref_matrix <- ref_matrix[, common, drop = FALSE]
            detectors <- common
            # Get labels only for common
            labels <- det_info$labels[match(common, det_info$names)]
        }
    } else {
        # Fallback: sort by laser group first, then detector number.
        # Desired order: UV, V, B, YG, R.
        detector_key <- toupper(gsub("\\s+", "", gsub("-A$", "", detectors)))
        laser_group <- vapply(detector_key, function(k) {
            if (grepl("^UV", k)) return(1L)
            if (grepl("^V", k)) return(2L)
            if (grepl("^B", k)) return(3L)
            if (grepl("^YG", k) || grepl("^Y", k) || grepl("^G", k)) return(4L)
            if (grepl("^R", k)) return(5L)
            return(99L)
        }, integer(1))

        det_num <- suppressWarnings(as.integer(sub("^[A-Z]+([0-9]+).*$", "\\1", detector_key)))
        det_num[!is.finite(det_num)] <- 999L

        ord <- order(laser_group, det_num, detectors)
        ref_matrix <- ref_matrix[, ord, drop = FALSE]
        detectors <- colnames(ref_matrix)
        labels <- detectors
    }

    long <- data.frame(
        Fluorophore = rep(rownames(ref_matrix), ncol(ref_matrix)),
        Detector = rep(detectors, each = nrow(ref_matrix)),
        Intensity = as.vector(ref_matrix)
    )
    
    long$Detector <- factor(long$Detector, levels = detectors, labels = labels)
    fluor_levels <- .natural_reference_row_levels(rownames(ref_matrix))
    long$Fluorophore <- factor(long$Fluorophore, levels = fluor_levels)

    n_fluor <- nrow(ref_matrix)
    legend_cols <- if (n_fluor > 16) 3L else if (n_fluor > 8) 2L else 1L
    all_af_rows <- all(grepl("^AF(?:_\\d+)?$", rownames(ref_matrix), ignore.case = TRUE, perl = TRUE))
    annotate_mode <- annotate_peaks
    if (is.logical(annotate_mode)) {
        do_annotate_peaks <- isTRUE(annotate_mode)
    } else {
        annotate_mode <- match.arg(as.character(annotate_mode)[1], c("auto", "always", "never"))
        do_annotate_peaks <- if (identical(annotate_mode, "always")) {
            TRUE
        } else if (identical(annotate_mode, "never")) {
            FALSE
        } else {
            n_fluor <= peak_label_max && !all_af_rows
        }
    }

    p <- ggplot2::ggplot(long, ggplot2::aes(Detector, Intensity, color = Fluorophore, group = Fluorophore)) +
        ggplot2::geom_line(linewidth = 0.7) +
        ggplot2::scale_y_continuous(
            breaks = seq(0, 1, by = 0.2),
            labels = sprintf("%.1f", seq(0, 1, by = 0.2))
        ) +
        ggplot2::theme_minimal(base_size = 13.75) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(size = 6.25, angle = 90, hjust = 1, vjust = 0.5),
            legend.text = ggplot2::element_text(size = 7.5),
            legend.key.size = ggplot2::unit(4, "mm"),
            legend.spacing.y = ggplot2::unit(0.5, "mm"),
            plot.margin = ggplot2::margin(8, 14, 6, 6)
        ) +
        ggplot2::guides(color = ggplot2::guide_legend(ncol = legend_cols)) +
        ggplot2::labs(
            title = "Reference Spectra Overlay",
            x = "Detector",
            y = "Normalized Intensity"
        )

    if (isTRUE(do_annotate_peaks)) {
        peak_df <- .build_reference_peak_label_df(
            ref_matrix = ref_matrix,
            detectors = detectors,
            labels = labels,
            fluor_levels = fluor_levels
        )
        if (nrow(peak_df) > 0L) {
            y_max <- max(c(long$Intensity, peak_df$LabelY), na.rm = TRUE)
            p <- p +
                ggplot2::geom_text(
                    data = peak_df,
                    ggplot2::aes(x = Detector, y = LabelY, label = Label),
                    inherit.aes = FALSE,
                    color = "black",
                    angle = 45,
                    hjust = 0,
                    vjust = 0.5,
                    size = peak_label_size
                ) +
                ggplot2::coord_cartesian(ylim = c(0, y_max + 0.04), clip = "off")
        }
    }

    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, p, width = width, height = height, units = unit, dpi = dpi)
    }
    return(p)
}


.extract_unmix_scatter_data_list <- function(unmixed_list) {
    if (!is.list(unmixed_list) || length(unmixed_list) == 0) {
        stop("unmixed_list must be a non-empty list from unmix_samples().")
    }

    extract_data <- function(x) {
        if (is.list(x) && "data" %in% names(x) && is.data.frame(x$data)) return(x$data)
        if (is.data.frame(x)) return(x)
        NULL
    }

    data_list <- lapply(unmixed_list, extract_data)
    valid <- !vapply(data_list, is.null, logical(1))
    data_list <- data_list[valid]
    if (length(data_list) == 0) stop("No valid data frames found in unmixed_list.")

    sample_ids <- names(data_list)
    if (is.null(sample_ids) || any(sample_ids == "")) {
        sample_ids <- paste0("sample_", seq_along(data_list))
    }
    names(data_list) <- sample_ids
    data_list
}

.normalize_unmix_scatter_mapping <- function(sample_ids, sample_to_marker = NULL) {
    normalize_id <- function(x) tools::file_path_sans_ext(basename(as.character(x)))
    sample_keys <- normalize_id(sample_ids)

    if (is.null(sample_to_marker)) {
        sample_to_marker <- stats::setNames(sample_ids, sample_ids)
    }
    marker_names_raw <- names(sample_to_marker)
    sample_to_marker <- as.character(sample_to_marker)
    names(sample_to_marker) <- marker_names_raw
    names(sample_to_marker) <- normalize_id(names(sample_to_marker))
    sample_to_marker <- sample_to_marker[!is.na(names(sample_to_marker)) & names(sample_to_marker) != ""]
    sample_to_marker <- sample_to_marker[!duplicated(names(sample_to_marker))]

    sample_stains <- sample_to_marker[sample_keys]
    names(sample_stains) <- sample_ids
    sample_stains
}

.resolve_unmix_scatter_markers <- function(sample_stains, markers = NULL, marker_display = NULL) {
    if (is.null(markers)) {
        markers <- unique(sample_stains[!is.na(sample_stains) & sample_stains != ""])
    } else {
        markers <- unique(as.character(markers))
    }
    markers <- markers[!is.na(markers) & markers != ""]
    if (length(markers) < 2) stop("Need at least two marker names for scatter matrix.")

    marker_labels <- stats::setNames(markers, markers)
    if (!is.null(marker_display)) {
        md <- as.character(marker_display)
        names(md) <- names(marker_display)
        md <- md[!is.na(names(md)) & names(md) != ""]
        hits <- intersect(names(md), markers)
        if (length(hits) > 0) {
            marker_labels[hits] <- md[hits]
        }
    }

    list(markers = markers, marker_labels = marker_labels)
}

.compute_unmix_scatter_limits <- function(v) {
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
        pad <- max(1e-6, abs(lo) * 0.1)
        lo <- lo - pad
        hi <- hi + pad
    } else {
        pad <- (hi - lo) * 0.05
        lo <- lo - pad
        hi <- hi + pad
    }
    c(lo, hi)
}

.add_unmix_scatter_density_colors <- function(plot_df) {
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

.compute_unmix_fixed_scatter_limits <- function(axis_limit,
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

.build_unmix_scatter_panel_data <- function(data_list,
                                            sample_stains,
                                            markers,
                                            max_points_per_sample = 3000,
                                            transform = c("none", "asinh"),
                                            asinh_cofactor = 150,
                                            axis_limit = NULL) {
    transform <- match.arg(transform)
    fixed_limits <- .compute_unmix_fixed_scatter_limits(
        axis_limit = axis_limit,
        transform = transform,
        asinh_cofactor = asinh_cofactor
    )
    panel_data <- list()
    panel_limits <- list()
    k <- 1
    lim_k <- 1

    for (i in seq_along(data_list)) {
        d <- data_list[[i]]
        row_stain <- sample_stains[[i]]
        if (is.na(row_stain) || row_stain == "" || !(row_stain %in% markers)) next
        if (!(row_stain %in% colnames(d))) next

        marker_idx <- match(row_stain, markers)
        if (is.na(marker_idx) || marker_idx <= 1) next
        x_markers <- markers[seq_len(marker_idx - 1)]
        x_markers <- x_markers[x_markers %in% colnames(d)]
        if (length(x_markers) == 0) next

        for (xm in x_markers) {
            d_pair <- d[, c(xm, row_stain), drop = FALSE]
            d_pair <- d_pair[stats::complete.cases(d_pair), , drop = FALSE]
            if (nrow(d_pair) == 0) next
            if (nrow(d_pair) > max_points_per_sample) {
                d_pair <- d_pair[sample.int(nrow(d_pair), max_points_per_sample), , drop = FALSE]
            }

            x_vals <- d_pair[[xm]]
            y_vals <- d_pair[[row_stain]]
            if (transform == "asinh") {
                x_vals <- asinh(x_vals / asinh_cofactor)
                y_vals <- asinh(y_vals / asinh_cofactor)
            }

            if (is.null(fixed_limits)) {
                x_lim <- .compute_unmix_scatter_limits(x_vals)
                y_lim <- .compute_unmix_scatter_limits(y_vals)
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
            k <- k + 1
            lim_k <- lim_k + 1
        }
    }

    if (length(panel_data) == 0) stop("No panel data available for scatter matrix plot.")

    plot_df <- do.call(rbind, panel_data)
    lim_df <- do.call(rbind, panel_limits)
    plot_df <- .add_unmix_scatter_density_colors(plot_df)
    plot_df$panel_col <- factor(plot_df$panel_col, levels = markers)
    plot_df$panel_row <- factor(plot_df$panel_row, levels = markers)
    lim_df$panel_col <- factor(lim_df$panel_col, levels = markers)
    lim_df$panel_row <- factor(lim_df$panel_row, levels = markers)

    list(plot_df = plot_df, lim_df = lim_df)
}

.build_unmix_scatter_plot <- function(plot_df, lim_df, markers, marker_labels) {
    panel_used <- unique(lim_df[, c("panel_row", "panel_col"), drop = FALSE])
    panel_grid <- expand.grid(
        panel_row = factor(markers, levels = markers),
        panel_col = factor(markers, levels = markers),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )

    ggplot2::ggplot() +
        ggplot2::geom_blank(data = panel_grid, ggplot2::aes(x = 0, y = 0)) +
        ggplot2::geom_blank(data = lim_df, ggplot2::aes(x = x_low, y = y_low)) +
        ggplot2::geom_blank(data = lim_df, ggplot2::aes(x = x_high, y = y_high)) +
        ggplot2::geom_hline(
            data = panel_used,
            ggplot2::aes(yintercept = 0),
            inherit.aes = FALSE,
            color = "grey45",
            linewidth = 0.1
        ) +
        ggplot2::geom_vline(
            data = panel_used,
            ggplot2::aes(xintercept = 0),
            inherit.aes = FALSE,
            color = "grey45",
            linewidth = 0.1
        ) +
        ggplot2::geom_point(
            data = plot_df,
            ggplot2::aes(x = x, y = y, color = color),
            alpha = 0.6,
            size = 0.25,
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
            title = "Unmixing Scatter Matrix",
            subtitle = "Good: row-stain events are high on Y and near zero on X (other markers).\nBad: large off-axis clouds indicate cross-talk, control mislabeling, or unstable unmixing."
        ) +
        ggplot2::theme_bw(base_size = 8.75) +
        ggplot2::theme(
            legend.position = "none",
            panel.grid = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(color = "grey60", linewidth = 0.25),
            strip.background = ggplot2::element_rect(fill = "grey95", color = "grey80"),
            strip.text = ggplot2::element_text(size = 7.5, face = "bold"),
            axis.title = ggplot2::element_blank(),
            axis.text = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            axis.line = ggplot2::element_blank(),
            panel.spacing = grid::unit(0.3, "mm"),
            plot.title = ggplot2::element_text(size = 12.5, face = "bold", hjust = 0.5),
            plot.subtitle = ggplot2::element_text(size = 8.4, hjust = 0.5, lineheight = 1.1)
        )
}

.resolve_unmix_scatter_plot_size <- function(markers, panel_size_mm = 30, max_size_mm = 1200) {
    n <- length(markers)
    panel_size_mm <- as.numeric(panel_size_mm)[1]
    if (!is.finite(panel_size_mm) || panel_size_mm <= 0) panel_size_mm <- 30
    max_size_mm <- as.numeric(max_size_mm)[1]
    if (!is.finite(max_size_mm) || max_size_mm <= 0) max_size_mm <- 1200

    width_mm <- max(220, panel_size_mm * n + 35)
    height_mm <- max(220, panel_size_mm * n + 30)
    scale <- min(1, max_size_mm / max(width_mm, height_mm))
    list(
        width_mm = width_mm * scale,
        height_mm = height_mm * scale
    )
}

#' Plot NxN Unmixing Scatter Matrix
#'
#' Creates a lower-triangle marker-vs-marker scatter matrix from unmixed SCC output.
#' Each populated panel uses events from only one single-stain file (row stain),
#' while the column stain is permuted across the lower triangle.
#'
#' @param unmixed_list Output list from `unmix_samples()`.
#' @param sample_to_marker Named character vector mapping sample IDs (basename
#'   without extension) to fluorophore/marker names.
#' @param markers Optional character vector of marker columns to include.
#' @param marker_display Optional named character vector mapping primary marker
#'   names to display labels (for example `"BUV395 / CD45RA"`).
#' @param output_file Optional path to save a PNG output. Set `NULL` to return the plot without writing a file.
#' @param max_points_per_sample Maximum events sampled per sample for plotting.
#' @param transform One of `"asinh"` or `"none"`.
#' @param asinh_cofactor Cofactor used when `transform = "asinh"`.
#' @param axis_limit Optional fixed symmetric scatter axis limit. The default
#'   `NULL` uses local per-panel ranges. Use `1e5` for `c(-1e5, 1e5)` on
#'   every panel.
#' @param panel_size_mm Size per matrix panel in millimeters.
#' @param seed Optional integer seed for deterministic point subsampling.
#' @return A `ggplot` object.
#' @export
#' @examples
#' M_demo <- rbind(
#'   FITC = c(1.00, 0.20, 0.05),
#'   PE = c(0.10, 1.00, 0.20),
#'   APC = c(0.05, 0.15, 1.00)
#' )
#' colnames(M_demo) <- c("B2-A", "YG1-A", "R1-A")
#'
#' simulate_sample <- function(dominant_marker, M, n_cells = 120) {
#'   markers <- rownames(M)
#'   marker_signal <- matrix(rexp(n_cells * length(markers), rate = 8), ncol = length(markers))
#'   colnames(marker_signal) <- markers
#'   marker_signal[, dominant_marker] <- rexp(n_cells, rate = 0.6) + 2
#'   raw_signal <- marker_signal %*% M + matrix(rnorm(n_cells * ncol(M), sd = 0.03), ncol = ncol(M))
#'   exprs_mat <- cbind(
#'     raw_signal,
#'     Time = seq_len(n_cells),
#'     "FSC-A" = rnorm(n_cells, mean = 90000, sd = 7000),
#'     "SSC-A" = rnorm(n_cells, mean = 45000, sd = 5000)
#'   )
#'   colnames(exprs_mat)[seq_len(ncol(M))] <- colnames(M)
#'   flowCore::flowFrame(exprs_mat)
#' }
#'
#' toy_fs <- flowCore::flowSet(list(
#'   FITC_sample = simulate_sample("FITC", M_demo),
#'   PE_sample = simulate_sample("PE", M_demo),
#'   APC_sample = simulate_sample("APC", M_demo)
#' ))
#' toy_unmixed <- unmix_samples(toy_fs, M = M_demo, method = "OLS", write_fcs = FALSE)
#'
#' p <- plot_unmixing_scatter_matrix(
#'   unmixed_list = toy_unmixed,
#'   sample_to_marker = c(FITC_sample = "FITC", PE_sample = "PE", APC_sample = "APC"),
#'   markers = rownames(M_demo),
#'   output_file = NULL,
#'   max_points_per_sample = 80,
#'   seed = 1
#' )
#' print(p)
plot_unmixing_scatter_matrix <- function(
    unmixed_list,
    sample_to_marker = NULL,
    markers = NULL,
    marker_display = NULL,
    output_file = NULL,
    max_points_per_sample = 1000,
    transform = c("none", "asinh"),
    asinh_cofactor = 150,
    axis_limit = NULL,
    panel_size_mm = 30,
    seed = NULL
) {
    transform <- match.arg(transform)
    .with_optional_seed(seed)

    data_list <- .extract_unmix_scatter_data_list(unmixed_list)
    sample_stains <- .normalize_unmix_scatter_mapping(names(data_list), sample_to_marker = sample_to_marker)
    marker_info <- .resolve_unmix_scatter_markers(sample_stains, markers = markers, marker_display = marker_display)
    panel_info <- .build_unmix_scatter_panel_data(
        data_list = data_list,
        sample_stains = sample_stains,
        markers = marker_info$markers,
        max_points_per_sample = max_points_per_sample,
        transform = transform,
        asinh_cofactor = asinh_cofactor,
        axis_limit = axis_limit
    )
    p <- .build_unmix_scatter_plot(
        plot_df = panel_info$plot_df,
        lim_df = panel_info$lim_df,
        markers = marker_info$markers,
        marker_labels = marker_info$marker_labels
    )

    plot_size <- .resolve_unmix_scatter_plot_size(marker_info$markers, panel_size_mm = panel_size_mm)
    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, p, width = plot_size$width_mm, height = plot_size$height_mm, units = "mm", dpi = 350)
    }

    return(p)
}
