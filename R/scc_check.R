# Internal helpers for SCC report generation.

.prepare_scc_report_plot_dir <- function(qc_plot_dir, save_qc_pngs = FALSE) {
    if (isTRUE(save_qc_pngs)) {
        plot_dir <- normalizePath(qc_plot_dir, mustWork = FALSE)
        dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
        return(list(plot_dir = plot_dir, cleanup_dir = NULL, retained_qc_plot_dir = plot_dir))
    }

    plot_dir <- tempfile("spectreasy_scc_report_plots_")
    dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
    list(plot_dir = plot_dir, cleanup_dir = plot_dir, retained_qc_plot_dir = NULL)
}

.read_report_png_grob <- function(path) {
    if (!file.exists(path)) {
        return(NULL)
    }
    img <- png::readPNG(path)
    grid::rasterGrob(img, interpolate = TRUE)
}

.draw_report_image_panel <- function(path, title, x, y, width, height) {
    grob <- .read_report_png_grob(path)
    grid::grid.text(
        title,
        x = x,
        y = y + height / 2 + grid::unit(4, "mm"),
        gp = grid::gpar(fontsize = 11, fontface = "bold")
    )
    if (is.null(grob)) {
        grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = "grey75", fill = NA))
        grid::grid.text("Plot not available", x = x, y = y, gp = grid::gpar(col = "grey40", fontsize = 10))
    } else {
        grid::grid.draw(
            grid::editGrob(
                grob,
                vp = grid::viewport(x = x, y = y, width = width, height = height)
            )
        )
    }
}

.format_scc_summary_lines <- function(df) {
    if (nrow(df) == 0) {
        return("No SCC summary rows available.")
    }
    if (!("marker" %in% colnames(df))) {
        df$marker <- ""
    }
    if (!("stain_index" %in% colnames(df))) {
        df$stain_index <- NA_real_
    }
    if (!("saturated" %in% colnames(df))) {
        df$saturated <- "OK"
    }
    df <- df[, c(
        "fluorophore",
        "marker",
        "sample",
        "type",
        "peak_channel",
        "n_total",
        "n_final",
        "scatter_gate_pct",
        "histogram_gate_pct",
        "stain_index",
        "saturated"
    ), drop = FALSE]

    headers <- c("Fluor", "Marker", "Sample", "Type", "Peak", "Events", "Final", "Scatter%", "Gate%", "SI", "Sat")
    widths <- c(14, 14, 20, 6, 8, 7, 7, 8, 7, 6, 4)

    trim_to <- function(x, w) {
        x <- as.character(x)
        x[is.na(x)] <- ""
        vapply(x, function(s) {
            if (nchar(s, type = "width") > w) {
                paste0(substr(s, 1, max(1, w - 3)), "...")
            } else {
                s
            }
        }, character(1))
    }

    fmt_num <- function(x, digits = 1) {
        out <- ifelse(is.na(x), "", format(round(as.numeric(x), digits), nsmall = digits, trim = TRUE))
        as.character(out)
    }

    vals <- list(
        trim_to(df$fluorophore, widths[1]),
        trim_to(df$marker, widths[2]),
        trim_to(df$sample, widths[3]),
        trim_to(df$type, widths[4]),
        trim_to(df$peak_channel, widths[5]),
        trim_to(df$n_total, widths[6]),
        trim_to(df$n_final, widths[7]),
        trim_to(fmt_num(df$scatter_gate_pct, 1), widths[8]),
        trim_to(fmt_num(df$histogram_gate_pct, 1), widths[9]),
        trim_to(fmt_num(df$stain_index, 1), widths[10]),
        trim_to(df$saturated, widths[11])
    )

    make_line <- function(parts) {
        paste(vapply(seq_along(parts), function(i) sprintf(paste0("%-", widths[i], "s"), parts[[i]]), character(1)), collapse = " ")
    }

    header_line <- make_line(as.list(headers))
    sep_line <- paste(vapply(widths, function(w) paste(rep("-", w), collapse = ""), character(1)), collapse = " ")
    row_lines <- vapply(seq_len(nrow(df)), function(i) make_line(lapply(vals, `[[`, i)), character(1))

    c(header_line, sep_line, row_lines)
}

.draw_report_ggplot_page <- function(p, square = FALSE, height_ratio = 1.0) {
    if (is.null(p)) {
        return(invisible(NULL))
    }
    grid::grid.newpage()
    grob <- ggplot2::ggplotGrob(p)
    if (!isTRUE(square) && height_ratio == 1.0) {
        grid::grid.draw(grob)
    } else {
        vp_h <- if (isTRUE(square)) 0.99 else height_ratio
        vp_w <- if (isTRUE(square)) 0.78 else 0.95
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
    }
    invisible(NULL)
}

.draw_scc_report_sample_page <- function(row, report_plot_dir, use_scatter_gating = TRUE) {
    sample_id <- row$sample[[1]]
    fluor <- row$fluorophore[[1]]
    marker <- if ("marker" %in% colnames(row)) trimws(as.character(row$marker[[1]])) else ""
    title <- if (!is.na(marker) && nzchar(marker) && tolower(marker) != tolower(fluor)) {
        paste0(fluor, " / ", marker, " (", sample_id, ")")
    } else {
        paste0(fluor, " (", sample_id, ")")
    }
    subtitle <- paste0(
        "Type: ", row$type[[1]],
        " | Peak channel: ", row$peak_channel[[1]],
        " | Events: ", row$n_total[[1]],
        " | Scatter gate: ", row$n_scatter_gated[[1]], " (", row$scatter_gate_pct[[1]], "%)",
        " | Final gate: ", row$n_final[[1]], " (", row$histogram_gate_pct[[1]], "%)"
    )

    grid::grid.newpage()
    grid::grid.text(title, x = 0.05, y = 0.96, just = c("left", "top"), gp = grid::gpar(fontsize = 16, fontface = "bold"))
    grid::grid.text(subtitle, x = 0.05, y = 0.92, just = c("left", "top"), gp = grid::gpar(fontsize = 10))

    .draw_report_image_panel(
        file.path(report_plot_dir, "fsc_ssc", paste0(sample_id, "_fsc_ssc.png")),
        "FSC/SSC Auto-Gate",
        x = grid::unit(0.27, "npc"),
        y = grid::unit(0.62, "npc"),
        width = grid::unit(0.46, "npc"),
        height = grid::unit(0.44, "npc")
    )
    gate_type <- if ("intensity_gate_type" %in% colnames(row)) {
        trimws(as.character(row$intensity_gate_type[[1]]))
    } else {
        ""
    }
    spectral_selection_path <- file.path(report_plot_dir, "spectral_selection", paste0(sample_id, "_spectral_selection.png"))
    if (identical(gate_type, "af_cosine") && file.exists(spectral_selection_path)) {
        gate_path <- spectral_selection_path
        gate_title <- "SCC/AF Spectral Selection"
    } else {
        gate_dir <- if (isTRUE(use_scatter_gating)) "intensity_scatter" else "histogram"
        gate_suffix <- if (isTRUE(use_scatter_gating)) "_intensity_scatter.png" else "_histogram.png"
        gate_path <- file.path(report_plot_dir, gate_dir, paste0(sample_id, gate_suffix))
        gate_title <- if (isTRUE(use_scatter_gating)) "Scatter/Spectral Event Gate" else "Peak-Channel Histogram Gate"
    }
    .draw_report_image_panel(
        gate_path,
        gate_title,
        x = grid::unit(0.74, "npc"),
        y = grid::unit(0.62, "npc"),
        width = grid::unit(0.42, "npc"),
        height = grid::unit(0.44, "npc")
    )
    .draw_report_image_panel(
        file.path(report_plot_dir, "spectrum", paste0(sample_id, "_spectrum.png")),
        "Per-Event Spectrum Distribution",
        x = grid::unit(0.5, "npc"),
        y = grid::unit(0.23, "npc"),
        width = grid::unit(0.9, "npc"),
        height = grid::unit(0.32, "npc")
    )

    invisible(NULL)
}

.format_af_bank_summary_lines <- function(af_bank_info) {
    if (is.null(af_bank_info) || is.null(af_bank_info$sources) || nrow(af_bank_info$sources) == 0) {
        return("No AF bank source metadata available.")
    }
    sources <- as.data.frame(af_bank_info$sources, stringsAsFactors = FALSE)
    headers <- c("File", "Source", "Events", "Gated", "Scatter%", "FSC", "SSC")
    widths <- c(30, 18, 8, 8, 9, 10, 10)

    trim_to <- function(x, w) {
        x <- as.character(x)
        x[is.na(x)] <- ""
        vapply(x, function(s) {
            if (nchar(s, type = "width") > w) {
                paste0(substr(s, 1, max(1, w - 3)), "...")
            } else {
                s
            }
        }, character(1))
    }
    fmt_num <- function(x, digits = 1) {
        out <- ifelse(is.na(x), "", format(round(as.numeric(x), digits), nsmall = digits, trim = TRUE))
        as.character(out)
    }
    make_line <- function(parts) {
        paste(vapply(seq_along(parts), function(i) sprintf(paste0("%-", widths[i], "s"), parts[[i]]), character(1)), collapse = " ")
    }

    vals <- list(
        trim_to(sources$file, widths[1]),
        trim_to(sources$source_type, widths[2]),
        trim_to(sources$n_total, widths[3]),
        trim_to(sources$n_scatter_gated, widths[4]),
        trim_to(fmt_num(sources$scatter_gate_pct, 1), widths[5]),
        trim_to(sources$fsc_channel, widths[6]),
        trim_to(sources$ssc_channel, widths[7])
    )

    c(
        make_line(as.list(headers)),
        paste(vapply(widths, function(w) paste(rep("-", w), collapse = ""), character(1)), collapse = " "),
        vapply(seq_len(nrow(sources)), function(i) make_line(lapply(vals, `[[`, i)), character(1))
    )
}

.draw_af_bank_qc_pages <- function(M, af_bank_info, pd = NULL) {
    if (is.null(af_bank_info)) {
        return(invisible(NULL))
    }
    source_count <- af_bank_info$source_count
    if (is.null(source_count)) source_count <- 1
    derived_bands <- af_bank_info$derived_bands
    if (is.null(derived_bands)) derived_bands <- 1

    if (source_count <= 1 && derived_bands <= 1) {
        return(invisible(NULL))
    }

    af_rows <- grepl("^AF($|_)", rownames(M), ignore.case = TRUE)
    if (any(af_rows)) {
        af_plot <- plot_spectra(M[af_rows, , drop = FALSE], pd = pd, output_file = NULL) +
            ggplot2::labs(title = "Autofluorescence Band Spectra Overlay") +
            ggplot2::theme(legend.position = "none")
        .draw_report_ggplot_page(af_plot, height_ratio = 0.72)
    }

    invisible(NULL)
}

.normalize_scc_report_summary <- function(qc_summary) {
    if (is.null(qc_summary)) {
        return(data.frame())
    }
    qc_summary <- as.data.frame(qc_summary, stringsAsFactors = FALSE)
    if (nrow(qc_summary) > 0 && all(c("fluorophore", "sample") %in% colnames(qc_summary))) {
        qc_summary <- qc_summary[order(qc_summary$fluorophore, qc_summary$sample), , drop = FALSE]
    }
    qc_summary
}

.resolve_scc_report_marker_mapping <- function(qc_summary) {
    sample_to_marker <- NULL
    marker_display <- NULL
    if (nrow(qc_summary) == 0) {
        return(list(sample_to_marker = sample_to_marker, marker_display = marker_display))
    }

    sample_keys <- as.character(qc_summary$sample)
    primary_vals <- trimws(as.character(qc_summary$fluorophore))
    secondary_vals <- if ("marker" %in% colnames(qc_summary)) {
        trimws(as.character(qc_summary$marker))
    } else {
        rep("", length(primary_vals))
    }
    secondary_vals[is.na(secondary_vals)] <- ""

    keep <- !is.na(sample_keys) & sample_keys != "" & !is.na(primary_vals) & primary_vals != ""
    if (any(keep)) {
        sample_to_marker <- stats::setNames(primary_vals[keep], sample_keys[keep])
        sample_to_marker <- sample_to_marker[!duplicated(names(sample_to_marker))]

        display_vals <- primary_vals
        show_secondary <- secondary_vals != "" &
            toupper(secondary_vals) != "AUTOFLUORESCENCE" &
            tolower(secondary_vals) != tolower(primary_vals)
        display_vals[show_secondary] <- paste0(primary_vals[show_secondary], " / ", secondary_vals[show_secondary])
        marker_display <- stats::setNames(display_vals[keep], primary_vals[keep])
        marker_display <- marker_display[names(marker_display) != ""]
        marker_display <- marker_display[!duplicated(names(marker_display))]
    }

    list(sample_to_marker = sample_to_marker, marker_display = marker_display)
}

.infer_scc_report_marker_mapping <- function(sample_ids, markers) {
    sample_ids <- as.character(sample_ids)
    markers <- as.character(markers)
    markers <- markers[!is.na(markers) & markers != ""]
    if (length(sample_ids) == 0 || length(markers) == 0) {
        return(NULL)
    }

    marker_key <- stats::setNames(markers, tolower(markers))
    sample_base <- tools::file_path_sans_ext(basename(sample_ids))
    sample_base <- trimws(sub("\\s*\\([^)]*\\)\\s*$", "", sample_base))
    hits <- marker_key[tolower(sample_base)]
    keep <- !is.na(hits) & hits != ""
    if (!any(keep)) {
        return(NULL)
    }
    stats::setNames(as.character(hits[keep]), sample_ids[keep])
}

.as_scc_report_unmixed_data <- function(x) {
    if (is.list(x) && "data" %in% names(x) && is.data.frame(x$data)) {
        return(as.data.frame(x$data, stringsAsFactors = FALSE, check.names = FALSE))
    }
    if (is.data.frame(x)) {
        return(as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE))
    }
    NULL
}

.combine_scc_report_unmixed_data <- function(data_list, names_to_keep) {
    names_to_keep <- intersect(names_to_keep, names(data_list))
    if (length(names_to_keep) == 0) {
        return(NULL)
    }
    parts <- data_list[names_to_keep]
    parts <- parts[!vapply(parts, is.null, logical(1))]
    if (length(parts) == 0) {
        return(NULL)
    }
    data.table::rbindlist(parts, fill = TRUE)
}

.scc_report_sample_key <- function(x) {
    tools::file_path_sans_ext(basename(trimws(as.character(x))))
}

.scc_report_find_sample_name <- function(sample_id, available_names) {
    if (length(available_names) == 0 || is.na(sample_id) || !nzchar(sample_id)) {
        return(NA_character_)
    }
    if (sample_id %in% available_names) {
        return(sample_id)
    }
    sample_key <- tolower(.scc_report_sample_key(sample_id))
    available_keys <- tolower(.scc_report_sample_key(available_names))
    hit <- which(available_keys == sample_key)
    if (length(hit) > 0) {
        return(available_names[hit[1]])
    }
    NA_character_
}

.scc_report_zero_negative_data <- function(n, markers) {
    n <- max(1L, as.integer(n))
    out <- as.data.frame(matrix(0, nrow = n, ncol = length(markers)))
    colnames(out) <- markers
    out
}

.scc_report_negative_reference <- function(sample_type,
                                           target,
                                           source_df,
                                           data_by_sample,
                                           cell_negative_names,
                                           bead_negative_names,
                                           markers) {
    if (identical(sample_type, "cells") && length(cell_negative_names) > 0) {
        neg <- .combine_scc_report_unmixed_data(data_by_sample, cell_negative_names)
        if (!is.null(neg)) return(list(data = neg, label = "unstained_cells"))
    }
    if (identical(sample_type, "beads") && length(bead_negative_names) > 0) {
        neg <- .combine_scc_report_unmixed_data(data_by_sample, bead_negative_names)
        if (!is.null(neg)) return(list(data = neg, label = "unstained_beads"))
    }
    if (identical(sample_type, "beads") && target %in% colnames(source_df)) {
        target_vals <- source_df[[target]]
        cutoff <- stats::quantile(target_vals, 0.20, na.rm = TRUE, names = FALSE)
        keep <- is.finite(target_vals) & target_vals <= cutoff
        if (sum(keep) >= 20L) {
            return(list(data = source_df[keep, , drop = FALSE], label = "low_target_beads"))
        }
    }
    list(data = .scc_report_zero_negative_data(nrow(source_df), markers), label = "zero_fallback")
}

.scc_report_safe_mad <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    stats::mad(x, na.rm = TRUE)
}

.scc_report_pair_slope <- function(target_vals, off_vals) {
    ok <- is.finite(target_vals) & is.finite(off_vals)
    if (sum(ok) < 3L) return(NA_real_)
    v <- stats::var(target_vals[ok], na.rm = TRUE)
    if (!is.finite(v) || v <= 1e-9) return(NA_real_)
    stats::cov(target_vals[ok], off_vals[ok], use = "complete.obs") / v
}

.scc_report_max_finite <- function(x, default = NA_real_) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(default)
    max(x, na.rm = TRUE)
}

.scc_report_metric_reason <- function(control_pairs, neg_ref_label) {
    if (identical(neg_ref_label, "zero_fallback")) {
        return(list(severity = "FAIL", reason = "No matched negative reference"))
    }
    nps <- .scc_report_max_finite(control_pairs$nps, default = 0)
    fpr <- .scc_report_max_finite(control_pairs$fpr, default = 0)
    bias_z <- .scc_report_max_finite(control_pairs$bias_z, default = 0)
    slope <- .scc_report_max_finite(abs(control_pairs$slope), default = 0)
    values <- c(NPS = nps / 3, FPR = fpr / 0.01, Bias = bias_z / 3, Slope = slope / 0.05)
    metric <- names(values)[which.max(values)]
    severity <- if (nps > 10 || fpr > 0.05 || bias_z > 10 || slope > 0.10) {
        "FAIL"
    } else if (nps > 3 || fpr > 0.01 || bias_z > 3 || slope > 0.05) {
        "WARN"
    } else {
        "OK"
    }
    reason <- switch(
        metric,
        NPS = sprintf("Worst NPS %.2f", nps),
        FPR = sprintf("Worst FPR %.2f%%", 100 * fpr),
        Bias = sprintf("Worst |bias|/MAD %.2f", bias_z),
        Slope = sprintf("Worst |slope| %.3f", slope),
        "Within limits"
    )
    if (identical(severity, "OK")) {
        reason <- "Within limits"
    }
    list(severity = severity, reason = reason)
}

.scc_report_floor <- function(x, floor = 1) {
    if (!is.finite(x) || x < floor) {
        return(floor)
    }
    x
}

.compute_scc_post_unmix_qc <- function(unmixed_list,
                                       qc_summary,
                                       markers) {
    if (is.null(unmixed_list) || length(unmixed_list) == 0 ||
        is.null(qc_summary) || nrow(qc_summary) == 0 || length(markers) < 2) {
        return(NULL)
    }

    data_by_sample <- lapply(unmixed_list, .as_scc_report_unmixed_data)
    data_by_sample <- data_by_sample[!vapply(data_by_sample, is.null, logical(1))]
    if (length(data_by_sample) == 0) {
        return(NULL)
    }

    sample_names <- names(data_by_sample)
    sample_files <- paste0(sample_names, ".fcs")
    bead_negative <- vapply(sample_files, .reference_is_bead_negative_file, logical(1))
    cell_negative <- !bead_negative & vapply(
        sample_files,
        function(x) .is_af_control_row(filename = x),
        logical(1)
    )
    cell_negative_names <- sample_names[cell_negative]
    bead_negative_names <- sample_names[bead_negative]

    controls <- qc_summary
    controls$sample_name <- vapply(
        as.character(controls$sample),
        .scc_report_find_sample_name,
        character(1),
        available_names = sample_names
    )
    controls <- controls[
        !is.na(controls$sample_name) &
            controls$fluorophore %in% markers &
            !grepl("^AF($|_)", controls$fluorophore, ignore.case = TRUE),
        ,
        drop = FALSE
    ]
    if (nrow(controls) == 0) {
        return(NULL)
    }

    pair_rows <- list()
    overview_rows <- list()
    for (i in seq_len(nrow(controls))) {
        target <- as.character(controls$fluorophore[i])
        source_name <- as.character(controls$sample_name[i])
        source_df <- data_by_sample[[source_name]]
        if (is.null(source_df) || !(target %in% colnames(source_df))) next
        source_type <- tolower(as.character(controls$type[i]))
        source_type <- if (source_type %in% c("cells", "beads")) source_type else "unknown"
        off_targets <- setdiff(intersect(markers, colnames(source_df)), target)
        if (length(off_targets) == 0) next

        neg_ref <- .scc_report_negative_reference(
            sample_type = source_type,
            target = target,
            source_df = source_df,
            data_by_sample = data_by_sample,
            cell_negative_names = cell_negative_names,
            bead_negative_names = bead_negative_names,
            markers = markers
        )
        neg_df <- neg_ref$data

        target_vals <- source_df[[target]]
        control_pairs <- lapply(off_targets, function(marker) {
            src <- source_df[[marker]]
            neg <- if (marker %in% colnames(neg_df)) neg_df[[marker]] else rep(0, max(1L, nrow(neg_df)))
            src_mad <- .scc_report_safe_mad(src)
            neg_mad <- .scc_report_safe_mad(neg)
            src_median <- stats::median(src, na.rm = TRUE)
            neg_median <- stats::median(neg, na.rm = TRUE)
            threshold <- stats::quantile(neg, 0.995, na.rm = TRUE, names = FALSE)
            if (!is.finite(threshold)) threshold <- 0
            bias <- src_median - neg_median
            neg_scale <- .scc_report_floor(neg_mad, 1)
            data.frame(
                control = source_name,
                type = source_type,
                target = target,
                marker = marker,
                n_events = nrow(source_df),
                negative_reference = neg_ref$label,
                nps = src_mad / neg_scale,
                source_mad = src_mad,
                negative_mad = neg_mad,
                bias = bias,
                bias_norm = bias / neg_scale,
                bias_z = abs(bias) / neg_scale,
                fpr_threshold = threshold,
                fpr = mean(src > threshold, na.rm = TRUE),
                slope = .scc_report_pair_slope(target_vals, src),
                stringsAsFactors = FALSE
            )
        })
        control_pairs <- do.call(rbind, control_pairs)
        pair_rows[[length(pair_rows) + 1L]] <- control_pairs

        worst_nps <- .scc_report_max_finite(control_pairs$nps)
        worst_fpr <- .scc_report_max_finite(control_pairs$fpr)
        worst_bias <- .scc_report_max_finite(abs(control_pairs$bias))
        worst_slope <- .scc_report_max_finite(abs(control_pairs$slope), default = 0)
        worst_bias_z <- .scc_report_max_finite(control_pairs$bias_z)
        flag_info <- .scc_report_metric_reason(control_pairs, neg_ref$label)
        overview_rows[[length(overview_rows) + 1L]] <- data.frame(
            control = source_name,
            type = source_type,
            target = target,
            n_events = nrow(source_df),
            worst_nps = worst_nps,
            worst_fpr = worst_fpr,
            worst_abs_bias = worst_bias,
            worst_abs_slope = worst_slope,
            worst_bias_mad = worst_bias_z,
            negative_reference = neg_ref$label,
            flag = flag_info$severity,
            reason = flag_info$reason,
            stringsAsFactors = FALSE
        )
    }

    if (length(pair_rows) == 0) {
        return(NULL)
    }
    list(
        overview = do.call(rbind, overview_rows),
        pairs = do.call(rbind, pair_rows)
    )
}

.draw_scc_post_unmix_overview_page <- function(metrics) {
    overview <- metrics$overview
    overview <- overview[order(match(overview$flag, c("FAIL", "WARN", "OK")), -overview$worst_fpr, -overview$worst_nps), , drop = FALSE]
    overview$control_label <- vapply(as.character(overview$control), function(x) {
        if (nchar(x, type = "width") > 28) paste0(substr(x, 1, 25), "...") else x
    }, character(1))
    table_df <- data.frame(
        row = rep(seq_len(nrow(overview)), 9),
        column = rep(seq_len(9), each = nrow(overview)),
        label = c(
            overview$control_label,
            overview$type,
            overview$target,
            as.character(overview$n_events),
            sprintf("%.2f", overview$worst_nps),
            sprintf("%.2f%%", 100 * overview$worst_fpr),
            sprintf("%.2f", overview$worst_bias_mad),
            sprintf("%.3f", overview$worst_abs_slope),
            overview$reason
        ),
        stringsAsFactors = FALSE
    )
    headers <- data.frame(
        row = 0,
        column = seq_len(10),
        label = c("Control", "Type", "Target", "Events", "NPS", "FPR", "abs bias/MAD", "abs slope", "Reason", "Status"),
        stringsAsFactors = FALSE
    )
    severity_df <- data.frame(
        row = seq_len(nrow(overview)),
        column = 10,
        label = overview$flag,
        stringsAsFactors = FALSE
    )
    severity_cols <- c(OK = "#2E7D32", WARN = "#F9A825", FAIL = "#C62828")
    table_df$fill <- ifelse(table_df$row %% 2 == 0, "#F7F9FC", "white")
    headers$fill <- "#263238"
    p <- ggplot2::ggplot() +
        ggplot2::geom_tile(
            data = table_df,
            ggplot2::aes(column, -row, fill = fill),
            color = "#E0E0E0",
            linewidth = 0.25
        ) +
        ggplot2::geom_tile(
            data = headers,
            ggplot2::aes(column, -row),
            fill = "#263238",
            color = "white",
            linewidth = 0.25
        ) +
        ggplot2::geom_tile(
            data = severity_df,
            ggplot2::aes(column, -row, fill = label),
            color = "white",
            linewidth = 0.25,
            width = 0.82,
            height = 0.62
        ) +
        ggplot2::geom_text(data = headers, ggplot2::aes(column, -row, label = label), color = "white", fontface = "bold", size = 2.7) +
        ggplot2::geom_text(data = table_df, ggplot2::aes(column, -row, label = label), color = "#1F2933", size = 2.45) +
        ggplot2::geom_text(data = severity_df, ggplot2::aes(column, -row, label = label), color = "white", fontface = "bold", size = 2.5) +
        ggplot2::scale_x_continuous(
            breaks = seq_len(10),
            labels = c("", "", "", "", "", "", "", "", "", ""),
            expand = ggplot2::expansion(mult = c(0.01, 0.01))
        ) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.03))) +
        ggplot2::scale_fill_manual(values = c(severity_cols, "white" = "white", "#F7F9FC" = "#F7F9FC"), guide = "none") +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::labs(
            title = "Post-unmixing control QC",
            subtitle = .scc_report_wrap_subtitle("Marker-space control diagnostics from already-unmixed SCC controls. Rows are sorted by severity and worst off-target behavior."),
            x = NULL,
            y = NULL
        ) +
        ggplot2::theme_void(base_size = 12) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 16),
            plot.subtitle = ggplot2::element_text(size = 9.5, lineheight = 1.05),
            plot.margin = ggplot2::margin(10, 10, 10, 10)
        )
    .draw_report_ggplot_page(p, height_ratio = 0.88)
    invisible(NULL)
}

.scc_report_clip <- function(x, cap) {
    if (is.null(cap) || !is.finite(cap) || cap <= 0) {
        return(x)
    }
    pmin(pmax(x, -cap), cap)
}

.scc_report_label_values <- function(x, digits, threshold = 0) {
    out <- ifelse(is.finite(x) & abs(x) >= threshold, sprintf(paste0("%.", digits, "f"), x), "")
    out[is.na(out)] <- ""
    out[grepl("^-?0(?:\\.0+)?$", out)] <- ""
    out
}

.scc_report_scale_cap <- function(x, requested_cap = NULL, diverging = FALSE) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(1)
    cap <- if (!is.null(requested_cap) && is.finite(requested_cap) && requested_cap > 0) requested_cap else max(abs(x), na.rm = TRUE)
    if (!is.finite(cap) || cap <= 0) cap <- 1
    if (isTRUE(diverging)) cap else max(cap, 1e-9)
}

.write_scc_report_csv <- function(x, path, row_id = NULL) {
    if (length(path) == 0 || is.null(path) || is.na(path[1]) || !nzchar(trimws(as.character(path[1])))) {
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

.write_scc_report_matrix_metric <- function(mat, path, row_id = "marker") {
    if (length(path) == 0 || is.null(path) || is.null(mat) || length(mat) == 0) {
        return(invisible(NULL))
    }
    .write_scc_report_csv(as.data.frame(mat, check.names = FALSE), path, row_id = row_id)
}

.write_scc_control_qc_summary <- function(qc_summary, qc_metrics_dir) {
    if (is.null(qc_metrics_dir) || length(qc_metrics_dir) == 0 ||
        is.na(qc_metrics_dir[1]) || !nzchar(trimws(as.character(qc_metrics_dir)[1])) ||
        is.null(qc_summary) || nrow(qc_summary) == 0) {
        return(invisible(NULL))
    }
    .write_scc_report_csv(qc_summary, file.path(qc_metrics_dir, "control_qc_summary.csv"))
}

.write_scc_post_unmix_metrics <- function(metrics, qc_metrics_dir) {
    if (is.null(qc_metrics_dir) || length(qc_metrics_dir) == 0 ||
        is.na(qc_metrics_dir[1]) || !nzchar(trimws(as.character(qc_metrics_dir)[1])) ||
        is.null(metrics)) {
        return(invisible(NULL))
    }
    if (!is.null(metrics$pairs) && nrow(metrics$pairs) > 0) {
        .write_scc_report_csv(
            metrics$pairs,
            file.path(qc_metrics_dir, "post_unmix_control_qc_pairs.csv")
        )
    }
    if (!is.null(metrics$overview) && nrow(metrics$overview) > 0) {
        .write_scc_report_csv(
            metrics$overview,
            file.path(qc_metrics_dir, "post_unmix_control_qc_overview.csv")
        )
    }
    invisible(NULL)
}

.scc_report_diag_marker_df <- function(marker_levels) {
    data.frame(
        target = marker_levels,
        marker = marker_levels,
        source = factor(marker_levels, levels = rev(marker_levels)),
        marker_f = factor(marker_levels, levels = marker_levels),
        stringsAsFactors = FALSE
    )
}

.scc_report_detector_order <- function(detectors) {
    detector_key <- toupper(gsub("\\s+", "", gsub("-A$", "", as.character(detectors))))
    laser_group <- vapply(detector_key, function(k) {
        if (grepl("^UV", k)) return(1L)
        if (grepl("^V", k)) return(2L)
        if (grepl("^B", k)) return(3L)
        if (grepl("^YG", k) || grepl("^Y", k) || grepl("^G", k)) return(4L)
        if (grepl("^R", k)) return(5L)
        99L
    }, integer(1))
    det_num <- suppressWarnings(as.integer(sub("^[A-Z]+([0-9]+).*$", "\\1", detector_key)))
    det_num[is.na(det_num)] <- 999L
    order(laser_group, det_num, detector_key, na.last = TRUE)
}

.scc_report_marker_levels <- function(markers, qc_summary = NULL) {
    markers <- unique(as.character(markers))
    markers <- markers[nzchar(markers) & !is.na(markers)]
    if (length(markers) == 0) {
        return(character())
    }
    if (!is.null(qc_summary) && all(c("fluorophore", "peak_channel") %in% colnames(qc_summary))) {
        q <- as.data.frame(qc_summary, stringsAsFactors = FALSE)
        q <- q[
            q$fluorophore %in% markers &
                !is.na(q$peak_channel) &
                nzchar(trimws(as.character(q$peak_channel))),
            ,
            drop = FALSE
        ]
        if (nrow(q) > 0) {
            detector_order <- .scc_report_detector_order(q$peak_channel)
            q$.detector_rank <- seq_along(detector_order)
            q$.detector_rank[detector_order] <- seq_along(detector_order)
            q <- q[order(q$.detector_rank, tolower(q$fluorophore)), , drop = FALSE]
            levels <- unique(as.character(q$fluorophore))
            return(c(levels, sort(setdiff(markers, levels))))
        }
    }
    sort(markers)
}

.scc_report_wrap_subtitle <- function(x, width = 72) {
    paste(strwrap(x, width = width), collapse = "\n")
}

.plot_scc_pair_heatmap <- function(pair_df,
                                   value_col,
                                   title,
                                   fill_label,
                                   subtitle,
                                   marker_levels,
                                   diverging = FALSE,
                                   percent = FALSE,
                                   label_digits = 2,
                                   color_cap = NULL,
                                   label_threshold = 0) {
    marker_levels <- unique(as.character(marker_levels))
    marker_levels <- marker_levels[marker_levels %in% unique(c(pair_df$target, pair_df$marker))]
    if (length(marker_levels) == 0) {
        return(NULL)
    }

    df <- pair_df[, c("target", "marker", value_col), drop = FALSE]
    names(df)[names(df) == value_col] <- "value"
    df <- df[is.finite(df$value), , drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    grid_df <- expand.grid(
        target = marker_levels,
        marker = marker_levels,
        stringsAsFactors = FALSE
    )
    df <- merge(grid_df, df, by = c("target", "marker"), all.x = TRUE, sort = FALSE)
    df$source <- factor(df$target, levels = rev(marker_levels))
    df$marker <- factor(df$marker, levels = marker_levels)
    if (isTRUE(percent)) {
        df$value <- 100 * df$value
    }
    label_df <- df[is.finite(df$value), , drop = FALSE]
    label_df$label <- .scc_report_label_values(label_df$value, digits = label_digits, threshold = label_threshold)
    label_df <- label_df[nzchar(label_df$label), , drop = FALSE]
    df$plot_value <- .scc_report_clip(df$value, color_cap)
    n_markers <- length(marker_levels)
    text_size <- max(1.7, min(3.6, 36 / max(1, n_markers)))
    diag_df <- .scc_report_diag_marker_df(marker_levels)
    p <- ggplot2::ggplot(df, ggplot2::aes(marker, source, fill = plot_value)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.25) +
        ggplot2::geom_tile(
            data = diag_df,
            ggplot2::aes(marker_f, source),
            inherit.aes = FALSE,
            fill = "#E6E8EB",
            color = "white",
            linewidth = 0.25
        ) +
        ggplot2::labs(
            title = title,
            subtitle = .scc_report_wrap_subtitle(subtitle),
            x = "Off-target unmixed marker channel",
            y = "Source SCC fluorophore",
            fill = fill_label
        ) +
        ggplot2::theme_minimal(base_size = 12.5) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
            axis.text.y = ggplot2::element_text(size = 8),
            plot.subtitle = ggplot2::element_text(size = 9.5, lineheight = 1.05),
            panel.grid = ggplot2::element_blank()
        )
    if (isTRUE(diverging)) {
        lim <- .scc_report_scale_cap(df$value, requested_cap = color_cap, diverging = TRUE)
        if (!is.finite(lim) || lim <= 0) {
            lim <- 1
        }
        label_df$text_color <- ifelse(abs(.scc_report_clip(label_df$value, lim)) > 0.62 * lim, "white", "black")
        p +
            ggplot2::geom_text(
                data = label_df,
                ggplot2::aes(marker, source, label = label, color = text_color),
                inherit.aes = FALSE,
                size = text_size,
                show.legend = FALSE
            ) +
            ggplot2::scale_color_identity() +
            ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-lim, lim), na.value = "grey92")
    } else {
        max_val <- .scc_report_scale_cap(df$value, requested_cap = color_cap, diverging = FALSE)
        if (!is.finite(max_val) || max_val <= 0) {
            max_val <- 1
        }
        label_df$text_color <- ifelse(.scc_report_clip(label_df$value, max_val) > 0.62 * max_val, "white", "black")
        p +
            ggplot2::geom_text(
                data = label_df,
                ggplot2::aes(marker, source, label = label, color = text_color),
                inherit.aes = FALSE,
                size = text_size,
                show.legend = FALSE
            ) +
            ggplot2::scale_color_identity() +
            ggplot2::scale_fill_gradient(low = "white", high = "#B2182B", na.value = "grey92")
    }
}

.plot_scc_nps_marker_summary <- function(pair_df, marker_levels) {
    df <- pair_df |>
        dplyr::group_by(marker) |>
        dplyr::summarise(NPS = max(nps, na.rm = TRUE), .groups = "drop")
    if (nrow(df) == 0) return(NULL)
    marker_levels <- unique(as.character(marker_levels))
    marker_levels <- marker_levels[marker_levels %in% df$marker]
    df <- df[match(marker_levels, df$marker), , drop = FALSE]
    df <- df[!is.na(df$marker), , drop = FALSE]
    df$marker <- factor(df$marker, levels = marker_levels)
    df$label <- sprintf("%.2f", df$NPS)
    y_max <- max(df$NPS, na.rm = TRUE)
    if (!is.finite(y_max) || y_max <= 0) y_max <- 1
    ggplot2::ggplot(df, ggplot2::aes(marker, NPS)) +
        ggplot2::geom_col(fill = "#C44E52", width = 0.7) +
        ggplot2::geom_text(ggplot2::aes(label = label), vjust = -0.25, size = 3) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.14))) +
        ggplot2::coord_cartesian(ylim = c(0, y_max * 1.14), clip = "off") +
        ggplot2::labs(
            title = "Off-target negative spread by marker",
            subtitle = .scc_report_wrap_subtitle("X-axis entries are off-target unmixed marker channels. Each bar is the worst NPS seen in that channel across all source SCC controls."),
            x = "Off-target unmixed marker channel",
            y = "Worst NPS (unitless MAD ratio)"
        ) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
            plot.subtitle = ggplot2::element_text(size = 9.5, lineheight = 1.05)
        )
}

.draw_scc_post_unmix_qc_pages <- function(unmixed_list, qc_summary, markers, qc_metrics_dir = NULL) {
    metrics <- .compute_scc_post_unmix_qc(
        unmixed_list = unmixed_list,
        qc_summary = qc_summary,
        markers = markers
    )
    if (is.null(metrics)) {
        return(invisible(NULL))
    }
    .write_scc_post_unmix_metrics(metrics, qc_metrics_dir)
    marker_levels <- .scc_report_marker_levels(markers, qc_summary = qc_summary)
    plot_pages <- list(
        .plot_scc_nps_marker_summary(metrics$pairs, marker_levels = marker_levels),
        .plot_scc_pair_heatmap(
            metrics$pairs,
            "nps",
            "Pairwise off-target negative spread",
            "NPS",
            subtitle = "Rows are source SCCs; columns are off-target unmixed marker channels. Color is capped at NPS = 10 so one extreme pair does not flatten the rest.",
            marker_levels = marker_levels,
            label_digits = 2,
            color_cap = 10
        ),
        .plot_scc_pair_heatmap(
            metrics$pairs,
            "fpr",
            "Pairwise false-positive rate",
            "FPR (%)",
            subtitle = "Percent of events in each source SCC whose off-target value exceeds the 99.5th percentile of the matched negative reference.",
            marker_levels = marker_levels,
            percent = TRUE,
            label_digits = 2
        ),
        .plot_scc_pair_heatmap(
            metrics$pairs,
            "bias_norm",
            "Pairwise normalized off-target bias",
            "Bias / MADneg",
            subtitle = "Signed median off-target offset divided by negative-reference MAD. Color is symmetrically clipped at +/-10.",
            marker_levels = marker_levels,
            diverging = TRUE,
            label_digits = 2,
            color_cap = 10,
            label_threshold = 0.2
        ),
        .plot_scc_pair_heatmap(
            metrics$pairs,
            "slope",
            "Pairwise target-driven off-target slope",
            "Slope",
            subtitle = "Covariance slope of off-target signal against the true target signal. Larger absolute values indicate target-dependent spillover artifacts.",
            marker_levels = marker_levels,
            diverging = TRUE,
            label_digits = 3,
            color_cap = 0.05,
            label_threshold = 0.005
        )
    )
    for (p in plot_pages) {
        if (!is.null(p)) {
            .draw_report_ggplot_page(p, height_ratio = 0.78)
        }
    }
    invisible(metrics)
}

.read_scc_report_unmixed_fcs <- function(unmixed_dir) {
    if (is.null(unmixed_dir) || !dir.exists(unmixed_dir)) {
        return(NULL)
    }
    fcs_files <- list.files(unmixed_dir, pattern = "_unmixed\\.fcs$", full.names = TRUE, ignore.case = TRUE)
    if (length(fcs_files) == 0) {
        return(NULL)
    }

    out <- lapply(fcs_files, function(path) {
        ff <- flowCore::read.FCS(path, transformation = FALSE, truncate_max_range = FALSE)
        data <- as.data.frame(flowCore::exprs(ff), stringsAsFactors = FALSE, check.names = FALSE)
        list(data = data, residuals = NULL)
    })
    names(out) <- sub("_unmixed\\.fcs$", "", basename(fcs_files), ignore.case = TRUE)
    class(out) <- c("spectreasy_unmixed_results", "list")
    out
}

.write_scc_qc_report <- function(M_built,
                                 M_report = NULL,
                                 qc_summary = NULL,
                                 report_plot_dir = NULL,
                                 unmixed_list = NULL,
                                 scc_dir = "scc",
                                 output_file = "spectreasy_outputs/unmix_controls/qc_controls_report.pdf",
                                 cytometer = "auto",
                                 method = "AutoSpectral",
                                 use_scatter_gating = TRUE,
                                 unmix_scatter_max_points = 1000,
                                 unmix_scatter_axis_limit = NULL,
                                 seed = NULL,
                                 qc_metrics_dir = NULL,
                                 retained_qc_plot_dir = NULL) {
    if (is.null(output_file) || !nzchar(trimws(as.character(output_file)[1]))) {
        stop("Please supply output_file to save the SCC PDF report.", call. = FALSE)
    }
    method <- .normalize_unmix_method(method)
    if (is.null(M_built) || nrow(M_built) == 0) {
        stop("No valid spectra found while generating the SCC report.")
    }
    if (is.null(M_report)) {
        M_report <- M_built
    } else {
        M_report <- .as_reference_matrix(M_report, "M_report")
    }
    qc_summary <- .normalize_scc_report_summary(qc_summary)
    if (is.null(report_plot_dir) || !dir.exists(report_plot_dir)) {
        report_plot_dir <- attr(M_built, "qc_plot_dir")
    }

    fcs_files <- if (dir.exists(scc_dir)) {
        list.files(scc_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
    } else {
        character()
    }
    pd <- attr(M_built, "detector_pd")
    if (length(fcs_files) > 0) {
        ff_meta <- flowCore::read.FCS(fcs_files[1], transformation = FALSE, truncate_max_range = FALSE)
        pd <- flowCore::pData(flowCore::parameters(ff_meta))
    }
    if (!is.data.frame(pd)) {
        pd <- NULL
    }

    if (!requireNamespace("png", quietly = TRUE)) {
        stop("Package 'png' is required to embed SCC QC plots in the report.")
    }

    out_dir <- dirname(output_file)
    if (!is.na(out_dir) && nzchar(out_dir) && out_dir != ".") {
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    }
    grDevices::pdf(output_file, width = 11, height = 8.5)
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
    .write_scc_control_qc_summary(qc_summary, qc_metrics_dir)

    scc_dir_line <- paste(strwrap(paste0("SCC directory: ", normalizePath(scc_dir, mustWork = FALSE)), width = 85), collapse = "\n")

    grid::grid.newpage()
    grid::grid.text("spectreasy: Single-Color Control Review", x = 0.5, y = 0.7, gp = grid::gpar(fontsize = 20, fontface = "bold"))
    grid::grid.text(
        paste0(
            "Generated on: ", Sys.time(), "\n",
            scc_dir_line, "\n",
            "Controls processed: ", nrow(qc_summary), "\n",
            "Unmixing method for scatter matrix: ", method, "\n",
            "Workflow intent: review SCC quality tied to the reported unmixing matrix."
        ),
        x = 0.5,
        y = 0.48,
        just = "center",
        gp = grid::gpar(fontsize = 11, lineheight = 1.3)
    )

    .draw_af_bank_qc_pages(M_built, attr(M_built, "af_bank_info"), pd = pd)

    keep_non_af <- !grepl("^AF($|_)", rownames(M_report), ignore.case = TRUE)
    M_no_af <- M_report[keep_non_af, , drop = FALSE]
    M_af_report <- M_report[!keep_non_af, , drop = FALSE]

    if (nrow(M_no_af) > 0) {
        .write_scc_report_matrix_metric(
            M_no_af,
            file.path(qc_metrics_dir, "reference_spectra.csv"),
            row_id = "fluorophore"
        )
        .draw_report_ggplot_page(plot_spectra(M_no_af, pd = pd, output_file = NULL), height_ratio = 0.6)
    }
    if (nrow(M_af_report) > 0) {
        .write_scc_report_matrix_metric(
            M_af_report,
            file.path(qc_metrics_dir, "af_bank_spectra.csv"),
            row_id = "af_band"
        )
    }

    if (nrow(M_no_af) > 1) {
        sim_mat <- calculate_similarity_matrix(M_no_af)
        .write_scc_report_matrix_metric(
            sim_mat,
            file.path(qc_metrics_dir, "fluorophore_spectral_similarity.csv"),
            row_id = "fluorophore"
        )
        .draw_report_ggplot_page(plot_similarity_matrix(sim_mat, output_file = NULL))
        ssm_method <- if (method %in% c("NNLS", "RWLS", "AutoSpectral")) {
            if (identical(method, "RWLS")) "WLS" else "OLS"
        } else {
            method
        }
        ssm_mat <- calculate_ssm(M_no_af, method = ssm_method)
        .write_scc_report_matrix_metric(
            ssm_mat,
            file.path(qc_metrics_dir, "spectral_spread_matrix.csv"),
            row_id = "spilling_marker"
        )
        .draw_report_ggplot_page(plot_ssm(ssm_mat, output_file = NULL))

        if (is.null(unmixed_list)) {
            unmixed_list <- unmix_samples(
                sample_dir = scc_dir,
                M = M_report,
                method = method,
                cytometer = cytometer,
                write_fcs = FALSE,
                save_report = FALSE,
                verbose = FALSE
            )
        }
        marker_mapping <- .resolve_scc_report_marker_mapping(qc_summary)
        if (is.null(marker_mapping$sample_to_marker) && !is.null(unmixed_list)) {
            marker_mapping$sample_to_marker <- .infer_scc_report_marker_mapping(
                sample_ids = names(unmixed_list),
                markers = rownames(M_no_af)
            )
        }
        p_scatter <- tryCatch(
            plot_unmixing_scatter_matrix(
                unmixed_list = unmixed_list,
                sample_to_marker = marker_mapping$sample_to_marker,
                markers = rownames(M_no_af),
                marker_display = NULL,
                output_file = NULL,
                max_points_per_sample = unmix_scatter_max_points,
                axis_limit = unmix_scatter_axis_limit,
                seed = seed
            ),
            error = function(e) NULL
        )
        if (!is.null(p_scatter)) {
            .draw_report_ggplot_page(p_scatter, square = TRUE)
        }

        .draw_scc_post_unmix_qc_pages(
            unmixed_list = unmixed_list,
            qc_summary = qc_summary,
            markers = rownames(M_no_af),
            qc_metrics_dir = qc_metrics_dir
        )
    }

    if (nrow(qc_summary) > 0 && !is.null(report_plot_dir) && dir.exists(report_plot_dir)) {
        for (i in seq_len(nrow(qc_summary))) {
            .draw_scc_report_sample_page(qc_summary[i, , drop = FALSE], report_plot_dir = report_plot_dir, use_scatter_gating = use_scatter_gating)
        }
    }

    message("SCC report saved to: ", output_file)
    invisible(list(M = M_built, qc_summary = qc_summary, qc_plot_dir = retained_qc_plot_dir, qc_metrics_dir = qc_metrics_dir, af_bank_info = attr(M_built, "af_bank_info"), method = method))
}

#' Generate SCC QC Report
#'
#' Builds a reference matrix from single-color controls, assembles a PDF report
#' for pre-unmix review, adds post-unmixing control QC from the already-unmixed
#' controls, and can optionally retain the intermediate QC PNG plots.
#'
#' @param results Optional result list returned by [unmix_controls()]. When
#'   supplied, its in-memory matrix and unmixed SCCs are reused instead of
#'   rerunning SCC unmixing.
#' @param M Optional precomputed reference matrix. If omitted, the report uses the
#'   matrix generated from the supplied SCC files.
#' @param unmixing_matrix_file Optional CSV path to a saved reference matrix.
#'   Used when `M` is not supplied. By default this points to the reference matrix
#'   produced by [unmix_controls()] (`"scc_reference_matrix.csv"`).
#' @param scc_dir Directory containing SCC FCS files.
#' @param output_file Path to save the PDF report. Defaults to `"spectreasy_outputs/unmix_controls/qc_controls_report.pdf"`.
#' @param unmix_controls_dir Optional directory from a previous `unmix_controls()`
#'   run. When supplied,
#'   defaults for `output_file`, `qc_plot_dir`, and `unmixing_matrix_file` are
#'   resolved relative to this directory.
#' @param control_file Control mapping CSV path.
#' @param cytometer Cytometer name passed to [build_reference_matrix()]. The
#'   default, `"auto"`, infers the cytometer from FCS detector names when possible.
#' @param method Unmixing method used for the control scatter matrix
#'   (`"AutoSpectral"`, `"OLS"`, `"WLS"`, `"RWLS"`, or `"NNLS"`).
#' @param qc_plot_dir Directory where FSC/SSC, intensity-gate, optional
#'   spectral-selection, and spectrum PNGs are written when
#'   `save_qc_pngs = TRUE`.
#' @param save_qc_pngs Logical; if `TRUE`, keep the intermediate QC PNG files in
#'   `qc_plot_dir`. If `FALSE` (default), PNGs are written to a temporary directory
#'   for report assembly and removed afterward.
#' @param use_scatter_gating Logical; if `TRUE` (default), use broad scatter
#'   cleanup plus the intensity-vs-FSC GMM/EM selector and show the scatter
#'   gate plot in the report. If `FALSE`, use and show the legacy
#'   one-dimensional histogram gate. Reports only show the spectral-selection
#'   plot for controls whose recorded gate type is `af_cosine`.
#' @param af_bands_per_file Deprecated compatibility argument. Multiple AF
#'   sources are pooled before SOM extraction; `af_n_bands`/`af_auto_max_bands`
#'   control the size of the one shared AF bank.
#' @param unmix_scatter_max_points Maximum events sampled per control for the
#'   SCC unmixing scatter matrix.
#' @param unmix_scatter_axis_limit Optional fixed symmetric axis limit for the
#'   SCC unmixing scatter matrix. The default `NULL` uses local per-panel
#'   ranges. Use `1e5` for `c(-1e5, 1e5)` on every panel.
#' @param seed Optional integer seed for deterministic subsampling/clustering.
#' @param ... Additional arguments forwarded to [build_reference_matrix()].
#' @return Invisibly returns a list with `M`, `qc_summary`, `qc_plot_dir`,
#'   `af_bank_info`, and `method`.
#'   `qc_plot_dir` is `NULL` unless `save_qc_pngs = TRUE`.
#' @examples
#' if (interactive()) {
#'   qc_controls(
#'     scc_dir = "scc",
#'     control_file = "fcs_mapping.csv"
#'   )
#' }
#' @export
qc_controls <- function(
    results = NULL,
    M = NULL,
    unmixing_matrix_file = file.path("spectreasy_outputs", "unmix_controls", "scc_reference_matrix.csv"),
    scc_dir = "scc",
    output_file = "spectreasy_outputs/unmix_controls/qc_controls_report.pdf",
    unmix_controls_dir = NULL,
    control_file = "fcs_mapping.csv",
    cytometer = "auto",
    method = "AutoSpectral",
    qc_plot_dir = file.path("spectreasy_outputs", "scc_report_plots"),
    save_qc_pngs = FALSE,
    use_scatter_gating = TRUE,
    af_bands_per_file = NULL,
    unmix_scatter_max_points = 1000,
    unmix_scatter_axis_limit = NULL,
    seed = NULL,
    ...
) {
    output_file_missing <- missing(output_file)
    qc_plot_dir_missing <- missing(qc_plot_dir)
    method_missing <- missing(method)
    if (!is.null(unmix_controls_dir)) {
        if (output_file_missing) {
            output_file <- file.path(unmix_controls_dir, "qc_controls_report.pdf")
        }
        if (qc_plot_dir_missing) {
            qc_plot_dir <- unmix_controls_dir
        }
        if (missing(unmixing_matrix_file)) {
            unmixing_matrix_file <- file.path(unmix_controls_dir, "scc_reference_matrix.csv")
        }
    }
    if (is.null(output_file) || !nzchar(trimws(as.character(output_file)[1]))) {
        stop("Please supply output_file to save the SCC PDF report.", call. = FALSE)
    }
    method <- .normalize_unmix_method(method)

    message("Generating SCC QC report...")
    if (!is.null(results)) {
        M_built <- .as_reference_matrix(results$M, "results$M")
        M_report <- if (!is.null(M)) .as_reference_matrix(M, "M") else M_built
        qc_summary <- if (!is.null(results$qc_summary)) results$qc_summary else attr(M_built, "qc_summary")
        report_plot_dir <- if (!is.null(results$qc_plot_dir)) results$qc_plot_dir else attr(M_built, "qc_plot_dir")
        retained_qc_plot_dir <- if (isTRUE(save_qc_pngs)) report_plot_dir else NULL
        method <- if (method_missing && !is.null(results$static_unmixing_matrix_method)) {
            results$static_unmixing_matrix_method
        } else {
            method
        }
        return(.write_scc_qc_report(
            M_built = M_built,
            M_report = M_report,
            qc_summary = qc_summary,
            report_plot_dir = report_plot_dir,
            unmixed_list = results$unmixed_list,
            scc_dir = scc_dir,
            output_file = output_file,
            cytometer = cytometer,
            method = method,
            use_scatter_gating = use_scatter_gating,
            unmix_scatter_max_points = unmix_scatter_max_points,
            unmix_scatter_axis_limit = unmix_scatter_axis_limit,
            seed = seed,
            retained_qc_plot_dir = retained_qc_plot_dir
        ))
    }

    if (!is.null(unmix_controls_dir) && file.exists(unmixing_matrix_file)) {
        M_built <- .read_unmixing_matrix_csv(unmixing_matrix_file)
        M_built <- .as_reference_matrix(M_built, "unmixing_matrix_file")
        M_report <- if (!is.null(M)) .as_reference_matrix(M, "M") else M_built
        report_plot_dir <- if (dir.exists(file.path(unmix_controls_dir, "fsc_ssc"))) unmix_controls_dir else NULL
        unmixed_list <- .read_scc_report_unmixed_fcs(file.path(unmix_controls_dir, "unmixed_fcs"))
        return(.write_scc_qc_report(
            M_built = M_built,
            M_report = M_report,
            qc_summary = attr(M_built, "qc_summary"),
            report_plot_dir = report_plot_dir,
            unmixed_list = unmixed_list,
            scc_dir = scc_dir,
            output_file = output_file,
            cytometer = cytometer,
            method = method,
            use_scatter_gating = use_scatter_gating,
            unmix_scatter_max_points = unmix_scatter_max_points,
            unmix_scatter_axis_limit = unmix_scatter_axis_limit,
            seed = seed,
            retained_qc_plot_dir = if (isTRUE(save_qc_pngs)) report_plot_dir else NULL
        ))
    }

    plot_dir_info <- .prepare_scc_report_plot_dir(qc_plot_dir = qc_plot_dir, save_qc_pngs = save_qc_pngs)
    if (!is.null(plot_dir_info$cleanup_dir)) {
        on.exit(unlink(plot_dir_info$cleanup_dir, recursive = TRUE, force = TRUE), add = TRUE)
    }

    control_resolved <- .resolve_control_file_path(control_file)
    control_input <- if (file.exists(control_resolved)) control_resolved else NULL

    M_built <- build_reference_matrix(
        input_folder = scc_dir,
        output_folder = plot_dir_info$plot_dir,
        save_qc_plots = TRUE,
        control_df = control_input,
        af_bands_per_file = af_bands_per_file,
        cytometer = cytometer,
        use_scatter_gating = use_scatter_gating,
        seed = seed,
        ...
    )
    if (is.null(M_built) || nrow(M_built) == 0) {
        stop("No valid spectra found while generating the SCC report.")
    }

    qc_summary <- attr(M_built, "qc_summary")
    if (is.null(qc_summary)) {
        qc_summary <- data.frame()
    } else {
        qc_summary <- as.data.frame(qc_summary, stringsAsFactors = FALSE)
        qc_summary <- qc_summary[order(qc_summary$fluorophore, qc_summary$sample), , drop = FALSE]
    }

    report_plot_dir <- attr(M_built, "qc_plot_dir")
    if (is.null(report_plot_dir) || !dir.exists(report_plot_dir)) {
        report_plot_dir <- plot_dir_info$plot_dir
    }
    retained_qc_plot_dir <- if (isTRUE(save_qc_pngs)) report_plot_dir else NULL

    M_report <- NULL
    if (!is.null(M)) {
        M_report <- .as_reference_matrix(M, "M")
    } else if (!is.null(unmixing_matrix_file)) {
        if (file.exists(unmixing_matrix_file)) {
            .stop_if_static_unmixing_matrix_path(unmixing_matrix_file, arg_name = "unmixing_matrix_file")
            M_report <- .read_unmixing_matrix_csv(unmixing_matrix_file)
            M_report <- .as_reference_matrix(M_report, "M")
        }
    }

    if (is.null(M_report)) {
        M_report <- M_built
    }

    .write_scc_qc_report(
        M_built = M_built,
        M_report = M_report,
        qc_summary = qc_summary,
        report_plot_dir = report_plot_dir,
        unmixed_list = NULL,
        scc_dir = scc_dir,
        output_file = output_file,
        cytometer = cytometer,
        method = method,
        use_scatter_gating = use_scatter_gating,
        unmix_scatter_max_points = unmix_scatter_max_points,
        unmix_scatter_axis_limit = unmix_scatter_axis_limit,
        seed = seed,
        retained_qc_plot_dir = retained_qc_plot_dir
    )
}
