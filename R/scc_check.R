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
    grid::grid.newpage()
    grid::grid.text("AF Bank QC", x = 0.05, y = 0.95, just = c("left", "top"), gp = grid::gpar(fontsize = 16, fontface = "bold"))
    grid::grid.text(
        paste0(
            "AF sources pooled: ", af_bank_info$source_count, "\n",
            "Pooled scatter-gated AF events: ", af_bank_info$pooled_events, "\n",
            "Bands requested: ", af_bank_info$requested_bands, "\n",
            "Bands derived: ", af_bank_info$derived_bands, "\n",
            "Mode: ", af_bank_info$mode
        ),
        x = 0.05,
        y = 0.88,
        just = c("left", "top"),
        gp = grid::gpar(fontsize = 11, lineheight = 1.25)
    )
    grid::grid.text(
        paste(.format_af_bank_summary_lines(af_bank_info), collapse = "\n"),
        x = 0.05,
        y = 0.68,
        just = c("left", "top"),
        gp = grid::gpar(fontsize = 8.5, fontfamily = "mono", lineheight = 1.15)
    )

    if (any(af_rows)) {
        af_plot <- plot_spectra(M[af_rows, , drop = FALSE], pd = pd, output_file = NULL) +
            ggplot2::labs(title = "Autofluorescence Band Spectra Overlay")
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
        flag <- if (identical(neg_ref$label, "zero_fallback")) {
            "CHECK"
        } else if (is.finite(worst_nps) && is.finite(worst_fpr) &&
            (worst_nps > 3 || worst_fpr > 0.01 || worst_bias_z > 3 || worst_slope > 0.05)) {
            "WARN"
        } else {
            "OK"
        }
        overview_rows[[length(overview_rows) + 1L]] <- data.frame(
            control = source_name,
            type = source_type,
            target = target,
            n_events = nrow(source_df),
            worst_nps = worst_nps,
            worst_fpr = worst_fpr,
            worst_abs_bias = worst_bias,
            worst_abs_slope = worst_slope,
            negative_reference = neg_ref$label,
            flag = flag,
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

.format_scc_post_unmix_overview_lines <- function(overview) {
    if (is.null(overview) || nrow(overview) == 0) {
        return("No post-unmixing control QC metrics available.")
    }
    headers <- c("Control", "Type", "Target", "Events", "Worst NPS", "Worst FPR", "Worst |bias|", "Worst |slope|", "Flag")
    widths <- c(24, 6, 12, 8, 10, 10, 13, 14, 6)
    trim_to <- function(x, w) {
        x <- as.character(x)
        x[is.na(x)] <- ""
        vapply(x, function(s) {
            if (nchar(s, type = "width") > w) paste0(substr(s, 1, max(1, w - 3)), "...") else s
        }, character(1))
    }
    fmt <- function(x, digits = 2) {
        x <- as.numeric(x)
        ifelse(is.finite(x), format(round(x, digits), nsmall = digits, trim = TRUE), "")
    }
    vals <- list(
        trim_to(overview$control, widths[1]),
        trim_to(overview$type, widths[2]),
        trim_to(overview$target, widths[3]),
        trim_to(overview$n_events, widths[4]),
        trim_to(fmt(overview$worst_nps, 2), widths[5]),
        trim_to(sprintf("%.2f%%", 100 * overview$worst_fpr), widths[6]),
        trim_to(fmt(overview$worst_abs_bias, 1), widths[7]),
        trim_to(fmt(overview$worst_abs_slope, 3), widths[8]),
        trim_to(overview$flag, widths[9])
    )
    make_line <- function(parts) {
        paste(vapply(seq_along(parts), function(i) sprintf(paste0("%-", widths[i], "s"), parts[[i]]), character(1)), collapse = " ")
    }
    c(
        make_line(as.list(headers)),
        paste(vapply(widths, function(w) paste(rep("-", w), collapse = ""), character(1)), collapse = " "),
        vapply(seq_len(nrow(overview)), function(i) make_line(lapply(vals, `[[`, i)), character(1))
    )
}

.draw_scc_post_unmix_overview_page <- function(metrics) {
    overview <- metrics$overview
    overview <- overview[order(match(overview$flag, c("CHECK", "WARN", "OK")), -overview$worst_fpr, -overview$worst_nps), , drop = FALSE]
    grid::grid.newpage()
    grid::grid.text("Post-unmixing control QC", x = 0.05, y = 0.95, just = c("left", "top"), gp = grid::gpar(fontsize = 16, fontface = "bold"))
    grid::grid.text(
        "Marker-space control diagnostics from already-unmixed SCC controls. Cell SCCs are compared to unstained cells; bead SCCs to bead-negative controls or low-target bead events.",
        x = 0.05,
        y = 0.90,
        just = c("left", "top"),
        gp = grid::gpar(fontsize = 10, lineheight = 1.15)
    )
    grid::grid.text(
        paste(.format_scc_post_unmix_overview_lines(overview), collapse = "\n"),
        x = 0.05,
        y = 0.82,
        just = c("left", "top"),
        gp = grid::gpar(fontsize = 8.2, fontfamily = "mono", lineheight = 1.14)
    )
}

.plot_scc_pair_heatmap <- function(pair_df,
                                   value_col,
                                   title,
                                   fill_label,
                                   diverging = FALSE,
                                   percent = FALSE) {
    df <- pair_df
    df$value <- df[[value_col]]
    df <- df[is.finite(df$value), , drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    df$source <- factor(df$target, levels = unique(df$target))
    df$marker <- factor(df$marker, levels = unique(df$marker))
    if (isTRUE(percent)) {
        df$value <- 100 * df$value
    }
    p <- ggplot2::ggplot(df, ggplot2::aes(marker, source, fill = value)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.25) +
        ggplot2::labs(title = title, x = "Off-target unmixed marker", y = "Source SCC fluorophore", fill = fill_label) +
        ggplot2::theme_minimal(base_size = 12.5) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
            axis.text.y = ggplot2::element_text(size = 8),
            panel.grid = ggplot2::element_blank()
        )
    if (isTRUE(diverging)) {
        lim <- max(abs(df$value), na.rm = TRUE)
        if (!is.finite(lim) || lim <= 0) {
            lim <- 1
        }
        p + ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-lim, lim))
    } else {
        p + ggplot2::scale_fill_gradient(low = "white", high = "#B2182B", na.value = "grey90")
    }
}

.plot_scc_nps_marker_summary <- function(pair_df) {
    df <- pair_df |>
        dplyr::group_by(marker) |>
        dplyr::summarise(NPS = max(nps, na.rm = TRUE), .groups = "drop")
    if (nrow(df) == 0) return(NULL)
    df <- df[order(df$NPS, decreasing = TRUE), , drop = FALSE]
    df$marker <- factor(df$marker, levels = df$marker)
    ggplot2::ggplot(df, ggplot2::aes(marker, NPS)) +
        ggplot2::geom_col(fill = "#C44E52", width = 0.7) +
        ggplot2::labs(
            title = "Off-target negative spread by marker",
            subtitle = "Bars show the worst source-SCC NPS for each off-target marker.",
            x = "Off-target unmixed marker",
            y = "Worst NPS"
        ) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))
}

.plot_scc_top_fpr_pairs <- function(pair_df, top_n = 15L) {
    df <- pair_df[is.finite(pair_df$fpr), , drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    df <- df[order(df$fpr, df$nps, decreasing = TRUE), , drop = FALSE]
    df <- df[seq_len(min(as.integer(top_n), nrow(df))), , drop = FALSE]
    df$pair <- paste0(df$target, " -> ", df$marker)
    df$pair <- factor(df$pair, levels = rev(df$pair))
    ggplot2::ggplot(df, ggplot2::aes(pair, 100 * fpr, fill = nps)) +
        ggplot2::geom_col(width = 0.72) +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_gradient(low = "#FEE8C8", high = "#B2182B", na.value = "grey80") +
        ggplot2::labs(
            title = "Top false-positive off-target pairs",
            subtitle = "FPR uses the 99.5th percentile of the matched negative reference as threshold.",
            x = "Source SCC -> off-target marker",
            y = "False-positive rate (%)",
            fill = "NPS"
        ) +
        ggplot2::theme_minimal(base_size = 12.5)
}

.draw_scc_post_unmix_qc_pages <- function(unmixed_list, qc_summary, markers) {
    metrics <- .compute_scc_post_unmix_qc(
        unmixed_list = unmixed_list,
        qc_summary = qc_summary,
        markers = markers
    )
    if (is.null(metrics)) {
        return(invisible(NULL))
    }
    .draw_scc_post_unmix_overview_page(metrics)
    plot_pages <- list(
        .plot_scc_nps_marker_summary(metrics$pairs),
        .plot_scc_pair_heatmap(metrics$pairs, "nps", "Pairwise off-target negative spread", "NPS"),
        .plot_scc_pair_heatmap(metrics$pairs, "fpr", "Pairwise false-positive rate", "FPR (%)", percent = TRUE),
        .plot_scc_top_fpr_pairs(metrics$pairs),
        .plot_scc_pair_heatmap(metrics$pairs, "bias", "Pairwise off-target bias", "Bias", diverging = TRUE),
        .plot_scc_pair_heatmap(metrics$pairs, "slope", "Pairwise target-driven off-target slope", "Slope", diverging = TRUE)
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

    if (nrow(qc_summary) > 0) {
        rows_per_page <- 24
        page_starts <- seq(1, nrow(qc_summary), by = rows_per_page)
        for (start_idx in page_starts) {
            end_idx <- min(start_idx + rows_per_page - 1, nrow(qc_summary))
            block <- qc_summary[start_idx:end_idx, , drop = FALSE]
            grid::grid.newpage()
            grid::grid.text("SCC Summary", x = 0.05, y = 0.95, just = c("left", "top"), gp = grid::gpar(fontsize = 16, fontface = "bold"))
            grid::grid.text(
                paste(.format_scc_summary_lines(block), collapse = "\n"),
                x = 0.05,
                y = 0.9,
                just = c("left", "top"),
                gp = grid::gpar(fontsize = 9, fontfamily = "mono", lineheight = 1.15)
            )
        }
    }

    .draw_af_bank_qc_pages(M_built, attr(M_built, "af_bank_info"), pd = pd)

    keep_non_af <- !grepl("^AF($|_)", rownames(M_report), ignore.case = TRUE)
    M_no_af <- M_report[keep_non_af, , drop = FALSE]

    if (nrow(M_no_af) > 0) {
        .draw_report_ggplot_page(plot_spectra(M_no_af, pd = pd, output_file = NULL), height_ratio = 0.6)
    }

    if (nrow(M_no_af) > 1) {
        sim_mat <- calculate_similarity_matrix(M_no_af)
        .draw_report_ggplot_page(plot_similarity_matrix(sim_mat, output_file = NULL))
        ssm_method <- if (method %in% c("NNLS", "RWLS", "AutoSpectral")) {
            if (identical(method, "RWLS")) "WLS" else "OLS"
        } else {
            method
        }
        .draw_report_ggplot_page(plot_ssm(calculate_ssm(M_no_af, method = ssm_method), output_file = NULL))

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
            markers = rownames(M_no_af)
        )
    }

    if (nrow(qc_summary) > 0 && !is.null(report_plot_dir) && dir.exists(report_plot_dir)) {
        for (i in seq_len(nrow(qc_summary))) {
            .draw_scc_report_sample_page(qc_summary[i, , drop = FALSE], report_plot_dir = report_plot_dir, use_scatter_gating = use_scatter_gating)
        }
    }

    message("SCC report saved to: ", output_file)
    invisible(list(M = M_built, qc_summary = qc_summary, qc_plot_dir = retained_qc_plot_dir, af_bank_info = attr(M_built, "af_bank_info"), method = method))
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
