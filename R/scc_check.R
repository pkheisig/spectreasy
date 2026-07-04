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
    gate_dir <- if (isTRUE(use_scatter_gating)) "intensity_scatter" else "histogram"
    gate_suffix <- if (isTRUE(use_scatter_gating)) "_intensity_scatter.png" else "_histogram.png"
    gate_title <- if (isTRUE(use_scatter_gating)) "Intensity-vs-FSC Scatter Gate" else "Peak-Channel Histogram Gate"
    .draw_report_image_panel(
        file.path(report_plot_dir, gate_dir, paste0(sample_id, gate_suffix)),
        gate_title,
        x = grid::unit(0.74, "npc"),
        y = grid::unit(0.71, "npc"),
        width = grid::unit(0.42, "npc"),
        height = grid::unit(0.26, "npc")
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

.scc_reference_overlay_matrix <- function(M) {
    M <- .as_reference_matrix(M, "M")
    af_rows <- grepl("^AF($|_)", rownames(M), ignore.case = TRUE)
    if (sum(af_rows, na.rm = TRUE) <= 1L) {
        return(M)
    }
    M[!af_rows, , drop = FALSE]
}

#' Generate SCC QC Report
#'
#' Builds a reference matrix from single-color controls, assembles a PDF report
#' for pre-unmix review, and can optionally retain the intermediate QC PNG plots.
#'
#' @param M Optional precomputed reference matrix. If omitted, the report uses the
#'   matrix generated from the supplied SCC files.
#' @param unmixing_matrix_file Optional CSV path to a saved reference matrix.
#'   Used when `M` is not supplied. By default this points to the reference matrix
#'   produced by [unmix_controls()] (`"scc_reference_matrix.csv"`).
#' @param scc_dir Directory containing SCC FCS files.
#' @param output_file Path to save the PDF report. Defaults to `"spectreasy_outputs/unmix_controls/qc_controls_report.pdf"`.
#' @param control_file Control mapping CSV path.
#' @param cytometer Cytometer name passed to [build_reference_matrix()]. The
#'   default, `"auto"`, infers the cytometer from FCS detector names when possible.
#' @param unmixing_method Unmixing method used for the control scatter matrix
#'   (`"WLS"`, `"RWLS"`, `"OLS"`, or `"NNLS"`).
#' @param qc_plot_dir Directory where FSC/SSC, intensity-gate, and spectrum PNGs are written
#'   when `save_qc_pngs = TRUE`.
#' @param save_qc_pngs Logical; if `TRUE`, keep the intermediate QC PNG files in
#'   `qc_plot_dir`. If `FALSE` (default), PNGs are written to a temporary directory
#'   for report assembly and removed afterward.
#' @param qc_metrics_dir Optional directory where plot-ready SCC QC metric
#'   CSVs are written alongside the PDF report.
#' @param unmixed_list Optional precomputed SCC unmixing results from
#'   [unmix_controls()]. When supplied with `M`, the report reuses these results
#'   instead of unmixing SCC files again.
#' @param qc_summary Optional precomputed SCC gating summary from
#'   `attr(M, "qc_summary")`.
#' @param report_plot_dir Optional directory containing precomputed SCC
#'   FSC/SSC, gate, and spectrum PNGs from [build_reference_matrix()].
#' @param pd Optional detector metadata (`flowCore::pData(parameters(ff))`).
#' @param af_bank_info Optional precomputed AF bank metadata from
#'   `attr(M, "af_bank_info")`.
#' @param cleanup_report_plot_dir Logical; if `TRUE`, remove
#'   `report_plot_dir` after the PDF is assembled. Intended for internal
#'   temporary plot directories created by [unmix_controls()].
#' @param use_scatter_gating Logical; if `TRUE` (default), use scatter/intensity
#'   gating and show scatter gate plots in the report. If `FALSE`, use and show
#'   the legacy histogram gate.
#' @param unmix_scatter_max_points Maximum events sampled per control for the
#'   SCC unmixing scatter matrix.
#' @param unmix_scatter_axis_limit Optional fixed symmetric axis limit for the
#'   SCC unmixing scatter matrix. The default `NULL` uses local per-panel
#'   ranges. Use `1e5` for `c(-1e5, 1e5)` on every panel.
#' @param seed Optional integer seed for deterministic subsampling/clustering.
#' @param ... Additional arguments forwarded to [build_reference_matrix()].
#' @return Invisibly returns a list with `M`, `qc_summary`, `qc_plot_dir`,
#'   `qc_metrics_dir`, `af_bank_info`, and `unmixing_method`.
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
    M = NULL,
    unmixing_matrix_file = file.path("spectreasy_outputs", "unmix_controls", "scc_reference_matrix.csv"),
    scc_dir = "scc",
    output_file = "spectreasy_outputs/unmix_controls/qc_controls_report.pdf",
    control_file = "fcs_mapping.csv",
    cytometer = "auto",
    unmixing_method = "WLS",
    qc_plot_dir = file.path("spectreasy_outputs", "scc_report_plots"),
    save_qc_pngs = FALSE,
    qc_metrics_dir = NULL,
    unmixed_list = NULL,
    qc_summary = NULL,
    report_plot_dir = NULL,
    pd = NULL,
    af_bank_info = NULL,
    cleanup_report_plot_dir = FALSE,
    use_scatter_gating = TRUE,
    unmix_scatter_max_points = 1000,
    unmix_scatter_axis_limit = NULL,
    seed = NULL,
    ...
) {
    if (is.null(output_file) || !nzchar(trimws(as.character(output_file)[1]))) {
        stop("Please supply output_file to save the SCC PDF report.", call. = FALSE)
    }
    extra_args <- list(...)
    unmixing_method <- .normalize_unmix_method(unmixing_method)

    message("Generating SCC QC report...")
    out_dir <- dirname(output_file)
    if (!is.na(out_dir) && nzchar(out_dir) && out_dir != ".") {
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    }

    qc_metrics_dir <- .prepare_qc_report_metrics_dir(qc_metrics_dir)

    can_reuse_reference <- !is.null(M) && !is.null(qc_summary) && !is.null(report_plot_dir)
    if (isTRUE(can_reuse_reference)) {
        M_built <- .as_reference_matrix(M, "M")
        report_plot_dir <- normalizePath(report_plot_dir, mustWork = FALSE)
        if (!dir.exists(report_plot_dir)) {
            stop("report_plot_dir not found: ", report_plot_dir, call. = FALSE)
        }
        if (isTRUE(cleanup_report_plot_dir)) {
            on.exit(unlink(report_plot_dir, recursive = TRUE, force = TRUE), add = TRUE)
        }
    } else {
        plot_dir_info <- .prepare_scc_report_plot_dir(qc_plot_dir = qc_plot_dir, save_qc_pngs = save_qc_pngs)
        if (!is.null(plot_dir_info$cleanup_dir)) {
            on.exit(unlink(plot_dir_info$cleanup_dir, recursive = TRUE, force = TRUE), add = TRUE)
        }

        control_resolved <- .resolve_control_file_path(control_file)
        control_input <- if (file.exists(control_resolved)) control_resolved else NULL

        build_args <- c(
            list(
                input_folder = scc_dir,
                output_folder = plot_dir_info$plot_dir,
                save_qc_plots = TRUE,
                control_df = control_input,
                cytometer = cytometer,
                use_scatter_gating = use_scatter_gating,
                seed = seed
            ),
            extra_args
        )
        M_built <- do.call(build_reference_matrix, build_args)
        if (is.null(M_built) || nrow(M_built) == 0) {
            stop("No valid spectra found while generating the SCC report.")
        }

        report_plot_dir <- attr(M_built, "qc_plot_dir")
        if (is.null(report_plot_dir) || !dir.exists(report_plot_dir)) {
            report_plot_dir <- plot_dir_info$plot_dir
        }
        af_bank_info <- attr(M_built, "af_bank_info")
    }

    if (is.null(qc_summary)) {
        qc_summary <- attr(M_built, "qc_summary")
    }
    if (is.null(qc_summary)) {
        qc_summary <- data.frame()
    } else {
        qc_summary <- as.data.frame(qc_summary, stringsAsFactors = FALSE)
        qc_summary <- qc_summary[order(qc_summary$fluorophore, qc_summary$sample), , drop = FALSE]
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

    if (is.null(pd)) {
        fcs_files <- list.files(scc_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
        if (length(fcs_files) == 0) {
            stop("No FCS files found in scc_dir: ", scc_dir)
        }
        ff_meta <- flowCore::read.FCS(fcs_files[1], transformation = FALSE, truncate_max_range = FALSE)
        pd <- flowCore::pData(flowCore::parameters(ff_meta))
    }

    if (!requireNamespace("png", quietly = TRUE)) {
        stop("Package 'png' is required to embed SCC QC plots in the report.")
    }

    grDevices::pdf(output_file, width = 11, height = 8.5)
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

    if (is.null(af_bank_info)) {
        af_bank_info <- attr(M_built, "af_bank_info")
    }
    .draw_af_bank_qc_pages(M_built, af_bank_info, pd = pd)

    af_rows <- grepl("^AF($|_)", rownames(M_report), ignore.case = TRUE)
    M_no_af <- M_report[!af_rows, , drop = FALSE]
    M_af <- M_report[af_rows, , drop = FALSE]
    M_reference_overlay <- .scc_reference_overlay_matrix(M_report)

    if (!is.null(qc_metrics_dir)) {
        if (nrow(M_no_af) > 0) {
            .write_qc_report_matrix_metric(
                M_no_af,
                file.path(qc_metrics_dir, "reference_spectra.csv"),
                row_id = "fluorophore"
            )
        }
        if (nrow(M_af) > 0) {
            .write_qc_report_matrix_metric(
                M_af,
                file.path(qc_metrics_dir, "af_bank_spectra.csv"),
                row_id = "af_band"
            )
        }
    }

    if (nrow(M_reference_overlay) > 0) {
        .draw_report_ggplot_page(plot_spectra(M_reference_overlay, pd = pd, output_file = NULL), height_ratio = 0.6)
    }

    if (nrow(M_no_af) > 1) {
        sim_mat <- calculate_similarity_matrix(M_no_af)
        .write_qc_report_matrix_metric(
            sim_mat,
            file.path(qc_metrics_dir, "fluorophore_spectral_similarity.csv"),
            row_id = "fluorophore"
        )
        .draw_report_ggplot_page(plot_similarity_matrix(sim_mat, output_file = NULL))
        if (is.null(unmixed_list)) {
            unmixed_list <- unmix_samples(
                sample_dir = scc_dir,
                M = M_report,
                unmixing_method = unmixing_method,
                write_fcs = FALSE,
                save_report = FALSE,
                verbose = FALSE
            )
        }
        if (!identical(unmixing_method, "NNLS")) {
            control_results_df <- .normalize_qc_report_results_df(unmixed_list)
            nps_scores <- calculate_nps(control_results_df)
            nps_scores <- nps_scores[!grepl("^AF($|_)", nps_scores$Marker, ignore.case = TRUE), , drop = FALSE]
            if (nrow(nps_scores) > 0) {
                .write_qc_report_csv(
                    nps_scores,
                    file.path(qc_metrics_dir, "negative_population_spread.csv")
                )
            }
        }

        detector_rms <- .compute_qc_report_detector_rms(unmixed_list, M = M_report, pd = pd, unmixing_method = unmixing_method)
        if (!is.null(detector_rms)) {
            .write_qc_report_csv(
                detector_rms,
                file.path(qc_metrics_dir, "rms_residual_per_detector.csv")
            )
        }
        sample_rms <- .compute_qc_report_sample_rms(unmixed_list, M = M_report, unmixing_method = unmixing_method)
        if (!is.null(sample_rms)) {
            .write_qc_report_csv(
                sample_rms,
                file.path(qc_metrics_dir, "detector_reconstruction_error.csv")
            )
        }

        sample_to_marker <- NULL
        marker_display <- NULL
        if (nrow(qc_summary) > 0) {
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
        }
        scatter_markers <- rownames(M_no_af)
        scatter_pages <- .build_qc_report_control_scatter_pages(
            unmixed_list = unmixed_list,
            sample_to_marker = sample_to_marker,
            markers = scatter_markers,
            marker_display = NULL,
            rows_per_page = 10,
            max_points_per_sample = unmix_scatter_max_points,
            axis_limit = unmix_scatter_axis_limit,
            seed = seed
        )
        for (p_scatter in scatter_pages) {
            .draw_report_ggplot_page(p_scatter, square = TRUE)
        }
    }

    if (nrow(qc_summary) > 0) {
        for (i in seq_len(nrow(qc_summary))) {
            .draw_scc_report_sample_page(qc_summary[i, , drop = FALSE], report_plot_dir = report_plot_dir, use_scatter_gating = use_scatter_gating)
        }
    }

    message("SCC report saved to: ", output_file)
    invisible(list(
        M = M_built,
        qc_summary = qc_summary,
        qc_plot_dir = retained_qc_plot_dir,
        qc_metrics_dir = qc_metrics_dir,
        af_bank_info = af_bank_info,
        unmixing_method = unmixing_method
    ))
}
