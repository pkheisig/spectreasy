# Internal helpers for SCC report generation.
.resolve_scc_report_control_input <- function(control_df = NULL, control_file = "fcs_mapping.csv") {
    control_file <- .resolve_control_file_path(control_file)
    control_input <- control_df
    if (is.null(control_input) && file.exists(control_file)) {
        control_input <- control_file
    }
    control_input
}

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
    df <- df[, c(
        "fluorophore",
        "marker",
        "sample",
        "type",
        "peak_channel",
        "n_total",
        "n_scatter_gated",
        "n_final",
        "scatter_gate_pct",
        "histogram_gate_pct"
    ), drop = FALSE]

    headers <- c("Fluor", "Marker", "Sample", "Type", "Peak", "Events", "Scatter", "Final", "Scatter%", "Hist%")
    widths <- c(14, 14, 24, 8, 8, 9, 9, 7, 9, 7)

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
        trim_to(df$n_scatter_gated, widths[7]),
        trim_to(df$n_final, widths[8]),
        trim_to(fmt_num(df$scatter_gate_pct, 1), widths[9]),
        trim_to(fmt_num(df$histogram_gate_pct, 1), widths[10])
    )

    make_line <- function(parts) {
        paste(vapply(seq_along(parts), function(i) sprintf(paste0("%-", widths[i], "s"), parts[[i]]), character(1)), collapse = " ")
    }

    header_line <- make_line(as.list(headers))
    sep_line <- paste(vapply(widths, function(w) paste(rep("-", w), collapse = ""), character(1)), collapse = " ")
    row_lines <- vapply(seq_len(nrow(df)), function(i) make_line(lapply(vals, `[[`, i)), character(1))

    c(header_line, sep_line, row_lines)
}

.draw_report_ggplot_page <- function(p) {
    if (is.null(p)) {
        return(invisible(NULL))
    }
    grid::grid.newpage()
    grid::grid.draw(ggplot2::ggplotGrob(p))
    invisible(NULL)
}

.draw_scc_report_sample_page <- function(row, report_plot_dir) {
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
        " | Final: ", row$n_final[[1]], " (", row$histogram_gate_pct[[1]], "%)"
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
    .draw_report_image_panel(
        file.path(report_plot_dir, "histogram", paste0(sample_id, "_histogram.png")),
        "Peak-Channel Histogram Gate",
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

#' Generate SCC QC Report
#'
#' Builds a reference matrix from single-color controls, assembles a PDF report
#' for pre-unmix review, and can optionally retain the intermediate QC PNG plots.
#'
#' @param M Optional precomputed reference matrix. If omitted, the report uses the
#'   matrix generated from the supplied SCC files.
#' @param scc_dir Directory containing SCC FCS files.
#' @param output_file Path to save the PDF report. Must be supplied explicitly.
#' @param control_df Optional control mapping data.frame or CSV path.
#' @param control_file Control mapping CSV path used when `control_df` is `NULL`.
#' @param cytometer Cytometer name passed to [build_reference_matrix()].
#' @param qc_plot_dir Directory where FSC/SSC, histogram, and spectrum PNGs are written
#'   when `save_qc_pngs = TRUE`.
#' @param save_qc_pngs Logical; if `TRUE`, keep the intermediate QC PNG files in
#'   `qc_plot_dir`. If `FALSE` (default), PNGs are written to a temporary directory
#'   for report assembly and removed afterward.
#' @param include_multi_af Logical; forward to [build_reference_matrix()].
#' @param af_dir AF directory forwarded to [build_reference_matrix()].
#' @param include_ssm Logical; include the spectral spread matrix page.
#' @param seed Optional integer seed for deterministic subsampling/clustering.
#' @param ... Additional arguments forwarded to [build_reference_matrix()].
#' @return Invisibly returns a list with `M`, `qc_summary`, and `qc_plot_dir`.
#'   `qc_plot_dir` is `NULL` unless `save_qc_pngs = TRUE`.
#' @examples
#' if (interactive()) {
#'   generate_scc_report(
#'     scc_dir = "scc",
#'     control_file = "fcs_mapping.csv",
#'     output_file = file.path("spectreasy_outputs", "SCC_QC_Report.pdf")
#'   )
#' }
#' @export
generate_scc_report <- function(
    M = NULL,
    scc_dir = "scc",
    output_file = NULL,
    control_df = NULL,
    control_file = "fcs_mapping.csv",
    cytometer = "Aurora",
    qc_plot_dir = file.path("spectreasy_outputs", "scc_report_plots"),
    save_qc_pngs = FALSE,
    include_multi_af = FALSE,
    af_dir = "af",
    include_ssm = TRUE,
    seed = NULL,
    ...
) {
    if (is.null(output_file) || !nzchar(trimws(as.character(output_file)[1]))) {
        stop("Please supply output_file to save the SCC PDF report.", call. = FALSE)
    }

    message("Generating SCC QC report...")
    out_dir <- dirname(output_file)
    if (!is.na(out_dir) && nzchar(out_dir) && out_dir != ".") {
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    }

    plot_dir_info <- .prepare_scc_report_plot_dir(qc_plot_dir = qc_plot_dir, save_qc_pngs = save_qc_pngs)
    if (!is.null(plot_dir_info$cleanup_dir)) {
        on.exit(unlink(plot_dir_info$cleanup_dir, recursive = TRUE, force = TRUE), add = TRUE)
    }

    control_input <- .resolve_scc_report_control_input(control_df = control_df, control_file = control_file)

    M_built <- build_reference_matrix(
        input_folder = scc_dir,
        output_folder = plot_dir_info$plot_dir,
        save_qc_plots = TRUE,
        control_df = control_input,
        include_multi_af = include_multi_af,
        af_dir = af_dir,
        cytometer = cytometer,
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

    M_report <- if (is.null(M)) M_built else .as_reference_matrix(M, "M")

    fcs_files <- list.files(scc_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
    if (length(fcs_files) == 0) {
        stop("No FCS files found in scc_dir: ", scc_dir)
    }
    ff_meta <- flowCore::read.FCS(fcs_files[1], transformation = FALSE, truncate_max_range = FALSE)
    pd <- flowCore::pData(flowCore::parameters(ff_meta))

    if (!requireNamespace("png", quietly = TRUE)) {
        stop("Package 'png' is required to embed SCC QC plots in the report.")
    }

    grDevices::pdf(output_file, width = 11, height = 8.5)
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

    grid::grid.newpage()
    grid::grid.text("spectreasy: Single-Color Control Review", x = 0.5, y = 0.7, gp = grid::gpar(fontsize = 20, fontface = "bold"))
    grid::grid.text(
        paste0(
            "Generated on: ", Sys.time(), "\n",
            "SCC directory: ", normalizePath(scc_dir, mustWork = FALSE), "\n",
            "Controls processed: ", nrow(qc_summary), "\n",
            "Workflow intent: review SCC quality before autounmix_controls()."
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

    keep_non_af <- !grepl("^AF($|_)", rownames(M_report), ignore.case = TRUE)
    M_no_af <- M_report[keep_non_af, , drop = FALSE]

    if (nrow(M_no_af) > 0) {
        .draw_report_ggplot_page(plot_spectra(M_no_af, pd = pd, output_file = NULL))
    }

    if (isTRUE(include_ssm) && nrow(M_no_af) > 1) {
        .draw_report_ggplot_page(plot_ssm(calculate_ssm(M_no_af), output_file = NULL))
    }

    if (nrow(qc_summary) > 0) {
        for (i in seq_len(nrow(qc_summary))) {
            .draw_scc_report_sample_page(qc_summary[i, , drop = FALSE], report_plot_dir = report_plot_dir)
        }
    }

    message("SCC report saved to: ", output_file)
    invisible(list(M = M_built, qc_summary = qc_summary, qc_plot_dir = retained_qc_plot_dir))
}
