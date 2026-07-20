.sample_pdf_draw_plot <- function(state, plot, square = FALSE, height_ratio = 1, width_ratio = 1) {
    .with_known_qc_plot_warnings_suppressed(
        .draw_qc_report_plot_page(
            plot,
            square = square,
            height_ratio = height_ratio,
            width_ratio = width_ratio,
            newpage = state$page_started
        )
    )
    state$page_started <- TRUE
    invisible(NULL)
}

.sample_pdf_draw_spectra <- function(state, plot) {
    .with_known_qc_plot_warnings_suppressed(
        .draw_qc_report_spectra_page(plot, newpage = state$page_started)
    )
    state$page_started <- TRUE
    invisible(NULL)
}

.sample_pdf_reference_pages <- function(M, pd, qc_metrics_dir, plot_dir, matrix_markers_per_page, state) {
    keep <- !grepl("^AF($|_)", rownames(M), ignore.case = TRUE)
    fluor <- M[keep, , drop = FALSE]
    af <- M[!keep, , drop = FALSE]
    if (!is.null(qc_metrics_dir)) {
        if (nrow(fluor) > 0) {
            .write_qc_report_matrix_metric(
                fluor, file.path(qc_metrics_dir, "reference_spectra.csv"), row_id = "fluorophore"
            )
        }
        if (nrow(af) > 0) {
            .write_qc_report_matrix_metric(
                af, file.path(qc_metrics_dir, "af_bank_spectra.csv"), row_id = "af_band"
            )
        }
    }
    .spectreasy_console_step("Spectra overlay")
    if (nrow(fluor) > 0) {
        spectra <- plot_spectra(fluor, pd = pd, output_file = NULL)
        .save_qc_report_png(spectra, plot_dir, "spectra_overlay.png")
        .sample_pdf_draw_spectra(state, spectra)
    }
    if (nrow(fluor) <= 1) return(invisible(NULL))
    .spectreasy_console_step("Similarity matrix")
    similarity <- calculate_similarity_matrix(fluor)
    .write_qc_report_matrix_metric(
        similarity,
        file.path(qc_metrics_dir, "fluorophore_spectral_similarity.csv"),
        row_id = "fluorophore"
    )
    pages <- .build_qc_report_matrix_pages(
        similarity,
        plot_fun = plot_similarity_matrix,
        max_markers_per_page = matrix_markers_per_page,
        item_label = "Markers"
    )
    for (i in seq_along(pages)) {
        .save_qc_report_png(pages[[i]], plot_dir, sprintf("similarity_matrix_%02d.png", i))
        .sample_pdf_draw_plot(state, pages[[i]])
    }
    invisible(NULL)
}

.sample_pdf_nps_pages <- function(results_df, method, qc_metrics_dir, plot_dir, max_files, state) {
    if (identical(method, "NNLS")) {
        .spectreasy_console_step("NPS diagnostics", "skipped for NNLS")
        return(invisible(NULL))
    }
    .spectreasy_console_step("NPS diagnostics")
    scores <- calculate_nps(results_df)
    scores <- scores[!grepl("^AF($|_)", scores$Marker, ignore.case = TRUE), , drop = FALSE]
    if (nrow(scores) == 0) return(invisible(NULL))
    .write_qc_report_csv(scores, file.path(qc_metrics_dir, "negative_population_spread.csv"))
    pages <- .build_qc_report_nps_pages(scores, max_files_per_page = max_files)
    for (i in seq_along(pages)) {
        .save_qc_report_png(pages[[i]], plot_dir, sprintf("negative_population_spread_%02d.png", i))
        .sample_pdf_draw_plot(state, pages[[i]])
    }
    invisible(NULL)
}

.sample_pdf_scatter_pages <- function(results_df, options, plot_dir, state) {
    .spectreasy_console_step("Sample NxN scatter")
    markers <- setdiff(colnames(results_df), .get_result_metadata_columns(colnames(results_df)))
    markers <- markers[!grepl("^AF($|_)", markers, ignore.case = TRUE)]
    pages <- .build_qc_report_sample_scatter_pages(
        results_df = results_df,
        markers = markers,
        rows_per_page = options$sample_nxn_rows_per_page,
        max_points = options$sample_nxn_max_points,
        transform = options$sample_nxn_transform,
        asinh_cofactor = options$sample_nxn_asinh_cofactor,
        axis_limit = options$sample_nxn_axis_limit,
        all_samples = options$nxn_all_samples
    )
    for (i in seq_along(pages)) {
        .save_qc_report_png(pages[[i]], plot_dir, sprintf("sample_nxn_scatter_%02d.png", i))
        .sample_pdf_draw_plot(state, pages[[i]], square = TRUE)
    }
    invisible(NULL)
}

.sample_pdf_residual_pages <- function(res_list, M, pd, method, qc_metrics_dir, plot_dir, max_files, state) {
    if (is.null(res_list)) return(invisible(NULL))
    .spectreasy_console_step("Detector RMS")
    detector_rms <- .compute_qc_report_detector_rms(res_list, M = M, pd = pd, unmixing_method = method)
    if (!is.null(detector_rms) && !is.null(qc_metrics_dir)) {
        .write_qc_report_csv(detector_rms, file.path(qc_metrics_dir, "rms_residual_per_detector.csv"))
    }
    plot <- plot_detector_rms_residuals(res_list, M = M, pd = pd, output_file = NULL, unmixing_method = method)
    if (!is.null(plot)) {
        .save_qc_report_png(plot, plot_dir, "rms_residual_per_detector.png")
        .sample_pdf_draw_plot(state, plot, height_ratio = 0.74, width_ratio = 0.90)
    }
    .spectreasy_console_step("Reconstruction")
    sample_rms <- .compute_qc_report_sample_rms(res_list, M = M, unmixing_method = method)
    if (!is.null(sample_rms) && !is.null(qc_metrics_dir)) {
        .write_qc_report_csv(sample_rms, file.path(qc_metrics_dir, "detector_reconstruction_error.csv"))
    }
    pages <- .build_qc_report_rms_pages(
        res_list, M = M, max_files_per_page = max_files, unmixing_method = method
    )
    for (i in seq_along(pages)) {
        .save_qc_report_png(pages[[i]], plot_dir, sprintf("detector_reconstruction_error_%02d.png", i))
        .sample_pdf_draw_plot(state, pages[[i]], height_ratio = 0.72, width_ratio = 0.72)
    }
    invisible(NULL)
}

.prepare_sample_pdf_results <- function(results, res_list, max_events_per_sample) {
    report_results <- if (is.data.frame(results)) {
        .cap_qc_report_results_df(results, max_events_per_sample = max_events_per_sample)
    } else {
        .cap_qc_report_results_list(results, max_events_per_sample = max_events_per_sample)
    }
    if (is.null(res_list) && is.list(report_results) && !is.data.frame(report_results)) {
        has_residuals <- any(vapply(report_results, function(x) is.list(x) && !is.null(x$residuals), logical(1)))
        if (has_residuals) res_list <- report_results
    } else if (!is.null(res_list)) {
        res_list <- .cap_qc_report_results_list(res_list, max_events_per_sample = max_events_per_sample)
    }
    list(data = .normalize_qc_report_results_df(report_results), residuals = res_list)
}

.qc_samples_pdf <- function(results, M, output_file, method, res_list, pd, options) {
    prepared <- .prepare_sample_pdf_results(results, res_list, options$max_events_per_sample)
    if (!"File" %in% colnames(prepared$data)) stop("results must contain a 'File' column.")
    output_dir <- dirname(output_file)
    if (!is.na(output_dir) && nzchar(output_dir) && output_dir != ".") {
        dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    }
    plot_dir <- .prepare_qc_report_png_dir(options$qc_plot_dir, options$save_qc_pngs, output_file)
    metrics_dir <- .prepare_qc_report_metrics_dir(options$qc_metrics_dir)
    grDevices::pdf(output_file, width = 11, height = 8.5)
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
    state <- new.env(parent = emptyenv())
    state$page_started <- FALSE

    file_counts <- .qc_report_file_counts(results)
    if (!is.null(options$max_events_per_sample) &&
        any(as.numeric(file_counts) > as.numeric(options$max_events_per_sample)[1], na.rm = TRUE)) {
        .spectreasy_console_step("Event cap", paste0(as.integer(options$max_events_per_sample[1]), " events per sample"))
    }
    .sample_pdf_reference_pages(M, pd, metrics_dir, plot_dir, options$matrix_markers_per_page, state)
    .sample_pdf_nps_pages(prepared$data, method, metrics_dir, plot_dir, options$overview_files_per_page, state)
    .sample_pdf_scatter_pages(prepared$data, options, plot_dir, state)
    .sample_pdf_residual_pages(
        prepared$residuals, M, pd, method, metrics_dir, plot_dir,
        options$overview_files_per_page, state
    )
    .spectreasy_console_field("Saved", .spectreasy_console_path(output_file))
    .spectreasy_console_footer(blank = FALSE)
    invisible(list(output_file = output_file, qc_plot_dir = plot_dir, qc_metrics_dir = metrics_dir))
}
