.prepare_control_report_reference <- function(M,
                                              output_format,
                                              qc_summary,
                                              report_plot_dir,
                                              cleanup_report_plot_dir,
                                              qc_plot_dir,
                                              save_qc_pngs,
                                              scc_dir,
                                              control_file,
                                              cytometer,
                                              seed,
                                              extra_args) {
    reuse <- !is.null(M) && (
        identical(output_format, "html") || (!is.null(qc_summary) && !is.null(report_plot_dir))
    )
    cleanup_dirs <- character()
    collector_plot_dir <- NULL
    af_bank_info <- NULL
    if (isTRUE(reuse)) {
        built <- .as_reference_matrix(M, "M")
        if (is.null(report_plot_dir) && identical(output_format, "html")) {
            report_plot_dir <- tempfile("spectreasy_control_html_gates_")
            dir.create(report_plot_dir, recursive = TRUE, showWarnings = FALSE)
        } else {
            report_plot_dir <- normalizePath(report_plot_dir, mustWork = FALSE)
            if (!dir.exists(report_plot_dir)) stop("report_plot_dir not found: ", report_plot_dir, call. = FALSE)
        }
        if (isTRUE(cleanup_report_plot_dir)) cleanup_dirs <- c(cleanup_dirs, report_plot_dir)
    } else {
        plot_info <- .prepare_scc_report_plot_dir(qc_plot_dir = qc_plot_dir, save_qc_pngs = save_qc_pngs)
        collector_plot_dir <- plot_info$plot_dir
        if (!is.null(plot_info$cleanup_dir)) cleanup_dirs <- c(cleanup_dirs, plot_info$cleanup_dir)
        control_resolved <- .resolve_control_file_path(control_file)
        control_input <- if (file.exists(control_resolved)) control_resolved else NULL
        build_args <- c(list(
            input_folder = scc_dir,
            output_folder = plot_info$plot_dir,
            save_qc_plots = TRUE,
            control_df = control_input,
            cytometer = cytometer,
            seed = seed
        ), extra_args)
        built <- do.call(build_reference_matrix, build_args)
        if (is.null(built) || nrow(built) == 0) stop("No valid spectra found while generating the SCC report.")
        report_plot_dir <- attr(built, "qc_plot_dir")
        if (is.null(report_plot_dir) || !dir.exists(report_plot_dir)) report_plot_dir <- plot_info$plot_dir
        af_bank_info <- attr(built, "af_bank_info")
    }
    list(
        M = built,
        report_plot_dir = report_plot_dir,
        collector_plot_dir = collector_plot_dir,
        cleanup_dirs = unique(cleanup_dirs),
        af_bank_info = af_bank_info
    )
}

.normalize_control_report_summary <- function(qc_summary, M) {
    if (is.null(qc_summary)) qc_summary <- attr(M, "qc_summary")
    if (is.null(qc_summary)) return(data.frame())
    qc_summary <- as.data.frame(qc_summary, stringsAsFactors = FALSE)
    qc_summary[order(qc_summary$fluorophore, qc_summary$sample), , drop = FALSE]
}

.resolve_control_report_matrix <- function(M, unmixing_matrix_file, fallback) {
    if (!is.null(M)) return(.as_reference_matrix(M, "M"))
    if (!is.null(unmixing_matrix_file) && file.exists(unmixing_matrix_file)) {
        .stop_if_static_unmixing_matrix_path(unmixing_matrix_file, arg_name = "unmixing_matrix_file")
        return(.as_reference_matrix(.read_unmixing_matrix_csv(unmixing_matrix_file), "M"))
    }
    fallback
}

.resolve_control_report_pd <- function(pd, M, output_format, scc_dir) {
    if (!is.null(pd)) return(pd)
    stored <- attr(M, "detector_pd")
    if (is.data.frame(stored)) return(stored)
    if (identical(output_format, "html") && !dir.exists(scc_dir)) return(NULL)
    if (!dir.exists(scc_dir)) .spectreasy_stop_missing_directory(scc_dir, label = "scc_dir")
    files <- list.files(scc_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
    if (length(files) == 0) {
        if (!identical(output_format, "html")) .spectreasy_stop_empty_fcs_directory(scc_dir, label = "scc_dir")
        return(NULL)
    }
    frame <- .spectreasy_read_fcs(files[1], label = "SCC FCS file")
    flowCore::pData(flowCore::parameters(frame))
}

.qc_controls_html <- function(M,
                              scc_dir,
                              control_file,
                              cytometer,
                              method,
                              unmixed_list,
                              qc_summary,
                              report_plot_dir,
                              pd,
                              af_bank_info,
                              matrix_source,
                              artifact_paths,
                              report_plots,
                              run_settings,
                              collector_plot_dir,
                              qc_plot_dir,
                              save_qc_pngs,
                              seed,
                              unmix_scatter_max_points,
                              unmix_scatter_axis_limit,
                              output_file,
                              overwrite,
                              project_path,
                              qc_metrics_dir = NULL) {
    cleanup_dir <- NULL
    if (is.null(collector_plot_dir)) {
        info <- .prepare_scc_report_plot_dir(qc_plot_dir = qc_plot_dir, save_qc_pngs = save_qc_pngs)
        collector_plot_dir <- info$plot_dir
        cleanup_dir <- info$cleanup_dir
    }
    if (!is.null(cleanup_dir)) on.exit(unlink(cleanup_dir, recursive = TRUE, force = TRUE), add = TRUE)
    report_data <- collect_control_report_data(
        M = M,
        scc_dir = scc_dir,
        control_file = control_file,
        cytometer = cytometer,
        unmixing_method = method,
        unmixed_list = unmixed_list,
        qc_summary = qc_summary,
        report_plot_dir = report_plot_dir,
        pd = pd,
        af_bank_info = af_bank_info,
        matrix_source = matrix_source,
        artifact_paths = artifact_paths,
        plots = report_plots,
        run_settings = c(list(
            unmix_scatter_max_points = unmix_scatter_max_points,
            unmix_scatter_axis_limit = unmix_scatter_axis_limit,
            seed = seed,
            save_qc_pngs = save_qc_pngs,
            report_format = "html"
        ), run_settings),
        plot_dir = collector_plot_dir,
        qc_metrics_dir = qc_metrics_dir,
        project_path = project_path
    )
    if (!isTRUE(save_qc_pngs)) report_data <- .report_embed_plot_manifest(report_data)
    render_qc_html_report(report_data, output_file, overwrite = overwrite)
}
