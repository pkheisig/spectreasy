.qc_samples_per_sample_pdf <- function(results, M, output_file, unmixing_method, res_list, pd, options) {
    sample_names <- .qc_report_sample_names(results)
    if (!length(sample_names)) stop("No sample names were available for per-sample PDF reports.", call. = FALSE)
    sample_slugs <- .report_filename_slug(sample_names, fallback = "sample")
    stem <- tools::file_path_sans_ext(basename(output_file))
    sample_files <- file.path(dirname(output_file), paste0(stem, "_", sample_slugs, ".pdf"))
    sample_results <- lapply(seq_along(sample_names), function(index) {
        qc_samples(
            results = .subset_qc_report_sample(results, sample_names[index]),
            M = M,
            output_file = sample_files[index],
            unmixing_method = unmixing_method,
            res_list = .subset_qc_report_sample(res_list, sample_names[index]),
            pd = pd,
            max_events_per_sample = options$max_events_per_sample,
            overview_files_per_page = options$overview_files_per_page,
            matrix_markers_per_page = options$matrix_markers_per_page,
            sample_nxn_rows_per_page = options$sample_nxn_rows_per_page,
            sample_nxn_max_points = options$sample_nxn_max_points,
            sample_nxn_transform = options$sample_nxn_transform,
            sample_nxn_asinh_cofactor = options$sample_nxn_asinh_cofactor,
            sample_nxn_axis_limit = options$sample_nxn_axis_limit,
            nxn_all_samples = TRUE,
            qc_plot_dir = if (is.null(options$qc_plot_dir)) NULL else file.path(options$qc_plot_dir, sample_slugs[index]),
            save_qc_pngs = options$save_qc_pngs,
            qc_metrics_dir = if (is.null(options$qc_metrics_dir)) NULL else file.path(options$qc_metrics_dir, sample_slugs[index]),
            report_format = "pdf",
            report_per_sample = FALSE,
            overwrite = options$overwrite,
            report_run_settings = c(list(report_per_sample = TRUE), options$report_run_settings),
            report_artifact_paths = options$report_artifact_paths,
            project_path = options$project_path
        )
    })
    sample_files <- stats::setNames(normalizePath(sample_files, mustWork = TRUE), sample_names)
    invisible(list(
        output_file = sample_files,
        qc_plot_dir = lapply(sample_results, `[[`, "qc_plot_dir"),
        qc_metrics_dir = lapply(sample_results, `[[`, "qc_metrics_dir"),
        ai_qc_prompt_path = lapply(sample_results, `[[`, "ai_qc_prompt_path"),
        ai_qc_data_paths = lapply(sample_results, `[[`, "ai_qc_data_paths")
    ))
}

.sample_report_collect_args <- function(results, M, unmixing_method, res_list, pd, matrix_source, options) {
    list(
        results = results,
        M = M,
        unmixing_method = unmixing_method,
        res_list = res_list,
        pd = pd,
        matrix_source = matrix_source,
        qc_metrics_dir = options$qc_metrics_dir,
        artifact_paths = options$report_artifact_paths,
        run_settings = c(
            list(
                max_events_per_sample = options$max_events_per_sample,
                overview_files_per_page = options$overview_files_per_page,
                matrix_markers_per_page = options$matrix_markers_per_page,
                sample_nxn_rows_per_page = options$sample_nxn_rows_per_page,
                sample_nxn_max_points = options$sample_nxn_max_points,
                sample_nxn_transform = options$sample_nxn_transform,
                sample_nxn_asinh_cofactor = options$sample_nxn_asinh_cofactor,
                sample_nxn_axis_limit = options$sample_nxn_axis_limit,
                nxn_all_samples = options$nxn_all_samples,
                report_per_sample = options$report_per_sample,
                save_qc_pngs = options$save_qc_pngs,
                report_format = "html"
            ),
            options$report_run_settings
        ),
        max_events_per_sample = options$max_events_per_sample,
        overview_files_per_page = options$overview_files_per_page,
        matrix_markers_per_page = options$matrix_markers_per_page,
        sample_nxn_rows_per_page = options$sample_nxn_rows_per_page,
        sample_nxn_max_points = options$sample_nxn_max_points,
        sample_nxn_transform = options$sample_nxn_transform,
        sample_nxn_asinh_cofactor = options$sample_nxn_asinh_cofactor,
        sample_nxn_axis_limit = options$sample_nxn_axis_limit,
        nxn_all_samples = options$nxn_all_samples,
        report_per_sample = options$report_per_sample,
        project_path = options$project_path
    )
}

.qc_samples_html <- function(results, M, output_file, unmixing_method, res_list, pd, matrix_source, options) {
    plot_dir <- if (isTRUE(options$save_qc_pngs)) {
        .prepare_qc_report_png_dir(
            qc_plot_dir = options$qc_plot_dir,
            save_qc_pngs = TRUE,
            output_file = output_file
        )
    } else {
        path <- tempfile("spectreasy_sample_html_plots_")
        dir.create(path, showWarnings = FALSE, recursive = TRUE)
        path
    }
    if (!isTRUE(options$save_qc_pngs)) on.exit(unlink(plot_dir, recursive = TRUE, force = TRUE), add = TRUE)
    collect_args <- .sample_report_collect_args(
        results, M, unmixing_method, res_list, pd, matrix_source, options
    )
    collect_args$plot_dir <- plot_dir
    report_data <- do.call(collect_sample_report_data, collect_args)

    if (isTRUE(options$report_per_sample)) {
        sample_names <- .qc_report_sample_names(results)
        sample_slugs <- .report_filename_slug(sample_names, fallback = "sample")
        sample_reports <- lapply(seq_along(sample_names), function(index) {
            sample_options <- options
            sample_options$qc_metrics_dir <- if (is.null(options$qc_metrics_dir)) NULL else file.path(options$qc_metrics_dir, sample_slugs[index])
            sample_options$report_run_settings <- c(list(report_per_sample = TRUE), options$report_run_settings)
            sample_options$nxn_all_samples <- FALSE
            sample_options$report_per_sample <- FALSE
            args <- .sample_report_collect_args(
                .subset_qc_report_sample(results, sample_names[index]),
                M,
                unmixing_method,
                .subset_qc_report_sample(res_list, sample_names[index]),
                pd,
                matrix_source,
                sample_options
            )
            args$plot_dir <- file.path(plot_dir, "samples", sample_slugs[index])
            sample_data <- do.call(collect_sample_report_data, args)
            if (!isTRUE(options$save_qc_pngs)) sample_data <- .report_embed_plot_manifest(sample_data)
            sample_data
        })
        names(sample_reports) <- sample_names
        report_data$sample_reports <- sample_reports
    }
    if (!isTRUE(options$save_qc_pngs)) report_data <- .report_embed_plot_manifest(report_data)
    render_qc_html_report(report_data, output_file, overwrite = options$overwrite)
}
