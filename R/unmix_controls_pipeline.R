.resolve_unmix_controls_gating <- function(gating_mode,
                                            manual_gate_file,
                                            manual_gate_file_explicit,
                                            scc_dir,
                                            control_file) {
    supplied <- length(manual_gate_file) > 0L &&
        !is.na(manual_gate_file[1]) &&
        nzchar(trimws(as.character(manual_gate_file)[1]))
    manual_gate_file <- if (supplied) normalizePath(as.character(manual_gate_file)[1], mustWork = FALSE) else NULL
    exists <- !is.null(manual_gate_file) && file.exists(manual_gate_file)

    if (identical(gating_mode, "interactive") && interactive()) {
        manual_gate_file <- gate_controls(
            scc_dir = scc_dir,
            control_file = control_file,
            gate_file = manual_gate_file,
            open_browser = TRUE
        )
        return(list(mode = gating_mode, file = manual_gate_file))
    }
    if (identical(gating_mode, "interactive")) {
        if (exists && manual_gate_file_explicit) {
            warning(
                "gating_mode = 'interactive' requires an interactive R session; reusing the existing manual_gate_file instead.",
                call. = FALSE
            )
            .spectreasy_console_field("Gate CSV", .spectreasy_console_path(manual_gate_file))
            return(list(mode = "reuse", file = manual_gate_file))
        }
        if (manual_gate_file_explicit) {
            stop(
                "gating_mode = 'interactive' requires an interactive R session and the supplied manual_gate_file does not exist: ",
                manual_gate_file,
                call. = FALSE
            )
        }
        stop(
            "gating_mode = 'interactive' requires an interactive R session. ",
            "Run interactively to review gates, explicitly provide manual_gate_file, ",
            "or set gating_mode = 'automatic'.",
            call. = FALSE
        )
    }
    if (identical(gating_mode, "reuse")) {
        if (!exists) {
            stop(
                "gating_mode = 'reuse' requires an existing manual_gate_file: ",
                if (is.null(manual_gate_file)) "<not supplied>" else manual_gate_file,
                call. = FALSE
            )
        }
        .spectreasy_console_field("Gate CSV", .spectreasy_console_path(manual_gate_file))
        return(list(mode = gating_mode, file = manual_gate_file))
    }
    .spectreasy_console_field("Gating", "automatic (GMM)")
    list(mode = gating_mode, file = NULL)
}

.unmix_controls_report_paths <- function(output_paths, save_report, report_format) {
    if (!isTRUE(save_report)) return(list(directory = NULL, output_file = NULL))
    directory <- .next_safe_output_dir(output_paths$qc_controls_dir)
    output_file <- file.path(directory, paste0("qc_controls_report.", report_format))
    .spectreasy_console_field("Report", .spectreasy_console_path(output_file))
    list(directory = directory, output_file = output_file)
}

.build_unmix_controls_reference <- function(scc_dir,
                                            output_dir,
                                            save_qc_png,
                                            save_report,
                                            control_df,
                                            cytometer,
                                            af_profile,
                                            af_n_bands,
                                            manual_gate_file,
                                            unmixing_method,
                                            scc_background_args,
                                            autospectral_n_candidates,
                                            autospectral_n_spectral,
                                            autospectral_min_events,
                                            autospectral_refine,
                                            n_threads,
                                            seed,
                                            extra_args) {
    save_plots <- isTRUE(save_qc_png) || isTRUE(save_report)
    plot_dir <- if (isTRUE(save_qc_png)) {
        output_dir
    } else if (isTRUE(save_report)) {
        tempfile("scc_report_plots_")
    } else {
        output_dir
    }
    args <- c(list(
        input_folder = scc_dir,
        output_folder = plot_dir,
        save_qc_plots = save_plots,
        control_df = control_df,
        cytometer = cytometer,
        af_profile = af_profile,
        af_n_bands = af_n_bands,
        manual_gate_file = manual_gate_file,
        unmixing_method = unmixing_method,
        scc_background_method = scc_background_args$method,
        scc_background_k = scc_background_args$k,
        autospectral_n_candidates = autospectral_n_candidates,
        autospectral_n_spectral = autospectral_n_spectral,
        autospectral_min_events = autospectral_min_events,
        refine = autospectral_refine,
        n_threads = n_threads,
        seed = seed
    ), extra_args)
    M <- do.call(build_reference_matrix, args)
    if (is.null(M) || nrow(M) == 0) stop("No valid spectra found while building reference matrix.")
    M
}

.learn_unmix_controls_variants <- function(M,
                                           enabled,
                                           output_file,
                                           som_nodes,
                                           cosine_threshold,
                                           max_variants,
                                           min_events,
                                           seed) {
    if (!isTRUE(enabled)) return(list(M = M, library = NULL, file = NULL))
    library <- tryCatch(
        .learn_spectral_variant_library(
            M = M,
            enabled = TRUE,
            som_nodes = som_nodes,
            cosine_threshold = cosine_threshold,
            max_variants = max_variants,
            min_events = min_events,
            seed = seed,
            warn = TRUE
        ),
        error = function(e) {
            warning(
                "Spectral-variant optimization could not build a variant library; using the base reference matrix. Reason: ",
                conditionMessage(e),
                call. = FALSE
            )
            NULL
        }
    )
    if (is.null(library)) return(list(M = M, library = NULL, file = NULL))
    attr(M, "spectral_variant_library") <- library
    .save_spectral_variant_library(library, output_file)
    list(M = M, library = library, file = output_file)
}

.plot_unmix_controls_spectra <- function(M, pd, save_qc_png, output_paths) {
    af_rows <- grepl("^AF($|_)", rownames(M), ignore.case = TRUE)
    fluor <- M[!af_rows, , drop = FALSE]
    af <- M[af_rows, , drop = FALSE]
    spectra <- if (nrow(fluor) > 0) {
        .run_optional_unmix_artifact(
            "SCC spectra plot",
            plot_spectra(fluor, pd = pd, output_file = if (isTRUE(save_qc_png)) output_paths$spectra_file else NULL)
        )
    } else NULL
    af_spectra <- if (nrow(af) > 0) {
        .run_optional_unmix_artifact(
            "AF spectra plot",
            plot_spectra(af, pd = pd, output_file = if (isTRUE(save_qc_png)) output_paths$af_spectra_file else NULL)
        )
    } else NULL
    list(spectra = spectra, af_spectra = af_spectra, af_matrix = af)
}

.unmix_control_files <- function(meta_info,
                                 M,
                                 control_file,
                                 unmixing_method,
                                 rwls_max_iter,
                                 n_threads,
                                 spectral_variant_library,
                                 spectral_variant_top_k,
                                 unmixed_controls_dir) {
    args <- list(
        sample_dir = meta_info$fcs_files,
        M = M,
        control_file = control_file,
        unmixing_method = unmixing_method,
        rwls_max_iter = rwls_max_iter,
        n_threads = n_threads,
        spectral_variant_library = spectral_variant_library,
        spectral_variant_top_k = spectral_variant_top_k,
        output_dir = .as_resolved_unmixed_fcs_dir(unmixed_controls_dir),
        write_fcs = TRUE,
        save_report = FALSE,
        verbose = FALSE
    )
    .spectreasy_console_field("Unmixing", paste0(length(meta_info$fcs_files), " control file(s)"))
    do.call(unmix_samples, args)
}

.unmix_controls_scatter_plot <- function(M,
                                         unmixed_list,
                                         control_df,
                                         save_qc_png,
                                         output_paths,
                                         panel_size,
                                         seed) {
    mapping <- .resolve_unmix_marker_mappings(control_df)
    markers <- rownames(M)
    markers <- markers[!grepl("^AF_", markers, ignore.case = TRUE)]
    .run_optional_unmix_artifact(
        "SCC unmixing scatter plot",
        plot_unmixing_scatter_matrix(
            unmixed_list = unmixed_list,
            sample_to_marker = mapping$sample_to_marker,
            markers = markers,
            marker_display = NULL,
            output_file = if (isTRUE(save_qc_png)) output_paths$unmixing_scatter_png else NULL,
            transform = "none",
            panel_size_mm = panel_size,
            seed = seed
        )
    )
}

.run_unmix_controls_report <- function(save_report,
                                       M,
                                       scc_dir,
                                       report_paths,
                                       control_file,
                                       cytometer,
                                       unmixing_method,
                                       save_qc_png,
                                       unmixed_list,
                                       meta_info,
                                       seed,
                                       output_paths,
                                       spectral_variant_library_file,
                                       manual_gate_file,
                                       report_format,
                                       report_plots,
                                       run_settings,
                                       project_path) {
    if (!isTRUE(save_report)) return(NULL)
    report <- .run_optional_unmix_artifact(
        "Automatic SCC QC report",
        qc_controls(
            M = M,
            scc_dir = scc_dir,
            output_file = report_paths$output_file,
            control_file = control_file,
            cytometer = cytometer,
            unmixing_method = unmixing_method,
            qc_plot_dir = report_paths$directory,
            save_qc_pngs = save_qc_png,
            qc_metrics_dir = report_paths$directory,
            unmixed_list = unmixed_list,
            qc_summary = attr(M, "qc_summary"),
            report_plot_dir = attr(M, "qc_plot_dir"),
            pd = meta_info$pd,
            af_bank_info = attr(M, "af_bank_info"),
            cleanup_report_plot_dir = !isTRUE(save_qc_png),
            unmix_scatter_max_points = 1000,
            seed = seed,
            unmixing_matrix_file = output_paths$reference_matrix_csv,
            report_format = report_format,
            overwrite = "overwrite",
            report_plots = report_plots,
            report_artifact_paths = list(
                reference_matrix = output_paths$reference_matrix_csv,
                detector_noise = output_paths$detector_noise_csv,
                spectral_variant_library = spectral_variant_library_file,
                unmixing_matrix = output_paths$unmixing_matrix_csv,
                gate_file = manual_gate_file
            ),
            report_run_settings = run_settings,
            project_path = project_path
        )
    )
    rendered <- if (is.list(report)) report$output_file else NULL
    if (is.null(rendered) || !file.exists(rendered)) {
        warning(
            "Automatic SCC QC report was requested but was not created at: ",
            report_paths$output_file,
            ". Continuing unmixing.",
            call. = FALSE
        )
    }
    report
}

.unmix_controls_result <- function(M,
                                   W,
                                   unmixed_list,
                                   output_paths,
                                   control_file,
                                   report,
                                   report_paths,
                                   save_report,
                                   save_qc_png,
                                   output_dir,
                                   unmixed_controls_dir,
                                   gating_mode,
                                   manual_gate_file,
                                   static_info,
                                   plots,
                                   spectral_variant_library,
                                   spectral_variant_library_file,
                                   autospectral_refine) {
    invisible(list(
        M = M,
        W = W,
        unmixed_list = unmixed_list,
        reference_matrix_file = output_paths$reference_matrix_csv,
        detector_noise_file = output_paths$detector_noise_csv,
        spectral_variant_library = spectral_variant_library,
        spectral_variant_library_file = spectral_variant_library_file,
        spectral_variant_info = if (!is.null(spectral_variant_library)) spectral_variant_library$info else NULL,
        unmixing_matrix_file = output_paths$unmixing_matrix_csv,
        control_file = control_file,
        control_mapping_file = output_paths$control_mapping_csv,
        qc_report_file = if (isTRUE(save_report) && is.list(report) && !is.null(report$output_file) && file.exists(report$output_file)) report$output_file else NULL,
        qc_nxn_files = if (isTRUE(save_report) && is.list(report)) report$companion_files else NULL,
        qc_controls_dir = if (isTRUE(save_report)) report_paths$directory else NULL,
        qc_metrics_dir = if (isTRUE(save_report)) report_paths$directory else NULL,
        qc_report = report,
        qc_report_data = if (is.list(report)) report$report_data else NULL,
        spectra_file = if (isTRUE(save_qc_png)) output_paths$spectra_file else NULL,
        af_spectra_file = if (isTRUE(save_qc_png) && nrow(plots$af_matrix) > 0) output_paths$af_spectra_file else NULL,
        unmixing_scatter_file = if (isTRUE(save_qc_png)) output_paths$unmixing_scatter_png else NULL,
        qc_plot_dir = if (isTRUE(save_qc_png)) output_dir else NULL,
        unmixed_fcs_dir = unmixed_controls_dir,
        gating_mode = gating_mode,
        manual_gate_file = manual_gate_file,
        static_unmixing_matrix_method = static_info$static_unmixing_matrix_method,
        spectra_plot = plots$spectra,
        af_spectra_plot = plots$af_spectra,
        unmixing_scatter_plot = plots$scatter,
        autospectral_refine = autospectral_refine
    ))
}
