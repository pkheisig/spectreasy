#' Unmix Experimental Samples
#'
#' @param sample_dir Directory containing experimental FCS files, a character
#'   vector of FCS file paths, a `flowCore::flowSet`, or a
#'   `SingleCellExperiment` for in-memory workflows.
#' @param samples_dir Alias for `sample_dir`, accepted for workflows and scripts
#'   that use the plural directory name.
#' @param M Optional reference matrix (Markers x Detectors). If supplied,
#'   unmixing is computed dynamically using this matrix. If not supplied,
#'   it is loaded from the CSV path provided in `unmixing_matrix_file`.
#' @param unmixing_matrix_file Optional CSV path to a saved reference matrix.
#'   Used when `M` is not supplied. By default this points to the reference matrix
#'   produced by [unmix_controls()] (`"scc_reference_matrix.csv"`).
#' @param detector_noise_file Optional CSV path to detector-specific WLS noise
#'   floors, as written by [unmix_controls()] (`"scc_detector_noise.csv"`). If
#'   omitted, `unmix_samples()` first looks beside `unmixing_matrix_file`, then
#'   falls back to the built-in scalar noise floor.
#' @param control_file Optional control mapping CSV used to assign unique
#'   fluorophore channel names (`$PnN`) and marker annotations (`$PnS`) in output
#'   FCS files. Marker annotations may be duplicated. An explicit value takes
#'   precedence. When omitted, the exact mapping carried by an in-memory matrix
#'   returned from [unmix_controls()] or the `fcs_mapping_used.csv` saved beside
#'   `unmixing_matrix_file` is used.
#' @param unmixing_method Unmixing method (`"WLS"`, `"RWLS"`, `"OLS"`,
#'   `"NNLS"`, or `"AutoSpectral"`). `AutoSpectral` uses
#'   per-event AF assignment with marker + selected-AF OLS, plus SCC-derived
#'   spectral-variant optimization when available.
#' @param rwls_max_iter Positive integer; number of robust reweighting
#'   iterations used when `unmixing_method = "RWLS"`. The default, 1, preserves the
#'   historical behavior.
#' @param n_threads Positive integer; number of threads for event-wise
#'   AutoSpectral AF assignment and NNLS/WLS/RWLS fitting. OLS uses
#'   vectorized matrix operations. The default, 1, keeps explicit event loops
#'   single-threaded.
#' @param spectral_variant_library Optional in-memory AutoSpectral
#'   spectral-variant library. Used only when `unmixing_method = "AutoSpectral"`.
#' @param spectral_variant_library_file Optional RDS path to an AutoSpectral
#'   spectral-variant library. When omitted, `unmix_samples()` looks for
#'   `scc_spectral_variants.rds` beside `unmixing_matrix_file`.
#' @param spectral_variant_top_k Number of best variant candidates to test per
#'   positive fluorophore for AutoSpectral.
#' @param spectral_variant_min_abundance Minimum unmixed abundance for a
#'   fluorophore to be eligible for AutoSpectral variant testing.
#' @param spectral_variant_positive_fraction Additional positivity threshold
#'   as a fraction of the event's strongest fluorophore abundance for
#'   AutoSpectral variant testing.
#' @param spectral_variant_min_improvement Minimum fractional residual
#'   improvement required before accepting a cell-specific AutoSpectral
#'   variant refit.
#' @param estimate_af Logical; if `TRUE`, estimate AF signatures directly from
#'   stained sample event-wise WLS residuals, select the best candidate model by
#'   held-out WLS residual score, and append the selected AF rows to the
#'   reference matrix before unmixing. This is intended for workflows where no
#'   unstained cell control is available. Default is `FALSE`.
#' @param output_dir Root output directory. Sample-stage artifacts are written
#'   under `output_dir/unmix_samples`; unmixed FCS files are saved in its
#'   `unmixed_fcs` subdirectory.
#' @param write_fcs Logical; if `TRUE`, write unmixed FCS files to
#'   `output_dir/unmix_samples/unmixed_fcs`.
#'   Defaults to `TRUE` so unmixed FCS files are written unless disabled explicitly.
#' @param save_report Logical; if `TRUE`, write a sample QC report and
#'   sample QC metric CSVs from the in-memory unmixing results without rerunning
#'   unmixing. Defaults to `TRUE`.
#' @param save_ai_qc Logical; when a QC report is written, also write its local
#'   paste-ready numerical interpretation prompt. Defaults to `save_report`.
#' @param ai_qc_detail Retained for development-branch call compatibility; the
#'   report prompt is automatically context-bounded.
#' @param ai_qc_privacy Retained for call compatibility. Prompts contain local
#'   report measurements, project basename, and no raw events.
#' @param ai_qc_reference Retained for call compatibility; no reference-derived
#'   labels are assigned.
#' @param report_format Report format, either `"html"` (default) or `"pdf"`.
#'   Only the selected format is written. Matching is case-insensitive.
#' @param report_per_sample Logical; if `TRUE`, PDF output writes one report
#'   per sample, while HTML output adds a sample selector to one report.
#'   Defaults to `FALSE`.
#' @param save_qc_plots Logical; if `TRUE`, save QC report plots as PNG files
#'   in `qc_plot_dir` while creating the report.
#' @param qc_plot_dir Directory for sample QC report PNG files when
#'   `save_qc_plots = TRUE`.
#' @param plot_n_events Optional integer; number of events per sample retained
#'   for automatic QC report plots and returned in-memory results. The full
#'   unsampled unmixed data will still be written to FCS files. Defaults to
#'   10000. Set `plot_n_events = NULL` to keep all events in the returned
#'   object and report input.
#' @param chunk_size Integer number of events unmixed at a time per sample.
#'   Defaults to 50000 to reduce peak memory use for large FCS files. Set
#'   `chunk_size = NULL` to process each sample in one chunk.
#' @param seed Optional integer seed for deterministic subsampling.
#' @param return_type Return format: `"list"` (default), `"flowSet"`, or
#'   `"SingleCellExperiment"`. When `"flowSet"`, detector residuals are attached
#'   as `attr(x, "spectreasy_residuals")`. When `"SingleCellExperiment"`, cell-level
#'   unmixed values are returned in assay `"unmixed"`, with detector residuals in
#'   `altExp(x, "detector_residuals")` when available.
#' @param verbose Logical; if `TRUE`, print progress messages while unmixing
#'   each sample.
#' @param project_path Project directory recorded in generated report metadata.
#' @return Either a named list with one element per sample, a `flowSet`, or a
#'   `SingleCellExperiment` depending on `return_type`. For `return_type = "list"`,
#'   the result has class `spectreasy_unmixed_results`; list elements contain
#'   `data` (unmixed abundances plus retained acquisition parameters) and
#'   `residuals` (detector residual matrix when available, otherwise `NULL`).
#'   The list can be passed directly to `qc_samples(results = ...)`
#'   or coerced with `as.data.frame()`. The return value is provided invisibly to
#'   avoid printing large result objects during interactive or Quarto execution.
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
#'
#' unmixed <- unmix_samples(toy_fs, M = M_demo, unmixing_method = "OLS", output_dir = tempdir())
#' names(unmixed)
#' @export
unmix_samples <- function(sample_dir = "samples",
                          samples_dir = NULL,
                          M = NULL,
                          unmixing_matrix_file = file.path("spectreasy_outputs", "unmix_controls", "scc_reference_matrix.csv"),
                          detector_noise_file = NULL,
                          control_file = NULL,
                          unmixing_method = "AutoSpectral",
                          rwls_max_iter = 1L,
                          n_threads = 1L,
                          spectral_variant_library = NULL,
                          spectral_variant_library_file = NULL,
                          spectral_variant_top_k = 3L,
                          spectral_variant_min_abundance = 1,
                          spectral_variant_positive_fraction = 0.02,
                          spectral_variant_min_improvement = 0.01,
                          estimate_af = FALSE,
                          output_dir = "spectreasy_outputs",
                          write_fcs = TRUE,
                          save_report = TRUE,
                          save_ai_qc = save_report,
                          ai_qc_detail = "standard",
                          ai_qc_privacy = "standard",
                          ai_qc_reference = "auto",
                          report_format = "html",
                          report_per_sample = FALSE,
                          save_qc_plots = FALSE,
                          qc_plot_dir = NULL,
                          plot_n_events = 10000L,
                          chunk_size = 50000L,
                          seed = NULL,
                          return_type = c("list", "flowSet", "SingleCellExperiment"),
                          verbose = TRUE,
                          project_path = getwd()) {
    unmixing_matrix_file_missing <- missing(unmixing_matrix_file)
    return_type <- .match_arg_ci(
        return_type,
        c("list", "flowSet", "SingleCellExperiment"),
        "return_type"
    )
    report_format <- .match_arg_ci(report_format, c("html", "pdf"), "report_format")
    write_fcs <- .normalize_scalar_logical(write_fcs, "write_fcs")
    save_report <- .normalize_scalar_logical(save_report, "save_report")
    save_ai_qc <- .normalize_scalar_logical(save_ai_qc, "save_ai_qc")
    ai_qc_detail <- .match_arg_ci(ai_qc_detail, c("compact", "standard", "full"), "ai_qc_detail")
    ai_qc_privacy <- .match_arg_ci(ai_qc_privacy, c("standard", "strict", "none"), "ai_qc_privacy")
    report_per_sample <- .normalize_scalar_logical(report_per_sample, "report_per_sample")
    save_qc_plots <- .normalize_scalar_logical(save_qc_plots, "save_qc_plots")
    verbose <- .normalize_scalar_logical(verbose, "verbose")
    .with_optional_seed(seed)
    scc_dir <- NULL
    if (!is.null(samples_dir)) {
        if (!identical(sample_dir, "samples") && !identical(sample_dir, samples_dir)) {
            stop("Use only one of sample_dir or samples_dir.", call. = FALSE)
        }
        sample_dir <- samples_dir
    }
    sample_entries <- .prepare_unmix_samples_input(sample_dir)
    if (!is.null(control_file)) {
        if (!is.character(control_file) || length(control_file) != 1L ||
            is.na(control_file) || !nzchar(trimws(control_file))) {
            stop("control_file must be NULL or a single non-empty CSV path.", call. = FALSE)
        }
        control_file <- .resolve_control_file_path(control_file)
        if (!file.exists(control_file)) {
            .spectreasy_stop_missing_file(control_file, label = "control_file")
        }
    }
    if (!is.null(unmixing_matrix_file) && !unmixing_matrix_file_missing && !file.exists(unmixing_matrix_file)) {
        .spectreasy_stop_missing_file(unmixing_matrix_file, label = "unmixing_matrix_file")
    }
    unmixing_matrix_file <- .resolve_unmixing_matrix_file_for_samples(
        unmixing_matrix_file = unmixing_matrix_file,
        output_dir = output_dir
    )

    M <- .load_unmix_samples_matrix(
        M = M,
        unmixing_matrix_file = unmixing_matrix_file,
        detector_noise_file = detector_noise_file,
        scc_dir = scc_dir
    )

    output_control_context <- .resolve_unmix_output_control_context(
        M = M,
        control_file = control_file,
        unmixing_matrix_file = unmixing_matrix_file
    )

    method <- .normalize_unmix_method(unmixing_method)
    rwls_max_iter <- .normalize_rwls_max_iter(rwls_max_iter)
    n_threads <- .normalize_n_threads(n_threads)
    estimate_af <- .normalize_scalar_logical(estimate_af, "estimate_af")
    chunk_size <- .normalize_unmix_chunk_size(chunk_size)
    plot_n_events <- .normalize_unmix_plot_n_events(plot_n_events)
    .validate_unmix_samples_output(output_dir, write_fcs, save_report)
    output_paths <- .unmix_samples_output_paths(output_dir)
    unmixed_output_dir <- if (isTRUE(write_fcs)) {
        .next_safe_output_dir(output_paths$unmixed_dir)
    } else {
        output_paths$unmixed_dir
    }

    .print_unmix_samples_header(
        verbose = verbose,
        sample_entries = sample_entries,
        method = method,
        M = M,
        chunk_size = chunk_size,
        plot_n_events = plot_n_events,
        write_fcs = write_fcs,
        unmixed_output_dir = unmixed_output_dir
    )

    if (estimate_af) {
        if (isTRUE(verbose)) {
            .spectreasy_console_step("Estimate AF", "from stained samples")
        }
        M <- .estimate_blind_af_reference(
            M = M,
            sample_entries = sample_entries,
            n_bands = 10L,
            candidate_quantile = 0.90,
            max_training_events = 20000L,
            max_evaluation_events = 5000L,
            seed = seed,
            verbose = verbose
        )
    }

    spectral_variant_library_resolved <- if (.is_autospectral_method(method)) {
        .resolve_spectral_variant_library_for_unmixing(
            M = M,
            spectral_variant_library = spectral_variant_library,
            spectral_variant_library_file = spectral_variant_library_file,
            unmixing_matrix_file = unmixing_matrix_file
        )
    } else {
        NULL
    }

    results <- list()
    fluorophore_source_all <- rownames(M)
    output_marker_map <- .resolve_output_marker_map(
        fluorophore_source_all,
        control_file = output_control_context$control_file,
        control_df = output_control_context$control_df
    )

    .ensure_unmix_output_directory(unmixed_output_dir, write_fcs)
    processing_context <- list(
        method = method,
        rwls_max_iter = rwls_max_iter,
        n_threads = n_threads,
        spectral_variant_library = spectral_variant_library_resolved,
        spectral_variant_top_k = spectral_variant_top_k,
        spectral_variant_min_abundance = spectral_variant_min_abundance,
        spectral_variant_positive_fraction = spectral_variant_positive_fraction,
        spectral_variant_min_improvement = spectral_variant_min_improvement,
        chunk_size = chunk_size,
        plot_n_events = plot_n_events,
        write_fcs = write_fcs,
        output_marker_map = output_marker_map,
        unmixed_output_dir = unmixed_output_dir,
        verbose = verbose
    )
    results <- list()
    for (entry in sample_entries) {
        results[[entry$sample_name]] <- .process_unmix_sample_entry(entry, M, processing_context)
    }
    results <- .initialize_unmixed_results(
        results = results,
        method = method,
        M = M,
        spectral_variant_library = spectral_variant_library_resolved,
        write_fcs = write_fcs,
        unmixed_output_dir = unmixed_output_dir
    )

    if (isTRUE(save_report)) {
        results <- .save_unmix_samples_qc(
            results = results,
            M = M,
            output_paths = output_paths,
            report_format = report_format,
            report_per_sample = report_per_sample,
            save_qc_plots = save_qc_plots,
            qc_plot_dir = qc_plot_dir,
            method = method,
            unmixing_matrix_file = unmixing_matrix_file,
            unmixing_matrix_file_missing = unmixing_matrix_file_missing,
            detector_noise_file = detector_noise_file,
            spectral_variant_library_file = spectral_variant_library_file,
            unmixed_output_dir = unmixed_output_dir,
            run_settings = list(
                rwls_max_iter = rwls_max_iter,
                n_threads = n_threads,
                spectral_variant_top_k = spectral_variant_top_k,
                spectral_variant_min_abundance = spectral_variant_min_abundance,
                spectral_variant_positive_fraction = spectral_variant_positive_fraction,
                spectral_variant_min_improvement = spectral_variant_min_improvement,
                estimate_af = estimate_af,
                write_fcs = write_fcs,
                save_qc_plots = save_qc_plots,
                report_per_sample = report_per_sample,
                plot_n_events = plot_n_events,
                chunk_size = chunk_size,
                save_ai_qc = save_ai_qc,
                seed = seed
            ),
            project_path = project_path
        )
    }

    if (isTRUE(save_ai_qc)) {
        report_file <- attr(results, "qc_report_file")
        report_data <- attr(results, "qc_report_data")
        if (!is.null(report_data) && !is.null(report_file) &&
            length(report_file) == 1L && file.exists(report_file)) {
            prompt_export <- .export_report_ai_qc(
                report_data,
                report_file,
                numeric_paths = report_data$qc_metric_paths %||% character()
            )
            attr(results, "ai_qc_prompt_path") <- prompt_export$prompt
            attr(results, "ai_qc_data_paths") <- prompt_export$numeric_sources
            attr(results, "ai_qc_paths") <- c(prompt = prompt_export$prompt)
        }
    }

    if (identical(return_type, "flowSet")) {
        return(invisible(.unmixed_results_to_flowset(results)))
    }
    if (identical(return_type, "SingleCellExperiment")) {
        return(invisible(.unmixed_results_to_sce(results, sample_entries = sample_entries)))
    }

    invisible(results)
}
