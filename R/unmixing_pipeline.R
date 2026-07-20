.load_unmix_samples_matrix <- function(M,
                                       unmixing_matrix_file,
                                       detector_noise_file,
                                       scc_dir = NULL) {
    if (!is.null(M)) {
        M <- .as_reference_matrix(M, "M")
        return(.load_detector_noise_for_unmixing(
            M,
            detector_noise_file = detector_noise_file,
            scc_dir = scc_dir
        ))
    }
    if (is.null(unmixing_matrix_file) || !file.exists(unmixing_matrix_file)) {
        stop(
            "No reference matrix provided. Run unmix_controls() first, then supply M or a valid unmixing_matrix_file.",
            call. = FALSE
        )
    }
    .stop_if_static_unmixing_matrix_path(unmixing_matrix_file, arg_name = "unmixing_matrix_file")
    M <- .as_reference_matrix(.read_unmixing_matrix_csv(unmixing_matrix_file), "M")
    .load_detector_noise_for_unmixing(
        M,
        detector_noise_file = detector_noise_file,
        unmixing_matrix_file = unmixing_matrix_file,
        scc_dir = scc_dir
    )
}

.validate_unmix_samples_output <- function(output_dir, write_fcs, save_report) {
    if (!isTRUE(write_fcs) && !isTRUE(save_report)) return(invisible(output_dir))
    if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) || !nzchar(trimws(output_dir))) {
        stop("output_dir must be a non-empty directory path when writing FCS files or a report.", call. = FALSE)
    }
    if (file.exists(output_dir) && !dir.exists(output_dir)) {
        stop("output_dir points to a file, not a directory: ", output_dir, call. = FALSE)
    }
    invisible(output_dir)
}

.print_unmix_samples_header <- function(verbose,
                                        sample_entries,
                                        method,
                                        M,
                                        chunk_size,
                                        plot_n_events,
                                        write_fcs,
                                        unmixed_output_dir) {
    if (!isTRUE(verbose)) return(invisible(NULL))
    .spectreasy_console_header("unmix_samples")
    .spectreasy_console_field("Samples", length(sample_entries))
    .spectreasy_console_field("Method", method)
    .spectreasy_console_field("AF bands", .reference_af_band_count(M))
    if (is.finite(chunk_size)) {
        .spectreasy_console_field("Chunk", paste0(format(chunk_size, big.mark = ","), " events"))
    }
    if (!is.null(plot_n_events)) {
        .spectreasy_console_field("Plot data", paste0(format(plot_n_events, big.mark = ","), " events/sample for plots"))
    }
    if (isTRUE(write_fcs)) {
        .spectreasy_console_field("Output", .spectreasy_console_path(unmixed_output_dir))
    }
    invisible(NULL)
}

.save_unmix_samples_qc <- function(results,
                                   M,
                                   output_paths,
                                   report_format,
                                   report_per_sample,
                                   save_qc_plots,
                                   qc_plot_dir,
                                   method,
                                   unmixing_matrix_file,
                                   unmixing_matrix_file_missing,
                                   detector_noise_file,
                                   spectral_variant_library_file,
                                   unmixed_output_dir,
                                   run_settings,
                                   project_path) {
    qc_samples_dir <- .next_safe_output_dir(output_paths$qc_samples_dir)
    output_file <- file.path(qc_samples_dir, paste0("qc_samples_report.", report_format))
    .spectreasy_console_field("Report", .spectreasy_console_path(output_file))
    report <- qc_samples(
        results = results,
        M = M,
        output_file = output_file,
        unmixing_method = method,
        qc_metrics_dir = qc_samples_dir,
        qc_plot_dir = qc_plot_dir,
        save_qc_pngs = save_qc_plots,
        report_format = report_format,
        report_per_sample = report_per_sample,
        overwrite = "overwrite",
        report_artifact_paths = list(
            matrix = if (!isTRUE(unmixing_matrix_file_missing)) unmixing_matrix_file else NULL,
            detector_noise = detector_noise_file,
            spectral_variant_library = spectral_variant_library_file,
            output_dir = unmixed_output_dir
        ),
        report_run_settings = run_settings,
        project_path = project_path
    )
    attr(results, "qc_report_file") <- report$output_file
    attr(results, "qc_samples_dir") <- qc_samples_dir
    attr(results, "qc_metrics_dir") <- qc_samples_dir
    attr(results, "qc_plot_dir") <- report$qc_plot_dir
    attr(results, "qc_report_data") <- report$report_data
    missing_files <- report$output_file[!file.exists(report$output_file)]
    if (length(missing_files)) {
        stop(
            "Automatic sample QC report was requested but was not created at: ",
            paste(missing_files, collapse = ", "),
            call. = FALSE
        )
    }
    results
}

.ensure_unmix_output_directory <- function(path, write_fcs) {
    if (!isTRUE(write_fcs)) return(invisible(path))
    if (!dir.exists(path)) {
        dir.create(path, showWarnings = FALSE, recursive = TRUE)
        if (!dir.exists(path)) {
            Sys.sleep(0.5)
            dir.create(path, showWarnings = FALSE, recursive = TRUE)
        }
    }
    if (!dir.exists(path)) {
        stop("Could not create unmixed FCS output directory: ", path, call. = FALSE)
    }
    invisible(path)
}

.unmix_chunk_calc_args <- function(flow_frame, M, context, sample_name) {
    args <- list(
        flow_frame = flow_frame,
        M = M,
        method = context$method,
        file_name = sample_name,
        rwls_max_iter = context$rwls_max_iter,
        n_threads = context$n_threads,
        spectral_variant_library = context$spectral_variant_library,
        spectral_variant_top_k = context$spectral_variant_top_k,
        spectral_variant_min_abundance = context$spectral_variant_min_abundance,
        spectral_variant_positive_fraction = context$spectral_variant_positive_fraction,
        spectral_variant_min_improvement = context$spectral_variant_min_improvement,
        return_residuals = TRUE
    )
    if (identical(context$method, "Spectreasy")) {
        args$spectreasy_weight_quantile <- context$spectreasy_weight_quantile
    }
    args
}

.unmix_chunk_write_data <- function(result, fluorophore_source, method) {
    data <- .aggregate_af_columns_for_output(result$data, fluorophore_source = fluorophore_source)
    fluorophores <- fluorophore_source[!grepl("^AF_", fluorophore_source, ignore.case = TRUE)]
    fluorophores <- intersect(colnames(data), fluorophores)
    autospectral <- if (.is_autospectral_style_method(method) && "AF Index" %in% colnames(data)) "AF Index" else character()
    passthrough <- .get_passthrough_parameter_names(colnames(data))
    list(
        matrix = as.matrix(data[, unique(c(fluorophores, autospectral, passthrough)), drop = FALSE]),
        fluorophores = fluorophores
    )
}

.store_unmix_write_chunk <- function(storage, chunk, event_idx, n_events, sample_name) {
    if (is.null(storage)) {
        storage <- matrix(
            NA_real_, nrow = n_events, ncol = ncol(chunk),
            dimnames = list(NULL, colnames(chunk))
        )
    } else if (!identical(colnames(chunk), colnames(storage))) {
        stop("Internal error: unmixed FCS columns changed between chunks for sample '", sample_name, "'.", call. = FALSE)
    }
    storage[event_idx, ] <- chunk
    storage
}

.write_unmixed_sample_fcs <- function(write_data,
                                      source_frame,
                                      fluorophores,
                                      marker_map,
                                      output_dir,
                                      sample_name,
                                      method,
                                      M) {
    storage.mode(write_data) <- "numeric"
    output_frame <- flowCore::flowFrame(write_data)
    output_frame <- .apply_output_fcs_feature_labels(
        target_ff = output_frame,
        source_ff = source_frame,
        fluorophore_cols = fluorophores,
        marker_label_map = marker_map
    )
    path <- file.path(output_dir, .unmixed_fcs_filename(sample_name, method, M))
    flowCore::write.FCS(output_frame, path)
    invisible(path)
}

.process_unmix_sample_entry <- function(entry, M, context) {
    sample_name <- entry$sample_name
    if (isTRUE(context$verbose)) .spectreasy_console_field("Sample", sample_name)
    source_frame <- if (inherits(entry$flow_frame, "flowFrame")) {
        entry$flow_frame
    } else {
        .spectreasy_read_fcs(entry$file_path, label = "sample FCS file")
    }
    n_events <- nrow(flowCore::exprs(source_frame))
    if (n_events == 0L) {
        stop("Sample '", sample_name, "' contains no events to unmix.", call. = FALSE)
    }

    chunk_indices <- .unmix_chunk_indices(n_events, chunk_size = context$chunk_size)
    keep_global <- .unmix_sample_keep_indices(n_events, plot_n_events = context$plot_n_events)
    fluorophore_source <- rownames(M)
    write_data <- NULL
    fluorophores_to_keep <- character()
    data_chunks <- list()
    residual_chunks <- list()
    variant_infos <- list()
    decoder_weights <- NULL

    for (chunk_i in seq_along(chunk_indices)) {
        event_idx <- chunk_indices[[chunk_i]]
        chunk_frame <- if (length(chunk_indices) == 1L) source_frame else source_frame[event_idx, ]
        result <- do.call(
            calc_residuals,
            .unmix_chunk_calc_args(chunk_frame, M, context, sample_name)
        )
        if (isTRUE(context$write_fcs)) {
            output <- .unmix_chunk_write_data(result, fluorophore_source, context$method)
            fluorophores_to_keep <- output$fluorophores
            write_data <- .store_unmix_write_chunk(write_data, output$matrix, event_idx, n_events, sample_name)
        }
        local_keep <- which(event_idx %in% keep_global)
        if (length(local_keep) > 0L) {
            data_chunks[[length(data_chunks) + 1L]] <- result$data[local_keep, , drop = FALSE]
            if (!is.null(result$residuals)) {
                residual_chunks[[length(residual_chunks) + 1L]] <- result$residuals[local_keep, , drop = FALSE]
            }
        }
        if (!is.null(result$spectral_variant_info)) {
            variant_infos[[length(variant_infos) + 1L]] <- result$spectral_variant_info
        }
        if (is.null(decoder_weights) && !is.null(result$spectreasy_decoder_weights)) {
            decoder_weights <- result$spectreasy_decoder_weights
        }
        if (chunk_i %% 5L == 0L) gc(verbose = FALSE)
    }

    if (isTRUE(context$write_fcs)) {
        .write_unmixed_sample_fcs(
            write_data = write_data,
            source_frame = source_frame,
            fluorophores = fluorophores_to_keep,
            marker_map = context$output_marker_map,
            output_dir = context$unmixed_output_dir,
            sample_name = sample_name,
            method = context$method,
            M = M
        )
    }
    result <- .combine_unmix_result_chunks(
        data_chunks = data_chunks,
        residual_chunks = residual_chunks,
        method = context$method,
        M = M,
        variant_infos = variant_infos,
        spectreasy_decoder_weights = decoder_weights
    )
    if (!is.null(entry$cell_ids)) {
        attr(result, "source_cell_ids") <- entry$cell_ids[keep_global]
    }
    gc(verbose = FALSE)
    result
}

.initialize_unmixed_results <- function(results,
                                        method,
                                        M,
                                        spectral_variant_library,
                                        write_fcs,
                                        unmixed_output_dir) {
    class(results) <- c("spectreasy_unmixed_results", "list")
    attr(results, "method") <- method
    attr(results, "reference_matrix") <- M
    attr(results, "spectral_variant_library") <- spectral_variant_library
    attr(results, "blind_af_info") <- attr(M, "blind_af_info")
    attr(results, "qc_report_file") <- NULL
    attr(results, "qc_samples_dir") <- NULL
    attr(results, "qc_metrics_dir") <- NULL
    attr(results, "qc_plot_dir") <- NULL
    attr(results, "qc_report_data") <- NULL
    attr(results, "unmixed_fcs_dir") <- if (isTRUE(write_fcs)) unmixed_output_dir else NULL
    results
}
