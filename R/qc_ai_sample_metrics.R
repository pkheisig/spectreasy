.ai_qc_results_list <- function(results) {
    if (is.null(results)) return(NULL)
    if (is.data.frame(results)) return(NULL)
    if (is.list(results) && length(results)) results else NULL
}

.ai_qc_sample_metrics <- function(results = NULL, M = NULL, sample_report_data = NULL, method = NULL) {
    report <- sample_report_data %||% list()
    if (is.null(M)) M <- report$matrix %||% attr(results, "M") %||% NULL
    method <- .normalize_unmix_method(method %||% report$unmixing_method %||% attr(results, "method") %||% "AutoSpectral")
    results_df <- tryCatch(
        if (is.null(results)) data.frame() else .normalize_qc_report_results_df(results),
        error = function(e) data.frame()
    )
    samples <- if (nrow(results_df) && "File" %in% names(results_df)) unique(as.character(results_df$File)) else character()
    residual_list <- .ai_qc_results_list(results)
    detector_rms <- report$detector_rms %||% if (!is.null(residual_list)) {
        tryCatch(.compute_qc_report_detector_rms(residual_list, M = M, unmixing_method = method), error = function(e) NULL)
    } else NULL
    reconstruction <- report$reconstruction_error %||% if (!is.null(residual_list)) {
        tryCatch(.compute_qc_report_sample_rms(residual_list, M = M, unmixing_method = method), error = function(e) NULL)
    } else NULL
    if (!is.null(reconstruction) && "sample" %in% names(reconstruction)) samples <- unique(c(samples, as.character(reconstruction$sample)))
    nps <- report$nps %||% if (nrow(results_df) && !identical(method, "NNLS")) {
        tryCatch(calculate_nps(results_df), error = function(e) NULL)
    } else NULL
    nps <- as.data.frame(nps %||% data.frame(), stringsAsFactors = FALSE)
    entities <- lapply(samples, function(sample) {
        idx <- if (nrow(results_df) && "File" %in% names(results_df)) which(as.character(results_df$File) == sample) else integer()
        row_reconstruction <- if (!is.null(reconstruction) && "sample" %in% names(reconstruction)) reconstruction[as.character(reconstruction$sample) == sample, , drop = FALSE] else data.frame()
        af_index <- if ("AF Index" %in% names(results_df)) suppressWarnings(as.numeric(results_df$`AF Index`[idx])) else numeric()
        list(
            id = sample,
            metrics = list(
                qc_events = .ai_qc_metric("QC-SAMPLE-EVENTS", sample, "samples", length(idx), "events", direction = "higher_better"),
                reconstruction_error = .ai_qc_metric(
                    "QC-SAMPLE-RECONSTRUCTION", sample, "samples",
                    if (nrow(row_reconstruction)) row_reconstruction$median_rms_percent_of_peak[1] else NULL,
                    "percent_of_peak", direction = "lower_better",
                    missing_reason = "Residual matrices were not retained in the supplied result or saved artifacts.",
                    metadata = list(residual_profile = if (.uses_wls_residual_metric(method)) "wls_weighted_rms" else "raw_rms", unmixing_method = method)
                ),
                af_index_median = .ai_qc_metric(
                    "QC-AF-INDEX-MEDIAN", sample, "samples",
                    if (any(is.finite(af_index))) stats::median(af_index, na.rm = TRUE) else NULL,
                    "band_index", direction = "descriptive",
                    missing_reason = "AF Index was not available for this sample."
                )
            )
        )
    })
    nps_metrics <- if (nrow(nps) && all(c("File", "Marker", "NPS") %in% names(nps))) {
        lapply(seq_len(nrow(nps)), function(i) .ai_qc_metric(
            "QC-SAMPLE-NEGATIVE-TAIL-MAD", paste(nps$File[i], nps$Marker[i], sep = ":"), "samples",
            suppressWarnings(as.numeric(nps$NPS[i])), "unmixed_intensity_mad",
            direction = "lower_better",
            missing_reason = "Negative-tail MAD proxy was unavailable.",
            metadata = list(negative_population = "ungated lower 20 percent", proxy = TRUE)
        ))
    } else list()
    negative_bias_metrics <- list()
    if (nrow(results_df) && "File" %in% names(results_df) && !identical(method, "NNLS")) {
        markers <- setdiff(names(results_df), .get_result_metadata_columns(names(results_df)))
        markers <- markers[vapply(results_df[markers], is.numeric, logical(1))]
        for (sample in unique(as.character(results_df$File))) {
            sample_rows <- results_df[as.character(results_df$File) == sample, , drop = FALSE]
            for (marker in markers) {
                values <- suppressWarnings(as.numeric(sample_rows[[marker]]))
                values <- values[is.finite(values)]
                lower <- if (length(values)) values[values <= stats::quantile(values, 0.2, names = FALSE)] else numeric()
                negative_bias_metrics[[length(negative_bias_metrics) + 1L]] <- .ai_qc_metric(
                    "QC-SAMPLE-NEGATIVE-BIAS", paste(sample, marker, sep = ":"), "samples",
                    if (length(lower)) stats::median(lower) else NULL,
                    "unmixed_intensity", direction = "closer_to_zero",
                    missing_reason = "The marker-aware lower-tail bias proxy was unavailable.",
                    metadata = list(negative_population = "ungated lower 20 percent", proxy = TRUE)
                )
            }
        }
    }
    detector_metrics <- if (!is.null(detector_rms) && nrow(detector_rms)) {
        lapply(seq_len(nrow(detector_rms)), function(i) .ai_qc_metric(
            "QC-DETECTOR-RMS", as.character(detector_rms$detector[i]), "samples",
            suppressWarnings(as.numeric(detector_rms$rms_residual[i])), "residual_rms",
            direction = "lower_better",
            metadata = list(
                residual_profile = as.character(detector_rms$residual_metric[i]),
                unmixing_method = as.character(detector_rms$unmixing_method[i]),
                laser = as.character(detector_rms$laser[i])
            )
        ))
    } else list(.ai_qc_metric("QC-DETECTOR-RMS", "detectors", "samples", missing_reason = "Detector residual matrices were not retained."))
    detector_metrics <- .ai_qc_apply_within_run(detector_metrics)
    list(
        entities = entities,
        negative_tail_mad = nps_metrics,
        negative_bias = negative_bias_metrics,
        detector_metrics = detector_metrics,
        residual_profile = if (.uses_wls_residual_metric(method)) "wls_weighted_rms" else "raw_rms",
        method = method,
        source = if (nrow(results_df)) "in_memory_results" else "unavailable"
    )
}
