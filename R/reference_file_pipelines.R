.reference_stain_metrics <- function(positive_values, negative_values, raw_data, detector_names) {
    mfi_pos <- stats::median(positive_values, na.rm = TRUE)
    mfi_neg <- stats::median(negative_values, na.rm = TRUE)
    sd_neg <- stats::mad(negative_values, na.rm = TRUE)
    if (is.na(sd_neg) || sd_neg == 0) sd_neg <- stats::sd(negative_values, na.rm = TRUE)
    if (is.na(sd_neg) || sd_neg == 0) sd_neg <- 1e-6
    list(
        stain_index = (mfi_pos - mfi_neg) / (2 * sd_neg),
        saturated = any(raw_data[, detector_names, drop = FALSE] >= 260000, na.rm = TRUE)
    )
}

.save_processed_reference_qc <- function(context, final_gated_data, hist_info) {
    if (!isTRUE(context$config$save_qc_plots)) return(invisible(NULL))
    .save_reference_qc_plots_safely(
        sn = context$sn,
        raw_data = context$raw_data,
        gated_data = context$scatter$gated_data,
        final_gated_data = final_gated_data,
        pd = context$pd,
        fsc = context$scatter$fsc,
        ssc = context$scatter$ssc,
        fsc_max = context$scatter$fsc_max,
        ssc_max = context$scatter$ssc_max,
        final_gate = context$scatter$final_gate,
        vals_log = hist_info$vals_log,
        peak_vals = context$peak_vals,
        gate_min = hist_info$gate_min,
        gate_max = hist_info$gate_max,
        peak_channel = context$peak_channel,
        detector_names = context$metadata$detector_names,
        detector_labels = context$metadata$detector_labels,
        det_info = context$metadata$det_info,
        out_path = context$config$out_path,
        manual_gate_info = context$scatter$manual_gate_info,
        hist_info = hist_info,
        gui_scatter_domains = context$config$gui_scatter_domains,
        histogram_domain = .reference_gui_channel_extent(
            context$raw_data, context$peak_channel, context$sn_ext, max_points = 100000L
        ),
        source_filename = context$sn_ext,
        sample_type = context$sample_type
    )
    invisible(NULL)
}

.reference_processed_file_result <- function(context,
                                             spectrum,
                                             final_count,
                                             histogram_gate_type,
                                             stain_metrics,
                                             scc_background_method = NULL) {
    qc_fields <- list(
        sample = context$sn,
        fluorophore = context$fluor_name,
        marker = context$marker_name,
        type = context$sample_type,
        peak_channel = context$peak_channel,
        fsc_channel = context$scatter$fsc,
        ssc_channel = context$scatter$ssc,
        n_total = nrow(context$raw_data),
        n_scatter_gated = nrow(context$scatter$gated_data),
        n_final = final_count,
        scatter_gate_pct = round(100 * nrow(context$scatter$gated_data) / max(nrow(context$raw_data), 1), 1),
        histogram_gate_pct = round(100 * final_count / max(nrow(context$scatter$gated_data), 1), 1),
        intensity_gate_type = histogram_gate_type,
        stain_index = round(stain_metrics$stain_index, 1),
        saturated = ifelse(stain_metrics$saturated, "YES", "OK")
    )
    if (!is.null(scc_background_method)) qc_fields$scc_background_method <- scc_background_method
    list(
        sample_name = context$sn,
        result = data.table::data.table(
            sample = context$sn,
            fluorophore = context$fluor_name,
            type = context$sample_type,
            n_total = nrow(context$raw_data),
            n_final = final_count,
            spectrum = list(spectrum)
        ),
        qc_summary = data.table::as.data.table(qc_fields)
    )
}

.process_reference_autospectral <- function(context,
                                            row_info,
                                            af_data_raw,
                                            universal_negatives,
                                            bead_negative,
                                            scc_background) {
    positive_gate <- .resolve_reference_positive_histogram_gate(
        gated_data = context$scatter$gated_data,
        peak_channel = context$peak_channel,
        filename = context$sn_ext,
        sample_type = context$sample_type,
        manual_gates = context$config$manual_gates,
        config = context$config,
        row_info = row_info
    )
    if (is.null(positive_gate)) return(NULL)
    clean_data <- positive_gate$final_gated_data
    internal_negative <- .reference_histogram_negative_source(
        gated_data = context$scatter$gated_data,
        peak_channel = context$peak_channel,
        filename = context$sn_ext,
        sample_type = context$sample_type,
        manual_gates = context$config$manual_gates,
        detector_names = context$metadata$detector_names,
        fsc = context$scatter$fsc,
        ssc = context$scatter$ssc,
        hist_info = positive_gate$hist_info
    )
    use_background <- isTRUE(.get_reference_config_value(context$config, "scc_background_enabled", FALSE))
    negative_source <- .resolve_reference_scc_negative_source(
        row_info = row_info,
        sample_type = context$sample_type,
        af_data_raw = af_data_raw,
        scc_background = scc_background,
        universal_negatives = universal_negatives,
        bead_negative = bead_negative
    )
    selected_background <- if (use_background) negative_source$background else NULL
    selected_negative <- if (use_background) negative_source$negative else NULL
    if (use_background && is.null(selected_negative) && !is.null(internal_negative)) {
        selected_negative <- internal_negative
        selected_background <- attr(internal_negative, "scc_background", exact = TRUE)
    }
    extraction <- .compute_reference_autospectral_scc(
        clean_data = clean_data,
        detector_names = context$metadata$detector_names,
        peak_channel = context$peak_channel,
        sample_type = context$sample_type,
        af_data_raw = selected_negative,
        scc_background = selected_background,
        n_candidates = .get_reference_config_value(context$config, "autospectral_n_candidates", 1000L),
        n_spectral = .get_reference_config_value(context$config, "autospectral_n_spectral", 200L),
        min_events = .get_reference_config_value(context$config, "autospectral_min_events", 10L),
        scc_background_k = .get_reference_config_value(context$config, "scc_background_k", 2L)
    )
    if (is.null(extraction)) return(NULL)

    peak_values <- clean_data[, context$peak_channel]
    selected_peak <- peak_values[extraction$positive_idx]
    negative_count <- max(1L, floor(0.10 * length(peak_values)))
    negative_values <- peak_values[order(peak_values)[seq_len(negative_count)]]
    metrics <- .reference_stain_metrics(
        selected_peak, negative_values, context$raw_data, context$metadata$detector_names
    )
    .save_processed_reference_qc(context, extraction$final_gated_data, positive_gate$hist_info)
    background <- attr(extraction$spectrum, "scc_background")
    .reference_processed_file_result(
        context = context,
        spectrum = extraction$spectrum,
        final_count = extraction$n_selected,
        histogram_gate_type = extraction$extraction_method,
        stain_metrics = metrics,
        scc_background_method = if (!is.null(background$method)) background$method else "none"
    )
}

.process_reference_standard <- function(context,
                                        row_info,
                                        af_data_raw,
                                        universal_negatives,
                                        bead_negative) {
    positive_gate <- .resolve_reference_positive_histogram_gate(
        gated_data = context$scatter$gated_data,
        peak_channel = context$peak_channel,
        filename = context$sn_ext,
        sample_type = context$sample_type,
        manual_gates = context$config$manual_gates,
        config = context$config,
        row_info = row_info
    )
    if (is.null(positive_gate)) return(NULL)
    final_data <- positive_gate$final_gated_data
    hist_info <- positive_gate$hist_info
    negative_idx <- which(context$peak_vals <= 10^stats::quantile(hist_info$vals_log, 0.15, na.rm = TRUE))
    metrics <- if (length(negative_idx) > 0) {
        .reference_stain_metrics(
            final_data[, context$peak_channel],
            context$peak_vals[negative_idx],
            context$raw_data,
            context$metadata$detector_names
        )
    } else {
        list(
            stain_index = NA_real_,
            saturated = any(context$raw_data[, context$metadata$detector_names, drop = FALSE] >= 260000, na.rm = TRUE)
        )
    }
    spectrum <- .compute_reference_spectrum(
        final_gated_data = final_data,
        gated_data = context$scatter$gated_data,
        peak_vals = context$peak_vals,
        vals_log = hist_info$vals_log,
        detector_names = context$metadata$detector_names,
        row_info = row_info,
        sample_type = context$sample_type,
        af_data_raw = af_data_raw,
        universal_negatives = universal_negatives,
        bead_negative = bead_negative
    )
    .save_processed_reference_qc(context, final_data, hist_info)
    gate_type <- attr(hist_info$vals_log, "gate_type")
    if (is.null(gate_type) || !nzchar(gate_type)) gate_type <- "histogram"
    .reference_processed_file_result(
        context = context,
        spectrum = spectrum,
        final_count = nrow(final_data),
        histogram_gate_type = gate_type,
        stain_metrics = metrics
    )
}
