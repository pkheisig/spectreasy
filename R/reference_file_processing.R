.process_reference_file <- function(fcs_file,
                                    control_df,
                                    sample_patterns,
                                    metadata,
                                    config,
                                    af_data_raw = NULL,
                                    universal_negatives = NULL,
                                    bead_negative = NULL,
                                    scc_background = NULL) {
    sn_ext <- basename(fcs_file)
    sn <- tools::file_path_sans_ext(sn_ext)

    is_extra_af <- FALSE

    row_info <- .get_control_rows_for_reference(control_df, c(sn_ext, sn))
    sample_info <- .resolve_reference_sample_type(
        filename = sn,
        row_info = row_info,
        patterns = sample_patterns,
        default = config$default_sample_type
    )

    fluor_name <- if (nrow(row_info) > 0 && !is.na(row_info$fluorophore[1])) row_info$fluorophore[1] else sample_info$pattern
    marker_name <- if (nrow(row_info) > 0 && "marker" %in% colnames(row_info) && !is.na(row_info$marker[1])) {
        trimws(as.character(row_info$marker[1]))
    } else {
        ""
    }

    if (is_extra_af) {
        return(NULL)
    }
    if (.is_af_control_row(fluorophore = fluor_name, marker = marker_name, filename = sn_ext)) {
        return(NULL)
    }

    ff <- .spectreasy_read_fcs(fcs_file, label = "SCC FCS file")
    pd <- flowCore::pData(flowCore::parameters(ff))
    raw_data <- flowCore::exprs(ff)

    if (!.validate_reference_raw_data(raw_data, sn, detector_names = metadata$detector_names)) {
        return(NULL)
    }

    scatter_info <- .apply_reference_manual_scatter_gates(
        raw_data = raw_data,
        sample_type = sample_info$type,
        filename = sn_ext,
        manual_gates = config$manual_gates
    )
    if (is.null(scatter_info)) {
        scatter_info <- .compute_reference_scatter_gate(
            raw_data = raw_data,
            pd = pd,
            sample_type = sample_info$type,
            outlier_percentile = config$outlier_percentile,
            debris_percentile = config$debris_percentile,
            subsample_n = config$subsample_n,
            max_clusters = config$max_clusters,
            min_cluster_proportion = config$min_cluster_proportion,
            gate_contour_beads = config$gate_contour_beads,
            gate_contour_cells = config$gate_contour_cells,
            bead_gate_scale = config$bead_gate_scale
        )
    }
    if (is.null(scatter_info)) {
        return(NULL)
    }

    peak_info <- .select_reference_peak_channel(
        gated_data = scatter_info$gated_data,
        detector_names = metadata$detector_names,
        row_info = row_info,
        channel_alias_map = metadata$channel_alias_map,
        sn_ext = sn_ext,
        sn = sn,
        cytometer = config$cytometer
    )
    peak_channel <- peak_info$peak_channel
    .spectreasy_console_field("SCC", paste0(fluor_name, " (", sn, ") -> ", peak_channel))

    peak_vals <- scatter_info$gated_data[, peak_channel]
    context <- list(
        sn = sn,
        sn_ext = sn_ext,
        fluor_name = fluor_name,
        marker_name = marker_name,
        sample_type = sample_info$type,
        raw_data = raw_data,
        pd = pd,
        scatter = scatter_info,
        peak_channel = peak_channel,
        peak_vals = peak_vals,
        metadata = metadata,
        config = config
    )
    if (isTRUE(.get_reference_config_value(config, "spectral_scc_pipeline", FALSE))) {
        return(.process_reference_autospectral(
            context = context,
            row_info = row_info,
            af_data_raw = af_data_raw,
            universal_negatives = universal_negatives,
            bead_negative = bead_negative,
            scc_background = scc_background
        ))
    }
    .process_reference_standard(
        context = context,
        row_info = row_info,
        af_data_raw = af_data_raw,
        universal_negatives = universal_negatives,
        bead_negative = bead_negative
    )

}

# Combines individual sample spectra and AF signatures into a single spillover matrix.
# Extracts spectra, cleans up row/column names, structures the metadata attributes
# (such as QC summary and parameter info), and performs basic sanity checks.
# Returns the finalized reference matrix.
.finalize_reference_matrix <- function(results_list,
                                       qc_summary_list,
                                       af_signatures_norm,
                                       af_data_raw,
                                       af_bank_info,
                                       detector_names,
                                       pd_meta,
                                       save_qc_plots = FALSE,
                                       out_path = NULL) {
    results_dt <- if (length(results_list) > 0) data.table::rbindlist(results_list) else data.table::data.table()
    spectra_list <- if (nrow(results_dt) > 0) results_dt$spectrum else list()
    if (nrow(results_dt) > 0) names(spectra_list) <- results_dt$fluorophore
    scc_positive_events <- lapply(spectra_list, function(x) attr(x, "scc_positive_events"))
    has_scc_positive_events <- vapply(
        scc_positive_events,
        function(x) is.matrix(x) && nrow(x) > 0L && ncol(x) == length(detector_names),
        logical(1)
    )

    if (is.null(af_signatures_norm) && !is.null(af_data_raw)) {
        af_vec <- pmax(af_data_raw, 0)
        af_max <- max(af_vec, na.rm = TRUE)
        if (is.finite(af_max) && af_max > 0) {
            af_signatures_norm <- matrix(af_vec / af_max, nrow = 1)
            rownames(af_signatures_norm) <- "AF"
            colnames(af_signatures_norm) <- detector_names
        }
    }

    if (!is.null(af_signatures_norm) && nrow(af_signatures_norm) > 0) {
        for (i in seq_len(nrow(af_signatures_norm))) {
            nm <- rownames(af_signatures_norm)[i]
            if (is.na(nm) || !nzchar(nm)) nm <- if (i == 1) "AF" else paste0("AF_", i)
            base_nm <- nm
            k <- 2L
            while (nm %in% names(spectra_list)) {
                nm <- paste0(base_nm, "_", k)
                k <- k + 1L
            }
            spectra_list[[nm]] <- as.numeric(af_signatures_norm[i, ])
        }
    }

    if (length(spectra_list) == 0) {
        warning("No valid spectra found!")
        return(NULL)
    }

    M <- do.call(rbind, spectra_list)
    colnames(M) <- detector_names

    if (length(qc_summary_list) > 0) {
        attr(M, "qc_summary") <- data.table::rbindlist(qc_summary_list)
    }
    if (isTRUE(save_qc_plots)) {
        attr(M, "qc_plot_dir") <- out_path
    }
    if (!is.null(af_bank_info)) {
        attr(M, "af_bank_info") <- af_bank_info
    }
    if (length(scc_positive_events) > 0L && any(has_scc_positive_events)) {
        attr(M, "scc_positive_events") <- scc_positive_events[has_scc_positive_events]
    }
    attr(M, "detector_pd") <- pd_meta
    M
}
