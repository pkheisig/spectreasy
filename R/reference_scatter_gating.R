.compute_reference_scatter_gate <- function(raw_data,
                                            pd,
                                            sample_type,
                                            outlier_percentile,
                                            debris_percentile,
                                            subsample_n,
                                            max_clusters,
                                            min_cluster_proportion,
                                            gate_contour_beads,
                                            gate_contour_cells,
                                            bead_gate_scale) {
    scatter <- .resolve_reference_scatter_channels(raw_data)
    if (is.null(scatter)) {
        return(NULL)
    }
    fsc <- scatter$fsc
    ssc <- scatter$ssc

    data_raw_scatter <- raw_data[, c(fsc, ssc)]
    fsc_max <- quantile(data_raw_scatter[, 1], 1 - outlier_percentile, na.rm = TRUE)
    ssc_max <- quantile(data_raw_scatter[, 2], 1 - outlier_percentile, na.rm = TRUE)
    fsc_lower_limit <- debris_percentile * fsc_max
    fsc_upper_limit <- .reference_cell_fsc_upper_fraction * fsc_max
    ssc_upper_limit <- .reference_cell_ssc_upper_fraction * ssc_max
    valid_idx <- which(data_raw_scatter[, 1] < fsc_max & data_raw_scatter[, 2] < ssc_max & data_raw_scatter[, 1] > 0 & data_raw_scatter[, 2] > 0)
    data_filtered <- data_raw_scatter[valid_idx, ]
    if (sample_type %in% c("cells", "unstained")) {
        debris_threshold <- fsc_lower_limit
        data_filtered <- data_filtered[
            data_filtered[, 1] >= fsc_lower_limit &
                data_filtered[, 1] <= fsc_upper_limit &
                data_filtered[, 2] <= ssc_upper_limit,
        ]
    } else {
        debris_threshold <- 0
        fsc_upper_limit <- fsc_max
        ssc_upper_limit <- ssc_max
    }
    if (nrow(data_filtered) < 100) return(NULL)

    if (!is.null(subsample_n) && nrow(data_filtered) > subsample_n) {
        data_fit <- data_filtered[sample(nrow(data_filtered), subsample_n), ]
    } else {
        data_fit <- data_filtered
    }

    gmm_result <- .fit_reference_gmm_populations(data_fit, max_k = max_clusters, min_prop = min_cluster_proportion)
    if (is.null(gmm_result)) return(NULL)

    if (sample_type == "beads") {
        selected_pops <- .select_reference_bead_population(gmm_result)$selected
        gate_level <- gate_contour_beads
    } else {
        selected_pops <- .select_reference_cell_populations(
            gmm_result,
            fsc_min = debris_threshold,
            fsc_max = fsc_upper_limit,
            ssc_max = ssc_upper_limit,
            ratio_max = if (sample_type == "unstained") .reference_unstained_ssc_fsc_ratio_max else Inf
        )$selected
        gate_level <- gate_contour_cells
        if (length(selected_pops) == 0) selected_pops <- .select_reference_bead_population(gmm_result)$selected
    }
    if (length(selected_pops) == 0) return(NULL)

    final_gate <- .create_reference_merged_gate(
        gmm_result,
        selected_pops,
        gate_level,
        scale = if (sample_type == "beads") bead_gate_scale else 1.0,
        clip_x = if (sample_type %in% c("cells", "unstained")) fsc_upper_limit else Inf,
        clip_y = if (sample_type %in% c("cells", "unstained")) ssc_upper_limit else Inf
    )
    if (is.null(final_gate) || nrow(final_gate) == 0) return(NULL)
    if (sample_type %in% c("cells", "unstained")) {
        final_gate$x <- pmax(fsc_lower_limit, final_gate$x)
        if (nrow(final_gate) > 0 && (final_gate$x[1] != final_gate$x[nrow(final_gate)] || final_gate$y[1] != final_gate$y[nrow(final_gate)])) {
            final_gate <- rbind(final_gate, final_gate[1, , drop = FALSE])
        }
    }

    inside_gate <- sp::point.in.polygon(raw_data[, fsc], raw_data[, ssc], final_gate$x, final_gate$y) > 0
    if (sample_type %in% c("cells", "unstained")) {
        inside_gate <- inside_gate &
            raw_data[, fsc] >= fsc_lower_limit &
            raw_data[, fsc] <= fsc_upper_limit &
            raw_data[, ssc] <= ssc_upper_limit
    }
    gated_data <- raw_data[inside_gate, ]
    if (nrow(gated_data) < 100) return(NULL)

    list(
        gated_data = gated_data,
        final_gate = final_gate,
        fsc = fsc,
        ssc = ssc,
        fsc_max = fsc_max,
        ssc_max = ssc_max
    )
}

# Identifies the peak/primary channel of high signal for single-color controls.
# For unstained controls, returns the channel with the highest median. For stained controls,
# infers the channel via the 99.9% quantile across all detectors, and cross-references/validates
# it against the target channel in metadata.
# Returns a list containing the resolved peak channel and 99.9% quantiles for all channels.
.select_reference_peak_channel <- function(gated_data, detector_names, row_info, channel_alias_map, sn_ext, sn, cytometer = "auto") {
    is_unstained <- grepl("unstained|autofluorescence|\\bAF\\b", paste(sn_ext, sn), ignore.case = TRUE)
    if (nrow(row_info) > 0) {
        is_unstained <- is_unstained || .is_af_control_row(
            fluorophore = if ("fluorophore" %in% colnames(row_info)) row_info$fluorophore[1] else "",
            marker = if ("marker" %in% colnames(row_info)) row_info$marker[1] else "",
            filename = sn_ext
        )
    }

    q999_by_channel <- apply(
        gated_data[, detector_names, drop = FALSE],
        2,
        function(x) stats::quantile(x, 0.999, na.rm = TRUE)
    )
    if (is_unstained) {
        med_by_channel <- apply(gated_data[, detector_names, drop = FALSE], 2, stats::median, na.rm = TRUE)
        return(list(peak_channel = names(which.max(med_by_channel)), q999_by_channel = q999_by_channel))
    }

    inferred_peak_channel <- detector_names[which.max(q999_by_channel)]
    peak_channel <- inferred_peak_channel

    expected_channel <- ""
    cytometer_id <- .resolve_cytometer_id(cytometer, allow_auto = TRUE, unknown_as_auto = TRUE)
    if (!identical(cytometer_id, "auto") && nrow(row_info) > 0 && "fluorophore" %in% colnames(row_info)) {
        expected_ref <- tryCatch(.load_control_file_shipped_reference(cytometer_id), error = function(e) NULL)
        if (!is.null(expected_ref)) {
            expected_raw <- .control_file_expected_channel_for_fluor(row_info$fluorophore[1], expected_ref$fluor_peak_channel_map)
            expected_channel <- .resolve_reference_control_channel(expected_raw, detector_names, channel_alias_map = channel_alias_map)
        }
    }
    if (nzchar(expected_channel)) {
        peak_channel <- expected_channel
        return(list(peak_channel = peak_channel, q999_by_channel = q999_by_channel))
    }

    if (nrow(row_info) > 0 && !is.na(row_info$channel[1]) && row_info$channel[1] != "") {
        resolved_channel <- .resolve_reference_control_channel(row_info$channel[1], detector_names, channel_alias_map = channel_alias_map)
        if (nzchar(resolved_channel)) {
            peak_channel <- resolved_channel
        } else {
            warning("Control channel '", row_info$channel[1], "' for ", sn_ext, " not found in file. Falling back to inferred channel ", inferred_peak_channel, ".")
        }
    }

    list(peak_channel = peak_channel, q999_by_channel = q999_by_channel)
}

# Gates positive and negative populations on the log-transformed peak channel.
# Uses a Gaussian Mixture Model (GMM) or density peaks to robustly separate background/negative
# events from stained/positive events, allowing flexible gating directions and percentiles.
# Returns a list of indices/logic indicating which events fall in the negative and positive gates.
