.compute_reference_autospectral_scc <- function(clean_data,
                                                detector_names,
                                                peak_channel,
                                                sample_type,
                                                af_data_raw = NULL,
                                                scc_background = NULL,
                                                n_candidates = 1000L,
                                                n_spectral = 200L,
                                                min_events = 10L,
                                                scc_background_k = 2L) {
    if (!all(c(detector_names, peak_channel) %in% colnames(clean_data))) {
        return(NULL)
    }
    event_mat <- as.matrix(clean_data[, detector_names, drop = FALSE])
    complete <- stats::complete.cases(event_mat)
    peak_vals <- clean_data[, peak_channel]
    valid <- complete & is.finite(peak_vals) & peak_vals > 0 & rowSums(event_mat, na.rm = TRUE) > 0
    if (sum(valid, na.rm = TRUE) < min_events) {
        return(NULL)
    }

    n_candidates <- suppressWarnings(as.integer(n_candidates[1]))
    if (!is.finite(n_candidates) || is.na(n_candidates) || n_candidates < 1L) {
        n_candidates <- 1000L
    }
    n_spectral <- suppressWarnings(as.integer(n_spectral[1]))
    if (!is.finite(n_spectral) || is.na(n_spectral) || n_spectral < 1L) {
        n_spectral <- 200L
    }
    scc_background_k <- suppressWarnings(as.integer(scc_background_k[1]))
    if (!is.finite(scc_background_k) || is.na(scc_background_k) || scc_background_k < 1L) {
        scc_background_k <- 2L
    }

    vals_log <- log10(pmax(peak_vals, 1))
    attr(vals_log, "positive_gate_present") <- TRUE
    attr(vals_log, "negative_gate_present") <- FALSE
    attr(vals_log, "mapped_peak_channel") <- peak_channel

    af_summary <- .reference_external_negative_summary(
        detector_names = detector_names,
        af_data_raw = af_data_raw,
        scc_background = scc_background
    )
    valid_idx <- which(valid)
    candidate_n <- min(n_candidates, length(valid_idx))
    candidate_idx <- valid_idx[order(peak_vals[valid_idx], decreasing = TRUE)[seq_len(candidate_n)]]

    if (!is.null(af_summary)) {
        af_mean <- af_summary$mean[detector_names]
        af_median <- af_summary$median[detector_names]
        af_norm <- sqrt(sum(af_mean^2, na.rm = TRUE))
        if (!is.finite(af_norm) || af_norm <= 0) {
            return(NULL)
        }
        event_complete <- pmax(event_mat[complete, , drop = FALSE], 0)
        af_unit <- af_mean / (af_norm + 1e-9)
        projection <- as.numeric(event_complete %*% af_unit)
        residual <- pmax(event_complete - projection %o% af_unit, 0)
        colnames(residual) <- detector_names
        empirical_peak <- names(which.max(colMeans(residual, na.rm = TRUE)))

        af_cosine_candidates <- .reference_cosine_to_vector(
            pmax(event_mat[candidate_idx, , drop = FALSE], 0),
            af_median
        )
        cosine_ok <- is.finite(af_cosine_candidates)
        if (sum(cosine_ok) < min_events) {
            return(NULL)
        }
        candidate_idx <- candidate_idx[cosine_ok]
        af_cosine_candidates <- af_cosine_candidates[cosine_ok]
        selected_n <- min(n_spectral, length(candidate_idx))
        selected_local <- order(af_cosine_candidates, decreasing = FALSE)[seq_len(selected_n)]
        selected_idx <- candidate_idx[selected_local]
        if (length(selected_idx) < min_events) {
            return(NULL)
        }

        selected_events <- clean_data[selected_idx, , drop = FALSE]
        matched_background <- .scc_background_match(
            events = selected_events,
            background = scc_background,
            k = scc_background_k
        )
        if (is.null(matched_background)) {
            matched_background <- matrix(
                af_median,
                nrow = nrow(selected_events),
                ncol = length(detector_names),
                byrow = TRUE,
                dimnames = list(NULL, detector_names)
            )
            background_method <- "external_median"
        } else {
            background_method <- "scatter_knn"
        }
        corrected_events <- as.matrix(selected_events[, detector_names, drop = FALSE]) - matched_background
        colnames(corrected_events) <- detector_names
        positive_events <- pmax(corrected_events, 0)
        pos_spectrum_raw <- apply(positive_events, 2, stats::median, na.rm = TRUE)
        sig_pure <- pmax(pos_spectrum_raw, 0)
        max_val <- max(sig_pure, na.rm = TRUE)
        if (!is.finite(max_val) || max_val <= 0) {
            return(NULL)
        }
        spectrum_norm <- sig_pure / max_val

        selected_cos <- af_cosine_candidates[selected_local]
        attr(vals_log, "gate_type") <- "autospectral_external"
        attr(vals_log, "gate_method") <- paste0(
            "AutoSpectral-style external-negative selector after FSC/SSC gate: top ",
            length(candidate_idx),
            " peak-bright candidate event(s), kept ",
            length(selected_idx),
            " least-AF-like event(s)"
        )
        attr(vals_log, "af_cosine_max_selected") <- max(selected_cos, na.rm = TRUE)
        attr(vals_log, "af_cosine_median_selected") <- stats::median(selected_cos, na.rm = TRUE)
        attr(vals_log, "af_score_valid_events") <- length(valid_idx)
        attr(vals_log, "empirical_peak_channel") <- empirical_peak
        attr(vals_log, "scc_background_method") <- background_method

        positive_idx <- rep(FALSE, nrow(clean_data))
        positive_idx[selected_idx] <- TRUE
        af_cosine_full <- rep(NA_real_, nrow(clean_data))
        af_cosine_full[candidate_idx] <- af_cosine_candidates
        attr(spectrum_norm, "variance") <- apply(positive_events, 2, stats::var, na.rm = TRUE) / (max_val^2)
        attr(spectrum_norm, "scc_background") <- list(
            method = background_method,
            k = if (identical(background_method, "scatter_knn")) scc_background_k else NA_integer_,
            matched_events = nrow(positive_events)
        )
        attr(spectrum_norm, "scc_positive_events") <- positive_events

        return(list(
            spectrum = spectrum_norm,
            final_gated_data = selected_events,
            positive_events = positive_events,
            positive_idx = positive_idx,
            vals_log = vals_log,
            gate_min = min(peak_vals[selected_idx], na.rm = TRUE),
            gate_max = max(peak_vals[selected_idx], na.rm = TRUE),
            peak_vals = peak_vals,
            spectral_gate_info = list(
                af_basis = af_summary$spectra,
                af_cosine = af_cosine_full
            ),
            extraction_method = "autospectral_external",
            n_candidates = length(candidate_idx),
            n_selected = length(selected_idx)
        ))
    }

    top_n <- min(100L, length(candidate_idx))
    selected_idx <- candidate_idx[seq_len(top_n)]
    bottom_n <- max(10L, floor(0.10 * length(valid_idx)))
    bottom_n <- min(bottom_n, length(valid_idx) - length(selected_idx))
    if (bottom_n < min_events || length(selected_idx) < min_events) {
        return(NULL)
    }
    negative_pool <- setdiff(valid_idx[order(peak_vals[valid_idx], decreasing = FALSE)], selected_idx)
    negative_idx <- head(negative_pool, bottom_n)
    neg_spectrum <- apply(event_mat[negative_idx, , drop = FALSE], 2, mean, na.rm = TRUE)
    corrected_events <- sweep(as.matrix(clean_data[selected_idx, detector_names, drop = FALSE]), 2, neg_spectrum, "-")
    colnames(corrected_events) <- detector_names
    positive_events <- pmax(corrected_events, 0)
    sig_pure <- pmax(colMeans(positive_events, na.rm = TRUE), 0)
    max_val <- max(sig_pure, na.rm = TRUE)
    if (!is.finite(max_val) || max_val <= 0) {
        return(NULL)
    }
    spectrum_norm <- sig_pure / max_val
    positive_idx <- rep(FALSE, nrow(clean_data))
    positive_idx[selected_idx] <- TRUE
    negative_bool <- rep(FALSE, nrow(clean_data))
    negative_bool[negative_idx] <- TRUE
    attr(vals_log, "gate_type") <- "autospectral_internal"
    attr(vals_log, "gate_method") <- paste0(
        "AutoSpectral-style internal-negative fallback after FSC/SSC gate: top ",
        length(selected_idx),
        " peak-bright event(s) minus bottom ",
        length(negative_idx),
        " event(s)"
    )
    attr(vals_log, "negative_idx") <- negative_bool
    attr(vals_log, "scc_background_method") <- "internal_negative"
    attr(spectrum_norm, "variance") <- apply(positive_events, 2, stats::var, na.rm = TRUE) / (max_val^2)
    attr(spectrum_norm, "scc_background") <- list(
        method = "internal_negative",
        k = NA_integer_,
        matched_events = 0L
    )
    attr(spectrum_norm, "scc_positive_events") <- positive_events

    list(
        spectrum = spectrum_norm,
        final_gated_data = clean_data[selected_idx, , drop = FALSE],
        positive_events = positive_events,
        positive_idx = positive_idx,
        vals_log = vals_log,
        gate_min = min(peak_vals[selected_idx], na.rm = TRUE),
        gate_max = max(peak_vals[selected_idx], na.rm = TRUE),
        peak_vals = peak_vals,
        spectral_gate_info = NULL,
        extraction_method = "autospectral_internal",
        n_candidates = length(candidate_idx),
        n_selected = length(selected_idx)
    )
}

# Computes the normalized spectral signature for a stained control sample.
# Calculates the median intensity in each detector for the positive and negative gates,
# subtracts the negative/background control, and normalizes the spectrum by the peak signal.
# Returns the normalized spectrum vector.
.compute_reference_spectrum <- function(final_gated_data,
                                        gated_data,
                                        peak_vals,
                                        vals_log,
                                        detector_names,
                                        row_info,
                                        sample_type = "beads",
                                        af_data_raw = NULL,
                                        universal_negatives = NULL,
                                        bead_negative = NULL) {
    pos_spectrum_raw <- apply(final_gated_data[, detector_names, drop = FALSE], 2, median, na.rm = TRUE)
    neg_log_min <- attr(vals_log, "neg_log_min")
    neg_log_max <- attr(vals_log, "neg_log_max")
    if (isTRUE(attr(vals_log, "negative_gate_present")) &&
        is.finite(neg_log_min) && is.finite(neg_log_max) && neg_log_max > neg_log_min) {
        neg_events <- gated_data[
            peak_vals >= 10^neg_log_min & peak_vals <= 10^neg_log_max,
            detector_names,
            drop = FALSE
        ]
    } else {
        neg_events <- gated_data[
            peak_vals <= 10^quantile(vals_log, 0.15, na.rm = TRUE),
            detector_names,
            drop = FALSE
        ]
    }
    if (nrow(neg_events) < 10) {
        neg_events <- gated_data[
            peak_vals <= 10^quantile(vals_log, 0.15, na.rm = TRUE),
            detector_names,
            drop = FALSE
        ]
    }
    neg_spectrum_raw <- apply(neg_events, 2, median, na.rm = TRUE)
    uv_val <- if (nrow(row_info) > 0 && "universal.negative" %in% colnames(row_info)) {
        trimws(as.character(row_info$universal.negative[1]))
    } else {
        ""
    }
    uv_upper <- toupper(uv_val)
    uv_key <- tools::file_path_sans_ext(basename(uv_val))
    use_af_negative <- uv_upper %in% c("TRUE", "AF")
    use_named_negative <- nzchar(uv_key) &&
        !uv_upper %in% c("FALSE", "TRUE", "AF") &&
        !is.null(universal_negatives) &&
        uv_key %in% names(universal_negatives)
    final_neg <- if (use_af_negative && !is.null(af_data_raw)) {
        af_data_raw
    } else if (use_named_negative) {
        universal_negatives[[uv_key]]
    } else if (identical(sample_type, "beads") && !is.null(bead_negative)) {
        bead_negative
    } else {
        neg_spectrum_raw
    }
    sig_pure <- pmax(pos_spectrum_raw - final_neg, 0)
    max_val <- max(sig_pure, na.rm = TRUE)
    if (max_val <= 0) max_val <- max(pos_spectrum_raw, na.rm = TRUE)
    res <- sig_pure / max_val

    res
}

# Standardizes a file path to its base name without file extension for matching.
# Returns the cleaned character string/vector.
