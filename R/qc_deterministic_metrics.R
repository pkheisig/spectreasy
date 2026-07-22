.qc_cosine_rows <- function(x, reference) {
    x <- as.matrix(x)
    reference <- as.numeric(reference)
    denom <- sqrt(rowSums(x^2, na.rm = TRUE)) * sqrt(sum(reference^2, na.rm = TRUE))
    out <- rep(NA_real_, nrow(x))
    keep <- is.finite(denom) & denom > 0
    if (any(keep)) out[keep] <- as.numeric(x[keep, , drop = FALSE] %*% reference) / denom[keep]
    out
}

.qc_entropy <- function(values) {
    values <- pmax(as.numeric(values), 0)
    total <- sum(values)
    if (!is.finite(total) || total <= 0) return(NA_real_)
    p <- values / total
    p <- p[p > 0]
    if (length(p) <= 1L) return(0)
    -sum(p * log(p)) / log(length(values))
}

.qc_metric_table <- function(values) {
    data.frame(
        metric_id = names(values),
        value = as.numeric(unlist(values, use.names = FALSE)),
        stringsAsFactors = FALSE
    )
}

#' Calculate deterministic event-level SCC quality metrics
#'
#' Uses the same negative subtraction and peak normalization convention as the
#' reference-spectrum workflow. This function is intended for already selected
#' positive and negative/background event pools; it does not choose biological
#' gates and never mutates its inputs.
#'
#' @param positive Numeric positive-event matrix (events by detectors).
#' @param negative Numeric negative/background-event matrix.
#' @param expected_peak Expected peak detector name.
#' @param detector_ranges Optional named upper acquisition ranges.
#' @param final_reference Optional final normalized reference spectrum.
#' @param off_target_references Optional matrix of compatible normalized spectra.
#' @param minimum_signal Minimum positive-minus-negative event amplitude used by
#'   spectral-consistency calculations.
#' @return A list with aggregate metrics, detector metrics, variability bands,
#'   event-level consistency values, and calculation provenance.
#' @export
calculate_control_qc_metrics <- function(
    positive, negative, expected_peak = NULL, detector_ranges = NULL,
    final_reference = NULL, off_target_references = NULL, minimum_signal = NULL
) {
    positive <- as.matrix(positive)
    negative <- as.matrix(negative)
    storage.mode(positive) <- "double"
    storage.mode(negative) <- "double"
    if (!nrow(positive) || !nrow(negative)) stop("Positive and negative event pools must be non-empty.", call. = FALSE)
    if (is.null(colnames(positive)) || is.null(colnames(negative))) stop("Event matrices require detector names.", call. = FALSE)
    if (!identical(colnames(positive), colnames(negative))) stop("Positive and negative detector order must match exactly.", call. = FALSE)
    detectors <- colnames(positive)
    background <- apply(negative, 2, stats::median, na.rm = TRUE)
    corrected <- sweep(positive, 2, background, "-")
    corrected_nonnegative <- pmax(corrected, 0)
    center <- apply(corrected_nonnegative, 2, stats::median, na.rm = TRUE)
    peak_value <- max(center, na.rm = TRUE)
    observed_peak <- if (is.finite(peak_value)) detectors[which.max(center)] else NA_character_
    reference <- if (is.null(final_reference)) {
        if (is.finite(peak_value) && peak_value > 0) center / peak_value else rep(NA_real_, length(center))
    } else {
        final_reference <- as.numeric(final_reference[detectors])
        if (any(!is.finite(final_reference))) stop("final_reference is missing one or more detector values.", call. = FALSE)
        final_reference / max(final_reference)
    }
    names(reference) <- detectors
    peak <- if (!is.null(expected_peak) && expected_peak %in% detectors) expected_peak else observed_peak
    pos_peak <- positive[, peak]
    neg_peak <- negative[, peak]
    pos_location <- stats::median(pos_peak, na.rm = TRUE)
    neg_location <- stats::median(neg_peak, na.rm = TRUE)
    raw_difference <- pos_location - neg_location
    robust_spread <- stats::mad(neg_peak, center = neg_location, constant = 1.4826, na.rm = TRUE)
    threshold <- (pos_location + neg_location) / 2
    overlap <- mean(c(pos_peak <= threshold, neg_peak > threshold), na.rm = TRUE)
    amplitude <- apply(corrected_nonnegative, 1, max, na.rm = TRUE)
    if (is.null(minimum_signal)) minimum_signal <- max(stats::median(apply(abs(sweep(negative, 2, background, "-")), 1, max), na.rm = TRUE) * 3, .Machine$double.eps)
    usable <- is.finite(amplitude) & amplitude >= minimum_signal
    event_normalized <- corrected_nonnegative[usable, , drop = FALSE]
    if (nrow(event_normalized)) event_normalized <- event_normalized / pmax(apply(event_normalized, 1, max), .Machine$double.eps)
    consistency <- if (nrow(event_normalized)) .qc_cosine_rows(event_normalized, reference) else numeric()
    q <- function(x, probability) if (length(x) && any(is.finite(x))) as.numeric(stats::quantile(x, probability, na.rm = TRUE, names = FALSE)) else NA_real_
    detector_ranges <- detector_ranges %||% stats::setNames(rep(NA_real_, length(detectors)), detectors)
    detector_ranges <- suppressWarnings(as.numeric(detector_ranges[detectors]))
    names(detector_ranges) <- detectors
    upper_hit <- sweep(positive, 2, detector_ranges, ">=")
    valid_ranges <- is.finite(detector_ranges) & detector_ranges > 0
    if (any(!valid_ranges)) {
        upper_hit[, !valid_ranges] <- NA
    }
    detector_metrics <- data.frame(
        detector = detectors,
        positive_median = apply(positive, 2, stats::median, na.rm = TRUE),
        negative_median = background,
        residual_background_median = apply(sweep(negative, 2, background, "-"), 2, stats::median, na.rm = TRUE),
        residual_background_mad = apply(sweep(negative, 2, background, "-"), 2, stats::mad, constant = 1.4826, na.rm = TRUE),
        upper_boundary_fraction = colMeans(upper_hit, na.rm = TRUE),
        stringsAsFactors = FALSE
    )
    any_fraction <- function(x) mean(rowSums(x, na.rm = TRUE) > 0, na.rm = TRUE)
    peak_index <- match(peak, detectors)
    expected_peak_upper <- if (is.finite(peak_index) && valid_ranges[[peak_index]]) mean(upper_hit[, peak_index], na.rm = TRUE) else NA_real_
    variability <- if (nrow(event_normalized)) data.frame(
        detector = detectors,
        q10 = apply(event_normalized, 2, q, probability = 0.10),
        q25 = apply(event_normalized, 2, q, probability = 0.25),
        median = apply(event_normalized, 2, stats::median, na.rm = TRUE),
        q75 = apply(event_normalized, 2, q, probability = 0.75),
        q90 = apply(event_normalized, 2, q, probability = 0.90),
        final_reference = as.numeric(reference),
        stringsAsFactors = FALSE
    ) else data.frame(detector = detectors, q10 = NA_real_, q25 = NA_real_, median = NA_real_, q75 = NA_real_, q90 = NA_real_, final_reference = as.numeric(reference))
    secondary <- list(max_projection = NA_real_, aligned_reference = NA_character_, intended_ratio = NA_real_, fraction_affected = NA_real_)
    if (!is.null(off_target_references) && nrow(event_normalized)) {
        refs <- as.matrix(off_target_references)
        if (!all(detectors %in% colnames(refs))) stop("off_target_references must use the same detector schema.", call. = FALSE)
        refs <- refs[, detectors, drop = FALSE]
        projections <- vapply(seq_len(nrow(refs)), function(i) stats::median(.qc_cosine_rows(event_normalized, refs[i, ]), na.rm = TRUE), numeric(1))
        affected <- vapply(seq_len(nrow(refs)), function(i) mean(.qc_cosine_rows(event_normalized, refs[i, ]) > 0.9, na.rm = TRUE), numeric(1))
        worst <- which.max(projections)
        intended <- stats::median(consistency, na.rm = TRUE)
        secondary <- list(
            max_projection = projections[worst],
            aligned_reference = rownames(refs)[worst] %||% as.character(worst),
            intended_ratio = projections[worst] / pmax(intended, .Machine$double.eps),
            fraction_affected = affected[worst]
        )
    }
    aggregate <- c(
        total_positive_events = nrow(positive),
        total_negative_events = nrow(negative),
        positive_peak_median = pos_location,
        negative_peak_median = neg_location,
        robust_separation = raw_difference / pmax(2 * robust_spread, .Machine$double.eps),
        raw_positive_minus_negative = raw_difference,
        overlap_misclassification_fraction = overlap,
        peak_agreement = as.numeric(!is.null(expected_peak) && identical(expected_peak, observed_peak)),
        upper_boundary_fraction = mean(upper_hit, na.rm = TRUE),
        expected_peak_upper_boundary_fraction = expected_peak_upper,
        any_detector_upper_boundary_fraction = any_fraction(upper_hit),
        spectral_cosine_median = if (length(consistency)) stats::median(consistency, na.rm = TRUE) else NA_real_,
        spectral_cosine_q10 = q(consistency, 0.10),
        spectral_events_evaluated = length(consistency),
        spectral_events_excluded_low_signal = sum(!usable)
    )
    list(
        metrics = .qc_metric_table(aggregate),
        detector_metrics = detector_metrics,
        spectrum_variability = variability,
        event_consistency = consistency,
        reference = reference,
        expected_peak = expected_peak,
        observed_peak = observed_peak,
        clipping_detectors = detectors[colSums(upper_hit, na.rm = TRUE) > 0],
        secondary_component = secondary,
        provenance = list(
            background = "negative detector-wise median",
            normalization = "event-wise maximum after non-negative background subtraction",
            minimum_usable_signal = minimum_signal,
            formulas_version = "1.0.0"
        )
    )
}

#' Calculate deterministic AF-bank diagnostics
#'
#' @param af_bank Numeric AF spectra matrix (bands by detectors).
#' @param assignments Optional event-wise band labels or one-based indices.
#'   Missing, unknown, non-integer, and out-of-range assignments are excluded
#'   from occupancy and reported explicitly in the aggregate metrics.
#' @param reconstruction_error Optional event-wise AF reconstruction errors.
#' @param requested_bands Optional requested number of AF bands.
#' @return A list containing aggregate metrics, occupancy, similarity, and
#'   detector variability.
#' @export
calculate_af_bank_qc_metrics <- function(af_bank, assignments = NULL, reconstruction_error = NULL, requested_bands = NULL) {
    af_bank <- as.matrix(af_bank)
    storage.mode(af_bank) <- "double"
    if (!nrow(af_bank) || !ncol(af_bank)) stop("af_bank must contain at least one band and detector.", call. = FALSE)
    if (is.null(rownames(af_bank))) rownames(af_bank) <- paste0("AF_", seq_len(nrow(af_bank)))
    normalized <- af_bank / pmax(sqrt(rowSums(af_bank^2)), .Machine$double.eps)
    similarity <- normalized %*% t(normalized)
    pairwise <- if (nrow(similarity) > 1L) similarity[upper.tri(similarity)] else numeric()
    sv <- svd(af_bank, nu = 0, nv = 0)$d
    effective_rank <- sum(sv > max(dim(af_bank)) * max(sv) * .Machine$double.eps)
    occupancy <- stats::setNames(rep(NA_real_, nrow(af_bank)), rownames(af_bank))
    assignment_count <- if (is.null(assignments)) 0L else length(assignments)
    valid_assignment_count <- 0L
    invalid_assignment_count <- 0L
    if (!is.null(assignments) && length(assignments)) {
        assignment_values <- trimws(as.character(assignments))
        labels <- rep(NA_character_, length(assignment_values))
        label_index <- match(assignment_values, rownames(af_bank))
        exact_label <- !is.na(label_index)
        labels[exact_label] <- rownames(af_bank)[label_index[exact_label]]

        unresolved <- !exact_label
        numeric_values <- suppressWarnings(as.numeric(assignment_values))
        valid_index <- unresolved & is.finite(numeric_values) &
            numeric_values == floor(numeric_values) &
            numeric_values >= 1 & numeric_values <= nrow(af_bank)
        labels[valid_index] <- rownames(af_bank)[as.integer(numeric_values[valid_index])]

        valid <- !is.na(labels)
        valid_assignment_count <- sum(valid)
        invalid_assignment_count <- length(labels) - valid_assignment_count
        if (valid_assignment_count > 0L) {
            tab <- table(factor(labels[valid], levels = rownames(af_bank)))
            occupancy <- as.numeric(tab) / valid_assignment_count
            names(occupancy) <- rownames(af_bank)
        }
    }
    detector_variability <- data.frame(
        detector = colnames(af_bank) %||% paste0("detector_", seq_len(ncol(af_bank))),
        median = apply(af_bank, 2, stats::median),
        mad = apply(af_bank, 2, stats::mad, constant = 1.4826),
        stringsAsFactors = FALSE
    )
    metrics <- c(
        requested_bands = requested_bands %||% NA_real_,
        returned_bands = nrow(af_bank),
        effective_rank = effective_rank,
        assignment_count = assignment_count,
        valid_assignment_count = valid_assignment_count,
        invalid_assignment_count = invalid_assignment_count,
        invalid_assignment_fraction = if (assignment_count > 0L) invalid_assignment_count / assignment_count else NA_real_,
        maximum_pairwise_similarity = if (length(pairwise)) max(pairwise) else NA_real_,
        median_pairwise_similarity = if (length(pairwise)) stats::median(pairwise) else NA_real_,
        reconstruction_error_median = if (!is.null(reconstruction_error)) stats::median(reconstruction_error, na.rm = TRUE) else NA_real_,
        reconstruction_error_q90 = if (!is.null(reconstruction_error)) as.numeric(stats::quantile(reconstruction_error, 0.9, na.rm = TRUE)) else NA_real_
    )
    list(
        metrics = .qc_metric_table(metrics),
        occupancy = data.frame(band = names(occupancy), fraction = as.numeric(occupancy), stringsAsFactors = FALSE),
        similarity = similarity,
        detector_variability = detector_variability
    )
}

.control_qc_numeric_tables <- function(M, qc_summary = NULL) {
    positive <- attr(M, "scc_qc_positive_events")
    negative <- attr(M, "scc_qc_negative_events")
    controls <- intersect(names(positive), names(negative))
    if (!length(controls)) return(list(summary = data.frame(), variability = data.frame()))
    qc_summary <- as.data.frame(qc_summary %||% attr(M, "qc_summary") %||% data.frame(), stringsAsFactors = FALSE)
    pd <- attr(M, "detector_pd")
    ranges <- NULL
    if (is.data.frame(pd) && all(c("name", "maxRange") %in% names(pd))) {
        ranges <- stats::setNames(suppressWarnings(as.numeric(pd$maxRange)), as.character(pd$name))
    }
    summary_rows <- list()
    variability_rows <- list()
    for (control in controls) {
        expected <- NULL
        if (nrow(qc_summary) && all(c("fluorophore", "peak_channel") %in% names(qc_summary))) {
            hit <- which(as.character(qc_summary$fluorophore) == control & nzchar(as.character(qc_summary$peak_channel)))
            if (length(hit)) expected <- as.character(qc_summary$peak_channel[hit[1]])
        }
        reference <- if (control %in% rownames(M)) M[control, , drop = TRUE] else NULL
        calculation <- tryCatch(calculate_control_qc_metrics(
            positive[[control]], negative[[control]], expected_peak = expected,
            detector_ranges = ranges, final_reference = reference
        ), error = function(e) NULL)
        if (is.null(calculation)) next
        values <- stats::setNames(calculation$metrics$value, calculation$metrics$metric_id)
        summary_rows[[length(summary_rows) + 1L]] <- data.frame(
            fluorophore = control,
            expected_peak = expected %||% "",
            observed_peak = calculation$observed_peak %||% "",
            positive_events = unname(values[["total_positive_events"]]),
            negative_events = unname(values[["total_negative_events"]]),
            positive_peak_median = unname(values[["positive_peak_median"]]),
            negative_peak_median = unname(values[["negative_peak_median"]]),
            robust_separation = unname(values[["robust_separation"]]),
            midpoint_overlap_fraction = unname(values[["overlap_misclassification_fraction"]]),
            spectral_cosine_median = unname(values[["spectral_cosine_median"]]),
            spectral_cosine_q10 = unname(values[["spectral_cosine_q10"]]),
            expected_peak_upper_boundary_fraction = unname(values[["expected_peak_upper_boundary_fraction"]]),
            any_detector_upper_boundary_fraction = unname(values[["any_detector_upper_boundary_fraction"]]),
            stringsAsFactors = FALSE
        )
        variability <- calculation$spectrum_variability
        variability$fluorophore <- control
        variability_rows[[length(variability_rows) + 1L]] <- variability[, c("fluorophore", "detector", "q10", "median", "q90", "final_reference"), drop = FALSE]
    }
    list(
        summary = if (length(summary_rows)) do.call(rbind, summary_rows) else data.frame(),
        variability = if (length(variability_rows)) do.call(rbind, variability_rows) else data.frame()
    )
}

.af_qc_summary_table <- function(M, info = NULL) {
    af_rows <- grepl("^AF($|_)", rownames(M), ignore.case = TRUE)
    if (!any(af_rows)) return(data.frame())
    result <- calculate_af_bank_qc_metrics(
        M[af_rows, , drop = FALSE],
        requested_bands = info$requested_bands %||% NULL
    )
    values <- stats::setNames(result$metrics$value, result$metrics$metric_id)
    data.frame(
        requested_bands = unname(values[["requested_bands"]]),
        returned_bands = unname(values[["returned_bands"]]),
        effective_rank = unname(values[["effective_rank"]]),
        maximum_pairwise_cosine_similarity = unname(values[["maximum_pairwise_similarity"]]),
        median_pairwise_cosine_similarity = unname(values[["median_pairwise_similarity"]]),
        stringsAsFactors = FALSE
    )
}

.af_band_usage_table <- function(results_df, af_band_names = NULL) {
    if (!is.data.frame(results_df) || !all(c("File", "AF Index") %in% names(results_df))) return(data.frame())
    samples <- unique(as.character(results_df$File))
    rows <- lapply(samples, function(sample) {
        values <- suppressWarnings(as.integer(results_df$`AF Index`[as.character(results_df$File) == sample]))
        values <- values[is.finite(values)]
        if (!length(values)) return(NULL)
        observed <- sort(unique(values))
        counts <- tabulate(match(values, observed), nbins = length(observed))
        labels <- if (!is.null(af_band_names) && all(observed >= 1L & observed <= length(af_band_names))) af_band_names[observed] else as.character(observed)
        data.frame(
            sample = sample, af_band_index = observed, af_band = labels,
            events = counts, fraction = counts / sum(counts), stringsAsFactors = FALSE
        )
    })
    rows <- Filter(Negate(is.null), rows)
    if (length(rows)) do.call(rbind, rows) else data.frame()
}
