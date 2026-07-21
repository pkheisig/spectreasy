.ai_qc_control_metrics <- function(M = NULL, control_report_data = NULL, controls = NULL) {
    report <- control_report_data %||% list()
    if (is.null(M)) M <- report$matrix %||% controls$M %||% NULL
    summary <- report$qc_summary %||% if (!is.null(M)) attr(M, "qc_summary") else NULL
    summary <- as.data.frame(summary %||% data.frame(), stringsAsFactors = FALSE)
    rows <- list()
    if (!nrow(summary)) {
        return(list(
            entities = list(),
            metrics = list(.ai_qc_metric(
                "QC-CTRL-AVAILABILITY", "controls", "controls", missing_reason =
                    "No in-memory control QC summary or compatible saved control metrics were supplied."
            )),
            source = "unavailable"
        ))
    }
    numeric_col <- function(name, i) {
        if (!name %in% names(summary)) return(NULL)
        suppressWarnings(as.numeric(summary[[name]][i]))
    }
    text_col <- function(name, i) {
        if (!name %in% names(summary)) return(NULL)
        as.character(summary[[name]][i])
    }
    for (i in seq_len(nrow(summary))) {
        entity <- text_col("fluorophore", i) %||% text_col("sample", i) %||% paste0("control_", i)
        total <- numeric_col("n_total", i)
        scatter <- numeric_col("n_scatter_gated", i)
        final <- numeric_col("n_final", i)
        fixed_target <- numeric_col("actual_spectral_events", i) %||% numeric_col("spectral_event_target", i)
        retention <- if (is.finite(total %||% NA_real_) && total > 0 && is.finite(final %||% NA_real_)) 100 * final / total else NULL
        fixed_selection <- is.finite(fixed_target %||% NA_real_) && is.finite(final %||% NA_real_) && final <= fixed_target
        event_metric <- .ai_qc_metric(
            "QC-CTRL-EVENTS-FINAL", entity, "controls", final, "events",
            direction = "higher_better",
            missing_reason = "Final control event count was not recorded.",
            metadata = list(intentional_fixed_selection = isTRUE(fixed_selection), configured_target = fixed_target)
        )
        if (isTRUE(event_metric$availability)) {
            event_metric <- .ai_qc_hard_grade(
                event_metric,
                failure = final <= 0,
                review = final < 10 && !isTRUE(fixed_selection),
                explanation = if (isTRUE(fixed_selection)) {
                    "The retained count reflects an intentional fixed-size spectral selection and is not penalized for percentage retention."
                } else if (final <= 0) "No final control events were available." else if (final < 10) "Fewer than ten final control events were available." else "A non-empty reliable final control pool was available."
            )
        }
        saturated <- toupper(text_col("saturated", i) %||% "") %in% c("YES", "TRUE", "1")
        saturation_metric <- .ai_qc_metric(
            "QC-CTRL-SATURATION", entity, "controls", saturated, "logical",
            direction = "lower_better", availability = "saturated" %in% names(summary),
            missing_reason = "Saturation status was not recorded."
        )
        if (isTRUE(saturation_metric$availability)) {
            saturation_metric <- .ai_qc_hard_grade(
                saturation_metric, review = saturated,
                explanation = if (saturated) "One or more detector measurements reached the configured saturation boundary." else "No saturation was recorded."
            )
        }
        rows[[length(rows) + 1L]] <- list(
            id = entity,
            sample = text_col("sample", i),
            fluorophore = text_col("fluorophore", i),
            marker = text_col("marker", i),
            control_type = text_col("type", i),
            peak_detector = text_col("peak_channel", i),
            gate_method = text_col("intensity_gate_type", i),
            background_method = text_col("scc_background_method", i),
            metrics = list(
                total_events = .ai_qc_metric("QC-CTRL-EVENTS-TOTAL", entity, "controls", total, "events", direction = "higher_better", missing_reason = "Total event count was not recorded."),
                scatter_events = .ai_qc_metric("QC-CTRL-EVENTS-SCATTER", entity, "controls", scatter, "events", direction = "higher_better", missing_reason = "Scatter-gated event count was not recorded."),
                final_events = event_metric,
                retention = .ai_qc_metric("QC-CTRL-RETENTION", entity, "controls", retention, "percent", direction = "higher_better", missing_reason = "Retention could not be calculated.", metadata = list(intentional_fixed_selection = isTRUE(fixed_selection))),
                stain_index = .ai_qc_metric("QC-CTRL-STAIN-INDEX", entity, "controls", numeric_col("stain_index", i), "ratio", direction = "higher_better", missing_reason = "Robust stain index was not recorded."),
                saturation = saturation_metric,
                saturation_fraction = .ai_qc_metric("QC-CTRL-SATURATION-FRACTION", entity, "controls", numeric_col("saturation_fraction", i), "fraction", direction = "lower_better", missing_reason = "A numerical saturation fraction was not retained."),
                positive_negative_separation = .ai_qc_metric("QC-CTRL-SEPARATION", entity, "controls", numeric_col("positive_negative_separation", i), "robust_separation", direction = "higher_better", missing_reason = "Positive/negative separation was not retained as a separate machine-readable value."),
                spectral_stability = .ai_qc_metric("QC-CTRL-SPECTRAL-STABILITY", entity, "controls", numeric_col("spectral_stability", i), "robust_similarity", direction = "higher_better", missing_reason = "Replicate spectral stability was not available."),
                background_subtraction = .ai_qc_metric(
                    "QC-CTRL-BACKGROUND-SUBTRACTION", entity, "controls",
                    if (!is.null(text_col("scc_background_method", i))) text_col("scc_background_method", i) else NULL,
                    "method", direction = "descriptive",
                    missing_reason = "Background-subtraction provenance was not recorded."
                )
            )
        )
    }
    list(
        entities = rows, metrics = list(), source = "control_qc_summary",
        post_unmixing_pairwise = .ai_qc_pairwise_metrics(controls$unmixed_list %||% NULL)
    )
}
