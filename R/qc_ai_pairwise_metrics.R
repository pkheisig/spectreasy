.ai_qc_pairwise_metrics <- function(results = NULL, max_pairs = 25L) {
    df <- tryCatch(if (is.null(results)) data.frame() else .normalize_qc_report_results_df(results), error = function(e) data.frame())
    if (!nrow(df) || !"File" %in% names(df)) return(list())
    markers <- setdiff(names(df), .get_result_metadata_columns(names(df)))
    markers <- markers[vapply(df[markers], is.numeric, logical(1))]
    if (length(markers) < 2L) return(list())
    pairs <- utils::combn(markers, 2L, simplify = FALSE)
    rows <- list()
    for (sample in unique(as.character(df$File))) {
        x <- df[df$File == sample, markers, drop = FALSE]
        for (pair in pairs) {
            keep <- is.finite(x[[pair[1]]]) & is.finite(x[[pair[2]]])
            if (sum(keep) < 3L) next
            covariance <- stats::cov(x[[pair[1]]][keep], x[[pair[2]]][keep])
            slope <- tryCatch(stats::coef(stats::lm(x[[pair[2]]][keep] ~ x[[pair[1]]][keep]))[[2]], error = function(e) NA_real_)
            rows[[length(rows) + 1L]] <- list(
                sample = sample, marker_a = pair[1], marker_b = pair[2],
                covariance = covariance, slope = slope,
                supporting_metric_ids = c("QC-PAIR-COVARIANCE", "QC-PAIR-SLOPE"),
                calibration_status = "descriptive_only",
                alternative_explanations = c("biological co-expression", "spillover or unmixing coupling", "sampling variation")
            )
        }
    }
    if (!length(rows)) return(list())
    strength <- vapply(rows, function(x) abs(x$covariance), numeric(1))
    rows[head(order(strength, decreasing = TRUE, na.last = NA), max_pairs)]
}

.ai_qc_cross_stage_findings <- function(controls, samples, pairwise) {
    findings <- list()
    control_metrics <- .ai_qc_walk_metrics(controls)
    sample_metrics <- .ai_qc_walk_metrics(samples)
    if (length(control_metrics) && length(sample_metrics)) {
        findings[[1]] <- list(
            finding_id = "QC-XSTAGE-001",
            statement = "Control-stage and sample-stage evidence are available for joint review.",
            supporting_metric_ids = unique(c(
                head(vapply(control_metrics, `[[`, character(1), "metric_id"), 3L),
                head(vapply(sample_metrics, `[[`, character(1), "metric_id"), 3L)
            )),
            strength = "descriptive",
            direction = "undetermined",
            calibration_status = "not_calibrated",
            alternative_explanations = c("panel complexity", "sample biology", "instrument state", "control quality"),
            causal_claim = FALSE
        )
    }
    findings
}
