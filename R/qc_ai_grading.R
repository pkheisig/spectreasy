.ai_qc_hard_grade <- function(metric, failure = FALSE, review = FALSE, explanation = NULL) {
    if (!.ai_qc_is_metric(metric)) return(metric)
    grade <- if (isTRUE(failure)) "poor" else if (isTRUE(review)) "review" else "good"
    metric$grade <- grade
    metric$grade_provenance <- .ai_qc_grade_provenance(
        basis = "hard_validation",
        direction = metric$direction,
        explanation = explanation %||% switch(
            grade,
            good = "Hard structural validation passed.",
            review = "Hard structural validation requires review.",
            poor = "Hard structural validation failed."
        )
    )
    metric
}

.ai_qc_threshold_grade <- function(
    metric, review_threshold, poor_threshold,
    basis = c("literature_rule", "user_reference"),
    assumptions_met = TRUE, event_reliable = TRUE,
    rule_name = NULL, rule_version = .ai_qc_metric_version,
    explanation = NULL
) {
    basis <- .match_arg_ci(basis, c("literature_rule", "user_reference"), "basis")
    if (!.ai_qc_is_metric(metric) || !isTRUE(metric$availability)) return(metric)
    if (!isTRUE(assumptions_met) || !isTRUE(event_reliable)) {
        metric$grade <- "not_graded"
        metric$grade_provenance <- .ai_qc_grade_provenance(
            basis = "not_graded", direction = metric$direction,
            threshold = c(review = review_threshold, poor = poor_threshold),
            profile_name = rule_name, profile_version = rule_version,
            explanation = "The rule was not applied because its assumptions or event-reliability requirement was not satisfied."
        )
        return(metric)
    }
    value <- suppressWarnings(as.numeric(metric$value)[1])
    if (!is.finite(value)) return(metric)
    adverse <- if (identical(metric$direction, "higher_better")) -value else value
    adverse_review <- if (identical(metric$direction, "higher_better")) -review_threshold else review_threshold
    adverse_poor <- if (identical(metric$direction, "higher_better")) -poor_threshold else poor_threshold
    metric$grade <- if (adverse >= adverse_poor) "poor" else if (adverse >= adverse_review) "review" else "good"
    metric$grade_provenance <- .ai_qc_grade_provenance(
        basis = basis,
        threshold = c(review = review_threshold, poor = poor_threshold),
        profile_name = rule_name, profile_version = rule_version,
        direction = metric$direction,
        explanation = explanation %||% paste("Applied", gsub("_", " ", basis), "with inclusive review and poor boundaries.")
    )
    metric
}

.ai_qc_profile_grade <- function(metric, reference, profile, minimum_n = 5L) {
    if (!.ai_qc_is_metric(metric) || !isTRUE(metric$availability) || is.null(reference)) return(metric)
    n <- as.integer(reference$n %||% 0L)
    if (n < minimum_n) {
        metric$grade <- "not_graded"
        metric$grade_provenance <- .ai_qc_grade_provenance(
            basis = "not_graded", profile_name = profile$name,
            profile_version = profile$version, reference_n = n,
            matching_strata = profile$strata,
            direction = metric$direction,
            interval = reference$interval %||% NULL,
            explanation = "Fewer than five matched reviewed clean runs; the range is descriptive only."
        )
        return(metric)
    }
    value <- suppressWarnings(as.numeric(metric$value)[1])
    med <- suppressWarnings(as.numeric(reference$median)[1])
    mad <- suppressWarnings(as.numeric(reference$mad)[1])
    if (!all(is.finite(c(value, med, mad))) || mad <= 0) return(metric)
    z <- (value - med) / (1.4826 * mad)
    adverse <- if (identical(metric$direction, "higher_better")) -z else z
    metric$grade <- if (adverse >= 4) "poor" else if (adverse >= 2.5) "review" else "good"
    metric$grade_provenance <- .ai_qc_grade_provenance(
        basis = "empirical_reference", threshold = c(review = 2.5, poor = 4),
        interval = reference$interval %||% NULL, profile_name = profile$name,
        profile_version = profile$version, reference_n = n,
        matching_strata = profile$strata, direction = metric$direction,
        explanation = sprintf("Context-matched robust z-score %.3f against reviewed clean runs.", z)
    )
    metric$metadata$robust_z <- z
    metric
}

.ai_qc_apply_profile <- function(x, profile = NULL) {
    if (is.null(profile)) return(x)
    apply_one <- function(value) {
        if (.ai_qc_is_metric(value)) {
            ref <- profile$metrics[[value$metric_id]] %||% NULL
            return(.ai_qc_profile_grade(value, ref, profile))
        }
        if (is.list(value) && !is.data.frame(value)) value <- lapply(value, apply_one)
        value
    }
    apply_one(x)
}

.ai_qc_apply_within_run <- function(metrics) {
    if (length(metrics) < 5L) return(metrics)
    vals <- vapply(metrics, function(m) if (.ai_qc_is_metric(m) && isTRUE(m$availability)) as.numeric(m$value)[1] else NA_real_, numeric(1))
    finite <- is.finite(vals)
    if (sum(finite) < 5L) return(metrics)
    center <- stats::median(vals[finite])
    spread <- stats::mad(vals[finite])
    if (!is.finite(spread) || spread <= 0) return(metrics)
    z <- abs((vals - center) / (1.4826 * spread))
    for (i in which(finite & z > 3.5)) {
        if (identical(metrics[[i]]$grade, "not_graded")) {
            metrics[[i]]$grade <- "review"
            metrics[[i]]$grade_provenance <- .ai_qc_grade_provenance(
                basis = "relative_within_run", direction = metrics[[i]]$direction,
                threshold = 3.5, reference_n = sum(finite),
                explanation = sprintf("Within-run robust outlier (absolute z-score %.3f); this cannot establish Poor without calibrated or hard evidence.", z[i])
            )
        }
    }
    metrics
}

.ai_qc_grade_summary <- function(x) {
    metrics <- .ai_qc_walk_metrics(x)
    grades <- factor(vapply(metrics, `[[`, character(1), "grade"), levels = .ai_qc_grades)
    counts <- as.list(stats::setNames(as.integer(table(grades)), .ai_qc_grades))
    list(total_metrics = length(metrics), counts = counts)
}
