.ai_qc_schema_name <- "spectreasy-ai-qc"
.ai_qc_schema_version <- "1.0.0"
.ai_qc_metric_version <- "1.0.0"
.ai_qc_grades <- c("good", "review", "poor", "not_graded")
.ai_qc_grade_bases <- c(
    "hard_validation", "literature_rule", "empirical_reference",
    "relative_within_run", "user_reference", "not_graded"
)
.ai_qc_top_level <- c(
    "schema", "provenance", "privacy", "scope", "experiment",
    "quality_reference", "overall_summary", "structural_validation",
    "controls", "reference_matrix", "af", "samples", "detectors",
    "pairwise_findings", "grade_summary", "limitations",
    "metric_definitions", "references"
)

.ai_qc_scalar <- function(x, default = NULL) {
    if (is.null(x) || !length(x) || is.na(x[[1]])) default else x[[1]]
}

.ai_qc_available <- function(value) {
    if (is.null(value) || !length(value)) return(FALSE)
    if (is.numeric(value)) return(any(is.finite(value)))
    any(!is.na(value))
}

.ai_qc_grade_provenance <- function(
    basis = "not_graded", threshold = NULL, interval = NULL,
    profile_name = NULL, profile_version = NULL, reference_n = 0L,
    matching_strata = list(), direction = "descriptive",
    explanation = "No applicable calibrated grading rule."
) {
    basis <- .match_arg_ci(basis, .ai_qc_grade_bases, "basis")
    list(
        basis = basis,
        threshold = threshold,
        interval = interval,
        profile_name = profile_name,
        profile_version = profile_version,
        reference_n = as.integer(reference_n %||% 0L),
        matching_strata = matching_strata %||% list(),
        direction = as.character(direction)[1],
        explanation = as.character(explanation)[1]
    )
}

.ai_qc_metric <- function(
    id, entity, stage, value = NULL, units = NULL, definition_ref = id,
    direction = "descriptive", availability = NULL, missing_reason = NULL,
    grade = "not_graded", grade_provenance = NULL, metadata = list()
) {
    grade <- .match_arg_ci(grade, .ai_qc_grades, "grade")
    if (is.null(availability)) availability <- .ai_qc_available(value)
    if (!isTRUE(availability)) value <- NULL
    if (isTRUE(availability)) missing_reason <- NULL
    if (is.null(grade_provenance)) {
        grade_provenance <- .ai_qc_grade_provenance(
            direction = direction,
            explanation = if (isTRUE(availability)) {
                "Metric is descriptive because no compatible grading rule or reviewed reference cohort was available."
            } else {
                missing_reason %||% "Metric was not available from the supplied artifacts."
            }
        )
    }
    list(
        metric_id = as.character(id)[1],
        entity = as.character(entity)[1],
        stage = as.character(stage)[1],
        value = value,
        units = units,
        definition_ref = as.character(definition_ref)[1],
        direction = as.character(direction)[1],
        availability = isTRUE(availability),
        missing_reason = missing_reason,
        grade = grade,
        grade_provenance = grade_provenance,
        metadata = metadata %||% list()
    )
}

.ai_qc_is_metric <- function(x) {
    is.list(x) && all(c("metric_id", "availability", "grade", "grade_provenance") %in% names(x))
}

.ai_qc_walk_metrics <- function(x) {
    out <- list()
    visit <- function(value) {
        if (.ai_qc_is_metric(value)) {
            out[[length(out) + 1L]] <<- value
        } else if (is.list(value) && !is.data.frame(value)) {
            lapply(value, visit)
        }
        invisible(NULL)
    }
    visit(x)
    out
}

.ai_qc_normalize_nonfinite <- function(x) {
    nonfinite <- list()
    visit <- function(value, path = "root") {
        if (is.data.frame(value)) {
            value[] <- lapply(seq_along(value), function(i) visit(value[[i]], paste0(path, ".", names(value)[i])))
            return(value)
        }
        if (is.list(value)) {
            nms <- names(value)
            for (i in seq_along(value)) {
                key <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else as.character(i)
                value[i] <- list(visit(value[[i]], paste0(path, ".", key)))
            }
            return(value)
        }
        if (is.numeric(value) && any(!is.finite(value), na.rm = TRUE)) {
            idx <- which(!is.finite(value))
            nonfinite[[length(nonfinite) + 1L]] <<- list(path = path, indices = as.integer(idx))
            value[idx] <- NA_real_
        }
        value
    }
    value <- visit(x)
    list(value = value, nonfinite = nonfinite)
}

.ai_qc_stable_order <- function(x) {
    if (is.data.frame(x)) {
        if (nrow(x) > 1L) {
            key_cols <- intersect(c("metric_id", "entity", "stage", "sample", "detector", "marker", "control"), names(x))
            if (length(key_cols)) {
                ord <- do.call(order, c(x[key_cols], list(na.last = TRUE, method = "radix")))
                x <- x[ord, , drop = FALSE]
            }
        }
        rownames(x) <- NULL
        return(x)
    }
    if (!is.list(x)) return(x)
    if (.ai_qc_is_metric(x)) return(lapply(x, .ai_qc_stable_order))
    nms <- names(x)
    if (!is.null(nms) && length(nms)) {
        preferred <- if (identical(sort(nms), sort(.ai_qc_top_level))) .ai_qc_top_level else character()
        rest <- sort(setdiff(nms, preferred), method = "radix")
        x <- x[c(preferred, rest)]
    }
    lapply(x, .ai_qc_stable_order)
}

#' Validate a Spectreasy AI-ready QC object
#'
#' @param x Object to validate.
#' @param error Whether to stop on validation failures.
#' @return Invisibly returns `TRUE`, or `FALSE` when `error = FALSE`.
#' @export
validate_ai_qc <- function(x, error = TRUE) {
    problems <- character()
    if (!inherits(x, "spectreasy_ai_qc")) problems <- c(problems, "Object does not inherit from spectreasy_ai_qc.")
    missing <- setdiff(.ai_qc_top_level, names(x))
    if (length(missing)) problems <- c(problems, paste0("Missing top-level components: ", paste(missing, collapse = ", "), "."))
    if (!identical(x$schema$name, .ai_qc_schema_name)) problems <- c(problems, "Unexpected schema name.")
    if (!identical(x$schema$version, .ai_qc_schema_version)) problems <- c(problems, "Unexpected schema version.")
    metrics <- .ai_qc_walk_metrics(x)
    bad_grade <- vapply(metrics, function(m) !m$grade %in% .ai_qc_grades, logical(1))
    bad_basis <- vapply(metrics, function(m) !m$grade_provenance$basis %in% .ai_qc_grade_bases, logical(1))
    if (any(bad_grade)) problems <- c(problems, "One or more metrics use an invalid grade.")
    if (any(bad_basis)) problems <- c(problems, "One or more metrics use an invalid grade basis.")
    if (length(problems) && isTRUE(error)) stop(paste(problems, collapse = " "), call. = FALSE)
    invisible(!length(problems))
}

.new_spectreasy_ai_qc <- function(components) {
    missing <- setdiff(.ai_qc_top_level, names(components))
    if (length(missing)) stop("Internal AI-QC schema error; missing: ", paste(missing, collapse = ", "), call. = FALSE)
    components <- components[.ai_qc_top_level]
    class(components) <- c("spectreasy_ai_qc", "list")
    validate_ai_qc(components)
    components
}
