.ai_qc_sections <- function(x, detail = c("standard", "compact", "full")) {
    detail <- .match_arg_ci(detail, c("compact", "standard", "full"), "detail")
    fmt <- function(value) {
        if (is.null(value) || !length(value)) return("Not available.")
        if (is.atomic(value) && length(value) <= 8L) return(paste(value, collapse = ", "))
        jsonlite::toJSON(.ai_qc_stable_order(value), auto_unbox = TRUE, null = "null", na = "null", digits = 6, pretty = detail == "full")
    }
    view <- x
    if (detail != "full") {
        pair_limit <- if (detail == "compact") 6L else 12L
        metric_limit <- if (detail == "compact") 8L else 16L
        view$reference_matrix$pairwise_similarity <- head(view$reference_matrix$pairwise_similarity %||% list(), pair_limit)
        view$af$marker_interactions <- head(view$af$marker_interactions %||% list(), pair_limit)
        view$pairwise_findings$within_stage <- head(view$pairwise_findings$within_stage %||% list(), pair_limit)
        negative_metrics <- view$samples$negative_tail_mad %||% list()
        priority <- Filter(function(metric) metric$grade %in% c("poor", "review"), negative_metrics)
        remainder <- Filter(function(metric) !metric$grade %in% c("poor", "review"), negative_metrics)
        view$samples$negative_tail_mad <- head(c(priority, remainder), metric_limit)
        view$samples$negative_tail_mad_summary <- list(
            total = length(negative_metrics), included = length(view$samples$negative_tail_mad),
            note = "The canonical JSON retains complete metric records."
        )
        view$samples$negative_bias <- head(view$samples$negative_bias %||% list(), metric_limit)
        detector_metrics <- view$detectors$metrics %||% list()
        detector_priority <- Filter(function(metric) metric$grade %in% c("poor", "review"), detector_metrics)
        detector_remainder <- Filter(function(metric) !metric$grade %in% c("poor", "review"), detector_metrics)
        view$detectors$metrics <- head(c(detector_priority, detector_remainder), if (detail == "compact") 10L else 24L)
    }
    findings <- view$overall_summary$top_findings %||% list()
    missing <- Filter(function(m) !isTRUE(m$availability), .ai_qc_walk_metrics(x))
    definitions <- x$metric_definitions
    if (detail == "standard") definitions <- lapply(definitions, function(definition) definition[intersect(c("name", "formula", "grading_method", "limitations"), names(definition))])
    if (detail == "compact") {
        findings <- head(findings, 4L)
        definitions <- list(note = "Definitions are retained in the canonical JSON artifact.")
    }
    list(
        "Local-only notice" = x$privacy$notice,
        "Scope" = fmt(view$scope),
        "Experiment" = fmt(view$experiment),
        "Quality reference" = fmt(view$quality_reference),
        "Overall summary" = fmt(view$overall_summary$status),
        "Structural validation" = fmt(view$structural_validation),
        "Controls" = fmt(view$controls),
        "Reference matrix" = fmt(view$reference_matrix),
        "Autofluorescence" = fmt(view$af),
        "Samples" = fmt(view$samples),
        "Detectors" = fmt(view$detectors),
        "Pairwise findings" = fmt(view$pairwise_findings$within_stage),
        "Cross-stage associations" = fmt(view$pairwise_findings$cross_stage),
        "Grade summary" = fmt(view$grade_summary),
        "Top findings" = fmt(findings),
        "Missing evidence" = fmt(lapply(missing, function(m) list(metric_id = m$metric_id, entity = m$entity, reason = m$missing_reason))),
        "Limitations" = paste(x$limitations, collapse = "\n"),
        "Metric definitions" = fmt(definitions),
        "Scientific references" = fmt(x$references),
        "Reproducibility" = fmt(list(schema = x$schema, provenance = x$provenance, privacy = x$privacy)),
        "Advisory interpretation" = "Any AI-assisted interpretation is advisory. Verify every claim against the canonical measurements, grading provenance, assumptions, and experimental context."
    )
}

.render_ai_qc_text <- function(x, detail = "standard", markdown = FALSE) {
    sections <- .ai_qc_sections(x, detail)
    stopifnot(length(sections) == 21L)
    header <- if (isTRUE(markdown)) "# Spectreasy AI-ready QC" else "SPECTREASY AI-READY QC"
    rendered <- unlist(Map(function(title, body) {
        heading <- if (isTRUE(markdown)) paste0("## ", title) else paste0(toupper(title), "\n", strrep("-", nchar(title)))
        c(heading, as.character(body), "")
    }, names(sections), sections), use.names = FALSE)
    paste(c(header, "", rendered), collapse = "\n")
}

.ai_qc_pdf_caption <- function(x, export_path = NULL) {
    if (!inherits(x, "spectreasy_ai_qc")) return(NULL)
    counts <- x$grade_summary$counts
    count_line <- paste(sprintf("%s: %s", gsub("_", " ", .ai_qc_grades), unlist(counts[.ai_qc_grades])), collapse = " | ")
    profile <- x$quality_reference$profile
    profile_line <- if (is.null(profile)) "profile: none compatible" else paste0("profile: ", profile$name, " ", profile$version)
    export_line <- if (!is.null(export_path)) paste0("export: ", export_path) else "export: use export_ai_qc() or save_ai_qc = TRUE"
    paste(
        "AI-ready QC - local deterministic summary; no model, provider, key, upload, or raw events",
        paste0("status: ", x$overall_summary$status, " | privacy: ", x$privacy$mode, " | scope: ", x$scope$value, " | ", profile_line),
        count_line,
        paste0(export_line, " | AI interpretation is advisory."),
        sep = "\n"
    )
}
