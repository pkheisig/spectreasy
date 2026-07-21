.report_relative_artifact_path <- function(path, from_file) {
    target <- strsplit(normalizePath(path, mustWork = FALSE), "/", fixed = TRUE)[[1]]
    base <- strsplit(normalizePath(dirname(from_file), mustWork = FALSE), "/", fixed = TRUE)[[1]]
    common <- 0L
    limit <- min(length(target), length(base))
    while (common < limit && identical(target[common + 1L], base[common + 1L])) common <- common + 1L
    parts <- c(rep("..", length(base) - common), target[(common + 1L):length(target)])
    paste(parts[nzchar(parts)], collapse = "/")
}

.report_ai_qc_html <- function(report_data, output_file) {
    qc <- report_data$ai_qc %||% NULL
    if (!inherits(qc, "spectreasy_ai_qc")) return("")
    counts <- qc$grade_summary$counts
    cards <- paste(vapply(.ai_qc_grades, function(grade) paste0(
        "<div class=\"ai-qc-grade ai-qc-grade-", grade, "\"><small>",
        .report_html_escape(gsub("_", " ", grade)), "</small><strong>",
        counts[[grade]] %||% 0L, "</strong></div>"
    ), character(1)), collapse = "")
    findings <- qc$overall_summary$top_findings %||% list()
    finding_html <- if (length(findings)) paste0("<ol>", paste(vapply(findings, function(item) paste0(
        "<li><strong>", .report_html_escape(item$grade), " &middot; ", .report_html_escape(item$metric_id),
        "</strong> &mdash; ", .report_html_escape(item$entity), ". ", .report_html_escape(item$explanation), "</li>"
    ), character(1)), collapse = ""), "</ol>") else "<p>No deterministic top findings were available.</p>"
    paths <- report_data$ai_qc_artifact_paths %||% character()
    paths <- paths[file.exists(paths)]
    link_html <- if (length(paths)) paste0("<div class=\"ai-qc-links\">", paste(vapply(names(paths), function(kind) paste0(
        "<a href=\"", .report_html_escape(.report_relative_artifact_path(paths[[kind]], output_file)), "\" download>",
        .report_html_escape(toupper(kind)), "</a>"
    ), character(1)), collapse = ""), "</div>") else "<p class=\"note\">AI-QC artifacts can be generated independently of this visual report.</p>"
    prompt <- build_ai_qc_prompt(qc, detail = "standard")
    profile <- qc$quality_reference$profile
    profile_text <- if (is.null(profile)) "No compatible reviewed profile" else paste(profile$name, profile$version)
    paste0(
        "<section id=\"ai-ready-qc\" class=\"ai-qc-section\"><h2><button type=\"button\">AI-ready QC</button></h2><div class=\"section-body\">",
        "<p class=\"note\">", .report_html_escape(qc$privacy$notice), "</p>",
        "<div class=\"ai-qc-meta\"><span>Status <strong>", .report_html_escape(qc$overall_summary$status), "</strong></span><span>Privacy <strong>", .report_html_escape(qc$privacy$mode), "</strong></span><span>Profile <strong>", .report_html_escape(profile_text), "</strong></span></div>",
        "<div class=\"ai-qc-grades\">", cards, "</div><h3>Deterministic top findings</h3>", finding_html,
        "<details><summary>How grades were assigned</summary><p>Hard validation is applied first; compatible literature or reviewed empirical profiles follow. Within-run outliers can reach Review but not Poor without calibrated or hard evidence. Uncalibrated metrics remain Not graded.</p></details>",
        link_html,
        "<details class=\"ai-qc-prompt\"><summary>Paste-ready evaluation prompt</summary><button class=\"ai-qc-copy\" type=\"button\">Copy prompt</button><pre>", .report_html_escape(prompt), "</pre></details>",
        "</div></section>"
    )
}

#' Render a Spectreasy QC report object as HTML
#' @param report_data Object returned by a Spectreasy report data collector.
#' @param output_file Target `.html` file.
#' @param overwrite Collision policy: `"version"`, `"overwrite"`, or `"error"`.
#' @return Structured report metadata including the final output path, source
#'   fingerprint, self-contained status, and cached report data.
#' @export
render_qc_html_report <- function(report_data, output_file, overwrite=c("version","overwrite","error")) {
    if(!inherits(report_data,"spectreasy_report_data")) stop("report_data must be created by a Spectreasy report data collector.",call.=FALSE)
    spec <- .report_output_spec(output_file,"html")
    include_nxn <- nrow(.report_nxn_entries(report_data)) > 0L
    output_file <- .report_resolve_bundle(spec$path, overwrite, report_data = report_data)
    dir.create(dirname(output_file),recursive=TRUE,showWarnings=FALSE)
    nxn_companion_files <- if (include_nxn) .report_nxn_companion_paths(output_file, report_data) else character()
    if (length(nxn_companion_files)) report_data$nxn_companion_files <- nxn_companion_files
    report_data$artifacts <- rbind(
        report_data$artifacts,
        data.frame(
            path = normalizePath(output_file, mustWork = FALSE),
            type = "HTML",
            role = "Rendered QC report",
            exists = TRUE,
            modified = format(report_data$created_at, "%Y-%m-%d %H:%M:%S %Z"),
            stringsAsFactors = FALSE
        )
    )
    template <- system.file("report_templates","qc_report.html",package="spectreasy")
    if(!nzchar(template)) template <- file.path("inst","report_templates","qc_report.html")
    if(!file.exists(template)) stop("Spectreasy HTML report template was not found.",call.=FALSE)
    html <- paste(readLines(template,warn=FALSE,encoding="UTF-8"),collapse="\n")
    sample_reports <- report_data$sample_reports %||% list()
    sections <- if (length(sample_reports)) .report_sections(sample_reports[[1]], nxn_companion_files = nxn_companion_files) else .report_sections(report_data, nxn_companion_files = nxn_companion_files)
    toc_titles <- vapply(sections,`[`,character(1),1)
    toc_labels <- vapply(seq_along(sections),function(i) .report_toc_label(names(sections)[i],toc_titles[i]),character(1))
    toc_prefix <- if (length(sample_reports)) "sample-1-" else ""
    toc <- paste0(if (inherits(report_data$ai_qc, "spectreasy_ai_qc")) "<a data-section=\"ai-ready-qc\" href=\"#ai-ready-qc\">AI-ready QC</a>" else "", paste0("<a data-section=\"", names(sections), "\" href=\"#",toc_prefix,names(sections),"\">",.report_html_escape(toc_labels),"</a>",collapse=""))
    selector <- ""
    if (length(sample_reports)) {
        sample_names <- names(sample_reports)
        if (is.null(sample_names) || any(!nzchar(sample_names))) sample_names <- paste("Sample", seq_along(sample_reports))
        selector <- paste0(
            "<div class=\"report-sample-selector\"><label for=\"report-sample\">Sample</label><select id=\"report-sample\">",
            paste0("<option value=\"", seq_along(sample_reports), "\">", .report_html_escape(sample_names), "</option>", collapse = ""),
            "</select></div>"
        )
        body <- paste(vapply(seq_along(sample_reports), function(sample_index) {
            sample_sections <- .report_sections(sample_reports[[sample_index]], nxn_companion_files = nxn_companion_files)
            sample_body <- paste(vapply(seq_along(sample_sections), function(i) {
                section_html <- .report_section(names(sample_sections)[i], sample_sections[[i]][1], sample_sections[[i]][2])
                sub(paste0("id=\"", names(sample_sections)[i], "\""), paste0("id=\"sample-", sample_index, "-", names(sample_sections)[i], "\""), section_html, fixed = TRUE)
            }, character(1)), collapse = "")
            paste0("<div class=\"sample-report-panel", if (sample_index == 1L) " is-active" else "", "\" data-sample-panel=\"", sample_index, "\">", sample_body, "</div>")
        }, character(1)), collapse = "")
    } else {
        body <- paste(vapply(seq_along(sections),function(i) .report_section(names(sections)[i],sections[[i]][1],sections[[i]][2]),character(1)),collapse="")
    }
    body <- paste0(.report_ai_qc_html(report_data, output_file), body)
    counts <- report_data$counts
    summary <- c(list(`Created`=format(report_data$created_at,"%Y-%m-%d %H:%M:%S %Z"),`Spectreasy`=report_data$version,`Method`=report_data$unmixing_method,`Cytometer`=report_data$cytometer %||% "Not recorded"),counts,.report_changed_run_settings(report_data))
    summary_html <- paste(vapply(names(summary),function(nm) paste0("<div><small>",.report_html_escape(nm),"</small><strong>",.report_html_escape(.report_scalar(summary[[nm]])),"</strong></div>"),character(1)),collapse="")
    header <- paste0("<header><span class=\"kicker\">Spectreasy &middot; ",.report_html_escape(report_data$report_type),"</span><h1>",.report_html_escape(report_data$report_type)," report</h1><p class=\"paths\">",.report_html_escape(report_data$project_path),"</p><div class=\"summary-grid\">",summary_html,"</div></header>")
    report_class <- if (inherits(report_data, "spectreasy_control_report_data")) "report-control" else "report-sample"
    replacements <- list("{{REPORT_TITLE}}"=paste("Spectreasy",report_data$report_type,"report"),"{{REPORT_CLASS}}"=report_class,"{{REPORT_TOC}}"=toc,"{{REPORT_HEADER}}"=header,"{{REPORT_SELECTOR}}"=selector,"{{REPORT_BODY}}"=body)
    for(key in names(replacements)) html <- gsub(key,replacements[[key]],html,fixed=TRUE)
    writeLines(html,output_file,useBytes=TRUE)
    companion_files <- character()
    if (length(nxn_companion_files)) companion_files <- .report_render_nxn_companions(report_data, nxn_companion_files)
    if (inherits(report_data, "spectreasy_control_report_data") && length(companion_files)) {
        report_data$plot_manifest$nxn <- companion_files
        report_data$nxn_companion_files <- companion_files
        report_data$artifacts <- rbind(
            report_data$artifacts,
            .report_artifacts(companion_files, roles = "High-resolution control NxN matrix")
        )
    }
    source_paths <- report_data$source_fingerprint$path %||% character()
    metadata <- list(output_file=normalizePath(output_file,mustWork=TRUE),companion_files=companion_files,assets_dir=NULL,self_contained=!include_nxn,format="html",report_type=report_data$report_type,created_at=report_data$created_at,source_fingerprint=report_data$source_fingerprint,stale=.report_is_stale(output_file,source_paths),report_data=report_data)
    class(metadata) <- c("spectreasy_report_result","list")
    invisible(metadata)
}
