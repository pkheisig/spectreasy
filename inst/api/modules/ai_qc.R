gui_ai_qc_output_dir <- function(root, output_root = "spectreasy_outputs") {
    gui_workflow_resolve_path(output_root, root)
}

gui_ai_qc_artifacts <- function(root, output_root = "spectreasy_outputs") {
    directory <- gui_ai_qc_output_dir(root, output_root)
    if (!dir.exists(directory)) return(character())
    files <- list.files(
        directory,
        pattern = "_ai_qc_prompt\\.txt$",
        recursive = TRUE,
        full.names = TRUE,
        ignore.case = TRUE
    )
    files[order(file.info(files)$mtime, decreasing = TRUE)]
}

gui_ai_qc_relative <- function(path, root) {
    normalized <- normalizePath(path, mustWork = FALSE)
    prefix <- paste0(normalizePath(root, mustWork = TRUE), .Platform$file.sep)
    if (!startsWith(normalized, prefix)) stop("QC prompt is outside the active project.", call. = FALSE)
    substring(normalized, nchar(prefix) + 1L)
}

gui_ai_qc_matrix_file <- function(root, output_root = "spectreasy_outputs") {
    candidates <- c(
        file.path(root, output_root, "unmix_controls", "scc_reference_matrix.csv"),
        file.path(root, "spectreasy_outputs", "unmix_controls", "scc_reference_matrix.csv")
    )
    hit <- candidates[file.exists(candidates)]
    if (length(hit)) hit[1] else NULL
}

gui_ai_qc_sources <- function(root, output_root = "spectreasy_outputs") {
    output <- gui_ai_qc_output_dir(root, output_root)
    candidates <- c(
        gui_ai_qc_matrix_file(root, output_root),
        file.path(root, "fcs_mapping.csv"),
        list.files(output, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    )
    unique(candidates[!is.na(candidates) & file.exists(candidates)])
}

gui_ai_qc_source_hashes <- function(files, root) {
    lapply(files, function(path) list(
        path = gui_ai_qc_relative(path, root),
        sha256 = spectreasy:::.ai_qc_file_hash(path)
    ))
}

gui_ai_qc_readiness <- function(root, output_root = "spectreasy_outputs") {
    prompts <- gui_ai_qc_artifacts(root, output_root)
    source_files <- gui_ai_qc_sources(root, output_root)
    latest <- if (length(prompts)) prompts[1] else NULL
    stale <- !is.null(latest) && length(source_files) && any(file.info(source_files)$mtime > file.info(latest)$mtime)
    status <- if (is.null(latest)) "not_generated" else if (stale) "stale" else "ready"
    list(
        status = status,
        available_scopes = c("control", "sample"),
        default_scope = "control",
        artifact_paths = lapply(prompts, gui_ai_qc_relative, root = root),
        source_paths = lapply(source_files, gui_ai_qc_relative, root = root),
        source_hashes = gui_ai_qc_source_hashes(source_files, root),
        stale = stale,
        generated_at = if (!is.null(latest)) format(file.info(latest)$mtime, "%Y-%m-%dT%H:%M:%S%z") else NULL,
        privacy = "local numerical report evidence",
        detail = "context-bounded"
    )
}

gui_ai_qc_preview <- function(root, output_root = "spectreasy_outputs", detail = "standard") {
    readiness <- gui_ai_qc_readiness(root, output_root)
    prompts <- gui_ai_qc_artifacts(root, output_root)
    prompt <- if (length(prompts)) paste(readLines(prompts[1], warn = FALSE, encoding = "UTF-8"), collapse = "\n") else ""
    c(readiness, list(
        scope = if (length(prompts) && grepl("qc_samples", prompts[1], ignore.case = TRUE)) "sample" else "control",
        prompt = prompt,
        prompt_characters = nchar(prompt),
        estimated_tokens = ceiling(nchar(prompt) / 4),
        measurements_only = TRUE,
        labels_assigned = FALSE,
        findings = list(),
        warnings = character(),
        missing_sections = character()
    ))
}

gui_ai_qc_profiles <- function(root) {
    list(profiles = list())
}

gui_ai_qc_generate <- function(body) {
    root <- gui_workflow_root(body)
    output_root <- gui_workflow_value(body, "output_root", "spectreasy_outputs")
    preview <- gui_ai_qc_preview(root, output_root)
    if (!nzchar(preview$prompt)) {
        stop("No saved QC prompt is available. Generate a control or sample QC report first.", call. = FALSE)
    }
    preview
}
