gui_ai_qc_output_dir <- function(root, output_root = "spectreasy_outputs") {
    output <- gui_workflow_resolve_path(output_root, root)
    file.path(output, "ai_qc")
}

gui_ai_qc_artifacts <- function(root, output_root = "spectreasy_outputs") {
    directory <- gui_ai_qc_output_dir(root, output_root)
    if (!dir.exists(directory)) return(character())
    files <- list.files(directory, pattern = "^spectreasy_ai_qc_.*\\.(json|txt|md)$", full.names = TRUE, ignore.case = TRUE)
    files[order(file.info(files)$mtime, decreasing = TRUE)]
}

gui_ai_qc_relative <- function(path, root) {
    normalized <- normalizePath(path, mustWork = FALSE)
    prefix <- paste0(normalizePath(root, mustWork = TRUE), .Platform$file.sep)
    if (!startsWith(normalized, prefix)) stop("AI-QC artifact is outside the active project.", call. = FALSE)
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
    candidates <- c(
        gui_ai_qc_matrix_file(root, output_root),
        file.path(root, "fcs_mapping.csv"),
        file.path(root, output_root, "unmix_controls", "scc_detector_noise.csv"),
        list.files(file.path(root, output_root, "unmix_samples", "qc_samples"), pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
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
    artifacts <- gui_ai_qc_artifacts(root, output_root)
    json_files <- artifacts[grepl("(?<!manifest)\\.json$", artifacts, perl = TRUE, ignore.case = TRUE)]
    source_files <- gui_ai_qc_sources(root, output_root)
    latest <- if (length(json_files)) json_files[1] else NULL
    source_hashes <- gui_ai_qc_source_hashes(source_files, root)
    stale <- FALSE
    if (!is.null(latest) && length(source_files)) {
        manifest_path <- sub("\\.json$", "_manifest.json", latest, ignore.case = TRUE)
        manifest <- if (file.exists(manifest_path)) tryCatch(jsonlite::fromJSON(manifest_path, simplifyVector = FALSE), error = function(e) NULL) else NULL
        previous <- if (is.null(manifest$source_artifacts)) list() else manifest$source_artifacts
        previous_by_name <- stats::setNames(
            vapply(previous, function(item) as.character(if (is.null(item$sha256)) "" else item$sha256), character(1)),
            vapply(previous, function(item) basename(as.character(if (is.null(item$path)) "" else item$path)), character(1))
        )
        current_by_name <- stats::setNames(
            vapply(source_hashes, `[[`, character(1), "sha256"),
            vapply(source_hashes, function(item) basename(item$path), character(1))
        )
        comparable <- intersect(names(previous_by_name)[nzchar(previous_by_name)], names(current_by_name))
        stale <- if (length(comparable)) {
            any(previous_by_name[comparable] != current_by_name[comparable])
        } else {
            any(file.info(source_files)$mtime > file.info(latest)$mtime)
        }
    }
    has_control <- !is.null(gui_ai_qc_matrix_file(root, output_root))
    has_sample <- length(list.files(file.path(root, output_root, "unmix_samples"), pattern = "\\.(fcs|csv)$", recursive = TRUE, ignore.case = TRUE)) > 0L
    status <- if (is.null(latest)) "not_generated" else if (stale) "stale" else "ready"
    list(
        status = status,
        available_scopes = c(if (has_control) "control", if (has_sample) "sample", if (has_control && has_sample) "combined"),
        default_scope = if (has_control && has_sample) "combined" else if (has_control) "control" else if (has_sample) "sample" else "combined",
        artifact_paths = lapply(artifacts, gui_ai_qc_relative, root = root),
        source_paths = lapply(source_files, gui_ai_qc_relative, root = root),
        source_hashes = source_hashes,
        stale = stale,
        generated_at = if (!is.null(latest)) format(file.info(latest)$mtime, "%Y-%m-%dT%H:%M:%S%z") else NULL,
        privacy = c("standard", "strict", "none"), detail = c("compact", "standard", "full")
    )
}

gui_ai_qc_latest_object <- function(root, output_root = "spectreasy_outputs") {
    artifacts <- gui_ai_qc_artifacts(root, output_root)
    json_files <- artifacts[grepl("(?<!manifest)\\.json$", artifacts, perl = TRUE, ignore.case = TRUE)]
    if (!length(json_files)) return(NULL)
    jsonlite::fromJSON(json_files[1], simplifyVector = FALSE)
}

gui_ai_qc_preview <- function(root, output_root = "spectreasy_outputs", detail = "standard") {
    readiness <- gui_ai_qc_readiness(root, output_root)
    object <- gui_ai_qc_latest_object(root, output_root)
    if (is.null(object)) return(c(readiness, list(summary = NULL, prompt = "")))
    class(object) <- c("spectreasy_ai_qc", "list")
    prompt <- spectreasy::build_ai_qc_prompt(object, detail = detail)
    c(readiness, list(
        schema = object$schema,
        profile = object$quality_reference$profile,
        scope = object$scope$value,
        privacy_mode = object$privacy$mode,
        grade_counts = object$grade_summary$counts,
        findings = object$overall_summary$top_findings,
        warnings = character(),
        missing_sections = unique(vapply(Filter(function(metric) !isTRUE(metric$availability), spectreasy:::.ai_qc_walk_metrics(object)), `[[`, character(1), "stage")),
        prompt = prompt,
        prompt_characters = nchar(prompt),
        estimated_tokens = ceiling(nchar(prompt) / 4)
    ))
}

gui_ai_qc_profiles <- function(root) {
    profiles <- spectreasy::list_qc_reference_profiles(root)
    list(profiles = profiles)
}

gui_ai_qc_generate <- function(body) {
    root <- gui_workflow_root(body)
    output_root <- gui_workflow_value(body, "output_root", "spectreasy_outputs")
    matrix_file <- gui_ai_qc_matrix_file(root, output_root)
    M <- if (!is.null(matrix_file)) spectreasy:::.read_unmixing_matrix_csv(matrix_file) else NULL
    if (is.null(M)) stop("No machine-readable reference matrix is available. Run the control stage first.", call. = FALSE)
    noise_file <- file.path(dirname(matrix_file), "scc_detector_noise.csv")
    if (file.exists(noise_file)) M <- spectreasy:::.attach_detector_noise(M, utils::read.csv(noise_file, check.names = FALSE), source = noise_file)
    scope <- gui_workflow_value(body, "scope", "control")
    result <- spectreasy::export_ai_qc(
        M = M,
        controls = if (scope %in% c("control", "combined")) list(M = M) else NULL,
        project_dir = root,
        existing_output_dir = file.path(root, output_root),
        output_dir = gui_ai_qc_output_dir(root, output_root),
        scope = scope,
        detail = gui_workflow_value(body, "detail", "standard"),
        privacy = gui_workflow_value(body, "privacy", "standard"),
        reference = gui_workflow_value(body, "reference", "auto"),
        context = gui_workflow_value(body, "context", NULL),
        overwrite = "overwrite"
    )
    preview <- gui_ai_qc_preview(root, output_root, detail = gui_workflow_value(body, "detail", "standard"))
    preview$status <- if (length(result$missing_sections)) "partial" else "ready"
    preview$artifact_paths <- lapply(result$paths, gui_ai_qc_relative, root = root)
    preview$content_sha256 <- result$content_sha256
    preview
}
