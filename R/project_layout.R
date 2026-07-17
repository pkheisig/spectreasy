# Internal project-layout helpers shared by the cockpit launcher and API.

.spectreasy_project_layout_defaults <- function() {
    list(control_input_dir = "scc", sample_input_dir = "samples")
}

.spectreasy_project_layout_marker_name <- function() ".spectreasy-input-role.json"

.spectreasy_normalize_project_input_dir <- function(value, root, label = "input directory") {
    value <- gsub("\\\\", "/", trimws(as.character(value)[1]))
    value <- sub("/+$", "", value)
    if (!nzchar(value)) stop(label, " cannot be empty.", call. = FALSE)
    if (grepl("^(/|[A-Za-z]:[/\\\\])", value)) stop(label, " must be a path inside the active project.", call. = FALSE)
    parts <- strsplit(value, "/", fixed = TRUE)[[1]]
    if (any(parts %in% c("", ".", ".."))) stop(label, " contains an invalid path component.", call. = FALSE)
    if (identical(parts[[1]], ".spectreasy")) stop(label, " cannot be stored inside the .spectreasy configuration directory.", call. = FALSE)
    root <- normalizePath(root, mustWork = TRUE)
    candidate <- normalizePath(file.path(root, do.call(file.path, as.list(parts))), mustWork = FALSE)
    if (!startsWith(candidate, paste0(root, .Platform$file.sep))) stop(label, " must remain inside the active project.", call. = FALSE)
    paste(parts, collapse = "/")
}

.spectreasy_project_layout_path <- function(root) {
    file.path(normalizePath(root, mustWork = TRUE), ".spectreasy", "project.json")
}

.spectreasy_write_project_layout <- function(root, layout) {
    root <- normalizePath(root, mustWork = TRUE)
    control_dir <- .spectreasy_normalize_project_input_dir(layout$control_input_dir, root, "Control input directory")
    sample_dir <- .spectreasy_normalize_project_input_dir(layout$sample_input_dir, root, "Sample input directory")
    if (identical(control_dir, sample_dir)) stop("Control and sample input directories must be different.", call. = FALSE)
    path <- .spectreasy_project_layout_path(root)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    payload <- list(version = 1L, control_input_dir = control_dir, sample_input_dir = sample_dir, updated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
    temporary <- tempfile("project_layout_", tmpdir = dirname(path), fileext = ".json")
    jsonlite::write_json(payload, temporary, auto_unbox = TRUE, pretty = TRUE)
    if (!file.rename(temporary, path)) {
        if (!file.copy(temporary, path, overwrite = TRUE)) stop("Could not save the project input-directory configuration.", call. = FALSE)
        unlink(temporary, force = TRUE)
    }
    payload
}

.spectreasy_write_project_input_marker <- function(directory, role) {
    if (!dir.exists(directory)) return(invisible(FALSE))
    jsonlite::write_json(
        list(role = role, version = 1L),
        file.path(directory, .spectreasy_project_layout_marker_name()),
        auto_unbox = TRUE,
        pretty = TRUE
    )
    invisible(TRUE)
}

.spectreasy_relative_project_path <- function(path, root) {
    root <- normalizePath(root, mustWork = TRUE)
    root_pattern <- gsub("([.|(){}+*?^$\\[\\]\\\\])", "\\\\\\1", root)
    gsub("\\\\", "/", sub(paste0("^", root_pattern, "[/\\\\]?"), "", normalizePath(path, mustWork = FALSE)))
}

.spectreasy_discover_project_input_markers <- function(root) {
    root <- normalizePath(root, mustWork = TRUE)
    marker_name <- .spectreasy_project_layout_marker_name()
    markers <- list.files(root, pattern = paste0("^", gsub("\\.", "\\\\.", marker_name), "$"), recursive = TRUE, full.names = TRUE, all.files = TRUE)
    result <- list(controls = character(), samples = character())
    for (marker in markers) {
        value <- tryCatch(jsonlite::fromJSON(marker, simplifyVector = TRUE), error = function(e) NULL)
        role <- if (is.null(value$role)) "" else trimws(as.character(value$role)[1])
        if (!role %in% names(result)) next
        relative <- tryCatch(.spectreasy_normalize_project_input_dir(.spectreasy_relative_project_path(dirname(marker), root), root), error = function(e) "")
        if (nzchar(relative)) result[[role]] <- unique(c(result[[role]], relative))
    }
    result
}

.spectreasy_project_layout <- function(root, reconcile = TRUE, ensure_markers = TRUE, persist = TRUE) {
    root <- normalizePath(root, mustWork = TRUE)
    defaults <- .spectreasy_project_layout_defaults()
    path <- .spectreasy_project_layout_path(root)
    stored <- if (file.exists(path)) tryCatch(jsonlite::fromJSON(path, simplifyVector = TRUE), error = function(e) list()) else list()
    layout <- list(
        control_input_dir = tryCatch(.spectreasy_normalize_project_input_dir(stored$control_input_dir %||% defaults$control_input_dir, root), error = function(e) defaults$control_input_dir),
        sample_input_dir = tryCatch(.spectreasy_normalize_project_input_dir(stored$sample_input_dir %||% defaults$sample_input_dir, root), error = function(e) defaults$sample_input_dir)
    )
    changed <- !file.exists(path)
    if (isTRUE(reconcile)) {
        discovered <- .spectreasy_discover_project_input_markers(root)
        role_fields <- c(controls = "control_input_dir", samples = "sample_input_dir")
        for (role in names(role_fields)) {
            field <- role_fields[[role]]
            if (!dir.exists(file.path(root, layout[[field]])) && length(discovered[[role]]) == 1L) {
                layout[[field]] <- discovered[[role]][[1]]
                changed <- TRUE
            }
        }
    }
    if (identical(layout$control_input_dir, layout$sample_input_dir)) {
        layout$sample_input_dir <- defaults$sample_input_dir
        changed <- TRUE
    }
    if (isTRUE(ensure_markers) && isTRUE(persist)) {
        .spectreasy_write_project_input_marker(file.path(root, layout$control_input_dir), "controls")
        .spectreasy_write_project_input_marker(file.path(root, layout$sample_input_dir), "samples")
    }
    if (changed && isTRUE(persist)) .spectreasy_write_project_layout(root, layout)
    layout
}

.spectreasy_project_input_path <- function(root, role) {
    root <- normalizePath(root, mustWork = TRUE)
    layout <- .spectreasy_project_layout(root)
    field <- switch(tolower(trimws(as.character(role)[1])), controls = "control_input_dir", control = "control_input_dir", samples = "sample_input_dir", sample = "sample_input_dir", stop("Input-directory role must be controls or samples.", call. = FALSE))
    normalizePath(file.path(root, layout[[field]]), mustWork = FALSE)
}

.spectreasy_update_project_input_dir <- function(root, role, value) {
    root <- normalizePath(root, mustWork = TRUE)
    role <- tolower(trimws(as.character(role)[1]))
    field <- switch(role, controls = "control_input_dir", samples = "sample_input_dir", NULL)
    if (is.null(field)) stop("Input-directory role must be controls or samples.", call. = FALSE)
    layout <- .spectreasy_project_layout(root)
    next_value <- .spectreasy_normalize_project_input_dir(value, root, paste(if (role == "controls") "Control" else "Sample", "input directory"))
    other_field <- if (field == "control_input_dir") "sample_input_dir" else "control_input_dir"
    if (identical(next_value, layout[[other_field]])) stop("Control and sample input directories must be different.", call. = FALSE)
    previous_value <- layout[[field]]
    previous_path <- file.path(root, previous_value)
    next_path <- file.path(root, next_value)
    if (!identical(previous_value, next_value)) {
        if (dir.exists(previous_path) && dir.exists(next_path)) stop("The requested input directory already exists. Remove it or choose a different name before renaming.", call. = FALSE)
        if (dir.exists(previous_path)) {
            dir.create(dirname(next_path), recursive = TRUE, showWarnings = FALSE)
            if (!file.rename(previous_path, next_path)) stop("Could not rename the input directory. Check filesystem permissions and try again.", call. = FALSE)
        } else if (!dir.exists(next_path) && !dir.create(next_path, recursive = TRUE, showWarnings = FALSE)) {
            stop("Could not create the requested input directory.", call. = FALSE)
        }
    }
    layout[[field]] <- next_value
    layout <- .spectreasy_write_project_layout(root, layout)
    .spectreasy_write_project_input_marker(next_path, role)
    layout
}
