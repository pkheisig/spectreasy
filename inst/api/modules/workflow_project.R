# Workflow wrappers keep the browser orchestration layer small. All scientific
# work remains delegated to the exported Spectreasy functions below.
gui_workflow_body <- function(req) {
    tryCatch(
        jsonlite::fromJSON(req$postBody, simplifyVector = FALSE),
        error = function(e) stop("Invalid JSON request body: ", conditionMessage(e), call. = FALSE)
    )
}

gui_workflow_value <- function(body, key, fallback = NULL) {
    value <- body[[key]]
    if (is.null(value) || length(value) == 0 || !nzchar(trimws(as.character(value[1])))) {
        return(fallback)
    }
    as.character(value[1])
}

gui_workflow_bool <- function(body, key, fallback = FALSE) {
    value <- body[[key]]
    if (is.null(value) || length(value) == 0 || is.na(value[1])) return(isTRUE(fallback))
    if (is.logical(value[1])) return(isTRUE(value[1]))
    normalized <- tolower(trimws(as.character(value[1])))
    if (normalized %in% c("true", "1", "yes", "on")) return(TRUE)
    if (normalized %in% c("false", "0", "no", "off")) return(FALSE)
    isTRUE(fallback)
}

gui_workflow_number <- function(body, key, fallback = 0, integer = FALSE, minimum = NULL, maximum = NULL) {
    supplied <- !is.null(body[[key]]) && length(body[[key]]) > 0L
    value <- suppressWarnings(as.numeric(gui_workflow_value(body, key, fallback)))
    if (length(value) == 0 || !is.finite(value[1])) {
        if (supplied) stop("Invalid numeric value for '", key, "'.", call. = FALSE)
        value <- fallback
    }
    value <- as.numeric(value[1])
    if (!is.null(minimum) && value < minimum) {
        stop("'", key, "' must be at least ", minimum, ".", call. = FALSE)
    }
    if (!is.null(maximum) && value > maximum) {
        stop("'", key, "' must be at most ", maximum, ".", call. = FALSE)
    }
    if (isTRUE(integer)) {
        if (!isTRUE(all.equal(value, round(value)))) stop("'", key, "' must be an integer.", call. = FALSE)
        value <- as.integer(value)
    }
    value
}

gui_workflow_path <- function(body, key, fallback = "", allow_empty = TRUE) {
    value <- gui_workflow_value(body, key, fallback)
    if (is.null(value) || is.na(value)) return(if (isTRUE(allow_empty)) "" else fallback)
    value <- trimws(as.character(value[1]))
    if (!nzchar(value) && isTRUE(allow_empty)) return("")
    value
}

gui_method_optional_args <- function(method, body) {
    resolved <- tryCatch(spectreasy:::.normalize_unmix_method(method), error = function(e) "")
    if (identical(resolved, "Spectreasy")) {
        return(list(spectreasy_weight_quantile = gui_workflow_number(body, "spectreasy_weight_quantile", 0.65, minimum = 0, maximum = 1)))
    }
    list()
}

gui_workflow_root <- function(body) {
    root <- gui_workflow_value(body, "projectPath", get_matrix_dir())
    if (!dir.exists(root)) stop("Project folder not found: ", root, call. = FALSE)
    normalizePath(root, mustWork = TRUE)
}

gui_workflow_resolve_path <- function(path, root, allow_empty = FALSE) {
    if (is.null(path) || length(path) == 0L || is.na(path[1])) {
        if (isTRUE(allow_empty)) return("")
        stop("Project-relative path is missing.", call. = FALSE)
    }
    path <- trimws(as.character(path[1]))
    if (!nzchar(path)) {
        if (isTRUE(allow_empty)) return("")
        stop("Project-relative path is empty.", call. = FALSE)
    }
    root <- normalizePath(root, mustWork = TRUE)
    expanded <- path.expand(path)
    candidate_input <- if (grepl("^(/|[A-Za-z]:[/\\\\])", expanded)) expanded else file.path(root, expanded)

    # Check the nearest existing ancestor too, so a project symlink cannot be
    # used to route a future output path outside the project.
    ancestor <- candidate_input
    suffix <- character()
    while (!file.exists(ancestor) && !dir.exists(ancestor)) {
        parent <- dirname(ancestor)
        if (identical(parent, ancestor)) break
        suffix <- c(basename(ancestor), suffix)
        ancestor <- parent
    }
    ancestor <- normalizePath(ancestor, mustWork = TRUE)
    candidate <- if (length(suffix) > 0L) {
        do.call(file.path, c(list(ancestor), as.list(suffix)))
    } else {
        ancestor
    }
    candidate <- normalizePath(candidate, mustWork = FALSE)
    root_prefix <- paste0(root, .Platform$file.sep)
    inside <- function(value) identical(value, root) || startsWith(value, root_prefix)
    if (!inside(candidate) || !inside(ancestor)) {
        stop("Workflow path is outside the active project: ", path, call. = FALSE)
    }
    candidate
}

gui_workflow_run <- function(body, action, expr) {
    root <- gui_workflow_root(body)
    messages <- character()
    warnings <- character()
    tryCatch(
        {
            result <- NULL
            output <- capture.output(
                result <- withCallingHandlers(
                    force(expr),
                    message = function(condition) {
                        messages <<- c(messages, conditionMessage(condition))
                    },
                    warning = function(condition) {
                        warnings <<- c(warnings, paste0("Warning: ", conditionMessage(condition)))
                    }
                ),
                type = "output"
            )
            if (length(output)) cat(paste0(output, collapse = "\n"), "\n")
            list(
                success = TRUE,
                action = action,
                project_path = root,
                finished_at = as.character(Sys.time()),
                logs = c(output, messages, warnings),
                result = result
            )
        },
        error = function(e) {
            list(
                success = FALSE,
                action = action,
                project_path = root,
                logs = c(messages, warnings),
                error = conditionMessage(e),
                suggested_next_step = "Review the project inputs and the full action log, then retry from the relevant workflow card."
            )
        }
    )
}

gui_workflow_file_or_null <- function(path, root) {
    if (is.null(path) || length(path) == 0 || is.na(path[1])) return(NULL)
    path <- trimws(as.character(path[1]))
    if (!nzchar(path)) return(NULL)
    candidate <- gui_workflow_resolve_path(path, root)
    if (!file.exists(candidate)) return(NULL)
    normalizePath(candidate, mustWork = TRUE)
}

gui_active_af_config_path <- function(root = get_matrix_dir()) {
    file.path(normalizePath(root, mustWork = FALSE), ".spectreasy", "active_af_profile.json")
}

gui_read_active_af_profile <- function(root = get_matrix_dir()) {
    path <- gui_active_af_config_path(root)
    if (!file.exists(path)) return("")
    value <- tryCatch(jsonlite::fromJSON(path, simplifyVector = TRUE), error = function(e) NULL)
    name <- if (is.list(value)) value$profile_name else NULL
    if (is.null(name) || length(name) == 0L || is.na(name[1])) return("")
    name <- trimws(as.character(name[1]))
    profile_exists <- nzchar(name) && tryCatch(
        file.exists(spectreasy:::.af_profile_file(name, create_dir = FALSE)),
        error = function(e) FALSE
    )
    if (!profile_exists) {
        unlink(path, force = TRUE)
        return("")
    }
    name
}

gui_write_active_af_profile <- function(name, root = get_matrix_dir()) {
    name <- trimws(as.character(name)[1])
    if (!nzchar(name)) stop("A saved AF profile name is required.", call. = FALSE)
    profile <- spectreasy::load_af_profile(name, show_plot = FALSE)
    scc_dir <- gui_project_input_path(root, "controls")
    fcs_files <- if (dir.exists(scc_dir)) list.files(scc_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE) else character()
    fcs_files <- fcs_files[!startsWith(basename(fcs_files), "._")]
    if (length(fcs_files) > 0L) {
        detector_names <- spectreasy:::.prepare_reference_detector_info(fcs_files[1])$detector_names
        profile_detectors <- colnames(profile$profile)
        if (!setequal(detector_names, profile_detectors)) {
            stop("Saved AF profile detectors do not match this dataset's SCC detector set.", call. = FALSE)
        }
    }
    path <- gui_active_af_config_path(root)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile("active_af_profile_", tmpdir = dirname(path), fileext = ".json")
    jsonlite::write_json(
        list(profile_name = name, updated = as.character(Sys.time())),
        tmp,
        auto_unbox = TRUE,
        pretty = TRUE
    )
    if (!file.rename(tmp, path)) {
        if (!file.copy(tmp, path, overwrite = TRUE)) stop("Could not save the active AF profile selection.", call. = FALSE)
        unlink(tmp)
    }
    name
}

gui_unlink_active_af_profile <- function(name, root = get_matrix_dir()) {
    name <- trimws(as.character(name)[1])
    if (!nzchar(name)) stop("A saved AF profile name is required.", call. = FALSE)
    active <- gui_read_active_af_profile(root)
    if (!nzchar(active)) return(name)
    if (!identical(active, name)) {
        stop(name, " is not linked to this dataset.", call. = FALSE)
    }
    path <- gui_active_af_config_path(root)
    if (file.exists(path)) {
        status <- unlink(path, force = TRUE)
        if (!identical(status, 0L)) stop("Could not unlink the active AF profile.", call. = FALSE)
    }
    name
}

gui_primary_unstained_rows <- function(rows) {
    if (is.null(rows) || !is.data.frame(rows) || nrow(rows) == 0L) return(logical(0))
    fluorophore <- if ("fluorophore" %in% colnames(rows)) as.character(rows$fluorophore) else rep("", nrow(rows))
    marker <- if ("marker" %in% colnames(rows)) as.character(rows$marker) else rep("", nrow(rows))
    filename <- if ("filename" %in% colnames(rows)) as.character(rows$filename) else rep("", nrow(rows))
    control_type <- if ("control.type" %in% colnames(rows)) tolower(trimws(as.character(rows$control.type))) else rep("cells", nrow(rows))
    is_af <- grepl("^AF($|_|\\b)", trimws(fluorophore), ignore.case = TRUE) |
        grepl("autofluorescence|unstained", trimws(marker), ignore.case = TRUE) |
        grepl("unstained", basename(filename), ignore.case = TRUE)
    is_dead <- grepl("dead", paste(fluorophore, marker, filename), ignore.case = TRUE)
    is_bead <- grepl("bead", control_type, ignore.case = TRUE) | grepl("bead", basename(filename), ignore.case = TRUE)
    is_af & !is_dead & !is_bead
}

gui_annotate_active_af_mapping <- function(rows, root = get_matrix_dir()) {
    rows <- as.data.frame(rows, stringsAsFactors = FALSE, check.names = FALSE)
    active <- gui_read_active_af_profile(root)
    ignored <- rep(FALSE, nrow(rows))
    if (nzchar(active)) ignored <- gui_primary_unstained_rows(rows)
    rows$ignored <- ignored
    rows$ignored_reason <- ifelse(
        ignored,
        paste0("Ignored because saved AF profile '", active, "' is used as the unstained cell control."),
        ""
    )
    rows
}

gui_filter_active_af_mapping <- function(rows, root = get_matrix_dir()) {
    if (!nzchar(gui_read_active_af_profile(root)) || is.null(rows) || nrow(rows) == 0L) return(rows)
    ignored <- gui_primary_unstained_rows(rows)
    ignored_files <- if ("filename" %in% colnames(rows)) as.character(rows$filename[ignored]) else character()
    rows <- rows[!ignored, , drop = FALSE]
    if ("universal.negative" %in% colnames(rows)) {
        refs <- trimws(as.character(rows$universal.negative))
        rows$universal.negative[refs %in% c("AF", ignored_files)] <- ""
    }
    rows
}

gui_filtered_control_file <- function(root = get_matrix_dir()) {
    source <- file.path(root, "fcs_mapping.csv")
    if (!file.exists(source)) return(source)
    rows <- utils::read.csv(source, stringsAsFactors = FALSE, check.names = FALSE)
    filtered <- gui_filter_active_af_mapping(rows, root = root)
    if (nrow(filtered) == nrow(rows)) return(source)
    dir.create(file.path(root, ".spectreasy"), recursive = TRUE, showWarnings = FALSE)
    target <- tempfile("fcs_mapping_active_af_", tmpdir = file.path(root, ".spectreasy"), fileext = ".csv")
    utils::write.csv(filtered, target, row.names = FALSE, quote = TRUE)
    target
}

gui_project_scan <- function(root) {
    root <- gui_request_project_root(root)
    layout <- gui_project_layout(root, persist = FALSE)
    files <- if (dir.exists(root)) list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE) else character()
    files <- normalizePath(files[file.exists(files)], mustWork = FALSE)
    relative <- if (length(files) > 0) {
        sub(paste0("^", gsub("([.|(){}+*?^$\\[\\]\\\\])", "\\\\\\1", root), "[/\\\\]?"), "", files)
    } else {
        character()
    }
    relative <- gsub("\\\\", "/", relative)
    count_matches <- function(pattern) sum(grepl(pattern, relative, ignore.case = TRUE, perl = TRUE))
    control_files <- if (dir.exists(file.path(root, layout$control_input_dir))) {
        list.files(file.path(root, layout$control_input_dir), pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
    } else character()
    sample_files <- if (dir.exists(file.path(root, layout$sample_input_dir))) {
        list.files(file.path(root, layout$sample_input_dir), pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
    } else character()
    controls <- length(control_files)
    samples <- length(sample_files)
    matrices <- count_matches("matrix.*\\.csv$|unmixing.*\\.csv$|detector_noise.*\\.csv$")
    reports <- count_matches("\\.(html?|pdf)$")
    gates <- count_matches("gate.*\\.csv$")
    qc_metrics <- count_matches("metric.*\\.csv$")
    spectral_variants <- count_matches("variant.*\\.(rds|csv)$")
    summary <- if (length(files) == 0) {
        "empty project"
    } else if (matrices > 0 && reports > 0 && samples > 0) {
        "mixed/partial project"
    } else if (matrices > 0) {
        "reference built"
    } else if (controls > 0) {
        "controls imported"
    } else {
        "project detected"
    }
    list(
        project_path = root,
        layout = layout,
        missing_input_dirs = gui_project_relative_path(gui_missing_project_input_dirs(root), root),
        files = relative,
        scan = list(
            controls = controls,
            samples = samples,
            matrices = matrices,
            reports = reports,
            gates = gates,
            qc_metrics = qc_metrics,
            spectral_variants = spectral_variants
        ),
        summary = summary,
        recommended_next_action = if (matrices == 0) "Review controls and build a reference matrix" else if (samples > 0) "Run sample unmixing" else "Import samples"
    )
}

gui_missing_project_input_dirs <- function(project_path) {
    project_path <- normalizePath(project_path, mustWork = TRUE)
    layout <- gui_project_layout(project_path, ensure_markers = FALSE, persist = FALSE)
    paths <- file.path(project_path, c(layout$control_input_dir, layout$sample_input_dir))
    paths[!dir.exists(paths)]
}

gui_ensure_project_input_dirs <- function(project_path) {
    project_path <- normalizePath(project_path, mustWork = TRUE)
    layout <- gui_project_layout(project_path, ensure_markers = FALSE)
    roles <- c(controls = layout$control_input_dir, samples = layout$sample_input_dir)
    paths <- file.path(project_path, unname(roles))
    for (index in seq_along(paths)) {
        path <- paths[[index]]
        if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
            stop("Could not create project input folder: ", path, call. = FALSE)
        }
        gui_write_project_input_marker(path, names(roles)[[index]])
    }
    gui_write_project_layout(project_path, layout)
    invisible(paths)
}

gui_project_file_location <- function(kind, filename = NULL, project_path = "") {
    kind <- tolower(trimws(as.character(kind)[1]))
    root <- gui_request_project_root(project_path)
    layout <- gui_project_layout(root, persist = FALSE)
    folder <- switch(kind, controls = layout$control_input_dir, samples = layout$sample_input_dir, NULL)
    if (is.null(folder)) stop("File kind must be 'controls' or 'samples'.", call. = FALSE)
    directory <- file.path(root, folder)
    if (is.null(filename)) return(list(kind = kind, folder = folder, directory = directory))
    candidate <- gsub("\\\\", "/", trimws(as.character(filename)[1]))
    safe_name <- basename(candidate)
    if (!nzchar(safe_name) || !identical(candidate, safe_name) || safe_name %in% c(".", "..")) {
        stop("Invalid project filename.", call. = FALSE)
    }
    if (!grepl("\\.fcs$", safe_name, ignore.case = TRUE)) {
        stop("Only FCS files can be managed here.", call. = FALSE)
    }
    list(kind = kind, folder = folder, directory = directory, filename = safe_name, path = file.path(directory, safe_name))
}

gui_project_file_rows <- function(kind, project_path = "") {
    location <- gui_project_file_location(kind, project_path = project_path)
    if (!dir.exists(location$directory)) return(data.frame())
    files <- list.files(location$directory, pattern = "\\.fcs$", ignore.case = TRUE, full.names = TRUE)
    if (!length(files)) return(data.frame())
    info <- file.info(files)
    rows <- data.frame(
        name = basename(files),
        size = as.numeric(info$size),
        modified = format(info$mtime, "%Y-%m-%dT%H:%M:%S%z"),
        modified_epoch = as.numeric(info$mtime),
        kind = location$kind,
        stringsAsFactors = FALSE
    )
    rows[order(gui_natural_sort_key(rows$name), tolower(rows$name)), , drop = FALSE]
}

.gui_project_uploads <- new.env(parent = emptyenv())

gui_project_upload_id <- function() {
    paste(sample(c(letters, LETTERS, 0:9), 32L, replace = TRUE), collapse = "")
}

gui_project_upload_session <- function(upload_id) {
    upload_id <- trimws(as.character(upload_id)[1])
    if (!nzchar(upload_id) || !exists(upload_id, envir = .gui_project_uploads, inherits = FALSE)) {
        stop("Upload session not found or expired.", call. = FALSE)
    }
    get(upload_id, envir = .gui_project_uploads, inherits = FALSE)
}

gui_project_upload_discard <- function(upload_id) {
    upload_id <- trimws(as.character(upload_id)[1])
    if (nzchar(upload_id) && exists(upload_id, envir = .gui_project_uploads, inherits = FALSE)) {
        session <- get(upload_id, envir = .gui_project_uploads, inherits = FALSE)
        unlink(session$temporary, force = TRUE)
        rm(list = upload_id, envir = .gui_project_uploads)
    }
    invisible(NULL)
}

gui_project_upload_discard_stale <- function(max_age_seconds = 3600) {
    upload_ids <- ls(envir = .gui_project_uploads, all.names = TRUE)
    now <- Sys.time()
    for (upload_id in upload_ids) {
        session <- get(upload_id, envir = .gui_project_uploads, inherits = FALSE)
        age <- suppressWarnings(as.numeric(difftime(now, session$created, units = "secs")))
        if (!is.finite(age) || age > max_age_seconds) gui_project_upload_discard(upload_id)
    }
    invisible(NULL)
}

gui_validate_fcs_upload <- function(path) {
    size <- file.info(path)$size
    if (!is.finite(size) || size < 6) stop("Uploaded FCS file is empty or truncated.", call. = FALSE)
    connection <- file(path, open = "rb")
    on.exit(close(connection), add = TRUE)
    header <- rawToChar(readBin(connection, what = "raw", n = 6L))
    if (!grepl("^FCS[0-9]\\.[0-9]$", header)) {
        stop("Uploaded file does not contain a valid FCS header.", call. = FALSE)
    }
    invisible(TRUE)
}

gui_natural_sort_key <- function(x) {
    vapply(as.character(x), function(value) {
        parts <- regmatches(value, gregexpr("[0-9]+|[^0-9]+", value, perl = TRUE))[[1]]
        paste(vapply(parts, function(part) {
            if (grepl("^[0-9]+$", part)) sprintf("%020.0f", as.numeric(part)) else tolower(part)
        }, character(1)), collapse = "")
    }, character(1))
}

gui_set_project_context <- function(project_path, reset_gate_cache = TRUE) {
    if (!dir.exists(project_path)) stop("Project folder not found: ", project_path, call. = FALSE)
    project_path <- normalizePath(project_path, mustWork = TRUE)
    layout <- gui_project_layout(project_path)
    project_options <- list(
        spectreasy.project_dir = project_path,
        spectreasy.matrix_dir = project_path,
        spectreasy.samples_dir = file.path(project_path, layout$sample_input_dir),
        spectreasy.gating_scc_dir = file.path(project_path, layout$control_input_dir),
        spectreasy.gating_control_file = file.path(project_path, "fcs_mapping.csv"),
        spectreasy.gating_gate_file = file.path(project_path, "ssc_gate_config.csv")
    )
    if (isTRUE(reset_gate_cache)) {
        project_options <- c(project_options, list(
            spectreasy.gating_payload_cache = list(),
            spectreasy.gating_spectrum_cache = list(),
            spectreasy.gating_detector_cache = list(),
            spectreasy.gating_state_cache = NULL
        ))
    }
    options(project_options)
    options(spectreasy.project_selected = TRUE)
    project_path
}

gui_pick_project_directory <- function(initial_dir = path.expand("~"), allow_create = FALSE) {
    initial_dir <- normalizePath(initial_dir, mustWork = FALSE)
    prompt <- if (isTRUE(allow_create)) "Create a Spectreasy project folder" else "Open a Spectreasy project folder"
    if (identical(Sys.info()[["sysname"]], "Darwin")) {
        script <- paste0(
            "POSIX path of (choose folder with prompt ", gate_applescript_quote(prompt), " default location POSIX file ",
            gate_applescript_quote(initial_dir), ")"
        )
        result <- suppressWarnings(system2("osascript", c("-e", shQuote(script)), stdout = TRUE, stderr = FALSE))
        if (!is.null(attr(result, "status")) || length(result) == 0L) return(NULL)
        return(sub("/$", "", trimws(result[[1]])))
    }
    if (.Platform$OS.type == "windows") {
        result <- utils::choose.dir(default = initial_dir, caption = prompt)
        if (is.na(result)) return(NULL)
        return(result)
    }
    picker <- Sys.which("zenity")
    if (!nzchar(picker)) picker <- Sys.which("kdialog")
    if (!nzchar(picker)) stop("No graphical folder picker is available. Install zenity or kdialog.", call. = FALSE)
    args <- if (grepl("zenity$", picker)) {
        c("--file-selection", "--directory", paste0("--title=", prompt), paste0("--filename=", initial_dir, "/"))
    } else {
        c("--getexistingdirectory", initial_dir, "--title", prompt)
    }
    result <- suppressWarnings(system2(picker, args, stdout = TRUE, stderr = FALSE))
    if (!is.null(attr(result, "status")) || length(result) == 0L) return(NULL)
    trimws(result[[1]])
}
