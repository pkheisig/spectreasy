# gui_api_adjust.R
# API for the spectreasy Interactive Tuner (Adjustment/Crosstalk Correction)

get_gui_default_project_dir <- function() {
    configured <- getOption("spectreasy.project_dir", "")
    if (!is.null(configured) && length(configured) > 0 && nzchar(trimws(as.character(configured[1])))) {
        return(normalizePath(as.character(configured[1]), mustWork = FALSE))
    }

    cwd <- normalizePath(getwd(), mustWork = FALSE)
    package_root <- normalizePath(file.path(cwd, "..", ".."), mustWork = FALSE)
    is_source_package <- file.exists(file.path(package_root, "DESCRIPTION")) &&
        file.exists(file.path(package_root, "inst", "api", "gui_api.R"))
    if (is_source_package) package_root else cwd
}

get_matrix_dir <- function() {
    configured <- getOption("spectreasy.matrix_dir", getOption("spectreasy.project_dir", ""))
    if (is.null(configured) || length(configured) == 0 || !nzchar(trimws(as.character(configured[1])))) {
        configured <- get_gui_default_project_dir()
    }
    normalizePath(as.character(configured[1]), mustWork = FALSE)
}

get_samples_dir <- function() {
    default_samples <- gui_project_input_path(get_matrix_dir(), "samples")
    normalizePath(getOption("spectreasy.samples_dir", default_samples), mustWork = FALSE)
}

get_unmixing_method <- function() {
    method <- getOption("spectreasy.unmixing_method", "AutoSpectral")
    spectreasy:::.normalize_unmix_method(method)
}

get_config_dir <- function() {
    get_user_gui_config_dir()
}

normalize_config_filename <- function(filename) {
    if (is.null(filename) || !nzchar(trimws(filename))) {
        return("gui_config.json")
    }
    out <- basename(trimws(filename))
    if (!grepl("\\.json$", out, ignore.case = TRUE)) {
        out <- paste0(out, ".json")
    }
    out
}

is_probably_matrix_csv <- function(path) {
    df <- tryCatch(
        utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, nrows = 200),
        error = function(e) NULL
    )
    if (is.null(df) || !is.data.frame(df) || ncol(df) < 2) {
        return(FALSE)
    }

    first_name <- tolower(trimws(colnames(df)[1]))
    first_col <- df[[1]]
    has_marker_col <- first_name %in% c("marker", "fluorophore", "file") || !is.numeric(first_col)
    mat_df <- if (has_marker_col) df[, -1, drop = FALSE] else df
    if (ncol(mat_df) == 0) {
        return(FALSE)
    }

    numeric_ok <- vapply(mat_df, function(col) {
        converted <- suppressWarnings(as.numeric(col))
        all(is.na(col) | !is.na(converted))
    }, logical(1))

    mean(numeric_ok) >= 0.8
}

matrix_filename_contains_matrix <- function(filename) {
    if (is.null(filename) || length(filename) == 0 || is.na(filename[1])) {
        return(FALSE)
    }
    normalized <- tolower(gsub("[^[:alnum:]]+", "", basename(as.character(filename)[1])))
    isTRUE(grepl("matrix", normalized, fixed = TRUE))
}

gui_request_project_root <- function(project_path = "", fallback = TRUE) {
    value <- if (is.null(project_path) || !length(project_path)) "" else trimws(as.character(project_path[1]))
    if (!nzchar(value) && isTRUE(fallback)) value <- get_matrix_dir()
    if (!nzchar(value)) stop("A project folder is required for this request.", call. = FALSE)
    if (!dir.exists(value)) stop("Project folder not found: ", value, call. = FALSE)
    normalizePath(value, mustWork = TRUE)
}

# Delegate the API surface to the package-level resolver so the launcher,
# cockpit, gating applets, and workflows share one layout implementation.
gui_project_layout_config_path <- function(project_path) spectreasy:::.spectreasy_project_layout_path(gui_request_project_root(project_path))
gui_project_input_marker_name <- function() spectreasy:::.spectreasy_project_layout_marker_name()
gui_normalize_project_input_dir <- function(value, root, label = "input directory") spectreasy:::.spectreasy_normalize_project_input_dir(value, gui_request_project_root(root), label)
gui_write_project_layout <- function(root, layout) spectreasy:::.spectreasy_write_project_layout(gui_request_project_root(root), layout)
gui_write_project_input_marker <- function(directory, role) spectreasy:::.spectreasy_write_project_input_marker(directory, role)
gui_discover_project_input_markers <- function(root) spectreasy:::.spectreasy_discover_project_input_markers(gui_request_project_root(root))
gui_project_layout <- function(root, reconcile = TRUE, ensure_markers = TRUE, persist = TRUE) {
    spectreasy:::.spectreasy_project_layout(
        gui_request_project_root(root),
        reconcile = reconcile,
        ensure_markers = ensure_markers,
        persist = persist
    )
}
gui_project_input_path <- function(root, role) spectreasy:::.spectreasy_project_input_path(gui_request_project_root(root), role)
gui_update_project_input_dir <- function(root, role, value) {
    root <- gui_request_project_root(root)
    layout <- spectreasy:::.spectreasy_update_project_input_dir(root, role, value)
    gui_set_project_context(root, reset_gate_cache = identical(tolower(trimws(as.character(role)[1])), "controls"))
    layout
}

gui_project_value <- function(body, fallback = TRUE) {
    value <- gui_workflow_value(body, "projectPath", gui_workflow_value(body, "project_path", ""))
    gui_request_project_root(value, fallback = fallback)
}

.gui_project_gate_sessions <- new.env(parent = emptyenv())

gui_reset_project_gate_session <- function(project_path) {
    key <- gui_request_project_root(project_path)
    if (exists(key, envir = .gui_project_gate_sessions, inherits = FALSE)) {
        rm(list = key, envir = .gui_project_gate_sessions)
    }
    invisible(NULL)
}

gui_with_project_context <- function(project_path, expr) {
    root <- gui_request_project_root(project_path)
    keys <- c(
        "spectreasy.project_dir", "spectreasy.matrix_dir", "spectreasy.samples_dir",
        "spectreasy.gating_scc_dir", "spectreasy.gating_control_file", "spectreasy.gating_gate_file",
        "spectreasy.gating_payload_cache", "spectreasy.gating_spectrum_cache",
        "spectreasy.gating_detector_cache", "spectreasy.gating_state_cache", "spectreasy.project_selected"
    )
    previous <- options()[keys]
    session_key <- normalizePath(root, mustWork = FALSE)
    session <- if (exists(session_key, envir = .gui_project_gate_sessions, inherits = FALSE)) {
        get(session_key, envir = .gui_project_gate_sessions, inherits = FALSE)
    } else list()
    gui_set_project_context(root, reset_gate_cache = FALSE)
    for (key in names(session)) options(stats::setNames(list(session[[key]]), key))
    on.exit({
        current <- options()[c(
            "spectreasy.gating_payload_cache", "spectreasy.gating_spectrum_cache",
            "spectreasy.gating_detector_cache", "spectreasy.gating_state_cache"
        )]
        assign(session_key, current, envir = .gui_project_gate_sessions)
        options(previous)
    }, add = TRUE)
    force(expr)
}

list_matrix_csv_files <- function(matrix_dir = get_matrix_dir()) {
    matrix_dir <- normalizePath(matrix_dir, mustWork = FALSE)
    if (!dir.exists(matrix_dir)) {
        return(character(0))
    }

    files <- list.files(
        matrix_dir,
        pattern = "\\.csv$",
        recursive = TRUE,
        full.names = FALSE,
        ignore.case = TRUE
    )
    files <- gsub("\\\\", "/", files)
    files <- files[vapply(files, matrix_filename_contains_matrix, logical(1))]
    sort(files)
}

matrix_path <- function(filename, matrix_dir = get_matrix_dir()) {
    if (is.null(filename) || length(filename) == 0 || is.na(filename[1])) {
        stop("Invalid matrix filename")
    }
    rel <- trimws(as.character(filename)[1])
    rel <- gsub("\\\\", "/", rel)
    rel <- sub("^/+", "", rel)
    parts <- strsplit(rel, "/", fixed = TRUE)[[1]]
    if (!isTRUE(nzchar(rel)) || any(parts %in% c("", ".", ".."))) {
        stop("Invalid matrix filename")
    }
    file.path(normalizePath(matrix_dir, mustWork = FALSE), do.call(file.path, as.list(parts)))
}

read_matrix_csv <- function(path) {
    df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    if (ncol(df) > 0 && colnames(df)[1] %in% c("V1", "")) {
        colnames(df)[1] <- "Marker"
    } else if (ncol(df) > 0 && colnames(df)[1] != "Marker" && is.character(df[[1]])) {
        colnames(df)[1] <- "Marker"
    }
    df
}

is_af_matrix_row <- function(df) {
    if (!is.data.frame(df) || ncol(df) == 0) return(logical(0))
    grepl("^AF($|_)", as.character(df[[1]]), ignore.case = TRUE)
}

merge_hidden_af_rows <- function(df, source_paths = character()) {
    if (!is.data.frame(df) || ncol(df) == 0) {
        return(df)
    }
    if (any(is_af_matrix_row(df))) {
        return(df)
    }

    source_paths <- unique(source_paths[file.exists(source_paths)])
    for (source_path in source_paths) {
        existing_df <- read_matrix_csv(source_path)
        af_rows <- existing_df[is_af_matrix_row(existing_df), , drop = FALSE]
        if (nrow(af_rows) == 0) next

        colnames(af_rows)[1] <- colnames(df)[1]
        af_rows_aligned <- af_rows[, intersect(colnames(df), colnames(af_rows)), drop = FALSE]
        missing_cols <- setdiff(colnames(df), colnames(af_rows_aligned))
        for (m_col in missing_cols) {
            af_rows_aligned[[m_col]] <- 0
        }
        af_rows_aligned <- af_rows_aligned[, colnames(df), drop = FALSE]
        return(rbind(df, af_rows_aligned))
    }

    df
}

raw_data_to_df <- function(raw_data_json) {
    if (is.data.frame(raw_data_json)) {
        return(as.data.frame(raw_data_json, check.names = FALSE))
    }
    if (is.list(raw_data_json) && length(raw_data_json) > 0 && all(vapply(raw_data_json, is.list, logical(1)))) {
        rows <- lapply(raw_data_json, function(row) {
            as.data.frame(row, check.names = FALSE, stringsAsFactors = FALSE)
        })
        return(do.call(rbind, rows))
    }
    as.data.frame(raw_data_json, check.names = FALSE)
}

get_user_gui_config_dir <- function() {
    cfg_dir <- file.path(tools::R_user_dir("spectreasy", which = "config"), "gui_configs")
    if (!dir.exists(cfg_dir)) {
        dir.create(cfg_dir, recursive = TRUE, showWarnings = FALSE)
    }
    cfg_dir
}

normalize_gui_module <- function(module) {
    module <- trimws(as.character(module)[1])
    if (is.na(module) || !nzchar(module)) module <- "matrix_tuner"
    module <- gsub("[^A-Za-z0-9_-]+", "_", module)
    module
}

user_gui_config_path <- function(module, project_path = "") {
    project_path <- if (is.null(project_path) || !length(project_path)) "" else trimws(as.character(project_path[1]))
    if (nzchar(project_path)) {
        root <- gui_request_project_root(project_path)
        directory <- file.path(root, ".spectreasy", "gui")
        dir.create(directory, recursive = TRUE, showWarnings = FALSE)
        return(file.path(directory, paste0(normalize_gui_module(module), ".json")))
    }
    file.path(get_user_gui_config_dir(), paste0(normalize_gui_module(module), ".json"))
}
