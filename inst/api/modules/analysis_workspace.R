# Project-scoped post-unmixing analysis helpers. This module deliberately keeps
# editable analysis state outside FCS files: FCS exports contain selected events
# and provenance, while the workspace JSON remains the hierarchy source of truth.

if (!exists("%||%", mode = "function", inherits = TRUE)) {
    `%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x
}

.gui_analysis_cache <- new.env(parent = emptyenv())
.gui_analysis_python_cache <- new.env(parent = emptyenv())
.gui_analysis_jobs <- new.env(parent = emptyenv())

gui_analysis_python_executable <- function() {
    configured <- c(
        as.character(getOption("spectreasy.analysis_python", ""))[1],
        Sys.getenv("SPECTREASY_ANALYSIS_PYTHON", unset = ""),
        file.path(path.expand("~/Library/Caches/spectreasy/analysis-python"), "bin", "python"),
        file.path(path.expand("~/.cache/spectreasy/analysis-python"), "bin", "python"),
        Sys.which("python3")
    )
    configured <- unique(configured[nzchar(configured)])
    hits <- configured[file.exists(configured) & file.access(configured, 1L) == 0L]
    if (!length(hits)) return("")
    path.expand(hits[[1]])
}

gui_analysis_python_bridge <- function() {
    frame_files <- unlist(lapply(sys.frames(), function(frame) {
        value <- frame$ofile
        if (is.null(value) || !length(value)) character() else as.character(value[[1]])
    }), use.names = FALSE)
    candidates <- c(
        if (length(frame_files)) file.path(dirname(dirname(dirname(frame_files))), "python", "analysis_bridge.py") else character(),
        file.path(getwd(), "inst", "python", "analysis_bridge.py"),
        file.path(getwd(), "..", "inst", "python", "analysis_bridge.py"),
        system.file("python", "analysis_bridge.py", package = "spectreasy")
    )
    candidates <- unique(candidates[nzchar(candidates)])
    hits <- candidates[file.exists(candidates)]
    if (!length(hits)) return("")
    normalizePath(hits[[1]], mustWork = TRUE)
}

gui_analysis_worker_script <- function() {
    candidates <- c(
        file.path(getwd(), "inst", "python", "analysis_worker.R"),
        file.path(getwd(), "..", "inst", "python", "analysis_worker.R"),
        system.file("python", "analysis_worker.R", package = "spectreasy")
    )
    hits <- candidates[file.exists(candidates)]
    if (!length(hits)) return("")
    normalizePath(hits[[1]], mustWork = TRUE)
}

gui_analysis_python_packages <- function(refresh = FALSE) {
    python <- gui_analysis_python_executable()
    bridge <- gui_analysis_python_bridge()
    cache_key <- paste(python, bridge, sep = "\u241f")
    if (!isTRUE(refresh) && exists(cache_key, envir = .gui_analysis_python_cache, inherits = FALSE)) {
        return(get(cache_key, envir = .gui_analysis_python_cache, inherits = FALSE))
    }
    package_names <- c("phate", "nptsne", "scanpy", "palantir", "spectreasy_builtin")
    empty <- stats::setNames(
        lapply(package_names, function(x) list(available = FALSE, version = "")),
        package_names
    )
    if (!nzchar(python) || !nzchar(bridge)) return(empty)
    output <- suppressWarnings(system2(python, shQuote(c(bridge, "probe")), stdout = TRUE, stderr = TRUE))
    if (!identical(attr(output, "status") %||% 0L, 0L)) return(empty)
    parsed <- tryCatch(jsonlite::fromJSON(paste(output, collapse = "\n"), simplifyVector = FALSE), error = function(e) NULL)
    if (!is.list(parsed)) return(empty)
    for (name in names(empty)) {
        if (is.list(parsed[[name]])) empty[[name]] <- parsed[[name]]
    }
    assign(cache_key, empty, envir = .gui_analysis_python_cache)
    empty
}

gui_analysis_python_run <- function(arguments) {
    python <- gui_analysis_python_executable()
    bridge <- gui_analysis_python_bridge()
    if (!nzchar(python) || !nzchar(bridge)) {
        stop("The managed Spectreasy Python analysis runtime is unavailable.", call. = FALSE)
    }
    output <- suppressWarnings(system2(python, shQuote(c(bridge, arguments)), stdout = TRUE, stderr = TRUE))
    status <- attr(output, "status") %||% 0L
    if (!identical(as.integer(status), 0L)) {
        detail <- trimws(paste(output, collapse = "\n"))
        if (!nzchar(detail)) detail <- paste("Python adapter exited with status", status)
        stop(detail, call. = FALSE)
    }
    invisible(output)
}

gui_analysis_relative_path <- function(path, root) {
    root <- normalizePath(root, mustWork = TRUE)
    path <- normalizePath(path, mustWork = FALSE)
    prefix <- paste0(root, .Platform$file.sep)
    if (identical(path, root)) return(".")
    if (!startsWith(path, prefix)) stop("Analysis path is outside the active project.", call. = FALSE)
    gsub("\\\\", "/", substring(path, nchar(prefix) + 1L))
}

gui_analysis_stable_integer <- function(...) {
    text <- enc2utf8(paste(..., collapse = "\u241f"))
    value <- 17
    for (code in utf8ToInt(text)) value <- (value * 131 + code) %% 2147483646
    as.integer(value + 1L)
}

gui_analysis_stable_id <- function(prefix, ...) {
    paste0(prefix, "-", sprintf("%08x", gui_analysis_stable_integer(...)))
}

gui_analysis_workspace_directory <- function(root, create = FALSE) {
    root <- gui_request_project_root(root)
    path <- file.path(root, ".spectreasy", "analysis-v2")
    if (isTRUE(create) && !dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
        stop("Could not create the analysis workspace directory.", call. = FALSE)
    }
    path
}

gui_analysis_workspace_path <- function(root, create = FALSE) {
    file.path(gui_analysis_workspace_directory(root, create = create), "workspace.json")
}

gui_analysis_default_workspace <- function() {
    list(
        schema_version = 2L,
        updated_at = NULL,
        source_path = "",
        selected_file = "",
        active_population_id = "root",
        seed = 20260723L,
        populations = list(list(
            id = "root", name = "All events", parent_id = NULL, type = "root",
            role = NULL, source_file = NULL, geometry = NULL
        )),
        plots = list(list(
            id = "plot-1", type = "scatter", population_id = "root",
            x = "FSC-A", y = "SSC-A", color_by = "density",
            x_transform = "linear", y_transform = "linear"
        )),
        annotations = list(),
        root_event_id = NULL,
        root_population_id = NULL,
        root_source_file = NULL
    )
}

gui_analysis_normalize_workspace <- function(value) {
    if (!is.list(value)) stop("Analysis workspace must be a JSON object.", call. = FALSE)
    out <- gui_analysis_default_workspace()
    for (name in names(value)) out[[name]] <- value[[name]]
    if (!is.list(out$populations) || length(out$populations) < 1L) {
        stop("Analysis workspace must contain a population hierarchy.", call. = FALSE)
    }
    if (length(out$populations) > 500L) stop("Analysis workspace has too many populations.", call. = FALSE)
    ids <- vapply(out$populations, function(node) trimws(as.character(node$id %||% "")[1]), character(1))
    if (any(!nzchar(ids)) || anyDuplicated(ids)) stop("Population IDs must be non-empty and unique.", call. = FALSE)
    if (!"root" %in% ids) stop("Analysis workspace is missing its root population.", call. = FALSE)
    out$schema_version <- 2L
    out$updated_at <- as.character(Sys.time())
    out
}

gui_analysis_read_workspace <- function(root) {
    path <- gui_analysis_workspace_path(root)
    if (!file.exists(path)) return(gui_analysis_default_workspace())
    value <- tryCatch(jsonlite::fromJSON(path, simplifyVector = FALSE), error = function(e) {
        stop("Could not read the analysis workspace: ", conditionMessage(e), call. = FALSE)
    })
    gui_analysis_normalize_workspace(value)
}

gui_analysis_write_workspace <- function(root, value) {
    value <- gui_analysis_normalize_workspace(value)
    path <- gui_analysis_workspace_path(root, create = TRUE)
    temporary <- tempfile("analysis-workspace-", tmpdir = dirname(path), fileext = ".json")
    jsonlite::write_json(value, temporary, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = 15)
    if (!file.rename(temporary, path)) {
        if (!file.copy(temporary, path, overwrite = TRUE)) stop("Could not save the analysis workspace.", call. = FALSE)
        unlink(temporary, force = TRUE)
    }
    value
}

gui_analysis_fcs_directories <- function(root) {
    root <- gui_request_project_root(root)
    layout <- gui_project_layout(root, persist = FALSE)
    candidates <- c(
        file.path(root, layout$sample_input_dir),
        file.path(root, layout$control_input_dir)
    )
    output_roots <- list.dirs(root, recursive = FALSE, full.names = TRUE)
    output_roots <- output_roots[grepl("^spectreasy_outputs", basename(output_roots), ignore.case = TRUE)]
    if (length(output_roots)) {
        derived <- unlist(lapply(output_roots, list.dirs, recursive = TRUE, full.names = TRUE), use.names = FALSE)
        candidates <- c(candidates, derived[grepl("(^|[/\\\\])(unmixed_fcs|scc_unmixed)$", derived, ignore.case = TRUE)])
    }
    candidates <- unique(normalizePath(candidates[dir.exists(candidates)], mustWork = TRUE))
    candidates[vapply(candidates, function(path) {
        length(list.files(path, pattern = "[.]fcs$", ignore.case = TRUE, full.names = FALSE)) > 0L
    }, logical(1))]
}

gui_analysis_header <- function(path) {
    header <- flowCore::read.FCSheader(path)[[1]]
    n_parameters <- suppressWarnings(as.integer(unname(header["$PAR"])))
    n_events <- suppressWarnings(as.integer(unname(header["$TOT"])))
    indexes <- if (is.finite(n_parameters) && n_parameters > 0L) seq_len(n_parameters) else integer()
    channel_names <- unname(header[paste0("$P", indexes, "N")])
    descriptions <- unname(header[paste0("$P", indexes, "S")])
    descriptions[is.na(descriptions)] <- ""
    list(
        event_count = n_events,
        parameter_count = n_parameters,
        channels = as.character(channel_names),
        descriptions = as.character(descriptions),
        written_by = if ("WRITTEN_BY" %in% names(header)) unname(header["WRITTEN_BY"]) else ""
    )
}

gui_analysis_source_inventory <- function(root) {
    root <- gui_request_project_root(root)
    layout <- gui_project_layout(root, persist = FALSE)
    sample_dir <- normalizePath(file.path(root, layout$sample_input_dir), mustWork = FALSE)
    control_dir <- normalizePath(file.path(root, layout$control_input_dir), mustWork = FALSE)
    directories <- gui_analysis_fcs_directories(root)
    lapply(directories, function(directory) {
        files <- list.files(directory, pattern = "[.]fcs$", full.names = TRUE, ignore.case = TRUE)
        files <- files[!startsWith(basename(files), "._")]
        info <- file.info(files)
        rows <- lapply(seq_along(files), function(index) {
            path <- files[[index]]
            metadata <- tryCatch(gui_analysis_header(path), error = function(e) list(
                event_count = NA_integer_, parameter_count = NA_integer_, channels = character(),
                descriptions = character(), written_by = "", error = conditionMessage(e)
            ))
            relative <- gui_analysis_relative_path(path, root)
            c(list(
                id = gui_analysis_stable_id("file", relative),
                name = basename(path),
                path = relative,
                size = as.numeric(info$size[[index]]),
                modified = format(info$mtime[[index]], "%Y-%m-%dT%H:%M:%S%z")
            ), metadata)
        })
        relative_directory <- gui_analysis_relative_path(directory, root)
        role <- if (identical(directory, sample_dir)) "samples" else if (identical(directory, control_dir)) "controls" else "unmixed"
        list(
            id = gui_analysis_stable_id("source", relative_directory),
            label = switch(role, samples = "Sample input", controls = "SCC controls", unmixed = paste("Unmixed", basename(dirname(directory)))),
            role = role,
            path = relative_directory,
            files = rows
        )
    })
}

gui_analysis_resolve_fcs <- function(root, relative_path) {
    root <- gui_request_project_root(root)
    path <- gui_workflow_resolve_path(relative_path, root)
    if (!file.exists(path) || !grepl("[.]fcs$", path, ignore.case = TRUE)) {
        stop("Analysis FCS file not found: ", relative_path, call. = FALSE)
    }
    normalizePath(path, mustWork = TRUE)
}

gui_analysis_flow_frame <- function(root, relative_path) {
    path <- gui_analysis_resolve_fcs(root, relative_path)
    info <- file.info(path)
    signature <- paste(path, info$size, as.numeric(info$mtime), sep = ":")
    key <- gui_analysis_stable_id("frame", path)
    cached <- if (exists(key, envir = .gui_analysis_cache, inherits = FALSE)) get(key, envir = .gui_analysis_cache) else NULL
    if (is.list(cached) && identical(cached$signature, signature)) return(cached$frame)
    frame <- flowCore::read.FCS(path, transformation = FALSE, truncate_max_range = FALSE)
    rm(list = ls(.gui_analysis_cache, all.names = TRUE), envir = .gui_analysis_cache)
    assign(key, list(signature = signature, frame = frame), envir = .gui_analysis_cache)
    frame
}

gui_analysis_population <- function(workspace, population_id) {
    matches <- Filter(function(node) identical(as.character(node$id)[1], as.character(population_id)[1]), workspace$populations)
    if (!length(matches)) stop("Population not found: ", population_id, call. = FALSE)
    matches[[1]]
}

gui_analysis_polygon_mask <- function(x, y, points) {
    if (!is.list(points) || length(points) < 3L) return(rep(FALSE, length(x)))
    coordinate <- function(point, name, index) {
        value <- if (is.list(point) && !is.null(point[[name]])) point[[name]] else if (length(point) >= index) point[[index]] else NA_real_
        as.numeric(value)[1]
    }
    px <- vapply(points, coordinate, numeric(1), name = "x", index = 1L)
    py <- vapply(points, coordinate, numeric(1), name = "y", index = 2L)
    valid <- is.finite(px) & is.finite(py)
    if (sum(valid) < 3L) return(rep(FALSE, length(x)))
    sp::point.in.polygon(x, y, px[valid], py[valid]) > 0L
}

gui_analysis_gate_mask <- function(data, workspace, population_id = "root", source_file = "") {
    memo <- new.env(parent = emptyenv())
    visiting <- character()
    evaluate <- function(id) {
        if (exists(id, envir = memo, inherits = FALSE)) return(get(id, envir = memo))
        if (id %in% visiting) stop("Population hierarchy contains a cycle.", call. = FALSE)
        visiting <<- c(visiting, id)
        on.exit(visiting <<- setdiff(visiting, id), add = TRUE)
        node <- gui_analysis_population(workspace, id)
        if (identical(id, "root") || identical(node$type, "root")) {
            mask <- rep(TRUE, nrow(data))
        } else {
            parent_id <- as.character(node$parent_id %||% "root")[1]
            mask <- evaluate(parent_id)
            scoped_file <- trimws(as.character(node$source_file %||% "")[1])
            if (nzchar(scoped_file) && !identical(scoped_file, source_file)) {
                mask[] <- FALSE
            } else {
                geometry <- node$geometry %||% list()
                x_name <- as.character(node$x %||% geometry$x %||% "")[1]
                y_name <- as.character(node$y %||% geometry$y %||% "")[1]
                if (!nzchar(x_name) || !x_name %in% colnames(data)) stop("Gate x channel is unavailable: ", x_name, call. = FALSE)
                x <- data[, x_name]
                type <- tolower(as.character(node$type %||% "polygon")[1])
                local <- switch(type,
                    rectangle = {
                        if (!nzchar(y_name) || !y_name %in% colnames(data)) stop("Gate y channel is unavailable: ", y_name, call. = FALSE)
                        x >= as.numeric(geometry$x_min) & x <= as.numeric(geometry$x_max) &
                            data[, y_name] >= as.numeric(geometry$y_min) & data[, y_name] <= as.numeric(geometry$y_max)
                    },
                    ellipse = {
                        if (!nzchar(y_name) || !y_name %in% colnames(data)) stop("Gate y channel is unavailable: ", y_name, call. = FALSE)
                        center_x <- as.numeric(geometry$center_x)
                        center_y <- as.numeric(geometry$center_y)
                        radius_x <- as.numeric(geometry$radius_x)
                        radius_y <- as.numeric(geometry$radius_y)
                        if (!all(is.finite(c(center_x, center_y, radius_x, radius_y))) || radius_x <= 0 || radius_y <= 0) {
                            stop("Ellipse gate geometry is invalid.", call. = FALSE)
                        }
                        ((x - center_x) / radius_x)^2 + ((data[, y_name] - center_y) / radius_y)^2 <= 1
                    },
                    range = x >= as.numeric(geometry$min) & x <= as.numeric(geometry$max),
                    polygon = {
                        if (!nzchar(y_name) || !y_name %in% colnames(data)) stop("Gate y channel is unavailable: ", y_name, call. = FALSE)
                        gui_analysis_polygon_mask(x, data[, y_name], geometry$points)
                    },
                    stop("Unsupported gate type: ", type, call. = FALSE)
                )
                local[is.na(local)] <- FALSE
                mask <- mask & local
            }
        }
        assign(id, mask, envir = memo)
        mask
    }
    evaluate(as.character(population_id)[1])
}

gui_analysis_sample_indices <- function(indices, n, seed, identity = "") {
    indices <- as.integer(indices)
    if (!length(indices)) return(indices)
    n <- as.integer(n)
    if (!is.finite(n) || n <= 0L || n >= length(indices)) return(indices)
    derived <- gui_analysis_stable_integer(seed, identity)
    sort(withr::with_seed(derived, sample(indices, size = n, replace = FALSE)))
}

gui_analysis_marker_labels <- function(frame) {
    parameters <- flowCore::parameters(frame)@data
    descriptions <- if ("desc" %in% colnames(parameters)) as.character(parameters$desc) else rep("", nrow(parameters))
    descriptions[is.na(descriptions)] <- ""
    data.frame(channel = colnames(flowCore::exprs(frame)), marker = descriptions, stringsAsFactors = FALSE)
}

gui_analysis_event_payload <- function(root, file, population_id = "root", x = "FSC-A", y = "SSC-A", color = "", max_points = 8000L, seed = 20260723L) {
    workspace <- gui_analysis_read_workspace(root)
    frame <- gui_analysis_flow_frame(root, file)
    data <- flowCore::exprs(frame)
    channels <- colnames(data)
    if (!x %in% channels) x <- channels[[1]]
    if (!y %in% channels) y <- channels[[min(2L, length(channels))]]
    if (!nzchar(color) || !color %in% channels) color <- ""
    mask <- gui_analysis_gate_mask(data, workspace, population_id, source_file = file)
    population_indices <- which(mask)
    selected <- gui_analysis_sample_indices(population_indices, max_points, seed, paste(file, population_id, "display"))
    events <- data.frame(
        event_id = selected,
        x = unname(data[selected, x]),
        y = unname(data[selected, y]),
        stringsAsFactors = FALSE
    )
    if (nzchar(color)) events$color <- unname(data[selected, color])
    parent_id <- as.character(gui_analysis_population(workspace, population_id)$parent_id %||% "root")[1]
    parent_count <- if (identical(population_id, "root")) nrow(data) else sum(gui_analysis_gate_mask(data, workspace, parent_id, source_file = file))
    list(
        file = file,
        population_id = population_id,
        x = x,
        y = y,
        color = color,
        events = events,
        population_count = length(population_indices),
        parent_count = parent_count,
        total_count = nrow(data),
        displayed_count = length(selected),
        channels = gui_analysis_marker_labels(frame)
    )
}

gui_analysis_population_path <- function(workspace, population_id) {
    names <- character()
    seen <- character()
    current <- population_id
    while (nzchar(current) && !current %in% seen) {
        seen <- c(seen, current)
        node <- gui_analysis_population(workspace, current)
        names <- c(as.character(node$name %||% current)[1], names)
        current <- as.character(node$parent_id %||% "")[1]
    }
    paste(names, collapse = "/")
}

gui_analysis_population_statistics <- function(root, file, population_id = "root", markers = character()) {
    workspace <- gui_analysis_read_workspace(root)
    frame <- gui_analysis_flow_frame(root, file)
    data <- flowCore::exprs(frame)
    mask <- gui_analysis_gate_mask(data, workspace, population_id, source_file = file)
    node <- gui_analysis_population(workspace, population_id)
    parent_id <- as.character(node$parent_id %||% "root")[1]
    parent_mask <- if (identical(population_id, "root")) rep(TRUE, nrow(data)) else gui_analysis_gate_mask(data, workspace, parent_id, source_file = file)
    markers <- intersect(as.character(unlist(markers, use.names = FALSE)), colnames(data))
    summaries <- lapply(markers, function(marker) {
        values <- data[mask, marker]
        data.frame(
            marker = marker,
            median = if (length(values)) stats::median(values, na.rm = TRUE) else NA_real_,
            mean = if (length(values)) mean(values, na.rm = TRUE) else NA_real_,
            robust_sd = if (length(values) > 1L) stats::mad(values, constant = 1.4826, na.rm = TRUE) else NA_real_,
            stringsAsFactors = FALSE
        )
    })
    list(
        population_id = population_id,
        population_name = as.character(node$name %||% population_id)[1],
        count = sum(mask),
        parent_count = sum(parent_mask),
        total_count = nrow(data),
        percent_parent = if (sum(parent_mask)) 100 * sum(mask) / sum(parent_mask) else 0,
        percent_total = if (nrow(data)) 100 * sum(mask) / nrow(data) else 0,
        markers = if (length(summaries)) do.call(rbind, summaries) else data.frame()
    )
}

gui_analysis_staining_index <- function(root, file, marker) {
    workspace <- gui_analysis_read_workspace(root)
    frame <- gui_analysis_flow_frame(root, file)
    data <- flowCore::exprs(frame)
    if (!marker %in% colnames(data)) stop("Staining-index marker is unavailable: ", marker, call. = FALSE)
    roles <- vapply(workspace$populations, function(node) tolower(as.character(node$role %||% "")[1]), character(1))
    positive <- workspace$populations[roles == "positive"]
    negative <- workspace$populations[roles == "negative"]
    if (length(positive) != 1L || length(negative) != 1L) {
        stop("Assign exactly one POS and one NEG population before calculating staining index.", call. = FALSE)
    }
    positive <- positive[[1]]
    negative <- negative[[1]]
    if (!identical(as.character(positive$parent_id %||% "root")[1], as.character(negative$parent_id %||% "root")[1])) {
        stop("POS and NEG populations must share the same parent population.", call. = FALSE)
    }
    positive_values <- data[gui_analysis_gate_mask(data, workspace, positive$id, source_file = file), marker]
    negative_values <- data[gui_analysis_gate_mask(data, workspace, negative$id, source_file = file), marker]
    if (length(positive_values) < 2L || length(negative_values) < 3L) stop("POS or NEG population has too few events.", call. = FALSE)
    mfi_positive <- stats::median(positive_values, na.rm = TRUE)
    mfi_negative <- stats::median(negative_values, na.rm = TRUE)
    robust_sd <- stats::mad(negative_values, center = mfi_negative, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(robust_sd) || robust_sd <= 0) stop("NEG robust SD is zero or unavailable.", call. = FALSE)
    list(
        marker = marker,
        positive_population_id = positive$id,
        negative_population_id = negative$id,
        positive_count = length(positive_values),
        negative_count = length(negative_values),
        mfi = "median",
        robust_sd_method = "1.4826 * MAD",
        mfi_positive = mfi_positive,
        mfi_negative = mfi_negative,
        robust_sd_negative = robust_sd,
        staining_index = (mfi_positive - mfi_negative) / (2 * robust_sd)
    )
}

gui_analysis_export <- function(root, body) {
    root <- gui_request_project_root(root)
    workspace <- gui_analysis_read_workspace(root)
    files <- as.character(unlist(body$files %||% body$file, use.names = FALSE))
    files <- files[nzchar(files)]
    if (!length(files)) stop("Select at least one FCS file to export.", call. = FALSE)
    population_id <- gui_workflow_value(body, "populationId", "root")
    format <- tolower(gui_workflow_value(body, "format", "fcs"))
    if (!format %in% c("fcs", "csv", "both")) stop("Export format must be fcs, csv, or both.", call. = FALSE)
    target_count <- gui_workflow_number(body, "maxEvents", 0, integer = TRUE, minimum = 0)
    seed <- gui_workflow_number(body, "seed", workspace$seed %||% 20260723L, integer = TRUE, minimum = 1)
    output_value <- gui_workflow_path(body, "outputFolder", "spectreasy_outputs/analysis/exports")
    output <- gui_workflow_resolve_path(output_value, root)
    if (!dir.exists(output) && !dir.create(output, recursive = TRUE, showWarnings = FALSE)) stop("Could not create the export folder.", call. = FALSE)
    population_path <- gui_analysis_population_path(workspace, population_id)
    records <- lapply(files, function(file) {
        frame <- gui_analysis_flow_frame(root, file)
        data <- flowCore::exprs(frame)
        mask <- gui_analysis_gate_mask(data, workspace, population_id, source_file = file)
        population_indices <- which(mask)
        selected <- gui_analysis_sample_indices(population_indices, target_count, seed, paste(file, population_id, "export"))
        safe_population <- gsub("[^[:alnum:]_.-]+", "_", gui_analysis_population(workspace, population_id)$name %||% population_id)
        stem <- tools::file_path_sans_ext(basename(file))
        paths <- character()
        if (format %in% c("fcs", "both")) {
            output_frame <- frame[selected, ]
            keywords <- flowCore::keyword(output_frame)
            keywords$WRITTEN_BY <- "Spectreasy v2"
            keywords$SpectreasyPopulationPath <- population_path
            keywords$SpectreasyPopulationId <- population_id
            keywords$SpectreasySeed <- as.character(seed)
            keywords$SpectreasySourceFile <- file
            flowCore::keyword(output_frame) <- keywords
            fcs_path <- file.path(output, paste0(stem, "_", safe_population, ".fcs"))
            flowCore::write.FCS(output_frame, fcs_path)
            paths <- c(paths, fcs_path)
        }
        if (format %in% c("csv", "both")) {
            csv_path <- file.path(output, paste0(stem, "_", safe_population, ".csv"))
            csv <- data.frame(
                sample_name = stem,
                sample_id = gui_analysis_stable_id("sample", file),
                event_id = selected,
                data[selected, , drop = FALSE],
                check.names = FALSE,
                stringsAsFactors = FALSE
            )
            data.table::fwrite(csv, csv_path)
            paths <- c(paths, csv_path)
        }
        lapply(paths, function(path) list(
            path = gui_analysis_relative_path(path, root),
            events = length(selected),
            source_events = length(population_indices),
            sha256 = gui_analysis_file_sha256(path)
        ))
    })
    list(
        population_id = population_id,
        population_path = population_path,
        seed = seed,
        max_events = target_count,
        files = unlist(records, recursive = FALSE)
    )
}

gui_analysis_file_sha256 <- function(path) {
    connection <- file(path, "rb")
    on.exit(close(connection), add = TRUE)
    unname(paste0(openssl::sha256(connection)))
}

gui_analysis_export_statistics <- function(root, body) {
    root <- gui_request_project_root(root)
    workspace <- gui_analysis_read_workspace(root)
    files <- as.character(unlist(body$files %||% body$file, use.names = FALSE))
    files <- files[nzchar(files)]
    if (!length(files)) stop("Select at least one FCS file for statistics export.", call. = FALSE)
    markers <- as.character(unlist(body$markers %||% character(), use.names = FALSE))
    population_ids <- as.character(unlist(body$populationIds %||% vapply(workspace$populations, `[[`, character(1), "id"), use.names = FALSE))
    rows <- unlist(lapply(files, function(file) {
        unlist(lapply(population_ids, function(population_id) {
            statistics <- gui_analysis_population_statistics(root, file, population_id, markers)
            marker_rows <- statistics$markers
            if (!nrow(marker_rows)) marker_rows <- data.frame(marker = "", median = NA_real_, mean = NA_real_, robust_sd = NA_real_)
            lapply(seq_len(nrow(marker_rows)), function(index) data.frame(
                source_file = file,
                population_id = population_id,
                population_path = gui_analysis_population_path(workspace, population_id),
                count = statistics$count,
                parent_count = statistics$parent_count,
                total_count = statistics$total_count,
                percent_parent = statistics$percent_parent,
                percent_total = statistics$percent_total,
                marker = marker_rows$marker[[index]],
                median = marker_rows$median[[index]],
                mean = marker_rows$mean[[index]],
                robust_sd = marker_rows$robust_sd[[index]],
                stringsAsFactors = FALSE
            ))
        }), recursive = FALSE)
    }), recursive = FALSE)
    table <- do.call(rbind, rows)
    output_value <- gui_workflow_path(body, "outputFolder", "spectreasy_outputs/analysis/exports")
    output <- gui_workflow_resolve_path(output_value, root)
    if (!dir.exists(output) && !dir.create(output, recursive = TRUE, showWarnings = FALSE)) stop("Could not create the statistics export folder.", call. = FALSE)
    path <- file.path(output, "population_statistics.csv")
    data.table::fwrite(table, path)
    list(
        path = gui_analysis_relative_path(path, root),
        rows = nrow(table), files = length(files), populations = length(population_ids),
        sha256 = gui_analysis_file_sha256(path)
    )
}

gui_analysis_parameter <- function(
    id, label, type, default, minimum = NULL, maximum = NULL, step = NULL,
    choices = NULL, description = ""
) {
    list(
        id = id,
        label = label,
        type = type,
        default = default,
        minimum = minimum,
        maximum = maximum,
        step = step,
        choices = choices,
        description = description
    )
}

gui_analysis_method_parameter_definitions <- function() {
    number <- function(id, label, default, minimum, maximum, step, description = "") {
        gui_analysis_parameter(id, label, "number", default, minimum, maximum, step, description = description)
    }
    integer <- function(id, label, default, minimum, maximum, step = 1L, description = "") {
        gui_analysis_parameter(id, label, "integer", default, minimum, maximum, step, description = description)
    }
    boolean <- function(id, label, default, description = "") {
        gui_analysis_parameter(id, label, "boolean", default, description = description)
    }
    select <- function(id, label, default, choices, description = "") {
        gui_analysis_parameter(
            id, label, "select", default,
            choices = lapply(names(choices), function(value) list(value = value, label = unname(choices[[value]]))),
            description = description
        )
    }
    list(
        pca = list(
            boolean("center", "Center markers", TRUE, "Subtract each marker mean before decomposition."),
            boolean("scale", "Scale markers", TRUE, "Scale each marker to unit variance before decomposition.")
        ),
        tsne = list(
            number("perplexity", "Perplexity", 30, 2, 200, 1, "Effective neighborhood size; automatically capped for small populations."),
            number("theta", "Barnes-Hut theta", 0.5, 0, 1, 0.05, "Zero is exact t-SNE; larger values are faster approximations."),
            integer("iterations", "Iterations", 1000L, 250L, 5000L, 50L),
            number("learning_rate", "Learning rate", 200, 10, 2000, 10),
            number("exaggeration", "Early exaggeration", 12, 1, 50, 0.5),
            boolean("normalize", "Normalize input", TRUE)
        ),
        umap = list(
            integer("neighbors", "Neighbors (k)", 15L, 2L, 200L, description = "Local neighborhood size."),
            number("min_dist", "Minimum distance", 0.01, 0, 0.99, 0.01),
            number("spread", "Spread", 1, 0.1, 10, 0.1),
            select("metric", "Distance metric", "euclidean", c(euclidean = "Euclidean", cosine = "Cosine", manhattan = "Manhattan", correlation = "Correlation")),
            integer("epochs", "Training epochs", 500L, 10L, 5000L, 10L),
            number("learning_rate", "Learning rate", 1, 0.01, 10, 0.05),
            number("repulsion", "Repulsion strength", 1, 0.1, 10, 0.1),
            integer("negative_samples", "Negative samples", 5L, 1L, 50L)
        ),
        `diffusion-map` = list(
            integer("neighbors", "Neighbors (k)", 15L, 2L, 200L),
            integer("eigenvectors", "Eigenvectors", 10L, 3L, 50L),
            select("distance", "Distance metric", "euclidean", c(euclidean = "Euclidean", cosine = "Cosine", rankcor = "Rank correlation", l2 = "L2")),
            boolean("density_normalization", "Density normalization", TRUE)
        ),
        phate = list(
            integer("neighbors", "Neighbors (k)", 5L, 2L, 200L),
            integer("decay", "Kernel decay", 40L, 1L, 100L),
            number("gamma", "Potential gamma", 1, -1, 1, 0.05),
            integer("landmarks", "Landmarks", 2000L, 20L, 10000L, 10L),
            integer("pca_components", "PCA components", 100L, 2L, 500L),
            select("mds", "MDS type", "metric", c(metric = "Metric", nonmetric = "Non-metric", classic = "Classical")),
            select("distance", "kNN distance", "euclidean", c(euclidean = "Euclidean", cosine = "Cosine", manhattan = "Manhattan"))
        ),
        hsne = list(
            integer("scales", "Hierarchy scales", 3L, 2L, 5L),
            integer("iterations", "Embedding iterations", 500L, 100L, 3000L, 50L)
        ),
        flowsom = list(
            integer("clusters", "Metaclusters", 10L, 2L, 100L),
            integer("xdim", "SOM grid width", 10L, 2L, 30L),
            integer("ydim", "SOM grid height", 10L, 2L, 30L),
            integer("rlen", "SOM training iterations", 10L, 1L, 100L),
            number("alpha_start", "Initial learning rate", 0.05, 0.001, 1, 0.005),
            number("alpha_end", "Final learning rate", 0.01, 0.0001, 1, 0.001)
        ),
        phenograph = list(
            integer("neighbors", "Neighbors (k)", 30L, 2L, 200L)
        ),
        dpt = list(
            integer("neighbors", "Neighbors (k)", 15L, 2L, 200L),
            integer("eigenvectors", "Diffusion eigenvectors", 10L, 3L, 50L),
            select("distance", "Distance metric", "euclidean", c(euclidean = "Euclidean", cosine = "Cosine", rankcor = "Rank correlation", l2 = "L2")),
            boolean("density_normalization", "Density normalization", TRUE)
        ),
        slingshot = list(
            select("distance", "Cluster distance", "simple", c(simple = "Euclidean centers (robust)", slingshot = "Slingshot covariance-aware")),
            number("shrink", "Shared-lineage shrinkage", 1, 0, 1, 0.05),
            select("extend", "Curve extension", "y", c(y = "Data range", n = "None", pc1 = "Terminal PC1")),
            boolean("reweight", "Reweight shared cells", TRUE),
            boolean("reassign", "Reassign curve membership", TRUE),
            integer("iterations", "Maximum iterations", 15L, 1L, 100L),
            number("stretch", "Curve stretch", 2, 0, 10, 0.25),
            integer("approximation_points", "Curve approximation points", 150L, 20L, 1000L, 10L),
            select("smoother", "Curve smoother", "smooth.spline", c(`smooth.spline` = "Smoothing spline", loess = "LOESS")),
            boolean("allow_breaks", "Allow diverging curves to separate", TRUE)
        ),
        tscan = list(),
        palantir = list(
            integer("neighbors", "Neighbors (k)", 30L, 2L, 200L),
            integer("diffusion_components", "Diffusion components", 10L, 3L, 50L),
            number("diffusion_alpha", "Diffusion alpha", 0, 0, 1, 0.05),
            integer("waypoints", "Waypoints", 1200L, 20L, 5000L, 10L),
            integer("iterations", "Maximum iterations", 25L, 1L, 200L),
            boolean("scale_components", "Scale diffusion components", TRUE)
        ),
        `paga-dpt` = list(
            integer("neighbors", "Neighbors (k)", 15L, 2L, 200L),
            select("metric", "Neighbor distance", "euclidean", c(euclidean = "Euclidean", cosine = "Cosine", manhattan = "Manhattan")),
            select("paga_model", "PAGA connectivity model", "v1.2", c(`v1.2` = "PAGA v1.2", `v1.0` = "PAGA v1.0")),
            integer("diffusion_components", "Diffusion components", 10L, 3L, 50L),
            integer("branchings", "DPT branchings", 0L, 0L, 10L),
            number("minimum_group_size", "Minimum group size", 0.01, 0.001, 0.5, 0.005)
        ),
        wanderlust = list(
            integer("neighbors", "Retained neighbors (k)", 15L, 3L, 100L, description = "Neighbors retained in each bootstrapped graph."),
            integer("candidate_neighbors", "Candidate neighbors (l)", 20L, 4L, 200L, description = "Must be greater than retained neighbors."),
            integer("graphs", "Graph ensemble size", 5L, 1L, 50L),
            integer("waypoints", "Waypoints", 64L, 10L, 400L, 10L),
            integer("iterations", "Alignment iterations", 15L, 2L, 100L),
            select("weighting", "Waypoint weighting", "exponential", c(exponential = "Exponential", linear = "Linear", uniform = "Uniform"))
        ),
        wishbone = list(
            integer("neighbors", "Retained neighbors (k)", 15L, 3L, 100L, description = "Neighbors retained in each bootstrapped graph."),
            integer("candidate_neighbors", "Candidate neighbors (l)", 20L, 4L, 200L, description = "Must be greater than retained neighbors."),
            integer("graphs", "Graph ensemble size", 5L, 1L, 50L),
            integer("waypoints", "Waypoints", 96L, 20L, 400L, 10L),
            integer("iterations", "Alignment iterations", 15L, 2L, 100L),
            integer("diffusion_components", "Diffusion components", 10L, 3L, 50L),
            number("branch_confidence", "Branch confidence", 0.65, 0.5, 0.95, 0.01)
        )
    )
}

gui_analysis_method_parameters <- function(method, supplied = NULL) {
    definitions <- method$parameters %||% list()
    supplied <- supplied %||% list()
    if (!is.list(supplied)) stop(method$name, " advanced settings must be an object.", call. = FALSE)
    known <- vapply(definitions, `[[`, character(1), "id")
    unknown <- setdiff(names(supplied), known)
    if (length(unknown)) stop(method$name, " has no advanced setting named ", unknown[[1]], ".", call. = FALSE)
    values <- lapply(definitions, function(definition) {
        id <- definition$id
        value <- supplied[[id]] %||% definition$default
        if (identical(definition$type, "boolean")) {
            if (!is.logical(value) || length(value) != 1L || is.na(value)) {
                stop(definition$label, " must be true or false.", call. = FALSE)
            }
            return(isTRUE(value))
        }
        if (identical(definition$type, "select")) {
            value <- as.character(value)[1]
            choices <- vapply(definition$choices, `[[`, character(1), "value")
            if (!value %in% choices) stop(definition$label, " has an unsupported value.", call. = FALSE)
            return(value)
        }
        value <- suppressWarnings(as.numeric(value)[1])
        if (!is.finite(value)) stop(definition$label, " must be numeric.", call. = FALSE)
        if (!is.null(definition$minimum) && value < definition$minimum) stop(definition$label, " is below its minimum.", call. = FALSE)
        if (!is.null(definition$maximum) && value > definition$maximum) stop(definition$label, " is above its maximum.", call. = FALSE)
        if (identical(definition$type, "integer")) value <- as.integer(round(value))
        value
    })
    resolved <- stats::setNames(values, known)
    if (identical(method$id, "flowsom") && resolved$alpha_end > resolved$alpha_start) {
        stop("FlowSOM final learning rate cannot exceed its initial learning rate.", call. = FALSE)
    }
    if (
        method$id %in% c("wanderlust", "wishbone") &&
        resolved$candidate_neighbors <= resolved$neighbors
    ) {
        stop(method$name, " candidate neighbors must exceed retained neighbors.", call. = FALSE)
    }
    resolved
}

gui_analysis_request_method_parameters <- function(body, method) {
    all_values <- body$advancedSettings %||% body$advanced_settings %||% list()
    supplied <- if (is.list(all_values)) all_values[[method$id]] %||% list() else list()
    gui_analysis_method_parameters(method, supplied)
}

gui_analysis_method_registry <- function() {
    python_packages <- gui_analysis_python_packages()
    parameter_definitions <- gui_analysis_method_parameter_definitions()
    method <- function(
        id, name, family, package, citation, doi, requirements = character(),
        research_only = FALSE, adapter_verified = TRUE,
        prerequisites = character(), outputs = character(), pipeline = id,
        user_prerequisites = character(), supports_3d = FALSE,
        python_package = "", visible = TRUE,
        additional_packages = character()
    ) {
        python_info <- if (nzchar(python_package)) python_packages[[python_package]] %||% list(available = FALSE, version = "") else NULL
        r_packages <- unique(c(package, additional_packages))
        r_packages <- r_packages[nzchar(r_packages)]
        package_available <- if (nzchar(python_package)) {
            isTRUE(python_info$available)
        } else if (!length(r_packages)) {
            TRUE
        } else {
            all(vapply(r_packages, requireNamespace, logical(1), quietly = TRUE))
        }
        available <- isTRUE(package_available) && isTRUE(adapter_verified) && !isTRUE(research_only)
        availability_state <- if (available) {
            "ready"
        } else if (!isTRUE(package_available) && isTRUE(adapter_verified) && !isTRUE(research_only)) {
            "setup-required"
        } else {
            "unavailable"
        }
        blocker <- paste(as.character(unlist(requirements, use.names = FALSE)), collapse = "; ")
        next_action <- if (identical(availability_state, "setup-required")) {
            if (nzchar(python_package)) {
                paste0("Install the maintained Spectreasy Python analysis runtime with ", python_package, ", then reopen this workspace.")
            } else {
                paste0("Install the ", paste(r_packages, collapse = " and "), " package", if (length(r_packages) > 1L) "s" else "", ", restart Spectreasy, and reopen this workspace.")
            }
        } else if (identical(availability_state, "unavailable")) {
            "This method is not available in this build. Choose an enabled method."
        } else {
            ""
        }
        list(
            id = id, name = name, family = family,
            package = if (nzchar(python_package)) python_package else package,
            runtime = if (nzchar(python_package)) "python" else "R",
            available = available,
            availability_state = availability_state,
            availability_label = switch(
                availability_state,
                ready = "Ready",
                `setup-required` = "Setup required",
                "Unavailable"
            ),
            installed = isTRUE(package_available),
            research_only = isTRUE(research_only),
            adapter_verified = isTRUE(adapter_verified),
            version = if (nzchar(python_package)) {
                as.character(python_info$version %||% "")
            } else if (package_available && nzchar(package)) {
                as.character(utils::packageVersion(package))
            } else "",
            citation = citation, doi = doi, requirements = requirements,
            blocker = blocker, next_action = next_action,
            prerequisites = prerequisites,
            automatic_prerequisites = prerequisites,
            user_prerequisites = user_prerequisites,
            outputs = outputs, pipeline = pipeline,
            supports_3d = isTRUE(supports_3d),
            visible = isTRUE(visible),
            parameters = parameter_definitions[[id]] %||% list()
        )
    }
    list(
        method("pca", "PCA", "reduction", "", "Pearson (1901), Philosophical Magazine", "10.1080/14786440109462720", outputs = "embedding", supports_3d = TRUE),
        method("tsne", "t-SNE", "reduction", "Rtsne", "van der Maaten and Hinton (2008), JMLR", "", outputs = "embedding", supports_3d = TRUE),
        method("umap", "UMAP", "reduction", "uwot", "McInnes, Healy and Melville (2018)", "10.48550/arXiv.1802.03426", outputs = "embedding", supports_3d = TRUE),
        method("diffusion-map", "Diffusion map", "reduction", "destiny", "Angerer et al. (2016), Bioinformatics", "10.1093/bioinformatics/btv715", outputs = c("diffusion_model", "embedding"), supports_3d = TRUE),
        method("phate", "PHATE", "reduction", "", "Moon et al. (2019), Nature Biotechnology", "10.1038/s41587-019-0336-3", "Python phate 2.0.0 or newer", outputs = "embedding", supports_3d = TRUE, python_package = "phate"),
        method("hsne", "HSNE", "reduction", "", "Pezzotti et al. (2016), Computer Graphics Forum", "10.1111/cgf.12878", "Python nptsne 2.0.3 or newer", outputs = c("hsne_hierarchy", "landmark_embedding"), python_package = "nptsne"),
        method("flowsom", "FlowSOM", "clustering", "FlowSOM", "Van Gassen et al. (2015), Cytometry A", "10.1002/cyto.a.22625", outputs = "cluster_labels"),
        method("phenograph", "PhenoGraph", "clustering", "Rphenograph", "Levine et al. (2015), Cell", "10.1016/j.cell.2015.05.047", outputs = "cluster_labels"),
        method("dpt", "Diffusion pseudotime", "trajectory", "destiny", "Haghverdi et al. (2016), Nature Methods", "10.1038/nmeth.3971", "Select a trajectory root event in the gating workspace", prerequisites = "diffusion-map", user_prerequisites = "trajectory-root", outputs = c("embedding", "pseudotime"), pipeline = c("diffusion-map", "dpt"), supports_3d = TRUE),
        method("slingshot", "Slingshot", "trajectory", "slingshot", "Street et al. (2018), BMC Genomics", "10.1186/s12864-018-4772-0", "Current Bioconductor slingshot with DelayedMatrixStats", prerequisites = c("reduction", "clustering"), user_prerequisites = "trajectory-root", outputs = c("lineages", "pseudotime"), pipeline = c("reduction", "clustering", "slingshot"), additional_packages = "DelayedMatrixStats"),
        method("tscan", "TSCAN", "trajectory", "TSCAN", "Ji and Ji (2016), Nucleic Acids Research", "10.1093/nar/gkw430", "Current Bioconductor TSCAN", prerequisites = c("reduction", "clustering"), user_prerequisites = "trajectory-root", outputs = c("minimum_spanning_tree", "pseudotime"), pipeline = c("reduction", "clustering", "tscan")),
        method("palantir", "Palantir", "trajectory", "", "Setty et al. (2019), Nature Biotechnology", "10.1038/s41587-019-0068-4", "Python palantir 1.4.5 or newer", prerequisites = "diffusion-map", user_prerequisites = "trajectory-root", outputs = c("pseudotime", "branch_probabilities"), supports_3d = TRUE, python_package = "palantir"),
        method("paga-dpt", "PAGA + DPT", "trajectory", "", "Wolf et al. (2019), Genome Biology", "10.1186/s13059-019-1663-x", "Python scanpy 1.12.2 or newer", prerequisites = c("neighbor_graph", "clustering"), user_prerequisites = "trajectory-root", outputs = c("graph", "pseudotime"), supports_3d = TRUE, python_package = "scanpy"),
        method("wanderlust", "Wanderlust", "trajectory", "", "Bendall et al. (2014), Cell", "10.1016/j.cell.2014.04.005", "Spectreasy independent maintained rewrite", prerequisites = "neighbor_graph", user_prerequisites = "trajectory-root", outputs = "pseudotime", supports_3d = TRUE, python_package = "spectreasy_builtin"),
        method("wishbone", "Wishbone", "trajectory", "", "Setty et al. (2016), Nature Biotechnology", "10.1038/nbt.3569", "Spectreasy independent maintained rewrite", prerequisites = "diffusion-map", user_prerequisites = "trajectory-root", outputs = c("pseudotime", "branch_labels"), supports_3d = TRUE, python_package = "spectreasy_builtin")
    )
}

gui_analysis_method <- function(id) {
    methods <- gui_analysis_method_registry()
    hits <- Filter(function(method) identical(method$id, id), methods)
    if (!length(hits)) stop("Unknown analysis method: ", id, call. = FALSE)
    hits[[1]]
}

gui_analysis_numeric_matrix <- function(data, markers, cofactor = 150) {
    markers <- intersect(as.character(unlist(markers, use.names = FALSE)), colnames(data))
    if (length(markers) < 2L) stop("Select at least two available markers.", call. = FALSE)
    matrix <- as.matrix(data[, markers, drop = FALSE])
    matrix <- asinh(matrix / cofactor)
    keep <- apply(matrix, 2, function(column) stats::sd(column, na.rm = TRUE) > 0)
    matrix <- matrix[, keep, drop = FALSE]
    if (ncol(matrix) < 2L) stop("Selected markers do not contain enough variation.", call. = FALSE)
    matrix
}

gui_analysis_marker_display_labels <- function(frame, channels) {
    channels <- as.character(channels)
    parameters <- tryCatch(flowCore::parameters(frame)@data, error = function(e) NULL)
    labels <- rep("", length(channels))
    if (is.data.frame(parameters) && all(c("name", "desc") %in% names(parameters))) {
        matched <- match(channels, as.character(parameters$name))
        descriptions <- as.character(parameters$desc)
        valid <- !is.na(matched)
        labels[valid] <- descriptions[matched[valid]]
    }
    labels[is.na(labels)] <- ""
    fallback <- sub("-(A|H|W)$", "", channels, ignore.case = TRUE)
    labels[!nzchar(trimws(labels))] <- fallback[!nzchar(trimws(labels))]
    make.unique(trimws(labels), sep = " ")
}

gui_analysis_result_plot <- function(result, method, output_dir, root) {
    plot_data <- result
    colour_name <- NULL
    if ("pseudotime" %in% names(plot_data)) colour_name <- "pseudotime"
    if ("cluster_id" %in% names(plot_data)) {
        plot_data$cluster <- factor(plot_data$cluster_id)
        colour_name <- "cluster"
    }
    mapping <- if (identical(colour_name, "cluster")) {
        ggplot2::aes(x = dimension_1, y = dimension_2, colour = cluster)
    } else if (identical(colour_name, "pseudotime")) {
        ggplot2::aes(x = dimension_1, y = dimension_2, colour = pseudotime)
    } else {
        ggplot2::aes(x = dimension_1, y = dimension_2)
    }
    plot <- ggplot2::ggplot(plot_data, mapping) +
        ggplot2::labs(
            title = method$name,
            subtitle = paste(format(nrow(plot_data), big.mark = ","), "events"),
            x = "Dimension 1", y = "Dimension 2"
        ) +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(face = "bold")
        )
    if (is.null(colour_name)) plot <- plot + ggplot2::geom_point(colour = "#247f9e", size = 0.45, alpha = 0.72) else plot <- plot + ggplot2::geom_point(size = 0.45, alpha = 0.72)
    if (identical(colour_name, "cluster")) plot <- plot + ggplot2::labs(colour = "Cluster")
    if (identical(colour_name, "pseudotime")) plot <- plot + ggplot2::labs(colour = "Pseudotime")
    if (identical(colour_name, "pseudotime")) plot <- plot + ggplot2::scale_colour_gradientn(colours = c("#315a86", "#19a899", "#dfc34f", "#d94b43"))

    png_path <- file.path(output_dir, "plot.png")
    pdf_path <- file.path(output_dir, "plot.pdf")
    svg_path <- file.path(output_dir, "plot.svg")
    png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
    ggplot2::ggsave(png_path, plot, width = 5.4, height = 5.4, dpi = 180, device = png_device, bg = "white")
    ggplot2::ggsave(pdf_path, plot, width = 5.4, height = 5.4, device = grDevices::pdf, bg = "white")
    if (requireNamespace("svglite", quietly = TRUE)) {
        ggplot2::ggsave(svg_path, plot, width = 5.4, height = 5.4, device = svglite::svglite, bg = "white")
    }
    html_path <- file.path(output_dir, "plot.html")
    graphic <- if (file.exists(svg_path)) paste(readLines(svg_path, warn = FALSE), collapse = "\n") else paste0("<img src=\"plot.png\" alt=\"", method$name, " plot\">")
    writeLines(c(
        "<!doctype html><html><head><meta charset=\"utf-8\">",
        paste0("<title>Spectreasy - ", method$name, "</title>"),
        "<style>body{margin:0;padding:24px;background:#f5f7f8;color:#17252d;font-family:system-ui,sans-serif}main{max-width:1100px;margin:auto;background:white;padding:20px;border:1px solid #ccd6db}svg,img{width:100%;height:auto}</style></head><body><main>",
        graphic,
        paste0("<p>Generated by Spectreasy v2. Method: ", method$citation, ".</p>"),
        "</main></body></html>"
    ), html_path, useBytes = TRUE)
    paths <- c(png_path, pdf_path, if (file.exists(svg_path)) svg_path, html_path)
    lapply(paths, function(path) list(
        format = sub("^.*[.]", "", path),
        path = gui_analysis_relative_path(path, root),
        sha256 = gui_analysis_file_sha256(path)
    ))
}

gui_analysis_manifest_id <- function(prefix, value) {
    json <- jsonlite::toJSON(value, auto_unbox = TRUE, null = "null", digits = 15)
    hash <- as.character(openssl::sha256(charToRaw(enc2utf8(json))))
    paste0(prefix, "-", substr(hash, 1L, 16L))
}

gui_analysis_object_directory <- function(root, kind, id) {
    path <- file.path(root, "spectreasy_outputs", "analysis", "objects", kind, id)
    if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
        stop("Could not create the analysis object directory.", call. = FALSE)
    }
    path
}

gui_analysis_artifact_file <- function(path, root) {
    list(
        path = gui_analysis_relative_path(path, root),
        sha256 = gui_analysis_file_sha256(path)
    )
}

gui_analysis_fit_clustering <- function(
    root, method, matrix, indices, input_id, seed,
    requested_neighbors, requested_clusters, advanced_parameters = NULL
) {
    if (is.null(method)) return(list(cluster = NULL, artifact = NULL, stage = NULL))
    parameters <- advanced_parameters %||% gui_analysis_method_parameters(method)
    parameters <- if (identical(method$id, "flowsom")) {
        default_side <- min(10L, max(2L, ceiling(sqrt(nrow(matrix) / 10))))
        parameters$clusters <- min(as.integer(parameters$clusters %||% requested_clusters), max(2L, floor(sqrt(nrow(matrix) / 2))))
        parameters$xdim <- as.integer(parameters$xdim %||% default_side)
        parameters$ydim <- as.integer(parameters$ydim %||% default_side)
        parameters$rlen <- as.integer(parameters$rlen %||% 10L)
        parameters$alpha_start <- as.numeric(parameters$alpha_start %||% 0.05)
        parameters$alpha_end <- as.numeric(parameters$alpha_end %||% 0.01)
        if (parameters$alpha_end > parameters$alpha_start) {
            stop("FlowSOM final learning rate cannot exceed its initial learning rate.", call. = FALSE)
        }
        parameters
    } else {
        parameters$neighbors <- min(as.integer(parameters$neighbors %||% requested_neighbors), nrow(matrix) - 1L)
        parameters
    }
    artifact_id <- gui_analysis_manifest_id("cluster", list(
        input_id = input_id, method = method$id, version = method$version,
        parameters = parameters, seed = seed
    ))
    directory <- gui_analysis_object_directory(root, "clustering", artifact_id)
    object_path <- file.path(directory, "object.rds")
    labels_path <- file.path(directory, "labels.csv")
    reused <- FALSE
    cluster <- NULL
    if (file.exists(object_path) && file.exists(labels_path)) {
        cached_fit <- tryCatch(readRDS(object_path), error = function(e) NULL)
        labels <- tryCatch(data.table::fread(labels_path), error = function(e) NULL)
        if (
            !is.null(cached_fit) &&
            is.data.frame(labels) &&
            identical(as.integer(labels$event_id), as.integer(indices)) &&
            length(labels$cluster_id) == length(indices)
        ) {
            cluster <- as.integer(labels$cluster_id)
            reused <- TRUE
        }
    }
    if (is.null(cluster)) {
        if (identical(method$id, "flowsom")) {
            fit <- FlowSOM::FlowSOM(
                matrix, transform = FALSE, scale = FALSE,
                colsToUse = seq_len(ncol(matrix)),
                xdim = parameters$xdim, ydim = parameters$ydim,
                nClus = min(parameters$clusters, parameters$xdim^2),
                rlen = parameters$rlen,
                alpha = c(parameters$alpha_start, parameters$alpha_end),
                seed = seed, silent = TRUE
            )
            cluster <- as.integer(fit$metaclustering[fit$map$mapping[, 1]])
        } else if (identical(method$id, "phenograph")) {
            fit <- withr::with_seed(
                seed,
                Rphenograph::Rphenograph(matrix, k = parameters$neighbors)
            )
            cluster <- as.integer(igraph::membership(fit[[2]]))
        } else {
            stop("Unsupported clustering method: ", method$id, call. = FALSE)
        }
        saveRDS(fit, object_path, compress = FALSE)
        data.table::fwrite(data.frame(event_id = indices, cluster_id = cluster), labels_path)
    }
    list(
        cluster = cluster,
        artifact = list(
            id = artifact_id,
            type = "clustering",
            method = method,
            reused = reused,
            parameters = parameters,
            object = gui_analysis_artifact_file(object_path, root),
            values = gui_analysis_artifact_file(labels_path, root)
        ),
        stage = list(
            id = method$id,
            name = method$name,
            inputs = "asinh_transformed_matrix",
            outputs = c("clustering_object", "cluster_labels"),
            artifact_id = artifact_id,
            reused = reused
        )
    )
}

gui_analysis_fit_reduction <- function(
    root, method, matrix, indices, input_id, seed,
    output_dimensions, requested_neighbors, requested_perplexity,
    advanced_parameters = NULL
) {
    parameters <- advanced_parameters %||% gui_analysis_method_parameters(method)
    if (identical(method$id, "tsne")) {
        parameters$perplexity <- min(
            as.numeric(parameters$perplexity %||% requested_perplexity),
            max(2, floor((nrow(matrix) - 1) / 3))
        )
    }
    if (method$id %in% c("umap", "diffusion-map", "phate")) {
        parameters$neighbors <- min(as.integer(parameters$neighbors %||% requested_neighbors), nrow(matrix) - 1L)
    }
    if (identical(method$id, "hsne")) {
        parameters$scales <- max(2L, min(as.integer(parameters$scales), 5L))
    }
    if (identical(method$id, "hsne")) output_dimensions <- 2L
    artifact_id <- gui_analysis_manifest_id("embedding", list(
        input_id = input_id, method = method$id, version = method$version,
        dimensions = output_dimensions,
        parameters = parameters, seed = seed
    ))
    directory <- gui_analysis_object_directory(root, "embeddings", artifact_id)
    object_path <- file.path(directory, "object.rds")
    coordinates_path <- file.path(directory, "coordinates.csv")
    reused <- FALSE
    coordinates <- NULL
    fit <- NULL
    if (file.exists(object_path) && file.exists(coordinates_path)) {
        fit <- tryCatch(readRDS(object_path), error = function(e) NULL)
        values <- tryCatch(data.table::fread(coordinates_path), error = function(e) NULL)
        coordinate_names <- if (is.data.frame(values)) grep("^dimension_[0-9]+$", names(values), value = TRUE) else character()
        if (
            !is.null(fit) &&
            is.data.frame(values) &&
            identical(as.integer(values$event_id), as.integer(indices)) &&
            length(coordinate_names) >= 2L
        ) {
            coordinates <- as.matrix(as.data.frame(values)[, coordinate_names, drop = FALSE])
            reused <- TRUE
        }
    }
    if (is.null(coordinates)) {
        native_model_path <- NULL
        if (identical(method$id, "pca")) {
            fit <- stats::prcomp(
                matrix,
                center = isTRUE(parameters$center),
                scale. = isTRUE(parameters$scale),
                rank. = output_dimensions
            )
            coordinates <- fit$x[, seq_len(min(output_dimensions, ncol(fit$x))), drop = FALSE]
        } else if (identical(method$id, "tsne")) {
            fit <- withr::with_seed(
                seed,
                Rtsne::Rtsne(
                    matrix, dims = output_dimensions,
                    perplexity = parameters$perplexity,
                    theta = parameters$theta,
                    max_iter = parameters$iterations,
                    eta = parameters$learning_rate,
                    exaggeration_factor = parameters$exaggeration,
                    normalize = isTRUE(parameters$normalize),
                    check_duplicates = FALSE, pca = TRUE, verbose = FALSE,
                    num_threads = 1L
                )
            )
            coordinates <- fit$Y
        } else if (identical(method$id, "umap")) {
            fit <- uwot::umap(
                matrix,
                n_neighbors = parameters$neighbors,
                n_components = output_dimensions,
                metric = parameters$metric,
                n_epochs = parameters$epochs,
                learning_rate = parameters$learning_rate,
                spread = parameters$spread,
                min_dist = parameters$min_dist,
                repulsion_strength = parameters$repulsion,
                negative_sample_rate = parameters$negative_samples,
                seed = seed, n_threads = 1L,
                ret_model = TRUE, verbose = FALSE
            )
            coordinates <- fit$embedding
        } else if (identical(method$id, "diffusion-map")) {
            fit <- withr::with_seed(
                seed,
                destiny::DiffusionMap(
                    matrix,
                    k = parameters$neighbors,
                    n_eigs = min(parameters$eigenvectors, nrow(matrix) - 2L),
                    distance = parameters$distance,
                    density_norm = isTRUE(parameters$density_normalization),
                    verbose = FALSE
                )
            )
            eigenvectors <- destiny::eigenvectors(fit)
            coordinates <- eigenvectors[, seq_len(min(output_dimensions, ncol(eigenvectors))), drop = FALSE]
        } else if (method$id %in% c("phate", "hsne")) {
            input_path <- tempfile("spectreasy-python-matrix-", tmpdir = directory, fileext = ".csv")
            output_path <- tempfile("spectreasy-python-coordinates-", tmpdir = directory, fileext = ".csv")
            parameters_path <- tempfile("spectreasy-python-parameters-", tmpdir = directory, fileext = ".json")
            on.exit(unlink(c(input_path, output_path, parameters_path), force = TRUE), add = TRUE)
            data.table::fwrite(as.data.frame(matrix), input_path, col.names = FALSE)
            jsonlite::write_json(parameters, parameters_path, auto_unbox = TRUE, null = "null")
            native_model_path <- file.path(directory, if (identical(method$id, "hsne")) "hierarchy.hsne" else "model.pkl")
            gui_analysis_python_run(c(
                "reduce",
                "--method", method$id,
                "--input", input_path,
                "--output", output_path,
                "--model", native_model_path,
                "--parameters", parameters_path,
                "--dimensions", as.character(output_dimensions),
                "--neighbors", as.character(parameters$neighbors %||% requested_neighbors),
                "--seed", as.character(seed)
            ))
            values <- data.table::fread(output_path)
            coordinate_names <- grep("^dimension_[0-9]+$", names(values), value = TRUE)
            if (length(coordinate_names) < 2L || nrow(values) != nrow(matrix)) {
                stop(method$name, " returned an invalid event-aligned embedding.", call. = FALSE)
            }
            coordinates <- as.matrix(as.data.frame(values)[, coordinate_names, drop = FALSE])
            fit <- list(
                method = method$id,
                runtime = "python",
                version = method$version,
                parameters = parameters,
                native_model = basename(native_model_path)
            )
        } else {
            stop("Unsupported dimensional-reduction method: ", method$id, call. = FALSE)
        }
        if (ncol(coordinates) < 2L) coordinates <- cbind(coordinates[, 1], 0)
        coordinate_table <- as.data.frame(coordinates, stringsAsFactors = FALSE)
        names(coordinate_table) <- paste0("dimension_", seq_len(ncol(coordinate_table)))
        saveRDS(fit, object_path, compress = FALSE)
        data.table::fwrite(cbind(data.frame(event_id = indices), coordinate_table), coordinates_path)
    }
    list(
        coordinates = coordinates,
        fit = fit,
        artifact = list(
            id = artifact_id,
            type = "embedding",
            method = method,
            reused = reused,
            parameters = parameters,
            object = gui_analysis_artifact_file(object_path, root),
            values = gui_analysis_artifact_file(coordinates_path, root)
        ),
        stage = list(
            id = method$id,
            name = method$name,
            inputs = "asinh_transformed_matrix",
            outputs = c("embedding_object", "embedding_coordinates"),
            artifact_id = artifact_id,
            reused = reused
        )
    )
}

gui_analysis_normalize_pseudotime <- function(value) {
    value <- as.numeric(value)
    finite <- is.finite(value)
    if (!any(finite)) stop("The trajectory method did not return finite pseudotime values.", call. = FALSE)
    if (!all(finite)) value[!finite] <- stats::median(value[finite])
    limits <- range(value)
    if (!is.finite(diff(limits)) || diff(limits) <= sqrt(.Machine$double.eps)) return(rep(0, length(value)))
    (value - limits[[1]]) / diff(limits)
}

gui_analysis_fit_trajectory <- function(
    root, method, matrix, indices, input_id, seed, root_event_id,
    clustering = NULL, reduction = NULL, requested_neighbors = 15L,
    output_dimensions = 3L, advanced_parameters = NULL
) {
    root_index <- match(root_event_id, indices)
    if (!is.finite(root_index)) stop(method$name, " requires a root event in the analyzed events.", call. = FALSE)
    cluster_artifact_id <- clustering$artifact$id %||% ""
    reduction_artifact_id <- reduction$artifact$id %||% ""
    parameters <- advanced_parameters %||% gui_analysis_method_parameters(method)
    if ("neighbors" %in% names(parameters)) {
        parameters$neighbors <- min(as.integer(parameters$neighbors), nrow(matrix) - 1L)
    }
    if ("candidate_neighbors" %in% names(parameters)) {
        parameters$candidate_neighbors <- min(
            max(as.integer(parameters$candidate_neighbors), as.integer(parameters$neighbors) + 1L),
            nrow(matrix) - 1L
        )
    }
    parameters$dimensions <- output_dimensions
    parameters$root_event_id <- root_event_id
    artifact_id <- gui_analysis_manifest_id("trajectory", list(
        input_id = input_id, method = method$id, version = method$version,
        cluster_artifact_id = cluster_artifact_id,
        reduction_artifact_id = reduction_artifact_id,
        parameters = parameters, seed = seed
    ))
    directory <- gui_analysis_object_directory(root, "trajectories", artifact_id)
    object_path <- file.path(directory, "object.rds")
    values_path <- file.path(directory, "values.csv")
    reused <- FALSE
    fit <- NULL
    values <- NULL
    if (file.exists(object_path) && file.exists(values_path)) {
        fit <- tryCatch(readRDS(object_path), error = function(e) NULL)
        cached <- tryCatch(data.table::fread(values_path), error = function(e) NULL)
        if (
            !is.null(fit) && is.data.frame(cached) &&
            identical(as.integer(cached$event_id), as.integer(indices)) &&
            "pseudotime" %in% names(cached)
        ) {
            values <- cached
            reused <- TRUE
        }
    }

    if (is.null(values)) {
        if (method$id %in% c("slingshot", "tscan")) {
            if (is.null(clustering$cluster) || is.null(reduction$coordinates)) {
                stop(method$name, " requires cached clustering and dimensional-reduction objects.", call. = FALSE)
            }
            clusters <- factor(clustering$cluster)
            start_cluster <- as.character(clusters[[root_index]])
            coordinates <- reduction$coordinates
            rownames(coordinates) <- as.character(indices)
            if (identical(method$id, "slingshot")) {
                fit <- withr::with_seed(seed, slingshot::slingshot(
                    coordinates,
                    clusterLabels = clusters,
                    start.clus = start_cluster,
                    dist.method = parameters$distance,
                    shrink = parameters$shrink,
                    extend = parameters$extend,
                    reweight = isTRUE(parameters$reweight),
                    reassign = isTRUE(parameters$reassign),
                    maxit = parameters$iterations,
                    stretch = parameters$stretch,
                    approx_points = parameters$approximation_points,
                    smoother = parameters$smoother,
                    allow.breaks = isTRUE(parameters$allow_breaks)
                ))
                pseudotime_matrix <- slingshot::slingPseudotime(fit, na = FALSE)
                pseudotime <- rowMeans(pseudotime_matrix, na.rm = TRUE)
                lineage_weights <- slingshot::slingCurveWeights(fit)
                branch <- if (is.null(dim(lineage_weights))) {
                    rep(1L, nrow(coordinates))
                } else {
                    max.col(lineage_weights, ties.method = "first")
                }
            } else {
                mst <- TSCAN::createClusterMST(coordinates, clusters = clusters)
                mapping <- TSCAN::mapCellsToEdges(coordinates, mst, clusters = clusters)
                ordering <- TSCAN::orderCells(mapping, mst, start = start_cluster)
                pseudotime <- TrajectoryUtils::averagePseudotime(ordering)
                fit <- list(mst = mst, mapping = mapping, ordering = ordering)
                branch <- rep(1L, nrow(coordinates))
            }
            coordinate_table <- as.data.frame(coordinates, stringsAsFactors = FALSE)
            names(coordinate_table) <- paste0("dimension_", seq_len(ncol(coordinate_table)))
            values <- cbind(
                data.frame(event_id = indices, stringsAsFactors = FALSE),
                coordinate_table,
                pseudotime = gui_analysis_normalize_pseudotime(pseudotime),
                trajectory_branch = as.integer(branch)
            )
        } else if (method$id %in% c("palantir", "paga-dpt", "wanderlust", "wishbone")) {
            if (identical(method$id, "paga-dpt") && is.null(clustering$cluster)) {
                stop("PAGA + DPT requires a cached clustering object.", call. = FALSE)
            }
            input_path <- tempfile("spectreasy-python-matrix-", tmpdir = directory, fileext = ".csv")
            output_path <- tempfile("spectreasy-python-trajectory-", tmpdir = directory, fileext = ".csv")
            cluster_path <- tempfile("spectreasy-python-clusters-", tmpdir = directory, fileext = ".csv")
            parameters_path <- tempfile("spectreasy-python-parameters-", tmpdir = directory, fileext = ".json")
            on.exit(unlink(c(input_path, output_path, cluster_path, parameters_path), force = TRUE), add = TRUE)
            data.table::fwrite(as.data.frame(matrix), input_path, col.names = FALSE)
            jsonlite::write_json(parameters, parameters_path, auto_unbox = TRUE, null = "null")
            native_model_path <- file.path(
                directory,
                if (identical(method$id, "paga-dpt")) "model.h5ad" else "model.pkl"
            )
            arguments <- c(
                "trajectory",
                "--method", method$id,
                "--input", input_path,
                "--output", output_path,
                "--model", native_model_path,
                "--parameters", parameters_path,
                "--dimensions", as.character(output_dimensions),
                "--neighbors", as.character(parameters$neighbors %||% requested_neighbors),
                "--root-index", as.character(root_index - 1L),
                "--seed", as.character(seed)
            )
            if (identical(method$id, "paga-dpt")) {
                data.table::fwrite(data.frame(cluster = clustering$cluster), cluster_path, col.names = FALSE)
                arguments <- c(arguments, "--clusters", cluster_path)
            }
            gui_analysis_python_run(arguments)
            python_values <- data.table::fread(output_path)
            coordinate_names <- grep("^dimension_[0-9]+$", names(python_values), value = TRUE)
            if (
                nrow(python_values) != nrow(matrix) ||
                length(coordinate_names) < 2L ||
                !"pseudotime" %in% names(python_values)
            ) {
                stop(method$name, " returned an invalid event-aligned trajectory.", call. = FALSE)
            }
            python_values$pseudotime <- gui_analysis_normalize_pseudotime(python_values$pseudotime)
            values <- cbind(data.frame(event_id = indices), python_values)
            fit <- list(
                method = method$id,
                runtime = "python",
                version = method$version,
                parameters = parameters,
                native_model = basename(native_model_path)
            )
        } else {
            stop("Unsupported trajectory method: ", method$id, call. = FALSE)
        }
        saveRDS(fit, object_path, compress = FALSE)
        data.table::fwrite(values, values_path)
    }

    coordinate_names <- grep("^dimension_[0-9]+$", names(values), value = TRUE)
    list(
        coordinates = as.matrix(as.data.frame(values)[, coordinate_names, drop = FALSE]),
        pseudotime = as.numeric(values$pseudotime),
        branch = if ("trajectory_branch" %in% names(values)) as.integer(values$trajectory_branch) else NULL,
        artifact = list(
            id = artifact_id,
            type = "trajectory",
            method = method,
            reused = reused,
            parameters = parameters,
            object = gui_analysis_artifact_file(object_path, root),
            values = gui_analysis_artifact_file(values_path, root)
        ),
        stage = list(
            id = method$id,
            name = method$name,
            inputs = unique(c(
                "root_event",
                if (nzchar(cluster_artifact_id)) "clustering_object" else character(),
                if (nzchar(reduction_artifact_id)) "embedding_object" else "asinh_transformed_matrix"
            )),
            outputs = method$outputs,
            artifact_id = artifact_id,
            reused = reused
        )
    )
}

gui_analysis_annotation_directory <- function(root, analysis_id, create = FALSE) {
    analysis_id <- trimws(as.character(analysis_id %||% "")[1])
    if (!nzchar(analysis_id) || !grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", analysis_id)) {
        stop("A valid analysis result is required for cell annotation.", call. = FALSE)
    }
    path <- file.path(gui_request_project_root(root), "spectreasy_outputs", "analysis", analysis_id)
    if (!dir.exists(path)) stop("The selected analysis result no longer exists.", call. = FALSE)
    if (isTRUE(create)) {
        path <- file.path(path, "identities")
        if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
            stop("Could not create the cell-identity output directory.", call. = FALSE)
        }
    }
    path
}

gui_analysis_identity_signatures <- function(value, available_markers) {
    signatures <- value %||% list()
    if (!is.list(signatures) || length(signatures) < 2L) {
        stop("Define at least two cell identities before annotation.", call. = FALSE)
    }
    if (length(signatures) > 30L) stop("At most 30 cell identities can be scored at once.", call. = FALSE)
    fallback_colors <- c("#197783", "#d06d32", "#7d58a6", "#4f8f45", "#bc4d65", "#416fae", "#9c7427", "#607d8b")
    normalized <- lapply(seq_along(signatures), function(index) {
        signature <- signatures[[index]]
        name <- trimws(as.character(signature$name %||% "")[1])
        if (!nzchar(name)) stop("Every cell identity needs a name.", call. = FALSE)
        positive <- unique(as.character(unlist(signature$positive_markers %||% signature$positiveMarkers %||% character(), use.names = FALSE)))
        negative <- unique(as.character(unlist(signature$negative_markers %||% signature$negativeMarkers %||% character(), use.names = FALSE)))
        positive <- positive[nzchar(positive)]
        negative <- negative[nzchar(negative)]
        if (!length(c(positive, negative))) {
            stop("Cell identity '", name, "' needs at least one positive or negative marker.", call. = FALSE)
        }
        overlap <- intersect(positive, negative)
        if (length(overlap)) {
            stop("Marker ", overlap[[1]], " cannot be both positive and negative for ", name, ".", call. = FALSE)
        }
        missing <- setdiff(c(positive, negative), available_markers)
        if (length(missing)) {
            stop("Cell identity '", name, "' uses unavailable marker ", missing[[1]], ".", call. = FALSE)
        }
        color <- as.character(signature$color %||% fallback_colors[[((index - 1L) %% length(fallback_colors)) + 1L]])[1]
        if (!grepl("^#[0-9A-Fa-f]{6}$", color)) color <- fallback_colors[[((index - 1L) %% length(fallback_colors)) + 1L]]
        list(name = name, color = color, positive_markers = positive, negative_markers = negative)
    })
    names <- vapply(normalized, `[[`, character(1), "name")
    if (anyDuplicated(tolower(names))) stop("Cell identity names must be unique.", call. = FALSE)
    normalized
}

gui_analysis_robust_scale <- function(matrix) {
    centers <- apply(matrix, 2L, stats::median, na.rm = TRUE)
    scales <- apply(matrix, 2L, stats::mad, constant = 1.4826, na.rm = TRUE)
    fallback <- apply(matrix, 2L, stats::IQR, na.rm = TRUE) / 1.349
    invalid <- !is.finite(scales) | scales <= sqrt(.Machine$double.eps)
    scales[invalid] <- fallback[invalid]
    scales[!is.finite(scales) | scales <= sqrt(.Machine$double.eps)] <- 1
    scaled <- sweep(sweep(matrix, 2L, centers, "-"), 2L, scales, "/")
    scaled[!is.finite(scaled)] <- 0
    scaled <- pmax(-8, pmin(8, scaled))
    list(values = scaled, centers = centers, scales = scales)
}

gui_analysis_score_identities <- function(
    matrix, signatures, min_score = 0.55, min_margin = 0.08,
    evidence_sensitivity = 1
) {
    matrix <- as.matrix(matrix)
    storage.mode(matrix) <- "double"
    if (!nrow(matrix) || !ncol(matrix)) stop("No marker values are available for annotation.", call. = FALSE)
    signatures <- gui_analysis_identity_signatures(signatures, colnames(matrix))
    min_score <- as.numeric(min_score)[1]
    min_margin <- as.numeric(min_margin)[1]
    evidence_sensitivity <- as.numeric(evidence_sensitivity)[1]
    if (!is.finite(min_score) || min_score < 0 || min_score > 1) stop("Minimum match score must be between 0 and 1.", call. = FALSE)
    if (!is.finite(min_margin) || min_margin < 0 || min_margin > 1) stop("Minimum score margin must be between 0 and 1.", call. = FALSE)
    if (!is.finite(evidence_sensitivity) || evidence_sensitivity < 0.25 || evidence_sensitivity > 4) {
        stop("Evidence sensitivity must be between 0.25 and 4.", call. = FALSE)
    }

    scaled <- gui_analysis_robust_scale(matrix)
    positive_evidence <- matrix(
        stats::plogis(scaled$values * evidence_sensitivity),
        nrow = nrow(matrix),
        ncol = ncol(matrix),
        dimnames = dimnames(matrix)
    )
    scores <- vapply(signatures, function(signature) {
        evidence <- matrix(numeric(), nrow = nrow(matrix), ncol = 0L)
        if (length(signature$positive_markers)) {
            evidence <- cbind(evidence, positive_evidence[, signature$positive_markers, drop = FALSE])
        }
        if (length(signature$negative_markers)) {
            evidence <- cbind(evidence, 1 - positive_evidence[, signature$negative_markers, drop = FALSE])
        }
        rowMeans(evidence)
    }, numeric(nrow(matrix)))
    if (is.null(dim(scores))) scores <- matrix(scores, ncol = length(signatures))
    colnames(scores) <- vapply(signatures, `[[`, character(1), "name")
    ordering <- t(apply(scores, 1L, order, decreasing = TRUE))
    best_index <- ordering[, 1L]
    runner_index <- ordering[, 2L]
    row_index <- seq_len(nrow(scores))
    best_score <- scores[cbind(row_index, best_index)]
    runner_score <- scores[cbind(row_index, runner_index)]
    margin <- best_score - runner_score
    assigned <- best_score >= min_score & margin >= min_margin
    labels <- colnames(scores)[best_index]
    labels[!assigned] <- "Unassigned"
    list(
        labels = labels,
        best_score = best_score,
        margin = margin,
        scores = scores,
        signatures = signatures,
        scaling = list(
            centers = stats::setNames(as.numeric(scaled$centers), colnames(matrix)),
            robust_sd = stats::setNames(as.numeric(scaled$scales), colnames(matrix))
        ),
        thresholds = list(
            min_score = min_score,
            min_margin = min_margin,
            evidence_sensitivity = evidence_sensitivity
        )
    )
}

gui_analysis_annotate_result <- function(root, body) {
    root <- gui_request_project_root(root)
    analysis_id <- gui_workflow_value(body, "analysisId", "")
    analysis_dir <- gui_analysis_annotation_directory(root, analysis_id)
    metadata_path <- file.path(analysis_dir, "metadata.json")
    events_path <- file.path(analysis_dir, "events.csv")
    if (!file.exists(metadata_path) || !file.exists(events_path)) {
        stop("The selected analysis result is incomplete.", call. = FALSE)
    }
    metadata <- jsonlite::fromJSON(metadata_path, simplifyVector = FALSE)
    marker_columns <- metadata$marker_columns %||% list()
    markers <- vapply(marker_columns, function(item) as.character(item$marker %||% "")[1], character(1))
    columns <- vapply(marker_columns, function(item) as.character(item$column %||% "")[1], character(1))
    if (!length(markers) || any(!nzchar(markers)) || any(!nzchar(columns))) {
        stop("This result predates marker-aware annotation. Rerun the dimensional-reduction pipeline.", call. = FALSE)
    }
    events <- data.table::fread(events_path, data.table = FALSE)
    if (!all(c("event_id", columns) %in% names(events))) {
        stop("The analysis event table is missing marker values.", call. = FALSE)
    }
    marker_matrix <- as.matrix(events[, columns, drop = FALSE])
    colnames(marker_matrix) <- markers
    scored <- gui_analysis_score_identities(
        marker_matrix,
        body$signatures,
        min_score = gui_workflow_number(body, "minScore", 0.55, minimum = 0, maximum = 1),
        min_margin = gui_workflow_number(body, "minMargin", 0.08, minimum = 0, maximum = 1),
        evidence_sensitivity = gui_workflow_number(body, "evidenceSensitivity", 1, minimum = 0.25, maximum = 4)
    )
    signature_payload <- lapply(scored$signatures, function(signature) {
        list(
            name = signature$name,
            color = signature$color,
            positive_markers = signature$positive_markers,
            negative_markers = signature$negative_markers
        )
    })
    identity_id <- gui_analysis_manifest_id("identity", list(
        analysis_id = analysis_id,
        signatures = signature_payload,
        thresholds = scored$thresholds,
        scaling = scored$scaling
    ))
    identity_root <- gui_analysis_annotation_directory(root, analysis_id, create = TRUE)
    output_dir <- file.path(identity_root, identity_id)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    score_columns <- paste0("identity_score_", seq_along(scored$signatures))
    score_table <- data.frame(
        event_id = events$event_id,
        predicted_identity = scored$labels,
        identity_score = scored$best_score,
        identity_margin = scored$margin,
        stringsAsFactors = FALSE
    )
    score_table[score_columns] <- as.data.frame(scored$scores, stringsAsFactors = FALSE)
    scores_path <- file.path(output_dir, "scores.csv")
    model_path <- file.path(output_dir, "model.rds")
    annotation_metadata_path <- file.path(output_dir, "metadata.json")
    data.table::fwrite(score_table, scores_path)
    saveRDS(list(
        method = "robust-signed-marker-score",
        signatures = signature_payload,
        marker_columns = marker_columns,
        scaling = scored$scaling,
        thresholds = scored$thresholds,
        score_columns = stats::setNames(score_columns, vapply(scored$signatures, `[[`, character(1), "name"))
    ), model_path)
    counts <- as.list(table(factor(scored$labels, levels = c(vapply(scored$signatures, `[[`, character(1), "name"), "Unassigned"))))
    names(counts) <- c(vapply(scored$signatures, `[[`, character(1), "name"), "Unassigned")
    annotation_metadata <- list(
        identity_id = identity_id,
        analysis_id = analysis_id,
        method = "robust-signed-marker-score",
        signatures = signature_payload,
        thresholds = scored$thresholds,
        counts = counts,
        event_count = nrow(score_table),
        assigned_count = sum(scored$labels != "Unassigned"),
        unassigned_count = sum(scored$labels == "Unassigned"),
        model = gui_analysis_artifact_file(model_path, root),
        scores = gui_analysis_artifact_file(scores_path, root),
        created_at = as.character(Sys.time())
    )
    jsonlite::write_json(annotation_metadata, annotation_metadata_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
    annotation_metadata$metadata <- gui_analysis_artifact_file(annotation_metadata_path, root)
    list(metadata = annotation_metadata, events = score_table[, 1:4, drop = FALSE])
}

gui_analysis_find_cluster_markers <- function(root, body) {
    root <- gui_request_project_root(root)
    analysis_id <- as.character(body$analysisId %||% body$analysis_id %||% "")[1]
    analysis_dir <- gui_analysis_annotation_directory(root, analysis_id)
    metadata_path <- file.path(analysis_dir, "metadata.json")
    events_path <- file.path(analysis_dir, "events.csv")
    metadata <- jsonlite::fromJSON(metadata_path, simplifyVector = FALSE)
    events <- data.table::fread(events_path, data.table = FALSE)
    if (!"cluster_id" %in% names(events)) {
        stop("Cluster marker discovery requires a completed clustering step.", call. = FALSE)
    }
    marker_columns <- metadata$marker_columns %||% list()
    markers <- vapply(marker_columns, function(item) as.character(item$marker %||% "")[1], character(1))
    channels <- vapply(marker_columns, function(item) as.character(item$channel %||% "")[1], character(1))
    columns <- vapply(marker_columns, function(item) as.character(item$column %||% "")[1], character(1))
    if (!length(columns) || any(!columns %in% names(events))) {
        stop("The analysis event table is missing marker values.", call. = FALSE)
    }
    top_n <- gui_workflow_number(body, "topN", 10L, integer = TRUE, minimum = 1L, maximum = 100L)
    minimum_auc <- gui_workflow_number(body, "minimumAuc", 0.55, minimum = 0.5, maximum = 1)
    clusters <- sort(unique(events$cluster_id[!is.na(events$cluster_id)]))
    if (length(clusters) < 2L) stop("At least two clusters are required to discover markers.", call. = FALSE)
    marker_matrix <- as.matrix(events[, columns, drop = FALSE])
    storage.mode(marker_matrix) <- "double"
    rows <- lapply(clusters, function(cluster) {
        inside <- events$cluster_id == cluster
        n_inside <- sum(inside)
        n_outside <- sum(!inside)
        if (!n_inside || !n_outside) return(NULL)
        pieces <- lapply(seq_along(columns), function(index) {
            positive <- marker_matrix[inside, index]
            negative <- marker_matrix[!inside, index]
            pooled <- c(positive, negative)
            ranks <- rank(pooled, ties.method = "average")
            rank_sum <- sum(ranks[seq_along(positive)])
            auc <- (rank_sum - n_inside * (n_inside + 1) / 2) / (n_inside * n_outside)
            standard_error <- sqrt(n_inside * n_outside * (n_inside + n_outside + 1) / 12)
            z <- if (standard_error > 0) (rank_sum - n_inside * (n_inside + n_outside + 1) / 2) / standard_error else 0
            threshold <- stats::median(pooled, na.rm = TRUE)
            scale <- stats::mad(pooled, constant = 1.4826, na.rm = TRUE)
            if (!is.finite(scale) || scale <= .Machine$double.eps) scale <- stats::sd(pooled, na.rm = TRUE)
            if (!is.finite(scale) || scale <= .Machine$double.eps) scale <- 1
            data.frame(
                cluster_id = cluster,
                marker = markers[[index]],
                channel = channels[[index]],
                median_cluster = stats::median(positive, na.rm = TRUE),
                median_other = stats::median(negative, na.rm = TRUE),
                robust_effect = (stats::median(positive, na.rm = TRUE) - stats::median(negative, na.rm = TRUE)) / scale,
                auc = auc,
                pct_cluster_above_global_median = mean(positive > threshold),
                pct_other_above_global_median = mean(negative > threshold),
                p_value = 2 * stats::pnorm(-abs(z)),
                stringsAsFactors = FALSE
            )
        })
        do.call(rbind, pieces)
    })
    table <- do.call(rbind, rows)
    table$p_adjusted <- stats::p.adjust(table$p_value, method = "BH")
    table <- table[order(table$cluster_id, -table$auc, -table$robust_effect, table$marker), , drop = FALSE]
    table$rank <- ave(seq_len(nrow(table)), table$cluster_id, FUN = seq_along)
    selected <- table[table$auc >= minimum_auc & table$rank <= top_n, , drop = FALSE]
    output_dir <- file.path(analysis_dir, "cluster-markers")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    full_path <- file.path(output_dir, "all-markers.csv")
    selected_path <- file.path(output_dir, "top-markers.csv")
    metadata_output_path <- file.path(output_dir, "metadata.json")
    data.table::fwrite(table, full_path)
    data.table::fwrite(selected, selected_path)
    output <- list(
        analysis_id = analysis_id,
        method = "cluster-vs-rest-rank-auc",
        cluster_count = length(clusters),
        marker_count = length(columns),
        top_n = top_n,
        minimum_auc = minimum_auc,
        full_table = gui_analysis_artifact_file(full_path, root),
        selected_table = gui_analysis_artifact_file(selected_path, root),
        created_at = as.character(Sys.time())
    )
    jsonlite::write_json(output, metadata_output_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
    output$metadata <- gui_analysis_artifact_file(metadata_output_path, root)
    list(metadata = output, markers = selected)
}

gui_analysis_job_directory <- function(root, job_id = NULL, create = FALSE) {
    directory <- file.path(gui_analysis_workspace_directory(root, create = create), "jobs")
    if (isTRUE(create)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    if (is.null(job_id)) return(directory)
    job_id <- as.character(job_id)[1]
    if (!grepl("^[a-z0-9-]{6,80}$", job_id)) stop("Invalid analysis job identifier.", call. = FALSE)
    file.path(directory, job_id)
}

gui_analysis_start_job <- function(root, body) {
    root <- gui_request_project_root(root)
    if (!requireNamespace("processx", quietly = TRUE)) {
        stop("Background analysis requires the suggested processx package.", call. = FALSE)
    }
    worker <- gui_analysis_worker_script()
    if (!nzchar(worker)) stop("The background analysis worker is unavailable.", call. = FALSE)
    job_id <- paste0(
        format(Sys.time(), "%Y%m%d-%H%M%S"), "-",
        substr(gui_analysis_stable_id("job", Sys.getpid(), runif(1), as.numeric(Sys.time())), 5L, 12L)
    )
    directory <- gui_analysis_job_directory(root, job_id, create = TRUE)
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    request_path <- file.path(directory, "request.json")
    result_path <- file.path(directory, "result.rds")
    status_path <- file.path(directory, "status.json")
    stdout_path <- file.path(directory, "stdout.log")
    stderr_path <- file.path(directory, "stderr.log")
    worker_repository <- normalizePath(file.path(dirname(worker), "..", ".."), mustWork = FALSE)
    repository <- if (file.exists(file.path(worker_repository, "DESCRIPTION"))) worker_repository else if (file.exists(file.path(getwd(), "DESCRIPTION"))) normalizePath(getwd(), mustWork = TRUE) else ""
    jsonlite::write_json(
        list(project = root, repository = repository, request = body),
        request_path, auto_unbox = TRUE, pretty = TRUE, null = "null"
    )
    jsonlite::write_json(
        list(state = "queued", message = "Queued", updated_at = as.character(Sys.time()), error = NULL),
        status_path, auto_unbox = TRUE, pretty = TRUE, null = "null"
    )
    process <- processx::process$new(
        command = file.path(R.home("bin"), "Rscript"),
        args = c(worker, request_path, result_path, status_path),
        stdout = stdout_path,
        stderr = stderr_path,
        cleanup = FALSE
    )
    assign(job_id, list(process = process, root = root, directory = directory), envir = .gui_analysis_jobs)
    list(job_id = job_id, state = "queued", message = "Queued")
}

gui_analysis_job_status <- function(root, job_id) {
    root <- gui_request_project_root(root)
    directory <- gui_analysis_job_directory(root, job_id)
    status_path <- file.path(directory, "status.json")
    if (!file.exists(status_path)) stop("Analysis job not found.", call. = FALSE)
    status <- jsonlite::fromJSON(status_path, simplifyVector = FALSE)
    entry <- if (exists(job_id, envir = .gui_analysis_jobs, inherits = FALSE)) get(job_id, envir = .gui_analysis_jobs) else NULL
    if (is.list(entry) && identical(status$state, "queued") && entry$process$is_alive()) {
        status$state <- "running"
        status$message <- "Starting analysis"
    }
    status$job_id <- job_id
    if (identical(status$state, "completed")) {
        result_path <- file.path(directory, "result.rds")
        if (!file.exists(result_path)) stop("Analysis job completed without a result.", call. = FALSE)
        status$result <- readRDS(result_path)
    }
    status
}

gui_analysis_cancel_job <- function(root, job_id) {
    root <- gui_request_project_root(root)
    directory <- gui_analysis_job_directory(root, job_id)
    status_path <- file.path(directory, "status.json")
    if (!file.exists(status_path)) stop("Analysis job not found.", call. = FALSE)
    entry <- if (exists(job_id, envir = .gui_analysis_jobs, inherits = FALSE)) get(job_id, envir = .gui_analysis_jobs) else NULL
    if (is.list(entry) && entry$process$is_alive()) entry$process$kill()
    status <- list(
        job_id = job_id,
        state = "cancelled",
        message = "Cancelled by user",
        updated_at = as.character(Sys.time()),
        error = NULL
    )
    jsonlite::write_json(status, status_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
    status
}

gui_analysis_run_method <- function(root, body) {
    root <- gui_request_project_root(root)
    workspace <- gui_analysis_read_workspace(root)
    requested_method_id <- gui_workflow_value(body, "method", "")
    cluster_method_id <- gui_workflow_value(body, "clusterMethod", "")
    reduction_method_id <- gui_workflow_value(body, "reductionMethod", "")
    requested_method <- if (nzchar(requested_method_id)) gui_analysis_method(requested_method_id) else NULL
    trajectory_request <- is.list(requested_method) && identical(requested_method$family, "trajectory")
    pipeline_request <- !trajectory_request && (nzchar(cluster_method_id) || nzchar(reduction_method_id))

    cluster_method <- NULL
    if ((pipeline_request || trajectory_request) && nzchar(cluster_method_id) && !identical(cluster_method_id, "none")) {
        cluster_method <- gui_analysis_method(cluster_method_id)
        if (!identical(cluster_method$family, "clustering")) {
            stop(cluster_method$name, " is not a clustering method.", call. = FALSE)
        }
        if (!isTRUE(cluster_method$available)) {
            requirement <- paste(as.character(unlist(cluster_method$requirements, use.names = FALSE)), collapse = "; ")
            stop(cluster_method$name, " is not executable in this runtime. ", requirement, call. = FALSE)
        }
    }

    if (trajectory_request) {
        method <- requested_method
        method_id <- method$id
        needs_clustering <- "clustering" %in% method$prerequisites
        needs_reduction <- "reduction" %in% method$prerequisites
        if (needs_clustering && is.null(cluster_method)) {
            stop(method$name, " requires a clustering method.", call. = FALSE)
        }
        if (needs_reduction && !nzchar(reduction_method_id)) {
            stop(method$name, " requires a dimensional-reduction method.", call. = FALSE)
        }
    } else if (pipeline_request) {
        if (!nzchar(reduction_method_id)) stop("Select a dimensional-reduction method.", call. = FALSE)
        method <- gui_analysis_method(reduction_method_id)
        if (!identical(method$family, "reduction")) {
            stop(method$name, " is not a dimensional-reduction method.", call. = FALSE)
        }
        method_id <- method$id
    } else {
        method_id <- if (nzchar(requested_method_id)) requested_method_id else "pca"
        method <- requested_method %||% gui_analysis_method(method_id)
        if (identical(method$family, "clustering")) {
            stop("Choose a dimensional-reduction method to visualize ", method$name, " clusters.", call. = FALSE)
        }
    }
    if (!isTRUE(method$available)) {
        requirement <- paste(as.character(unlist(method$requirements, use.names = FALSE)), collapse = "; ")
        stop(method$name, " is not executable in this runtime. ", requirement, call. = FALSE)
    }

    files <- unique(as.character(unlist(body$files %||% body$file %||% "", use.names = FALSE)))
    files <- files[nzchar(trimws(files))]
    if (!length(files)) stop("Select at least one FCS file.", call. = FALSE)
    file <- files[[1]]
    population_id <- gui_workflow_value(body, "populationId", "root")
    seed <- gui_workflow_number(body, "seed", workspace$seed %||% 20260723L, integer = TRUE, minimum = 1)
    max_events <- gui_workflow_number(body, "maxEvents", 20000L, integer = TRUE, minimum = 100L, maximum = 200000L)
    cofactor <- gui_workflow_number(body, "cofactor", 150, minimum = .Machine$double.eps)
    requested_neighbors <- gui_workflow_number(body, "neighbors", 15L, integer = TRUE, minimum = 2L, maximum = 200L)
    requested_clusters <- gui_workflow_number(body, "clusters", 10L, integer = TRUE, minimum = 2L, maximum = 100L)
    requested_perplexity <- gui_workflow_number(body, "perplexity", 30, minimum = 2, maximum = 200)
    requested_markers <- as.character(unlist(body$markers %||% character(), use.names = FALSE))
    frames <- lapply(files, function(candidate) gui_analysis_flow_frame(root, candidate))
    common_channels <- Reduce(intersect, lapply(frames, function(frame) colnames(flowCore::exprs(frame))))
    missing_markers <- setdiff(requested_markers, common_channels)
    if (length(missing_markers)) {
        stop(
            "Selected markers are not present in every pooled file: ",
            paste(missing_markers, collapse = ", "),
            call. = FALSE
        )
    }
    pooled <- lapply(seq_along(files), function(index) {
        candidate <- files[[index]]
        values <- flowCore::exprs(frames[[index]])[, common_channels, drop = FALSE]
        mask <- gui_analysis_gate_mask(values, workspace, population_id, source_file = candidate)
        selected <- which(mask)
        list(
            data = values[selected, , drop = FALSE],
            source_file = rep(candidate, length(selected)),
            source_event_id = selected
        )
    })
    data <- do.call(rbind, lapply(pooled, `[[`, "data"))
    source_files <- unlist(lapply(pooled, `[[`, "source_file"), use.names = FALSE)
    source_event_ids <- as.integer(unlist(lapply(pooled, `[[`, "source_event_id"), use.names = FALSE))
    if (!nrow(data)) stop("The selected population has no events in the selected files.", call. = FALSE)
    sampling_identity <- paste(files, collapse = "|")
    sampling_identity <- paste(sampling_identity, population_id, "analysis-input", paste(sort(requested_markers), collapse = ","), cofactor)
    indices <- gui_analysis_sample_indices(seq_len(nrow(data)), max_events, seed, sampling_identity)
    root_event_id <- NULL
    root_source_event_id <- NULL
    root_source_file <- NULL
    if ("trajectory-root" %in% method$user_prerequisites) {
        root_source_event_id <- suppressWarnings(as.integer(gui_workflow_value(body, "rootEventId", workspace$root_event_id %||% NA_integer_)))
        root_source_file <- gui_workflow_value(body, "rootSourceFile", workspace$root_source_file %||% file)
        root_matches <- which(source_files == root_source_file & source_event_ids == root_source_event_id)
        root_event_id <- if (length(root_matches)) root_matches[[1]] else NA_integer_
        if (!is.finite(root_event_id)) {
            stop(method$name, " requires a root event inside the selected population.", call. = FALSE)
        }
        if (!root_event_id %in% indices) {
            indices <- sort(c(if (length(indices) >= max_events) head(indices, -1L) else indices, root_event_id))
        }
    }
    if (length(indices) < 20L) stop("The selected population has too few events for this method.", call. = FALSE)
    matrix <- gui_analysis_numeric_matrix(data[indices, , drop = FALSE], requested_markers, cofactor = cofactor)
    marker_labels <- gui_analysis_marker_display_labels(frames[[1]], colnames(matrix))
    output_dimensions <- min(3L, ncol(matrix), nrow(matrix) - 1L)
    source_manifest <- lapply(files, function(candidate) {
        source_path <- gui_analysis_resolve_fcs(root, candidate)
        source_info <- file.info(source_path)
        list(
            source_file = candidate,
            source_size = unname(source_info$size),
            source_modified = as.numeric(source_info$mtime)
        )
    })
    input_id <- gui_analysis_manifest_id("input", list(
        sources = source_manifest,
        population_id = population_id,
        populations = workspace$populations,
        markers = colnames(matrix),
        transform = list(name = "asinh", cofactor = cofactor),
        max_events = max_events,
        seed = seed,
        events = lapply(indices, function(index) list(
            source_file = source_files[[index]],
            source_event_id = source_event_ids[[index]]
        ))
    ))
    started <- proc.time()[[3]]
    pseudotime <- NULL
    intermediate_tables <- list()
    pipeline <- list(list(
        id = "input", name = "Gated event preparation",
        inputs = c("selected_population", "selected_markers"),
        outputs = "asinh_transformed_matrix"
    ))
    method_parameters <- gui_analysis_request_method_parameters(body, method)
    cluster_parameters <- if (!is.null(cluster_method)) {
        gui_analysis_request_method_parameters(body, cluster_method)
    } else {
        NULL
    }

    clustering <- gui_analysis_fit_clustering(
        root, cluster_method, matrix, indices, input_id, seed,
        requested_neighbors, requested_clusters,
        advanced_parameters = cluster_parameters
    )
    if (!is.null(clustering$stage)) pipeline <- c(pipeline, list(clustering$stage))

    trajectory <- NULL
    if (identical(method_id, "dpt")) {
        reduction_method <- gui_analysis_method("diffusion-map")
        diffusion_parameters <- list(
            neighbors = method_parameters$neighbors,
            eigenvectors = method_parameters$eigenvectors,
            distance = method_parameters$distance,
            density_normalization = method_parameters$density_normalization
        )
        reduction <- gui_analysis_fit_reduction(
            root, reduction_method, matrix, indices, input_id, seed,
            output_dimensions, requested_neighbors, requested_perplexity,
            advanced_parameters = diffusion_parameters
        )
        coordinates <- reduction$coordinates
        root_index <- match(root_event_id, indices)
        # DPT1 is the accumulated diffusion distance from the selected tip.
        # Constructing the full DPT branching object invokes an unrelated
        # automatic branch heuristic that can fail on valid linear or compact
        # populations. Directly evaluating the public DPT distance accessor is
        # equivalent to DPT(..., tips = root)$DPT1 and avoids that failure mode.
        ordering <- methods::new(
            "DPT",
            branch = matrix(numeric(), nrow = 0L, ncol = 0L),
            tips = matrix(numeric(), nrow = 0L, ncol = 0L),
            dm = reduction$fit
        )
        pseudotime <- gui_analysis_normalize_pseudotime(ordering[root_index, seq_len(nrow(matrix))])
        diffusion_table <- as.data.frame(coordinates, stringsAsFactors = FALSE)
        names(diffusion_table) <- paste0("diffusion_", seq_len(ncol(diffusion_table)))
        intermediate_tables$diffusion_map <- cbind(
            data.frame(event_id = indices, stringsAsFactors = FALSE),
            diffusion_table
        )
        pipeline <- c(pipeline, list(
            reduction$stage,
            list(
                id = "dpt", name = "Diffusion pseudotime",
                inputs = c("embedding_object", "root_event"),
                outputs = "pseudotime"
            )
        ))
    } else if (trajectory_request) {
        if (method_id %in% c("slingshot", "tscan")) {
            reduction_method <- gui_analysis_method(reduction_method_id)
            if (!identical(reduction_method$family, "reduction") || !isTRUE(reduction_method$available)) {
                stop(method$name, " requires an enabled dimensional-reduction method.", call. = FALSE)
            }
            reduction_parameters <- gui_analysis_request_method_parameters(body, reduction_method)
            reduction <- gui_analysis_fit_reduction(
                root, reduction_method, matrix, indices, input_id, seed,
                output_dimensions, requested_neighbors, requested_perplexity,
                advanced_parameters = reduction_parameters
            )
            pipeline <- c(pipeline, list(reduction$stage))
        } else {
            reduction_method <- if (identical(method_id, "palantir")) gui_analysis_method("diffusion-map") else NULL
            reduction <- list(artifact = NULL)
        }
        trajectory <- gui_analysis_fit_trajectory(
            root, method, matrix, indices, input_id, seed, root_event_id,
            clustering = clustering,
            reduction = reduction,
            requested_neighbors = requested_neighbors,
            output_dimensions = output_dimensions,
            advanced_parameters = method_parameters
        )
        coordinates <- trajectory$coordinates
        pseudotime <- trajectory$pseudotime
        pipeline <- c(pipeline, list(trajectory$stage))
    } else {
        reduction_method <- method
        reduction <- gui_analysis_fit_reduction(
            root, reduction_method, matrix, indices, input_id, seed,
            output_dimensions, requested_neighbors, requested_perplexity,
            advanced_parameters = method_parameters
        )
        coordinates <- reduction$coordinates
        pipeline <- c(pipeline, list(reduction$stage))
    }

    runtime <- proc.time()[[3]] - started
    coordinate_table <- as.data.frame(coordinates, stringsAsFactors = FALSE)
    names(coordinate_table) <- paste0("dimension_", seq_len(ncol(coordinate_table)))
    result <- cbind(
        data.frame(
            event_id = indices,
            source_file = source_files[indices],
            source_event_id = source_event_ids[indices],
            sample_id = match(source_files[indices], files),
            stringsAsFactors = FALSE
        ),
        coordinate_table
    )
    marker_columns <- paste0("marker_", seq_len(ncol(matrix)))
    marker_table <- as.data.frame(matrix, stringsAsFactors = FALSE)
    names(marker_table) <- marker_columns
    result <- cbind(result, marker_table)
    if (!is.null(clustering$cluster)) result$cluster_id <- clustering$cluster
    if (!is.null(pseudotime)) result$pseudotime <- pseudotime
    if (!is.null(trajectory$branch)) result$trajectory_branch <- trajectory$branch
    rownames(result) <- NULL

    display_name <- if (trajectory_request || identical(method_id, "dpt")) {
        method$name
    } else if (!is.null(cluster_method)) {
        paste(cluster_method$name, "\u2192", reduction_method$name)
    } else {
        method$name
    }
    analysis_slug <- if (!trajectory_request && !identical(method_id, "dpt") && !is.null(cluster_method)) {
        paste(cluster_method$id, reduction_method$id, sep = "-")
    } else {
        method_id
    }
    analysis_id <- paste0(
        format(Sys.time(), "%Y%m%d-%H%M%S"), "-", analysis_slug, "-",
        substr(gui_analysis_stable_id("run", file, population_id, seed, as.numeric(Sys.time())), 5L, 12L)
    )
    output_dir <- file.path(root, "spectreasy_outputs", "analysis", analysis_id)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    intermediate_files <- lapply(names(intermediate_tables), function(name) {
        path <- file.path(output_dir, paste0(name, ".csv"))
        data.table::fwrite(intermediate_tables[[name]], path)
        list(id = name, path = gui_analysis_relative_path(path, root), sha256 = gui_analysis_file_sha256(path))
    })
    csv_path <- file.path(output_dir, "events.csv")
    metadata_path <- file.path(output_dir, "metadata.json")
    data.table::fwrite(result, csv_path)
    plot_method <- method
    plot_method$name <- gsub("\u2192", "+", display_name, fixed = TRUE)
    plot_files <- gui_analysis_result_plot(result, plot_method, output_dir, root)
    coordinate_basis_id <- if (
        trajectory_request &&
        method_id %in% c("palantir", "paga-dpt", "wanderlust", "wishbone")
    ) method_id else reduction_method$id
    coordinate_labels <- switch(
        coordinate_basis_id,
        pca = paste("PC", seq_len(ncol(coordinate_table))),
        tsne = paste("t-SNE", seq_len(ncol(coordinate_table))),
        umap = paste("UMAP", seq_len(ncol(coordinate_table))),
        `diffusion-map` = paste("DC", seq_len(ncol(coordinate_table))),
        phate = paste("PHATE", seq_len(ncol(coordinate_table))),
        hsne = paste("HSNE", seq_len(ncol(coordinate_table))),
        palantir = paste("Palantir DC", seq_len(ncol(coordinate_table))),
        `paga-dpt` = paste("PAGA DC", seq_len(ncol(coordinate_table))),
        wanderlust = paste("Wanderlust graph", seq_len(ncol(coordinate_table))),
        wishbone = paste("Wishbone DC", seq_len(ncol(coordinate_table))),
        paste("Dimension", seq_len(ncol(coordinate_table)))
    )
    metadata <- list(
        analysis_id = analysis_id,
        display_name = display_name,
        method = method,
        cluster_method = cluster_method,
        reduction_method = reduction_method,
        source_file = file,
        source_files = files,
        population_id = population_id,
        population_path = gui_analysis_population_path(workspace, population_id),
        markers = marker_labels,
        marker_channels = colnames(matrix),
        marker_columns = lapply(seq_along(marker_columns), function(index) list(
            marker = marker_labels[[index]],
            channel = colnames(matrix)[[index]],
            column = marker_columns[[index]]
        )),
        coordinate_count = ncol(coordinate_table),
        coordinate_labels = coordinate_labels,
        transform = list(name = "asinh", cofactor = cofactor),
        input_id = input_id,
        pipeline = pipeline,
        artifacts = list(
            clustering = clustering$artifact,
            embedding = reduction$artifact,
            trajectory = trajectory$artifact
        ),
        intermediate_files = intermediate_files,
        parameters = list(
            max_events = max_events,
            cofactor = cofactor,
            method = method_parameters,
            clustering = cluster_parameters,
            reduction = if (exists("reduction_parameters", inherits = FALSE)) reduction_parameters else NULL
        ),
        seed = seed,
        event_count = nrow(result),
        root_event_id = root_source_event_id,
        root_source_file = root_source_file,
        runtime_seconds = unname(runtime),
        package_version = method$version,
        plot_files = plot_files,
        created_at = as.character(Sys.time()),
        events_file = gui_analysis_relative_path(csv_path, root)
    )
    jsonlite::write_json(metadata, metadata_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
    list(metadata = metadata, events = result, metadata_file = gui_analysis_relative_path(metadata_path, root))
}
