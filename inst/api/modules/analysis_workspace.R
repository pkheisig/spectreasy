# Project-scoped post-unmixing analysis helpers. This module deliberately keeps
# editable analysis state outside FCS files: FCS exports contain selected events
# and provenance, while the workspace JSON remains the hierarchy source of truth.

if (!exists("%||%", mode = "function", inherits = TRUE)) {
    `%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x
}

.gui_analysis_cache <- new.env(parent = emptyenv())

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
        schema_version = 1L,
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
        root_event_id = NULL
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
    out$schema_version <- 1L
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

gui_analysis_method_registry <- function() {
    method <- function(id, name, family, package, citation, doi, requirements = character(), research_only = FALSE, adapter_verified = TRUE) {
        package_available <- if (!nzchar(package)) TRUE else requireNamespace(package, quietly = TRUE)
        list(
            id = id, name = name, family = family, package = package,
            available = isTRUE(package_available) && isTRUE(adapter_verified) && !isTRUE(research_only),
            installed = isTRUE(package_available),
            research_only = isTRUE(research_only),
            adapter_verified = isTRUE(adapter_verified),
            version = if (package_available && nzchar(package)) as.character(utils::packageVersion(package)) else "",
            citation = citation, doi = doi, requirements = requirements
        )
    }
    list(
        method("pca", "PCA", "reduction", "", "Pearson (1901), Philosophical Magazine", "10.1080/14786440109462720"),
        method("tsne", "t-SNE", "reduction", "Rtsne", "van der Maaten and Hinton (2008), JMLR", ""),
        method("umap", "UMAP", "reduction", "uwot", "McInnes, Healy and Melville (2018)", "10.48550/arXiv.1802.03426"),
        method("diffusion-map", "Diffusion map", "reduction", "destiny", "Angerer et al. (2016), Bioinformatics", "10.1093/bioinformatics/btv715"),
        method("phate", "PHATE", "reduction", "", "Moon et al. (2019), Nature Biotechnology", "10.1038/s41587-019-0336-3", "Python package phate", TRUE),
        method("hsne", "HSNE", "reduction", "", "Pezzotti et al. (2016), Computer Graphics Forum", "10.1111/cgf.12878", "Validated HSNE adapter", TRUE),
        method("flowsom", "FlowSOM", "clustering", "FlowSOM", "Van Gassen et al. (2015), Cytometry A", "10.1002/cyto.a.22625"),
        method("phenograph", "PhenoGraph", "clustering", "Rphenograph", "Levine et al. (2015), Cell", "10.1016/j.cell.2015.05.047"),
        method("dpt", "Diffusion pseudotime", "trajectory", "destiny", "Haghverdi et al. (2016), Nature Methods", "10.1038/nmeth.3971", "A selected root event"),
        method("slingshot", "Slingshot", "trajectory", "slingshot", "Street et al. (2018), BMC Genomics", "10.1186/s12864-018-4772-0", "Install slingshot, then validate its cluster-to-lineage adapter", adapter_verified = FALSE),
        method("tscan", "TSCAN", "trajectory", "TSCAN", "Ji and Ji (2016), Nucleic Acids Research", "10.1093/nar/gkw430", "Install TSCAN, then validate its cluster-to-lineage adapter", adapter_verified = FALSE),
        method("palantir", "Palantir", "trajectory", "", "Setty et al. (2019), Nature Biotechnology", "10.1038/s41587-019-0068-4", "Validated Python palantir environment", TRUE),
        method("paga-dpt", "PAGA + DPT", "trajectory", "", "Wolf et al. (2019), Genome Biology", "10.1186/s13059-019-1663-x", "Validated Python scanpy environment", TRUE),
        method("wanderlust", "Wanderlust", "trajectory", "", "Bendall et al. (2014), Cell", "10.1016/j.cell.2014.04.005", "Validated compatibility adapter", TRUE),
        method("wishbone", "Wishbone", "trajectory", "", "Setty et al. (2016), Nature Biotechnology", "10.1038/nbt.3569", "Validated compatibility adapter", TRUE)
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
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), plot.title = ggplot2::element_text(face = "bold"))
    if (is.null(colour_name)) plot <- plot + ggplot2::geom_point(colour = "#247f9e", size = 0.45, alpha = 0.72) else plot <- plot + ggplot2::geom_point(size = 0.45, alpha = 0.72)
    if (identical(colour_name, "cluster")) plot <- plot + ggplot2::labs(colour = "Cluster")
    if (identical(colour_name, "pseudotime")) plot <- plot + ggplot2::labs(colour = "Pseudotime")
    if (identical(colour_name, "pseudotime")) plot <- plot + ggplot2::scale_colour_gradientn(colours = c("#315a86", "#19a899", "#dfc34f", "#d94b43"))

    png_path <- file.path(output_dir, "plot.png")
    pdf_path <- file.path(output_dir, "plot.pdf")
    svg_path <- file.path(output_dir, "plot.svg")
    png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
    ggplot2::ggsave(png_path, plot, width = 7.2, height = 5.4, dpi = 180, device = png_device, bg = "white")
    ggplot2::ggsave(pdf_path, plot, width = 7.2, height = 5.4, device = grDevices::pdf, bg = "white")
    if (requireNamespace("svglite", quietly = TRUE)) {
        ggplot2::ggsave(svg_path, plot, width = 7.2, height = 5.4, device = svglite::svglite, bg = "white")
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

gui_analysis_run_method <- function(root, body) {
    root <- gui_request_project_root(root)
    workspace <- gui_analysis_read_workspace(root)
    method_id <- gui_workflow_value(body, "method", "pca")
    method <- gui_analysis_method(method_id)
    if (!isTRUE(method$available)) {
        requirement <- paste(as.character(unlist(method$requirements, use.names = FALSE)), collapse = "; ")
        stop(method$name, " is not executable in this runtime. ", requirement, call. = FALSE)
    }
    file <- gui_workflow_value(body, "file", "")
    population_id <- gui_workflow_value(body, "populationId", "root")
    seed <- gui_workflow_number(body, "seed", workspace$seed %||% 20260723L, integer = TRUE, minimum = 1)
    max_events <- gui_workflow_number(body, "maxEvents", 20000L, integer = TRUE, minimum = 100L, maximum = 200000L)
    cofactor <- gui_workflow_number(body, "cofactor", 150, minimum = .Machine$double.eps)
    requested_neighbors <- gui_workflow_number(body, "neighbors", 15L, integer = TRUE, minimum = 2L, maximum = 200L)
    requested_clusters <- gui_workflow_number(body, "clusters", 10L, integer = TRUE, minimum = 2L, maximum = 100L)
    requested_perplexity <- gui_workflow_number(body, "perplexity", 30, minimum = 2, maximum = 200)
    frame <- gui_analysis_flow_frame(root, file)
    data <- flowCore::exprs(frame)
    mask <- gui_analysis_gate_mask(data, workspace, population_id, source_file = file)
    indices <- gui_analysis_sample_indices(which(mask), max_events, seed, paste(file, population_id, method_id))
    root_event_id <- NULL
    if (identical(method_id, "dpt")) {
        root_event_id <- suppressWarnings(as.integer(gui_workflow_value(body, "rootEventId", workspace$root_event_id %||% NA_integer_)))
        if (!is.finite(root_event_id) || root_event_id < 1L || root_event_id > nrow(data) || !isTRUE(mask[[root_event_id]])) {
            stop("DPT requires a root event inside the selected population.", call. = FALSE)
        }
        if (!root_event_id %in% indices) {
            indices <- sort(c(if (length(indices) >= max_events) head(indices, -1L) else indices, root_event_id))
        }
    }
    if (length(indices) < 20L) stop("The selected population has too few events for this method.", call. = FALSE)
    matrix <- gui_analysis_numeric_matrix(data[indices, , drop = FALSE], body$markers, cofactor = cofactor)
    started <- proc.time()[[3]]
    coordinates <- NULL
    cluster <- NULL
    pseudotime <- NULL
    if (identical(method_id, "pca")) {
        fit <- stats::prcomp(matrix, center = TRUE, scale. = TRUE, rank. = 2L)
        coordinates <- fit$x[, seq_len(min(2L, ncol(fit$x))), drop = FALSE]
    } else if (identical(method_id, "tsne")) {
        perplexity <- min(requested_perplexity, max(2, floor((nrow(matrix) - 1) / 3)))
        fit <- withr::with_seed(seed, Rtsne::Rtsne(matrix, dims = 2L, perplexity = perplexity, check_duplicates = FALSE, pca = TRUE, verbose = FALSE))
        coordinates <- fit$Y
    } else if (identical(method_id, "umap")) {
        neighbors <- min(requested_neighbors, nrow(matrix) - 1L)
        coordinates <- uwot::umap(matrix, n_neighbors = neighbors, n_components = 2L, seed = seed, n_threads = 1L, verbose = FALSE)
    } else if (method_id %in% c("diffusion-map", "dpt")) {
        k <- min(requested_neighbors, nrow(matrix) - 1L)
        fit <- destiny::DiffusionMap(matrix, k = k, n_eigs = min(10L, nrow(matrix) - 2L), verbose = FALSE)
        eigenvectors <- destiny::eigenvectors(fit)
        coordinates <- eigenvectors[, seq_len(min(2L, ncol(eigenvectors))), drop = FALSE]
        if (identical(method_id, "dpt")) {
            root_index <- match(root_event_id, indices)
            ordering <- destiny::DPT(fit, tips = root_index)
            pseudotime <- as.numeric(ordering$DPT1)
        }
    } else if (identical(method_id, "flowsom")) {
        n_clusters <- min(requested_clusters, max(2L, floor(sqrt(nrow(matrix) / 2))))
        fit <- FlowSOM::FlowSOM(matrix, transform = FALSE, scale = FALSE, colsToUse = seq_len(ncol(matrix)), nClus = n_clusters, seed = seed, silent = TRUE)
        cluster <- as.integer(fit$metaclustering[fit$map$mapping[, 1]])
        coordinates <- stats::prcomp(matrix, center = TRUE, scale. = TRUE, rank. = 2L)$x[, 1:2, drop = FALSE]
    } else if (identical(method_id, "phenograph")) {
        k <- min(requested_neighbors, nrow(matrix) - 1L)
        fit <- Rphenograph::Rphenograph(matrix, k = k)
        cluster <- as.integer(igraph::membership(fit[[2]]))
        coordinates <- stats::prcomp(matrix, center = TRUE, scale. = TRUE, rank. = 2L)$x[, 1:2, drop = FALSE]
    } else {
        stop("The installed adapter for ", method$name, " has not passed the v2 execution contract yet.", call. = FALSE)
    }
    runtime <- proc.time()[[3]] - started
    if (ncol(coordinates) < 2L) coordinates <- cbind(coordinates[, 1], 0)
    result <- data.frame(
        event_id = indices,
        dimension_1 = as.numeric(coordinates[, 1]),
        dimension_2 = as.numeric(coordinates[, 2]),
        stringsAsFactors = FALSE
    )
    if (!is.null(cluster)) result$cluster_id <- cluster
    if (!is.null(pseudotime)) result$pseudotime <- pseudotime
    analysis_id <- paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), "-", method_id, "-", substr(gui_analysis_stable_id("run", file, population_id, seed, as.numeric(Sys.time())), 5L, 12L))
    output_dir <- file.path(root, "spectreasy_outputs", "analysis", analysis_id)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    csv_path <- file.path(output_dir, "events.csv")
    metadata_path <- file.path(output_dir, "metadata.json")
    data.table::fwrite(result, csv_path)
    plot_files <- gui_analysis_result_plot(result, method, output_dir, root)
    metadata <- list(
        analysis_id = analysis_id,
        method = method,
        source_file = file,
        population_id = population_id,
        population_path = gui_analysis_population_path(workspace, population_id),
        markers = colnames(matrix),
        transform = list(name = "asinh", cofactor = cofactor),
        parameters = list(
            max_events = max_events,
            neighbors = if (exists("neighbors", inherits = FALSE)) neighbors else if (exists("k", inherits = FALSE)) k else NULL,
            clusters = if (exists("n_clusters", inherits = FALSE)) n_clusters else NULL,
            perplexity = if (exists("perplexity", inherits = FALSE)) perplexity else NULL
        ),
        seed = seed,
        event_count = nrow(result),
        root_event_id = root_event_id,
        runtime_seconds = unname(runtime),
        package_version = method$version,
        plot_files = plot_files,
        created_at = as.character(Sys.time()),
        events_file = gui_analysis_relative_path(csv_path, root)
    )
    jsonlite::write_json(metadata, metadata_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
    list(metadata = metadata, events = result, metadata_file = gui_analysis_relative_path(metadata_path, root))
}
