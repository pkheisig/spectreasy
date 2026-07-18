# Project-scoped benchmark configuration and sequential workflow orchestration.

.benchmark_id_pattern <- "^[A-Za-z0-9][A-Za-z0-9._-]*$"

.benchmark_generated_id <- function(prefix) {
    token <- paste0(
        sprintf("%08x", as.integer(Sys.time()) %% .Machine$integer.max),
        paste(sample(c(0:9, letters[1:6]), 24L, replace = TRUE), collapse = "")
    )
    paste0(prefix, "-", substr(token, 1, 8), "-", substr(token, 9, 12), "-", substr(token, 13, 16), "-", substr(token, 17, 20), "-", substr(token, 21, 32))
}

.validate_benchmark_id <- function(x, label = "benchmark_id") {
    x <- trimws(as.character(x)[1])
    if (!nzchar(x) || !grepl(.benchmark_id_pattern, x)) {
        stop(label, " must contain only letters, numbers, '.', '_', and '-'.", call. = FALSE)
    }
    x
}

.benchmark_root <- function(project_path, benchmark_id = NULL, create = FALSE) {
    project_path <- normalizePath(project_path, mustWork = TRUE)
    base <- file.path(project_path, "spectreasy_outputs", "benchmarks")
    path <- if (is.null(benchmark_id)) base else file.path(base, .validate_benchmark_id(benchmark_id))
    if (isTRUE(create)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    path
}

.benchmark_status <- function(state = "not_run", message = "", updated_at = NULL) {
    list(state = state, message = message, updated_at = updated_at)
}

.benchmark_empty_control_outputs <- function() {
    list(matrix_file = "", detector_noise_file = "", variant_library_file = "", report_file = "", output_dir = "")
}

.benchmark_empty_sample_outputs <- function() {
    list(report_file = "", output_dir = "")
}

.normalize_benchmark_entry <- function(entry, index = 1L) {
    entry_id <- entry$entry_id %||% entry$entryId %||% .benchmark_generated_id("entry")
    entry_id <- .validate_benchmark_id(entry_id, "entry_id")
    method <- .normalize_unmix_method(entry$method %||% "AutoSpectral")
    normalize_status <- function(x) {
        x <- x %||% list()
        state <- as.character(x$state %||% "not_run")[1]
        allowed <- c("not_run", "queued", "running", "complete", "failed", "stale")
        if (!state %in% allowed) state <- "not_run"
        .benchmark_status(state, as.character(x$message %||% "")[1], x$updated_at %||% NULL)
    }
    list(
        entry_id = entry_id,
        label = as.character(entry$label %||% paste("Entry", index))[1],
        enabled = if (is.null(entry$enabled)) TRUE else isTRUE(entry$enabled),
        method = method,
        control_settings = entry$control_settings %||% entry$controlSettings %||% list(),
        sample_settings = entry$sample_settings %||% entry$sampleSettings %||% list(),
        control_status = normalize_status(entry$control_status %||% entry$controlStatus),
        sample_status = normalize_status(entry$sample_status %||% entry$sampleStatus),
        control_outputs = utils::modifyList(.benchmark_empty_control_outputs(), entry$control_outputs %||% entry$controlOutputs %||% list()),
        sample_outputs = utils::modifyList(.benchmark_empty_sample_outputs(), entry$sample_outputs %||% entry$sampleOutputs %||% list())
    )
}

.new_benchmark_config <- function(name = "Method comparison", benchmark_id = .benchmark_generated_id("benchmark")) {
    now <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    list(
        schema_version = 1L,
        benchmark_id = .validate_benchmark_id(benchmark_id),
        name = trimws(as.character(name)[1]),
        created_at = now,
        updated_at = now,
        entries = list()
    )
}

validate_benchmark_config <- function(config, require_runnable = FALSE) {
    if (!is.list(config)) stop("Benchmark configuration must be a list.", call. = FALSE)
    config$schema_version <- as.integer(config$schema_version %||% config$schemaVersion %||% 1L)
    if (!identical(config$schema_version, 1L)) stop("Unsupported benchmark schema version.", call. = FALSE)
    config$benchmark_id <- .validate_benchmark_id(config$benchmark_id %||% config$benchmarkId)
    config$name <- trimws(as.character(config$name %||% "Method comparison")[1])
    if (!nzchar(config$name)) stop("Benchmark name must not be empty.", call. = FALSE)
    entries <- config$entries %||% list()
    config$entries <- lapply(seq_along(entries), function(index) .normalize_benchmark_entry(entries[[index]], index))
    ids <- vapply(config$entries, `[[`, character(1), "entry_id")
    if (anyDuplicated(ids)) stop("Benchmark entry IDs must be unique.", call. = FALSE)
    if (isTRUE(require_runnable)) {
        enabled <- vapply(config$entries, `[[`, logical(1), "enabled")
        if (sum(enabled) < 2L) stop("At least two enabled benchmark entries are required.", call. = FALSE)
    }
    config$created_at <- as.character(config$created_at %||% config$createdAt %||% format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))[1]
    config$updated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    config
}

save_benchmark_config <- function(config, project_path) {
    config <- validate_benchmark_config(config)
    root <- .benchmark_root(project_path, config$benchmark_id, create = TRUE)
    path <- file.path(root, "benchmark_config.json")
    temporary <- tempfile("benchmark_config_", tmpdir = root, fileext = ".json")
    on.exit(unlink(temporary, force = TRUE), add = TRUE)
    jsonlite::write_json(config, temporary, auto_unbox = TRUE, pretty = TRUE, null = "null")
    parsed <- jsonlite::fromJSON(temporary, simplifyVector = FALSE)
    validate_benchmark_config(parsed)
    if (!file.rename(temporary, path) && !file.copy(temporary, path, overwrite = TRUE)) {
        stop("Could not save benchmark configuration.", call. = FALSE)
    }
    invisible(normalizePath(path, mustWork = TRUE))
}

load_benchmark_config <- function(benchmark_id, project_path) {
    path <- file.path(.benchmark_root(project_path, benchmark_id), "benchmark_config.json")
    if (!file.exists(path)) stop("Benchmark configuration not found: ", benchmark_id, call. = FALSE)
    validate_benchmark_config(jsonlite::fromJSON(path, simplifyVector = FALSE))
}

list_benchmark_configs <- function(project_path) {
    root <- .benchmark_root(project_path)
    if (!dir.exists(root)) return(list())
    paths <- list.files(root, pattern = "^benchmark_config\\.json$", recursive = TRUE, full.names = TRUE)
    configs <- lapply(paths, function(path) tryCatch(validate_benchmark_config(jsonlite::fromJSON(path, simplifyVector = FALSE)), error = function(e) NULL))
    configs[!vapply(configs, is.null, logical(1))]
}

.benchmark_relative_path <- function(path, project_path) {
    if (is.null(path) || !length(path) || is.na(path[1]) || !nzchar(as.character(path)[1])) return("")
    gui_root <- normalizePath(project_path, mustWork = TRUE)
    normalized <- normalizePath(as.character(path)[1], mustWork = FALSE)
    prefix <- paste0(gui_root, .Platform$file.sep)
    if (!startsWith(normalized, prefix)) stop("Benchmark output escaped the active project.", call. = FALSE)
    gsub("\\\\", "/", substring(normalized, nchar(prefix) + 1L))
}

.benchmark_runner_args <- function(settings, defaults, runner) {
    args <- utils::modifyList(defaults, settings %||% list())
    aliases <- c(method = "unmixing_method", threads = "n_threads", refine_af_bank = "autospectral_refine")
    for (old in names(aliases)) {
        if (!is.null(args[[old]]) && is.null(args[[aliases[[old]]]])) args[[aliases[[old]]]] <- args[[old]]
        args[[old]] <- NULL
    }
    accepted <- names(formals(runner))
    if (!"..." %in% accepted) args <- args[intersect(names(args), accepted)]
    args
}

run_benchmark_controls <- function(benchmark_id, project_path, runner = unmix_controls) {
    config <- validate_benchmark_config(load_benchmark_config(benchmark_id, project_path), require_runnable = TRUE)
    layout <- .spectreasy_project_layout(project_path)
    gate_candidates <- c(file.path(project_path, "ssc_gate_config.csv"), file.path(project_path, ".spectreasy", "gui", "control_gating.csv"))
    gate_file <- gate_candidates[file.exists(gate_candidates)][1]
    if (is.na(gate_file)) gate_file <- ""
    logs <- character()
    response <- list()
    for (index in seq_along(config$entries)) {
        entry <- config$entries[[index]]
        if (!isTRUE(entry$enabled)) next
        entry$control_status <- .benchmark_status("running", "", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
        if (identical(entry$sample_status$state, "complete")) entry$sample_status <- .benchmark_status("stale", "Control output was rerun.", entry$control_status$updated_at)
        config$entries[[index]] <- entry
        save_benchmark_config(config, project_path)
        output_root <- file.path(.benchmark_root(project_path, benchmark_id, create = TRUE), "controls", entry$entry_id)
        defaults <- list(
            scc_dir = file.path(project_path, layout$control_input_dir),
            control_file = file.path(project_path, "fcs_mapping.csv"),
            output_dir = output_root,
            unmixing_method = entry$method,
            gating_mode = if (nzchar(gate_file)) "reuse" else "automatic",
            manual_gate_file = if (nzchar(gate_file)) gate_file else NULL,
            save_report = TRUE,
            save_qc_png = TRUE,
            verbose = FALSE,
            project_path = project_path
        )
        result <- tryCatch(do.call(runner, .benchmark_runner_args(entry$control_settings, defaults, runner)), error = function(e) e)
        now <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
        if (inherits(result, "error")) {
            entry$control_status <- .benchmark_status("failed", conditionMessage(result), now)
            entry$sample_status <- .benchmark_status(
                "failed",
                "Blocked because the linked control run failed.",
                now
            )
            logs <- c(logs, paste0(entry$label, ": ", conditionMessage(result)))
        } else {
            entry$control_status <- .benchmark_status("complete", "", now)
            entry$control_outputs <- list(
                matrix_file = .benchmark_relative_path(result$reference_matrix_file, project_path),
                detector_noise_file = .benchmark_relative_path(result$detector_noise_file, project_path),
                variant_library_file = .benchmark_relative_path(result$spectral_variant_library_file, project_path),
                report_file = .benchmark_relative_path(result$qc_report_file, project_path),
                output_dir = .benchmark_relative_path(file.path(output_root, "unmix_controls"), project_path)
            )
        }
        config$entries[[index]] <- entry
        save_benchmark_config(config, project_path)
        response[[length(response) + 1L]] <- list(entry_id = entry$entry_id, state = entry$control_status$state, message = entry$control_status$message, report_path = entry$control_outputs$report_file, output_paths = entry$control_outputs)
    }
    list(success = all(vapply(response, function(x) identical(x$state, "complete"), logical(1))), message = "Benchmark control queue finished.", entries = response, logs = logs, benchmark = config)
}

run_benchmark_samples <- function(benchmark_id, project_path, runner = unmix_samples) {
    config <- validate_benchmark_config(load_benchmark_config(benchmark_id, project_path), require_runnable = TRUE)
    layout <- .spectreasy_project_layout(project_path)
    logs <- character()
    response <- list()
    for (index in seq_along(config$entries)) {
        entry <- config$entries[[index]]
        if (!isTRUE(entry$enabled)) next
        matrix_file <- file.path(project_path, entry$control_outputs$matrix_file %||% "")
        if (!identical(entry$control_status$state, "complete") || !file.exists(matrix_file)) {
            entry$sample_status <- .benchmark_status("failed", "Linked control output is not available.", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
            config$entries[[index]] <- entry
            save_benchmark_config(config, project_path)
            response[[length(response) + 1L]] <- list(entry_id = entry$entry_id, state = "failed", message = entry$sample_status$message)
            next
        }
        entry$sample_status <- .benchmark_status("running", "", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
        config$entries[[index]] <- entry
        save_benchmark_config(config, project_path)
        output_root <- file.path(.benchmark_root(project_path, benchmark_id, create = TRUE), "samples", entry$entry_id)
        variant_path <- entry$control_outputs$variant_library_file %||% ""
        defaults <- list(
            sample_dir = file.path(project_path, layout$sample_input_dir),
            unmixing_matrix_file = matrix_file,
            detector_noise_file = if (nzchar(entry$control_outputs$detector_noise_file %||% "")) file.path(project_path, entry$control_outputs$detector_noise_file) else NULL,
            spectral_variant_library_file = if (nzchar(variant_path)) file.path(project_path, variant_path) else NULL,
            output_dir = output_root,
            unmixing_method = entry$method,
            write_fcs = TRUE,
            save_report = TRUE,
            save_qc_plots = TRUE,
            verbose = FALSE,
            project_path = project_path
        )
        result <- tryCatch(do.call(runner, .benchmark_runner_args(entry$sample_settings, defaults, runner)), error = function(e) e)
        now <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
        if (inherits(result, "error")) {
            entry$sample_status <- .benchmark_status("failed", conditionMessage(result), now)
            logs <- c(logs, paste0(entry$label, ": ", conditionMessage(result)))
        } else {
            entry$sample_status <- .benchmark_status("complete", "", now)
            entry$sample_outputs <- list(
                report_file = .benchmark_relative_path(attr(result, "qc_report_file"), project_path),
                output_dir = .benchmark_relative_path(file.path(output_root, "unmix_samples"), project_path)
            )
        }
        config$entries[[index]] <- entry
        save_benchmark_config(config, project_path)
        response[[length(response) + 1L]] <- list(entry_id = entry$entry_id, state = entry$sample_status$state, message = entry$sample_status$message, report_path = entry$sample_outputs$report_file, output_paths = entry$sample_outputs)
    }
    list(success = all(vapply(response, function(x) identical(x$state, "complete"), logical(1))), message = "Benchmark sample queue finished.", entries = response, logs = logs, benchmark = config)
}
