# Load implementation helpers from small, responsibility-specific modules.
# Route annotations remain in this entry file so plumber can discover them.
.gui_api_module_directory <- function() {
    frame_files <- unlist(lapply(sys.frames(), function(frame) {
        value <- frame$ofile
        if (is.null(value) || !length(value)) character() else as.character(value[[1]])
    }), use.names = FALSE)
    frame_files <- frame_files[nzchar(frame_files)]
    candidates <- c(
        if (length(frame_files)) file.path(dirname(frame_files), "modules") else character(),
        file.path(getwd(), "inst", "api", "modules"),
        system.file("api", "modules", package = "spectreasy")
    )
    candidates <- unique(candidates[nzchar(candidates)])
    required <- "project_matrix.R"
    hit <- candidates[file.exists(file.path(candidates, required))]
    if (!length(hit)) stop("Could not locate Spectreasy GUI API modules.", call. = FALSE)
    normalizePath(hit[[1]], mustWork = TRUE)
}

.gui_api_modules <- c(
    "project_matrix.R",
    "gating_io.R",
    "gating_selection.R",
    "gating_spectrum.R",
    "workflow_project.R",
    "ai_qc.R"
)
.gui_api_module_dir <- .gui_api_module_directory()
for (.gui_api_module in .gui_api_modules) {
    sys.source(file.path(.gui_api_module_dir, .gui_api_module), envir = environment())
}
rm(.gui_api_module, .gui_api_module_dir, .gui_api_modules, .gui_api_module_directory)

#* @filter logger
function(req) {
    if (isTRUE(getOption("spectreasy.gui_request_log", FALSE))) {
        cat(as.character(Sys.time()), "-", req$REQUEST_METHOD, req$PATH_INFO, "\n")
    }
    plumber::forward()
}

gui_request_origin_allowed <- function(req) {
    origin <- req$HTTP_ORIGIN
    if (is.null(origin) || !nzchar(trimws(origin))) return(TRUE)
    allowed <- trimws(as.character(getOption("spectreasy.gui_allowed_origins", character())))
    allowed <- allowed[nzchar(allowed)]
    trimws(origin) %in% allowed
}

gui_api_token_value_allowed <- function(token) {
    expected <- as.character(getOption("spectreasy.gui_api_token", ""))[1]
    if (is.null(token)) token <- ""
    nzchar(expected) && identical(as.character(token)[1], expected)
}

gui_api_token_allowed <- function(req) {
    gui_api_token_value_allowed(req$HTTP_X_SPECTREASY_TOKEN)
}

#* @filter cors
function(req, res) {
    origin <- if (is.null(req$HTTP_ORIGIN)) "" else trimws(req$HTTP_ORIGIN)
    if (!gui_request_origin_allowed(req)) {
        res$status <- 403
        return(list(error = "This GUI origin is not authorized for the local Spectreasy session."))
    }
    mutating_method <- toupper(req$REQUEST_METHOD) %in% c("POST", "PUT", "PATCH", "DELETE")
    supplied_token <- req$HTTP_X_SPECTREASY_TOKEN
    validates_session <- identical(req$PATH_INFO, "/status") &&
        !is.null(supplied_token) && nzchar(trimws(as.character(supplied_token)[1]))
    if ((mutating_method || validates_session) && !gui_api_token_allowed(req)) {
        res$status <- 403
        return(list(error = "This request is not authorized for the active local Spectreasy session."))
    }
    if (nzchar(origin)) {
        res$setHeader("Access-Control-Allow-Origin", origin)
        res$setHeader("Vary", "Origin")
    }
    res$setHeader("Access-Control-Allow-Methods", "GET, POST, DELETE, OPTIONS")
    res$setHeader("Access-Control-Allow-Headers", "Content-Type, X-Spectreasy-Token")
    res$setHeader("Access-Control-Allow-Private-Network", "true")
    res$setHeader("Cache-Control", "no-store, no-cache, must-revalidate, max-age=0")
    res$setHeader("Pragma", "no-cache")
    res$setHeader("Expires", "0")
    if (req$REQUEST_METHOD == "OPTIONS") {
        res$status <- 200
        return(list())
    }
    plumber::forward()
}

#* Health Check
#* @get /status
function() {
    project_selected <- isTRUE(getOption("spectreasy.project_selected", TRUE))
    project_path <- if (project_selected) get_matrix_dir() else ""
    layout <- if (nzchar(project_path) && dir.exists(project_path)) gui_project_layout(project_path, persist = FALSE) else list()
    return(list(
        status = "ok",
        time = Sys.time(),
        wd = getwd(),
        matrix_dir = get_matrix_dir(),
        samples_dir = get_samples_dir(),
        project_layout = layout,
        unmixing_method = get_unmixing_method(),
        gui_mode = getOption("spectreasy.gui_mode", "tuner"),
        panel_cytometer = getOption("spectreasy.panel_cytometer", "aurora"),
        project_selected = project_selected,
        project_name = if (nzchar(project_path)) basename(project_path) else ""
    ))
}

#* Load persistent GUI state for one GUI module
#* @get /gui_state
#* @param module
function(module = "matrix_tuner", project_path = "") {
    path <- user_gui_config_path(module, project_path)
    if (!file.exists(path)) {
        return(list(module = normalize_gui_module(module), path = path, config = list()))
    }
    cfg <- tryCatch(jsonlite::fromJSON(path, simplifyVector = TRUE), error = function(e) list())
    list(module = normalize_gui_module(module), path = path, config = cfg)
}

#* Save persistent GUI state for one GUI module
#* @post /gui_state
function(req) {
    body <- jsonlite::fromJSON(req$postBody, simplifyVector = FALSE)
    module <- if (!is.null(body$module)) body$module else "matrix_tuner"
    cfg <- if (!is.null(body$config_json)) body$config_json else list()
    project_path <- if (!is.null(body$projectPath)) body$projectPath else if (!is.null(body$project_path)) body$project_path else ""
    path <- user_gui_config_path(module, project_path)
    jsonlite::write_json(cfg, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
    list(success = TRUE, module = normalize_gui_module(module), path = path)
}

#* List SCC control files for manual gating
#* @get /gate_files
function(project_path = "") {
    gui_with_project_context(project_path, {
        df <- gate_read_mapping()
        first <- df$filename[df$file_exists][1]
        meta <- if (!is.na(first)) gate_fcs_metadata(file.path(get_gate_scc_dir(), first)) else list()
        list(files = df, metadata = meta, gate_file = get_gate_file())
    })
}

#* Read the saved control mapping without synthesizing placeholder rows
#* @get /control_mapping
function(project_path = "") {
    root <- gui_request_project_root(project_path)
    path <- file.path(root, "fcs_mapping.csv")
    if (!file.exists(path)) return(list(rows = data.frame(), exists = FALSE, path = path))
    rows <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    if (is.null(rows)) return(list(rows = data.frame(), exists = TRUE, path = path, error = "The existing fcs_mapping.csv could not be read."))
    rows <- gui_annotate_active_af_mapping(rows, root)
    list(rows = rows, exists = TRUE, path = path)
}

#* Create fcs_mapping.csv from the active project's configured control input folder
#* @post /control_mapping/create
function(req) {
    body <- gui_workflow_body(req)
    tryCatch({
        cytometer <- getOption("spectreasy.panel_cytometer", "auto")
        if (is.null(cytometer) || length(cytometer) == 0L || is.na(cytometer[1]) || !nzchar(trimws(as.character(cytometer[1])))) cytometer <- "auto"
        root <- gui_project_value(body)
        rows <- spectreasy::create_control_file(
            input_folder = gui_project_input_path(root, "controls"),
            cytometer = cytometer,
            unknown_fluor_policy = "by_channel",
            output_file = file.path(root, "fcs_mapping.csv")
        )
        list(success = TRUE, path = file.path(root, "fcs_mapping.csv"), rows = rows)
    }, error = function(e) list(success = FALSE, error = conditionMessage(e)))
}

#* Save the control mapping edited in the cockpit
#* @post /control_mapping
function(req) {
    body <- jsonlite::fromJSON(req$postBody, simplifyVector = FALSE)
    rows <- body$rows
    if (is.null(rows) || length(rows) == 0) {
        return(list(success = FALSE, error = "No control mapping rows were supplied."))
    }
    row_value <- function(row, keys, fallback = "") {
        for (key in keys) {
            value <- row[[key]]
            if (!is.null(value) && length(value) > 0 && !is.na(value[1]) && nzchar(as.character(value[1]))) {
                return(as.character(value[1]))
            }
        }
        fallback
    }
    mapping <- do.call(rbind, lapply(rows, function(row) {
        data.frame(
            filename = row_value(row, c("file", "filename")),
            fluorophore = row_value(row, "fluorophore"),
            marker = row_value(row, "marker"),
            channel = row_value(row, "channel"),
            control.type = if (grepl("bead", row_value(row, c("controlType", "control.type"), "cell"), ignore.case = TRUE)) "beads" else "cells",
            universal.negative = row_value(row, c("universalNegative", "universal.negative")),
            is.viability = row_value(row, c("isViability", "is.viability"), "FALSE"),
            stringsAsFactors = FALSE
        )
    }))
    root <- gui_project_value(body)
    path <- file.path(root, "fcs_mapping.csv")
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(mapping, path, row.names = FALSE, quote = TRUE)
    list(success = TRUE, path = path, rows = nrow(mapping))
}

#* Load one downsampled SCC control payload
#* @get /gate_events
#* @param filename
#* @param max_points
function(filename, max_points = 3000, project_path = "") {
    gui_with_project_context(project_path, gate_payload_for_file(filename = filename, max_points = max_points))
}

#* Preload all downsampled SCC control payloads
#* @get /gate_preload
#* @param max_points
function(max_points = 3000, project_path = "") {
    gui_with_project_context(project_path, {
    df <- gate_read_mapping()
    payloads <- lapply(df$filename, function(filename) {
        tryCatch(gate_payload_for_file(filename = filename, max_points = max_points), error = function(e) {
            list(error = conditionMessage(e), filename = filename)
        })
    })
    list(payloads = payloads, max_points = as.integer(max_points))
    })
}

#* Preload all SCC controls in a compact float32 representation
#* @get /gate_preload_compact
#* @param max_points
function(max_points = 3000, project_path = "") {
    gui_with_project_context(project_path, {
    df <- gate_read_mapping()
    payloads <- lapply(df$filename, function(filename) {
        tryCatch({
            payload <- gate_payload_for_file(
                filename = filename,
                max_points = max_points,
                cache_result = FALSE
            )
            compact <- gate_compact_payload(payload)
            rm(payload)
            compact
        }, error = function(e) {
            list(error = conditionMessage(e), filename = filename)
        })
    })
    list(payloads = payloads, max_points = as.integer(max_points))
    })
}

#* Auto-generate the required histogram gates for all non-AF controls
#* @post /gate_histogram_autogate
function(req) {
    tryCatch({
        body <- jsonlite::fromJSON(req$postBody, simplifyVector = FALSE)
        result <- gui_with_project_context(gui_project_value(body), gate_autogenerate_histograms(body$gates))
        list(
            success = TRUE,
            gates = result$gates,
            spectra = result$spectra,
            files_processed = result$files_processed,
            files_succeeded = result$files_succeeded,
            files_failed = result$files_failed,
            gates_generated = result$gates_generated,
            gates_preserved = result$gates_preserved,
            failures = result$failures
        )
    }, error = function(e) {
        list(success = FALSE, error = conditionMessage(e))
    })
}

#* Compute selected SCC spectrum bins for the current gate cache
#* @get /gate_spectrum
#* @param filename
#* @param dark
function(filename, dark = "false", project_path = "") {
    tryCatch(
        list(spectrum = gui_with_project_context(project_path, gate_spectrum_for_file(filename = filename))),
        error = function(e) list(error = conditionMessage(e), spectrum = NULL)
    )
}

#* Compute and cache spectrum bins for one or more SCC controls using explicit gates
#* @post /gate_spectra
function(req) {
    tryCatch({
        body <- jsonlite::fromJSON(req$postBody, simplifyVector = FALSE)
        filenames <- unlist(body$filenames, use.names = FALSE)
        filenames <- as.character(filenames[nzchar(as.character(filenames))])
        gates <- body$gates
        spectra <- gui_with_project_context(gui_project_value(body), stats::setNames(lapply(filenames, function(filename) {
            gate_spectrum_for_file(
                filename = filename,
                gates = gates
            )
        }), filenames))
        list(success = TRUE, spectra = spectra)
    }, error = function(e) {
        list(success = FALSE, error = conditionMessage(e), spectra = list())
    })
}

#* List gate CSV files
#* @get /gate_configs
function(project_path = "") {
    gui_with_project_context(project_path, {
        dir.create(dirname(get_gate_file()), recursive = TRUE, showWarnings = FALSE)
        files <- list.files(dirname(get_gate_file()), pattern = "\\.csv$", full.names = FALSE, ignore.case = TRUE)
        list(configs = sort(files), config_dir = dirname(get_gate_file()), active = basename(get_gate_file()))
    })
}

#* Load gate CSV
#* @get /gate_config
#* @param filename
function(filename = "", project_path = "") {
    gui_with_project_context(project_path, {
    path <- if (is.null(filename) || !nzchar(trimws(as.character(filename)[1]))) {
        get_gate_file()
    } else {
        file.path(dirname(get_gate_file()), basename(as.character(filename)[1]))
    }
    if (!file.exists(path)) {
        return(list(path = path, rows = gate_empty_config()))
    }
    rows <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    list(path = path, rows = rows)
    })
}

#* Save gate CSV
#* @post /gate_config
function(req) {
    body <- jsonlite::fromJSON(req$postBody, simplifyVector = TRUE)
    written <- gui_with_project_context(gui_project_value(body), gate_write_config_csv(body$rows, get_gate_file()))
    list(success = TRUE, path = written$path, rows = written$rows)
}

#* Save gate CSV with a system file picker
#* @post /gate_config_save_dialog
function(req) {
    body <- jsonlite::fromJSON(req$postBody, simplifyVector = TRUE)
    root <- gui_project_value(body)
    if (!gate_has_system_file_picker()) {
        return(list(success = FALSE, cancelled = FALSE, fallback = TRUE, message = "No backend system file picker is available."))
    }
    path <- gui_with_project_context(root, gate_pick_csv_file(mode = "save"))
    if (identical(path, "CANCEL")) {
        return(list(success = FALSE, cancelled = TRUE, fallback = FALSE, message = "Save cancelled."))
    }
    if (!nzchar(path)) {
        return(list(success = FALSE, cancelled = FALSE, fallback = TRUE, message = "System file picker failed. Falling back to browser picker."))
    }
    written <- gate_write_config_csv(body$rows, path)
    list(success = TRUE, cancelled = FALSE, path = written$path, rows = written$rows)
}

#* Load gate CSV with a system file picker
#* @get /gate_config_load_dialog
function(project_path = "") {
    if (!gate_has_system_file_picker()) {
        return(list(success = FALSE, cancelled = FALSE, fallback = TRUE, message = "No backend system file picker is available."))
    }
    path <- gui_with_project_context(project_path, gate_pick_csv_file(mode = "open"))
    if (identical(path, "CANCEL")) {
        return(list(success = FALSE, cancelled = TRUE, fallback = FALSE, message = "Load cancelled."))
    }
    if (!nzchar(path)) {
        return(list(success = FALSE, cancelled = FALSE, fallback = TRUE, message = "System file picker failed. Falling back to browser picker."))
    }
    if (!file.exists(path)) {
        return(list(success = FALSE, cancelled = FALSE, fallback = FALSE, message = paste("File does not exist:", path)))
    }
    rows <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    list(success = TRUE, cancelled = FALSE, path = path, rows = rows)
}

#* Load in-memory gate cache
#* @get /gate_cache
function(project_path = "") {
    gui_with_project_context(project_path, {
    cache <- getOption("spectreasy.gating_state_cache")
    if (is.null(cache)) {
        return(list(gates = list(), pointSize = 1.5, maxPoints = 50000, histogramBins = 100, histogramTransform = "auto", histogramTransformVersion = 2, viewSettings = list(), eventCountVersion = 2))
    }
    cache
    })
}

#* Save in-memory gate cache
#* @post /gate_cache
function(req) {
    body <- jsonlite::fromJSON(req$postBody, simplifyVector = FALSE)
    gui_with_project_context(gui_project_value(body), options(spectreasy.gating_state_cache = list(
        gates = body$gates,
        pointSize = body$pointSize,
        maxPoints = body$maxPoints,
        histogramBins = body$histogramBins,
        histogramTransform = body$histogramTransform,
        histogramTransformVersion = body$histogramTransformVersion,
        viewSettings = body$viewSettings,
        eventCountVersion = body$eventCountVersion
    )))
    list(success = TRUE)
}

#* Shut down manual gating GUI
#* @post /gate_shutdown
function(req) {
    options(spectreasy.gui_shutdown_requested = TRUE)
    server <- getOption("spectreasy.gui_server", NULL)
    options(spectreasy.gui_server = NULL)
    later::later(function() {
        if (!is.null(server)) {
            try(httpuv::stopServer(server), silent = TRUE)
        }
    }, delay = 0.1)
    list(success = TRUE, message = "Gate config saved. Manual gating GUI is shutting down.")
}

#* Spectral panel builder metadata and current selection
#* @get /spectral_panel
#* @param cytometer
#* @param configuration
function(cytometer = "", configuration = "") {
    selected_cytometer <- if (is.null(cytometer) || !nzchar(trimws(as.character(cytometer)[1]))) {
        getOption("spectreasy.panel_cytometer", "aurora")
    } else {
        cytometer
    }
    selected_configuration <- if (is.null(configuration) || !nzchar(trimws(as.character(configuration)[1]))) NULL else configuration
    tryCatch(
        spectreasy:::.spectral_panel_payload(
            cytometer = selected_cytometer,
            fluorophores = character(),
            configuration = selected_configuration
        ),
        error = function(e) list(error = conditionMessage(e))
    )
}

#* CORS preflight for spectral_panel_metrics
#* @options /spectral_panel_metrics
function(res) {
    return("")
}

#* Recalculate spectral panel metrics for selected fluorophores
#* @post /spectral_panel_metrics
function(req) {
    body <- jsonlite::fromJSON(req$postBody, simplifyVector = TRUE)
    cytometer <- if (!is.null(body$cytometer)) body$cytometer else getOption("spectreasy.panel_cytometer", "aurora")
    configuration <- if (!is.null(body$configuration)) body$configuration else NULL
    fluorophores <- if (!is.null(body$fluorophores)) body$fluorophores else character()
    tryCatch(
        spectreasy:::.spectral_panel_payload(
            cytometer = cytometer,
            fluorophores = fluorophores,
            configuration = configuration
        ),
        error = function(e) list(error = conditionMessage(e))
    )
}

#* CORS preflight for export_spectral_panel_overview
#* @options /export_spectral_panel_overview
function(res) {
    return("")
}

#* Export spectral panel overview PDF
#* @post /export_spectral_panel_overview
function(req) {
    body <- jsonlite::fromJSON(req$postBody, simplifyVector = TRUE)
    cytometer <- if (!is.null(body$cytometer)) body$cytometer else getOption("spectreasy.panel_cytometer", "aurora")
    configuration <- if (!is.null(body$configuration)) body$configuration else NULL
    fluorophores <- if (!is.null(body$fluorophores)) body$fluorophores else character()
    markers <- if (!is.null(body$markers)) body$markers else character()

    tryCatch({
        output_file <- tempfile("spectreasy_panel_overview_", fileext = ".pdf")
        on.exit(unlink(output_file), add = TRUE)
        spectreasy:::.write_spectral_panel_overview_pdf(
            cytometer = cytometer,
            configuration = configuration,
            fluorophores = fluorophores,
            markers = markers,
            output_file = output_file
        )
        payload <- readBin(output_file, what = "raw", n = file.info(output_file)$size)
        list(
            filename = paste0("spectreasy_", cytometer, "_", ifelse(is.null(configuration), "panel", configuration), "_overview.pdf"),
            content_type = "application/pdf",
            content_base64 = jsonlite::base64_enc(payload)
        )
    }, error = function(e) list(error = conditionMessage(e)))
}

#* List saved AF profiles
#* @get /af_profiles
function(project_path = "") {
    root <- gui_request_project_root(project_path)
    profiles <- tryCatch(spectreasy::list_af_profiles(), error = function(e) data.frame())
    if (nrow(profiles) > 0L) profiles$active <- profiles$name == gui_read_active_af_profile(root)
    list(profiles = profiles)
}

#* Return the detector-wise spectra for one saved AF profile
#* @get /af_profiles/data
#* @param name Profile name
function(name = "") {
    name <- trimws(as.character(name)[1])
    if (!nzchar(name)) return(list(error = "A profile name is required."))
    tryCatch({
        profile <- spectreasy::load_af_profile(name, show_plot = FALSE)$profile
        list(
            name = name,
            detectors = colnames(profile),
            spectra = lapply(seq_len(nrow(profile)), function(index) list(
                name = rownames(profile)[index],
                values = as.numeric(profile[index, , drop = TRUE])
            ))
        )
    }, error = function(e) list(error = conditionMessage(e)))
}

gui_pick_af_source_file <- function(initial_dir = gui_project_input_path(get_matrix_dir(), "controls")) {
    initial_dir <- normalizePath(initial_dir, mustWork = FALSE)
    sysname <- Sys.info()[["sysname"]]
    if (identical(sysname, "Darwin")) {
        script <- paste0(
            "POSIX path of (choose file with prompt \"Select unstained FCS file\" default location POSIX file ",
            gate_applescript_quote(initial_dir), ")"
        )
        result <- suppressWarnings(system2("osascript", c("-e", shQuote(script)), stdout = TRUE, stderr = FALSE))
        if (!is.null(attr(result, "status")) || length(result) == 0L) return(NULL)
        return(normalizePath(trimws(result[[1]]), mustWork = TRUE))
    }
    if (.Platform$OS.type == "windows") {
        selected <- tryCatch(file.choose(new = FALSE), error = function(e) "")
        if (!nzchar(selected)) return(NULL)
        return(normalizePath(selected, mustWork = TRUE))
    }
    picker <- Sys.which("zenity")
    if (!nzchar(picker)) picker <- Sys.which("kdialog")
    if (!nzchar(picker)) stop("No graphical file picker is available. Install zenity or kdialog.", call. = FALSE)
    args <- if (grepl("zenity$", picker)) c("--file-selection", "--title=Select unstained FCS file", paste0("--filename=", initial_dir, "/"), "--file-filter=FCS files | *.fcs *.FCS") else c("--getopenfilename", initial_dir, "FCS files (*.fcs *.FCS)")
    result <- suppressWarnings(system2(picker, args, stdout = TRUE, stderr = FALSE))
    if (!is.null(attr(result, "status")) || length(result) == 0L) return(NULL)
    normalizePath(trimws(result[[1]]), mustWork = TRUE)
}

#* Open the native FCS picker for standalone AF extraction
#* @post /af_profiles/select-source
function(req) {
    body <- gui_workflow_body(req)
    tryCatch({
        selected <- gui_pick_af_source_file(gui_project_input_path(gui_project_value(body), "controls"))
        if (is.null(selected) || !nzchar(selected)) return(list(success = FALSE, cancelled = TRUE))
        if (!grepl("\\.fcs$", selected, ignore.case = TRUE)) return(list(success = FALSE, cancelled = FALSE, error = "Select an FCS file."))
        list(success = TRUE, cancelled = FALSE, path = selected)
    }, error = function(e) list(success = FALSE, cancelled = FALSE, error = conditionMessage(e)))
}

#* Use a saved AF profile as the active dataset unstained control
#* @post /af_profiles/activate
function(req) {
    body <- jsonlite::fromJSON(req$postBody, simplifyVector = TRUE)
    root <- gui_project_value(body)
    name <- if (!is.null(body$profile_name)) trimws(as.character(body$profile_name)[1]) else ""
    if (!nzchar(name)) return(list(success = FALSE, error = "A profile name is required."))
    tryCatch({
        active <- gui_write_active_af_profile(name, root)
        gui_reset_project_gate_session(root)
        list(success = TRUE, profile_name = active)
    }, error = function(e) list(success = FALSE, error = conditionMessage(e)))
}

#* Stop using a saved AF profile for the active dataset
#* @post /af_profiles/deactivate
function(req) {
    body <- jsonlite::fromJSON(req$postBody, simplifyVector = TRUE)
    root <- gui_project_value(body)
    name <- if (!is.null(body$profile_name)) trimws(as.character(body$profile_name)[1]) else ""
    if (!nzchar(name)) return(list(success = FALSE, error = "A profile name is required."))
    tryCatch({
        removed <- gui_unlink_active_af_profile(name, root)
        gui_reset_project_gate_session(root)
        list(success = TRUE, profile_name = removed)
    }, error = function(e) list(success = FALSE, error = conditionMessage(e)))
}

#* Delete a saved AF profile
#* @delete /af_profiles/delete
#* @param name Profile name
function(name = "", project_path = "") {
    if (is.null(name) || !nzchar(trimws(as.character(name)[1]))) {
        return(list(success = FALSE, error = "A profile name is required."))
    }
    tryCatch({
        spectreasy::delete_af_profile(as.character(name)[1])
        root <- gui_request_project_root(project_path)
        if (identical(gui_read_active_af_profile(root), as.character(name)[1])) {
            unlink(gui_active_af_config_path(root), force = TRUE)
        }
        list(success = TRUE, name = as.character(name)[1])
    }, error = function(e) list(success = FALSE, error = conditionMessage(e)))
}

#* Rename a saved AF profile and keep an active project link synchronized
#* @post /af_profiles/rename
function(req) {
    body <- gui_workflow_body(req)
    root <- gui_project_value(body)
    name <- trimws(as.character(gui_workflow_value(body, "profile_name", ""))[1])
    new_name <- trimws(as.character(gui_workflow_value(body, "new_name", ""))[1])
    if (!nzchar(name) || !nzchar(new_name)) return(list(success = FALSE, error = "Current and new profile names are required."))
    tryCatch({
        was_active <- identical(gui_read_active_af_profile(root), name)
        spectreasy::rename_af_profile(name, new_name)
        if (was_active) {
            linked <- tryCatch(gui_write_active_af_profile(new_name, root), error = function(e) e)
            if (inherits(linked, "error")) {
                try(spectreasy::rename_af_profile(new_name, name), silent = TRUE)
                stop(conditionMessage(linked), call. = FALSE)
            }
            gui_reset_project_gate_session(root)
        }
        list(success = TRUE, profile_name = name, new_name = new_name, active = was_active)
    }, error = function(e) list(success = FALSE, error = conditionMessage(e)))
}

#* Apply a saved AF profile to a matrix
#* @post /af_profiles/apply
function(req) {
    body <- jsonlite::fromJSON(req$postBody, simplifyVector = TRUE)
    body_value <- function(name, fallback = "") {
        value <- body[[name]]
        if (is.null(value) || length(value) == 0 || is.na(value[1])) fallback else value[1]
    }
    matrix_filename <- trimws(as.character(body_value("matrix_filename"))[1])
    profile_name <- trimws(as.character(body_value("profile_name"))[1])
    output_filename <- trimws(as.character(body_value("output_filename", matrix_filename))[1])
    if (!nzchar(matrix_filename) || !nzchar(profile_name) || !nzchar(output_filename)) {
        return(list(success = FALSE, error = "Matrix, profile, and output names are required."))
    }
    tryCatch({
        root <- gui_project_value(body)
        matrix_file <- matrix_path(matrix_filename, root)
        output_file <- matrix_path(output_filename, root)
        if (!file.exists(matrix_file)) stop("Matrix file not found: ", matrix_filename, call. = FALSE)
        matrix_data <- read_matrix_csv(matrix_file)
        profile <- spectreasy::load_af_profile(profile_name, show_plot = FALSE)
        adjusted <- spectreasy::add_af_profile(matrix_data, profile, replace_existing = TRUE)
        dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
        utils::write.csv(adjusted, output_file, row.names = FALSE, quote = TRUE)
        list(success = TRUE, path = output_file, filename = output_filename)
    }, error = function(e) list(success = FALSE, error = conditionMessage(e)))
}

#* List available matrices
#* @get /matrices
function(project_path = "") {
    root <- gui_request_project_root(project_path)
    files <- list_matrix_csv_files(root)
    if (length(files) == 0) {
        return(character(0))
    }
    paths <- file.path(root, files)
    keep <- vapply(paths, is_probably_matrix_csv, logical(1))
    return(as.character(files[keep]))
}

#* List available sample files
#* @get /samples
function(project_path = "") {
    samples_dir <- gui_project_input_path(gui_request_project_root(project_path), "samples")
    if (!dir.exists(samples_dir)) return(character(0))
    files <- list.files(samples_dir, pattern = "\\.fcs$", ignore.case = TRUE)
    return(as.character(sort(files)))
}

#* List GUI config presets
#* @get /configs
function() {
    cfg_dir <- get_config_dir()
    files <- list.files(cfg_dir, pattern = "\\.json$", ignore.case = TRUE)
    return(as.character(sort(files)))
}

#* Load GUI config preset
#* @get /load_config
#* @param filename
function(filename) {
    cfg_name <- normalize_config_filename(filename)
    path <- file.path(get_config_dir(), cfg_name)
    if (!file.exists(path)) {
        return(list(error = paste("Config file not found:", path)))
    }
    cfg <- jsonlite::fromJSON(path, simplifyVector = TRUE)
    return(cfg)
}

#* Save GUI config preset
#* @post /save_config
function(req) {
    body <- jsonlite::fromJSON(req$postBody)
    cfg_name <- normalize_config_filename(body$filename)
    cfg <- body$config_json
    if (is.null(cfg)) {
        return(list(error = "Missing config_json in request body"))
    }
    path <- file.path(get_config_dir(), cfg_name)
    jsonlite::write_json(cfg, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
    return(list(success = TRUE, filename = cfg_name, path = path))
}

#* Load a specific matrix
#* @get /load_matrix
#* @param filename
function(filename, project_path = "") {
    path <- matrix_path(filename, gui_request_project_root(project_path))
    if (!file.exists(path)) {
        return(list(error = paste("File not found:", path)))
    }

    df <- read_matrix_csv(path)

    # Filter out AF rows (matching ^AF($|_) case-insensitively)
    df <- df[!is_af_matrix_row(df), , drop = FALSE]

    return(df)
}

#* Save the adjusted matrix
#* @post /save_matrix
function(req) {
    body <- jsonlite::fromJSON(req$postBody)
    root <- gui_project_value(body)
    filename <- body$filename
    source_filename <- body$source_filename
    matrix_data <- body$matrix_json
    df <- as.data.frame(matrix_data, check.names = FALSE)
    path <- matrix_path(filename, root)
    source_path <- if (!is.null(source_filename) && length(source_filename) &&
        !is.na(source_filename[1]) && nzchar(trimws(as.character(source_filename[1])))) {
        matrix_path(as.character(source_filename[1]), root)
    } else {
        character()
    }
    if (nrow(df) < 1L || ncol(df) < 2L) stop("Matrix must contain at least one row and one detector column.", call. = FALSE)
    if (!"Marker" %in% colnames(df)) stop("Matrix must contain a Marker column.", call. = FALSE)
    markers <- trimws(as.character(df$Marker))
    if (any(is.na(markers) | !nzchar(markers)) || anyDuplicated(markers)) stop("Matrix row names must be non-empty and unique.", call. = FALSE)
    detector_names <- setdiff(colnames(df), "Marker")
    if (any(!nzchar(trimws(detector_names))) || anyDuplicated(detector_names)) stop("Matrix detector names must be non-empty and unique.", call. = FALSE)
    numeric_columns <- lapply(df[detector_names], function(values) suppressWarnings(as.numeric(values)))
    if (any(vapply(seq_along(numeric_columns), function(index) any(is.na(numeric_columns[[index]]) & !is.na(df[[detector_names[index]]])), logical(1)))) {
        stop("Matrix detector values must be numeric.", call. = FALSE)
    }
    df[detector_names] <- numeric_columns
    df$Marker <- markers
    df <- merge_hidden_af_rows(df, source_paths = c(path, source_path))

    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    temporary <- tempfile("matrix_save_", tmpdir = dirname(path), fileext = ".csv")
    on.exit(unlink(temporary, force = TRUE), add = TRUE)
    utils::write.csv(df, temporary, row.names = FALSE, quote = TRUE)
    verified <- read_matrix_csv(temporary)
    if (!identical(dim(verified), dim(df)) || !identical(colnames(verified), colnames(df))) {
        stop("Temporary matrix verification failed; the original file was not replaced.", call. = FALSE)
    }
    if (!file.rename(temporary, path)) {
        stop("Could not atomically replace the selected matrix file.", call. = FALSE)
    }
    return(list(success = TRUE, path = path))
}

#* CORS preflight for import_matrix
#* @options /import_matrix
function(res) {
    return("")
}

#* Import a matrix from an external path
#* @post /import_matrix
#* @param path Absolute path to the CSV file
function(path, project_path = "") {
    if (!file.exists(path)) {
        return(list(error = paste("File not found:", path)))
    }
    filename <- basename(path)
    dest <- file.path(gui_request_project_root(project_path), filename)
    file.copy(path, dest, overwrite = TRUE)
    return(list(success = TRUE, filename = filename))
}

#* CORS preflight for import_matrix_content
#* @options /import_matrix_content
function(res) {
    return("")
}

#* Import a matrix from uploaded content
#* @post /import_matrix_content
#* @param filename The filename
#* @param content The CSV content as text
function(filename, content, project_path = "") {
    dest <- matrix_path(filename, gui_request_project_root(project_path))
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    writeLines(content, dest)
    return(list(success = TRUE, filename = filename))
}

#* Get unmixed data for a subset of cells from a sample
#* @get /data
#* @param sample_name The name of the sample to load
function(sample_name = "", project_path = "") {
    # If no sample is provided, use the first file in the configured sample input directory.
    samples_dir <- gui_project_input_path(gui_request_project_root(project_path), "samples")
    files <- sort(list.files(samples_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE))
    if (length(files) == 0) {
        return(list(error = paste0("No FCS files found in samples directory: ", samples_dir)))
    }

    if (is.null(sample_name) || !nzchar(trimws(as.character(sample_name)))) {
        sample_path <- files[1]
    } else {
        sn <- trimws(as.character(sample_name))
        candidates <- c(
            file.path(samples_dir, sn),
            file.path(samples_dir, paste0(sn, ".fcs"))
        )
        sample_path <- candidates[file.exists(candidates)][1]
        if (is.na(sample_path) || !nzchar(sample_path)) {
            target <- tolower(tools::file_path_sans_ext(basename(sn)))
            file_ids <- tolower(tools::file_path_sans_ext(basename(files)))
            idx <- which(file_ids == target)
            if (length(idx) > 0) {
                sample_path <- files[idx[1]]
            } else {
                return(list(error = paste("Sample not found:", sn)))
            }
        }
    }

    ff <- gate_read_fcs(sample_path)
    raw_data <- flowCore::exprs(ff)

    # Subsample for speed - smaller for fast interactive updates
    n_sub <- 2000
    if (nrow(raw_data) > n_sub) {
        set.seed(123)
        raw_data <- raw_data[sample(nrow(raw_data), n_sub), ]
    }

    pd <- flowCore::pData(flowCore::parameters(ff))
    # Helper to get sorted detectors (copying logic from spectreasy if not exported)
    # Assuming spectreasy is loaded or we implement basic logic
    det_info <- tryCatch(
        {
            spectreasy::get_sorted_detectors(pd)
        },
        error = function(e) {
            # Fallback if function not accessible
            fl_cols <- grep("FL", pd$name, value = TRUE)
            list(names = fl_cols, labels = pd$desc[match(fl_cols, pd$name)])
        }
    )

    return(list(
        raw_data = as.data.frame(raw_data),
        sample_name = basename(sample_path),
        detector_names = det_info$names,
        detector_labels = det_info$labels
    ))
}

#* CORS preflight for unmix
#* @options /unmix
function(res) {
    return("")
}

#* Run unmixing (On-demand unmixing endpoint)
#* @post /unmix
#* @param matrix_json The matrix (M or W)
#* @param raw_data_json The raw data
#* @param type "reference" (M) or "unmixing" (W)
function(matrix_json, raw_data_json, type = "reference", matrix_filename = "", method = "", project_path = "") {
    # matrix_json format: {MarkerName: {det1: val, det2: val, ...}, ...}
    markers <- names(matrix_json)
    detectors <- names(matrix_json[[1]])

    mat <- matrix(0, nrow = length(markers), ncol = length(detectors))
    for (i in seq_along(markers)) {
        mat[i, ] <- as.numeric(unlist(matrix_json[[markers[i]]][detectors]))
    }
    rownames(mat) <- markers
    colnames(mat) <- detectors

    Y <- as.matrix(raw_data_to_df(raw_data_json))
    suppressWarnings(storage.mode(Y) <- "numeric")

    if (!grepl("unmixing", tolower(type)) &&
        !is.null(matrix_filename) &&
        nzchar(trimws(as.character(matrix_filename)[1]))) {
        mat_df <- as.data.frame(mat, check.names = FALSE)
        mat_df$Marker <- rownames(mat)
        mat_df <- mat_df[, c("Marker", colnames(mat)), drop = FALSE]
        mat_df <- merge_hidden_af_rows(mat_df, source_paths = matrix_path(matrix_filename, gui_request_project_root(project_path)))
        markers <- as.character(mat_df[[1]])
        mat <- as.matrix(mat_df[, -1, drop = FALSE])
        storage.mode(mat) <- "numeric"
        rownames(mat) <- markers
    }

    # Matching columns
    common_dets <- intersect(colnames(Y), colnames(mat))
    if (length(common_dets) == 0) {
        marker_cols <- intersect(colnames(Y), markers)
        if (length(marker_cols) > 0) {
            return(as.data.frame(Y[, marker_cols, drop = FALSE]))
        }
        return(list(error = "No matching detectors found between data and matrix"))
    }

    Y_sub <- Y[, common_dets, drop = FALSE]
    mat_sub <- mat[, common_dets, drop = FALSE]

    # Perform Unmixing
    if (grepl("unmixing", tolower(type))) {
        # Provided matrix IS the unmixing matrix (W)
        # Unmixed = Raw * t(W)
        unmixed <- Y_sub %*% t(mat_sub)
    } else {
        method_resolved <- tryCatch(
            spectreasy:::.normalize_unmix_method(if (!nzchar(trimws(as.character(method)[1]))) get_unmixing_method() else method),
            error = function(e) e
        )
        if (inherits(method_resolved, "error")) {
            return(list(error = conditionMessage(method_resolved)))
        }
        ff <- flowCore::flowFrame(Y_sub)
        res <- tryCatch(
            spectreasy::calc_residuals(ff, mat_sub, method = method_resolved),
            error = function(e) {
                return(list(error = conditionMessage(e)))
            }
        )
        if (is.list(res) && !is.data.frame(res) && !is.null(res$error)) {
            return(res)
        }
        unmixed <- as.matrix(res[, rownames(mat_sub), drop = FALSE])
    }

    return(as.data.frame(unmixed))
}


#* Project scan and workflow prerequisites
#* @get /project/status
#* @param project_path Optional project directory
function(project_path = "") {
    if (!isTRUE(getOption("spectreasy.project_selected", TRUE)) &&
        (is.null(project_path) || !nzchar(trimws(as.character(project_path)[1])))) {
        return(list(
            project_path = "",
            files = character(),
            scan = list(controls = 0, samples = 0, matrices = 0, reports = 0, gates = 0, qc_metrics = 0, spectral_variants = 0),
            summary = "no project selected",
            recommended_next_action = "Choose a project folder"
        ))
    }
    root <- if (is.null(project_path) || !nzchar(trimws(as.character(project_path)[1]))) get_matrix_dir() else as.character(project_path)[1]
    root <- normalizePath(root, mustWork = FALSE)
    gui_project_scan(root)
}

#* Read the active project's control and sample input-directory configuration
#* @get /project/layout
function(project_path = "") {
    tryCatch({
        root <- gui_request_project_root(project_path)
        list(success = TRUE, layout = gui_project_layout(root, persist = FALSE), project = gui_project_scan(root))
    }, error = function(e) list(success = FALSE, error = conditionMessage(e)))
}

#* Rename or adopt one project input directory and persist the new layout
#* @post /project/layout
function(req) {
    body <- gui_workflow_body(req)
    tryCatch({
        root <- gui_project_value(body)
        layout <- gui_update_project_input_dir(
            root,
            gui_workflow_value(body, "role", ""),
            gui_workflow_value(body, "path", "")
        )
        list(success = TRUE, layout = layout, project = gui_project_scan(root))
    }, error = function(e) list(success = FALSE, error = conditionMessage(e)))
}

gui_project_relative_path <- function(path, root) {
    root_pattern <- gsub("([.|(){}+*?^$\\[\\]\\\\])", "\\\\\\1", normalizePath(root, mustWork = FALSE))
    gsub("\\\\", "/", sub(paste0("^", root_pattern, "[/\\\\]?"), "", path))
}

gui_project_report_files <- function(root, files = NULL, output_root = "", report_type = "") {
    root <- normalizePath(root, mustWork = FALSE)
    output_root <- trimws(as.character(output_root)[1])
    report_type <- tolower(trimws(as.character(report_type)[1]))
    if (nzchar(output_root) && report_type %in% c("control", "sample")) {
        output_base <- if (grepl("^(/|[A-Za-z]:[/\\\\])", output_root)) output_root else file.path(root, output_root)
        output_base <- normalizePath(output_base, mustWork = FALSE)
        root_prefix <- paste0(root, .Platform$file.sep)
        if (!identical(output_base, root) && !startsWith(output_base, root_prefix)) {
            stop("Report output folder is outside the active project.", call. = FALSE)
        }
        stage_dir <- file.path(output_base, if (identical(report_type, "sample")) "unmix_samples" else "unmix_controls")
        report_dir_name <- if (identical(report_type, "sample")) "qc_samples" else "qc_controls"
        report_dirs <- if (dir.exists(stage_dir)) {
            candidates <- list.dirs(stage_dir, recursive = FALSE, full.names = TRUE)
            candidates[grepl(paste0("^", report_dir_name, "(?:_(?:[2-9]|[1-9][0-9]+))?$"), basename(candidates), perl = TRUE)]
        } else character()
        report_name <- if (identical(report_type, "sample")) "qc_samples_report" else "qc_controls_report"
        legacy_reports <- unlist(lapply(report_dirs, function(report_dir) {
            list.files(
                report_dir,
                pattern = paste0("^", report_name, "\\.(html?|pdf)$"),
                full.names = TRUE,
                ignore.case = TRUE
            )
        }), use.names = FALSE)
        return(unique(legacy_reports))
    }
    if (is.null(files)) {
        top_level <- if (dir.exists(root)) list.dirs(root, recursive = FALSE, full.names = TRUE) else character()
        report_roots <- top_level[
            grepl("^(reports|spectreasy_outputs(?:_[^/]+)?)$", basename(top_level), ignore.case = TRUE, perl = TRUE)
        ]
        files <- unique(unlist(lapply(report_roots, function(directory) {
            list.files(directory, recursive = TRUE, full.names = TRUE, all.files = FALSE)
        }), use.names = FALSE))
    }
    relative_files <- gui_project_relative_path(files, root)
    files[
        grepl("\\.(html?|pdf)$", relative_files, ignore.case = TRUE) &
            grepl("(^|/)(reports|spectreasy_outputs(?:_[^/]+)?)/", relative_files, ignore.case = TRUE, perl = TRUE)
    ]
}

gui_project_report_source_files <- function(root, sample = FALSE, output_root = "") {
    input_dir <- gui_project_input_path(root, if (isTRUE(sample)) "samples" else "controls")
    input_files <- if (dir.exists(input_dir)) {
        list.files(input_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
    } else {
        character()
    }
    shared <- file.path(root, c("fcs_mapping.csv", "ssc_gate_config.csv"))
    shared <- shared[file.exists(shared)]
    output_root <- trimws(as.character(output_root)[1])
    output_roots <- if (nzchar(output_root)) {
        candidate <- if (grepl("^(/|[A-Za-z]:[/\\\\])", output_root)) output_root else file.path(root, output_root)
        candidate <- normalizePath(candidate, mustWork = FALSE)
        if (dir.exists(candidate)) candidate else character()
    } else {
        candidates <- list.dirs(root, recursive = FALSE, full.names = TRUE)
        candidates[grepl("^spectreasy_outputs(?:_[^/]+)?$", basename(candidates), ignore.case = TRUE, perl = TRUE)]
    }
    artifact_pattern <- if (isTRUE(sample)) {
        "(reference_matrix|detector_noise|spectral_variant).*\\.(csv|rds)$"
    } else {
        "(reference_matrix|detector_noise).*\\.csv$"
    }
    artifacts <- unique(unlist(lapply(output_roots, function(directory) {
        list.files(directory, pattern = artifact_pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    }), use.names = FALSE))
    unique(c(input_files, shared, artifacts))
}

gui_project_report_type <- function(relative_path) {
    relative_path <- gsub("\\\\", "/", as.character(relative_path)[1])
    if (grepl("panel", relative_path, ignore.case = TRUE)) return("Panel overview")
    if (grepl("(^|/)(unmix_samples|qc_samples)(/|$)|qc_samples_report", relative_path, ignore.case = TRUE, perl = TRUE)) {
        return("Sample QC")
    }
    "Control QC"
}

gui_project_report_prompt_file <- function(report) {
    report <- normalizePath(report, mustWork = FALSE)
    candidate <- file.path(
        dirname(report),
        paste0(tools::file_path_sans_ext(basename(report)), "_ai_qc_prompt.txt")
    )
    if (file.exists(candidate) && !dir.exists(candidate)) return(normalizePath(candidate, mustWork = TRUE))

    ""
}

gui_resolve_project_report <- function(path, root = get_matrix_dir()) {
    root <- normalizePath(root, mustWork = FALSE)
    relative <- gsub("\\\\", "/", trimws(as.character(path)[1]))
    parts <- strsplit(relative, "/", fixed = TRUE)[[1]]
    if (!nzchar(relative) || any(parts %in% c("", ".", ".."))) {
        stop("Invalid project report path.", call. = FALSE)
    }
    target <- normalizePath(file.path(root, do.call(file.path, as.list(parts))), mustWork = FALSE)
    root_prefix <- paste0(root, .Platform$file.sep)
    if (!identical(target, root) && !startsWith(target, root_prefix)) {
        stop("Project report is outside the active project.", call. = FALSE)
    }
    if (!file.exists(target) || dir.exists(target)) {
        stop("Project report not found.", call. = FALSE)
    }
    if (!tolower(tools::file_ext(target)) %in% c("html", "htm")) {
        stop("Only an existing HTML report can be exported to PDF.", call. = FALSE)
    }
    target
}

gui_find_chromium <- function() {
    candidates <- unique(c(
        unname(Sys.which(c("google-chrome", "google-chrome-stable", "chromium", "chromium-browser", "chrome"))),
        "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
        "/Applications/Chromium.app/Contents/MacOS/Chromium",
        file.path(Sys.getenv("PROGRAMFILES"), "Google", "Chrome", "Application", "chrome.exe"),
        file.path(Sys.getenv("PROGRAMFILES(X86)"), "Google", "Chrome", "Application", "chrome.exe")
    ))
    candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
    if (length(candidates)) candidates[[1]] else ""
}

gui_export_html_report_pdf <- function(html_file, output_file = tempfile(fileext = ".pdf")) {
    browser <- gui_find_chromium()
    if (!nzchar(browser)) {
        stop("Chrome or Chromium is required to export the existing HTML report as PDF.", call. = FALSE)
    }
    if (!requireNamespace("chromote", quietly = TRUE)) {
        stop("The optional 'chromote' package is required to export the existing HTML report as PDF.", call. = FALSE)
    }
    html_file <- normalizePath(html_file, mustWork = TRUE)
    output_file <- normalizePath(output_file, mustWork = FALSE)
    previous_browser <- Sys.getenv("CHROMOTE_CHROME", unset = NA_character_)
    Sys.setenv(CHROMOTE_CHROME = browser)
    on.exit({
        if (is.na(previous_browser)) Sys.unsetenv("CHROMOTE_CHROME") else Sys.setenv(CHROMOTE_CHROME = previous_browser)
    }, add = TRUE)
    session <- chromote::ChromoteSession$new()
    on.exit(try(session$close(), silent = TRUE), add = TRUE)
    session$Page$navigate(utils::URLencode(paste0("file://", html_file), reserved = FALSE))
    for (attempt in seq_len(40L)) {
        state <- tryCatch(
            session$Runtime$evaluate("document.readyState", returnByValue = TRUE)$result$value,
            error = function(e) ""
        )
        if (identical(state, "complete")) break
        Sys.sleep(0.1)
    }
    Sys.sleep(0.4)
    pdf <- session$Page$printToPDF(printBackground = TRUE, preferCSSPageSize = TRUE)
    base64enc::base64decode(pdf$data, output = output_file)
    if (!file.exists(output_file) || file.info(output_file)$size < 5L) {
        stop("Chrome could not export the HTML report as PDF.", call. = FALSE)
    }
    output_file
}

#* Discover report artifacts and compare them with upstream project files
#* @get /project/reports
function(project_path = "", output_root = "", report_type = "") {
    root <- if (is.null(project_path) || !nzchar(trimws(as.character(project_path)[1]))) get_matrix_dir() else as.character(project_path)[1]
    root <- normalizePath(root, mustWork = FALSE)
    reports <- gui_project_report_files(root, output_root = output_root, report_type = report_type)
    if (length(reports) == 0L) return(list(reports = data.frame()))
    relative <- function(path) gui_project_relative_path(path, root)
    control_sources <- gui_project_report_source_files(root, sample = FALSE, output_root = output_root)
    sample_sources <- gui_project_report_source_files(root, sample = TRUE, output_root = output_root)
    rows <- lapply(reports, function(report) {
        relative_report <- relative(report)
        report_type <- gui_project_report_type(relative_report)
        is_sample <- identical(report_type, "Sample QC")
        source_pattern <- if (is_sample) {
            "\\.fcs$|reference_matrix.*\\.csv$|detector_noise.*\\.csv$|spectral_variant.*\\.rds$"
        } else {
            "\\.fcs$|fcs_mapping\\.csv$|gate.*\\.csv$|reference_matrix.*\\.csv$|detector_noise.*\\.csv$"
        }
        source_files <- if (is_sample) sample_sources else control_sources
        source_files <- source_files[grepl(source_pattern, relative(source_files), ignore.case = TRUE, perl = TRUE)]
        source_files <- setdiff(source_files, report)
        stale <- length(source_files) > 0L && any(file.info(source_files)$mtime > file.info(report)$mtime)
        prompt_file <- gui_project_report_prompt_file(report)
        data.frame(
            path = relative_report,
            report_type = report_type,
            format = if (grepl("\\.pdf$", report, ignore.case = TRUE)) "PDF" else "HTML",
            created = format(file.info(report)$mtime, "%Y-%m-%dT%H:%M:%S%z"),
            created_epoch = as.numeric(file.info(report)$mtime),
            status = if (stale) "stale" else "current",
            prompt_path = if (nzchar(prompt_file)) relative(prompt_file) else "",
            prompt_bytes = if (nzchar(prompt_file)) as.numeric(file.info(prompt_file)$size) else 0,
            stringsAsFactors = FALSE
        )
    })
    list(reports = do.call(rbind, rows))
}

#* Stream a project artifact through the local R backend
#* @get /project/file
#* @param path Relative path inside the active project
#* @param project_path Active cockpit project directory
function(path = "", project_path = "", req, res) {
    if (!gui_api_token_allowed(req)) {
        res$status <- 403
        return(list(error = "This project artifact belongs to a different or inactive Spectreasy session."))
    }
    root <- if (is.null(project_path) || !nzchar(trimws(as.character(project_path)[1]))) get_matrix_dir() else as.character(project_path)[1]
    root <- normalizePath(root, mustWork = FALSE)
    relative <- gsub("\\\\", "/", trimws(as.character(path)[1]))
    parts <- strsplit(relative, "/", fixed = TRUE)[[1]]
    if (!nzchar(relative) || any(parts %in% c("", ".", ".."))) {
        res$status <- 400
        return(list(error = "Invalid project artifact path."))
    }
    target <- normalizePath(file.path(root, do.call(file.path, as.list(parts))), mustWork = FALSE)
    root_prefix <- paste0(root, .Platform$file.sep)
    if (!identical(target, root) && !startsWith(target, root_prefix)) {
        res$status <- 403
        return(list(error = "Artifact path is outside the active project."))
    }
    if (!file.exists(target) || dir.exists(target)) {
        res$status <- 404
        return(list(error = "Project artifact not found."))
    }
    extension <- tolower(tools::file_ext(target))
    content_type <- switch(
        extension,
        html = "text/html; charset=utf-8",
        htm = "text/html; charset=utf-8",
        pdf = "application/pdf",
        csv = "text/csv; charset=utf-8",
        json = "application/json; charset=utf-8",
        md = "text/markdown; charset=utf-8",
        txt = "text/plain; charset=utf-8",
        "application/octet-stream"
    )
    res$setHeader("Content-Type", content_type)
    res$body <- readBin(target, what = "raw", n = file.info(target)$size)
    res
}

#* Export an existing HTML QC report to PDF without rerunning QC
#* @post /project/report/export-pdf
function(req) {
    body <- gui_workflow_body(req)
    tryCatch({
        root <- gui_workflow_root(body)
        report_file <- gui_resolve_project_report(gui_workflow_value(body, "path", ""), root = root)
        pdf_file <- gui_export_html_report_pdf(report_file)
        on.exit(unlink(pdf_file, force = TRUE), add = TRUE)
        list(
            success = TRUE,
            filename = paste0(tools::file_path_sans_ext(basename(report_file)), ".pdf"),
            content_type = "application/pdf",
            content_base64 = base64enc::base64encode(pdf_file)
        )
    }, error = function(e) list(success = FALSE, error = conditionMessage(e)))
}

#* Change the active project folder from the cockpit
#* @post /project/context
function(req) {
    body <- gui_workflow_body(req)
    project_path <- gui_workflow_path(body, "projectPath", "", allow_empty = FALSE)
    tryCatch({
        project_path <- gui_request_project_root(project_path, fallback = FALSE)
        list(success = TRUE, project = gui_project_scan(project_path))
    }, error = function(e) list(success = FALSE, error = conditionMessage(e)))
}

#* Create the standard input folders after explicit cockpit confirmation
#* @post /project/initialize
function(req) {
    body <- gui_workflow_body(req)
    tryCatch({
        root <- gui_workflow_root(body)
        created <- gui_project_relative_path(gui_missing_project_input_dirs(root), root)
        gui_ensure_project_input_dirs(root)
        list(success = TRUE, created = created, project = gui_project_scan(root))
    }, error = function(e) list(success = FALSE, error = conditionMessage(e)))
}

#* Open the native folder picker and activate the selected project
#* @post /project/select
function(req) {
    tryCatch({
        selected <- gui_pick_project_directory(get_matrix_dir(), allow_create = FALSE)
        if (is.null(selected) || !nzchar(selected)) return(list(success = FALSE, cancelled = TRUE))
        selected <- gui_request_project_root(selected, fallback = FALSE)
        list(success = TRUE, cancelled = FALSE, project = gui_project_scan(selected))
    }, error = function(e) list(success = FALSE, cancelled = FALSE, error = conditionMessage(e)))
}

#* Create and activate a project folder with the native folder picker
#* @post /project/create
function(req) {
    tryCatch({
        selected <- gui_pick_project_directory(get_matrix_dir(), allow_create = TRUE)
        if (is.null(selected) || !nzchar(selected)) return(list(success = FALSE, cancelled = TRUE))
        selected <- gui_request_project_root(selected, fallback = FALSE)
        list(success = TRUE, cancelled = FALSE, project = gui_project_scan(selected))
    }, error = function(e) list(success = FALSE, cancelled = FALSE, error = conditionMessage(e)))
}

#* List FCS files in the active project's controls or samples folder
#* @get /project/files
#* @param kind Either controls or samples
#* @param project_path Active cockpit project directory
function(kind = "controls", project_path = "") {
    tryCatch(
        list(success = TRUE, files = gui_project_file_rows(kind, project_path)),
        error = function(e) list(success = FALSE, error = conditionMessage(e), files = data.frame())
    )
}

#* Begin a bounded-memory FCS upload
#* @post /project/upload-start
function(req) {
    body <- gui_workflow_body(req)
    tryCatch({
        gui_project_upload_discard_stale()
        location <- gui_project_file_location(
            gui_workflow_value(body, "kind", ""),
            gui_workflow_value(body, "filename", ""),
            gui_project_value(body)
        )
        if (!dir.exists(location$directory)) {
            stop("Create the project's configured control and sample input folders before adding files.", call. = FALSE)
        }
        size <- gui_workflow_number(body, "size", 0, integer = TRUE, minimum = 1)
        maximum <- as.numeric(getOption("spectreasy.gui_max_upload_bytes", 50 * 1024^3))
        if (is.finite(maximum) && size > maximum) {
            stop("Uploaded file exceeds the configured size limit.", call. = FALSE)
        }
        if (file.exists(location$path)) stop("A file named '", location$filename, "' already exists.", call. = FALSE)
        upload_id <- gui_project_upload_id()
        temporary <- tempfile("spectreasy-upload-", tmpdir = location$directory)
        if (!file.create(temporary)) stop("Could not prepare the upload destination.", call. = FALSE)
        assign(upload_id, list(
            location = location,
            temporary = temporary,
            size = size,
            written = 0,
            created = Sys.time()
        ), envir = .gui_project_uploads)
        list(success = TRUE, upload_id = upload_id)
    }, error = function(e) list(success = FALSE, error = conditionMessage(e)))
}

#* Append one encoded chunk to an FCS upload
#* @post /project/upload-chunk
function(req) {
    body <- gui_workflow_body(req)
    upload_id <- gui_workflow_value(body, "upload_id", "")
    tryCatch({
        session <- gui_project_upload_session(upload_id)
        offset <- gui_workflow_number(body, "offset", 0, integer = TRUE, minimum = 0)
        if (!identical(as.numeric(offset), as.numeric(session$written))) {
            stop("Upload chunk is out of sequence; restart the upload.", call. = FALSE)
        }
        content <- gui_workflow_value(body, "content_base64", "")
        payload <- tryCatch(jsonlite::base64_dec(content), error = function(e) NULL)
        if (is.null(payload) || !length(payload)) stop("Upload chunk is empty or invalid.", call. = FALSE)
        if (session$written + length(payload) > session$size) stop("Upload exceeds the declared file size.", call. = FALSE)
        connection <- file(session$temporary, open = "ab")
        on.exit(try(close(connection), silent = TRUE), add = TRUE)
        writeBin(payload, connection)
        close(connection)
        session$written <- session$written + length(payload)
        assign(upload_id, session, envir = .gui_project_uploads)
        list(success = TRUE, written = session$written)
    }, error = function(e) {
        gui_project_upload_discard(upload_id)
        list(success = FALSE, error = conditionMessage(e))
    })
}

#* Validate and commit an FCS upload
#* @post /project/upload-finish
function(req) {
    body <- gui_workflow_body(req)
    upload_id <- gui_workflow_value(body, "upload_id", "")
    tryCatch({
        session <- gui_project_upload_session(upload_id)
        if (!identical(as.numeric(session$written), as.numeric(session$size))) {
            stop("Upload is incomplete.", call. = FALSE)
        }
        gui_validate_fcs_upload(session$temporary)
        if (file.exists(session$location$path)) stop("The destination file now exists; the upload was not overwritten.", call. = FALSE)
        if (!file.rename(session$temporary, session$location$path)) stop("Could not save uploaded file.", call. = FALSE)
        rm(list = upload_id, envir = .gui_project_uploads)
        rows <- gui_project_file_rows(session$location$kind, dirname(session$location$directory))
        row <- rows[rows$name == session$location$filename, , drop = FALSE]
        list(success = TRUE, file = row)
    }, error = function(e) {
        gui_project_upload_discard(upload_id)
        list(success = FALSE, error = conditionMessage(e))
    })
}

#* Cancel an incomplete FCS upload
#* @post /project/upload-abort
function(req) {
    body <- gui_workflow_body(req)
    upload_id <- gui_workflow_value(body, "upload_id", "")
    gui_project_upload_discard(upload_id)
    list(success = TRUE)
}

#* Delete one FCS file from the active project's controls or samples folder
#* @delete /project/files
#* @param kind Either controls or samples
#* @param filename FCS filename inside that folder
#* @param project_path Active cockpit project directory
function(kind = "", filename = "", project_path = "") {
    tryCatch({
        location <- gui_project_file_location(kind, filename, project_path)
        if (!file.exists(location$path) || dir.exists(location$path)) stop("Project file not found.", call. = FALSE)
        removed <- unlink(location$path, force = TRUE)
        if (!identical(removed, 0L)) stop("Could not delete project file.", call. = FALSE)
        list(success = TRUE, filename = location$filename)
    }, error = function(e) list(success = FALSE, error = conditionMessage(e)))
}

#* Delete all FCS files from one active project input folder
#* @delete /project/files/all
#* @param kind Either controls or samples
#* @param project_path Active cockpit project directory
function(kind = "", project_path = "") {
    tryCatch({
        location <- gui_project_file_location(kind, project_path = project_path)
        if (!dir.exists(location$directory)) {
            return(list(success = TRUE, deleted = 0L))
        }
        files <- list.files(location$directory, pattern = "\\.fcs$", ignore.case = TRUE, full.names = TRUE)
        if (length(files)) {
            removed <- unlink(files, force = TRUE)
            if (!identical(removed, 0L)) stop("Could not delete all project files.", call. = FALSE)
        }
        list(success = TRUE, deleted = length(files))
    }, error = function(e) list(success = FALSE, error = conditionMessage(e)))
}

#* CORS preflight for workflow_control
#* @options /workflow/control
function(res) {
    return("")
}

#* Run the complete control-stage workflow through Spectreasy R
#* @post /workflow/control
function(req) {
    body <- gui_workflow_body(req)
    root <- gui_workflow_root(body)
    layout <- gui_project_layout(root)
    scc_dir <- gui_workflow_resolve_path(gui_workflow_value(body, "scc_dir", layout$control_input_dir), root)
    control_file <- gui_workflow_resolve_path(gui_workflow_value(body, "control_file", "fcs_mapping.csv"), root)
    active_af_profile <- gui_read_active_af_profile(root)
    if (nzchar(active_af_profile)) {
        control_file <- gui_filtered_control_file(root)
        if (!identical(normalizePath(control_file, mustWork = FALSE), normalizePath(file.path(root, "fcs_mapping.csv"), mustWork = FALSE))) {
            on.exit(unlink(control_file, force = TRUE), add = TRUE)
        }
    }
    output_dir <- gui_workflow_resolve_path(gui_workflow_value(body, "output_dir", "spectreasy_outputs"), root)
    method <- gui_workflow_value(body, "method", get_unmixing_method())
    cytometer <- gui_workflow_value(body, "cytometer", "auto")
    gate_file <- gui_workflow_file_or_null(
        gui_workflow_value(body, "manual_gate_file", ""),
        root = gui_workflow_root(body)
    )
    control_args <- c(
        list(
            scc_dir = scc_dir,
            control_file = control_file,
            auto_create_mapping = gui_workflow_bool(body, "auto_create_mapping", TRUE),
            cytometer = cytometer,
            auto_unknown_fluor_policy = gui_workflow_value(body, "auto_unknown_fluor_policy", "by_channel"),
            output_dir = output_dir,
            unmixing_method = method,
            unmix_scatter_panel_size_mm = gui_workflow_number(body, "unmix_scatter_panel_size_mm", 30, minimum = 1),
            seed = gui_workflow_number(body, "seed", 1, integer = TRUE, minimum = 1),
            af_n_bands = gui_workflow_number(body, "af_n_bands", 100, integer = TRUE, minimum = 1),
            af_max_cells = gui_workflow_number(body, "af_max_cells", 50000, integer = TRUE, minimum = 1),
            default_sample_type = gui_workflow_value(body, "default_sample_type", "beads"),
            histogram_pct_beads = gui_workflow_number(body, "histogram_pct_beads", 0.98, minimum = 0, maximum = 1),
            histogram_direction_beads = gui_workflow_value(body, "histogram_direction_beads", "right"),
            histogram_pct_cells = gui_workflow_number(body, "histogram_pct_cells", 0.35, minimum = 0, maximum = 1),
            histogram_direction_cells = gui_workflow_value(body, "histogram_direction_cells", "right"),
            outlier_percentile = gui_workflow_number(body, "outlier_percentile", 0.02, minimum = 0, maximum = 1),
            debris_percentile = gui_workflow_number(body, "debris_percentile", 0.08, minimum = 0, maximum = 1),
            bead_gate_scale = gui_workflow_number(body, "bead_gate_scale", 1.3, minimum = 0),
            max_clusters = gui_workflow_number(body, "max_clusters", 10, integer = TRUE, minimum = 1),
            min_cluster_proportion = gui_workflow_number(body, "min_cluster_proportion", 0.03, minimum = 0, maximum = 1),
            gate_contour_beads = gui_workflow_number(body, "gate_contour_beads", 0.95, minimum = 0, maximum = 1),
            gate_contour_cells = gui_workflow_number(body, "gate_contour_cells", 0.9, minimum = 0, maximum = 1),
            subsample_n = gui_workflow_number(body, "subsample_n", 5000, integer = TRUE, minimum = 1),
            rwls_max_iter = gui_workflow_number(body, "rwls_max_iter", 1, integer = TRUE, minimum = 1),
            n_threads = gui_workflow_number(body, "n_threads", 1, integer = TRUE, minimum = 1),
            save_qc_png = gui_workflow_bool(body, "save_qc_png", TRUE),
            save_report = gui_workflow_bool(body, "save_report", TRUE),
            save_ai_qc = gui_workflow_bool(body, "save_ai_qc", gui_workflow_bool(body, "save_report", TRUE)),
            ai_qc_detail = gui_workflow_value(body, "ai_qc_detail", "standard"),
            ai_qc_privacy = gui_workflow_value(body, "ai_qc_privacy", "standard"),
            ai_qc_reference = gui_workflow_value(body, "ai_qc_reference", "auto"),
            report_format = gui_workflow_value(body, "report_format", "html"),
            gating_mode = gui_workflow_value(
                body,
                "gating_mode",
                "interactive"
            ),
            manual_gate_file = gate_file,
            scc_background_method = gui_workflow_value(body, "scc_background_method", "scatter_knn"),
            scc_background_k = gui_workflow_number(body, "scc_background_k", 2, integer = TRUE, minimum = 1),
            spectral_variant_som_nodes = gui_workflow_number(body, "spectral_variant_som_nodes", 16, integer = TRUE, minimum = 1),
            spectral_variant_top_k = gui_workflow_number(body, "spectral_variant_top_k", 3, integer = TRUE, minimum = 1),
            spectral_variant_cosine_threshold = gui_workflow_number(body, "spectral_variant_cosine_threshold", 0.98, minimum = 0, maximum = 1),
            spectral_variant_max_variants = gui_workflow_number(body, "spectral_variant_max_variants", 8, integer = TRUE, minimum = 1),
            spectral_variant_min_events = gui_workflow_number(body, "spectral_variant_min_events", 50, integer = TRUE, minimum = 1),
            autospectral_n_candidates = gui_workflow_number(body, "autospectral_n_candidates", 1000, integer = TRUE, minimum = 1),
            autospectral_n_spectral = gui_workflow_number(body, "autospectral_n_spectral", 200, integer = TRUE, minimum = 1),
            autospectral_min_events = gui_workflow_number(body, "autospectral_min_events", 10, integer = TRUE, minimum = 1),
            autospectral_refine = gui_workflow_bool(body, "autospectral_refine", FALSE),
            project_path = root
        ),
        if (nzchar(active_af_profile)) list(af_profile = active_af_profile) else list()
    )
    run <- gui_workflow_run(
        body,
        "control",
        do.call(spectreasy::unmix_controls, control_args)
    )
    if (!isTRUE(run$success)) return(run)
    result <- run$result
    run$result <- list(
        reference_matrix_file = result$reference_matrix_file,
        detector_noise_file = result$detector_noise_file,
        unmixing_matrix_file = result$unmixing_matrix_file,
        spectral_variant_library_file = result$spectral_variant_library_file,
        qc_report_file = result$qc_report_file,
        spectra_file = result$spectra_file,
        unmixing_scatter_file = result$unmixing_scatter_file,
        ai_qc_paths = result$ai_qc_paths,
        ai_qc_prompt_path = result$ai_qc_prompt_path,
        ai_qc_data_paths = result$ai_qc_data_paths
    )
    run
}

#* CORS preflight for workflow_sample
#* @options /workflow/sample
function(res) {
    return("")
}

#* Run sample unmixing and sample QC through Spectreasy R
#* @post /workflow/sample
function(req) {
    body <- gui_workflow_body(req)
    root <- gui_workflow_root(body)
    layout <- gui_project_layout(root)
    sample_dir <- gui_workflow_resolve_path(gui_workflow_value(body, "sample_dir", layout$sample_input_dir), root)
    matrix_file <- gui_workflow_resolve_path(gui_workflow_value(body, "matrix_file", file.path("spectreasy_outputs", "unmix_controls", "scc_reference_matrix.csv")), root)
    noise_file <- gui_workflow_file_or_null(gui_workflow_value(body, "detector_noise_file", ""), root = root)
    output_dir <- gui_workflow_resolve_path(gui_workflow_value(body, "output_dir", "spectreasy_outputs"), root)
    method <- gui_workflow_value(body, "method", get_unmixing_method())
    sample_args <- c(
        list(
            sample_dir = sample_dir,
            unmixing_matrix_file = matrix_file,
            detector_noise_file = noise_file,
            unmixing_method = method,
            rwls_max_iter = gui_workflow_number(body, "rwls_max_iter", 1, integer = TRUE, minimum = 1),
            n_threads = gui_workflow_number(body, "n_threads", 1, integer = TRUE, minimum = 1),
            spectral_variant_top_k = gui_workflow_number(body, "spectral_variant_top_k", 3, integer = TRUE, minimum = 1),
            spectral_variant_min_abundance = gui_workflow_number(body, "spectral_variant_min_abundance", 1, minimum = 0),
            spectral_variant_positive_fraction = gui_workflow_number(body, "spectral_variant_positive_fraction", 0.02, minimum = 0, maximum = 1),
            spectral_variant_min_improvement = gui_workflow_number(body, "spectral_variant_min_improvement", 0.01, minimum = 0),
            spectral_variant_library_file = gui_workflow_file_or_null(gui_workflow_value(body, "spectral_variant_library_file", ""), root = root),
            estimate_af = gui_workflow_bool(body, "estimate_af", FALSE),
            output_dir = output_dir,
            write_fcs = gui_workflow_bool(body, "write_fcs", TRUE),
            save_report = gui_workflow_bool(body, "save_report", TRUE),
            save_ai_qc = gui_workflow_bool(body, "save_ai_qc", gui_workflow_bool(body, "save_report", TRUE)),
            ai_qc_detail = gui_workflow_value(body, "ai_qc_detail", "standard"),
            ai_qc_privacy = gui_workflow_value(body, "ai_qc_privacy", "standard"),
            ai_qc_reference = gui_workflow_value(body, "ai_qc_reference", "auto"),
            report_format = gui_workflow_value(body, "report_format", "html"),
            report_per_sample = gui_workflow_bool(body, "report_per_sample", FALSE),
            save_qc_plots = gui_workflow_bool(body, "save_qc_plots", TRUE),
            plot_n_events = gui_workflow_number(body, "plot_n_events", 10000, integer = TRUE, minimum = 1),
            chunk_size = gui_workflow_number(body, "chunk_size", 50000, integer = TRUE, minimum = 1),
            seed = gui_workflow_number(body, "seed", 1, integer = TRUE, minimum = 1),
            return_type = gui_workflow_value(body, "return_type", "list"),
            verbose = FALSE,
            project_path = root
        )
    )
    run <- gui_workflow_run(
        body,
        "sample",
        do.call(spectreasy::unmix_samples, sample_args)
    )
    if (!isTRUE(run$success)) return(run)
    result <- run$result
    run$result <- list(
        qc_report_file = attr(result, "qc_report_file"),
        qc_samples_dir = attr(result, "qc_samples_dir"),
        qc_metrics_dir = attr(result, "qc_metrics_dir"),
        output_dir = output_dir,
        ai_qc_paths = attr(result, "ai_qc_paths"),
        ai_qc_prompt_path = attr(result, "ai_qc_prompt_path"),
        ai_qc_data_paths = attr(result, "ai_qc_data_paths")
    )
    run
}

#* Inspect AI-ready QC readiness and stale state
#* @get /ai-qc/readiness
function(project_path = "", output_root = "spectreasy_outputs") {
    root <- gui_request_project_root(project_path)
    gui_ai_qc_readiness(root, output_root = output_root)
}

#* Preview the canonical AI-ready QC summary and paste-ready prompt
#* @get /ai-qc/preview
function(project_path = "", output_root = "spectreasy_outputs", detail = "standard") {
    root <- gui_request_project_root(project_path)
    tryCatch(gui_ai_qc_preview(root, output_root = output_root, detail = detail), error = function(e) list(status = "failed", error = conditionMessage(e)))
}

#* List compatible local AI-ready QC reference profiles
#* @get /ai-qc/profiles
function(project_path = "") {
    root <- gui_request_project_root(project_path)
    gui_ai_qc_profiles(root)
}

#* Generate or refresh local AI-ready QC artifacts
#* @post /ai-qc/generate
function(req) {
    body <- gui_workflow_body(req)
    tryCatch(gui_ai_qc_generate(body), error = function(e) list(status = "failed", error = conditionMessage(e)))
}

#* CORS preflight for workflow_report
#* @options /workflow/report
function(res) {
    return("")
}

#* Compare selected unmixing methods on one active sample
#* @post /workflow/compare
function(req) {
    body <- gui_workflow_body(req)
    root <- gui_workflow_root(body)
    raw_methods <- body$methods
    methods <- unique(as.character(unlist(raw_methods, recursive = TRUE, use.names = FALSE)))
    methods <- methods[nzchar(methods)]
    if (length(methods) == 0) methods <- c("AutoSpectral", "OLS", "WLS", "NNLS")
    matrix_input <- gui_workflow_value(body, "matrix_file", "")
    sample_dir_input <- gui_workflow_value(body, "sample_dir", gui_project_layout(root)$sample_input_dir)
    matrix_file <- gui_workflow_file_or_null(matrix_input, root = root)
    sample_dir <- if (grepl("^(/|[A-Za-z]:[/\\\\])", sample_dir_input)) sample_dir_input else file.path(root, sample_dir_input)
    if (is.null(matrix_file) || !file.exists(matrix_file)) {
        return(list(success = FALSE, error = "Select a readable reference matrix before comparing methods."))
    }
    sample_files <- sort(list.files(sample_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE))
    if (length(sample_files) == 0) {
        return(list(success = FALSE, error = paste("No FCS sample files found in", sample_dir)))
    }
    matrix_df <- read_matrix_csv(matrix_file)
    if (ncol(matrix_df) < 2) return(list(success = FALSE, error = "The selected matrix has no detector columns."))
    marker_names <- as.character(matrix_df[[1]])
    matrix_values <- as.matrix(matrix_df[, -1, drop = FALSE])
    rownames(matrix_values) <- marker_names
    sample_frame <- gate_read_fcs(sample_files[1])
    events <- flowCore::exprs(sample_frame)
    if (nrow(events) > 2000) {
        set.seed(gui_workflow_number(body, "seed", 1, integer = TRUE, minimum = 1))
        events <- events[sort(sample.int(nrow(events), 2000)), , drop = FALSE]
    }
    common_detectors <- intersect(colnames(events), colnames(matrix_values))
    if (length(common_detectors) == 0) return(list(success = FALSE, error = "The sample and matrix share no detector channels."))
    comparison <- lapply(methods, function(method) {
        resolved <- tryCatch(spectreasy:::.normalize_unmix_method(method), error = function(e) e)
        if (inherits(resolved, "error")) return(data.frame(method = method, status = "error", residual_rms = NA_real_, message = conditionMessage(resolved), stringsAsFactors = FALSE))
        result <- tryCatch(
            spectreasy::calc_residuals(
                flowCore::flowFrame(events[, common_detectors, drop = FALSE]),
                matrix_values[, common_detectors, drop = FALSE],
                method = resolved
            ),
            error = function(e) e
        )
        if (inherits(result, "error")) return(data.frame(method = method, status = "error", residual_rms = NA_real_, message = conditionMessage(result), stringsAsFactors = FALSE))
        numeric_result <- suppressWarnings(as.matrix(result))
        data.frame(method = method, status = "complete", residual_rms = sqrt(mean(numeric_result^2, na.rm = TRUE)), message = "", stringsAsFactors = FALSE)
    })
    comparison_df <- do.call(rbind, comparison)
    output_dir <- file.path(root, "spectreasy_outputs", "method_comparison")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    output_file <- file.path(output_dir, "method_comparison.csv")
    utils::write.csv(comparison_df, output_file, row.names = FALSE, quote = TRUE)
    list(success = TRUE, result = list(sample = basename(sample_files[1]), matrix = matrix_file, output_file = output_file, rows = comparison_df))
}

#* Render a report from existing Spectreasy outputs
#* @post /workflow/report
function(req) {
    body <- gui_workflow_body(req)
    report_type <- tolower(gui_workflow_value(body, "report_type", "control"))
    report_format <- tolower(gui_workflow_value(body, "report_format", "html"))
    if (!report_format %in% c("html", "pdf")) {
        return(list(success = FALSE, error = "Report format must be 'html' or 'pdf'."))
    }
    overwrite <- tolower(gui_workflow_value(body, "overwrite", "overwrite"))
    if (!overwrite %in% c("version", "overwrite", "error")) {
        return(list(success = FALSE, error = "Overwrite behavior must be 'version', 'overwrite', or 'error'."))
    }
    root <- gui_workflow_root(body)
    run <- gui_workflow_run(
        body,
        "report",
        if (identical(report_type, "sample")) {
            stop("Sample report rendering requires the current sample results object. Run sample unmixing from the Samples workspace first.")
        } else {
            matrix_file <- file.path(root, "spectreasy_outputs", "unmix_controls", "scc_reference_matrix.csv")
            if (!file.exists(matrix_file)) stop("Reference matrix not found: ", matrix_file)
            M <- spectreasy:::.read_unmixing_matrix_csv(matrix_file)
            report_dir <- file.path(root, "spectreasy_outputs", "unmix_controls", "qc_controls")
            dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
            report_file <- file.path(report_dir, paste0("qc_controls_report.", report_format))
            spectreasy::qc_controls(
                M = M,
                unmixing_matrix_file = matrix_file,
                scc_dir = gui_project_input_path(root, "controls"),
                control_file = file.path(root, "fcs_mapping.csv"),
                output_file = report_file,
                report_format = report_format,
                overwrite = overwrite
            )
        }
    )
    if (!isTRUE(run$success)) return(run)
    report_path <- if (is.list(run$result) && !is.null(run$result$output_file)) run$result$output_file else report_file
    run$result <- list(report_file = report_path, report_format = report_format)
    run
}

#* CORS preflight for workflow_af
#* @options /workflow/af
function(res) {
    return("")
}

gui_default_af_profile_name <- function(fcs_file, existing_names = character()) {
    stem <- tools::file_path_sans_ext(basename(as.character(fcs_file)[1]))
    stem <- gsub("[^A-Za-z0-9._-]+", "_", trimws(stem))
    stem <- gsub("^_+|_+$", "", stem)
    if (!nzchar(stem)) stem <- "af_profile"
    candidate <- stem
    suffix <- 2L
    while (candidate %in% existing_names) {
        candidate <- paste0(stem, "_", suffix)
        suffix <- suffix + 1L
    }
    candidate
}

#* Extract one AF profile through Spectreasy R
#* @post /workflow/af
function(req) {
    body <- gui_workflow_body(req)
    root <- gui_workflow_root(body)
    fcs_file <- gui_workflow_resolve_path(gui_workflow_value(body, "fcs_file", file.path(gui_project_layout(root)$control_input_dir, "unstained_cells.fcs")), root)
    bands <- gui_workflow_number(body, "af_n_bands", 100, integer = TRUE, minimum = 1)
    run <- gui_workflow_run(
        body,
        "af",
        spectreasy::extract_af_profile(
            fcs_file = fcs_file,
            af_n_bands = bands,
            af_max_cells = gui_workflow_number(body, "af_max_cells", 50000, integer = TRUE, minimum = 1),
            seed = gui_workflow_number(body, "seed", 1, integer = TRUE, minimum = 1),
            show_plot = FALSE,
            verbose = FALSE
        )
    )
    if (!isTRUE(run$success)) return(run)
    profile <- run$result
    save_name <- trimws(as.character(gui_workflow_value(body, "save_name", ""))[1])
    save_overwrite <- isTRUE(gui_workflow_bool(body, "save_overwrite", FALSE))
    if (!nzchar(save_name)) {
        existing_profiles <- tryCatch(spectreasy::list_af_profiles(), error = function(e) data.frame())
        existing_names <- if (nrow(existing_profiles) > 0L) existing_profiles$name else character()
        save_name <- gui_default_af_profile_name(fcs_file, existing_names)
    }
    saved_path <- tryCatch(
        spectreasy::save_af_profile(save_name, profile, overwrite = save_overwrite),
        error = function(e) e
    )
    if (inherits(saved_path, "error")) {
        run$success <- FALSE
        run$error <- conditionMessage(saved_path)
        return(run)
    }
    saved_path <- as.character(saved_path)
    run$result <- list(
        bands = nrow(profile$profile),
        detectors = ncol(profile$profile),
        source = fcs_file,
        profile_name = save_name,
        profile_path = saved_path
    )
    run
}
