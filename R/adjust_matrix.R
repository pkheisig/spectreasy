# Internal helpers for GUI launch.
.prepare_gui_paths <- function() {
    api_path <- system.file("api", "gui_api.R", package = "spectreasy")
    gui_path <- system.file("gui", package = "spectreasy")

    if (!nzchar(api_path) || !file.exists(api_path)) {
        stop(
            "Could not find the bundled GUI API script in the installed spectreasy package.",
            call. = FALSE
        )
    }
    if (!nzchar(gui_path) || !dir.exists(gui_path)) {
        stop(
            "Could not find the bundled GUI assets in the installed spectreasy package.",
            call. = FALSE
        )
    }

    list(
        api_path = api_path,
        gui_path = gui_path,
        dist_path = file.path(gui_path, "dist")
    )
}

.normalize_gui_dirs <- function(matrix_dir, samples_dir = NULL) {
    if (!dir.exists(matrix_dir)) {
        .spectreasy_stop_missing_directory(matrix_dir, label = "matrix_dir")
    }
    matrix_dir <- normalizePath(matrix_dir, mustWork = FALSE)
    if (is.null(samples_dir)) {
        samples_dir <- .default_adjust_matrix_samples_dir(matrix_dir = matrix_dir)
        if (is.null(samples_dir)) {
            samples_dir <- .spectreasy_project_input_path(matrix_dir, "samples")
        }
    }
    samples_dir <- normalizePath(samples_dir, mustWork = FALSE)

    list(matrix_dir = matrix_dir, samples_dir = samples_dir)
}

.default_adjust_matrix_matrix_dir <- function() {
    unmix_controls_dir <- file.path(getwd(), "spectreasy_outputs", "unmix_controls")
    if (dir.exists(unmix_controls_dir)) {
        return(unmix_controls_dir)
    }
    getwd()
}

.infer_project_dir_from_matrix_dir <- function(matrix_dir) {
    if (is.null(matrix_dir) || !nzchar(matrix_dir)) {
        return(NULL)
    }
    matrix_dir <- normalizePath(matrix_dir, mustWork = FALSE)
    if (basename(matrix_dir) == "unmix_controls" &&
        basename(dirname(matrix_dir)) == "spectreasy_outputs") {
        return(dirname(dirname(matrix_dir)))
    }
    dirname(matrix_dir)
}

.gui_project_dir <- function(matrix_dir, mode) {
    if (identical(mode, "cockpit")) {
        return(normalizePath(matrix_dir, mustWork = FALSE))
    }
    .infer_project_dir_from_matrix_dir(matrix_dir)
}

.default_adjust_matrix_samples_dir <- function(matrix_dir = NULL) {
    project_dirs <- unique(c(getwd(), .infer_project_dir_from_matrix_dir(matrix_dir)))
    project_dirs <- project_dirs[nzchar(project_dirs)]

    raw_sample_dirs <- vapply(project_dirs, function(project_dir) {
        if (!dir.exists(project_dir)) return("")
        layout <- .spectreasy_project_layout(project_dir, ensure_markers = FALSE, persist = FALSE)
        file.path(project_dir, layout$sample_input_dir)
    }, character(1))
    raw_sample_dirs <- unname(raw_sample_dirs[nzchar(raw_sample_dirs)])
    raw_hit <- raw_sample_dirs[dir.exists(raw_sample_dirs)][1]
    if (!is.na(raw_hit) && nzchar(raw_hit)) {
        return(raw_hit)
    }

    unmixed_fcs_dirs <- file.path(project_dirs, "spectreasy_outputs", "unmix_samples", "unmixed_fcs")
    unmixed_fcs_hit <- unmixed_fcs_dirs[dir.exists(unmixed_fcs_dirs)][1]
    if (!is.na(unmixed_fcs_hit) && nzchar(unmixed_fcs_hit)) {
        return(unmixed_fcs_hit)
    }

    unmix_samples_dirs <- file.path(project_dirs, "spectreasy_outputs", "unmix_samples")
    for (unmix_samples_dir in unmix_samples_dirs[dir.exists(unmix_samples_dirs)]) {
        if (length(list.files(unmix_samples_dir, pattern = "\\.fcs$", ignore.case = TRUE)) > 0) {
            return(unmix_samples_dir)
        }
    }
    NULL
}

.default_gate_controls_scc_dir <- function() {
    layout <- .spectreasy_project_layout(getwd(), ensure_markers = FALSE, persist = FALSE)
    file.path(getwd(), layout$control_input_dir)
}

.normalize_gate_controls_paths <- function(scc_dir = .default_gate_controls_scc_dir(),
                                           control_file = "fcs_mapping.csv",
                                           gate_file = "ssc_gate_config.csv") {
    if (!dir.exists(scc_dir)) {
        .spectreasy_stop_missing_directory(scc_dir, label = "scc_dir")
    }
    scc_dir <- normalizePath(scc_dir, mustWork = FALSE)
    control_file <- as.character(control_file)[1]
    if (!is.na(control_file) && nzchar(trimws(control_file)) && !file.exists(control_file)) {
        candidate <- file.path(dirname(scc_dir), basename(control_file))
        if (file.exists(candidate)) {
            control_file <- candidate
        }
    }
    control_file <- normalizePath(control_file, mustWork = FALSE)

    gate_file <- as.character(gate_file)[1]
    if (is.na(gate_file) || !nzchar(trimws(gate_file))) {
        gate_file <- "ssc_gate_config.csv"
    }
    if (!grepl("\\.csv$", gate_file, ignore.case = TRUE)) {
        gate_file <- paste0(gate_file, ".csv")
    }
    if (!grepl("^(/|[A-Za-z]:)", gate_file)) {
        gate_file <- file.path(getwd(), gate_file)
    }
    gate_file <- normalizePath(gate_file, mustWork = FALSE)

    list(scc_dir = scc_dir, control_file = control_file, gate_file = gate_file)
}

.normalize_gui_port <- function(port, arg_name = "port") {
    value <- suppressWarnings(as.integer(port))
    if (length(port) != 1L || is.na(value) || value < 1L || value > 65535L || !identical(as.numeric(port), as.numeric(value))) {
        stop(arg_name, " must be one integer between 1 and 65535.", call. = FALSE)
    }
    value
}

.gui_tcp_port_open <- function(port, host = "127.0.0.1") {
    connection <- tryCatch(
        suppressWarnings(socketConnection(
            host = host,
            port = .normalize_gui_port(port),
            open = "r+b",
            blocking = TRUE,
            timeout = 0.2
        )),
        error = function(e) NULL
    )
    if (is.null(connection)) return(FALSE)
    close(connection)
    TRUE
}

.gui_port_is_available <- function(port) {
    if (.gui_tcp_port_open(port)) return(FALSE)
    socket <- tryCatch(serverSocket(.normalize_gui_port(port)), error = function(e) NULL)
    if (is.null(socket)) return(FALSE)
    close(socket)
    TRUE
}

.resolve_gui_dev_port <- function(dev_frontend_port = NULL, preferred = 5174L) {
    if (!is.null(dev_frontend_port)) {
        port <- .normalize_gui_port(dev_frontend_port, "dev_frontend_port")
        if (!.gui_port_is_available(port)) {
            stop("The requested Vite frontend port is already in use: ", port, call. = FALSE)
        }
        return(port)
    }

    candidates <- seq.int(.normalize_gui_port(preferred), min(65535L, preferred + 100L))
    available <- candidates[vapply(candidates, .gui_port_is_available, logical(1))]
    if (length(available) == 0L) {
        stop("Could not find an available Vite frontend port between ", preferred, " and ", max(candidates), ".", call. = FALSE)
    }
    available[[1]]
}

.resolve_gui_frontend <- function(gui_path,
                                  dist_path,
                                  port,
                                  dev_mode = FALSE,
                                  npm_bin = Sys.which("npm"),
                                  dev_frontend_port = NULL) {
    frontend_url <- paste0("http://127.0.0.1:", port)

    if (isTRUE(dev_mode)) {
        if (!nzchar(npm_bin)) {
            stop(
                "Developer mode requires npm, but it was not found on PATH. ",
                "Install Node.js/npm or run adjust_matrix(dev_mode = FALSE) to use bundled assets.",
                call. = FALSE
            )
        }

        node_modules <- file.path(gui_path, "node_modules")
        if (!dir.exists(node_modules)) {
            stop(
                "Developer mode requires GUI dependencies in: ", gui_path, "\n",
                "Run `npm install` in that directory before calling adjust_matrix(dev_mode = TRUE).",
                call. = FALSE
            )
        }

        frontend_port <- .resolve_gui_dev_port(dev_frontend_port)
        return(list(
            frontend_url = paste0("http://127.0.0.1:", frontend_port),
            frontend_port = frontend_port,
            npm_bin = npm_bin,
            mode = "dev"
        ))
    }

    if (!dir.exists(dist_path) || !file.exists(file.path(dist_path, "index.html"))) {
        stop(
            "Bundled GUI assets not found at: ", dist_path, "\n",
            "Reinstall spectreasy from a build that includes inst/gui/dist.",
            call. = FALSE
        )
    }

    list(frontend_url = frontend_url, npm_bin = npm_bin, mode = "bundled")
}

.spawn_gui_dev_process <- function(gui_path, api_port, frontend_port, npm_bin) {
    if (!requireNamespace("processx", quietly = TRUE)) {
        stop(
            "Developer mode requires package 'processx' to manage the Vite server safely.",
            call. = FALSE
        )
    }
    child_env <- Sys.getenv()
    child_env[["VITE_API_BASE"]] <- paste0("http://127.0.0.1:", api_port)
    processx::process$new(
        npm_bin,
        c(
            "run", "dev", "--",
            "--host", "127.0.0.1",
            "--port", as.character(frontend_port),
            "--strictPort"
        ),
        wd = gui_path,
        env = child_env,
        stdout = "|",
        stderr = "2>&1",
        cleanup = TRUE,
        cleanup_tree = TRUE
    )
}

.wait_for_gui_frontend <- function(process, frontend_port, timeout = 15) {
    deadline <- Sys.time() + timeout
    repeat {
        if (!isTRUE(process$is_alive())) {
            output <- paste(process$read_all_output_lines(), collapse = "\n")
            stop(
                "The Vite frontend stopped before becoming ready.",
                if (nzchar(output)) paste0("\n", output) else "",
                call. = FALSE
            )
        }
        if (.gui_tcp_port_open(frontend_port)) return(invisible(TRUE))
        if (Sys.time() >= deadline) {
            process$kill_tree()
            stop("Timed out waiting for the Vite frontend on port ", frontend_port, ".", call. = FALSE)
        }
        Sys.sleep(0.1)
    }
}

.start_gui_dev_server <- function(gui_path,
                                  api_port,
                                  frontend_port,
                                  npm_bin,
                                  spawn = .spawn_gui_dev_process,
                                  wait_until_ready = .wait_for_gui_frontend) {
    .spectreasy_console_field("Frontend", paste0("starting npm dev server on port ", frontend_port))
    process <- spawn(gui_path, api_port, frontend_port, npm_bin)
    tryCatch(
        wait_until_ready(process, frontend_port),
        error = function(e) {
            if (isTRUE(process$is_alive())) process$kill_tree()
            stop(conditionMessage(e), call. = FALSE)
        }
    )
    process
}

.spectreasy_gui_mode_label <- function(mode) {
    if (identical(mode, "cockpit")) return("Cockpit")
    if (identical(mode, "control-gating")) return("Control gating GUI")
    if (identical(mode, "panel-builder")) return("Spectral panel builder GUI")
    "Matrix adjustment GUI"
}

.spectreasy_gui_display_path <- function(path) {
    .spectreasy_console_path(path)
}

.message_spectreasy_gui_startup <- function(mode,
                                            port,
                                            frontend_url,
                                            asset_mode = "bundled",
                                            gate_file = NULL) {
    .spectreasy_console_header(.spectreasy_gui_mode_label(mode))
    .spectreasy_console_field("Frontend", frontend_url)
    .spectreasy_console_field("API port", port)
    .spectreasy_console_field("Assets", asset_mode)
    gate_file_display <- .spectreasy_gui_display_path(gate_file)
    if (nzchar(gate_file_display)) {
        .spectreasy_console_field("Gate CSV", gate_file_display)
    }
    .spectreasy_console_field("INFO", "Press Ctrl + C in this R console to terminate the GUI app.")
    .spectreasy_console_footer()
}

.run_gui_until_shutdown <- function(pr, port, host = "127.0.0.1", quiet = FALSE, announce = TRUE) {
    options(spectreasy.gui_shutdown_requested = FALSE)
    server <- httpuv::startServer(host = host, port = port, app = pr, quiet = quiet)
    options(spectreasy.gui_server = server)
    on.exit({
        try(httpuv::stopServer(server), silent = TRUE)
        options(
            spectreasy.gui_server = NULL,
            spectreasy.gui_shutdown_requested = NULL
        )
    }, add = TRUE)

    if (isTRUE(announce)) {
        .spectreasy_console_field("API", paste0("http://", host, ":", port))
    }
    while (!isTRUE(getOption("spectreasy.gui_shutdown_requested", FALSE))) {
        httpuv::service(timeout = 250)
    }
    invisible(NULL)
}

.gui_session_token <- function(n = 32L) {
    bytes <- if (requireNamespace("openssl", quietly = TRUE)) {
        openssl::rand_bytes(n)
    } else {
        tryCatch(readBin("/dev/urandom", what = "raw", n = n), error = function(e) raw())
    }
    if (length(bytes) != n) {
        bytes <- as.raw(sample.int(256L, n, replace = TRUE) - 1L)
    }
    paste(sprintf("%02x", as.integer(bytes)), collapse = "")
}

.gui_url_origin <- function(url) {
    origin <- sub("^(https?://[^/]+).*$", "\\1", url, ignore.case = TRUE)
    if (!grepl("^https?://", origin, ignore.case = TRUE)) "" else origin
}

.launch_spectreasy_gui <- function(matrix_dir = NULL,
                                   samples_dir = NULL,
                                   port = 8000,
                                   open_browser = TRUE,
                                   dev_mode = FALSE,
                                   dev_frontend_port = NULL,
                                   unmixing_method = "AutoSpectral",
                                   mode = "tuner",
                                   panel_cytometer = NULL,
                                   gating_scc_dir = NULL,
                                   gating_control_file = NULL,
                                   gating_gate_file = NULL,
                                   hosted_frontend_url = NULL,
                                   initial_project_selected = TRUE) {
    if (!requireNamespace("plumber", quietly = TRUE)) {
        stop(
            "Package 'plumber' is required for the spectreasy GUI. ",
            "Please install the suggested dependency to use the GUI.",
            call. = FALSE
        )
    }
    if (!requireNamespace("httpuv", quietly = TRUE) || !requireNamespace("later", quietly = TRUE)) {
        stop(
            "Packages 'httpuv' and 'later' are required for the spectreasy GUI. ",
            "Please install the suggested GUI dependencies.",
            call. = FALSE
        )
    }
    unmixing_method <- .normalize_unmix_method(unmixing_method)

    paths <- .prepare_gui_paths()
    if (identical(mode, "control-gating")) {
        gate_paths <- .normalize_gate_controls_paths(
            scc_dir = if (is.null(gating_scc_dir)) .default_gate_controls_scc_dir() else gating_scc_dir,
            control_file = if (is.null(gating_control_file)) "fcs_mapping.csv" else gating_control_file,
            gate_file = if (is.null(gating_gate_file)) "ssc_gate_config.csv" else gating_gate_file
        )
        dirs <- list(matrix_dir = getwd(), samples_dir = gate_paths$scc_dir)
    } else if (is.null(matrix_dir)) {
        matrix_dir <- .default_adjust_matrix_matrix_dir()
        dirs <- .normalize_gui_dirs(matrix_dir = matrix_dir, samples_dir = samples_dir)
    } else {
        dirs <- .normalize_gui_dirs(matrix_dir = matrix_dir, samples_dir = samples_dir)
    }

    options(
        spectreasy.project_dir = dirs$matrix_dir,
        spectreasy.matrix_dir = dirs$matrix_dir,
        spectreasy.samples_dir = dirs$samples_dir,
        spectreasy.unmixing_method = unmixing_method,
        spectreasy.gui_mode = mode,
        spectreasy.panel_cytometer = panel_cytometer,
        spectreasy.project_selected = isTRUE(initial_project_selected)
    )
    if (!identical(mode, "control-gating")) {
        project_dir <- .gui_project_dir(dirs$matrix_dir, mode = mode)
        if (is.null(project_dir) || !dir.exists(project_dir)) project_dir <- dirs$matrix_dir
        project_layout <- .spectreasy_project_layout(project_dir)
        options(
            spectreasy.samples_dir = file.path(project_dir, project_layout$sample_input_dir),
            spectreasy.gating_scc_dir = file.path(project_dir, project_layout$control_input_dir),
            spectreasy.gating_control_file = file.path(project_dir, "fcs_mapping.csv"),
            spectreasy.gating_gate_file = file.path(project_dir, "ssc_gate_config.csv")
        )
    }
    if (identical(mode, "control-gating")) {
        options(
            spectreasy.gating_scc_dir = gate_paths$scc_dir,
            spectreasy.gating_control_file = gate_paths$control_file,
            spectreasy.gating_gate_file = gate_paths$gate_file,
            spectreasy.gating_payload_cache = list(),
            spectreasy.gating_spectrum_cache = list(),
            spectreasy.gating_detector_cache = list()
        )
    }

    frontend <- if (!is.null(hosted_frontend_url) && nzchar(trimws(hosted_frontend_url))) {
        list(mode = "hosted", frontend_url = sub("/+$", "", trimws(hosted_frontend_url)))
    } else {
        .resolve_gui_frontend(
            gui_path = paths$gui_path,
            dist_path = paths$dist_path,
            port = port,
            dev_mode = dev_mode,
            dev_frontend_port = dev_frontend_port
        )
    }

    dev_server <- NULL
    if (identical(frontend$mode, "dev")) {
        dev_server <- .start_gui_dev_server(
            gui_path = paths$gui_path,
            api_port = port,
            frontend_port = frontend$frontend_port,
            npm_bin = frontend$npm_bin
        )
        on.exit({
            if (isTRUE(dev_server$is_alive())) dev_server$kill_tree()
        }, add = TRUE)
    }
    api_token <- .gui_session_token()
    options(
        spectreasy.gui_api_token = api_token,
        spectreasy.gui_allowed_origins = .gui_url_origin(frontend$frontend_url)
    )
    query_separator <- if (grepl("?", frontend$frontend_url, fixed = TRUE)) "&" else "?"
    frontend_url <- paste0(
        frontend$frontend_url,
        query_separator,
        "mode=", utils::URLencode(mode, reserved = TRUE),
        "&api=", utils::URLencode(paste0("http://127.0.0.1:", port), reserved = TRUE),
        "#token=", utils::URLencode(api_token, reserved = TRUE)
    )

    .message_spectreasy_gui_startup(
        mode = mode,
        port = port,
        frontend_url = sub("#token=[^&]*", "", frontend_url),
        asset_mode = if (identical(frontend$mode, "dev")) "Vite dev server" else if (identical(frontend$mode, "hosted")) "GitHub Pages cockpit" else "bundled package assets",
        gate_file = if (identical(mode, "control-gating")) gate_paths$gate_file else NULL
    )

    pr <- plumber::plumb(paths$api_path)
    pr <- plumber::pr_filter(pr, "nocache", function(req, res) {
        res$setHeader("Cache-Control", "no-store, no-cache, must-revalidate, max-age=0")
        res$setHeader("Pragma", "no-cache")
        res$setHeader("Expires", "0")
        plumber::forward()
    })
    if (!isTRUE(dev_mode)) {
        plumber::pr_static(pr, "/", paths$dist_path)
    }
    if (isTRUE(open_browser)) {
        local({
            browser_url <- frontend_url
            later::later(function() utils::browseURL(browser_url), delay = 0.15)
        })
    }
    if (identical(mode, "control-gating")) {
        .run_gui_until_shutdown(pr, port = port, host = "127.0.0.1", announce = FALSE)
    } else {
        pr$run(port = port, host = "127.0.0.1", docs = FALSE, quiet = TRUE)
    }

    invisible(NULL)
}

#' Adjust an Unmixing Matrix Interactively
#'
#' Starts the backend Plumber API for the interactive matrix adjustment interface.
#' By default, the frontend is served from bundled package assets.
#'
#' @param matrix_dir Directory containing matrix CSV files. If `NULL`, defaults to
#'   `spectreasy_outputs/unmix_controls` when it exists under the current working
#'   directory, otherwise the current working directory.
#' @param samples_dir Directory containing FCS sample files. If `NULL`, defaults to
#'   the `samples` folder under the current working directory or inferred project
#'   directory when it exists, then falls back to
#'   `spectreasy_outputs/unmix_samples/unmixed_fcs` or
#'   `spectreasy_outputs/unmix_samples` when available. Otherwise, uses the
#'   `samples` subfolder of `matrix_dir`.
#' @param port API port (default: 8000)
#' @param open_browser Logical. Open browser automatically? (default: TRUE)
#' @param dev_mode Logical. If `FALSE` (default), serves bundled GUI assets from the
#'   installed package and requires no Node.js/npm on user machines. If `TRUE`, starts
#'   the Vite dev server (`npm run dev`) from the packaged GUI source folder.
#' @param dev_frontend_port Optional Vite frontend port used only with
#'   `dev_mode = TRUE`. When `NULL`, the first free port from 5174 onward is used.
#' @param unmixing_method Method used by the preview unmixing backend. Defaults
#'   to `"AutoSpectral"` and accepts the same values as [unmix_samples()].
#' @return Invisibly returns NULL. This function blocks while the API is running.
#' @export
#' @examples
#' if (interactive()) {
#'   adjust_matrix(open_browser = FALSE)
#'   adjust_matrix(matrix_dir = "/path/to/my/matrices", open_browser = FALSE)
#'   adjust_matrix(dev_mode = TRUE, open_browser = FALSE)
#' }
adjust_matrix <- function(matrix_dir = NULL,
                          samples_dir = NULL,
                          port = 8000,
                          open_browser = TRUE,
                          dev_mode = FALSE,
                          dev_frontend_port = NULL,
                          unmixing_method = "AutoSpectral") {
    .launch_spectreasy_gui(
        matrix_dir = matrix_dir,
        samples_dir = samples_dir,
        port = port,
        open_browser = open_browser,
        dev_mode = dev_mode,
        dev_frontend_port = dev_frontend_port,
        unmixing_method = unmixing_method,
        mode = "tuner"
    )
    invisible(NULL)
}

#' Launch the Spectreasy Cockpit
#'
#' Starts the local Spectreasy R backend and serves the bundled cockpit from
#' that same localhost process. The public GitHub Pages site remains the
#' installation and onboarding entry point; using localhost here avoids browser
#' local-network permission and cross-origin restrictions. The
#' current R working directory is opened as the initial project. Use the
#' project button in the header to create or open a different project.
#'
#' This is the recommended entry point for interactive use:
#' `library(spectreasy); spectreasy_gui()`.
#'
#' @param port Local API and cockpit port. Defaults to `8000`. Use another
#'   available port, for example `8001`, to run a separate test instance.
#' @return Invisibly returns `NULL`. This function blocks while the local GUI
#' application is running.
#' @export
#' @examples
#' if (interactive()) {
#'   spectreasy_gui()
#'   spectreasy_gui(port = 8001)
#' }
spectreasy_gui <- function(port = 8000) {
    project_dir <- normalizePath(getwd(), mustWork = TRUE)
    samples_dir <- .spectreasy_project_input_path(project_dir, "samples")
    .launch_spectreasy_gui(
        matrix_dir = project_dir,
        samples_dir = samples_dir,
        port = port,
        mode = "cockpit",
        hosted_frontend_url = NULL,
        initial_project_selected = TRUE
    )
    invisible(NULL)
}
