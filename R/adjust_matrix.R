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
            samples_dir <- file.path(matrix_dir, "samples")
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

    raw_sample_dirs <- file.path(project_dirs, "samples")
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
    file.path(getwd(), "scc")
}

.normalize_gate_controls_paths <- function(scc_dir = "scc",
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

.resolve_gui_frontend <- function(gui_path, dist_path, port, dev_mode = FALSE, npm_bin = Sys.which("npm")) {
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

        return(list(frontend_url = "http://127.0.0.1:5174", npm_bin = npm_bin, mode = "dev"))
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

.start_gui_dev_server <- function(gui_path, port, npm_bin) {
    .spectreasy_console_field("Frontend", "starting npm dev server")
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(gui_path)
    system2(
        npm_bin,
        args = c("run", "dev"),
        env = paste0("VITE_API_BASE=http://127.0.0.1:", port),
        wait = FALSE,
        stdout = FALSE,
        stderr = FALSE
    )
    Sys.sleep(2)

    invisible(NULL)
}

.spectreasy_gui_mode_label <- function(mode) {
    if (identical(mode, "cockpit")) return("Spectreasy Cockpit")
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
                                   unmixing_method = "Spectreasy",
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
        options(
            spectreasy.gating_scc_dir = file.path(project_dir, "scc"),
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
            dev_mode = dev_mode
        )
    }

    if (identical(frontend$mode, "dev")) {
        .start_gui_dev_server(
            gui_path = paths$gui_path,
            port = port,
            npm_bin = frontend$npm_bin
        )
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
        "&token=", utils::URLencode(api_token, reserved = TRUE)
    )

    .message_spectreasy_gui_startup(
        mode = mode,
        port = port,
        frontend_url = sub("&token=[^&]*", "", frontend_url),
        asset_mode = if (identical(frontend$mode, "dev")) "Vite dev server" else if (identical(frontend$mode, "hosted")) "GitHub Pages cockpit" else "bundled package assets",
        gate_file = if (identical(mode, "control-gating")) gate_paths$gate_file else NULL
    )

    if (isTRUE(open_browser)) utils::browseURL(frontend_url)

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
#' @param unmixing_method Method used by the preview unmixing backend. Defaults
#'   to `"Spectreasy"` and accepts the same values as [unmix_samples()].
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
                          unmixing_method = "Spectreasy") {
    .launch_spectreasy_gui(
        matrix_dir = matrix_dir,
        samples_dir = samples_dir,
        port = port,
        open_browser = open_browser,
        dev_mode = dev_mode,
        unmixing_method = unmixing_method,
        mode = "tuner"
    )
    invisible(NULL)
}

#' Launch the Spectreasy Cockpit
#'
#' Starts the local Spectreasy R backend and opens the hosted cockpit. Choose
#' the project folder from the Active project button in the header; the R
#' working directory is not used as an implicit project.
#'
#' This is the recommended entry point for interactive use:
#' `library(spectreasy); spectreasy_gui()`.
#'
#' @return Invisibly returns `NULL`. This function blocks while the local GUI
#' application is running.
#' @export
#' @examples
#' if (interactive()) {
#'   spectreasy_gui()
#' }
spectreasy_gui <- function() {
    project_dir <- file.path(tempdir(), "spectreasy-cockpit-session")
    dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
    samples_dir <- file.path(project_dir, "samples")
    .launch_spectreasy_gui(
        matrix_dir = project_dir,
        samples_dir = samples_dir,
        mode = "cockpit",
        hosted_frontend_url = getOption(
            "spectreasy.cockpit_url",
            "https://pkheisig.github.io/spectreasy/"
        ),
        initial_project_selected = FALSE
    )
    invisible(NULL)
}
