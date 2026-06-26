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
    matrix_dir <- normalizePath(matrix_dir, mustWork = TRUE)
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
    message("Starting frontend (npm run dev)...")
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

.launch_spectreasy_gui <- function(matrix_dir = NULL,
                                   samples_dir = NULL,
                                   port = 8000,
                                   open_browser = TRUE,
                                   dev_mode = FALSE,
                                   mode = "tuner",
                                   panel_cytometer = NULL) {
    if (!requireNamespace("plumber", quietly = TRUE)) {
        stop(
            "Package 'plumber' is required for the spectreasy GUI. ",
            "Please install the suggested dependency to use the GUI.",
            call. = FALSE
        )
    }

    paths <- .prepare_gui_paths()
    if (is.null(matrix_dir)) {
        matrix_dir <- .default_adjust_matrix_matrix_dir()
    }
    dirs <- .normalize_gui_dirs(matrix_dir = matrix_dir, samples_dir = samples_dir)

    options(
        spectreasy.matrix_dir = dirs$matrix_dir,
        spectreasy.samples_dir = dirs$samples_dir,
        spectreasy.gui_mode = mode,
        spectreasy.panel_cytometer = panel_cytometer
    )

    frontend <- .resolve_gui_frontend(
        gui_path = paths$gui_path,
        dist_path = paths$dist_path,
        port = port,
        dev_mode = dev_mode
    )

    if (identical(frontend$mode, "dev")) {
        .start_gui_dev_server(
            gui_path = paths$gui_path,
            port = port,
            npm_bin = frontend$npm_bin
        )
    } else {
        message("Using bundled GUI assets from: ", paths$dist_path)
    }

    message("Starting spectreasy API on port ", port)
    message("Matrix directory: ", dirs$matrix_dir)
    message("Samples directory: ", dirs$samples_dir)
    frontend_url <- paste0(frontend$frontend_url, "?mode=", utils::URLencode(mode, reserved = TRUE))

    message("Frontend: ", frontend_url)

    if (isTRUE(open_browser)) utils::browseURL(frontend_url)

    pr <- plumber::plumb(paths$api_path)
    if (!isTRUE(dev_mode)) {
        plumber::pr_static(pr, "/", paths$dist_path)
    }
    pr$run(port = port, host = "127.0.0.1")

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
#' @return Invisibly returns NULL. This function blocks while the API is running.
#' @export
#' @examples
#' if (interactive()) {
#'   adjust_matrix(open_browser = FALSE)
#'   adjust_matrix(matrix_dir = "/path/to/my/matrices", open_browser = FALSE)
#'   adjust_matrix(dev_mode = TRUE, open_browser = FALSE)
#' }
adjust_matrix <- function(matrix_dir = NULL, samples_dir = NULL, port = 8000, open_browser = TRUE, dev_mode = FALSE) {
    .launch_spectreasy_gui(
        matrix_dir = matrix_dir,
        samples_dir = samples_dir,
        port = port,
        open_browser = open_browser,
        dev_mode = dev_mode,
        mode = "tuner"
    )
    invisible(NULL)
}
