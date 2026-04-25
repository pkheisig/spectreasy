# Internal helpers for GUI launch.
.prepare_launch_gui_paths <- function() {
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
        samples_dir <- file.path(matrix_dir, "samples")
    }
    samples_dir <- normalizePath(samples_dir, mustWork = FALSE)

    list(matrix_dir = matrix_dir, samples_dir = samples_dir)
}

.resolve_launch_gui_frontend <- function(gui_path, dist_path, port, dev_mode = FALSE, npm_bin = Sys.which("npm")) {
    frontend_url <- paste0("http://127.0.0.1:", port)

    if (isTRUE(dev_mode)) {
        if (!nzchar(npm_bin)) {
            stop(
                "Developer mode requires npm, but it was not found on PATH. ",
                "Install Node.js/npm or run launch_gui(dev_mode = FALSE) to use bundled assets.",
                call. = FALSE
            )
        }

        node_modules <- file.path(gui_path, "node_modules")
        if (!dir.exists(node_modules)) {
            stop(
                "Developer mode requires GUI dependencies in: ", gui_path, "\n",
                "Run `npm install` in that directory before calling launch_gui(dev_mode = TRUE).",
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

.start_launch_gui_dev_server <- function(gui_path, port, npm_bin) {
    message("Starting frontend (npm run dev)...")
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(gui_path)
    system2(
        npm_bin,
        args = c("run", "dev"),
        env = c(VITE_API_BASE = paste0("http://127.0.0.1:", port)),
        wait = FALSE,
        stdout = FALSE,
        stderr = FALSE
    )
    Sys.sleep(2)

    invisible(NULL)
}

#' Launch spectreasy Interactive Adjustment GUI
#'
#' Starts the backend Plumber API for the interactive matrix adjustment interface.
#' By default, the frontend is served from bundled package assets.
#'
#' @param matrix_dir Directory containing matrix CSV files (default: current working directory)
#' @param samples_dir Directory containing FCS sample files (default: "samples" subfolder of matrix_dir)
#' @param port API port (default: 8000)
#' @param open_browser Logical. Open browser automatically? (default: TRUE)
#' @param dev_mode Logical. If `FALSE` (default), serves bundled GUI assets from the
#'   installed package and requires no Node.js/npm on user machines. If `TRUE`, starts
#'   the Vite dev server (`npm run dev`) from the packaged GUI source folder.
#' @return Invisibly returns NULL. This function blocks while the API is running.
#' @export
#' @examples
#' if (interactive()) {
#'   launch_gui(open_browser = FALSE)
#'   launch_gui(matrix_dir = "/path/to/my/matrices", open_browser = FALSE)
#'   launch_gui(dev_mode = TRUE, open_browser = FALSE)
#' }
launch_gui <- function(matrix_dir = getwd(), samples_dir = NULL, port = 8000, open_browser = TRUE, dev_mode = FALSE) {
    if (!requireNamespace("plumber", quietly = TRUE)) {
        stop(
            "Package 'plumber' is required for launch_gui(). ",
            "Please install the suggested dependency to use the GUI.",
            call. = FALSE
        )
    }

    paths <- .prepare_launch_gui_paths()
    dirs <- .normalize_gui_dirs(matrix_dir = matrix_dir, samples_dir = samples_dir)

    options(
        spectreasy.matrix_dir = dirs$matrix_dir,
        spectreasy.samples_dir = dirs$samples_dir
    )

    frontend <- .resolve_launch_gui_frontend(
        gui_path = paths$gui_path,
        dist_path = paths$dist_path,
        port = port,
        dev_mode = dev_mode
    )

    if (identical(frontend$mode, "dev")) {
        .start_launch_gui_dev_server(
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
    message("Frontend: ", frontend$frontend_url)

    if (isTRUE(open_browser)) utils::browseURL(frontend$frontend_url)

    pr <- plumber::plumb(paths$api_path)
    if (!isTRUE(dev_mode)) {
        plumber::pr_static(pr, "/", paths$dist_path)
    }
    pr$run(port = port, host = "127.0.0.1")

    invisible(NULL)
}
