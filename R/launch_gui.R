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
#'   package and requires no Node.js/npm on user machines. If `TRUE`, starts the Vite
#'   dev server (`npm run dev`) from the GUI source folder.
#' @return Invisibly returns NULL. This function blocks while the API is running.
#' @export
#' @examples
#' \dontrun{
#' # Start API + bundled GUI (recommended for users; no npm required)
#' launch_gui()
#'
#' # Start with custom matrix directory
#' launch_gui(matrix_dir = "/path/to/my/matrices")
#'
#' # Developer mode (uses Vite)
#' launch_gui(dev_mode = TRUE)
#' }
launch_gui <- function(matrix_dir = getwd(), samples_dir = NULL, port = 8000, open_browser = TRUE, dev_mode = FALSE) {
    api_path <- system.file("api/gui_api.R", package = "spectreasy")
    gui_path <- system.file("gui", package = "spectreasy")

    if (api_path == "") {
        api_path <- file.path(getwd(), "inst", "api", "gui_api.R")
    }
    if (gui_path == "" || !dir.exists(gui_path)) {
        gui_path <- file.path(getwd(), "inst", "gui")
    }
    if (!dir.exists(gui_path)) {
        gui_path <- file.path(getwd(), "gui")
    }

    if (!file.exists(api_path)) stop("Could not find gui_api.R")
    if (!dir.exists(gui_path)) stop("Could not find gui folder")

    dist_path <- file.path(gui_path, "dist")

    matrix_dir <- normalizePath(matrix_dir, mustWork = TRUE)
    if (is.null(samples_dir)) samples_dir <- file.path(matrix_dir, "samples")
    samples_dir <- normalizePath(samples_dir, mustWork = FALSE)

    options(
        spectreasy.matrix_dir = matrix_dir,
        spectreasy.samples_dir = samples_dir
    )

    frontend_url <- paste0("http://localhost:", port)
    if (isTRUE(dev_mode)) {
        node_modules <- file.path(gui_path, "node_modules")
        if (!dir.exists(node_modules)) {
            stop("node_modules not found. Run this once in terminal:\n  cd ", gui_path, " && npm install")
        }
        message("Starting frontend (npm run dev)...")
        old_wd <- getwd()
        on.exit(setwd(old_wd), add = TRUE)
        setwd(gui_path)
        system2(
            "npm",
            args = c("run", "dev"),
            env = c(VITE_API_BASE = paste0("http://localhost:", port)),
            wait = FALSE,
            stdout = FALSE,
            stderr = FALSE
        )
        Sys.sleep(2)
        frontend_url <- "http://localhost:5174"
    } else {
        if (!dir.exists(dist_path) || !file.exists(file.path(dist_path, "index.html"))) {
            stop(
                "Bundled GUI assets not found at: ", dist_path, "\n",
                "Reinstall spectreasy from a build that includes inst/gui/dist."
            )
        }
        message("Using bundled GUI assets from: ", dist_path)
    }

    message("Starting spectreasy API on port ", port)
    message("Matrix directory: ", matrix_dir)
    message("Samples directory: ", samples_dir)
    message("Frontend: ", frontend_url)

    if (open_browser) utils::browseURL(frontend_url)

    pr <- plumber::plumb(api_path)
    if (!isTRUE(dev_mode)) {
        plumber::pr_static(pr, "/", dist_path)
    }
    pr$run(port = port, host = "0.0.0.0")

    invisible(NULL)
}
