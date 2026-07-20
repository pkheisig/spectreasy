#' Launch Spectral Panel Builder
#'
#' Opens an interactive browser-based spectral panel builder for packaged
#' theoretical spectra. The panel builder currently supports Aurora, Discover,
#' ID7000, and Xenith libraries.
#'
#' @param port API port (default: 8000).
#' @param open_browser Logical. Open browser automatically?
#' @param dev_mode Logical. Use the Vite dev server instead of bundled assets.
#' @param dev_frontend_port Optional Vite frontend port used only with
#'   `dev_mode = TRUE`. When `NULL`, the first free port from 5174 onward is used.
#' @return Invisibly returns NULL. This function blocks while the API is running.
#' @export
#' @examples
#' if (interactive()) {
#'   build_panel(open_browser = FALSE)
#' }
build_panel <- function(port = 8000,
                        open_browser = TRUE,
                        dev_mode = FALSE,
                        dev_frontend_port = NULL) {
    .launch_spectreasy_gui(
        matrix_dir = NULL,
        samples_dir = NULL,
        port = port,
        open_browser = open_browser,
        dev_mode = dev_mode,
        dev_frontend_port = dev_frontend_port,
        mode = "panel-builder",
        panel_cytometer = "aurora"
    )
    invisible(NULL)
}
