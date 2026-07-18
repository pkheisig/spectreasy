#' Gate Single-Color Controls Interactively
#'
#' Launches the browser-based manual gating interface for single-color controls.
#' The gate CSV is written to the current working directory by default so it can
#' be reused by [unmix_controls()] or inspected directly.
#'
#' @param scc_dir Directory containing SCC `.fcs` files. Defaults to `"scc"`.
#' @param control_file Control mapping CSV. Defaults to `"fcs_mapping.csv"`.
#' @param gate_file Output gate CSV. Relative paths are resolved from the current
#'   working directory. Defaults to `"ssc_gate_config.csv"`.
#' @param port API port. Defaults to `8000`.
#' @param open_browser Logical. Open browser automatically?
#' @param dev_mode Logical. Use the Vite dev server instead of bundled assets.
#' @param dev_frontend_port Optional Vite frontend port used only with
#'   `dev_mode = TRUE`. When `NULL`, the first free port from 5174 onward is used.
#' @return Invisibly returns the normalized gate CSV path. The function blocks
#'   while the GUI is running.
#' @export
#' @examples
#' if (interactive()) {
#'   gate_controls(scc_dir = "scc", open_browser = FALSE)
#' }
gate_controls <- function(scc_dir = "scc",
                          control_file = "fcs_mapping.csv",
                          gate_file = "ssc_gate_config.csv",
                          port = 8000,
                          open_browser = TRUE,
                          dev_mode = FALSE,
                          dev_frontend_port = NULL) {
    paths <- .normalize_gate_controls_paths(
        scc_dir = scc_dir,
        control_file = control_file,
        gate_file = gate_file
    )
    .launch_spectreasy_gui(
        matrix_dir = getwd(),
        samples_dir = paths$scc_dir,
        port = port,
        open_browser = open_browser,
        dev_mode = dev_mode,
        dev_frontend_port = dev_frontend_port,
        mode = "control-gating",
        gating_scc_dir = paths$scc_dir,
        gating_control_file = paths$control_file,
        gating_gate_file = paths$gate_file
    )
    invisible(paths$gate_file)
}
