#' Create spectreasy Control File
#'
#' Generates a spectreasy-compatible control CSV.
#'
#' @param input_folder Directory containing single-stained control FCS files.
#' @param cytometer Cytometer name used as a channel-mapping hint. The default,
#'   `"auto"`, infers the cytometer from FCS detector names when possible.
#' @param unknown_fluor_policy How to fill unresolved fluorophores:
#'   `"empty"` (recommended), `"by_channel"` (best-effort guess), `"filename"`.
#' @param output_file Path where the CSV will be saved (default: "fcs_mapping.csv").
#' @param custom_fluorophores Optional named vector to map filenames to fluorophore names.
#' @return A data frame containing the control file information.
#'
#' The generated file contains only columns that usually need review:
#' `filename`, `fluorophore`, `marker`, `channel`, `control.type`, and
#' `is.viability`. It auto-detects `control.type` from filename tokens
#' (for example `"beads"` or `"cells"`). Unstained/negative bead files are
#' denoted as `AF_beads` with `control.type = "beads"` so they can be used
#' as bead-background negatives. Ordinary unstained cell controls are numbered
#' `AF`, `AF_2`, `AF_3`, and so on, and are pooled as AF-bank sources.
#' Empty `universal.negative` values are filled with the exact matching negative
#' control filename: the primary AF file for ordinary cell SCCs, the `AF_dead`
#' file for viability SCCs, and the `AF_beads` file for bead SCCs. Legacy `AF`
#' references are replaced by the primary AF filename; explicit filename
#' mappings are preserved.
#' @export
#' @examples
#' make_example_ff <- function(main, n = 250) {
#'   exprs <- cbind(
#'     "B2-A" = pmax(rnorm(n, main[1], 40), 1),
#'     "YG1-A" = pmax(rnorm(n, main[2], 15), 1),
#'     "R1-A" = pmax(rnorm(n, main[3], 10), 1),
#'     "FSC-A" = rnorm(n, 90000, 7000),
#'     "SSC-A" = rnorm(n, 45000, 5000),
#'     Time = seq_len(n)
#'   )
#'   flowCore::flowFrame(exprs)
#' }
#'
#' td <- tempfile("spectreasy-")
#' dir.create(td)
#' scc_dir <- file.path(td, "scc")
#' dir.create(scc_dir)
#' flowCore::write.FCS(make_example_ff(c(800, 80, 50)), file.path(scc_dir, "FITC_cells.fcs"))
#' flowCore::write.FCS(make_example_ff(c(80, 820, 60)), file.path(scc_dir, "PE_cells.fcs"))
#'
#' control_df <- create_control_file(
#'   input_folder = scc_dir,
#'   cytometer = "auto",
#'   output_file = file.path(td, "fcs_mapping.csv")
#' )
#' head(control_df)
create_control_file <- function(input_folder = "scc",
                                cytometer = "auto",
                                unknown_fluor_policy = c("empty", "by_channel", "filename"),
                                output_file = "fcs_mapping.csv",
                                custom_fluorophores = NULL) {
    unknown_fluor_policy <- .match_arg_ci(
        unknown_fluor_policy,
        c("empty", "by_channel", "filename"),
        "unknown_fluor_policy"
    )
    if (!dir.exists(input_folder)) {
        .spectreasy_stop_missing_directory(input_folder, label = "input_folder")
    }
    scc_files <- list.files(input_folder, pattern = "\\.fcs$", full.names = FALSE, ignore.case = TRUE)

    if (length(scc_files) == 0) {
        .spectreasy_stop_empty_fcs_directory(input_folder, label = "input_folder")
    }

    cytometer_resolved <- .resolve_cytometer_from_files(
        cytometer,
        files = file.path(input_folder, scc_files)
    )
    ref <- .prepare_control_file_reference(cytometer_resolved)
    df <- .build_control_file_scc_df(scc_files, ref = ref, custom_fluorophores = custom_fluorophores)

    df <- .annotate_control_file_rows(
        df = df,
        input_folder = input_folder,
        cytometer = cytometer_resolved,
        unknown_fluor_policy = unknown_fluor_policy,
        ref = ref
    )
    df <- .finalize_control_file_df(
        df = df,
        scc_files = scc_files
    )

    utils::write.csv(df, output_file, row.names = FALSE, quote = TRUE)
    df
}

#' Get Spectra via Internal Backend
#'
#' Extracts SCC signatures using internal logic. AF/unstained controls should be
#' mapped as files inside `control_dir`.
#'
#' @param flow_frame A `flowFrame` used to determine detector ordering.
#' @param control_file Path to spectreasy-compatible control CSV.
#' @param control_dir Directory containing SCC FCS files.
#' @param cytometer Cytometer name used as a channel-mapping hint. The default,
#'   `"auto"`, infers the cytometer from FCS detector names when possible.
#' @return Expanded reference matrix aligned to the detectors in `flow_frame`.
#' @export
#' @examples
#' if (interactive()) {
#'   ff <- flowCore::read.FCS("samples/Sample1.fcs", transformation = FALSE)
#'   M_ctrl <- get_control_spectra(
#'     flow_frame = ff,
#'     control_file = "fcs_mapping.csv",
#'     control_dir = "scc",
#'     cytometer = "auto"
#'   )
#'   dim(M_ctrl)
#' }
get_control_spectra <- function(flow_frame,
                                control_file = "fcs_mapping.csv",
                                control_dir = "scc",
                                cytometer = "auto") {
    # 1. Get detector info
    pd <- flowCore::pData(flowCore::parameters(flow_frame))
    det_info <- get_sorted_detectors(pd)
    detector_names <- det_info$names
    cytometer_resolved <- .resolve_cytometer_from_pd(cytometer, pd)

    control_file <- .resolve_control_file_path(control_file)

    .spectreasy_console_step("Reference matrix", "from single-color controls")
    control_df <- utils::read.csv(control_file, stringsAsFactors = FALSE, check.names = FALSE)
    M_scc <- build_reference_matrix(input_folder = control_dir, control_df = control_df, cytometer = cytometer_resolved)
    M_scc[, detector_names, drop = FALSE]
}
