#' Auto-Unmix Single-Color Controls
#'
#' SCC workflow helper that:
#' 1) validates/creates the control file,
#' 2) builds the reference matrix from SCCs,
#' 3) unmixed SCC files,
#' 4) saves reference/unmixing matrices,
#' 5) plots spectra, unmixing matrix, and SCC unmixing scatter matrix.
#'
#' This function is intended as the control-stage step before GUI adjustment
#' and downstream sample unmixing.
#'
#' @param scc_dir Directory containing SCC FCS files.
#' @param control_df Optional control mapping data.frame, or a path to control CSV.
#' @param control_file Path to control mapping CSV used when `control_df` is `NULL`.
#' @param auto_create_control Logical; auto-generate control file when missing.
#' @param cytometer Cytometer name (for example `"Aurora"`).
#' @param auto_default_control_type Deprecated and ignored.
#' @param auto_unknown_fluor_policy Auto-fill policy for unresolved fluorophores
#'   when creating controls (`"by_channel"`, `"empty"`, `"filename"`).
#' @param output_dir Output directory for SCC workflow artifacts.
#' @param unmix_method SCC unmixing method (`"WLS"`, `"OLS"`, `"NNLS"`, `"AutoSpectral"`).
#' @param build_qc_plots Logical; keep detailed build_reference_matrix plots.
#' @param unmix_scatter_panel_size_mm Panel size for SCC unmixing scatter matrix plot.
#' @param ... Additional arguments forwarded to [build_reference_matrix()].
#'
#' @return List with `M`, `W`, `unmixed_list`, and key output file paths.
#' @export
#' @examples
#' \dontrun{
#' ctrl <- autounmix_controls(
#'   scc_dir = "scc",
#'   control_file = "fcs_control_file.csv",
#'   auto_create_control = TRUE,
#'   cytometer = "Aurora",
#'   output_dir = "spectreasy_outputs/autounmix_controls"
#' )
#' ctrl$unmixing_matrix_file
#' }
autounmix_controls <- function(
    scc_dir = "scc",
    control_df = NULL,
    control_file = "fcs_control_file.csv",
    auto_create_control = TRUE,
    cytometer = "Aurora",
    auto_default_control_type = "beads",
    auto_unknown_fluor_policy = c("by_channel", "empty", "filename"),
    output_dir = "spectreasy_outputs/autounmix_controls",
    unmix_method = "WLS",
    build_qc_plots = FALSE,
    unmix_scatter_panel_size_mm = 30,
    ...
) {
    auto_unknown_fluor_policy <- match.arg(auto_unknown_fluor_policy)
    user_supplied_control_df <- !is.null(control_df)
    created_control_file <- FALSE

    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    if (!dir.exists(scc_dir)) stop("scc_dir not found: ", scc_dir)
    required_cols <- c("filename", "fluorophore", "channel")

    if (is.character(control_df) && length(control_df) == 1 && !is.na(control_df)) {
        control_file <- control_df
        control_df <- NULL
    }
    if (!is.null(control_df) && !is.data.frame(control_df)) {
        stop("control_df must be either a data.frame or a single CSV path.")
    }

    if (is.null(control_df)) {
        if (!file.exists(control_file)) {
            if (!isTRUE(auto_create_control)) {
                stop(
                    "Control file not found: ", control_file, "\n",
                    "Set auto_create_control = TRUE to auto-generate it from SCC files."
                )
            }
            message("Control file not found: ", control_file)
            message("Auto-generating control file from SCC filenames and peak channels...")
            create_autospectral_control_file(
                input_folder = scc_dir,
                include_af_folder = FALSE,
                cytometer = cytometer,
                default_control_type = auto_default_control_type,
                unknown_fluor_policy = auto_unknown_fluor_policy,
                output_file = control_file
            )
            created_control_file <- TRUE
            message("Auto-generated control file: ", control_file, " (please review marker/fluorophore/channel columns).")
        }
        control_df <- tryCatch(
            utils::read.csv(control_file, stringsAsFactors = FALSE, check.names = FALSE),
            error = function(e) NULL
        )
        if (is.null(control_df)) {
            stop("Could not read control file: ", control_file)
        }
    }

    missing_cols <- setdiff(required_cols, colnames(control_df))
    if (length(missing_cols) > 0) {
        stop(
            "Control mapping is missing required columns: ", paste(missing_cols, collapse = ", "), "\n",
            "Required columns: ", paste(required_cols, collapse = ", "), "\n",
            "Provide at least filename, fluorophore, and channel."
        )
    }
    if (!("universal.negative" %in% colnames(control_df))) {
        control_df$universal.negative <- ""
    }

    confirm_created_control_file <- function(path) {
        msg <- paste(
            c(
                "A new control file was created:",
                paste0(" - ", path),
                "Please review it before unmixing.",
                "Check at least these columns:",
                " - fluorophore / marker / channel mappings",
                " - control.type: fill non-AF rows manually with 'beads' or 'cells'",
                " - universal.negative: leave empty unless you explicitly use it"
            ),
            collapse = "\n"
        )

        if (!interactive()) {
            stop(
                msg, "\n",
                "R session is non-interactive, so confirmation prompt is not possible.\n",
                "After reviewing the file, rerun autounmix_controls().",
                call. = FALSE
            )
        }

        message(msg)
        repeat {
            ans <- tolower(trimws(readline("Proceed with autounmix_controls now? [y/n]: ")))
            if (ans %in% c("y", "yes")) return(invisible(TRUE))
            if (ans %in% c("n", "no", "")) {
                stop("Stopped after control-file creation. Review and rerun autounmix_controls() when ready.", call. = FALSE)
            }
            message("Please answer 'y' or 'n'.")
        }
    }

    if (created_control_file) {
        confirm_created_control_file(control_file)
    }

    run_preflight <- function(df) {
        validate_control_file_mapping(
            control_df = df,
            scc_dir = scc_dir,
            include_multi_af = FALSE,
            af_dir = "af",
            require_all_scc_mapped = TRUE,
            require_channels = TRUE,
            stop_on_error = FALSE
        )
    }

    preflight <- run_preflight(control_df)
    if (!preflight$ok && isTRUE(auto_create_control) && !user_supplied_control_df) {
        message("Control preflight failed; attempting automatic control-file regeneration...")
        create_autospectral_control_file(
            input_folder = scc_dir,
            include_af_folder = FALSE,
            cytometer = cytometer,
            default_control_type = auto_default_control_type,
            unknown_fluor_policy = auto_unknown_fluor_policy,
            output_file = control_file
        )
        control_df <- utils::read.csv(control_file, stringsAsFactors = FALSE, check.names = FALSE)
        preflight <- run_preflight(control_df)
    }

    if (!preflight$ok) {
        hint <- NULL
        if (any(grepl("^Empty fluorophore", preflight$errors)) && identical(auto_unknown_fluor_policy, "empty")) {
            hint <- "Tip: rerun with auto_unknown_fluor_policy = \"by_channel\" to auto-fill common fluorophore names from detected channels."
        }
        stop(
            paste(
                c(
                    "autounmix_controls preflight failed:",
                    paste0(" - ", preflight$errors),
                    if (length(preflight$warnings) > 0) c("Warnings:", paste0(" - ", preflight$warnings)) else NULL,
                    hint,
                    "Fix the control file and rerun autounmix_controls()."
                ),
                collapse = "\n"
            ),
            call. = FALSE
        )
    }

    build_plots_dir <- file.path(output_dir, "build_reference_plots")
    unmixed_dir <- file.path(output_dir, "scc_unmixed")
    spectra_file <- file.path(output_dir, "scc_spectra.png")
    reference_matrix_csv <- file.path(output_dir, "scc_reference_matrix.csv")
    unmixing_matrix_csv <- file.path(output_dir, "scc_unmixing_matrix.csv")
    unmixing_matrix_png <- file.path(output_dir, "scc_unmixing_matrix.png")
    unmixing_scatter_png <- file.path(output_dir, "scc_unmixing_scatter_matrix.png")

    M <- build_reference_matrix(
        input_folder = scc_dir,
        output_folder = build_plots_dir,
        save_qc_plots = build_qc_plots,
        control_df = control_df,
        cytometer = cytometer,
        ...
    )
    if (is.null(M) || nrow(M) == 0) stop("No valid spectra found while building reference matrix.")

    M_df <- as.data.frame(M, check.names = FALSE)
    M_df$Marker <- rownames(M)
    M_df <- M_df[, c("Marker", setdiff(colnames(M_df), "Marker")), drop = FALSE]
    utils::write.csv(M_df, reference_matrix_csv, row.names = FALSE, quote = TRUE)

    fcs_files <- list.files(scc_dir, pattern = "\\.fcs$", full.names = TRUE)
    ff_meta <- flowCore::read.FCS(fcs_files[1], transformation = FALSE, truncate_max_range = FALSE)
    pd <- flowCore::pData(flowCore::parameters(ff_meta))

    p_spectra <- plot_spectra(M, pd = pd, output_file = spectra_file)

    unmixed_list <- unmix_samples(
        sample_dir = scc_dir,
        M = M,
        method = unmix_method,
        cytometer = cytometer,
        output_dir = unmixed_dir
    )

    W <- derive_unmixing_matrix(M, method = "OLS")
    save_unmixing_matrix(W, unmixing_matrix_csv)
    p_unmix <- plot_unmixing_matrix(W, pd = pd)
    ggplot2::ggsave(unmixing_matrix_png, p_unmix, width = 200, height = 150, units = "mm")

    sample_to_marker <- NULL
    if (is.data.frame(control_df) && all(c("filename", "fluorophore") %in% colnames(control_df))) {
        sample_keys <- tools::file_path_sans_ext(basename(as.character(control_df$filename)))
        sample_vals <- trimws(as.character(control_df$fluorophore))
        keep <- !is.na(sample_keys) & sample_keys != "" & !is.na(sample_vals) & sample_vals != ""
        if (any(keep)) {
            sample_to_marker <- stats::setNames(sample_vals[keep], sample_keys[keep])
            sample_to_marker <- sample_to_marker[!duplicated(names(sample_to_marker))]
        }
    }

    plot_unmixing_scatter_matrix(
        unmixed_list = unmixed_list,
        sample_to_marker = sample_to_marker,
        markers = rownames(M),
        output_file = unmixing_scatter_png,
        transform = "none",
        panel_size_mm = unmix_scatter_panel_size_mm
    )

    invisible(list(
        M = M,
        W = W,
        unmixed_list = unmixed_list,
        reference_matrix_file = reference_matrix_csv,
        unmixing_matrix_file = unmixing_matrix_csv,
        spectra_file = spectra_file,
        unmixing_matrix_plot = unmixing_matrix_png,
        unmixing_scatter_file = unmixing_scatter_png,
        spectra_plot = p_spectra,
        unmixing_plot = p_unmix
    ))
}
