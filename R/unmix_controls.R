# Internal helpers for unmix_controls().
.unmix_confirm_created_control_file <- function(path) {
    msg <- paste(
        c(
            "A new control file was created:",
            paste0(" - ", path),
            "Please review it before unmixing.",
            "Check at least these columns:",
            " - fluorophore / marker / channel mappings",
            " - control.type: auto-detected from filename tokens ('beads'/'cells'); fix unknown rows as needed"
        ),
        collapse = "\n"
    )

    if (!interactive()) {
        stop(
            msg, "\n",
            "R session is non-interactive, so confirmation prompt is not possible.\n",
            "After reviewing the file, rerun unmix_controls().",
            call. = FALSE
        )
    }

    message(msg)
    repeat {
        ans <- tolower(trimws(readline("Proceed with unmix_controls now? [y/n]: ")))
        if (ans %in% c("y", "yes")) return(invisible(TRUE))
        if (ans %in% c("n", "no", "")) {
            stop("Stopped after control-file creation. Review and rerun unmix_controls() when ready.", call. = FALSE)
        }
        message("Please answer 'y' or 'n'.")
    }
}

.unmix_read_control_df <- function(control_file) {
    tryCatch(
        utils::read.csv(control_file, stringsAsFactors = FALSE, check.names = FALSE),
        error = function(e) NULL
    )
}

.unmix_generate_control_file <- function(scc_dir,
                                         control_file,
                                         cytometer,
                                         auto_unknown_fluor_policy) {
    message("Control file not found: ", control_file)
    message("Auto-generating control file from SCC filenames and peak channels...")
    create_control_file(
        input_folder = scc_dir,
        cytometer = cytometer,
        unknown_fluor_policy = auto_unknown_fluor_policy,
        output_file = control_file
    )
    message("Auto-generated control file: ", control_file, " (please review marker/fluorophore/channel columns).")
    invisible(control_file)
}

.prepare_unmix_control_df <- function(control_file,
                                      scc_dir,
                                      auto_create_mapping,
                                      cytometer,
                                      auto_unknown_fluor_policy) {
    created_control_file <- FALSE

    if (!file.exists(control_file)) {
        if (!isTRUE(auto_create_mapping)) {
            stop(
                "Control file not found: ", control_file, "\n",
                "Set auto_create_mapping = TRUE to auto-generate it from SCC files."
            )
        }
        .unmix_generate_control_file(
            scc_dir = scc_dir,
            control_file = control_file,
            cytometer = cytometer,
            auto_unknown_fluor_policy = auto_unknown_fluor_policy
        )
        created_control_file <- TRUE
    }
    control_df <- .unmix_read_control_df(control_file)
    if (is.null(control_df)) {
        stop("Could not read control file: ", control_file)
    }

    list(control_df = control_df, control_file = control_file, created_control_file = created_control_file)
}

.normalize_unmix_control_df <- function(control_df, required_cols = c("filename", "fluorophore", "channel")) {
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
    control_df
}

.drop_missing_unmix_af_rows <- function(control_df, scc_dir, emit_message = TRUE) {
    if (is.null(control_df) || nrow(control_df) == 0 || !("filename" %in% colnames(control_df))) {
        return(control_df)
    }

    af_rows <- .is_af_control_row(
        fluorophore = if ("fluorophore" %in% colnames(control_df)) control_df$fluorophore else NULL,
        marker = if ("marker" %in% colnames(control_df)) control_df$marker else NULL,
        filename = control_df$filename
    )
    if (!any(af_rows)) {
        return(control_df)
    }

    af_files <- file.path(scc_dir, as.character(control_df$filename[af_rows]))
    missing_af <- af_rows
    missing_af[af_rows] <- !file.exists(af_files)
    if (!any(missing_af)) {
        return(control_df)
    }

    if (isTRUE(emit_message)) {
        message(
            "AF/unstained control file(s) are missing; building a marker-only reference matrix. ",
            "Use unmix_samples(estimate_af = TRUE) to estimate AF from stained samples."
        )
    }
    control_df[!missing_af, , drop = FALSE]
}

.run_unmix_preflight <- function(control_df, scc_dir) {
    validate_control_file_mapping(
        control_df = control_df,
        scc_dir = scc_dir,
        require_all_scc_mapped = TRUE,
        require_channels = TRUE,
        stop_on_error = FALSE
    )
}

.stop_unmix_preflight <- function(preflight, auto_unknown_fluor_policy) {
    hint <- NULL
    if (any(grepl("^Empty fluorophore", preflight$errors)) && identical(auto_unknown_fluor_policy, "empty")) {
        hint <- "Tip: rerun with auto_unknown_fluor_policy = \"by_channel\" to auto-fill common fluorophore names from detected channels."
    }
    stop(
        paste(
            c(
                "unmix_controls preflight failed:",
                paste0(" - ", preflight$errors),
                if (length(preflight$warnings) > 0) c("Preflight notes:", paste0(" - ", preflight$warnings)) else NULL,
                hint,
                "Fix the control file and rerun unmix_controls()."
            ),
            collapse = "\n"
        )
    )
}

.normalize_unmix_controls_output_dir <- function(output_dir) {
    output_dir <- as.character(output_dir)[1]
    if (identical(basename(normalizePath(output_dir, mustWork = FALSE)), "unmixed_fcs")) {
        return(dirname(output_dir))
    }
    output_dir
}

.unmix_output_paths <- function(output_dir) {
    output_dir <- .normalize_unmix_controls_output_dir(output_dir)
    list(
        unmixed_dir = file.path(output_dir, "unmixed_fcs"),
        qc_controls_dir = file.path(output_dir, "qc_controls"),
        spectra_file = file.path(output_dir, "scc_spectra.png"),
        af_spectra_file = file.path(output_dir, "scc_af_spectra.png"),
        reference_matrix_csv = file.path(output_dir, "scc_reference_matrix.csv"),
        detector_noise_csv = file.path(output_dir, "scc_detector_noise.csv"),
        unmixing_matrix_csv = file.path(output_dir, "scc_unmixing_matrix.csv"),
        unmixing_scatter_png = file.path(output_dir, "scc_unmixing_scatter_matrix.png"),
        qc_report_pdf = file.path(output_dir, "qc_controls", "qc_controls_report.pdf")
    )
}

.save_reference_matrix_csv <- function(M, path) {
    M_df <- as.data.frame(M, check.names = FALSE)
    M_df$Marker <- rownames(M)
    M_df <- M_df[, c("Marker", setdiff(colnames(M_df), "Marker")), drop = FALSE]
    utils::write.csv(M_df, path, row.names = FALSE, quote = TRUE)
    invisible(path)
}

.read_unmix_metadata_pd <- function(scc_dir) {
    fcs_files <- list.files(scc_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
    if (length(fcs_files) == 0) {
        stop("No FCS files found in scc_dir: ", scc_dir)
    }
    ff_meta <- flowCore::read.FCS(fcs_files[1], transformation = FALSE, truncate_max_range = FALSE)
    list(fcs_files = fcs_files, pd = flowCore::pData(flowCore::parameters(ff_meta)))
}

.derive_unmix_static_matrix <- function(M, fcs_files, unmixing_method) {
    # Static exported matrices cannot represent per-event AF selection reliably.
    # Drop all AF rows and keep AF modeling in dynamic unmixing results.
    marker_names <- rownames(M)
    af_rows <- grepl("^AF($|_)", marker_names, ignore.case = TRUE)
    if (any(af_rows)) {
        M_static <- M[!af_rows, , drop = FALSE]
        detector_noise <- attr(M, "detector_noise")
        if (!is.null(detector_noise)) {
            attr(M_static, "detector_noise") <- detector_noise
        }
    } else {
        M_static <- M
    }

    static_unmixing_matrix_method <- toupper(unmixing_method)
    if (static_unmixing_matrix_method %in% c("WLS", "RWLS")) {
        W <- tryCatch(
            derive_unmixing_matrix(M_static, method = "WLS"),
            error = function(e) {
                if (!grepl("singular", conditionMessage(e), ignore.case = TRUE)) {
                    stop(e)
                }
                .derive_regularized_static_wls_matrix(M_static)
            }
        )
    } else {
        W <- derive_unmixing_matrix(M_static, method = static_unmixing_matrix_method)
    }

    list(W = W, static_unmixing_matrix_method = static_unmixing_matrix_method)
}

.derive_regularized_static_wls_matrix <- function(M) {
    detector_weights <- .wls_static_detector_weights(
        M = M,
        background_noise = .default_wls_background_noise(),
        signal_scale = .default_wls_signal_scale(),
        max_weight_ratio = .default_wls_max_weight_ratio()
    )
    V_inv <- diag(detector_weights)
    Mt <- t(M)
    MVMt <- M %*% V_inv %*% Mt
    ridge <- max(c(diag(MVMt), 1), na.rm = TRUE) * 1e-8
    W <- solve(MVMt + diag(ridge, nrow(MVMt)), M %*% V_inv)
    rownames(W) <- rownames(M)
    colnames(W) <- colnames(M)
    W
}

.resolve_unmix_marker_mappings <- function(control_df) {
    sample_to_marker <- NULL
    marker_display <- NULL

    if (is.data.frame(control_df) && all(c("filename", "fluorophore") %in% colnames(control_df))) {
        sample_keys <- tools::file_path_sans_ext(basename(as.character(control_df$filename)))
        primary_vals <- trimws(as.character(control_df$fluorophore))
        secondary_vals <- if ("marker" %in% colnames(control_df)) {
            trimws(as.character(control_df$marker))
        } else {
            rep("", length(primary_vals))
        }
        secondary_vals[is.na(secondary_vals)] <- ""
        keep <- !is.na(sample_keys) & sample_keys != "" & !is.na(primary_vals) & primary_vals != ""
        if (any(keep)) {
            sample_to_marker <- stats::setNames(primary_vals[keep], sample_keys[keep])
            sample_to_marker <- sample_to_marker[!duplicated(names(sample_to_marker))]

            display_vals <- primary_vals
            show_secondary <- secondary_vals != "" &
                toupper(secondary_vals) != "AUTOFLUORESCENCE" &
                tolower(secondary_vals) != tolower(primary_vals)
            display_vals[show_secondary] <- paste0(primary_vals[show_secondary], " / ", secondary_vals[show_secondary])
            marker_display <- stats::setNames(display_vals[keep], primary_vals[keep])
            marker_display <- marker_display[names(marker_display) != ""]
            marker_display <- marker_display[!duplicated(names(marker_display))]
        }
    }

    list(sample_to_marker = sample_to_marker, marker_display = marker_display)
}

#' Unmix Single-Color Controls
#'
#' SCC workflow helper that:
#' 1) validates/creates the control file,
#' 2) builds the reference matrix from SCCs,
#' 3) unmixed SCC files,
#' 4) saves reference/unmixing matrices and SCC-derived WLS detector noise floors,
#' 5) saves standalone spectra, unmixing-matrix, and SCC scatter PNGs, and
#'    optionally writes per-control QC PNGs.
#'
#' This function is intended as the control-stage step before GUI adjustment
#' and downstream sample unmixing.
#'
#' @param scc_dir Directory containing SCC FCS files.
#' @param control_file Path to control mapping CSV.
#' @param auto_create_mapping Logical; auto-generate the FCS mapping CSV when
#'   `control_file` is missing.
#' @param cytometer Cytometer name used as a channel-mapping hint. The default,
#'   `"auto"`, infers the cytometer from FCS detector names when possible.
#' @param auto_unknown_fluor_policy Auto-fill policy for unresolved fluorophores
#'   when creating controls (`"by_channel"`, `"empty"`, `"filename"`).
#' @param output_dir Output directory for SCC workflow artifacts.
#' @param unmixing_method SCC unmixing method (`"WLS"`, `"RWLS"`, `"OLS"`, `"NNLS"`).
#' @param unmix_scatter_panel_size_mm Panel size for SCC unmixing scatter matrix plot.
#' @param seed Optional integer seed for deterministic subsampling and plotting.
#' @param af_n_bands Number of AF basis signatures to extract from pooled
#'   unstained/AF control events. The default, `10`, is a conservative fixed
#'   AF bank size.
#' @param af_min_cluster_events Minimum number of AF events required to keep a
#'   k-means AF cluster. Used together with `af_min_cluster_proportion`.
#' @param af_min_cluster_proportion Minimum fraction of modeled scatter-gated AF
#'   events required to keep a k-means AF cluster. Default is 0.005.
#' @param rwls_max_iter Positive integer; number of robust reweighting
#'   iterations used when `unmixing_method = "RWLS"`. The default, 1, preserves the
#'   historical behavior.
#' @param unmix_threads Positive integer; number of threads to use for event-wise
#'   multi-AF WLS/RWLS SCC unmixing. The default, 1, keeps execution single-threaded.
#' @param save_qc_plots Logical; whether to write per-control FSC/SSC,
#'   intensity-gate, and spectrum PNGs under `output_dir`.
#' @param save_report Logical; if `TRUE`, write the SCC QC PDF report
#'   automatically after controls are unmixed. Defaults to `TRUE`.
#' @param use_scatter_gating Logical; if `TRUE` (default), use the intensity-vs-FSC
#'   scatter gate for final positive/negative population selection. If `FALSE`,
#'   use the legacy one-dimensional histogram gate.
#' @param ... Additional arguments forwarded to [build_reference_matrix()].
#'
#' @return List with `M`, `W`, `unmixed_list`, and key output file paths,
#'   including `detector_noise_file` for the SCC-derived WLS noise floors.
#' @export
#' @examples
#' if (interactive()) {
#'   ctrl <- unmix_controls(
#'     scc_dir = "scc",
#'     control_file = "fcs_mapping.csv",
#'     auto_create_mapping = TRUE,
#'     cytometer = "auto",
#'     output_dir = "spectreasy_outputs/unmix_controls"
#'   )
#'   ctrl$unmixing_matrix_file
#' }
unmix_controls <- function(
    scc_dir = "scc",
    control_file = "fcs_mapping.csv",
    auto_create_mapping = TRUE,
    cytometer = "auto",
    auto_unknown_fluor_policy = c("by_channel", "empty", "filename"),
    output_dir = "spectreasy_outputs/unmix_controls",
    unmixing_method = "WLS",
    unmix_scatter_panel_size_mm = 30,
    seed = NULL,
    af_n_bands = 10,
    af_min_cluster_events = 20,
    af_min_cluster_proportion = 0.005,
    rwls_max_iter = 1L,
    unmix_threads = 1L,
    save_qc_plots = FALSE,
    save_report = TRUE,
    use_scatter_gating = TRUE,
    ...
) {
    auto_unknown_fluor_policy <- match.arg(auto_unknown_fluor_policy)
    .with_optional_seed(seed)
    extra_args <- list(...)

    control_file <- .resolve_control_file_path(control_file)
    output_dir <- .normalize_unmix_controls_output_dir(output_dir)

    if (!dir.exists(output_dir)) {
        dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
        if (!dir.exists(output_dir)) {
            Sys.sleep(0.5)
            dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
        }
    }
    if (!dir.exists(scc_dir)) stop("scc_dir not found: ", scc_dir)

    control_info <- .prepare_unmix_control_df(
        control_file = control_file,
        scc_dir = scc_dir,
        auto_create_mapping = auto_create_mapping,
        cytometer = cytometer,
        auto_unknown_fluor_policy = auto_unknown_fluor_policy
    )
    control_df <- .normalize_unmix_control_df(control_info$control_df)
    control_file <- control_info$control_file

    control_df <- .drop_missing_unmix_af_rows(control_df, scc_dir = scc_dir)

    if (isTRUE(control_info$created_control_file)) {
        .unmix_confirm_created_control_file(control_file)
    }

    preflight <- .run_unmix_preflight(control_df, scc_dir = scc_dir)
    if (!preflight$ok) {
        .stop_unmix_preflight(preflight, auto_unknown_fluor_policy = auto_unknown_fluor_policy)
    }

    output_paths <- .unmix_output_paths(output_dir)
    qc_controls_dir <- NULL
    output_file <- NULL
    if (isTRUE(save_report)) {
        qc_controls_dir <- .next_safe_output_dir(output_paths$qc_controls_dir)
        output_file <- file.path(qc_controls_dir, "qc_controls_report.pdf")
        message("Automatic SCC QC report enabled: ", output_file)
    }
    build_qc_plots <- isTRUE(save_qc_plots) || isTRUE(save_report)
    build_output_folder <- if (isTRUE(save_qc_plots)) {
        output_dir
    } else if (isTRUE(save_report)) {
        tempfile("scc_report_plots_")
    } else {
        output_dir
    }
    build_args <- c(list(
        input_folder = scc_dir,
        output_folder = build_output_folder,
        save_qc_plots = build_qc_plots,
        control_df = control_df,
        cytometer = cytometer,
        af_n_bands = af_n_bands,
        af_min_cluster_events = af_min_cluster_events,
        af_min_cluster_proportion = af_min_cluster_proportion,
        use_scatter_gating = use_scatter_gating,
        seed = seed
    ), extra_args)
    M <- do.call(build_reference_matrix, build_args)
    if (is.null(M) || nrow(M) == 0) stop("No valid spectra found while building reference matrix.")

    .save_reference_matrix_csv(M, output_paths$reference_matrix_csv)
    .save_detector_noise_csv(attr(M, "detector_noise"), output_paths$detector_noise_csv)
    meta_info <- .read_unmix_metadata_pd(scc_dir)

    # Split M into fluorophore rows and AF/autofluorescence rows to plot separately
    af_rows <- grepl("^AF($|_)", rownames(M), ignore.case = TRUE)
    M_fluor <- M[!af_rows, , drop = FALSE]
    M_af <- M[af_rows, , drop = FALSE]

    p_spectra <- if (nrow(M_fluor) > 0) {
        plot_spectra(M_fluor, pd = meta_info$pd, output_file = if (isTRUE(save_qc_plots)) output_paths$spectra_file else NULL)
    } else {
        NULL
    }

    p_af_spectra <- if (nrow(M_af) > 0) {
        plot_spectra(M_af, pd = meta_info$pd, output_file = if (isTRUE(save_qc_plots)) output_paths$af_spectra_file else NULL)
    } else {
        NULL
    }

    unmixed_list <- unmix_samples(
        sample_dir = scc_dir,
        M = M,
        unmixing_method = unmixing_method,
        rwls_max_iter = rwls_max_iter,
        n_threads = unmix_threads,
        output_dir = output_paths$unmixed_dir,
        write_fcs = TRUE,
        save_report = FALSE
    )

    static_info <- .derive_unmix_static_matrix(M, fcs_files = meta_info$fcs_files, unmixing_method = unmixing_method)
    W <- static_info$W
    save_unmixing_matrix(W, output_paths$unmixing_matrix_csv)

    marker_mapping <- .resolve_unmix_marker_mappings(control_df)
    scatter_markers <- rownames(M)
    extra_af_rows <- grepl("^AF_", scatter_markers, ignore.case = TRUE)
    scatter_markers <- scatter_markers[!extra_af_rows]

    p_scatter <- plot_unmixing_scatter_matrix(
        unmixed_list = unmixed_list,
        sample_to_marker = marker_mapping$sample_to_marker,
        markers = scatter_markers,
        marker_display = NULL,
        output_file = if (isTRUE(save_qc_plots)) output_paths$unmixing_scatter_png else NULL,
        transform = "none",
        panel_size_mm = unmix_scatter_panel_size_mm,
        seed = seed
    )

    qc_report <- NULL
    if (isTRUE(save_report)) {
        qc_report <- qc_controls(
            M = M,
            scc_dir = scc_dir,
            output_file = output_file,
            control_file = control_file,
            cytometer = cytometer,
            unmixing_method = unmixing_method,
            qc_plot_dir = qc_controls_dir,
            save_qc_pngs = save_qc_plots,
            qc_metrics_dir = qc_controls_dir,
            unmixed_list = unmixed_list,
            qc_summary = attr(M, "qc_summary"),
            report_plot_dir = attr(M, "qc_plot_dir"),
            pd = meta_info$pd,
            af_bank_info = attr(M, "af_bank_info"),
            cleanup_report_plot_dir = !isTRUE(save_qc_plots),
            use_scatter_gating = use_scatter_gating,
            unmix_scatter_max_points = 1000,
            seed = seed,
            unmixing_matrix_file = output_paths$reference_matrix_csv
        )
        if (!file.exists(output_file)) {
            stop("Automatic SCC QC report was requested but was not created at: ", output_file, call. = FALSE)
        }
    }

    invisible(list(
        M = M,
        W = W,
        unmixed_list = unmixed_list,
        reference_matrix_file = output_paths$reference_matrix_csv,
        detector_noise_file = output_paths$detector_noise_csv,
        unmixing_matrix_file = output_paths$unmixing_matrix_csv,
        qc_report_file = if (isTRUE(save_report)) output_file else NULL,
        qc_controls_dir = if (isTRUE(save_report)) qc_controls_dir else NULL,
        qc_metrics_dir = if (isTRUE(save_report)) qc_controls_dir else NULL,
        qc_report = qc_report,
        spectra_file = if (isTRUE(save_qc_plots)) output_paths$spectra_file else NULL,
        af_spectra_file = if (isTRUE(save_qc_plots) && nrow(M_af) > 0) output_paths$af_spectra_file else NULL,
        unmixing_scatter_file = if (isTRUE(save_qc_plots)) output_paths$unmixing_scatter_png else NULL,
        qc_plot_dir = if (isTRUE(save_qc_plots)) output_dir else NULL,
        static_unmixing_matrix_method = static_info$static_unmixing_matrix_method,
        spectra_plot = p_spectra,
        af_spectra_plot = p_af_spectra,
        unmixing_scatter_plot = p_scatter
    ))
}
