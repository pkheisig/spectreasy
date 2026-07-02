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
                                         auto_default_control_type,
                                         auto_unknown_fluor_policy) {
    message("Control file not found: ", control_file)
    message("Auto-generating control file from SCC filenames and peak channels...")
    create_control_file(
        input_folder = scc_dir,
        include_af_folder = FALSE,
        cytometer = cytometer,
        default_control_type = auto_default_control_type,
        unknown_fluor_policy = auto_unknown_fluor_policy,
        output_file = control_file
    )
    message("Auto-generated control file: ", control_file, " (please review marker/fluorophore/channel columns).")
    invisible(control_file)
}

.prepare_unmix_control_df <- function(control_file,
                                      scc_dir,
                                      auto_create_control,
                                      cytometer,
                                      auto_default_control_type,
                                      auto_unknown_fluor_policy) {
    created_control_file <- FALSE

    if (!file.exists(control_file)) {
        if (!isTRUE(auto_create_control)) {
            stop(
                "Control file not found: ", control_file, "\n",
                "Set auto_create_control = TRUE to auto-generate it from SCC files."
            )
        }
        .unmix_generate_control_file(
            scc_dir = scc_dir,
            control_file = control_file,
            cytometer = cytometer,
            auto_default_control_type = auto_default_control_type,
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

.filter_unmix_af_rows <- function(control_df, exclude_af = FALSE, emit_message = TRUE) {
    excluded_af_filenames <- character()
    if (!isTRUE(exclude_af)) {
        return(list(control_df = control_df, excluded_af_filenames = excluded_af_filenames))
    }

    af_rows <- .is_af_control_row(
        fluorophore = control_df$fluorophore,
        marker = if ("marker" %in% colnames(control_df)) control_df$marker else NULL,
        filename = control_df$filename
    )
    if (any(af_rows)) {
        excluded_af_filenames <- unique(as.character(control_df$filename[af_rows]))
        if (isTRUE(emit_message)) {
            message("exclude_af = TRUE: ignoring ", sum(af_rows), " AF/unstained control row(s) from control mapping.")
        }
        control_df <- control_df[!af_rows, , drop = FALSE]
    }

    list(control_df = control_df, excluded_af_filenames = excluded_af_filenames)
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

.run_unmix_preflight <- function(control_df, scc_dir, exclude_af = FALSE, include_multi_af = FALSE) {
    validate_control_file_mapping(
        control_df = control_df,
        scc_dir = scc_dir,
        include_multi_af = include_multi_af,
        exclude_af = exclude_af,
        af_dir = "af",
        require_all_scc_mapped = TRUE,
        require_channels = TRUE,
        stop_on_error = FALSE
    )
}

.unmix_regenerate_control_df <- function(control_file,
                                         scc_dir,
                                         cytometer,
                                         auto_default_control_type,
                                         auto_unknown_fluor_policy,
                                         exclude_af = FALSE) {
    message("Control preflight failed; attempting automatic control-file regeneration...")
    create_control_file(
        input_folder = scc_dir,
        include_af_folder = FALSE,
        cytometer = cytometer,
        default_control_type = auto_default_control_type,
        unknown_fluor_policy = auto_unknown_fluor_policy,
        output_file = control_file
    )
    control_df <- utils::read.csv(control_file, stringsAsFactors = FALSE, check.names = FALSE)
    filtered <- .filter_unmix_af_rows(control_df, exclude_af = exclude_af, emit_message = FALSE)
    list(control_df = filtered$control_df, excluded_af_filenames = filtered$excluded_af_filenames)
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

.unmix_output_paths <- function(output_dir) {
    list(
        unmixed_dir = file.path(output_dir, "unmixed_fcs"),
        qc_controls_dir = file.path(output_dir, "qc_controls"),
        spectra_file = file.path(output_dir, "scc_spectra.png"),
        af_spectra_file = file.path(output_dir, "scc_af_spectra.png"),
        reference_matrix_csv = file.path(output_dir, "scc_reference_matrix.csv"),
        detector_noise_csv = file.path(output_dir, "scc_detector_noise.csv"),
        unmixing_matrix_csv = file.path(output_dir, "scc_unmixing_matrix.csv"),
        unmixing_scatter_png = file.path(output_dir, "scc_unmixing_scatter_matrix.png"),
        variances_csv = file.path(output_dir, "scc_variances.csv"),
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

.filter_unmix_excluded_af_outputs <- function(unmixed_list, scc_dir, excluded_af_filenames, unmixed_dir, exclude_af = FALSE) {
    if (!isTRUE(exclude_af)) {
        return(unmixed_list)
    }

    scc_basenames <- list.files(scc_dir, pattern = "\\.fcs$", full.names = FALSE, ignore.case = TRUE)
    excluded_af_filenames <- unique(c(excluded_af_filenames, scc_basenames[.is_af_filename(scc_basenames)]))
    excluded_af_keys <- tools::file_path_sans_ext(basename(excluded_af_filenames))
    excluded_af_keys <- excluded_af_keys[!is.na(excluded_af_keys) & excluded_af_keys != ""]
    if (length(excluded_af_keys) == 0) {
        return(unmixed_list)
    }

    unmixed_list <- unmixed_list[!(names(unmixed_list) %in% excluded_af_keys)]
    for (sample_key in excluded_af_keys) {
        out_file <- file.path(unmixed_dir, paste0(sample_key, "_unmixed.fcs"))
        if (file.exists(out_file)) file.remove(out_file)
    }

    unmixed_list
}

.derive_unmix_static_matrix <- function(M, fcs_files, unmix_method) {
    # If M has multiple AF bands, only use the base "AF" band (and non-AF markers) for the static unmixing matrix
    # to avoid singularity (since we cannot statically unmix more markers than detectors).
    marker_names <- rownames(M)
    extra_af_rows <- grepl("^AF_", marker_names, ignore.case = TRUE)
    if (any(extra_af_rows)) {
        M_static <- M[!extra_af_rows, , drop = FALSE]
        vars <- attr(M, "variances")
        if (!is.null(vars)) {
            attr(M_static, "variances") <- vars[!extra_af_rows, , drop = FALSE]
        }
        detector_noise <- attr(M, "detector_noise")
        if (!is.null(detector_noise)) {
            attr(M_static, "detector_noise") <- detector_noise
        }
    } else {
        M_static <- M
    }

    static_unmixing_matrix_method <- toupper(unmix_method)
    if (static_unmixing_matrix_method %in% c("WLS", "RWLS")) {
        W <- derive_unmixing_matrix(M_static, method = "WLS")
    } else {
        W <- derive_unmixing_matrix(M_static, method = static_unmixing_matrix_method)
    }

    list(W = W, static_unmixing_matrix_method = static_unmixing_matrix_method)
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
#' @param auto_create_control Logical; auto-generate control file when missing.
#' @param cytometer Cytometer name used as a channel-mapping hint. The default,
#'   `"auto"`, infers the cytometer from FCS detector names when possible.
#' @param auto_default_control_type Deprecated and ignored.
#' @param auto_unknown_fluor_policy Auto-fill policy for unresolved fluorophores
#'   when creating controls (`"by_channel"`, `"empty"`, `"filename"`).
#' @param output_dir Output directory for SCC workflow artifacts.
#' @param exclude_af Logical; if `TRUE`, ignore AF/unstained controls even when
#'   they are present in the SCC folder or control mapping.
#' @param unmix_method SCC unmixing method (`"WLS"`, `"RWLS"`, `"OLS"`, `"NNLS"`).
#' @param unmix_scatter_panel_size_mm Panel size for SCC unmixing scatter matrix plot.
#' @param seed Optional integer seed for deterministic subsampling and plotting.
#' @param af_n_bands Number of AF bands to extract from the unstained control
#'   when only one AF source is available. Use `"auto"` to select the count
#'   from AF event shapes and prune near-duplicate AF signatures. Default is
#'   `"auto"`.
#' @param af_bands_per_file Number of AF bands requested per AF file when
#'   multiple AF sources are pooled. Default is 5.
#' @param af_auto_max_bands Maximum AF bands that `"auto"` may test/select.
#'   Default is 20.
#' @param af_min_cluster_events Minimum number of AF events required to keep a
#'   k-means AF cluster. Used together with `af_min_cluster_proportion`.
#' @param af_min_cluster_proportion Minimum fraction of modeled scatter-gated AF
#'   events required to keep a k-means AF cluster. Default is 0.005.
#' @param af_n_bands_sensitivity Normalized sensitivity for adding AF bands
#'   when `af_n_bands = "auto"`. Lower values allow more bands; higher values
#'   select fewer bands before near-duplicate AF signatures are pruned. Default
#'   is `1.5`.
#' @param include_multi_af Logical; whether to include additional AF files from `af_dir`. Default is FALSE.
#' @param rwls_max_iter Positive integer; number of robust reweighting
#'   iterations used when `unmix_method = "RWLS"`. The default, 1, preserves the
#'   historical behavior.
#' @param unmix_threads Positive integer; number of threads to use for event-wise
#'   multi-AF WLS/RWLS SCC unmixing. The default, 1, keeps execution single-threaded.
#' @param save_qc_plots Logical; whether to write per-control FSC/SSC,
#'   intensity-gate, and spectrum PNGs under `output_dir`.
#' @param save_report Logical; if `TRUE`, write the SCC QC PDF report
#'   automatically after controls are unmixed. Defaults to `TRUE`.
#' @param output_file Optional output path for the SCC QC PDF report. When this
#'   is `NULL`, each automatic report is written to a fresh `qc_controls`,
#'   `qc_controls_2`, ... folder under `output_dir`.
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
#'     auto_create_control = TRUE,
#'     cytometer = "auto",
#'     output_dir = "spectreasy_outputs/unmix_controls"
#'   )
#'   ctrl$unmixing_matrix_file
#' }
unmix_controls <- function(
    scc_dir = "scc",
    control_file = "fcs_mapping.csv",
    auto_create_control = TRUE,
    cytometer = "auto",
    auto_default_control_type = "beads",
    auto_unknown_fluor_policy = c("by_channel", "empty", "filename"),
    output_dir = "spectreasy_outputs/unmix_controls",
    exclude_af = FALSE,
    unmix_method = "WLS",
    unmix_scatter_panel_size_mm = 30,
    seed = NULL,
    af_n_bands = "auto",
    af_bands_per_file = 5,
    af_auto_max_bands = 20,
    af_min_cluster_events = 20,
    af_min_cluster_proportion = 0.005,
    af_n_bands_sensitivity = 1.5,
    include_multi_af = FALSE,
    rwls_max_iter = 1L,
    unmix_threads = 1L,
    save_qc_plots = FALSE,
    save_report = TRUE,
    output_file = NULL,
    use_scatter_gating = TRUE,
    ...
) {
    auto_unknown_fluor_policy <- match.arg(auto_unknown_fluor_policy)
    .with_optional_seed(seed)

    control_file <- .resolve_control_file_path(control_file)

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
        auto_create_control = auto_create_control,
        cytometer = cytometer,
        auto_default_control_type = auto_default_control_type,
        auto_unknown_fluor_policy = auto_unknown_fluor_policy
    )
    control_df <- .normalize_unmix_control_df(control_info$control_df)
    control_file <- control_info$control_file

    control_df <- .drop_missing_unmix_af_rows(control_df, scc_dir = scc_dir)
    filtered <- .filter_unmix_af_rows(control_df, exclude_af = exclude_af)
    control_df <- filtered$control_df
    excluded_af_filenames <- filtered$excluded_af_filenames

    if (isTRUE(control_info$created_control_file)) {
        .unmix_confirm_created_control_file(control_file)
    }

    preflight <- .run_unmix_preflight(control_df, scc_dir = scc_dir, exclude_af = exclude_af, include_multi_af = include_multi_af)
    if (!preflight$ok) {
        .stop_unmix_preflight(preflight, auto_unknown_fluor_policy = auto_unknown_fluor_policy)
    }

    output_paths <- .unmix_output_paths(output_dir)
    qc_controls_dir <- NULL
    if (isTRUE(save_report)) {
        if (is.null(output_file)) {
            qc_controls_dir <- .next_safe_output_dir(output_paths$qc_controls_dir)
            output_file <- file.path(qc_controls_dir, "qc_controls_report.pdf")
        } else {
            qc_controls_dir <- dirname(output_file)
        }
        message("Automatic SCC QC report enabled: ", output_file)
    }
    M <- build_reference_matrix(
        input_folder = scc_dir,
        output_folder = output_dir,
        save_qc_plots = save_qc_plots,
        control_df = control_df,
        cytometer = cytometer,
        exclude_af = exclude_af,
        af_n_bands = af_n_bands,
        af_bands_per_file = af_bands_per_file,
        af_auto_max_bands = af_auto_max_bands,
        af_min_cluster_events = af_min_cluster_events,
        af_min_cluster_proportion = af_min_cluster_proportion,
        af_n_bands_sensitivity = af_n_bands_sensitivity,
        include_multi_af = include_multi_af,
        use_scatter_gating = use_scatter_gating,
        seed = seed,
        ...
    )
    if (is.null(M) || nrow(M) == 0) stop("No valid spectra found while building reference matrix.")

    .save_reference_matrix_csv(M, output_paths$reference_matrix_csv)
    if (!is.null(attr(M, "variances"))) {
        .save_reference_matrix_csv(attr(M, "variances"), output_paths$variances_csv)
    }
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
        method = unmix_method,
        rwls_max_iter = rwls_max_iter,
        n_threads = unmix_threads,
        cytometer = cytometer,
        output_dir = output_paths$unmixed_dir,
        write_fcs = TRUE,
        save_report = FALSE
    )
    unmixed_list <- .filter_unmix_excluded_af_outputs(
        unmixed_list = unmixed_list,
        scc_dir = scc_dir,
        excluded_af_filenames = excluded_af_filenames,
        unmixed_dir = output_paths$unmixed_dir,
        exclude_af = exclude_af
    )

    static_info <- .derive_unmix_static_matrix(M, fcs_files = meta_info$fcs_files, unmix_method = unmix_method)
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
            method = unmix_method,
            qc_plot_dir = qc_controls_dir,
            save_qc_pngs = save_qc_plots,
            use_scatter_gating = use_scatter_gating,
            include_multi_af = include_multi_af,
            exclude_af = exclude_af,
            af_bands_per_file = af_bands_per_file,
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
        variances_file = output_paths$variances_csv,
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
