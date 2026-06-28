# Internal helpers for unmix_controls().
.unmix_confirm_created_control_file <- function(path) {
    msg <- paste(
        c(
            "A new control file was created:",
            paste0(" - ", path),
            "Please review it before unmixing.",
            "Check at least these columns:",
            " - fluorophore / marker / channel mappings",
            " - control.type: auto-detected from filename tokens ('beads'/'cells'); fix unknown rows as needed",
            " - universal.negative: leave empty unless you explicitly use it"
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
        spectra_file = file.path(output_dir, "scc_spectra.png"),
        af_spectra_file = file.path(output_dir, "scc_af_spectra.png"),
        reference_matrix_csv = file.path(output_dir, "scc_reference_matrix.csv"),
        detector_noise_csv = file.path(output_dir, "scc_detector_noise.csv"),
        unmixing_matrix_csv = file.path(output_dir, "scc_unmixing_matrix.csv"),
        spectral_variant_library_rds = file.path(output_dir, "scc_spectral_variants.rds"),
        unmixing_scatter_png = file.path(output_dir, "scc_unmixing_scatter_matrix.png"),
        variances_csv = file.path(output_dir, "scc_variances.csv"),
        qc_report_pdf = file.path(output_dir, "qc_controls_report.pdf")
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
#' @param af_n_bands Number of AF basis signatures to extract from the pooled
#'   unstained/AF control events. Use `"auto"` to build the default SOM AF bank.
#'   Default is `"auto"`.
#' @param af_bands_per_file Deprecated compatibility argument. Multiple AF
#'   sources are pooled before SOM extraction; `af_n_bands`/`af_auto_max_bands`
#'   control the size of the one shared AF bank.
#' @param af_auto_max_bands Maximum SOM nodes that `"auto"` may create before
#'   prepending the mean AF row.
#'   Default is 100.
#' @param af_min_cluster_events Compatibility argument retained for older
#'   workflows.
#' @param af_min_cluster_proportion Compatibility argument retained for older
#'   workflows.
#' @param af_n_bands_sensitivity Compatibility argument retained for older
#'   workflows.
#' @param af_refine Logical; if `TRUE`, run the optional second-pass AF
#'   refinement while building the reference matrix.
#' @param af_refine_problem_quantile Quantile used to choose high-error
#'   unstained cells for AF refinement.
#' @param include_multi_af Logical; whether to include additional AF files from `af_dir`. Default is FALSE.
#' @param rwls_max_iter Positive integer; number of robust reweighting
#'   iterations used when `unmix_method = "RWLS"`. The default, 1, preserves the
#'   historical behavior.
#' @param multithreading Logical; if `TRUE`, allow event-wise multi-AF WLS/RWLS
#'   SCC unmixing to use multiple threads. The default, `FALSE`, keeps
#'   execution single-threaded.
#' @param n_threads `"auto"` or positive integer; thread count to use when
#'   `multithreading = TRUE`. `"auto"` uses `RcppParallel::defaultNumThreads()`.
#'   Integers larger than the available thread count are clipped to the
#'   available count.
#' @param save_qc_plots Logical; whether to write per-control FSC/SSC,
#'   intensity-gate, and spectrum PNGs under `output_dir`.
#' @param save_report Logical; if `TRUE`, write the SCC QC PDF report from the
#'   matrix and unmixed controls produced by this call, without rerunning SCC
#'   unmixing.
#' @param output_file Optional output path for the SCC QC PDF report.
#' @param use_scatter_gating Logical; if `TRUE` (default), use the intensity-vs-FSC
#'   scatter gate for final positive/negative population selection. If `FALSE`,
#'   use the legacy one-dimensional histogram gate.
#' @param clean_scc_with_unstained Logical; if `TRUE`, clean cell SCC spectra
#'   with scatter-matched unstained/AF events before deriving spectra and
#'   spectral variants.
#' @param scc_background_method Background method for cell SCC cleaning.
#'   `"scatter_knn"` matches stained cells to unstained cells by FSC/SSC.
#' @param scc_background_k Number of nearest unstained cells averaged for
#'   scatter-matched background subtraction.
#' @param optimize_spectral_variants Logical; if `TRUE`, learn conservative
#'   per-fluorophore spectral variants from SCC controls and use them
#'   automatically while unmixing controls and later samples.
#' @param spectral_variant_som_nodes Number of SOM nodes used per fluorophore
#'   when learning spectral variants.
#' @param spectral_variant_top_k Number of best variant candidates to test per
#'   positive fluorophore during event-wise optimization.
#' @param spectral_variant_cosine_threshold Minimum cosine similarity to the
#'   base fluorophore spectrum required for a variant to be retained.
#' @param spectral_variant_max_variants Maximum retained variants per
#'   fluorophore.
#' @param spectral_variant_min_events Minimum cleaned positive events required
#'   to learn variants for one fluorophore.
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
    af_bands_per_file = NULL,
    af_auto_max_bands = 100,
    af_min_cluster_events = 20,
    af_min_cluster_proportion = 0.005,
    af_n_bands_sensitivity = 1.5,
    af_refine = FALSE,
    af_refine_problem_quantile = 0.99,
    include_multi_af = FALSE,
    rwls_max_iter = 1L,
    multithreading = FALSE,
    n_threads = "auto",
    save_qc_plots = FALSE,
    save_report = TRUE,
    output_file = NULL,
    use_scatter_gating = TRUE,
    clean_scc_with_unstained = TRUE,
    scc_background_method = c("scatter_knn", "none"),
    scc_background_k = 3L,
    optimize_spectral_variants = TRUE,
    spectral_variant_som_nodes = 16L,
    spectral_variant_top_k = 3L,
    spectral_variant_cosine_threshold = 0.98,
    spectral_variant_max_variants = 8L,
    spectral_variant_min_events = 50L,
    ...
) {
    auto_unknown_fluor_policy <- match.arg(auto_unknown_fluor_policy)
    scc_background_args <- .validate_scc_background_args(
        clean_scc_with_unstained = clean_scc_with_unstained,
        scc_background_method = scc_background_method,
        scc_background_k = scc_background_k
    )
    clean_scc_with_unstained <- scc_background_args$enabled
    scc_background_method <- scc_background_args$method
    scc_background_k <- scc_background_args$k
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
    report_plot_output_dir <- output_dir
    cleanup_report_plot_dir <- NULL
    build_report_plots <- isTRUE(save_qc_plots) || isTRUE(save_report)
    if (isTRUE(save_report) && !isTRUE(save_qc_plots)) {
        report_plot_output_dir <- tempfile("spectreasy_controls_report_plots_")
        cleanup_report_plot_dir <- report_plot_output_dir
        on.exit(unlink(cleanup_report_plot_dir, recursive = TRUE, force = TRUE), add = TRUE)
    }

    M <- build_reference_matrix(
        input_folder = scc_dir,
        output_folder = report_plot_output_dir,
        save_qc_plots = build_report_plots,
        control_df = control_df,
        cytometer = cytometer,
        exclude_af = exclude_af,
        af_n_bands = af_n_bands,
        af_bands_per_file = af_bands_per_file,
        af_auto_max_bands = af_auto_max_bands,
        af_min_cluster_events = af_min_cluster_events,
        af_min_cluster_proportion = af_min_cluster_proportion,
        af_n_bands_sensitivity = af_n_bands_sensitivity,
        af_refine = af_refine,
        af_refine_problem_quantile = af_refine_problem_quantile,
        include_multi_af = include_multi_af,
        use_scatter_gating = use_scatter_gating,
        clean_scc_with_unstained = clean_scc_with_unstained,
        scc_background_method = scc_background_method,
        scc_background_k = scc_background_k,
        seed = seed,
        ...
    )
    if (is.null(M) || nrow(M) == 0) stop("No valid spectra found while building reference matrix.")

    spectral_variant_library <- NULL
    spectral_variant_library_file <- NULL
    if (isTRUE(optimize_spectral_variants)) {
        spectral_variant_library <- tryCatch(
            .learn_spectral_variant_library(
                scc_dir = scc_dir,
                control_df = control_df,
                M = M,
                enabled = TRUE,
                som_nodes = spectral_variant_som_nodes,
                cosine_threshold = spectral_variant_cosine_threshold,
                max_variants = spectral_variant_max_variants,
                min_events = spectral_variant_min_events,
                clean_scc_with_unstained = clean_scc_with_unstained,
                scc_background_method = scc_background_method,
                scc_background_k = scc_background_k,
                include_multi_af = include_multi_af,
                exclude_af = exclude_af,
                seed = seed,
                warn = TRUE
            ),
            error = function(e) {
                warning(
                    "Spectral-variant optimization could not build a variant library; using the base reference matrix. Reason: ",
                    conditionMessage(e),
                    call. = FALSE
                )
                NULL
            }
        )
        if (!is.null(spectral_variant_library)) {
            attr(M, "spectral_variant_library") <- spectral_variant_library
            .save_spectral_variant_library(spectral_variant_library, output_paths$spectral_variant_library_rds)
            spectral_variant_library_file <- output_paths$spectral_variant_library_rds
        }
    }

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
        multithreading = multithreading,
        n_threads = n_threads,
        cytometer = cytometer,
        output_dir = output_paths$unmixed_dir,
        write_fcs = TRUE,
        save_report = FALSE,
        optimize_spectral_variants = isTRUE(optimize_spectral_variants),
        spectral_variant_library = spectral_variant_library,
        spectral_variant_top_k = spectral_variant_top_k
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

    if (is.null(output_file)) {
        output_file <- output_paths$qc_report_pdf
    }
    qc_report <- NULL
    if (isTRUE(save_report)) {
        qc_report <- .write_scc_qc_report(
            M_built = M,
            M_report = M,
            qc_summary = attr(M, "qc_summary"),
            report_plot_dir = attr(M, "qc_plot_dir"),
            unmixed_list = unmixed_list,
            scc_dir = scc_dir,
            output_file = output_file,
            cytometer = cytometer,
            method = unmix_method,
            use_scatter_gating = use_scatter_gating,
            seed = seed,
            retained_qc_plot_dir = if (isTRUE(save_qc_plots)) output_dir else NULL
        )
    }

    M_return <- M
    if (!isTRUE(save_qc_plots)) {
        attr(M_return, "qc_plot_dir") <- NULL
    }

    invisible(list(
        M = M_return,
        W = W,
        unmixed_list = unmixed_list,
        reference_matrix_file = output_paths$reference_matrix_csv,
        detector_noise_file = output_paths$detector_noise_csv,
        spectral_variant_library = spectral_variant_library,
        spectral_variant_library_file = spectral_variant_library_file,
        spectral_variant_info = if (!is.null(spectral_variant_library)) spectral_variant_library$info else NULL,
        unmixing_matrix_file = output_paths$unmixing_matrix_csv,
        variances_file = output_paths$variances_csv,
        qc_report_file = if (isTRUE(save_report)) output_file else NULL,
        qc_report = qc_report,
        qc_summary = attr(M, "qc_summary"),
        af_bank_info = attr(M, "af_bank_info"),
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
