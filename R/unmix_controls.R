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

    .spectreasy_console_header("control file created")
    .spectreasy_console_field("File", .spectreasy_console_path(path))
    .spectreasy_console_field("Review", "marker, fluorophore, channel, and control.type columns")
    .spectreasy_console_footer(blank = FALSE)
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
    .spectreasy_console_field("Controls", paste0(.spectreasy_console_path(control_file), " not found"))
    .spectreasy_console_step("Create mapping", "from SCC filenames and peak channels")
    create_control_file(
        input_folder = scc_dir,
        cytometer = cytometer,
        unknown_fluor_policy = auto_unknown_fluor_policy,
        output_file = control_file
    )
    .spectreasy_console_field("Controls", paste0(.spectreasy_console_path(control_file), " created"))
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
        .spectreasy_console_step("AF controls", "missing; building marker-only reference matrix")
        .spectreasy_console_step("Next option", "use unmix_samples(estimate_af = TRUE) for stained-sample AF")
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
        spectral_variant_library_rds = file.path(output_dir, "scc_spectral_variants.rds"),
        unmixing_scatter_png = file.path(output_dir, "scc_unmixing_scatter_matrix.png"),
        qc_report_pdf = file.path(output_dir, "qc_controls", "qc_controls_report.pdf")
    )
}

.run_optional_unmix_artifact <- function(label, expr) {
    tryCatch(
        force(expr),
        error = function(e) {
            warning(
                label, " could not be created: ",
                conditionMessage(e),
                ". Continuing unmixing.",
                call. = FALSE
            )
            NULL
        }
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
    if (!dir.exists(scc_dir)) {
        .spectreasy_stop_missing_directory(scc_dir, label = "scc_dir")
    }
    fcs_files <- list.files(scc_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
    if (length(fcs_files) == 0) {
        .spectreasy_stop_empty_fcs_directory(scc_dir, label = "scc_dir")
    }
    ff_meta <- .spectreasy_read_fcs(fcs_files[1], label = "SCC FCS file")
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

    static_unmixing_matrix_method <- .normalize_unmix_method(unmixing_method)
    static_solver_method <- .solver_method_for_unmix(static_unmixing_matrix_method)
    if (static_solver_method %in% c("WLS", "RWLS")) {
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
        W <- derive_unmixing_matrix(M_static, method = static_solver_method)
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
#' @param unmixing_method SCC unmixing method (`"WLS"`, `"RWLS"`,
#'   `"OLS"`, `"NNLS"`, `"AutoSpectral"`, or `"Spectreasy"`).
#'   `AutoSpectral` keeps the k-means AF bank controlled by `af_n_bands`,
#'   then uses OLS after AF matching plus AutoSpectral-style SCC cleanup and
#'   spectral variants. `Spectreasy` uses the same cleanup and variant
#'   machinery, then blends the AutoSpectral-style OLS marker fit with a
#'   marker-only OLS anchor using decoder-projected AF impact weights.
#' @param unmix_scatter_panel_size_mm Panel size for SCC unmixing scatter matrix plot.
#' @param seed Optional integer seed for deterministic subsampling and plotting.
#' @param af_n_bands Number of AF basis signatures to extract from pooled
#'   unstained/AF control events. The default, `100`, builds a broad fixed
#'   AF bank for Spectreasy unmixing.
#' @param af_min_cluster_events Minimum number of AF events required to keep a
#'   k-means AF cluster. Used together with `af_min_cluster_proportion`.
#' @param af_min_cluster_proportion Minimum fraction of modeled scatter-gated AF
#'   events required to keep a k-means AF cluster. Default is 0.005.
#' @param rwls_max_iter Positive integer; number of robust reweighting
#'   iterations used when `unmixing_method = "RWLS"`. The default, 1, preserves the
#'   historical behavior.
#' @param n_threads Positive integer; number of threads to use for event-wise
#'   multi-AF WLS/RWLS SCC unmixing. The default, 1, keeps execution single-threaded.
#' @param save_qc_plots Logical; whether to write per-control FSC/SSC,
#'   intensity-gate, and spectrum PNGs under `output_dir`.
#' @param save_report Logical; if `TRUE`, write the SCC QC PDF report
#'   automatically after controls are unmixed. Defaults to `TRUE`.
#' @param use_scatter_gating Logical; if `TRUE` (default), use the intensity-vs-FSC
#'   scatter gate for final positive/negative population selection. If `FALSE`,
#'   use the legacy one-dimensional histogram gate. Ignored when
#'   `unmixing_method = "AutoSpectral"` or `"Spectreasy"`; these methods use
#'   their own post-FSC/SSC SCC spectral selector.
#' @param manual_gating Logical; if `TRUE` (default), launch the manual control
#'   gating GUI before building the SCC reference matrix. If `FALSE`, the GUI is
#'   skipped; an existing `gating_file`/`manual_gate_file` is still reused when
#'   supplied.
#' @param manual_gate_file Gate CSV written by [gate_controls()]. Relative paths
#'   are resolved from the current working directory. Kept for backwards
#'   compatibility; prefer `gating_file` in new code.
#' @param gating_file Optional gate CSV to reuse without launching the GUI, for
#'   example `file.path(getwd(), "ssc_gate_config.csv")`. If supplied explicitly,
#'   the file must exist. Defaults to `manual_gate_file`.
#' @param clean_scc_with_unstained Logical; when `unmixing_method =
#'   `"AutoSpectral"` or `"Spectreasy"`, subtract matching unstained/negative
#'   background events before calculating SCC spectra.
#' @param scc_background_method Background subtraction method for AutoSpectral
#'   SCC cleanup (`"scatter_knn"` or `"none"`).
#' @param scc_background_k Number of nearest unstained/negative events averaged
#'   for scatter-matched SCC background subtraction.
#' @param spectral_variant_som_nodes Number of nodes used per fluorophore when
#'   learning AutoSpectral spectral variants.
#' @param spectral_variant_top_k Number of best variant candidates to test per
#'   positive fluorophore during AutoSpectral unmixing.
#' @param spectral_variant_cosine_threshold Minimum cosine similarity to the
#'   base fluorophore spectrum required for a variant to be retained.
#' @param spectral_variant_max_variants Maximum retained variants per
#'   fluorophore.
#' @param spectral_variant_min_events Minimum cleaned positive SCC events
#'   required to learn variants for one fluorophore.
#' @param spectreasy_weight_quantile Numeric in `[0, 1]`; only accepted when
#'   `unmixing_method = "Spectreasy"`. Controls the quantile of
#'   decoder-projected AF impacts used as the soft-saturation scale for
#'   marker-specific AutoSpectral mixing. The default is `0.9`.
#' @param autospectral_n_candidates Number of peak-bright SCC candidate events
#'   considered by the AutoSpectral-style selector.
#' @param autospectral_n_spectral Number of least-background-like SCC events
#'   kept for spectrum calculation by the AutoSpectral-style selector.
#' @param autospectral_min_events Minimum event count required by the
#'   AutoSpectral-style SCC selector.
#' @param refine Logical; when `unmixing_method = "AutoSpectral"`, refine the
#'   fixed-size k-means AF bank with native AutoSpectral-style unstained residual
#'   modulation. `TRUE` is rejected for all other unmixing methods.
#' @param unmix_threads Deprecated compatibility alias for `n_threads`. New code
#'   should use `n_threads`.
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
    unmixing_method = "Spectreasy",
    unmix_scatter_panel_size_mm = 30,
    seed = NULL,
    af_n_bands = 100,
    af_min_cluster_events = 20,
    af_min_cluster_proportion = 0.005,
    rwls_max_iter = 1L,
    n_threads = 1L,
    save_qc_plots = FALSE,
    save_report = TRUE,
    use_scatter_gating = TRUE,
    manual_gating = TRUE,
    manual_gate_file = "ssc_gate_config.csv",
    gating_file = manual_gate_file,
    clean_scc_with_unstained = TRUE,
    scc_background_method = c("scatter_knn", "none"),
    scc_background_k = 2L,
    spectral_variant_som_nodes = 16L,
    spectral_variant_top_k = 3L,
    spectral_variant_cosine_threshold = 0.98,
    spectral_variant_max_variants = 8L,
    spectral_variant_min_events = 50L,
    spectreasy_weight_quantile = 0.9,
    autospectral_n_candidates = 1000L,
    autospectral_n_spectral = 200L,
    autospectral_min_events = 10L,
    refine = FALSE,
    unmix_threads = NULL,
    ...
) {
    spectreasy_weight_quantile_missing <- missing(spectreasy_weight_quantile)
    manual_gate_file_missing <- missing(manual_gate_file)
    gating_file_missing <- missing(gating_file)
    n_threads_missing <- missing(n_threads)
    if (!is.null(unmix_threads)) {
        if (!n_threads_missing && !identical(as.integer(n_threads[1]), as.integer(unmix_threads[1]))) {
            stop("n_threads and unmix_threads must match when both are supplied.", call. = FALSE)
        }
        n_threads <- unmix_threads
    }
    n_threads <- .normalize_unmix_threads(n_threads)
    auto_unknown_fluor_policy <- match.arg(auto_unknown_fluor_policy)
    unmixing_method <- .normalize_unmix_method(unmixing_method)
    use_autospectral <- .is_autospectral_style_method(unmixing_method)
    use_refine <- identical(unmixing_method, "AutoSpectral")
    refine <- .validate_reference_refine_arg(refine)
    if (isTRUE(refine) && !use_refine) {
        stop("refine = TRUE is only accepted when unmixing_method = \"AutoSpectral\".", call. = FALSE)
    }
    if (!identical(unmixing_method, "Spectreasy") && !spectreasy_weight_quantile_missing) {
        stop("spectreasy_weight_quantile is only accepted when unmixing_method = \"Spectreasy\".", call. = FALSE)
    }
    if (identical(unmixing_method, "Spectreasy")) {
        spectreasy_weight_quantile <- .normalize_spectreasy_weight_quantile(spectreasy_weight_quantile)
    }
    scc_background_args <- .validate_scc_background_args(
        clean_scc_with_unstained = use_autospectral && isTRUE(clean_scc_with_unstained),
        scc_background_method = scc_background_method,
        scc_background_k = scc_background_k
    )
    .with_optional_seed(seed)
    extra_args <- list(...)

    control_file <- .resolve_control_file_path(control_file)
    output_dir <- .normalize_unmix_controls_output_dir(output_dir)
    if (!gating_file_missing) {
        manual_gate_file <- gating_file
    }
    manual_gate_file_explicit <- !gating_file_missing || !manual_gate_file_missing

    if (!dir.exists(output_dir)) {
        dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
        if (!dir.exists(output_dir)) {
            Sys.sleep(0.5)
            dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
        }
    }
    if (!dir.exists(scc_dir)) .spectreasy_stop_missing_directory(scc_dir, label = "scc_dir")

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

    .spectreasy_console_header("unmix_controls")
    .spectreasy_console_field("SCC dir", .spectreasy_console_path(scc_dir))
    .spectreasy_console_field("Controls", .spectreasy_console_path(control_file))
    .spectreasy_console_field("Output", .spectreasy_console_path(output_dir))
    .spectreasy_console_field("Method", unmixing_method)
    .spectreasy_console_field("AF bands", af_n_bands)

    if (isTRUE(manual_gating)) {
        if (!interactive()) {
            if (manual_gate_file_explicit && length(manual_gate_file) > 0 && !is.na(manual_gate_file[1]) && nzchar(trimws(as.character(manual_gate_file)[1]))) {
                manual_gate_file <- normalizePath(as.character(manual_gate_file)[1], mustWork = FALSE)
                if (!file.exists(manual_gate_file)) {
                    stop("manual_gating = TRUE requires an interactive R session and the supplied gating_file does not exist: ", manual_gate_file, call. = FALSE)
                }
                warning("manual_gating = TRUE requires an interactive R session; reusing supplied gating_file instead.", call. = FALSE)
            } else {
                warning("manual_gating = TRUE requires an interactive R session; continuing with automatic gating.", call. = FALSE)
                manual_gate_file <- NULL
            }
        } else {
            manual_gate_file <- gate_controls(
                scc_dir = scc_dir,
                control_file = control_file,
                gate_file = manual_gate_file,
                open_browser = TRUE
            )
        }
    } else {
        if (length(manual_gate_file) == 0 || is.na(manual_gate_file[1]) || !nzchar(trimws(as.character(manual_gate_file)[1]))) {
            manual_gate_file <- NULL
        } else {
            manual_gate_file <- normalizePath(as.character(manual_gate_file)[1], mustWork = FALSE)
            if (!file.exists(manual_gate_file)) {
                if (manual_gate_file_explicit) {
                    stop("gating_file not found: ", manual_gate_file, call. = FALSE)
                }
                manual_gate_file <- NULL
            } else {
                .spectreasy_console_field("Gate CSV", .spectreasy_console_path(manual_gate_file))
            }
        }
    }

    output_paths <- .unmix_output_paths(output_dir)
    qc_controls_dir <- NULL
    output_file <- NULL
    if (isTRUE(save_report)) {
        qc_controls_dir <- .next_safe_output_dir(output_paths$qc_controls_dir)
        output_file <- file.path(qc_controls_dir, "qc_controls_report.pdf")
        .spectreasy_console_field("Report", .spectreasy_console_path(output_file))
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
        manual_gate_file = manual_gate_file,
        autospectral_scc_cleanup = use_autospectral,
        clean_scc_with_unstained = scc_background_args$enabled,
        scc_background_method = scc_background_args$method,
        scc_background_k = scc_background_args$k,
        autospectral_n_candidates = autospectral_n_candidates,
        autospectral_n_spectral = autospectral_n_spectral,
        autospectral_min_events = autospectral_min_events,
        refine = refine,
        seed = seed
    ), extra_args)
    M <- do.call(build_reference_matrix, build_args)
    if (is.null(M) || nrow(M) == 0) stop("No valid spectra found while building reference matrix.")

    spectral_variant_library <- NULL
    spectral_variant_library_file <- NULL
    if (use_autospectral) {
        spectral_variant_library <- tryCatch(
            .learn_spectral_variant_library(
                M = M,
                enabled = TRUE,
                som_nodes = spectral_variant_som_nodes,
                cosine_threshold = spectral_variant_cosine_threshold,
                max_variants = spectral_variant_max_variants,
                min_events = spectral_variant_min_events,
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
    .save_detector_noise_csv(attr(M, "detector_noise"), output_paths$detector_noise_csv)
    meta_info <- .read_unmix_metadata_pd(scc_dir)

    # Split M into fluorophore rows and AF/autofluorescence rows to plot separately
    af_rows <- grepl("^AF($|_)", rownames(M), ignore.case = TRUE)
    M_fluor <- M[!af_rows, , drop = FALSE]
    M_af <- M[af_rows, , drop = FALSE]

    p_spectra <- if (nrow(M_fluor) > 0) {
        .run_optional_unmix_artifact(
            "SCC spectra plot",
            plot_spectra(M_fluor, pd = meta_info$pd, output_file = if (isTRUE(save_qc_plots)) output_paths$spectra_file else NULL)
        )
    } else {
        NULL
    }

    p_af_spectra <- if (nrow(M_af) > 0) {
        .run_optional_unmix_artifact(
            "AF spectra plot",
            plot_spectra(M_af, pd = meta_info$pd, output_file = if (isTRUE(save_qc_plots)) output_paths$af_spectra_file else NULL)
        )
    } else {
        NULL
    }

    unmix_sample_args <- list(
        sample_dir = scc_dir,
        M = M,
        unmixing_method = unmixing_method,
        rwls_max_iter = rwls_max_iter,
        n_threads = n_threads,
        spectral_variant_library = spectral_variant_library,
        spectral_variant_top_k = spectral_variant_top_k,
        output_dir = output_paths$unmixed_dir,
        write_fcs = TRUE,
        save_report = FALSE
    )
    if (identical(unmixing_method, "Spectreasy")) {
        unmix_sample_args$spectreasy_weight_quantile <- spectreasy_weight_quantile
    }
    .spectreasy_console_field("Unmixing", paste0(length(meta_info$fcs_files), " control file(s)"))
    unmix_sample_args$verbose <- FALSE
    unmixed_list <- do.call(unmix_samples, unmix_sample_args)

    static_info <- .derive_unmix_static_matrix(M, fcs_files = meta_info$fcs_files, unmixing_method = unmixing_method)
    W <- static_info$W
    save_unmixing_matrix(W, output_paths$unmixing_matrix_csv)

    marker_mapping <- .resolve_unmix_marker_mappings(control_df)
    scatter_markers <- rownames(M)
    extra_af_rows <- grepl("^AF_", scatter_markers, ignore.case = TRUE)
    scatter_markers <- scatter_markers[!extra_af_rows]

    p_scatter <- .run_optional_unmix_artifact(
        "SCC unmixing scatter plot",
        plot_unmixing_scatter_matrix(
            unmixed_list = unmixed_list,
            sample_to_marker = marker_mapping$sample_to_marker,
            markers = scatter_markers,
            marker_display = NULL,
            output_file = if (isTRUE(save_qc_plots)) output_paths$unmixing_scatter_png else NULL,
            transform = "none",
            panel_size_mm = unmix_scatter_panel_size_mm,
            seed = seed
        )
    )

    qc_report <- NULL
    if (isTRUE(save_report)) {
        qc_report <- .run_optional_unmix_artifact(
            "Automatic SCC QC report",
            qc_controls(
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
        )
        if (!file.exists(output_file)) {
            warning("Automatic SCC QC report was requested but was not created at: ", output_file, ". Continuing unmixing.", call. = FALSE)
        }
    }

    invisible(list(
        M = M,
        W = W,
        unmixed_list = unmixed_list,
        reference_matrix_file = output_paths$reference_matrix_csv,
        detector_noise_file = output_paths$detector_noise_csv,
        spectral_variant_library = spectral_variant_library,
        spectral_variant_library_file = spectral_variant_library_file,
        spectral_variant_info = if (!is.null(spectral_variant_library)) spectral_variant_library$info else NULL,
        unmixing_matrix_file = output_paths$unmixing_matrix_csv,
        qc_report_file = if (isTRUE(save_report) && !is.null(output_file) && file.exists(output_file)) output_file else NULL,
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
        unmixing_scatter_plot = p_scatter,
        refine = refine
    ))
}
