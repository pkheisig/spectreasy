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
#' @param output_dir Root output directory. Control-stage artifacts are written
#'   under `output_dir/unmix_controls`, including the automatic QC report and
#'   its HTML NxN companion.
#' @param unmixing_method SCC unmixing method (`"WLS"`, `"RWLS"`,
#'   `"OLS"`, `"NNLS"`, `"AutoSpectral"`, or `"Spectreasy"`).
#'   `AutoSpectral` keeps the k-means AF bank controlled by `af_n_bands`,
#'   then uses OLS after AF matching plus AutoSpectral-style SCC cleanup and
#'   spectral variants. Both spectral methods first restrict stained SCC events
#'   with the saved positive histogram gate or the automatic histogram fallback.
#'   `Spectreasy` uses the same cleanup and variant
#'   machinery, then blends the AutoSpectral-style OLS marker fit with a
#'   marker-only OLS anchor using decoder-projected AF impact weights.
#' @param unmix_scatter_panel_size_mm Panel size for SCC unmixing scatter matrix plot.
#' @param seed Optional integer seed for deterministic subsampling and plotting.
#' @param af_profile Optional saved AF profile name, `spectreasy_af_profile`, or
#'   AF-only matrix forwarded explicitly to [build_reference_matrix()].
#' @param af_n_bands Exact number of AF basis signatures to extract from pooled
#'   unstained/AF control events. The default is `100`. If the events cannot
#'   support that many distinct, non-empty bands, extraction stops with an error
#'   instead of returning a smaller bank.
#' @param rwls_max_iter Positive integer; number of robust reweighting
#'   iterations used when `unmixing_method = "RWLS"`. The default, 1, preserves the
#'   historical behavior.
#' @param n_threads Positive integer; number of threads for event-wise
#'   AutoSpectral/Spectreasy AF assignment and NNLS/WLS/RWLS SCC fitting. OLS
#'   uses vectorized matrix operations. The default, 1, keeps explicit event
#'   loops single-threaded.
#' @param save_qc_png Logical; whether to write per-control FSC/SSC,
#'   intensity-gate, and spectrum PNGs under `output_dir/unmix_controls`.
#' @param save_report Logical; if `TRUE`, write the SCC QC report, its HTML NxN
#'   companion, and supporting QC metrics under
#'   `output_dir/unmix_controls/qc_controls` after controls are unmixed.
#'   Defaults to `TRUE`.
#' @param report_format Report format, either `"html"` (default) or `"pdf"`.
#'   Only the selected format is written. Matching is case-insensitive.
#' @param gating_mode Control-gating mode. `"interactive"` (default) opens the
#'   gating GUI and preloads `manual_gate_file` when it exists; `"reuse"` skips
#'   the GUI and requires an existing gate CSV; `"automatic"` ignores any gate
#'   CSV and uses automatic GMM gating. Matching is case-insensitive.
#' @param manual_gate_file Gate CSV written by [gate_controls()]. Relative paths
#'   are resolved from the current working directory. In `"interactive"` mode
#'   the file is preloaded when it exists; `"reuse"` requires it to exist.
#' @param scc_background_method Background subtraction method for Spectreasy/
#'   AutoSpectral SCC cleanup (`"scatter_knn"` or `"none"`). The default enables
#'   scatter-matched subtraction from the resolved external or internal negative
#'   source; use `"none"` only to opt out explicitly.
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
#'   marker-specific AutoSpectral mixing. The default is `0.65`.
#' @param autospectral_n_candidates Number of peak-bright SCC candidate events
#'   considered by the AutoSpectral-style selector.
#' @param autospectral_n_spectral Number of least-background-like SCC events
#'   kept for spectrum calculation by the AutoSpectral-style selector.
#' @param autospectral_min_events Minimum event count required by the
#'   AutoSpectral-style SCC selector.
#' @param autospectral_refine Logical; when `unmixing_method = "AutoSpectral"`, refine the
#'   fixed-size k-means AF bank with native AutoSpectral-style unstained residual
#'   modulation. `TRUE` is rejected for all other unmixing methods.
#' @param project_path Project directory recorded in generated report metadata.
#' @param ... Additional arguments forwarded to [build_reference_matrix()].
#'
#' @return List with `M`, `W`, `unmixed_list`, the effective `gating_mode`, the
#'   reused or created `manual_gate_file` (when applicable), and key output paths,
#'   including `detector_noise_file` for the SCC-derived WLS noise floors.
#' @export
#' @examples
#' if (interactive()) {
#'   ctrl <- unmix_controls(
#'     scc_dir = "scc",
#'     control_file = "fcs_mapping.csv",
#'     auto_create_mapping = TRUE,
#'     cytometer = "auto",
#'     output_dir = "spectreasy_outputs"
#'   )
#'   ctrl$unmixing_matrix_file
#' }
unmix_controls <- function(
    scc_dir = "scc",
    control_file = "fcs_mapping.csv",
    auto_create_mapping = TRUE,
    cytometer = "auto",
    auto_unknown_fluor_policy = c("by_channel", "empty", "filename"),
    output_dir = "spectreasy_outputs",
    unmixing_method = "AutoSpectral",
    unmix_scatter_panel_size_mm = 30,
    seed = NULL,
    af_profile = NULL,
    af_n_bands = 100,
    rwls_max_iter = 1L,
    n_threads = 1L,
    save_qc_png = FALSE,
    save_report = TRUE,
    report_format = "html",
    gating_mode = "interactive",
    manual_gate_file = "ssc_gate_config.csv",
    scc_background_method = c("scatter_knn", "none"),
    scc_background_k = 2L,
    spectral_variant_som_nodes = 16L,
    spectral_variant_top_k = 3L,
    spectral_variant_cosine_threshold = 0.98,
    spectral_variant_max_variants = 8L,
    spectral_variant_min_events = 50L,
    spectreasy_weight_quantile = 0.65,
    autospectral_n_candidates = 1000L,
    autospectral_n_spectral = 200L,
    autospectral_min_events = 10L,
    autospectral_refine = FALSE,
    project_path = getwd(),
    ...
) {
    spectreasy_weight_quantile_missing <- missing(spectreasy_weight_quantile)
    manual_gate_file_missing <- missing(manual_gate_file)
    report_format <- .match_arg_ci(report_format, c("html", "pdf"), "report_format")
    gating_mode <- .match_arg_ci(
        gating_mode,
        c("interactive", "reuse", "automatic"),
        "gating_mode"
    )
    n_threads <- .normalize_n_threads(n_threads)
    save_qc_png <- .normalize_scalar_logical(save_qc_png, "save_qc_png")
    auto_unknown_fluor_policy <- .match_arg_ci(
        auto_unknown_fluor_policy,
        c("by_channel", "empty", "filename"),
        "auto_unknown_fluor_policy"
    )
    unmixing_method <- .normalize_unmix_method(unmixing_method)
    use_autospectral <- .is_autospectral_style_method(unmixing_method)
    use_refine <- identical(unmixing_method, "AutoSpectral")
    autospectral_refine <- .validate_reference_refine_arg(autospectral_refine)
    if (isTRUE(autospectral_refine) && !use_refine) {
        stop("autospectral_refine = TRUE is only accepted when unmixing_method = \"AutoSpectral\".", call. = FALSE)
    }
    if (!identical(unmixing_method, "Spectreasy") && !spectreasy_weight_quantile_missing) {
        stop("spectreasy_weight_quantile is only accepted when unmixing_method = \"Spectreasy\".", call. = FALSE)
    }
    if (identical(unmixing_method, "Spectreasy")) {
        spectreasy_weight_quantile <- .normalize_spectreasy_weight_quantile(spectreasy_weight_quantile)
    }
    scc_background_args <- .validate_scc_background_args(
        scc_background_method = scc_background_method,
        scc_background_k = scc_background_k,
        enabled = use_autospectral
    )
    .with_optional_seed(seed)
    extra_args <- list(...)
    removed_args <- intersect(names(extra_args), c("unmix_threads", "save_qc_plots", "refine"))
    if (length(removed_args) > 0L) {
        stop(
            "unmix_controls() no longer accepts: ", paste(removed_args, collapse = ", "),
            ". Use n_threads, save_qc_png, and autospectral_refine.",
            call. = FALSE
        )
    }
    if ("manual_gating" %in% names(extra_args)) {
        stop(
            "manual_gating has been replaced by gating_mode. Use ",
            "gating_mode = 'interactive', 'reuse', or 'automatic'.",
            call. = FALSE
        )
    }

    control_file <- .resolve_control_file_path(control_file)
    if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) || !nzchar(trimws(output_dir))) {
        stop("output_dir must be a non-empty root directory path.", call. = FALSE)
    }
    if (file.exists(output_dir) && !dir.exists(output_dir)) {
        stop("output_dir points to a file, not a directory: ", output_dir, call. = FALSE)
    }
    output_dir <- .unmix_controls_stage_dir(output_dir)
    manual_gate_file_explicit <- !manual_gate_file_missing

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

    gating <- .resolve_unmix_controls_gating(
        gating_mode = gating_mode,
        manual_gate_file = manual_gate_file,
        manual_gate_file_explicit = manual_gate_file_explicit,
        scc_dir = scc_dir,
        control_file = control_file
    )
    gating_mode_used <- gating$mode
    manual_gate_file <- gating$file

    output_paths <- .unmix_output_paths(output_dir)
    tryCatch(
        utils::write.csv(
            control_df,
            output_paths$control_mapping_csv,
            row.names = FALSE,
            quote = TRUE,
            na = ""
        ),
        error = function(e) {
            stop("Could not save the control mapping used for this run: ", conditionMessage(e), call. = FALSE)
        }
    )
    report_paths <- .unmix_controls_report_paths(output_paths, save_report, report_format)
    M <- .build_unmix_controls_reference(
        scc_dir = scc_dir,
        output_dir = output_dir,
        save_qc_png = save_qc_png,
        save_report = save_report,
        control_df = control_df,
        cytometer = cytometer,
        af_profile = af_profile,
        af_n_bands = af_n_bands,
        manual_gate_file = manual_gate_file,
        unmixing_method = unmixing_method,
        scc_background_args = scc_background_args,
        autospectral_n_candidates = autospectral_n_candidates,
        autospectral_n_spectral = autospectral_n_spectral,
        autospectral_min_events = autospectral_min_events,
        autospectral_refine = autospectral_refine,
        n_threads = n_threads,
        seed = seed,
        extra_args = extra_args
    )
    attr(M, "spectreasy_control_file") <- normalizePath(control_file, mustWork = TRUE)
    attr(M, "spectreasy_control_df") <- control_df

    variants <- .learn_unmix_controls_variants(
        M = M,
        enabled = use_autospectral,
        output_file = output_paths$spectral_variant_library_rds,
        som_nodes = spectral_variant_som_nodes,
        cosine_threshold = spectral_variant_cosine_threshold,
        max_variants = spectral_variant_max_variants,
        min_events = spectral_variant_min_events,
        seed = seed
    )
    M <- variants$M
    spectral_variant_library <- variants$library
    spectral_variant_library_file <- variants$file

    .save_reference_matrix_csv(M, output_paths$reference_matrix_csv)
    .save_detector_noise_csv(attr(M, "detector_noise"), output_paths$detector_noise_csv)
    meta_info <- .read_unmix_metadata_pd(scc_dir, control_df = control_df)

    plots <- .plot_unmix_controls_spectra(M, meta_info$pd, save_qc_png, output_paths)
    unmixed_controls_dir <- .next_safe_output_dir(output_paths$unmixed_dir)
    unmixed_list <- .unmix_control_files(
        meta_info = meta_info,
        M = M,
        control_file = control_file,
        unmixing_method = unmixing_method,
        rwls_max_iter = rwls_max_iter,
        n_threads = n_threads,
        spectral_variant_library = spectral_variant_library,
        spectral_variant_top_k = spectral_variant_top_k,
        spectreasy_weight_quantile = spectreasy_weight_quantile,
        unmixed_controls_dir = unmixed_controls_dir
    )
    static_info <- .derive_unmix_static_matrix(
        M,
        fcs_files = meta_info$fcs_files,
        unmixing_method = unmixing_method
    )
    W <- static_info$W
    save_unmixing_matrix(W, output_paths$unmixing_matrix_csv)
    plots$scatter <- .unmix_controls_scatter_plot(
        M = M,
        unmixed_list = unmixed_list,
        control_df = control_df,
        save_qc_png = save_qc_png,
        output_paths = output_paths,
        panel_size = unmix_scatter_panel_size_mm,
        seed = seed
    )
    qc_report <- .run_unmix_controls_report(
        save_report = save_report,
        M = M,
        scc_dir = scc_dir,
        report_paths = report_paths,
        control_file = control_file,
        cytometer = cytometer,
        unmixing_method = unmixing_method,
        save_qc_png = save_qc_png,
        unmixed_list = unmixed_list,
        meta_info = meta_info,
        seed = seed,
        output_paths = output_paths,
        spectral_variant_library_file = spectral_variant_library_file,
        manual_gate_file = manual_gate_file,
        report_format = report_format,
        report_plots = list(
            spectra = plots$spectra,
            af = plots$af_spectra,
            unmixing_scatter = plots$scatter
        ),
        run_settings = list(
            auto_create_mapping = auto_create_mapping,
            auto_unknown_fluor_policy = auto_unknown_fluor_policy,
            unmix_scatter_panel_size_mm = unmix_scatter_panel_size_mm,
            seed = seed,
            af_n_bands = af_n_bands,
            rwls_max_iter = rwls_max_iter,
            n_threads = n_threads,
            save_qc_png = save_qc_png,
            gating_mode = gating_mode_used,
            scc_background_method = scc_background_args$method,
            scc_background_k = scc_background_args$k,
            spectral_variant_som_nodes = spectral_variant_som_nodes,
            spectral_variant_top_k = spectral_variant_top_k,
            spectral_variant_cosine_threshold = spectral_variant_cosine_threshold,
            spectral_variant_max_variants = spectral_variant_max_variants,
            spectral_variant_min_events = spectral_variant_min_events,
            spectreasy_weight_quantile = spectreasy_weight_quantile,
            autospectral_n_candidates = autospectral_n_candidates,
            autospectral_n_spectral = autospectral_n_spectral,
            autospectral_min_events = autospectral_min_events,
            autospectral_refine = autospectral_refine
        ),
        project_path = project_path
    )
    .unmix_controls_result(
        M = M,
        W = W,
        unmixed_list = unmixed_list,
        output_paths = output_paths,
        control_file = control_file,
        report = qc_report,
        report_paths = report_paths,
        save_report = save_report,
        save_qc_png = save_qc_png,
        output_dir = output_dir,
        unmixed_controls_dir = unmixed_controls_dir,
        gating_mode = gating_mode_used,
        manual_gate_file = manual_gate_file,
        static_info = static_info,
        plots = plots,
        spectral_variant_library = spectral_variant_library,
        spectral_variant_library_file = spectral_variant_library_file,
        autospectral_refine = autospectral_refine
    )

}
