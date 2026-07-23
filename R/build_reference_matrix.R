#' Build a Reference Matrix from Single-Color Controls
#'
#' Reads SCC FCS files, performs FSC/SSC gating and positive-peak intensity gating,
#' then computes one normalized spectrum per fluorophore.
#'
#' This is the core matrix-construction step used before unmixing.
#'
#' @param input_folder Directory containing SCC `.fcs` files.
#' @param output_folder Directory where gating/spectrum plots are written.
#' @param save_qc_plots Logical; if `TRUE`, write FSC/SSC, intensity-gate, and spectrum plots.
#'   When `FALSE` (default), the function returns the matrix without writing QC files.
#' @param control_df Optional control mapping as a data.frame or CSV path.
#'   Expected columns: `filename`, `fluorophore`, `channel`; `universal.negative` is optional.
#' @param af_profile Optional saved AF profile name, `spectreasy_af_profile`, or
#'   AF-only matrix. When supplied, its spectra replace extraction from mapped
#'   unstained cell controls.
#' @param af_n_bands Exact number of AF basis signatures to extract from pooled
#'   unstained/AF control events. The default is `100`. If the events cannot
#'   support that many distinct, non-empty bands, extraction stops with an error
#'   instead of returning a smaller bank.
#' @param af_max_cells Maximum number of scatter-gated AF events used when
#'   deriving AF basis signatures.
#' @param seed Optional integer seed for deterministic subsampling/clustering.
#' @param n_threads Positive integer; number of threads used for event-wise
#'   AutoSpectral AF assignment during AF-bank refinement.
#' @param default_sample_type Fallback type when filename heuristics are ambiguous (`"beads"` or `"cells"`).
#' @param cytometer Cytometer name used as a channel-mapping hint. The default,
#'   `"auto"`, infers the cytometer from FCS detector names when possible.
#' @param manual_gate_file Optional gate CSV from [gate_controls()]. When
#'   provided, manual cell/singlet gates are used before automatic SCC spectrum
#'   extraction, and manual positive gates are used for standard SCC extraction.
#' @param histogram_pct_beads Quantile width for the bead histogram gate.
#' @param histogram_direction_beads Histogram gate direction for beads: `"right"` starts at the median,
#'   `"both"` centers on the median, and `"left"` ends at the median.
#' @param histogram_pct_cells Quantile width for the cell histogram gate.
#' @param histogram_direction_cells Histogram gate direction for cells: `"right"` starts at the median,
#'   `"both"` centers on the median, and `"left"` ends at the median.
#' @param outlier_percentile Upper-tail FSC/SSC filtering percentile.
#' @param debris_percentile Lower FSC range fraction for cell/unstained controls.
#'   For example, `0.08` excludes events below 8\% of the per-file high FSC scale.
#' @param bead_gate_scale Ellipse scaling factor for bead FSC/SSC gate.
#' @param max_clusters Maximum number of GMM components tested.
#' @param min_cluster_proportion Minimum population proportion kept from GMM fit.
#' @param gate_contour_beads Contour probability for bead gating ellipse/hull.
#' @param gate_contour_cells Contour probability for cell gating ellipse/hull.
#' @param subsample_n Maximum number of events used for GMM fitting per file.
#' @param unmixing_method SCC processing method. `"AutoSpectral"` (default)
#'   first applies the saved positive histogram gate or the same automatic
#'   histogram fallback used by the other methods, then enables spectral SCC
#'   selection, resolved-negative background subtraction, and spectral-variant
#'   learning. Other methods stop after calculating the conventional
#'   histogram-gated reference spectrum.
#' @param scc_background_method Background subtraction method for AutoSpectral
#'   SCC cleanup (`"scatter_knn"` or `"none"`).
#' @param scc_background_k Number of nearest unstained/negative events averaged
#'   for scatter-matched SCC background subtraction.
#' @param autospectral_n_candidates Number of peak-bright SCC candidate events
#'   considered by the AutoSpectral-style selector.
#' @param autospectral_n_spectral Number of least-background-like SCC events
#'   kept for spectrum calculation by the AutoSpectral-style selector.
#' @param autospectral_min_events Minimum event count required by the
#'   AutoSpectral-style SCC selector.
#' @param refine Logical; if `TRUE`, refine the fixed-size k-means AF bank with
#'   native AutoSpectral-style unstained residual modulation. This is only
#'   supported with `unmixing_method = "AutoSpectral"`.
#'
#' @return Numeric matrix with rows = fluorophores and columns = detectors
#'   (normalized spectra). The matrix carries SCC-derived detector noise floors
#'   in `attr(M, "detector_noise")` for WLS/RWLS unmixing.
#' @export
#' @examples
#' if (interactive()) {
#'   M <- build_reference_matrix(
#'     input_folder = "scc",
#'     output_folder = "spectreasy_outputs/build_reference_plots",
#'     save_qc_plots = TRUE,
#'     control_df = "fcs_mapping.csv",
#'     cytometer = "auto"
#'   )
#'
#'   M_fast <- build_reference_matrix(
#'     input_folder = "scc",
#'     save_qc_plots = FALSE
#'   )
#' }
build_reference_matrix <- function(
  input_folder = "scc",
  output_folder = "gating_and_spectrum_plots",
  save_qc_plots = FALSE,
  control_df = NULL,
  af_profile = NULL,
  af_n_bands = 100,
  af_max_cells = 50000,
  seed = NULL,
  default_sample_type = "beads",
  cytometer = "auto",
  manual_gate_file = NULL,
  histogram_pct_beads = 0.98,
  histogram_direction_beads = "right",
  histogram_pct_cells = 0.35,
  histogram_direction_cells = "right",
  outlier_percentile = 0.02,
  debris_percentile = 0.08,
  bead_gate_scale = 1.3,
  max_clusters = 10,
  min_cluster_proportion = 0.03,
  gate_contour_beads = 0.95,
  gate_contour_cells = 0.90,
  subsample_n = 5000,
  unmixing_method = "AutoSpectral",
  scc_background_method = c("scatter_knn", "none"),
  scc_background_k = 2L,
  autospectral_n_candidates = 1000L,
  autospectral_n_spectral = 200L,
  autospectral_min_events = 10L,
  refine = FALSE,
  n_threads = 1L
) {
    control_df <- .normalize_build_reference_control_df(control_df)
    unmixing_method <- .normalize_unmix_method(unmixing_method)
    n_threads <- .normalize_n_threads(n_threads)
    default_sample_type <- .match_arg_ci(default_sample_type, c("beads", "cells"), "default_sample_type")
    histogram_direction_beads <- .match_arg_ci(
        histogram_direction_beads, c("right", "both", "left"), "histogram_direction_beads"
    )
    histogram_direction_cells <- .match_arg_ci(
        histogram_direction_cells, c("right", "both", "left"), "histogram_direction_cells"
    )
    histogram_pct_beads <- .normalize_unit_interval(histogram_pct_beads, "histogram_pct_beads")
    histogram_pct_cells <- .normalize_unit_interval(histogram_pct_cells, "histogram_pct_cells")
    outlier_percentile <- .normalize_unit_interval(
        outlier_percentile, "outlier_percentile", upper_inclusive = FALSE
    )
    debris_percentile <- .normalize_unit_interval(debris_percentile, "debris_percentile")
    min_cluster_proportion <- .normalize_unit_interval(
        min_cluster_proportion, "min_cluster_proportion", lower_inclusive = FALSE
    )
    gate_contour_beads <- .normalize_unit_interval(
        gate_contour_beads, "gate_contour_beads", lower_inclusive = FALSE, upper_inclusive = FALSE
    )
    gate_contour_cells <- .normalize_unit_interval(
        gate_contour_cells, "gate_contour_cells", lower_inclusive = FALSE, upper_inclusive = FALSE
    )
    bead_gate_scale <- .normalize_positive_number(bead_gate_scale, "bead_gate_scale")
    max_clusters <- .normalize_positive_integer(max_clusters, "max_clusters")
    subsample_n <- .normalize_positive_integer(subsample_n, "subsample_n")
    use_autospectral <- .is_autospectral_method(unmixing_method)
    refine <- .validate_reference_refine_arg(refine)
    if (isTRUE(refine) && !identical(unmixing_method, "AutoSpectral")) {
        stop("refine = TRUE is only supported with unmixing_method = \"AutoSpectral\".", call. = FALSE)
    }
    af_args <- .validate_build_reference_af_args(
        af_n_bands = af_n_bands,
        af_max_cells = af_max_cells
    )
    af_n_bands <- af_args$af_n_bands
    af_max_cells <- af_args$af_max_cells
    scc_background_args <- .validate_scc_background_args(
        scc_background_method = scc_background_method,
        scc_background_k = scc_background_k,
        enabled = use_autospectral
    )
    autospectral_n_candidates <- .normalize_positive_integer(autospectral_n_candidates, "autospectral_n_candidates")
    autospectral_n_spectral <- .normalize_positive_integer(autospectral_n_spectral, "autospectral_n_spectral")
    autospectral_min_events <- .normalize_positive_integer(autospectral_min_events, "autospectral_min_events")

    .with_optional_seed(seed)
    manual_gates <- .read_reference_manual_gates(manual_gate_file)

    sample_patterns <- get_fluorophore_patterns()
    file_info <- .prepare_reference_file_set(input_folder = input_folder, control_df = control_df)
    out_path <- .prepare_reference_output_path(output_folder = output_folder, save_qc_plots = save_qc_plots)
    metadata <- .prepare_reference_detector_info(file_info$fcs_files[1])
    cytometer <- .resolve_cytometer_from_pd(cytometer, metadata$pd_meta)
    gui_scatter_domains <- if (isTRUE(save_qc_plots)) {
        .reference_compute_gui_scatter_domains(file_info$fcs_files_all, max_points = 100000L)
    } else {
        NULL
    }

    config <- list(
        default_sample_type = default_sample_type,
        outlier_percentile = outlier_percentile,
        debris_percentile = debris_percentile,
        subsample_n = subsample_n,
        max_clusters = max_clusters,
        min_cluster_proportion = min_cluster_proportion,
        gate_contour_beads = gate_contour_beads,
        gate_contour_cells = gate_contour_cells,
        bead_gate_scale = bead_gate_scale,
        manual_gates = manual_gates,
        histogram_pct_beads = histogram_pct_beads,
        histogram_direction_beads = histogram_direction_beads,
        histogram_pct_cells = histogram_pct_cells,
        histogram_direction_cells = histogram_direction_cells,
        save_qc_plots = save_qc_plots,
        out_path = out_path,
        cytometer = cytometer,
        spectral_scc_pipeline = use_autospectral,
        scc_background_enabled = scc_background_args$enabled,
        scc_background_method = scc_background_args$method,
        scc_background_k = scc_background_args$k,
        autospectral_n_candidates = autospectral_n_candidates,
        autospectral_n_spectral = autospectral_n_spectral,
        autospectral_min_events = autospectral_min_events,
        refine = refine,
        gui_scatter_domains = gui_scatter_domains
    )

    .spectreasy_console_field("Detectors", paste0(length(metadata$detector_names), " spectral channel(s), sorted by laser"))
    .validate_reference_detector_consistency(
        fcs_files = file_info$fcs_files_all,
        detector_names = metadata$detector_names
    )

    af_profiles <- if (!is.null(af_profile)) {
        profile_name <- if (is.character(af_profile) && length(af_profile) == 1L) af_profile else "saved AF profile"
        profile_object <- if (is.character(af_profile) && length(af_profile) == 1L) {
            load_af_profile(af_profile, show_plot = FALSE)
        } else {
            af_profile
        }
        profile_matrix <- .coerce_af_profile_matrix(profile_object, arg_name = "af_profile")
        if (!setequal(colnames(profile_matrix), metadata$detector_names)) {
            missing_detectors <- setdiff(metadata$detector_names, colnames(profile_matrix))
            extra_detectors <- setdiff(colnames(profile_matrix), metadata$detector_names)
            stop(
                "Saved AF profile detectors do not match the SCC detector set.",
                if (length(missing_detectors) > 0L) paste0(" Missing: ", paste(missing_detectors, collapse = ", "), ".") else "",
                if (length(extra_detectors) > 0L) paste0(" Extra: ", paste(extra_detectors, collapse = ", "), ".") else "",
                call. = FALSE
            )
        }
        profile_matrix <- profile_matrix[, metadata$detector_names, drop = FALSE]
        profile_background <- if (.is_af_profile_object(profile_object)) profile_object$scc_background else NULL
        profile_raw_median <- if (.is_af_profile_object(profile_object)) profile_object$raw_median else NULL
        saved_profiles <- list(
            af_data_raw = profile_raw_median,
            af_signatures_norm = profile_matrix,
            af_bank_info = list(
                source_count = 0L,
                sources = data.frame(),
                pooled_events = 0L,
                requested_bands = nrow(profile_matrix),
                derived_bands = nrow(profile_matrix),
                mode = "saved_profile",
                profile_name = profile_name
            ),
            scc_background = profile_background,
            af_events = if (!is.null(profile_background$spectra)) profile_background$spectra else NULL
        )
        .spectreasy_console_field("AF bank", paste0(nrow(profile_matrix), " saved signature(s) from ", profile_name))
        saved_profiles
    } else {
        .collect_reference_af_profiles(
            control_df = control_df,
            fcs_files = file_info$fcs_files,
            fcs_files_all = file_info$fcs_files_all,
            detector_names = metadata$detector_names,
            af_n_bands = af_n_bands,
            af_max_cells = af_max_cells,
            config = config
        )
    }

    universal_negatives <- .collect_reference_universal_negatives(
        control_df = control_df,
        fcs_files = file_info$fcs_files_all,
        detector_names = metadata$detector_names,
        sample_patterns = sample_patterns,
        config = config
    )
    bead_negative <- .collect_reference_unstained_bead_negative(
        control_df = control_df,
        fcs_files = file_info$fcs_files_all,
        detector_names = metadata$detector_names,
        config = config
    )

    results_list <- list()
    qc_summary_list <- list()
    for (fcs_file in file_info$fcs_files_all) {
        processed <- .process_reference_file(
            fcs_file = fcs_file,
            control_df = control_df,
            sample_patterns = sample_patterns,
            metadata = metadata,
            config = config,
            af_data_raw = af_profiles$af_data_raw,
            universal_negatives = universal_negatives,
            bead_negative = bead_negative,
            scc_background = if (isTRUE(config$spectral_scc_pipeline) && isTRUE(config$scc_background_enabled)) {
                af_profiles$scc_background
            } else {
                NULL
            }
        )
        if (is.null(processed)) {
            next
        }
        results_list[[processed$sample_name]] <- processed$result
        qc_summary_list[[processed$sample_name]] <- processed$qc_summary
    }
    .validate_reference_complete_controls(
        control_df = control_df,
        fcs_files_all = file_info$fcs_files_all,
        processed_results = results_list
    )

    M <- .finalize_reference_matrix(
        results_list = results_list,
        qc_summary_list = qc_summary_list,
        af_signatures_norm = af_profiles$af_signatures_norm,
        af_data_raw = af_profiles$af_data_raw,
        af_bank_info = af_profiles$af_bank_info,
        detector_names = metadata$detector_names,
        pd_meta = metadata$pd_meta,
        save_qc_plots = save_qc_plots,
        out_path = out_path
    )
    if (is.null(M)) {
        return(M)
    }
    if (isTRUE(refine)) {
        refined_af <- .reference_refine_af_bank(
            M = M,
            af_events = af_profiles$af_events,
            af_n_bands = af_n_bands,
            af_max_cells = af_max_cells,
            n_threads = n_threads,
            seed = seed,
            verbose = TRUE
        )
        if (!is.null(refined_af) && !is.null(refined_af$signatures) && nrow(refined_af$signatures) > 0L) {
            af_bank_info <- attr(M, "af_bank_info")
            if (is.null(af_bank_info)) af_bank_info <- list()
            af_bank_info$requested_bands <- af_n_bands
            af_bank_info$derived_bands <- nrow(refined_af$signatures)
            af_bank_info$selection <- refined_af$selection
            af_bank_info$refine <- TRUE
            M <- .replace_reference_af_rows(
                M = M,
                af_signatures_norm = refined_af$signatures,
                af_bank_info = af_bank_info
            )
            .spectreasy_console_field("AF refine", paste0(nrow(refined_af$signatures), " fixed signature(s)"))
        } else {
            af_bank_info <- attr(M, "af_bank_info")
            if (!is.null(af_bank_info)) {
                af_bank_info$refine <- FALSE
                af_bank_info$refine_reason <- "no_refined_af_candidates"
                attr(M, "af_bank_info") <- af_bank_info
            }
            .spectreasy_console_field("AF refine", "kept base AF bank")
        }
    }
    .attach_estimated_wls_detector_noise(
        M = M,
        scc_dir = input_folder,
        fcs_files = file_info$fcs_files_all,
        fallback = .default_wls_background_noise(),
        warn = FALSE
    )
}
