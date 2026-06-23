# spectreasy user function reference
#
# This file lists the user-facing functions most people touch in a normal
# spectreasy workflow, with every default argument written out.
#
# It is wrapped in if (FALSE) so sourcing this file is safe. Copy the call you
# need, then replace placeholder objects such as flow_frame, M, results, and W.

if (FALSE) {
    library(spectreasy)

    # -------------------------------------------------------------------------
    # Typical simple workflow
    # -------------------------------------------------------------------------

    control_file <- create_control_file(
        input_folder = "scc",
        af_folder = "af",
        include_af_folder = TRUE,
        cytometer = "Aurora",
        default_control_type = "cells",
        unknown_fluor_policy = c("empty", "by_channel", "filename"),
        output_file = "fcs_mapping.csv",
        custom_fluorophores = NULL
    )

    control_check <- validate_control_file_mapping(
        control_df = control_file,
        scc_dir = "scc",
        include_multi_af = FALSE,
        exclude_af = FALSE,
        af_dir = "af",
        require_all_scc_mapped = TRUE,
        require_channels = TRUE,
        stop_on_error = FALSE
    )

    gate_settings <- gating_options(
        use_scatter_gating = TRUE,
        histogram_pct_beads = 0.98,
        histogram_direction_beads = "right",
        histogram_pct_cells = 0.35,
        histogram_direction_cells = "right"
    )

    ctrl <- unmix_controls(
        scc_dir = "scc",
        control_file = "fcs_mapping.csv",
        auto_create_control = TRUE,
        cytometer = "Aurora",
        auto_default_control_type = "beads",
        auto_unknown_fluor_policy = c("by_channel", "empty", "filename"),
        output_dir = "spectreasy_outputs/unmix_controls",
        exclude_af = FALSE,
        unmix_method = "WLS",
        unmix_scatter_panel_size_mm = 30,
        seed = NULL,
        af_n_bands = 10,
        af_bands_per_file = 5,
        include_multi_af = FALSE,
        rwls_max_iter = 1L,
        unmix_threads = 1L,
        save_qc_plots = FALSE,
        use_scatter_gating = TRUE
        # Additional build_reference_matrix() arguments can be passed here.
    )

    qc_controls_report <- qc_controls(
        M = NULL,
        unmixing_matrix_file = file.path(
            "spectreasy_outputs",
            "unmix_controls",
            "scc_reference_matrix.csv"
        ),
        scc_dir = "scc",
        output_file = "spectreasy_outputs/unmix_controls/qc_controls_report.pdf",
        control_file = "fcs_mapping.csv",
        cytometer = "Aurora",
        method = "WLS",
        qc_plot_dir = file.path("spectreasy_outputs", "scc_report_plots"),
        save_qc_pngs = FALSE,
        use_scatter_gating = TRUE,
        include_multi_af = FALSE,
        af_dir = "af",
        af_bands_per_file = 5,
        unmix_scatter_max_points = 1000,
        unmix_scatter_axis_limit = NULL,
        seed = NULL
        # Additional build_reference_matrix() arguments can be passed here.
    )

    launch_gui(
        matrix_dir = NULL,
        samples_dir = NULL,
        port = 8000,
        open_browser = TRUE,
        dev_mode = FALSE
    )

    unmixed <- unmix_samples(
        sample_dir = "samples",
        M = NULL,
        unmixing_matrix_file = file.path(
            "spectreasy_outputs",
            "unmix_controls",
            "scc_reference_matrix.csv"
        ),
        variances_file = file.path(
            "spectreasy_outputs",
            "unmix_controls",
            "scc_variances.csv"
        ),
        detector_noise_file = NULL,
        method = "WLS",
        rwls_max_iter = 1L,
        n_threads = 1L,
        cytometer = "Aurora",
        scc_dir = NULL,
        control_file = NULL,
        af_n_bands = 10,
        af_bands_per_file = 5,
        exclude_af = FALSE,
        include_multi_af = FALSE,
        output_dir = file.path(
            "spectreasy_outputs",
            "unmix_samples",
            "unmixed_fcs"
        ),
        write_fcs = TRUE,
        subsample_n = NULL,
        seed = NULL,
        return_type = c("list", "flowSet", "SingleCellExperiment"),
        verbose = TRUE
    )

    qc_samples_report <- qc_samples(
        results = unmixed,
        M = NULL,
        unmixing_matrix_file = file.path(
            "spectreasy_outputs",
            "unmix_controls",
            "scc_reference_matrix.csv"
        ),
        output_file = "spectreasy_outputs/unmix_samples/qc_samples_report.pdf",
        method = NULL,
        res_list = NULL,
        png_dir = NULL,
        pd = NULL,
        max_events_per_sample = 1000,
        overview_files_per_page = 15,
        matrix_markers_per_page = 15,
        sample_nxn_rows_per_page = 10,
        sample_nxn_max_points = 1000,
        sample_nxn_transform = c("none", "asinh"),
        sample_nxn_asinh_cofactor = 150,
        sample_nxn_axis_limit = NULL,
        nxn_all_samples = FALSE
    )

    # -------------------------------------------------------------------------
    # Autofluorescence profile library workflow
    # -------------------------------------------------------------------------

    profile <- extract_af_profile(
        fcs_file = "af/unstained.fcs",
        af_n_bands = "auto",
        af_max_cells = 50000,
        seed = NULL,
        show_plot = TRUE,
        verbose = TRUE
    )

    save_af_profile(
        name = "my_af_profile",
        x = profile,
        overwrite = FALSE
    )

    profiles <- list_af_profiles()

    profile_directory <- af_profile_dir(
        create = TRUE
    )

    loaded_profile <- load_af_profile(
        name = "my_af_profile",
        show_plot = FALSE
    )

    profile_plot <- plot_af_profile(
        x = loaded_profile,
        show = TRUE
    )

    M_with_profile <- add_af_profile(
        M = ctrl$M,
        profile = loaded_profile,
        replace_existing = TRUE
    )

    delete_af_profile(
        name = "my_af_profile"
    )

    # -------------------------------------------------------------------------
    # Lower-level reference matrix and unmixing helpers
    # -------------------------------------------------------------------------

    M <- build_reference_matrix(
        input_folder = "scc",
        output_folder = "gating_and_spectrum_plots",
        save_qc_plots = FALSE,
        control_df = NULL,
        include_multi_af = FALSE,
        exclude_af = FALSE,
        af_dir = "af",
        af_n_bands = 10,
        af_bands_per_file = 5,
        af_max_cells = 50000,
        seed = NULL,
        default_sample_type = "beads",
        cytometer = "Aurora",
        use_scatter_gating = TRUE,
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
        gate_contour_cells = 0.9,
        subsample_n = 5000
    )

    residuals <- calc_residuals(
        flow_frame = flow_frame,
        M = M,
        file_name = NULL,
        method = "WLS",
        return_residuals = FALSE,
        background_noise = 125,
        wls_signal_scale = 1,
        wls_max_weight_ratio = 1600,
        rwls_max_iter = 1L,
        n_threads = 1L
    )

    W <- derive_unmixing_matrix(
        M = M,
        method = "OLS",
        variances = NULL,
        background_noise = 125,
        wls_signal_scale = 1,
        wls_max_weight_ratio = 1600
    )

    save_unmixing_matrix(
        W = W,
        file = "unmixing_matrix.csv"
    )

    spectra_plot <- plot_spectra(
        ref_matrix = M,
        pd = NULL,
        output_file = NULL,
        width = 250,
        height = 100,
        unit = "mm",
        dpi = 600,
        theme_custom = NULL
    )

    unmixing_matrix_plot <- plot_unmixing_matrix(
        W = W,
        pd = NULL
    )

    scatter_plot <- plot_unmixing_scatter_matrix(
        unmixed_list = unmixed,
        sample_to_marker = NULL,
        markers = NULL,
        marker_display = NULL,
        output_file = NULL,
        max_points_per_sample = 1000,
        transform = c("none", "asinh"),
        asinh_cofactor = 150,
        axis_limit = NULL,
        panel_size_mm = 30,
        seed = NULL
    )

    detector_residual_plot <- plot_detector_residuals(
        res_list = residuals,
        M = M,
        top_n = 50,
        output_file = NULL,
        width = 250,
        height = 120,
        pd = NULL
    )

    nps <- calculate_nps(
        data = as.data.frame(unmixed),
        markers = NULL
    )

    nps_plot <- plot_nps(
        nps_results = nps,
        output_file = NULL,
        width = 200
    )

    ssm <- calculate_ssm(
        M = M,
        method = "OLS"
    )

    ssm_plot <- plot_ssm(
        SSM = ssm,
        output_file = NULL,
        width = 200,
        height = 180
    )

    spectra <- get_control_spectra(
        flow_frame = flow_frame,
        control_file = "fcs_mapping.csv",
        control_dir = "scc",
        af_dir = "af",
        method = "WLS",
        cytometer = "Aurora"
    )

    # -------------------------------------------------------------------------
    # Small utility helpers users may inspect or call
    # -------------------------------------------------------------------------

    positive_cells <- gate_positive_cells(
        mat = matrix_or_data_frame,
        histogram_pct = 0.98,
        histogram_direction = "right"
    )

    fluorophore_patterns <- get_fluorophore_patterns()

    sorted_detectors <- get_sorted_detectors(
        pd = parameter_data
    )

    example_data <- spectreasy_example_data(
        asset = NULL,
        cache_dir = file.path(
            tools::R_user_dir("spectreasy", which = "cache"),
            "example-data"
        ),
        dest_dir = NULL,
        force = FALSE,
        quiet = FALSE
    )
}
