    # -------------------------------------------------------------------------
    # Typical simple workflow
    # -------------------------------------------------------------------------

    control_file <- create_control_file(
        input_folder = "scc",
        cytometer = "auto",
        unknown_fluor_policy = c("empty", "by_channel", "filename"),
        output_file = "fcs_mapping.csv",
        custom_fluorophores = NULL
    )

    control_check <- validate_control_file_mapping(
        control_df = control_file,
        scc_dir = "scc",
        require_all_scc_mapped = TRUE,
        require_channels = TRUE,
        stop_on_error = FALSE
    )

    gate_settings <- gating_options(
        histogram_pct_beads = 0.98,
        histogram_direction_beads = "right",
        histogram_pct_cells = 0.35,
        histogram_direction_cells = "right"
    )

    ctrl <- unmix_controls(
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
        manual_gating = TRUE,
        manual_gate_file = "ssc_gate_config.csv",
        gating_file = manual_gate_file,
        rwls_max_iter = 1L,
        n_threads = 1L,
        save_qc_plots = FALSE,
        save_report = TRUE,
        scc_background_method = c("scatter_knn", "none"),
        scc_background_k = 2L,
        spectral_variant_som_nodes = 16L,
        spectral_variant_top_k = 3L,
        spectral_variant_cosine_threshold = 0.98,
        spectral_variant_max_variants = 8L,
        spectral_variant_min_events = 50L,
        # Spectreasy only:
        spectreasy_weight_quantile = 0.9,
        autospectral_n_candidates = 1000L,
        autospectral_n_spectral = 200L,
        autospectral_min_events = 10L,
        refine = FALSE
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
        cytometer = "auto",
        unmixing_method = "WLS",
        qc_plot_dir = file.path("spectreasy_outputs", "scc_report_plots"),
        save_qc_pngs = FALSE,
        qc_metrics_dir = NULL,
        unmix_scatter_max_points = 1000,
        unmix_scatter_axis_limit = NULL,
        seed = NULL
        # Additional build_reference_matrix() arguments can be passed here.
    )

    unmixed <- unmix_samples(
        sample_dir = "samples",
        M = NULL,
        unmixing_matrix_file = file.path(
            "spectreasy_outputs",
            "unmix_controls",
            "scc_reference_matrix.csv"
        ),
        detector_noise_file = NULL,
        unmixing_method = "Spectreasy",
        rwls_max_iter = 1L,
        n_threads = 1L,
        spectral_variant_library = NULL,
        spectral_variant_library_file = NULL,
        spectral_variant_top_k = 3L,
        spectral_variant_min_abundance = 1,
        spectral_variant_positive_fraction = 0.02,
        spectral_variant_min_improvement = 0.01,
        # Spectreasy only:
        spectreasy_weight_quantile = 0.9,
        estimate_af = FALSE,
        output_dir = file.path(
            "spectreasy_outputs",
            "unmix_samples",
            "unmixed_fcs"
        ),
        write_fcs = TRUE,
        save_report = TRUE,
        save_qc_plots = FALSE,
        qc_plot_dir = NULL,
        subsample_n = NULL,
        seed = NULL,
        return_type = c("list", "flowSet", "SingleCellExperiment"),
        verbose = TRUE
    )

    # Spectreasy is the default unmixing method. It builds on the AutoSpectral
    # approach with AF-band assignment and spectral variants, then reweights
    # each marker between an AF-aware fit and a marker-only OLS anchor.
    # Lower spectreasy_weight_quantile values increase the marker blend weights;
    # higher values make the blend more conservative.
    ctrl_spectreasy <- unmix_controls(
        scc_dir = "scc",
        control_file = "fcs_mapping.csv",
        output_dir = "spectreasy_outputs/unmix_controls_spectreasy",
        unmixing_method = "Spectreasy",
        af_n_bands = 100,
        manual_gating = TRUE,
        gating_file = file.path(getwd(), "ssc_gate_config.csv"),
        spectreasy_weight_quantile = 0.9
    )

    unmixed_spectreasy <- unmix_samples(
        sample_dir = "samples",
        M = ctrl_spectreasy$M,
        unmixing_method = "Spectreasy",
        spectral_variant_library = ctrl_spectreasy$spectral_variant_library,
        output_dir = file.path(
            "spectreasy_outputs",
            "unmix_samples_spectreasy",
            "unmixed_fcs"
        ),
        spectreasy_weight_quantile = 0.9
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
        unmixing_method = NULL,
        res_list = NULL,
        pd = NULL,
        max_events_per_sample = 1000,
        overview_files_per_page = 15,
        matrix_markers_per_page = 20,
        sample_nxn_rows_per_page = 10,
        sample_nxn_max_points = max_events_per_sample,
        sample_nxn_transform = c("none", "asinh"),
        sample_nxn_asinh_cofactor = 150,
        sample_nxn_axis_limit = NULL,
        nxn_all_samples = FALSE,
        qc_plot_dir = NULL,
        save_qc_pngs = FALSE,
        qc_metrics_dir = NULL
    )

    # -------------------------------------------------------------------------
    # Interactive tools
    # -------------------------------------------------------------------------

    spectreasy_gui()

    adjust_matrix(
        matrix_dir = NULL,
        samples_dir = NULL,
        port = 8000,
        open_browser = TRUE,
        dev_mode = FALSE
    )

    build_panel(
        port = 8000,
        open_browser = TRUE,
        dev_mode = FALSE
    )

    gate_controls(
        scc_dir = "scc",
        control_file = "fcs_mapping.csv",
        gate_file = "ssc_gate_config.csv",
        port = 8000,
        open_browser = TRUE,
        dev_mode = FALSE
    )

    # -------------------------------------------------------------------------
    # Autofluorescence profile library workflow
    # -------------------------------------------------------------------------

    profile <- extract_af_profile(
        fcs_file = "scc/Unstained.fcs",
        af_n_bands = 100,
        af_max_cells = 50000,
        af_min_cluster_events = 20,
        af_min_cluster_proportion = 0.005,
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
        af_n_bands = 100,
        af_max_cells = 50000,
        af_min_cluster_events = 20,
        af_min_cluster_proportion = 0.005,
        seed = NULL,
        default_sample_type = "beads",
        cytometer = "auto",
        unmixing_method = "Spectreasy",
        scc_background_method = c("scatter_knn", "none"),
        scc_background_k = 2L,
        autospectral_n_candidates = 1000L,
        autospectral_n_spectral = 200L,
        autospectral_min_events = 10L,
        refine = FALSE,
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
        n_threads = 1L,
        spectral_variant_library = NULL,
        spectral_variant_top_k = 3L,
        spectral_variant_min_abundance = 1,
        spectral_variant_positive_fraction = 0.02,
        spectral_variant_min_improvement = 0.01
        # Spectreasy only:
        # , spectreasy_weight_quantile = 0.9
    )

    W <- derive_unmixing_matrix(
        M = M,
        method = "OLS",
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
        theme_custom = NULL,
        annotate_peaks = "auto",
        peak_label_max = 40L,
        peak_label_size = 2.6
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

    nps <- calculate_nps(
        data = as.data.frame(unmixed),
        markers = NULL
    )

    nps_plot <- plot_nps(
        nps_results = nps,
        output_file = NULL,
        width = 200
    )


    spectra <- get_control_spectra(
        flow_frame = flow_frame,
        control_file = "fcs_mapping.csv",
        control_dir = "scc",
        cytometer = "auto"
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

    cytometers <- supported_cytometers(
        include_auto = FALSE
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
