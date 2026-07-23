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
        output_dir = "spectreasy_outputs",
        output_collision = c("version", "overwrite", "error"),
        unmixing_method = "AutoSpectral",
        unmix_scatter_panel_size_mm = 30,
        seed = NULL,
        af_profile = NULL,
        af_n_bands = 100,
        gating_mode = "interactive",
        manual_gate_file = "ssc_gate_config.csv",
        rwls_max_iter = 1L,
        n_threads = 1L,
        save_qc_png = FALSE,
        save_report = TRUE,
        save_ai_qc = TRUE,
        ai_qc_detail = "standard",
        ai_qc_privacy = "standard",
        ai_qc_reference = "auto",
        report_format = "html",
        scc_background_method = c("scatter_knn", "none"),
        scc_background_k = 2L,
        spectral_variant_som_nodes = 16L,
        spectral_variant_top_k = 3L,
        spectral_variant_cosine_threshold = 0.98,
        spectral_variant_max_variants = 8L,
        spectral_variant_min_events = 50L,
        autospectral_n_candidates = 1000L,
        autospectral_n_spectral = 200L,
        autospectral_min_events = 10L,
        autospectral_refine = FALSE,
        project_path = getwd()
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
        output_file = "spectreasy_outputs/unmix_controls/qc_controls_report.html",
        control_file = "fcs_mapping.csv",
        cytometer = "auto",
        unmixing_method = "WLS",
        qc_plot_dir = file.path("spectreasy_outputs", "scc_report_plots"),
        save_qc_pngs = FALSE,
        qc_metrics_dir = NULL,
        unmixed_list = NULL,
        qc_summary = NULL,
        report_plot_dir = NULL,
        pd = NULL,
        af_bank_info = NULL,
        cleanup_report_plot_dir = FALSE,
        unmix_scatter_max_points = 1000,
        unmix_scatter_axis_limit = NULL,
        seed = NULL,
        report_format = "html",
        overwrite = "version",
        report_plots = list(),
        report_run_settings = list(),
        report_artifact_paths = list(),
        project_path = getwd()
        # Additional build_reference_matrix() arguments can be passed here.
    )

    unmixed <- unmix_samples(
        sample_dir = "samples",
        samples_dir = NULL,
        M = NULL,
        unmixing_matrix_file = file.path(
            "spectreasy_outputs",
            "unmix_controls",
            "scc_reference_matrix.csv"
        ),
        detector_noise_file = NULL,
        control_file = NULL,
        unmixing_method = "AutoSpectral",
        rwls_max_iter = 1L,
        n_threads = 1L,
        spectral_variant_library = NULL,
        spectral_variant_library_file = NULL,
        spectral_variant_top_k = 3L,
        spectral_variant_min_abundance = 1,
        spectral_variant_positive_fraction = 0.02,
        spectral_variant_min_improvement = 0.01,
        estimate_af = FALSE,
        output_dir = "spectreasy_outputs",
        write_fcs = TRUE,
        save_report = TRUE,
        save_ai_qc = TRUE,
        ai_qc_detail = "standard",
        ai_qc_privacy = "standard",
        ai_qc_reference = "auto",
        report_format = "html",
        report_per_sample = FALSE,
        save_qc_plots = FALSE,
        qc_plot_dir = NULL,
        plot_n_events = 10000L,
        chunk_size = 50000L,
        seed = NULL,
        return_type = c("list", "flowSet", "SingleCellExperiment"),
        verbose = TRUE,
        project_path = getwd()
    )

    qc_samples_report <- qc_samples(
        results = unmixed,
        M = NULL,
        unmixing_matrix_file = file.path(
            "spectreasy_outputs",
            "unmix_controls",
            "scc_reference_matrix.csv"
        ),
        output_file = "spectreasy_outputs/unmix_samples/qc_samples_report.html",
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
        qc_metrics_dir = NULL,
        report_format = "html",
        report_per_sample = FALSE,
        overwrite = "version",
        report_run_settings = list(),
        report_artifact_paths = list(),
        project_path = getwd()
    )

    # Lower-level reusable HTML report workflow. Most users can call
    # qc_controls(..., report_format = "html") or
    # qc_samples(..., report_format = "html") directly instead.
    control_report_data <- collect_control_report_data(
        M = ctrl$M,
        scc_dir = "scc",
        control_file = "fcs_mapping.csv",
        cytometer = "auto",
        unmixing_method = "AutoSpectral",
        unmixed_list = ctrl$unmixed_list,
        qc_summary = NULL,
        report_plot_dir = NULL,
        pd = NULL,
        af_bank_info = NULL,
        matrix_source = ctrl$reference_matrix_file,
        artifact_paths = list(),
        plots = list(),
        warnings = character(),
        run_settings = list(),
        project_path = getwd(),
        plot_dir = NULL,
        qc_metrics_dir = NULL
    )

    sample_report_data <- collect_sample_report_data(
        results = unmixed,
        M = ctrl$M,
        unmixing_method = "AutoSpectral",
        res_list = NULL,
        pd = NULL,
        matrix_source = ctrl$reference_matrix_file,
        detector_noise_file = ctrl$detector_noise_file,
        spectral_variant_library_file = NULL,
        output_dir = file.path(
            "spectreasy_outputs",
            "unmix_samples",
            "unmixed_fcs"
        ),
        qc_metrics_dir = NULL,
        artifact_paths = list(),
        warnings = character(),
        run_settings = list(),
        project_path = getwd(),
        max_events_per_sample = 1000,
        overview_files_per_page = 15,
        matrix_markers_per_page = 20,
        sample_nxn_rows_per_page = 10,
        sample_nxn_max_points = max_events_per_sample,
        sample_nxn_transform = "none",
        sample_nxn_asinh_cofactor = 150,
        sample_nxn_axis_limit = NULL,
        nxn_all_samples = FALSE,
        report_per_sample = FALSE,
        plot_dir = NULL
    )

    html_report <- render_qc_html_report(
        report_data = sample_report_data,
        output_file = "spectreasy_outputs/unmix_samples/qc_samples_report.html",
        overwrite = "version"
    )

    control_metric_details <- calculate_control_qc_metrics(
        positive = matrix(numeric(), nrow = 0),
        negative = matrix(numeric(), nrow = 0),
        expected_peak = NULL,
        detector_ranges = NULL,
        final_reference = NULL,
        off_target_references = NULL,
        minimum_signal = NULL
    )

    af_metric_details <- calculate_af_bank_qc_metrics(
        af_bank = matrix(numeric(), nrow = 0),
        assignments = NULL,
        reconstruction_error = NULL,
        requested_bands = NULL
    )

    # -------------------------------------------------------------------------
    # Post-unmixing population analysis
    # -------------------------------------------------------------------------

    available_analysis_methods <- analysis_methods(
        refresh = FALSE
    )

    analysis_palettes <- population_analysis_palettes()

    population_analysis <- analyze_population(
        project_path = getwd(),
        file = "samples/sample.fcs",
        markers = c("CD3-A", "CD19-A", "CD56-A"),
        population_id = "root",
        clustering = "flowsom",
        reduction = "umap",
        trajectory = NULL,
        max_events = 20000L,
        cofactor = 150,
        seed = 20260723L,
        root_event_id = NULL,
        cluster_settings = list(clusters = 10L),
        reduction_settings = list(neighbors = 15L, min_dist = 0.01),
        trajectory_settings = list(),
        advanced_settings = list(),
        neighbors = 15L,
        clusters = 10L,
        perplexity = 30
    )

    saved_population_analysis <- load_population_analysis(
        project_path = getwd(),
        analysis_id = population_analysis$metadata$analysis_id,
        identity_id = NULL
    )

    identity_templates <- population_identity_templates(
        markers = c("CD3-A", "CD4-A", "CD8-A", "CD19-A", "CD56-A", "CD14-A")
    )

    population_annotation <- annotate_population(
        analysis = population_analysis,
        project_path = getwd(),
        signatures = identity_templates,
        min_score = 0.55,
        min_margin = 0.08,
        evidence_sensitivity = 1
    )

    population_plot <- plot_population_analysis(
        analysis = population_analysis,
        x = 1L,
        y = 2L,
        z = NULL,
        color_by = "CD3",
        palette = "viridis",
        point_size = 0.7,
        alpha = 0.82,
        title = NULL,
        background = "#f8f7f3"
    )

    population_export <- export_population_analysis(
        analysis = population_analysis,
        file = "spectreasy_outputs/analysis/plot.svg",
        format = "svg",
        width = 6,
        height = 6,
        dpi = 300
    )

    # -------------------------------------------------------------------------
    # Interactive tools
    # -------------------------------------------------------------------------

    spectreasy_gui(
        port = 8000
    )

    adjust_matrix(
        matrix_dir = NULL,
        samples_dir = NULL,
        port = 8000,
        open_browser = TRUE,
        dev_mode = FALSE,
        dev_frontend_port = NULL,
        unmixing_method = "AutoSpectral"
    )

    build_panel(
        port = 8000,
        open_browser = TRUE,
        dev_mode = FALSE,
        dev_frontend_port = NULL
    )

    gate_controls(
        scc_dir = "scc",
        control_file = "fcs_mapping.csv",
        gate_file = "ssc_gate_config.csv",
        port = 8000,
        open_browser = TRUE,
        dev_mode = FALSE,
        dev_frontend_port = NULL
    )

    # -------------------------------------------------------------------------
    # Population analysis workspace and gated exports
    # -------------------------------------------------------------------------

    workspace <- analysis_workspace(
        project_path = "."
    )

    save_analysis_workspace(
        project_path = ".",
        workspace = workspace
    )

    exported_gates <- export_population_gates(
        project_path = ".",
        path = "spectreasy_outputs/analysis/population-gates.csv"
    )

    imported_gates <- import_population_gates(
        project_path = ".",
        path = exported_gates[["path"]]
    )

    add_population_gate(
        project_path = ".",
        name = "CD3 positive",
        parent_id = "root",
        type = "rectangle",
        x = "CD3-A",
        y = "CD4-A",
        geometry = list(x_min = 0, x_max = 1, y_min = 0, y_max = 1),
        role = NULL,
        source_file = NULL,
        id = NULL
    )

    update_population_gate(
        project_path = ".",
        population_id = "population-1"
    )

    delete_population_gate(
        project_path = ".",
        population_id = "population-1"
    )

    statistics <- population_statistics(
        project_path = ".",
        file = "unmixed_fcs/sample.fcs",
        population_id = "root",
        markers = NULL
    )

    index <- staining_index(
        project_path = ".",
        file = "unmixed_fcs/sample.fcs",
        marker = "CD3-A"
    )

    exported_population <- export_gated_population(
        project_path = ".",
        files = "unmixed_fcs/sample.fcs",
        population_id = "root",
        format = "fcs",
        max_events = 0,
        seed = 1L,
        output_folder = "spectreasy_outputs/analysis/exports"
    )

    markers <- find_population_markers(
        analysis = analysis,
        project_path = ".",
        top_n = 10L,
        minimum_auc = 0.6
    )

    runtime <- analysis_runtime_status()

    install_analysis_runtime(
        python = NULL,
        path = NULL,
        rebuild = FALSE
    )

    install_analysis_dependencies(
        include_python = TRUE,
        python = Sys.which("python3"),
        ask = FALSE,
        reinstall = FALSE
    )

    # -------------------------------------------------------------------------
    # Autofluorescence profile library workflow
    # -------------------------------------------------------------------------

    profile <- extract_af_profile(
        fcs_file = "scc/Unstained.fcs",
        af_n_bands = 100,
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

    rename_af_profile(
        name = "my_af_profile",
        new_name = "renamed_af_profile"
    )

    delete_af_profile(
        name = "renamed_af_profile"
    )

    # -------------------------------------------------------------------------
    # Lower-level reference matrix and unmixing helpers
    # -------------------------------------------------------------------------

    M <- build_reference_matrix(
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
        unmixing_method = "AutoSpectral",
        scc_background_method = c("scatter_knn", "none"),
        scc_background_k = 2L,
        autospectral_n_candidates = 1000L,
        autospectral_n_spectral = 200L,
        autospectral_min_events = 10L,
        refine = FALSE,
        n_threads = 1L,
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

    similarity_matrix <- calculate_similarity_matrix(
        M = M
    )

    similarity_plot <- plot_similarity_matrix(
        similarity_matrix = similarity_matrix,
        output_file = NULL,
        width = 180,
        height = 160
    )

    sample_rms_plot <- plot_sample_rms_residuals(
        results = unmixed,
        M = M,
        output_file = NULL,
        width = 225,
        height = 125,
        unmixing_method = NULL
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
