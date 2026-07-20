.write_spectral_panel_overview_pdf <- function(cytometer = "aurora",
                                               configuration = NULL,
                                               fluorophores,
                                               markers = NULL,
                                               output_file) {
    if (missing(output_file) || is.null(output_file) || !nzchar(trimws(as.character(output_file)[1]))) {
        stop("output_file is required.", call. = FALSE)
    }

    requested_fluorophores <- unique(trimws(as.character(fluorophores)))
    requested_fluorophores <- requested_fluorophores[nzchar(requested_fluorophores)]
    requested_markers <- if (is.null(markers)) rep("", length(requested_fluorophores)) else trimws(as.character(markers))
    length(requested_markers) <- length(requested_fluorophores)
    requested_markers[is.na(requested_markers)] <- ""
    marker_map <- stats::setNames(requested_markers, requested_fluorophores)

    panel <- .build_spectral_panel_data(
        fluorophores = requested_fluorophores,
        cytometer = cytometer,
        configuration = configuration,
        strict = FALSE
    )
    spectra <- panel$spectra
    if (nrow(spectra) == 0) {
        stop("Select at least one fluorophore before exporting a panel overview.", call. = FALSE)
    }

    markers <- unname(marker_map[rownames(spectra)])
    markers[is.na(markers)] <- ""

    id <- .resolve_spectral_panel_cytometer(cytometer)
    config <- .resolve_spectral_panel_configuration(id, configuration)
    libs <- .spectral_panel_libraries()
    cyt_label <- libs$label[match(id, libs$id)]
    if (is.na(cyt_label) || !nzchar(cyt_label)) cyt_label <- id
    configs <- .spectral_panel_configurations(id)
    config_label <- configs$label[match(config, configs$id)]
    if (is.na(config_label) || !nzchar(config_label)) config_label <- config

    signature_pages <- .plot_spectral_panel_bands(spectra, cytometer = id, markers = markers)
    sim <- panel$similarity_matrix

    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    grDevices::pdf(output_file, width = 11, height = 8.5)
    on.exit(grDevices::dev.off(), add = TRUE)

    .draw_spectral_panel_similarity_pages(
        similarity_matrix = sim,
        complexity_index = panel$complexity_index,
        cytometer_label = cyt_label,
        configuration_label = config_label,
        fluorophore_count = nrow(spectra)
    )
    .draw_spectral_panel_signature_pages(signature_pages, plots_per_page = 2L)

    invisible(output_file)
}

.spectral_panel_payload <- function(cytometer = "aurora",
                                    fluorophores = character(),
                                    configuration = NULL) {
    id <- .resolve_spectral_panel_cytometer(cytometer)
    config <- .resolve_spectral_panel_configuration(id, configuration)
    full_library <- .load_spectral_library(id, renormalize = FALSE)
    config_detectors <- .spectral_panel_configuration_detectors(
        cytometer = id,
        configuration = config,
        detectors = colnames(full_library)
    )
    retained_signal <- apply(abs(full_library[, config_detectors, drop = FALSE]), 1, max, na.rm = TRUE)
    retained_signal[!is.finite(retained_signal)] <- 0
    available <- retained_signal >= 0.02
    library <- .normalize_spectral_rows(full_library[available, config_detectors, drop = FALSE])
    detector_info <- .spectral_detector_metadata(id, colnames(library))
    fluorophores <- unique(trimws(as.character(fluorophores)))
    fluorophores <- fluorophores[nzchar(fluorophores)]

    selected <- .spectral_panel_configuration_spectra(
        cytometer = id,
        configuration = config,
        fluorophores = fluorophores,
        strict = FALSE
    )
    selected_names <- rownames(selected)
    if (is.null(selected_names)) selected_names <- character()
    sim <- if (length(selected_names) > 0) calculate_similarity_matrix(selected) else matrix(numeric(0), nrow = 0, ncol = 0)
    complexity <- if (length(selected_names) > 0) .calculate_panel_complexity(selected) else NA_real_
    peaks <- if (length(selected_names) > 0) {
        apply(selected, 1, function(x) colnames(selected)[which.max(x)])
    } else {
        character()
    }

    fluor_table <- data.frame(
        fluorophore = rownames(library),
        peak_detector = apply(library, 1, function(x) colnames(library)[which.max(x)]),
        retained_signal = as.numeric(retained_signal[rownames(library)]),
        stringsAsFactors = FALSE
    )
    fluor_table$peak_laser <- detector_info$laser[match(fluor_table$peak_detector, detector_info$detector)]
    fluor_table$peak_color <- detector_info$color[match(fluor_table$peak_detector, detector_info$detector)]
    fluor_table <- fluor_table[order(fluor_table$peak_laser, fluor_table$fluorophore), , drop = FALSE]

    list(
        cytometer = id,
        configuration = config,
        libraries = .spectral_panel_libraries(),
        configurations = .spectral_panel_configurations(id),
        detectors = detector_info,
        fluorophores = fluor_table,
        selected = selected_names,
        spectra = .matrix_to_named_rows(selected),
        similarity = .matrix_to_named_rows(sim),
        complexity_index = complexity,
        peak_detectors = unname(peaks)
    )
}
