.calculate_panel_complexity <- function(spectra) {
    spectra <- .as_reference_matrix(spectra, "spectra")
    if (nrow(spectra) < 2) return(1)
    sv <- svd(spectra)$d
    sv <- sv[is.finite(sv) & sv > 0]
    if (length(sv) == 0) return(NA_real_)
    round(max(sv) / min(sv), 2)
}

.matrix_to_named_rows <- function(M, row_name = "fluorophore") {
    if (length(M) == 0 || nrow(M) == 0) return(data.frame())
    out <- as.data.frame(M, check.names = FALSE)
    out[[row_name]] <- rownames(M)
    out[, c(row_name, setdiff(colnames(out), row_name)), drop = FALSE]
}

.build_spectral_panel_data <- function(fluorophores,
                                       cytometer = "aurora",
                                       configuration = NULL,
                                       strict = TRUE) {
    spectra <- .spectral_panel_configuration_spectra(
        cytometer = cytometer,
        configuration = configuration,
        fluorophores = fluorophores,
        strict = strict
    )
    peak_detectors <- apply(spectra, 1, function(x) colnames(spectra)[which.max(x)])
    list(
        spectra = spectra,
        similarity_matrix = calculate_similarity_matrix(spectra),
        complexity_index = .calculate_panel_complexity(spectra),
        peak_detectors = peak_detectors
    )
}

.spectral_panel_export_table <- function(fluorophores, markers = NULL) {
    fluorophores <- trimws(as.character(fluorophores))
    fluorophores <- fluorophores[nzchar(fluorophores)]
    markers <- if (is.null(markers)) rep("", length(fluorophores)) else trimws(as.character(markers))
    length(markers) <- length(fluorophores)
    markers[is.na(markers)] <- ""
    data.frame(
        Marker = markers,
        Fluorophore = fluorophores,
        stringsAsFactors = FALSE
    )
}

.spectral_panel_label_rows <- function(fluorophores, markers = NULL) {
    export_df <- .spectral_panel_export_table(fluorophores, markers)
    labels <- ifelse(
        nzchar(export_df$Marker),
        paste0(export_df$Marker, " / ", export_df$Fluorophore),
        export_df$Fluorophore
    )
    make.unique(labels, sep = " ")
}

.spectral_panel_signature_bins <- function(signature,
                                           detector_labels,
                                           y_power = 1.5,
                                           max_log_intensity = 6,
                                           n_bins = 37L) {
    signature <- as.numeric(signature)
    signature[!is.finite(signature)] <- 0
    signature <- pmax(0, pmin(1, signature))

    offsets <- seq(-0.42, 0.42, length.out = n_bins)
    center_weight <- 1 - abs(seq(-1, 1, length.out = n_bins))^1.6
    center_weight <- pmax(0.08, center_weight)

    rows <- lapply(seq_along(signature), function(i) {
        center_log <- 0.35 + (signature[i]^0.72) * (max_log_intensity - 0.35)
        y_orig <- pmax(0.05, pmin(max_log_intensity, center_log + offsets))
        data.frame(
            ch_idx = i,
            Detector = detector_labels[i],
            y_orig = y_orig,
            y = y_orig^y_power,
            fill = pmin(1, center_weight * pmax(0.12, sqrt(signature[i])) * 1.35),
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    out$Detector <- factor(out$Detector, levels = detector_labels)
    out
}

.plot_spectral_panel_signature <- function(signature,
                                           detector_labels,
                                           title,
                                           y_power = 1.5,
                                           max_log_intensity = 6) {
    dt_c <- .spectral_panel_signature_bins(
        signature = signature,
        detector_labels = detector_labels,
        y_power = y_power,
        max_log_intensity = max_log_intensity
    )
    y_breaks_orig <- 0:max_log_intensity
    y_breaks_trans <- y_breaks_orig^y_power
    y_labels <- vapply(y_breaks_orig, function(x) paste0("10^", x), character(1))

    ggplot2::ggplot(dt_c, ggplot2::aes(ch_idx, y, fill = fill)) +
        ggplot2::geom_tile(width = 0.76, height = 0.055) +
        ggplot2::scale_fill_gradientn(
            colors = c("#0000FF", "#00FFFF", "#00FF00", "#FFFF00", "#FF0000"),
            limits = c(0, 1),
            oob = scales::squish,
            guide = "none"
        ) +
        ggplot2::scale_x_continuous(
            breaks = seq_along(detector_labels),
            labels = detector_labels
        ) +
        ggplot2::scale_y_continuous(
            limits = c(-0.1, (max_log_intensity + 0.35)^y_power),
            breaks = y_breaks_trans,
            labels = y_labels
        ) +
        ggplot2::coord_cartesian(expand = FALSE) +
        ggplot2::labs(title = title, x = NULL, y = "Intensity") +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5),
            axis.text.y = ggplot2::element_text(size = 8),
            axis.title.y = ggplot2::element_text(size = 9, face = "bold"),
            panel.grid.major = ggplot2::element_line(color = "#eeeeee", linewidth = 0.25),
            panel.grid.minor = ggplot2::element_line(color = "#f4f4f4", linewidth = 0.18),
            panel.background = ggplot2::element_rect(fill = "white", color = NA),
            plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
            plot.margin = ggplot2::margin(4, 8, 4, 4)
        )
}

.plot_spectral_panel_bands <- function(spectra, cytometer = "aurora", markers = NULL) {
    spectra <- .as_reference_matrix(spectra, "spectra")
    labels <- .spectral_panel_label_rows(rownames(spectra), markers)
    detector_info <- .spectral_detector_metadata(cytometer, colnames(spectra))
    detector_labels <- detector_info$label

    pages <- vector("list", nrow(spectra))
    for (i in seq_len(nrow(spectra))) {
        pages[[i]] <- .plot_spectral_panel_signature(
            signature = spectra[i, ],
            detector_labels = detector_labels,
            title = labels[i]
        )
    }
    pages
}

.draw_spectral_panel_signature_pages <- function(signature_pages, plots_per_page = 2L) {
    if (length(signature_pages) == 0) return(invisible(NULL))
    plots_per_page <- max(1L, as.integer(plots_per_page[1]))
    starts <- seq(1L, length(signature_pages), by = plots_per_page)
    for (start_idx in starts) {
        page_plots <- signature_pages[start_idx:min(start_idx + plots_per_page - 1L, length(signature_pages))]
        grid::grid.newpage()
        if (length(page_plots) == 1L) {
            grid::grid.draw(grid::editGrob(
                ggplot2::ggplotGrob(page_plots[[1]]),
                vp = grid::viewport(x = 0.5, y = 0.50, width = 0.94, height = 0.72)
            ))
        } else {
            grid::grid.draw(grid::editGrob(
                ggplot2::ggplotGrob(page_plots[[1]]),
                vp = grid::viewport(x = 0.5, y = 0.72, width = 0.94, height = 0.43)
            ))
            grid::grid.draw(grid::editGrob(
                ggplot2::ggplotGrob(page_plots[[2]]),
                vp = grid::viewport(x = 0.5, y = 0.27, width = 0.94, height = 0.43)
            ))
        }
    }
    invisible(NULL)
}

.draw_spectral_panel_similarity_pages <- function(similarity_matrix,
                                                  complexity_index,
                                                  cytometer_label = "",
                                                  configuration_label = "",
                                                  fluorophore_count = NULL,
                                                  max_markers_per_page = NULL) {
    if (is.null(similarity_matrix) || nrow(similarity_matrix) < 2L) return(invisible(NULL))
    n_markers <- nrow(similarity_matrix)
    if (is.null(max_markers_per_page)) {
        max_markers_per_page <- if (n_markers <= 15L) 15L else 8L
    }
    pages <- .build_qc_report_matrix_pages(
        similarity_matrix,
        plot_fun = plot_similarity_matrix,
        max_markers_per_page = max_markers_per_page,
        item_label = "Fluorophores"
    )
    if (length(pages) == 0) return(invisible(NULL))

    panel_label_parts <- c(
        cytometer_label,
        configuration_label,
        if (!is.null(fluorophore_count)) paste0(fluorophore_count, " fluorophore(s)") else NULL
    )
    panel_label <- paste(panel_label_parts[nzchar(panel_label_parts)], collapse = " | ")

    for (p in pages) {
        p <- p +
            ggplot2::theme(
                plot.title = ggplot2::element_blank(),
                plot.subtitle = ggplot2::element_blank(),
                axis.text.x = ggplot2::element_text(size = 8, angle = 45, hjust = 1),
                axis.text.y = ggplot2::element_text(size = 8),
                legend.position = "right",
                plot.margin = ggplot2::margin(8, 8, 8, 8)
            )
        grid::grid.newpage()
        grid::grid.text(
            "Fluorophore Spectral Similarity",
            x = 0.04,
            y = 0.965,
            just = c("left", "top"),
            gp = grid::gpar(fontsize = 18, fontface = "bold")
        )
        if (nzchar(panel_label)) {
            grid::grid.text(
                panel_label,
                x = 0.04,
                y = 0.925,
                just = c("left", "top"),
                gp = grid::gpar(fontsize = 10, col = "#475569")
            )
        }
        grid::grid.text(
            paste0("Complexity Index: ", format(round(complexity_index, 2), nsmall = 2)),
            x = 0.965,
            y = 0.955,
            just = c("right", "top"),
            gp = grid::gpar(fontsize = 12, fontface = "bold")
        )
        grid::grid.draw(grid::editGrob(
            ggplot2::ggplotGrob(p),
            vp = grid::viewport(x = 0.50, y = 0.46, width = 0.94, height = 0.78)
        ))
    }
    invisible(NULL)
}
