# SCC positive-event spectral variation summaries and plots.

.collect_scc_variation_data <- function(M,
                                        lower_quantile = 0.10,
                                        upper_quantile = 0.90,
                                        min_events = 20L) {
    M <- .as_reference_matrix(M, "M")
    lower_quantile <- as.numeric(lower_quantile)[1]
    upper_quantile <- as.numeric(upper_quantile)[1]
    min_events <- as.integer(min_events)[1]
    if (!is.finite(lower_quantile) || !is.finite(upper_quantile) ||
        lower_quantile < 0 || upper_quantile > 1 || lower_quantile >= upper_quantile) {
        stop("lower_quantile and upper_quantile must satisfy 0 <= lower < upper <= 1.", call. = FALSE)
    }
    if (!is.finite(min_events) || min_events < 1L) {
        stop("min_events must be a positive integer.", call. = FALSE)
    }

    detector_order <- colnames(M)[.residual_detector_channel_order(colnames(M))]
    fluorophores <- rownames(M)[!grepl("^AF($|_)", rownames(M), ignore.case = TRUE)]
    positive_events <- attr(M, "scc_positive_events")
    data <- list()
    metadata <- lapply(fluorophores, function(fluorophore) {
        raw <- if (is.list(positive_events)) positive_events[[fluorophore]] else NULL
        raw_count <- if (is.matrix(raw) || is.data.frame(raw)) nrow(raw) else 0L
        status <- "missing_positive_events"
        used_count <- 0L

        if (!is.null(raw)) {
            raw <- as.matrix(raw)
            if (is.null(colnames(raw)) || length(setdiff(detector_order, colnames(raw))) > 0L) {
                status <- "detector_mismatch"
            } else {
                raw <- raw[, detector_order, drop = FALSE]
                storage.mode(raw) <- "double"
                raw[!is.finite(raw)] <- 0
                raw[raw < 0] <- 0
                maxima <- apply(raw, 1L, max)
                valid <- is.finite(maxima) & maxima > 0
                normalized <- raw[valid, , drop = FALSE]
                if (nrow(normalized)) normalized <- normalized / maxima[valid]
                used_count <- nrow(normalized)
                if (used_count == 0L) {
                    status <- "invalid_event_spectra"
                } else if (used_count < min_events) {
                    status <- "insufficient_events"
                } else {
                    reference <- as.numeric(M[fluorophore, detector_order, drop = TRUE])
                    reference[!is.finite(reference)] <- 0
                    reference[reference < 0] <- 0
                    reference_max <- max(reference)
                    if (!is.finite(reference_max) || reference_max <= 0) {
                        status <- "invalid_event_spectra"
                    } else {
                        status <- "available"
                        data[[fluorophore]] <<- data.frame(
                            detector = detector_order,
                            lower = apply(normalized, 2L, stats::quantile, probs = lower_quantile, names = FALSE, type = 7),
                            upper = apply(normalized, 2L, stats::quantile, probs = upper_quantile, names = FALSE, type = 7),
                            reference = reference / reference_max,
                            stringsAsFactors = FALSE
                        )
                    }
                }
            }
        }
        data.frame(
            fluorophore = fluorophore,
            event_count_raw = as.integer(raw_count),
            event_count_used = as.integer(used_count),
            status = status,
            stringsAsFactors = FALSE
        )
    })
    metadata <- if (length(metadata)) do.call(rbind, metadata) else data.frame(
        fluorophore = character(), event_count_raw = integer(),
        event_count_used = integer(), status = character(), stringsAsFactors = FALSE
    )
    list(data = data, metadata = metadata, detector_order = detector_order)
}

.plot_scc_variation_band <- function(variation_data, fluorophore, pd = NULL) {
    if (is.list(variation_data) && !is.data.frame(variation_data)) {
        variation_data <- variation_data$data[[fluorophore]]
    }
    if (is.null(variation_data) || !is.data.frame(variation_data) || !nrow(variation_data)) {
        stop("No SCC variation data are available for fluorophore: ", fluorophore, call. = FALSE)
    }
    variation_data$x <- seq_len(nrow(variation_data))
    ggplot2::ggplot(variation_data, ggplot2::aes(x = .data$x)) +
        ggplot2::geom_ribbon(
            ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
            fill = "#4f94c7", alpha = 0.30
        ) +
        ggplot2::geom_line(ggplot2::aes(y = .data$reference), colour = "black", linewidth = 1.2) +
        ggplot2::scale_x_continuous(
            breaks = variation_data$x,
            labels = variation_data$detector,
            expand = ggplot2::expansion(mult = c(0.01, 0.01))
        ) +
        ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
        ggplot2::labs(title = fluorophore, x = "Detector", y = "Normalized intensity") +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            plot.title = ggplot2::element_text(face = "bold", size = 12)
        )
}

.save_scc_variation_plots <- function(variation, plot_dir, pd = NULL) {
    if (!length(variation$data)) return(character())
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    fluorophores <- names(variation$data)
    paths <- vapply(fluorophores, function(fluorophore) {
        slug <- .report_filename_slug(fluorophore, fallback = "fluorophore")
        .report_plot_file(
            .plot_scc_variation_band(variation, fluorophore, pd = pd),
            file.path(plot_dir, paste0(slug, ".png")),
            width = 7.2, height = 4.3, dpi = 220
        )
    }, character(1))
    stats::setNames(paths, fluorophores)
}
