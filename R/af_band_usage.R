# AF-band assignment summaries used by sample workflow reporting.

.empty_af_band_usage <- function() {
    data.frame(
        sample = character(),
        af_band = character(),
        assigned_events = integer(),
        eligible_events = integer(),
        usage_fraction = numeric(),
        stringsAsFactors = FALSE
    )
}

.af_band_usage_rows <- function(sample, af_bands, assignment_counts, eligible_events) {
    if (!length(af_bands)) return(.empty_af_band_usage())
    assignment_counts <- as.integer(rep_len(assignment_counts, length(af_bands)))
    eligible_events <- as.integer(eligible_events)[1]
    data.frame(
        sample = rep(as.character(sample)[1], length(af_bands)),
        af_band = as.character(af_bands),
        assigned_events = assignment_counts,
        eligible_events = rep(eligible_events, length(af_bands)),
        usage_fraction = if (eligible_events > 0L) assignment_counts / eligible_events else rep(0, length(af_bands)),
        stringsAsFactors = FALSE
    )
}

.normalize_af_band_usage <- function(x) {
    required <- c("sample", "af_band", "assigned_events", "eligible_events", "usage_fraction")
    if (is.null(x) || !is.data.frame(x) || !all(required %in% colnames(x))) return(.empty_af_band_usage())
    out <- x[, required, drop = FALSE]
    out$sample <- as.character(out$sample)
    out$af_band <- as.character(out$af_band)
    out$assigned_events <- as.integer(out$assigned_events)
    out$eligible_events <- as.integer(out$eligible_events)
    out$usage_fraction <- as.numeric(out$usage_fraction)
    out$usage_fraction[!is.finite(out$usage_fraction)] <- 0
    out
}

.build_af_band_usage_pages <- function(usage, max_samples_per_page = 30L) {
    usage <- .normalize_af_band_usage(usage)
    if (!nrow(usage)) return(list())
    max_samples_per_page <- max(1L, as.integer(max_samples_per_page)[1])
    sample_order <- unique(usage$sample)
    band_order <- unique(usage$af_band)
    max_usage <- max(usage$usage_fraction, na.rm = TRUE)
    if (!is.finite(max_usage) || max_usage <= 0) max_usage <- 1
    groups <- split(sample_order, ceiling(seq_along(sample_order) / max_samples_per_page))
    show_text <- length(sample_order) <= 20L && length(band_order) <= 20L
    lapply(groups, function(samples) {
        page <- usage[usage$sample %in% samples, , drop = FALSE]
        page$sample <- factor(page$sample, levels = rev(samples))
        page$af_band <- factor(page$af_band, levels = band_order)
        p <- ggplot2::ggplot(page, ggplot2::aes(x = .data$af_band, y = .data$sample, fill = .data$usage_fraction)) +
            ggplot2::geom_tile(colour = "white", linewidth = 0.35) +
            ggplot2::scale_fill_viridis_c(
                limits = c(0, max_usage), option = "D", name = "Fraction of events"
            ) +
            ggplot2::labs(x = "AF bands", y = "Samples") +
            ggplot2::theme_minimal(base_size = 11) +
            ggplot2::theme(
                panel.grid = ggplot2::element_blank(),
                axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
            )
        if (show_text) {
            p <- p + ggplot2::geom_text(
                ggplot2::aes(label = sprintf("%.1f%%", 100 * .data$usage_fraction)),
                size = 3
            )
        }
        attr(p, "spectreasy_af_usage_limits") <- c(0, max_usage)
        p
    })
}

.save_af_band_usage_plots <- function(usage, plot_dir) {
    pages <- .build_af_band_usage_pages(usage)
    if (!length(pages)) return(character())
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    vapply(seq_along(pages), function(index) {
        .report_plot_file(
            pages[[index]],
            file.path(plot_dir, sprintf("af_band_usage_%02d.png", index)),
            width = 8.5,
            height = max(4, min(10, length(unique(usage$sample)) * 0.30 + 2.2)),
            dpi = 220
        )
    }, character(1))
}
