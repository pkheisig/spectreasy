#' Plot Detector-Level Residuals
#' 
#' Identifies which detectors contribute most to the unmixing mismatch for the
#' highest-residual cells. Overlays the reference signatures to help identify
#' the source of the mismatch.
#' 
#' @param res_list List returned by calc_residuals with return_residuals=TRUE
#' @param M Reference matrix
#' @param top_n Number of top high-error cells to analyze
#' @param output_file Path to save the plot
#' @param width Plot width.
#' @param height Plot height.
#' @param pd Optional pData for descriptive labels
#' @return A `ggplot` object or `NULL` when residuals are unavailable.
#' @examples
#' \dontrun{
#' p <- plot_detector_residuals(res_list[[1]], M, top_n = 50)
#' print(p)
#' }
#' @export
plot_detector_residuals <- function(res_list, M, top_n = 50, output_file = "detector_residuals.png", width = 250, height = 120, pd = NULL) {
    M <- .as_reference_matrix(M, "M")
    data <- res_list$data
    residuals <- res_list$residuals
    
    if (is.null(residuals) || nrow(residuals) == 0) {
        warning("No residuals provided to plot_detector_residuals. Skipping plot.")
        return(NULL)
    }

    # Identify highest-residual cells by detector residual energy
    residual_score <- sqrt(rowMeans(residuals^2, na.rm = TRUE))
    residual_score[!is.finite(residual_score)] <- -Inf
    idx <- order(residual_score, decreasing = TRUE)[1:min(top_n, nrow(data))]
    R_sub <- residuals[idx, , drop = FALSE]
    
    # Convert to long format for plotting
    long <- as.data.frame(R_sub)
    long$Cell <- seq_len(nrow(R_sub))
    long <- tidyr::pivot_longer(long, cols = -Cell, names_to = "Detector", values_to = "Residual")
    
    # Sort and label detectors
    det_names <- colnames(M)
    if (!is.null(pd)) {
        det_info <- get_sorted_detectors(pd)
        common <- intersect(det_info$names, det_names)
        levels_sorted <- common
        labels_sorted <- det_info$labels[match(common, det_info$names)]
    } else {
        nums <- as.numeric(gsub("[^0-9]", "", det_names))
        levels_sorted <- det_names[order(nums)]
        labels_sorted <- levels_sorted
    }
    
    long$Detector <- factor(long$Detector, levels = levels_sorted, labels = labels_sorted)
    
    # Prepare M for overlay
    # Scale M to the max absolute residual for visual comparison
    max_res <- max(abs(long$Residual), na.rm = TRUE)
    M_overlay <- as.data.frame(M[, levels_sorted, drop = FALSE])
    M_overlay$Fluorophore <- rownames(M)
    M_long <- tidyr::pivot_longer(M_overlay, cols = -Fluorophore, names_to = "Detector", values_to = "Signature")
    M_long$Signature <- M_long$Signature * max_res
    M_long$Detector <- factor(M_long$Detector, levels = levels_sorted, labels = labels_sorted)

    # Use signed pseudo-log scaling: expand around 0, compress large magnitudes.
    y_all <- c(long$Residual, M_long$Signature)
    y_all <- y_all[is.finite(y_all)]
    y_breaks <- NULL
    y_sigma <- 5
    if (length(y_all) > 0) {
        max_abs <- max(abs(y_all), na.rm = TRUE)
        if (is.finite(max_abs) && max_abs > 0) {
            max_pow <- max(1, ceiling(log10(max_abs)))
            pos_breaks <- 10^(1:max_pow)
            pos_breaks <- pos_breaks[pos_breaks <= (max_abs * 1.05)]
            if (length(pos_breaks) == 0) pos_breaks <- 10

            has_negative <- any(y_all < 0, na.rm = TRUE)
            if (has_negative) {
                neg_breaks <- -rev(pos_breaks)
                y_breaks <- c(neg_breaks, 0, pos_breaks)
            } else {
                y_breaks <- c(0, pos_breaks)
            }
            y_breaks <- sort(unique(y_breaks))
        }
    }
    
    median_df <- long |>
        dplyr::group_by(Detector) |>
        dplyr::summarise(MedianResidual = stats::median(Residual, na.rm = TRUE), .groups = "drop")

    p <- ggplot2::ggplot() +
        # Spectra overlay (background)
        ggplot2::geom_line(data = M_long, ggplot2::aes(Detector, Signature, group = Fluorophore, color = Fluorophore), 
                           alpha = 0.5, linewidth = 0.5) +
        # Median residual trend (foreground)
        ggplot2::geom_line(
            data = median_df,
            ggplot2::aes(Detector, MedianResidual, group = 1),
            color = "black",
            linewidth = 0.8
        ) +
        ggplot2::geom_point(
            data = median_df,
            ggplot2::aes(Detector, MedianResidual),
            color = "black",
            size = 1.1
        ) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        ggplot2::labs(title = paste("Residual Contributions for Top", top_n, "Highest-Residual Cells"),
                      subtitle = "Good: median residual (black) stays near 0 across detectors. Bad: consistent detector-specific shifts (positive/negative) indicate missing signatures, matrix mismatch, or calibration drift.",
                      x = "Detector", y = "Residual Value / Scaled Signature") +
        ggplot2::scale_y_continuous(
            trans = scales::pseudo_log_trans(base = 10, sigma = y_sigma),
            breaks = y_breaks,
            labels = scales::label_number(accuracy = 1, big.mark = ",", trim = TRUE)
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6))
    
    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, p, width = width, height = height, units = "mm", dpi = 300)
    }
    return(p)
}

#' Calculate Negative Population Spread (NPS)
#' 
#' Quantifies the spreading of negative populations after unmixing. 
#' High spread indicates poor unmixing or high spillover noise.
#' 
#' @param data Unmixed data frame
#' @param markers Vector of marker names to analyze
#' @return A data frame with NPS values (MAD) per marker per file
#' @examples
#' demo_df <- data.frame(
#'   File = rep(c("A", "B"), each = 50),
#'   FITC = rnorm(100, sd = 0.3),
#'   PE = rnorm(100, sd = 0.4)
#' )
#' nps <- calculate_nps(demo_df)
#' head(nps)
#' @export
calculate_nps <- function(data, markers = NULL) {
    if (is.null(markers)) {
        exclude <- .get_result_metadata_columns(colnames(data))
        markers <- setdiff(colnames(data), exclude)
        markers <- markers[!grepl("^AF($|_)", markers, ignore.case = TRUE)]
    }
    
    # For each marker, isolate the negative population (intensity < threshold)
    # We use a robust estimate (Median Absolute Deviation) of the spread
    nps_results <- data |> 
        dplyr::group_by(File) |> 
        dplyr::summarize(dplyr::across(dplyr::all_of(markers), function(x) {
            # Heuristic: negative population is around 0. 
            # We take values between -2SD and +2SD to avoid true positives
            # But simpler: just take the MAD of all values < quantile(0.2)
            neg_subset <- x[x < quantile(x, 0.2)]
            stats::mad(neg_subset, na.rm = TRUE)
        })) |> 
        tidyr::pivot_longer(cols = -File, names_to = "Marker", values_to = "NPS")
    
    return(nps_results)
}

#' Plot Negative Population Spread
#' @param nps_results Output of [calculate_nps()].
#' @param output_file Path to save the plot (set `NULL` to skip saving).
#' @param width Plot width.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' nps <- calculate_nps(results_df)
#' p <- plot_nps(nps, output_file = "nps_plot.png")
#' print(p)
#' }
#' @export
plot_nps <- function(nps_results, output_file = "nps_plot.png", width = 200) {
    p <- ggplot2::ggplot(nps_results, ggplot2::aes(Marker, NPS, fill = File)) +
        ggplot2::geom_bar(stat = "identity", position = "dodge") +
        ggplot2::labs(title = "Negative Population Spread (Unmixing Noise Floor)",
                      subtitle = "Good: low, similar MAD across files and markers. Bad: isolated high bars indicate broad negative-population spread, often from spillover into dim channels.",
                      y = "Spread (MAD)", x = "Unmixed Marker") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    
    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, p, width = width, height = 120, units = "mm", dpi = 300)
    }
    return(p)
}
