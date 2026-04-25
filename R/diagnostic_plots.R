.extract_top_detector_residuals <- function(res_list, top_n = 50) {
    residuals <- res_list$residuals
    if (is.null(residuals) || nrow(residuals) == 0) {
        warning("No residuals provided to plot_detector_residuals. Skipping plot.")
        return(NULL)
    }

    residual_score <- sqrt(rowMeans(residuals^2, na.rm = TRUE))
    residual_score[!is.finite(residual_score)] <- -Inf
    idx <- order(residual_score, decreasing = TRUE)[seq_len(min(top_n, nrow(res_list$data)))]
    residuals[idx, , drop = FALSE]
}

.resolve_detector_residual_labels <- function(det_names, pd = NULL) {
    if (!is.null(pd)) {
        det_info <- get_sorted_detectors(pd)
        common <- intersect(det_info$names, det_names)
        return(list(levels_sorted = common, labels_sorted = det_info$labels[match(common, det_info$names)]))
    }

    nums <- as.numeric(gsub("[^0-9]", "", det_names))
    levels_sorted <- det_names[order(nums)]
    list(levels_sorted = levels_sorted, labels_sorted = levels_sorted)
}

.prepare_detector_residual_long_df <- function(R_sub, levels_sorted, labels_sorted) {
    long <- as.data.frame(R_sub)
    long$Cell <- seq_len(nrow(R_sub))
    long <- tidyr::pivot_longer(long, cols = -Cell, names_to = "Detector", values_to = "Residual")
    long$Detector <- factor(long$Detector, levels = levels_sorted, labels = labels_sorted)
    long
}

.prepare_detector_overlay_df <- function(M, levels_sorted, labels_sorted, max_res) {
    M_overlay <- as.data.frame(M[, levels_sorted, drop = FALSE])
    M_overlay$Fluorophore <- rownames(M)
    M_long <- tidyr::pivot_longer(M_overlay, cols = -Fluorophore, names_to = "Detector", values_to = "Signature")
    M_long$Signature <- M_long$Signature * max_res
    M_long$Detector <- factor(M_long$Detector, levels = levels_sorted, labels = labels_sorted)
    M_long
}

.compute_detector_residual_breaks <- function(long, M_long) {
    y_all <- c(long$Residual, M_long$Signature)
    y_all <- y_all[is.finite(y_all)]
    if (length(y_all) == 0) return(NULL)

    max_abs <- max(abs(y_all), na.rm = TRUE)
    if (!is.finite(max_abs) || max_abs <= 0) return(NULL)

    max_pow <- max(1, ceiling(log10(max_abs)))
    pos_breaks <- 10^seq_len(max_pow)
    pos_breaks <- pos_breaks[pos_breaks <= (max_abs * 1.05)]
    if (length(pos_breaks) == 0) pos_breaks <- 10

    has_negative <- any(y_all < 0, na.rm = TRUE)
    y_breaks <- if (has_negative) c(-rev(pos_breaks), 0, pos_breaks) else c(0, pos_breaks)
    sort(unique(y_breaks))
}

.build_detector_residual_plot <- function(long, M_long, top_n = 50, y_breaks = NULL) {
    median_df <- long |>
        dplyr::group_by(Detector) |>
        dplyr::summarise(MedianResidual = stats::median(Residual, na.rm = TRUE), .groups = "drop")

    ggplot2::ggplot() +
        ggplot2::geom_line(
            data = M_long,
            ggplot2::aes(Detector, Signature, group = Fluorophore, color = Fluorophore),
            alpha = 0.5,
            linewidth = 0.5
        ) +
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
        ggplot2::labs(
            title = paste("Residual Contributions for Top", top_n, "Highest-Residual Cells"),
            subtitle = "Good: median residual (black) stays near 0 across detectors. Bad: consistent detector-specific shifts (positive/negative) indicate missing signatures, matrix mismatch, or calibration drift.",
            x = "Detector",
            y = "Residual Value / Scaled Signature"
        ) +
        ggplot2::scale_y_continuous(
            trans = scales::pseudo_log_trans(base = 10, sigma = 5),
            breaks = y_breaks,
            labels = scales::label_number(accuracy = 1, big.mark = ",", trim = TRUE)
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6))
}

#' Plot Detector-Level Residuals
#' 
#' Identifies which detectors contribute most to the unmixing mismatch for the
#' highest-residual cells. Overlays the reference signatures to help identify
#' the source of the mismatch.
#' 
#' @param res_list List returned by calc_residuals with return_residuals=TRUE
#' @param M Reference matrix
#' @param top_n Number of top high-error cells to analyze
#' @param output_file Optional path to save the plot. Set `NULL` to return the plot without writing a file.
#' @param width Plot width.
#' @param height Plot height.
#' @param pd Optional pData for descriptive labels
#' @return A `ggplot` object or `NULL` when residuals are unavailable.
#' @examples
#' M_demo <- rbind(
#'   FITC = c(1.00, 0.20, 0.05),
#'   PE = c(0.10, 1.00, 0.20),
#'   APC = c(0.05, 0.15, 1.00)
#' )
#' colnames(M_demo) <- c("B2-A", "YG1-A", "R1-A")
#'
#' marker_signal <- matrix(rexp(250 * nrow(M_demo), rate = 8), ncol = nrow(M_demo))
#' colnames(marker_signal) <- rownames(M_demo)
#' marker_signal[, "FITC"] <- rexp(250, rate = 0.6) + 2
#' raw_signal <- marker_signal %*% M_demo +
#'   matrix(rnorm(250 * ncol(M_demo), sd = 0.03), ncol = ncol(M_demo))
#' exprs_mat <- cbind(
#'   raw_signal,
#'   Time = seq_len(250),
#'   "FSC-A" = rnorm(250, mean = 90000, sd = 7000),
#'   "SSC-A" = rnorm(250, mean = 45000, sd = 5000)
#' )
#' colnames(exprs_mat)[seq_len(ncol(M_demo))] <- colnames(M_demo)
#' ff <- flowCore::flowFrame(exprs_mat)
#' res <- calc_residuals(ff, M_demo, return_residuals = TRUE)
#'
#' p <- plot_detector_residuals(res, M_demo, top_n = 20, output_file = NULL)
#' print(p)
#' @export
plot_detector_residuals <- function(res_list, M, top_n = 50, output_file = NULL, width = 250, height = 120, pd = NULL) {
    M <- .as_reference_matrix(M, "M")
    R_sub <- .extract_top_detector_residuals(res_list, top_n = top_n)
    if (is.null(R_sub)) {
        return(NULL)
    }

    label_info <- .resolve_detector_residual_labels(colnames(M), pd = pd)
    long <- .prepare_detector_residual_long_df(R_sub, label_info$levels_sorted, label_info$labels_sorted)
    max_res <- max(abs(long$Residual), na.rm = TRUE)
    M_long <- .prepare_detector_overlay_df(M, label_info$levels_sorted, label_info$labels_sorted, max_res = max_res)
    y_breaks <- .compute_detector_residual_breaks(long, M_long)
    p <- .build_detector_residual_plot(long, M_long, top_n = top_n, y_breaks = y_breaks)

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
#' demo_df <- data.frame(
#'   File = rep(c("sample_a", "sample_b"), each = 60),
#'   FITC = c(rnorm(60, 0.2, 0.2), rnorm(60, 0.1, 0.15)),
#'   PE = c(rnorm(60, 0.3, 0.25), rnorm(60, 0.2, 0.2))
#' )
#' nps <- calculate_nps(demo_df)
#' p <- plot_nps(nps, output_file = NULL)
#' print(p)
#' @export
plot_nps <- function(nps_results, output_file = NULL, width = 200) {
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
