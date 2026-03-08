#' Generate a Full QC PDF Report
#'
#' Creates a multi-page report summarizing unmixing quality, including spectra,
#' detector residuals, spread matrix, and NPS.
#'
#' @param results_df Combined unmixed data frame (typically `rbind` of sample results).
#' @param M Reference matrix used for unmixing.
#' @param output_file Output PDF file path.
#' @param res_list Optional residual object/list from `calc_residuals(..., return_residuals = TRUE)`.
#' @param png_dir Deprecated and ignored (kept for backward compatibility).
#' @param pd Optional detector metadata (`flowCore::pData(parameters(ff))`) for axis labels.
#'
#' @return Invisibly returns `NULL`; writes report to disk.
#' @export
#' @examples
#' \dontrun{
#' # results_df <- data.table::rbindlist(lapply(unmixed_list, `[[`, "data"))
#' generate_qc_report(
#'   results_df = results_df,
#'   M = M,
#'   output_file = file.path("spectreasy_outputs", "Sample_QC_Report.pdf")
#' )
#' }
generate_qc_report <- function(results_df, M, output_file = file.path("spectreasy_outputs", "Sample_QC_Report.pdf"), res_list = NULL, png_dir = NULL, pd = NULL) {
    message("Generating spectreasy Summary Report...")
    if (!is.null(png_dir)) {
        warning("png_dir is deprecated and ignored; report output is PDF-only.")
    }
    results_df <- as.data.frame(results_df, stringsAsFactors = FALSE)
    if (!("File" %in% colnames(results_df))) {
        stop("results_df must contain a 'File' column.")
    }
    out_dir <- dirname(output_file)
    if (!is.na(out_dir) && nzchar(out_dir) && out_dir != ".") {
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    }
    grDevices::pdf(output_file, width = 11, height = 8.5)
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
    
    grid::grid.newpage()
    file_counts <- table(results_df$File)
    median_events_txt <- if (length(file_counts) > 0) {
        as.character(stats::median(as.numeric(file_counts)))
    } else {
        "NA"
    }
    summary_txt <- paste0(
        "spectreasy: Spectral Unmixing Quality Control Report\n",
        "Generated on: ", Sys.time(), "\n\n",
        "Files Processed: ", length(unique(results_df$File)), "\n",
        "Total Markers: ", nrow(M), "\n",
        "Median Events per File: ", median_events_txt
    )
    grid::grid.text(summary_txt, x = 0.5, y = 0.6, just = "center", gp = grid::gpar(fontsize = 15))
    
    message("  - Adding spectra overlay...")
    p_spectra <- plot_spectra(M, pd = pd, output_file = NULL)
    if (!is.null(p_spectra)) print(p_spectra)
    
    if (!is.null(res_list)) {
        message("  - Adding detector-level residual diagnostics...")
        rep_res <- if (!is.null(res_list$residuals)) res_list else res_list[[1]]
        p_det <- plot_detector_residuals(rep_res, M = M, top_n = 50, output_file = NULL)
        if (!is.null(p_det)) print(p_det)
    }
    
    message("  - Adding Spread Matrix...")
    ssm <- calculate_ssm(M)
    p_ssm <- plot_ssm(ssm, output_file = NULL)
    if (!is.null(p_ssm)) print(p_ssm)
    
    message("  - Adding NPS diagnostics...")
    nps_scores <- calculate_nps(results_df)
    p_nps <- plot_nps(nps_scores, output_file = NULL)
    if (!is.null(p_nps)) print(p_nps)
    
    message("  - Adding Recommendations page...")
    grid::grid.newpage()
    high_spread <- which(ssm > 10, arr.ind = TRUE)
    spread_msgs <- if(nrow(high_spread) > 0) sapply(1:nrow(high_spread), function(i) paste0("- ", rownames(ssm)[high_spread[i,1]], " spreads noise heavily into ", colnames(ssm)[high_spread[i,2]])) else "- No extreme noise spread detected between markers."
    nps_summary <- nps_scores |>
        dplyr::group_by(Marker) |>
        dplyr::summarise(Max_NPS = max(NPS, na.rm = TRUE), .groups = "drop")
    nps_summary <- nps_summary[is.finite(nps_summary$Max_NPS), , drop = FALSE]
    high_nps <- nps_summary |>
        dplyr::arrange(dplyr::desc(Max_NPS)) |>
        dplyr::slice_head(n = min(5, nrow(nps_summary)))
    nps_msgs <- if (nrow(high_nps) > 0) {
        paste0("- Highest negative-population spread observed in: ", paste(high_nps$Marker, collapse = ", "))
    } else {
        "- Negative-population spread could not be ranked from the provided data."
    }
    wrap_lines <- function(x, width = 80) paste(strwrap(x, width = width), collapse = "\n")
    rec_txt <- paste0("spectreasy: Conclusions & Recommendations\n\n", "Spectral Spread Analysis:\n", wrap_lines(paste(spread_msgs, collapse = "\n")), "\n\n", "Negative Population Spread:\n", wrap_lines(nps_msgs), "\n\n", "General Recommendations:\n", "1. Review SCC gating and histogram selection for controls that retain unusually few events.\n", "2. Markers with high spread should not be used to resolve dim co-expressed populations.\n", "3. If detector residuals show laser-specific patterns, check instrument calibration or missing fluorophores.")
    grid::grid.text(rec_txt, x = 0.1, y = 0.9, just = c("left", "top"), gp = grid::gpar(fontsize = 11, lineheight = 1.2))
    
    grDevices::dev.off()
    message("Report saved to: ", output_file)
}
