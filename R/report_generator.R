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
#'   If omitted, `attr(M, "detector_pd")` is used when available.
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

    M <- .as_reference_matrix(M, "M")
    if (is.null(pd)) {
        pd_attr <- attr(M, "detector_pd")
        if (is.data.frame(pd_attr)) {
            pd <- pd_attr
        }
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

    draw_plot_page <- function(p) {
        if (is.null(p)) return(invisible(NULL))
        grid::grid.newpage()
        grid::grid.draw(ggplot2::ggplotGrob(p))
        invisible(NULL)
    }

    grid::grid.newpage()
    file_counts <- table(results_df$File)
    median_events_txt <- if (length(file_counts) > 0) {
        as.character(stats::median(as.numeric(file_counts)))
    } else {
        "NA"
    }

    keep_non_af <- !grepl("^AF($|_)", rownames(M), ignore.case = TRUE)
    M_no_af <- M[keep_non_af, , drop = FALSE]

    summary_txt <- paste0(
        "spectreasy: Spectral Unmixing Quality Control Report\n",
        "Generated on: ", Sys.time(), "\n\n",
        "Files Processed: ", length(unique(results_df$File)), "\n",
        "Total Markers (non-AF): ", nrow(M_no_af), "\n",
        "Median Events per File: ", median_events_txt
    )
    grid::grid.text(summary_txt, x = 0.5, y = 0.6, just = "center", gp = grid::gpar(fontsize = 15))

    message("  - Adding spectra overlay...")
    if (nrow(M_no_af) > 0) {
        p_spectra <- plot_spectra(M_no_af, pd = pd, output_file = NULL)
        draw_plot_page(p_spectra)
    }

    if (!is.null(res_list)) {
        message("  - Adding detector-level residual diagnostics...")
        rep_res <- if (!is.null(res_list$residuals)) res_list else res_list[[1]]
        p_det <- plot_detector_residuals(rep_res, M = if (nrow(M_no_af) > 0) M_no_af else M, top_n = 50, output_file = NULL)
        draw_plot_page(p_det)
    }

    message("  - Adding Spread Matrix...")
    if (nrow(M_no_af) > 1) {
        ssm <- calculate_ssm(M_no_af)
        p_ssm <- plot_ssm(ssm, output_file = NULL)
        draw_plot_page(p_ssm)
    }

    message("  - Adding NPS diagnostics...")
    nps_scores <- calculate_nps(results_df)
    nps_scores <- nps_scores[!grepl("^AF($|_)", nps_scores$Marker, ignore.case = TRUE), , drop = FALSE]
    if (nrow(nps_scores) > 0) {
        p_nps <- plot_nps(nps_scores, output_file = NULL)
        draw_plot_page(p_nps)
    }

    message("Report saved to: ", output_file)
}
