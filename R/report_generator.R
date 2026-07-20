#' Generate a Full Sample PDF Report
#'
#' Creates a multi-page report summarizing unmixing quality, including spectra,
#' RMS residuals per detector, matrix diagnostics, NPS, per-sample NxN marker
#' scatter pages, and detector reconstruction error.
#'
#' `qc_samples()` expects a combined data frame or the raw list returned
#' by [unmix_samples()]. In the usual workflow, pass the unmixed results object
#' directly.
#'
#' `M` should be the reference matrix used for the unmixing context, supplied as
#' a numeric matrix or detector-column data frame. In the usual
#' `unmix_controls()` workflow, pass `ctrl$M` or load
#' `scc_reference_matrix.csv`. Do not pass the path to
#' `scc_unmixing_matrix.csv` here.
#'
#' @param results Combined unmixed data frame, or the list returned by
#'   [unmix_samples()]. `qc_samples()` will automatically bind
#'   per-sample `$data` elements when needed.
#' @param M Reference matrix used for report context. Must be a numeric matrix or
#'   a data frame with detector columns, for example `ctrl$M` from
#'   `unmix_controls()` or `utils::read.csv("scc_reference_matrix.csv",
#'   check.names = FALSE)`.
#' @param unmixing_matrix_file Optional CSV path to a saved reference matrix.
#'   Used when `M` is not supplied. By default this points to the reference matrix
#'   produced by [unmix_controls()] (`"scc_reference_matrix.csv"`).
#' @param output_file Output report path. Defaults to
#'   `"spectreasy_outputs/unmix_samples/qc_samples_report.html"`.
#' @param unmixing_method Unmixing method used to create `results`
#'   (`"AutoSpectral"`, `"Spectreasy"`, `"OLS"`, `"WLS"`, `"RWLS"`, or
#'   `"NNLS"`). When `"NNLS"`, the negative population spread page is skipped
#'   because constrained NNLS results are non-negative by construction.
#' @param res_list Optional residual object/list from `calc_residuals(..., return_residuals = TRUE)`.
#' @param pd Optional detector metadata (`flowCore::pData(parameters(ff))`) for axis labels.
#'   If omitted, `attr(M, "detector_pd")` is used when available.
#' @param max_events_per_sample Maximum events per sample used for report-wide
#'   plots and diagnostics. Defaults to 1000 to keep large FCS reports
#'   responsive. Set `NULL` to use all events.
#' @param overview_files_per_page Maximum files shown on each overall detector
#'   reconstruction error and NPS overview page.
#' @param matrix_markers_per_page Marker cutoff for similarity and spread
#'   matrix pages. Up to this many markers are shown on one page; larger panels
#'   are split into two balanced pages up to twice this cutoff, and three
#'   balanced pages above that.
#' @param sample_nxn_rows_per_page Number of marker rows and columns to show per per-sample NxN page block.
#'   Defaults to 10, which standardizes geometry across pages and samples.
#' @param sample_nxn_max_points Maximum cells sampled per sample for each NxN page.
#' @param sample_nxn_transform One of `"none"` or `"asinh"` for per-sample NxN pages.
#' @param sample_nxn_asinh_cofactor Cofactor used when `sample_nxn_transform = "asinh"`.
#' @param sample_nxn_axis_limit Optional fixed symmetric NxN scatter axis limit.
#'   The default `NULL` uses local per-panel ranges. Use `1e5` for
#'   `c(-1e5, 1e5)` on every NxN panel.
#' @param nxn_all_samples Logical; if `TRUE`, include per-sample NxN pages for all samples.
#'   If `FALSE` (default), only include NxN pages for the first sample in `results`.
#' @param qc_plot_dir Directory where report PNG plots are written when
#'   `save_qc_pngs = TRUE`.
#' @param save_qc_pngs Logical; if `TRUE`, save report plot pages as PNG files
#'   alongside the report.
#' @param qc_metrics_dir Optional directory where plot-ready QC metric
#'   CSVs are written alongside the report.
#' @param report_format Report format, `"html"` (default) or `"pdf"`.
#'   Matching is case-insensitive.
#' @param report_per_sample Logical; when `TRUE`, PDF output writes one report
#'   per sample and HTML output provides a sample selector. Default is `FALSE`.
#' @param overwrite HTML collision policy: create a versioned filename
#'   (default), overwrite, or error. Existing PDF behavior is unchanged.
#' @param report_run_settings Additional workflow settings recorded in HTML.
#' @param report_artifact_paths Additional input/output paths recorded in HTML.
#' @param project_path Project directory recorded in generated report metadata.
#'
#' @return Invisibly returns a list with `output_file`, `qc_plot_dir`, and
#'   `qc_metrics_dir`; writes report artifacts to disk.
#' @export
#' @examples
#' M_demo <- rbind(
#'   FITC = c(1.00, 0.20, 0.05),
#'   PE = c(0.10, 1.00, 0.20),
#'   APC = c(0.05, 0.15, 1.00)
#' )
#' colnames(M_demo) <- c("B2-A", "YG1-A", "R1-A")
#'
#' results <- data.frame(
#'   File = rep(c("sample_a", "sample_b"), each = 120),
#'   FITC = c(rnorm(120, 2, 0.4), rnorm(120, 0.1, 0.2)),
#'   PE = c(rnorm(120, 0.2, 0.2), rnorm(120, 2.5, 0.5)),
#'   APC = rnorm(240, 0.3, 0.3)
#' )
#'
#' # Typical workflow after unmix_samples():
#' # qc_samples(
#' #   results = unmixed,
#' #   M = ctrl$M
#' # )
#'
#' pdf_file <- tempfile(fileext = ".pdf")
#' qc_samples(results = results, M = M_demo, output_file = pdf_file)
#' file.exists(pdf_file)
qc_samples <- function(results,
                       M = NULL,
                       unmixing_matrix_file = file.path("spectreasy_outputs", "unmix_controls", "scc_reference_matrix.csv"),
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
                       overwrite = c("version", "overwrite", "error"),
                       report_run_settings = list(),
                       report_artifact_paths = list(),
                       project_path = getwd()) {
    output_file_missing <- missing(output_file)
    report_format_missing <- missing(report_format)
    output_spec <- .report_output_spec(
        output_file,
        if (report_format_missing) NULL else report_format,
        default_format = "html",
        output_missing = output_file_missing
    )
    output_file <- output_spec$path
    report_per_sample <- .normalize_scalar_logical(report_per_sample, "report_per_sample")
    if (is.null(output_file) || !nzchar(trimws(as.character(output_file)[1]))) {
        stop("Please supply output_file to save the QC report.", call. = FALSE)
    }
    sample_nxn_transform <- .match_arg_ci(sample_nxn_transform, c("none", "asinh"), "sample_nxn_transform")
    method_attr <- attr(results, "method")
    if (is.null(unmixing_method)) {
        unmixing_method <- if (!is.null(method_attr)) method_attr else "AutoSpectral"
    }
    unmixing_method <- .normalize_unmix_method(unmixing_method)

    .spectreasy_console_header("sample QC report")
    .spectreasy_console_field("Report", .spectreasy_console_path(output_file))

    if (!is.null(M)) {
        M <- .as_reference_matrix(M, "M")
    } else if (!is.null(unmixing_matrix_file)) {
        if (!file.exists(unmixing_matrix_file)) {
            .spectreasy_stop_missing_file(unmixing_matrix_file, label = "unmixing_matrix_file")
        }
        .stop_if_static_unmixing_matrix_path(unmixing_matrix_file, arg_name = "unmixing_matrix_file")
        M <- .read_unmixing_matrix_csv(unmixing_matrix_file)
        M <- .as_reference_matrix(M, "M")
    }

    if (is.null(M)) {
        stop("No reference matrix provided. Supply either M or a valid unmixing_matrix_file.")
    }

    if (is.null(pd)) {
        pd_attr <- attr(M, "detector_pd")
        if (is.data.frame(pd_attr)) {
            pd <- pd_attr
        }
    }

    entry_options <- list(
        max_events_per_sample = max_events_per_sample,
        overview_files_per_page = overview_files_per_page,
        matrix_markers_per_page = matrix_markers_per_page,
        sample_nxn_rows_per_page = sample_nxn_rows_per_page,
        sample_nxn_max_points = sample_nxn_max_points,
        sample_nxn_transform = sample_nxn_transform,
        sample_nxn_asinh_cofactor = sample_nxn_asinh_cofactor,
        sample_nxn_axis_limit = sample_nxn_axis_limit,
        nxn_all_samples = nxn_all_samples,
        qc_plot_dir = qc_plot_dir,
        save_qc_pngs = save_qc_pngs,
        qc_metrics_dir = qc_metrics_dir,
        report_per_sample = report_per_sample,
        overwrite = overwrite,
        report_run_settings = report_run_settings,
        report_artifact_paths = report_artifact_paths,
        project_path = project_path
    )
    if (identical(output_spec$format, "pdf") && isTRUE(report_per_sample)) {
        return(.qc_samples_per_sample_pdf(
            results, M, output_file, unmixing_method, res_list, pd, entry_options
        ))
    }
    if (identical(output_spec$format, "html")) {
        return(.qc_samples_html(
            results, M, output_file, unmixing_method, res_list, pd,
            unmixing_matrix_file, entry_options
        ))
    }

    .qc_samples_pdf(
        results = results,
        M = M,
        output_file = output_file,
        method = unmixing_method,
        res_list = res_list,
        pd = pd,
        options = entry_options
    )

}
