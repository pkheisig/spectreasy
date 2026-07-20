#' Collect control QC report data
#'
#' Builds a reusable, structured report object without rerunning control
#' unmixing. The object can be passed to [render_qc_html_report()].
#' @param M Reference matrix used by the control workflow.
#' @param scc_dir Directory containing control FCS files.
#' @param control_file Control mapping CSV.
#' @param cytometer Cytometer identifier or label.
#' @param unmixing_method Method used for control unmixing.
#' @param unmixed_list Optional precomputed control unmixing results.
#' @param qc_summary Optional precomputed control gating summary.
#' @param report_plot_dir Optional directory of existing control QC PNG files.
#' @param pd Optional detector parameter metadata.
#' @param af_bank_info Optional AF bank metadata.
#' @param matrix_source Optional source path for `M`.
#' @param artifact_paths Named list of additional artifact paths.
#' @param plots Named list of reusable `spectra`, `af`, and
#'   `unmixing_scatter` plot objects.
#' @param warnings Character vector of captured report warnings.
#' @param run_settings Named list of workflow and report settings.
#' @param project_path Project path displayed in the report.
#' @param plot_dir Directory used for cached report PNG files.
#' @return A `spectreasy_control_report_data` object.
#' @export
collect_control_report_data <- function(M, scc_dir="scc", control_file="fcs_mapping.csv", cytometer="auto",
                                        unmixing_method="WLS", unmixed_list=NULL, qc_summary=NULL,
                                        report_plot_dir=NULL, pd=NULL, af_bank_info=NULL,
                                        matrix_source=NULL, artifact_paths=list(), plots=list(),
                                        warnings=character(), run_settings=list(), project_path=getwd(),
                                        plot_dir=NULL) {
    if (is.null(M)) stop("No reference matrix provided for the control HTML report.", call.=FALSE)
    M <- .as_reference_matrix(M, "M")
    method <- .normalize_unmix_method(unmixing_method)
    plot_dir <- plot_dir %||% tempfile("spectreasy_control_html_plots_")
    dir.create(plot_dir, recursive=TRUE, showWarnings=FALSE)
    af_rows <- grepl("^AF($|_)", rownames(M), ignore.case=TRUE)
    mapping <- .report_control_mapping(control_file, scc_dir)
    similarity <- if (sum(!af_rows) > 1L) calculate_similarity_matrix(M[!af_rows,,drop=FALSE]) else NULL
    af_info <- af_bank_info %||% attr(M,"af_bank_info")
    af_source_count <- suppressWarnings(as.numeric(af_info$source_count %||% 1)[1])
    af_derived_bands <- suppressWarnings(as.numeric(af_info$derived_bands %||% sum(af_rows))[1])
    show_af_page <- any(af_rows) && (isTRUE(af_source_count > 1) || isTRUE(af_derived_bands > 1))
    reference_matrix <- .scc_reference_overlay_matrix(M)
    spectra <- if (nrow(reference_matrix)) plot_spectra(reference_matrix,pd=pd,output_file=NULL) else NULL
    af_plot <- if (show_af_page) {
        af_matrix <- M[af_rows,,drop=FALSE]
        .style_af_bank_plot(
            plots$af %||% plot_spectra(af_matrix,pd=pd,output_file=NULL),
            rownames(af_matrix)
        ) + ggplot2::labs(title = "Autofluorescence Band Spectra Overlay")
    } else NULL
    similarity_plot <- if (!is.null(similarity)) plot_similarity_matrix(similarity,output_file=NULL) else NULL
    reference_file <- .report_plot_file(spectra,file.path(plot_dir,"reference_spectra.png"))
    af_file <- .report_plot_file(af_plot,file.path(plot_dir,"af_bank.png"))
    similarity_file <- .report_plot_file(similarity_plot,file.path(plot_dir,"spectral_similarity.png"))
    scatter_pages <- list()
    if (!is.null(unmixed_list) && sum(!af_rows) > 1L) {
        sample_to_marker <- NULL
        if (nrow(mapping)) {
            keys <- tools::file_path_sans_ext(as.character(mapping$filename))
            values <- as.character(mapping$fluorophore)
            keep <- nzchar(keys) & nzchar(values)
            sample_to_marker <- stats::setNames(values[keep], keys[keep])
        }
        scatter_pages <- tryCatch(.build_qc_report_control_scatter_pages(
            unmixed_list = unmixed_list,
            sample_to_marker = sample_to_marker,
            markers = rownames(M)[!af_rows],
            rows_per_page = max(2L, sum(!af_rows)),
            max_points_per_sample = 1000,
            point_size = 0.35,
            point_alpha = 0.82
        ), error = function(e) list())
    }
    if (!length(scatter_pages) && !is.null(plots$unmixing_scatter)) scatter_pages <- list(plots$unmixing_scatter)
    scatter_files <- if (length(scatter_pages)) {
        nxn_width <- max(12, min(30, sum(!af_rows) * 0.9))
        .report_plot_file(
            scatter_pages[[1]],
            file.path(plot_dir, "control_nxn_full.png"),
            width = nxn_width,
            height = nxn_width,
            dpi = 300
        )
    } else character()
    if (length(scatter_files)) names(scatter_files) <- "Controls"
    control_panels <- .report_control_panels(report_plot_dir, qc_summary)
    plot_files <- .report_existing_paths(c(reference_file, af_file, similarity_file, scatter_files, unlist(lapply(control_panels, `[[`, "paths"), use.names = FALSE)))
    control_paths <- if (nrow(mapping)) file.path(scc_dir,mapping$filename) else character()
    source_paths <- unique(c(matrix_source, tryCatch(.resolve_control_file_path(control_file),error=function(e) control_file), control_paths, artifact_paths$gate_file))
    artifacts <- .report_artifacts(c(source_paths, unlist(artifact_paths,recursive=TRUE,use.names=FALSE), plot_files))
    out <- list(
        report_type="Control QC", project_path=normalizePath(project_path,mustWork=FALSE), created_at=Sys.time(), version=as.character(utils::packageVersion("spectreasy")),
        unmixing_method=method, cytometer=cytometer, matrix=M, matrix_source=matrix_source, matrix_preview=.report_matrix_preview(M), peak_detectors=.report_peak_detectors(M),
        detector_metadata=pd, detector_noise=attr(M,"detector_noise"), af_bank_info=af_info, spectral_variant_metadata=attr(M,"spectral_variant_library"),
        mapping=mapping, qc_summary=as.data.frame(qc_summary %||% attr(M,"qc_summary") %||% data.frame(),stringsAsFactors=FALSE), gating_summary=as.data.frame(qc_summary %||% data.frame(),stringsAsFactors=FALSE),
        control_files=control_paths, warnings=unique(as.character(warnings)), similarity=similarity, nps=NULL, detector_rms=NULL, reconstruction_error=NULL,
        plots=plot_files, plot_manifest=list(af=af_file, reference=reference_file, similarity=similarity_file, nxn=scatter_files, controls=control_panels), artifacts=artifacts, source_fingerprint=.report_source_fingerprint(source_paths), run_settings=run_settings,
        counts=list(controls=if(nrow(mapping)) nrow(mapping) else length(control_paths), markers=sum(!af_rows), detectors=ncol(M), af_bands=sum(af_rows)),
        input_status=if(isTRUE(attr(M,"adjusted"))) "Adjusted" else if(isTRUE(attr(M,"synthetic"))) "Synthetic" else "Measured"
    )
    class(out) <- c("spectreasy_control_report_data","spectreasy_report_data","list")
    out
}
