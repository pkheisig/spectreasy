.report_sample_nxn_interactive_data <- function(results_df,
                                                markers,
                                                max_points = 1000,
                                                transform = c("none", "asinh"),
                                                asinh_cofactor = 150,
                                                axis_limit = NULL,
                                                all_samples = FALSE) {
    transform <- .match_arg_ci(transform, c("none", "asinh"), "transform")
    markers <- intersect(as.character(markers), colnames(results_df))
    if (!("File" %in% colnames(results_df)) || length(markers) < 2L) return(NULL)
    sample_names <- unique(as.character(results_df$File))
    sample_names <- sample_names[!is.na(sample_names) & nzchar(sample_names)]
    if (!isTRUE(all_samples) && length(sample_names)) sample_names <- sample_names[1]
    fixed_limits <- .compute_qc_report_fixed_scatter_limits(axis_limit, transform, asinh_cofactor)
    samples <- lapply(sample_names, function(sample_name) {
        values <- suppressWarnings(as.matrix(results_df[results_df$File == sample_name, markers, drop = FALSE]))
        storage.mode(values) <- "double"
        if (nrow(values) > max_points) {
            keep <- unique(round(seq(1, nrow(values), length.out = max_points)))
            values <- values[keep, , drop = FALSE]
        }
        if (identical(transform, "asinh")) values <- asinh(values / asinh_cofactor)
        limits <- lapply(seq_along(markers), function(index) {
            if (!is.null(fixed_limits)) fixed_limits else .compute_qc_report_scatter_limits(values[, index])
        })
        list(name = sample_name, values = values, limits = limits)
    })
    list(
        markers = markers,
        transform = transform,
        asinh_cofactor = asinh_cofactor,
        samples = samples
    )
}

#' Collect sample QC report data
#' @param results Unmixed sample results from [unmix_samples()] or a compatible
#'   combined data frame with a `File` column.
#' @param M Reference matrix used for sample unmixing.
#' @param unmixing_method Method used for sample unmixing.
#' @param res_list Optional precomputed residual result list.
#' @param pd Optional detector parameter metadata.
#' @param matrix_source Optional source path for `M`.
#' @param detector_noise_file Optional detector noise CSV path.
#' @param spectral_variant_library_file Optional spectral variant library path.
#' @param output_dir Optional directory containing unmixed FCS outputs.
#' @param qc_metrics_dir Optional directory containing QC metric CSV files.
#' @param artifact_paths Named list of additional artifact paths.
#' @param warnings Character vector of captured report warnings.
#' @param run_settings Named list of workflow and report settings.
#' @param project_path Project path displayed in the report.
#' @param max_events_per_sample Maximum retained events per sample.
#' @param overview_files_per_page Maximum samples combined in PDF-equivalent
#'   NPS and reconstruction overview plots.
#' @param matrix_markers_per_page Marker cutoff used by the PDF-equivalent
#'   similarity plot pages.
#' @param sample_nxn_rows_per_page Marker rows/columns per PDF NxN plot page.
#'   HTML output instead renders one coherent matrix per selected sample in a
#'   separate companion HTML file.
#' @param sample_nxn_max_points Maximum plotted events per sample.
#' @param sample_nxn_transform NxN transform, `"none"` or `"asinh"`.
#' @param sample_nxn_asinh_cofactor Asinh cofactor for NxN plots.
#' @param sample_nxn_axis_limit Optional symmetric NxN axis limit.
#' @param nxn_all_samples Whether to include NxN pages for every sample.
#' @param report_per_sample Whether the HTML report should expose a sample
#'   selector and sample-specific report content.
#' @param plot_dir Directory used for cached report PNG files.
#' @return A `spectreasy_sample_report_data` object.
#' @export
collect_sample_report_data <- function(results, M, unmixing_method=NULL, res_list=NULL, pd=NULL,
                                       matrix_source=NULL, detector_noise_file=NULL, spectral_variant_library_file=NULL,
                                       output_dir=NULL, qc_metrics_dir=NULL, artifact_paths=list(), warnings=character(),
                                       run_settings=list(), project_path=getwd(), max_events_per_sample=1000,
                                       overview_files_per_page=15, matrix_markers_per_page=20,
                                       sample_nxn_rows_per_page=10, sample_nxn_max_points=max_events_per_sample,
                                       sample_nxn_transform=c("none","asinh"), sample_nxn_asinh_cofactor=150,
                                       sample_nxn_axis_limit=NULL, nxn_all_samples=FALSE,
                                       report_per_sample=FALSE, plot_dir=NULL) {
    if (is.null(results)) stop("No unmixed results provided for the sample HTML report.",call.=FALSE)
    if (is.null(M)) stop("No reference matrix provided for the sample HTML report.",call.=FALSE)
    M <- .as_reference_matrix(M,"M")
    method <- .normalize_unmix_method(unmixing_method %||% attr(results,"method") %||% "AutoSpectral")
    sample_nxn_transform <- .match_arg_ci(sample_nxn_transform, c("none", "asinh"), "sample_nxn_transform")
    report_results <- if(is.data.frame(results)) .cap_qc_report_results_df(results,max_events_per_sample) else .cap_qc_report_results_list(results,max_events_per_sample)
    results_df <- .normalize_qc_report_results_df(report_results)
    if(!"File" %in% colnames(results_df)) stop("results must contain a 'File' column.",call.=FALSE)
    if(is.null(res_list) && is.list(report_results) && !is.data.frame(report_results)) res_list <- report_results
    af_rows <- grepl("^AF($|_)",rownames(M),ignore.case=TRUE)
    markers <- setdiff(colnames(results_df),.get_result_metadata_columns(colnames(results_df)))
    samples <- unique(as.character(results_df$File))
    full_counts <- .qc_report_file_counts(results)
    event_counts <- as.integer(full_counts[samples])
    event_counts[is.na(event_counts)] <- as.integer(table(factor(results_df$File,levels=samples)))[is.na(event_counts)]
    sample_table <- data.frame(sample=samples,event_count=event_counts,detector_compatible=TRUE,unmixed=TRUE,stringsAsFactors=FALSE)
    similarity <- if(sum(!af_rows)>1L) calculate_similarity_matrix(M[!af_rows,,drop=FALSE]) else NULL
    spread <- NULL
    nps <- if(!identical(method,"NNLS")) tryCatch(calculate_nps(results_df),error=function(e) NULL) else NULL
    if (.report_has_rows(nps) && "Marker" %in% colnames(nps)) {
        nps <- nps[!grepl("^AF($|_)", nps$Marker, ignore.case = TRUE), , drop = FALSE]
    }
    detector_rms <- if(!is.null(res_list)) tryCatch(.compute_qc_report_detector_rms(res_list,M=M,pd=pd,unmixing_method=method),error=function(e) NULL) else NULL
    reconstruction <- if(!is.null(res_list)) tryCatch(.compute_qc_report_sample_rms(res_list,M=M,unmixing_method=method),error=function(e) NULL) else NULL
    qc_metrics_dir <- .prepare_qc_report_metrics_dir(qc_metrics_dir)
    if (!is.null(qc_metrics_dir)) {
        if (any(!af_rows)) {
            .write_qc_report_matrix_metric(
                M[!af_rows, , drop = FALSE],
                file.path(qc_metrics_dir, "reference_spectra.csv"),
                row_id = "fluorophore"
            )
        }
        if (any(af_rows)) {
            .write_qc_report_matrix_metric(
                M[af_rows, , drop = FALSE],
                file.path(qc_metrics_dir, "af_bank_spectra.csv"),
                row_id = "af_band"
            )
        }
        if (!is.null(similarity)) {
            .write_qc_report_matrix_metric(
                similarity,
                file.path(qc_metrics_dir, "fluorophore_spectral_similarity.csv"),
                row_id = "fluorophore"
            )
        }
        if (.report_has_rows(nps)) {
            .write_qc_report_csv(
                nps,
                file.path(qc_metrics_dir, "negative_population_spread.csv")
            )
        }
        if (.report_has_rows(detector_rms)) {
            .write_qc_report_csv(
                detector_rms,
                file.path(qc_metrics_dir, "rms_residual_per_detector.csv")
            )
        }
        if (.report_has_rows(reconstruction)) {
            .write_qc_report_csv(
                reconstruction,
                file.path(qc_metrics_dir, "detector_reconstruction_error.csv")
            )
        }
    }
    plot_dir <- plot_dir %||% tempfile("spectreasy_sample_html_plots_")
    dir.create(plot_dir,recursive=TRUE,showWarnings=FALSE)
    reference_file <- if(any(!af_rows)) .report_plot_file(plot_spectra(M[!af_rows,,drop=FALSE],pd=pd,output_file=NULL),file.path(plot_dir,"reference_spectra.png")) else NULL
    similarity_files <- character()
    if(!is.null(similarity)) {
        similarity_pages <- .build_qc_report_matrix_pages(
            similarity,
            plot_fun = plot_similarity_matrix,
            max_markers_per_page = matrix_markers_per_page,
            item_label = "Markers"
        )
        if(length(similarity_pages)) similarity_files <- vapply(seq_along(similarity_pages), function(i) {
            .report_plot_file(similarity_pages[[i]], file.path(plot_dir, sprintf("similarity_matrix_%02d.png", i)))
        }, character(1))
    }
    nps_files <- character()
    if(.report_has_rows(nps)) {
        nps_pages <- .build_qc_report_nps_pages(nps, max_files_per_page = overview_files_per_page)
        if(length(nps_pages)) nps_files <- vapply(seq_along(nps_pages), function(i) {
            .report_plot_file(nps_pages[[i]],file.path(plot_dir,sprintf("negative_population_spread_%02d.png",i)))
        }, character(1))
    }
    nxn_markers <- markers[!grepl("^AF($|_)",markers,ignore.case=TRUE)]
    include_all_nxn <- isTRUE(nxn_all_samples) || isTRUE(report_per_sample)
    nxn <- .build_qc_report_sample_scatter_pages(
        results_df,
        markers = nxn_markers,
        rows_per_page = max(2L, length(nxn_markers)),
        max_points = sample_nxn_max_points,
        transform = sample_nxn_transform,
        asinh_cofactor = sample_nxn_asinh_cofactor,
        axis_limit = sample_nxn_axis_limit,
        all_samples = include_all_nxn
    )
    selected_samples <- samples
    if (!include_all_nxn && length(selected_samples)) selected_samples <- selected_samples[1]
    if (length(nxn)) names(nxn) <- rep_len(selected_samples, length(nxn))
    nxn_files <- character()
    if(length(nxn)) {
        nxn_width <- max(12, min(24, length(nxn_markers) * 0.9))
        nxn_files <- vapply(seq_along(nxn), function(i) {
            .report_plot_file(nxn[[i]],file.path(plot_dir,sprintf("sample_nxn_full_%02d.png",i)),width=nxn_width,height=nxn_width)
        }, character(1))
        names(nxn_files) <- names(nxn)
    }
    detector_rms_file <- NULL
    reconstruction_files <- character()
    if(!is.null(res_list)) {
        p <- tryCatch(plot_detector_rms_residuals(res_list,M=M,pd=pd,output_file=NULL,unmixing_method=method),error=function(e) NULL)
        detector_rms_file <- .report_plot_file(p,file.path(plot_dir,"detector_rms.png"))
        reconstruction_pages <- tryCatch(.build_qc_report_rms_pages(res_list,M=M,max_files_per_page=overview_files_per_page,unmixing_method=method),error=function(e) list())
        if(length(reconstruction_pages)) reconstruction_files <- vapply(seq_along(reconstruction_pages), function(i) {
            .report_plot_file(reconstruction_pages[[i]],file.path(plot_dir,sprintf("detector_reconstruction_%02d.png",i)))
        }, character(1))
    }
    plot_files <- .report_existing_paths(c(reference_file, similarity_files, nps_files, nxn_files, detector_rms_file, reconstruction_files))
    fcs_paths <- if(!is.null(output_dir) && dir.exists(output_dir)) list.files(output_dir,pattern="\\.fcs$",full.names=TRUE,ignore.case=TRUE) else character()
    metric_paths <- if(!is.null(qc_metrics_dir) && dir.exists(qc_metrics_dir)) list.files(qc_metrics_dir,pattern="\\.csv$",full.names=TRUE,ignore.case=TRUE) else character()
    source_paths <- unique(c(matrix_source,detector_noise_file,spectral_variant_library_file,unlist(artifact_paths,recursive=TRUE,use.names=FALSE)))
    nxn_interactive <- .report_sample_nxn_interactive_data(
        results_df = results_df,
        markers = nxn_markers,
        max_points = sample_nxn_max_points,
        transform = sample_nxn_transform,
        asinh_cofactor = sample_nxn_asinh_cofactor,
        axis_limit = sample_nxn_axis_limit,
        all_samples = include_all_nxn
    )
    out <- list(report_type="Sample QC",project_path=normalizePath(project_path,mustWork=FALSE),created_at=Sys.time(),version=as.character(utils::packageVersion("spectreasy")),unmixing_method=method,
        matrix=M,matrix_source=matrix_source,matrix_preview=.report_matrix_preview(M),peak_detectors=.report_peak_detectors(M),detector_metadata=pd %||% attr(M,"detector_pd"),detector_noise_file=detector_noise_file,
        residual_metadata=list(metric=if(.uses_wls_residual_metric(method)) "wls_weighted_rms" else "raw_rms"),af_metadata=attr(results,"blind_af_info") %||% attr(M,"blind_af_info"),spectral_variant_metadata=attr(results,"spectral_variant_library"),spectral_variant_library_file=spectral_variant_library_file,
        samples=sample_table,markers=markers,detectors=colnames(M),warnings=unique(as.character(warnings)),nps=nps,nps_note=if(identical(method,"NNLS")) "Skipped for NNLS because constrained outputs are non-negative by construction." else NULL,
        similarity=similarity,spread=spread,detector_rms=detector_rms,reconstruction_error=reconstruction,nxn_interactive=nxn_interactive,plots=unique(plot_files),plot_manifest=list(reference=reference_file,similarity=similarity_files,nps=nps_files,nxn=nxn_files,detector_rms=detector_rms_file,reconstruction=reconstruction_files),artifacts=.report_artifacts(c(source_paths,fcs_paths,metric_paths,plot_files)),source_fingerprint=.report_source_fingerprint(source_paths),run_settings=c(list(report_per_sample=isTRUE(report_per_sample)),run_settings),
        counts=list(samples=length(samples),markers=sum(!af_rows),detectors=ncol(M),af_bands=sum(af_rows)),input_status=if(isTRUE(attr(M,"adjusted"))) "Adjusted" else if(isTRUE(attr(M,"synthetic"))) "Synthetic" else "Measured")
    class(out) <- c("spectreasy_sample_report_data","spectreasy_report_data","list")
    out$ai_qc <- collect_ai_qc(
        samples = results, sample_report_data = out, M = M,
        scope = "sample", privacy = "standard", reference = "none",
        project_dir = project_path, generated_at = out$created_at
    )
    out$ai_qc_summary <- list(
        status = out$ai_qc$overall_summary$status,
        grade_counts = out$ai_qc$grade_summary$counts,
        profile = out$ai_qc$quality_reference$profile,
        privacy = out$ai_qc$privacy$mode,
        top_findings = out$ai_qc$overall_summary$top_findings
    )
    out
}
