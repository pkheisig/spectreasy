# Shared data collection and template rendering for HTML QC reports.

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

.report_scalar <- function(x, default = "Missing") {
    if (is.null(x) || length(x) == 0L || is.na(x[1]) || !nzchar(trimws(as.character(x[1])))) default else as.character(x[1])
}

.report_html_escape <- function(x) {
    x <- as.character(x)
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x <- gsub("\"", "&quot;", x, fixed = TRUE)
    x
}

.report_output_spec <- function(output_file, output_format = NULL, default_format = "pdf", output_missing = FALSE) {
    if (is.null(output_file) || length(output_file) != 1L || is.na(output_file) ||
        !is.character(output_file) || !nzchar(trimws(output_file))) {
        stop("output_file must be a single non-empty file path.", call. = FALSE)
    }
    ext <- tolower(tools::file_ext(.report_scalar(output_file, "")))
    inferred <- if (ext %in% c("html", "htm")) "html" else if (identical(ext, "pdf")) "pdf" else NULL
    if (is.null(output_format) || length(output_format) == 0L || is.na(output_format[1]) || !nzchar(output_format[1])) {
        output_format <- inferred %||% default_format
    }
    output_format <- match.arg(tolower(output_format[1]), c("html", "pdf"))
    if (isTRUE(output_missing)) {
        output_file <- sub("\\.(html?|pdf)$", paste0(".", output_format), output_file, ignore.case = TRUE)
        ext <- output_format
    }
    if (!ext %in% c("html", "htm", "pdf")) {
        output_file <- paste0(output_file, ".", output_format)
        ext <- output_format
    }
    if ((identical(output_format, "html") && !ext %in% c("html", "htm")) ||
        (identical(output_format, "pdf") && !identical(ext, "pdf"))) {
        stop("output_format = '", output_format, "' conflicts with output_file extension '.", ext, "'.", call. = FALSE)
    }
    list(format = output_format, path = as.character(output_file)[1])
}

.report_versioned_path <- function(path) {
    if (!file.exists(path)) return(path)
    stem <- tools::file_path_sans_ext(basename(path))
    ext <- tools::file_ext(path)
    for (i in seq_len(9999L)) {
        candidate <- file.path(dirname(path), sprintf("%s_v%03d.%s", stem, i, ext))
        if (!file.exists(candidate)) return(candidate)
    }
    stop("Could not find a free versioned report filename for: ", path, call. = FALSE)
}

.report_resolve_overwrite <- function(path, overwrite = c("version", "overwrite", "error")) {
    overwrite <- match.arg(overwrite)
    if (!file.exists(path) || identical(overwrite, "overwrite")) return(path)
    if (identical(overwrite, "version")) return(.report_versioned_path(path))
    stop("Report already exists: ", path, ". Use overwrite = 'overwrite' or 'version'.", call. = FALSE)
}

.report_artifacts <- function(paths = character(), roles = NULL) {
    paths <- unique(as.character(unlist(paths, recursive = TRUE, use.names = FALSE)))
    paths <- paths[!is.na(paths) & nzchar(paths)]
    if (!length(paths)) return(data.frame(path=character(), type=character(), role=character(), exists=logical(), modified=character(), stringsAsFactors=FALSE))
    if (is.null(roles)) {
        roles <- ifelse(
            grepl("reference_matrix|unmixing_matrix", paths, ignore.case=TRUE), "Reference or unmixing matrix",
            ifelse(grepl("detector_noise", paths, ignore.case=TRUE), "Detector noise metadata",
            ifelse(grepl("spectral_variant", paths, ignore.case=TRUE), "Spectral variant library",
            ifelse(grepl("gate", paths, ignore.case=TRUE), "Gating input or plot",
            ifelse(grepl("\\.fcs$", paths, ignore.case=TRUE), "FCS input or unmixed output",
            ifelse(grepl("\\.csv$", paths, ignore.case=TRUE), "QC metric or mapping CSV", "Report input/output"))))))
    }
    roles <- rep_len(as.character(roles), length(paths))
    exists <- file.exists(paths)
    modified <- rep(NA_character_, length(paths))
    modified[exists] <- format(file.info(paths[exists])$mtime, "%Y-%m-%d %H:%M:%S %Z")
    data.frame(path=paths, type=toupper(tools::file_ext(paths)), role=roles, exists=exists, modified=modified, stringsAsFactors=FALSE)
}

.report_source_fingerprint <- function(paths) {
    paths <- unique(as.character(paths))
    paths <- paths[!is.na(paths) & nzchar(paths)]
    info <- file.info(paths)
    data.frame(path=paths, exists=file.exists(paths), size=ifelse(file.exists(paths), info$size, NA_real_), mtime=ifelse(file.exists(paths), as.numeric(info$mtime), NA_real_), stringsAsFactors=FALSE)
}

.report_is_stale <- function(report_file, source_paths) {
    if (!file.exists(report_file)) return(TRUE)
    sources <- source_paths[file.exists(source_paths)]
    length(sources) > 0L && any(file.info(sources)$mtime > file.info(report_file)$mtime)
}

.report_matrix_preview <- function(M, rows = 12L, cols = 12L) {
    if (is.null(M)) return(data.frame())
    out <- as.data.frame(M[seq_len(min(nrow(M), rows)), seq_len(min(ncol(M), cols)), drop=FALSE], check.names=FALSE)
    out <- data.frame(marker=rownames(M)[seq_len(nrow(out))], out, check.names=FALSE)
    out
}

.report_peak_detectors <- function(M) {
    if (is.null(M) || !nrow(M) || !ncol(M)) return(data.frame())
    idx <- max.col(as.matrix(M), ties.method="first")
    data.frame(marker=rownames(M), peak_detector=colnames(M)[idx], peak_value=as.numeric(M[cbind(seq_len(nrow(M)), idx)]), stringsAsFactors=FALSE)
}

.report_control_mapping <- function(control_file, scc_dir) {
    path <- tryCatch(.resolve_control_file_path(control_file), error=function(e) control_file)
    if (!file.exists(path)) return(data.frame())
    x <- utils::read.csv(path, stringsAsFactors=FALSE, check.names=FALSE)
    wanted <- c("filename","fluorophore","marker","channel","control.type","is.viability","universal.negative")
    for (nm in setdiff(wanted, colnames(x))) x[[nm]] <- ""
    full <- file.path(scc_dir, x$filename)
    x$file.exists <- file.exists(full)
    x$validation.warning <- ""
    x$validation.warning[!x$file.exists] <- "File missing"
    x$validation.warning[!nzchar(trimws(x$fluorophore))] <- paste(x$validation.warning[!nzchar(trimws(x$fluorophore))], "Missing fluorophore")
    universal <- toupper(trimws(as.character(x$universal.negative))) %in% c("TRUE", "YES", "1")
    viability <- toupper(trimws(as.character(x$is.viability))) %in% c("TRUE", "YES", "1")
    bead_negative <- tolower(x$control.type)=="beads" & grepl("negative|unstained|blank", x$filename, ignore.case=TRUE)
    x$control.class <- ifelse(
        grepl("^AF($|_)", x$fluorophore, ignore.case=TRUE), "unstained / AF",
        ifelse(viability, "viability control",
        ifelse(universal, "universal negative",
        ifelse(bead_negative, "bead negative",
        ifelse(tolower(x$control.type)=="beads", "bead control", ifelse(tolower(x$control.type)=="cells", "cell control", "unspecified")))))
    )
    x[, c(wanted, "control.class", "file.exists", "validation.warning"), drop=FALSE]
}

.report_plot_file <- function(plot, path, width=11, height=7) {
    if (is.null(plot)) return(NULL)
    dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
    .with_known_qc_plot_warnings_suppressed(ggplot2::ggsave(path, plot=plot, width=width, height=height, units="in", dpi=150))
    path
}

.report_collect_existing_plots <- function(path) {
    if (is.null(path) || !dir.exists(path)) return(character())
    sort(list.files(path, pattern="\\.png$", full.names=TRUE, recursive=TRUE, ignore.case=TRUE))
}

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
    nps <- if (!is.null(unmixed_list) && !identical(method,"NNLS")) tryCatch(calculate_nps(.normalize_qc_report_results_df(unmixed_list)), error=function(e) NULL) else NULL
    detector_rms <- if (!is.null(unmixed_list)) tryCatch(.compute_qc_report_detector_rms(unmixed_list,M=M,pd=pd,unmixing_method=method),error=function(e) NULL) else NULL
    reconstruction <- if (!is.null(unmixed_list)) tryCatch(.compute_qc_report_sample_rms(unmixed_list,M=M,unmixing_method=method),error=function(e) NULL) else NULL
    spectra <- plots$spectra %||% if (any(!af_rows)) plot_spectra(M[!af_rows,,drop=FALSE],pd=pd,output_file=NULL) else NULL
    af_plot <- plots$af %||% if (any(af_rows)) plot_spectra(M[af_rows,,drop=FALSE],pd=pd,output_file=NULL) else NULL
    similarity_plot <- if (!is.null(similarity)) plot_similarity_matrix(similarity,output_file=NULL) else NULL
    plot_files <- c(
        .report_plot_file(spectra,file.path(plot_dir,"control_spectra.png")),
        .report_plot_file(af_plot,file.path(plot_dir,"af_bank.png")),
        .report_plot_file(similarity_plot,file.path(plot_dir,"spectral_similarity.png")),
        .report_plot_file(plots$unmixing_scatter %||% NULL,file.path(plot_dir,"scc_unmixing_scatter.png"))
    )
    plot_files <- unique(c(plot_files, .report_collect_existing_plots(report_plot_dir)))
    control_paths <- if (nrow(mapping)) file.path(scc_dir,mapping$filename) else character()
    source_paths <- unique(c(matrix_source, tryCatch(.resolve_control_file_path(control_file),error=function(e) control_file), control_paths, artifact_paths$gate_file))
    artifacts <- .report_artifacts(c(source_paths, unlist(artifact_paths,recursive=TRUE,use.names=FALSE), plot_files))
    out <- list(
        report_type="Control QC", project_path=normalizePath(project_path,mustWork=FALSE), created_at=Sys.time(), version=as.character(utils::packageVersion("spectreasy")),
        unmixing_method=method, cytometer=cytometer, matrix=M, matrix_source=matrix_source, matrix_preview=.report_matrix_preview(M), peak_detectors=.report_peak_detectors(M),
        detector_metadata=pd, detector_noise=attr(M,"detector_noise"), af_bank_info=af_bank_info %||% attr(M,"af_bank_info"), spectral_variant_metadata=attr(M,"spectral_variant_library"),
        mapping=mapping, qc_summary=as.data.frame(qc_summary %||% attr(M,"qc_summary") %||% data.frame(),stringsAsFactors=FALSE), gating_summary=as.data.frame(qc_summary %||% data.frame(),stringsAsFactors=FALSE),
        control_files=control_paths, warnings=unique(as.character(warnings)), similarity=similarity, nps=nps, detector_rms=detector_rms, reconstruction_error=reconstruction,
        plots=plot_files, artifacts=artifacts, source_fingerprint=.report_source_fingerprint(source_paths), run_settings=run_settings,
        counts=list(controls=if(nrow(mapping)) nrow(mapping) else length(control_paths), markers=sum(!af_rows), detectors=ncol(M), af_bands=sum(af_rows)),
        input_status=if(isTRUE(attr(M,"adjusted"))) "Adjusted" else if(isTRUE(attr(M,"synthetic"))) "Synthetic" else "Measured"
    )
    class(out) <- c("spectreasy_control_report_data","spectreasy_report_data","list")
    out
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
#' @param sample_nxn_rows_per_page Marker rows/columns per NxN plot page.
#' @param sample_nxn_max_points Maximum plotted events per sample.
#' @param sample_nxn_transform NxN transform, `"none"` or `"asinh"`.
#' @param sample_nxn_asinh_cofactor Asinh cofactor for NxN plots.
#' @param sample_nxn_axis_limit Optional symmetric NxN axis limit.
#' @param nxn_all_samples Whether to include NxN pages for every sample.
#' @param plot_dir Directory used for cached report PNG files.
#' @return A `spectreasy_sample_report_data` object.
#' @export
collect_sample_report_data <- function(results, M, unmixing_method=NULL, res_list=NULL, pd=NULL,
                                       matrix_source=NULL, detector_noise_file=NULL, spectral_variant_library_file=NULL,
                                       output_dir=NULL, qc_metrics_dir=NULL, artifact_paths=list(), warnings=character(),
                                       run_settings=list(), project_path=getwd(), max_events_per_sample=1000,
                                       sample_nxn_rows_per_page=10, sample_nxn_max_points=max_events_per_sample,
                                       sample_nxn_transform=c("none","asinh"), sample_nxn_asinh_cofactor=150,
                                       sample_nxn_axis_limit=NULL, nxn_all_samples=FALSE, plot_dir=NULL) {
    if (is.null(results)) stop("No unmixed results provided for the sample HTML report.",call.=FALSE)
    if (is.null(M)) stop("No reference matrix provided for the sample HTML report.",call.=FALSE)
    M <- .as_reference_matrix(M,"M")
    method <- .normalize_unmix_method(unmixing_method %||% attr(results,"method") %||% "Spectreasy")
    sample_nxn_transform <- match.arg(sample_nxn_transform)
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
    detector_rms <- if(!is.null(res_list)) tryCatch(.compute_qc_report_detector_rms(res_list,M=M,pd=pd,unmixing_method=method),error=function(e) NULL) else NULL
    reconstruction <- if(!is.null(res_list)) tryCatch(.compute_qc_report_sample_rms(res_list,M=M,unmixing_method=method),error=function(e) NULL) else NULL
    plot_dir <- plot_dir %||% tempfile("spectreasy_sample_html_plots_")
    dir.create(plot_dir,recursive=TRUE,showWarnings=FALSE)
    plot_files <- character()
    if(any(!af_rows)) plot_files <- c(plot_files,.report_plot_file(plot_spectra(M[!af_rows,,drop=FALSE],pd=pd,output_file=NULL),file.path(plot_dir,"reference_spectra.png")))
    if(!is.null(similarity)) plot_files <- c(plot_files,.report_plot_file(plot_similarity_matrix(similarity,output_file=NULL),file.path(plot_dir,"similarity_matrix.png")))
    nxn <- .build_qc_report_sample_scatter_pages(results_df,markers=markers[!grepl("^AF($|_)",markers,ignore.case=TRUE)],rows_per_page=sample_nxn_rows_per_page,max_points=sample_nxn_max_points,transform=sample_nxn_transform,asinh_cofactor=sample_nxn_asinh_cofactor,axis_limit=sample_nxn_axis_limit,all_samples=nxn_all_samples)
    if(length(nxn)) for(i in seq_along(nxn)) plot_files <- c(plot_files,.report_plot_file(nxn[[i]],file.path(plot_dir,sprintf("sample_nxn_%02d.png",i)),width=10,height=10))
    if(!is.null(res_list)) {
        p <- tryCatch(plot_detector_rms_residuals(res_list,M=M,pd=pd,output_file=NULL,unmixing_method=method),error=function(e) NULL)
        plot_files <- c(plot_files,.report_plot_file(p,file.path(plot_dir,"detector_rms.png")))
    }
    fcs_paths <- if(!is.null(output_dir) && dir.exists(output_dir)) list.files(output_dir,pattern="\\.fcs$",full.names=TRUE,ignore.case=TRUE) else character()
    metric_paths <- if(!is.null(qc_metrics_dir) && dir.exists(qc_metrics_dir)) list.files(qc_metrics_dir,pattern="\\.csv$",full.names=TRUE,ignore.case=TRUE) else character()
    source_paths <- unique(c(matrix_source,detector_noise_file,spectral_variant_library_file,unlist(artifact_paths,recursive=TRUE,use.names=FALSE)))
    out <- list(report_type="Sample QC",project_path=normalizePath(project_path,mustWork=FALSE),created_at=Sys.time(),version=as.character(utils::packageVersion("spectreasy")),unmixing_method=method,
        matrix=M,matrix_source=matrix_source,matrix_preview=.report_matrix_preview(M),peak_detectors=.report_peak_detectors(M),detector_metadata=pd %||% attr(M,"detector_pd"),detector_noise_file=detector_noise_file,
        residual_metadata=list(metric=if(.uses_wls_residual_metric(method)) "wls_weighted_rms" else "raw_rms"),af_metadata=attr(results,"blind_af_info") %||% attr(M,"blind_af_info"),spectral_variant_metadata=attr(results,"spectral_variant_library"),spectral_variant_library_file=spectral_variant_library_file,
        samples=sample_table,markers=markers,detectors=colnames(M),warnings=unique(as.character(warnings)),nps=nps,nps_note=if(identical(method,"NNLS")) "Skipped for NNLS because constrained outputs are non-negative by construction." else NULL,
        similarity=similarity,spread=spread,detector_rms=detector_rms,reconstruction_error=reconstruction,plots=unique(plot_files),artifacts=.report_artifacts(c(source_paths,fcs_paths,metric_paths,plot_files)),source_fingerprint=.report_source_fingerprint(source_paths),run_settings=run_settings,
        counts=list(samples=length(samples),markers=sum(!af_rows),detectors=ncol(M),af_bands=sum(af_rows)),input_status=if(isTRUE(attr(M,"adjusted"))) "Adjusted" else if(isTRUE(attr(M,"synthetic"))) "Synthetic" else "Measured")
    class(out) <- c("spectreasy_sample_report_data","spectreasy_report_data","list")
    out
}

.report_table_html <- function(x, empty="Not available for this run.") {
    if(is.null(x) || !is.data.frame(x) || !nrow(x)) return(paste0("<p class=\"note\">",.report_html_escape(empty),"</p>"))
    cells <- lapply(x,function(v){v<-as.character(v);v[is.na(v)|!nzchar(v)]<-"Missing";.report_html_escape(v)})
    header <- paste0("<th>",.report_html_escape(colnames(x)),"</th>",collapse="")
    rows <- vapply(seq_len(nrow(x)),function(i){vals<-vapply(cells,`[`,character(1),i);classes<-ifelse(vals=="Missing"," class=\"missing\"","");paste0("<tr>",paste0("<td",classes,">",vals,"</td>",collapse=""),"</tr>")},character(1))
    paste0("<div class=\"table-wrap\"><table><thead><tr>",header,"</tr></thead><tbody>",paste(rows,collapse=""),"</tbody></table></div>")
}

.report_list_html <- function(x) {
    if(is.null(x) || !length(x)) return("<p class=\"note\">Not available for this run.</p>")
    flat <- unlist(x,recursive=TRUE,use.names=TRUE)
    .report_table_html(data.frame(setting=names(flat),value=as.character(flat),stringsAsFactors=FALSE))
}

.report_image_uri <- function(path) {
    if(!file.exists(path)) return(NULL)
    mime <- if(tolower(tools::file_ext(path)) %in% c("jpg","jpeg")) "image/jpeg" else "image/png"
    paste0("data:",mime,";base64,",base64enc::base64encode(path))
}

.report_plots_html <- function(paths, empty="No plot artifact was available.") {
    paths <- unique(paths[file.exists(paths)])
    if(!length(paths)) return(paste0("<p class=\"note\">",empty,"</p>"))
    paste(vapply(paths,function(path) paste0("<h3>",.report_html_escape(gsub("[_-]+"," ",tools::file_path_sans_ext(basename(path)))),"</h3><div class=\"plot-wrap\"><img loading=\"lazy\" alt=\"",.report_html_escape(basename(path)),"\" src=\"",.report_image_uri(path),"\"></div>"),character(1)),collapse="")
}

.report_section <- function(id,title,body) paste0("<section id=\"",id,"\"><h2><button type=\"button\">",.report_html_escape(title),"</button></h2><div class=\"section-body\">",body,"</div></section>")

.report_qc_summary <- function(x) {
    badges <- c(if(length(x$warnings)) "Warning" else "OK",x$input_status)
    missing <- sum(!x$artifacts$exists)
    if(missing>0) badges <- c(badges,"Missing")
    badge_html <- paste0("<span class=\"badge ",tolower(badges),"\">",badges,"</span>",collapse="")
    warning_html <- if(length(x$warnings)) paste0("<ul class=\"warning-list\"><li>",paste(.report_html_escape(x$warnings),collapse="</li><li>"),"</li></ul>") else "<p class=\"note\">No report-level warnings were captured.</p>"
    paste0("<div class=\"badges\">",badge_html,"</div>",warning_html)
}

.report_has_rows <- function(x) is.data.frame(x) && nrow(x) > 0L

.report_has_values <- function(x) {
    if (is.null(x) || length(x) == 0L) return(FALSE)
    flat <- unlist(x, recursive = TRUE, use.names = FALSE)
    length(flat) > 0L && any(!is.na(flat) & nzchar(trimws(as.character(flat))))
}

.report_metric_html <- function(x, plots = character()) {
    plots <- plots[file.exists(plots)]
    if (length(plots)) return(.report_plots_html(plots))
    .report_table_html(as.data.frame(x %||% data.frame(), check.names = FALSE))
}

.report_sections <- function(x) {
    if(inherits(x,"spectreasy_control_report_data")) {
        gate_plots <- x$plots[grepl("gate|singlet|fsc|histogram", x$plots, ignore.case = TRUE)]
        spectra_plots <- x$plots[grepl("control_spectra|spectrum", x$plots, ignore.case = TRUE)]
        af_plots <- x$plots[grepl("af_bank|af_", x$plots, ignore.case = TRUE)]
        scc_plots <- x$plots[grepl("unmix|scatter", x$plots, ignore.case = TRUE)]
        similarity_plots <- x$plots[grepl("similar", x$plots, ignore.case = TRUE)]
        sections <- Filter(Negate(is.null), list(
            overview=c("Overview",.report_list_html(c(x$counts,list(input_status=x$input_status)))),
            "input-files"=if (.report_has_rows(x$mapping)) c("Input files and mapping",.report_table_html(x$mapping)) else NULL,
            "mapping-validation"=if (.report_has_rows(x$mapping)) c("Control mapping validation",.report_table_html(x$mapping[,intersect(c("filename","file.exists","validation.warning"),colnames(x$mapping)),drop=FALSE])) else NULL,
            gating=if (.report_has_rows(x$gating_summary) || length(gate_plots)) c("Gating summary",paste0(.report_table_html(x$gating_summary),.report_plots_html(gate_plots))) else NULL,
            reference=c("Reference matrix summary",paste0(.report_list_html(list(rows=nrow(x$matrix),detectors=ncol(x$matrix),source=x$matrix_source,status=x$input_status)),.report_table_html(x$peak_detectors))),
            spectra=if (length(spectra_plots)) c("Control spectra",.report_plots_html(spectra_plots)) else NULL,
            af=if (.report_has_values(x$af_bank_info) || length(af_plots)) c("AF bank summary",paste0(.report_list_html(x$af_bank_info),.report_plots_html(af_plots))) else NULL,
            noise=if (.report_has_values(x$detector_noise)) c("Detector noise summary", .report_table_html(as.data.frame(x$detector_noise, check.names = FALSE))) else NULL,
            scc=if (.report_has_rows(x$detector_rms) || .report_has_rows(x$reconstruction_error) || length(scc_plots)) c("SCC unmixing diagnostics", paste0(.report_table_html(x$detector_rms), .report_table_html(x$reconstruction_error), .report_plots_html(scc_plots))) else NULL,
            similarity=if (.report_has_values(x$similarity) || length(similarity_plots)) c("Spectral similarity or library comparison", .report_metric_html(x$similarity, similarity_plots)) else NULL,
            spread=if (.report_has_rows(x$nps)) c("Negative population spread and spread-related diagnostics", .report_table_html(x$nps)) else NULL,
            artifacts=c("Generated artifacts",.report_table_html(x$artifacts)),
            warnings=if (length(x$warnings)) c("Warnings and logs", .report_qc_summary(x)) else NULL,
            settings=if (.report_has_values(x$run_settings)) c("Settings",.report_list_html(x$run_settings)) else NULL
        ))
        sections
    } else {
        reference_plots <- x$plots[grepl("reference_spectra", x$plots, ignore.case = TRUE)]
        rms_plots <- x$plots[grepl("rms", x$plots, ignore.case = TRUE)]
        matrix_plots <- x$plots[grepl("similarity", x$plots, ignore.case = TRUE)]
        nxn_plots <- x$plots[grepl("sample_nxn", x$plots, ignore.case = TRUE)]
        sections <- Filter(Negate(is.null), list(
            overview=c("Overview",.report_list_html(c(x$counts,list(input_status=x$input_status)))),
            samples=c("Input samples",.report_table_html(x$samples)),reference=c("Reference matrix used",paste0(.report_list_html(list(source=x$matrix_source,rows=nrow(x$matrix),detectors=ncol(x$matrix),af_bands=x$counts$af_bands,status=x$input_status,detector_noise_file=x$detector_noise_file,spectral_variant_library=x$spectral_variant_library_file)),.report_plots_html(reference_plots))),
            unmixing=if (.report_has_values(x$run_settings)) c("Unmixing settings",.report_list_html(x$run_settings)) else NULL,
            outputs=c("Sample output summary",.report_table_html(x$artifacts[grepl("FCS|CSV|HTML|PDF",x$artifacts$type),,drop=FALSE])),
            residual=if (.report_has_rows(x$detector_rms) || length(rms_plots)) c("Detector residual summary", .report_metric_html(x$detector_rms, rms_plots)) else NULL,
            nps=if (!is.null(x$nps_note) || .report_has_rows(x$nps)) c("Negative population spread", paste0(if(!is.null(x$nps_note)) paste0("<p class=\"note\">",x$nps_note,"</p>") else "",.report_table_html(x$nps))) else NULL,
            matrix=if (.report_has_values(x$similarity) || .report_has_values(x$spread) || length(matrix_plots)) c("Matrix diagnostics", paste0("<h3>Similarity</h3>", .report_metric_html(x$similarity, matrix_plots), if (.report_has_values(x$spread)) paste0("<h3>Spread</h3>", .report_table_html(as.data.frame(x$spread, check.names = FALSE))) else "")) else NULL,
            nxn=if (length(nxn_plots)) c("Per-sample marker plots", .report_plots_html(nxn_plots)) else NULL,
            reconstruction=if (.report_has_rows(x$reconstruction_error)) c("Detector reconstruction error", .report_table_html(x$reconstruction_error)) else NULL,
            af=if (.report_has_values(x$af_metadata) || .report_has_values(x$spectral_variant_metadata)) c("AF and spectral variant diagnostics", .report_list_html(list(af=x$af_metadata,spectral_variants=x$spectral_variant_metadata))) else NULL,
            artifacts=c("Generated artifacts",.report_table_html(x$artifacts)),
            warnings=if (length(x$warnings)) c("Warnings and logs", .report_qc_summary(x)) else NULL
        ))
        sections
    }
}

#' Render a Spectreasy QC report object as HTML
#' @param report_data Object returned by a Spectreasy report data collector.
#' @param output_file Target `.html` file.
#' @param overwrite Collision policy: `"version"`, `"overwrite"`, or `"error"`.
#' @return Structured report metadata including the final output path, source
#'   fingerprint, self-contained status, and cached report data.
#' @export
render_qc_html_report <- function(report_data, output_file, overwrite=c("version","overwrite","error")) {
    if(!inherits(report_data,"spectreasy_report_data")) stop("report_data must be created by a Spectreasy report data collector.",call.=FALSE)
    spec <- .report_output_spec(output_file,"html")
    output_file <- .report_resolve_overwrite(spec$path,overwrite)
    dir.create(dirname(output_file),recursive=TRUE,showWarnings=FALSE)
    report_data$artifacts <- rbind(
        report_data$artifacts,
        data.frame(
            path = normalizePath(output_file, mustWork = FALSE),
            type = "HTML",
            role = "Rendered QC report",
            exists = TRUE,
            modified = format(report_data$created_at, "%Y-%m-%d %H:%M:%S %Z"),
            stringsAsFactors = FALSE
        )
    )
    template <- system.file("report_templates","qc_report.html",package="spectreasy")
    if(!nzchar(template)) template <- file.path("inst","report_templates","qc_report.html")
    if(!file.exists(template)) stop("Spectreasy HTML report template was not found.",call.=FALSE)
    html <- paste(readLines(template,warn=FALSE,encoding="UTF-8"),collapse="\n")
    sections <- .report_sections(report_data)
    toc <- paste0("<a href=\"#",names(sections),"\">",vapply(sections,`[`,character(1),1),"</a>",collapse="")
    body <- paste(vapply(seq_along(sections),function(i) .report_section(names(sections)[i],sections[[i]][1],sections[[i]][2]),character(1)),collapse="")
    counts <- report_data$counts
    summary <- c(list(`Created`=format(report_data$created_at,"%Y-%m-%d %H:%M:%S %Z"),`Spectreasy`=report_data$version,`Method`=report_data$unmixing_method,`Cytometer`=report_data$cytometer %||% "Not recorded",`Matrix source`=report_data$matrix_source %||% "In-memory",`Report output`=normalizePath(output_file,mustWork=FALSE),`Input status`=report_data$input_status),counts)
    summary_html <- paste(vapply(names(summary),function(nm) paste0("<div><small>",.report_html_escape(nm),"</small><strong>",.report_html_escape(.report_scalar(summary[[nm]])),"</strong></div>"),character(1)),collapse="")
    header <- paste0("<header><span class=\"kicker\">Spectreasy &middot; ",.report_html_escape(report_data$report_type),"</span><h1>",.report_html_escape(report_data$report_type)," report</h1><p class=\"paths\">",.report_html_escape(report_data$project_path),"</p><div class=\"summary-grid\">",summary_html,"</div></header>")
    replacements <- list("{{REPORT_TITLE}}"=paste("Spectreasy",report_data$report_type,"report"),"{{REPORT_TOC}}"=toc,"{{REPORT_HEADER}}"=header,"{{REPORT_BODY}}"=body)
    for(key in names(replacements)) html <- gsub(key,replacements[[key]],html,fixed=TRUE)
    writeLines(html,output_file,useBytes=TRUE)
    source_paths <- report_data$source_fingerprint$path %||% character()
    metadata <- list(output_file=normalizePath(output_file,mustWork=TRUE),assets_dir=NULL,self_contained=TRUE,format="html",report_type=report_data$report_type,created_at=report_data$created_at,source_fingerprint=report_data$source_fingerprint,stale=.report_is_stale(output_file,source_paths),report_data=report_data)
    class(metadata) <- c("spectreasy_report_result","list")
    invisible(metadata)
}
