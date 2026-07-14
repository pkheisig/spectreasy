# Shared data collection and template rendering for HTML QC reports.

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

.report_scalar <- function(x, default = "Missing") {
    if (is.null(x) || length(x) == 0L || is.na(x[1]) || !nzchar(trimws(as.character(x[1])))) default else as.character(x[1])
}

.report_setting_equal <- function(value, default) {
    isTRUE(all.equal(value, default, check.attributes = FALSE))
}

.report_setting_value <- function(value) {
    if (is.logical(value) && length(value) == 1L && !is.na(value)) return(if (value) "Yes" else "No")
    .report_scalar(value)
}

.report_run_setting_schema <- function(report_data) {
    if (inherits(report_data, "spectreasy_control_report_data")) {
        defaults <- list(
            auto_create_mapping = TRUE,
            auto_unknown_fluor_policy = "by_channel",
            unmix_scatter_panel_size_mm = 30,
            seed = NULL,
            af_min_cluster_events = 20,
            af_min_cluster_proportion = 0.005,
            rwls_max_iter = 1L,
            n_threads = 1L,
            save_qc_png = FALSE,
            gating_mode = "interactive",
            scc_background_method = "scatter_knn",
            scc_background_k = 2L,
            spectral_variant_som_nodes = 16L,
            spectral_variant_top_k = 3L,
            spectral_variant_cosine_threshold = 0.98,
            spectral_variant_max_variants = 8L,
            spectral_variant_min_events = 50L,
            spectreasy_weight_quantile = 0.9,
            autospectral_n_candidates = 1000L,
            autospectral_n_spectral = 200L,
            autospectral_min_events = 10L,
            autospectral_refine = FALSE,
            unmix_scatter_max_points = 1000,
            unmix_scatter_axis_limit = NULL,
            save_qc_pngs = FALSE
        )
        labels <- c(
            auto_create_mapping = "Auto-create mapping", auto_unknown_fluor_policy = "Unknown fluor policy",
            unmix_scatter_panel_size_mm = "Scatter panel size (mm)", seed = "Seed",
            af_min_cluster_events = "Minimum AF cluster events", af_min_cluster_proportion = "Minimum AF cluster proportion",
            rwls_max_iter = "RWLS iterations", n_threads = "Threads", save_qc_png = "Save QC PNGs",
            gating_mode = "Gating mode", scc_background_method = "SCC background method",
            scc_background_k = "SCC background neighbours", spectral_variant_som_nodes = "Spectral variant SOM nodes",
            spectral_variant_top_k = "Spectral variant top K", spectral_variant_cosine_threshold = "Variant cosine threshold",
            spectral_variant_max_variants = "Maximum spectral variants", spectral_variant_min_events = "Minimum variant events",
            spectreasy_weight_quantile = "Spectreasy weight quantile", autospectral_n_candidates = "AutoSpectral candidates",
            autospectral_n_spectral = "AutoSpectral spectra", autospectral_min_events = "AutoSpectral minimum events",
            autospectral_refine = "AutoSpectral refinement", unmix_scatter_max_points = "Control NxN events",
            unmix_scatter_axis_limit = "Control NxN axis limit", save_qc_pngs = "Save report PNGs"
        )
    } else {
        defaults <- list(
            rwls_max_iter = 1L,
            n_threads = 1L,
            spectral_variant_top_k = 3L,
            spectral_variant_min_abundance = 1,
            spectral_variant_positive_fraction = 0.02,
            spectral_variant_min_improvement = 0.01,
            spectreasy_weight_quantile = 0.9,
            estimate_af = FALSE,
            write_fcs = TRUE,
            save_qc_plots = FALSE,
            plot_n_events = 10000L,
            chunk_size = 50000L,
            seed = NULL,
            max_events_per_sample = 1000,
            overview_files_per_page = 15,
            matrix_markers_per_page = 20,
            sample_nxn_rows_per_page = 10,
            sample_nxn_max_points = 1000,
            sample_nxn_transform = "none",
            sample_nxn_asinh_cofactor = 150,
            sample_nxn_axis_limit = NULL,
            nxn_all_samples = FALSE,
            save_qc_pngs = FALSE
        )
        labels <- c(
            rwls_max_iter = "RWLS iterations", n_threads = "Threads", spectral_variant_top_k = "Spectral variant top K",
            spectral_variant_min_abundance = "Minimum variant abundance", spectral_variant_positive_fraction = "Variant positive fraction",
            spectral_variant_min_improvement = "Minimum variant improvement", spectreasy_weight_quantile = "Spectreasy weight quantile",
            estimate_af = "Estimate sample AF", write_fcs = "Write FCS", save_qc_plots = "Save QC plots",
            plot_n_events = "Plot events per sample", chunk_size = "Chunk size", seed = "Seed",
            max_events_per_sample = "Report events per sample", overview_files_per_page = "Overview samples per page",
            matrix_markers_per_page = "Matrix markers per page", sample_nxn_rows_per_page = "NxN rows per page",
            sample_nxn_max_points = "NxN events per sample", sample_nxn_transform = "NxN transform",
            sample_nxn_asinh_cofactor = "NxN asinh cofactor", sample_nxn_axis_limit = "NxN axis limit",
            nxn_all_samples = "NxN for all samples", save_qc_pngs = "Save report PNGs"
        )
    }
    list(defaults = defaults, labels = labels)
}

.report_changed_run_settings <- function(report_data) {
    settings <- report_data$run_settings %||% list()
    if (!length(settings) || is.null(names(settings))) return(list())
    settings <- settings[!duplicated(names(settings), fromLast = TRUE)]
    schema <- .report_run_setting_schema(report_data)
    if (inherits(report_data, "spectreasy_sample_report_data") && "max_events_per_sample" %in% names(settings)) {
        schema$defaults$sample_nxn_max_points <- settings[["max_events_per_sample"]]
    }
    candidates <- intersect(names(schema$defaults), names(settings))
    changed <- candidates[!vapply(candidates, function(name) {
        .report_setting_equal(settings[[name]], schema$defaults[[name]])
    }, logical(1))]
    values <- lapply(changed, function(name) .report_setting_value(settings[[name]]))
    names(values) <- unname(schema$labels[changed])
    values
}

.report_html_escape <- function(x) {
    x <- as.character(x)
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x <- gsub("\"", "&quot;", x, fixed = TRUE)
    x
}

.report_output_spec <- function(output_file, report_format = NULL, default_format = "html", output_missing = FALSE) {
    if (is.null(output_file) || length(output_file) != 1L || is.na(output_file) ||
        !is.character(output_file) || !nzchar(trimws(output_file))) {
        stop("output_file must be a single non-empty file path.", call. = FALSE)
    }
    ext <- tolower(tools::file_ext(.report_scalar(output_file, "")))
    inferred <- if (ext %in% c("html", "htm")) "html" else if (identical(ext, "pdf")) "pdf" else NULL
    if (is.null(report_format) || length(report_format) == 0L || is.na(report_format[1]) || !nzchar(report_format[1])) {
        report_format <- inferred %||% default_format
    }
    report_format <- .match_arg_ci(report_format, c("html", "pdf"), "report_format")
    if (isTRUE(output_missing)) {
        output_file <- sub("\\.(html?|pdf)$", paste0(".", report_format), output_file, ignore.case = TRUE)
        ext <- report_format
    }
    if (!ext %in% c("html", "htm", "pdf")) {
        output_file <- paste0(output_file, ".", report_format)
        ext <- report_format
    }
    if ((identical(report_format, "html") && !ext %in% c("html", "htm")) ||
        (identical(report_format, "pdf") && !identical(ext, "pdf"))) {
        stop("report_format = '", report_format, "' conflicts with output_file extension '.", ext, "'.", call. = FALSE)
    }
    list(format = report_format, path = as.character(output_file)[1])
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
    overwrite <- .match_arg_ci(overwrite, c("version", "overwrite", "error"), "overwrite")
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

.report_plot_file <- function(plot, path, width=11, height=7, dpi=150) {
    if (is.null(plot)) return(NULL)
    dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
    .with_known_qc_plot_warnings_suppressed(ggplot2::ggsave(path, plot=plot, width=width, height=height, units="in", dpi=dpi))
    path
}

.report_nxn_plot_size <- function(plot, width = 10) {
    n_rows <- suppressWarnings(as.numeric(attr(plot, "spectreasy_nxn_rows")))
    n_cols <- suppressWarnings(as.numeric(attr(plot, "spectreasy_nxn_cols")))
    if (length(n_rows) != 1L || length(n_cols) != 1L || !is.finite(n_rows) || !is.finite(n_cols) || n_rows <= 0 || n_cols <= 0) {
        return(c(width = width, height = width))
    }
    c(width = width, height = max(2.5, width * n_rows / n_cols))
}

.report_collect_existing_plots <- function(path) {
    if (is.null(path) || !dir.exists(path)) return(character())
    sort(list.files(path, pattern="\\.png$", full.names=TRUE, recursive=TRUE, ignore.case=TRUE))
}

.report_is_embedded_image <- function(path) {
    !is.na(path) & (startsWith(path, "data:image/png;base64,") | startsWith(path, "data:image/jpeg;base64,"))
}

.report_existing_paths <- function(paths) {
    paths <- unique(as.character(paths))
    nonempty <- !is.na(paths) & nzchar(paths)
    embedded <- nonempty & .report_is_embedded_image(paths)
    exists <- embedded
    exists[nonempty & !embedded] <- file.exists(paths[nonempty & !embedded])
    paths[exists]
}

.report_control_caption <- function(row) {
    if (is.null(row) || !nrow(row)) return(NULL)
    field <- function(name, default = NULL) {
        if (!name %in% colnames(row)) return(default)
        value <- row[[name]][1]
        if (length(value) == 0L || is.na(value) || !nzchar(trimws(as.character(value)))) default else as.character(value)
    }
    parts <- c(
        if (!is.null(field("type"))) paste0("Type: ", field("type")) else NULL,
        if (!is.null(field("peak_channel"))) paste0("Peak channel: ", field("peak_channel")) else NULL,
        if (!is.null(field("n_total"))) paste0("Events: ", field("n_total")) else NULL,
        if (!is.null(field("n_scatter_gated"))) {
            pct <- field("scatter_gate_pct")
            paste0("Scatter gate: ", field("n_scatter_gated"), if (!is.null(pct)) paste0(" (", pct, "%)") else "")
        } else NULL,
        if (!is.null(field("n_final"))) {
            pct <- field("histogram_gate_pct")
            paste0("Final gate: ", field("n_final"), if (!is.null(pct)) paste0(" (", pct, "%)") else "")
        } else NULL
    )
    if (length(parts)) paste(parts, collapse = " | ") else NULL
}

.report_control_panels <- function(report_plot_dir, qc_summary = data.frame()) {
    if (is.null(report_plot_dir) || !dir.exists(report_plot_dir)) return(list())
    qc_summary <- as.data.frame(qc_summary %||% data.frame(), stringsAsFactors = FALSE)
    samples <- if (nrow(qc_summary) && "sample" %in% colnames(qc_summary)) {
        as.character(qc_summary$sample)
    } else {
        sub("_spectrum\\.png$", "", basename(list.files(file.path(report_plot_dir, "spectrum"), pattern = "_spectrum\\.png$", full.names = TRUE, ignore.case = TRUE)))
    }
    samples <- unique(samples[!is.na(samples) & nzchar(samples)])
    panels <- lapply(samples, function(sample_id) {
        row <- if (nrow(qc_summary) && "sample" %in% colnames(qc_summary)) {
            qc_summary[as.character(qc_summary$sample) == sample_id, , drop = FALSE][1, , drop = FALSE]
        } else data.frame()
        fluor <- if (nrow(row) && "fluorophore" %in% colnames(row)) trimws(as.character(row$fluorophore[1])) else ""
        marker <- if (nrow(row) && "marker" %in% colnames(row)) trimws(as.character(row$marker[1])) else ""
        title <- if (nzchar(fluor)) {
            if (nzchar(marker) && !identical(tolower(marker), tolower(fluor))) paste0(fluor, " / ", marker, " (", sample_id, ")") else paste0(fluor, " (", sample_id, ")")
        } else sample_id
        singlet <- file.path(report_plot_dir, "singlet", paste0(sample_id, "_singlet.png"))
        has_singlet <- file.exists(singlet)
        gate_title <- if (has_singlet) "Histogram" else "Peak-channel histogram gate"
        image_paths <- c(
            file.path(report_plot_dir, "fsc_ssc", paste0(sample_id, "_fsc_ssc.png")),
            if (has_singlet) singlet else NULL,
            file.path(report_plot_dir, "histogram", paste0(sample_id, "_histogram.png")),
            file.path(report_plot_dir, "spectrum", paste0(sample_id, "_spectrum.png"))
        )
        image_titles <- c(
            if (has_singlet) "Cell gate" else "FSC/SSC gate",
            if (has_singlet) "Singlet gate" else NULL,
            gate_title,
            "Per-event spectrum distribution"
        )
        keep <- file.exists(image_paths)
        list(
            id = sample_id,
            title = title,
            caption = .report_control_caption(row),
            paths = image_paths[keep],
            image_titles = image_titles[keep]
        )
    })
    Filter(function(panel) length(panel$paths) > 0L, panels)
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
    af_info <- af_bank_info %||% attr(M,"af_bank_info")
    af_source_count <- suppressWarnings(as.numeric(af_info$source_count %||% 1)[1])
    af_derived_bands <- suppressWarnings(as.numeric(af_info$derived_bands %||% sum(af_rows))[1])
    show_af_page <- any(af_rows) && (isTRUE(af_source_count > 1) || isTRUE(af_derived_bands > 1))
    reference_matrix <- .scc_reference_overlay_matrix(M)
    spectra <- if (nrow(reference_matrix)) plot_spectra(reference_matrix,pd=pd,output_file=NULL) else NULL
    af_plot <- if (show_af_page) {
        (plots$af %||% plot_spectra(M[af_rows,,drop=FALSE],pd=pd,output_file=NULL)) +
            ggplot2::labs(title = "Autofluorescence Band Spectra Overlay") +
            ggplot2::theme(legend.position = "none")
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
                                       sample_nxn_axis_limit=NULL, nxn_all_samples=FALSE, plot_dir=NULL) {
    if (is.null(results)) stop("No unmixed results provided for the sample HTML report.",call.=FALSE)
    if (is.null(M)) stop("No reference matrix provided for the sample HTML report.",call.=FALSE)
    M <- .as_reference_matrix(M,"M")
    method <- .normalize_unmix_method(unmixing_method %||% attr(results,"method") %||% "Spectreasy")
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
    nxn <- .build_qc_report_sample_scatter_pages(
        results_df,
        markers = nxn_markers,
        rows_per_page = max(2L, length(nxn_markers)),
        max_points = sample_nxn_max_points,
        transform = sample_nxn_transform,
        asinh_cofactor = sample_nxn_asinh_cofactor,
        axis_limit = sample_nxn_axis_limit,
        all_samples = nxn_all_samples
    )
    selected_samples <- samples
    if (!isTRUE(nxn_all_samples) && length(selected_samples)) selected_samples <- selected_samples[1]
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
    out <- list(report_type="Sample QC",project_path=normalizePath(project_path,mustWork=FALSE),created_at=Sys.time(),version=as.character(utils::packageVersion("spectreasy")),unmixing_method=method,
        matrix=M,matrix_source=matrix_source,matrix_preview=.report_matrix_preview(M),peak_detectors=.report_peak_detectors(M),detector_metadata=pd %||% attr(M,"detector_pd"),detector_noise_file=detector_noise_file,
        residual_metadata=list(metric=if(.uses_wls_residual_metric(method)) "wls_weighted_rms" else "raw_rms"),af_metadata=attr(results,"blind_af_info") %||% attr(M,"blind_af_info"),spectral_variant_metadata=attr(results,"spectral_variant_library"),spectral_variant_library_file=spectral_variant_library_file,
        samples=sample_table,markers=markers,detectors=colnames(M),warnings=unique(as.character(warnings)),nps=nps,nps_note=if(identical(method,"NNLS")) "Skipped for NNLS because constrained outputs are non-negative by construction." else NULL,
        similarity=similarity,spread=spread,detector_rms=detector_rms,reconstruction_error=reconstruction,plots=unique(plot_files),plot_manifest=list(reference=reference_file,similarity=similarity_files,nps=nps_files,nxn=nxn_files,detector_rms=detector_rms_file,reconstruction=reconstruction_files),artifacts=.report_artifacts(c(source_paths,fcs_paths,metric_paths,plot_files)),source_fingerprint=.report_source_fingerprint(source_paths),run_settings=run_settings,
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
    if (length(path) == 1L && .report_is_embedded_image(path)) return(path)
    if(!file.exists(path)) return(NULL)
    mime <- if(tolower(tools::file_ext(path)) %in% c("jpg","jpeg")) "image/jpeg" else "image/png"
    paste0("data:",mime,";base64,",base64enc::base64encode(path))
}

.report_embed_plot_manifest <- function(report_data) {
    manifest <- report_data$plot_manifest %||% list()
    embed_paths <- function(paths) {
        if (is.null(paths) || !length(paths)) return(paths)
        embedded <- vapply(as.character(paths), function(path) {
            uri <- .report_image_uri(path)
            if (is.null(uri)) path else uri
        }, character(1), USE.NAMES = FALSE)
        names(embedded) <- names(paths)
        embedded
    }
    fields <- intersect(c("af", "reference", "similarity", "nps", "nxn", "detector_rms", "reconstruction"), names(manifest))
    if (inherits(report_data, "spectreasy_control_report_data")) {
        fields <- setdiff(fields, "nxn")
    }
    for (field in fields) {
        manifest[[field]] <- embed_paths(manifest[[field]])
    }
    if (length(manifest$controls)) {
        manifest$controls <- lapply(manifest$controls, function(panel) {
            panel$paths <- embed_paths(panel$paths)
            panel
        })
    }
    report_data$plot_manifest <- manifest
    report_data
}

.report_plots_html <- function(paths, titles = NULL, empty="No plot artifact was available.", image_class = NULL) {
    paths <- .report_existing_paths(paths)
    if(!length(paths)) return(paste0("<p class=\"note\">",empty,"</p>"))
    if (is.null(titles)) titles <- rep("", length(paths))
    titles <- rep_len(as.character(titles), length(paths))
    class_attr <- if (!is.null(image_class) && nzchar(image_class)) paste0(" class=\"", .report_html_escape(image_class), "\"") else ""
    paste(vapply(seq_along(paths),function(i) {
        paste0("<figure class=\"plot-wrap\"><img",class_attr," alt=\"",.report_html_escape(if (nzchar(titles[i])) titles[i] else basename(paths[i])),"\" src=\"",.report_image_uri(paths[i]),"\"></figure>")
    },character(1)),collapse="")
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
    plots <- .report_existing_paths(plots)
    if (length(plots)) return(.report_plots_html(plots))
    "<p class=\"note\">The corresponding PDF-style plot was not available for this run.</p>"
}

.report_control_panels_html <- function(panels) {
    if (!length(panels)) return("<p class=\"note\">No per-control PDF diagnostic panels were available.</p>")
    paste(vapply(panels, function(panel) {
        caption <- if (!is.null(panel$caption) && nzchar(panel$caption)) paste0("<p class=\"control-caption\">", .report_html_escape(panel$caption), "</p>") else ""
        paste0(
            "<article class=\"control-review\"><h3>", .report_html_escape(panel$title), "</h3>", caption,
            "<div class=\"control-plot-grid\">",
            .report_plots_html(panel$paths, titles = panel$image_titles),
            "</div></article>"
        )
    }, character(1)), collapse = "")
}

.report_nxn_link_html <- function(paths, report_type = c("sample", "control")) {
    report_type <- .match_arg_ci(report_type, c("sample", "control"), "report_type")
    if (is.null(paths) || !length(paths)) return("<p class=\"note\">No NxN scatter matrix was available for this run.</p>")
    labels <- names(paths)
    if (is.null(labels) || length(labels) != length(paths) || any(!nzchar(labels))) {
        labels <- if (identical(report_type, "control")) "Controls" else paste("Sample", seq_along(paths))
    }
    links <- paste(vapply(seq_along(paths), function(i) {
        basename_path <- basename(paths[i])
        label <- if (identical(report_type, "control")) "Open high-resolution control NxN matrix" else paste0("Open ", labels[i], " NxN matrix")
        paste0(
            "<a class=\"companion-link\" data-companion=\"", .report_html_escape(basename_path),
            "\" href=\"", .report_html_escape(basename_path), "\">", .report_html_escape(label), " &rarr;</a>"
        )
    }, character(1)), collapse = "")
    description <- if (identical(report_type, "control")) {
        "The complete control NxN scatter matrix is saved separately as a high-resolution PNG for detailed review and zooming."
    } else {
        "Each selected sample has its own compact one-page NxN viewer containing one complete matrix with zoom controls and two-dimensional scrolling."
    }
    paste0("<p>", description, "</p><div class=\"companion-links\">", links, "</div>")
}

.report_toc_label <- function(id, title) {
    labels <- c(
        af = "AF bank",
        spectra = "Spectra",
        similarity = "Similarity",
        gating = "Gating",
        residual = "Residuals",
        reconstruction = "Reconstruction",
        nps = "NPS",
        nxn = "NxN matrix"
    )
    label <- unname(labels[id])
    if (length(label) && !is.na(label)) label else title
}

.report_filename_slug <- function(x, fallback = "matrix") {
    slug <- gsub("[^A-Za-z0-9._-]+", "_", trimws(as.character(x)))
    slug <- gsub("^[_\\.]+|[_\\.]+$", "", slug)
    slug[!nzchar(slug)] <- fallback
    make.unique(slug, sep = "_")
}

.report_nxn_entries <- function(report_data) {
    raw_paths <- report_data$plot_manifest$nxn %||% character()
    keep <- !is.na(raw_paths) & nzchar(raw_paths)
    embedded <- keep & .report_is_embedded_image(raw_paths)
    keep[keep & !embedded] <- file.exists(raw_paths[keep & !embedded])
    paths <- as.character(raw_paths[keep])
    if (!length(paths)) return(data.frame(label = character(), path = character(), stringsAsFactors = FALSE))
    labels <- names(raw_paths)[keep]
    is_control <- inherits(report_data, "spectreasy_control_report_data")
    if (is.null(labels) || length(labels) != length(paths)) labels <- rep("", length(paths))
    empty <- !nzchar(labels)
    labels[empty] <- if (is_control) "Controls" else paste("Sample", which(empty))
    if (is_control && length(paths) != 1L) {
        stop("Control HTML reports require one coherent NxN matrix plot, not PDF-style split pages. Recollect the report data with the current collector.", call. = FALSE)
    }
    data.frame(label = labels, path = paths, stringsAsFactors = FALSE)
}

.report_nxn_companion_paths <- function(output_file, report_data) {
    entries <- .report_nxn_entries(report_data)
    if (!nrow(entries)) return(character())
    stem <- tools::file_path_sans_ext(basename(output_file))
    paths <- if (inherits(report_data, "spectreasy_control_report_data")) {
        file.path(dirname(output_file), paste0(stem, "_nxn.png"))
    } else {
        slugs <- .report_filename_slug(entries$label, fallback = "sample")
        file.path(dirname(output_file), paste0(stem, "_nxn_", slugs, ".html"))
    }
    stats::setNames(paths, entries$label)
}

.report_resolve_bundle <- function(path, overwrite, report_data = NULL) {
    overwrite <- .match_arg_ci(overwrite, c("version", "overwrite", "error"), "overwrite")
    paths_exist <- function(candidate) {
        companion_paths <- if (is.null(report_data)) character() else .report_nxn_companion_paths(candidate, report_data)
        any(file.exists(c(candidate, companion_paths)))
    }
    if (identical(overwrite, "overwrite") || !paths_exist(path)) return(path)
    if (identical(overwrite, "error")) {
        stop("Report already exists: ", path, ". Use overwrite = 'overwrite' or 'version'.", call. = FALSE)
    }
    stem <- tools::file_path_sans_ext(basename(path))
    ext <- tools::file_ext(path)
    for (i in seq_len(9999L)) {
        candidate <- file.path(dirname(path), sprintf("%s_v%03d.%s", stem, i, ext))
        if (!paths_exist(candidate)) return(candidate)
    }
    stop("Could not find a free versioned report filename for: ", path, call. = FALSE)
}

.report_render_nxn_companions <- function(report_data, output_files) {
    entries <- .report_nxn_entries(report_data)
    if (!nrow(entries)) return(character())
    if (length(output_files) != nrow(entries)) stop("NxN companion path count does not match the matrix plot count.", call. = FALSE)
    is_control <- inherits(report_data, "spectreasy_control_report_data")
    if (is_control) {
        rendered <- vapply(seq_len(nrow(entries)), function(i) {
            source <- entries$path[i]
            target <- output_files[i]
            if (.report_is_embedded_image(source)) {
                payload <- sub("^data:image/(png|jpeg);base64,", "", source)
                base64enc::base64decode(payload, output = target)
            } else if (identical(normalizePath(source, mustWork = FALSE), normalizePath(target, mustWork = FALSE))) {
                if (!file.exists(target)) stop("Control NxN PNG not found: ", target, call. = FALSE)
            } else if (!file.copy(source, target, overwrite = TRUE)) {
                stop("Could not write the control NxN PNG: ", target, call. = FALSE)
            }
            normalizePath(target, mustWork = TRUE)
        }, character(1))
        return(stats::setNames(rendered, entries$label))
    }
    template <- system.file("report_templates", "nxn_report.html", package = "spectreasy")
    if (!nzchar(template)) template <- file.path("inst", "report_templates", "nxn_report.html")
    if (!file.exists(template)) stop("Spectreasy NxN HTML report template was not found.", call. = FALSE)
    rendered <- vapply(seq_len(nrow(entries)), function(i) {
        html <- paste(readLines(template, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
        report_label <- paste0("Sample NxN scatter matrix: ", entries$label[i])
        replacements <- list(
            "{{REPORT_TITLE}}" = paste("Spectreasy", report_label),
            "{{REPORT_CLASS}}" = "report-sample",
            "{{REPORT_HEADING}}" = report_label,
            "{{REPORT_META}}" = "One complete matrix - use the zoom controls and two-dimensional scrolling for detailed review",
            "{{MATRIX_LABEL}}" = .report_html_escape(entries$label[i]),
            "{{MATRIX_IMAGE}}" = .report_image_uri(entries$path[i])
        )
        for (key in names(replacements)) html <- gsub(key, replacements[[key]], html, fixed = TRUE)
        writeLines(html, output_files[i], useBytes = TRUE)
        normalizePath(output_files[i], mustWork = TRUE)
    }, character(1))
    stats::setNames(rendered, entries$label)
}

.report_sections <- function(x, nxn_companion_files = NULL) {
    if(inherits(x,"spectreasy_control_report_data")) {
        manifest <- x$plot_manifest %||% list()
        companions <- nxn_companion_files %||% x$nxn_companion_files
        spectra_paths <- c(.report_existing_paths(manifest$af), .report_existing_paths(manifest$reference))
        spectra_titles <- c(
            if (length(.report_existing_paths(manifest$af))) "Autofluorescence band spectra overlay" else NULL,
            if (length(.report_existing_paths(manifest$reference))) "Reference spectra overlay" else NULL
        )
        sections <- Filter(Negate(is.null), list(
            spectra=if (length(spectra_paths)) c("Spectra",.report_plots_html(spectra_paths, titles = spectra_titles)) else NULL,
            similarity=if (length(.report_existing_paths(manifest$similarity))) c("Fluorophore similarity",.report_plots_html(manifest$similarity, titles = "Fluorophore spectral similarity")) else NULL,
            nxn=if (length(.report_existing_paths(manifest$nxn))) c("Control NxN scatter matrix",.report_nxn_link_html(companions, report_type = "control")) else NULL,
            gating=if (length(manifest$controls)) c("Gating",.report_control_panels_html(manifest$controls)) else NULL
        ))
        sections
    } else {
        manifest <- x$plot_manifest %||% list()
        companions <- nxn_companion_files %||% x$nxn_companion_files
        sections <- Filter(Negate(is.null), list(
            spectra=if (length(.report_existing_paths(manifest$reference))) c("Reference spectra",.report_plots_html(manifest$reference, titles = "Reference spectra overlay")) else NULL,
            similarity=if (length(.report_existing_paths(manifest$similarity))) c("Fluorophore similarity",.report_plots_html(manifest$similarity, titles = "Fluorophore spectral similarity")) else NULL,
            nps=if (!is.null(x$nps_note) || length(.report_existing_paths(manifest$nps))) c("Negative population spread", if (!is.null(x$nps_note)) paste0("<p class=\"note\">", .report_html_escape(x$nps_note), "</p>") else .report_plots_html(manifest$nps, titles = "Negative population spread")) else NULL,
            nxn=if (length(.report_existing_paths(manifest$nxn))) c("Sample NxN scatter matrices",.report_nxn_link_html(companions, report_type = "sample")) else NULL,
            residual=if (length(.report_existing_paths(manifest$detector_rms))) c("Detector RMS residual",.report_plots_html(manifest$detector_rms, titles = "RMS residual per detector")) else NULL,
            reconstruction=if (length(.report_existing_paths(manifest$reconstruction))) c("Reconstruction error",.report_plots_html(manifest$reconstruction, titles = "Reconstruction error per sample")) else NULL
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
    include_nxn <- nrow(.report_nxn_entries(report_data)) > 0L
    output_file <- .report_resolve_bundle(spec$path, overwrite, report_data = report_data)
    dir.create(dirname(output_file),recursive=TRUE,showWarnings=FALSE)
    nxn_companion_files <- if (include_nxn) .report_nxn_companion_paths(output_file, report_data) else character()
    if (length(nxn_companion_files)) report_data$nxn_companion_files <- nxn_companion_files
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
    sections <- .report_sections(report_data, nxn_companion_files = nxn_companion_files)
    toc_titles <- vapply(sections,`[`,character(1),1)
    toc_labels <- vapply(seq_along(sections),function(i) .report_toc_label(names(sections)[i],toc_titles[i]),character(1))
    toc <- paste0("<a href=\"#",names(sections),"\">",.report_html_escape(toc_labels),"</a>",collapse="")
    body <- paste(vapply(seq_along(sections),function(i) .report_section(names(sections)[i],sections[[i]][1],sections[[i]][2]),character(1)),collapse="")
    counts <- report_data$counts
    summary <- c(list(`Created`=format(report_data$created_at,"%Y-%m-%d %H:%M:%S %Z"),`Spectreasy`=report_data$version,`Method`=report_data$unmixing_method,`Cytometer`=report_data$cytometer %||% "Not recorded"),counts,.report_changed_run_settings(report_data))
    summary_html <- paste(vapply(names(summary),function(nm) paste0("<div><small>",.report_html_escape(nm),"</small><strong>",.report_html_escape(.report_scalar(summary[[nm]])),"</strong></div>"),character(1)),collapse="")
    header <- paste0("<header><span class=\"kicker\">Spectreasy &middot; ",.report_html_escape(report_data$report_type),"</span><h1>",.report_html_escape(report_data$report_type)," report</h1><p class=\"paths\">",.report_html_escape(report_data$project_path),"</p><div class=\"summary-grid\">",summary_html,"</div></header>")
    report_class <- if (inherits(report_data, "spectreasy_control_report_data")) "report-control" else "report-sample"
    replacements <- list("{{REPORT_TITLE}}"=paste("Spectreasy",report_data$report_type,"report"),"{{REPORT_CLASS}}"=report_class,"{{REPORT_TOC}}"=toc,"{{REPORT_HEADER}}"=header,"{{REPORT_BODY}}"=body)
    for(key in names(replacements)) html <- gsub(key,replacements[[key]],html,fixed=TRUE)
    writeLines(html,output_file,useBytes=TRUE)
    companion_files <- character()
    if (length(nxn_companion_files)) companion_files <- .report_render_nxn_companions(report_data, nxn_companion_files)
    if (inherits(report_data, "spectreasy_control_report_data") && length(companion_files)) {
        report_data$plot_manifest$nxn <- companion_files
        report_data$nxn_companion_files <- companion_files
        report_data$artifacts <- rbind(
            report_data$artifacts,
            .report_artifacts(companion_files, roles = "High-resolution control NxN matrix")
        )
    }
    source_paths <- report_data$source_fingerprint$path %||% character()
    metadata <- list(output_file=normalizePath(output_file,mustWork=TRUE),companion_files=companion_files,assets_dir=NULL,self_contained=!include_nxn,format="html",report_type=report_data$report_type,created_at=report_data$created_at,source_fingerprint=report_data$source_fingerprint,stale=.report_is_stale(output_file,source_paths),report_data=report_data)
    class(metadata) <- c("spectreasy_report_result","list")
    invisible(metadata)
}
