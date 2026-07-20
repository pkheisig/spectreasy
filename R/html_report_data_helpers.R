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
            spectreasy_weight_quantile = 0.65,
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
            spectreasy_weight_quantile = 0.65,
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
            report_per_sample = FALSE,
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
            nxn_all_samples = "NxN for all samples", report_per_sample = "Report per sample", save_qc_pngs = "Save report PNGs"
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
