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
        label <- if (identical(report_type, "control")) "Open high-resolution control NxN matrix" else "Open interactive sample NxN matrix"
        paste0(
            "<a class=\"companion-link\" data-companion=\"", .report_html_escape(basename_path),
            "\" href=\"", .report_html_escape(basename_path), "\">", .report_html_escape(label), " &rarr;</a>"
        )
    }, character(1)), collapse = "")
    description <- if (identical(report_type, "control")) {
        "The complete control NxN scatter matrix is saved separately as a high-resolution PNG for detailed review and zooming."
    } else {
        "The interactive NxN viewer provides sample selection, color gradients, plot types, point sizing, and adjustable matrix cells."
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
        file.path(dirname(output_file), paste0(stem, "_nxn.html"))
    }
    stats::setNames(paths, if (inherits(report_data, "spectreasy_control_report_data")) entries$label else "Samples")
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
    is_control <- inherits(report_data, "spectreasy_control_report_data")
    if (is_control) {
        if (length(output_files) != nrow(entries)) stop("NxN companion path count does not match the matrix plot count.", call. = FALSE)
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
    if (length(output_files) != 1L) stop("Sample HTML reports require one interactive NxN companion.", call. = FALSE)
    matrix_data <- report_data$nxn_interactive
    if (is.null(matrix_data) || !length(matrix_data$samples)) stop("Interactive sample NxN data was not available.", call. = FALSE)
    matrix_json <- jsonlite::toJSON(matrix_data, auto_unbox = TRUE, dataframe = "rows", matrix = "rowmajor", na = "null", digits = 9)
    matrix_json <- gsub("</", "<\\/", matrix_json, fixed = TRUE)
    html <- paste(readLines(template, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    replacements <- list(
        "{{REPORT_TITLE}}" = "Spectreasy sample NxN scatter matrices",
        "{{REPORT_CLASS}}" = "report-sample",
        "{{REPORT_HEADING}}" = "Sample NxN scatter matrix",
        "{{REPORT_META}}" = "Interactive sample review",
        "{{MATRIX_DATA}}" = matrix_json
    )
    for (key in names(replacements)) html <- gsub(key, replacements[[key]], html, fixed = TRUE)
    writeLines(html, output_files[1], useBytes = TRUE)
    stats::setNames(normalizePath(output_files[1], mustWork = TRUE), "Samples")
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
