.ai_qc_file_hash <- function(path) {
    con <- file(path, "rb")
    on.exit(close(con), add = TRUE)
    bytes <- readBin(con, what = "raw", n = file.info(path)$size)
    paste0(openssl::sha256(bytes))
}

.ai_qc_write_text <- function(text, path) {
    con <- file(path, "wb")
    on.exit(close(con), add = TRUE)
    writeBin(charToRaw(enc2utf8(text)), con)
    path
}

.ai_qc_report_companion_path <- function(report_file) {
    report_file <- normalizePath(as.character(report_file)[[1]], mustWork = FALSE)
    file.path(
        dirname(report_file),
        paste0(tools::file_path_sans_ext(basename(report_file)), "_ai_qc_prompt.txt")
    )
}

.ai_qc_report_csv_label <- function(path, report_file) {
    path <- normalizePath(path, mustWork = FALSE)
    report_dir <- normalizePath(dirname(report_file), mustWork = FALSE)
    prefix <- paste0(report_dir, .Platform$file.sep)
    if (startsWith(path, prefix)) {
        return(gsub("\\\\", "/", substring(path, nchar(prefix) + 1L)))
    }
    basename(path)
}

.qc_prompt_drop_empty <- function(x) {
    x <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
    populated <- function(value) {
        if (is.numeric(value)) any(is.finite(value)) else any(!is.na(value) & nzchar(trimws(as.character(value))))
    }
    if (ncol(x)) x <- x[, vapply(x, populated, logical(1)), drop = FALSE]
    if (nrow(x) && ncol(x)) {
        keep <- apply(x, 1, function(row) any(!is.na(row) & nzchar(trimws(as.character(row)))))
        x <- x[keep, , drop = FALSE]
    }
    rownames(x) <- NULL
    x
}

.qc_prompt_csv_text <- function(x) {
    lines <- character()
    con <- textConnection("lines", "w", local = TRUE)
    utils::write.csv(x, con, row.names = FALSE, na = "")
    close(con)
    paste(lines, collapse = "\n")
}

.qc_prompt_numeric_summary <- function(x) {
    numeric_names <- names(x)[vapply(x, is.numeric, logical(1))]
    rows <- lapply(numeric_names, function(name) {
        values <- as.numeric(x[[name]])
        finite <- values[is.finite(values)]
        if (!length(finite)) return(NULL)
        data.frame(
            variable = name, n = length(finite), missing = sum(!is.finite(values)),
            minimum = min(finite), q25 = as.numeric(stats::quantile(finite, 0.25, names = FALSE)),
            median = stats::median(finite), q75 = as.numeric(stats::quantile(finite, 0.75, names = FALSE)),
            maximum = max(finite), stringsAsFactors = FALSE
        )
    })
    rows <- Filter(Negate(is.null), rows)
    if (length(rows)) do.call(rbind, rows) else data.frame()
}

.qc_prompt_table_block <- function(path, report_file, max_rows = 250L) {
    table <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    if (is.null(table)) return(NULL)
    table <- .qc_prompt_drop_empty(table)
    label <- .ai_qc_report_csv_label(path, report_file)
    dimensions <- paste0(nrow(table), " rows x ", ncol(table), " columns")
    if (!nrow(table) || !ncol(table)) return(paste0("--- ", label, " (", dimensions, ") ---\nNo measured values."))
    if (nrow(table) <= max_rows) {
        return(paste0("--- ", label, " (", dimensions, ") ---\n", .qc_prompt_csv_text(table)))
    }
    summary <- .qc_prompt_numeric_summary(table)
    paste0(
        "--- ", label, " (", dimensions, ") ---\n",
        "The table is too long to paste without obscuring the rest of the evidence. Numerical distribution summary:\n",
        if (nrow(summary)) .qc_prompt_csv_text(summary) else "No finite numeric columns."
    )
}

.build_report_ai_qc_prompt <- function(report_data, report_file, numeric_paths = character(), context = NULL) {
    if (!inherits(report_data, "spectreasy_report_data")) stop("report_data must be created by a Spectreasy report data collector.", call. = FALSE)
    numeric_paths <- unique(as.character(unlist(c(numeric_paths, report_data$qc_metric_paths), recursive = TRUE, use.names = FALSE)))
    numeric_paths <- numeric_paths[!is.na(numeric_paths) & nzchar(numeric_paths) & file.exists(numeric_paths) & !dir.exists(numeric_paths)]
    numeric_paths <- numeric_paths[grepl("\\.csv$", numeric_paths, ignore.case = TRUE)]
    large_tables <- c("reference_spectra.csv", "af_bank_spectra.csv", "control_spectrum_variability.csv")
    embedded <- numeric_paths[!tolower(basename(numeric_paths)) %in% large_tables]
    embedded <- embedded[order(basename(embedded), method = "radix")]
    blocks <- Filter(Negate(is.null), lapply(embedded, .qc_prompt_table_block, report_file = report_file))
    inventory <- lapply(numeric_paths[tolower(basename(numeric_paths)) %in% large_tables], function(path) {
        header <- tryCatch(utils::read.csv(path, nrows = 1L, check.names = FALSE), error = function(e) NULL)
        rows <- tryCatch(length(readLines(path, warn = FALSE)) - 1L, error = function(e) NA_integer_)
        paste0(.ai_qc_report_csv_label(path, report_file), ": ", rows, " data rows x ", if (is.null(header)) NA_integer_ else ncol(header), " columns; retained as CSV and represented here by derived summaries.")
    })
    metadata <- c(
        paste0("Report type: ", report_data$report_type),
        paste0("Project: ", basename(normalizePath(report_data$project_path, mustWork = FALSE))),
        paste0("Generated: ", format(report_data$created_at, "%Y-%m-%d %H:%M:%S %Z")),
        paste0("Spectreasy: ", report_data$version),
        paste0("Unmixing method: ", report_data$unmixing_method),
        paste0("Cytometer: ", report_data$cytometer %||% "not recorded"),
        paste0(names(report_data$counts), ": ", unlist(report_data$counts, use.names = FALSE))
    )
    runtime_notes <- unique(as.character(report_data$warnings %||% character()))
    runtime_notes <- runtime_notes[!is.na(runtime_notes) & nzchar(trimws(runtime_notes))]
    definitions <- c(
        "control_signal_metrics.csv: robust_separation = (positive peak median - negative peak median) / (2 x negative peak MAD); midpoint_overlap_fraction = fraction of positive and negative peak-channel events falling on the opposite side of their median midpoint; spectral_cosine_* = cosine similarity of background-subtracted, event-normalized spectra to the final reference; upper_boundary_fraction = observed fraction at the acquisition maximum recorded in detector metadata.",
        "control_spectrum_variability.csv: q10, median, and q90 are detector-wise quantiles across background-subtracted, event-normalized positive-control spectra; retained as a separate CSV to avoid bloating this prompt.",
        "af_bank_summary.csv: effective_rank is the numerical matrix rank; pairwise values are cosine similarities between AF basis spectra.",
        "af_band_usage.csv: event counts and fractions assigned to each AF basis, grouped by sample.",
        "negative_population_spread.csv, residual, reconstruction, and similarity tables retain the definitions used by the standard Spectreasy QC report."
    )
    paste(
        "SPECTREASY NUMERICAL QC INTERPRETATION PROMPT", "",
        "You are assisting a scientist with interpretation of measured spectral-flow cytometry QC output. Separate observations from interpretation. Do not invent thresholds, categorical quality labels, or unavailable measurements. Do not infer values from plots that are not included. Explain uncertainty and alternative explanations. Compare residual quantities only when their metric definitions and unmixing methods match. Spectral similarity describes panel complexity and is not by itself evidence of a defective control. Do not recommend uploading raw FCS events.", "",
        "Return these sections: Measurement overview; Control and gating observations; Reference matrix and spectral overlap; Autofluorescence observations; Sample and residual observations; Plausible technical explanations; Additional measurements that would reduce uncertainty.", "",
        "RUN AND METHOD", paste(metadata, collapse = "\n"), "",
        "ANALYST CONTEXT", context %||% "No additional analyst context was supplied.", "",
        "METRIC DEFINITIONS", paste(definitions, collapse = "\n"), "",
        "NUMERICAL TABLES", if (length(blocks)) paste(blocks, collapse = "\n\n") else "No machine-readable QC tables were available.", "",
        "LARGE TABLES RETAINED OUTSIDE THE PROMPT", if (length(inventory)) paste(unlist(inventory), collapse = "\n") else "None.", "",
        "RUNTIME NOTES", if (length(runtime_notes)) paste(runtime_notes, collapse = "\n") else "None recorded.", "",
        "LIMITATIONS", "The prompt contains summaries and persisted QC tables, not raw events. Biological interpretation requires panel design, sample preparation, gating, instrument settings, and experimental context. Associations within one run are descriptive and do not establish causation.",
        sep = "\n"
    )
}

.export_report_ai_qc <- function(report_data, report_file, numeric_paths = character(), context = NULL) {
    if (!inherits(report_data, "spectreasy_report_data")) return(NULL)
    report_file <- normalizePath(as.character(report_file)[[1]], mustWork = FALSE)
    prompt <- .build_report_ai_qc_prompt(report_data, report_file, numeric_paths = numeric_paths, context = context)
    prompt_path <- .ai_qc_report_companion_path(report_file)
    .ai_qc_write_text(prompt, prompt_path)
    list(
        prompt = prompt_path,
        numeric_sources = unique(as.character(unlist(c(numeric_paths, report_data$qc_metric_paths), recursive = TRUE, use.names = FALSE))),
        bytes = unname(file.info(prompt_path)$size), estimated_tokens = ceiling(nchar(prompt) / 4)
    )
}
