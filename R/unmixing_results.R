.next_safe_output_dir <- function(path) {
    if (!dir.exists(path)) {
        return(path)
    }

    i <- 2L
    repeat {
        candidate <- paste0(path, "_", i)
        if (!dir.exists(candidate)) {
            return(candidate)
        }
        i <- i + 1L
    }
}

.unmix_method_file_label <- function(method) {
    method <- .normalize_unmix_method(method)
    gsub("[^A-Za-z0-9]+", "", method)
}

.reference_af_band_count <- function(M) {
    M <- .as_reference_matrix(M, "M")
    sum(grepl("^AF($|_)", rownames(M), ignore.case = TRUE))
}

.unmixed_fcs_filename <- function(sample_name, method, M) {
    sample_name <- trimws(as.character(sample_name)[1])
    sample_name <- gsub("[[:cntrl:]/\\\\:]", "_", sample_name)
    sample_name <- gsub("[. ]+$", "", sample_name)
    if (!nzchar(sample_name) || sample_name %in% c(".", "..")) sample_name <- "sample"
    sample_name <- substr(sample_name, 1L, 180L)
    paste0(
        sample_name,
        "_",
        .unmix_method_file_label(method),
        "-",
        .reference_af_band_count(M),
        "AF.fcs"
    )
}

.normalize_unmix_chunk_size <- function(chunk_size) {
    .normalize_positive_integer(chunk_size, "chunk_size", allow_null = TRUE)
}

.normalize_unmix_plot_n_events <- function(plot_n_events) {
    if (is.null(plot_n_events)) {
        return(NULL)
    }
    .normalize_positive_integer(plot_n_events, "plot_n_events")
}

.unmix_chunk_indices <- function(n_events, chunk_size) {
    if (!is.finite(chunk_size) || n_events <= chunk_size) {
        return(list(seq_len(n_events)))
    }
    starts <- seq.int(1L, n_events, by = chunk_size)
    lapply(starts, function(start) seq.int(start, min(n_events, start + chunk_size - 1L)))
}

.unmix_sample_keep_indices <- function(n_events, plot_n_events = NULL) {
    if (is.null(plot_n_events) || n_events <= plot_n_events) {
        return(seq_len(n_events))
    }
    sort(sample.int(n_events, plot_n_events))
}

.merge_unmix_variant_info <- function(infos) {
    infos <- infos[!vapply(infos, is.null, logical(1))]
    if (length(infos) == 0) {
        return(NULL)
    }
    changed_events <- sum(vapply(infos, function(x) {
        val <- x$changed_events
        if (is.null(val) || !is.finite(val)) 0 else as.integer(val)
    }, integer(1)))
    event_counts <- vapply(infos, function(x) {
        val <- x$event_count
        if (!is.null(val) && length(val) == 1L && is.finite(val) && val >= 0) {
            return(as.integer(val))
        }
        if (is.list(x$selected)) return(length(x$selected))
        fraction <- x$changed_fraction
        changed <- x$changed_events
        if (!is.null(fraction) && is.finite(fraction) && fraction > 0 &&
            !is.null(changed) && is.finite(changed)) {
            return(as.integer(round(changed / fraction)))
        }
        0L
    }, integer(1))
    event_count <- sum(event_counts)
    changed_fraction <- if (event_count > 0L) {
        changed_events / event_count
    } else {
        mean(vapply(infos, function(x) {
            val <- x$changed_fraction
            if (is.null(val) || !is.finite(val)) 0 else as.numeric(val)
        }, numeric(1)), na.rm = TRUE)
    }
    selected <- do.call(c, lapply(seq_along(infos), function(index) {
        value <- infos[[index]]$selected
        if (is.list(value)) return(value)
        vector("list", event_counts[[index]])
    }))
    out <- infos[[1]]
    out$changed_events <- changed_events
    out$changed_fraction <- changed_fraction
    out$event_count <- event_count
    out$selected <- selected
    out
}

.combine_unmix_result_chunks <- function(data_chunks,
                                         residual_chunks,
                                         method,
                                         M,
                                         variant_infos = list(),
                                         spectreasy_decoder_weights = NULL) {
    data_chunks <- data_chunks[!vapply(data_chunks, is.null, logical(1))]
    residual_chunks <- residual_chunks[!vapply(residual_chunks, is.null, logical(1))]
    out <- list(
        data = if (length(data_chunks) > 0) do.call(rbind, data_chunks) else data.frame(),
        residuals = if (length(residual_chunks) > 0) do.call(rbind, residual_chunks) else NULL
    )
    variant_info <- .merge_unmix_variant_info(variant_infos)
    if (!is.null(variant_info)) {
        out$spectral_variant_info <- variant_info
    }
    if (!is.null(spectreasy_decoder_weights)) {
        out$spectreasy_decoder_weights <- spectreasy_decoder_weights
    }
    attr(out, "method") <- method
    attr(out, "reference_matrix") <- M
    out
}

.as_resolved_unmixed_fcs_dir <- function(path) {
    structure(as.character(path)[1], class = c("spectreasy_resolved_unmixed_fcs_dir", "character"))
}

.unmix_samples_output_paths <- function(output_dir) {
    is_resolved <- inherits(output_dir, "spectreasy_resolved_unmixed_fcs_dir")
    output_dir <- as.character(output_dir)[1]
    if (is_resolved) {
        sample_dir <- dirname(output_dir)
        return(list(
            sample_dir = sample_dir,
            unmixed_dir = unclass(output_dir),
            qc_samples_dir = file.path(sample_dir, "qc_samples")
        ))
    }
    sample_dir <- file.path(output_dir, "unmix_samples")
    list(
        sample_dir = sample_dir,
        unmixed_dir = file.path(sample_dir, "unmixed_fcs"),
        qc_samples_dir = file.path(sample_dir, "qc_samples")
    )
}

.default_unmix_samples_report_dir <- function(output_dir) {
    .unmix_samples_output_paths(output_dir)$qc_samples_dir
}

.as_unmixed_results_data_frame <- function(x, arg_name = "x") {
    if (!is.list(x) || length(x) == 0) {
        stop(arg_name, " must be a non-empty list returned by unmix_samples().")
    }

    sample_names <- names(x)
    if (is.null(sample_names)) {
        sample_names <- rep("", length(x))
    }

    data_parts <- vector("list", length(x))
    for (i in seq_along(x)) {
        res_obj <- x[[i]]
        if (!is.list(res_obj) || !is.data.frame(res_obj$data)) {
            stop(arg_name, " must contain one data frame per sample in $data.")
        }

        data_df <- as.data.frame(res_obj$data, stringsAsFactors = FALSE, check.names = FALSE)
        if (!("File" %in% colnames(data_df))) {
            sample_name <- if (nzchar(sample_names[[i]])) sample_names[[i]] else paste0("sample_", i)
            data_df$File <- sample_name
        }
        data_parts[[i]] <- data_df
    }

    do.call(rbind, data_parts)
}

#' @export
as.data.frame.spectreasy_unmixed_results <- function(x, row.names = NULL, optional = FALSE, ...) {
    .as_unmixed_results_data_frame(x, arg_name = "x")
}

.unmixed_results_to_flowset <- function(results) {
    if (!is.list(results) || length(results) == 0) {
        stop("results must be a non-empty list returned by unmix_samples().")
    }

    sample_names <- names(results)
    if (is.null(sample_names) || any(sample_names == "")) {
        sample_names <- paste0("sample_", seq_along(results))
        names(results) <- sample_names
    }

    fs_list <- vector("list", length(results))
    names(fs_list) <- sample_names
    residuals_attr <- vector("list", length(results))
    names(residuals_attr) <- sample_names

    for (i in seq_along(results)) {
        sample_name <- sample_names[[i]]
        res_obj <- results[[i]]
        if (!is.list(res_obj) || !is.data.frame(res_obj$data)) {
            stop("Each result element must contain a data frame in $data.")
        }

        expr_cols <- setdiff(colnames(res_obj$data), "File")
        expr_mat <- as.matrix(res_obj$data[, expr_cols, drop = FALSE])
        storage.mode(expr_mat) <- "numeric"
        ff <- flowCore::flowFrame(expr_mat)

        pd <- flowCore::pData(flowCore::parameters(ff))
        if (all(c("name", "desc") %in% colnames(pd))) {
            pd$desc <- as.character(pd$name)
            flowCore::parameters(ff) <- methods::new("AnnotatedDataFrame", data = pd)
        }

        fs_list[[i]] <- ff
        residuals_attr[[i]] <- if ("residuals" %in% names(res_obj)) res_obj$residuals else NULL
    }

    pheno <- methods::new(
        "AnnotatedDataFrame",
        data = data.frame(sample_id = sample_names, row.names = sample_names, stringsAsFactors = FALSE)
    )
    fs <- flowCore::flowSet(fs_list, phenoData = pheno)
    attr(fs, "spectreasy_residuals") <- residuals_attr
    fs
}

.unmixed_results_to_sce <- function(results, sample_entries) {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE) ||
        !requireNamespace("SummarizedExperiment", quietly = TRUE) ||
        !requireNamespace("S4Vectors", quietly = TRUE)) {
        stop(
            "Packages 'SingleCellExperiment', 'SummarizedExperiment', and 'S4Vectors' are required for return_type = 'SingleCellExperiment'.",
            call. = FALSE
        )
    }
    if (!is.list(results) || length(results) == 0) {
        stop("results must be a non-empty list returned by unmix_samples().")
    }

    sample_names <- names(results)
    if (is.null(sample_names) || any(sample_names == "")) {
        sample_names <- vapply(sample_entries, function(x) x$sample_name, character(1))
    }

    expr_cols <- unique(unlist(lapply(results, function(res_obj) setdiff(colnames(res_obj$data), "File"))))
    if (length(expr_cols) == 0) {
        stop("No numeric result columns available to construct SingleCellExperiment output.")
    }

    expr_parts <- list()
    coldata_parts <- list()
    residual_parts <- list()

    for (i in seq_along(results)) {
        res_obj <- results[[i]]
        entry <- sample_entries[[i]]
        sample_name <- if (!is.null(names(results)[i]) && nzchar(names(results)[i])) names(results)[i] else entry$sample_name
        data_df <- res_obj$data

        missing_cols <- setdiff(expr_cols, colnames(data_df))
        if (length(missing_cols) > 0) {
            for (mc in missing_cols) data_df[[mc]] <- NA_real_
        }
        data_df <- data_df[, expr_cols, drop = FALSE]
        expr_mat <- t(as.matrix(data_df))
        storage.mode(expr_mat) <- "numeric"

        source_cell_ids <- attr(res_obj, "source_cell_ids", exact = TRUE)
        if (is.null(source_cell_ids) && !is.null(entry$cell_ids) && length(entry$cell_ids) == nrow(res_obj$data)) {
            source_cell_ids <- entry$cell_ids
        }
        if (!is.null(source_cell_ids) && length(source_cell_ids) == nrow(res_obj$data)) {
            cell_ids <- paste0(sample_name, "__", make.unique(as.character(source_cell_ids)))
        } else {
            cell_ids <- paste0(sample_name, "_", seq_len(nrow(res_obj$data)))
        }
        colnames(expr_mat) <- cell_ids
        expr_parts[[i]] <- expr_mat

        source_name <- if (!is.null(entry$source_name) && nzchar(entry$source_name)) entry$source_name else sample_name
        coldata_parts[[i]] <- data.frame(
            sample_id = rep(sample_name, length(cell_ids)),
            source_name = rep(source_name, length(cell_ids)),
            row.names = cell_ids,
            stringsAsFactors = FALSE
        )

        residuals <- if ("residuals" %in% names(res_obj)) res_obj$residuals else NULL
        if (!is.null(residuals)) {
            residual_mat <- t(as.matrix(residuals))
            storage.mode(residual_mat) <- "numeric"
            colnames(residual_mat) <- cell_ids
            residual_parts[[i]] <- residual_mat
        } else {
            residual_parts[[i]] <- NULL
        }
    }

    assay_mat <- do.call(cbind, expr_parts)
    col_data <- do.call(rbind, coldata_parts)
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(unmixed = assay_mat),
        colData = S4Vectors::DataFrame(col_data)
    )

    non_null_residuals <- residual_parts[!vapply(residual_parts, is.null, logical(1))]
    if (length(non_null_residuals) > 0) {
        residual_features <- unique(unlist(lapply(non_null_residuals, rownames)))
        residual_aligned <- lapply(seq_along(residual_parts), function(i) {
            mat <- residual_parts[[i]]
            cell_ids <- colnames(expr_parts[[i]])
            out <- matrix(
                NA_real_,
                nrow = length(residual_features),
                ncol = length(cell_ids),
                dimnames = list(residual_features, cell_ids)
            )
            if (!is.null(mat)) {
                out[rownames(mat), colnames(mat)] <- mat
            }
            out
        })
        residual_mat <- do.call(cbind, residual_aligned)
        residual_sce <- SingleCellExperiment::SingleCellExperiment(
            assays = list(residuals = residual_mat)
        )
        SingleCellExperiment::altExp(sce, "detector_residuals") <- residual_sce
    }

    sce
}
