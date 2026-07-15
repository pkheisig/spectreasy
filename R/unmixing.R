# Internal helpers for sample unmixing.
.prepare_unmix_samples_input_from_flowset <- function(sample_input) {
    sample_names <- flowCore::sampleNames(sample_input)
    if (length(sample_names) == 0) {
        sample_names <- paste0("sample_", seq_len(length(sample_input)))
    }
    sample_keys <- make.unique(tools::file_path_sans_ext(basename(sample_names)))
    lapply(seq_along(sample_names), function(i) {
        list(
            sample_name = sample_keys[[i]],
            source_name = sample_names[[i]],
            flow_frame = sample_input[[i]],
            cell_ids = NULL
        )
    })
}

.select_sce_unmix_assay <- function(sample_input) {
    assay_names <- SummarizedExperiment::assayNames(sample_input)
    if (length(assay_names) == 0) {
        stop("SingleCellExperiment input must contain at least one assay.")
    }
    preferred <- c("exprs", "counts", "normcounts", "logcounts")
    hit <- intersect(preferred, assay_names)
    if (length(hit) > 0) hit[[1]] else assay_names[[1]]
}

.resolve_sce_sample_column <- function(sample_input) {
    col_data <- as.data.frame(SummarizedExperiment::colData(sample_input), stringsAsFactors = FALSE)
    candidate_cols <- c("sample_id", "sample", "sample_name", "File", "file")
    hit <- intersect(candidate_cols, colnames(col_data))
    if (length(hit) == 0) {
        return(rep("sample_1", ncol(sample_input)))
    }
    sample_ids <- trimws(as.character(col_data[[hit[[1]]]]))
    sample_ids[is.na(sample_ids) | sample_ids == ""] <- "sample_1"
    sample_ids
}

.prepare_unmix_samples_input_from_sce <- function(sample_input) {
    assay_name <- .select_sce_unmix_assay(sample_input)
    assay_mat <- SummarizedExperiment::assay(sample_input, assay_name)
    if (is.null(rownames(assay_mat)) || any(rownames(assay_mat) == "")) {
        stop(
            "SingleCellExperiment input must have non-empty feature names in rownames() ",
            "so detector and acquisition channels can be matched."
        )
    }

    sample_ids <- .resolve_sce_sample_column(sample_input)
    sample_levels <- unique(sample_ids)
    sample_keys <- make.unique(tools::file_path_sans_ext(basename(sample_levels)))
    cell_names <- colnames(sample_input)
    if (is.null(cell_names) || any(cell_names == "")) {
        cell_names <- paste0("cell_", seq_len(ncol(sample_input)))
    }

    lapply(seq_along(sample_levels), function(i) {
        sample_id <- sample_levels[[i]]
        idx <- which(sample_ids == sample_id)
        expr_mat <- t(as.matrix(assay_mat[, idx, drop = FALSE]))
        storage.mode(expr_mat) <- "numeric"
        colnames(expr_mat) <- rownames(assay_mat)
        ff <- flowCore::flowFrame(expr_mat)
        list(
            sample_name = sample_keys[[i]],
            source_name = sample_id,
            flow_frame = ff,
            cell_ids = cell_names[idx]
        )
    })
}

.prepare_unmix_samples_input_from_directory <- function(sample_input) {
    if (!dir.exists(sample_input)) {
        .spectreasy_stop_missing_directory(sample_input, label = "sample_dir")
    }
    sample_dir_resolved <- normalizePath(sample_input, mustWork = FALSE)
    fcs_files <- list.files(sample_dir_resolved, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
    if (length(fcs_files) == 0) {
        .spectreasy_stop_empty_fcs_directory(sample_dir_resolved, label = "sample_dir")
    }

    sample_keys <- make.unique(tools::file_path_sans_ext(basename(fcs_files)))
    lapply(seq_along(fcs_files), function(i) {
        f <- fcs_files[[i]]
        list(
            sample_name = sample_keys[[i]],
            source_name = basename(f),
            file_path = f,
            flow_frame = NULL,
            cell_ids = NULL
        )
    })
}

.prepare_unmix_samples_input_from_files <- function(fcs_files) {
    fcs_files <- as.character(fcs_files)
    if (length(fcs_files) == 0L) {
        stop("sample_dir must contain at least one FCS file path.", call. = FALSE)
    }
    missing_files <- fcs_files[is.na(fcs_files) | !file.exists(fcs_files) | dir.exists(fcs_files)]
    if (length(missing_files) > 0) {
        stop(
            "sample_dir contains missing or invalid FCS file paths: ",
            paste(missing_files, collapse = ", "),
            call. = FALSE
        )
    }
    invalid_ext <- fcs_files[!grepl("\\.fcs$", fcs_files, ignore.case = TRUE)]
    if (length(invalid_ext) > 0) {
        stop("sample_dir file paths must point to .fcs files: ", paste(invalid_ext, collapse = ", "), call. = FALSE)
    }

    normalized_files <- normalizePath(fcs_files, mustWork = TRUE)
    duplicate_files <- unique(normalized_files[duplicated(normalized_files)])
    if (length(duplicate_files) > 0L) {
        stop("sample_dir contains duplicate FCS file paths: ", paste(duplicate_files, collapse = ", "), call. = FALSE)
    }
    sample_keys <- make.unique(tools::file_path_sans_ext(basename(fcs_files)))
    lapply(seq_along(fcs_files), function(i) {
        f <- fcs_files[[i]]
        list(
            sample_name = sample_keys[[i]],
            source_name = basename(f),
            file_path = normalized_files[[i]],
            flow_frame = NULL,
            cell_ids = NULL
        )
    })
}

.prepare_unmix_samples_input <- function(sample_input) {
    if (inherits(sample_input, "flowSet")) {
        return(.prepare_unmix_samples_input_from_flowset(sample_input))
    }

    if (inherits(sample_input, "SingleCellExperiment")) {
        if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
            stop(
                "Package 'SummarizedExperiment' is required to unmix SingleCellExperiment inputs.",
                call. = FALSE
            )
        }
        return(.prepare_unmix_samples_input_from_sce(sample_input))
    }

    if (is.character(sample_input) &&
        (length(sample_input) > 1L ||
            (length(sample_input) == 1L && !is.na(sample_input) && grepl("\\.fcs$", sample_input, ignore.case = TRUE)))) {
        return(.prepare_unmix_samples_input_from_files(sample_input))
    }

    if (!is.character(sample_input) || length(sample_input) != 1 || is.na(sample_input)) {
        stop("sample_dir must be a directory path, one or more FCS file paths, a flowCore::flowSet, or a SingleCellExperiment.", call. = FALSE)
    }

    .prepare_unmix_samples_input_from_directory(sample_input)
}

.aggregate_af_columns_for_output <- function(data, fluorophore_source) {
    af_source <- fluorophore_source[grepl("^AF($|_)", fluorophore_source, ignore.case = TRUE)]
    af_source <- intersect(af_source, colnames(data))
    if ("AF" %in% fluorophore_source && length(af_source) > 1L) {
        data[["AF"]] <- rowSums(as.matrix(data[, af_source, drop = FALSE]))
    }
    data
}

.read_unmixing_matrix_csv <- function(path) {
    if (!file.exists(path)) .spectreasy_stop_missing_file(path, label = "unmixing_matrix_file")
    df <- tryCatch(
        utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
        error = function(e) {
            stop("Could not read unmixing_matrix_file: ", path, call. = FALSE)
        }
    )
    if (ncol(df) < 2) stop("Invalid unmixing_matrix_file: expected at least 2 columns in ", path, call. = FALSE)

    first_name <- tolower(trimws(colnames(df)[1]))
    first_col <- df[[1]]
    has_marker_col <- first_name %in% c("marker", "fluorophore", "file") || !is.numeric(first_col)
    if (has_marker_col) {
        marker_names <- trimws(as.character(first_col))
        mat_df <- df[, -1, drop = FALSE]
    } else {
        marker_names <- rownames(df)
        mat_df <- df
    }
    W_mat <- as.matrix(mat_df)
    storage.mode(W_mat) <- "numeric"
    if (is.null(marker_names) || all(marker_names == "")) {
        marker_names <- paste0("marker_", seq_len(nrow(W_mat)))
    }
    rownames(W_mat) <- marker_names
    colnames(W_mat) <- colnames(mat_df)
    W_mat
}

.default_unmixing_matrix_file <- function() {
    file.path("spectreasy_outputs", "unmix_controls", "scc_reference_matrix.csv")
}

.candidate_unmixing_matrix_files <- function(unmixing_matrix_file, output_dir = NULL) {
    candidates <- character()

    if (!is.null(unmixing_matrix_file) && length(unmixing_matrix_file) > 0) {
        matrix_file <- as.character(unmixing_matrix_file)[1]
        if (!is.na(matrix_file) && nzchar(trimws(matrix_file))) {
            candidates <- c(candidates, matrix_file)
            if (identical(basename(matrix_file), "scc_reference_matrix.csv")) {
                matrix_dir <- dirname(matrix_file)
                candidates <- c(candidates, file.path(matrix_dir, "unmixed_fcs", "scc_reference_matrix.csv"))
                if (identical(basename(normalizePath(matrix_dir, mustWork = FALSE)), "unmixed_fcs")) {
                    candidates <- c(candidates, file.path(dirname(matrix_dir), "scc_reference_matrix.csv"))
                }
            }
        }
    }

    if (!is.null(output_dir) && length(output_dir) > 0) {
        sample_output_dir <- as.character(output_dir)[1]
        if (!is.na(sample_output_dir) && nzchar(trimws(sample_output_dir))) {
            sample_output_dir <- normalizePath(sample_output_dir, mustWork = FALSE)
            project_dir <- if (inherits(output_dir, "spectreasy_resolved_unmixed_fcs_dir")) {
                dirname(dirname(sample_output_dir))
            } else {
                sample_output_dir
            }
            candidates <- c(
                candidates,
                file.path(project_dir, "unmix_controls", "scc_reference_matrix.csv")
            )
        }
    }

    unique(candidates[!is.na(candidates) & nzchar(trimws(candidates))])
}

.resolve_unmixing_matrix_file_for_samples <- function(unmixing_matrix_file, output_dir = NULL) {
    candidates <- .candidate_unmixing_matrix_files(
        unmixing_matrix_file = unmixing_matrix_file,
        output_dir = output_dir
    )
    existing <- candidates[file.exists(candidates)]
    if (length(existing) > 0) {
        return(existing[[1]])
    }
    unmixing_matrix_file
}

.resolve_detector_noise_file_for_unmixing <- function(unmixing_matrix_file,
                                                     detector_noise_file = NULL) {
    if (!is.null(detector_noise_file) && file.exists(detector_noise_file)) {
        return(detector_noise_file)
    }
    if (!is.null(unmixing_matrix_file) && file.exists(unmixing_matrix_file)) {
        sibling <- file.path(dirname(unmixing_matrix_file), "scc_detector_noise.csv")
        if (file.exists(sibling)) {
            return(sibling)
        }
    }
    detector_noise_file
}

.load_detector_noise_for_unmixing <- function(M,
                                             detector_noise_file = NULL,
                                             unmixing_matrix_file = NULL,
                                             scc_dir = NULL) {
    M <- .as_reference_matrix(M, "M")
    if (!is.null(attr(M, "detector_noise"))) {
        return(M)
    }

    detector_noise_file <- .resolve_detector_noise_file_for_unmixing(
        unmixing_matrix_file = unmixing_matrix_file,
        detector_noise_file = detector_noise_file
    )
    if (!is.null(detector_noise_file) && file.exists(detector_noise_file)) {
        detector_noise <- utils::read.csv(detector_noise_file, stringsAsFactors = FALSE, check.names = FALSE)
        return(.attach_detector_noise(M, detector_noise, source = detector_noise_file))
    }

    if (!is.null(scc_dir) && dir.exists(scc_dir)) {
        return(.attach_estimated_wls_detector_noise(
            M = M,
            scc_dir = scc_dir,
            fallback = .default_wls_background_noise(),
            warn = TRUE
        ))
    }

    M
}

.resolve_output_marker_map <- function(fluorophore_names, control_file = NULL, control_df = NULL) {
    fluorophore_names <- trimws(as.character(fluorophore_names))
    labels <- stats::setNames(fluorophore_names, fluorophore_names)
    if (length(fluorophore_names) == 0) return(labels)
    mapping_supplied <- !is.null(control_file) || !is.null(control_df)

    if (!is.null(control_file)) {
        control_df <- tryCatch(
            utils::read.csv(control_file, stringsAsFactors = FALSE, check.names = FALSE),
            error = function(e) {
                stop(
                    "Could not read control_file for output FCS labels: ",
                    conditionMessage(e),
                    call. = FALSE
                )
            }
        )
    }
    if (is.null(control_df) || !is.data.frame(control_df) || nrow(control_df) == 0) {
        if (mapping_supplied) {
            stop("The control mapping used for output FCS labels is empty or invalid.", call. = FALSE)
        }
        return(labels)
    }

    lower_names <- tolower(colnames(control_df))
    fluor_idx <- match("fluorophore", lower_names)
    marker_idx <- match("marker", lower_names)
    if (is.na(fluor_idx) || is.na(marker_idx)) {
        stop(
            "The control mapping used for output FCS labels must contain fluorophore and marker columns.",
            call. = FALSE
        )
    }

    fluor_vals <- trimws(as.character(control_df[[fluor_idx]]))
    marker_vals <- trimws(as.character(control_df[[marker_idx]]))
    keep <- !is.na(fluor_vals) & fluor_vals != "" & !is.na(marker_vals) & marker_vals != ""
    if (!any(keep)) return(labels)

    key <- tolower(fluor_vals[keep])
    val <- marker_vals[keep]
    map_ci <- stats::setNames(val, key)
    map_ci <- map_ci[!duplicated(names(map_ci))]

    for (nm in names(labels)) {
        k <- tolower(trimws(nm))
        if (k %in% names(map_ci) && nzchar(map_ci[[k]])) {
            labels[[nm]] <- map_ci[[k]]
        }
    }
    labels
}

.resolve_unmix_output_control_context <- function(M,
                                                  control_file = NULL,
                                                  unmixing_matrix_file = NULL) {
    if (!is.null(control_file)) {
        return(list(control_file = control_file, control_df = NULL))
    }

    control_df <- attr(M, "spectreasy_control_df", exact = TRUE)
    if (is.data.frame(control_df) && nrow(control_df) > 0L) {
        return(list(control_file = NULL, control_df = control_df))
    }

    matrix_control_file <- attr(M, "spectreasy_control_file", exact = TRUE)
    if (is.character(matrix_control_file) && length(matrix_control_file) == 1L &&
        !is.na(matrix_control_file) && file.exists(matrix_control_file)) {
        return(list(control_file = matrix_control_file, control_df = NULL))
    }

    if (!is.null(unmixing_matrix_file) && file.exists(unmixing_matrix_file)) {
        mapping_used_file <- file.path(dirname(unmixing_matrix_file), "fcs_mapping_used.csv")
        if (file.exists(mapping_used_file)) {
            return(list(control_file = mapping_used_file, control_df = NULL))
        }
    }

    list(control_file = NULL, control_df = NULL)
}

.apply_output_fcs_feature_labels <- function(target_ff, source_ff, fluorophore_cols, marker_label_map) {
    pd_new <- flowCore::pData(flowCore::parameters(target_ff))
    if (!all(c("name", "desc") %in% colnames(pd_new))) {
        return(target_ff)
    }

    param_names <- as.character(pd_new$name)
    primary_vals <- param_names
    desc_vals <- param_names

    pd_src <- tryCatch(flowCore::pData(flowCore::parameters(source_ff)), error = function(e) NULL)
    if (!is.null(pd_src) && all(c("name", "desc") %in% colnames(pd_src))) {
        src_name <- as.character(pd_src$name)
        src_desc <- as.character(pd_src$desc)
        src_desc[is.na(src_desc)] <- ""
        src_map <- stats::setNames(src_desc, src_name)
        matched <- param_names %in% names(src_map)
        desc_vals[matched] <- ifelse(
            nzchar(src_map[param_names[matched]]),
            src_map[param_names[matched]],
            param_names[matched]
        )
    }

    fluorophore_cols <- as.character(fluorophore_cols)
    for (fluorophore in fluorophore_cols) {
        idx <- which(param_names == fluorophore)
        if (length(idx) == 0) next
        marker <- if (fluorophore %in% names(marker_label_map)) {
            marker_label_map[[fluorophore]]
        } else {
            fluorophore
        }
        if (is.na(marker) || !nzchar(trimws(marker))) marker <- fluorophore
        primary_vals[idx] <- as.character(marker)
        desc_vals[idx] <- fluorophore
    }

    primary_vals <- make.unique(primary_vals, sep = "_")
    output_exprs <- flowCore::exprs(target_ff)
    colnames(output_exprs) <- primary_vals
    output_ff <- flowCore::flowFrame(output_exprs)
    output_pd <- flowCore::pData(flowCore::parameters(output_ff))
    output_pd$desc <- desc_vals
    flowCore::parameters(output_ff) <- methods::new("AnnotatedDataFrame", data = output_pd)
    output_ff
}

.next_safe_output_path <- function(path) {
    if (!file.exists(path)) {
        return(path)
    }

    dir_path <- dirname(path)
    ext <- tools::file_ext(path)
    ext_suffix <- if (nzchar(ext)) paste0(".", ext) else ""
    stem <- tools::file_path_sans_ext(basename(path))

    i <- 2L
    repeat {
        candidate <- file.path(dir_path, paste0(stem, "_", i, ext_suffix))
        if (!file.exists(candidate)) {
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
    changed_fraction <- mean(vapply(infos, function(x) {
        val <- x$changed_fraction
        if (is.null(val) || !is.finite(val)) 0 else as.numeric(val)
    }, numeric(1)), na.rm = TRUE)
    out <- infos[[1]]
    out$changed_events <- changed_events
    out$changed_fraction <- changed_fraction
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

        if (!is.null(entry$cell_ids) && length(entry$cell_ids) == nrow(res_obj$data)) {
            cell_ids <- paste0(sample_name, "__", make.unique(as.character(entry$cell_ids)))
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

#' Unmix Experimental Samples
#' 
#' @param sample_dir Directory containing experimental FCS files, a character
#'   vector of FCS file paths, a `flowCore::flowSet`, or a
#'   `SingleCellExperiment` for in-memory workflows.
#' @param samples_dir Alias for `sample_dir`, accepted for workflows and scripts
#'   that use the plural directory name.
#' @param M Optional reference matrix (Markers x Detectors). If supplied,
#'   unmixing is computed dynamically using this matrix. If not supplied,
#'   it is loaded from the CSV path provided in `unmixing_matrix_file`.
#' @param unmixing_matrix_file Optional CSV path to a saved reference matrix.
#'   Used when `M` is not supplied. By default this points to the reference matrix
#'   produced by [unmix_controls()] (`"scc_reference_matrix.csv"`).
#' @param detector_noise_file Optional CSV path to detector-specific WLS noise
#'   floors, as written by [unmix_controls()] (`"scc_detector_noise.csv"`). If
#'   omitted, `unmix_samples()` first looks beside `unmixing_matrix_file`, then
#'   falls back to the built-in scalar noise floor.
#' @param control_file Optional control mapping CSV used to assign marker names
#'   as primary FCS parameter names and fluorophores as secondary names. An
#'   explicit value takes precedence. When omitted, the exact mapping carried by
#'   an in-memory matrix returned from [unmix_controls()] or the
#'   `fcs_mapping_used.csv` saved beside `unmixing_matrix_file` is used.
#' @param unmixing_method Unmixing method (`"WLS"`, `"RWLS"`, `"OLS"`,
#'   `"NNLS"`, `"AutoSpectral"`, or `"Spectreasy"`). `AutoSpectral` uses
#'   per-event AF assignment with marker + selected-AF OLS, plus SCC-derived
#'   spectral-variant optimization when available. `Spectreasy` uses the same
#'   AutoSpectral-style OLS fit, then blends marker abundances with a marker-only
#'   OLS anchor using decoder-projected AF impact weights.
#' @param rwls_max_iter Positive integer; number of robust reweighting
#'   iterations used when `unmixing_method = "RWLS"`. The default, 1, preserves the
#'   historical behavior.
#' @param n_threads Positive integer; number of threads for event-wise
#'   AutoSpectral/Spectreasy AF assignment and NNLS/WLS/RWLS fitting. OLS uses
#'   vectorized matrix operations. The default, 1, keeps explicit event loops
#'   single-threaded.
#' @param spectral_variant_library Optional in-memory AutoSpectral
#'   spectral-variant library. Used only when `unmixing_method = "AutoSpectral"`
#'   or `"Spectreasy"`.
#' @param spectral_variant_library_file Optional RDS path to an AutoSpectral
#'   spectral-variant library. When omitted, `unmix_samples()` looks for
#'   `scc_spectral_variants.rds` beside `unmixing_matrix_file`.
#' @param spectral_variant_top_k Number of best variant candidates to test per
#'   positive fluorophore for AutoSpectral.
#' @param spectral_variant_min_abundance Minimum unmixed abundance for a
#'   fluorophore to be eligible for AutoSpectral variant testing.
#' @param spectral_variant_positive_fraction Additional positivity threshold
#'   as a fraction of the event's strongest fluorophore abundance for
#'   AutoSpectral variant testing.
#' @param spectral_variant_min_improvement Minimum fractional residual
#'   improvement required before accepting a cell-specific AutoSpectral
#'   variant refit.
#' @param spectreasy_weight_quantile Numeric in `[0, 1]`; only accepted when
#'   `unmixing_method = "Spectreasy"`. Controls the quantile of
#'   decoder-projected AF impacts used as the soft-saturation scale for
#'   marker-specific AutoSpectral mixing. The default is `0.9`.
#' @param estimate_af Logical; if `TRUE`, estimate AF signatures directly from
#'   stained sample event-wise WLS residuals, select the best candidate model by
#'   held-out WLS residual score, and append the selected AF rows to the
#'   reference matrix before unmixing. This is intended for workflows where no
#'   unstained cell control is available. Default is `FALSE`.
#' @param output_dir Root output directory. Sample-stage artifacts are written
#'   under `output_dir/unmix_samples`; unmixed FCS files are saved in its
#'   `unmixed_fcs` subdirectory.
#' @param write_fcs Logical; if `TRUE`, write unmixed FCS files to
#'   `output_dir/unmix_samples/unmixed_fcs`.
#'   Defaults to `TRUE` so unmixed FCS files are written unless disabled explicitly.
#' @param save_report Logical; if `TRUE`, write a sample QC report and
#'   sample QC metric CSVs from the in-memory unmixing results without rerunning
#'   unmixing. Defaults to `TRUE`.
#' @param report_format Report format, either `"html"` (default) or `"pdf"`.
#'   Only the selected format is written. Matching is case-insensitive.
#' @param save_qc_plots Logical; if `TRUE`, save QC report plots as PNG files
#'   in `qc_plot_dir` while creating the report.
#' @param qc_plot_dir Directory for sample QC report PNG files when
#'   `save_qc_plots = TRUE`.
#' @param plot_n_events Optional integer; number of events per sample retained
#'   for automatic QC report plots and returned in-memory results. The full
#'   unsampled unmixed data will still be written to FCS files. Defaults to
#'   10000. Set `plot_n_events = NULL` to keep all events in the returned
#'   object and report input.
#' @param chunk_size Integer number of events unmixed at a time per sample.
#'   Defaults to 50000 to reduce peak memory use for large FCS files. Set
#'   `chunk_size = NULL` to process each sample in one chunk.
#' @param seed Optional integer seed for deterministic subsampling.
#' @param return_type Return format: `"list"` (default), `"flowSet"`, or
#'   `"SingleCellExperiment"`. When `"flowSet"`, detector residuals are attached
#'   as `attr(x, "spectreasy_residuals")`. When `"SingleCellExperiment"`, cell-level
#'   unmixed values are returned in assay `"unmixed"`, with detector residuals in
#'   `altExp(x, "detector_residuals")` when available.
#' @param verbose Logical; if `TRUE`, print progress messages while unmixing
#'   each sample.
#' @return Either a named list with one element per sample, a `flowSet`, or a
#'   `SingleCellExperiment` depending on `return_type`. For `return_type = "list"`,
#'   the result has class `spectreasy_unmixed_results`; list elements contain
#'   `data` (unmixed abundances plus retained acquisition parameters) and
#'   `residuals` (detector residual matrix when available, otherwise `NULL`).
#'   The list can be passed directly to `qc_samples(results = ...)`
#'   or coerced with `as.data.frame()`. The return value is provided invisibly to
#'   avoid printing large result objects during interactive or Quarto execution.
#' @examples
#' M_demo <- rbind(
#'   FITC = c(1.00, 0.20, 0.05),
#'   PE = c(0.10, 1.00, 0.20),
#'   APC = c(0.05, 0.15, 1.00)
#' )
#' colnames(M_demo) <- c("B2-A", "YG1-A", "R1-A")
#'
#' simulate_sample <- function(dominant_marker, M, n_cells = 120) {
#'   markers <- rownames(M)
#'   marker_signal <- matrix(rexp(n_cells * length(markers), rate = 8), ncol = length(markers))
#'   colnames(marker_signal) <- markers
#'   marker_signal[, dominant_marker] <- rexp(n_cells, rate = 0.6) + 2
#'   raw_signal <- marker_signal %*% M + matrix(rnorm(n_cells * ncol(M), sd = 0.03), ncol = ncol(M))
#'   exprs_mat <- cbind(
#'     raw_signal,
#'     Time = seq_len(n_cells),
#'     "FSC-A" = rnorm(n_cells, mean = 90000, sd = 7000),
#'     "SSC-A" = rnorm(n_cells, mean = 45000, sd = 5000)
#'   )
#'   colnames(exprs_mat)[seq_len(ncol(M))] <- colnames(M)
#'   flowCore::flowFrame(exprs_mat)
#' }
#'
#' toy_fs <- flowCore::flowSet(list(
#'   FITC_sample = simulate_sample("FITC", M_demo),
#'   PE_sample = simulate_sample("PE", M_demo),
#'   APC_sample = simulate_sample("APC", M_demo)
#' ))
#'
#' unmixed <- unmix_samples(toy_fs, M = M_demo, unmixing_method = "OLS", output_dir = tempdir())
#' names(unmixed)
#' @export
unmix_samples <- function(sample_dir = "samples",
                          samples_dir = NULL,
                          M = NULL, 
                          unmixing_matrix_file = file.path("spectreasy_outputs", "unmix_controls", "scc_reference_matrix.csv"),
                          detector_noise_file = NULL,
                          control_file = NULL,
                          unmixing_method = "Spectreasy", 
                          rwls_max_iter = 1L,
                          n_threads = 1L,
                          spectral_variant_library = NULL,
                          spectral_variant_library_file = NULL,
                          spectral_variant_top_k = 3L,
                          spectral_variant_min_abundance = 1,
                          spectral_variant_positive_fraction = 0.02,
                          spectral_variant_min_improvement = 0.01,
                          spectreasy_weight_quantile = 0.9,
                          estimate_af = FALSE,
                          output_dir = "spectreasy_outputs",
                          write_fcs = TRUE,
                          save_report = TRUE,
                          report_format = "html",
                          save_qc_plots = FALSE,
                          qc_plot_dir = NULL,
                          plot_n_events = 10000L,
                          chunk_size = 50000L,
                          seed = NULL,
                          return_type = c("list", "flowSet", "SingleCellExperiment"),
                          verbose = TRUE) {
    spectreasy_weight_quantile_missing <- missing(spectreasy_weight_quantile)
    unmixing_matrix_file_missing <- missing(unmixing_matrix_file)
    return_type <- .match_arg_ci(
        return_type,
        c("list", "flowSet", "SingleCellExperiment"),
        "return_type"
    )
    report_format <- .match_arg_ci(report_format, c("html", "pdf"), "report_format")
    write_fcs <- .normalize_scalar_logical(write_fcs, "write_fcs")
    save_report <- .normalize_scalar_logical(save_report, "save_report")
    save_qc_plots <- .normalize_scalar_logical(save_qc_plots, "save_qc_plots")
    verbose <- .normalize_scalar_logical(verbose, "verbose")
    .with_optional_seed(seed)
    scc_dir <- NULL
    if (!is.null(samples_dir)) {
        if (!identical(sample_dir, "samples") && !identical(sample_dir, samples_dir)) {
            stop("Use only one of sample_dir or samples_dir.", call. = FALSE)
        }
        sample_dir <- samples_dir
    }
    sample_entries <- .prepare_unmix_samples_input(sample_dir)
    if (!is.null(control_file)) {
        if (!is.character(control_file) || length(control_file) != 1L ||
            is.na(control_file) || !nzchar(trimws(control_file))) {
            stop("control_file must be NULL or a single non-empty CSV path.", call. = FALSE)
        }
        control_file <- .resolve_control_file_path(control_file)
        if (!file.exists(control_file)) {
            .spectreasy_stop_missing_file(control_file, label = "control_file")
        }
    }
    if (!is.null(unmixing_matrix_file) && !unmixing_matrix_file_missing && !file.exists(unmixing_matrix_file)) {
        .spectreasy_stop_missing_file(unmixing_matrix_file, label = "unmixing_matrix_file")
    }
    unmixing_matrix_file <- .resolve_unmixing_matrix_file_for_samples(
        unmixing_matrix_file = unmixing_matrix_file,
        output_dir = output_dir
    )

    if (!is.null(M)) {
        M <- .as_reference_matrix(M, "M")
        M <- .load_detector_noise_for_unmixing(
            M,
            detector_noise_file = detector_noise_file,
            scc_dir = scc_dir
        )
    } else if (!is.null(unmixing_matrix_file) && file.exists(unmixing_matrix_file)) {
        .stop_if_static_unmixing_matrix_path(unmixing_matrix_file, arg_name = "unmixing_matrix_file")
        M <- .read_unmixing_matrix_csv(unmixing_matrix_file)
        M <- .as_reference_matrix(M, "M")
        M <- .load_detector_noise_for_unmixing(
            M,
            detector_noise_file = detector_noise_file,
            unmixing_matrix_file = unmixing_matrix_file,
            scc_dir = scc_dir
        )
    } else {
        stop(
            "No reference matrix provided. Run unmix_controls() first, then supply M or a valid unmixing_matrix_file.",
            call. = FALSE
        )
    }

    output_control_context <- .resolve_unmix_output_control_context(
        M = M,
        control_file = control_file,
        unmixing_matrix_file = unmixing_matrix_file
    )

    method <- .normalize_unmix_method(unmixing_method)
    rwls_max_iter <- .normalize_rwls_max_iter(rwls_max_iter)
    n_threads <- .normalize_n_threads(n_threads)
    if (!identical(method, "Spectreasy") && !spectreasy_weight_quantile_missing) {
        stop("spectreasy_weight_quantile is only accepted when unmixing_method = \"Spectreasy\".", call. = FALSE)
    }
    if (identical(method, "Spectreasy")) {
        spectreasy_weight_quantile <- .normalize_spectreasy_weight_quantile(spectreasy_weight_quantile)
    }
    estimate_af <- .normalize_scalar_logical(estimate_af, "estimate_af")
    chunk_size <- .normalize_unmix_chunk_size(chunk_size)
    plot_n_events <- .normalize_unmix_plot_n_events(plot_n_events)
    if (isTRUE(write_fcs) || isTRUE(save_report)) {
        if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) || !nzchar(trimws(output_dir))) {
            stop("output_dir must be a non-empty directory path when writing FCS files or a report.", call. = FALSE)
        }
        if (file.exists(output_dir) && !dir.exists(output_dir)) {
            stop("output_dir points to a file, not a directory: ", output_dir, call. = FALSE)
        }
    }
    output_paths <- .unmix_samples_output_paths(output_dir)
    unmixed_output_dir <- output_paths$unmixed_dir

    if (isTRUE(verbose)) {
        .spectreasy_console_header("unmix_samples")
        .spectreasy_console_field("Samples", length(sample_entries))
        .spectreasy_console_field("Method", method)
        .spectreasy_console_field("AF bands", .reference_af_band_count(M))
        if (is.finite(chunk_size)) {
            .spectreasy_console_field("Chunk", paste0(format(chunk_size, big.mark = ","), " events"))
        }
        if (!is.null(plot_n_events)) {
            label <- paste0(format(plot_n_events, big.mark = ","), " events/sample for plots")
            .spectreasy_console_field("Plot data", label)
        }
        if (isTRUE(write_fcs)) {
            .spectreasy_console_field("Output", .spectreasy_console_path(unmixed_output_dir))
        }
    }

    if (estimate_af) {
        if (isTRUE(verbose)) {
            .spectreasy_console_step("Estimate AF", "from stained samples")
        }
        M <- .estimate_blind_af_reference(
            M = M,
            sample_entries = sample_entries,
            n_bands = 10L,
            candidate_quantile = 0.90,
            max_training_events = 20000L,
            max_evaluation_events = 5000L,
            seed = seed,
            verbose = verbose
        )
    }

    spectral_variant_library_resolved <- if (.is_autospectral_style_method(method)) {
        .resolve_spectral_variant_library_for_unmixing(
            M = M,
            spectral_variant_library = spectral_variant_library,
            spectral_variant_library_file = spectral_variant_library_file,
            unmixing_matrix_file = unmixing_matrix_file
        )
    } else {
        NULL
    }

    results <- list()
    fluorophore_source_all <- rownames(M)
    output_marker_map <- .resolve_output_marker_map(
        fluorophore_source_all,
        control_file = output_control_context$control_file,
        control_df = output_control_context$control_df
    )

    if (isTRUE(write_fcs)) {
        if (!dir.exists(unmixed_output_dir)) {
            dir.create(unmixed_output_dir, showWarnings = FALSE, recursive = TRUE)
            if (!dir.exists(unmixed_output_dir)) {
                Sys.sleep(0.5)
                dir.create(unmixed_output_dir, showWarnings = FALSE, recursive = TRUE)
            }
        }
        if (!dir.exists(unmixed_output_dir)) {
            stop("Could not create unmixed FCS output directory: ", unmixed_output_dir, call. = FALSE)
        }
    }

    for (entry in sample_entries) {
        sn <- entry$sample_name
        if (isTRUE(verbose)) .spectreasy_console_field("Sample", sn)
        ff <- if (inherits(entry$flow_frame, "flowFrame")) {
            entry$flow_frame
        } else {
            .spectreasy_read_fcs(entry$file_path, label = "sample FCS file")
        }
        n_events <- nrow(flowCore::exprs(ff))
        if (n_events == 0L) {
            stop("Sample '", sn, "' contains no events to unmix.", call. = FALSE)
        }
        chunk_indices <- .unmix_chunk_indices(n_events, chunk_size = chunk_size)
        keep_global <- .unmix_sample_keep_indices(n_events, plot_n_events = plot_n_events)
        fluorophore_source_all <- rownames(M)
        fluorophore_source_for_output <- fluorophore_source_all[
            !grepl("^AF_", fluorophore_source_all, ignore.case = TRUE)
        ]

        write_data_all <- NULL
        data_chunks <- list()
        residual_chunks <- list()
        variant_infos <- list()
        spectreasy_decoder_weights <- NULL
        fluorophores_to_keep <- character()

        for (chunk_i in seq_along(chunk_indices)) {
            event_idx <- chunk_indices[[chunk_i]]
            chunk_ff <- if (length(chunk_indices) == 1L) ff else ff[event_idx, ]
            calc_args <- list(
                flow_frame = chunk_ff,
                M = M,
                method = method,
                file_name = sn,
                rwls_max_iter = rwls_max_iter,
                n_threads = n_threads,
                spectral_variant_library = spectral_variant_library_resolved,
                spectral_variant_top_k = spectral_variant_top_k,
                spectral_variant_min_abundance = spectral_variant_min_abundance,
                spectral_variant_positive_fraction = spectral_variant_positive_fraction,
                spectral_variant_min_improvement = spectral_variant_min_improvement,
                return_residuals = TRUE
            )
            if (identical(method, "Spectreasy")) {
                calc_args$spectreasy_weight_quantile <- spectreasy_weight_quantile
            }
            res_chunk <- do.call(calc_residuals, calc_args)

            if (isTRUE(write_fcs)) {
                write_data <- .aggregate_af_columns_for_output(
                    res_chunk$data,
                    fluorophore_source = fluorophore_source_all
                )
                fluorophores_to_keep <- intersect(colnames(write_data), fluorophore_source_for_output)
                autospectral_cols <- if (.is_autospectral_style_method(method) && "AF Index" %in% colnames(write_data)) "AF Index" else character()
                passthrough_cols <- .get_passthrough_parameter_names(colnames(write_data))
                cols_to_write <- unique(c(fluorophores_to_keep, autospectral_cols, passthrough_cols))
                write_chunk <- as.matrix(write_data[, cols_to_write, drop = FALSE])
                if (is.null(write_data_all)) {
                    write_data_all <- matrix(
                        NA_real_, nrow = n_events, ncol = ncol(write_chunk),
                        dimnames = list(NULL, colnames(write_chunk))
                    )
                } else if (!identical(colnames(write_chunk), colnames(write_data_all))) {
                    stop("Internal error: unmixed FCS columns changed between chunks for sample '", sn, "'.", call. = FALSE)
                }
                write_data_all[event_idx, ] <- write_chunk
            }

            local_keep <- which(event_idx %in% keep_global)
            if (length(local_keep) > 0L) {
                data_chunks[[length(data_chunks) + 1L]] <- res_chunk$data[local_keep, , drop = FALSE]
                if (!is.null(res_chunk$residuals)) {
                    residual_chunks[[length(residual_chunks) + 1L]] <- res_chunk$residuals[local_keep, , drop = FALSE]
                }
            }
            if (!is.null(res_chunk$spectral_variant_info)) {
                variant_infos[[length(variant_infos) + 1L]] <- res_chunk$spectral_variant_info
            }
            if (is.null(spectreasy_decoder_weights) && !is.null(res_chunk$spectreasy_decoder_weights)) {
                spectreasy_decoder_weights <- res_chunk$spectreasy_decoder_weights
            }

            rm(res_chunk, chunk_ff)
            if (chunk_i %% 5L == 0L) {
                gc(verbose = FALSE)
            }
        }

        if (isTRUE(write_fcs)) {
            unmixed_exprs <- write_data_all
            storage.mode(unmixed_exprs) <- "numeric"
            new_ff <- flowCore::flowFrame(unmixed_exprs)
            new_ff <- .apply_output_fcs_feature_labels(
                target_ff = new_ff,
                source_ff = ff,
                fluorophore_cols = fluorophores_to_keep,
                marker_label_map = output_marker_map
            )

            target_output_path <- file.path(unmixed_output_dir, .unmixed_fcs_filename(sn, method, M))
            output_path <- .next_safe_output_path(target_output_path)
            if (!identical(output_path, target_output_path)) {
                .spectreasy_console_step("Safe output", basename(output_path))
            }
            flowCore::write.FCS(new_ff, output_path)
            rm(write_data_all, unmixed_exprs, new_ff)
        }

        results[[sn]] <- .combine_unmix_result_chunks(
            data_chunks = data_chunks,
            residual_chunks = residual_chunks,
            method = method,
            M = M,
            variant_infos = variant_infos,
            spectreasy_decoder_weights = spectreasy_decoder_weights
        )
        rm(ff, data_chunks, residual_chunks, variant_infos)
        gc(verbose = FALSE)
    }

    class(results) <- c("spectreasy_unmixed_results", "list")
    attr(results, "method") <- method
    attr(results, "reference_matrix") <- M
    attr(results, "spectral_variant_library") <- spectral_variant_library_resolved
    attr(results, "blind_af_info") <- attr(M, "blind_af_info")
    attr(results, "qc_report_file") <- NULL
    attr(results, "qc_samples_dir") <- NULL
    attr(results, "qc_metrics_dir") <- NULL
    attr(results, "qc_plot_dir") <- NULL
    attr(results, "qc_report_data") <- NULL

    if (isTRUE(save_report)) {
        qc_samples_dir <- .next_safe_output_dir(output_paths$qc_samples_dir)
        output_file <- file.path(qc_samples_dir, paste0("qc_samples_report.", report_format))
        .spectreasy_console_field("Report", .spectreasy_console_path(output_file))
        report_res <- qc_samples(
            results = results,
            M = M,
            output_file = output_file,
            unmixing_method = method,
            qc_metrics_dir = qc_samples_dir,
            qc_plot_dir = qc_plot_dir,
            save_qc_pngs = save_qc_plots,
            report_format = report_format,
            overwrite = "overwrite",
            report_artifact_paths = list(
                matrix = if (!isTRUE(unmixing_matrix_file_missing)) unmixing_matrix_file else NULL,
                detector_noise = detector_noise_file,
                spectral_variant_library = spectral_variant_library_file,
                output_dir = unmixed_output_dir
            ),
            report_run_settings = list(
                rwls_max_iter = rwls_max_iter,
                n_threads = n_threads,
                spectral_variant_top_k = spectral_variant_top_k,
                spectral_variant_min_abundance = spectral_variant_min_abundance,
                spectral_variant_positive_fraction = spectral_variant_positive_fraction,
                spectral_variant_min_improvement = spectral_variant_min_improvement,
                spectreasy_weight_quantile = spectreasy_weight_quantile,
                estimate_af = estimate_af,
                write_fcs = write_fcs,
                save_qc_plots = save_qc_plots,
                plot_n_events = plot_n_events,
                chunk_size = chunk_size,
                seed = seed
            )
        )
        attr(results, "qc_report_file") <- report_res$output_file
        attr(results, "qc_samples_dir") <- qc_samples_dir
        attr(results, "qc_metrics_dir") <- qc_samples_dir
        attr(results, "qc_plot_dir") <- report_res$qc_plot_dir
        attr(results, "qc_report_data") <- report_res$report_data
        if (!file.exists(report_res$output_file)) {
            stop("Automatic sample QC report was requested but was not created at: ", report_res$output_file, call. = FALSE)
        }
    }

    if (identical(return_type, "flowSet")) {
        return(invisible(.unmixed_results_to_flowset(results)))
    }
    if (identical(return_type, "SingleCellExperiment")) {
        return(invisible(.unmixed_results_to_sce(results, sample_entries = sample_entries)))
    }

    invisible(results)
}
