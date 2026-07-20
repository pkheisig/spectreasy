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
