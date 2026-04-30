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
    cell_names <- colnames(sample_input)
    if (is.null(cell_names) || any(cell_names == "")) {
        cell_names <- paste0("cell_", seq_len(ncol(sample_input)))
    }

    lapply(sample_levels, function(sample_id) {
        idx <- which(sample_ids == sample_id)
        expr_mat <- t(as.matrix(assay_mat[, idx, drop = FALSE]))
        storage.mode(expr_mat) <- "numeric"
        colnames(expr_mat) <- rownames(assay_mat)
        ff <- flowCore::flowFrame(expr_mat)
        list(
            sample_name = make.unique(tools::file_path_sans_ext(basename(sample_id)))[[1]],
            source_name = sample_id,
            flow_frame = ff,
            cell_ids = cell_names[idx]
        )
    })
}

.prepare_unmix_samples_input_from_directory <- function(sample_input) {
    sample_dir_resolved <- normalizePath(sample_input, mustWork = TRUE)
    fcs_files <- list.files(sample_dir_resolved, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
    if (length(fcs_files) == 0) stop("No FCS files found in ", sample_dir_resolved)

    lapply(fcs_files, function(f) {
        list(
            sample_name = tools::file_path_sans_ext(basename(f)),
            source_name = basename(f),
            file_path = f,
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

    if (!is.character(sample_input) || length(sample_input) != 1 || is.na(sample_input)) {
        stop("sample_dir must be either a directory path, a flowCore::flowSet, or a SingleCellExperiment.")
    }

    .prepare_unmix_samples_input_from_directory(sample_input)
}

.build_result_from_static_unmix <- function(flow_frame, W_sub, file_name) {
    full_data <- flowCore::exprs(flow_frame)
    detectors <- colnames(W_sub)
    Y <- full_data[, detectors, drop = FALSE]
    abundances <- Y %*% t(W_sub)
    colnames(abundances) <- rownames(W_sub)
    residuals <- NULL

    out <- as.data.frame(abundances)
    out <- .append_passthrough_parameters(out, full_data, detector_names = detectors)
    out$File <- file_name

    list(data = out, residuals = residuals)
}

.read_unmixing_matrix_csv <- function(path) {
    if (!file.exists(path)) stop("unmixing_matrix_file not found: ", path)
    df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    if (ncol(df) < 2) stop("Invalid unmixing matrix CSV: expected at least 2 columns in ", path)

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

.resolve_secondary_label_map <- function(primary_names, sample_dir) {
    primary_names <- trimws(as.character(primary_names))
    labels <- stats::setNames(primary_names, primary_names)
    if (length(primary_names) == 0) return(labels)

    opt_control <- getOption("spectreasy.control_file", "")
    sample_parent <- if (is.character(sample_dir) && length(sample_dir) == 1 && !is.na(sample_dir)) {
        dirname(normalizePath(sample_dir, mustWork = FALSE))
    } else {
        getwd()
    }
    candidates <- unique(c(
        as.character(opt_control),
        .resolve_control_file_path("fcs_mapping.csv"),
        file.path(sample_parent, "fcs_mapping.csv")
    ))
    candidates <- candidates[!is.na(candidates) & nzchar(trimws(candidates))]
    existing <- candidates[file.exists(candidates)]
    if (length(existing) == 0) return(labels)

    control_df <- tryCatch(
        utils::read.csv(existing[1], stringsAsFactors = FALSE, check.names = FALSE),
        error = function(e) NULL
    )
    if (is.null(control_df) || nrow(control_df) == 0) return(labels)

    lower_names <- tolower(colnames(control_df))
    fluor_idx <- match("fluorophore", lower_names)
    marker_idx <- match("marker", lower_names)
    if (is.na(fluor_idx) || is.na(marker_idx)) return(labels)

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

.apply_feature_secondary_labels <- function(target_ff, source_ff, marker_cols, secondary_label_map) {
    pd_new <- flowCore::pData(flowCore::parameters(target_ff))
    if (!all(c("name", "desc") %in% colnames(pd_new))) {
        return(target_ff)
    }

    param_names <- as.character(pd_new$name)
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

    marker_cols <- as.character(marker_cols)
    for (mk in marker_cols) {
        idx <- which(param_names == mk)
        if (length(idx) == 0) next
        lbl <- if (mk %in% names(secondary_label_map)) secondary_label_map[[mk]] else mk
        if (is.na(lbl) || !nzchar(trimws(lbl))) lbl <- mk
        desc_vals[idx] <- as.character(lbl)
    }

    pd_new$desc <- desc_vals
    flowCore::parameters(target_ff) <- methods::new("AnnotatedDataFrame", data = pd_new)
    target_ff
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
#' @param sample_dir Directory containing experimental FCS files, a
#'   `flowCore::flowSet`, or a `SingleCellExperiment` for in-memory workflows.
#' @param M Optional reference matrix (Markers x Detectors). Required for
#'   dynamic unmixing methods (`"WLS"`, `"OLS"`, `"NNLS"`).
#' @param W Optional static unmixing matrix (Markers x Detectors). If supplied,
#'   unmixing is performed using matrix multiplication with the transposed W.
#' @param unmixing_matrix_file Optional CSV path to a saved unmixing matrix.
#'   Used when `W` is not supplied. By default this points to the matrix produced
#'   by [autounmix_controls()].
#' @param method Unmixing method (`"WLS"`, `"OLS"`, or `"NNLS"`).
#' @param cytometer Reserved for compatibility with older workflows.
#' @param output_dir Directory to save unmixed FCS files when `write_fcs = TRUE`.
#' @param write_fcs Logical; if `TRUE`, write unmixed FCS files to `output_dir`.
#'   Defaults to `TRUE` so unmixed FCS files are written unless disabled explicitly.
#' @param return_type Return format: `"list"` (default), `"flowSet"`, or
#'   `"SingleCellExperiment"`. When `"flowSet"`, detector residuals are attached
#'   as `attr(x, "spectreasy_residuals")`. When `"SingleCellExperiment"`, cell-level
#'   unmixed values are returned in assay `"unmixed"`, with detector residuals in
#'   `altExp(x, "detector_residuals")` when available.
#' @return Either a named list with one element per sample, a `flowSet`, or a
#'   `SingleCellExperiment` depending on `return_type`. List elements contain
#'   `data` (unmixed abundances plus retained acquisition parameters) and
#'   `residuals` (detector residual matrix when available, otherwise `NULL`).
#'   The return value is provided invisibly to avoid printing large result objects
#'   during interactive or Quarto execution.
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
#' unmixed <- unmix_samples(toy_fs, M = M_demo, method = "OLS", output_dir = tempdir())
#' names(unmixed)
#' @export
unmix_samples <- function(sample_dir = "samples", 
                          M = NULL, 
                          W = NULL,
                          unmixing_matrix_file = file.path("spectreasy_outputs", "autounmix_controls", "scc_unmixing_matrix.csv"),
                          method = "WLS", 
                          cytometer = "Aurora",
                          output_dir = file.path("spectreasy_outputs", "unmix_samples"),
                          write_fcs = TRUE,
                          return_type = c("list", "flowSet", "SingleCellExperiment")) {
    return_type <- match.arg(return_type)

    if (!is.null(M)) {
        M <- .as_reference_matrix(M, "M")
    }

    sample_entries <- .prepare_unmix_samples_input(sample_dir)

    using_static_W <- FALSE
    W_use <- NULL
    if (!is.null(W)) {
        W_use <- as.matrix(W)
        using_static_W <- TRUE
    } else if (is.null(M) && !is.null(unmixing_matrix_file)) {
        W_use <- .read_unmixing_matrix_csv(unmixing_matrix_file)
        using_static_W <- TRUE
    }

    if (is.null(M) && !using_static_W) {
        stop(
            "No unmixing input provided. Supply either:\n",
            " - M (reference matrix), or\n",
            " - W / unmixing_matrix_file (static unmixing matrix)."
        )
    }
    if (!is.null(M) && using_static_W) {
        message("Both M and W/unmixing_matrix_file provided. Using static unmixing matrix.")
    }

    method_upper <- toupper(method)
    allowed_methods <- c("WLS", "OLS", "NNLS")
    if (!using_static_W && !(method_upper %in% allowed_methods)) {
        stop("method must be one of: ", paste(allowed_methods, collapse = ", "))
    }

    results <- list()
    marker_source_all <- if (using_static_W) rownames(W_use) else rownames(M)
    secondary_label_map <- .resolve_secondary_label_map(marker_source_all, sample_dir = sample_dir)

    if (isTRUE(write_fcs)) {
        dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    }

    for (entry in sample_entries) {
        sn <- entry$sample_name
        message("  Unmixing sample: ", sn)
        ff <- if (inherits(entry$flow_frame, "flowFrame")) {
            entry$flow_frame
        } else {
            flowCore::read.FCS(entry$file_path, transformation = FALSE, truncate_max_range = FALSE)
        }

        if (using_static_W) {
            raw_data <- flowCore::exprs(ff)
            detectors <- colnames(W_use)
            missing <- setdiff(detectors, colnames(raw_data))
            if (length(missing) > 0) {
                stop("Detectors in unmixing matrix not found in sample '", sn, "': ", paste(missing, collapse = ", "))
            }
            res_obj <- .build_result_from_static_unmix(ff, W_use, sn)
        } else {
            res_obj <- calc_residuals(ff, M, method = method_upper, file_name = sn, return_residuals = TRUE)
        }
        
        if (isTRUE(write_fcs)) {
            marker_source <- if (using_static_W) rownames(W_use) else rownames(M)
            markers_to_keep <- intersect(colnames(res_obj$data), marker_source)
            passthrough_cols <- .get_passthrough_parameter_names(colnames(res_obj$data))
            cols_to_write <- unique(c(markers_to_keep, passthrough_cols))
            unmixed_exprs <- as.matrix(res_obj$data[, cols_to_write, drop = FALSE])
            
            new_ff <- flowCore::flowFrame(unmixed_exprs)
            new_ff <- .apply_feature_secondary_labels(
                target_ff = new_ff,
                source_ff = ff,
                marker_cols = markers_to_keep,
                secondary_label_map = secondary_label_map
            )

            output_path <- .next_safe_output_path(file.path(output_dir, paste0(sn, "_unmixed.fcs")))
            if (!identical(output_path, file.path(output_dir, paste0(sn, "_unmixed.fcs")))) {
                message("    Existing output detected; writing to safe path: ", basename(output_path))
            }
            flowCore::write.FCS(new_ff, output_path)
        }
        results[[sn]] <- res_obj
    }
    
    if (identical(return_type, "flowSet")) {
        return(invisible(.unmixed_results_to_flowset(results)))
    }
    if (identical(return_type, "SingleCellExperiment")) {
        return(invisible(.unmixed_results_to_sce(results, sample_entries = sample_entries)))
    }

    invisible(results)
}
