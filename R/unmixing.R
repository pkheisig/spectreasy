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

.start_unmix_samples_progress <- function(total, verbose) {
    if (!isTRUE(verbose) || !interactive() || total < 1L) {
        return(NULL)
    }
    utils::txtProgressBar(min = 0, max = total, initial = 0, style = 3)
}

.tick_unmix_samples_progress <- function(progress_bar, value) {
    if (is.null(progress_bar)) {
        return(invisible(FALSE))
    }
    utils::setTxtProgressBar(progress_bar, value)
    invisible(TRUE)
}

.close_unmix_samples_progress <- function(progress_bar) {
    if (!is.null(progress_bar)) {
        close(progress_bar)
    }
    invisible(NULL)
}

.default_unmix_samples_report_file <- function(output_dir) {
    output_dir <- as.character(output_dir)[1]
    if (is.na(output_dir) || !nzchar(output_dir)) {
        output_dir <- file.path("spectreasy_outputs", "unmix_samples", "unmixed_fcs")
    }
    report_dir <- if (identical(basename(normalizePath(output_dir, mustWork = FALSE)), "unmixed_fcs")) {
        dirname(output_dir)
    } else {
        output_dir
    }
    file.path(report_dir, "qc_samples_report.pdf")
}

.default_unmix_samples_report_dir <- function(output_dir) {
    output_dir <- as.character(output_dir)[1]
    if (is.na(output_dir) || !nzchar(output_dir)) {
        output_dir <- file.path("spectreasy_outputs", "unmix_samples", "unmixed_fcs")
    }
    report_parent <- if (identical(basename(normalizePath(output_dir, mustWork = FALSE)), "unmixed_fcs")) {
        dirname(output_dir)
    } else {
        output_dir
    }
    file.path(report_parent, "qc_samples")
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

.attach_matching_variances <- function(M, V, source = "variances") {
    M <- .as_reference_matrix(M, "M")
    V <- .as_reference_matrix(V, source)

    missing_rows <- setdiff(rownames(M), rownames(V))
    missing_cols <- setdiff(colnames(M), colnames(V))
    if (length(missing_rows) > 0 || length(missing_cols) > 0) {
        stop(
            "WLS variances from ", source, " do not match the reference matrix.\n",
            if (length(missing_rows) > 0) paste0("Missing marker rows: ", paste(missing_rows, collapse = ", "), "\n") else "",
            if (length(missing_cols) > 0) paste0("Missing detector columns: ", paste(missing_cols, collapse = ", ")) else ""
        )
    }

    V <- V[rownames(M), colnames(M), drop = FALSE]
    attr(M, "variances") <- V
    M
}

.load_variances_for_unmixing <- function(M, variances_file = NULL) {
    if (is.null(variances_file) || !file.exists(variances_file)) {
        return(M)
    }

    V <- .read_unmixing_matrix_csv(variances_file)
    .attach_matching_variances(M, V, source = variances_file)
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

.resolve_variances_file_for_unmixing <- function(unmixing_matrix_file,
                                                variances_file = NULL,
                                                prefer_sibling = FALSE) {
    if (is.null(variances_file)) {
        return(variances_file)
    }

    if (!is.null(unmixing_matrix_file) && file.exists(unmixing_matrix_file)) {
        sibling_variances_file <- file.path(dirname(unmixing_matrix_file), "scc_variances.csv")
        if (isTRUE(prefer_sibling) && file.exists(sibling_variances_file)) {
            return(sibling_variances_file)
        }
    }

    if (file.exists(variances_file)) {
        return(variances_file)
    }

    default_variances_file <- file.path("spectreasy_outputs", "unmix_controls", "scc_variances.csv")
    if (!identical(normalizePath(variances_file, mustWork = FALSE), normalizePath(default_variances_file, mustWork = FALSE))) {
        return(variances_file)
    }

    if (is.null(unmixing_matrix_file) || !file.exists(unmixing_matrix_file)) {
        return(variances_file)
    }

    sibling_variances_file <- file.path(dirname(unmixing_matrix_file), "scc_variances.csv")
    if (file.exists(sibling_variances_file)) {
        return(sibling_variances_file)
    }

    variances_file
}

.ensure_wls_variances <- function(M,
                                  method,
                                  variances_file = NULL,
                                  scc_dir = NULL,
                                  control_file = NULL,
                                  af_n_bands = "auto",
                                  af_auto_max_bands = 100,
                                  af_min_cluster_events = 20,
                                  af_min_cluster_proportion = 0.005,
                                  af_n_bands_sensitivity = 1.5,
                                  af_refine = FALSE,
                                  af_refine_problem_quantile = 0.99,
                                  exclude_af = FALSE,
                                  cytometer = "auto",
                                  seed = NULL) {
    if (!(toupper(method) %in% c("WLS", "RWLS"))) {
        return(M)
    }

    .load_variances_for_unmixing(M, variances_file = variances_file)
}

.resolve_unmix_samples_matrix <- function(M,
                                          unmixing_matrix_file,
                                          variances_file,
                                          variances_file_was_missing,
                                          detector_noise_file = NULL,
                                          scc_dir = NULL) {
    if (!is.null(M)) {
        M <- .as_reference_matrix(M, "M")
        M <- .load_variances_for_unmixing(M, variances_file = variances_file)
        M <- .load_detector_noise_for_unmixing(
            M,
            detector_noise_file = detector_noise_file,
            scc_dir = scc_dir
        )
        return(list(M = M, variances_file = variances_file))
    }

    if (!is.null(unmixing_matrix_file) && file.exists(unmixing_matrix_file)) {
        .stop_if_static_unmixing_matrix_path(unmixing_matrix_file, arg_name = "unmixing_matrix_file")
        M <- .read_unmixing_matrix_csv(unmixing_matrix_file)
        M <- .as_reference_matrix(M, "M")
        variances_file <- .resolve_variances_file_for_unmixing(
            unmixing_matrix_file = unmixing_matrix_file,
            variances_file = variances_file,
            prefer_sibling = variances_file_was_missing
        )
        M <- .load_variances_for_unmixing(M, variances_file = variances_file)
        M <- .load_detector_noise_for_unmixing(
            M,
            detector_noise_file = detector_noise_file,
            unmixing_matrix_file = unmixing_matrix_file,
            scc_dir = scc_dir
        )
        return(list(M = M, variances_file = variances_file))
    }

    if (!is.null(unmixing_matrix_file)) {
        stop(
            "Reference matrix not found: ", unmixing_matrix_file, ".\n",
            "Run unmix_controls() first, or pass an explicit reference matrix with M.",
            call. = FALSE
        )
    }

    stop(
        "No reference matrix was provided.\n",
        "Run unmix_controls() first, or pass an explicit reference matrix with M.",
        call. = FALSE
    )
}

.normalize_unmix_samples_options <- function(method,
                                             rwls_max_iter,
                                             multithreading,
                                             n_threads) {
    method_label <- .normalize_unmix_method(method)

    list(
        method_label = method_label,
        solver_method = .solver_method_for_unmix(method_label),
        rwls_max_iter = .normalize_rwls_max_iter(rwls_max_iter),
        n_threads = .normalize_unmix_threads(multithreading = multithreading, n_threads = n_threads)
    )
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

.ensure_unmix_output_dir <- function(output_dir, write_fcs) {
    if (!isTRUE(write_fcs)) {
        return(invisible(NULL))
    }
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
        if (!dir.exists(output_dir)) {
            Sys.sleep(0.5)
            dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
        }
    }
    if (!dir.exists(output_dir)) {
        stop("Could not create output_dir: ", output_dir, call. = FALSE)
    }
    invisible(NULL)
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

.next_safe_output_dir <- function(path) {
    if (!dir.exists(path) && !file.exists(path)) {
        return(path)
    }

    parent <- dirname(path)
    stem <- basename(path)
    i <- 2L
    repeat {
        candidate <- file.path(parent, paste0(stem, "_", i))
        if (!dir.exists(candidate) && !file.exists(candidate)) {
            return(candidate)
        }
        i <- i + 1L
    }
}

.read_unmix_sample_flow_frame <- function(entry) {
    if (inherits(entry$flow_frame, "flowFrame")) {
        return(entry$flow_frame)
    }
    flowCore::read.FCS(entry$file_path, transformation = FALSE, truncate_max_range = FALSE)
}

.make_unmixed_output_frame <- function(res_obj, M, source_ff, secondary_label_map) {
    marker_source <- rownames(M)
    marker_source <- marker_source[!grepl("^AF_", marker_source, ignore.case = TRUE)]
    markers_to_keep <- intersect(colnames(res_obj$data), marker_source)
    passthrough_cols <- .get_passthrough_parameter_names(colnames(res_obj$data))
    cols_to_write <- unique(c(markers_to_keep, passthrough_cols))
    unmixed_exprs <- as.matrix(res_obj$data[, cols_to_write, drop = FALSE])

    new_ff <- flowCore::flowFrame(unmixed_exprs)
    .apply_feature_secondary_labels(
        target_ff = new_ff,
        source_ff = source_ff,
        marker_cols = markers_to_keep,
        secondary_label_map = secondary_label_map
    )
}

.write_unmixed_sample_fcs <- function(res_obj, M, source_ff, output_dir, sample_name, secondary_label_map) {
    new_ff <- .make_unmixed_output_frame(
        res_obj = res_obj,
        M = M,
        source_ff = source_ff,
        secondary_label_map = secondary_label_map
    )

    default_path <- file.path(output_dir, paste0(sample_name, "_unmixed.fcs"))
    output_path <- .next_safe_output_path(default_path)
    if (!identical(output_path, default_path)) {
        message("    Existing output detected; writing to safe path: ", basename(output_path))
    }
    flowCore::write.FCS(new_ff, output_path)
    output_path
}

.subsample_unmix_result <- function(res_obj, subsample_n) {
    if (is.null(subsample_n)) {
        return(res_obj)
    }

    n_events <- nrow(res_obj$data)
    if (n_events > subsample_n) {
        idx <- sample.int(n_events, subsample_n)
        res_obj$data <- res_obj$data[idx, , drop = FALSE]
        if (!is.null(res_obj$residuals)) {
            res_obj$residuals <- res_obj$residuals[idx, , drop = FALSE]
        }
    }
    res_obj
}

.unmix_one_sample <- function(entry,
                              M,
                              method_label,
                              rwls_max_iter,
                              multithreading,
                              n_threads,
                              spectral_variant_library,
                              spectral_variant_top_k,
                              spectral_variant_min_abundance,
                              write_fcs,
                              output_dir,
                              secondary_label_map,
                              subsample_n,
                              verbose) {
    sample_name <- entry$sample_name
    if (isTRUE(verbose)) message("  Unmixing sample: ", sample_name)

    ff <- .read_unmix_sample_flow_frame(entry)
    res_obj <- calc_residuals(
        ff,
        M,
        method = method_label,
        file_name = sample_name,
        rwls_max_iter = rwls_max_iter,
        multithreading = multithreading,
        n_threads = n_threads,
        spectral_variant_library = spectral_variant_library,
        spectral_variant_top_k = spectral_variant_top_k,
        spectral_variant_min_abundance = spectral_variant_min_abundance,
        return_residuals = TRUE
    )

    if (isTRUE(write_fcs)) {
        .write_unmixed_sample_fcs(
            res_obj = res_obj,
            M = M,
            source_ff = ff,
            output_dir = output_dir,
            sample_name = sample_name,
            secondary_label_map = secondary_label_map
        )
    }

    .subsample_unmix_result(res_obj, subsample_n = subsample_n)
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

    all_cols <- unique(unlist(lapply(data_parts, colnames), use.names = FALSE))
    data_parts <- lapply(data_parts, function(data_df) {
        missing_cols <- setdiff(all_cols, colnames(data_df))
        for (col in missing_cols) {
            data_df[[col]] <- NA
        }
        data_df[, all_cols, drop = FALSE]
    })

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
#' @param sample_dir Directory containing experimental FCS files, a
#'   `flowCore::flowSet`, or a `SingleCellExperiment` for in-memory workflows.
#' @param M Optional reference matrix (Markers x Detectors). If supplied,
#'   unmixing is computed dynamically using this matrix. If not supplied,
#'   it is loaded from the CSV path provided in `unmixing_matrix_file`.
#' @param unmixing_matrix_file Optional CSV path to a saved reference matrix.
#'   Used when `M` is not supplied. By default this points to the reference matrix
#'   produced by [unmix_controls()] (`"scc_reference_matrix.csv"`).
#' @param variances_file Optional CSV path to a saved variances matrix. SCC
#'   variances are loaded as reference QC metadata when available, but the
#'   default WLS model does not use them as detector weights.
#' @param detector_noise_file Optional CSV path to detector-specific WLS noise
#'   floors, as written by [unmix_controls()] (`"scc_detector_noise.csv"`). If
#'   omitted, `unmix_samples()` first looks beside `unmixing_matrix_file`, then
#'   estimates the floors from `scc_dir` when available, and otherwise falls
#'   back to the built-in scalar noise floor.
#' @param method Unmixing method (`"AutoSpectral"`, `"OLS"`, `"WLS"`,
#'   `"RWLS"`, or `"NNLS"`). `AutoSpectral` is the default and combines
#'   per-cell SCC/AF band matching with an OLS refit.
#' @param rwls_max_iter Positive integer; number of robust reweighting
#'   iterations used when `method = "RWLS"`. The default, 1, preserves the
#'   historical behavior.
#' @param multithreading Logical; if `TRUE`, allow event-wise multi-AF WLS/RWLS
#'   unmixing to use multiple threads. The default, `FALSE`, keeps execution
#'   single-threaded.
#' @param n_threads `"auto"` or positive integer; thread count to use when
#'   `multithreading = TRUE`. `"auto"` uses `RcppParallel::defaultNumThreads()`.
#'   Integers larger than the available thread count are clipped to the
#'   available count.
#' @param cytometer Legacy compatibility argument. Reference-matrix construction
#'   now belongs in [unmix_controls()]; `unmix_samples()` requires either `M` or
#'   a valid `unmixing_matrix_file`.
#' @param scc_dir Optional SCC directory used only to locate saved detector-noise
#'   metadata when available. It is not used to rebuild a missing reference
#'   matrix.
#' @param control_file Legacy compatibility argument. To create a reference
#'   matrix from a control mapping, run [unmix_controls()] first.
#' @param af_n_bands Number of k-means AF basis signatures to extract from
#'   pooled unstained/AF control events, or `"auto"` to keep distinct signatures
#'   from up to `af_auto_max_bands` k-means centers. Legacy compatibility
#'   argument; build AF references with [unmix_controls()].
#' @param af_auto_max_bands Maximum k-means centers that `"auto"` may score.
#'   Legacy compatibility argument; build AF
#'   references with [unmix_controls()].
#' @param af_min_cluster_events Minimum k-means cluster size retained in
#'   `"auto"` mode.
#' @param af_min_cluster_proportion Minimum k-means cluster proportion retained
#'   in `"auto"` mode.
#' @param af_n_bands_sensitivity Compatibility argument retained for older
#'   workflows.
#' @param af_refine Logical; if `TRUE`, run the optional second-pass AF
#'   refinement. Legacy compatibility argument; build refined AF references with
#'   [unmix_controls()].
#' @param af_refine_problem_quantile Quantile used to choose high-error
#'   unstained cells for AF refinement. Legacy compatibility argument.
#' @param af_deduplicate Logical; if `TRUE`, remove near-identical AF spectra
#'   when building a reference matrix. Legacy compatibility argument.
#' @param af_deduplication_threshold Cosine similarity threshold used when
#'   `af_deduplicate = TRUE`.
#' @param af_contaminant_threshold Cosine similarity threshold used to remove
#'   AF candidates that are too similar to known fluorophore spectra. Legacy
#'   compatibility argument.
#' @param spectral_variant_library Optional in-memory spectral-variant library,
#'   usually returned by [unmix_controls()].
#' @param spectral_variant_library_file Optional `.rds` path to a saved
#'   spectral-variant library. If omitted, `unmix_samples()` looks for
#'   `scc_spectral_variants.rds` beside `unmixing_matrix_file`.
#' @param spectral_variant_top_k Number of best variant candidates to test per
#'   positive fluorophore.
#' @param spectral_variant_min_abundance Minimum unmixed abundance for a
#'   fluorophore to be eligible for variant testing.
#' @param exclude_af Legacy compatibility argument. Exclude AF while building
#'   the matrix in [unmix_controls()], not in `unmix_samples()`.
#' @param estimate_af Logical; if `TRUE`, estimate AF signatures directly from
#'   stained sample event-wise WLS residuals, select the best candidate model by
#'   held-out WLS residual score, and append the selected AF rows to the
#'   reference matrix before unmixing. This is intended for workflows where no
#'   unstained cell control is available. Default is `FALSE`.
#' @param output_dir Directory to save unmixed FCS files when `write_fcs = TRUE`.
#' @param write_fcs Logical; if `TRUE`, write unmixed FCS files to `output_dir`.
#'   Defaults to `TRUE`. Existing FCS files are not overwritten; a numeric suffix
#'   is added when needed.
#' @param save_report Logical; if `TRUE`, write a sample QC PDF report and
#'   sample QC metric CSVs from the in-memory unmixing results without rerunning
#'   unmixing. Defaults to `TRUE`.
#' @param output_file Optional output path for the sample QC PDF report. When
#'   this is `NULL`, each `unmix_samples()` report run is written to a fresh
#'   `qc_samples`, `qc_samples_2`, ... folder beside the unmixed FCS output.
#' @param save_qc_plots Logical; if `TRUE`, save QC report plots as PNG files
#'   in `qc_plot_dir` while creating the PDF report.
#' @param qc_plot_dir Directory for sample QC report PNG files when
#'   `save_qc_plots = TRUE`.
#' @param subsample_n Optional integer; if provided, subsample the returned in-memory
#'   results to at most `subsample_n` events per sample. The full unsampled unmixed
#'   data will still be written to the FCS files when `write_fcs = TRUE`.
#' @param seed Optional integer seed for deterministic subsampling.
#' @param return_type Return format: `"list"` (default), `"flowSet"`, or
#'   `"SingleCellExperiment"`. When `"flowSet"`, detector residuals are attached
#'   as `attr(x, "spectreasy_residuals")`. When `"SingleCellExperiment"`, cell-level
#'   unmixed values are returned in assay `"unmixed"`, with detector residuals in
#'   `altExp(x, "detector_residuals")` when available.
#' @param verbose Logical; if `TRUE`, print progress messages while unmixing
#'   each sample and show an interactive console progress bar when available.
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
#' unmixed <- unmix_samples(toy_fs, M = M_demo, method = "OLS", output_dir = tempdir())
#' names(unmixed)
#' @export
unmix_samples <- function(sample_dir = "samples", 
                          M = NULL, 
                          unmixing_matrix_file = file.path("spectreasy_outputs", "unmix_controls", "scc_reference_matrix.csv"),
                          variances_file = file.path("spectreasy_outputs", "unmix_controls", "scc_variances.csv"),
                          detector_noise_file = NULL,
                          method = "AutoSpectral",
                          rwls_max_iter = 1L,
                          multithreading = FALSE,
                          n_threads = "auto",
                          cytometer = "auto",
                          scc_dir = NULL,
                          control_file = NULL,
                          af_n_bands = "auto",
                          af_auto_max_bands = 100,
                          af_min_cluster_events = 20,
                          af_min_cluster_proportion = 0.005,
                          af_n_bands_sensitivity = 1.5,
                          af_refine = FALSE,
                          af_refine_problem_quantile = 0.99,
                          af_deduplicate = FALSE,
                          af_deduplication_threshold = 0.99,
                          af_contaminant_threshold = 0.99,
                          spectral_variant_library = NULL,
                          spectral_variant_library_file = NULL,
                          spectral_variant_top_k = 3L,
                          spectral_variant_min_abundance = 1,
                          exclude_af = FALSE,
                          estimate_af = FALSE,
                          output_dir = file.path("spectreasy_outputs", "unmix_samples", "unmixed_fcs"),
                          write_fcs = TRUE,
                          save_report = TRUE,
                          output_file = NULL,
                          save_qc_plots = FALSE,
                          qc_plot_dir = NULL,
                          subsample_n = NULL,
                          seed = NULL,
                          return_type = c("list", "flowSet", "SingleCellExperiment"),
                          verbose = TRUE) {
    return_type <- match.arg(return_type)
    .with_optional_seed(seed)
    variances_file_was_missing <- missing(variances_file)

    matrix_res <- .resolve_unmix_samples_matrix(
        M = M,
        unmixing_matrix_file = unmixing_matrix_file,
        variances_file = variances_file,
        variances_file_was_missing = variances_file_was_missing,
        detector_noise_file = detector_noise_file,
        scc_dir = scc_dir
    )
    M <- matrix_res$M
    variances_file <- matrix_res$variances_file

    run_options <- .normalize_unmix_samples_options(
        method = method,
        rwls_max_iter = rwls_max_iter,
        multithreading = multithreading,
        n_threads = n_threads
    )
    method_label <- run_options$method_label
    solver_method <- run_options$solver_method
    rwls_max_iter <- run_options$rwls_max_iter
    n_threads <- run_options$n_threads
    estimate_af <- isTRUE(estimate_af)
    spectral_variant_library <- .resolve_spectral_variant_library_for_unmixing(
        M = M,
        spectral_variant_library = spectral_variant_library,
        spectral_variant_library_file = spectral_variant_library_file,
        unmixing_matrix_file = unmixing_matrix_file
    )
    M <- .ensure_wls_variances(
        M = M,
        method = solver_method,
        variances_file = variances_file,
        scc_dir = scc_dir,
        control_file = control_file,
        af_n_bands = af_n_bands,
        af_auto_max_bands = af_auto_max_bands,
        af_min_cluster_events = af_min_cluster_events,
        af_min_cluster_proportion = af_min_cluster_proportion,
        af_n_bands_sensitivity = af_n_bands_sensitivity,
        exclude_af = exclude_af,
        cytometer = cytometer,
        seed = seed
    )

    sample_entries <- .prepare_unmix_samples_input(sample_dir)

    if (estimate_af) {
        if (isTRUE(verbose)) {
            message("  estimate_af = TRUE: estimating AF from stained samples.")
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

    results <- list()
    marker_source_all <- rownames(M)
    secondary_label_map <- .resolve_secondary_label_map(marker_source_all, sample_dir = sample_dir)
    .ensure_unmix_output_dir(output_dir = output_dir, write_fcs = write_fcs)

    progress_bar <- .start_unmix_samples_progress(length(sample_entries), verbose = verbose)
    on.exit(.close_unmix_samples_progress(progress_bar), add = TRUE)

    for (entry_i in seq_along(sample_entries)) {
        entry <- sample_entries[[entry_i]]
        results[[entry$sample_name]] <- .unmix_one_sample(
            entry = entry,
            M = M,
            method_label = method_label,
            rwls_max_iter = rwls_max_iter,
            multithreading = multithreading,
            n_threads = n_threads,
            spectral_variant_library = spectral_variant_library,
            spectral_variant_top_k = spectral_variant_top_k,
            spectral_variant_min_abundance = spectral_variant_min_abundance,
            write_fcs = write_fcs,
            output_dir = output_dir,
            secondary_label_map = secondary_label_map,
            subsample_n = subsample_n,
            verbose = verbose
        )
        .tick_unmix_samples_progress(progress_bar, entry_i)
    }

    class(results) <- c("spectreasy_unmixed_results", "list")
    attr(results, "method") <- method_label
    attr(results, "reference_matrix") <- M
    attr(results, "blind_af_info") <- attr(M, "blind_af_info")
    attr(results, "spectral_variant_library") <- spectral_variant_library
    attr(results, "qc_report_file") <- NULL
    attr(results, "qc_samples_dir") <- NULL
    attr(results, "qc_metrics_dir") <- NULL
    attr(results, "qc_plot_dir") <- NULL

    if (isTRUE(save_report)) {
        qc_samples_dir <- NULL
        if (is.null(output_file)) {
            qc_samples_dir <- .next_safe_output_dir(.default_unmix_samples_report_dir(output_dir))
            output_file <- file.path(qc_samples_dir, "qc_samples_report.pdf")
        } else {
            qc_samples_dir <- dirname(output_file)
        }
        report_res <- qc_samples(
            results = results,
            M = M,
            output_file = output_file,
            method = method_label,
            qc_metrics_dir = qc_samples_dir,
            qc_plot_dir = qc_plot_dir,
            save_qc_pngs = save_qc_plots
        )
        attr(results, "qc_report_file") <- report_res$output_file
        attr(results, "qc_samples_dir") <- qc_samples_dir
        attr(results, "qc_metrics_dir") <- qc_samples_dir
        attr(results, "qc_plot_dir") <- report_res$qc_plot_dir
    }

    if (identical(return_type, "flowSet")) {
        return(invisible(.unmixed_results_to_flowset(results)))
    }
    if (identical(return_type, "SingleCellExperiment")) {
        return(invisible(.unmixed_results_to_sce(results, sample_entries = sample_entries)))
    }

    invisible(results)
}
