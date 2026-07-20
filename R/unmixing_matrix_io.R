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
        desc_vals[idx] <- as.character(marker)
    }

    output_ff <- target_ff
    output_pd <- pd_new
    output_pd$desc <- desc_vals
    flowCore::parameters(output_ff) <- methods::new("AnnotatedDataFrame", data = output_pd)
    output_ff
}
