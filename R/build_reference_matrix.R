# Normalizes the control data.frame mapping.
# Reads the CSV file if a path is provided, validates that the required "filename" column is present,
# fills in missing columns with defaults, and cleans up string values by trimming whitespace.
# Returns a standardized data.frame or NULL/empty data.frame.
.normalize_build_reference_control_df <- function(control_df) {
    if (is.character(control_df) && length(control_df) == 1 && !is.na(control_df)) {
        if (!file.exists(control_df)) stop("control_df file not found: ", control_df)
        control_df <- utils::read.csv(control_df, stringsAsFactors = FALSE, check.names = FALSE)
    } else if (!is.null(control_df) && !is.data.frame(control_df)) {
        stop("control_df must be either a data.frame or a single CSV path.")
    }

    if (is.data.frame(control_df)) {
        if (!("filename" %in% colnames(control_df))) {
            stop("control_df is missing required column: filename")
        }
        if (!("fluorophore" %in% colnames(control_df))) control_df$fluorophore <- ""
        if (!("channel" %in% colnames(control_df))) control_df$channel <- ""
        if (!("control.type" %in% colnames(control_df))) control_df$control.type <- ""
        if (!("universal.negative" %in% colnames(control_df))) control_df$universal.negative <- ""

        control_df$filename <- trimws(as.character(control_df$filename))
        control_df$fluorophore <- trimws(as.character(control_df$fluorophore))
        control_df$channel <- trimws(as.character(control_df$channel))
        control_df$control.type <- tolower(trimws(as.character(control_df$control.type)))
        control_df$universal.negative <- trimws(as.character(control_df$universal.negative))
    }

    control_df
}

# Validates the Autofluorescence (AF) modeling parameters.
# Ensures that AF bank-size parameters and maximum events are valid.
# Returns a validated list containing these parameters.
.validate_build_reference_af_args <- function(af_n_bands,
                                              af_max_cells,
                                              af_bands_per_file = NULL,
                                              af_auto_max_bands = 100,
                                              af_min_cluster_events = 20,
                                              af_min_cluster_proportion = 0.005,
                                              af_n_bands_sensitivity = 1.5) {
    af_n_bands_raw <- af_n_bands[1]
    af_n_bands <- if (is.character(af_n_bands_raw) && identical(tolower(trimws(af_n_bands_raw)), "auto")) {
        "auto"
    } else {
        as.integer(af_n_bands_raw)
    }
    if (!identical(af_n_bands, "auto") && (!is.finite(af_n_bands) || is.na(af_n_bands) || af_n_bands < 1)) {
        stop("af_n_bands must be an integer >= 1 or \"auto\".")
    }

    if (is.null(af_bands_per_file) || length(af_bands_per_file) == 0 || is.na(af_bands_per_file[1])) {
        af_bands_per_file <- NA_integer_
    } else {
        af_bands_per_file <- as.integer(af_bands_per_file[1])
        if (!is.finite(af_bands_per_file) || is.na(af_bands_per_file) || af_bands_per_file < 1) {
            stop("af_bands_per_file must be an integer >= 1 when supplied.")
        }
    }

    af_max_cells <- as.integer(af_max_cells[1])
    if (!is.finite(af_max_cells) || is.na(af_max_cells) || af_max_cells < 100) {
        stop("af_max_cells must be an integer >= 100.")
    }

    af_auto_max_bands <- as.integer(af_auto_max_bands[1])
    if (!is.finite(af_auto_max_bands) || is.na(af_auto_max_bands) || af_auto_max_bands < 1) {
        stop("af_auto_max_bands must be an integer >= 1.")
    }

    af_min_cluster_events <- as.integer(af_min_cluster_events[1])
    if (!is.finite(af_min_cluster_events) || is.na(af_min_cluster_events) || af_min_cluster_events < 1) {
        stop("af_min_cluster_events must be an integer >= 1.")
    }

    af_min_cluster_proportion <- as.numeric(af_min_cluster_proportion[1])
    if (!is.finite(af_min_cluster_proportion) || is.na(af_min_cluster_proportion) ||
        af_min_cluster_proportion < 0 || af_min_cluster_proportion > 1) {
        stop("af_min_cluster_proportion must be a number between 0 and 1.")
    }

    if (is.null(af_n_bands_sensitivity)) {
        stop("af_n_bands_sensitivity must be a number between 0.1 and 5.")
    }
    af_n_bands_sensitivity <- as.numeric(af_n_bands_sensitivity[1])
    if (!is.finite(af_n_bands_sensitivity) || is.na(af_n_bands_sensitivity) ||
        af_n_bands_sensitivity < 0.1 || af_n_bands_sensitivity > 5) {
        stop("af_n_bands_sensitivity must be a number between 0.1 and 5.")
    }

    list(
        af_n_bands = af_n_bands,
        af_bands_per_file = af_bands_per_file,
        af_max_cells = af_max_cells,
        af_auto_max_bands = af_auto_max_bands,
        af_min_cluster_events = af_min_cluster_events,
        af_min_cluster_proportion = af_min_cluster_proportion,
        af_n_bands_sensitivity = af_n_bands_sensitivity
    )
}

.reference_min_af_cluster_size <- function(n_events,
                                           min_cluster_events = 20,
                                           min_cluster_proportion = 0.005) {
    max(
        as.integer(min_cluster_events),
        as.integer(ceiling(as.numeric(min_cluster_proportion) * n_events))
    )
}

.reference_som_grid_dims <- function(n_nodes) {
    n_nodes <- max(1L, as.integer(n_nodes))
    xdim <- max(1L, floor(sqrt(n_nodes)))
    ydim <- ceiling(n_nodes / xdim)
    c(xdim = as.integer(xdim), ydim = as.integer(ydim))
}

.reference_som_af_centers <- function(af_shape, n_nodes) {
    if (!requireNamespace("FlowSOM", quietly = TRUE)) {
        stop("FlowSOM is required for AF SOM extraction.", call. = FALSE)
    }
    dims <- .reference_som_grid_dims(n_nodes)
    map <- FlowSOM::SOM(
        as.matrix(af_shape),
        xdim = dims[["xdim"]],
        ydim = dims[["ydim"]],
        silent = TRUE
    )
    centers <- as.matrix(map$codes)
    cluster_sizes <- tabulate(map$mapping[, 1], nbins = nrow(centers))
    if (nrow(centers) > n_nodes) {
        centers <- centers[seq_len(n_nodes), , drop = FALSE]
        cluster_sizes <- cluster_sizes[seq_len(n_nodes)]
    }
    list(
        centers = centers,
        cluster_sizes = cluster_sizes,
        xdim = dims[["xdim"]],
        ydim = dims[["ydim"]]
    )
}

.reference_cosine_matrix <- function(A, B) {
    A <- as.matrix(A)
    B <- as.matrix(B)
    if (nrow(A) == 0 || nrow(B) == 0) {
        return(matrix(numeric(), nrow = nrow(A), ncol = nrow(B)))
    }
    A_norm <- sqrt(rowSums(A^2)) + 1e-9
    B_norm <- sqrt(rowSums(B^2)) + 1e-9
    (A %*% t(B)) / (A_norm %o% B_norm)
}

.reference_qc_af_spectra <- function(af_spectra,
                                     spectra,
                                     threshold = 0.99,
                                     remove = TRUE) {
    af_spectra <- as.matrix(af_spectra)
    spectra <- as.matrix(spectra)
    if (nrow(af_spectra) == 0 || nrow(spectra) == 0 || !isTRUE(remove)) {
        return(af_spectra)
    }
    sim <- .reference_cosine_matrix(af_spectra, spectra)
    max_sim <- apply(sim, 1, max, na.rm = TRUE)
    keep <- !is.finite(max_sim) | max_sim < threshold
    n_removed <- sum(!keep)
    if (n_removed > 0) {
        warning(
            n_removed,
            " autofluorescence spectrum/spectra were removed because their cosine similarity to a fluorophore spectrum was >= ",
            threshold,
            ".",
            call. = FALSE
        )
    }
    af_spectra[keep, , drop = FALSE]
}

.reference_deduplicate_af_spectra <- function(af_spectra,
                                              threshold = 0.99,
                                              deduplicate = FALSE) {
    af_spectra <- as.matrix(af_spectra)
    if (!isTRUE(deduplicate) || nrow(af_spectra) <= 1L) {
        return(af_spectra)
    }
    threshold <- as.numeric(threshold[1])
    if (!is.finite(threshold) || is.na(threshold) || threshold <= 0 || threshold > 1) {
        stop("af_deduplication_threshold must be a number > 0 and <= 1.", call. = FALSE)
    }

    keep <- integer()
    for (i in seq_len(nrow(af_spectra))) {
        if (!length(keep)) {
            keep <- i
            next
        }
        sim <- .reference_cosine_matrix(af_spectra[i, , drop = FALSE], af_spectra[keep, , drop = FALSE])
        if (!any(is.finite(sim) & sim >= threshold)) {
            keep <- c(keep, i)
        }
    }
    af_spectra[keep, , drop = FALSE]
}

.reference_filter_contaminant_events <- function(event_mat, spectra_mat, threshold = 0.99) {
    event_mat <- as.matrix(event_mat)
    spectra_mat <- as.matrix(spectra_mat)
    if (nrow(event_mat) == 0 || nrow(spectra_mat) == 0) {
        return(rep(TRUE, nrow(event_mat)))
    }
    event_norms <- sqrt(rowSums(event_mat^2)) + 1e-9
    keep <- rep(TRUE, nrow(event_mat))
    for (i in seq_len(nrow(spectra_mat))) {
        spec <- spectra_mat[i, ]
        spec_norm <- sqrt(sum(spec^2)) + 1e-9
        cosine <- as.numeric(event_mat %*% spec) / (event_norms * spec_norm)
        keep <- keep & (!is.finite(cosine) | cosine < threshold)
        if (!any(keep)) {
            break
        }
    }
    keep
}

.reference_unmix_no_af <- function(events, spectra) {
    spectra <- as.matrix(spectra)
    unmixing_matrix <- .reference_ols_unmixing_matrix(spectra)
    as.matrix(events) %*% t(unmixing_matrix)
}

.reference_assign_af_fluorophores <- function(events,
                                             spectra,
                                             af_spectra,
                                             af_assignment = "projection",
                                             return_details = FALSE) {
    .assign_af_candidates(
        Y = events,
        marker_M = spectra,
        af_M = af_spectra,
        af_assignment = af_assignment,
        return_details = return_details
    )
}

.reference_unmix_selected_af <- function(events, spectra, af_spectra, af_assignments) {
    events <- as.matrix(events)
    spectra <- as.matrix(spectra)
    af_spectra <- as.matrix(af_spectra)
    af_assignments <- as.integer(af_assignments)

    n_events <- nrow(events)
    n_fluors <- nrow(spectra)
    coeffs <- matrix(0, nrow = n_events, ncol = n_fluors + 1L)
    colnames(coeffs) <- c("AF", rownames(spectra))
    residuals <- matrix(0, nrow = n_events, ncol = ncol(spectra), dimnames = list(NULL, colnames(spectra)))
    proj_fluor <- matrix(0, nrow = n_events, ncol = ncol(spectra), dimnames = list(NULL, colnames(spectra)))

    for (af_id in sort(unique(af_assignments))) {
        idx <- which(af_assignments == af_id)
        if (!length(idx) || is.na(af_id) || af_id < 1L || af_id > nrow(af_spectra)) {
            next
        }
        combined <- rbind(af_spectra[af_id, , drop = FALSE], spectra)
        unmixing_matrix <- solve(combined %*% t(combined), combined)
        coeff_i <- events[idx, , drop = FALSE] %*% t(unmixing_matrix)
        coeffs[idx, ] <- coeff_i
        residuals[idx, ] <- events[idx, , drop = FALSE] - (coeff_i %*% combined)
        proj_fluor[idx, ] <- coeff_i[, seq.int(2L, n_fluors + 1L), drop = FALSE] %*% spectra
    }

    list(coefficients = coeffs, residuals = residuals, projected_fluor = proj_fluor)
}

.reference_refine_af_spectra <- function(events,
                                         spectra,
                                         af_spectra,
                                         problem_quantile = 0.99,
                                         af_assignment = "projection",
                                         remove_contaminants = TRUE,
                                         contaminant_threshold = 0.99) {
    events <- as.matrix(events)
    spectra <- as.matrix(spectra)
    af_spectra <- as.matrix(af_spectra)
    if (nrow(events) < 20 || nrow(spectra) == 0 || nrow(af_spectra) == 0) {
        return(list(signatures = af_spectra, info = list(enabled = TRUE, added = 0L, reason = "insufficient_input")))
    }
    if (!requireNamespace("FlowSOM", quietly = TRUE)) {
        stop("FlowSOM is required for AF refinement. Install FlowSOM or set af_refine = FALSE.", call. = FALSE)
    }

    problem_quantile <- as.numeric(problem_quantile[1])
    if (!is.finite(problem_quantile) || is.na(problem_quantile) || problem_quantile <= 0 || problem_quantile >= 1) {
        stop("af_refine_problem_quantile must be a number between 0 and 1.", call. = FALSE)
    }

    af_assignments <- .reference_assign_af_fluorophores(
        events = events,
        spectra = spectra,
        af_spectra = af_spectra,
        af_assignment = af_assignment
    )
    fit <- .reference_unmix_selected_af(
        events = events,
        spectra = spectra,
        af_spectra = af_spectra,
        af_assignments = af_assignments
    )

    error <- fit$residuals + fit$projected_fluor
    fluor_coeffs <- fit$coefficients[, seq.int(2L, ncol(fit$coefficients)), drop = FALSE]
    error_magnitude <- if (ncol(fluor_coeffs) > 1) {
        sqrt(rowSums(fluor_coeffs^2))
    } else {
        abs(fluor_coeffs[, 1])
    }

    q <- problem_quantile
    repeat {
        threshold <- stats::quantile(error_magnitude, q, na.rm = TRUE, names = FALSE)
        problem_idx <- which(is.finite(error_magnitude) & error_magnitude > threshold)
        if (length(problem_idx) >= 500 || q < 0.5) {
            break
        }
        q <- q - 0.05
    }

    if (length(problem_idx) <= 10) {
        return(list(
            signatures = af_spectra,
            info = list(enabled = TRUE, added = 0L, problem_cells = length(problem_idx), problem_quantile = q, reason = "too_few_problem_cells")
        ))
    }

    af_abundance <- fit$coefficients[problem_idx, 1]
    af_abundance[!is.finite(af_abundance) | af_abundance == 0] <- 1e-6
    spill_ratios <- sweep(error[problem_idx, , drop = FALSE], 1, af_abundance, "/")
    complete_spill <- stats::complete.cases(spill_ratios)
    complete_problem_idx <- problem_idx[complete_spill]
    spill_ratios <- spill_ratios[complete_spill, , drop = FALSE]
    if (nrow(spill_ratios) <= 10) {
        return(list(
            signatures = af_spectra,
            info = list(enabled = TRUE, added = 0L, problem_cells = length(problem_idx), problem_quantile = q, reason = "too_few_complete_spill_ratios")
        ))
    }

    som_dim_error <- max(2L, floor(sqrt(nrow(spill_ratios) / 3)))
    map_error <- FlowSOM::SOM(
        spill_ratios,
        xdim = som_dim_error,
        ydim = som_dim_error,
        silent = TRUE
    )
    cluster_ids <- unique(map_error$mapping[, 1])

    modulated <- lapply(cluster_ids, function(cluster_id) {
        cluster_sub_idx <- which(map_error$mapping[, 1] == cluster_id)
        global_idx <- complete_problem_idx[cluster_sub_idx]
        if (!length(global_idx)) {
            return(NULL)
        }
        median_ratio <- apply(spill_ratios[cluster_sub_idx, , drop = FALSE], 2, stats::median, na.rm = TRUE)
        contributing_af <- unique(af_assignments[global_idx])
        out <- lapply(contributing_af, function(af_id) {
            if (is.na(af_id) || af_id < 1L || af_id > nrow(af_spectra)) {
                return(NULL)
            }
            updated <- af_spectra[af_id, ] * (1 + median_ratio)
            peak <- max(abs(updated), na.rm = TRUE)
            if (is.finite(peak) && peak > 1e-12) {
                updated <- updated / peak
            }
            updated
        })
        do.call(rbind, out[!vapply(out, is.null, logical(1))])
    })
    modulated <- do.call(rbind, modulated[!vapply(modulated, is.null, logical(1))])
    if (is.null(modulated) || nrow(modulated) == 0) {
        return(list(
            signatures = af_spectra,
            info = list(enabled = TRUE, added = 0L, problem_cells = length(problem_idx), problem_quantile = q, reason = "no_modulated_spectra")
        ))
    }
    modulated <- as.matrix(stats::na.omit(modulated))
    if (nrow(modulated) == 0) {
        return(list(
            signatures = af_spectra,
            info = list(enabled = TRUE, added = 0L, problem_cells = length(problem_idx), problem_quantile = q, reason = "no_complete_modulated_spectra")
        ))
    }

    combined <- rbind(af_spectra, modulated)
    combined <- .reference_qc_af_spectra(
        af_spectra = combined,
        spectra = spectra,
        threshold = contaminant_threshold,
        remove = remove_contaminants
    )
    rownames(combined) <- c("AF", if (nrow(combined) > 1L) paste0("AF_", seq.int(2L, nrow(combined))) else NULL)
    colnames(combined) <- colnames(af_spectra)

    list(
        signatures = combined,
        info = list(
            enabled = TRUE,
            added = max(0L, nrow(combined) - nrow(af_spectra)),
            problem_cells = length(problem_idx),
            problem_quantile = q,
            som_dim_error = som_dim_error,
            total_bands = nrow(combined)
        )
    )
}

.reference_cell_fsc_upper_fraction <- 0.98
.reference_cell_ssc_upper_fraction <- 0.98
.reference_unstained_ssc_fsc_ratio_max <- 1.25

# Helper to retrieve rows in the control mapping that match the given filenames.
# Used to query metadata (like fluorophore name or target channel) for specific FCS files.
# Returns a subset of the control mapping data.frame.
.get_control_rows_for_reference <- function(df, filenames) {
    if (is.null(df) || !("filename" %in% colnames(df))) {
        return(data.frame())
    }
    fn <- as.character(df$filename)
    df[fn %in% filenames, ]
}

# Filters out Autofluorescence (AF) or unstained control files from the main FCS files set.
# Used when exclude_af is TRUE to ensure only single-color controls (and not unstained files)
# are processed during reference matrix construction.
# Returns the filtered FCS file list, stopping if no non-AF files remain.
.filter_reference_af_files <- function(fcs_files, input_folder, exclude_af = FALSE) {
    if (!isTRUE(exclude_af)) {
        return(fcs_files)
    }

    keep_non_af_files <- !.is_af_filename(fcs_files)
    n_excluded_files <- sum(!keep_non_af_files)
    if (n_excluded_files > 0) {
        message("Excluding ", n_excluded_files, " AF/unstained SCC file(s) from reference-matrix construction.")
    }
    fcs_files <- fcs_files[keep_non_af_files]
    if (length(fcs_files) == 0) {
        stop("No non-AF FCS files found in ", input_folder, " after exclude_af filtering.")
    }

    fcs_files
}

# Prepares the complete list of FCS files to be processed.
# Locates FCS files in the input folder and applies AF filtering,
# and optionally appends additional AF files from an external directory (af_dir).
# Returns a list containing 'fcs_files' (standard controls) and 'fcs_files_all' (all controls including AF).
.prepare_reference_file_set <- function(input_folder, include_multi_af = FALSE, af_dir = "af", exclude_af = FALSE) {
    fcs_files <- list.files(input_folder, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
    if (length(fcs_files) == 0) stop("No FCS files found in ", input_folder)

    fcs_files <- .filter_reference_af_files(fcs_files, input_folder = input_folder, exclude_af = exclude_af)

    fcs_files_all <- fcs_files
    if (!isTRUE(exclude_af) && isTRUE(include_multi_af) && dir.exists(af_dir)) {
        af_files <- list.files(af_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
        message("Found ", length(af_files), " extra AF files in '", af_dir, "'")
        fcs_files_all <- c(af_files, fcs_files)
    }

    list(fcs_files = fcs_files, fcs_files_all = fcs_files_all)
}

# Prepares and returns the normalized path of the output directory.
# If save_qc_plots is TRUE, recursively creates subdirectories for FSC/SSC plots,
# intensity gates, and spectra.
# Returns the absolute path string.
.prepare_reference_output_path <- function(output_folder, save_qc_plots = FALSE) {
    out_path <- normalizePath(output_folder, mustWork = FALSE)
    if (isTRUE(save_qc_plots)) {
        dir.create(file.path(out_path, "fsc_ssc"), showWarnings = FALSE, recursive = TRUE)
        dir.create(file.path(out_path, "histogram"), showWarnings = FALSE, recursive = TRUE)
        dir.create(file.path(out_path, "intensity_scatter"), showWarnings = FALSE, recursive = TRUE)
        dir.create(file.path(out_path, "spectrum"), showWarnings = FALSE, recursive = TRUE)
    }
    out_path
}

# Helper to normalize a channel name by converting to uppercase, removing whitespace,
# and handling NA values.
# Returns the normalized character vector.
.normalize_reference_channel <- function(x) {
    out <- toupper(gsub("\\s+", "", trimws(as.character(x))))
    out <- gsub("([A-Z]+)-([0-9])", "\\1\\2", out, perl = TRUE)
    out[is.na(out)] <- ""
    out
}

# Resolves a fluorophore's target channel name to an actual detector in the FCS file.
# Normalizes string cases/whitespace and tries matching the channel name directly, with/without
# standard suffix/prefix options (e.g., "-A"), or using a cytometer-specific channel alias map.
# Returns the matched detector name string or an empty string if unresolved.
.resolve_reference_control_channel <- function(channel_value, det_names, channel_alias_map = character()) {
    if (is.null(channel_value) || is.na(channel_value) || trimws(channel_value) == "") return("")
    if (channel_value %in% det_names) return(channel_value)

    det_norm <- .normalize_reference_channel(det_names)
    name_by_norm <- stats::setNames(det_names, det_norm)
    key <- .normalize_reference_channel(channel_value)
    candidates <- unique(c(
        key,
        gsub("-A$", "", key),
        paste0(key, "-A"),
        if (length(channel_alias_map) > 0 && key %in% names(channel_alias_map)) channel_alias_map[[key]] else NULL
    ))
    candidates <- candidates[nzchar(candidates)]
    for (cand in candidates) {
        if (cand %in% det_norm) return(name_by_norm[[cand]])
    }
    ""
}

# Extracts detector information from the first FCS file in the set.
# Reads the file metadata to determine detector names, labels, and channel mappings.
# Returns a list containing the parameter data table (pd_meta), sorted detector info, and alias map.
.prepare_reference_detector_info <- function(first_fcs_file) {
    ff_meta <- flowCore::read.FCS(first_fcs_file, transformation = FALSE, truncate_max_range = FALSE)
    pd_meta <- flowCore::pData(flowCore::parameters(ff_meta))
    det_info <- get_sorted_detectors(pd_meta)

    list(
        pd_meta = pd_meta,
        det_info = det_info,
        detector_names = det_info$names,
        detector_labels = det_info$labels,
        channel_alias_map = .build_channel_alias_map_from_pd(pd_meta)
    )
}

# Automatically identifies the primary Forward Scatter (FSC) and Side Scatter (SSC) channels.
# Matches names in the raw data columns using standard patterns for FSC and SSC.
# Returns a list with 'fsc' and 'ssc' channel names, or NULL if they cannot be resolved.
.resolve_reference_scatter_channels <- function(raw_data) {
    primary <- .get_primary_scatter_channels(colnames(raw_data))

    if (is.na(primary$fsc) || is.na(primary$ssc) || !all(c(primary$fsc, primary$ssc) %in% colnames(raw_data))) {
        return(NULL)
    }

    primary
}

# Computes the 2D ellipse coordinates at a given confidence level.
# Used for gating populations in FSC/SSC scatter plots based on GMM component mean and variance.
# Returns a data.table containing the ellipse points.
.get_reference_ellipse <- function(mean, sigma, level = 0.95, n = 100, scale = 1.0) {
    chi2_val <- qchisq(level, df = 2)
    eig <- eigen(sigma)
    a <- sqrt(eig$values[1] * chi2_val) * scale
    b <- sqrt(eig$values[2] * chi2_val) * scale
    angle <- atan2(eig$vectors[2, 1], eig$vectors[1, 1])
    theta <- seq(0, 2 * pi, length.out = n)
    ellipse_x <- a * cos(theta)
    ellipse_y <- b * sin(theta)
    rot_x <- ellipse_x * cos(angle) - ellipse_y * sin(angle)
    rot_y <- ellipse_x * sin(angle) + ellipse_y * cos(angle)
    data.table::data.table(x = rot_x + mean[1], y = rot_y + mean[2])
}

# Fits a Gaussian Mixture Model (GMM) on scatter data using mclust.
# Tests components from 1 to max_k to automatically cluster cells or beads based on FSC/SSC density.
# Returns a list with the model fit, mixture proportions, cluster means, covariances, and main components.
.fit_reference_gmm_populations <- function(data, max_k = 5, min_prop = 0.05) {
    mclustBIC <- get("mclustBIC", envir = asNamespace("mclust"))
    fit <- mclust::Mclust(data, G = seq_len(max_k), verbose = FALSE)
    if (is.null(fit)) {
        return(NULL)
    }
    sigmas <- lapply(seq_len(fit$G), function(k) {
        if (fit$modelName %in% c("EII", "VII")) {
            diag(fit$parameters$variance$sigmasq[k], nrow = 2)
        } else {
            fit$parameters$variance$sigma[, , k]
        }
    })
    list(
        fit = fit,
        proportions = fit$parameters$pro,
        means = fit$parameters$mean,
        sigmas = sigmas,
        main_populations = which(fit$parameters$pro >= min_prop)
    )
}

# Resolves sample type (cells/beads/unstained) from the FCS filename based on pre-defined regex patterns.
# Returns a list with the detected type and matching pattern.
.get_reference_sample_type <- function(filename, patterns, default) {
    for (type in names(patterns)) {
        pats <- patterns[[type]][order(-nchar(patterns[[type]]))]
        for (p in pats) {
            if (grepl(p, filename, fixed = FALSE, ignore.case = TRUE)) {
                return(list(type = type, pattern = p))
            }
        }
    }
    list(type = default, pattern = "default")
}

# Resolves control sample type, giving priority to control_df metadata over filename heuristics.
# Returns the final sample type configuration.
.resolve_reference_sample_type <- function(filename, row_info, patterns, default) {
    sample_info <- .get_reference_sample_type(filename, patterns, default)
    if (nrow(row_info) > 0 && "control.type" %in% colnames(row_info)) {
        control_type <- tolower(trimws(as.character(row_info$control.type[1])))
        if (control_type %in% c("beads", "cells")) {
            sample_info$type <- control_type
        }
    }
    sample_info
}

# Selects the most likely bead population from GMM components.
# Uses local density in a normalized FSC/SSC space to identify the actual peak/bead population
# rather than debris or background.
# Returns a list with the index of the selected population.
.select_reference_bead_population <- function(gmm_result) {
    if (length(gmm_result$main_populations) == 0) {
        return(NULL)
    }
    
    # Calculate local density of raw data around each GMM mean to ensure we select the actual peak
    data_fit <- gmm_result$fit$data
    means <- gmm_result$means
    
    sd_fsc <- stats::sd(data_fit[, 1])
    sd_ssc <- stats::sd(data_fit[, 2])
    
    # Avoid division by zero
    if (is.na(sd_fsc) || sd_fsc == 0) sd_fsc <- 1
    if (is.na(sd_ssc) || sd_ssc == 0) sd_ssc <- 1
    
    # Normalized data
    data_norm <- data_fit
    data_norm[, 1] <- data_fit[, 1] / sd_fsc
    data_norm[, 2] <- data_fit[, 2] / sd_ssc
    
    # Compute local density within r = 0.10 in normalized space
    local_densities <- sapply(gmm_result$main_populations, function(k) {
        mu_norm <- c(means[1, k] / sd_fsc, means[2, k] / sd_ssc)
        dists <- sqrt((data_norm[, 1] - mu_norm[1])^2 + (data_norm[, 2] - mu_norm[2])^2)
        sum(dists <= 0.10)
    })
    
    best <- gmm_result$main_populations[which.max(local_densities)]
    list(selected = best)
}

# Selects valid cell populations from GMM components.
# Filters components using minimum and maximum FSC/SSC thresholds and ratio boundaries
# to exclude debris and keep only the intact single-cell population.
# Returns a list with the indices of all valid cell populations.
.select_reference_cell_populations <- function(gmm_result,
                                               fsc_min,
                                               fsc_max = Inf,
                                               ssc_max = Inf,
                                               ratio_max = Inf) {
    valid <- c()
    means <- gmm_result$means
    fsc <- means[1, ]
    ssc <- means[2, ]
    ratio <- ssc / pmax(fsc, 1)
    if (!is.finite(fsc_max)) fsc_max <- Inf
    if (!is.finite(ssc_max)) ssc_max <- Inf
    if (!is.finite(ratio_max)) ratio_max <- Inf

    for (k in gmm_result$main_populations) {
        if (fsc[k] >= fsc_min &&
            fsc[k] <= fsc_max &&
            ssc[k] <= ssc_max &&
            ratio[k] <= ratio_max) {
            valid <- c(valid, k)
        }
    }
    list(selected = valid)
}

# Creates a polygon gate enclosing the selected GMM populations.
# Computes confidence ellipses at the specified level. If multiple populations are selected,
# computes the convex hull around the union of their ellipses. Optionally clips the gate boundaries.
# Returns a data.table of (x, y) coordinates representing the gate polygon.
.create_reference_merged_gate <- function(gmm_result, populations, level, scale = 1.0, clip_x = Inf, clip_y = Inf) {
    if (length(populations) == 0) {
        return(NULL)
    }
    if (length(populations) == 1) {
        ell <- .get_reference_ellipse(gmm_result$means[, populations], gmm_result$sigmas[[populations]], level, scale = scale)
    } else {
        all_pts <- data.table::rbindlist(
            lapply(populations, function(k) {
                .get_reference_ellipse(gmm_result$means[, k], gmm_result$sigmas[[k]], level, scale = scale)
            })
        )
        hull_idx <- grDevices::chull(all_pts$x, all_pts$y)
        ell <- all_pts[hull_idx, ]
    }
    ell$x <- pmax(0, pmin(ell$x, clip_x))
    ell$y <- pmax(0, pmin(ell$y, clip_y))
    if (nrow(ell) > 0 && (ell$x[1] != ell$x[nrow(ell)] || ell$y[1] != ell$y[nrow(ell)])) {
        ell <- rbind(ell, ell[1, , drop = FALSE])
    }
    ell
}

# Isolates the peak population from a 1D density distribution of log-transformed values.
# Detects peaks and troughs in the density curve, selects the rightmost peak exceeding
# the height cutoff, and extracts values bounded by the adjacent troughs.
# Returns the vector of values corresponding to this peak population.
.select_reference_hist_peak_population <- function(vals_log, rel_height_cutoff = 0.15, min_points = 20) {
    vals_log <- vals_log[is.finite(vals_log)]
    if (length(vals_log) < min_points) return(vals_log)

    d <- density(vals_log, n = 1024)
    x <- d$x
    y <- d$y

    peak_idx <- which(diff(sign(diff(y))) == -2) + 1
    trough_idx <- which(diff(sign(diff(y))) == 2) + 1
    if (length(peak_idx) == 0) return(vals_log)

    h_cut <- max(y, na.rm = TRUE) * rel_height_cutoff
    keep_peaks <- peak_idx[y[peak_idx] >= h_cut]
    if (length(keep_peaks) == 0) keep_peaks <- peak_idx
    sel_peak <- keep_peaks[which.max(x[keep_peaks])]

    left_trough <- trough_idx[trough_idx < sel_peak]
    right_trough <- trough_idx[trough_idx > sel_peak]
    left_x <- if (length(left_trough) > 0) x[max(left_trough)] else min(x, na.rm = TRUE)
    right_x <- if (length(right_trough) > 0) x[min(right_trough)] else max(x, na.rm = TRUE)

    peak_vals <- vals_log[vals_log >= left_x & vals_log <= right_x]
    if (length(peak_vals) < min_points) return(vals_log)
    peak_vals
}

# Extracts normalized Autofluorescence (AF) profiles from gated unstained events.
# Uses a FlowSOM grid on scale-free spectral shapes to build the multi-AF bank.
# Each SOM code is normalized by its peak detector to yield a relative signature from 0 to 1.
# Returns a list with the raw median spectrum and a matrix of normalized AF basis signatures.
.extract_reference_af_profiles <- function(ff_af = NULL,
                                           detector_names,
                                           n_bands = "auto",
                                           max_cells = 50000,
                                           af_events = NULL,
                                           auto_max_bands = 100,
                                           fluor_spectra = NULL,
                                           remove_contaminants = TRUE,
                                           contaminant_threshold = 0.99,
                                           af_deduplicate = FALSE,
                                           af_deduplication_threshold = 0.99,
                                           af_assignment = "projection",
                                           refine = FALSE,
                                           refine_problem_quantile = 0.99) {
    if (is.null(af_events)) {
        if (is.null(ff_af)) {
            stop("Either ff_af or af_events must be provided.")
        }
        raw <- flowCore::exprs(ff_af)
        af_events <- raw[, detector_names, drop = FALSE]
    }

    af_events <- af_events[stats::complete.cases(af_events), , drop = FALSE]
    if (nrow(af_events) == 0) {
        return(list(raw_median = NULL, signatures = NULL))
    }
    af_assignment <- .normalize_af_assignment(af_assignment, choices = c("projection", "residual_alignment"))

    if (nrow(af_events) > max_cells) {
        af_events <- af_events[sample.int(nrow(af_events), max_cells), , drop = FALSE]
    }

    raw_median <- apply(af_events, 2, stats::median, na.rm = TRUE)

    fluor_spectra <- if (!is.null(fluor_spectra)) {
        fluor_spectra <- as.matrix(fluor_spectra)
        fluor_spectra <- fluor_spectra[, detector_names, drop = FALSE]
        fluor_spectra[!grepl("^AF($|_)", rownames(fluor_spectra), ignore.case = TRUE), , drop = FALSE]
    } else {
        NULL
    }
    if (!is.null(fluor_spectra) && nrow(fluor_spectra) > 0 && isTRUE(remove_contaminants)) {
        af_centered <- sweep(af_events, 2, colMeans(af_events, na.rm = TRUE), "-")
        keep_events <- .reference_filter_contaminant_events(
            event_mat = af_centered,
            spectra_mat = fluor_spectra,
            threshold = contaminant_threshold
        )
        if (sum(keep_events, na.rm = TRUE) >= 100) {
            af_events <- af_events[keep_events, , drop = FALSE]
        }
    }

    af_pos <- pmax(af_events, 0)
    row_scale <- apply(af_pos, 1, max, na.rm = TRUE)
    keep <- is.finite(row_scale) & row_scale > 0
    if (!any(keep)) {
        stop("AF events do not contain positive spectral signal for SOM extraction.", call. = FALSE)
    }

    af_shape <- af_pos[keep, , drop = FALSE] / row_scale[keep]
    som_input <- af_shape
    if (!is.null(fluor_spectra) && nrow(fluor_spectra) > 0) {
        no_af <- .reference_unmix_no_af(af_events[keep, , drop = FALSE], fluor_spectra)
        colnames(no_af) <- rownames(fluor_spectra)
        som_input <- cbind(af_events[keep, detector_names, drop = FALSE], no_af)
    }
    selection <- NULL
    if (is.character(n_bands) && identical(tolower(trimws(n_bands[1])), "auto")) {
        n_bands <- auto_max_bands
        selection <- list(
            n_bands = as.integer(n_bands),
            method = "som_grid",
            max_bands = as.integer(auto_max_bands),
            hit_max_bands = TRUE
        )
    }
    n_eff <- min(as.integer(n_bands), nrow(af_shape))
    if (n_eff < 1) n_eff <- 1

    base_center_count <- NULL

    som <- .reference_som_af_centers(som_input, n_eff)
    centers <- som$centers[, detector_names, drop = FALSE]
    centers <- t(apply(centers, 1, function(v) {
        vmax <- max(abs(v), na.rm = TRUE)
        if (!is.finite(vmax) || vmax <= 0) {
            return(rep(0, length(v)))
        }
        v / vmax
    }))
    centers <- as.matrix(stats::na.omit(centers))
    if (!is.null(selection)) {
        selection$xdim <- som$xdim
        selection$ydim <- som$ydim
    }
    mean_af <- colMeans(centers, na.rm = TRUE)
    centers <- rbind(mean_af, centers)
    if (!is.null(fluor_spectra) && nrow(fluor_spectra) > 0) {
        centers <- .reference_qc_af_spectra(
            af_spectra = centers,
            spectra = fluor_spectra,
            threshold = contaminant_threshold,
            remove = remove_contaminants
        )
    }
    centers <- .reference_deduplicate_af_spectra(
        centers,
        threshold = af_deduplication_threshold,
        deduplicate = af_deduplicate
    )
    base_center_count <- nrow(centers)
    if (isTRUE(refine) && !is.null(fluor_spectra) && nrow(fluor_spectra) > 0) {
        refined <- .reference_refine_af_spectra(
            events = af_events[keep, detector_names, drop = FALSE],
            spectra = fluor_spectra,
            af_spectra = centers,
            problem_quantile = refine_problem_quantile,
            af_assignment = af_assignment,
            remove_contaminants = remove_contaminants,
            contaminant_threshold = contaminant_threshold
        )
        centers <- refined$signatures
        centers <- .reference_deduplicate_af_spectra(
            centers,
            threshold = af_deduplication_threshold,
            deduplicate = af_deduplicate
        )
        if (!is.null(selection)) {
            selection$refine <- refined$info
        }
    }

    if (is.null(dim(centers))) {
        centers <- matrix(centers, nrow = 1)
    }

    if (!is.null(selection)) {
        selection$raw_selected_bands <- selection$n_bands
        selection$n_bands <- nrow(centers)
        selection$final_bands <- nrow(centers)
        selection$raw_center_count <- if (!is.null(base_center_count)) base_center_count else nrow(centers)
        if (is.null(selection$refine)) {
            selection$refine <- list(enabled = isTRUE(refine), added = 0L)
        }
    }
    rownames(centers) <- c("AF", if (nrow(centers) > 1) paste0("AF_", seq.int(2, nrow(centers))) else NULL)
    colnames(centers) <- detector_names

    list(raw_median = raw_median, signatures = centers, selection = selection)
}

# Determines the name of the autofluorescence (unstained) control file.
# Looks it up in the control mapping or falls back to identifying files matching AF file naming heuristics.
# Returns the AF file basename (without extension) or NULL.
.resolve_reference_af_name <- function(control_df, fcs_files, exclude_af = FALSE) {
    if (isTRUE(exclude_af)) {
        return(NULL)
    }

    af_fn <- NULL
    if (!is.null(control_df)) {
        af_rows <- if ("fluorophore" %in% colnames(control_df)) {
            control_df[.is_af_control_row(
                fluorophore = control_df$fluorophore,
                marker = if ("marker" %in% colnames(control_df)) control_df$marker else NULL,
                filename = control_df$filename
            ), ]
        } else {
            data.frame()
        }
        if (nrow(af_rows) > 0) af_fn <- tools::file_path_sans_ext(basename(af_rows$filename[1]))
    }
    if (is.null(af_fn)) {
        af_idx_tmp <- which(.is_af_filename(fcs_files))
        if (length(af_idx_tmp) > 0) af_fn <- tools::file_path_sans_ext(basename(fcs_files[af_idx_tmp[1]]))
    }

    af_fn
}

# Safe retriever for key values from the configuration list.
# Returns the configured value if present, otherwise returns the specified default.
.get_reference_config_value <- function(config, name, default) {
    if (!is.null(config) && name %in% names(config)) {
        return(config[[name]])
    }
    default
}

# Reads and gates an AF / unstained control FCS file.
# Applies scatter-gating on FSC-SSC to isolate the main cell/bead population from debris.
# Returns a list containing the gated event data across all spectral detectors and gating metadata.
.extract_reference_af_gated_events <- function(fcs_file, detector_names, config) {
    ff_af <- tryCatch(
        flowCore::read.FCS(fcs_file, transformation = FALSE, truncate_max_range = FALSE),
        error = function(e) NULL
    )
    if (is.null(ff_af)) {
        warning("Could not read AF control file: ", fcs_file)
        return(NULL)
    }

    pd_af <- flowCore::pData(flowCore::parameters(ff_af))
    raw_af <- flowCore::exprs(ff_af)
    sn <- tools::file_path_sans_ext(basename(fcs_file))
    if (!.validate_reference_raw_data(raw_af, sn)) {
        return(NULL)
    }

    scatter_info <- .compute_reference_scatter_gate(
        raw_data = raw_af,
        pd = pd_af,
        sample_type = "unstained",
        outlier_percentile = .get_reference_config_value(config, "outlier_percentile", 0.02),
        debris_percentile = .get_reference_config_value(config, "debris_percentile", 0.08),
        subsample_n = .get_reference_config_value(config, "subsample_n", 5000),
        max_clusters = .get_reference_config_value(config, "max_clusters", 10),
        min_cluster_proportion = .get_reference_config_value(config, "min_cluster_proportion", 0.03),
        gate_contour_beads = .get_reference_config_value(config, "gate_contour_beads", 0.95),
        gate_contour_cells = .get_reference_config_value(config, "gate_contour_cells", 0.90),
        bead_gate_scale = .get_reference_config_value(config, "bead_gate_scale", 1.3)
    )
    if (is.null(scatter_info)) {
        message("  Skipping AF control ", basename(fcs_file), ": scatter gate found too few valid cells.")
        return(NULL)
    }

    list(
        events = scatter_info$gated_data[, detector_names, drop = FALSE],
        scatter = scatter_info$gated_data[, c(scatter_info$fsc, scatter_info$ssc), drop = FALSE],
        scatter_names = c(scatter_info$fsc, scatter_info$ssc),
        source = data.table::data.table(
            file = basename(fcs_file),
            path = normalizePath(fcs_file, mustWork = FALSE),
            n_total = nrow(raw_af),
            n_scatter_gated = nrow(scatter_info$gated_data),
            scatter_gate_pct = round(100 * nrow(scatter_info$gated_data) / max(nrow(raw_af), 1), 1),
            fsc_channel = scatter_info$fsc,
            ssc_channel = scatter_info$ssc
        )
    )
}

# Gathers all autofluorescence and unstained control files and extracts their signatures.
# Locates primary and extra AF files, loads and gates them, pools events, and computes
# the specified number of AF basis signatures.
# Returns a list containing raw medians, normalized basis matrices, and banking metadata.
.collect_reference_af_profiles <- function(control_df,
                                           fcs_files,
                                           detector_names,
                                           af_n_bands,
                                           af_bands_per_file,
                                           af_max_cells,
                                           af_auto_max_bands,
                                           af_deduplicate = FALSE,
                                           af_deduplication_threshold = 0.99,
                                           af_contaminant_threshold = 0.99,
                                           af_assignment = "projection",
                                           af_refine = FALSE,
                                           af_refine_problem_quantile = 0.99,
                                           exclude_af = FALSE,
                                           fcs_files_all = fcs_files,
                                           config = NULL,
                                           fluor_spectra = NULL,
                                           extract_signatures = TRUE) {
    af_data_raw <- NULL
    af_signatures_norm <- NULL
    af_bank_info <- NULL
    scc_background <- NULL
    af_fn <- .resolve_reference_af_name(control_df = control_df, fcs_files = fcs_files, exclude_af = exclude_af)

    if (isTRUE(exclude_af)) {
        return(list(
            af_data_raw = af_data_raw,
            af_signatures_norm = af_signatures_norm,
            af_bank_info = af_bank_info,
            scc_background = scc_background
        ))
    }

    af_paths <- character()
    af_source_types <- character()
    if (!is.null(af_fn)) {
        af_path <- fcs_files[grep(af_fn, fcs_files, fixed = TRUE)]
        if (length(af_path) > 0) {
            af_paths <- c(af_paths, af_path[1])
            af_source_types <- c(af_source_types, "primary_unstained")
        }
    }
    extra_af_paths <- fcs_files_all[vapply(
        fcs_files_all,
        .is_reference_extra_af_file,
        logical(1),
        include_multi_af = .get_reference_config_value(config, "include_multi_af", FALSE),
        af_dir = .get_reference_config_value(config, "af_dir", "af")
    )]
    af_paths <- c(af_paths, extra_af_paths)
    af_source_types <- c(af_source_types, rep("extra_af", length(extra_af_paths)))
    keep_unique <- !duplicated(normalizePath(af_paths, mustWork = FALSE))
    af_paths <- af_paths[keep_unique]
    af_source_types <- af_source_types[keep_unique]

    if (length(af_paths) > 0) {
        af_gated_list <- lapply(af_paths, .extract_reference_af_gated_events, detector_names = detector_names, config = config)
        keep_gated <- vapply(af_gated_list, function(x) !is.null(x) && !is.null(x$events) && nrow(x$events) > 0, logical(1))
        af_gated_list <- af_gated_list[keep_gated]
        af_source_types <- af_source_types[keep_gated]
        if (length(af_gated_list) == 0) {
            return(list(
                af_data_raw = af_data_raw,
                af_signatures_norm = af_signatures_norm,
                af_bank_info = af_bank_info,
                scc_background = scc_background
            ))
        }

        af_events <- do.call(rbind, lapply(af_gated_list, `[[`, "events"))
        scc_background <- .scc_background_from_gated_af_list(
            af_gated_list = af_gated_list,
            detector_names = detector_names
        )
        n_af_sources <- length(af_gated_list)
        requested_bands <- af_n_bands
        af_data_raw <- apply(af_events, 2, stats::median, na.rm = TRUE)
        af_profiles <- list(raw_median = af_data_raw, signatures = NULL, selection = NULL)
        if (isTRUE(extract_signatures)) {
            af_profiles <- .extract_reference_af_profiles(
                detector_names = detector_names,
                n_bands = requested_bands,
                max_cells = af_max_cells,
                af_events = af_events,
                auto_max_bands = af_auto_max_bands,
                fluor_spectra = fluor_spectra,
                contaminant_threshold = af_contaminant_threshold,
                af_deduplicate = af_deduplicate,
                af_deduplication_threshold = af_deduplication_threshold,
                af_assignment = af_assignment,
                refine = af_refine,
                refine_problem_quantile = af_refine_problem_quantile
            )
            af_data_raw <- af_profiles$raw_median
            af_signatures_norm <- af_profiles$signatures
        }
        af_sources <- data.table::rbindlist(lapply(af_gated_list, `[[`, "source"))
        af_sources$source_type <- af_source_types
        data.table::setcolorder(af_sources, c("file", "source_type", "n_total", "n_scatter_gated", "scatter_gate_pct", "fsc_channel", "ssc_channel", "path"))
        af_bank_info <- list(
            source_count = n_af_sources,
            sources = af_sources,
            pooled_events = nrow(af_events),
            requested_bands = requested_bands,
            af_bands_per_file = NA_integer_,
            derived_bands = if (!is.null(af_signatures_norm)) nrow(af_signatures_norm) else 0L,
            af_auto_max_bands = if (identical(requested_bands, "auto")) af_auto_max_bands else NA_integer_,
            af_deduplicate = isTRUE(af_deduplicate),
            af_deduplication_threshold = if (isTRUE(af_deduplicate)) af_deduplication_threshold else NA_real_,
            af_contaminant_threshold = af_contaminant_threshold,
            af_assignment = af_assignment,
            af_refine = isTRUE(af_refine),
            af_refine_problem_quantile = if (isTRUE(af_refine)) af_refine_problem_quantile else NA_real_,
            auto_selection = af_profiles$selection,
            mode = if (n_af_sources > 1) "pooled_af_sources" else "single_af"
        )
        if (!is.null(af_signatures_norm)) {
            msg <- if (n_af_sources == 1) {
                "primary unstained control"
            } else {
                paste0(n_af_sources, " pooled AF control files")
            }
            auto_msg <- if (!is.null(af_profiles$selection)) {
                paste0(" (auto-selected from ", af_profiles$selection$method, ")")
            } else {
                ""
            }
            message("Derived ", nrow(af_signatures_norm), " AF basis signature(s) from ", msg, auto_msg, ".")
            if (!is.null(af_profiles$selection) && isTRUE(af_profiles$selection$hit_max_bands)) {
                message(
                    "  Auto AF selection reached af_auto_max_bands = ",
                    af_profiles$selection$max_bands,
                    "; consider increasing it if QC improves with more AF bands."
                )
            }
        }
    }

    list(
        af_data_raw = af_data_raw,
        af_signatures_norm = af_signatures_norm,
        af_bank_info = af_bank_info,
        scc_background = scc_background
    )
}

# Checks if a file path points to an extra autofluorescence file.
# Specifically checks if the file resides in the designated AF directory and include_multi_af is enabled.
# Returns TRUE if it is an extra AF file, otherwise FALSE.
.is_reference_extra_af_file <- function(fcs_file, include_multi_af = FALSE, af_dir = "af") {
    if (!isTRUE(include_multi_af)) {
        return(FALSE)
    }
    grepl(
        normalizePath(af_dir, mustWork = FALSE),
        normalizePath(fcs_file, mustWork = FALSE),
        fixed = TRUE
    )
}

# Validates raw FCS data for corruption or extreme values.
# Ensures the dataset contains no infinite values, fewer than 10% NA values, and no extreme
# values exceeding 1e9 which indicate file corruption.
# Returns TRUE if valid, otherwise FALSE.
.validate_reference_raw_data <- function(raw_data, sn) {
    if (any(is.infinite(raw_data))) {
        message("  Skipping ", sn, ": Contains infinite values.")
        return(FALSE)
    }
    na_prop <- sum(is.na(raw_data)) / length(raw_data)
    if (na_prop > 0.1) {
        message("  Skipping ", sn, ": Too many NAs (", round(na_prop * 100, 1), "%).")
        return(FALSE)
    }
    max_val <- max(raw_data, na.rm = TRUE)
    if (max_val > 1e9) {
        message("  Skipping ", sn, ": Extreme values detected (max > 1e9). File may be corrupted.")
        return(FALSE)
    }
    TRUE
}

.validate_reference_detector_consistency <- function(fcs_files, detector_names) {
    if (length(fcs_files) == 0 || length(detector_names) == 0) {
        return(invisible(TRUE))
    }

    mismatches <- character()
    for (fcs_file in fcs_files) {
        ff <- tryCatch(
            flowCore::read.FCS(fcs_file, transformation = FALSE, truncate_max_range = FALSE),
            error = function(e) {
                stop("Could not read FCS file while checking detector consistency: ", fcs_file, call. = FALSE)
            }
        )
        pd <- flowCore::pData(flowCore::parameters(ff))
        current_detectors <- tryCatch(
            get_sorted_detectors(pd)$names,
            error = function(e) character()
        )
        if (length(current_detectors) == 0) {
            current_detectors <- intersect(detector_names, colnames(flowCore::exprs(ff)))
        }
        missing <- setdiff(detector_names, current_detectors)
        extra <- setdiff(current_detectors, detector_names)
        if (length(missing) > 0 || length(extra) > 0) {
            mismatches <- c(
                mismatches,
                paste0(
                    basename(fcs_file),
                    if (length(missing) > 0) paste0(" is missing: ", paste(missing, collapse = ", ")) else "",
                    if (length(missing) > 0 && length(extra) > 0) "; " else "",
                    if (length(extra) > 0) paste0("has extra detectors: ", paste(extra, collapse = ", ")) else ""
                )
            )
        }
    }

    if (length(mismatches) > 0) {
        stop(
            paste(
                c(
                    "Detector set mismatch across SCC/AF files.",
                    "All files used to build a reference matrix must contain the detectors found in the first SCC file.",
                    paste0(" - ", mismatches)
                ),
                collapse = "\n"
            ),
            call. = FALSE
        )
    }

    invisible(TRUE)
}

.active_reference_control_rows <- function(control_df, fcs_files_all, exclude_af = FALSE) {
    if (is.null(control_df) || !is.data.frame(control_df) || nrow(control_df) == 0 ||
        !all(c("filename", "fluorophore") %in% colnames(control_df))) {
        return(data.frame())
    }

    known_files <- basename(fcs_files_all)
    known_keys <- tools::file_path_sans_ext(known_files)
    control_keys <- tools::file_path_sans_ext(basename(as.character(control_df$filename)))
    active <- as.character(control_df$filename) %in% known_files | control_keys %in% known_keys
    is_af <- .is_af_control_row(
        fluorophore = control_df$fluorophore,
        marker = if ("marker" %in% colnames(control_df)) control_df$marker else NULL,
        filename = control_df$filename
    )
    if (isTRUE(exclude_af)) {
        active <- active & !is_af
    } else {
        active <- active & !is_af
    }

    control_df[active, , drop = FALSE]
}

.validate_reference_complete_controls <- function(control_df, fcs_files_all, processed_results, exclude_af = FALSE) {
    active_rows <- .active_reference_control_rows(
        control_df = control_df,
        fcs_files_all = fcs_files_all,
        exclude_af = exclude_af
    )
    if (nrow(active_rows) == 0) {
        return(invisible(TRUE))
    }

    processed_fluors <- if (length(processed_results) > 0) {
        trimws(as.character(vapply(processed_results, function(x) x$fluorophore[[1]], character(1))))
    } else {
        character()
    }
    expected_fluors <- trimws(as.character(active_rows$fluorophore))
    expected_fluors <- expected_fluors[nzchar(expected_fluors)]
    missing <- setdiff(unique(expected_fluors), unique(processed_fluors))
    if (length(missing) == 0) {
        return(invisible(TRUE))
    }

    missing_rows <- active_rows[trimws(as.character(active_rows$fluorophore)) %in% missing, , drop = FALSE]
    details <- paste0(
        as.character(missing_rows$filename),
        " -> ",
        as.character(missing_rows$fluorophore)
    )
    stop(
        paste(
            c(
                "Reference matrix construction did not produce spectra for all mapped non-AF controls.",
                "This usually means one or more controls failed reading, scatter gating, or histogram gating.",
                paste0(" - ", unique(details)),
                "Fix the listed controls or remove them from the control file before continuing."
            ),
            collapse = "\n"
        ),
        call. = FALSE
    )
}

# Computes a 2D scatter gate in FSC-SSC space for a sample.
# Filters out extreme outliers and debris, fits a Gaussian Mixture Model (GMM) to find clusters,
# selects the appropriate cell/bead population component(s), and builds a polygon gating boundary.
# Returns a list containing gated event data, gate coordinates, and scatter parameter bounds.
.compute_reference_scatter_gate <- function(raw_data,
                                            pd,
                                            sample_type,
                                            outlier_percentile,
                                            debris_percentile,
                                            subsample_n,
                                            max_clusters,
                                            min_cluster_proportion,
                                            gate_contour_beads,
                                            gate_contour_cells,
                                            bead_gate_scale) {
    scatter <- .resolve_reference_scatter_channels(raw_data)
    if (is.null(scatter)) {
        return(NULL)
    }
    fsc <- scatter$fsc
    ssc <- scatter$ssc

    data_raw_scatter <- raw_data[, c(fsc, ssc)]
    fsc_max <- quantile(data_raw_scatter[, 1], 1 - outlier_percentile, na.rm = TRUE)
    ssc_max <- quantile(data_raw_scatter[, 2], 1 - outlier_percentile, na.rm = TRUE)
    fsc_lower_limit <- debris_percentile * fsc_max
    fsc_upper_limit <- .reference_cell_fsc_upper_fraction * fsc_max
    ssc_upper_limit <- .reference_cell_ssc_upper_fraction * ssc_max
    valid_idx <- which(data_raw_scatter[, 1] < fsc_max & data_raw_scatter[, 2] < ssc_max & data_raw_scatter[, 1] > 0 & data_raw_scatter[, 2] > 0)
    data_filtered <- data_raw_scatter[valid_idx, ]
    if (sample_type %in% c("cells", "unstained")) {
        debris_threshold <- fsc_lower_limit
        data_filtered <- data_filtered[
            data_filtered[, 1] >= fsc_lower_limit &
                data_filtered[, 1] <= fsc_upper_limit &
                data_filtered[, 2] <= ssc_upper_limit,
        ]
    } else {
        debris_threshold <- 0
        fsc_upper_limit <- fsc_max
        ssc_upper_limit <- ssc_max
    }
    if (nrow(data_filtered) < 100) return(NULL)

    if (!is.null(subsample_n) && nrow(data_filtered) > subsample_n) {
        data_fit <- data_filtered[sample(nrow(data_filtered), subsample_n), ]
    } else {
        data_fit <- data_filtered
    }

    gmm_result <- .fit_reference_gmm_populations(data_fit, max_k = max_clusters, min_prop = min_cluster_proportion)
    if (is.null(gmm_result)) return(NULL)

    if (sample_type == "beads") {
        selected_pops <- .select_reference_bead_population(gmm_result)$selected
        gate_level <- gate_contour_beads
    } else {
        selected_pops <- .select_reference_cell_populations(
            gmm_result,
            fsc_min = debris_threshold,
            fsc_max = fsc_upper_limit,
            ssc_max = ssc_upper_limit,
            ratio_max = if (sample_type == "unstained") .reference_unstained_ssc_fsc_ratio_max else Inf
        )$selected
        gate_level <- gate_contour_cells
        if (length(selected_pops) == 0) selected_pops <- .select_reference_bead_population(gmm_result)$selected
    }
    if (length(selected_pops) == 0) return(NULL)

    final_gate <- .create_reference_merged_gate(
        gmm_result,
        selected_pops,
        gate_level,
        scale = if (sample_type == "beads") bead_gate_scale else 1.0,
        clip_x = if (sample_type %in% c("cells", "unstained")) fsc_upper_limit else Inf,
        clip_y = if (sample_type %in% c("cells", "unstained")) ssc_upper_limit else Inf
    )
    if (sample_type %in% c("cells", "unstained")) {
        final_gate$x <- pmax(fsc_lower_limit, final_gate$x)
        if (nrow(final_gate) > 0 && (final_gate$x[1] != final_gate$x[nrow(final_gate)] || final_gate$y[1] != final_gate$y[nrow(final_gate)])) {
            final_gate <- rbind(final_gate, final_gate[1, , drop = FALSE])
        }
    }

    inside_gate <- sp::point.in.polygon(raw_data[, fsc], raw_data[, ssc], final_gate$x, final_gate$y) > 0
    if (sample_type %in% c("cells", "unstained")) {
        inside_gate <- inside_gate &
            raw_data[, fsc] >= fsc_lower_limit &
            raw_data[, fsc] <= fsc_upper_limit &
            raw_data[, ssc] <= ssc_upper_limit
    }
    gated_data <- raw_data[inside_gate, ]
    if (nrow(gated_data) < 100) return(NULL)

    list(
        gated_data = gated_data,
        final_gate = final_gate,
        fsc = fsc,
        ssc = ssc,
        fsc_max = fsc_max,
        ssc_max = ssc_max
    )
}

# Identifies the peak/primary channel of high signal for single-color controls.
# For unstained controls, returns the channel with the highest median. For stained controls,
# infers the channel via the 99.9% quantile across all detectors, and cross-references/validates
# it against the target channel in metadata.
# Returns a list containing the resolved peak channel and 99.9% quantiles for all channels.
.select_reference_peak_channel <- function(gated_data,
                                           detector_names,
                                           row_info,
                                           channel_alias_map,
                                           sn_ext,
                                           sn,
                                           af_data_raw = NULL,
                                           clean_with_af = FALSE) {
    is_unstained <- grepl("unstained|autofluorescence|\\bAF\\b", paste(sn_ext, sn), ignore.case = TRUE)
    if (nrow(row_info) > 0) {
        is_unstained <- is_unstained || .is_af_control_row(
            fluorophore = if ("fluorophore" %in% colnames(row_info)) row_info$fluorophore[1] else "",
            marker = if ("marker" %in% colnames(row_info)) row_info$marker[1] else "",
            filename = sn_ext
        )
    }

    peak_matrix <- if (isTRUE(clean_with_af) && !is_unstained) {
        .scc_background_clean_peak_matrix(
            gated_data = gated_data,
            detector_names = detector_names,
            af_data_raw = af_data_raw
        )
    } else {
        gated_data[, detector_names, drop = FALSE]
    }

    q999_by_channel <- apply(
        peak_matrix,
        2,
        function(x) stats::quantile(x, 0.999, na.rm = TRUE)
    )
    if (is_unstained) {
        med_by_channel <- apply(gated_data[, detector_names, drop = FALSE], 2, stats::median, na.rm = TRUE)
        return(list(peak_channel = names(which.max(med_by_channel)), q999_by_channel = q999_by_channel))
    }

    inferred_peak_channel <- detector_names[which.max(q999_by_channel)]
    peak_channel <- inferred_peak_channel

    if (nrow(row_info) > 0 && !is.na(row_info$channel[1]) && row_info$channel[1] != "") {
        resolved_channel <- .resolve_reference_control_channel(row_info$channel[1], detector_names, channel_alias_map = channel_alias_map)
        if (nzchar(resolved_channel)) {
            inferred_val <- as.numeric(q999_by_channel[inferred_peak_channel])
            resolved_val <- as.numeric(q999_by_channel[resolved_channel])
            ranked <- names(sort(q999_by_channel, decreasing = TRUE))
            resolved_rank <- match(resolved_channel, ranked)
            use_inferred <- is.finite(inferred_val) &&
                is.finite(resolved_val) &&
                inferred_val > 0 &&
                !is.na(resolved_rank) &&
                resolved_rank > 8 &&
                (resolved_val / inferred_val) < 0.25

            if (use_inferred) {
                warning(
                    "Control channel '", row_info$channel[1], "' for ", sn_ext,
                    " appears inconsistent with signal profile (rank ",
                    resolved_rank, ", 99.9% ratio ",
                    round(resolved_val / inferred_val, 3),
                    "). Falling back to inferred channel ",
                    inferred_peak_channel, "."
                )
                peak_channel <- inferred_peak_channel
            } else {
                peak_channel <- resolved_channel
            }
        } else {
            warning("Control channel '", row_info$channel[1], "' for ", sn_ext, " not found in file. Falling back to inferred channel ", inferred_peak_channel, ".")
        }
    }

    list(peak_channel = peak_channel, q999_by_channel = q999_by_channel)
}

# Gates positive and negative populations on the log-transformed peak channel.
# Uses a Gaussian Mixture Model (GMM) or density peaks to robustly separate background/negative
# events from stained/positive events, allowing flexible gating directions and percentiles.
# Returns a list of indices/logic indicating which events fall in the negative and positive gates.
.compute_reference_histogram_gate <- function(peak_vals,
                                              sample_type,
                                              histogram_pct_beads,
                                              histogram_direction_beads,
                                              histogram_pct_cells,
                                              histogram_direction_cells,
                                              is_viability = FALSE) {
    vals_log <- log10(pmax(peak_vals, 1))
    vals_log <- vals_log[is.finite(vals_log)]
    positive_gate_present <- !(sample_type %in% c("unstained"))
    vals_for_gate <- if (sample_type %in% c("unstained", "cells")) vals_log else .select_reference_hist_peak_population(vals_log)
    pct <- if (sample_type %in% c("unstained", "cells")) histogram_pct_cells else histogram_pct_beads
    dir <- if (sample_type %in% c("unstained", "cells")) histogram_direction_cells else histogram_direction_beads
    pct <- min(max(pct, 0.01), 0.999)

    neg_log_min <- min(vals_log, na.rm = TRUE)
    neg_log_max <- as.numeric(stats::quantile(vals_log, 0.15, na.rm = TRUE))
    neg_gate_method <- "density lower-tail fallback"
    negative_gate_present <- positive_gate_present
    negative_components <- NULL
    density_gate <- NULL

    posterior_boundary <- function(fit, left_k, right_k, left_mean, right_mean) {
        if (!is.finite(left_mean) || !is.finite(right_mean) || right_mean <= left_mean) {
            return(NA_real_)
        }
        grid <- seq(left_mean, right_mean, length.out = 512)
        pred <- tryCatch(stats::predict(fit, newdata = grid)$z, error = function(e) NULL)
        if (is.null(pred) || ncol(pred) < max(left_k, right_k)) {
            return(mean(c(left_mean, right_mean)))
        }
        diff_post <- pred[, left_k] - pred[, right_k]
        cross <- which(diff_post <= 0)
        if (length(cross) == 0) {
            return(mean(c(left_mean, right_mean)))
        }
        grid[min(cross)]
    }

    component_stats <- function(x, cls) {
        comp <- sort(unique(cls))
        do.call(rbind, lapply(comp, function(k) {
            xk <- x[cls == k]
            data.frame(
                component = k,
                n = length(xk),
                prop = length(xk) / length(x),
                mean = mean(xk, na.rm = TRUE),
                sd = stats::sd(xk, na.rm = TRUE),
                q005 = as.numeric(stats::quantile(xk, 0.005, na.rm = TRUE)),
                q995 = as.numeric(stats::quantile(xk, 0.995, na.rm = TRUE)),
                stringsAsFactors = FALSE
            )
        }))
    }

    if (positive_gate_present && length(vals_log) >= 80) {
        d <- stats::density(vals_log, n = 2048)
        peak_idx <- which(diff(sign(diff(d$y))) == -2) + 1
        trough_idx <- which(diff(sign(diff(d$y))) == 2) + 1
        if (length(peak_idx) > 0) {
            if (sample_type %in% c("cells") && isTRUE(is_viability)) {
                sig_peaks <- peak_idx[
                    d$y[peak_idx] >= max(d$y, na.rm = TRUE) * 0.03 &
                        d$x[peak_idx] > 0.75
                ]
                sig_peaks <- sig_peaks[order(d$x[sig_peaks])]
                if (length(sig_peaks) >= 2) {
                    neg_peak <- sig_peaks[length(sig_peaks) - 1L]
                    pos_peak <- sig_peaks[length(sig_peaks)]
                    boundary_candidates <- trough_idx[trough_idx > neg_peak & trough_idx < pos_peak]
                    boundary <- if (length(boundary_candidates) > 0) {
                        d$x[boundary_candidates[which.min(d$y[boundary_candidates])]]
                    } else {
                        mean(c(d$x[neg_peak], d$x[pos_peak]))
                    }
                    neg_left <- trough_idx[trough_idx < neg_peak]
                    neg_log_min <- if (length(neg_left) > 0) {
                        d$x[max(neg_left)]
                    } else {
                        as.numeric(stats::quantile(vals_log[vals_log <= boundary], 0.005, na.rm = TRUE))
                    }
                    neg_log_max <- boundary
                    pos_lower <- boundary
                    pos_upper <- as.numeric(stats::quantile(vals_log[vals_log >= boundary], 0.999, na.rm = TRUE))
                    if (is.finite(pos_lower) && is.finite(pos_upper) && pos_upper > pos_lower &&
                        is.finite(neg_log_min) && is.finite(neg_log_max) && neg_log_max > neg_log_min) {
                        neg_gate_method <- paste0(
                            "viability gate: negative low mode at ", round(d$x[neg_peak], 2),
                            "; positive high mode at ", round(d$x[pos_peak], 2)
                        )
                        density_gate <- list(
                            pos_lower = pos_lower,
                            pos_upper = pos_upper,
                            pos_peak = d$x[pos_peak],
                            neg_peak = d$x[neg_peak],
                            component_means = d$x[sig_peaks]
                        )
                    }
                }
            }

            if (sample_type %in% c("cells") && is.null(density_gate)) {
                dominant_peak <- peak_idx[which.max(d$y[peak_idx])]
                left <- trough_idx[trough_idx < dominant_peak]
                right <- trough_idx[trough_idx > dominant_peak]
                bright_dominant <- d$x[dominant_peak] >= stats::median(vals_log, na.rm = TRUE) &&
                    d$x[dominant_peak] >= max(vals_log, na.rm = TRUE) - 0.75
                if (bright_dominant) {
                    sig_left_peaks <- peak_idx[
                        peak_idx < dominant_peak &
                            d$y[peak_idx] >= max(d$y, na.rm = TRUE) * 0.03 &
                            d$x[peak_idx] > 0.75 &
                            d$x[peak_idx] <= d$x[dominant_peak] - 0.75
                    ]
                    neg_peak <- if (length(sig_left_peaks) > 0) sig_left_peaks[which.max(d$x[sig_left_peaks])] else NA_integer_
                    if (!is.na(neg_peak)) {
                        neg_group <- neg_peak
                        prev_peaks <- rev(sig_left_peaks[sig_left_peaks < min(neg_group)])
                        for (prev_peak in prev_peaks) {
                            between_troughs <- trough_idx[trough_idx > prev_peak & trough_idx < min(neg_group)]
                            if (length(between_troughs) == 0) break
                            valley_y <- min(d$y[between_troughs], na.rm = TRUE)
                            adjacent_peak_y <- min(d$y[c(prev_peak, min(neg_group))], na.rm = TRUE)
                            if (is.finite(valley_y) && is.finite(adjacent_peak_y) && valley_y >= adjacent_peak_y * 0.60) {
                                neg_group <- c(prev_peak, neg_group)
                            } else {
                                break
                            }
                        }
                        neg_group <- sort(neg_group)
                        neg_left <- trough_idx[trough_idx < min(neg_group)]
                        neg_right <- trough_idx[trough_idx > max(neg_group) & trough_idx < dominant_peak]
                        neg_log_max <- if (length(neg_right) > 0) d$x[min(neg_right)] else d$x[max(left)]
                        density_floor <- max(d$y[neg_group], na.rm = TRUE) * 0.10
                        left_drop <- which(d$x < d$x[min(neg_group)] & d$y <= density_floor)
                        neg_log_min <- if (length(left_drop) > 0) {
                            d$x[max(left_drop)]
                        } else if (length(neg_left) > 0) {
                            d$x[max(neg_left)]
                        } else {
                            as.numeric(stats::quantile(vals_log[vals_log <= neg_log_max], 0.005, na.rm = TRUE))
                        }
                        neg_gate_method <- paste0("density middle negative cell mode at ", paste(round(d$x[neg_group], 2), collapse = "/"))
                    }
                    pos_lower <- if (length(left) > 0) d$x[max(left)] else as.numeric(stats::quantile(vals_log, 0.005, na.rm = TRUE))
                    pos_upper <- as.numeric(stats::quantile(vals_log, 0.999, na.rm = TRUE))
                    if (is.na(neg_peak)) {
                        neg_gate_method <- paste0("density dominant bright cell peak at ", round(d$x[dominant_peak], 2))
                    }
                    density_gate <- list(
                        pos_lower = pos_lower,
                        pos_upper = pos_upper,
                        pos_peak = d$x[dominant_peak],
                        neg_peak = if (!is.na(neg_peak)) d$x[neg_peak] else NA_real_,
                        component_means = d$x[peak_idx[d$y[peak_idx] >= max(d$y, na.rm = TRUE) * 0.005]]
                    )
                } else {
                    neg_log_min <- if (length(left) > 0) d$x[max(left)] else as.numeric(stats::quantile(vals_log, 0.005, na.rm = TRUE))
                    neg_log_max <- if (length(right) > 0) d$x[min(right)] else as.numeric(stats::quantile(vals_log, 0.85, na.rm = TRUE))
                    pos_lower <- max(
                        neg_log_max,
                        as.numeric(stats::quantile(vals_log, 0.95, na.rm = TRUE)),
                        na.rm = TRUE
                    )
                    pos_upper <- as.numeric(stats::quantile(vals_log, 0.999, na.rm = TRUE))
                    if (is.finite(pos_lower) && is.finite(pos_upper) && pos_upper > pos_lower) {
                        neg_gate_method <- paste0("density dominant negative peak at ", round(d$x[dominant_peak], 2))
                        density_gate <- list(
                            pos_lower = pos_lower,
                            pos_upper = pos_upper,
                            pos_peak = d$x[peak_idx[which.max(d$x[peak_idx])]],
                            neg_peak = d$x[dominant_peak],
                            component_means = d$x[peak_idx[d$y[peak_idx] >= max(d$y, na.rm = TRUE) * 0.005]]
                        )
                    }
                }
            }

            keep <- d$y[peak_idx] >= max(d$y, na.rm = TRUE) * 0.03
            keep <- keep & d$x[peak_idx] > 0.75
            sig_peaks <- peak_idx[keep]
            if (is.null(density_gate) && length(sig_peaks) > 1) {
                grouped <- list()
                current <- sig_peaks[1]
                for (idx in sig_peaks[-1]) {
                    if ((d$x[idx] - d$x[current[length(current)]]) <= 0.35) {
                        current <- c(current, idx)
                    } else {
                        grouped[[length(grouped) + 1L]] <- current
                        current <- idx
                    }
                }
                grouped[[length(grouped) + 1L]] <- current
                sig_peaks <- as.integer(vapply(grouped, function(g) g[which.max(d$y[g])], numeric(1)))
            }
            if (is.null(density_gate) && length(sig_peaks) >= 1) {
                sig_peaks <- sig_peaks[order(d$x[sig_peaks])]
                pos_peak <- sig_peaks[length(sig_peaks)]
                neg_peak <- if (length(sig_peaks) >= 2) sig_peaks[length(sig_peaks) - 1L] else NA_integer_
                nearest_left_trough <- function(peak) {
                    left <- trough_idx[trough_idx < peak]
                    if (length(left) > 0) d$x[max(left)] else min(vals_log, na.rm = TRUE)
                }
                nearest_right_trough <- function(peak) {
                    right <- trough_idx[trough_idx > peak]
                    if (length(right) > 0) d$x[min(right)] else max(vals_log, na.rm = TRUE)
                }
                pos_vals <- vals_log[vals_log >= nearest_left_trough(pos_peak)]
                pos_upper <- as.numeric(stats::quantile(pos_vals[pos_vals >= d$x[pos_peak]], 0.999, na.rm = TRUE))
                pos_right_trough <- nearest_right_trough(pos_peak)
                if (is.finite(pos_right_trough) && pos_right_trough > d$x[pos_peak] && pos_right_trough < max(vals_log, na.rm = TRUE)) {
                    pos_upper <- min(pos_upper, pos_right_trough, na.rm = TRUE)
                }
                if (!is.na(neg_peak)) {
                    neg_log_min <- nearest_left_trough(neg_peak)
                    neg_log_max <- nearest_right_trough(neg_peak)
                    neg_gate_method <- paste0("density negative peak at ", round(d$x[neg_peak], 2))
                }
                pos_peak_val <- d$x[pos_peak]
                left_trough_val <- nearest_left_trough(pos_peak)
                max_val_log <- max(vals_log, na.rm = TRUE)
                is_cut_off <- (max_val_log - pos_peak_val) < 0.25

                if (dir == "right" && !is_cut_off) {
                    pos_lower <- pos_peak_val
                } else {
                    pos_lower <- left_trough_val
                }

                density_gate <- list(
                    pos_lower = pos_lower,
                    pos_upper = pos_upper,
                    pos_peak = pos_peak_val,
                    neg_peak = if (!is.na(neg_peak)) d$x[neg_peak] else NA_real_,
                    component_means = d$x[sig_peaks]
                )
            }
        }
    }

    if (positive_gate_present && is.null(density_gate) && length(vals_log) >= 80 && requireNamespace("mclust", quietly = TRUE)) {
        neg_fit <- tryCatch(
            mclust::Mclust(vals_log, G = seq_len(min(6, max(1, floor(length(vals_log) / 40)))), verbose = FALSE),
            error = function(e) NULL
        )
        if (!is.null(neg_fit) && !is.null(neg_fit$classification)) {
            neg_stats <- component_stats(vals_log, neg_fit$classification)
            negative_components <- neg_stats
            min_real_prop <- 0.02
            floor_component <- neg_stats$mean <= 0.1 | neg_stats$q995 <= 0.25
            high_artifact <- neg_stats$prop < min_real_prop & neg_stats$mean > stats::median(vals_log, na.rm = TRUE)
            real_idx <- which(!floor_component & !high_artifact & neg_stats$prop >= min_real_prop)
            if (length(real_idx) >= 2) {
                positive_row <- real_idx[which.max(neg_stats$mean[real_idx])]
                negative_candidates <- real_idx[neg_stats$mean[real_idx] < neg_stats$mean[positive_row]]
                if (length(negative_candidates) > 0) {
                    negative_row <- negative_candidates[which.max(neg_stats$mean[negative_candidates])]
                    negative_k <- neg_stats$component[negative_row]
                    positive_k <- neg_stats$component[positive_row]
                    left_rows <- which(neg_stats$mean < neg_stats$mean[negative_row])
                    left_bound <- neg_stats$q005[negative_row]
                    if (length(left_rows) > 0) {
                        left_row <- left_rows[which.max(neg_stats$mean[left_rows])]
                        left_cross <- posterior_boundary(
                            neg_fit,
                            neg_stats$component[left_row],
                            negative_k,
                            neg_stats$mean[left_row],
                            neg_stats$mean[negative_row]
                        )
                        if (is.finite(left_cross)) left_bound <- max(left_bound, left_cross)
                    }
                    right_bound <- neg_stats$q995[negative_row]
                    right_cross <- posterior_boundary(
                        neg_fit,
                        negative_k,
                        positive_k,
                        neg_stats$mean[negative_row],
                        neg_stats$mean[positive_row]
                    )
                    if (is.finite(right_cross)) right_bound <- min(right_bound, right_cross)
                    if (is.finite(left_bound) && is.finite(right_bound) && right_bound > left_bound) {
                        neg_log_min <- left_bound
                        neg_log_max <- right_bound
                        neg_gate_method <- paste0("GMM negative component ", negative_k, " left of positive component ", positive_k)
                    }
                }
            }
        }
    }

    fit_vals <- vals_for_gate[is.finite(vals_for_gate)]
    model_info <- list(method = "quantile fallback", components = negative_components)
    gmm_upper <- NA_real_
    if (positive_gate_present && length(fit_vals) >= 80 && requireNamespace("mclust", quietly = TRUE)) {
        fit <- tryCatch(
            mclust::Mclust(fit_vals, G = seq_len(min(6, max(1, floor(length(fit_vals) / 40)))), verbose = FALSE),
            error = function(e) NULL
        )
        if (!is.null(fit) && !is.null(fit$classification)) {
            stats_df <- component_stats(fit_vals, fit$classification)
            min_real_prop <- 0.02
            q995_fit <- as.numeric(stats::quantile(fit_vals, 0.995, na.rm = TRUE))
            q999_fit <- as.numeric(stats::quantile(fit_vals, 0.999, na.rm = TRUE))
            high_artifact <- (stats_df$prop < min_real_prop & stats_df$mean > stats::median(fit_vals, na.rm = TRUE)) |
                (stats_df$prop < 0.05 & stats_df$mean >= q995_fit) |
                (is.finite(stats_df$sd) & stats_df$sd < 0.005 & stats_df$mean >= q999_fit)
            real_idx <- which(!high_artifact & stats_df$prop >= min_real_prop)
            if (length(real_idx) == 0) {
                real_idx <- which(!high_artifact)
            }
            if (length(real_idx) > 0) {
                bright_row <- real_idx[which.max(stats_df$mean[real_idx])]
                artifact_rows <- which(stats_df$mean > stats_df$mean[bright_row] & high_artifact)
                if (length(artifact_rows) > 0) {
                    gmm_upper <- q995_fit
                }
                model_info <- list(
                    method = paste0(
                        "GMM bright component ", stats_df$component[bright_row],
                        if (length(artifact_rows) > 0) "; high-end artifact components excluded" else "; no separated high-end artifact detected",
                        "; negative: ", neg_gate_method
                    ),
                    components = if (!is.null(negative_components)) negative_components else stats_df
                )
            }
        }
    }

    if (dir == "right") {
        lq <- 0.5
        uq <- min(1, 0.5 + pct)
    } else if (dir == "left") {
        lq <- max(0, 0.5 - pct)
        uq <- 0.5
    } else {
        lq <- 0.5 - pct / 2
        uq <- 0.5 + pct / 2
    }
    q_gate <- c(
        lower = as.numeric(stats::quantile(vals_for_gate, max(0, lq), na.rm = TRUE)),
        upper = as.numeric(stats::quantile(vals_for_gate, min(1, uq), na.rm = TRUE))
    )
    if (!positive_gate_present) {
        lower_log <- min(vals_log, na.rm = TRUE)
        upper_log <- max(vals_log, na.rm = TRUE)
        negative_gate_present <- FALSE
        model_info <- list(
            method = "AF/unstained control: scatter-gated only; no positive histogram gate",
            components = negative_components
        )
    } else if (!is.null(density_gate) && is.finite(density_gate$pos_lower) && is.finite(density_gate$pos_upper) && density_gate$pos_upper > density_gate$pos_lower) {
        lower_log <- density_gate$pos_lower
        upper_log <- density_gate$pos_upper
        component_df <- data.frame(
            component = seq_along(density_gate$component_means),
            n = NA_integer_,
            prop = NA_real_,
            mean = density_gate$component_means,
            sd = NA_real_,
            q005 = NA_real_,
            q995 = NA_real_
        )
        model_info <- list(
            method = paste0(
                "density/GMM mode gate: positive peak at ", round(density_gate$pos_peak, 2),
                "; negative: ", neg_gate_method
            ),
            components = component_df
        )
    } else if (dir == "right") {
        max_val_log <- max(vals_for_gate, na.rm = TRUE)
        med_val_log <- stats::median(vals_for_gate, na.rm = TRUE)
        is_cut_off <- (max_val_log - med_val_log) < 0.25
        if (is_cut_off) {
            lower_log <- as.numeric(stats::quantile(vals_for_gate, max(0, 0.5 - pct / 2), na.rm = TRUE))
            upper_log <- min(as.numeric(stats::quantile(vals_for_gate, min(1, 0.5 + pct / 2), na.rm = TRUE)), if (is.finite(gmm_upper)) gmm_upper else Inf, na.rm = TRUE)
        } else {
            lower_log <- med_val_log
            upper_log <- if (is.finite(gmm_upper)) gmm_upper else q_gate[["upper"]]
        }
    } else if (dir == "left") {
        lower_log <- q_gate[["lower"]]
        upper_log <- stats::median(vals_for_gate, na.rm = TRUE)
    } else {
        lower_log <- q_gate[["lower"]]
        upper_log <- min(q_gate[["upper"]], if (is.finite(gmm_upper)) gmm_upper else Inf, na.rm = TRUE)
    }

    retained_fraction <- mean(peak_vals >= 10^lower_log & peak_vals <= 10^upper_log, na.rm = TRUE)
    if (is.null(density_gate) && dir == "right" && (!is.finite(retained_fraction) || retained_fraction < 0.15)) {
        upper_log <- min(
            as.numeric(stats::quantile(fit_vals, 0.995, na.rm = TRUE)),
            q_gate[["upper"]],
            na.rm = TRUE
        )
    }
    if (!is.finite(lower_log) || !is.finite(upper_log) || upper_log <= lower_log) {
        lower_log <- q_gate[["lower"]]
        upper_log <- q_gate[["upper"]]
    }

    attr(vals_log, "neg_log_min") <- neg_log_min
    attr(vals_log, "neg_log_max") <- neg_log_max
    attr(vals_log, "negative_gate_present") <- negative_gate_present
    attr(vals_log, "positive_gate_present") <- positive_gate_present
    attr(vals_log, "gate_method") <- model_info$method
    attr(vals_log, "gmm_components") <- model_info$components
    attr(vals_log, "gate_type") <- "histogram"

    list(vals_log = vals_log, gate_min = 10^lower_log, gate_max = 10^upper_log)
}

.group_reference_density_peaks <- function(x, y, peak_idx, merge_distance = 0.35) {
    if (length(peak_idx) == 0) return(integer())
    peak_idx <- peak_idx[order(x[peak_idx])]
    groups <- list()
    current <- peak_idx[[1]]
    if (length(peak_idx) > 1) {
        for (idx in peak_idx[-1]) {
            if ((x[idx] - x[current[length(current)]]) <= merge_distance) {
                current <- c(current, idx)
            } else {
                groups[[length(groups) + 1L]] <- current
                current <- idx
            }
        }
    }
    groups[[length(groups) + 1L]] <- current
    as.integer(vapply(groups, function(g) g[which.max(y[g])], numeric(1)))
}

.compute_reference_density_scatter_gate <- function(vals_log,
                                                    min_prop = 0.02,
                                                    min_mean_gap = 0.35,
                                                    min_signal_log = 0.5,
                                                    peak_height_fraction = 0.03,
                                                    merge_distance = 0.35) {
    vals_log <- vals_log[is.finite(vals_log)]
    if (length(vals_log) < 80) return(NULL)
    gate_vals <- vals_log[vals_log > min_signal_log]
    if (length(gate_vals) < 80) return(NULL)

    d <- stats::density(gate_vals, n = 2048)
    peak_idx <- which(diff(sign(diff(d$y))) == -2) + 1
    if (length(peak_idx) == 0) return(NULL)

    signal_region <- d$x > min_signal_log
    reference_height <- if (any(signal_region)) max(d$y[signal_region], na.rm = TRUE) else max(d$y, na.rm = TRUE)
    signal_peaks <- peak_idx[
        d$x[peak_idx] > min_signal_log &
            d$y[peak_idx] >= reference_height * peak_height_fraction
    ]
    signal_peaks <- .group_reference_density_peaks(d$x, d$y, signal_peaks, merge_distance = merge_distance)
    if (length(signal_peaks) < 2) return(NULL)

    centers <- sort(d$x[signal_peaks])
    boundaries <- (centers[-1] + centers[-length(centers)]) / 2
    assigned_gate <- findInterval(gate_vals, boundaries) + 1L

    stats_df <- do.call(rbind, lapply(seq_along(centers), function(i) {
        xi <- gate_vals[assigned_gate == i]
        if (length(xi) == 0) return(NULL)
        spread <- stats::mad(xi, constant = 1.4826, na.rm = TRUE)
        if (!is.finite(spread) || spread <= 0) spread <- stats::sd(xi, na.rm = TRUE)
        data.frame(
            component = i,
            n = length(xi),
            prop = length(xi) / length(vals_log),
            mean = mean(xi, na.rm = TRUE),
            peak = centers[[i]],
            sd = stats::sd(xi, na.rm = TRUE),
            spread = spread,
            q005 = as.numeric(stats::quantile(xi, 0.005, na.rm = TRUE, names = FALSE)),
            q01 = as.numeric(stats::quantile(xi, 0.01, na.rm = TRUE, names = FALSE)),
            q15 = as.numeric(stats::quantile(xi, 0.15, na.rm = TRUE, names = FALSE)),
            q20 = as.numeric(stats::quantile(xi, 0.20, na.rm = TRUE, names = FALSE)),
            q85 = as.numeric(stats::quantile(xi, 0.85, na.rm = TRUE, names = FALSE)),
            q90 = as.numeric(stats::quantile(xi, 0.90, na.rm = TRUE, names = FALSE)),
            q99 = as.numeric(stats::quantile(xi, 0.99, na.rm = TRUE, names = FALSE)),
            q995 = as.numeric(stats::quantile(xi, 0.995, na.rm = TRUE, names = FALSE)),
            stringsAsFactors = FALSE
        )
    }))
    if (is.null(stats_df) || nrow(stats_df) < 2) return(NULL)

    real_idx <- which(stats_df$prop >= min_prop)
    if (length(real_idx) < 2) return(NULL)

    pos_idx <- real_idx[which.max(stats_df$peak[real_idx])]
    neg_candidates <- real_idx[stats_df$peak[real_idx] < (stats_df$peak[pos_idx] - min_mean_gap)]
    if (length(neg_candidates) == 0) return(NULL)
    neg_idx <- neg_candidates[which.max(stats_df$peak[neg_candidates])]

    window_for <- function(i, kind = c("negative", "positive")) {
        kind <- match.arg(kind)
        xi <- gate_vals[assigned_gate == i]
        row <- stats_df[i, ]
        if (identical(kind, "negative")) {
            lo <- row$q20
            hi <- row$q85
            spread_mult <- 2.0
        } else {
            lo <- row$q01
            hi <- row$q99
            spread_mult <- 2.75
        }
        if (is.finite(row$spread) && row$spread > 0) {
            lo <- max(lo, row$peak - spread_mult * row$spread)
            hi <- min(hi, row$peak + spread_mult * row$spread)
        }
        if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
            lo <- min(xi, na.rm = TRUE)
            hi <- max(xi, na.rm = TRUE)
        }
        c(lower = lo, upper = hi)
    }

    neg_window <- window_for(neg_idx, kind = "negative")
    pos_window <- window_for(pos_idx, kind = "positive")
    if ((pos_window[["lower"]] - neg_window[["upper"]]) < 0) {
        split <- mean(c(stats_df$peak[neg_idx], stats_df$peak[pos_idx]))
        neg_window[["upper"]] <- min(neg_window[["upper"]], split)
        pos_window[["lower"]] <- max(pos_window[["lower"]], split)
    }

    positive <- vals_log >= pos_window[["lower"]] & vals_log <= pos_window[["upper"]]
    negative <- vals_log >= neg_window[["lower"]] & vals_log <= neg_window[["upper"]]
    if (sum(positive, na.rm = TRUE) < 20 || sum(negative, na.rm = TRUE) < 20) {
        return(NULL)
    }

    list(
        mode = "density_peak_window_pos_neg",
        positive = positive,
        negative = negative,
        pos_lower = pos_window[["lower"]],
        pos_upper = pos_window[["upper"]],
        neg_lower = neg_window[["lower"]],
        neg_upper = neg_window[["upper"]],
        components = stats_df[, c("component", "n", "prop", "mean", "sd", "q005", "q995")]
    )
}

.compute_reference_gmm_scatter_gate <- function(vals_log,
                                                min_prop = 0.02,
                                                min_mean_gap = 0.35,
                                                max_components = 6L) {
    if (!requireNamespace("mclust", quietly = TRUE)) return(NULL)
    max_g <- min(as.integer(max_components), max(1L, floor(length(vals_log) / 40)))
    fit <- tryCatch(
        mclust::Mclust(vals_log, G = seq_len(max_g), verbose = FALSE),
        error = function(e) NULL
    )
    if (is.null(fit) || is.null(fit$classification)) return(NULL)

    component_stats <- function(x, cls) {
        comp <- sort(unique(cls))
        do.call(rbind, lapply(comp, function(k) {
            xk <- x[cls == k]
            data.frame(
                component = k,
                n = length(xk),
                prop = length(xk) / length(x),
                mean = mean(xk, na.rm = TRUE),
                sd = stats::sd(xk, na.rm = TRUE),
                q005 = as.numeric(stats::quantile(xk, 0.005, na.rm = TRUE)),
                q995 = as.numeric(stats::quantile(xk, 0.995, na.rm = TRUE)),
                stringsAsFactors = FALSE
            )
        }))
    }

    posterior_boundary <- function(left_k, right_k, left_mean, right_mean) {
        grid <- seq(left_mean, right_mean, length.out = 1024)
        pred <- tryCatch(stats::predict(fit, newdata = grid)$z, error = function(e) NULL)
        if (is.null(pred) || ncol(pred) < max(left_k, right_k)) return(mean(c(left_mean, right_mean)))
        diff_post <- pred[, left_k] - pred[, right_k]
        cross <- which(diff_post <= 0)
        if (length(cross) == 0) return(mean(c(left_mean, right_mean)))
        grid[min(cross)]
    }

    stats_df <- component_stats(vals_log, fit$classification)
    if (nrow(stats_df) < 2) return(NULL)

    floor_component <- stats_df$mean <= 0.1 | stats_df$q995 <= 0.25
    high_artifact <- stats_df$prop < min_prop & stats_df$mean > stats::median(vals_log, na.rm = TRUE)
    real_idx <- which(!floor_component & !high_artifact & stats_df$prop >= min_prop)
    if (length(real_idx) == 0) real_idx <- which(!floor_component & !high_artifact)
    if (length(real_idx) < 2) return(NULL)

    pos_idx <- real_idx[which.max(stats_df$mean[real_idx])]
    neg_candidates <- real_idx[stats_df$mean[real_idx] < stats_df$mean[pos_idx]]
    if (length(neg_candidates) == 0) return(NULL)
    neg_idx <- neg_candidates[which.max(stats_df$mean[neg_candidates])]
    neg_row <- stats_df[neg_idx, ]
    pos_row <- stats_df[pos_idx, ]
    if ((pos_row$mean - neg_row$mean) < min_mean_gap) return(NULL)

    boundary <- posterior_boundary(
        left_k = neg_row$component,
        right_k = pos_row$component,
        left_mean = neg_row$mean,
        right_mean = pos_row$mean
    )
    if (!is.finite(boundary) || boundary <= neg_row$mean || boundary >= pos_row$mean) {
        boundary <- mean(c(neg_row$mean, pos_row$mean))
    }

    left_bound <- neg_row$q005
    left_candidates <- which(stats_df$mean < neg_row$mean)
    if (length(left_candidates) > 0) {
        left_idx <- left_candidates[which.max(stats_df$mean[left_candidates])]
        left_cross <- posterior_boundary(
            left_k = stats_df$component[left_idx],
            right_k = neg_row$component,
            left_mean = stats_df$mean[left_idx],
            right_mean = neg_row$mean
        )
        if (is.finite(left_cross)) left_bound <- max(left_bound, left_cross)
    }

    list(
        mode = if (nrow(stats_df) == 2) "two_component_pos_neg" else "multi_component_nearest_negative_pos",
        positive = vals_log >= boundary & vals_log <= pos_row$q995,
        negative = vals_log >= left_bound & vals_log < boundary,
        pos_lower = boundary,
        pos_upper = pos_row$q995,
        neg_lower = left_bound,
        neg_upper = boundary,
        components = stats_df
    )
}

.compute_reference_one_cluster_scatter_gate <- function(vals_log) {
    signal_vals <- vals_log[is.finite(vals_log) & vals_log > 0.5]
    if (length(signal_vals) < 20) signal_vals <- vals_log[is.finite(vals_log)]
    pos_lower <- as.numeric(stats::quantile(signal_vals, 0.005, names = FALSE, na.rm = TRUE))
    pos_upper <- as.numeric(stats::quantile(signal_vals, 0.995, names = FALSE, na.rm = TRUE))
    positive <- vals_log >= pos_lower & vals_log <= pos_upper
    list(
        mode = "one_component_positive",
        positive = positive,
        negative = vals_log < pos_lower,
        pos_lower = pos_lower,
        pos_upper = pos_upper,
        neg_lower = min(vals_log, na.rm = TRUE),
        neg_upper = pos_lower,
        components = data.frame(
            component = 1L,
            n = length(vals_log),
            prop = 1,
            mean = mean(vals_log, na.rm = TRUE),
            sd = stats::sd(vals_log, na.rm = TRUE),
            q005 = as.numeric(stats::quantile(vals_log, 0.005, names = FALSE, na.rm = TRUE)),
            q995 = as.numeric(stats::quantile(vals_log, 0.995, names = FALSE, na.rm = TRUE))
        )
    )
}

.compute_reference_scatter_intensity_gate <- function(peak_vals,
                                                      sample_type,
                                                      histogram_pct_beads,
                                                      histogram_direction_beads,
                                                      histogram_pct_cells,
                                                      histogram_direction_cells,
                                                      is_viability = FALSE) {
    vals_log <- log10(pmax(peak_vals, 1))
    vals_log <- vals_log[is.finite(vals_log)]
    positive_gate_present <- !(sample_type %in% c("unstained"))
    if (!positive_gate_present) {
        out <- .compute_reference_histogram_gate(
            peak_vals = peak_vals,
            sample_type = sample_type,
            histogram_pct_beads = histogram_pct_beads,
            histogram_direction_beads = histogram_direction_beads,
            histogram_pct_cells = histogram_pct_cells,
            histogram_direction_cells = histogram_direction_cells,
            is_viability = is_viability
        )
        attr(out$vals_log, "gate_type") <- "scatter"
        return(out)
    }

    gate <- .compute_reference_density_scatter_gate(vals_log)
    if (is.null(gate)) {
        gate <- .compute_reference_gmm_scatter_gate(vals_log)
    }
    if (is.null(gate)) {
        gate <- .compute_reference_one_cluster_scatter_gate(vals_log)
    }

    if (!is.finite(gate$pos_lower) || !is.finite(gate$pos_upper) || gate$pos_upper <= gate$pos_lower) {
        out <- .compute_reference_histogram_gate(
            peak_vals = peak_vals,
            sample_type = sample_type,
            histogram_pct_beads = histogram_pct_beads,
            histogram_direction_beads = histogram_direction_beads,
            histogram_pct_cells = histogram_pct_cells,
            histogram_direction_cells = histogram_direction_cells,
            is_viability = is_viability
        )
        attr(out$vals_log, "gate_type") <- "histogram_fallback"
        return(out)
    }

    attr(vals_log, "neg_log_min") <- gate$neg_lower
    attr(vals_log, "neg_log_max") <- gate$neg_upper
    attr(vals_log, "negative_gate_present") <- isTRUE(sum(gate$negative, na.rm = TRUE) >= 10)
    attr(vals_log, "positive_gate_present") <- TRUE
    attr(vals_log, "gate_method") <- paste0("intensity scatter/GMM gate: ", gate$mode)
    attr(vals_log, "gmm_components") <- gate$components
    attr(vals_log, "gate_type") <- "scatter"

    list(vals_log = vals_log, gate_min = 10^gate$pos_lower, gate_max = 10^gate$pos_upper)
}

# Computes the normalized spectral signature for a stained control sample.
# Calculates the median intensity in each detector for the positive and negative gates,
# subtracts the negative/background control, normalizes the spectrum by the peak signal,
# and estimates the per-channel variance.
# Returns a list with the normalized spectrum vector, raw positive/negative medians, and variance vector.
.compute_reference_spectrum <- function(final_gated_data,
                                        gated_data,
                                        peak_vals,
                                        vals_log,
                                        detector_names,
                                        row_info,
                                        sample_type = "beads",
                                        af_data_raw = NULL,
                                        universal_negatives = NULL,
                                        scc_background = NULL,
                                        scc_background_method = "none",
                                        scc_background_k = 3L) {
    use_matched_background <- identical(scc_background_method, "scatter_knn") &&
        identical(sample_type, "cells") &&
        !is.null(scc_background)

    clean_pos <- if (use_matched_background) {
        .scc_background_clean_events(
            events = final_gated_data,
            detector_names = detector_names,
            background = scc_background,
            k = scc_background_k
        )
    } else {
        NULL
    }
    if (!is.null(clean_pos)) {
        pos_spectrum_raw <- apply(pmax(clean_pos, 0), 2, median, na.rm = TRUE)
        final_neg <- rep(0, length(detector_names))
        names(final_neg) <- detector_names
        neg_events <- matrix(0, nrow = 0, ncol = length(detector_names), dimnames = list(NULL, detector_names))
    } else {
        pos_spectrum_raw <- apply(final_gated_data[, detector_names, drop = FALSE], 2, median, na.rm = TRUE)
    }
    neg_log_min <- attr(vals_log, "neg_log_min")
    neg_log_max <- attr(vals_log, "neg_log_max")

    if (is.null(clean_pos)) {
        if (isTRUE(attr(vals_log, "negative_gate_present")) &&
            is.finite(neg_log_min) && is.finite(neg_log_max) && neg_log_max > neg_log_min) {
            neg_events <- gated_data[
                peak_vals >= 10^neg_log_min & peak_vals <= 10^neg_log_max,
                detector_names,
                drop = FALSE
            ]
        } else {
            neg_events <- gated_data[
                peak_vals <= 10^quantile(vals_log, 0.15, na.rm = TRUE),
                detector_names,
                drop = FALSE
            ]
        }
        if (nrow(neg_events) < 10) {
            neg_events <- gated_data[
                peak_vals <= 10^quantile(vals_log, 0.15, na.rm = TRUE),
                detector_names,
                drop = FALSE
            ]
        }
        neg_spectrum_raw <- apply(neg_events, 2, median, na.rm = TRUE)
        uv_val <- if (nrow(row_info) > 0 && "universal.negative" %in% colnames(row_info)) {
            trimws(as.character(row_info$universal.negative[1]))
        } else {
            ""
        }
        uv_upper <- toupper(uv_val)
        uv_key <- tools::file_path_sans_ext(basename(uv_val))
        use_af_negative <- uv_upper %in% c("TRUE", "AF")
        use_named_negative <- nzchar(uv_key) &&
            !uv_upper %in% c("FALSE", "TRUE", "AF") &&
            !is.null(universal_negatives) &&
            uv_key %in% names(universal_negatives)
        final_neg <- if (use_af_negative && !is.null(af_data_raw)) {
            af_data_raw
        } else if (use_named_negative) {
            universal_negatives[[uv_key]]
        } else {
            neg_spectrum_raw
        }
    }
    sig_pure <- pmax(pos_spectrum_raw - final_neg, 0)
    max_val <- max(sig_pure, na.rm = TRUE)
    if (max_val <= 0) max_val <- max(pos_spectrum_raw, na.rm = TRUE)
    res <- sig_pure / max_val
    positive_events <- if (!is.null(clean_pos)) {
        pmax(clean_pos, 0)
    } else {
        detector_mat <- as.matrix(final_gated_data[, detector_names, drop = FALSE])
        pmax(sweep(detector_mat, 2, final_neg, "-"), 0)
    }
    colnames(positive_events) <- detector_names

    # Keep positive/negative population spread as reference QC metadata. These
    # values are not used as default WLS detector-error weights.
    pos_var <- if (!is.null(clean_pos)) {
        apply(pmax(clean_pos, 0), 2, stats::var, na.rm = TRUE)
    } else {
        apply(final_gated_data[, detector_names, drop = FALSE], 2, stats::var, na.rm = TRUE)
    }
    neg_var <- if (nrow(neg_events) > 1L) apply(neg_events, 2, stats::var, na.rm = TRUE) else rep(0, length(detector_names))
    tot_var <- pos_var + neg_var
    tot_var[is.na(tot_var) | tot_var <= 0] <- 0
    if (max_val > 0) {
        tot_var <- tot_var / (max_val^2)
    }
    attr(res, "variance") <- tot_var
    attr(res, "scc_background") <- list(
        method = if (!is.null(clean_pos)) "scatter_knn" else "none",
        k = if (!is.null(clean_pos)) as.integer(scc_background_k) else NA_integer_,
        matched_events = if (!is.null(clean_pos)) nrow(clean_pos) else 0L
    )
    attr(res, "scc_positive_events") <- positive_events

    res
}

# Standardizes a file path to its base name without file extension for matching.
# Returns the cleaned character string/vector.
.reference_negative_key <- function(x) {
    tools::file_path_sans_ext(basename(trimws(as.character(x))))
}

# Identifies and loads the universal/marker-specific negative control profiles.
# If specified in the control mapping, reads the designated negative control files and computes
# their median spectral profiles to be subtracted from target samples.
# Returns a list of numeric vectors mapping fluorophore names to their negative control profiles.
.collect_reference_universal_negatives <- function(control_df,
                                                   fcs_files,
                                                   detector_names,
                                                   sample_patterns,
                                                   config) {
    out <- list()
    if (is.null(control_df) || !is.data.frame(control_df) || !("universal.negative" %in% colnames(control_df))) {
        return(out)
    }

    uv_vals <- trimws(as.character(control_df$universal.negative))
    uv_vals[is.na(uv_vals)] <- ""
    uv_vals <- unique(uv_vals[nzchar(uv_vals) & !toupper(uv_vals) %in% c("FALSE", "TRUE", "AF")])
    if (length(uv_vals) == 0) {
        return(out)
    }

    file_keys <- .reference_negative_key(fcs_files)
    names(fcs_files) <- file_keys

    for (uv in uv_vals) {
        key <- .reference_negative_key(uv)
        if (!nzchar(key) || !key %in% names(fcs_files)) {
            next
        }

        fcs_file <- unname(fcs_files[[key]])
        sn_ext <- basename(fcs_file)
        sn <- tools::file_path_sans_ext(sn_ext)
        row_info <- .get_control_rows_for_reference(control_df, c(sn_ext, sn))
        sample_info <- .resolve_reference_sample_type(
            filename = sn,
            row_info = row_info,
            patterns = sample_patterns,
            default = config$default_sample_type
        )

        ff <- tryCatch(
            flowCore::read.FCS(fcs_file, transformation = FALSE, truncate_max_range = FALSE),
            error = function(e) NULL
        )
        if (is.null(ff)) {
            warning("Could not read universal.negative file: ", fcs_file)
            next
        }

        pd <- flowCore::pData(flowCore::parameters(ff))
        raw_data <- flowCore::exprs(ff)
        scatter_info <- .compute_reference_scatter_gate(
            raw_data = raw_data,
            pd = pd,
            sample_type = sample_info$type,
            outlier_percentile = config$outlier_percentile,
            debris_percentile = config$debris_percentile,
            subsample_n = config$subsample_n,
            max_clusters = config$max_clusters,
            min_cluster_proportion = config$min_cluster_proportion,
            gate_contour_beads = config$gate_contour_beads,
            gate_contour_cells = config$gate_contour_cells,
            bead_gate_scale = config$bead_gate_scale
        )

        neg_data <- if (!is.null(scatter_info)) scatter_info$gated_data else raw_data
        out[[key]] <- apply(neg_data[, detector_names, drop = FALSE], 2, stats::median, na.rm = TRUE)
        message("Using universal negative file for SCC subtraction: ", basename(fcs_file))
    }

    out
}

# Helper to extract the description attribute of a channel for plotting.
# Falls back to the channel name if no description column or value is found.
# Returns the description string.
.get_reference_axis_label <- function(ch_name, pd_tbl) {
    idx <- match(ch_name, as.character(pd_tbl$name))
    if (!is.na(idx) && "desc" %in% colnames(pd_tbl)) {
        desc_val <- trimws(as.character(pd_tbl$desc[idx]))
        if (!is.na(desc_val) && nzchar(desc_val)) {
            return(desc_val)
        }
    }
    ch_name
}

# Generates and saves PDF quality control (QC) plots for a processed sample.
# Produces three plots: a 2D density plot of FSC-SSC gating, a 1D density histogram of positive/negative
# peak gating, and a spectral profile plot of the normalized signature.
# Writes files to the output folder.
.save_reference_qc_plots <- function(sn,
                                     raw_data,
                                     gated_data,
                                     final_gated_data,
                                     pd,
                                     fsc,
                                     ssc,
                                     fsc_max,
                                     ssc_max,
                                     final_gate,
                                     vals_log,
                                     peak_vals,
                                     gate_min,
                                     gate_max,
                                     peak_channel,
                                     detector_names,
                                     detector_labels,
                                     det_info,
                                     out_path,
                                     use_scatter_gating = TRUE) {
    fsc_desc <- .get_reference_axis_label(fsc, pd)
    ssc_desc <- .get_reference_axis_label(ssc, pd)
    dt_plot_gating <- data.table::data.table(x = raw_data[, fsc], y = raw_data[, ssc])
    x_breaks <- seq(0, max(fsc_max, ssc_max) * 1.05, length.out = 201)
    y_breaks <- x_breaks
    x_bin <- findInterval(dt_plot_gating$x, x_breaks)
    y_bin <- findInterval(dt_plot_gating$y, y_breaks)
    keep_bin <- x_bin >= 1 & x_bin <= 200 & y_bin >= 1 & y_bin <= 200
    dt2d <- data.table::as.data.table(as.data.frame(table(
        x_bin = x_bin[keep_bin],
        y_bin = y_bin[keep_bin]
    )))
    data.table::setnames(dt2d, c("x_bin", "y_bin", "count"))
    dt2d$x_bin <- as.integer(as.character(dt2d$x_bin))
    dt2d$y_bin <- as.integer(as.character(dt2d$y_bin))
    dt2d$count <- as.integer(as.character(dt2d$count))
    dt2d <- dt2d[dt2d$count > 0, ]
    dt2d$x <- x_breaks[dt2d$x_bin]
    dt2d$y <- y_breaks[dt2d$y_bin]
    dt2d$fill <- log10(dt2d$count + 1)
    p1 <- ggplot2::ggplot(dt2d, ggplot2::aes(x, y, fill = fill)) +
        ggplot2::geom_tile(width = diff(x_breaks)[1], height = diff(y_breaks)[1]) +
        ggplot2::scale_fill_gradientn(colors = c("#0000FF", "#00FFFF", "#00FF00", "#FFFF00", "#FF0000"), guide = "none") +
        ggplot2::geom_path(data = final_gate, ggplot2::aes(x, y), inherit.aes = FALSE, color = "red", linewidth = 1) +
        ggplot2::labs(title = paste0(sn, " - FSC/SSC"), subtitle = paste0(round(100 * nrow(gated_data) / nrow(raw_data), 1), "% gated"), x = fsc_desc, y = ssc_desc) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none", panel.grid = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = "white", color = NA), plot.subtitle = ggplot2::element_text(size = 10.6)) +
        ggplot2::coord_cartesian(xlim = c(0, max(fsc_max, ssc_max) * 1.05), ylim = c(0, max(fsc_max, ssc_max) * 1.05))
    ggplot2::ggsave(file.path(out_path, "fsc_ssc", paste0(sn, "_fsc_ssc.png")), p1, width = 5, height = 5, dpi = 300)

    neg_log_min <- attr(vals_log, "neg_log_min")
    neg_log_max <- attr(vals_log, "neg_log_max")
    gate_method <- attr(vals_log, "gate_method")
    negative_gate_present <- isTRUE(attr(vals_log, "negative_gate_present"))
    positive_gate_present <- isTRUE(attr(vals_log, "positive_gate_present"))
    comp <- attr(vals_log, "gmm_components")
    comp_means <- if (!is.null(comp) && nrow(comp) > 0) comp$mean else numeric()

    hist_subtitle_lines <- if (positive_gate_present) {
        c(
            paste0(
                round(100 * nrow(final_gated_data) / nrow(gated_data), 1),
                "% positive gated | blue = negative gate | red = bright gate"
            ),
            gate_method
        )
    } else {
        c("AF/unstained control: scatter-gated only; no positive histogram gate", gate_method)
    }
    hist_subtitle <- paste(
        unlist(lapply(hist_subtitle_lines, strwrap, width = 95), use.names = FALSE),
        collapse = "\n"
    )

    p2 <- ggplot2::ggplot(data.table::data.table(x = vals_log), ggplot2::aes(x)) +
        ggplot2::geom_density(fill = "grey80", color = "grey40") +
        ggplot2::labs(
            title = paste0(sn, " - ", peak_channel),
            subtitle = hist_subtitle,
            x = paste0("log10(", peak_channel, ")")
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none", plot.subtitle = ggplot2::element_text(size = 8.4, lineheight = 1.05))
    if (negative_gate_present && is.finite(neg_log_min) && is.finite(neg_log_max) && neg_log_max > neg_log_min) {
        p2 <- p2 +
            ggplot2::annotate("rect", xmin = neg_log_min, xmax = neg_log_max, ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "#2C7BE5") +
            ggplot2::geom_vline(xintercept = c(neg_log_min, neg_log_max), color = "#2C7BE5", linewidth = 0.8)
    }
    if (positive_gate_present) {
        p2 <- p2 +
            ggplot2::annotate("rect", xmin = log10(gate_min), xmax = log10(gate_max), ymin = -Inf, ymax = Inf, alpha = 0.18, fill = "#D62728") +
            ggplot2::geom_vline(xintercept = c(log10(gate_min), log10(gate_max)), color = "#D62728", linewidth = 1)
    }
    if (length(comp_means) > 0) {
        p2 <- p2 + ggplot2::geom_vline(xintercept = comp_means, color = "black", linetype = "dashed", alpha = 0.35, linewidth = 0.5)
    }
    if (isTRUE(use_scatter_gating)) {
        scatter_x <- log10(pmax(peak_vals, 1))
        scatter_class <- ifelse(
            peak_vals >= gate_min & peak_vals <= gate_max,
            "positive",
            ifelse(
                negative_gate_present &
                    peak_vals >= 10^neg_log_min &
                    peak_vals <= 10^neg_log_max,
                "negative",
                "other"
            )
        )
        scatter_dt <- data.table::data.table(
            x = scatter_x,
            y = gated_data[, fsc],
            class = factor(scatter_class, levels = c("negative", "other", "positive"))
        )
        p2_scatter <- ggplot2::ggplot(scatter_dt, ggplot2::aes(x, y, color = class)) +
            ggplot2::annotate("rect", xmin = neg_log_min, xmax = neg_log_max, ymin = -Inf, ymax = Inf, alpha = 0.12, fill = "#2C7BE5") +
            ggplot2::annotate("rect", xmin = log10(gate_min), xmax = log10(gate_max), ymin = -Inf, ymax = Inf, alpha = 0.12, fill = "#D62728") +
            ggplot2::geom_point(size = 0.45, alpha = 0.45) +
            ggplot2::scale_color_manual(values = c(negative = "#2C7BE5", other = "grey55", positive = "#D62728"), drop = FALSE) +
            ggplot2::labs(
                title = paste0(sn, " - ", peak_channel, " corrected scatter gate"),
                subtitle = paste(unlist(lapply(c(attr(vals_log, "gate_method"), paste0(round(100 * nrow(final_gated_data) / nrow(gated_data), 1), "% positive gated | blue = negative gate | red = bright gate")), strwrap, width = 95), use.names = FALSE), collapse = "\n"),
                x = paste0("log10(", peak_channel, ")"),
                y = fsc_desc,
                color = "class"
            ) +
            ggplot2::theme_minimal(base_size = 11) +
            ggplot2::theme(legend.position = "bottom", plot.subtitle = ggplot2::element_text(size = 8.4, lineheight = 1.05))
        ggplot2::ggsave(file.path(out_path, "intensity_scatter", paste0(sn, "_intensity_scatter.png")), p2_scatter, width = 6.5, height = 4, dpi = 300)
    } else {
        ggplot2::ggsave(file.path(out_path, "histogram", paste0(sn, "_histogram.png")), p2, width = 6.5, height = 4, dpi = 300)
    }

    log_mat <- log10(pmax(final_gated_data[, detector_names, drop = FALSE], 1e-3))
    min_y <- floor(min(log_mat, na.rm = TRUE))
    max_y <- ceiling(max(log_mat, na.rm = TRUE))
    breaks <- seq(min_y, max_y, length.out = 151)
    bin_mid <- (breaks[-1] + breaks[-length(breaks)]) / 2
    bin_height <- breaks[2] - breaks[1]
    counts_mat <- vapply(seq_len(ncol(log_mat)), function(j) as.numeric(graphics::hist(log_mat[, j], breaks = breaks, plot = FALSE)$counts), numeric(length(bin_mid)))
    rownames(counts_mat) <- as.character(seq_along(bin_mid))
    colnames(counts_mat) <- as.character(seq_len(ncol(log_mat)))
    dt_c <- data.table::as.data.table(as.table(counts_mat))
    data.table::setnames(dt_c, c("bin_idx", "ch_idx", "count"))
    dt_c$bin_idx <- as.integer(as.character(dt_c$bin_idx))
    dt_c$ch_idx <- as.integer(as.character(dt_c$ch_idx))
    dt_c$y_orig <- bin_mid[dt_c$bin_idx]
    dt_c$fill <- log10(dt_c$count + 1)
    min_bin_count <- 3
    dt_c <- dt_c[dt_c$count >= min_bin_count, ]
    y_power <- 1.5
    dt_c$y <- dt_c$y_orig^y_power

    if (nrow(dt_c) == 0) {
        message("  Note: dt_c is empty for ", sn, ". Max count: ", max(counts_mat, na.rm = TRUE))
    } else {
        message("  Plotting ", nrow(dt_c), " bins for ", sn)
    }

    vlines <- which(diff(det_info$laser_nm) != 0) + 0.5
    fill_lo <- min(dt_c$fill, na.rm = TRUE)
    fill_hi <- quantile(dt_c$fill, 0.96, na.rm = TRUE)
    y_breaks_orig <- 0:ceiling(max_y)
    y_breaks_trans <- y_breaks_orig^y_power
    y_labels <- vapply(y_breaks_orig, function(x) paste0("10^", x), character(1))

    p3 <- ggplot2::ggplot(dt_c, ggplot2::aes(ch_idx, y, fill = fill)) +
        ggplot2::geom_tile(width = 0.7, height = bin_height * 3) +
        ggplot2::scale_fill_gradientn(
            colors = c("#0000FF", "#00FFFF", "#00FF00", "#FFFF00", "#FF0000"),
            limits = c(fill_lo, fill_hi), oob = scales::squish
        ) +
        ggplot2::scale_x_continuous(breaks = seq_along(detector_names), labels = detector_labels) +
        ggplot2::scale_y_continuous(limits = c(0, (max_y + 0.5)^y_power), breaks = y_breaks_trans, labels = y_labels) +
        ggplot2::coord_cartesian(expand = FALSE) +
        ggplot2::labs(title = paste0(sn, " - Spectrum"), x = NULL, y = "Intensity") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6), legend.position = "none", panel.background = ggplot2::element_rect(fill = "white", color = NA))
    ggplot2::ggsave(file.path(out_path, "spectrum", paste0(sn, "_spectrum.png")), p3, width = 300, height = 120, units = "mm", dpi = 600)

    invisible(NULL)
}

# Processes a single FCS file through the entire reference construction pipeline.
# Reads the file, performs scatter gating, identifies the peak channel, runs histogram gating
# to separate positive and negative events, calculates the normalized spectrum, and saves QC plots.
# Returns a list containing the processed spectrum, QC summary row, and sample name metadata.
.process_reference_file <- function(fcs_file,
                                    control_df,
                                    sample_patterns,
                                    metadata,
                                    config,
                                    af_data_raw = NULL,
                                    universal_negatives = NULL,
                                    scc_background = NULL) {
    sn_ext <- basename(fcs_file)
    sn <- tools::file_path_sans_ext(sn_ext)

    is_extra_af <- .is_reference_extra_af_file(
        fcs_file = fcs_file,
        include_multi_af = config$include_multi_af,
        af_dir = config$af_dir
    )

    row_info <- .get_control_rows_for_reference(control_df, c(sn_ext, sn))
    sample_info <- .resolve_reference_sample_type(
        filename = sn,
        row_info = row_info,
        patterns = sample_patterns,
        default = config$default_sample_type
    )

    fluor_name <- if (nrow(row_info) > 0 && !is.na(row_info$fluorophore[1])) row_info$fluorophore[1] else sample_info$pattern
    marker_name <- if (nrow(row_info) > 0 && "marker" %in% colnames(row_info) && !is.na(row_info$marker[1])) {
        trimws(as.character(row_info$marker[1]))
    } else {
        ""
    }

    if (isTRUE(config$exclude_af) && (.is_af_control_row(fluorophore = fluor_name, marker = marker_name, filename = sn_ext) || is_extra_af)) {
        return(NULL)
    }
    if (is_extra_af) {
        return(NULL)
    }
    if (.is_af_control_row(fluorophore = fluor_name, marker = marker_name, filename = sn_ext)) {
        return(NULL)
    }

    message("Processing SCC: ", fluor_name, " (", sn, ")")
    ff <- flowCore::read.FCS(fcs_file, transformation = FALSE, truncate_max_range = FALSE)
    pd <- flowCore::pData(flowCore::parameters(ff))
    raw_data <- flowCore::exprs(ff)

    if (!.validate_reference_raw_data(raw_data, sn)) {
        return(NULL)
    }

    scatter_info <- .compute_reference_scatter_gate(
        raw_data = raw_data,
        pd = pd,
        sample_type = sample_info$type,
        outlier_percentile = config$outlier_percentile,
        debris_percentile = config$debris_percentile,
        subsample_n = config$subsample_n,
        max_clusters = config$max_clusters,
        min_cluster_proportion = config$min_cluster_proportion,
        gate_contour_beads = config$gate_contour_beads,
        gate_contour_cells = config$gate_contour_cells,
        bead_gate_scale = config$bead_gate_scale
    )
    if (is.null(scatter_info)) {
        return(NULL)
    }

    peak_info <- .select_reference_peak_channel(
        gated_data = scatter_info$gated_data,
        detector_names = metadata$detector_names,
        row_info = row_info,
        channel_alias_map = metadata$channel_alias_map,
        sn_ext = sn_ext,
        sn = sn,
        af_data_raw = af_data_raw,
        clean_with_af = isTRUE(config$clean_scc_with_unstained) && identical(sample_info$type, "cells")
    )
    peak_channel <- peak_info$peak_channel
    message("  Peak channel: ", peak_channel)

    peak_vals <- scatter_info$gated_data[, peak_channel]
    gate_fun <- if (isTRUE(config$use_scatter_gating)) {
        .compute_reference_scatter_intensity_gate
    } else {
        .compute_reference_histogram_gate
    }
    hist_info <- gate_fun(
        peak_vals = peak_vals,
        sample_type = sample_info$type,
        histogram_pct_beads = config$histogram_pct_beads,
        histogram_direction_beads = config$histogram_direction_beads,
        histogram_pct_cells = config$histogram_pct_cells,
        histogram_direction_cells = config$histogram_direction_cells,
        is_viability = nrow(row_info) > 0 &&
            "is.viability" %in% colnames(row_info) &&
            toupper(trimws(as.character(row_info$is.viability[1]))) == "TRUE"
    )

    final_gated_data <- scatter_info$gated_data[peak_vals >= hist_info$gate_min & peak_vals <= hist_info$gate_max, ]
    if (nrow(final_gated_data) < 10) {
        return(NULL)
    }

    # Stain Index Calculation
    mfi_pos <- stats::median(final_gated_data[, peak_channel], na.rm = TRUE)
    neg_idx <- which(peak_vals <= 10^quantile(hist_info$vals_log, 0.15, na.rm = TRUE))
    if (length(neg_idx) > 0) {
        neg_vals <- peak_vals[neg_idx]
        mfi_neg <- stats::median(neg_vals, na.rm = TRUE)
        sd_neg <- stats::mad(neg_vals, na.rm = TRUE)
        if (is.na(sd_neg) || sd_neg == 0) {
            sd_neg <- stats::sd(neg_vals, na.rm = TRUE)
        }
        if (is.na(sd_neg) || sd_neg == 0) {
            sd_neg <- 1e-6
        }
        stain_index <- (mfi_pos - mfi_neg) / (2 * sd_neg)
    } else {
        stain_index <- NA_real_
    }

    # Detector Saturation Check
    any_sat <- any(raw_data[, metadata$detector_names, drop = FALSE] >= 260000, na.rm = TRUE)

    spectrum_norm <- .compute_reference_spectrum(
        final_gated_data = final_gated_data,
        gated_data = scatter_info$gated_data,
        peak_vals = peak_vals,
        vals_log = hist_info$vals_log,
        detector_names = metadata$detector_names,
        row_info = row_info,
        sample_type = sample_info$type,
        af_data_raw = af_data_raw,
        universal_negatives = universal_negatives,
        scc_background = scc_background,
        scc_background_method = if (isTRUE(config$clean_scc_with_unstained)) config$scc_background_method else "none",
        scc_background_k = config$scc_background_k
    )

    if (isTRUE(config$save_qc_plots)) {
        .save_reference_qc_plots(
            sn = sn,
            raw_data = raw_data,
            gated_data = scatter_info$gated_data,
            final_gated_data = final_gated_data,
            pd = pd,
            fsc = scatter_info$fsc,
            ssc = scatter_info$ssc,
            fsc_max = scatter_info$fsc_max,
            ssc_max = scatter_info$ssc_max,
            final_gate = scatter_info$final_gate,
            vals_log = hist_info$vals_log,
            peak_vals = peak_vals,
            gate_min = hist_info$gate_min,
            gate_max = hist_info$gate_max,
            peak_channel = peak_channel,
            detector_names = metadata$detector_names,
            detector_labels = metadata$detector_labels,
            det_info = metadata$det_info,
            out_path = config$out_path,
            use_scatter_gating = config$use_scatter_gating
        )
    }

    list(
        sample_name = sn,
        result = data.table::data.table(
            sample = sn,
            fluorophore = fluor_name,
            type = sample_info$type,
            n_total = nrow(raw_data),
            n_final = nrow(final_gated_data),
            spectrum = list(spectrum_norm)
        ),
        qc_summary = data.table::data.table(
            sample = sn,
            fluorophore = fluor_name,
            marker = marker_name,
            type = sample_info$type,
            peak_channel = peak_channel,
            fsc_channel = scatter_info$fsc,
            ssc_channel = scatter_info$ssc,
            n_total = nrow(raw_data),
            n_scatter_gated = nrow(scatter_info$gated_data),
            n_final = nrow(final_gated_data),
            scatter_gate_pct = round(100 * nrow(scatter_info$gated_data) / max(nrow(raw_data), 1), 1),
            histogram_gate_pct = round(100 * nrow(final_gated_data) / max(nrow(scatter_info$gated_data), 1), 1),
            intensity_gate_type = {
                gate_type <- attr(hist_info$vals_log, "gate_type")
                if (!is.null(gate_type) && nzchar(gate_type)) gate_type else if (isTRUE(config$use_scatter_gating)) "scatter" else "histogram"
            },
            scc_background_method = {
                bg_info <- attr(spectrum_norm, "scc_background")
                if (!is.null(bg_info) && !is.null(bg_info$method)) bg_info$method else "none"
            },
            stain_index = round(stain_index, 1),
            saturated = ifelse(any_sat, "YES", "OK")
        )
    )
}

# Combines individual sample spectra and AF signatures into a single spillover matrix.
# Extracts spectra and variances, cleans up row/column names, structures the metadata attributes
# (such as variances, QC summary, and parameter info), and performs basic sanity checks.
# Returns the finalized reference matrix.
.finalize_reference_matrix <- function(results_list,
                                       qc_summary_list,
                                       af_signatures_norm,
                                       af_data_raw,
                                       af_bank_info,
                                       detector_names,
                                       pd_meta,
                                       save_qc_plots = FALSE,
                                       out_path = NULL) {
    results_dt <- if (length(results_list) > 0) data.table::rbindlist(results_list) else data.table::data.table()
    spectra_list <- if (nrow(results_dt) > 0) results_dt$spectrum else list()
    if (nrow(results_dt) > 0) names(spectra_list) <- results_dt$fluorophore
    scc_positive_events <- lapply(spectra_list, function(x) attr(x, "scc_positive_events"))
    has_scc_positive_events <- vapply(
        scc_positive_events,
        function(x) is.matrix(x) && nrow(x) > 0L && ncol(x) == length(detector_names),
        logical(1)
    )

    if (is.null(af_signatures_norm) && !is.null(af_data_raw)) {
        af_vec <- pmax(af_data_raw, 0)
        af_max <- max(af_vec, na.rm = TRUE)
        if (is.finite(af_max) && af_max > 0) {
            af_signatures_norm <- matrix(af_vec / af_max, nrow = 1)
            rownames(af_signatures_norm) <- "AF"
            colnames(af_signatures_norm) <- detector_names
        }
    }

    if (!is.null(af_signatures_norm) && nrow(af_signatures_norm) > 0) {
        for (i in seq_len(nrow(af_signatures_norm))) {
            nm <- rownames(af_signatures_norm)[i]
            if (is.na(nm) || !nzchar(nm)) nm <- if (i == 1) "AF" else paste0("AF_", i)
            base_nm <- nm
            k <- 2L
            while (nm %in% names(spectra_list)) {
                nm <- paste0(base_nm, "_", k)
                k <- k + 1L
            }
            spectra_list[[nm]] <- as.numeric(af_signatures_norm[i, ])
        }
    }

    if (length(spectra_list) == 0) {
        warning("No valid spectra found!")
        return(NULL)
    }

    M <- do.call(rbind, spectra_list)
    colnames(M) <- detector_names
    V <- do.call(rbind, lapply(spectra_list, function(x) {
        v <- attr(x, "variance")
        if (is.null(v)) rep(0, ncol(M)) else v
    }))
    rownames(V) <- rownames(M)
    colnames(V) <- colnames(M)
    attr(M, "variances") <- V

    if (length(qc_summary_list) > 0) {
        attr(M, "qc_summary") <- data.table::rbindlist(qc_summary_list)
    }
    if (isTRUE(save_qc_plots)) {
        attr(M, "qc_plot_dir") <- out_path
    }
    if (!is.null(af_bank_info)) {
        attr(M, "af_bank_info") <- af_bank_info
    }
    if (length(scc_positive_events) > 0L && any(has_scc_positive_events)) {
        attr(M, "scc_positive_events") <- scc_positive_events[has_scc_positive_events]
    }
    attr(M, "detector_pd") <- pd_meta
    M
}

#' Build a Reference Matrix from Single-Color Controls
#'
#' Reads SCC FCS files, performs FSC/SSC gating and positive-peak intensity gating,
#' then computes one normalized spectrum per fluorophore.
#'
#' This is the core matrix-construction step used before unmixing.
#'
#' @param input_folder Directory containing SCC `.fcs` files.
#' @param output_folder Directory where gating/spectrum plots are written.
#' @param save_qc_plots Logical; if `TRUE`, write FSC/SSC, intensity-gate, and spectrum plots.
#'   When `FALSE` (default), the function returns the matrix without writing QC files.
#' @param control_df Optional control mapping as a data.frame or CSV path.
#'   Expected columns: `filename`, `fluorophore`, `channel`; `universal.negative` is optional.
#' @param include_multi_af Logical; if `TRUE`, include additional AF controls from `af_dir`.
#' @param exclude_af Logical; if `TRUE`, ignore AF/unstained controls entirely,
#'   even when they are present in `control_df`, the SCC folder, or `af_dir`.
#' @param af_dir Directory with extra AF controls when `include_multi_af = TRUE`.
#' @param af_n_bands Number of SOM nodes used to extract AF basis signatures
#'   from pooled unstained/AF control events, or `"auto"` to use
#'   `af_auto_max_bands`. A mean AF row is prepended to the SOM nodes, so the
#'   default creates 100 SOM rows plus the mean row before QC/removal steps.
#' @param af_bands_per_file Deprecated compatibility argument. Multiple AF
#'   sources are pooled before SOM extraction; `af_n_bands`/`af_auto_max_bands`
#'   control the size of the one shared AF bank.
#' @param af_max_cells Maximum number of scatter-gated AF events used when
#'   deriving AF basis signatures.
#' @param af_auto_max_bands Maximum SOM nodes that `"auto"` may create before
#'   prepending the mean AF row. Default is 100.
#' @param af_min_cluster_events Compatibility argument retained for older
#'   workflows.
#' @param af_min_cluster_proportion Compatibility argument retained for older
#'   workflows.
#' @param af_n_bands_sensitivity Compatibility argument retained for older
#'   workflows.
#' @param af_refine Logical; if `TRUE`, run a second AF-library refinement pass
#'   that modulates SOM AF spectra from high-error unstained cells.
#' @param af_refine_problem_quantile Quantile used to choose high-error
#'   unstained cells for AF refinement. Default is 0.99.
#' @param af_deduplicate Logical; if `TRUE`, remove near-identical AF spectra
#'   after SOM extraction/refinement using cosine similarity. Default is
#'   `FALSE` for backward compatibility.
#' @param af_deduplication_threshold Cosine similarity threshold used when
#'   `af_deduplicate = TRUE`. Default is 0.99.
#' @param af_contaminant_threshold Cosine similarity threshold used to remove
#'   AF candidates that are too similar to known fluorophore spectra. Default is
#'   0.99.
#' @param af_assignment AF candidate assignment used during AF refinement.
#'   `"projection"` (default) chooses the AF band whose expected fluorophore
#'   spillover best brings apparent fluorophore signal toward zero.
#'   `"residual_alignment"` keeps the older detector-residual alignment score
#'   as an explicit alternative.
#' @param seed Optional integer seed for deterministic subsampling/clustering.
#' @param default_sample_type Fallback type when filename heuristics are ambiguous (`"beads"` or `"cells"`).
#' @param cytometer Cytometer name used as a channel-mapping hint. The default,
#'   `"auto"`, infers the cytometer from FCS detector names when possible.
#' @param use_scatter_gating Logical; if `TRUE` (default), use the intensity-vs-FSC
#'   scatter gate for final positive/negative population selection. If `FALSE`,
#'   use the legacy one-dimensional histogram gate.
#' @param clean_scc_with_unstained Logical; if `TRUE`, cell SCC spectra are
#'   cleaned with scatter-matched unstained/AF events before spectrum
#'   derivation. Bead controls keep their existing background subtraction.
#' @param scc_background_method Background method for cell SCC cleaning.
#'   `"scatter_knn"` matches each stained event to unstained events by FSC/SSC.
#'   Use `"none"` to disable matched background subtraction.
#' @param scc_background_k Number of nearest unstained events averaged for
#'   `"scatter_knn"` cell SCC background subtraction.
#' @param histogram_pct_beads Quantile width for the bead histogram gate.
#' @param histogram_direction_beads Histogram gate direction for beads: `"right"` starts at the median,
#'   `"both"` centers on the median, and `"left"` ends at the median.
#' @param histogram_pct_cells Quantile width for the cell histogram gate.
#' @param histogram_direction_cells Histogram gate direction for cells: `"right"` starts at the median,
#'   `"both"` centers on the median, and `"left"` ends at the median.
#' @param outlier_percentile Upper-tail FSC/SSC filtering percentile.
#' @param debris_percentile Lower FSC range fraction for cell/unstained controls.
#'   For example, `0.08` excludes events below 8\% of the per-file high FSC scale.
#' @param bead_gate_scale Ellipse scaling factor for bead FSC/SSC gate.
#' @param max_clusters Maximum number of GMM components tested.
#' @param min_cluster_proportion Minimum population proportion kept from GMM fit.
#' @param gate_contour_beads Contour probability for bead gating ellipse/hull.
#' @param gate_contour_cells Contour probability for cell gating ellipse/hull.
#' @param subsample_n Maximum number of events used for GMM fitting per file.
#'
#' @return Numeric matrix with rows = fluorophores and columns = detectors
#'   (normalized spectra). The matrix carries SCC-derived detector noise floors
#'   in `attr(M, "detector_noise")` for default WLS unmixing.
#' @export
#' @examples
#' if (interactive()) {
#'   M <- build_reference_matrix(
#'     input_folder = "scc",
#'     output_folder = "spectreasy_outputs/build_reference_plots",
#'     save_qc_plots = TRUE,
#'     control_df = "fcs_mapping.csv",
#'     cytometer = "auto"
#'   )
#'
#'   M_fast <- build_reference_matrix(
#'     input_folder = "scc",
#'     save_qc_plots = FALSE
#'   )
#' }
build_reference_matrix <- function(
  input_folder = "scc",
  output_folder = "gating_and_spectrum_plots",
  save_qc_plots = FALSE,
  control_df = NULL,
  include_multi_af = FALSE,
  exclude_af = FALSE,
  af_dir = "af",
  af_n_bands = "auto",
  af_bands_per_file = NULL,
  af_max_cells = 50000,
  af_auto_max_bands = 100,
  af_min_cluster_events = 20,
  af_min_cluster_proportion = 0.005,
  af_n_bands_sensitivity = 1.5,
  af_refine = FALSE,
  af_refine_problem_quantile = 0.99,
  af_deduplicate = FALSE,
  af_deduplication_threshold = 0.99,
  af_contaminant_threshold = 0.99,
  af_assignment = "projection",
  seed = NULL,
  default_sample_type = "beads",
  cytometer = "auto",
  use_scatter_gating = TRUE,
  clean_scc_with_unstained = TRUE,
  scc_background_method = c("scatter_knn", "none"),
  scc_background_k = 3L,
  histogram_pct_beads = 0.98,
  histogram_direction_beads = "right",
  histogram_pct_cells = 0.35,
  histogram_direction_cells = "right",
  outlier_percentile = 0.02,
  debris_percentile = 0.08,
  bead_gate_scale = 1.3,
  max_clusters = 10,
  min_cluster_proportion = 0.03,
  gate_contour_beads = 0.95,
  gate_contour_cells = 0.90,
  subsample_n = 5000
) {
    control_df <- .normalize_build_reference_control_df(control_df)
    af_args <- .validate_build_reference_af_args(
        af_n_bands = af_n_bands,
        af_max_cells = af_max_cells,
        af_bands_per_file = af_bands_per_file,
        af_auto_max_bands = af_auto_max_bands,
        af_min_cluster_events = af_min_cluster_events,
        af_min_cluster_proportion = af_min_cluster_proportion,
        af_n_bands_sensitivity = af_n_bands_sensitivity
    )
    af_n_bands <- af_args$af_n_bands
    af_bands_per_file <- af_args$af_bands_per_file
    af_max_cells <- af_args$af_max_cells
    af_auto_max_bands <- af_args$af_auto_max_bands
    af_min_cluster_events <- af_args$af_min_cluster_events
    af_min_cluster_proportion <- af_args$af_min_cluster_proportion
    af_n_bands_sensitivity <- af_args$af_n_bands_sensitivity
    af_refine <- isTRUE(af_refine)
    af_deduplicate <- isTRUE(af_deduplicate)
    af_deduplication_threshold <- as.numeric(af_deduplication_threshold[1])
    if (!is.finite(af_deduplication_threshold) || is.na(af_deduplication_threshold) ||
        af_deduplication_threshold <= 0 || af_deduplication_threshold > 1) {
        stop("af_deduplication_threshold must be a number > 0 and <= 1.")
    }
    af_contaminant_threshold <- as.numeric(af_contaminant_threshold[1])
    if (!is.finite(af_contaminant_threshold) || is.na(af_contaminant_threshold) ||
        af_contaminant_threshold <= 0 || af_contaminant_threshold > 1) {
        stop("af_contaminant_threshold must be a number > 0 and <= 1.")
    }
    af_assignment <- .normalize_af_assignment(af_assignment, choices = c("projection", "residual_alignment"))
    scc_background_args <- .validate_scc_background_args(
        clean_scc_with_unstained = clean_scc_with_unstained,
        scc_background_method = scc_background_method,
        scc_background_k = scc_background_k
    )
    clean_scc_with_unstained <- scc_background_args$enabled
    scc_background_method <- scc_background_args$method
    scc_background_k <- scc_background_args$k
    af_refine_problem_quantile <- as.numeric(af_refine_problem_quantile[1])
    if (!is.finite(af_refine_problem_quantile) || is.na(af_refine_problem_quantile) ||
        af_refine_problem_quantile <= 0 || af_refine_problem_quantile >= 1) {
        stop("af_refine_problem_quantile must be a number between 0 and 1.")
    }

    .with_optional_seed(seed)

    sample_patterns <- get_fluorophore_patterns()
    file_info <- .prepare_reference_file_set(
        input_folder = input_folder,
        include_multi_af = include_multi_af,
        af_dir = af_dir,
        exclude_af = exclude_af
    )
    out_path <- .prepare_reference_output_path(output_folder = output_folder, save_qc_plots = save_qc_plots)
    metadata <- .prepare_reference_detector_info(file_info$fcs_files[1])
    cytometer <- .resolve_cytometer_from_pd(cytometer, metadata$pd_meta)

    config <- list(
        include_multi_af = include_multi_af,
        af_dir = af_dir,
        exclude_af = exclude_af,
        default_sample_type = default_sample_type,
        outlier_percentile = outlier_percentile,
        debris_percentile = debris_percentile,
        subsample_n = subsample_n,
        max_clusters = max_clusters,
        min_cluster_proportion = min_cluster_proportion,
        gate_contour_beads = gate_contour_beads,
        gate_contour_cells = gate_contour_cells,
        bead_gate_scale = bead_gate_scale,
        use_scatter_gating = use_scatter_gating,
        clean_scc_with_unstained = clean_scc_with_unstained,
        scc_background_method = scc_background_method,
        scc_background_k = scc_background_k,
        histogram_pct_beads = histogram_pct_beads,
        histogram_direction_beads = histogram_direction_beads,
        histogram_pct_cells = histogram_pct_cells,
        histogram_direction_cells = histogram_direction_cells,
        save_qc_plots = save_qc_plots,
        out_path = out_path,
        cytometer = cytometer,
        af_auto_max_bands = af_auto_max_bands,
        af_deduplicate = af_deduplicate,
        af_deduplication_threshold = af_deduplication_threshold,
        af_contaminant_threshold = af_contaminant_threshold,
        af_assignment = af_assignment,
        af_refine = af_refine,
        af_refine_problem_quantile = af_refine_problem_quantile
    )

    message("Found ", length(metadata$detector_names), " spectral detectors. Sorting by laser...")
    .validate_reference_detector_consistency(
        fcs_files = file_info$fcs_files_all,
        detector_names = metadata$detector_names
    )

    af_profiles_raw <- .collect_reference_af_profiles(
        control_df = control_df,
        fcs_files = file_info$fcs_files,
        fcs_files_all = file_info$fcs_files_all,
        detector_names = metadata$detector_names,
        af_n_bands = af_n_bands,
        af_bands_per_file = af_bands_per_file,
        af_max_cells = af_max_cells,
        af_auto_max_bands = af_auto_max_bands,
        af_deduplicate = af_deduplicate,
        af_deduplication_threshold = af_deduplication_threshold,
        af_contaminant_threshold = af_contaminant_threshold,
        af_assignment = af_assignment,
        af_refine = FALSE,
        af_refine_problem_quantile = af_refine_problem_quantile,
        exclude_af = exclude_af,
        config = config,
        extract_signatures = FALSE
    )

    universal_negatives <- .collect_reference_universal_negatives(
        control_df = control_df,
        fcs_files = file_info$fcs_files_all,
        detector_names = metadata$detector_names,
        sample_patterns = sample_patterns,
        config = config
    )

    results_list <- list()
    qc_summary_list <- list()
    for (fcs_file in file_info$fcs_files_all) {
        processed <- .process_reference_file(
            fcs_file = fcs_file,
            control_df = control_df,
            sample_patterns = sample_patterns,
            metadata = metadata,
            config = config,
            af_data_raw = af_profiles_raw$af_data_raw,
            universal_negatives = universal_negatives,
            scc_background = if (isTRUE(clean_scc_with_unstained)) af_profiles_raw$scc_background else NULL
        )
        if (is.null(processed)) {
            next
        }
        results_list[[processed$sample_name]] <- processed$result
        qc_summary_list[[processed$sample_name]] <- processed$qc_summary
    }
    .validate_reference_complete_controls(
        control_df = control_df,
        fcs_files_all = file_info$fcs_files_all,
        processed_results = results_list,
        exclude_af = exclude_af
    )

    marker_spectra <- NULL
    if (length(results_list) > 0) {
        marker_dt <- data.table::rbindlist(results_list)
        marker_spectra <- do.call(rbind, marker_dt$spectrum)
        rownames(marker_spectra) <- marker_dt$fluorophore
        colnames(marker_spectra) <- metadata$detector_names
    }

    af_profiles <- .collect_reference_af_profiles(
        control_df = control_df,
        fcs_files = file_info$fcs_files,
        fcs_files_all = file_info$fcs_files_all,
        detector_names = metadata$detector_names,
        af_n_bands = af_n_bands,
        af_bands_per_file = af_bands_per_file,
        af_max_cells = af_max_cells,
        af_auto_max_bands = af_auto_max_bands,
        af_deduplicate = af_deduplicate,
        af_deduplication_threshold = af_deduplication_threshold,
        af_contaminant_threshold = af_contaminant_threshold,
        af_assignment = af_assignment,
        af_refine = af_refine,
        af_refine_problem_quantile = af_refine_problem_quantile,
        exclude_af = exclude_af,
        config = config,
        fluor_spectra = marker_spectra,
        extract_signatures = TRUE
    )

    M <- .finalize_reference_matrix(
        results_list = results_list,
        qc_summary_list = qc_summary_list,
        af_signatures_norm = af_profiles$af_signatures_norm,
        af_data_raw = af_profiles$af_data_raw,
        af_bank_info = af_profiles$af_bank_info,
        detector_names = metadata$detector_names,
        pd_meta = metadata$pd_meta,
        save_qc_plots = save_qc_plots,
        out_path = out_path
    )
    if (is.null(M)) {
        return(M)
    }
    .attach_estimated_wls_detector_noise(
        M = M,
        scc_dir = input_folder,
        fallback = .default_wls_background_noise(),
        warn = FALSE
    )
}
