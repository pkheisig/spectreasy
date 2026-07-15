# Normalizes the control data.frame mapping.
# Reads the CSV file if a path is provided, validates that the required "filename" column is present,
# fills in missing columns with defaults, and cleans up string values by trimming whitespace.
# Returns a standardized data.frame or NULL/empty data.frame.
.normalize_build_reference_control_df <- function(control_df) {
    if (is.character(control_df) && length(control_df) == 1 && !is.na(control_df)) {
        if (!file.exists(control_df)) .spectreasy_stop_missing_file(control_df, label = "control_df file")
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
        if (!("is.viability" %in% colnames(control_df))) control_df$is.viability <- ""

        control_df$filename <- trimws(as.character(control_df$filename))
        control_df$fluorophore <- trimws(as.character(control_df$fluorophore))
        control_df$channel <- trimws(as.character(control_df$channel))
        control_df$control.type <- tolower(trimws(as.character(control_df$control.type)))
        control_df$universal.negative <- trimws(as.character(control_df$universal.negative))
        control_df$is.viability <- trimws(as.character(control_df$is.viability))
        if ("marker" %in% colnames(control_df)) {
            control_df$marker <- trimws(as.character(control_df$marker))
        }
        for (col in c("filename", "fluorophore", "channel", "control.type", "universal.negative", "is.viability")) {
            control_df[[col]][is.na(control_df[[col]])] <- ""
        }
        if ("marker" %in% colnames(control_df)) {
            control_df$marker[is.na(control_df$marker)] <- ""
        }
        control_df <- .assign_reference_automatic_negatives(control_df)
        control_df <- .canonicalize_primary_af_labels(control_df)
    }

    control_df
}

.assign_reference_automatic_negatives <- function(control_df) {
    if (is.null(control_df) || !is.data.frame(control_df) || nrow(control_df) == 0) {
        return(control_df)
    }
    if (!("universal.negative" %in% colnames(control_df))) control_df$universal.negative <- ""
    if (!("is.viability" %in% colnames(control_df))) control_df$is.viability <- ""

    dead_rows <- .is_dead_af_control_row(
        fluorophore = if ("fluorophore" %in% colnames(control_df)) control_df$fluorophore else NULL,
        marker = if ("marker" %in% colnames(control_df)) control_df$marker else NULL,
        filename = control_df$filename,
        control_type = if ("control.type" %in% colnames(control_df)) control_df$control.type else NULL
    )
    if (any(dead_rows)) {
        control_df$fluorophore[dead_rows] <- "AF_dead"
        if ("marker" %in% colnames(control_df)) {
            control_df$marker[dead_rows] <- "Dead cell background"
        }
    }

    viability_rows <- toupper(trimws(as.character(control_df$is.viability))) %in% c("TRUE", "T", "1", "YES", "Y")
    control_type <- if ("control.type" %in% colnames(control_df)) {
        tolower(trimws(as.character(control_df$control.type)))
    } else {
        rep("", nrow(control_df))
    }
    primary_af_rows <- .is_primary_af_control_row(
        fluorophore = if ("fluorophore" %in% colnames(control_df)) control_df$fluorophore else NULL,
        marker = if ("marker" %in% colnames(control_df)) control_df$marker else NULL,
        filename = control_df$filename
    ) & !dead_rows & control_type != "beads"
    bead_af_rows <- control_type == "beads" & (
        grepl("^AF_BEADS$", trimws(as.character(control_df$fluorophore)), ignore.case = TRUE) |
        vapply(control_df$filename, .reference_is_bead_negative_file, logical(1))
    )
    af_files <- control_df$filename[primary_af_rows & nzchar(control_df$filename)]
    dead_files <- control_df$filename[dead_rows & nzchar(control_df$filename)]
    bead_files <- control_df$filename[bead_af_rows & nzchar(control_df$filename)]

    empty_negative <- !nzchar(trimws(as.character(control_df$universal.negative)))
    source_rows <- primary_af_rows | dead_rows | bead_af_rows
    ordinary_cells <- control_type != "beads" & !viability_rows & !source_rows & empty_negative
    viability_targets <- control_type != "beads" & viability_rows & !source_rows & empty_negative
    bead_targets <- control_type == "beads" & !source_rows & empty_negative
    if (length(af_files) > 0L) control_df$universal.negative[ordinary_cells] <- "AF"
    if (length(dead_files) > 0L) control_df$universal.negative[viability_targets] <- dead_files[1]
    if (length(bead_files) > 0L) control_df$universal.negative[bead_targets] <- bead_files[1]
    control_df
}

# Validates the Autofluorescence (AF) modeling parameters.
# Ensures that the number of AF bands, maximum events per file, and bands per file parameters
# are positive integers, throwing an error if they are invalid.
# Returns a validated list containing these parameters.
.validate_build_reference_af_args <- function(af_n_bands,
                                              af_max_cells,
                                              af_min_cluster_events = 20,
                                              af_min_cluster_proportion = 0.005) {
    af_n_bands <- .normalize_positive_integer(af_n_bands, "af_n_bands")

    af_max_cells <- .normalize_positive_integer(af_max_cells, "af_max_cells")
    if (af_max_cells < 100L) {
        stop("af_max_cells must be an integer >= 100.")
    }

    af_min_cluster_events <- .normalize_positive_integer(af_min_cluster_events, "af_min_cluster_events")

    af_min_cluster_proportion <- .normalize_numeric_scalar(
        af_min_cluster_proportion, "af_min_cluster_proportion", lower = 0, upper = 1
    )

    list(
        af_n_bands = af_n_bands,
        af_max_cells = af_max_cells,
        af_min_cluster_events = af_min_cluster_events,
        af_min_cluster_proportion = af_min_cluster_proportion
    )
}

.validate_reference_refine_arg <- function(refine) {
    if (!is.logical(refine) || length(refine) != 1L || is.na(refine)) {
        stop("refine must be TRUE or FALSE.", call. = FALSE)
    }
    isTRUE(refine)
}

.reference_min_af_cluster_size <- function(n_events,
                                           min_cluster_events = 20,
                                           min_cluster_proportion = 0.005) {
    max(
        as.integer(min_cluster_events),
        as.integer(ceiling(as.numeric(min_cluster_proportion) * n_events))
    )
}

.reference_normalize_af_centers <- function(centers) {
    centers <- as.matrix(centers)
    centers[!is.finite(centers)] <- 0
    centers <- pmax(centers, 0)
    center_scale <- apply(centers, 1, max, na.rm = TRUE)
    ok <- is.finite(center_scale) & center_scale > 0
    if (any(ok)) {
        centers[ok, ] <- sweep(centers[ok, , drop = FALSE], 1, center_scale[ok], "/")
    }
    if (any(!ok)) {
        centers[!ok, ] <- 0
    }
    centers
}

.reference_af_row_mask <- function(M) {
    grepl("^AF($|_)", rownames(M), ignore.case = TRUE)
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

.reference_af_candidate_qc <- function(af_spectra,
                                       marker_spectra,
                                       contaminant_threshold = 0.99) {
    af_spectra <- as.matrix(af_spectra)
    marker_spectra <- as.matrix(marker_spectra)
    if (nrow(af_spectra) == 0L || nrow(marker_spectra) == 0L) {
        return(af_spectra)
    }
    sims <- .reference_cosine_matrix(af_spectra, marker_spectra)
    max_sim <- apply(sims, 1, max, na.rm = TRUE)
    keep <- !is.finite(max_sim) | max_sim < contaminant_threshold
    af_spectra[keep, , drop = FALSE]
}

.reference_filter_af_contaminant_events <- function(af_events,
                                                   marker_spectra,
                                                   contaminant_threshold = 0.99) {
    af_events <- as.matrix(af_events)
    marker_spectra <- as.matrix(marker_spectra)
    if (nrow(af_events) == 0L || nrow(marker_spectra) == 0L) {
        return(af_events)
    }
    centered <- sweep(af_events, 2, colMeans(af_events, na.rm = TRUE), "-")
    sims <- .reference_cosine_matrix(centered, marker_spectra)
    max_sim <- apply(sims, 1, max, na.rm = TRUE)
    keep <- !is.finite(max_sim) | max_sim < contaminant_threshold
    af_events[keep, , drop = FALSE]
}

.reference_cluster_error_ratios <- function(spill_ratios,
                                            n_centers,
                                            seed = NULL) {
    spill_ratios <- as.matrix(spill_ratios)
    spill_ratios[!is.finite(spill_ratios)] <- 0
    unique_count <- nrow(unique(as.data.frame(round(spill_ratios, digits = 8))))
    n_eff <- min(as.integer(n_centers), nrow(spill_ratios), unique_count)
    if (!is.finite(n_eff) || is.na(n_eff) || n_eff < 1L) {
        n_eff <- 1L
    }
    if (n_eff == 1L) {
        cluster <- rep(1L, nrow(spill_ratios))
        centers <- matrix(apply(spill_ratios, 2, stats::median, na.rm = TRUE), nrow = 1L)
        colnames(centers) <- colnames(spill_ratios)
        return(list(cluster = cluster, centers = centers))
    }
    if (!is.null(seed)) set.seed(seed)
    km <- withCallingHandlers(
        stats::kmeans(spill_ratios, centers = n_eff, nstart = 8, iter.max = 100),
        warning = function(w) {
            if (grepl("Quick-TRANSfer stage steps exceeded maximum", conditionMessage(w), fixed = TRUE)) {
                invokeRestart("muffleWarning")
            }
        }
    )
    list(cluster = km$cluster, centers = as.matrix(km$centers))
}

.reference_refine_af_bank <- function(M,
                                      af_events,
                                      af_n_bands,
                                      af_max_cells,
                                      af_min_cluster_events,
                                      af_min_cluster_proportion,
                                      n_threads = 1L,
                                      seed = NULL,
                                      problem_quantile = 0.99,
                                      contaminant_threshold = 0.99,
                                      verbose = TRUE) {
    M <- .as_reference_matrix(M, "M")
    af_rows <- .reference_af_row_mask(M)
    marker_rows <- !af_rows
    if (!any(af_rows) || !any(marker_rows) || is.null(af_events)) {
        return(NULL)
    }

    detector_names <- colnames(M)
    af_events <- as.matrix(af_events)
    if (!all(detector_names %in% colnames(af_events))) {
        return(NULL)
    }
    af_events <- af_events[, detector_names, drop = FALSE]
    af_events <- af_events[stats::complete.cases(af_events), , drop = FALSE]
    if (nrow(af_events) < 100L) {
        return(NULL)
    }
    if (nrow(af_events) > af_max_cells) {
        af_events <- af_events[sample.int(nrow(af_events), af_max_cells), , drop = FALSE]
    }

    marker_spectra <- M[marker_rows, , drop = FALSE]
    base_af <- M[af_rows, , drop = FALSE]
    af_events <- .reference_filter_af_contaminant_events(
        af_events = af_events,
        marker_spectra = marker_spectra,
        contaminant_threshold = contaminant_threshold
    )
    if (nrow(af_events) < 100L) {
        return(NULL)
    }

    if (isTRUE(verbose)) {
        .spectreasy_console_step("AF refine", "AutoSpectral-style unstained residual modulation")
    }
    assignments <- .autospectral_assign_af_fluorophores(
        Y = af_events,
        marker_spectra = marker_spectra,
        af_spectra = base_af,
        n_threads = n_threads
    )
    assignments[!is.finite(assignments) | assignments < 1L | assignments > nrow(base_af)] <- 1L

    marker_idx <- seq_len(nrow(marker_spectra)) + 1L
    unmixed_markers <- .unmix_with_fixed_reference(
        Y = af_events,
        M = marker_spectra,
        method = "OLS",
        wls_noise = list(noise_floor = numeric(0), signal_scale = numeric(0), max_weight_ratio = .default_wls_max_weight_ratio()),
        rwls_max_iter = 1L
    )
    unmixed <- cbind(AF = rep(0, nrow(af_events)), unmixed_markers)
    residuals <- matrix(0, nrow = nrow(af_events), ncol = ncol(M), dimnames = list(NULL, detector_names))
    projected_markers <- matrix(0, nrow = nrow(af_events), ncol = ncol(M), dimnames = list(NULL, detector_names))

    for (af_i in unique(assignments)) {
        event_idx <- which(assignments == af_i)
        X <- rbind(AF = base_af[af_i, ], marker_spectra)
        coeff <- .unmix_with_fixed_reference(
            Y = af_events[event_idx, , drop = FALSE],
            M = X,
            method = "OLS",
            wls_noise = list(noise_floor = numeric(0), signal_scale = numeric(0), max_weight_ratio = .default_wls_max_weight_ratio()),
            rwls_max_iter = 1L
        )
        unmixed[event_idx, ] <- coeff
        residuals[event_idx, ] <- af_events[event_idx, , drop = FALSE] - (coeff %*% X)
        projected_markers[event_idx, ] <- coeff[, marker_idx, drop = FALSE] %*% marker_spectra
    }

    error <- residuals + projected_markers
    marker_abundance <- unmixed[, marker_idx, drop = FALSE]
    if (ncol(marker_abundance) > 1L) {
        error_magnitude <- sqrt(rowSums(marker_abundance^2))
    } else {
        error_magnitude <- abs(marker_abundance[, 1])
    }
    error_magnitude[!is.finite(error_magnitude)] <- 0

    q <- as.numeric(problem_quantile[1])
    if (!is.finite(q) || is.na(q) || q <= 0 || q >= 1) q <- 0.99
    repeat {
        threshold <- stats::quantile(error_magnitude, q, na.rm = TRUE, names = FALSE)
        problem_idx <- which(error_magnitude > threshold)
        if (length(problem_idx) >= 500L || q < 0.5) break
        q <- q - 0.05
    }
    if (length(problem_idx) <= 10L) {
        return(NULL)
    }

    af_abundance <- unmixed[problem_idx, 1]
    af_abundance[!is.finite(af_abundance) | af_abundance == 0] <- 1e-6
    spill_ratios <- sweep(error[problem_idx, , drop = FALSE], 1, af_abundance, "/")
    clusters <- .reference_cluster_error_ratios(
        spill_ratios = spill_ratios,
        n_centers = af_n_bands,
        seed = seed
    )

    modulated <- list()
    for (cl in sort(unique(clusters$cluster))) {
        local_idx <- which(clusters$cluster == cl)
        global_idx <- problem_idx[local_idx]
        median_ratio <- apply(spill_ratios[local_idx, , drop = FALSE], 2, stats::median, na.rm = TRUE)
        median_ratio[!is.finite(median_ratio)] <- 0
        contributing_af <- unique(assignments[global_idx])
        for (af_i in contributing_af) {
            updated <- base_af[af_i, ] * (1 + median_ratio)
            updated <- pmax(updated, 0)
            peak <- max(updated, na.rm = TRUE)
            if (is.finite(peak) && peak > 1e-12) {
                updated <- updated / peak
                modulated[[length(modulated) + 1L]] <- updated
            }
        }
    }
    if (!length(modulated)) {
        return(NULL)
    }

    candidates <- rbind(base_af, do.call(rbind, modulated))
    colnames(candidates) <- detector_names
    candidates <- .reference_normalize_af_centers(candidates)
    candidates <- .reference_af_candidate_qc(
        af_spectra = candidates,
        marker_spectra = marker_spectra,
        contaminant_threshold = contaminant_threshold
    )
    if (nrow(candidates) == 0L) {
        candidates <- base_af
    }

    km <- .reference_kmeans_af_centers(
        af_shape = candidates,
        n_centers = af_n_bands,
        min_cluster_events = 1L,
        min_cluster_proportion = 0
    )
    refined <- .reference_normalize_af_centers(km$centers)
    rownames(refined) <- c("AF", if (nrow(refined) > 1L) paste0("AF_", seq.int(2L, nrow(refined))) else NULL)
    colnames(refined) <- detector_names

    list(
        signatures = refined,
        selection = list(
            method = if (identical(km$center_method, "median")) "median_refined_fixed" else "kmeans_refined_fixed",
            n_bands = nrow(refined),
            requested_bands = as.integer(af_n_bands),
            raw_center_count = nrow(km$centers),
            final_bands = nrow(refined),
            cluster_sizes = km$cluster_sizes,
            min_cluster_size = km$min_cluster_size,
            base_bands = nrow(base_af),
            candidate_bands = nrow(candidates),
            modulated_candidates = length(modulated),
            problem_cells = length(problem_idx),
            problem_quantile = q,
            problem_threshold = threshold
        )
    )
}

.replace_reference_af_rows <- function(M, af_signatures_norm, af_bank_info = NULL) {
    attrs <- attributes(M)
    af_rows <- .reference_af_row_mask(M)
    marker_M <- M[!af_rows, , drop = FALSE]
    out <- rbind(marker_M, af_signatures_norm)
    colnames(out) <- colnames(M)
    for (nm in setdiff(names(attrs), c("dim", "dimnames", "names"))) {
        attr(out, nm) <- attrs[[nm]]
    }
    if (!is.null(af_bank_info)) {
        attr(out, "af_bank_info") <- af_bank_info
    }
    out
}

.reference_kmeans_af_centers <- function(af_shape,
                                         n_centers,
                                         min_cluster_events = 20,
                                         min_cluster_proportion = 0.005,
                                         nstart = 10,
                                         iter.max = 100) {
    af_shape <- as.matrix(af_shape)
    if (nrow(af_shape) == 0) {
        return(list(centers = af_shape, cluster_sizes = integer(), requested_centers = as.integer(n_centers)))
    }
    unique_count <- nrow(unique(as.data.frame(round(af_shape, digits = 8))))
    n_eff <- min(as.integer(n_centers), nrow(af_shape), unique_count)
    if (!is.finite(n_eff) || is.na(n_eff) || n_eff < 1L) {
        n_eff <- 1L
    }

    if (n_eff == 1L) {
        centers <- matrix(apply(af_shape, 2, stats::median, na.rm = TRUE), nrow = 1L)
        colnames(centers) <- colnames(af_shape)
        cluster_sizes <- nrow(af_shape)
    } else {
        km <- withCallingHandlers(
            stats::kmeans(af_shape, centers = n_eff, nstart = nstart, iter.max = iter.max),
            warning = function(w) {
                if (grepl("Quick-TRANSfer stage steps exceeded maximum", conditionMessage(w), fixed = TRUE)) {
                    invokeRestart("muffleWarning")
                }
            }
        )
        centers <- as.matrix(km$centers)
        cluster_sizes <- tabulate(km$cluster, nbins = nrow(centers))
    }

    min_cluster_size <- .reference_min_af_cluster_size(
        n_events = nrow(af_shape),
        min_cluster_events = min_cluster_events,
        min_cluster_proportion = min_cluster_proportion
    )
    ord <- order(cluster_sizes, decreasing = TRUE)
    centers <- centers[ord, , drop = FALSE]
    cluster_sizes <- cluster_sizes[ord]

    # Small k-means components are usually rare event-level artefacts rather
    # than reusable AF phenotypes. The minimum was previously reported in the
    # selection metadata but never applied, allowing tiny, nonphysiological
    # spectra to enter the AF bank.
    keep <- cluster_sizes >= min_cluster_size
    if (!any(keep)) keep[1L] <- TRUE
    centers <- centers[keep, , drop = FALSE]
    cluster_sizes <- cluster_sizes[keep]

    list(
        centers = centers,
        cluster_sizes = cluster_sizes,
        min_cluster_size = min_cluster_size,
        requested_centers = as.integer(n_centers),
        effective_centers = as.integer(n_eff),
        center_method = if (n_eff == 1L) "median" else "kmeans"
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

# Prepares the complete list of FCS files to be processed.
# Locates FCS files in the input folder.
# Returns a list containing 'fcs_files' and 'fcs_files_all'.
.prepare_reference_file_set <- function(input_folder, control_df = NULL) {
    if (!dir.exists(input_folder)) {
        .spectreasy_stop_missing_directory(input_folder, label = "input_folder")
    }
    fcs_files_all <- list.files(input_folder, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
    if (length(fcs_files_all) == 0) {
        .spectreasy_stop_empty_fcs_directory(input_folder, label = "input_folder")
    }
    fcs_files_all <- .control_validation_select_scc_files(control_df, fcs_files_all)
    if (length(fcs_files_all) == 0) {
        stop("None of the FCS files listed in control_df were found in input_folder.", call. = FALSE)
    }

    list(fcs_files = fcs_files_all, fcs_files_all = fcs_files_all)
}

# Prepares and returns the normalized path of the output directory.
# If save_qc_plots is TRUE, recursively creates subdirectories for FSC/SSC plots,
# intensity gates, and spectra.
# Returns the absolute path string.
.prepare_reference_output_path <- function(output_folder, save_qc_plots = FALSE) {
    out_path <- normalizePath(output_folder, mustWork = FALSE)
    if (isTRUE(save_qc_plots)) {
        dir.create(file.path(out_path, "fsc_ssc"), showWarnings = FALSE, recursive = TRUE)
        dir.create(file.path(out_path, "singlet"), showWarnings = FALSE, recursive = TRUE)
        dir.create(file.path(out_path, "histogram"), showWarnings = FALSE, recursive = TRUE)
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
    ff_meta <- .spectreasy_read_fcs(first_fcs_file, label = "SCC FCS file")
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

.read_reference_manual_gates <- function(manual_gate_file = NULL) {
    if (is.null(manual_gate_file) || length(manual_gate_file) == 0 || is.na(manual_gate_file[1])) return(NULL)
    path <- normalizePath(as.character(manual_gate_file)[1], mustWork = FALSE)
    if (!file.exists(path)) return(NULL)
    df <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    required <- c("gate_type", "scope", "filename", "x_channel", "y_channel", "plot_mode", "vertex_index", "x", "y")
    if (is.null(df) || nrow(df) == 0) {
        stop("Could not read a valid manual gate CSV: ", path, call. = FALSE)
    }
    missing_cols <- setdiff(required, colnames(df))
    if (length(missing_cols) > 0) {
        stop("Manual gate CSV is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    settings <- list()
    setting_rows <- df[df$gate_type == "setting", , drop = FALSE]
    if (nrow(setting_rows) > 0) {
        for (i in seq_len(nrow(setting_rows))) {
            key <- trimws(as.character(setting_rows$x_channel[i]))
            value <- trimws(as.character(setting_rows$x[i]))
            if (!nzchar(key)) next
            if (identical(key, "point_size")) {
                settings$point_size <- suppressWarnings(as.numeric(value))
            } else if (identical(key, "max_points")) {
                settings$max_points <- suppressWarnings(as.integer(value))
            } else if (identical(key, "histogram_bins")) {
                settings$histogram_bins <- suppressWarnings(as.integer(value))
            } else if (identical(key, "histogram_transform")) {
                settings$histogram_transform <- tolower(value)
            }
        }
    }
    keep <- df$gate_type != "setting" &
        (df$plot_mode == "blocked" | (df$plot_mode != "missing" & suppressWarnings(as.integer(df$vertex_index)) > 0))
    df <- df[keep, , drop = FALSE]
    if (nrow(df) == 0) {
        out <- list()
    } else {
        out <- split(df, paste(df$gate_type, df$scope, df$filename, sep = "\r"))
    }
    attr(out, "settings") <- settings
    out
}

.reference_manual_gate_settings <- function(manual_gates) {
    settings <- if (!is.null(manual_gates)) attr(manual_gates, "settings") else NULL
    if (is.null(settings)) settings <- list()
    bins_value <- if (is.null(settings$histogram_bins)) 100L else settings$histogram_bins
    bins <- suppressWarnings(as.integer(bins_value))
    if (!is.finite(bins) || is.na(bins)) bins <- 100L
    bins <- min(max(bins, 5L), 500L)
    transform_value <- if (is.null(settings$histogram_transform)) "asinh" else settings$histogram_transform
    transform <- tolower(trimws(as.character(transform_value)))
    if (!transform %in% c("asinh", "linear", "log10", "biexponential")) transform <- "asinh"
    max_points_value <- if (is.null(settings$max_points)) 50000L else settings$max_points
    max_points <- suppressWarnings(as.integer(max_points_value))
    if (!is.finite(max_points) || is.na(max_points) || max_points <= 0L) max_points <- 50000L
    point_size_value <- if (is.null(settings$point_size)) 1.5 else settings$point_size
    point_size <- suppressWarnings(as.numeric(point_size_value))
    if (!is.finite(point_size) || is.na(point_size) || point_size <= 0) point_size <- 1.5
    list(
        histogram_bins = bins,
        histogram_transform = transform,
        max_points = max_points,
        point_size = point_size
    )
}

.reference_manual_gate_key <- function(gate_type, scope, filename = "") {
    paste(gate_type, scope, filename, sep = "\r")
}

.reference_manual_gate <- function(manual_gates, gate_type, filename, sample_type) {
    if (is.null(manual_gates)) return(NULL)
    keys <- c(
        .reference_manual_gate_key(gate_type, "file", filename),
        .reference_manual_gate_key(gate_type, sample_type, ""),
        .reference_manual_gate_key(gate_type, "cells", "")
    )
    hit <- keys[keys %in% names(manual_gates)][1]
    if (is.na(hit) || !nzchar(hit)) return(NULL)
    gate <- manual_gates[[hit]]
    gate[order(suppressWarnings(as.integer(gate$vertex_index))), , drop = FALSE]
}

.reference_gate_vertices <- function(gate) {
    if (is.null(gate) || nrow(gate) == 0) return(NULL)
    out <- data.frame(
        x = suppressWarnings(as.numeric(gate$x)),
        y = suppressWarnings(as.numeric(gate$y))
    )
    out <- out[stats::complete.cases(out), , drop = FALSE]
    if (nrow(out) == 0) return(NULL)
    out
}

.reference_gate_channel <- function(gate, column, fallback, raw_data) {
    channel <- if (!is.null(gate) && column %in% colnames(gate)) as.character(gate[[column]][1]) else fallback
    if (is.na(channel) || !nzchar(channel) || !channel %in% colnames(raw_data)) fallback else channel
}

.apply_reference_manual_scatter_gates <- function(raw_data, sample_type, filename, manual_gates) {
    scatter <- .resolve_reference_scatter_channels(raw_data)
    if (is.null(scatter)) return(NULL)
    fsc <- scatter$fsc
    ssc <- scatter$ssc
    fsc_h <- sub("-A$", "-H", fsc, ignore.case = TRUE)
    if (!fsc_h %in% colnames(raw_data)) {
        fsc_h <- grep("^FSC.*-H$", colnames(raw_data), value = TRUE, ignore.case = TRUE)[1]
    }

    keep <- rep(TRUE, nrow(raw_data))
    cell_keep <- rep(TRUE, nrow(raw_data))
    singlet_keep <- rep(TRUE, nrow(raw_data))
    final_gate <- NULL
    cell_gate <- .reference_manual_gate(manual_gates, "cell", filename, sample_type)
    cell_vertices <- .reference_gate_vertices(cell_gate)
    cell_x <- fsc
    cell_y <- ssc
    if (!is.null(cell_vertices) && nrow(cell_vertices) >= 3) {
        cell_x <- .reference_gate_channel(cell_gate, "x_channel", fsc, raw_data)
        cell_y <- .reference_gate_channel(cell_gate, "y_channel", ssc, raw_data)
        cell_keep <- sp::point.in.polygon(raw_data[, cell_x], raw_data[, cell_y], cell_vertices$x, cell_vertices$y) > 0
        keep <- keep & cell_keep
        if (identical(cell_x, fsc) && identical(cell_y, ssc)) {
            final_gate <- cell_vertices
        }
    }
    singlet_gate <- .reference_manual_gate(manual_gates, "singlet", filename, sample_type)
    singlet_vertices <- .reference_gate_vertices(singlet_gate)
    manual_attempted <- !is.null(cell_gate) || !is.null(singlet_gate)
    singlet_x <- fsc_h
    singlet_y <- fsc
    if (!is.null(singlet_vertices) && nrow(singlet_vertices) >= 3 && !is.na(fsc_h) && fsc_h %in% colnames(raw_data)) {
        singlet_x <- .reference_gate_channel(singlet_gate, "x_channel", fsc_h, raw_data)
        singlet_y <- .reference_gate_channel(singlet_gate, "y_channel", fsc, raw_data)
        singlet_keep <- sp::point.in.polygon(raw_data[, singlet_x], raw_data[, singlet_y], singlet_vertices$x, singlet_vertices$y) > 0
        keep <- keep & singlet_keep
    }
    if (sum(keep, na.rm = TRUE) < 10) return(NULL)
    if (all(keep) && !manual_attempted) return(NULL)
    gated_data <- raw_data[keep, , drop = FALSE]
    if (is.null(final_gate)) {
        x_min <- min(gated_data[, fsc], na.rm = TRUE)
        x_max <- max(gated_data[, fsc], na.rm = TRUE)
        y_min <- min(gated_data[, ssc], na.rm = TRUE)
        y_max <- max(gated_data[, ssc], na.rm = TRUE)
        final_gate <- data.frame(
            x = c(x_min, x_max, x_max, x_min, x_min),
            y = c(y_min, y_min, y_max, y_max, y_min)
        )
    }
    list(
        gated_data = gated_data,
        final_gate = final_gate,
        fsc = fsc,
        ssc = ssc,
        fsc_max = stats::quantile(raw_data[, fsc], 0.98, na.rm = TRUE),
        ssc_max = stats::quantile(raw_data[, ssc], 0.98, na.rm = TRUE),
        manual = TRUE,
        manual_gate_info = list(
            cell = list(vertices = cell_vertices, x_channel = cell_x, y_channel = cell_y, keep_count = sum(cell_keep, na.rm = TRUE)),
            singlet = list(
                vertices = singlet_vertices,
                x_channel = singlet_x,
                y_channel = singlet_y,
                source_count = sum(cell_keep, na.rm = TRUE),
                keep_count = sum(cell_keep & singlet_keep, na.rm = TRUE)
            ),
            settings = .reference_manual_gate_settings(manual_gates)
        )
    )
}

.apply_reference_manual_positive_gate <- function(gated_data, peak_channel, filename, sample_type, manual_gates) {
    gate <- .reference_manual_gate(manual_gates, "positive", filename, sample_type)
    verts <- .reference_gate_vertices(gate)
    if (is.null(gate) || is.null(verts) || nrow(verts) == 0) return(NULL)
    mode <- as.character(gate$plot_mode[1])
    peak_vals <- gated_data[, peak_channel]
    keep <- NULL
    if (identical(mode, "separator")) {
        keep <- peak_vals >= verts$x[1]
    } else if (identical(mode, "positive_1d") && nrow(verts) >= 2) {
        lim <- range(verts$x, na.rm = TRUE)
        keep <- peak_vals >= lim[1] & peak_vals <= lim[2]
    } else if (nrow(verts) >= 3) {
        x_channel <- as.character(gate$x_channel[1])
        y_channel <- as.character(gate$y_channel[1])
        if (!x_channel %in% colnames(gated_data)) x_channel <- peak_channel
        if (!y_channel %in% colnames(gated_data)) {
            scatter <- .resolve_reference_scatter_channels(gated_data)
            y_channel <- if (!is.null(scatter)) scatter$fsc else peak_channel
        }
        keep <- sp::point.in.polygon(gated_data[, x_channel], gated_data[, y_channel], verts$x, verts$y) > 0
    }
    if (is.null(keep) || sum(keep, na.rm = TRUE) < 10) return(NULL)
    vals_log <- log10(pmax(peak_vals, 1))
    attr(vals_log, "gate_type") <- paste0("manual_", mode)
    attr(vals_log, "gate_method") <- paste0("manual_", mode)
    attr(vals_log, "positive_gate_present") <- TRUE
    if (identical(mode, "positive_1d") && nrow(verts) >= 2) {
        positive_lim <- range(verts$x, na.rm = TRUE)
        if (all(is.finite(positive_lim)) && positive_lim[2] > positive_lim[1]) {
            attr(vals_log, "pos_raw_min") <- positive_lim[1]
            attr(vals_log, "pos_raw_max") <- positive_lim[2]
        }
    }
    negative_gate <- .reference_manual_gate(manual_gates, "negative", filename, sample_type)
    negative_verts <- .reference_gate_vertices(negative_gate)
    if (!is.null(negative_gate) && !is.null(negative_verts) && nrow(negative_verts) >= 2) {
        negative_mode <- as.character(negative_gate$plot_mode[1])
        if (identical(negative_mode, "negative_1d")) {
            negative_lim <- range(negative_verts$x, na.rm = TRUE)
            if (all(is.finite(negative_lim)) && negative_lim[2] > negative_lim[1]) {
                attr(vals_log, "negative_gate_present") <- TRUE
                attr(vals_log, "neg_log_min") <- log10(max(negative_lim[1], 1))
                attr(vals_log, "neg_log_max") <- log10(max(negative_lim[2], 1))
                attr(vals_log, "neg_raw_min") <- negative_lim[1]
                attr(vals_log, "neg_raw_max") <- negative_lim[2]
            }
        }
    }
    list(
        final_gated_data = gated_data[keep, , drop = FALSE],
        hist_info = list(
            vals_log = vals_log,
            gate_min = min(peak_vals[keep], na.rm = TRUE),
            gate_max = max(peak_vals[keep], na.rm = TRUE),
            positive_raw_min = attr(vals_log, "pos_raw_min"),
            positive_raw_max = attr(vals_log, "pos_raw_max"),
            negative_raw_min = attr(vals_log, "neg_raw_min"),
            negative_raw_max = attr(vals_log, "neg_raw_max")
        )
    )
}

.resolve_reference_positive_histogram_gate <- function(gated_data,
                                                       peak_channel,
                                                       filename,
                                                       sample_type,
                                                       manual_gates,
                                                       config,
                                                       row_info) {
    manual <- .apply_reference_manual_positive_gate(
        gated_data = gated_data,
        peak_channel = peak_channel,
        filename = filename,
        sample_type = sample_type,
        manual_gates = manual_gates
    )
    if (!is.null(manual)) return(manual)

    peak_vals <- gated_data[, peak_channel]
    hist_info <- .compute_reference_histogram_gate(
        peak_vals = peak_vals,
        sample_type = sample_type,
        histogram_pct_beads = config$histogram_pct_beads,
        histogram_direction_beads = config$histogram_direction_beads,
        histogram_pct_cells = config$histogram_pct_cells,
        histogram_direction_cells = config$histogram_direction_cells,
        is_viability = nrow(row_info) > 0 &&
            "is.viability" %in% colnames(row_info) &&
            toupper(trimws(as.character(row_info$is.viability[1]))) == "TRUE"
    )
    keep <- is.finite(peak_vals) &
        peak_vals >= hist_info$gate_min & peak_vals <= hist_info$gate_max
    keep[is.na(keep)] <- FALSE
    if (sum(keep) < 10L) return(NULL)
    list(
        final_gated_data = gated_data[keep, , drop = FALSE],
        hist_info = hist_info
    )
}

.reference_histogram_negative_source <- function(gated_data,
                                                  peak_channel,
                                                  filename,
                                                  sample_type,
                                                  manual_gates,
                                                  detector_names,
                                                  fsc,
                                                  ssc,
                                                  hist_info = NULL) {
    if (!all(c(peak_channel, detector_names, fsc, ssc) %in% colnames(gated_data))) return(NULL)

    gate <- .reference_manual_gate(manual_gates, "negative", filename, sample_type)
    verts <- .reference_gate_vertices(gate)
    keep <- NULL
    source <- "lower_tail"
    if (!is.null(gate) && !is.null(verts) && nrow(verts) >= 2L &&
        identical(as.character(gate$plot_mode[1]), "negative_1d")) {
        limits <- range(verts$x, na.rm = TRUE)
        if (all(is.finite(limits)) && limits[2] > limits[1]) {
            keep <- gated_data[, peak_channel] >= limits[1] & gated_data[, peak_channel] <= limits[2]
            source <- "manual_negative_gate"
        }
    }
    if (is.null(keep) && !is.null(hist_info) && !is.null(hist_info$vals_log)) {
        neg_min <- as.numeric(attr(hist_info$vals_log, "neg_raw_min", exact = TRUE))[1]
        neg_max <- as.numeric(attr(hist_info$vals_log, "neg_raw_max", exact = TRUE))[1]
        if (!is.finite(neg_min) || !is.finite(neg_max) || neg_max <= neg_min) {
            neg_min_log <- as.numeric(attr(hist_info$vals_log, "neg_log_min", exact = TRUE))[1]
            neg_max_log <- as.numeric(attr(hist_info$vals_log, "neg_log_max", exact = TRUE))[1]
            if (is.finite(neg_min_log) && is.finite(neg_max_log) && neg_max_log > neg_min_log) {
                neg_min <- 10^neg_min_log
                neg_max <- 10^neg_max_log
            }
        }
        if (is.finite(neg_min) && is.finite(neg_max) && neg_max > neg_min) {
            keep <- gated_data[, peak_channel] >= neg_min & gated_data[, peak_channel] <= neg_max
            source <- "automatic_negative_gate"
        }
    }
    if (is.null(keep) || sum(keep, na.rm = TRUE) < 10L) {
        peak_vals <- gated_data[, peak_channel]
        finite <- is.finite(peak_vals)
        if (sum(finite) < 10L) return(NULL)
        cutoff <- as.numeric(stats::quantile(peak_vals[finite], 0.15, na.rm = TRUE))
        keep <- finite & peak_vals <= cutoff
        source <- "lower_tail"
    }
    keep[is.na(keep)] <- FALSE
    negative_events <- gated_data[keep, , drop = FALSE]
    if (nrow(negative_events) < 10L) return(NULL)

    negative <- apply(negative_events[, detector_names, drop = FALSE], 2, stats::median, na.rm = TRUE)
    background <- .scc_background_from_gated_af_list(
        af_gated_list = list(list(
            events = negative_events[, detector_names, drop = FALSE],
            scatter = negative_events[, c(fsc, ssc), drop = FALSE],
            scatter_names = c(fsc, ssc)
        )),
        detector_names = detector_names
    )
    if (!is.null(background)) attr(negative, "scc_background") <- background
    attr(negative, "source") <- source
    negative
}

# Computes the 2D ellipse coordinates at a given confidence level.
# Used for gating populations in FSC/SSC scatter plots based on GMM component mean and variance.
# Returns a data.table containing the ellipse points.
.get_reference_ellipse <- function(mean, sigma, level = 0.95, n = 100, scale = 1.0) {
    if (length(mean) < 2 || any(!is.finite(mean)) ||
        !is.matrix(sigma) || any(dim(sigma) < 2) || any(!is.finite(sigma))) {
        return(NULL)
    }
    chi2_val <- qchisq(level, df = 2)
    eig <- eigen(sigma)
    eig_values <- pmax(Re(eig$values), 0)
    if (length(eig_values) < 2 || any(!is.finite(eig_values)) || sum(eig_values > 0) == 0) {
        return(NULL)
    }
    a <- sqrt(eig_values[1] * chi2_val) * scale
    b <- sqrt(eig_values[2] * chi2_val) * scale
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
        ellipses <- lapply(populations, function(k) {
            .get_reference_ellipse(gmm_result$means[, k], gmm_result$sigmas[[k]], level, scale = scale)
        })
        ellipses <- ellipses[!vapply(ellipses, is.null, logical(1))]
        if (length(ellipses) == 0) {
            return(NULL)
        }
        all_pts <- data.table::rbindlist(
            ellipses
        )
        hull_idx <- grDevices::chull(all_pts$x, all_pts$y)
        ell <- all_pts[hull_idx, ]
    }
    if (is.null(ell) || nrow(ell) == 0) {
        return(NULL)
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
# Uses the median normalized shape for one AF band and k-means on scale-free
# spectral shapes for multi-band AF banks.
# Each center is normalized by its peak detector to yield a relative signature from 0 to 1.
# Returns a list with the raw median spectrum and a matrix of normalized AF basis signatures.
.extract_reference_af_profiles <- function(ff_af = NULL,
                                           detector_names,
                                           n_bands = 10,
                                           max_cells = 50000,
                                           af_events = NULL,
                                           min_cluster_events = 20,
                                           min_cluster_proportion = 0.005) {
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

    if (nrow(af_events) > max_cells) {
        af_events <- af_events[sample.int(nrow(af_events), max_cells), , drop = FALSE]
    }

    raw_median <- apply(af_events, 2, stats::median, na.rm = TRUE)

    af_pos <- pmax(af_events, 0)
    row_scale <- apply(af_pos, 1, max, na.rm = TRUE)
    keep <- is.finite(row_scale) & row_scale > 0
    if (!any(keep)) {
        sig <- pmax(raw_median, 0)
        sig_max <- max(sig, na.rm = TRUE)
        if (is.finite(sig_max) && sig_max > 0) {
            sig <- sig / sig_max
        } else {
            sig <- rep(0, length(sig))
            names(sig) <- detector_names
        }
        sig_mat <- matrix(sig, nrow = 1)
        rownames(sig_mat) <- "AF"
        colnames(sig_mat) <- detector_names
        return(list(raw_median = raw_median, signatures = sig_mat))
    }

    af_shape <- af_pos[keep, , drop = FALSE] / row_scale[keep]
    n_bands <- suppressWarnings(as.integer(n_bands[1]))
    if (!is.finite(n_bands) || is.na(n_bands) || n_bands < 1) {
        stop("n_bands must be an integer >= 1.")
    }
    km <- .reference_kmeans_af_centers(
        af_shape = af_shape,
        n_centers = n_bands,
        min_cluster_events = min_cluster_events,
        min_cluster_proportion = min_cluster_proportion
    )
    centers <- .reference_normalize_af_centers(km$centers)
    if (ncol(centers) != length(detector_names) && nrow(centers) == length(detector_names)) {
        centers <- t(centers)
    }
    if (ncol(centers) != length(detector_names)) {
        stop(
            "AF center detector count mismatch: expected ",
            length(detector_names),
            " detector(s), got ",
            ncol(centers),
            ".",
            call. = FALSE
        )
    }
    selection <- list(
        method = if (identical(km$center_method, "median")) "median_fixed" else "kmeans_fixed",
        n_bands = nrow(centers),
        requested_bands = n_bands,
        raw_center_count = nrow(km$centers),
        final_bands = nrow(centers),
        cluster_sizes = km$cluster_sizes,
        min_cluster_size = km$min_cluster_size
    )
    rownames(centers) <- c("AF", if (nrow(centers) > 1) paste0("AF_", seq.int(2, nrow(centers))) else NULL)
    colnames(centers) <- detector_names

    list(raw_median = raw_median, signatures = centers, selection = selection)
}

# Determines the name of the autofluorescence (unstained) control file.
# Looks it up in the control mapping or falls back to identifying files matching AF file naming heuristics.
# Returns the AF file basename (without extension) or NULL.
.resolve_reference_af_name <- function(control_df, fcs_files) {
    af_fn <- NULL
    if (!is.null(control_df)) {
        af_rows <- if ("fluorophore" %in% colnames(control_df)) {
            control_df[.is_primary_af_control_row(
                fluorophore = control_df$fluorophore,
                marker = if ("marker" %in% colnames(control_df)) control_df$marker else NULL,
                filename = control_df$filename
            ), ]
        } else {
            data.frame()
        }
        if (nrow(af_rows) > 0) {
            control_type <- if ("control.type" %in% colnames(af_rows)) {
                tolower(trimws(as.character(af_rows$control.type)))
            } else {
                rep("", nrow(af_rows))
            }
            bead_negative <- vapply(af_rows$filename, .reference_is_bead_negative_file, logical(1))
            af_rows <- af_rows[control_type != "beads" & !bead_negative, , drop = FALSE]
        }
        if (nrow(af_rows) > 0) af_fn <- tools::file_path_sans_ext(basename(af_rows$filename[1]))
    }
    if (is.null(af_fn)) {
        af_idx_tmp <- which(.is_af_filename(fcs_files) & !vapply(fcs_files, .reference_is_bead_negative_file, logical(1)))
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

.reference_external_negative_summary <- function(detector_names,
                                                 af_data_raw = NULL,
                                                 scc_background = NULL) {
    bg <- if (!is.null(scc_background) && !is.null(scc_background$spectra)) {
        as.matrix(scc_background$spectra[, detector_names, drop = FALSE])
    } else {
        NULL
    }
    if (!is.null(bg) && nrow(bg) > 0L) {
        bg <- bg[stats::complete.cases(bg), , drop = FALSE]
        bg <- pmax(bg, 0)
        bg <- bg[rowSums(bg, na.rm = TRUE) > 0, , drop = FALSE]
        if (nrow(bg) > 0L) {
            return(list(
                mean = apply(bg, 2, mean, na.rm = TRUE),
                median = apply(bg, 2, stats::median, na.rm = TRUE),
                spectra = bg
            ))
        }
    }

    if (!is.null(af_data_raw) && all(detector_names %in% names(af_data_raw))) {
        af_vec <- pmax(as.numeric(af_data_raw[detector_names]), 0)
        names(af_vec) <- detector_names
        if (sum(af_vec^2, na.rm = TRUE) > 0) {
            return(list(
                mean = af_vec,
                median = af_vec,
                spectra = matrix(af_vec, nrow = 1L, dimnames = list(NULL, detector_names))
            ))
        }
    }

    NULL
}

.reference_cosine_to_vector <- function(mat, vec) {
    mat <- pmax(as.matrix(mat), 0)
    vec <- pmax(as.numeric(vec), 0)
    denom <- sqrt(rowSums(mat^2, na.rm = TRUE)) * sqrt(sum(vec^2, na.rm = TRUE))
    out <- rep(NA_real_, nrow(mat))
    ok <- is.finite(denom) & denom > 0
    if (any(ok)) {
        out[ok] <- as.numeric(mat[ok, , drop = FALSE] %*% vec) / denom[ok]
    }
    out
}

.reference_background_scatter_gate <- function(raw_data,
                                               pd,
                                               sample_type,
                                               filename,
                                               config,
                                               auto_sample_type = sample_type) {
    scatter_info <- .apply_reference_manual_scatter_gates(
        raw_data = raw_data,
        sample_type = sample_type,
        filename = filename,
        manual_gates = .get_reference_config_value(config, "manual_gates", NULL)
    )
    if (!is.null(scatter_info)) return(scatter_info)

    .compute_reference_scatter_gate(
        raw_data = raw_data,
        pd = pd,
        sample_type = auto_sample_type,
        outlier_percentile = .get_reference_config_value(config, "outlier_percentile", 0.02),
        debris_percentile = .get_reference_config_value(config, "debris_percentile", 0.08),
        subsample_n = .get_reference_config_value(config, "subsample_n", 5000),
        max_clusters = .get_reference_config_value(config, "max_clusters", 10),
        min_cluster_proportion = .get_reference_config_value(config, "min_cluster_proportion", 0.03),
        gate_contour_beads = .get_reference_config_value(config, "gate_contour_beads", 0.95),
        gate_contour_cells = .get_reference_config_value(config, "gate_contour_cells", 0.90),
        bead_gate_scale = .get_reference_config_value(config, "bead_gate_scale", 1.3)
    )
}

# Reads and gates an AF / unstained control FCS file.
# Applies scatter-gating on FSC-SSC to isolate the main cell/bead population from debris.
# Returns a list containing the gated event data across all spectral detectors and gating metadata.
.extract_reference_af_gated_events <- function(fcs_file,
                                               detector_names,
                                               config,
                                               sample_type = "cells") {
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

    sample_type <- if (identical(tolower(trimws(as.character(sample_type)[1])), "beads")) "beads" else "cells"
    scatter_info <- .reference_background_scatter_gate(
        raw_data = raw_af,
        pd = pd_af,
        sample_type = sample_type,
        filename = basename(fcs_file),
        config = config,
        auto_sample_type = if (identical(sample_type, "beads")) "beads" else "unstained"
    )
    if (is.null(scatter_info)) {
        .spectreasy_console_step("Skip AF", paste0(basename(fcs_file), " has too few scatter-gated cells"))
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
                                           af_max_cells,
                                           af_min_cluster_events,
                                           af_min_cluster_proportion,
                                           fcs_files_all = fcs_files,
                                           config = NULL) {
    af_data_raw <- NULL
    af_signatures_norm <- NULL
    af_bank_info <- NULL
    scc_background <- NULL
    af_events <- NULL
    af_fn <- .resolve_reference_af_name(control_df = control_df, fcs_files = fcs_files)

    af_paths <- character()
    af_source_types <- character()
    fcs_keys <- tools::file_path_sans_ext(basename(fcs_files_all))
    if (!is.null(control_df) && nrow(control_df) > 0) {
        af_rows <- .is_primary_af_control_row(
            fluorophore = if ("fluorophore" %in% colnames(control_df)) control_df$fluorophore else NULL,
            marker = if ("marker" %in% colnames(control_df)) control_df$marker else NULL,
            filename = control_df$filename
        )
        if (any(af_rows)) {
            af_df <- control_df[af_rows, , drop = FALSE]
            control_type <- if ("control.type" %in% colnames(af_df)) {
                tolower(trimws(as.character(af_df$control.type)))
            } else {
                rep("", nrow(af_df))
            }
            bead_negative <- vapply(af_df$filename, .reference_is_bead_negative_file, logical(1))
            af_df <- af_df[control_type != "beads" & !bead_negative, , drop = FALSE]
            if (nrow(af_df) > 0) {
                af_file_keys <- tools::file_path_sans_ext(basename(as.character(af_df$filename)))
                hits <- match(af_file_keys, fcs_keys)
                hits <- hits[!is.na(hits)]
                if (length(hits) > 0) {
                    af_paths <- c(af_paths, fcs_files_all[hits])
                    af_source_types <- c(af_source_types, rep("mapped_unstained", length(hits)))
                }
            }
        }
    }
    if (length(af_paths) == 0 && !is.null(af_fn)) {
        af_path <- fcs_files_all[match(af_fn, fcs_keys)]
        af_path <- af_path[!is.na(af_path)]
        if (length(af_path) > 0) {
            af_paths <- c(af_paths, af_path[1])
            af_source_types <- c(af_source_types, "primary_unstained")
        }
    }
    keep_unique <- !duplicated(normalizePath(af_paths, mustWork = FALSE))
    af_paths <- af_paths[keep_unique]
    af_source_types <- af_source_types[keep_unique]

    if (length(af_paths) > 0) {
        af_gated_list <- lapply(
            af_paths,
            .extract_reference_af_gated_events,
            detector_names = detector_names,
            config = config,
            sample_type = "cells"
        )
        keep_gated <- vapply(af_gated_list, function(x) !is.null(x) && !is.null(x$events) && nrow(x$events) > 0, logical(1))
        af_gated_list <- af_gated_list[keep_gated]
        af_source_types <- af_source_types[keep_gated]
        if (length(af_gated_list) == 0) {
            return(list(
                af_data_raw = af_data_raw,
                af_signatures_norm = af_signatures_norm,
                af_bank_info = af_bank_info,
                scc_background = scc_background,
                af_events = af_events
            ))
        }

        af_events <- do.call(rbind, lapply(af_gated_list, `[[`, "events"))
        scc_background <- .scc_background_from_gated_af_list(
            af_gated_list = af_gated_list,
            detector_names = detector_names
        )
        n_af_sources <- length(af_gated_list)
        requested_bands <- af_n_bands
        af_profiles <- .extract_reference_af_profiles(
            detector_names = detector_names,
            n_bands = requested_bands,
            max_cells = af_max_cells,
            af_events = af_events,
            min_cluster_events = af_min_cluster_events,
            min_cluster_proportion = af_min_cluster_proportion
        )
        af_data_raw <- af_profiles$raw_median
        af_signatures_norm <- af_profiles$signatures
        af_sources <- data.table::rbindlist(lapply(af_gated_list, `[[`, "source"))
        af_sources$source_type <- af_source_types
        data.table::setcolorder(af_sources, c("file", "source_type", "n_total", "n_scatter_gated", "scatter_gate_pct", "fsc_channel", "ssc_channel", "path"))
        af_bank_info <- list(
            source_count = n_af_sources,
            sources = af_sources,
            pooled_events = nrow(af_events),
            requested_bands = requested_bands,
            derived_bands = if (!is.null(af_signatures_norm)) nrow(af_signatures_norm) else 0L,
            af_min_cluster_events = af_min_cluster_events,
            af_min_cluster_proportion = af_min_cluster_proportion,
            selection = af_profiles$selection,
            mode = if (n_af_sources > 1) "pooled_af_sources" else "single_af"
        )
        if (!is.null(af_signatures_norm)) {
            msg <- if (n_af_sources == 1) {
                "primary unstained control"
            } else {
                paste0(n_af_sources, " pooled AF control files")
            }
            selection_msg <- if (!is.null(af_profiles$selection)) {
                paste0(" (selected by ", af_profiles$selection$method, ")")
            } else {
                ""
            }
            .spectreasy_console_field(
                "AF bank",
                paste0(nrow(af_signatures_norm), " signature(s) from ", msg, selection_msg)
            )
        }
    }

    list(
        af_data_raw = af_data_raw,
        af_signatures_norm = af_signatures_norm,
        af_bank_info = af_bank_info,
        scc_background = scc_background,
        af_events = af_events
    )
}

# Validates raw FCS data for corruption or extreme values.
# Ensures the dataset contains no infinite values, fewer than 10% NA values, and no extreme
# values exceeding 1e9 which indicate file corruption.
# Returns TRUE if valid, otherwise FALSE.
.validate_reference_raw_data <- function(raw_data, sn) {
    if (any(is.infinite(raw_data))) {
        .spectreasy_console_step("Skip SCC", paste0(sn, " has infinite values"))
        return(FALSE)
    }
    na_prop <- sum(is.na(raw_data)) / length(raw_data)
    if (na_prop > 0.1) {
        .spectreasy_console_step("Skip SCC", paste0(sn, " has too many NAs (", round(na_prop * 100, 1), "%)"))
        return(FALSE)
    }
    max_val <- max(raw_data, na.rm = TRUE)
    if (max_val > 1e9) {
        .spectreasy_console_step("Skip SCC", paste0(sn, " has extreme values (max > 1e9)"))
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

.active_reference_control_rows <- function(control_df, fcs_files_all) {
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
    active <- active & !is_af

    control_df[active, , drop = FALSE]
}

.validate_reference_complete_controls <- function(control_df, fcs_files_all, processed_results) {
    active_rows <- .active_reference_control_rows(
        control_df = control_df,
        fcs_files_all = fcs_files_all
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
    if (is.null(final_gate) || nrow(final_gate) == 0) return(NULL)
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
.select_reference_peak_channel <- function(gated_data, detector_names, row_info, channel_alias_map, sn_ext, sn, cytometer = "auto") {
    is_unstained <- grepl("unstained|autofluorescence|\\bAF\\b", paste(sn_ext, sn), ignore.case = TRUE)
    if (nrow(row_info) > 0) {
        is_unstained <- is_unstained || .is_af_control_row(
            fluorophore = if ("fluorophore" %in% colnames(row_info)) row_info$fluorophore[1] else "",
            marker = if ("marker" %in% colnames(row_info)) row_info$marker[1] else "",
            filename = sn_ext
        )
    }

    q999_by_channel <- apply(
        gated_data[, detector_names, drop = FALSE],
        2,
        function(x) stats::quantile(x, 0.999, na.rm = TRUE)
    )
    if (is_unstained) {
        med_by_channel <- apply(gated_data[, detector_names, drop = FALSE], 2, stats::median, na.rm = TRUE)
        return(list(peak_channel = names(which.max(med_by_channel)), q999_by_channel = q999_by_channel))
    }

    inferred_peak_channel <- detector_names[which.max(q999_by_channel)]
    peak_channel <- inferred_peak_channel

    expected_channel <- ""
    cytometer_id <- .resolve_cytometer_id(cytometer, allow_auto = TRUE, unknown_as_auto = TRUE)
    if (!identical(cytometer_id, "auto") && nrow(row_info) > 0 && "fluorophore" %in% colnames(row_info)) {
        expected_ref <- tryCatch(.load_control_file_shipped_reference(cytometer_id), error = function(e) NULL)
        if (!is.null(expected_ref)) {
            expected_raw <- .control_file_expected_channel_for_fluor(row_info$fluorophore[1], expected_ref$fluor_peak_channel_map)
            expected_channel <- .resolve_reference_control_channel(expected_raw, detector_names, channel_alias_map = channel_alias_map)
        }
    }
    if (nzchar(expected_channel)) {
        peak_channel <- expected_channel
        return(list(peak_channel = peak_channel, q999_by_channel = q999_by_channel))
    }

    if (nrow(row_info) > 0 && !is.na(row_info$channel[1]) && row_info$channel[1] != "") {
        resolved_channel <- .resolve_reference_control_channel(row_info$channel[1], detector_names, channel_alias_map = channel_alias_map)
        if (nzchar(resolved_channel)) {
            peak_channel <- resolved_channel
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

.compute_reference_autospectral_scc <- function(clean_data,
                                                detector_names,
                                                peak_channel,
                                                sample_type,
                                                af_data_raw = NULL,
                                                scc_background = NULL,
                                                n_candidates = 1000L,
                                                n_spectral = 200L,
                                                min_events = 10L,
                                                scc_background_k = 2L) {
    if (!all(c(detector_names, peak_channel) %in% colnames(clean_data))) {
        return(NULL)
    }
    event_mat <- as.matrix(clean_data[, detector_names, drop = FALSE])
    complete <- stats::complete.cases(event_mat)
    peak_vals <- clean_data[, peak_channel]
    valid <- complete & is.finite(peak_vals) & peak_vals > 0 & rowSums(event_mat, na.rm = TRUE) > 0
    if (sum(valid, na.rm = TRUE) < min_events) {
        return(NULL)
    }

    n_candidates <- suppressWarnings(as.integer(n_candidates[1]))
    if (!is.finite(n_candidates) || is.na(n_candidates) || n_candidates < 1L) {
        n_candidates <- 1000L
    }
    n_spectral <- suppressWarnings(as.integer(n_spectral[1]))
    if (!is.finite(n_spectral) || is.na(n_spectral) || n_spectral < 1L) {
        n_spectral <- 200L
    }
    scc_background_k <- suppressWarnings(as.integer(scc_background_k[1]))
    if (!is.finite(scc_background_k) || is.na(scc_background_k) || scc_background_k < 1L) {
        scc_background_k <- 2L
    }

    vals_log <- log10(pmax(peak_vals, 1))
    attr(vals_log, "positive_gate_present") <- TRUE
    attr(vals_log, "negative_gate_present") <- FALSE
    attr(vals_log, "mapped_peak_channel") <- peak_channel

    af_summary <- .reference_external_negative_summary(
        detector_names = detector_names,
        af_data_raw = af_data_raw,
        scc_background = scc_background
    )
    valid_idx <- which(valid)
    candidate_n <- min(n_candidates, length(valid_idx))
    candidate_idx <- valid_idx[order(peak_vals[valid_idx], decreasing = TRUE)[seq_len(candidate_n)]]

    if (!is.null(af_summary)) {
        af_mean <- af_summary$mean[detector_names]
        af_median <- af_summary$median[detector_names]
        af_norm <- sqrt(sum(af_mean^2, na.rm = TRUE))
        if (!is.finite(af_norm) || af_norm <= 0) {
            return(NULL)
        }
        event_complete <- pmax(event_mat[complete, , drop = FALSE], 0)
        af_unit <- af_mean / (af_norm + 1e-9)
        projection <- as.numeric(event_complete %*% af_unit)
        residual <- pmax(event_complete - projection %o% af_unit, 0)
        colnames(residual) <- detector_names
        empirical_peak <- names(which.max(colMeans(residual, na.rm = TRUE)))

        af_cosine_candidates <- .reference_cosine_to_vector(
            pmax(event_mat[candidate_idx, , drop = FALSE], 0),
            af_median
        )
        cosine_ok <- is.finite(af_cosine_candidates)
        if (sum(cosine_ok) < min_events) {
            return(NULL)
        }
        candidate_idx <- candidate_idx[cosine_ok]
        af_cosine_candidates <- af_cosine_candidates[cosine_ok]
        selected_n <- min(n_spectral, length(candidate_idx))
        selected_local <- order(af_cosine_candidates, decreasing = FALSE)[seq_len(selected_n)]
        selected_idx <- candidate_idx[selected_local]
        if (length(selected_idx) < min_events) {
            return(NULL)
        }

        selected_events <- clean_data[selected_idx, , drop = FALSE]
        matched_background <- .scc_background_match(
            events = selected_events,
            background = scc_background,
            k = scc_background_k
        )
        if (is.null(matched_background)) {
            matched_background <- matrix(
                af_median,
                nrow = nrow(selected_events),
                ncol = length(detector_names),
                byrow = TRUE,
                dimnames = list(NULL, detector_names)
            )
            background_method <- "external_median"
        } else {
            background_method <- "scatter_knn"
        }
        corrected_events <- as.matrix(selected_events[, detector_names, drop = FALSE]) - matched_background
        colnames(corrected_events) <- detector_names
        positive_events <- pmax(corrected_events, 0)
        pos_spectrum_raw <- apply(positive_events, 2, stats::median, na.rm = TRUE)
        sig_pure <- pmax(pos_spectrum_raw, 0)
        max_val <- max(sig_pure, na.rm = TRUE)
        if (!is.finite(max_val) || max_val <= 0) {
            return(NULL)
        }
        spectrum_norm <- sig_pure / max_val

        selected_cos <- af_cosine_candidates[selected_local]
        attr(vals_log, "gate_type") <- "autospectral_external"
        attr(vals_log, "gate_method") <- paste0(
            "AutoSpectral-style external-negative selector after FSC/SSC gate: top ",
            length(candidate_idx),
            " peak-bright candidate event(s), kept ",
            length(selected_idx),
            " least-AF-like event(s)"
        )
        attr(vals_log, "af_cosine_max_selected") <- max(selected_cos, na.rm = TRUE)
        attr(vals_log, "af_cosine_median_selected") <- stats::median(selected_cos, na.rm = TRUE)
        attr(vals_log, "af_score_valid_events") <- length(valid_idx)
        attr(vals_log, "empirical_peak_channel") <- empirical_peak
        attr(vals_log, "scc_background_method") <- background_method

        positive_idx <- rep(FALSE, nrow(clean_data))
        positive_idx[selected_idx] <- TRUE
        af_cosine_full <- rep(NA_real_, nrow(clean_data))
        af_cosine_full[candidate_idx] <- af_cosine_candidates
        attr(spectrum_norm, "variance") <- apply(positive_events, 2, stats::var, na.rm = TRUE) / (max_val^2)
        attr(spectrum_norm, "scc_background") <- list(
            method = background_method,
            k = if (identical(background_method, "scatter_knn")) scc_background_k else NA_integer_,
            matched_events = nrow(positive_events)
        )
        attr(spectrum_norm, "scc_positive_events") <- positive_events

        return(list(
            spectrum = spectrum_norm,
            final_gated_data = selected_events,
            positive_events = positive_events,
            positive_idx = positive_idx,
            vals_log = vals_log,
            gate_min = min(peak_vals[selected_idx], na.rm = TRUE),
            gate_max = max(peak_vals[selected_idx], na.rm = TRUE),
            peak_vals = peak_vals,
            spectral_gate_info = list(
                af_basis = af_summary$spectra,
                af_cosine = af_cosine_full
            ),
            extraction_method = "autospectral_external",
            n_candidates = length(candidate_idx),
            n_selected = length(selected_idx)
        ))
    }

    top_n <- min(100L, length(candidate_idx))
    selected_idx <- candidate_idx[seq_len(top_n)]
    bottom_n <- max(10L, floor(0.10 * length(valid_idx)))
    bottom_n <- min(bottom_n, length(valid_idx) - length(selected_idx))
    if (bottom_n < min_events || length(selected_idx) < min_events) {
        return(NULL)
    }
    negative_pool <- setdiff(valid_idx[order(peak_vals[valid_idx], decreasing = FALSE)], selected_idx)
    negative_idx <- head(negative_pool, bottom_n)
    neg_spectrum <- apply(event_mat[negative_idx, , drop = FALSE], 2, mean, na.rm = TRUE)
    corrected_events <- sweep(as.matrix(clean_data[selected_idx, detector_names, drop = FALSE]), 2, neg_spectrum, "-")
    colnames(corrected_events) <- detector_names
    positive_events <- pmax(corrected_events, 0)
    sig_pure <- pmax(colMeans(positive_events, na.rm = TRUE), 0)
    max_val <- max(sig_pure, na.rm = TRUE)
    if (!is.finite(max_val) || max_val <= 0) {
        return(NULL)
    }
    spectrum_norm <- sig_pure / max_val
    positive_idx <- rep(FALSE, nrow(clean_data))
    positive_idx[selected_idx] <- TRUE
    negative_bool <- rep(FALSE, nrow(clean_data))
    negative_bool[negative_idx] <- TRUE
    attr(vals_log, "gate_type") <- "autospectral_internal"
    attr(vals_log, "gate_method") <- paste0(
        "AutoSpectral-style internal-negative fallback after FSC/SSC gate: top ",
        length(selected_idx),
        " peak-bright event(s) minus bottom ",
        length(negative_idx),
        " event(s)"
    )
    attr(vals_log, "negative_idx") <- negative_bool
    attr(vals_log, "scc_background_method") <- "internal_negative"
    attr(spectrum_norm, "variance") <- apply(positive_events, 2, stats::var, na.rm = TRUE) / (max_val^2)
    attr(spectrum_norm, "scc_background") <- list(
        method = "internal_negative",
        k = NA_integer_,
        matched_events = 0L
    )
    attr(spectrum_norm, "scc_positive_events") <- positive_events

    list(
        spectrum = spectrum_norm,
        final_gated_data = clean_data[selected_idx, , drop = FALSE],
        positive_events = positive_events,
        positive_idx = positive_idx,
        vals_log = vals_log,
        gate_min = min(peak_vals[selected_idx], na.rm = TRUE),
        gate_max = max(peak_vals[selected_idx], na.rm = TRUE),
        peak_vals = peak_vals,
        spectral_gate_info = NULL,
        extraction_method = "autospectral_internal",
        n_candidates = length(candidate_idx),
        n_selected = length(selected_idx)
    )
}

# Computes the normalized spectral signature for a stained control sample.
# Calculates the median intensity in each detector for the positive and negative gates,
# subtracts the negative/background control, and normalizes the spectrum by the peak signal.
# Returns the normalized spectrum vector.
.compute_reference_spectrum <- function(final_gated_data,
                                        gated_data,
                                        peak_vals,
                                        vals_log,
                                        detector_names,
                                        row_info,
                                        sample_type = "beads",
                                        af_data_raw = NULL,
                                        universal_negatives = NULL,
                                        bead_negative = NULL) {
    pos_spectrum_raw <- apply(final_gated_data[, detector_names, drop = FALSE], 2, median, na.rm = TRUE)
    neg_log_min <- attr(vals_log, "neg_log_min")
    neg_log_max <- attr(vals_log, "neg_log_max")
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
    } else if (identical(sample_type, "beads") && !is.null(bead_negative)) {
        bead_negative
    } else {
        neg_spectrum_raw
    }
    sig_pure <- pmax(pos_spectrum_raw - final_neg, 0)
    max_val <- max(sig_pure, na.rm = TRUE)
    if (max_val <= 0) max_val <- max(pos_spectrum_raw, na.rm = TRUE)
    res <- sig_pure / max_val

    res
}

# Standardizes a file path to its base name without file extension for matching.
# Returns the cleaned character string/vector.
.reference_negative_key <- function(x) {
    tools::file_path_sans_ext(basename(trimws(as.character(x))))
}

.reference_is_bead_negative_file <- function(x) {
    stem <- tools::file_path_sans_ext(basename(trimws(as.character(x))))
    stem_lower <- tolower(stem)
    stem_norm <- gsub("[^[:alnum:]]+", "", stem_lower)
    tokens <- unlist(strsplit(stem_lower, "[^[:alnum:]]+"))
    tokens <- tokens[nzchar(tokens)]

    has_bead <- any(tokens %in% c("bead", "beads", "compbead", "compbeads")) ||
        grepl("compbeads?", stem_norm) ||
        grepl("beads?", stem_norm)
    has_negative <- any(tokens %in% c(
        "unstained", "unstain", "us", "ut", "usut", "usut1", "neg",
        "negative", "background", "bg", "blank", "minus"
    )) ||
        grepl("us[_ -]?ut", stem, ignore.case = TRUE) ||
        grepl("unstained|negative|background|blank", stem, ignore.case = TRUE) ||
        grepl("(^|[^[:alnum:]])(?:us|neg|bg)(?:[^[:alnum:]]|$)", stem, ignore.case = TRUE, perl = TRUE)

    has_bead && has_negative
}

.collect_reference_unstained_bead_negative <- function(control_df,
                                                       fcs_files,
                                                       detector_names,
                                                       config) {
    if (length(fcs_files) == 0) {
        return(NULL)
    }

    file_keys <- .reference_negative_key(fcs_files)
    names(fcs_files) <- file_keys
    bead_keys <- character()

    if (!is.null(control_df) && is.data.frame(control_df) && nrow(control_df) > 0 && "filename" %in% colnames(control_df)) {
        is_af <- .is_af_control_row(
            fluorophore = if ("fluorophore" %in% colnames(control_df)) control_df$fluorophore else NULL,
            marker = if ("marker" %in% colnames(control_df)) control_df$marker else NULL,
            filename = control_df$filename
        )
        control_type <- if ("control.type" %in% colnames(control_df)) {
            tolower(trimws(as.character(control_df$control.type)))
        } else {
            rep("", nrow(control_df))
        }
        file_is_bead <- grepl("beads?", basename(as.character(control_df$filename)), ignore.case = TRUE)
        file_is_bead_negative <- vapply(control_df$filename, .reference_is_bead_negative_file, logical(1))
        idx <- which((is_af & (control_type == "beads" | file_is_bead)) | file_is_bead_negative)
        if (length(idx) > 0) {
            bead_keys <- .reference_negative_key(control_df$filename[idx])
        }
    }

    if (length(bead_keys) == 0) {
        bead_keys <- file_keys[vapply(fcs_files, .reference_is_bead_negative_file, logical(1))]
    }
    bead_keys <- unique(bead_keys[nzchar(bead_keys) & bead_keys %in% names(fcs_files)])
    if (length(bead_keys) == 0) {
        return(NULL)
    }

    bead_events <- list()
    bead_gated_list <- list()
    for (key in bead_keys) {
        fcs_file <- unname(fcs_files[[key]])

        ff <- tryCatch(
            flowCore::read.FCS(fcs_file, transformation = FALSE, truncate_max_range = FALSE),
            error = function(e) NULL
        )
        if (is.null(ff)) {
            warning("Could not read unstained bead negative file: ", fcs_file)
            next
        }

        pd <- flowCore::pData(flowCore::parameters(ff))
        raw_data <- flowCore::exprs(ff)
        scatter_info <- .reference_background_scatter_gate(
            raw_data = raw_data,
            pd = pd,
            sample_type = "beads",
            filename = basename(fcs_file),
            config = config
        )

        neg_data <- if (!is.null(scatter_info)) scatter_info$gated_data else raw_data
        if (nrow(neg_data) > 0) {
            bead_events[[length(bead_events) + 1L]] <- neg_data[, detector_names, drop = FALSE]
            if (!is.null(scatter_info)) {
                bead_gated_list[[length(bead_gated_list) + 1L]] <- list(
                    events = neg_data[, detector_names, drop = FALSE],
                    scatter = neg_data[, c(scatter_info$fsc, scatter_info$ssc), drop = FALSE],
                    scatter_names = c(scatter_info$fsc, scatter_info$ssc)
                )
            }
            .spectreasy_console_field("Bead neg", basename(fcs_file))
        }
    }

    if (length(bead_events) == 0) {
        return(NULL)
    }

    bead_mat <- do.call(rbind, bead_events)
    bead_negative <- apply(bead_mat[, detector_names, drop = FALSE], 2, stats::median, na.rm = TRUE)
    bead_background <- .scc_background_from_gated_af_list(
        af_gated_list = bead_gated_list,
        detector_names = detector_names
    )
    if (!is.null(bead_background)) {
        attr(bead_negative, "scc_background") <- bead_background
    }
    bead_negative
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
        scatter_info <- .reference_background_scatter_gate(
            raw_data = raw_data,
            pd = pd,
            sample_type = sample_info$type,
            filename = basename(fcs_file),
            config = config
        )

        neg_data <- if (!is.null(scatter_info)) scatter_info$gated_data else raw_data
        negative <- apply(neg_data[, detector_names, drop = FALSE], 2, stats::median, na.rm = TRUE)
        if (!is.null(scatter_info)) {
            background <- .scc_background_from_gated_af_list(
                af_gated_list = list(list(
                    events = neg_data[, detector_names, drop = FALSE],
                    scatter = neg_data[, c(scatter_info$fsc, scatter_info$ssc), drop = FALSE],
                    scatter_names = c(scatter_info$fsc, scatter_info$ssc)
                )),
                detector_names = detector_names
            )
            if (!is.null(background)) attr(negative, "scc_background") <- background
        }
        out[[key]] <- negative
        .spectreasy_console_field("Negative", basename(fcs_file))
    }

    out
}

.resolve_reference_scc_negative_source <- function(row_info,
                                                   sample_type,
                                                   af_data_raw = NULL,
                                                   scc_background = NULL,
                                                   universal_negatives = NULL,
                                                   bead_negative = NULL) {
    if (identical(sample_type, "beads")) {
        return(list(
            negative = bead_negative,
            background = attr(bead_negative, "scc_background", exact = TRUE),
            source = "AF_beads"
        ))
    }

    uv_val <- if (nrow(row_info) > 0L && "universal.negative" %in% colnames(row_info)) {
        trimws(as.character(row_info$universal.negative[1]))
    } else {
        ""
    }
    uv_upper <- toupper(uv_val)
    uv_key <- .reference_negative_key(uv_val)
    named_negative <- if (
        nzchar(uv_key) &&
        !uv_upper %in% c("FALSE", "TRUE", "AF") &&
        !is.null(universal_negatives) &&
        uv_key %in% names(universal_negatives)
    ) {
        universal_negatives[[uv_key]]
    } else {
        NULL
    }

    if (!is.null(named_negative)) {
        return(list(
            negative = named_negative,
            background = attr(named_negative, "scc_background", exact = TRUE),
            source = uv_val
        ))
    }

    list(
        negative = af_data_raw,
        background = scc_background,
        source = "AF"
    )
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

.warn_reference_qc_plot_failure <- function(sn, plot_type, condition) {
    warning(
        "SCC QC plot skipped for ", sn, " (", plot_type, "): ",
        conditionMessage(condition),
        ". Continuing unmixing.",
        call. = FALSE
    )
    invisible(FALSE)
}

.save_reference_ggsave <- function(filename, plot, sn, plot_type, ...) {
    tryCatch(
        {
            ggplot2::ggsave(filename, plot, ...)
            invisible(TRUE)
        },
        error = function(e) .warn_reference_qc_plot_failure(sn, plot_type, e)
    )
}

.save_reference_qc_plots_safely <- function(...) {
    args <- list(...)
    sn <- if (!is.null(args$sn)) args$sn else "control"
    tryCatch(
        do.call(.save_reference_qc_plots, args),
        error = function(e) .warn_reference_qc_plot_failure(sn, "plot bundle", e)
    )
}

.reference_even_indices <- function(n, max_points = 50000L) {
    n <- as.integer(n)
    max_points <- suppressWarnings(as.integer(max_points))
    if (!is.finite(max_points) || is.na(max_points) || max_points <= 0L) max_points <- 50000L
    if (n <= max_points) return(seq_len(n))
    pmax(1L, pmin(n, floor((seq_len(max_points) - 1L) * n / max_points) + 1L))
}

.reference_density_palette <- function(size = 64L) {
    vapply(seq_len(size), function(i) {
        v <- (i - 1) / max(size - 1, 1)
        r <- max(0, min(1, min(4 * v - 1.5, -4 * v + 4.5)))
        g <- max(0, min(1, min(4 * v - 0.5, -4 * v + 3.5)))
        b <- max(0, min(1, min(4 * v + 0.5, -4 * v + 2.5)))
        grDevices::rgb(r, g, b, alpha = 0.70)
    }, character(1))
}

.reference_point_density <- function(x, y) {
    n <- length(x)
    if (n == 0) return(numeric())
    min_x <- min(x, na.rm = TRUE)
    max_x <- max(x, na.rm = TRUE)
    min_y <- min(y, na.rm = TRUE)
    max_y <- max(y, na.rm = TRUE)
    rx <- max(max_x - min_x, 1)
    ry <- max(max_y - min_y, 1)
    num_bins <- 160L
    grid_side <- num_bins + 1L
    bx <- pmax(0L, pmin(num_bins, floor(((x - min_x) / rx) * num_bins)))
    by <- pmax(0L, pmin(num_bins, floor(((y - min_y) / ry) * num_bins)))
    grid <- tabulate(by * grid_side + bx + 1L, nbins = grid_side * grid_side)
    radius <- 5L
    sigma <- 2
    offsets <- expand.grid(dx = -radius:radius, dy = -radius:radius)
    offsets$weight <- exp(-(offsets$dx^2 + offsets$dy^2) / (2 * sigma^2))
    densities <- numeric(n)
    for (i in seq_len(n)) {
        nx <- bx[i] + offsets$dx
        ny <- by[i] + offsets$dy
        ok <- nx >= 0L & nx <= num_bins & ny >= 0L & ny <= num_bins
        grid_idx <- ny[ok] * grid_side + nx[ok] + 1L
        densities[i] <- sum(grid[grid_idx] * offsets$weight[ok])
    }
    densities
}

.reference_pretty_k_label <- function(x) {
    x <- as.numeric(x)
    out <- ifelse(is.finite(x), ifelse(abs(x) >= 1000, paste0(round(x / 1000), "K"), as.character(round(x))), "")
    million <- is.finite(x) & abs(x) >= 1000000
    out[million] <- paste0(sub("\\.0$", "", sprintf("%.1f", x[million] / 1000000)), "M")
    out
}

.reference_gui_plot_aspect <- function() {
    (520 - 54 - 18) / (420 - 18 - 46)
}

.reference_gui_coord_ratio <- function(x_lim, y_lim) {
    x_span <- diff(as.numeric(x_lim[1:2]))
    y_span <- diff(as.numeric(y_lim[1:2]))
    if (!is.finite(x_span) || !is.finite(y_span) || x_span <= 0 || y_span <= 0) return(1)
    x_span / (y_span * .reference_gui_plot_aspect())
}

.reference_qc_scatter_density_plot <- function(data,
                                               x_channel,
                                               y_channel,
                                               pd,
                                               gate_vertices = NULL,
                                               title,
                                               subtitle = NULL,
                                               max_points = 50000L,
                                               point_size = 1.5,
                                               x_domain = NULL,
                                               y_domain = NULL,
                                               sample_type = NULL) {
    x_desc <- .get_reference_axis_label(x_channel, pd)
    y_desc <- .get_reference_axis_label(y_channel, pd)
    x_vals <- data[, x_channel]
    y_vals <- data[, y_channel]
    keep <- is.finite(x_vals) & is.finite(y_vals)
    x_vals <- x_vals[keep]
    y_vals <- y_vals[keep]
    idx <- .reference_even_indices(length(x_vals), max_points = max_points)
    plot_df <- data.frame(x = x_vals[idx], y = y_vals[idx])
    density <- .reference_point_density(plot_df$x, plot_df$y)
    palette <- .reference_density_palette()
    density[!is.finite(density)] <- 0
    max_density <- if (length(density) > 0) max(density, na.rm = TRUE) else 0
    bucket <- if (length(density) > 0 && is.finite(max_density) && max_density > 0) {
        pmax(1L, pmin(length(palette), floor((density / max_density) * (length(palette) - 1L)) + 1L))
    } else {
        rep(1L, nrow(plot_df))
    }
    bucket[!is.finite(bucket) | is.na(bucket)] <- 1L
    plot_df$color <- palette[bucket]
    x_lim <- .reference_gui_extent(x_vals)
    y_lim <- .reference_gui_extent(y_vals)
    if (!is.null(x_domain) && length(x_domain) >= 2L && all(is.finite(x_domain[1:2])) && x_domain[2] > x_domain[1]) {
        x_lim <- as.numeric(x_domain[1:2])
    }
    if (!is.null(y_domain) && length(y_domain) >= 2L && all(is.finite(y_domain[1:2])) && y_domain[2] > y_domain[1]) {
        y_lim <- as.numeric(y_domain[1:2])
    }
    coord_ratio <- .reference_gui_coord_ratio(x_lim, y_lim)
    point_size <- max(0.3, min(1.8, as.numeric(point_size) * 0.55))
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x, y)) +
        ggplot2::geom_point(color = plot_df$color, size = point_size, alpha = 0.95, stroke = 0) +
        ggplot2::labs(title = title, subtitle = subtitle, x = x_desc, y = y_desc) +
        ggplot2::scale_x_continuous(labels = .reference_pretty_k_label) +
        ggplot2::scale_y_continuous(labels = .reference_pretty_k_label) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            legend.position = "none",
            panel.grid = ggplot2::element_blank(),
            panel.background = ggplot2::element_rect(fill = "#fbfaf6", color = "#ded9cf", linewidth = 0.35),
            axis.line = ggplot2::element_line(color = "#6c746f", linewidth = 0.35),
            axis.ticks = ggplot2::element_line(color = "#6c746f", linewidth = 0.3),
            axis.title = ggplot2::element_text(color = "#626a65", face = "bold"),
            axis.text = ggplot2::element_text(color = "#626a65", face = "bold"),
            plot.title = ggplot2::element_text(color = "#17201c", face = "plain", size = 14),
            plot.subtitle = ggplot2::element_text(color = "#69716b", size = 10.6)
        ) +
        ggplot2::coord_fixed(ratio = coord_ratio, xlim = x_lim, ylim = y_lim, expand = FALSE)
    if (!is.null(gate_vertices) && nrow(gate_vertices) >= 3) {
        is_bead_gate <- identical(tolower(trimws(as.character(sample_type)[1])), "beads")
        gate_color <- if (is_bead_gate) "#56b4e9" else "#d65238"
        gate_alpha <- if (is_bead_gate) 0.16 else 0.16
        gate_closed <- rbind(gate_vertices, gate_vertices[1, , drop = FALSE])
        p <- p +
            ggplot2::geom_polygon(data = gate_vertices, ggplot2::aes(x, y), inherit.aes = FALSE, fill = gate_color, alpha = gate_alpha) +
            ggplot2::geom_path(data = gate_closed, ggplot2::aes(x, y), inherit.aes = FALSE, color = gate_color, linewidth = 1.05) +
            ggplot2::geom_point(data = gate_vertices, ggplot2::aes(x, y), inherit.aes = FALSE, shape = 21, size = 2.4, stroke = 0.9, color = "#d65238", fill = "#fbfaf6")
    }
    p
}

.reference_histogram_transform_values <- function(values, transform = "asinh") {
    transform <- tolower(trimws(as.character(transform)[1]))
    values <- as.numeric(values)
    if (identical(transform, "log10")) {
        return(log10(pmax(values, 1)))
    }
    if (identical(transform, "asinh")) {
        return(asinh(values / 150))
    }
    if (identical(transform, "biexponential")) {
        return(sign(values) * log10(1 + abs(values) / 50))
    }
    values
}

.reference_histogram_inverse_values <- function(values, transform = "asinh") {
    transform <- tolower(trimws(as.character(transform)[1]))
    values <- as.numeric(values)
    if (identical(transform, "log10")) {
        return(10^values)
    }
    if (identical(transform, "asinh")) {
        return(sinh(values) * 150)
    }
    if (identical(transform, "biexponential")) {
        return(sign(values) * 50 * (10^abs(values) - 1))
    }
    values
}

.reference_gui_extent <- function(values) {
    finite <- as.numeric(values[is.finite(values)])
    if (length(finite) == 0L) return(c(0, 1))
    val_max <- as.numeric(stats::quantile(finite, 0.998, na.rm = TRUE, names = FALSE))
    if (!is.finite(val_max) || val_max <= 0) val_max <- max(finite, na.rm = TRUE)
    if (!is.finite(val_max) || val_max <= 0) val_max <- 1
    c(0, val_max + val_max * 0.04)
}

.reference_gui_ticks <- function(domain, count = 5L) {
    domain <- as.numeric(domain[1:2])
    if (!all(is.finite(domain))) return(c(0, 1))
    if (domain[1] == domain[2]) return(domain[1])
    seq(domain[1], domain[2], length.out = count)
}

.reference_gui_domain_for_channel <- function(domains, channel, fallback_values = NULL) {
    if (!is.null(domains) && !is.null(channel) && nzchar(channel) && !is.null(domains[[channel]])) {
        domain <- as.numeric(domains[[channel]][1:2])
        if (length(domain) >= 2L && all(is.finite(domain)) && domain[2] > domain[1]) return(domain)
    }
    .reference_gui_extent(fallback_values)
}

.reference_gui_channel_extent <- function(raw_data, channel, filename, max_points = 100000L) {
    if (is.null(channel) || !nzchar(channel) || !channel %in% colnames(raw_data)) return(c(0, 1))
    idx <- .reference_gui_payload_indices(filename, nrow(raw_data), max_points = max_points)
    .reference_gui_extent(as.numeric(raw_data[idx, channel]))
}

.reference_gui_payload_indices <- function(filename, n, max_points = 100000L) {
    if (!is.finite(n) || n <= 0L) return(integer())
    max_points <- suppressWarnings(as.integer(max_points[1]))
    if (!is.finite(max_points) || is.na(max_points) || max_points <= 0L || n <= max_points) return(seq_len(n))
    seed <- sum(utf8ToInt(basename(filename))) %% .Machine$integer.max
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
        if (is.null(old_seed)) {
            if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
        } else {
            assign(".Random.seed", old_seed, envir = .GlobalEnv)
        }
    }, add = TRUE)
    set.seed(seed)
    sort(sample.int(n, max_points))
}

.reference_gui_display_data <- function(raw_data,
                                        filename,
                                        channels,
                                        display_max_points,
                                        preload_max_points = 100000L) {
    if (is.null(raw_data) || nrow(raw_data) == 0L) return(raw_data[0, , drop = FALSE])
    idx <- .reference_gui_payload_indices(filename, nrow(raw_data), max_points = preload_max_points)
    payload <- raw_data[idx, , drop = FALSE]
    channels <- unique(stats::na.omit(as.character(channels)))
    channels <- intersect(channels[nzchar(channels)], colnames(payload))
    if (length(channels) > 0L) {
        payload <- payload[stats::complete.cases(payload[, channels, drop = FALSE]), , drop = FALSE]
    }
    display_idx <- .reference_even_indices(nrow(payload), max_points = display_max_points)
    payload[display_idx, , drop = FALSE]
}

.reference_scatter_channels <- function(col_names) {
    unique(grep("^(FSC|SSC).*-(A|H|W)$", col_names, value = TRUE, ignore.case = TRUE))
}

.reference_compute_gui_scatter_domains <- function(fcs_files, max_points = 100000L) {
    domains <- list()
    for (fcs_file in fcs_files) {
        ff <- tryCatch(.spectreasy_read_fcs(fcs_file, label = "SCC FCS file"), error = function(e) NULL)
        if (is.null(ff)) next
        expr <- flowCore::exprs(ff)
        idx <- .reference_gui_payload_indices(basename(fcs_file), nrow(expr), max_points = max_points)
        channels <- .reference_scatter_channels(colnames(expr))
        for (channel in channels) {
            if (!channel %in% colnames(expr)) next
            values <- as.numeric(expr[idx, channel])
            domain <- .reference_gui_extent(values[is.finite(values)])
            if (is.null(domains[[channel]])) {
                domains[[channel]] <- domain
            } else {
                domains[[channel]] <- c(
                    min(domains[[channel]][1], domain[1], na.rm = TRUE),
                    max(domains[[channel]][2], domain[2], na.rm = TRUE)
                )
            }
        }
    }
    domains
}

.reference_histogram_density_curve <- function(values, domain, bins = 100L) {
    bins <- suppressWarnings(as.integer(bins))
    if (!is.finite(bins) || is.na(bins)) bins <- 100L
    bins <- min(max(bins, 5L), 500L)
    values <- values[is.finite(values)]
    if (length(values) == 0 || !all(is.finite(domain)) || domain[2] <= domain[1]) {
        return(data.frame(x = numeric(), y = numeric()))
    }
    idx <- pmax(0L, pmin(bins, floor(((values - domain[1]) / (domain[2] - domain[1])) * bins)))
    counts <- tabulate(idx + 1L, nbins = bins + 1L)
    smooth_radius <- max(1L, min(3L, floor(bins / 40L)))
    y <- numeric(bins + 1L)
    for (i in 0:bins) {
        js <- (-smooth_radius):smooth_radius
        bin_idx <- i + js
        ok <- bin_idx >= 0L & bin_idx <= bins
        weights <- exp(-(js[ok] * js[ok]) / 4)
        y[i + 1L] <- sum(counts[bin_idx[ok] + 1L] * weights) / max(sum(weights), 1)
    }
    data.frame(x = seq(domain[1], domain[2], length.out = bins + 1L), y = y)
}

.reference_qc_histogram_gui_plot <- function(peak_vals,
                                             pd,
                                             peak_channel,
                                             vals_log,
                                             gate_min,
                                             gate_max,
                                             settings,
                                             hist_info = NULL,
                                             x_domain = NULL) {
    transform <- settings$histogram_transform
    bins <- settings$histogram_bins
    values_t <- .reference_histogram_transform_values(peak_vals, transform = transform)
    pos_min <- if (!is.null(hist_info) && !is.null(hist_info$positive_raw_min)) hist_info$positive_raw_min else attr(vals_log, "pos_raw_min")
    pos_max <- if (!is.null(hist_info) && !is.null(hist_info$positive_raw_max)) hist_info$positive_raw_max else attr(vals_log, "pos_raw_max")
    scalar_finite <- function(x) length(x) == 1L && is.finite(x)
    if (!scalar_finite(pos_min)) pos_min <- gate_min
    if (!scalar_finite(pos_max)) pos_max <- gate_max
    neg_min <- if (!is.null(hist_info) && !is.null(hist_info$negative_raw_min)) hist_info$negative_raw_min else attr(vals_log, "neg_raw_min")
    neg_max <- if (!is.null(hist_info) && !is.null(hist_info$negative_raw_max)) hist_info$negative_raw_max else attr(vals_log, "neg_raw_max")
    neg_log_min <- attr(vals_log, "neg_log_min")
    neg_log_max <- attr(vals_log, "neg_log_max")
    if (!scalar_finite(neg_min) && scalar_finite(neg_log_min)) neg_min <- 10^neg_log_min
    if (!scalar_finite(neg_max) && scalar_finite(neg_log_max)) neg_max <- 10^neg_log_max
    raw_domain <- if (!is.null(x_domain) && length(x_domain) >= 2L && all(is.finite(x_domain[1:2])) && x_domain[2] > x_domain[1]) {
        as.numeric(x_domain[1:2])
    } else {
        .reference_gui_extent(peak_vals)
    }
    domain <- .reference_histogram_transform_values(raw_domain, transform = transform)
    if (!all(is.finite(domain)) || domain[2] <= domain[1]) domain <- c(0, 1)
    x_breaks <- .reference_gui_ticks(domain, 5L)
    x_labels <- .reference_pretty_k_label(.reference_histogram_inverse_values(x_breaks, transform = transform))
    curve <- .reference_histogram_density_curve(values_t, domain, bins = bins)
    ymax <- max(curve$y, 1, na.rm = TRUE)
    y_lim <- c(0, ymax * 1.12)
    coord_ratio <- .reference_gui_coord_ratio(domain, y_lim)
    x_desc <- .get_reference_axis_label(peak_channel, pd)
    p <- ggplot2::ggplot(curve, ggplot2::aes(x, y)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = y), fill = "#263f73", alpha = 0.28) +
        ggplot2::geom_line(color = "#263f73", linewidth = 0.85) +
        ggplot2::labs(title = "Histogram", x = x_desc, y = "Events per bin") +
        ggplot2::scale_x_continuous(breaks = x_breaks, labels = x_labels) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            legend.position = "none",
            panel.grid = ggplot2::element_blank(),
            panel.background = ggplot2::element_rect(fill = "#fbfaf6", color = "#ded9cf", linewidth = 0.35),
            axis.line = ggplot2::element_line(color = "#6c746f", linewidth = 0.35),
            axis.ticks = ggplot2::element_line(color = "#6c746f", linewidth = 0.3),
            axis.title = ggplot2::element_text(color = "#626a65", face = "bold"),
            axis.text = ggplot2::element_text(color = "#626a65", face = "bold"),
            plot.title = ggplot2::element_text(color = "#17201c", face = "plain", size = 14),
            plot.margin = ggplot2::margin(5.5, 24, 5.5, 5.5)
        ) +
        ggplot2::coord_fixed(ratio = coord_ratio, xlim = domain, ylim = y_lim, expand = FALSE, clip = "off")
    if (isTRUE(attr(vals_log, "negative_gate_present")) && scalar_finite(neg_min) && scalar_finite(neg_max) && neg_max > neg_min) {
        neg_t <- .reference_histogram_transform_values(c(neg_min, neg_max), transform = transform)
        neg_handles <- data.frame(x = neg_t, y = rep(ymax * 0.5, length(neg_t)))
        p <- p +
            ggplot2::annotate("rect", xmin = min(neg_t), xmax = max(neg_t), ymin = -Inf, ymax = Inf, alpha = 0.13, fill = "#263f73") +
            ggplot2::geom_vline(xintercept = neg_t, color = "#263f73", linewidth = 1) +
            ggplot2::geom_point(data = neg_handles, ggplot2::aes(x, y), inherit.aes = FALSE, shape = 21, size = 2.4, stroke = 0.9, color = "#263f73", fill = "#fbfaf6")
    }
    if (isTRUE(attr(vals_log, "positive_gate_present")) && scalar_finite(pos_min) && scalar_finite(pos_max) && pos_max > pos_min) {
        pos_t <- .reference_histogram_transform_values(c(pos_min, pos_max), transform = transform)
        pos_handles <- data.frame(x = pos_t, y = rep(ymax * 0.5, length(pos_t)))
        p <- p +
            ggplot2::annotate("rect", xmin = min(pos_t), xmax = max(pos_t), ymin = -Inf, ymax = Inf, alpha = 0.13, fill = "#d65238") +
            ggplot2::geom_vline(xintercept = pos_t, color = "#d65238", linewidth = 1) +
            ggplot2::geom_point(data = pos_handles, ggplot2::aes(x, y), inherit.aes = FALSE, shape = 21, size = 2.4, stroke = 0.9, color = "#d65238", fill = "#fbfaf6")
    }
    p
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
                                     manual_gate_info = NULL,
                                     hist_info = NULL,
                                     gui_scatter_domains = NULL,
                                     histogram_domain = NULL,
                                     source_filename = NULL,
                                     sample_type = NULL) {
    gate_settings <- if (!is.null(manual_gate_info) && !is.null(manual_gate_info$settings)) {
        manual_gate_info$settings
    } else {
        .reference_manual_gate_settings(NULL)
    }
    if (is.null(source_filename) || !nzchar(as.character(source_filename)[1])) {
        source_filename <- paste0(sn, ".fcs")
    }
    cell_info <- if (!is.null(manual_gate_info) && !is.null(manual_gate_info$cell)) manual_gate_info$cell else NULL
    singlet_info <- if (!is.null(manual_gate_info) && !is.null(manual_gate_info$singlet)) manual_gate_info$singlet else NULL
    cell_x <- if (!is.null(cell_info) && !is.null(cell_info$x_channel) && cell_info$x_channel %in% colnames(raw_data)) cell_info$x_channel else fsc
    cell_y <- if (!is.null(cell_info) && !is.null(cell_info$y_channel) && cell_info$y_channel %in% colnames(raw_data)) cell_info$y_channel else ssc
    cell_vertices <- if (!is.null(cell_info) && !is.null(cell_info$vertices)) cell_info$vertices else final_gate
    singlet_vertices <- if (!is.null(singlet_info) && !is.null(singlet_info$vertices)) singlet_info$vertices else NULL
    singlet_x <- if (!is.null(singlet_info) && !is.null(singlet_info$x_channel) && singlet_info$x_channel %in% colnames(raw_data)) singlet_info$x_channel else NA_character_
    singlet_y <- if (!is.null(singlet_info) && !is.null(singlet_info$y_channel) && singlet_info$y_channel %in% colnames(raw_data)) singlet_info$y_channel else NA_character_
    display_data <- .reference_gui_display_data(
        raw_data = raw_data,
        filename = source_filename,
        channels = c(cell_x, cell_y, singlet_x, singlet_y, peak_channel),
        display_max_points = gate_settings$max_points,
        preload_max_points = 100000L
    )
    cell_keep_count <- if (!is.null(cell_info) && !is.null(cell_info$keep_count) && is.finite(cell_info$keep_count)) {
        as.numeric(cell_info$keep_count)
    } else {
        nrow(gated_data)
    }
    p1 <- .reference_qc_scatter_density_plot(
        data = display_data,
        x_channel = cell_x,
        y_channel = cell_y,
        pd = pd,
        gate_vertices = cell_vertices,
        title = paste0(sn, " - Cell Gate"),
        subtitle = paste0(round(100 * cell_keep_count / nrow(raw_data), 1), "% gated"),
        max_points = nrow(display_data),
        point_size = gate_settings$point_size,
        x_domain = .reference_gui_domain_for_channel(gui_scatter_domains, cell_x, raw_data[, cell_x]),
        y_domain = .reference_gui_domain_for_channel(gui_scatter_domains, cell_y, raw_data[, cell_y]),
        sample_type = sample_type
    )
    .save_reference_ggsave(file.path(out_path, "fsc_ssc", paste0(sn, "_fsc_ssc.png")), p1, sn, "FSC/SSC", width = 5.2, height = 4.2, dpi = 300)
    fsc_desc <- .get_reference_axis_label(fsc, pd)

    cell_keep_display <- if (!is.null(cell_vertices) && nrow(cell_vertices) >= 3 && cell_x %in% colnames(display_data) && cell_y %in% colnames(display_data)) {
        sp::point.in.polygon(display_data[, cell_x], display_data[, cell_y], cell_vertices$x, cell_vertices$y) > 0
    } else {
        rep(TRUE, nrow(display_data))
    }
    histogram_source <- display_data[cell_keep_display, , drop = FALSE]
    if (!is.null(singlet_vertices) && nrow(singlet_vertices) >= 3 && !is.na(singlet_x) && !is.na(singlet_y)) {
        singlet_source <- histogram_source
        singlet_keep_count <- if (!is.null(singlet_info$keep_count) && is.finite(singlet_info$keep_count)) {
            as.numeric(singlet_info$keep_count)
        } else {
            nrow(gated_data)
        }
        singlet_source_count <- if (!is.null(singlet_info$source_count) && is.finite(singlet_info$source_count)) {
            as.numeric(singlet_info$source_count)
        } else {
            nrow(singlet_source)
        }
        p_singlet <- .reference_qc_scatter_density_plot(
            data = singlet_source,
            x_channel = singlet_x,
            y_channel = singlet_y,
            pd = pd,
            gate_vertices = singlet_vertices,
            title = paste0(sn, " - Singlet Gate"),
            subtitle = paste0(round(100 * singlet_keep_count / max(singlet_source_count, 1), 1), "% gated"),
            max_points = nrow(singlet_source),
            point_size = gate_settings$point_size,
            x_domain = .reference_gui_domain_for_channel(gui_scatter_domains, singlet_x, raw_data[, singlet_x]),
            y_domain = .reference_gui_domain_for_channel(gui_scatter_domains, singlet_y, raw_data[, singlet_y]),
            sample_type = sample_type
        )
        .save_reference_ggsave(file.path(out_path, "singlet", paste0(sn, "_singlet.png")), p_singlet, sn, "singlet", width = 5.2, height = 4.2, dpi = 300)
        if (singlet_x %in% colnames(histogram_source) && singlet_y %in% colnames(histogram_source)) {
            singlet_keep_display <- sp::point.in.polygon(histogram_source[, singlet_x], histogram_source[, singlet_y], singlet_vertices$x, singlet_vertices$y) > 0
            histogram_source <- histogram_source[singlet_keep_display, , drop = FALSE]
        } else {
            histogram_source <- histogram_source[0, , drop = FALSE]
        }
    }

    neg_log_min <- attr(vals_log, "neg_log_min")
    neg_log_max <- attr(vals_log, "neg_log_max")
    negative_gate_present <- isTRUE(attr(vals_log, "negative_gate_present"))
    positive_gate_present <- isTRUE(attr(vals_log, "positive_gate_present"))
    negative_gate_valid <- negative_gate_present &&
        is.finite(neg_log_min) &&
        is.finite(neg_log_max) &&
        neg_log_max > neg_log_min
    positive_gate_valid <- positive_gate_present &&
        is.finite(gate_min) &&
        is.finite(gate_max) &&
        gate_max > gate_min
    histogram_peak_vals <- if (peak_channel %in% colnames(histogram_source)) {
        histogram_source[, peak_channel]
    } else {
        numeric()
    }
    spectrum_plot_data <- if (positive_gate_valid && peak_channel %in% colnames(histogram_source)) {
        histogram_source[histogram_peak_vals >= gate_min & histogram_peak_vals <= gate_max, , drop = FALSE]
    } else {
        final_gated_data
    }
    p2 <- .reference_qc_histogram_gui_plot(
        peak_vals = histogram_peak_vals,
        pd = pd,
        peak_channel = peak_channel,
        vals_log = vals_log,
        gate_min = gate_min,
        gate_max = gate_max,
        settings = gate_settings,
        hist_info = hist_info,
        x_domain = histogram_domain
    )
    .save_reference_ggsave(file.path(out_path, "histogram", paste0(sn, "_histogram.png")), p2, sn, "histogram", width = 6.5, height = 4, dpi = 300)
    if (nrow(spectrum_plot_data) > 0L && all(detector_names %in% colnames(spectrum_plot_data))) {
        log_mat <- log10(pmax(spectrum_plot_data[, detector_names, drop = FALSE], 1e-3))
    } else {
        log_mat <- matrix(numeric(), nrow = 0L, ncol = length(detector_names), dimnames = list(NULL, detector_names))
    }
    finite_log <- log_mat[is.finite(log_mat)]
    if (length(finite_log) == 0L) {
        finite_log <- c(0, 1)
    }
    min_y <- floor(min(finite_log, na.rm = TRUE))
    max_y <- ceiling(max(finite_log, na.rm = TRUE))
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
        .spectreasy_console_step("QC plot", paste0(sn, " spectrum histogram is empty"))
    }

    vlines <- which(diff(det_info$laser_nm) != 0) + 0.5
    fill_lo <- min(dt_c$fill, na.rm = TRUE)
    fill_hi <- quantile(dt_c$fill, 0.96, na.rm = TRUE)
    if (!is.finite(fill_lo) || !is.finite(fill_hi) || fill_hi <= fill_lo) {
        fill_lo <- 0
        fill_hi <- 1
    }
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
    .save_reference_ggsave(file.path(out_path, "spectrum", paste0(sn, "_spectrum.png")), p3, sn, "spectrum", width = 300, height = 120, units = "mm", dpi = 600)

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
                                    bead_negative = NULL,
                                    scc_background = NULL) {
    sn_ext <- basename(fcs_file)
    sn <- tools::file_path_sans_ext(sn_ext)

    is_extra_af <- FALSE

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

    if (is_extra_af) {
        return(NULL)
    }
    if (.is_af_control_row(fluorophore = fluor_name, marker = marker_name, filename = sn_ext)) {
        return(NULL)
    }

    ff <- .spectreasy_read_fcs(fcs_file, label = "SCC FCS file")
    pd <- flowCore::pData(flowCore::parameters(ff))
    raw_data <- flowCore::exprs(ff)

    if (!.validate_reference_raw_data(raw_data, sn)) {
        return(NULL)
    }

    scatter_info <- .apply_reference_manual_scatter_gates(
        raw_data = raw_data,
        sample_type = sample_info$type,
        filename = sn_ext,
        manual_gates = config$manual_gates
    )
    if (is.null(scatter_info)) {
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
    }
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
        cytometer = config$cytometer
    )
    peak_channel <- peak_info$peak_channel
    .spectreasy_console_field("SCC", paste0(fluor_name, " (", sn, ") -> ", peak_channel))

    peak_vals <- scatter_info$gated_data[, peak_channel]
    if (isTRUE(.get_reference_config_value(config, "spectral_scc_pipeline", FALSE))) {
        positive_gate <- .resolve_reference_positive_histogram_gate(
            gated_data = scatter_info$gated_data,
            peak_channel = peak_channel,
            filename = sn_ext,
            sample_type = sample_info$type,
            manual_gates = config$manual_gates,
            config = config,
            row_info = row_info
        )
        if (is.null(positive_gate)) return(NULL)
        autospectral_clean_data <- positive_gate$final_gated_data
        internal_negative <- .reference_histogram_negative_source(
            gated_data = scatter_info$gated_data,
            peak_channel = peak_channel,
            filename = sn_ext,
            sample_type = sample_info$type,
            manual_gates = config$manual_gates,
            detector_names = metadata$detector_names,
            fsc = scatter_info$fsc,
            ssc = scatter_info$ssc,
            hist_info = positive_gate$hist_info
        )
        autospectral_peak_vals <- autospectral_clean_data[, peak_channel]
        use_background_subtraction <- isTRUE(.get_reference_config_value(config, "scc_background_enabled", FALSE))
        negative_source <- .resolve_reference_scc_negative_source(
            row_info = row_info,
            sample_type = sample_info$type,
            af_data_raw = af_data_raw,
            scc_background = scc_background,
            universal_negatives = universal_negatives,
            bead_negative = bead_negative
        )
        selected_background <- if (isTRUE(use_background_subtraction)) negative_source$background else NULL
        selected_negative <- if (isTRUE(use_background_subtraction)) negative_source$negative else NULL
        if (isTRUE(use_background_subtraction) && is.null(selected_negative) && !is.null(internal_negative)) {
            selected_negative <- internal_negative
            selected_background <- attr(internal_negative, "scc_background", exact = TRUE)
        }
        extraction <- .compute_reference_autospectral_scc(
            clean_data = autospectral_clean_data,
            detector_names = metadata$detector_names,
            peak_channel = peak_channel,
            sample_type = sample_info$type,
            af_data_raw = selected_negative,
            scc_background = selected_background,
            n_candidates = .get_reference_config_value(config, "autospectral_n_candidates", 1000L),
            n_spectral = .get_reference_config_value(config, "autospectral_n_spectral", 200L),
            min_events = .get_reference_config_value(config, "autospectral_min_events", 10L),
            scc_background_k = .get_reference_config_value(config, "scc_background_k", 2L)
        )
        if (is.null(extraction)) {
            return(NULL)
        }

        spectrum_norm <- extraction$spectrum
        final_gated_data <- extraction$final_gated_data
        hist_info <- list(
            vals_log = extraction$vals_log,
            gate_min = extraction$gate_min,
            gate_max = extraction$gate_max
        )
        qc_hist_info <- positive_gate$hist_info

        selected_peak <- autospectral_peak_vals[extraction$positive_idx]
        neg_vals <- autospectral_peak_vals[order(autospectral_peak_vals, decreasing = FALSE)[seq_len(max(1L, floor(0.10 * length(autospectral_peak_vals))))]]
        mfi_pos <- stats::median(selected_peak, na.rm = TRUE)
        mfi_neg <- stats::median(neg_vals, na.rm = TRUE)
        sd_neg <- stats::mad(neg_vals, na.rm = TRUE)
        if (is.na(sd_neg) || sd_neg == 0) sd_neg <- stats::sd(neg_vals, na.rm = TRUE)
        if (is.na(sd_neg) || sd_neg == 0) sd_neg <- 1e-6
        stain_index <- (mfi_pos - mfi_neg) / (2 * sd_neg)
        any_sat <- any(raw_data[, metadata$detector_names, drop = FALSE] >= 260000, na.rm = TRUE)

        if (isTRUE(config$save_qc_plots)) {
            .save_reference_qc_plots_safely(
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
                vals_log = qc_hist_info$vals_log,
                peak_vals = peak_vals,
                gate_min = qc_hist_info$gate_min,
                gate_max = qc_hist_info$gate_max,
                peak_channel = peak_channel,
                detector_names = metadata$detector_names,
                detector_labels = metadata$detector_labels,
                det_info = metadata$det_info,
                out_path = config$out_path,
                manual_gate_info = scatter_info$manual_gate_info,
                hist_info = qc_hist_info,
                gui_scatter_domains = config$gui_scatter_domains,
                histogram_domain = .reference_gui_channel_extent(raw_data, peak_channel, sn_ext, max_points = 100000L),
                source_filename = sn_ext,
                sample_type = sample_info$type
            )
        }

        return(list(
            sample_name = sn,
            result = data.table::data.table(
                sample = sn,
                fluorophore = fluor_name,
                type = sample_info$type,
                n_total = nrow(raw_data),
                n_final = extraction$n_selected,
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
                n_final = extraction$n_selected,
                scatter_gate_pct = round(100 * nrow(scatter_info$gated_data) / max(nrow(raw_data), 1), 1),
                histogram_gate_pct = round(100 * extraction$n_selected / max(nrow(scatter_info$gated_data), 1), 1),
                intensity_gate_type = extraction$extraction_method,
                scc_background_method = {
                    bg_info <- attr(spectrum_norm, "scc_background")
                    if (!is.null(bg_info) && !is.null(bg_info$method)) bg_info$method else "none"
                },
                stain_index = round(stain_index, 1),
                saturated = ifelse(any_sat, "YES", "OK")
            )
        ))
    }

    positive_gate <- .resolve_reference_positive_histogram_gate(
        gated_data = scatter_info$gated_data,
        peak_channel = peak_channel,
        filename = sn_ext,
        sample_type = sample_info$type,
        manual_gates = config$manual_gates,
        config = config,
        row_info = row_info
    )
    if (is.null(positive_gate)) return(NULL)
    hist_info <- positive_gate$hist_info
    final_gated_data <- positive_gate$final_gated_data

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
        bead_negative = bead_negative
    )

    if (isTRUE(config$save_qc_plots)) {
        .save_reference_qc_plots_safely(
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
            manual_gate_info = scatter_info$manual_gate_info,
            hist_info = hist_info,
            gui_scatter_domains = config$gui_scatter_domains,
            histogram_domain = .reference_gui_channel_extent(raw_data, peak_channel, sn_ext, max_points = 100000L),
            source_filename = sn_ext,
            sample_type = sample_info$type
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
                if (!is.null(gate_type) && nzchar(gate_type)) gate_type else "histogram"
            },
            stain_index = round(stain_index, 1),
            saturated = ifelse(any_sat, "YES", "OK")
        )
    )
}

# Combines individual sample spectra and AF signatures into a single spillover matrix.
# Extracts spectra, cleans up row/column names, structures the metadata attributes
# (such as QC summary and parameter info), and performs basic sanity checks.
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
#' @param af_profile Optional saved AF profile name, `spectreasy_af_profile`, or
#'   AF-only matrix. When supplied, its spectra replace extraction from mapped
#'   unstained cell controls.
#' @param af_n_bands Number of AF basis signatures to extract from pooled
#'   unstained/AF control events. The default, `100`, builds a broad fixed
#'   AF bank for Spectreasy unmixing.
#' @param af_max_cells Maximum number of scatter-gated AF events used when
#'   deriving AF basis signatures.
#' @param af_min_cluster_events Minimum number of AF events required to keep a
#'   k-means AF cluster. Used together with `af_min_cluster_proportion`; the
#'   larger threshold is applied.
#' @param af_min_cluster_proportion Minimum fraction of modeled scatter-gated AF
#'   events required to keep a k-means AF cluster. The default `0.005` means
#'   0.5\% of the AF events used for extraction.
#' @param seed Optional integer seed for deterministic subsampling/clustering.
#' @param n_threads Positive integer; number of threads used for event-wise
#'   AutoSpectral AF assignment during AF-bank refinement.
#' @param default_sample_type Fallback type when filename heuristics are ambiguous (`"beads"` or `"cells"`).
#' @param cytometer Cytometer name used as a channel-mapping hint. The default,
#'   `"auto"`, infers the cytometer from FCS detector names when possible.
#' @param manual_gate_file Optional gate CSV from [gate_controls()]. When
#'   provided, manual cell/singlet gates are used before automatic SCC spectrum
#'   extraction, and manual positive gates are used for standard SCC extraction.
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
#' @param unmixing_method SCC processing method. `"Spectreasy"` and
#'   `"AutoSpectral"` first apply the saved positive histogram gate or the same
#'   automatic histogram fallback used by legacy methods, then enable spectral
#'   SCC selection, resolved-negative background subtraction, and
#'   spectral-variant learning. Legacy methods stop after calculating the
#'   conventional histogram-gated reference spectrum.
#' @param scc_background_method Background subtraction method for Spectreasy/
#'   AutoSpectral SCC cleanup (`"scatter_knn"` or `"none"`).
#' @param scc_background_k Number of nearest unstained/negative events averaged
#'   for scatter-matched SCC background subtraction.
#' @param autospectral_n_candidates Number of peak-bright SCC candidate events
#'   considered by the AutoSpectral-style selector.
#' @param autospectral_n_spectral Number of least-background-like SCC events
#'   kept for spectrum calculation by the AutoSpectral-style selector.
#' @param autospectral_min_events Minimum event count required by the
#'   AutoSpectral-style SCC selector.
#' @param refine Logical; if `TRUE`, refine the fixed-size k-means AF bank with
#'   native AutoSpectral-style unstained residual modulation. This is only
#'   supported with `unmixing_method = "AutoSpectral"`.
#'
#' @return Numeric matrix with rows = fluorophores and columns = detectors
#'   (normalized spectra). The matrix carries SCC-derived detector noise floors
#'   in `attr(M, "detector_noise")` for WLS/RWLS unmixing.
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
  af_profile = NULL,
  af_n_bands = 100,
  af_max_cells = 50000,
  af_min_cluster_events = 20,
  af_min_cluster_proportion = 0.005,
  seed = NULL,
  default_sample_type = "beads",
  cytometer = "auto",
  manual_gate_file = NULL,
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
  subsample_n = 5000,
  unmixing_method = "Spectreasy",
  scc_background_method = c("scatter_knn", "none"),
  scc_background_k = 2L,
  autospectral_n_candidates = 1000L,
  autospectral_n_spectral = 200L,
  autospectral_min_events = 10L,
  refine = FALSE,
  n_threads = 1L
) {
    control_df <- .normalize_build_reference_control_df(control_df)
    unmixing_method <- .normalize_unmix_method(unmixing_method)
    n_threads <- .normalize_n_threads(n_threads)
    default_sample_type <- .match_arg_ci(default_sample_type, c("beads", "cells"), "default_sample_type")
    histogram_direction_beads <- .match_arg_ci(
        histogram_direction_beads, c("right", "both", "left"), "histogram_direction_beads"
    )
    histogram_direction_cells <- .match_arg_ci(
        histogram_direction_cells, c("right", "both", "left"), "histogram_direction_cells"
    )
    use_autospectral <- .is_autospectral_style_method(unmixing_method)
    refine <- .validate_reference_refine_arg(refine)
    if (isTRUE(refine) && !identical(unmixing_method, "AutoSpectral")) {
        stop("refine = TRUE is only supported with unmixing_method = \"AutoSpectral\".", call. = FALSE)
    }
    af_args <- .validate_build_reference_af_args(
        af_n_bands = af_n_bands,
        af_max_cells = af_max_cells,
        af_min_cluster_events = af_min_cluster_events,
        af_min_cluster_proportion = af_min_cluster_proportion
    )
    af_n_bands <- af_args$af_n_bands
    af_max_cells <- af_args$af_max_cells
    af_min_cluster_events <- af_args$af_min_cluster_events
    af_min_cluster_proportion <- af_args$af_min_cluster_proportion
    scc_background_args <- .validate_scc_background_args(
        scc_background_method = scc_background_method,
        scc_background_k = scc_background_k,
        enabled = use_autospectral
    )
    autospectral_n_candidates <- .normalize_positive_integer(autospectral_n_candidates, "autospectral_n_candidates")
    autospectral_n_spectral <- .normalize_positive_integer(autospectral_n_spectral, "autospectral_n_spectral")
    autospectral_min_events <- .normalize_positive_integer(autospectral_min_events, "autospectral_min_events")

    .with_optional_seed(seed)
    manual_gates <- .read_reference_manual_gates(manual_gate_file)

    sample_patterns <- get_fluorophore_patterns()
    file_info <- .prepare_reference_file_set(input_folder = input_folder, control_df = control_df)
    out_path <- .prepare_reference_output_path(output_folder = output_folder, save_qc_plots = save_qc_plots)
    metadata <- .prepare_reference_detector_info(file_info$fcs_files[1])
    cytometer <- .resolve_cytometer_from_pd(cytometer, metadata$pd_meta)
    gui_scatter_domains <- if (isTRUE(save_qc_plots)) {
        .reference_compute_gui_scatter_domains(file_info$fcs_files_all, max_points = 100000L)
    } else {
        NULL
    }

    config <- list(
        default_sample_type = default_sample_type,
        outlier_percentile = outlier_percentile,
        debris_percentile = debris_percentile,
        subsample_n = subsample_n,
        max_clusters = max_clusters,
        min_cluster_proportion = min_cluster_proportion,
        gate_contour_beads = gate_contour_beads,
        gate_contour_cells = gate_contour_cells,
        bead_gate_scale = bead_gate_scale,
        manual_gates = manual_gates,
        histogram_pct_beads = histogram_pct_beads,
        histogram_direction_beads = histogram_direction_beads,
        histogram_pct_cells = histogram_pct_cells,
        histogram_direction_cells = histogram_direction_cells,
        save_qc_plots = save_qc_plots,
        out_path = out_path,
        cytometer = cytometer,
        af_min_cluster_events = af_min_cluster_events,
        af_min_cluster_proportion = af_min_cluster_proportion,
        spectral_scc_pipeline = use_autospectral,
        scc_background_enabled = scc_background_args$enabled,
        scc_background_method = scc_background_args$method,
        scc_background_k = scc_background_args$k,
        autospectral_n_candidates = autospectral_n_candidates,
        autospectral_n_spectral = autospectral_n_spectral,
        autospectral_min_events = autospectral_min_events,
        refine = refine,
        gui_scatter_domains = gui_scatter_domains
    )

    .spectreasy_console_field("Detectors", paste0(length(metadata$detector_names), " spectral channel(s), sorted by laser"))
    .validate_reference_detector_consistency(
        fcs_files = file_info$fcs_files_all,
        detector_names = metadata$detector_names
    )

    af_profiles <- if (!is.null(af_profile)) {
        profile_name <- if (is.character(af_profile) && length(af_profile) == 1L) af_profile else "saved AF profile"
        profile_object <- if (is.character(af_profile) && length(af_profile) == 1L) {
            load_af_profile(af_profile, show_plot = FALSE)
        } else {
            af_profile
        }
        profile_matrix <- .coerce_af_profile_matrix(profile_object, arg_name = "af_profile")
        if (!setequal(colnames(profile_matrix), metadata$detector_names)) {
            missing_detectors <- setdiff(metadata$detector_names, colnames(profile_matrix))
            extra_detectors <- setdiff(colnames(profile_matrix), metadata$detector_names)
            stop(
                "Saved AF profile detectors do not match the SCC detector set.",
                if (length(missing_detectors) > 0L) paste0(" Missing: ", paste(missing_detectors, collapse = ", "), ".") else "",
                if (length(extra_detectors) > 0L) paste0(" Extra: ", paste(extra_detectors, collapse = ", "), ".") else "",
                call. = FALSE
            )
        }
        profile_matrix <- profile_matrix[, metadata$detector_names, drop = FALSE]
        profile_background <- if (.is_af_profile_object(profile_object)) profile_object$scc_background else NULL
        profile_raw_median <- if (.is_af_profile_object(profile_object)) profile_object$raw_median else NULL
        saved_profiles <- list(
            af_data_raw = profile_raw_median,
            af_signatures_norm = profile_matrix,
            af_bank_info = list(
                source_count = 0L,
                sources = data.frame(),
                pooled_events = 0L,
                requested_bands = nrow(profile_matrix),
                derived_bands = nrow(profile_matrix),
                mode = "saved_profile",
                profile_name = profile_name
            ),
            scc_background = profile_background,
            af_events = if (!is.null(profile_background$spectra)) profile_background$spectra else NULL
        )
        .spectreasy_console_field("AF bank", paste0(nrow(profile_matrix), " saved signature(s) from ", profile_name))
        saved_profiles
    } else {
        .collect_reference_af_profiles(
            control_df = control_df,
            fcs_files = file_info$fcs_files,
            fcs_files_all = file_info$fcs_files_all,
            detector_names = metadata$detector_names,
            af_n_bands = af_n_bands,
            af_max_cells = af_max_cells,
            af_min_cluster_events = af_min_cluster_events,
            af_min_cluster_proportion = af_min_cluster_proportion,
            config = config
        )
    }

    universal_negatives <- .collect_reference_universal_negatives(
        control_df = control_df,
        fcs_files = file_info$fcs_files_all,
        detector_names = metadata$detector_names,
        sample_patterns = sample_patterns,
        config = config
    )
    bead_negative <- .collect_reference_unstained_bead_negative(
        control_df = control_df,
        fcs_files = file_info$fcs_files_all,
        detector_names = metadata$detector_names,
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
            af_data_raw = af_profiles$af_data_raw,
            universal_negatives = universal_negatives,
            bead_negative = bead_negative,
            scc_background = if (isTRUE(config$spectral_scc_pipeline) && isTRUE(config$scc_background_enabled)) {
                af_profiles$scc_background
            } else {
                NULL
            }
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
        processed_results = results_list
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
    if (isTRUE(refine)) {
        refined_af <- .reference_refine_af_bank(
            M = M,
            af_events = af_profiles$af_events,
            af_n_bands = af_n_bands,
            af_max_cells = af_max_cells,
            af_min_cluster_events = af_min_cluster_events,
            af_min_cluster_proportion = af_min_cluster_proportion,
            n_threads = n_threads,
            seed = seed,
            verbose = TRUE
        )
        if (!is.null(refined_af) && !is.null(refined_af$signatures) && nrow(refined_af$signatures) > 0L) {
            af_bank_info <- attr(M, "af_bank_info")
            if (is.null(af_bank_info)) af_bank_info <- list()
            af_bank_info$requested_bands <- af_n_bands
            af_bank_info$derived_bands <- nrow(refined_af$signatures)
            af_bank_info$selection <- refined_af$selection
            af_bank_info$refine <- TRUE
            M <- .replace_reference_af_rows(
                M = M,
                af_signatures_norm = refined_af$signatures,
                af_bank_info = af_bank_info
            )
            .spectreasy_console_field("AF refine", paste0(nrow(refined_af$signatures), " fixed signature(s)"))
        } else {
            af_bank_info <- attr(M, "af_bank_info")
            if (!is.null(af_bank_info)) {
                af_bank_info$refine <- FALSE
                af_bank_info$refine_reason <- "no_refined_af_candidates"
                attr(M, "af_bank_info") <- af_bank_info
            }
            .spectreasy_console_field("AF refine", "kept base AF bank")
        }
    }
    .attach_estimated_wls_detector_noise(
        M = M,
        scc_dir = input_folder,
        fcs_files = file_info$fcs_files_all,
        fallback = .default_wls_background_noise(),
        warn = FALSE
    )
}
