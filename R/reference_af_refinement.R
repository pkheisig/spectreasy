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

    universal_negative <- trimws(as.character(control_df$universal.negative))
    universal_negative[is.na(universal_negative)] <- ""
    empty_negative <- !nzchar(universal_negative)
    legacy_af_negative <- toupper(universal_negative) == "AF"
    source_rows <- primary_af_rows | dead_rows | bead_af_rows
    ordinary_cells <- control_type != "beads" & !viability_rows & !source_rows &
        (empty_negative | legacy_af_negative)
    viability_targets <- control_type != "beads" & viability_rows & !source_rows & empty_negative
    bead_targets <- control_type == "beads" & !source_rows & empty_negative
    if (length(af_files) > 0L) {
        control_df$universal.negative[ordinary_cells] <- af_files[1]
    }
    if (length(dead_files) > 0L) control_df$universal.negative[viability_targets] <- dead_files[1]
    if (length(bead_files) > 0L) control_df$universal.negative[bead_targets] <- bead_files[1]
    control_df
}

# Validates the Autofluorescence (AF) modeling parameters.
# Ensures that the number of AF bands and maximum events per file are positive
# integers, throwing an error if they are invalid.
# Returns a validated list containing these parameters.
.validate_build_reference_af_args <- function(af_n_bands,
                                              af_max_cells) {
    af_n_bands <- .normalize_positive_integer(af_n_bands, "af_n_bands")

    af_max_cells <- .normalize_positive_integer(af_max_cells, "af_max_cells")
    if (af_max_cells < 100L) {
        stop("af_max_cells must be an integer >= 100.")
    }

    list(
        af_n_bands = af_n_bands,
        af_max_cells = af_max_cells
    )
}

.validate_reference_refine_arg <- function(refine) {
    if (!is.logical(refine) || length(refine) != 1L || is.na(refine)) {
        stop("refine must be TRUE or FALSE.", call. = FALSE)
    }
    isTRUE(refine)
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

.prepare_reference_af_refinement <- function(M, af_events, af_max_cells, contaminant_threshold) {
    M <- .as_reference_matrix(M, "M")
    af_rows <- .reference_af_row_mask(M)
    if (!any(af_rows) || !any(!af_rows) || is.null(af_events)) return(NULL)
    detector_names <- colnames(M)
    af_events <- as.matrix(af_events)
    if (!all(detector_names %in% colnames(af_events))) return(NULL)
    af_events <- af_events[, detector_names, drop = FALSE]
    af_events <- af_events[stats::complete.cases(af_events), , drop = FALSE]
    if (nrow(af_events) < 100L) return(NULL)
    if (nrow(af_events) > af_max_cells) {
        af_events <- af_events[sample.int(nrow(af_events), af_max_cells), , drop = FALSE]
    }
    markers <- M[!af_rows, , drop = FALSE]
    af_events <- .reference_filter_af_contaminant_events(af_events, markers, contaminant_threshold)
    if (nrow(af_events) < 100L) return(NULL)
    list(
        M = M,
        detector_names = detector_names,
        events = af_events,
        marker_spectra = markers,
        base_af = M[af_rows, , drop = FALSE]
    )
}

.fit_reference_af_events <- function(events, marker_spectra, base_af, assignments, detector_names) {
    marker_idx <- seq_len(nrow(marker_spectra)) + 1L
    markers <- .unmix_with_fixed_reference(
        Y = events,
        M = marker_spectra,
        method = "OLS",
        wls_noise = list(noise_floor = numeric(0), signal_scale = numeric(0), max_weight_ratio = .default_wls_max_weight_ratio()),
        rwls_max_iter = 1L
    )
    unmixed <- cbind(AF = rep(0, nrow(events)), markers)
    residuals <- matrix(0, nrow = nrow(events), ncol = length(detector_names), dimnames = list(NULL, detector_names))
    projected <- residuals
    for (af_i in unique(assignments)) {
        event_idx <- which(assignments == af_i)
        X <- rbind(AF = base_af[af_i, ], marker_spectra)
        coefficients <- .unmix_with_fixed_reference(
            Y = events[event_idx, , drop = FALSE],
            M = X,
            method = "OLS",
            wls_noise = list(noise_floor = numeric(0), signal_scale = numeric(0), max_weight_ratio = .default_wls_max_weight_ratio()),
            rwls_max_iter = 1L
        )
        unmixed[event_idx, ] <- coefficients
        residuals[event_idx, ] <- events[event_idx, , drop = FALSE] - coefficients %*% X
        projected[event_idx, ] <- coefficients[, marker_idx, drop = FALSE] %*% marker_spectra
    }
    list(unmixed = unmixed, error = residuals + projected, marker_idx = marker_idx)
}

.reference_af_problem_events <- function(marker_abundance, problem_quantile) {
    magnitude <- if (ncol(marker_abundance) > 1L) sqrt(rowSums(marker_abundance^2)) else abs(marker_abundance[, 1])
    magnitude[!is.finite(magnitude)] <- 0
    quantile <- as.numeric(problem_quantile[1])
    if (!is.finite(quantile) || is.na(quantile) || quantile <= 0 || quantile >= 1) quantile <- 0.99
    repeat {
        threshold <- stats::quantile(magnitude, quantile, na.rm = TRUE, names = FALSE)
        indices <- which(magnitude > threshold)
        if (length(indices) >= 500L || quantile < 0.5) break
        quantile <- quantile - 0.05
    }
    list(indices = indices, quantile = quantile, threshold = threshold)
}

.reference_modulated_af_candidates <- function(base_af, assignments, problem_idx, spill_ratios, clusters) {
    modulated <- list()
    for (cluster in sort(unique(clusters$cluster))) {
        local_idx <- which(clusters$cluster == cluster)
        global_idx <- problem_idx[local_idx]
        ratio <- apply(spill_ratios[local_idx, , drop = FALSE], 2, stats::median, na.rm = TRUE)
        ratio[!is.finite(ratio)] <- 0
        for (af_i in unique(assignments[global_idx])) {
            updated <- pmax(base_af[af_i, ] * (1 + ratio), 0)
            peak <- max(updated, na.rm = TRUE)
            if (is.finite(peak) && peak > 1e-12) modulated[[length(modulated) + 1L]] <- updated / peak
        }
    }
    modulated
}

.reference_refine_af_bank <- function(M,
                                      af_events,
                                      af_n_bands,
                                      af_max_cells,
                                      n_threads = 1L,
                                      seed = NULL,
                                      problem_quantile = 0.99,
                                      contaminant_threshold = 0.99,
                                      verbose = TRUE) {
    prepared <- .prepare_reference_af_refinement(
        M, af_events, af_max_cells, contaminant_threshold
    )
    if (is.null(prepared)) return(NULL)
    M <- prepared$M
    detector_names <- prepared$detector_names
    af_events <- prepared$events
    marker_spectra <- prepared$marker_spectra
    base_af <- prepared$base_af

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

    fit <- .fit_reference_af_events(
        af_events, marker_spectra, base_af, assignments, detector_names
    )
    problem <- .reference_af_problem_events(
        fit$unmixed[, fit$marker_idx, drop = FALSE],
        problem_quantile
    )
    problem_idx <- problem$indices
    q <- problem$quantile
    threshold <- problem$threshold

    if (length(problem_idx) <= 10L) {
        return(NULL)
    }

    af_abundance <- fit$unmixed[problem_idx, 1]
    af_abundance[!is.finite(af_abundance) | af_abundance == 0] <- 1e-6
    spill_ratios <- sweep(fit$error[problem_idx, , drop = FALSE], 1, af_abundance, "/")
    clusters <- .reference_cluster_error_ratios(
        spill_ratios = spill_ratios,
        n_centers = af_n_bands,
        seed = seed
    )

    modulated <- .reference_modulated_af_candidates(
        base_af, assignments, problem_idx, spill_ratios, clusters
    )
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
        n_centers = af_n_bands
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
                                         nstart = 10,
                                         iter.max = 100) {
    af_shape <- as.matrix(af_shape)
    if (nrow(af_shape) == 0) {
        return(list(centers = af_shape, cluster_sizes = integer(), requested_centers = as.integer(n_centers)))
    }
    n_eff <- suppressWarnings(as.integer(n_centers[1]))
    if (!is.finite(n_eff) || is.na(n_eff) || n_eff < 1L) {
        stop("n_centers must be an integer >= 1.", call. = FALSE)
    }
    unique_count <- nrow(unique(as.data.frame(round(af_shape, digits = 8))))
    if (nrow(af_shape) < n_eff || unique_count < n_eff) {
        stop(
            "Cannot derive exactly ", n_eff, " AF bands from ", nrow(af_shape),
            " event(s) containing ", unique_count, " distinct spectral shape(s). ",
            "Reduce af_n_bands or provide more diverse AF events.",
            call. = FALSE
        )
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

    if (nrow(centers) != n_eff || length(cluster_sizes) != n_eff || any(cluster_sizes < 1L)) {
        stop(
            "AF clustering failed to return exactly ", n_eff,
            " non-empty bands. Try a different seed or provide more diverse AF events.",
            call. = FALSE
        )
    }

    ord <- order(cluster_sizes, decreasing = TRUE)
    centers <- centers[ord, , drop = FALSE]
    cluster_sizes <- cluster_sizes[ord]

    list(
        centers = centers,
        cluster_sizes = cluster_sizes,
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
