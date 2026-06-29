# Internal helpers for SCC background subtraction.

.validate_scc_background_args <- function(clean_scc_with_unstained,
                                          scc_background_method,
                                          scc_background_k,
                                          require_for_af_cosine = FALSE) {
    method <- match.arg(scc_background_method, c("scatter_knn", "none"))
    k <- suppressWarnings(as.integer(scc_background_k[1]))
    if (!is.finite(k) || is.na(k) || k < 1L) {
        stop("scc_background_k must be an integer >= 1.", call. = FALSE)
    }
    enabled <- isTRUE(clean_scc_with_unstained) && method != "none"
    if (isTRUE(require_for_af_cosine) && !enabled) {
        warning(
            "use_af_cosine_scc_selection = TRUE requires scatter-matched ",
            "unstained SCC cleaning; enabling clean_scc_with_unstained = TRUE ",
            "with scc_background_method = 'scatter_knn'.",
            call. = FALSE
        )
        enabled <- TRUE
        method <- "scatter_knn"
    }
    list(
        enabled = enabled,
        method = method,
        k = k
    )
}

.scc_background_scale_scatter <- function(x, center = NULL, scale = NULL) {
    x <- as.matrix(x)
    if (is.null(center)) center <- apply(x, 2, stats::median, na.rm = TRUE)
    if (is.null(scale)) {
        scale <- apply(x, 2, stats::mad, na.rm = TRUE)
        bad <- !is.finite(scale) | scale <= 0
        if (any(bad)) {
            scale[bad] <- apply(x[, bad, drop = FALSE], 2, stats::sd, na.rm = TRUE)
        }
        bad <- !is.finite(scale) | scale <= 0
        if (any(bad)) scale[bad] <- 1
    }
    list(
        x = sweep(sweep(x, 2, center, "-"), 2, scale, "/"),
        center = center,
        scale = scale
    )
}

.scc_background_knn_indices <- function(query_scatter, reference_scatter, k, chunk_size = 500L) {
    k <- min(as.integer(k), nrow(reference_scatter))
    if (k < 1L || nrow(query_scatter) == 0 || nrow(reference_scatter) == 0) {
        return(matrix(integer(), nrow = nrow(query_scatter), ncol = 0))
    }

    if (requireNamespace("FNN", quietly = TRUE)) {
        get_knnx <- getExportedValue("FNN", "get.knnx")
        return(get_knnx(reference_scatter, query_scatter, k = k)$nn.index)
    }

    idx <- matrix(NA_integer_, nrow = nrow(query_scatter), ncol = k)
    starts <- seq.int(1L, nrow(query_scatter), by = chunk_size)
    for (start in starts) {
        stop <- min(start + chunk_size - 1L, nrow(query_scatter))
        for (i in seq.int(start, stop)) {
            d2 <- rowSums(sweep(reference_scatter, 2, query_scatter[i, ], "-")^2)
            idx[i, ] <- head(order(d2), k)
        }
    }
    idx
}

.scc_background_match <- function(events, background, k = 3L) {
    if (is.null(background) || is.null(background$scatter) || is.null(background$spectra)) {
        return(NULL)
    }
    scatter_names <- background$scatter_names
    detector_names <- background$detector_names
    if (!all(scatter_names %in% colnames(events)) || !all(detector_names %in% colnames(events))) {
        return(NULL)
    }

    complete_query <- stats::complete.cases(events[, c(scatter_names, detector_names), drop = FALSE])
    matched <- matrix(0, nrow = nrow(events), ncol = length(detector_names), dimnames = list(NULL, detector_names))
    if (!any(complete_query)) {
        return(matched)
    }

    ref_scaled <- .scc_background_scale_scatter(background$scatter)
    query_scaled <- .scc_background_scale_scatter(
        events[complete_query, scatter_names, drop = FALSE],
        center = ref_scaled$center,
        scale = ref_scaled$scale
    )$x

    idx <- .scc_background_knn_indices(
        query_scatter = query_scaled,
        reference_scatter = ref_scaled$x,
        k = k
    )
    if (ncol(idx) == 0) {
        return(NULL)
    }

    bg <- background$spectra
    matched[complete_query, ] <- t(apply(idx, 1, function(i) {
        colMeans(bg[i, , drop = FALSE], na.rm = TRUE)
    }))
    matched
}

.scc_background_clean_events <- function(events, detector_names, background, k = 3L) {
    matched <- .scc_background_match(events = events, background = background, k = k)
    if (is.null(matched)) {
        return(NULL)
    }
    clean <- as.matrix(events[, detector_names, drop = FALSE]) - matched
    colnames(clean) <- detector_names
    attr(clean, "matched_background") <- matched
    clean
}

.scc_background_clean_peak_matrix <- function(gated_data, detector_names, af_data_raw = NULL) {
    Y <- as.matrix(gated_data[, detector_names, drop = FALSE])
    if (is.null(af_data_raw) || !all(detector_names %in% names(af_data_raw))) {
        return(Y)
    }

    af_vec <- pmax(as.numeric(af_data_raw[detector_names]), 0)
    norm <- sqrt(sum(af_vec^2, na.rm = TRUE))
    if (!is.finite(norm) || norm <= 0) {
        return(Y)
    }

    af_unit <- af_vec / norm
    projection <- as.numeric(Y %*% af_unit)
    clean <- Y - projection %o% af_unit
    clean <- pmax(clean, 0)
    colnames(clean) <- detector_names
    clean
}

.scc_background_from_gated_af_list <- function(af_gated_list, detector_names) {
    if (!length(af_gated_list)) {
        return(NULL)
    }
    has_scatter <- vapply(
        af_gated_list,
        function(x) !is.null(x$scatter) && !is.null(x$events) && nrow(x$scatter) == nrow(x$events),
        logical(1)
    )
    af_gated_list <- af_gated_list[has_scatter]
    if (!length(af_gated_list)) {
        return(NULL)
    }

    scatter_names <- af_gated_list[[1]]$scatter_names
    scatter <- do.call(rbind, lapply(af_gated_list, `[[`, "scatter"))
    spectra <- do.call(rbind, lapply(af_gated_list, `[[`, "events"))
    keep <- stats::complete.cases(scatter, spectra)
    scatter <- scatter[keep, , drop = FALSE]
    spectra <- spectra[keep, detector_names, drop = FALSE]
    if (nrow(scatter) < 10L || nrow(spectra) < 10L) {
        return(NULL)
    }

    list(
        method = "scatter_knn",
        detector_names = detector_names,
        scatter_names = scatter_names,
        scatter = as.matrix(scatter),
        spectra = as.matrix(spectra),
        n_events = nrow(scatter)
    )
}

.scc_background_default_config <- function() {
    list(
        outlier_percentile = 0.02,
        debris_percentile = 0.08,
        subsample_n = 5000,
        max_clusters = 10,
        min_cluster_proportion = 0.03,
        gate_contour_beads = 0.95,
        gate_contour_cells = 0.90,
        bead_gate_scale = 1.3
    )
}

.collect_scc_background_from_controls <- function(scc_dir,
                                                  control_df,
                                                  detector_names,
                                                  k = 3L,
                                                  exclude_af = FALSE,
                                                  config = NULL) {
    if (isTRUE(exclude_af) || is.null(control_df) || !dir.exists(scc_dir)) {
        return(NULL)
    }
    fcs_files <- list.files(scc_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
    if (!length(fcs_files)) {
        return(NULL)
    }
    config <- utils::modifyList(.scc_background_default_config(), if (is.null(config)) list() else config)

    af_paths <- .resolve_reference_af_paths(control_df = control_df, fcs_files = fcs_files, exclude_af = exclude_af)$path
    if (!length(af_paths)) {
        return(NULL)
    }

    af_gated_list <- lapply(af_paths, .extract_reference_af_gated_events, detector_names = detector_names, config = config)
    keep <- vapply(af_gated_list, function(x) !is.null(x) && !is.null(x$events) && nrow(x$events) > 0, logical(1))
    background <- .scc_background_from_gated_af_list(af_gated_list[keep], detector_names = detector_names)
    if (!is.null(background)) {
        background$k <- as.integer(k)
    }
    background
}
