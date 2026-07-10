# Internal helpers for scatter-matched SCC background subtraction.

.validate_scc_background_args <- function(clean_scc_with_unstained,
                                          scc_background_method,
                                          scc_background_k) {
    method <- match.arg(scc_background_method, c("scatter_knn", "none"))
    clean_scc_with_unstained <- .normalize_scalar_logical(clean_scc_with_unstained, "clean_scc_with_unstained")
    k <- .normalize_positive_integer(scc_background_k, "scc_background_k")
    list(
        enabled = isTRUE(clean_scc_with_unstained) && method != "none",
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

.scc_background_match <- function(events, background, k = 2L) {
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

.scc_background_clean_events <- function(events, detector_names, background, k = 2L) {
    matched <- .scc_background_match(events = events, background = background, k = k)
    if (is.null(matched)) {
        return(NULL)
    }
    clean <- as.matrix(events[, detector_names, drop = FALSE]) - matched
    colnames(clean) <- detector_names
    attr(clean, "matched_background") <- matched
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
