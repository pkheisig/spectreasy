.spectral_variant_row_mask <- function(M) {
    !grepl("^AF($|_)", rownames(M), ignore.case = TRUE)
}

.spectral_variant_cosine <- function(A, B) {
    A <- as.matrix(A)
    B <- as.matrix(B)
    if (nrow(A) == 0 || nrow(B) == 0) {
        return(matrix(numeric(), nrow = nrow(A), ncol = nrow(B)))
    }
    denom <- sqrt(rowSums(A^2)) %o% sqrt(rowSums(B^2))
    out <- (A %*% t(B)) / denom
    out[!is.finite(out)] <- NA_real_
    out
}

.normalize_spectral_variant_shapes <- function(x) {
    x <- as.matrix(x)
    x[!is.finite(x)] <- 0
    x <- pmax(x, 0)
    scale <- apply(x, 1, max, na.rm = TRUE)
    keep <- is.finite(scale) & scale > 0
    if (!any(keep)) {
        return(x[FALSE, , drop = FALSE])
    }
    sweep(x[keep, , drop = FALSE], 1, scale[keep], "/")
}

.deduplicate_spectral_variant_shapes <- function(x, threshold = 0.995) {
    x <- as.matrix(x)
    if (nrow(x) <= 1L) return(x)
    keep <- integer()
    for (i in seq_len(nrow(x))) {
        if (!length(keep)) {
            keep <- i
            next
        }
        sim <- .spectral_variant_cosine(x[i, , drop = FALSE], x[keep, , drop = FALSE])
        if (!any(is.finite(sim) & sim >= threshold)) {
            keep <- c(keep, i)
        }
    }
    x[keep, , drop = FALSE]
}

.cluster_spectral_variant_shapes <- function(shape_mat, n_nodes = 16L, seed = NULL) {
    shape_mat <- .normalize_spectral_variant_shapes(shape_mat)
    if (nrow(shape_mat) == 0) return(shape_mat)

    n_nodes <- suppressWarnings(as.integer(n_nodes[1]))
    if (!is.finite(n_nodes) || is.na(n_nodes) || n_nodes < 1L) n_nodes <- 16L
    n_eff <- min(n_nodes, nrow(shape_mat))
    if (n_eff <= 1L) {
        out <- matrix(colMeans(shape_mat, na.rm = TRUE), nrow = 1)
        colnames(out) <- colnames(shape_mat)
        return(.normalize_spectral_variant_shapes(out))
    }

    rounded <- unique(as.data.frame(round(shape_mat, digits = 6)))
    n_eff <- min(n_eff, nrow(rounded))
    if (n_eff <= 1L) {
        out <- matrix(colMeans(shape_mat, na.rm = TRUE), nrow = 1)
        colnames(out) <- colnames(shape_mat)
        return(.normalize_spectral_variant_shapes(out))
    }

    if (!is.null(seed)) set.seed(seed)
    centers <- NULL
    if (requireNamespace("FlowSOM", quietly = TRUE)) {
        dims <- .reference_som_grid_dims(n_eff)
        centers <- tryCatch(
            {
                som <- FlowSOM::SOM(
                    shape_mat,
                    xdim = dims[["xdim"]],
                    ydim = dims[["ydim"]],
                    silent = TRUE
                )
                as.matrix(som$codes)
            },
            error = function(e) NULL
        )
    }
    if (is.null(centers)) {
        km <- stats::kmeans(shape_mat, centers = n_eff, nstart = 8, iter.max = 100)
        centers <- as.matrix(km$centers)
    }
    if (nrow(centers) > n_eff) centers <- centers[seq_len(n_eff), , drop = FALSE]
    colnames(centers) <- colnames(shape_mat)
    .normalize_spectral_variant_shapes(centers)
}

.finalize_spectral_variant_shapes <- function(shape_mat,
                                              fluorophore,
                                              M,
                                              event_count,
                                              som_nodes = 16L,
                                              cosine_threshold = 0.98,
                                              max_variants = 8L,
                                              min_events = 50L,
                                              seed = NULL) {
    if (is.null(shape_mat) || nrow(shape_mat) < min_events) {
        return(list(variants = NULL, info = list(reason = "insufficient_events", event_count = event_count)))
    }

    detectors <- colnames(M)
    base <- M[fluorophore, , drop = FALSE]
    centers <- .cluster_spectral_variant_shapes(shape_mat, n_nodes = som_nodes, seed = seed)
    sim_to_base <- as.numeric(.spectral_variant_cosine(centers, base))
    keep <- is.finite(sim_to_base) & sim_to_base >= cosine_threshold
    centers <- centers[keep, , drop = FALSE]
    sim_to_base <- sim_to_base[keep]
    if (nrow(centers) == 0) {
        return(list(variants = NULL, info = list(reason = "no_plausible_variants", event_count = event_count)))
    }

    delta_norm <- sqrt(rowSums((centers - matrix(base[1, ], nrow = nrow(centers), ncol = ncol(centers), byrow = TRUE))^2))
    keep_distinct <- is.finite(delta_norm) & delta_norm > 1e-4
    centers <- centers[keep_distinct, , drop = FALSE]
    sim_to_base <- sim_to_base[keep_distinct]
    delta_norm <- delta_norm[keep_distinct]
    if (nrow(centers) == 0) {
        return(list(variants = NULL, info = list(reason = "only_base_like_variants", event_count = event_count)))
    }

    centers <- centers[order(delta_norm, decreasing = TRUE), , drop = FALSE]
    centers <- .deduplicate_spectral_variant_shapes(centers, threshold = 0.995)
    max_variants <- suppressWarnings(as.integer(max_variants[1]))
    if (!is.finite(max_variants) || is.na(max_variants) || max_variants < 1L) max_variants <- 8L
    if (nrow(centers) > max_variants) {
        centers <- centers[seq_len(max_variants), , drop = FALSE]
    }
    rownames(centers) <- paste0(fluorophore, "_variant_", seq_len(nrow(centers)))
    colnames(centers) <- detectors

    list(
        variants = centers,
        info = list(
            reason = "ok",
            event_count = event_count,
            retained = nrow(centers),
            cosine_min = min(as.numeric(.spectral_variant_cosine(centers, base)), na.rm = TRUE),
            cosine_max = max(as.numeric(.spectral_variant_cosine(centers, base)), na.rm = TRUE)
        )
    )
}

.learn_one_spectral_variant_set_from_events <- function(events,
                                                        fluorophore,
                                                        M,
                                                        som_nodes = 16L,
                                                        cosine_threshold = 0.98,
                                                        max_variants = 8L,
                                                        min_events = 50L,
                                                        seed = NULL) {
    if (is.null(events) || !(fluorophore %in% rownames(M))) {
        return(NULL)
    }
    detectors <- colnames(M)
    events <- as.matrix(events)
    if (!all(detectors %in% colnames(events))) {
        return(list(variants = NULL, info = list(reason = "insufficient_events", event_count = 0L)))
    }
    events <- events[, detectors, drop = FALSE]
    events <- events[stats::complete.cases(events), , drop = FALSE]
    event_count <- nrow(events)
    shapes <- .normalize_spectral_variant_shapes(events)
    .finalize_spectral_variant_shapes(
        shape_mat = shapes,
        fluorophore = fluorophore,
        M = M,
        event_count = event_count,
        som_nodes = som_nodes,
        cosine_threshold = cosine_threshold,
        max_variants = max_variants,
        min_events = min_events,
        seed = seed
    )
}

.learn_spectral_variant_library <- function(scc_dir,
                                            control_df,
                                            M,
                                            enabled = TRUE,
                                            som_nodes = 16L,
                                            cosine_threshold = 0.98,
                                            max_variants = 8L,
                                            min_events = 50L,
                                            clean_scc_with_unstained = TRUE,
                                            scc_background_method = "scatter_knn",
                                            scc_background_k = 3L,
                                            exclude_af = FALSE,
                                            seed = NULL,
                                            warn = TRUE) {
    M <- .as_reference_matrix(M, "M")
    out <- list(
        enabled = isTRUE(enabled),
        detector_names = colnames(M),
        variants = list(),
        info = data.frame(),
        settings = list(
            som_nodes = as.integer(som_nodes[1]),
            cosine_threshold = as.numeric(cosine_threshold[1]),
            max_variants = as.integer(max_variants[1]),
            min_events = as.integer(min_events[1])
        )
    )
    class(out) <- c("spectreasy_spectral_variant_library", "list")
    if (!isTRUE(enabled)) return(out)
    fluorophores <- rownames(M)[.spectral_variant_row_mask(M)]
    reference_positive_events <- attr(M, "scc_positive_events")
    if (!is.list(reference_positive_events) || length(reference_positive_events) == 0L) {
        out$enabled <- FALSE
        if (isTRUE(warn)) {
            warning(
                "Spectral variants require SCC positive events stored on the reference matrix. ",
                "No `scc_positive_events` attribute was found, so variant optimization is disabled. ",
                "Re-run build_reference_matrix()/unmix_controls() with the current version.",
                call. = FALSE
            )
        }
        return(out)
    }
    bg_args <- .validate_scc_background_args(
        clean_scc_with_unstained = clean_scc_with_unstained,
        scc_background_method = scc_background_method,
        scc_background_k = scc_background_k
    )
    out$settings$clean_scc_with_unstained <- isTRUE(bg_args$enabled)
    out$settings$scc_background_method <- bg_args$method
    out$settings$scc_background_k <- bg_args$k

    rows <- list()
    for (fluor in fluorophores) {
        stored_events <- if (!is.null(reference_positive_events) && fluor %in% names(reference_positive_events)) {
            reference_positive_events[[fluor]]
        } else {
            NULL
        }
        learned <- if (!is.null(stored_events)) {
            .learn_one_spectral_variant_set_from_events(
                events = stored_events,
                fluorophore = fluor,
                M = M,
                som_nodes = som_nodes,
                cosine_threshold = cosine_threshold,
                max_variants = max_variants,
                min_events = min_events,
                seed = seed
            )
        } else {
            list(variants = NULL, info = list(reason = "missing_scc_positive_events", event_count = 0L))
        }
        if (is.null(learned)) {
            rows[[length(rows) + 1L]] <- data.frame(fluorophore = fluor, variants = 0L, reason = "missing_scc_positive_events", stringsAsFactors = FALSE)
            next
        }
        if (!is.null(learned$variants) && nrow(learned$variants) > 0) {
            out$variants[[fluor]] <- learned$variants
            rows[[length(rows) + 1L]] <- data.frame(
                fluorophore = fluor,
                variants = nrow(learned$variants),
                reason = learned$info$reason,
                event_count = learned$info$event_count,
                cosine_min = learned$info$cosine_min,
                cosine_max = learned$info$cosine_max,
                stringsAsFactors = FALSE
            )
        } else {
            rows[[length(rows) + 1L]] <- data.frame(
                fluorophore = fluor,
                variants = 0L,
                reason = learned$info$reason,
                event_count = learned$info$event_count,
                stringsAsFactors = FALSE
            )
        }
    }
    out$info <- if (length(rows)) data.table::rbindlist(rows, fill = TRUE) else data.frame()
    out
}

.save_spectral_variant_library <- function(library, path) {
    if (is.null(library)) return(invisible(NULL))
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(library, path, version = 3)
    invisible(path)
}

.load_spectral_variant_library <- function(x) {
    if (is.null(x) || (is.character(x) && (!length(x) || is.na(x[1]) || !nzchar(x[1])))) return(NULL)
    if (inherits(x, "spectreasy_spectral_variant_library")) return(x)
    if (is.character(x) && length(x) == 1L) {
        if (!file.exists(x)) return(NULL)
        obj <- tryCatch(readRDS(x), error = function(e) NULL)
        if (inherits(obj, "spectreasy_spectral_variant_library")) return(obj)
    }
    NULL
}

.spectral_variant_library_has_variants <- function(library) {
    inherits(library, "spectreasy_spectral_variant_library") &&
        length(library$variants) > 0 &&
        any(vapply(library$variants, function(x) is.matrix(x) && nrow(x) > 0, logical(1)))
}

.resolve_spectral_variant_library_for_unmixing <- function(M,
                                                           spectral_variant_library = NULL,
                                                           spectral_variant_library_file = NULL,
                                                           unmixing_matrix_file = NULL) {
    lib <- .load_spectral_variant_library(spectral_variant_library)
    if (!is.null(lib)) return(lib)
    attr_lib <- attr(M, "spectral_variant_library")
    lib <- .load_spectral_variant_library(attr_lib)
    if (!is.null(lib)) return(lib)
    lib <- .load_spectral_variant_library(spectral_variant_library_file)
    if (!is.null(lib)) return(lib)
    if (!is.null(unmixing_matrix_file) && file.exists(unmixing_matrix_file)) {
        sibling <- file.path(dirname(unmixing_matrix_file), "scc_spectral_variants.rds")
        lib <- .load_spectral_variant_library(sibling)
        if (!is.null(lib)) return(lib)
    }
    NULL
}

.spectral_variant_event_rows <- function(a_row, M) {
    non_af <- which(.spectral_variant_row_mask(M))
    af_idx <- which(!.spectral_variant_row_mask(M))
    rows <- non_af
    if (length(af_idx) == 1L) {
        rows <- c(rows, af_idx)
    } else if (length(af_idx) > 1L) {
        af_coeff <- abs(as.numeric(a_row[af_idx]))
        if (any(is.finite(af_coeff) & af_coeff > 0)) {
            rows <- c(rows, af_idx[[which.max(af_coeff)]])
        }
    }
    rows
}

.apply_spectral_variant_optimization <- function(Y,
                                                 M,
                                                 A,
                                                 method,
                                                 wls_noise,
                                                 rwls_max_iter,
                                                 spectral_variant_library,
                                                 top_k = 3L,
                                                 min_abundance = 1,
                                                 positive_fraction = 0.02,
                                                 min_improvement_fraction = 0.01) {
    if (!.spectral_variant_library_has_variants(spectral_variant_library)) {
        return(list(A = A, fitted = A %*% M, info = list(enabled = FALSE, changed_events = 0L)))
    }

    detectors <- colnames(M)
    if (!identical(detectors, spectral_variant_library$detector_names)) {
        warning("Spectral-variant library detectors do not match M; using the base reference matrix.", call. = FALSE)
        return(list(A = A, fitted = A %*% M, info = list(enabled = FALSE, changed_events = 0L, reason = "detector_mismatch")))
    }

    top_k <- max(1L, suppressWarnings(as.integer(top_k[1])))
    min_abundance <- suppressWarnings(as.numeric(min_abundance[1]))
    if (!is.finite(min_abundance) || is.na(min_abundance) || min_abundance < 0) min_abundance <- 1
    positive_fraction <- suppressWarnings(as.numeric(positive_fraction[1]))
    if (!is.finite(positive_fraction) || is.na(positive_fraction) || positive_fraction < 0) positive_fraction <- 0.02
    min_improvement_fraction <- suppressWarnings(as.numeric(min_improvement_fraction[1]))
    if (!is.finite(min_improvement_fraction) || is.na(min_improvement_fraction) || min_improvement_fraction < 0) {
        min_improvement_fraction <- 0.01
    }

    fitted <- A %*% M
    residuals <- Y - fitted
    A_out <- A
    fitted_out <- fitted
    changed <- logical(nrow(Y))
    selected <- vector("list", nrow(Y))
    fluor_rows <- rownames(M)[.spectral_variant_row_mask(M)]
    variant_names <- intersect(names(spectral_variant_library$variants), fluor_rows)
    if (!length(variant_names)) {
        return(list(A = A, fitted = fitted, info = list(enabled = FALSE, changed_events = 0L, reason = "no_matching_fluorophores")))
    }

    for (i in seq_len(nrow(Y))) {
        a_row <- A_out[i, ]
        marker_coeffs <- as.numeric(a_row[fluor_rows])
        names(marker_coeffs) <- fluor_rows
        finite_marker_coeffs <- marker_coeffs[is.finite(marker_coeffs)]
        if (!length(finite_marker_coeffs)) next
        max_signal <- max(abs(finite_marker_coeffs))
        if (!is.finite(max_signal) || max_signal < min_abundance) next
        threshold <- max(min_abundance, positive_fraction * max_signal)
        eligible <- names(marker_coeffs)[is.finite(marker_coeffs) & marker_coeffs > threshold]
        eligible <- intersect(eligible, variant_names)
        if (!length(eligible)) next

        base_r <- as.numeric(residuals[i, ])
        base_sse <- sum(base_r^2)
        if (!is.finite(base_sse) || base_sse <= 1e-12) next
        working_r <- base_r
        chosen <- list()

        for (fluor in eligible) {
            vars <- spectral_variant_library$variants[[fluor]]
            if (is.null(vars) || nrow(vars) == 0) next
            coeff <- as.numeric(a_row[[fluor]])
            if (!is.finite(coeff) || coeff <= threshold) next
            delta <- sweep(vars, 2, M[fluor, ], "-")
            predicted_sse <- rowSums((matrix(working_r, nrow = nrow(delta), ncol = ncol(delta), byrow = TRUE) - coeff * delta)^2)
            improvement <- sum(working_r^2) - predicted_sse
            top <- head(order(improvement, decreasing = TRUE), top_k)
            top <- top[is.finite(improvement[top]) & improvement[top] > max(1e-8, min_improvement_fraction * sum(working_r^2))]
            if (!length(top)) next
            best <- top[[which.max(improvement[top])]]
            chosen[[fluor]] <- vars[best, , drop = FALSE]
            working_r <- working_r - coeff * delta[best, ]
        }

        if (!length(chosen)) next
        M_event <- M
        for (fluor in names(chosen)) {
            M_event[fluor, ] <- chosen[[fluor]][1, ]
        }
        rows <- .spectral_variant_event_rows(a_row, M_event)
        X <- M_event[rows, , drop = FALSE]
        coeff_new <- tryCatch(
            .unmix_with_fixed_reference(
                Y = Y[i, , drop = FALSE],
                M = X,
                method = method,
                wls_noise = wls_noise,
                rwls_max_iter = rwls_max_iter
            ),
            error = function(e) NULL
        )
        if (is.null(coeff_new) || any(!is.finite(coeff_new))) next
        fitted_new <- coeff_new %*% X
        r_new <- as.numeric(Y[i, ] - fitted_new)
        sse_new <- sum(r_new^2)
        if (!is.finite(sse_new) || sse_new > base_sse * (1 - min_improvement_fraction)) next
        A_out[i, ] <- 0
        A_out[i, rows] <- coeff_new[1, ]
        fitted_out[i, ] <- fitted_new[1, ]
        changed[[i]] <- TRUE
        selected[[i]] <- names(chosen)
    }

    list(
        A = A_out,
        fitted = fitted_out,
        info = list(
            enabled = TRUE,
            changed_events = sum(changed),
            changed_fraction = mean(changed),
            selected = selected
        )
    )
}
