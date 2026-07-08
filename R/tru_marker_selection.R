.normalize_tru_gamma_mode <- function(tru_gamma_mode) {
    match.arg(tru_gamma_mode, choices = c("global", "event_specific"))
}

.normalize_tru_display <- function(tru_display) {
    match.arg(tru_display, choices = "ucm")
}

.normalize_tru_cutoff_quantile <- function(tru_cutoff_quantile, default = 0.995) {
    q <- if (is.null(tru_cutoff_quantile)) default else suppressWarnings(as.numeric(tru_cutoff_quantile[1]))
    if (!is.finite(q) || is.na(q) || q < 0 || q > 1) {
        stop("tru_cutoff_quantile must be a numeric value between 0 and 1.", call. = FALSE)
    }
    q
}

.tru_marker_rows <- function(M) {
    !grepl("^AF($|_)", rownames(M), ignore.case = TRUE)
}

.tru_ols <- function(Y, X, tol = 1e-10, ridge = 1e-8) {
    Y <- as.matrix(Y)
    X <- as.matrix(X)
    if (nrow(X) == 0L) {
        return(matrix(0, nrow = nrow(Y), ncol = 0L))
    }
    gram <- X %*% t(X)
    inv <- tryCatch(
        {
            if (rcond(gram) < tol) stop("ill-conditioned", call. = FALSE)
            solve(gram)
        },
        error = function(e) {
            scale <- max(c(diag(gram), 1), na.rm = TRUE)
            tryCatch(solve(gram + diag(ridge * scale, nrow(gram))), error = function(e2) NULL)
        }
    )
    if (is.null(inv)) {
        return(NULL)
    }
    Y %*% t(X) %*% inv
}

.autospectral_adjusted_marker_abundances <- function(Y, marker_spectra, af_spectra, tol = 1e-10) {
    Y <- as.matrix(Y)
    marker_spectra <- as.matrix(marker_spectra)
    af_spectra <- as.matrix(af_spectra)
    if (nrow(marker_spectra) == 0L || nrow(af_spectra) == 0L) {
        return(NULL)
    }
    if (!identical(colnames(marker_spectra), colnames(af_spectra))) {
        stop("AutoSpectral marker and AF spectra detector columns do not match.", call. = FALSE)
    }
    matrix_to_invert <- marker_spectra %*% t(marker_spectra)
    if (rcond(matrix_to_invert) < tol) {
        stop("Marker-only reference matrix is singular; TRU AF assignment cannot proceed.", call. = FALSE)
    }
    unmixing_matrix <- solve(matrix_to_invert, marker_spectra)
    marker_t <- t(marker_spectra)
    v_library <- unmixing_matrix %*% t(af_spectra)
    r_library <- t(af_spectra) - marker_t %*% v_library
    denominator <- colSums(r_library^2)
    valid <- is.finite(denominator) & denominator > tol
    if (!any(valid)) {
        adjusted <- array(0, dim = c(nrow(Y), nrow(marker_spectra), nrow(af_spectra)))
        dimnames(adjusted) <- list(NULL, rownames(marker_spectra), rownames(af_spectra))
        return(list(adjusted = adjusted, current_leakage = matrix(Inf, nrow(Y), nrow(af_spectra))))
    }

    numerator <- Y %*% r_library
    k_matrix <- matrix(0, nrow = nrow(Y), ncol = nrow(af_spectra))
    k_matrix[, valid] <- sweep(numerator[, valid, drop = FALSE], 2, denominator[valid], "/")
    unmixed_markers <- Y %*% t(unmixing_matrix)
    adjusted <- array(0, dim = c(nrow(Y), nrow(marker_spectra), nrow(af_spectra)))
    dimnames(adjusted) <- list(NULL, rownames(marker_spectra), rownames(af_spectra))
    current_leakage <- matrix(Inf, nrow = nrow(Y), ncol = nrow(af_spectra))
    for (i in which(valid)) {
        adj <- unmixed_markers - k_matrix[, i, drop = FALSE] %*% t(v_library[, i, drop = FALSE])
        adjusted[, , i] <- adj
        current_leakage[, i] <- rowSums(abs(adj))
    }
    list(adjusted = adjusted, current_leakage = current_leakage)
}

.autospectral_assign_af_fluorophores_tru <- function(Y,
                                                     marker_spectra,
                                                     af_spectra,
                                                     marker_cutoffs,
                                                     score_mode = c("thresholded_excess", "current"),
                                                     score_power = 2,
                                                     positive_only = TRUE,
                                                     tiebreak = c("current_leakage", "none"),
                                                     active_margin = 1.0,
                                                     tol = 1e-10,
                                                     return_details = FALSE) {
    score_mode <- match.arg(score_mode)
    tiebreak <- match.arg(tiebreak)
    details <- .autospectral_adjusted_marker_abundances(Y, marker_spectra, af_spectra, tol = tol)
    if (is.null(details)) {
        selected <- rep(NA_integer_, nrow(as.matrix(Y)))
        return(if (isTRUE(return_details)) list(selected_af = selected) else selected)
    }
    cutoffs <- as.numeric(marker_cutoffs[rownames(marker_spectra)])
    cutoffs[!is.finite(cutoffs)] <- 0
    margin_cutoffs <- active_margin * pmax(cutoffs, 0)
    score_power <- suppressWarnings(as.numeric(score_power[1]))
    if (!is.finite(score_power) || score_power <= 0) score_power <- 2

    n_events <- dim(details$adjusted)[1]
    n_af <- dim(details$adjusted)[3]
    score <- matrix(Inf, n_events, n_af)
    if (identical(score_mode, "current")) {
        score <- details$current_leakage
    } else {
        for (k in seq_len(n_af)) {
            adj <- matrix(
                details$adjusted[, , k, drop = FALSE],
                nrow = n_events,
                ncol = length(margin_cutoffs),
                dimnames = list(NULL, rownames(marker_spectra))
            )
            excess <- if (isTRUE(positive_only)) {
                sweep(adj, 2, margin_cutoffs, "-")
            } else {
                sweep(abs(adj), 2, margin_cutoffs, "-")
            }
            excess <- pmax(excess, 0)
            score[, k] <- rowSums(excess^score_power)
        }
    }

    selected <- max.col(-score, ties.method = "first")
    if (identical(tiebreak, "current_leakage") && n_af > 1L) {
        for (i in seq_len(n_events)) {
            best_score <- score[i, selected[[i]]]
            tied <- which(abs(score[i, ] - best_score) <= tol)
            if (length(tied) > 1L) {
                selected[[i]] <- tied[[which.min(details$current_leakage[i, tied])]]
            }
        }
    }

    if (!isTRUE(return_details)) {
        return(selected)
    }
    selected_adj <- matrix(0, nrow = n_events, ncol = nrow(marker_spectra), dimnames = list(NULL, rownames(marker_spectra)))
    for (i in seq_len(n_events)) {
        selected_adj[i, ] <- details$adjusted[i, , selected[[i]]]
    }
    list(
        selected_af = selected,
        selected_af_score = score[cbind(seq_len(n_events), selected)],
        current_leakage_score = details$current_leakage[cbind(seq_len(n_events), selected)],
        adjusted_marker_abundances_for_selected_af = selected_adj,
        score_matrix = score
    )
}

.recompute_tru_cutoffs <- function(tru_calibration, cutoff_quantile = NULL) {
    .validate_tru_calibration(tru_calibration)
    q <- .normalize_tru_cutoff_quantile(cutoff_quantile, default = tru_calibration$cutoff_quantile)
    null_abundances <- tru_calibration$null_abundances
    if (!is.null(null_abundances)) {
        null_abundances <- as.matrix(null_abundances)
        cutoffs <- apply(null_abundances, 2, function(x) max(0, stats::quantile(x, probs = q, na.rm = TRUE, names = FALSE, type = 7)))
        return(stats::setNames(as.numeric(cutoffs), colnames(null_abundances)))
    }
    if (!is.null(cutoff_quantile) && !identical(q, tru_calibration$cutoff_quantile)) {
        warning("TRU null distributions were not stored; using saved TRU cutoffs.", call. = FALSE)
    }
    tru_calibration$cutoffs
}

.build_tru_calibration <- function(M,
                                   unstained_events,
                                   cutoff_quantile = 0.995,
                                   min_null_events = 500L,
                                   max_events = 50000L,
                                   ns_correction = TRUE,
                                   store_null_distributions = TRUE) {
    M <- .as_reference_matrix(M, "M")
    marker_rows <- .tru_marker_rows(M)
    af_rows <- !marker_rows
    if (!any(marker_rows) || !any(af_rows) || is.null(unstained_events)) {
        return(NULL)
    }
    U <- as.matrix(unstained_events)
    detectors <- colnames(M)
    if (!all(detectors %in% colnames(U))) {
        return(NULL)
    }
    U <- U[, detectors, drop = FALSE]
    U <- U[stats::complete.cases(U), , drop = FALSE]
    min_null_events <- suppressWarnings(as.integer(min_null_events[1]))
    if (!is.finite(min_null_events) || min_null_events < 1L) min_null_events <- 500L
    if (nrow(U) < min_null_events) {
        return(NULL)
    }
    max_events <- suppressWarnings(as.integer(max_events[1]))
    if (!is.finite(max_events) || max_events < 1L) max_events <- 50000L
    if (nrow(U) > max_events) {
        U <- U[sample.int(nrow(U), max_events), , drop = FALSE]
    }

    null_fit <- tryCatch(
        .autospectral_unmix_legacy(Y = U, M = M),
        error = function(e) NULL
    )
    if (is.null(null_fit) || is.null(null_fit$A)) {
        return(NULL)
    }
    A_null <- null_fit$A[, rownames(M)[marker_rows], drop = FALSE]
    cutoff_quantile <- .normalize_tru_cutoff_quantile(cutoff_quantile)
    cutoffs <- apply(A_null, 2, function(x) max(0, stats::quantile(x, probs = cutoff_quantile, na.rm = TRUE, names = FALSE, type = 7)))
    null_center <- apply(A_null, 2, stats::median, na.rm = TRUE)
    null_mad <- apply(A_null, 2, stats::mad, na.rm = TRUE)
    q_probs <- c(q950 = 0.95, q975 = 0.975, q990 = 0.99, q995 = 0.995, q999 = 0.999)
    quant_mat <- t(vapply(seq_len(ncol(A_null)), function(j) {
        stats::quantile(A_null[, j], probs = q_probs, na.rm = TRUE, names = FALSE, type = 7)
    }, numeric(length(q_probs))))
    null_quantiles <- data.frame(marker = colnames(A_null), quant_mat, row.names = NULL, check.names = FALSE)
    o_ns <- as.numeric(colMeans(A_null, na.rm = TRUE) %*% M[marker_rows, , drop = FALSE])
    names(o_ns) <- detectors

    out <- list(
        version = 1L,
        detector_names = detectors,
        marker_names = rownames(M)[marker_rows],
        af_names = rownames(M)[af_rows],
        cutoff_quantile = cutoff_quantile,
        cutoffs = stats::setNames(as.numeric(cutoffs), colnames(A_null)),
        null_center = stats::setNames(as.numeric(null_center), colnames(A_null)),
        null_mad = stats::setNames(as.numeric(null_mad), colnames(A_null)),
        null_quantiles = null_quantiles,
        o_ns = o_ns,
        ns_correction = isTRUE(ns_correction),
        null_abundances = if (isTRUE(store_null_distributions)) A_null else NULL,
        n_null_events = as.integer(nrow(U)),
        created_from = list(
            method = "Spectreasy/TRU calibration",
            source = "unstained control",
            date = Sys.time()
        )
    )
    class(out) <- c("spectreasy_tru_calibration", "list")
    out
}

.validate_tru_calibration <- function(tru_calibration) {
    if (!inherits(tru_calibration, "spectreasy_tru_calibration")) {
        stop("Invalid TRU calibration object.", call. = FALSE)
    }
    required <- c("detector_names", "marker_names", "cutoffs")
    missing <- required[!vapply(required, function(nm) !is.null(tru_calibration[[nm]]), logical(1))]
    if (length(missing)) {
        stop("Invalid TRU calibration object; missing: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    invisible(TRUE)
}

.save_tru_calibration <- function(tru_calibration, path) {
    if (is.null(tru_calibration)) return(invisible(NULL))
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(tru_calibration, path, version = 3)
    invisible(path)
}

.write_tru_cutoffs_csv <- function(tru_calibration, path) {
    if (is.null(tru_calibration)) return(invisible(NULL))
    df <- data.frame(
        marker = names(tru_calibration$cutoffs),
        cutoff = as.numeric(tru_calibration$cutoffs),
        null_center = as.numeric(tru_calibration$null_center[names(tru_calibration$cutoffs)]),
        null_mad = as.numeric(tru_calibration$null_mad[names(tru_calibration$cutoffs)]),
        stringsAsFactors = FALSE
    )
    if (!is.null(tru_calibration$null_quantiles)) {
        df <- merge(df, tru_calibration$null_quantiles, by = "marker", all.x = TRUE, sort = FALSE)
    }
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(df, path, row.names = FALSE, quote = TRUE)
    invisible(path)
}

.load_tru_calibration <- function(x) {
    if (is.null(x) || (is.character(x) && (!length(x) || is.na(x[1]) || !nzchar(x[1])))) return(NULL)
    if (inherits(x, "spectreasy_tru_calibration")) return(x)
    if (is.character(x) && length(x) == 1L && file.exists(x)) {
        obj <- tryCatch(readRDS(x), error = function(e) NULL)
        if (inherits(obj, "spectreasy_tru_calibration")) return(obj)
    }
    NULL
}

.resolve_tru_calibration_for_unmixing <- function(M,
                                                  tru_calibration_file = NULL,
                                                  unmixing_matrix_file = NULL) {
    cal <- .load_tru_calibration(tru_calibration_file)
    if (!is.null(cal)) return(cal)
    cal <- .load_tru_calibration(attr(M, "tru_calibration"))
    if (!is.null(cal)) return(cal)
    if (!is.null(unmixing_matrix_file) && file.exists(unmixing_matrix_file)) {
        cal <- .load_tru_calibration(file.path(dirname(unmixing_matrix_file), "scc_tru_calibration.rds"))
        if (!is.null(cal)) return(cal)
    }
    stop(
        "method = \"TRU_Spectreasy\" requires a TRU calibration from an unstained/AF control. ",
        "Run unmix_controls() with an unstained control or supply tru_calibration_file.",
        call. = FALSE
    )
}

.spectreasy_tru_refit_event <- function(y,
                                        M,
                                        selected_af_row,
                                        active_marker_rows,
                                        method = "OLS",
                                        wls_noise = NULL,
                                        rwls_max_iter = 1L,
                                        cutoffs = NULL,
                                        active_margin = 1.0,
                                        prune_iterative = TRUE,
                                        max_iter = 5L,
                                        protected_marker_rows = integer(),
                                        tol = 1e-10) {
    y <- matrix(as.numeric(y), nrow = 1)
    active_marker_rows <- unique(as.integer(active_marker_rows))
    active_marker_rows <- active_marker_rows[active_marker_rows %in% seq_len(nrow(M))]
    protected_marker_rows <- unique(as.integer(protected_marker_rows))
    max_iter <- max(1L, suppressWarnings(as.integer(max_iter[1])))
    iterations <- 0L
    fallback <- FALSE

    repeat {
        iterations <- iterations + 1L
        rows_af <- c(active_marker_rows, selected_af_row)
        coeff_af <- .tru_ols(y, M[rows_af, , drop = FALSE], tol = tol)
        if (is.null(coeff_af) || any(!is.finite(coeff_af))) {
            fallback <- TRUE
            break
        }
        names_af <- rownames(M)[rows_af]
        af_marker_coeff <- stats::setNames(rep(0, length(active_marker_rows)), rownames(M)[active_marker_rows])
        if (length(active_marker_rows)) {
            af_marker_coeff[] <- coeff_af[1, seq_along(active_marker_rows)]
        }
        if (!isTRUE(prune_iterative) || length(active_marker_rows) == 0L || iterations >= max_iter || is.null(cutoffs)) {
            break
        }
        marker_names <- rownames(M)[active_marker_rows]
        thresholds <- active_margin * pmax(as.numeric(cutoffs[marker_names]), 0)
        thresholds[!is.finite(thresholds)] <- 0
        drop <- af_marker_coeff <= thresholds
        drop[active_marker_rows %in% protected_marker_rows] <- FALSE
        if (!any(drop)) break
        active_marker_rows <- active_marker_rows[!drop]
    }

    if (isTRUE(fallback)) {
        return(list(fallback = TRUE, iterations = iterations, active_marker_rows = active_marker_rows))
    }

    baseline_coeff <- if (length(active_marker_rows)) {
        .tru_ols(y, M[active_marker_rows, , drop = FALSE], tol = tol)
    } else {
        matrix(0, nrow = 1, ncol = 0L)
    }
    if (is.null(baseline_coeff) || any(!is.finite(baseline_coeff))) {
        return(list(fallback = TRUE, iterations = iterations, active_marker_rows = active_marker_rows))
    }

    A_af <- stats::setNames(rep(0, nrow(M)), rownames(M))
    A_base <- stats::setNames(rep(0, nrow(M)), rownames(M))
    rows_af <- c(active_marker_rows, selected_af_row)
    A_af[rows_af] <- as.numeric(coeff_af[1, seq_along(rows_af)])
    if (length(active_marker_rows)) {
        A_base[active_marker_rows] <- as.numeric(baseline_coeff[1, seq_along(active_marker_rows)])
    }
    list(
        A_af = A_af,
        A_base = A_base,
        active_marker_rows = active_marker_rows,
        iterations = iterations,
        fallback = FALSE
    )
}

.spectreasy_tru_event_gamma <- function(M, active_marker_rows, selected_af_row, tau, ridge = 1e-8, tol = 1e-10) {
    if (!length(active_marker_rows) || !is.finite(tau) || tau <= tol) {
        return(stats::setNames(rep(0, length(active_marker_rows)), rownames(M)[active_marker_rows]))
    }
    F <- M[active_marker_rows, , drop = FALSE]
    gram <- F %*% t(F) + diag(ridge, nrow(F))
    decoder <- tryCatch(solve(gram, F), error = function(e) NULL)
    if (is.null(decoder)) {
        return(stats::setNames(rep(0, length(active_marker_rows)), rownames(M)[active_marker_rows]))
    }
    impact <- as.numeric(decoder %*% M[selected_af_row, ])
    gamma <- impact^2 / (impact^2 + tau^2)
    gamma[!is.finite(gamma)] <- 0
    stats::setNames(gamma, rownames(M)[active_marker_rows])
}

.spectreasy_tru_global_gamma_blend <- function(A_base, A_af, gamma) {
    A_base + gamma * (A_af - A_base)
}

.tru_ucm_null_abundances <- function(tru_calibration, marker_names) {
    null_abundances <- tru_calibration$null_abundances
    if (is.null(null_abundances)) {
        stop(
            "method = \"TRU_Spectreasy\" now uses TRU-OLS + UCM and requires stored null abundances. ",
            "Run unmix_controls(..., tru_store_null_distributions = TRUE) to rebuild the TRU calibration.",
            call. = FALSE
        )
    }
    null_abundances <- as.matrix(null_abundances)
    missing <- setdiff(marker_names, colnames(null_abundances))
    if (length(missing)) {
        stop("TRU UCM calibration is missing null abundance column(s): ", paste(missing, collapse = ", "), call. = FALSE)
    }
    null_abundances[, marker_names, drop = FALSE]
}

.tru_ucm_values_for_event <- function(null_abundances, event_index, marker_names) {
    if (!length(marker_names)) {
        return(stats::setNames(numeric(0), character(0)))
    }
    null_index <- ((as.integer(event_index) - 1L) %% nrow(null_abundances)) + 1L
    out <- as.numeric(null_abundances[null_index, marker_names, drop = TRUE])
    out[!is.finite(out)] <- 0
    stats::setNames(out, marker_names)
}

.summarize_tru_diagnostics <- function(marker_names,
                                       cutoffs,
                                       active_mask,
                                       before,
                                       after,
                                       active_count,
                                       selected_af,
                                       refit_iterations,
                                       fallback) {
    marker_summary <- data.frame(
        marker = marker_names,
        cutoff = as.numeric(cutoffs[marker_names]),
        active_fraction = colMeans(active_mask),
        dropped_fraction = 1 - colMeans(active_mask),
        mean_abundance_before_tru = colMeans(before[, marker_names, drop = FALSE], na.rm = TRUE),
        mean_abundance_after_tru = colMeans(after[, marker_names, drop = FALSE], na.rm = TRUE),
        median_abundance_before_tru = apply(before[, marker_names, drop = FALSE], 2, stats::median, na.rm = TRUE),
        median_abundance_after_tru = apply(after[, marker_names, drop = FALSE], 2, stats::median, na.rm = TRUE),
        stringsAsFactors = FALSE,
        row.names = NULL
    )
    list(
        active_count = as.integer(active_count),
        dropped_count = as.integer(length(marker_names) - active_count),
        selected_af = as.integer(selected_af),
        refit_iterations = as.integer(refit_iterations),
        fallback = as.logical(fallback),
        marker_summary = marker_summary
    )
}

.spectreasy_tru_unmix <- function(Y,
                                  M,
                                  tru_calibration,
                                  spectreasy_weight_quantile = 0.9,
                                  tru_cutoff_quantile = NULL,
                                  tru_af_score_mode = c("thresholded_excess", "current"),
                                  tru_af_score_power = 2,
                                  tru_af_score_positive_only = TRUE,
                                  tru_af_score_tiebreak = c("current_leakage", "none"),
                                  tru_prune_iterative = TRUE,
                                  tru_max_iter = 5L,
                                  tru_active_margin = 1.0,
                                  tru_min_active_markers = 0L,
                                  tru_protected_markers = character(0),
                                  tru_gamma_mode = c("global", "event_specific"),
                                  tru_display = "ucm",
                                  write_tru_diagnostics = FALSE,
                                  tol = 1e-10) {
    M <- .as_reference_matrix(M, "M")
    Y <- as.matrix(Y)
    .validate_tru_calibration(tru_calibration)
    tru_af_score_mode <- match.arg(tru_af_score_mode)
    tru_af_score_tiebreak <- match.arg(tru_af_score_tiebreak)
    tru_gamma_mode <- .normalize_tru_gamma_mode(tru_gamma_mode)
    tru_display <- .normalize_tru_display(tru_display)
    marker_rows_logical <- .tru_marker_rows(M)
    marker_idx <- which(marker_rows_logical)
    af_idx <- which(!marker_rows_logical)
    if (!length(marker_idx) || !length(af_idx)) {
        stop("method = \"TRU_Spectreasy\" requires at least one marker row and one AF row.", call. = FALSE)
    }
    if (!identical(colnames(M), tru_calibration$detector_names)) {
        stop("TRU calibration detector names do not match the reference matrix.", call. = FALSE)
    }

    cutoffs <- .recompute_tru_cutoffs(tru_calibration, tru_cutoff_quantile)
    marker_names <- rownames(M)[marker_idx]
    missing_cutoffs <- setdiff(marker_names, names(cutoffs))
    if (length(missing_cutoffs)) {
        stop("TRU calibration is missing marker cutoff(s): ", paste(missing_cutoffs, collapse = ", "), call. = FALSE)
    }
    ucm_null_abundances <- .tru_ucm_null_abundances(tru_calibration, marker_names)

    o_ns <- if (isTRUE(tru_calibration$ns_correction) && !is.null(tru_calibration$o_ns)) {
        as.numeric(tru_calibration$o_ns[colnames(M)])
    } else {
        rep(0, ncol(M))
    }
    o_ns[!is.finite(o_ns)] <- 0
    Y_specific <- sweep(Y, 2, o_ns, "-")

    marker_spectra <- M[marker_idx, , drop = FALSE]
    af_spectra <- M[af_idx, , drop = FALSE]
    af_assignment_spectra <- af_spectra
    af_assignment_map <- seq_along(af_idx)
    if (nrow(af_assignment_spectra) == 1L) {
        af_assignment_spectra <- rbind(af_assignment_spectra, af_assignment_spectra)
        rownames(af_assignment_spectra) <- c("AF1", "AF2")
        af_assignment_map <- c(1L, 1L)
    }
    assignment <- .autospectral_assign_af_fluorophores_tru(
        Y = Y_specific,
        marker_spectra = marker_spectra,
        af_spectra = af_assignment_spectra,
        marker_cutoffs = cutoffs,
        score_mode = tru_af_score_mode,
        score_power = tru_af_score_power,
        positive_only = tru_af_score_positive_only,
        tiebreak = tru_af_score_tiebreak,
        active_margin = tru_active_margin,
        tol = tol,
        return_details = TRUE
    )
    selected_af <- af_assignment_map[assignment$selected_af]
    selected_af[!is.finite(selected_af)] <- 1L

    A <- matrix(0, nrow = nrow(Y), ncol = nrow(M), dimnames = list(NULL, rownames(M)))
    A_af_before <- A
    fitted <- matrix(0, nrow = nrow(Y), ncol = ncol(M), dimnames = list(NULL, colnames(M)))
    active_mask <- matrix(FALSE, nrow = nrow(Y), ncol = length(marker_idx), dimnames = list(NULL, marker_names))
    active_count <- integer(nrow(Y))
    iterations <- integer(nrow(Y))
    fallback <- logical(nrow(Y))
    protected_rows <- marker_idx[marker_names %in% tru_protected_markers]
    min_active <- max(0L, suppressWarnings(as.integer(tru_min_active_markers[1])))
    active_margin <- suppressWarnings(as.numeric(tru_active_margin[1]))
    if (!is.finite(active_margin) || active_margin < 0) active_margin <- 1
    global_gamma <- .decoder_projected_af_marker_weights(M, spectreasy_weight_quantile)
    tau <- attr(global_gamma, "tau")

    for (i in seq_len(nrow(Y))) {
        adj <- assignment$adjusted_marker_abundances_for_selected_af[i, marker_names]
        thresholds <- active_margin * pmax(as.numeric(cutoffs[marker_names]), 0)
        thresholds[!is.finite(thresholds)] <- 0
        active <- which(is.finite(adj) & adj > thresholds)
        if (length(protected_rows)) {
            active <- union(active, match(rownames(M)[protected_rows], marker_names))
        }
        if (min_active > 0L && length(active) < min_active) {
            order_idx <- order(adj, decreasing = TRUE, na.last = NA)
            active <- union(active, head(order_idx, min_active))
        }
        active <- active[!is.na(active)]
        active_rows <- marker_idx[active]
        selected_af_row <- af_idx[[selected_af[[i]]]]
        refit <- .spectreasy_tru_refit_event(
            y = Y_specific[i, ],
            M = M,
            selected_af_row = selected_af_row,
            active_marker_rows = active_rows,
            cutoffs = cutoffs,
            active_margin = active_margin,
            prune_iterative = tru_prune_iterative,
            max_iter = tru_max_iter,
            protected_marker_rows = protected_rows,
            tol = tol
        )
        if (isTRUE(refit$fallback)) {
            fallback[[i]] <- TRUE
            fallback_fit <- .autospectral_unmix_legacy(Y = Y_specific[i, , drop = FALSE], M = M)
            A_af_before[i, ] <- fallback_fit$A[1, ]
            A[i, ] <- fallback_fit$A[1, ]
            fitted[i, ] <- A[i, ] %*% M + o_ns
            next
        }
        active_rows <- refit$active_marker_rows
        active_names <- rownames(M)[active_rows]
        active_mask[i, active_names] <- TRUE
        active_count[[i]] <- length(active_rows)
        iterations[[i]] <- refit$iterations
        A_af_before[i, ] <- refit$A_af
        gamma <- if (identical(tru_gamma_mode, "event_specific")) {
            .spectreasy_tru_event_gamma(M, active_rows, selected_af_row, tau = tau)
        } else {
            global_gamma[active_names]
        }
        blended <- refit$A_base
        if (length(active_names)) {
            blended[active_names] <- .spectreasy_tru_global_gamma_blend(
                A_base = refit$A_base[active_names],
                A_af = refit$A_af[active_names],
                gamma = gamma[active_names]
            )
        }
        blended[selected_af_row] <- refit$A_af[[rownames(M)[selected_af_row]]]
        inactive_names <- setdiff(marker_names, active_names)
        blended[inactive_names] <- .tru_ucm_values_for_event(ucm_null_abundances, i, inactive_names)
        A[i, ] <- blended
        fitted[i, ] <- A[i, ] %*% M + o_ns
    }

    tru_info <- .summarize_tru_diagnostics(
        marker_names = marker_names,
        cutoffs = cutoffs,
        active_mask = active_mask,
        before = A_af_before,
        after = A,
        active_count = active_count,
        selected_af = selected_af,
        refit_iterations = iterations,
        fallback = fallback
    )
    autospectral_fit <- list(
        A = A_af_before,
        fitted = A_af_before %*% M,
        af_index = selected_af,
        selected_af_row = rownames(M)[af_idx[selected_af]]
    )
    list(
        A = A,
        fitted = fitted,
        autospectral_fit = autospectral_fit,
        spectreasy_decoder_weights = global_gamma,
        tru_info = tru_info
    )
}
