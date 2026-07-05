.normalize_rwls_max_iter <- function(rwls_max_iter) {
    rwls_max_iter <- suppressWarnings(as.integer(rwls_max_iter[1]))
    if (!is.finite(rwls_max_iter) || rwls_max_iter < 1L) {
        stop("rwls_max_iter must be an integer >= 1.", call. = FALSE)
    }
    rwls_max_iter
}

.normalize_unmix_threads <- function(n_threads) {
    n_threads <- suppressWarnings(as.integer(n_threads[1]))
    if (!is.finite(n_threads) || n_threads < 1L) {
        stop("n_threads must be an integer >= 1.", call. = FALSE)
    }
    n_threads
}

.normalize_unmix_method <- function(method,
                                    choices = c("AutoSpectral", "Spectreasy", "OLS", "NNLS", "WLS", "RWLS")) {
    method_raw <- toupper(trimws(as.character(method[1])))
    if (length(method_raw) == 0 || is.na(method_raw) || !nzchar(method_raw)) {
        method_raw <- "AUTOSPECTRAL"
    }
    method_raw <- gsub("[_-]", "", method_raw)
    lookup <- c(
        AUTOSPECTRAL = "AutoSpectral",
        SPECTREASY = "Spectreasy",
        OLS = "OLS",
        NNLS = "NNLS",
        WLS = "WLS",
        RWLS = "RWLS"
    )
    out <- unname(lookup[method_raw])
    allowed_lookup <- unname(lookup[gsub("[_-]", "", toupper(choices))])
    if (is.na(out) || !(out %in% allowed_lookup)) {
        stop("method must be one of: ", paste(choices, collapse = ", "), call. = FALSE)
    }
    out
}

.solver_method_for_unmix <- function(method) {
    if (method %in% c("AutoSpectral", "Spectreasy")) "OLS" else method
}

.unmix_with_fixed_reference <- function(Y, M, method, wls_noise, rwls_max_iter) {
    if (method == "OLS") {
        matrix_to_invert <- M %*% t(M)
        if (rcond(matrix_to_invert) < 1e-10) {
            stop("Reference Matrix is singular. You likely have collinear spectra.", call. = FALSE)
        }
        return(Y %*% t(M) %*% solve(matrix_to_invert))
    }
    if (method == "NNLS") {
        return(spectreasy_nnls_unmix_cpp(Y = Y, M = M))
    }
    if (method == "WLS") {
        return(spectreasy_wls_unmix_cpp(
            Y = Y,
            M = M,
            noise_floor = wls_noise$noise_floor,
            signal_scale = wls_noise$signal_scale,
            max_weight_ratio = wls_noise$max_weight_ratio
        ))
    }
    if (method == "RWLS") {
        return(spectreasy_rwls_unmix_cpp(
            Y = Y,
            M = M,
            noise_floor = wls_noise$noise_floor,
            signal_scale = wls_noise$signal_scale,
            max_weight_ratio = wls_noise$max_weight_ratio,
            max_iter = rwls_max_iter
        ))
    }
    stop("method must be one of: AutoSpectral, Spectreasy, OLS, NNLS, WLS, RWLS", call. = FALSE)
}

.is_autospectral_style_method <- function(method) {
    method %in% c("AutoSpectral", "Spectreasy")
}

.normalize_spectreasy_weight_quantile <- function(spectreasy_weight_quantile) {
    q <- suppressWarnings(as.numeric(spectreasy_weight_quantile[1]))
    if (!is.finite(q) || is.na(q) || q < 0 || q > 1) {
        stop("spectreasy_weight_quantile must be a numeric value between 0 and 1.", call. = FALSE)
    }
    q
}

.decoder_projected_af_marker_weights <- function(M,
                                                 spectreasy_weight_quantile = 0.9,
                                                 ridge = 1e-8,
                                                 tol = 1e-12) {
    M <- .as_reference_matrix(M, "M")
    spectreasy_weight_quantile <- .normalize_spectreasy_weight_quantile(spectreasy_weight_quantile)
    af_rows <- grepl("^AF($|_)", rownames(M), ignore.case = TRUE)
    marker_names <- rownames(M)[!af_rows]
    weights <- stats::setNames(rep(0, length(marker_names)), marker_names)
    if (!any(af_rows) || length(marker_names) == 0L) {
        return(weights)
    }

    marker_spectra <- M[marker_names, , drop = FALSE]
    af_spectra <- M[af_rows, , drop = FALSE]
    S <- t(marker_spectra)
    gram <- crossprod(S) + diag(ridge, ncol(S))
    decoder <- tryCatch(
        solve(gram, t(S)),
        error = function(e) NULL
    )
    if (is.null(decoder)) {
        return(weights)
    }
    rownames(decoder) <- marker_names

    impact_by_af <- decoder %*% t(af_spectra)
    impact <- sqrt(rowMeans(impact_by_af^2))
    tau <- stats::quantile(impact, probs = spectreasy_weight_quantile, na.rm = TRUE, names = FALSE, type = 7)
    if (!is.finite(tau) || tau <= tol) {
        return(weights)
    }

    weights <- impact^2 / (impact^2 + tau^2)
    weights[!is.finite(weights)] <- 0
    weights <- stats::setNames(as.numeric(weights), names(impact))
    attr(weights, "impact") <- stats::setNames(as.numeric(impact), names(impact))
    attr(weights, "tau") <- tau
    attr(weights, "quantile") <- spectreasy_weight_quantile
    weights
}

.spectreasy_marker_only_baseline_fit <- function(Y, M) {
    af_match <- grepl("^AF($|_)", rownames(M), ignore.case = TRUE)
    marker_M <- M[!af_match, , drop = FALSE]
    if (nrow(marker_M) == 0L) {
        stop("Spectreasy unmixing requires at least one non-AF reference row.", call. = FALSE)
    }
    baseline <- .unmix_with_fixed_reference(
        Y = Y,
        M = marker_M,
        method = "OLS",
        wls_noise = list(noise_floor = numeric(0), signal_scale = numeric(0), max_weight_ratio = .default_wls_max_weight_ratio()),
        rwls_max_iter = 1L
    )
    colnames(baseline) <- rownames(marker_M)
    baseline
}

.apply_spectreasy_decoder_blend <- function(Y, M, autospectral_A, spectreasy_weight_quantile = 0.9) {
    baseline_A <- .spectreasy_marker_only_baseline_fit(Y = Y, M = M)
    out_A <- autospectral_A
    marker_weights <- .decoder_projected_af_marker_weights(
        M,
        spectreasy_weight_quantile = spectreasy_weight_quantile
    )
    marker_cols <- intersect(names(marker_weights), colnames(out_A))
    for (marker in marker_cols) {
        weight <- marker_weights[[marker]]
        out_A[, marker] <- baseline_A[, marker] + weight * (autospectral_A[, marker] - baseline_A[, marker])
    }
    attr(out_A, "spectreasy_decoder_weights") <- marker_weights
    out_A
}

.autospectral_assign_af_fluorophores <- function(Y, marker_spectra, af_spectra, tol = 1e-10) {
    Y <- as.matrix(Y)
    marker_spectra <- as.matrix(marker_spectra)
    af_spectra <- as.matrix(af_spectra)
    if (nrow(marker_spectra) == 0L || nrow(af_spectra) == 0L) {
        return(rep(NA_integer_, nrow(Y)))
    }
    if (!identical(colnames(marker_spectra), colnames(af_spectra))) {
        stop("AutoSpectral marker and AF spectra detector columns do not match.", call. = FALSE)
    }

    matrix_to_invert <- marker_spectra %*% t(marker_spectra)
    if (rcond(matrix_to_invert) < tol) {
        stop("Marker-only reference matrix is singular; AutoSpectral AF assignment cannot proceed.", call. = FALSE)
    }
    unmixing_matrix <- solve(matrix_to_invert, marker_spectra)
    marker_t <- t(marker_spectra)
    v_library <- unmixing_matrix %*% t(af_spectra)
    r_library <- t(af_spectra) - marker_t %*% v_library
    denominator <- colSums(r_library^2)
    valid <- is.finite(denominator) & denominator > tol
    if (!any(valid)) {
        return(rep(1L, nrow(Y)))
    }

    numerator <- Y %*% r_library
    k_matrix <- matrix(0, nrow = nrow(Y), ncol = nrow(af_spectra))
    k_matrix[, valid] <- sweep(numerator[, valid, drop = FALSE], 2, denominator[valid], "/")
    unmixed_markers <- Y %*% t(unmixing_matrix)
    error_matrix <- matrix(Inf, nrow = nrow(Y), ncol = nrow(af_spectra))
    for (i in which(valid)) {
        adjusted <- unmixed_markers - k_matrix[, i, drop = FALSE] %*% t(v_library[, i, drop = FALSE])
        error_matrix[, i] <- rowSums(abs(adjusted))
    }
    max.col(-error_matrix, ties.method = "first")
}

.autospectral_unmix_legacy <- function(Y, M, tol = 1e-10) {
    Y <- as.matrix(Y)
    M <- .as_reference_matrix(M, "M")
    af_match <- grepl("^AF($|_)", rownames(M), ignore.case = TRUE)
    marker_idx <- which(!af_match)
    af_idx <- which(af_match)
    if (length(af_idx) == 0L) {
        A <- .unmix_with_fixed_reference(
            Y = Y,
            M = M,
            method = "OLS",
            wls_noise = list(noise_floor = numeric(0), signal_scale = numeric(0), max_weight_ratio = .default_wls_max_weight_ratio()),
            rwls_max_iter = 1L
        )
        return(list(
            A = A,
            fitted = A %*% M,
            af_index = rep(NA_integer_, nrow(Y)),
            selected_af_row = rep(NA_character_, nrow(Y))
        ))
    }
    if (length(marker_idx) == 0L) {
        stop("AutoSpectral unmixing requires at least one non-AF reference row.", call. = FALSE)
    }

    marker_spectra <- M[marker_idx, , drop = FALSE]
    af_spectra <- M[af_idx, , drop = FALSE]
    af_assignment_spectra <- af_spectra
    af_assignment_map <- seq_along(af_idx)
    if (nrow(af_assignment_spectra) == 1L) {
        # AutoSpectral expects at least two AF rows. Duplicating one row keeps
        # one-band AF behavior equivalent to a marker+AF OLS refit.
        af_assignment_spectra <- rbind(af_assignment_spectra, af_assignment_spectra)
        rownames(af_assignment_spectra) <- c("AF1", "AF2")
        af_assignment_map <- c(1L, 1L)
    }

    assignments_raw <- .autospectral_assign_af_fluorophores(
        Y = Y,
        marker_spectra = marker_spectra,
        af_spectra = af_assignment_spectra,
        tol = tol
    )
    selected_af <- af_assignment_map[assignments_raw]
    selected_af[!is.finite(selected_af)] <- 1L

    A <- matrix(0, nrow = nrow(Y), ncol = nrow(M), dimnames = list(NULL, rownames(M)))
    fitted <- matrix(0, nrow = nrow(Y), ncol = ncol(M), dimnames = list(NULL, colnames(M)))
    for (af_i in unique(selected_af)) {
        rows <- c(marker_idx, af_idx[[af_i]])
        event_idx <- which(selected_af == af_i)
        X <- M[rows, , drop = FALSE]
        coeff <- .unmix_with_fixed_reference(
            Y = Y[event_idx, , drop = FALSE],
            M = X,
            method = "OLS",
            wls_noise = list(noise_floor = numeric(0), signal_scale = numeric(0), max_weight_ratio = .default_wls_max_weight_ratio()),
            rwls_max_iter = 1L
        )
        A[event_idx, rows] <- coeff
        fitted[event_idx, ] <- coeff %*% X
    }

    list(
        A = A,
        fitted = fitted,
        af_index = selected_af,
        selected_af_row = rownames(M)[af_idx[selected_af]]
    )
}

#' Calculate unmixing residuals
#'
#' @param flow_frame A flowFrame object with raw fluorescence data
#' @param M Reference matrix (fluorophores x detectors)
#' @param file_name Optional file name to add to output
#' @param method Unmixing method: `"WLS"` (default), `"RWLS"`, `"OLS"`,
#'   `"NNLS"`, `"AutoSpectral"`, or `"Spectreasy"`. `AutoSpectral` assigns the best AF
#'   spectrum per event with the AutoSpectral fluorophore-leakage score, refits
#'   marker + selected-AF OLS, and applies spectral-variant optimization when a
#'   variant library is supplied. `Spectreasy` uses the same AutoSpectral-style
#'   fit, then blends marker abundances with a marker-only OLS anchor using
#'   decoder-projected AF-impact weights derived from the reference matrix.
#' @param return_residuals Logical. If TRUE, returns a list containing the unmixed
#'   data and the detector residual matrix.
#' @param background_noise Scalar or detector-length vector used as the WLS noise
#'   floor when `M` does not carry `attr(M, "detector_noise")`; the built-in
#'   fallback is 125 raw detector units.
#' @param wls_signal_scale Scalar or detector-length vector multiplying the
#'   non-negative event signal in the WLS variance model.
#' @param wls_max_weight_ratio Maximum detector weight ratio allowed per event.
#' @param rwls_max_iter Positive integer; number of robust reweighting
#'   iterations used when `method = "RWLS"`. The default, 1, preserves the
#'   historical behavior.
#' @param n_threads Positive integer; number of threads to use for event-wise
#'   multi-AF WLS/RWLS unmixing. The default, 1, keeps execution single-threaded.
#' @param spectral_variant_library Optional per-fluorophore spectral-variant
#'   library learned from single-color controls. Used only with
#'   `method = "AutoSpectral"` or `"Spectreasy"`.
#' @param spectral_variant_top_k Number of best variant candidates to test per
#'   positive fluorophore.
#' @param spectral_variant_min_abundance Minimum unmixed abundance for a
#'   fluorophore to be eligible for variant testing.
#' @param spectral_variant_positive_fraction Additional positivity threshold
#'   as a fraction of the event's strongest fluorophore abundance.
#' @param spectral_variant_min_improvement Minimum fractional residual
#'   improvement required before accepting a cell-specific variant refit.
#' @param spectreasy_weight_quantile Numeric in `[0, 1]`; only accepted when
#'   `method = "Spectreasy"`. Controls the quantile of decoder-projected AF
#'   impacts used as the soft-saturation scale for marker-specific AutoSpectral
#'   mixing. The default, `0.9`, uses the 90th percentile of marker AF impacts.
#' @return Data frame with unmixed abundances and retained acquisition parameters
#'   (`Time` plus all `FSC*`/`SSC*` columns, when available).
#'         If return_residuals=TRUE, returns a list with [[data]] and [[residuals]].
#' @importFrom RcppParallel defaultNumThreads
#' @examples
#' M_demo <- rbind(
#'   FITC = c(1.00, 0.20, 0.05),
#'   PE = c(0.10, 1.00, 0.20),
#'   APC = c(0.05, 0.15, 1.00)
#' )
#' colnames(M_demo) <- c("B2-A", "YG1-A", "R1-A")
#'
#' marker_signal <- matrix(rexp(300 * nrow(M_demo), rate = 8), ncol = nrow(M_demo))
#' colnames(marker_signal) <- rownames(M_demo)
#' marker_signal[, "FITC"] <- rexp(300, rate = 0.6) + 2
#' raw_signal <- marker_signal %*% M_demo +
#'   matrix(rnorm(300 * ncol(M_demo), sd = 0.03), ncol = ncol(M_demo))
#'
#' exprs_mat <- cbind(
#'   raw_signal,
#'   Time = seq_len(300),
#'   "FSC-A" = rnorm(300, mean = 90000, sd = 7000),
#'   "SSC-A" = rnorm(300, mean = 45000, sd = 5000)
#' )
#' colnames(exprs_mat)[seq_len(ncol(M_demo))] <- colnames(M_demo)
#' ff <- flowCore::flowFrame(exprs_mat)
#'
#' res <- calc_residuals(ff, M_demo, method = "OLS", return_residuals = TRUE)
#' head(res$data)
#' @export
calc_residuals <- function(flow_frame,
                           M,
                           file_name = NULL,
                           method = "WLS",
                           return_residuals = FALSE,
                           background_noise = .default_wls_background_noise(),
                           wls_signal_scale = .default_wls_signal_scale(),
                           wls_max_weight_ratio = .default_wls_max_weight_ratio(),
                           rwls_max_iter = 1L,
                           n_threads = 1L,
                           spectral_variant_library = NULL,
                           spectral_variant_top_k = 3L,
                           spectral_variant_min_abundance = 1,
                           spectral_variant_positive_fraction = 0.02,
                           spectral_variant_min_improvement = 0.01,
                           spectreasy_weight_quantile = 0.9) {
    spectreasy_weight_quantile_missing <- missing(spectreasy_weight_quantile)
    M <- .as_reference_matrix(M, "M")
    rwls_max_iter <- .normalize_rwls_max_iter(rwls_max_iter)
    n_threads <- .normalize_unmix_threads(n_threads)
    method <- .normalize_unmix_method(method)
    if (!identical(method, "Spectreasy") && !spectreasy_weight_quantile_missing) {
        stop("spectreasy_weight_quantile is only accepted when method = \"Spectreasy\".", call. = FALSE)
    }
    if (identical(method, "Spectreasy")) {
        spectreasy_weight_quantile <- .normalize_spectreasy_weight_quantile(spectreasy_weight_quantile)
    }
    solver_method <- .solver_method_for_unmix(method)
    use_spectral_variants <- .is_autospectral_style_method(method) &&
        .spectral_variant_library_has_variants(spectral_variant_library)
    full_data <- flowCore::exprs(flow_frame)
    detectors <- colnames(M)
    
    # Check for missing detectors
    missing <- setdiff(detectors, colnames(full_data))
    if (length(missing) > 0) {
        stop("Detectors in reference matrix not found in flow_frame: ", paste(missing, collapse = ", "))
    }
    
    Y <- full_data[, detectors, drop = FALSE]

    Mt <- t(M)

    af_match <- grepl("^AF($|_)", rownames(M), ignore.case = TRUE)
    wls_noise <- if (solver_method %in% c("WLS", "RWLS")) {
        .resolve_wls_noise_parameters(
            M = M,
            background_noise = background_noise,
            signal_scale = wls_signal_scale,
            max_weight_ratio = wls_max_weight_ratio
        )
    } else {
        list(noise_floor = numeric(0), signal_scale = numeric(0), max_weight_ratio = .default_wls_max_weight_ratio())
    }

    autospectral_fit <- NULL
    if (.is_autospectral_style_method(method)) {
        autospectral_fit <- .autospectral_unmix_legacy(Y = Y, M = M)
        A <- autospectral_fit$A
        Fitted <- autospectral_fit$fitted
    } else if (sum(af_match) > 1) {
        fluor_idx <- which(!af_match) - 1
        af_idx <- which(af_match) - 1
        A <- spectreasy_unmix_best_af_cpp(
            Y = Y,
            M = M,
            fluor_idx = fluor_idx,
            af_idx = af_idx,
            method = solver_method,
            noise_floor = wls_noise$noise_floor,
            signal_scale = wls_noise$signal_scale,
            max_weight_ratio = wls_noise$max_weight_ratio,
            rwls_max_iter = rwls_max_iter,
            n_threads = n_threads
        )
    } else {
        if (solver_method == "OLS") {
            # Ordinary Least Squares: A = Y * M^T * (M * M^T)^-1
            matrix_to_invert <- M %*% Mt
            if (rcond(matrix_to_invert) < 1e-10) {
                stop("Reference Matrix is singular. You likely have collinear spectra.")
            }
            A <- Y %*% Mt %*% solve(matrix_to_invert)
        } else if (solver_method == "NNLS") {
            A <- spectreasy_nnls_unmix_cpp(Y = Y, M = M)
        } else if (solver_method == "WLS") {
            A <- spectreasy_wls_unmix_cpp(
                Y = Y,
                M = M,
                noise_floor = wls_noise$noise_floor,
                signal_scale = wls_noise$signal_scale,
                max_weight_ratio = wls_noise$max_weight_ratio
            )
        } else if (solver_method == "RWLS") {
            A <- spectreasy_rwls_unmix_cpp(
                Y = Y,
                M = M,
                noise_floor = wls_noise$noise_floor,
                signal_scale = wls_noise$signal_scale,
                max_weight_ratio = wls_noise$max_weight_ratio,
                max_iter = rwls_max_iter
            )
        }
    }

    variant_info <- NULL
    if (use_spectral_variants) {
        variant_fit <- .apply_spectral_variant_optimization(
            Y = Y,
            M = M,
            A = A,
            method = solver_method,
            wls_noise = wls_noise,
            rwls_max_iter = rwls_max_iter,
            spectral_variant_library = spectral_variant_library,
            top_k = spectral_variant_top_k,
            min_abundance = spectral_variant_min_abundance,
            positive_fraction = spectral_variant_positive_fraction,
            min_improvement_fraction = spectral_variant_min_improvement
        )
        A <- variant_fit$A
        Fitted <- variant_fit$fitted
        variant_info <- variant_fit$info
    }
    spectreasy_decoder_weights <- NULL
    if (identical(method, "Spectreasy")) {
        A <- .apply_spectreasy_decoder_blend(
            Y = Y,
            M = M,
            autospectral_A = A,
            spectreasy_weight_quantile = spectreasy_weight_quantile
        )
        spectreasy_decoder_weights <- attr(A, "spectreasy_decoder_weights")
        Fitted <- A %*% M
    } else if (is.null(autospectral_fit)) {
        Fitted <- A %*% M
    }
    R <- Y - Fitted
    out <- as.data.frame(A)
    colnames(out) <- rownames(M)
    if (!is.null(autospectral_fit) && any(is.finite(autospectral_fit$af_index))) {
        out[["AF Index"]] <- autospectral_fit$af_index
    }

    out <- .append_passthrough_parameters(out, full_data, detector_names = detectors)

    if (!is.null(file_name)) out$File <- file_name

    if (return_residuals) {
        res <- list(data = out, residuals = R)
        if (!is.null(variant_info)) res$spectral_variant_info <- variant_info
        if (!is.null(spectreasy_decoder_weights)) res$spectreasy_decoder_weights <- spectreasy_decoder_weights
        attr(res, "method") <- method
        attr(res, "reference_matrix") <- M
        return(res)
    } else {
        return(out)
    }
}
