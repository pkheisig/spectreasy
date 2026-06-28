.normalize_rwls_max_iter <- function(rwls_max_iter) {
    rwls_max_iter <- suppressWarnings(as.integer(rwls_max_iter[1]))
    if (!is.finite(rwls_max_iter) || rwls_max_iter < 1L) {
        stop("rwls_max_iter must be an integer >= 1.", call. = FALSE)
    }
    rwls_max_iter
}

.available_unmix_threads <- function() {
    max(1L, as.integer(RcppParallel::defaultNumThreads()))
}

.normalize_unmix_threads <- function(multithreading = FALSE, n_threads = "auto") {
    multithreading <- isTRUE(multithreading[1])
    if (!multithreading) {
        return(1L)
    }

    available_threads <- .available_unmix_threads()
    if (is.character(n_threads) && length(n_threads) > 0 && identical(tolower(trimws(n_threads[1])), "auto")) {
        return(available_threads)
    }

    n_threads <- suppressWarnings(as.integer(n_threads[1]))
    if (!is.finite(n_threads) || n_threads < 1L) {
        stop("n_threads must be \"auto\" or an integer >= 1.", call. = FALSE)
    }
    min(n_threads, available_threads)
}

.ols_candidate_model_ok <- function(M, tol = 1e-10) {
    M <- as.matrix(M)
    cross <- M %*% t(M)
    is.finite(rcond(cross)) && rcond(cross) >= tol
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
    stop("method must be 'OLS', 'NNLS', 'WLS', or 'RWLS'", call. = FALSE)
}

.unmix_selected_af_groups <- function(Y,
                                      M,
                                      af_match,
                                      method,
                                      wls_noise,
                                      rwls_max_iter,
                                      af_assignment = "projection") {
    fluor_idx <- which(!af_match)
    af_idx <- which(af_match)
    marker_M <- M[fluor_idx, , drop = FALSE]
    af_M <- M[af_idx, , drop = FALSE]

    if (method == "OLS") {
        candidate_ok <- vapply(seq_along(af_idx), function(k) {
            .ols_candidate_model_ok(M[c(fluor_idx, af_idx[[k]]), , drop = FALSE])
        }, logical(1))
        if (!any(candidate_ok)) {
            stop(
                "No usable AF candidate model for OLS unmixing; all candidate matrices are singular or ill-conditioned.",
                call. = FALSE
            )
        }
    }

    assignments <- .assign_af_candidates(
        Y = Y,
        marker_M = marker_M,
        af_M = af_M,
        af_assignment = af_assignment
    )

    A <- matrix(0, nrow = nrow(Y), ncol = nrow(M), dimnames = list(NULL, rownames(M)))
    for (k in sort(unique(assignments))) {
        idx <- which(assignments == k)
        if (!length(idx) || is.na(k) || k < 1L || k > length(af_idx)) {
            next
        }
        rows <- c(fluor_idx, af_idx[[k]])
        X <- M[rows, , drop = FALSE]
        if (method == "OLS" && !.ols_candidate_model_ok(X)) {
            stop(
                "The selected AF candidate model for event ",
                idx[[1]],
                " is singular or ill-conditioned.",
                call. = FALSE
            )
        }
        coeff_i <- .unmix_with_fixed_reference(
            Y = Y[idx, , drop = FALSE],
            M = X,
            method = method,
            wls_noise = wls_noise,
            rwls_max_iter = rwls_max_iter
        )
        A[idx, fluor_idx] <- coeff_i[, seq_along(fluor_idx), drop = FALSE]
        A[idx, af_idx[[k]]] <- coeff_i[, length(rows), drop = TRUE]
    }
    A
}

#' Calculate unmixing residuals
#'
#' @param flow_frame A flowFrame object with raw fluorescence data
#' @param M Reference matrix (fluorophores x detectors)
#' @param file_name Optional file name to add to output
#' @param method Unmixing method: "WLS" (default), "RWLS", "OLS", or "NNLS".
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
#' @param multithreading Logical; if `TRUE`, allow event-wise multi-AF WLS/RWLS
#'   unmixing to use multiple threads. The default, `FALSE`, keeps execution
#'   single-threaded.
#' @param n_threads `"auto"` or positive integer; thread count to use when
#'   `multithreading = TRUE`. `"auto"` uses `RcppParallel::defaultNumThreads()`.
#'   Integers larger than the available thread count are clipped to the
#'   available count.
#' @param af_assignment How to choose one AF row per event when `M` contains
#'   multiple AF rows. `"projection"` (default) projects each AF candidate into
#'   fluorophore space and picks the candidate that brings apparent fluorophore
#'   abundance closest to zero. `"residual_alignment"` uses detector residual
#'   alignment as an explicit alternative. `"legacy"` uses the previous C++
#'   residual/covariance selector.
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
                           multithreading = FALSE,
                           n_threads = "auto",
                           af_assignment = "projection") {
    M <- .as_reference_matrix(M, "M")
    rwls_max_iter <- .normalize_rwls_max_iter(rwls_max_iter)
    n_threads <- .normalize_unmix_threads(multithreading = multithreading, n_threads = n_threads)
    af_assignment <- .normalize_af_assignment(af_assignment, choices = c("projection", "residual_alignment", "legacy"))
    full_data <- flowCore::exprs(flow_frame)
    detectors <- colnames(M)
    
    # Check for missing detectors
    missing <- setdiff(detectors, colnames(full_data))
    if (length(missing) > 0) {
        stop("Detectors in reference matrix not found in flow_frame: ", paste(missing, collapse = ", "))
    }
    
    Y <- full_data[, detectors, drop = FALSE]

    Mt <- t(M)
    method <- toupper(method)
    if (!method %in% c("OLS", "NNLS", "WLS", "RWLS")) {
        stop("method must be 'OLS', 'NNLS', 'WLS', or 'RWLS'")
    }

    af_match <- grepl("^AF($|_)", rownames(M), ignore.case = TRUE)
    wls_noise <- if (method %in% c("WLS", "RWLS")) {
        .resolve_wls_noise_parameters(
            M = M,
            background_noise = background_noise,
            signal_scale = wls_signal_scale,
            max_weight_ratio = wls_max_weight_ratio
        )
    } else {
        list(noise_floor = numeric(0), signal_scale = numeric(0), max_weight_ratio = .default_wls_max_weight_ratio())
    }

    use_eventwise_af_solver <- sum(af_match) > 1 && sum(!af_match) > 0

    if (use_eventwise_af_solver) {
        if (identical(af_assignment, "legacy")) {
            fluor_idx <- which(!af_match) - 1
            af_idx <- which(af_match) - 1
            A <- spectreasy_unmix_best_af_cpp(
                Y = Y,
                M = M,
                fluor_idx = fluor_idx,
                af_idx = af_idx,
                method = method,
                noise_floor = wls_noise$noise_floor,
                signal_scale = wls_noise$signal_scale,
                max_weight_ratio = wls_noise$max_weight_ratio,
                rwls_max_iter = rwls_max_iter,
                n_threads = n_threads
            )
        } else {
            A <- .unmix_selected_af_groups(
                Y = Y,
                M = M,
                af_match = af_match,
                method = method,
                wls_noise = wls_noise,
                rwls_max_iter = rwls_max_iter,
                af_assignment = af_assignment
            )
        }
    } else {
        A <- .unmix_with_fixed_reference(
            Y = Y,
            M = M,
            method = method,
            wls_noise = wls_noise,
            rwls_max_iter = rwls_max_iter
        )
    }

    Fitted <- A %*% M
    R <- Y - Fitted
    out <- as.data.frame(A)
    colnames(out) <- rownames(M)

    out <- .append_passthrough_parameters(out, full_data, detector_names = detectors)

    if (!is.null(file_name)) out$File <- file_name

    if (return_residuals) {
        return(list(data = out, residuals = R))
    } else {
        return(out)
    }
}
