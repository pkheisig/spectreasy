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
                                    choices = c("AutoSpectral", "OLS", "NNLS", "WLS", "RWLS")) {
    method_raw <- toupper(trimws(as.character(method[1])))
    if (length(method_raw) == 0 || is.na(method_raw) || !nzchar(method_raw)) {
        method_raw <- "AUTOSPECTRAL"
    }
    method_raw <- gsub("[_-]", "", method_raw)
    lookup <- c(
        AUTOSPECTRAL = "AutoSpectral",
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
#' @param n_threads Positive integer; number of threads to use for event-wise
#'   multi-AF WLS/RWLS unmixing. The default, 1, keeps execution single-threaded.
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
                           n_threads = 1L) {
    M <- .as_reference_matrix(M, "M")
    rwls_max_iter <- .normalize_rwls_max_iter(rwls_max_iter)
    n_threads <- .normalize_unmix_threads(n_threads)
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

    if (sum(af_match) > 1) {
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
        if (method == "OLS") {
            # Ordinary Least Squares: A = Y * M^T * (M * M^T)^-1
            matrix_to_invert <- M %*% Mt
            if (rcond(matrix_to_invert) < 1e-10) {
                stop("Reference Matrix is singular. You likely have collinear spectra.")
            }
            A <- Y %*% Mt %*% solve(matrix_to_invert)
        } else if (method == "NNLS") {
            A <- spectreasy_nnls_unmix_cpp(Y = Y, M = M)
        } else if (method == "WLS") {
            A <- spectreasy_wls_unmix_cpp(
                Y = Y,
                M = M,
                noise_floor = wls_noise$noise_floor,
                signal_scale = wls_noise$signal_scale,
                max_weight_ratio = wls_noise$max_weight_ratio
            )
        } else if (method == "RWLS") {
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

    Fitted <- A %*% M
    R <- Y - Fitted
    out <- as.data.frame(A)
    colnames(out) <- rownames(M)

    out <- .append_passthrough_parameters(out, full_data, detector_names = detectors)

    if (!is.null(file_name)) out$File <- file_name

    if (return_residuals) {
        res <- list(data = out, residuals = R)
        attr(res, "method") <- method
        attr(res, "reference_matrix") <- M
        return(res)
    } else {
        return(out)
    }
}
