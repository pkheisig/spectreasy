#' Calculate unmixing residuals
#'
#' @param flow_frame A flowFrame object with raw fluorescence data
#' @param M Reference matrix (fluorophores x detectors)
#' @param file_name Optional file name to add to output
#' @param method Unmixing method: "OLS" (default), "NNLS", or "WLS".
#' @param return_residuals Logical. If TRUE, returns a list containing the unmixed
#'   data and the detector residual matrix.
#' @return Data frame with unmixed abundances and retained acquisition parameters
#'   (`Time` plus all `FSC*`/`SSC*` columns, when available).
#'         If return_residuals=TRUE, returns a list with [[data]] and [[residuals]].
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
calc_residuals <- function(flow_frame, M, file_name = NULL, method = "OLS", 
                          return_residuals = FALSE) {
    M <- .as_reference_matrix(M, "M")
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
    if (!method %in% c("OLS", "NNLS", "WLS")) {
        stop("method must be 'OLS', 'NNLS', or 'WLS'")
    }

    af_match <- grepl("^AF($|_)", rownames(M), ignore.case = TRUE)
    variances <- attr(M, "variances")
    wls_detector_weights <- if (method == "WLS" && !is.null(variances)) {
        .wls_weights_from_variances(variances, n_detectors = ncol(M))
    } else {
        numeric(0)
    }
    if (method == "WLS" && length(wls_detector_weights) != ncol(M)) {
        stop(
            "WLS requires SCC-derived detector variances. ",
            "Run unmix_controls() first, pass M returned by build_reference_matrix(), ",
            "or use unmix_samples() with a matching variances_file."
        )
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
            detector_weights = wls_detector_weights
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
            w_scale <- sqrt(wls_detector_weights)
            M_scaled <- M %*% diag(w_scale)
            colnames(M_scaled) <- colnames(M)
            rownames(M_scaled) <- rownames(M)
            Y_scaled <- Y %*% diag(w_scale)
            colnames(Y_scaled) <- colnames(Y)
            Mt_scaled <- t(M_scaled)
            matrix_to_invert <- M_scaled %*% Mt_scaled
            if (rcond(matrix_to_invert) < 1e-10) {
                stop("Scaled Reference Matrix is singular.")
            }
            A <- Y_scaled %*% Mt_scaled %*% solve(matrix_to_invert)
        }
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
