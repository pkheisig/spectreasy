#' Derive Unmixing Matrix from Reference Matrix
#' 
#' Calculates a static unmixing matrix from a reference spectral matrix (M).
#' The output matrix (W) has the same dimensions as M (Markers x Detectors).
#' To unmix data manually: Unmixed_Data = Raw_Data %*% t(W)
#' 
#' @param M Reference matrix (Markers x Detectors)
#' @param method Unmixing method ("OLS", "WLS", or "NNLS")
#' @param variances Optional SCC-derived variance matrix with the same detector
#'   columns as `M`. If `NULL`, `attr(M, "variances")` is used.
#' @return A matrix of unmixing coefficients (Markers x Detectors)
#' @examples
#' M <- matrix(c(1, 0.2, 0.1, 1), nrow = 2, byrow = TRUE)
#' rownames(M) <- c("FITC", "PE")
#' colnames(M) <- c("B2-A", "YG1-A")
#' W <- derive_unmixing_matrix(M, method = "OLS")
#' W
#' @export
derive_unmixing_matrix <- function(M, method = "OLS", variances = NULL) {
    M <- .as_reference_matrix(M, "M")
    if (!is.null(variances)) {
        attr(M, "variances") <- .as_reference_matrix(variances, "variances")
    }
    # M is Markers (m) x Detectors (d)
    Mt <- t(M) # Detectors x Markers

    method_upper <- toupper(trimws(method))

    if (method_upper == "OLS") {
        # Standard analytical solution: W = (M %*% M^T)^-1 %*% M
        MMt <- M %*% Mt
        if (rcond(MMt) < 1e-10) stop("Reference Matrix is singular (collinear spectra).")
        W <- solve(MMt) %*% M

    } else if (method_upper == "WLS") {
        variances <- attr(M, "variances")
        if (is.null(variances)) {
            stop(
                "WLS requires SCC-derived detector variances. ",
                "Pass variances = <scc_variances matrix> or use an M returned by build_reference_matrix()."
            )
        }
        detector_weights <- .wls_weights_from_variances(variances, n_detectors = ncol(M))

        detector_weights <- as.numeric(detector_weights)
        if (length(detector_weights) != ncol(M)) {
            stop("WLS detector weights length must match number of detectors (", ncol(M), ").")
        }
        if (any(!is.finite(detector_weights)) || any(detector_weights <= 0)) {
            stop("WLS detector weights must be finite and > 0.")
        }

        V_inv <- diag(detector_weights)
        # W = (M V^-1 M^T)^-1 M V^-1
        MVMt <- M %*% V_inv %*% Mt
        if (rcond(MVMt) < 1e-10) stop("Weighted matrix is singular.")
        W <- solve(MVMt) %*% M %*% V_inv

    } else if (method_upper == "NNLS") {
        # A single static matrix cannot exactly reproduce per-cell NNLS
        # (piecewise-linear constraint). We export a deterministic linear proxy
        # by solving NNLS for each detector basis vector.
        d <- ncol(M)
        I_det <- diag(d)
        W <- t(spectreasy_nnls_unmix_cpp(Y = I_det, M = M))

        warning(
            "Static NNLS matrix is a linear proxy and may differ from per-cell NNLS solutions. ",
            "Use dynamic method = 'NNLS' in unmix_samples() for exact constrained per-cell fits."
        )

    } else {
        stop("method must be one of: 'OLS', 'WLS', 'NNLS'")
    }

    rownames(W) <- rownames(M)
    colnames(W) <- colnames(M)

    return(W)
}

#' Save Unmixing Matrix to CSV
#' @param W Unmixing matrix
#' @param file Path to save
#' @return Invisibly returns `NULL`; writes CSV to disk.
#' @examples
#' tmp <- tempfile(fileext = ".csv")
#' M <- matrix(c(1, 0.2, 0.1, 1), nrow = 2, byrow = TRUE)
#' rownames(M) <- c("FITC", "PE")
#' colnames(M) <- c("B2-A", "YG1-A")
#' save_unmixing_matrix(M, tmp)
#' @export
save_unmixing_matrix <- function(W, file = "unmixing_matrix.csv") {
    W_df <- as.data.frame(W)
    W_df$Marker <- rownames(W)
    W_df <- W_df[, c("Marker", setdiff(colnames(W_df), "Marker"))]
    utils::write.csv(W_df, file, row.names = FALSE, quote = TRUE)
    message("Unmixing matrix (", nrow(W), "x", ncol(W), ") saved to: ", file)
}

#' Plot Unmixing Matrix
#' @param W Unmixing matrix
#' @param pd Optional pData for descriptive labels
#' @return ggplot object
#' @examples
#' M <- matrix(c(1, 0.2, 0.1, 1), nrow = 2, byrow = TRUE)
#' rownames(M) <- c("FITC", "PE")
#' colnames(M) <- c("B2-A", "YG1-A")
#' p <- plot_unmixing_matrix(M)
#' print(p)
#' @export
plot_unmixing_matrix <- function(W, pd = NULL) {
    long <- as.data.frame(W)
    long$Marker <- rownames(W)
    long <- tidyr::pivot_longer(long, cols = -Marker, names_to = "Detector", values_to = "Coefficient")
    
    # Sort and label detectors
    det_names <- colnames(W)
    if (!is.null(pd)) {
        det_info <- get_sorted_detectors(pd)
        common <- intersect(det_info$names, det_names)
        levels_sorted <- common
        labels_sorted <- det_info$labels[match(common, det_info$names)]
    } else {
        nums <- as.numeric(gsub("[^0-9]", "", det_names))
        levels_sorted <- det_names[order(nums)]
        labels_sorted <- levels_sorted
    }
    
    long$Detector <- factor(long$Detector, levels = levels_sorted, labels = labels_sorted)
    
    # Cap values for color scale to avoid outliers dominating
    cap_val <- quantile(abs(long$Coefficient), 0.95, na.rm=TRUE)
    if (is.na(cap_val) || cap_val == 0) cap_val <- 1
    long$Color_Val <- pmin(pmax(long$Coefficient, -cap_val), cap_val)
    
    p <- ggplot2::ggplot(long, ggplot2::aes(Detector, Marker, fill = Color_Val)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", name = "Coefficient") +
        # Text color based on value
        ggplot2::geom_text(ggplot2::aes(label = round(Coefficient, 3), 
                                       color = abs(Color_Val) > cap_val/2), 
                           size = 1.5, show.legend = FALSE, angle = 90) +
        ggplot2::scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black")) +
        ggplot2::labs(title = "Unmixing Matrix Coefficients",
                      subtitle = "Good: dominant, stable coefficients align with expected detector-markers. Bad: widespread large-magnitude off-target coefficients suggest collinearity, poor controls, or unstable inversion.") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 5),
            plot.subtitle = ggplot2::element_text(size = 10.6, lineheight = 1.1)
        )
    
    return(p)
}
