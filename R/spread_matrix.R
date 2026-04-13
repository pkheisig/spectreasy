#' Calculate Spectral Spread Matrix (SSM)
#' 
#' Quantifies the spreading error introduced by unmixing between fluorophores.
#' High values indicate that marker i significantly increases the noise in marker j.
#' 
#' @param M Reference matrix (Markers x Detectors)
#' @param method Unmixing method ("OLS" or "WLS")
#' @param background_noise Background noise used for WLS
#' @return A matrix (Markers x Markers) representing unmixing spread
#' @examples
#' M <- matrix(c(1, 0.2, 0.1, 1), nrow = 2, byrow = TRUE)
#' rownames(M) <- c("FITC", "PE")
#' colnames(M) <- c("B2-A", "YG1-A")
#' ssm <- calculate_ssm(M)
#' ssm
#' @export
calculate_ssm <- function(M, method = "OLS", background_noise = 100) {
    M <- .as_reference_matrix(M, "M")
    # This function estimates unmixing-induced spread analytically.
    # Spreading error is proportional to the square root of signal intensity.
    # The analytical estimate for spread in unmixed space is:
    # Var(Unmixed) = W %*% diag(Poisson_Variance) %*% t(W)
    
    # M: Markers x Detectors
    n_markers <- nrow(M)
    marker_names <- rownames(M)
    
    # Derive unmixing matrix W (Markers x Detectors)
    W <- derive_unmixing_matrix(M, method = method)
    
    # Initialize SSM
    SSM <- matrix(0, nrow = n_markers, ncol = n_markers)
    rownames(SSM) <- marker_names
    colnames(SSM) <- marker_names
    
    # For each marker i (the "spilling" marker)
    for (i in seq_len(n_markers)) {
        # Assume a unit signal intensity for marker i
        # The expected photon counts in each detector for marker i is M[i, ]
        signal_i <- M[i, ] 
        
        # Variance in each detector due to marker i (Poisson)
        # We assume signal is large enough that Shot Noise dominates
        # Var_detectors = signal_i
        var_d <- diag(signal_i)
        
        # Propagation of variance to unmixed space: W %*% Var_d %*% W^T
        # Result is n_markers x n_markers
        unmixed_variance <- W %*% var_d %*% t(W)
        
        # The unmixed standard deviation (spread) in all other markers j
        unmixed_sd <- sqrt(pmax(diag(unmixed_variance), 0))
        
        # SSM[i, j] = SD_j / sqrt(Signal_i)
        SSM[i, ] <- unmixed_sd
    }
    
    # Set diagonal to 0 as we only care about spread into OTHER channels
    diag(SSM) <- 0
    
    return(SSM)
}

#' Plot Spectral Spread Matrix
#' @param SSM Matrix returned by calculate_ssm
#' @param output_file Path to save the plot
#' @param width Width of plot in mm
#' @param height Height of plot in mm
#' @return A `ggplot` object.
#' @examples
#' M <- matrix(c(1, 0.2, 0.1, 1), nrow = 2, byrow = TRUE)
#' rownames(M) <- c("FITC", "PE")
#' colnames(M) <- c("B2-A", "YG1-A")
#' ssm <- calculate_ssm(M)
#' p <- plot_ssm(ssm, output_file = NULL)
#' print(p)
#' @export
plot_ssm <- function(SSM, output_file = "spectral_spread_matrix.png", width = 200, height = 180) {
    long <- as.data.frame(SSM)
    long$Spilling_Marker <- rownames(SSM)
    long <- tidyr::pivot_longer(long, cols = -Spilling_Marker, names_to = "Receiving_Marker", values_to = "Spread")

    fmt_spread <- function(x) {
        if (!is.finite(x)) return("")
        ax <- abs(x)
        if (ax >= 1000) return(formatC(x, format = "e", digits = 1))
        if (ax >= 100) return(formatC(x, format = "f", digits = 0))
        if (ax >= 10) return(formatC(x, format = "f", digits = 1))
        formatC(x, format = "f", digits = 2)
    }

    n_markers <- max(nrow(SSM), ncol(SSM))
    text_size <- max(1.6, min(3.2, 24 / max(1, n_markers)))
    spread_cut <- stats::quantile(long$Spread, 0.7, na.rm = TRUE)

    p <- ggplot2::ggplot(long, ggplot2::aes(Receiving_Marker, Spilling_Marker, fill = Spread)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_viridis_c(option = "magma", name = "Spread Factor") +
        ggplot2::geom_text(
            ggplot2::aes(label = vapply(Spread, fmt_spread, character(1)), color = Spread > spread_cut),
            size = text_size,
            show.legend = FALSE
        ) +
        ggplot2::scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white")) +
        ggplot2::labs(title = "Spectral Spread Matrix",
                      subtitle = "Rows = noise source, columns = noise destination. Good: low off-diagonal spread (darker cells). Bad: bright off-diagonal cells mark problematic marker pairs for dim co-expression.",
                      x = "Receiving Marker (Noise Destination)", y = "Spilling Marker (Noise Source)") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, p, width = width, height = height, units = "mm", dpi = 300)
    }
    return(p)
}
