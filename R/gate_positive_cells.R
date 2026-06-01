#' Gate Positive Events from a Single-Stain Matrix
#'
#' Performs one-dimensional histogram gating on the channel with highest variance.
#' The function returns a logical vector indicating events inside the positive gate.
#'
#' @param mat Numeric matrix/data.frame of events x detectors.
#' @param histogram_pct Quantile width used for the histogram gate.
#' @param histogram_direction Gate direction: `"right"` starts at the median,
#'   `"both"` centers on the median, and `"left"` ends at the median.
#'
#' @return Logical vector of length `nrow(mat)`, `TRUE` for gated-in events.
#' @export
#' @examples
#' set.seed(1)
#' mat <- cbind(
#'   ch1 = c(rlnorm(500, 2, 0.4), rlnorm(500, 5, 0.3)),
#'   ch2 = rlnorm(1000, 2.5, 0.4)
#' )
#' idx <- gate_positive_cells(mat, histogram_pct = 0.4, histogram_direction = "right")
#' mean(idx)
gate_positive_cells <- function(mat,
                                histogram_pct = 0.98,
                                histogram_direction = "right") {
    # Identify peak channel by variance
    peak_channel <- which.max(apply(mat, 2, var))
    peak_vals <- mat[, peak_channel]
    vals_log <- log10(pmax(peak_vals, 1))

    # Calculate gate thresholds based on direction
    if (histogram_direction == "right") {
        lower_q <- 0.5
        upper_q <- 0.5 + histogram_pct
    } else if (histogram_direction == "left") {
        lower_q <- 0.5 - histogram_pct
        upper_q <- 0.5
    } else { # "both"
        lower_q <- 0.5 - histogram_pct / 2
        upper_q <- 0.5 + histogram_pct / 2
    }
    lower_q <- max(0, lower_q)
    upper_q <- min(1, upper_q)

    # Apply thresholds on log scale, return linear
    gate_min <- 10^quantile(vals_log, lower_q)
    gate_max <- 10^quantile(vals_log, upper_q)

    peak_vals >= gate_min & peak_vals <= gate_max
}
