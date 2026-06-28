#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(flowCore)
})

if (!requireNamespace("spectreasy", quietly = TRUE)) {
    stop("Install or load spectreasy before running this benchmark.", call. = FALSE)
}

set.seed(20260628)

marker_M <- rbind(
    F1 = c(1.00, 0.20, 0.00),
    F2 = c(0.00, 0.20, 1.00)
)
af_M <- rbind(
    AF = c(0.80, 1.00, 0.20),
    AF_2 = c(0.10, 1.00, 0.90)
)
colnames(marker_M) <- colnames(af_M) <- c("D1", "D2", "D3")
M <- rbind(marker_M, af_M)

simulate_events <- function(n = 2000, fluor_overlay = TRUE, noise_sd = 1) {
    true_af <- sample.int(nrow(af_M), n, replace = TRUE)
    af_signal <- runif(n, 40, 140)
    fluor_signal <- matrix(0, nrow = n, ncol = nrow(marker_M), dimnames = list(NULL, rownames(marker_M)))
    if (isTRUE(fluor_overlay)) {
        stained_idx <- sample.int(n, floor(n / 2))
        fluor_signal[stained_idx, "F1"] <- rgamma(length(stained_idx), shape = 1.4, rate = 0.05)
        fluor_signal[stained_idx, "F2"] <- rgamma(length(stained_idx), shape = 1.2, rate = 0.06)
    }
    Y <- fluor_signal %*% marker_M +
        af_signal * af_M[true_af, ] +
        matrix(rnorm(n * ncol(M), sd = noise_sd), ncol = ncol(M))
    colnames(Y) <- colnames(M)
    list(Y = Y, true_af = true_af, fluor_signal = fluor_signal)
}

cosine <- function(a, b) {
    sum(a * b) / sqrt(sum(a^2) * sum(b^2))
}

fit_by_assignment <- function(Y, assignments) {
    out <- matrix(0, nrow = nrow(Y), ncol = nrow(M), dimnames = list(NULL, rownames(M)))
    for (k in sort(unique(assignments))) {
        idx <- which(assignments == k)
        X <- rbind(marker_M, af_M[k, , drop = FALSE])
        coeff <- Y[idx, , drop = FALSE] %*% t(X) %*% solve(X %*% t(X))
        out[idx, rownames(marker_M)] <- coeff[, seq_len(nrow(marker_M)), drop = FALSE]
        out[idx, rownames(af_M)[[k]]] <- coeff[, ncol(coeff)]
    }
    out
}

score_method <- function(sim, method) {
    Y <- sim$Y
    if (identical(method, "projection")) {
        assignments <- spectreasy:::.assign_af_by_projection(Y, marker_M, af_M)
    } else if (identical(method, "residual_alignment")) {
        assignments <- spectreasy:::.assign_af_by_residual_alignment(Y, marker_M, af_M)
    } else if (identical(method, "legacy_residual")) {
        res <- spectreasy::calc_residuals(flowCore::flowFrame(Y), M, method = "OLS", af_assignment = "legacy", return_residuals = TRUE)
        assigned_names <- apply(as.matrix(res$data[, rownames(af_M), drop = FALSE]) != 0, 1, function(x) {
            hit <- which(x)
            if (length(hit)) hit[[1]] else 1L
        })
        assignments <- as.integer(assigned_names)
    } else if (identical(method, "random_baseline")) {
        assignments <- sample.int(nrow(af_M), nrow(Y), replace = TRUE)
    } else {
        stop("Unknown method: ", method, call. = FALSE)
    }

    coeff <- fit_by_assignment(Y, assignments)
    fitted <- coeff %*% M
    selected_cosine <- vapply(seq_along(assignments), function(i) {
        cosine(af_M[assignments[[i]], ], af_M[sim$true_af[[i]], ])
    }, numeric(1))

    fluor_est <- coeff[, rownames(marker_M), drop = FALSE]
    off_target <- fluor_est
    off_target[sim$fluor_signal > 0] <- 0
    unstained <- rowSums(sim$fluor_signal) == 0
    data.frame(
        method = method,
        af_accuracy = mean(assignments == sim$true_af),
        selected_af_cosine = mean(selected_cosine),
        fluor_mae = mean(abs(fluor_est - sim$fluor_signal)),
        false_fluor_unstained = mean(rowSums(abs(fluor_est[unstained, , drop = FALSE]))),
        residual_rmse = sqrt(mean((Y - fitted)^2)),
        off_target_correction = mean(abs(off_target)),
        stringsAsFactors = FALSE
    )
}

sim <- simulate_events()
methods <- c("projection", "residual_alignment", "legacy_residual", "random_baseline")
results <- do.call(rbind, lapply(methods, score_method, sim = sim))
print(results, row.names = FALSE)
