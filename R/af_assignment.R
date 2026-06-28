.normalize_af_assignment <- function(af_assignment,
                                     choices = c("projection", "residual_alignment", "legacy", "legacy_residual")) {
    af_assignment <- as.character(af_assignment[1])
    if (!nzchar(af_assignment) || is.na(af_assignment)) {
        af_assignment <- "projection"
    }
    match.arg(af_assignment, choices = choices)
}

.reference_ols_unmixing_matrix <- function(spectra, ridge = 1e-8) {
    spectra <- as.matrix(spectra)
    cross <- spectra %*% t(spectra)
    out <- tryCatch(
        qr.solve(cross, spectra),
        error = function(e) NULL
    )
    if (!is.null(out) && all(is.finite(out))) {
        return(out)
    }

    scale <- mean(diag(cross), na.rm = TRUE)
    if (!is.finite(scale) || scale <= 0) scale <- 1
    qr.solve(cross + diag(ridge * scale, nrow(cross)), spectra)
}

.project_af_into_fluor_space <- function(marker_M, af_M, eps = 1e-10) {
    marker_M <- as.matrix(marker_M)
    af_M <- as.matrix(af_M)
    if (nrow(marker_M) == 0) {
        stop("Projection AF assignment requires at least one non-AF fluorophore row.", call. = FALSE)
    }
    if (nrow(af_M) == 0) {
        stop("Projection AF assignment requires at least one AF candidate row.", call. = FALSE)
    }
    if (ncol(marker_M) != ncol(af_M) || !identical(colnames(marker_M), colnames(af_M))) {
        stop("marker_M and af_M must use the same detector columns.", call. = FALSE)
    }

    unmixing_matrix <- .reference_ols_unmixing_matrix(marker_M)
    v_library <- unmixing_matrix %*% t(af_M)
    r_library <- t(af_M) - (t(marker_M) %*% v_library)
    denominator <- colSums(r_library^2)
    denominator[!is.finite(denominator) | denominator <= eps] <- eps

    rownames(unmixing_matrix) <- rownames(marker_M)
    colnames(unmixing_matrix) <- colnames(marker_M)
    rownames(v_library) <- rownames(marker_M)
    colnames(v_library) <- rownames(af_M)
    rownames(r_library) <- colnames(marker_M)
    colnames(r_library) <- rownames(af_M)

    list(
        unmixing_matrix = unmixing_matrix,
        v_library = v_library,
        r_library = r_library,
        denominator = denominator
    )
}

.estimate_af_abundance_from_residual_component <- function(Y, r_library, denominator) {
    Y <- as.matrix(Y)
    r_library <- as.matrix(r_library)
    if (ncol(Y) != nrow(r_library) || !identical(colnames(Y), rownames(r_library))) {
        stop("Y and r_library must use matching detector columns.", call. = FALSE)
    }
    k_matrix <- sweep(Y %*% r_library, 2, denominator, "/")
    colnames(k_matrix) <- colnames(r_library)
    k_matrix
}

.selected_af_details <- function(assignments, k_matrix, v_library, score_matrix = NULL,
                                 r_library = NULL, return_details = FALSE) {
    assignments <- as.integer(assignments)
    n <- length(assignments)
    selected <- cbind(seq_len(n), assignments)
    af_abundance <- k_matrix[selected]
    projected_af_fluor <- matrix(0, nrow = n, ncol = nrow(v_library))
    colnames(projected_af_fluor) <- rownames(v_library)
    for (k in sort(unique(assignments))) {
        idx <- which(assignments == k)
        if (!length(idx) || is.na(k) || k < 1L || k > ncol(v_library)) next
        projected_af_fluor[idx, ] <- k_matrix[idx, k, drop = FALSE] %*% t(v_library[, k, drop = FALSE])
    }

    out <- list(
        assignments = assignments,
        af_abundance = as.numeric(af_abundance),
        projected_af_fluor = projected_af_fluor
    )
    if (isTRUE(return_details)) {
        out$score_matrix <- score_matrix
        out$v_library <- v_library
        out$r_library <- r_library
        out$af_abundance_matrix <- k_matrix
    }
    out
}

.assign_af_by_projection <- function(Y,
                                     marker_M,
                                     af_M,
                                     return_details = FALSE,
                                     eps = 1e-10) {
    Y <- as.matrix(Y)
    marker_M <- as.matrix(marker_M)
    af_M <- as.matrix(af_M)
    if (!isTRUE(return_details)) {
        return(as.integer(spectreasy_assign_af_projection_cpp(
            Y = Y,
            F = marker_M,
            AF = af_M,
            tol = eps
        )))
    }
    projection <- .project_af_into_fluor_space(marker_M, af_M, eps = eps)
    k_matrix <- .estimate_af_abundance_from_residual_component(
        Y = Y,
        r_library = projection$r_library,
        denominator = projection$denominator
    )
    f0 <- Y %*% t(projection$unmixing_matrix)
    colnames(f0) <- rownames(marker_M)

    score_matrix <- matrix(0, nrow = nrow(Y), ncol = nrow(af_M), dimnames = list(NULL, rownames(af_M)))
    for (k in seq_len(nrow(af_M))) {
        adjusted <- f0 - (k_matrix[, k, drop = FALSE] %*% t(projection$v_library[, k, drop = FALSE]))
        score_matrix[, k] <- rowSums(abs(adjusted))
    }
    assignments <- max.col(-score_matrix, ties.method = "first")

    details <- .selected_af_details(
        assignments = assignments,
        k_matrix = k_matrix,
        v_library = projection$v_library,
        score_matrix = score_matrix,
        r_library = projection$r_library,
        return_details = return_details
    )
    if (isTRUE(return_details)) details else details$assignments
}

.assign_af_by_residual_alignment <- function(Y,
                                             marker_M,
                                             af_M,
                                             return_details = FALSE,
                                             eps = 1e-10) {
    Y <- as.matrix(Y)
    marker_M <- as.matrix(marker_M)
    af_M <- as.matrix(af_M)
    projection <- .project_af_into_fluor_space(marker_M, af_M, eps = eps)
    f0 <- Y %*% t(projection$unmixing_matrix)
    residual <- Y - (f0 %*% marker_M)
    score_matrix <- residual %*% t(af_M)
    colnames(score_matrix) <- rownames(af_M)
    assignments <- max.col(score_matrix, ties.method = "first")
    k_matrix <- .estimate_af_abundance_from_residual_component(
        Y = Y,
        r_library = projection$r_library,
        denominator = projection$denominator
    )

    details <- .selected_af_details(
        assignments = assignments,
        k_matrix = k_matrix,
        v_library = projection$v_library,
        score_matrix = score_matrix,
        r_library = projection$r_library,
        return_details = return_details
    )
    if (isTRUE(return_details)) details else details$assignments
}

.assign_af_candidates <- function(Y,
                                  marker_M,
                                  af_M,
                                  af_assignment = c("projection", "residual_alignment"),
                                  return_details = FALSE,
                                  eps = 1e-10) {
    af_assignment <- .normalize_af_assignment(af_assignment, choices = c("projection", "residual_alignment"))
    if (nrow(af_M) <= 1L) {
        if (!isTRUE(return_details)) {
            return(rep(1L, nrow(Y)))
        }
    }
    if (identical(af_assignment, "projection")) {
        return(.assign_af_by_projection(Y, marker_M, af_M, return_details = return_details, eps = eps))
    }
    .assign_af_by_residual_alignment(Y, marker_M, af_M, return_details = return_details, eps = eps)
}
