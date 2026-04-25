test_that("gating_options returns named list", {
    opts <- spectreasy::gating_options(histogram_pct_beads = 0.9, histogram_pct_cells = 0.3)
    expect_type(opts, "list")
    expect_true(all(c(
        "histogram_pct_beads",
        "histogram_direction_beads",
        "histogram_pct_cells",
        "histogram_direction_cells"
    ) %in% names(opts)))
    expect_equal(opts$histogram_pct_beads, 0.9)
    expect_equal(opts$histogram_pct_cells, 0.3)
})

test_that("derive_unmixing_matrix returns finite matrix with expected dims", {
    M <- matrix(c(
        1, 0.2, 0.1,
        0.1, 1, 0.3
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B2-A", "YG1-A", "R1-A")

    W <- spectreasy::derive_unmixing_matrix(M, method = "OLS")
    expect_equal(dim(W), dim(M))
    expect_true(all(is.finite(W)))
})

test_that("derive_unmixing_matrix supports WLS fallback and NNLS proxy", {
    M <- matrix(c(
        1, 0.2, 0.1,
        0.1, 1, 0.3
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B2-A", "YG1-A", "R1-A")

    expect_warning(
        W_wls <- spectreasy::derive_unmixing_matrix(M, method = "WLS"),
        regexp = "global_weights"
    )
    expect_equal(dim(W_wls), dim(M))
    expect_true(all(is.finite(W_wls)))

    expect_warning(
        W_nnls <- spectreasy::derive_unmixing_matrix(M, method = "NNLS"),
        regexp = "linear proxy"
    )
    expect_equal(dim(W_nnls), dim(M))
    expect_true(all(is.finite(W_nnls)))
    expect_true(all(W_nnls >= -1e-12))
})

test_that("calc_residuals NNLS matches per-cell constrained reference fits", {
    M <- matrix(c(
        1, 0.2,
        0.1, 1
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B1-A", "YG1-A")

    exprs <- matrix(c(
        100,  20,
         10, 120,
         50,  40
    ), nrow = 3, byrow = TRUE)
    colnames(exprs) <- c("B1-A", "YG1-A")

    ff <- flowCore::flowFrame(exprs)
    res <- spectreasy::calc_residuals(ff, M, method = "NNLS")

    nnls_reference <- function(A, b, tol = 1e-10, max_outer = 500, max_inner = 500) {
        p <- ncol(A)
        x <- numeric(p)
        passive <- rep(FALSE, p)
        w <- drop(crossprod(A, b - A %*% x))
        outer_iter <- 0

        while (outer_iter < max_outer) {
            active_idx <- which(!passive & w > tol)
            if (length(active_idx) == 0) break

            t_idx <- active_idx[which.max(w[active_idx])]
            passive[t_idx] <- TRUE

            inner_iter <- 0
            repeat {
                s <- numeric(p)
                p_idx <- which(passive)
                if (length(p_idx) > 0) {
                    AP <- A[, p_idx, drop = FALSE]
                    normal_mat <- crossprod(AP)
                    rhs <- crossprod(AP, b)
                    z <- tryCatch(
                        solve(normal_mat, rhs),
                        error = function(e) {
                            sv <- svd(normal_mat)
                            keep <- sv$d > (max(sv$d) * 1e-10)
                            if (!any(keep)) {
                                return(matrix(0, nrow = ncol(normal_mat), ncol = 1))
                            }
                            sv$v[, keep, drop = FALSE] %*%
                                (diag(1 / sv$d[keep], nrow = sum(keep)) %*% t(sv$u[, keep, drop = FALSE]) %*% rhs)
                        }
                    )
                    s[p_idx] <- drop(z)
                }

                nonpos_idx <- which(passive & s <= tol)
                if (length(nonpos_idx) == 0 || inner_iter >= max_inner) {
                    x <- s
                    break
                }

                alpha <- 1
                for (idx in nonpos_idx) {
                    denom <- x[idx] - s[idx]
                    if (denom > 0) {
                        alpha <- min(alpha, x[idx] / denom)
                    }
                }

                x <- x + alpha * (s - x)
                reset_idx <- which(passive & x <= tol)
                x[reset_idx] <- 0
                passive[reset_idx] <- FALSE
                inner_iter <- inner_iter + 1
            }

            w <- drop(crossprod(A, b - A %*% x))
            outer_iter <- outer_iter + 1
        }

        pmax(x, 0)
    }

    expected <- t(apply(exprs, 1, function(y) {
        nnls_reference(t(M), y)
    }))
    colnames(expected) <- rownames(M)

    expect_equal(as.matrix(res[, rownames(M)]), expected, tolerance = 1e-6)
    expect_true(all(as.matrix(res[, rownames(M)]) >= -1e-10))
})

test_that("calc_residuals WLS matches small-matrix reference implementation", {
    M <- matrix(c(
        1, 0.2,
        0.1, 1
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B1-A", "YG1-A")

    exprs <- matrix(c(
        100,  20,
         10, 120,
         50,  40
    ), nrow = 3, byrow = TRUE)
    colnames(exprs) <- c("B1-A", "YG1-A")

    ff <- flowCore::flowFrame(exprs)
    res <- spectreasy::calc_residuals(ff, M, method = "WLS", background_noise = 25)

    safe_inverse <- function(mat, tol = 1e-10) {
        if (rcond(mat) > tol) {
            return(solve(mat))
        }
        sv <- svd(mat)
        if (length(sv$d) == 0 || max(sv$d) == 0) {
            return(matrix(0, nrow = nrow(mat), ncol = ncol(mat)))
        }
        keep <- sv$d > (max(sv$d) * tol)
        if (!any(keep)) {
            return(matrix(0, nrow = nrow(mat), ncol = ncol(mat)))
        }
        sv$v[, keep, drop = FALSE] %*%
            (diag(1 / sv$d[keep], nrow = sum(keep)) %*% t(sv$u[, keep, drop = FALSE]))
    }

    Mt <- t(M)
    MMt_inv <- safe_inverse(M %*% Mt)
    expected <- matrix(0, nrow = nrow(exprs), ncol = nrow(M))
    for (i in seq_len(nrow(exprs))) {
        weights_i <- 1 / (pmax(exprs[i, ], 0) + 25)
        Wi <- diag(weights_i)
        MWMt_i <- M %*% Wi %*% Mt
        if (rcond(MWMt_i) > 1e-10) {
            expected[i, ] <- exprs[i, , drop = FALSE] %*% Wi %*% Mt %*% solve(MWMt_i)
        } else {
            expected[i, ] <- exprs[i, , drop = FALSE] %*% Mt %*% MMt_inv
        }
    }
    colnames(expected) <- rownames(M)

    expect_equal(as.matrix(res[, rownames(M)]), expected, tolerance = 1e-6)
})

test_that("calc_residuals retains Time and all FSC/SSC parameters but not raw detectors", {
    M <- matrix(c(
        1, 0.2,
        0.1, 1
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B1-A", "YG1-A")

    exprs <- matrix(c(
        100,  20,  1, 1000, 1200, 500, 50,
         90,  30,  2, 1100, 1300, 550, 55,
         80,  40,  3,  900, 1250, 530, 60
    ), nrow = 3, byrow = TRUE)
    colnames(exprs) <- c("B1-A", "YG1-A", "Time", "FSC-A", "FSC-H", "SSC-A", "SSC-W")

    ff <- flowCore::flowFrame(exprs)
    res <- spectreasy::calc_residuals(ff, M, method = "OLS")

    expect_setequal(colnames(res), c("FITC", "PE", "Time", "FSC-A", "FSC-H", "SSC-A", "SSC-W"))
    expect_false(any(c("B1-A", "YG1-A") %in% colnames(res)))
})

test_that("unmix_samples writes unmixed FCS with passthrough acquisition parameters only", {
    M <- matrix(c(
        1, 0.2,
        0.1, 1
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B1-A", "YG1-A")

    exprs <- matrix(c(
        100,  20,  1, 1000, 1200, 500, 50,
         90,  30,  2, 1100, 1300, 550, 55,
         80,  40,  3,  900, 1250, 530, 60
    ), nrow = 3, byrow = TRUE)
    colnames(exprs) <- c("B1-A", "YG1-A", "Time", "FSC-A", "FSC-H", "SSC-A", "SSC-W")

    ff <- flowCore::flowFrame(exprs)
    sample_dir <- tempfile("spectreasy_test_samples_")
    output_dir <- tempfile("spectreasy_test_unmixed_")
    dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    sample_file <- file.path(sample_dir, "sample1.fcs")
    flowCore::write.FCS(ff, sample_file)

    unmixed <- spectreasy::unmix_samples(
        sample_dir = sample_dir,
        M = M,
        method = "OLS",
        output_dir = output_dir,
        write_fcs = TRUE
    )
    expect_setequal(colnames(unmixed$sample1$data), c("FITC", "PE", "Time", "FSC-A", "FSC-H", "SSC-A", "SSC-W", "File"))

    unmixed_ff <- flowCore::read.FCS(
        file.path(output_dir, "sample1_unmixed.fcs"),
        transformation = FALSE,
        truncate_max_range = FALSE
    )
    out_cols <- colnames(flowCore::exprs(unmixed_ff))

    expect_true(all(c("FITC", "PE", "Time", "FSC-A", "FSC-H", "SSC-A", "SSC-W") %in% out_cols))
    expect_false(any(c("B1-A", "YG1-A") %in% out_cols))
})

test_that("unmix_samples uses safe output filenames and supports flowSet return", {
    M <- matrix(c(
        1, 0.2,
        0.1, 1
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B1-A", "YG1-A")

    exprs <- matrix(c(
        100, 20, 1, 1000, 500,
         90, 30, 2, 1100, 550,
         80, 40, 3,  900, 530
    ), nrow = 3, byrow = TRUE)
    colnames(exprs) <- c("B1-A", "YG1-A", "Time", "FSC-A", "SSC-A")

    ff <- flowCore::flowFrame(exprs)
    fs <- flowCore::flowSet(list(sample1 = ff, sample2 = ff))
    output_dir <- tempfile("spectreasy_test_safe_unmixed_")
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    flowCore::write.FCS(ff, file.path(output_dir, "sample1_unmixed.fcs"))

    unmixed_fs <- spectreasy::unmix_samples(
        sample_dir = fs,
        M = M,
        method = "OLS",
        output_dir = output_dir,
        write_fcs = TRUE,
        return_type = "flowSet"
    )

    expect_s4_class(unmixed_fs, "flowSet")
    expect_setequal(flowCore::sampleNames(unmixed_fs), c("sample1", "sample2"))
    expect_true(file.exists(file.path(output_dir, "sample1_unmixed.fcs")))
    expect_true(file.exists(file.path(output_dir, "sample1_unmixed_2.fcs")))
    expect_true(file.exists(file.path(output_dir, "sample2_unmixed.fcs")))

    residuals_attr <- attr(unmixed_fs, "spectreasy_residuals")
    expect_true(is.list(residuals_attr))
    expect_setequal(names(residuals_attr), c("sample1", "sample2"))
})

test_that("unmix_samples supports SingleCellExperiment input and output", {
    skip_if_not_installed("SingleCellExperiment")
    skip_if_not_installed("SummarizedExperiment")

    M <- matrix(c(
        1, 0.2,
        0.1, 1
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B1-A", "YG1-A")

    sample_a <- t(matrix(c(
        100, 20, 1, 1000, 500,
         90, 30, 2, 1100, 550,
         80, 40, 3,  900, 530
    ), nrow = 3, byrow = TRUE))
    sample_b <- t(matrix(c(
        120, 15, 4, 1050, 520,
        110, 25, 5, 1120, 540,
         95, 35, 6,  980, 510
    ), nrow = 3, byrow = TRUE))
    rownames(sample_a) <- rownames(sample_b) <- c("B1-A", "YG1-A", "Time", "FSC-A", "SSC-A")
    toy_sce <- suppressWarnings(
        SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = cbind(sample_a, sample_b)),
            colData = S4Vectors::DataFrame(
                sample_id = rep(c("sampleA", "sampleB"), each = ncol(sample_a)),
                row.names = paste0("cell_", seq_len(ncol(sample_a) + ncol(sample_b)))
            )
        )
    )

    unmixed_sce <- suppressWarnings(
        spectreasy::unmix_samples(
            sample_dir = toy_sce,
            M = M,
            method = "OLS",
            write_fcs = FALSE,
            return_type = "SingleCellExperiment"
        )
    )

    expect_s4_class(unmixed_sce, "SingleCellExperiment")
    expect_true("unmixed" %in% SummarizedExperiment::assayNames(unmixed_sce))
    expect_true("detector_residuals" %in% SingleCellExperiment::altExpNames(unmixed_sce))
    expect_equal(ncol(unmixed_sce), ncol(toy_sce))
    expect_setequal(unique(as.character(SummarizedExperiment::colData(unmixed_sce)$sample_id)), c("sampleA", "sampleB"))
    expect_true(all(c("FITC", "PE", "Time", "FSC-A", "SSC-A") %in% rownames(unmixed_sce)))

    residual_alt <- SingleCellExperiment::altExp(unmixed_sce, "detector_residuals")
    expect_setequal(rownames(residual_alt), c("B1-A", "YG1-A"))
})
