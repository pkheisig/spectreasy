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

test_that("derive_unmixing_matrix requires variances for WLS and supports NNLS proxy", {
    M <- matrix(c(
        1, 0.2, 0.1,
        0.1, 1, 0.3
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B2-A", "YG1-A", "R1-A")

    expect_error(
        spectreasy::derive_unmixing_matrix(M, method = "WLS"),
        regexp = "SCC-derived detector variances"
    )

    V <- matrix(c(
        10, 20, 30,
        15, 25, 35
    ), nrow = 2, byrow = TRUE, dimnames = dimnames(M))
    W_wls <- spectreasy::derive_unmixing_matrix(M, method = "WLS", variances = V)
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

    V <- matrix(c(
        50, 10,
        5, 60
    ), nrow = 2, byrow = TRUE)
    rownames(V) <- rownames(M)
    colnames(V) <- colnames(M)
    attr(M, "variances") <- V

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

    V <- matrix(c(
        50, 10,
        5, 60
    ), nrow = 2, byrow = TRUE)
    rownames(V) <- rownames(M)
    colnames(V) <- colnames(M)
    attr(M, "variances") <- V

    exprs <- matrix(c(
        100,  20,
         10, 120,
         50,  40
    ), nrow = 3, byrow = TRUE)
    colnames(exprs) <- c("B1-A", "YG1-A")

    ff <- flowCore::flowFrame(exprs)
    res <- spectreasy::calc_residuals(ff, M, method = "WLS")

    Mt <- t(M)
    detector_weights <- spectreasy:::.wls_weights_from_variances(V, n_detectors = ncol(M))
    Wi <- diag(detector_weights)
    expected <- matrix(0, nrow = nrow(exprs), ncol = nrow(M))
    for (i in seq_len(nrow(exprs))) {
        expected[i, ] <- exprs[i, , drop = FALSE] %*% Wi %*% Mt %*% solve(M %*% Wi %*% Mt)
    }
    colnames(expected) <- rownames(M)

    expect_equal(as.matrix(res[, rownames(M)]), expected, tolerance = 1e-6)
})

test_that("calc_residuals WLS errors without SCC-derived variances", {
    M <- matrix(c(
        1, 0.2,
        0.1, 1
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B1-A", "YG1-A")

    exprs <- matrix(c(100, 20), nrow = 1)
    colnames(exprs) <- c("B1-A", "YG1-A")
    ff <- flowCore::flowFrame(exprs)

    expect_error(
        spectreasy::calc_residuals(ff, M, method = "WLS"),
        regexp = "SCC-derived detector variances"
    )
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

test_that(".compute_reference_spectrum computes and attaches variance attribute", {
    detector_names <- c("B1-A", "YG1-A")
    set.seed(42)
    final_gated_data <- matrix(rnorm(100, mean = 1000, sd = 50), ncol = 2)
    colnames(final_gated_data) <- detector_names
    gated_data <- matrix(rnorm(200, mean = 200, sd = 10), ncol = 2)
    colnames(gated_data) <- detector_names
    
    peak_vals <- gated_data[, 1]
    vals_log <- log10(pmax(peak_vals, 1))
    row_info <- data.frame(universal.negative = "FALSE", stringsAsFactors = FALSE)
    
    res <- spectreasy:::.compute_reference_spectrum(
        final_gated_data = final_gated_data,
        gated_data = gated_data,
        peak_vals = peak_vals,
        vals_log = vals_log,
        detector_names = detector_names,
        row_info = row_info
    )
    
    expect_type(res, "double")
    expect_equal(length(res), 2)
    expect_equal(names(res), detector_names)
    
    var_attr <- attr(res, "variance")
    expect_type(var_attr, "double")
    expect_equal(length(var_attr), 2)
    expect_true(all(var_attr >= 0))
})

test_that(".compute_reference_spectrum honors named universal negative files", {
    detector_names <- c("B1-A", "YG1-A")
    final_gated_data <- matrix(
        rep(c(100, 60), each = 20),
        ncol = 2,
        dimnames = list(NULL, detector_names)
    )
    gated_data <- matrix(
        rep(c(10, 10), each = 40),
        ncol = 2,
        dimnames = list(NULL, detector_names)
    )
    peak_vals <- gated_data[, 1]
    vals_log <- log10(pmax(peak_vals, 1))
    row_info <- data.frame(universal.negative = "unstained_cells.fcs", stringsAsFactors = FALSE)
    universal_negatives <- list(unstained_cells = c("B1-A" = 30, "YG1-A" = 10))

    res <- spectreasy:::.compute_reference_spectrum(
        final_gated_data = final_gated_data,
        gated_data = gated_data,
        peak_vals = peak_vals,
        vals_log = vals_log,
        detector_names = detector_names,
        row_info = row_info,
        universal_negatives = universal_negatives
    )

    expect_equal(as.numeric(res), c(1, 50 / 70), tolerance = 1e-6)
})

test_that("control.type from control file overrides filename fallback", {
    row_info <- data.frame(
        filename = "Ambiguous Control.fcs",
        fluorophore = "FITC",
        channel = "B1-A",
        control.type = "cells",
        stringsAsFactors = FALSE
    )

    sample_info <- spectreasy:::.resolve_reference_sample_type(
        filename = "Ambiguous Control",
        row_info = row_info,
        patterns = list(beads = "Beads", cells = "Cells"),
        default = "beads"
    )

    expect_equal(sample_info$type, "cells")
})

test_that("derive_unmixing_matrix uses variances attribute for WLS", {
    M <- matrix(c(
        1.0, 0.2, 0.05,
        0.1, 1.0, 0.25
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B1-A", "YG1-A", "R1-A")
    
    V <- matrix(c(
        100, 20, 5,
        15, 120, 50
    ), nrow = 2, byrow = TRUE)
    rownames(V) <- rownames(M)
    colnames(V) <- colnames(M)
    attr(M, "variances") <- V
    
    expect_silent(
        W_wls <- spectreasy::derive_unmixing_matrix(M, method = "WLS")
    )
    
    expect_equal(dim(W_wls), dim(M))

    W_wls_explicit <- spectreasy::derive_unmixing_matrix(M, method = "WLS", variances = V)
    expect_equal(W_wls, W_wls_explicit)
    
    M_no_attr <- M
    attr(M_no_attr, "variances") <- NULL
    expect_error(
        spectreasy::derive_unmixing_matrix(M_no_attr, method = "WLS"),
        regexp = "SCC-derived detector variances"
    )
})

test_that("calc_residuals WLS performs global scaling unmixing when variances attribute is present", {
    M <- matrix(c(
        1, 0.2,
        0.1, 1
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B1-A", "YG1-A")

    V <- matrix(c(
        50, 10,
        5, 60
    ), nrow = 2, byrow = TRUE)
    rownames(V) <- rownames(M)
    colnames(V) <- colnames(M)
    attr(M, "variances") <- V

    exprs <- matrix(c(
        100,  20,
         10, 120
    ), nrow = 2, byrow = TRUE)
    colnames(exprs) <- c("B1-A", "YG1-A")
    ff <- flowCore::flowFrame(exprs)

    res <- spectreasy::calc_residuals(ff, M, method = "WLS")
    expect_equal(dim(res), c(2, 2))
    expect_true(all(is.finite(as.matrix(res))))
})

test_that("calc_residuals multi-AF WLS uses SCC-derived detector weights", {
    M <- matrix(c(
        1.0, 0.2, 0.1, 0.05,
        0.1, 1.0, 0.2, 0.05,
        0.2, 0.2, 1.0, 0.05,
        0.05, 0.05, 0.05, 1.0
    ), nrow = 4, byrow = TRUE)
    rownames(M) <- c("FITC", "PE", "AF", "AF_2")
    colnames(M) <- c("B1-A", "YG1-A", "R1-A", "V1-A")

    V <- matrix(c(
        1, 100, 10, 50,
        2, 80, 12, 40,
        0, 0, 0, 0,
        0, 0, 0, 0
    ), nrow = 4, byrow = TRUE)
    rownames(V) <- rownames(M)
    colnames(V) <- colnames(M)
    attr(M, "variances") <- V

    exprs <- matrix(c(
        100, 10, 50, 0,
        10, 100, 0, 50
    ), nrow = 2, byrow = TRUE)
    colnames(exprs) <- colnames(M)
    ff <- flowCore::flowFrame(exprs)

    res <- spectreasy::calc_residuals(ff, M, method = "WLS")

    detector_weights <- spectreasy:::.wls_weights_from_variances(V, n_detectors = ncol(M))
    fluor_idx <- which(!grepl("^AF($|_)", rownames(M), ignore.case = TRUE))
    af_idx <- which(grepl("^AF($|_)", rownames(M), ignore.case = TRUE))
    expected <- matrix(0, nrow = nrow(exprs), ncol = nrow(M), dimnames = list(NULL, rownames(M)))

    for (i in seq_len(nrow(exprs))) {
        y <- exprs[i, ]
        best_rss <- Inf
        best_coeffs <- NULL
        best_af <- NA_integer_
        for (k in af_idx) {
            rows <- c(fluor_idx, k)
            X <- M[rows, , drop = FALSE]
            Xw <- sweep(X, 2, detector_weights, `*`)
            coeffs <- solve(Xw %*% t(X), X %*% (detector_weights * y))
            resid <- y - drop(as.numeric(coeffs) %*% X)
            rss <- sum(detector_weights * resid^2)
            if (rss < best_rss) {
                best_rss <- rss
                best_coeffs <- coeffs
                best_af <- k
            }
        }
        expected[i, fluor_idx] <- best_coeffs[seq_along(fluor_idx)]
        expected[i, best_af] <- best_coeffs[length(best_coeffs)]
    }

    expect_equal(as.matrix(res[, rownames(M)]), expected, tolerance = 1e-6)
})

test_that("unmix_samples integrates with variances CSV file", {
    M <- matrix(c(
        1, 0.2,
        0.1, 1
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B1-A", "YG1-A")

    V <- matrix(c(
        50, 10,
        5, 60
    ), nrow = 2, byrow = TRUE)
    rownames(V) <- rownames(M)
    colnames(V) <- colnames(M)
    
    tmp_ref <- tempfile("ref_", fileext = ".csv")
    tmp_var <- tempfile("var_", fileext = ".csv")
    
    ref_df <- as.data.frame(M)
    ref_df$Marker <- rownames(M)
    ref_df <- ref_df[, c("Marker", "B1-A", "YG1-A")]
    write.csv(ref_df, tmp_ref, row.names = FALSE)
    
    var_df <- as.data.frame(V)
    var_df$Marker <- rownames(V)
    var_df <- var_df[, c("Marker", "B1-A", "YG1-A")]
    write.csv(var_df, tmp_var, row.names = FALSE)
    
    exprs <- matrix(c(
        100,  20,
         10, 120
    ), nrow = 2, byrow = TRUE)
    colnames(exprs) <- c("B1-A", "YG1-A")
    ff <- flowCore::flowFrame(exprs)
    
    sample_dir <- tempfile("samples_")
    dir.create(sample_dir, showWarnings = FALSE)
    flowCore::write.FCS(ff, file.path(sample_dir, "sample.fcs"))
    
    output_dir <- tempfile("unmixed_")
    
    unmixed <- spectreasy::unmix_samples(
        sample_dir = sample_dir,
        M = NULL,
        unmixing_matrix_file = tmp_ref,
        variances_file = tmp_var,
        method = "WLS",
        output_dir = output_dir,
        write_fcs = TRUE
    )
    
    expect_true(file.exists(file.path(output_dir, "sample_unmixed.fcs")))
    expect_setequal(names(unmixed), "sample")
})
