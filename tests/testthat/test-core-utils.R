test_that("gating_options returns named list", {
    opts <- spectreasy::gating_options(histogram_pct_beads = 0.9, histogram_pct_cells = 0.3)
    expect_type(opts, "list")
    expect_true(all(c(
        "use_scatter_gating",
        "histogram_pct_beads",
        "histogram_direction_beads",
        "histogram_pct_cells",
        "histogram_direction_cells"
    ) %in% names(opts)))
    expect_true(opts$use_scatter_gating)
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

test_that("reference matrices reject non-finite values and duplicated names", {
    M_bad_value <- matrix(c(1, NA_real_, 0.1, 1), nrow = 2, byrow = TRUE)
    rownames(M_bad_value) <- c("FITC", "PE")
    colnames(M_bad_value) <- c("B2-A", "YG1-A")

    expect_error(
        spectreasy::derive_unmixing_matrix(M_bad_value, method = "OLS"),
        regexp = "non-finite"
    )

    M_dup <- matrix(c(1, 0.2, 0.1, 1), nrow = 2, byrow = TRUE)
    rownames(M_dup) <- c("FITC", "FITC")
    colnames(M_dup) <- c("B2-A", "YG1-A")

    expect_error(
        spectreasy::derive_unmixing_matrix(M_dup, method = "OLS"),
        regexp = "duplicated marker"
    )
})

test_that("derive_unmixing_matrix supports static WLS approximation and NNLS proxy", {
    M <- matrix(c(
        1, 0.2, 0.1,
        0.1, 1, 0.3
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B2-A", "YG1-A", "R1-A")

    V <- matrix(c(
        10, 20, 30,
        15, 25, 35
    ), nrow = 2, byrow = TRUE, dimnames = dimnames(M))
    W_wls_default <- spectreasy::derive_unmixing_matrix(M, method = "WLS")
    expect_equal(dim(W_wls_default), dim(M))
    expect_true(all(is.finite(W_wls_default)))

    W_wls <- spectreasy::derive_unmixing_matrix(M, method = "WLS", variances = V)
    expect_equal(dim(W_wls), dim(M))
    expect_true(all(is.finite(W_wls)))
    expect_equal(W_wls, W_wls_default)

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

test_that("calc_residuals multi-AF NNLS selects AF band by constrained residual", {
    M <- matrix(c(
        1.00000000, 0.46357026, 0.56448788,
        0.20952040, 1.00000000, 0.41805212,
        1.00000000, 0.31610837, 0.64491387
    ), nrow = 3, byrow = TRUE)
    rownames(M) <- c("FITC", "AF", "AF_2")
    colnames(M) <- c("D1-A", "D2-A", "D3-A")

    exprs <- matrix(c(1.19896866, 1.83872848, 0.58991772), nrow = 1)
    colnames(exprs) <- colnames(M)
    ff <- flowCore::flowFrame(exprs)

    res <- spectreasy::calc_residuals(ff, M, method = "NNLS", return_residuals = TRUE)

    expect_gt(res$data$AF, 0)
    expect_equal(res$data$AF_2, 0)
    expect_lt(sum(res$residuals^2), 0.21)
})

test_that("calc_residuals multi-AF OLS fails clearly when all candidates are ill-conditioned", {
    M <- matrix(c(
        1, 0,
        2, 0,
        3, 0
    ), nrow = 3, byrow = TRUE)
    rownames(M) <- c("FITC", "AF", "AF_2")
    colnames(M) <- c("D1-A", "D2-A")

    exprs <- matrix(c(10, 0), nrow = 1)
    colnames(exprs) <- colnames(M)
    ff <- flowCore::flowFrame(exprs)

    expect_error(
        spectreasy::calc_residuals(ff, M, method = "OLS"),
        regexp = "No usable AF candidate model"
    )
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
    res <- spectreasy::calc_residuals(ff, M, method = "WLS")

    Mt <- t(M)
    expected <- matrix(0, nrow = nrow(exprs), ncol = nrow(M))
    for (i in seq_len(nrow(exprs))) {
        detector_weights <- spectreasy:::.wls_event_weights(
            exprs[i, ],
            noise_floor = rep(spectreasy:::.default_wls_background_noise(), ncol(M)),
            signal_scale = rep(1, ncol(M)),
            max_weight_ratio = spectreasy:::.default_wls_max_weight_ratio()
        )
        Wi <- diag(detector_weights)
        expected[i, ] <- exprs[i, , drop = FALSE] %*% Wi %*% Mt %*% solve(M %*% Wi %*% Mt)
    }
    colnames(expected) <- rownames(M)

    expect_equal(as.matrix(res[, rownames(M)]), expected, tolerance = 1e-6)
})

test_that("calc_residuals WLS works without SCC-derived variances", {
    M <- matrix(c(
        1, 0.2,
        0.1, 1
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B1-A", "YG1-A")

    exprs <- matrix(c(100, 20), nrow = 1)
    colnames(exprs) <- c("B1-A", "YG1-A")
    ff <- flowCore::flowFrame(exprs)

    res <- spectreasy::calc_residuals(ff, M, method = "WLS")
    expect_equal(dim(res), c(1, 2))
    expect_true(all(is.finite(as.matrix(res))))
})

test_that("WLS detector noise can be estimated from SCC low-signal tails", {
    scc_dir <- tempfile("scc_noise_")
    dir.create(scc_dir, showWarnings = FALSE)

    exprs1 <- cbind(
        "B1-A" = c(seq(10, 90, length.out = 30), seq(500, 2000, length.out = 70)),
        "YG1-A" = c(seq(20, 120, length.out = 30), seq(600, 2200, length.out = 70)),
        "FSC-A" = seq(1000, 2000, length.out = 100),
        "SSC-A" = seq(500, 900, length.out = 100)
    )
    exprs2 <- cbind(
        "B1-A" = c(seq(15, 100, length.out = 30), seq(550, 2100, length.out = 70)),
        "YG1-A" = c(seq(25, 150, length.out = 30), seq(650, 2300, length.out = 70)),
        "FSC-A" = seq(1100, 2100, length.out = 100),
        "SSC-A" = seq(550, 950, length.out = 100)
    )
    flowCore::write.FCS(flowCore::flowFrame(exprs1), file.path(scc_dir, "FITC.fcs"))
    flowCore::write.FCS(flowCore::flowFrame(exprs2), file.path(scc_dir, "PE.fcs"))

    noise <- spectreasy:::.estimate_wls_detector_noise(
        scc_dir = scc_dir,
        detectors = c("B1-A", "YG1-A")
    )

    expect_equal(noise$detector, c("B1-A", "YG1-A"))
    expect_true(all(is.finite(noise$noise_floor)))
    expect_true(all(noise$noise_floor > 0))
    expect_equal(noise$signal_scale, c(1, 1))
})

test_that("unmix_samples loads sibling SCC detector-noise file for WLS", {
    M <- matrix(c(
        1.0, 0.2, 0.1,
        0.1, 1.0, 0.3
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B1-A", "YG1-A", "R1-A")

    matrix_dir <- tempfile("saved_matrix_noise_")
    dir.create(matrix_dir, showWarnings = FALSE)
    tmp_ref <- file.path(matrix_dir, "scc_reference_matrix.csv")
    ref_df <- as.data.frame(M)
    ref_df$Marker <- rownames(M)
    ref_df <- ref_df[, c("Marker", colnames(M))]
    write.csv(ref_df, tmp_ref, row.names = FALSE)

    detector_noise <- data.frame(
        detector = colnames(M),
        noise_floor = c(125, 400, 125),
        signal_scale = c(1, 1, 1)
    )
    write.csv(detector_noise, file.path(matrix_dir, "scc_detector_noise.csv"), row.names = FALSE)

    exprs <- matrix(c(
        100, 20, 40,
        10, 120, 50
    ), nrow = 2, byrow = TRUE)
    colnames(exprs) <- colnames(M)
    ff <- flowCore::flowFrame(exprs)

    sample_dir <- tempfile("samples_noise_")
    dir.create(sample_dir, showWarnings = FALSE)
    flowCore::write.FCS(ff, file.path(sample_dir, "sample.fcs"))

    output_dir <- tempfile("unmixed_noise_")
    unmixed <- spectreasy::unmix_samples(
        sample_dir = sample_dir,
        unmixing_matrix_file = tmp_ref,
        method = "WLS",
        output_dir = output_dir,
        write_fcs = FALSE
    )

    M_with_noise <- M
    attr(M_with_noise, "detector_noise") <- detector_noise
    expected <- spectreasy::calc_residuals(ff, M_with_noise, method = "WLS")

    expect_equal(
        as.matrix(unmixed$sample$data[, rownames(M)]),
        as.matrix(expected[, rownames(M)]),
        tolerance = 1e-6
    )
})

test_that("saved static unmixing matrices are rejected where reference matrices are expected", {
    W <- matrix(c(
        1.0, -0.2,
        -0.1, 1.0
    ), nrow = 2, byrow = TRUE)
    rownames(W) <- c("FITC", "PE")
    colnames(W) <- c("B1-A", "YG1-A")

    matrix_dir <- tempfile("saved_static_matrix_")
    dir.create(matrix_dir, showWarnings = FALSE)
    static_file <- file.path(matrix_dir, "scc_unmixing_matrix.csv")
    W_df <- as.data.frame(W, check.names = FALSE)
    W_df$Marker <- rownames(W)
    W_df <- W_df[, c("Marker", colnames(W)), drop = FALSE]
    write.csv(W_df, static_file, row.names = FALSE)

    exprs <- matrix(c(100, 20), nrow = 1)
    colnames(exprs) <- colnames(W)
    sample_dir <- tempfile("samples_static_matrix_")
    dir.create(sample_dir, showWarnings = FALSE)
    flowCore::write.FCS(flowCore::flowFrame(exprs), file.path(sample_dir, "sample.fcs"))

    expect_error(
        spectreasy::unmix_samples(
            sample_dir = sample_dir,
            unmixing_matrix_file = static_file,
            write_fcs = FALSE
        ),
        regexp = "static unmixing matrix"
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

test_that("cell population selection applies SSC/FSC ratio only when requested", {
    gmm_result <- list(
        main_populations = 1:3,
        means = matrix(
            c(
                100, 100,
                300, 600,
                700, 500
            ),
            nrow = 2,
            dimnames = list(c("FSC-A", "SSC-A"), NULL)
        )
    )

    with_ratio <- spectreasy:::.select_reference_cell_populations(
        gmm_result,
        fsc_min = 80,
        fsc_max = 800,
        ssc_max = 800,
        ratio_max = 1.25
    )$selected
    without_ratio <- spectreasy:::.select_reference_cell_populations(
        gmm_result,
        fsc_min = 80,
        fsc_max = 800,
        ssc_max = 800,
        ratio_max = Inf
    )$selected

    expect_equal(with_ratio, c(1, 3))
    expect_equal(without_ratio, 1:3)
})

test_that("reference scatter gate accepts non-area scatter channel names", {
    set.seed(1)
    raw_data <- cbind(
        "B1-A" = stats::rlnorm(250, log(100), 0.2),
        "FS-H" = stats::rnorm(250, 100000, 3000),
        "Side Scatter-H" = stats::rnorm(250, 45000, 2500)
    )
    pd <- data.frame(name = colnames(raw_data), stringsAsFactors = FALSE)

    scatter_info <- spectreasy:::.compute_reference_scatter_gate(
        raw_data = raw_data,
        pd = pd,
        sample_type = "beads",
        outlier_percentile = 0.02,
        debris_percentile = 0.08,
        subsample_n = 250,
        max_clusters = 3,
        min_cluster_proportion = 0.03,
        gate_contour_beads = 0.95,
        gate_contour_cells = 0.90,
        bead_gate_scale = 1.3
    )

    expect_false(is.null(scatter_info))
    expect_equal(scatter_info$fsc, "FS-H")
    expect_equal(scatter_info$ssc, "Side Scatter-H")
    expect_gt(nrow(scatter_info$gated_data), 100)
})

test_that("AF profile extraction clusters pooled AF phenotypes", {
    detector_names <- c("B1-A", "YG1-A", "V1-A")
    af_events <- rbind(
        matrix(rep(c(100, 15, 5), 120), ncol = 3, byrow = TRUE),
        matrix(rep(c(10, 90, 20), 120), ncol = 3, byrow = TRUE)
    )
    colnames(af_events) <- detector_names

    profiles <- spectreasy:::.extract_reference_af_profiles(
        detector_names = detector_names,
        n_bands = 2,
        max_cells = 500,
        af_events = af_events
    )

    expect_equal(nrow(profiles$signatures), 2)
    expect_equal(colnames(profiles$signatures), detector_names)
    expect_equal(rownames(profiles$signatures), c("AF", "AF_2"))

    expected_shapes <- rbind(
        c("B1-A" = 1, "YG1-A" = 0.15, "V1-A" = 0.05),
        c("B1-A" = 10 / 90, "YG1-A" = 1, "V1-A" = 20 / 90)
    )
    matched_shapes <- profiles$signatures[
        c(
            which.max(profiles$signatures[, "B1-A"]),
            which.max(profiles$signatures[, "YG1-A"])
        ),
        ,
        drop = FALSE
    ]
    expect_equal(unname(matched_shapes), unname(expected_shapes), tolerance = 1e-6)
    expect_equal(profiles$raw_median, c("B1-A" = 55, "YG1-A" = 52.5, "V1-A" = 12.5), tolerance = 1e-6)
})

test_that("AF profile extraction handles empty and all-zero AF events", {
    detector_names <- c("B1-A", "YG1-A")

    empty_profiles <- spectreasy:::.extract_reference_af_profiles(
        detector_names = detector_names,
        af_events = matrix(numeric(), nrow = 0, ncol = 2, dimnames = list(NULL, detector_names))
    )
    expect_null(empty_profiles$raw_median)
    expect_null(empty_profiles$signatures)

    zero_profiles <- spectreasy:::.extract_reference_af_profiles(
        detector_names = detector_names,
        af_events = matrix(0, nrow = 5, ncol = 2, dimnames = list(NULL, detector_names))
    )
    expect_equal(zero_profiles$raw_median, c("B1-A" = 0, "YG1-A" = 0))
    expect_equal(dim(zero_profiles$signatures), c(1, 2))
    expect_equal(rownames(zero_profiles$signatures), "AF")
    expect_equal(zero_profiles$signatures[1, ], c("B1-A" = 0, "YG1-A" = 0))
})

test_that("AF argument validation supports multi-AF bands per file", {
    args <- spectreasy:::.validate_build_reference_af_args(
        af_n_bands = 2,
        af_bands_per_file = 5,
        af_max_cells = 500
    )

    expect_equal(args$af_n_bands, 2L)
    expect_equal(args$af_bands_per_file, 5L)
    expect_equal(args$af_max_cells, 500L)
    expect_error(
        spectreasy:::.validate_build_reference_af_args(2, 500, af_bands_per_file = 0),
        "af_bands_per_file"
    )
    expect_error(
        spectreasy:::.validate_build_reference_af_args(0, 500, af_bands_per_file = 1),
        "af_n_bands"
    )
    expect_error(
        spectreasy:::.validate_build_reference_af_args(2, 99, af_bands_per_file = 1),
        "af_max_cells"
    )
})

test_that("extra AF files are banked centrally, not processed as one averaged SCC row", {
    af_file <- file.path(tempdir(), "af", "mixed_af.fcs")
    config <- list(include_multi_af = TRUE, af_dir = dirname(af_file), exclude_af = FALSE)

    processed <- spectreasy:::.process_reference_file(
        fcs_file = af_file,
        control_df = NULL,
        sample_patterns = spectreasy::get_fluorophore_patterns(),
        metadata = list(detector_names = c("B1-A")),
        config = config
    )

    expect_null(processed)
})

test_that("cell histogram gating keeps full middle negative mode for bright controls", {
    set.seed(1)
    vals_log <- c(
        stats::rnorm(500, 3.12, 0.18),
        stats::rnorm(500, 3.35, 0.16),
        stats::rnorm(4000, 6.0, 0.12)
    )

    gate <- spectreasy:::.compute_reference_histogram_gate(
        peak_vals = 10^vals_log,
        sample_type = "cells",
        histogram_pct_beads = 0.98,
        histogram_direction_beads = "right",
        histogram_pct_cells = 0.35,
        histogram_direction_cells = "right"
    )

    expect_true(isTRUE(attr(gate$vals_log, "negative_gate_present")))
    expect_lt(attr(gate$vals_log, "neg_log_min"), 3.0)
    expect_gt(attr(gate$vals_log, "neg_log_max"), 3.6)
    expect_gt(log10(gate$gate_min), 5.0)
})

test_that("scatter intensity gating separates nearest negative and bright positive modes", {
    set.seed(7)
    vals_log <- c(
        rep(0, 100),
        stats::rnorm(600, 2.45, 0.18),
        stats::rnorm(800, 5.2, 0.10)
    )

    gate <- spectreasy:::.compute_reference_scatter_intensity_gate(
        peak_vals = 10^vals_log,
        sample_type = "beads",
        histogram_pct_beads = 0.98,
        histogram_direction_beads = "right",
        histogram_pct_cells = 0.35,
        histogram_direction_cells = "right"
    )

    expect_equal(attr(gate$vals_log, "gate_type"), "scatter")
    expect_true(isTRUE(attr(gate$vals_log, "negative_gate_present")))
    expect_gt(attr(gate$vals_log, "neg_log_min"), 1.8)
    expect_lt(attr(gate$vals_log, "neg_log_max"), 3.1)
    expect_gt(log10(gate$gate_min), 4.8)
    expect_lt(log10(gate$gate_max), 5.6)
})

test_that("histogram gating cutoff detection extends right gate leftwards", {
    # 1. Test .compute_reference_histogram_gate density-based cutoff
    set.seed(42)
    # Generate a peak cut-off on the right: positive peak truncated at 5.4
    vals_log <- c(
        stats::rnorm(100, 2.0, 0.2), # negative peak
        pmin(stats::rnorm(1000, 5.35, 0.15), 5.4) # positive peak, truncated at 5.4
    )
    
    gate_right <- spectreasy:::.compute_reference_histogram_gate(
        peak_vals = 10^vals_log,
        sample_type = "cells",
        histogram_pct_beads = 0.98,
        histogram_direction_beads = "right",
        histogram_pct_cells = 0.35,
        histogram_direction_cells = "right"
    )
    
    # Without the cutoff fix, a right-gate would start near the peak (~5.35).
    # With the fix, it should extend to the left trough (around 3.0 - 4.5).
    expect_lt(log10(gate_right$gate_min), 4.8)
    
    # 2. Test gate_positive_cells utility function
    mat <- matrix(pmin(stats::rnorm(1000, 5.35, 0.15), 5.4), ncol = 1)
    
    idx_right <- gate_positive_cells(mat, histogram_pct = 0.35, histogram_direction = "right")
    # If the cutoff fix kicked in, it should behave like "both", which sets lower_q = 0.5 - pct/2 = 0.325.
    # Since the distribution is truncated on the right, upper_q (0.675) lands on the maximum value (5.4).
    # Therefore, the gated fraction is approximately 1 - 0.325 = 0.675 (67.5%).
    expect_gt(mean(idx_right), 0.60)
    expect_lt(mean(idx_right), 0.75)
})

test_that("reference spectrum uses histogram negative gate attributes", {
    detector_names <- c("B1-A", "YG1-A")
    final_gated_data <- matrix(
        rep(c(10000, 200), each = 20),
        ncol = 2,
        dimnames = list(NULL, detector_names)
    )
    gated_data <- rbind(
        matrix(rep(c(10, 5), each = 20), ncol = 2),
        matrix(rep(c(1000, 20), each = 20), ncol = 2)
    )
    colnames(gated_data) <- detector_names
    peak_vals <- gated_data[, "B1-A"]
    vals_log <- log10(pmax(peak_vals, 1))
    attr(vals_log, "negative_gate_present") <- TRUE
    attr(vals_log, "neg_log_min") <- 2.9
    attr(vals_log, "neg_log_max") <- 3.1

    res <- spectreasy:::.compute_reference_spectrum(
        final_gated_data = final_gated_data,
        gated_data = gated_data,
        peak_vals = peak_vals,
        vals_log = vals_log,
        detector_names = detector_names,
        row_info = data.frame(universal.negative = "FALSE", stringsAsFactors = FALSE)
    )

    expect_equal(as.numeric(res), c(1, 180 / 9000), tolerance = 1e-6)
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

test_that("derive_unmixing_matrix WLS ignores SCC variances and honors detector noise", {
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
    expect_equal(spectreasy::derive_unmixing_matrix(M_no_attr, method = "WLS"), W_wls)

    attr(M_no_attr, "detector_noise") <- data.frame(
        detector = colnames(M_no_attr),
        noise_floor = c(125, 200, 125),
        signal_scale = c(1, 1, 1)
    )
    W_noise <- spectreasy::derive_unmixing_matrix(M_no_attr, method = "WLS")
    expect_false(isTRUE(all.equal(W_noise, W_wls)))
})

test_that("calc_residuals WLS performs event-wise noise-model unmixing", {
    M <- matrix(c(
        1, 0.2,
        0.1, 1
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B1-A", "YG1-A")

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

test_that("calc_residuals multi-AF WLS uses event-wise detector weights", {
    M <- matrix(c(
        1.0, 0.2, 0.1, 0.05,
        0.1, 1.0, 0.2, 0.05,
        0.2, 0.2, 1.0, 0.05,
        0.05, 0.05, 0.05, 1.0
    ), nrow = 4, byrow = TRUE)
    rownames(M) <- c("FITC", "PE", "AF", "AF_2")
    colnames(M) <- c("B1-A", "YG1-A", "R1-A", "V1-A")

    exprs <- matrix(c(
        100, 10, 50, 0,
        10, 100, 0, 50
    ), nrow = 2, byrow = TRUE)
    colnames(exprs) <- colnames(M)
    ff <- flowCore::flowFrame(exprs)

    res <- spectreasy::calc_residuals(ff, M, method = "WLS")

    fluor_idx <- which(!grepl("^AF($|_)", rownames(M), ignore.case = TRUE))
    af_idx <- which(grepl("^AF($|_)", rownames(M), ignore.case = TRUE))
    expected <- matrix(0, nrow = nrow(exprs), ncol = nrow(M), dimnames = list(NULL, rownames(M)))

    for (i in seq_len(nrow(exprs))) {
        y <- exprs[i, ]
        detector_weights <- spectreasy:::.wls_event_weights(
            y,
            noise_floor = rep(spectreasy:::.default_wls_background_noise(), ncol(M)),
            signal_scale = rep(1, ncol(M)),
            max_weight_ratio = spectreasy:::.default_wls_max_weight_ratio()
        )
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

test_that("unmix_samples finds sibling variances for saved reference matrix", {
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

    matrix_dir <- tempfile("saved_matrix_")
    dir.create(matrix_dir, showWarnings = FALSE)
    tmp_ref <- file.path(matrix_dir, "scc_reference_matrix.csv")
    tmp_var <- file.path(matrix_dir, "scc_variances.csv")

    ref_df <- as.data.frame(M)
    ref_df$Marker <- rownames(M)
    ref_df <- ref_df[, c("Marker", "B1-A", "YG1-A")]
    write.csv(ref_df, tmp_ref, row.names = FALSE)

    var_df <- as.data.frame(V)
    var_df$Marker <- rownames(V)
    var_df <- var_df[, c("Marker", "B1-A", "YG1-A")]
    write.csv(var_df, tmp_var, row.names = FALSE)

    exprs <- matrix(c(
        100, 20,
        10, 120
    ), nrow = 2, byrow = TRUE)
    colnames(exprs) <- c("B1-A", "YG1-A")
    ff <- flowCore::flowFrame(exprs)

    sample_dir <- tempfile("samples_")
    dir.create(sample_dir, showWarnings = FALSE)
    flowCore::write.FCS(ff, file.path(sample_dir, "sample.fcs"))

    output_dir <- tempfile("unmixed_")

    work_dir <- tempfile("working_dir_")
    dir.create(file.path(work_dir, "spectreasy_outputs", "unmix_controls"), recursive = TRUE)
    wrong_var <- matrix(1, nrow = 1, ncol = 2, dimnames = list("Wrong", colnames(M)))
    wrong_var_df <- as.data.frame(wrong_var)
    wrong_var_df$Marker <- rownames(wrong_var)
    wrong_var_df <- wrong_var_df[, c("Marker", "B1-A", "YG1-A")]
    write.csv(
        wrong_var_df,
        file.path(work_dir, "spectreasy_outputs", "unmix_controls", "scc_variances.csv"),
        row.names = FALSE
    )
    withr::local_dir(work_dir)

    unmixed <- spectreasy::unmix_samples(
        sample_dir = sample_dir,
        unmixing_matrix_file = tmp_ref,
        method = "WLS",
        output_dir = output_dir,
        write_fcs = TRUE
    )

    expect_true(file.exists(file.path(output_dir, "sample_unmixed.fcs")))
    expect_setequal(names(unmixed), "sample")
})
