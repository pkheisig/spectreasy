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
    expect_error(spectreasy::gating_options(histogram_pct_beads = NA_real_), "finite numeric")
    expect_error(spectreasy::gating_options(histogram_pct_cells = -0.1), "must be in")
    expect_error(spectreasy::gating_options(histogram_pct_cells = 1.1), "must be in")
})

test_that("reference gating fractions reject missing and impossible values", {
    missing_input <- tempfile("missing_scc_")
    expect_error(
        spectreasy::build_reference_matrix(input_folder = missing_input, outlier_percentile = NA_real_),
        "outlier_percentile must be one finite numeric value"
    )
    expect_error(
        spectreasy::build_reference_matrix(input_folder = missing_input, gate_contour_beads = 1),
        "gate_contour_beads must be in"
    )
    expect_error(
        spectreasy::build_reference_matrix(input_folder = missing_input, bead_gate_scale = -1),
        "bead_gate_scale must be one finite number greater than 0"
    )
})

test_that("obsolete scatter-intensity gating is absent from the public and internal API", {
    public_functions <- list(
        spectreasy::gating_options,
        spectreasy::build_reference_matrix,
        spectreasy::unmix_controls,
        spectreasy::qc_controls,
        spectreasy::collect_control_report_data
    )
    expect_true(all(vapply(
        public_functions,
        function(fun) !"use_scatter_gating" %in% names(formals(fun)),
        logical(1)
    )))
    expect_false(exists(
        ".compute_reference_scatter_intensity_gate",
        envir = asNamespace("spectreasy"),
        inherits = FALSE
    ))
})

test_that("SCC background cleanup is method-driven without overlapping booleans", {
    expect_false("clean_scc_with_unstained" %in% names(formals(spectreasy::unmix_controls)))
    expect_false("clean_scc_with_unstained" %in% names(formals(spectreasy::build_reference_matrix)))
    expect_false("autospectral_scc_cleanup" %in% names(formals(spectreasy::build_reference_matrix)))
    expect_identical(formals(spectreasy::build_reference_matrix)$unmixing_method, "AutoSpectral")
    expect_true("scc_background_method" %in% names(formals(spectreasy::unmix_controls)))
    expect_true("scc_background_method" %in% names(formals(spectreasy::build_reference_matrix)))
})

test_that("all workflow-facing unmixing defaults use AutoSpectral", {
    expect_identical(formals(spectreasy::calc_residuals)$method, "AutoSpectral")
    expect_identical(formals(spectreasy::build_reference_matrix)$unmixing_method, "AutoSpectral")
    expect_identical(formals(spectreasy::unmix_controls)$unmixing_method, "AutoSpectral")
    expect_identical(formals(spectreasy::unmix_samples)$unmixing_method, "AutoSpectral")
    expect_identical(formals(spectreasy::adjust_matrix)$unmixing_method, "AutoSpectral")
})

test_that("documented diagnostic helpers are exported", {
    expect_true(all(c(
        "calculate_similarity_matrix",
        "plot_similarity_matrix",
        "plot_sample_rms_residuals"
    ) %in% getNamespaceExports("spectreasy")))
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

    W_wls_default <- spectreasy::derive_unmixing_matrix(M, method = "WLS")
    expect_equal(dim(W_wls_default), dim(M))
    expect_true(all(is.finite(W_wls_default)))

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

test_that("calc_residuals RWLS down-weights detector-level outliers", {
    M <- rbind(
        FITC = c(1, 0.2, 0.05, 0.1),
        PE = c(0.1, 1, 0.1, 0.05)
    )
    colnames(M) <- c("D1", "D2", "D3", "D4")

    true_abundance <- c(FITC = 100, PE = 20)
    y <- as.numeric(true_abundance %*% M)
    y[3] <- y[3] + 500
    ff <- flowCore::flowFrame(matrix(y, nrow = 1, dimnames = list(NULL, colnames(M))))

    wls <- spectreasy::calc_residuals(ff, M, method = "WLS")[1, names(true_abundance)]
    rwls <- spectreasy::calc_residuals(ff, M, method = "RWLS")[1, names(true_abundance)]

    wls_error <- sum((as.numeric(wls) - true_abundance)^2)
    rwls_error <- sum((as.numeric(rwls) - true_abundance)^2)
    expect_lt(rwls_error, wls_error * 0.1)
})

test_that("rwls_max_iter is exposed through the unmixing APIs", {
    expect_true("spectreasy_gui" %in% getNamespaceExports("spectreasy"))
    expect_identical(formals(spectreasy::spectreasy_gui)$port, 8000)
    expect_true("rwls_max_iter" %in% names(formals(spectreasy::calc_residuals)))
    expect_true("rwls_max_iter" %in% names(formals(spectreasy::unmix_samples)))
    expect_true("samples_dir" %in% names(formals(spectreasy::unmix_samples)))
    expect_true("plot_n_events" %in% names(formals(spectreasy::unmix_samples)))
    expect_false("subsample_n" %in% names(formals(spectreasy::unmix_samples)))
    expect_true("rwls_max_iter" %in% names(formals(spectreasy::unmix_controls)))
    expect_true("n_threads" %in% names(formals(spectreasy::unmix_controls)))
    expect_false("unmix_threads" %in% names(formals(spectreasy::unmix_controls)))
    expect_true("save_qc_png" %in% names(formals(spectreasy::unmix_controls)))
    expect_true("autospectral_refine" %in% names(formals(spectreasy::unmix_controls)))
    expect_false("save_qc_plots" %in% names(formals(spectreasy::unmix_controls)))
    expect_false("refine" %in% names(formals(spectreasy::unmix_controls)))
    expect_true("estimate_af" %in% names(formals(spectreasy::unmix_samples)))
    expect_true("unmixing_method" %in% names(formals(spectreasy::adjust_matrix)))
    expect_false(formals(spectreasy::unmix_samples)$estimate_af)
    expect_equal(formals(spectreasy::calc_residuals)$method, "AutoSpectral")
    expect_equal(formals(spectreasy::unmix_controls)$unmixing_method, "AutoSpectral")
    expect_equal(formals(spectreasy::unmix_samples)$unmixing_method, "AutoSpectral")
    expect_equal(formals(spectreasy::unmix_samples)$plot_n_events, 10000L)
    expect_equal(formals(spectreasy::adjust_matrix)$unmixing_method, "AutoSpectral")
})

test_that("unmix_controls validates n_threads", {
    expect_error(
        spectreasy::unmix_controls(scc_dir = tempfile(), n_threads = 0),
        "n_threads must be an integer >= 1"
    )
})

test_that("unmix_controls rejects removed parameter names", {
    for (removed_arg in c("unmix_threads", "save_qc_plots", "refine")) {
        args <- list(scc_dir = tempfile())
        args[[removed_arg]] <- FALSE
        expect_error(
            do.call(spectreasy::unmix_controls, args),
            paste0("no longer accepts: ", removed_arg),
            fixed = TRUE
        )
    }
})

test_that("calc_residuals multi-AF RWLS honors rwls_max_iter", {
    M <- rbind(
        FITC = c(1.00, 0.15, 0.10, 0.05, 0.02, 0.01),
        PE = c(0.10, 1.00, 0.15, 0.05, 0.02, 0.01),
        AF = c(0.20, 0.20, 0.80, 0.20, 0.10, 0.05),
        AF_2 = c(0.05, 0.10, 0.20, 0.80, 0.20, 0.10)
    )
    colnames(M) <- paste0("D", seq_len(ncol(M)))

    true_abundance <- c(FITC = 80, PE = 25, AF = 0, AF_2 = 120)
    y <- as.numeric(true_abundance %*% M)
    y[3] <- y[3] + 1000
    y[6] <- y[6] - 250
    ff <- flowCore::flowFrame(matrix(y, nrow = 1, dimnames = list(NULL, colnames(M))))

    iter_1 <- spectreasy::calc_residuals(ff, M, method = "RWLS", rwls_max_iter = 1)
    iter_5 <- spectreasy::calc_residuals(ff, M, method = "RWLS", rwls_max_iter = 5)

    expect_false(isTRUE(all.equal(
        as.matrix(iter_1[, rownames(M)]),
        as.matrix(iter_5[, rownames(M)]),
        tolerance = 1e-8
    )))
    expect_error(
        spectreasy::calc_residuals(ff, M, method = "RWLS", rwls_max_iter = 0),
        "rwls_max_iter must be an integer >= 1"
    )
})

test_that("calc_residuals WLS works without detector-noise metadata", {
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
        noise_floor = c(125, 400, 125)
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
        unmixing_method = "WLS",
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

test_that("saved SCC detector noise contains only canonical persisted columns", {
    detector_noise <- data.frame(
        detector = c("B1-A", "YG1-A"),
        noise_floor = c(125, 400),
        signal_scale = c(1, 1),
        transient_internal_value = c(2, 3)
    )
    path <- tempfile(fileext = ".csv")

    spectreasy:::.save_detector_noise_csv(detector_noise, path)

    saved <- utils::read.csv(path, check.names = FALSE)
    expect_identical(colnames(saved), c("detector", "noise_floor"))
    expect_equal(saved$detector, detector_noise$detector)
    expect_equal(saved$noise_floor, detector_noise$noise_floor)
})

test_that("unmix_samples can estimate missing AF from stained samples", {
    set.seed(42)
    M_marker <- rbind(
        FITC = c(1.00, 0.12, 0.04),
        PE = c(0.08, 1.00, 0.10)
    )
    colnames(M_marker) <- c("B1-A", "YG1-A", "R1-A")
    af_signature <- c("B1-A" = 0.20, "YG1-A" = 0.25, "R1-A" = 1.00)

    n <- 320
    abund <- cbind(
        FITC = rgamma(n, shape = 2, rate = 0.03),
        PE = rgamma(n, shape = 2, rate = 0.04)
    )
    af_abund <- rgamma(n, shape = 2, rate = 0.04)
    Y <- abund %*% M_marker +
        af_abund %*% matrix(af_signature, nrow = 1) +
        matrix(rnorm(n * ncol(M_marker), sd = 0.15), nrow = n)
    colnames(Y) <- colnames(M_marker)
    exprs_mat <- cbind(
        Y,
        "FSC-A" = rnorm(n, 100000, 4000),
        "SSC-A" = rnorm(n, 50000, 3000)
    )
    ff <- flowCore::flowFrame(exprs_mat)
    fs <- flowCore::flowSet(list(stained = ff))

    without_af <- spectreasy::unmix_samples(
        sample_dir = fs,
        M = M_marker,
        unmixing_method = "WLS",
        write_fcs = FALSE,
        verbose = FALSE
    )
    with_af <- spectreasy::unmix_samples(
        sample_dir = fs,
        M = M_marker,
        unmixing_method = "WLS",
        estimate_af = TRUE,
        write_fcs = FALSE,
        verbose = FALSE,
        seed = 11
    )

    M_est <- attr(with_af, "reference_matrix")
    expect_true(any(grepl("^AF($|_)", rownames(M_est), ignore.case = TRUE)))
    expect_true("AF" %in% colnames(with_af$stained$data))

    rmse_without <- sqrt(mean(without_af$stained$residuals^2))
    rmse_with <- sqrt(mean(with_af$stained$residuals^2))
    expect_lt(rmse_with, rmse_without * 0.8)

    blind_info <- attr(with_af, "blind_af_info")
    expect_false(is.null(blind_info))
    expect_gt(blind_info$derived_bands, 0)
    expect_true(blind_info$model_id %in% c("residual_wls_q90_k10", "residual_wls_low_marker_q90_k10"))
    expect_equal(blind_info$candidate_quantile, 0.90)
    expect_equal(blind_info$requested_bands, 10)
    expect_equal(blind_info$selected_by, "heldout_event_wise_wls")
    expect_s3_class(blind_info$heldout_scores, "data.frame")
    expect_setequal(
        blind_info$heldout_scores$model_id,
        c("residual_wls_q90_k10", "residual_wls_low_marker_q90_k10")
    )
})

test_that("multi-AF WLS threading matches single-threaded output", {
    M <- matrix(c(
        1.00, 0.20, 0.05, 0.01,
        0.10, 1.00, 0.20, 0.05,
        0.05, 0.10, 1.00, 0.20,
        0.30, 0.25, 0.15, 0.10,
        0.08, 0.18, 0.28, 0.38
    ), nrow = 5, byrow = TRUE)
    rownames(M) <- c("FITC", "PE", "APC", "AF", "AF_2")
    colnames(M) <- c("B1-A", "YG1-A", "R1-A", "V1-A")

    coeffs <- matrix(c(
        120, 15, 8, 25, 0,
        20, 110, 12, 0, 35,
        15, 10, 95, 18, 0,
        65, 55, 30, 0, 22
    ), nrow = 4, byrow = TRUE)
    Y <- coeffs %*% M
    colnames(Y) <- colnames(M)
    ff <- flowCore::flowFrame(Y)

    single <- spectreasy::calc_residuals(ff, M, method = "WLS", n_threads = 1)
    threaded <- spectreasy::calc_residuals(ff, M, method = "WLS", n_threads = 2)

    expect_equal(
        as.matrix(threaded[, rownames(M)]),
        as.matrix(single[, rownames(M)]),
        tolerance = 1e-8
    )
    expect_error(
        spectreasy::calc_residuals(ff, M, method = "WLS", n_threads = 0),
        "n_threads must be an integer >= 1"
    )
})

test_that("event-wise solvers match across thread counts", {
    set.seed(20260714)
    M <- rbind(
        FITC = c(1.00, 0.20, 0.05, 0.02),
        PE = c(0.08, 1.00, 0.18, 0.04),
        APC = c(0.03, 0.10, 1.00, 0.20),
        AF = c(0.30, 0.24, 0.15, 0.10),
        AF_2 = c(0.07, 0.17, 0.29, 0.40)
    )
    colnames(M) <- c("B1-A", "YG1-A", "R1-A", "V1-A")
    n <- 320L
    marker_coeffs <- matrix(rexp(n * 3L, rate = 0.02), nrow = n)
    af_choice <- sample.int(2L, n, replace = TRUE)
    af_coeffs <- matrix(0, nrow = n, ncol = 2L)
    af_coeffs[cbind(seq_len(n), af_choice)] <- rexp(n, rate = 0.025)
    Y <- cbind(marker_coeffs, af_coeffs) %*% M
    colnames(Y) <- colnames(M)
    ff <- flowCore::flowFrame(Y)

    for (method in c("AutoSpectral", "OLS", "NNLS", "WLS", "RWLS")) {
        args <- list(flow_frame = ff, M = M, method = method)
        if (identical(method, "RWLS")) args$rwls_max_iter <- 2L
        single <- do.call(spectreasy::calc_residuals, c(args, list(n_threads = 1L)))
        threaded <- do.call(spectreasy::calc_residuals, c(args, list(n_threads = 2L)))
        expect_equal(
            as.matrix(threaded[, rownames(M)]),
            as.matrix(single[, rownames(M)]),
            tolerance = 1e-8,
            info = method
        )
    }
})

test_that("threaded AutoSpectral AF assignment preserves the reference calculation", {
    set.seed(20260715)
    marker_spectra <- rbind(
        FITC = c(1.00, 0.12, 0.03, 0.01),
        PE = c(0.06, 1.00, 0.14, 0.02)
    )
    af_spectra <- rbind(
        AF = c(0.35, 0.28, 0.16, 0.10),
        AF_2 = c(0.08, 0.17, 0.31, 0.42),
        AF_3 = c(0.18, 0.35, 0.22, 0.16)
    )
    colnames(marker_spectra) <- colnames(af_spectra) <- paste0("D", 1:4)
    Y <- matrix(rexp(320L * 4L, rate = 0.01), ncol = 4L)
    colnames(Y) <- colnames(marker_spectra)

    unmixing_matrix <- solve(marker_spectra %*% t(marker_spectra), marker_spectra)
    marker_t <- t(marker_spectra)
    v_library <- unmixing_matrix %*% t(af_spectra)
    r_library <- t(af_spectra) - marker_t %*% v_library
    denominator <- colSums(r_library^2)
    numerator <- Y %*% r_library
    k_matrix <- sweep(numerator, 2, denominator, "/")
    unmixed_markers <- Y %*% t(unmixing_matrix)
    error_matrix <- vapply(seq_len(nrow(af_spectra)), function(i) {
        rowSums(abs(unmixed_markers - k_matrix[, i, drop = FALSE] %*% t(v_library[, i, drop = FALSE])))
    }, numeric(nrow(Y)))
    expected <- max.col(-error_matrix, ties.method = "first")

    single <- spectreasy:::.autospectral_assign_af_fluorophores(
        Y, marker_spectra, af_spectra, n_threads = 1L
    )
    threaded <- spectreasy:::.autospectral_assign_af_fluorophores(
        Y, marker_spectra, af_spectra, n_threads = 2L
    )
    expect_identical(single, expected)
    expect_identical(threaded, expected)
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
        unmixing_method = "OLS",
        output_dir = output_dir,
        write_fcs = TRUE
    )
    expect_setequal(colnames(unmixed$sample1$data), c("FITC", "PE", "Time", "FSC-A", "FSC-H", "SSC-A", "SSC-W", "File"))

    unmixed_ff <- flowCore::read.FCS(
        file.path(output_dir, "unmix_samples", "unmixed_fcs", "sample1_OLS-0AF.fcs"),
        transformation = FALSE,
        truncate_max_range = FALSE
    )
    out_cols <- colnames(flowCore::exprs(unmixed_ff))

    expect_true(all(c("FITC", "PE", "Time", "FSC-A", "FSC-H", "SSC-A", "SSC-W") %in% out_cols))
    expect_false(any(c("B1-A", "YG1-A") %in% out_cols))
})

test_that("unmix_samples versions the output directory and supports flowSet return", {
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

    unmixed_dir <- file.path(output_dir, "unmix_samples", "unmixed_fcs")
    dir.create(unmixed_dir, showWarnings = FALSE, recursive = TRUE)
    flowCore::write.FCS(ff, file.path(unmixed_dir, "sample1_OLS-0AF.fcs"))

    unmixed_fs <- spectreasy::unmix_samples(
        sample_dir = fs,
        M = M,
        unmixing_method = "OLS",
        output_dir = output_dir,
        write_fcs = TRUE,
        return_type = "flowSet"
    )

    expect_s4_class(unmixed_fs, "flowSet")
    expect_setequal(flowCore::sampleNames(unmixed_fs), c("sample1", "sample2"))
    expect_true(file.exists(file.path(unmixed_dir, "sample1_OLS-0AF.fcs")))
    expect_true(file.exists(file.path(paste0(unmixed_dir, "_2"), "sample1_OLS-0AF.fcs")))
    expect_true(file.exists(file.path(paste0(unmixed_dir, "_2"), "sample2_OLS-0AF.fcs")))

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
            unmixing_method = "OLS",
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

test_that("downsampled SingleCellExperiment output preserves source cell identities", {
    skip_if_not_installed("SingleCellExperiment")
    skip_if_not_installed("SummarizedExperiment")

    M <- diag(2)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B1-A", "YG1-A")
    assay_mat <- matrix(
        seq_len(20),
        nrow = 2,
        dimnames = list(colnames(M), paste0("original_cell_", seq_len(10)))
    )
    toy_sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = assay_mat),
        colData = S4Vectors::DataFrame(sample_id = rep("sampleA", 10))
    )

    unmixed_sce <- spectreasy::unmix_samples(
        sample_dir = toy_sce,
        M = M,
        unmixing_method = "OLS",
        write_fcs = FALSE,
        save_report = FALSE,
        plot_n_events = 3L,
        seed = 7L,
        return_type = "SingleCellExperiment",
        verbose = FALSE
    )

    source_ids <- sub("^sampleA__", "", colnames(unmixed_sce))
    expect_length(source_ids, 3L)
    expect_true(all(source_ids %in% colnames(toy_sce)))
    expect_false(any(grepl("^sampleA_[0-9]+$", colnames(unmixed_sce))))
})

test_that("SingleCellExperiment sample basenames remain unique", {
    skip_if_not_installed("SingleCellExperiment")
    skip_if_not_installed("SummarizedExperiment")

    M <- matrix(c(
        1, 0.2,
        0.1, 1
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B1-A", "YG1-A")

    assay_mat <- matrix(c(
        100, 20, 90, 30,
        120, 15, 110, 25
    ), nrow = 2, byrow = TRUE)
    rownames(assay_mat) <- c("B1-A", "YG1-A")
    colnames(assay_mat) <- paste0("cell_", seq_len(ncol(assay_mat)))
    toy_sce <- suppressWarnings(
        SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = assay_mat),
            colData = S4Vectors::DataFrame(
                sample_id = rep(c("run1/sample.fcs", "run2/sample.fcs"), each = 2),
                row.names = colnames(assay_mat)
            )
        )
    )

    out <- suppressWarnings(
        spectreasy::unmix_samples(
            sample_dir = toy_sce,
            M = M,
            unmixing_method = "OLS",
            write_fcs = FALSE,
            save_report = FALSE
        )
    )

    expect_equal(names(out), c("sample", "sample.1"))
    expect_length(out, 2)
})

test_that(".compute_reference_spectrum computes detector spectrum", {
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
    expect_null(attr(res, "variance"))
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

test_that(".compute_reference_spectrum uses unstained bead negative for bead SCCs", {
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
    row_info <- data.frame(universal.negative = "", stringsAsFactors = FALSE)
    bead_negative <- c("B1-A" = 30, "YG1-A" = 20)

    res <- spectreasy:::.compute_reference_spectrum(
        final_gated_data = final_gated_data,
        gated_data = gated_data,
        peak_vals = peak_vals,
        vals_log = vals_log,
        detector_names = detector_names,
        row_info = row_info,
        sample_type = "beads",
        bead_negative = bead_negative
    )

    expect_equal(as.numeric(res), c(1, 40 / 70), tolerance = 1e-6)
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

test_that("SCC report reference overlay includes only single-band AF", {
    M_one_af <- rbind(
        FITC = c(1, 0.2),
        PE = c(0.2, 1),
        AF = c(0.1, 0.3)
    )
    colnames(M_one_af) <- c("B1-A", "YG1-A")

    M_multi_af <- rbind(
        M_one_af,
        AF_2 = c(0.4, 0.1)
    )

    expect_equal(rownames(spectreasy:::.scc_reference_overlay_matrix(M_one_af)), c("FITC", "PE", "AF"))
    expect_equal(rownames(spectreasy:::.scc_reference_overlay_matrix(M_multi_af)), c("FITC", "PE"))
})

test_that("single-band AF uses robust median normalized shape", {
    detector_names <- c("B1-A", "YG1-A")
    af_events <- rbind(
        matrix(rep(c(100, 10), 5), ncol = 2, byrow = TRUE),
        c(1, 100)
    )
    colnames(af_events) <- detector_names

    profiles <- spectreasy:::.extract_reference_af_profiles(
        detector_names = detector_names,
        n_bands = 1,
        max_cells = 100,
        af_events = af_events
    )

    expect_equal(profiles$selection$method, "median_fixed")
    expect_equal(as.numeric(profiles$signatures[1, ]), c(1, 0.1), tolerance = 1e-6)
})

test_that("fixed AF bank size is the default for AF extraction APIs", {
    expect_equal(formals(spectreasy::build_reference_matrix)$af_n_bands, 100)
    expect_equal(formals(spectreasy::unmix_controls)$af_n_bands, 100)
    expect_true("af_profile" %in% names(formals(spectreasy::unmix_controls)))
    expect_null(formals(spectreasy::unmix_controls)$af_profile)
    expect_equal(formals(spectreasy::extract_af_profile)$af_n_bands, 100)
    expect_false("af_n_bands" %in% names(formals(spectreasy::unmix_samples)))
    for (fun in list(
        spectreasy::build_reference_matrix,
        spectreasy::unmix_controls,
        spectreasy::extract_af_profile
    )) {
        expect_false(any(c("af_min_cluster_events", "af_min_cluster_proportion") %in% names(formals(fun))))
    }
})

test_that("fixed AF k-means returns the exact requested band count", {
    detector_names <- c("B1-A", "YG1-A", "V1-A")
    af_events <- rbind(
        matrix(rep(c(100, 15, 5), 4978), ncol = 3, byrow = TRUE),
        matrix(rep(c(10, 90, 20), 22), ncol = 3, byrow = TRUE)
    )
    colnames(af_events) <- detector_names

    profiles <- spectreasy:::.extract_reference_af_profiles(
        detector_names = detector_names,
        n_bands = 2,
        max_cells = 5000,
        af_events = af_events
    )

    expect_equal(nrow(profiles$signatures), 2)
    expect_equal(profiles$selection$method, "kmeans_fixed")
    expect_equal(profiles$selection$requested_bands, 2)
    expect_equal(profiles$selection$final_bands, 2)
    expect_equal(profiles$selection$cluster_sizes, c(4978L, 22L))
})

test_that("fixed AF k-means fails instead of returning fewer bands", {
    shape <- matrix(rep(c(1, 0.2, 0.1), 100), ncol = 3, byrow = TRUE)
    colnames(shape) <- c("B1-A", "YG1-A", "V1-A")

    expect_error(
        spectreasy:::.reference_kmeans_af_centers(shape, n_centers = 2),
        "Cannot derive exactly 2 AF bands"
    )
})

test_that("a request for 100 AF bands returns exactly 100 bands", {
    x <- seq(0, 1, length.out = 500)
    shapes <- cbind(
        B1.A = 1,
        YG1.A = 0.05 + 0.9 * x,
        V1.A = 0.05 + 0.9 * x^2,
        R1.A = 0.05 + 0.9 * sqrt(x)
    )

    result <- spectreasy:::.reference_kmeans_af_centers(
        af_shape = shapes,
        n_centers = 100,
        nstart = 1
    )

    expect_equal(nrow(result$centers), 100)
    expect_equal(length(result$cluster_sizes), 100)
    expect_equal(sum(result$cluster_sizes), nrow(shapes))
    expect_true(all(result$cluster_sizes > 0L))
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
        n_bands = 1,
        af_events = matrix(0, nrow = 5, ncol = 2, dimnames = list(NULL, detector_names))
    )
    expect_equal(zero_profiles$raw_median, c("B1-A" = 0, "YG1-A" = 0))
    expect_equal(dim(zero_profiles$signatures), c(1, 2))
    expect_equal(rownames(zero_profiles$signatures), "AF")
    expect_equal(zero_profiles$signatures[1, ], c("B1-A" = 0, "YG1-A" = 0))
})

test_that("AF argument validation requires fixed numeric bands", {
    args <- spectreasy:::.validate_build_reference_af_args(
        af_n_bands = 2,
        af_max_cells = 500
    )

    expect_equal(args$af_n_bands, 2L)
    expect_equal(args$af_max_cells, 500L)
    expect_error(
        spectreasy:::.validate_build_reference_af_args(0, 500),
        "af_n_bands"
    )
    expect_error(
        spectreasy:::.validate_build_reference_af_args("auto", 500),
        "af_n_bands"
    )
    expect_error(
        spectreasy:::.validate_build_reference_af_args(2, 99),
        "af_max_cells"
    )
})

test_that("mapped SCC AF files are pooled into one AF bank size request", {
    detector_names <- c("B1-A", "YG1-A")
    fcs_files <- file.path(tempdir(), paste0(c("Unstained_1", "Unstained_2", "Unstained_3", "UnstainedDead"), ".fcs"))
    control_df <- data.frame(
        filename = basename(fcs_files),
        fluorophore = c("AF", "AF_2", "AF_3", "AF_Internal"),
        marker = "Autofluorescence",
        control.type = "cells",
        stringsAsFactors = FALSE
    )
    captured <- new.env(parent = emptyenv())
    testthat::local_mocked_bindings(
        .extract_reference_af_gated_events = function(fcs_file, detector_names, config = NULL, sample_type = "cells") {
            i <- match(normalizePath(fcs_file, mustWork = FALSE), normalizePath(fcs_files, mustWork = FALSE))
            events <- matrix(i, nrow = 3, ncol = length(detector_names), dimnames = list(NULL, detector_names))
            list(
                events = events,
                source = data.table::data.table(
                    file = basename(fcs_file),
                    path = fcs_file,
                    n_total = 3L,
                    n_scatter_gated = 3L,
                    scatter_gate_pct = 100,
                    fsc_channel = "FSC-A",
                    ssc_channel = "SSC-A"
                )
            )
        },
        .extract_reference_af_profiles = function(detector_names, n_bands, max_cells, af_events) {
            captured$n_bands <- n_bands
            captured$event_count <- nrow(af_events)
            n_bands <- max(1L, as.integer(n_bands))
            signatures <- matrix(0.1, nrow = n_bands, ncol = length(detector_names))
            for (i in seq_len(n_bands)) {
                signatures[i, ((i - 1L) %% length(detector_names)) + 1L] <- 1
            }
            rownames(signatures) <- c("AF", if (n_bands > 1L) paste0("AF_", seq.int(2L, n_bands)) else NULL)
            colnames(signatures) <- detector_names
            list(
                raw_median = stats::setNames(rep(2, length(detector_names)), detector_names),
                signatures = signatures,
                selection = list(method = "kmeans_fixed", n_bands = n_bands)
            )
        },
        .env = asNamespace("spectreasy")
    )

    out <- spectreasy:::.collect_reference_af_profiles(
        control_df = control_df,
        fcs_files = fcs_files,
        detector_names = detector_names,
        af_n_bands = 3L,
        af_max_cells = 500L,
        fcs_files_all = fcs_files
    )

    expect_equal(captured$n_bands, 3L)
    expect_equal(captured$event_count, 9L)
    expect_equal(out$af_bank_info$source_count, 3L)
    expect_equal(out$af_bank_info$requested_bands, 3L)
    expect_equal(out$af_bank_info$mode, "pooled_af_sources")
    expect_equal(out$af_bank_info$sources$source_type, rep("mapped_unstained", 3))
    expect_false("UnstainedDead.fcs" %in% out$af_bank_info$sources$file)
})

test_that("legacy ordinary AF_Internal mappings are canonicalized as AF bank sources", {
    control_df <- data.frame(
        filename = c(
            "Reference Group_Unstained (Cells).fcs",
            "Reference Group_Unstained_Fbs (Cells).fcs",
            "Reference Group_Unstained_Tumor (Cells).fcs"
        ),
        fluorophore = c("AF", "AF_Internal", "AF_Internal"),
        marker = "Autofluorescence",
        channel = c("V7-A", "V7-A", "UV7-A"),
        control.type = "cells",
        stringsAsFactors = FALSE
    )

    normalized <- spectreasy:::.normalize_build_reference_control_df(control_df)

    expect_equal(normalized$fluorophore, c("AF", "AF_2", "AF_3"))
    expect_true(all(spectreasy:::.is_primary_af_control_row(
        fluorophore = normalized$fluorophore,
        marker = normalized$marker,
        filename = normalized$filename
    )))
})

test_that("viability controls auto-use dead cell AF negatives", {
    control_df <- data.frame(
        filename = c(
            "scc_cells_AF_Unstained.fcs",
            "scc_cells_AF_UnstainedDead.fcs",
            "scc_cells_eFluor780_LiveDead.fcs",
            "scc_cells_FITC_CD4.fcs"
        ),
        fluorophore = c("AF", "AF_Internal", "eFluor 780", "FITC"),
        marker = c("Autofluorescence", "Autofluorescence", "Live", "CD4"),
        channel = c("FL13-A", "FL11-A", "FL44-A", "FL37-A"),
        control.type = c("cells", "cells", "cells", "cells"),
        is.viability = c("", "", "TRUE", ""),
        stringsAsFactors = FALSE
    )

    normalized <- spectreasy:::.normalize_build_reference_control_df(control_df)

    expect_equal(normalized$universal.negative[normalized$filename == "scc_cells_eFluor780_LiveDead.fcs"], "scc_cells_AF_UnstainedDead.fcs")
    expect_equal(
        normalized$universal.negative[normalized$filename == "scc_cells_FITC_CD4.fcs"],
        "scc_cells_AF_Unstained.fcs"
    )
})

test_that("automatic negatives follow cell, viability, and bead source types", {
    control_df <- data.frame(
        filename = c("af.fcs", "dead.fcs", "bead-neg.fcs", "fitc.fcs", "live.fcs", "pe-beads.fcs"),
        fluorophore = c("AF", "AF_dead", "AF_beads", "FITC", "Zombie NIR", "PE"),
        marker = c("Autofluorescence", "Dead cell background", "Bead background", "CD4", "Live", "CD8"),
        channel = c("V1-A", "R7-A", "V1-A", "B1-A", "R7-A", "YG1-A"),
        control.type = c("cells", "cells", "beads", "cells", "cells", "beads"),
        universal.negative = "",
        is.viability = c("", "", "", "", "TRUE", ""),
        stringsAsFactors = FALSE
    )

    normalized <- spectreasy:::.normalize_build_reference_control_df(control_df)
    by_file <- split(normalized, normalized$filename)

    expect_equal(by_file[["fitc.fcs"]]$universal.negative, "af.fcs")
    expect_equal(by_file[["live.fcs"]]$universal.negative, "dead.fcs")
    expect_equal(by_file[["pe-beads.fcs"]]$universal.negative, "bead-neg.fcs")
    expect_equal(by_file[["af.fcs"]]$universal.negative, "")
    expect_equal(by_file[["dead.fcs"]]$universal.negative, "")
    expect_equal(by_file[["bead-neg.fcs"]]$universal.negative, "")
})

test_that("legacy AF negative labels are replaced by the exact AF filename", {
    control_df <- data.frame(
        filename = c("Unstained_A2_Without (Cells).fcs", "BUV661 (Cells).fcs"),
        fluorophore = c("AF", "BUV661"),
        marker = c("Autofluorescence", "TCR"),
        channel = c("V7-A", "UV11-A"),
        control.type = "cells",
        universal.negative = c("", "AF"),
        is.viability = "",
        stringsAsFactors = FALSE
    )

    normalized <- spectreasy:::.normalize_build_reference_control_df(control_df)

    expect_equal(
        normalized$universal.negative[normalized$filename == "BUV661 (Cells).fcs"],
        "Unstained_A2_Without (Cells).fcs"
    )
})

test_that("SCC negative resolver selects AF_dead, AF_beads, and AF backgrounds", {
    af_background <- list(method = "scatter_knn", n_events = 101L)
    dead_background <- list(method = "scatter_knn", n_events = 51L)
    bead_background <- list(method = "scatter_knn", n_events = 81L)
    af_negative <- stats::setNames(c(1, 2), c("A", "B"))
    dead_negative <- stats::setNames(c(3, 4), c("A", "B"))
    bead_negative <- stats::setNames(c(5, 6), c("A", "B"))
    attr(dead_negative, "scc_background") <- dead_background
    attr(bead_negative, "scc_background") <- bead_background
    negatives <- list(dead = dead_negative)

    ordinary <- spectreasy:::.resolve_reference_scc_negative_source(
        row_info = data.frame(universal.negative = ""),
        sample_type = "cells",
        af_data_raw = af_negative,
        scc_background = af_background,
        universal_negatives = negatives,
        bead_negative = bead_negative
    )
    viability <- spectreasy:::.resolve_reference_scc_negative_source(
        row_info = data.frame(universal.negative = "dead.fcs"),
        sample_type = "cells",
        af_data_raw = af_negative,
        scc_background = af_background,
        universal_negatives = negatives,
        bead_negative = bead_negative
    )
    beads <- spectreasy:::.resolve_reference_scc_negative_source(
        row_info = data.frame(universal.negative = "bead-neg.fcs"),
        sample_type = "beads",
        af_data_raw = af_negative,
        scc_background = af_background,
        universal_negatives = negatives,
        bead_negative = bead_negative
    )

    expect_identical(ordinary$negative, af_negative)
    expect_identical(ordinary$background, af_background)
    expect_identical(viability$negative, dead_negative)
    expect_identical(viability$background, dead_background)
    expect_identical(beads$negative, bead_negative)
    expect_identical(beads$background, bead_background)
})

test_that("AutoSpectral internal background prefers a manual histogram negative gate", {
    detector_names <- c("B1-A", "YG1-A")
    gated <- data.frame(
        `FSC-A` = seq_len(100),
        `SSC-A` = seq_len(100) * 2,
        `B1-A` = c(rep(10, 20), rep(1000, 80)),
        `YG1-A` = c(rep(5, 20), rep(100, 80)),
        check.names = FALSE
    )
    manual_gates <- data.frame(
        gate_type = rep("negative", 2), scope = rep("file", 2),
        filename = rep("beads.fcs", 2), x_channel = rep("B1-A", 2),
        y_channel = rep("", 2), plot_mode = rep("negative_1d", 2),
        vertex_index = 1:2, x = c(5, 20), y = c(0, 0),
        stringsAsFactors = FALSE
    )
    manual_gates <- split(
        manual_gates,
        paste(manual_gates$gate_type, manual_gates$scope, manual_gates$filename, sep = "\r")
    )

    negative <- spectreasy:::.reference_histogram_negative_source(
        gated_data = as.matrix(gated), peak_channel = "B1-A",
        filename = "beads.fcs", sample_type = "beads",
        manual_gates = manual_gates, detector_names = detector_names,
        fsc = "FSC-A", ssc = "SSC-A"
    )

    expect_equal(as.numeric(negative), c(10, 5))
    expect_identical(attr(negative, "source"), "manual_negative_gate")
    expect_equal(attr(negative, "scc_background")$spectra[, "B1-A"], rep(10, 20))
})

test_that("spectral SCC processing resolves automatic positive and negative histogram gates", {
    set.seed(12)
    peak <- 10^c(stats::rnorm(500, 2.5, 0.12), stats::rnorm(1000, 5.5, 0.10))
    gated <- data.frame(
        `FSC-A` = seq_along(peak),
        `SSC-A` = seq_along(peak) * 2,
        `B1-A` = peak,
        `YG1-A` = peak * 0.2,
        check.names = FALSE
    )
    config <- list(
        histogram_pct_beads = 0.98,
        histogram_direction_beads = "right",
        histogram_pct_cells = 0.35,
        histogram_direction_cells = "right"
    )

    positive <- spectreasy:::.resolve_reference_positive_histogram_gate(
        gated_data = as.matrix(gated), peak_channel = "B1-A",
        filename = "cells.fcs", sample_type = "cells",
        manual_gates = list(), config = config, row_info = data.frame()
    )
    negative <- spectreasy:::.reference_histogram_negative_source(
        gated_data = as.matrix(gated), peak_channel = "B1-A",
        filename = "cells.fcs", sample_type = "cells",
        manual_gates = list(), detector_names = c("B1-A", "YG1-A"),
        fsc = "FSC-A", ssc = "SSC-A", hist_info = positive$hist_info
    )

    expect_lt(nrow(positive$final_gated_data), nrow(gated))
    expect_true(all(positive$final_gated_data[, "B1-A"] >= positive$hist_info$gate_min))
    expect_true(all(positive$final_gated_data[, "B1-A"] <= positive$hist_info$gate_max))
    expect_identical(attr(negative, "source"), "automatic_negative_gate")
    expect_lt(unname(negative[["B1-A"]]), positive$hist_info$gate_min)
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

test_that("manual positive gates use raw detector coordinates", {
    gate_csv <- tempfile(fileext = ".csv")
    rows <- data.frame(
        gate_type = c("setting", "positive", "positive", "negative", "negative"),
        scope = c("global", "file", "file", "file", "file"),
        filename = c("", "sample.fcs", "sample.fcs", "sample.fcs", "sample.fcs"),
        x_channel = c("histogram_transform", "Peak-A", "Peak-A", "Peak-A", "Peak-A"),
        y_channel = c("", "", "", "", ""),
        plot_mode = c("setting", "positive_1d", "positive_1d", "negative_1d", "negative_1d"),
        vertex_index = c(0, 1, 2, 1, 2),
        x = c("asinh", "100", "200", "1", "20"),
        y = c("", "0", "0", "0", "0"),
        stringsAsFactors = FALSE
    )
    utils::write.csv(rows, gate_csv, row.names = FALSE, quote = TRUE)

    manual_gates <- spectreasy:::.read_reference_manual_gates(gate_csv)
    gated_data <- data.frame(`Peak-A` = 1:300, check.names = FALSE)
    out <- spectreasy:::.apply_reference_manual_positive_gate(
        gated_data = gated_data,
        peak_channel = "Peak-A",
        filename = "sample.fcs",
        sample_type = "beads",
        manual_gates = manual_gates
    )

    expect_equal(range(out$final_gated_data[["Peak-A"]]), c(100L, 200L))
    expect_equal(nrow(out$final_gated_data), 101L)
    expect_equal(out$hist_info$positive_raw_min, 100)
    expect_equal(out$hist_info$positive_raw_max, 200)
    expect_equal(out$hist_info$negative_raw_min, 1)
    expect_equal(out$hist_info$negative_raw_max, 20)
})

test_that("QC histogram ticks are drawn in transform space but labeled as raw values", {
    raw_domain <- c(0, 277000)
    plot_domain <- spectreasy:::.reference_histogram_transform_values(raw_domain, transform = "asinh")
    breaks <- spectreasy:::.reference_gui_ticks(plot_domain, 5L)
    labels <- spectreasy:::.reference_pretty_k_label(
        spectreasy:::.reference_histogram_inverse_values(breaks, transform = "asinh")
    )

    expect_equal(labels[1], "0")
    expect_equal(labels[2], "575")
    expect_true(grepl("K$", labels[5]))
    expect_false(any(labels %in% c("2", "4", "6", "8")))
})

test_that("automatic histogram transforms calibrate independently and preserve raw coordinates", {
    narrow <- c(-100, -60, -20, 0, 200, 1000, 5000, 10000)
    wide <- c(-1000, -600, -200, 0, 20000, 100000, 500000, 1000000)
    narrow_cofactor <- spectreasy:::.reference_histogram_auto_cofactor(narrow)
    wide_cofactor <- spectreasy:::.reference_histogram_auto_cofactor(wide)

    expect_gt(wide_cofactor, narrow_cofactor)
    transformed <- spectreasy:::.reference_histogram_transform_values(
        wide, transform = "auto", cofactor = wide_cofactor
    )
    restored <- spectreasy:::.reference_histogram_inverse_values(
        transformed, transform = "auto", cofactor = wide_cofactor
    )
    expect_equal(restored, wide, tolerance = 1e-6)
    expect_lt(diff(range(transformed)), 20)
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

test_that("derive_unmixing_matrix WLS honors detector noise", {
    M <- matrix(c(
        1.0, 0.2, 0.05,
        0.1, 1.0, 0.25
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B1-A", "YG1-A", "R1-A")

    expect_silent(
        W_wls <- spectreasy::derive_unmixing_matrix(M, method = "WLS")
    )

    expect_equal(dim(W_wls), dim(M))

    M_no_attr <- M
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

test_that("SCC QC plot failures warn without aborting unmixing helpers", {
    bad_plot <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
        ggplot2::annotate("rect", ymin = -Inf, ymax = Inf, alpha = 0.1)

    expect_warning(
        ok <- spectreasy:::.save_reference_ggsave(
            tempfile(fileext = ".png"),
            bad_plot,
            sn = "bad_control",
            plot_type = "bad plot",
            width = 1,
            height = 1,
            dpi = 72
        ),
        "Continuing unmixing"
    )
    expect_false(ok)
})
