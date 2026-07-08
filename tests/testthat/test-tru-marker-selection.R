make_tru_test_matrix <- function() {
    M <- rbind(
        FITC = c(1.00, 0.08, 0.02),
        PE = c(0.10, 1.00, 0.04),
        AF = c(0.35, 0.20, 1.00),
        AF_2 = c(0.08, 0.45, 0.90)
    )
    colnames(M) <- c("B1-A", "YG1-A", "R1-A")
    M
}

make_tru_null_events <- function(M, n = 700) {
    af_names <- rownames(M)[grepl("^AF($|_)", rownames(M), ignore.case = TRUE)]
    coeff <- matrix(0, nrow = n, ncol = nrow(M), dimnames = list(NULL, rownames(M)))
    sampled_af <- sample(af_names, n, replace = TRUE)
    coeff[cbind(seq_len(n), match(sampled_af, rownames(M)))] <- stats::runif(n, 20, 80)
    out <- coeff %*% M + matrix(stats::rnorm(n * ncol(M), sd = 0.02), ncol = ncol(M))
    colnames(out) <- colnames(M)
    out
}

test_that("TRU_Spectreasy maps inactive markers with UCM", {
    set.seed(801)
    M <- make_tru_test_matrix()
    null_events <- make_tru_null_events(M, n = 700)
    cal <- spectreasy:::.build_tru_calibration(
        M = M,
        unstained_events = null_events,
        min_null_events = 100L,
        ns_correction = FALSE
    )
    expect_s3_class(cal, "spectreasy_tru_calibration")

    exprs <- rbind(
        c(FITC = 300, PE = 0, AF = 50, AF_2 = 0) %*% M,
        c(FITC = 0, PE = 260, AF = 0, AF_2 = 55) %*% M,
        c(FITC = 0, PE = 0, AF = 60, AF_2 = 0) %*% M
    )
    colnames(exprs) <- colnames(M)
    ff <- flowCore::flowFrame(exprs)

    res <- spectreasy::calc_residuals(
        ff,
        M,
        method = "TRU_Spectreasy",
        tru_calibration = cal,
        return_residuals = TRUE,
        write_tru_diagnostics = TRUE
    )

    expect_true("tru_info" %in% names(res))
    expect_true("TRU Active Count" %in% colnames(res$data))
    null_abundances <- as.matrix(cal$null_abundances)
    expect_equal(res$data$PE[1], unname(null_abundances[1, "PE"]), tolerance = 1e-6)
    expect_equal(res$data$FITC[2], unname(null_abundances[2, "FITC"]), tolerance = 1e-6)
    expect_equal(res$data$FITC[3], unname(null_abundances[3, "FITC"]), tolerance = 1e-6)
    expect_equal(res$data$PE[3], unname(null_abundances[3, "PE"]), tolerance = 1e-6)
    expect_true(res$data$FITC[1] > 100)
    expect_true(res$data$PE[2] > 100)
})

test_that("unmix_samples resolves sibling TRU calibration file", {
    set.seed(802)
    M <- make_tru_test_matrix()
    null_events <- make_tru_null_events(M, n = 700)
    cal <- spectreasy:::.build_tru_calibration(
        M = M,
        unstained_events = null_events,
        min_null_events = 100L,
        ns_correction = FALSE
    )
    matrix_dir <- tempfile("spectreasy_tru_matrix_")
    sample_dir <- tempfile("spectreasy_tru_samples_")
    dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    matrix_file <- file.path(matrix_dir, "scc_reference_matrix.csv")
    spectreasy:::.save_reference_matrix_csv(M, matrix_file)
    saveRDS(cal, file.path(matrix_dir, "scc_tru_calibration.rds"), version = 3)

    exprs <- c(FITC = 220, PE = 0, AF = 40, AF_2 = 0) %*% M
    colnames(exprs) <- colnames(M)
    flowCore::write.FCS(flowCore::flowFrame(exprs), file.path(sample_dir, "sample.fcs"))

    res <- spectreasy::unmix_samples(
        sample_dir = sample_dir,
        unmixing_matrix_file = matrix_file,
        unmixing_method = "TRU_Spectreasy",
        write_fcs = FALSE,
        save_report = FALSE,
        write_tru_diagnostics = TRUE,
        verbose = FALSE
    )

    expect_s3_class(res, "spectreasy_unmixed_results")
    expect_s3_class(attr(res, "tru_calibration"), "spectreasy_tru_calibration")
    expect_true("tru_info" %in% names(res$sample))
    expect_equal(res$sample$data$PE[1], unname(as.matrix(cal$null_abundances)[1, "PE"]), tolerance = 1e-6)
})
