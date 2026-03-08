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
