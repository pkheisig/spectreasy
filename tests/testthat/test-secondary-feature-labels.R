testthat::test_that("unmix_samples writes secondary feature names from control marker mapping", {
    testthat::skip_if_not_installed("spectreasy")

    tmp <- tempfile("spectreasy_secondary_labels_")
    dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(tmp)

    control_df <- data.frame(
        filename = c("sample1.fcs", "sample2.fcs"),
        fluorophore = c("FITC", "PE"),
        marker = c("CD45RA", "CD3"),
        channel = c("B1-A", "YG1-A"),
        stringsAsFactors = FALSE
    )
    utils::write.csv(control_df, "fcs_mapping.csv", row.names = FALSE, quote = TRUE)

    sample_dir <- file.path(tmp, "samples")
    output_dir <- file.path(tmp, "out")
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    exprs <- matrix(c(
        100, 20, 1, 1000, 1200, 500, 50,
         90, 30, 2, 1100, 1300, 550, 55,
         80, 40, 3,  900, 1250, 530, 60
    ), nrow = 3, byrow = TRUE)
    colnames(exprs) <- c("B1-A", "YG1-A", "Time", "FSC-A", "FSC-H", "SSC-A", "SSC-W")
    ff <- flowCore::flowFrame(exprs)
    flowCore::write.FCS(ff, file.path(sample_dir, "sample1.fcs"))

    M <- matrix(c(
        1, 0.2,
        0.1, 1
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B1-A", "YG1-A")

    spectreasy::unmix_samples(
        sample_dir = sample_dir,
        M = M,
        unmixing_method = "OLS",
        output_dir = output_dir,
        write_fcs = TRUE
    )

    out_ff <- flowCore::read.FCS(
        file.path(output_dir, "sample1_OLS-0AF.fcs"),
        transformation = FALSE,
        truncate_max_range = FALSE
    )
    pd <- flowCore::pData(flowCore::parameters(out_ff))
    desc_map <- stats::setNames(as.character(pd$desc), as.character(pd$name))

    testthat::expect_equal(desc_map[["FITC"]], "CD45RA")
    testthat::expect_equal(desc_map[["PE"]], "CD3")
    testthat::expect_equal(desc_map[["FSC-A"]], "FSC-A")
    testthat::expect_equal(desc_map[["Time"]], "Time")
})
