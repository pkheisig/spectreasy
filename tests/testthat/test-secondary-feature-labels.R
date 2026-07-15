testthat::test_that("unmix_samples writes markers as primary and fluorophores as secondary names", {
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
        file.path(output_dir, "unmix_samples", "unmixed_fcs", "sample1_OLS-0AF.fcs"),
        transformation = FALSE,
        truncate_max_range = FALSE
    )
    pd <- flowCore::pData(flowCore::parameters(out_ff))
    secondary_by_primary <- stats::setNames(as.character(pd$desc), as.character(pd$name))

    testthat::expect_equal(secondary_by_primary[["CD45RA"]], "FITC")
    testthat::expect_equal(secondary_by_primary[["CD3"]], "PE")
    testthat::expect_equal(secondary_by_primary[["FSC-A"]], "FSC-A")
    testthat::expect_equal(secondary_by_primary[["Time"]], "Time")
    testthat::expect_false(any(c("FITC", "PE") %in% pd$name))
})

testthat::test_that("duplicate marker names are disambiguated without losing fluorophore labels", {
    exprs <- matrix(
        seq_len(9),
        nrow = 3,
        dimnames = list(NULL, c("APC", "BUV661", "Time"))
    )
    target_ff <- flowCore::flowFrame(exprs)
    marker_map <- c(APC = "TCR", BUV661 = "TCR")

    output_ff <- spectreasy:::.apply_output_fcs_feature_labels(
        target_ff = target_ff,
        source_ff = target_ff,
        fluorophore_cols = c("APC", "BUV661"),
        marker_label_map = marker_map
    )
    pd <- flowCore::pData(flowCore::parameters(output_ff))

    testthat::expect_equal(unname(as.character(pd$name)), c("TCR", "TCR_1", "Time"))
    testthat::expect_equal(unname(as.character(pd$desc)), c("APC", "BUV661", "Time"))
})

testthat::test_that("sample-project mapping takes precedence over a working-directory mapping", {
    tmp <- tempfile("spectreasy_mapping_precedence_")
    project_dir <- file.path(tmp, "project")
    sample_dir <- file.path(project_dir, "samples")
    dir.create(sample_dir, recursive = TRUE)

    wrong_mapping <- data.frame(fluorophore = "FITC", marker = "Wrong marker")
    project_mapping <- data.frame(fluorophore = "FITC", marker = "CD8")
    utils::write.csv(wrong_mapping, file.path(tmp, "fcs_mapping.csv"), row.names = FALSE)
    utils::write.csv(project_mapping, file.path(project_dir, "fcs_mapping.csv"), row.names = FALSE)

    old_wd <- setwd(tmp)
    on.exit(setwd(old_wd), add = TRUE)
    labels <- spectreasy:::.resolve_output_marker_map("FITC", sample_dir = sample_dir)

    testthat::expect_identical(unname(labels[["FITC"]]), "CD8")
})
