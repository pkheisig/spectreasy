testthat::test_that("create_control_file recognizes fluor and control type from filename variants", {
    testthat::skip_if_not_installed("spectreasy")

    scc_dir <- tempfile("spectreasy_scc_")
    dir.create(scc_dir, recursive = TRUE, showWarnings = FALSE)

    files <- c(
        "LIVE DEAD NIR (Cells).fcs",
        "Alexa 647 (Beads).fcs",
        "AF594 (Beads).fcs",
        "AF488 CD4 (Beads).fcs",
        "BB515 CD123 (Beads).fcs",
        "BB700 HLA-DR (Beads).fcs",
        "APC-H7 CD279 (Beads).fcs",
        "APC PD-1 (Beads).fcs",
        "PE PD1 (Beads).fcs",
        "FITC CD39 (Beads).fcs",
        "FITC ENTPD1 (Beads).fcs",
        "BV510 CD1a (Beads).fcs",
        "BV421 CD371 (Beads).fcs",
        "FVD780 CD8 (Cells).fcs",
        "Zombie Aqua CCR7 (Cells).fcs",
        "PE-CF594 (Beads).fcs",
        "PECF594 (beads).fcs",
        "pe cy7 (bEaDs).fcs",
        "PE-Fire 700 (BEADS).fcs",
        "PE (Beads).fcs"
    )
    created <- file.create(file.path(scc_dir, files))
    testthat::expect_true(all(created))

    out_csv <- tempfile(fileext = ".csv")
    df <- spectreasy::create_control_file(
        input_folder = scc_dir,
        include_af_folder = FALSE,
        output_file = out_csv
    )

    by_file <- split(df, df$filename)

    testthat::expect_equal(by_file[["LIVE DEAD NIR (Cells).fcs"]]$fluorophore[[1]], "LIVE/DEAD NIR")
    testthat::expect_equal(by_file[["LIVE DEAD NIR (Cells).fcs"]]$control.type[[1]], "cells")
    testthat::expect_equal(by_file[["LIVE DEAD NIR (Cells).fcs"]]$is.viability[[1]], "TRUE")

    testthat::expect_equal(by_file[["Alexa 647 (Beads).fcs"]]$fluorophore[[1]], "Alexa Fluor 647")
    testthat::expect_equal(by_file[["AF594 (Beads).fcs"]]$fluorophore[[1]], "Alexa Fluor 594")

    testthat::expect_equal(by_file[["AF488 CD4 (Beads).fcs"]]$fluorophore[[1]], "Alexa Fluor 488")
    testthat::expect_equal(by_file[["AF488 CD4 (Beads).fcs"]]$marker[[1]], "CD4")
    testthat::expect_equal(by_file[["AF488 CD4 (Beads).fcs"]]$control.type[[1]], "beads")

    testthat::expect_equal(by_file[["BB515 CD123 (Beads).fcs"]]$fluorophore[[1]], "BB515")
    testthat::expect_equal(by_file[["BB515 CD123 (Beads).fcs"]]$marker[[1]], "CD123")

    testthat::expect_equal(by_file[["BB700 HLA-DR (Beads).fcs"]]$fluorophore[[1]], "BB700")
    testthat::expect_equal(by_file[["BB700 HLA-DR (Beads).fcs"]]$marker[[1]], "HLA-DR")

    testthat::expect_equal(by_file[["APC-H7 CD279 (Beads).fcs"]]$fluorophore[[1]], "APC-H7")
    testthat::expect_equal(by_file[["APC-H7 CD279 (Beads).fcs"]]$marker[[1]], "CD279")

    testthat::expect_equal(by_file[["APC PD-1 (Beads).fcs"]]$marker[[1]], "PD-1")
    testthat::expect_equal(by_file[["PE PD1 (Beads).fcs"]]$marker[[1]], "PD1")
    testthat::expect_equal(by_file[["FITC CD39 (Beads).fcs"]]$marker[[1]], "CD39")
    testthat::expect_equal(by_file[["FITC ENTPD1 (Beads).fcs"]]$marker[[1]], "ENTPD1")
    testthat::expect_equal(by_file[["BV510 CD1a (Beads).fcs"]]$marker[[1]], "CD1a")

    testthat::expect_equal(by_file[["BV421 CD371 (Beads).fcs"]]$marker[[1]], "CD371")

    testthat::expect_equal(by_file[["FVD780 CD8 (Cells).fcs"]]$fluorophore[[1]], "FVD780")
    testthat::expect_equal(by_file[["FVD780 CD8 (Cells).fcs"]]$marker[[1]], "CD8")
    testthat::expect_equal(by_file[["FVD780 CD8 (Cells).fcs"]]$is.viability[[1]], "TRUE")

    testthat::expect_equal(by_file[["Zombie Aqua CCR7 (Cells).fcs"]]$fluorophore[[1]], "Zombie Aqua")
    testthat::expect_equal(by_file[["Zombie Aqua CCR7 (Cells).fcs"]]$marker[[1]], "CCR7")
    testthat::expect_equal(by_file[["Zombie Aqua CCR7 (Cells).fcs"]]$is.viability[[1]], "TRUE")

    testthat::expect_equal(by_file[["PE-CF594 (Beads).fcs"]]$fluorophore[[1]], "PE-CF594")
    testthat::expect_equal(by_file[["PE-CF594 (Beads).fcs"]]$control.type[[1]], "beads")

    testthat::expect_equal(by_file[["PECF594 (beads).fcs"]]$fluorophore[[1]], "PE-CF594")
    testthat::expect_equal(by_file[["PECF594 (beads).fcs"]]$control.type[[1]], "beads")

    testthat::expect_equal(by_file[["pe cy7 (bEaDs).fcs"]]$fluorophore[[1]], "PE-Cy7")
    testthat::expect_equal(by_file[["pe cy7 (bEaDs).fcs"]]$control.type[[1]], "beads")

    testthat::expect_equal(by_file[["PE-Fire 700 (BEADS).fcs"]]$fluorophore[[1]], "PE-Fire 700")
    testthat::expect_equal(by_file[["PE-Fire 700 (BEADS).fcs"]]$control.type[[1]], "beads")

    testthat::expect_equal(by_file[["PE (Beads).fcs"]]$fluorophore[[1]], "PE")
    testthat::expect_equal(by_file[["PE (Beads).fcs"]]$control.type[[1]], "beads")
    testthat::expect_equal(by_file[["PE (Beads).fcs"]]$is.viability[[1]], "")
})

testthat::test_that("detector fallback matches whole detector codes", {
    testthat::skip_if_not_installed("spectreasy")

    infer <- spectreasy:::.infer_fluor_from_detector

    testthat::expect_equal(infer("V1-A", fluor_channel_map = character()), "BV421")
    testthat::expect_equal(infer("V10-A", fluor_channel_map = character()), "BV605")
    testthat::expect_equal(infer("V11-A", fluor_channel_map = character()), "BV650")
    testthat::expect_equal(infer("V13-A", fluor_channel_map = character()), "BV711")
    testthat::expect_equal(infer("UV10-A", fluor_channel_map = character()), "BUV615")
    testthat::expect_equal(infer("UV11-A", fluor_channel_map = character()), "BUV661")
    testthat::expect_equal(infer("UV17-A", fluor_channel_map = character()), "BUV395")
    testthat::expect_equal(infer("R10-A", fluor_channel_map = character()), "APC")
})

testthat::test_that("custom fluorophore overrides accept common filename forms", {
    testthat::skip_if_not_installed("spectreasy")

    scc_dir <- tempfile("spectreasy_scc_")
    dir.create(scc_dir, recursive = TRUE, showWarnings = FALSE)
    testthat::expect_true(file.create(file.path(scc_dir, "odd-control-name.fcs")))

    out_csv <- tempfile(fileext = ".csv")
    df <- spectreasy::create_control_file(
        input_folder = scc_dir,
        include_af_folder = FALSE,
        unknown_fluor_policy = "empty",
        output_file = out_csv,
        custom_fluorophores = c("odd-control-name" = "BUV737")
    )

    testthat::expect_equal(df$fluorophore[[1]], "BUV737")
})
