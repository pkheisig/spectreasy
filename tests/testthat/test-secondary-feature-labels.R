testthat::test_that("unmix_samples writes fluorophore channels with mapped marker annotations", {
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
    control_file <- file.path(tmp, "project_specific_bead_controls.csv")
    utils::write.csv(control_df, control_file, row.names = FALSE, quote = TRUE)

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
        control_file = control_file,
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
    marker_by_fluorophore <- stats::setNames(as.character(pd$desc), as.character(pd$name))

    testthat::expect_equal(marker_by_fluorophore[["FITC"]], "CD45RA")
    testthat::expect_equal(marker_by_fluorophore[["PE"]], "CD3")
    testthat::expect_equal(marker_by_fluorophore[["FSC-A"]], "FSC-A")
    testthat::expect_equal(marker_by_fluorophore[["Time"]], "Time")
    testthat::expect_true(all(c("FITC", "PE") %in% pd$name))
})

testthat::test_that("duplicate marker names are preserved exactly", {
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

    testthat::expect_equal(unname(as.character(pd$name)), c("APC", "BUV661", "Time"))
    testthat::expect_equal(unname(as.character(pd$desc)), c("TCR", "TCR", "Time"))

    output_file <- tempfile(fileext = ".fcs")
    flowCore::write.FCS(output_ff, output_file)
    written_ff <- suppressWarnings(flowCore::read.FCS(
        output_file,
        transformation = FALSE,
        truncate_max_range = FALSE,
        alter.names = FALSE
    ))
    written_keywords <- flowCore::keyword(written_ff)

    testthat::expect_identical(written_keywords[["$P1N"]], "APC")
    testthat::expect_identical(written_keywords[["$P2N"]], "BUV661")
    testthat::expect_identical(written_keywords[["$P1S"]], "TCR")
    testthat::expect_identical(written_keywords[["$P2S"]], "TCR")
})

testthat::test_that("output labels use only the control mapping passed by the workflow", {
    tmp <- tempfile("spectreasy_exact_mapping_")
    dir.create(tmp)
    exact_mapping <- file.path(tmp, "arbitrary_control_filename.csv")
    utils::write.csv(
        data.frame(fluorophore = "FITC", marker = "CD8"),
        exact_mapping,
        row.names = FALSE
    )

    labels <- spectreasy:::.resolve_output_marker_map(
        "FITC",
        control_file = exact_mapping
    )
    fallback <- spectreasy:::.resolve_output_marker_map("FITC")

    testthat::expect_identical(unname(labels[["FITC"]]), "CD8")
    testthat::expect_identical(unname(fallback[["FITC"]]), "FITC")
})

testthat::test_that("an invalid supplied control mapping never silently becomes fluorophore labels", {
    invalid_mapping <- data.frame(color = "FITC", target = "CD8")

    testthat::expect_error(
        spectreasy:::.resolve_output_marker_map("FITC", control_df = invalid_mapping),
        "must contain fluorophore and marker columns",
        fixed = TRUE
    )
})

testthat::test_that("control mapping context survives in-memory and saved-matrix workflows", {
    tmp <- tempfile("spectreasy_mapping_context_")
    dir.create(tmp)
    mapping <- data.frame(fluorophore = "FITC", marker = "CD8")
    mapping_used_file <- file.path(tmp, "fcs_mapping_used.csv")
    matrix_file <- file.path(tmp, "scc_reference_matrix.csv")
    utils::write.csv(mapping, mapping_used_file, row.names = FALSE)
    utils::write.csv(
        data.frame(fluorophore = "FITC", `B1-A` = 1, check.names = FALSE),
        matrix_file,
        row.names = FALSE
    )

    M <- matrix(1, nrow = 1, dimnames = list("FITC", "B1-A"))
    attr(M, "spectreasy_control_df") <- mapping
    in_memory <- spectreasy:::.resolve_unmix_output_control_context(M)

    plain_M <- matrix(1, nrow = 1, dimnames = list("FITC", "B1-A"))
    from_disk <- spectreasy:::.resolve_unmix_output_control_context(
        plain_M,
        unmixing_matrix_file = matrix_file
    )

    testthat::expect_identical(in_memory$control_df, mapping)
    testthat::expect_null(in_memory$control_file)
    testthat::expect_identical(from_disk$control_file, mapping_used_file)
    testthat::expect_null(from_disk$control_df)
})
