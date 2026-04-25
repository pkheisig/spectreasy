test_that("validate_control_file_mapping can ignore unstained files when exclude_af = TRUE", {
    src_fitc <- testthat::test_path("../../scc/FITC (Beads).fcs")
    src_unstained <- testthat::test_path("../../scc/Unstained (Cells).fcs")

    testthat::skip_if_not(file.exists(src_fitc))
    testthat::skip_if_not(file.exists(src_unstained))

    scc_dir <- tempfile("spectreasy_af_exclude_")
    dir.create(scc_dir, recursive = TRUE, showWarnings = FALSE)
    file.copy(src_fitc, file.path(scc_dir, "FITC (Beads).fcs"))
    file.copy(src_unstained, file.path(scc_dir, "Unstained (Cells).fcs"))

    control_df <- data.frame(
        filename = "FITC (Beads).fcs",
        fluorophore = "FITC",
        channel = "",
        stringsAsFactors = FALSE
    )

    preflight_default <- spectreasy::validate_control_file_mapping(
        control_df = control_df,
        scc_dir = scc_dir,
        require_channels = FALSE,
        exclude_af = FALSE
    )
    testthat::expect_false(preflight_default$ok)
    testthat::expect_true(any(grepl("SCC files missing from control file", preflight_default$errors)))

    preflight_excluding_af <- spectreasy::validate_control_file_mapping(
        control_df = control_df,
        scc_dir = scc_dir,
        require_channels = FALSE,
        exclude_af = TRUE
    )
    testthat::expect_true(preflight_excluding_af$ok)
})

test_that("build_reference_matrix omits AF basis rows when exclude_af = TRUE", {
    src_fitc <- testthat::test_path("../../scc/FITC (Beads).fcs")
    src_unstained <- testthat::test_path("../../scc/Unstained (Cells).fcs")

    testthat::skip_if_not(file.exists(src_fitc))
    testthat::skip_if_not(file.exists(src_unstained))

    scc_dir <- tempfile("spectreasy_build_af_exclude_")
    out_dir <- tempfile("spectreasy_build_af_exclude_out_")
    dir.create(scc_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    file.copy(src_fitc, file.path(scc_dir, "FITC (Beads).fcs"))
    file.copy(src_unstained, file.path(scc_dir, "Unstained (Cells).fcs"))

    control_df <- data.frame(
        filename = c("FITC (Beads).fcs", "Unstained (Cells).fcs"),
        fluorophore = c("FITC", "AF"),
        channel = c("", ""),
        marker = c("", "Autofluorescence"),
        stringsAsFactors = FALSE
    )

    M <- spectreasy::build_reference_matrix(
        input_folder = scc_dir,
        output_folder = out_dir,
        save_qc_plots = FALSE,
        control_df = control_df,
        exclude_af = TRUE,
        seed = 1,
        subsample_n = 1000
    )

    testthat::expect_true("FITC" %in% rownames(M))
    testthat::expect_false(any(grepl("^AF($|_)", rownames(M), ignore.case = TRUE)))
})

test_that("autounmix_controls exposes exclude_af", {
    testthat::expect_true("exclude_af" %in% names(formals(spectreasy::autounmix_controls)))
})
