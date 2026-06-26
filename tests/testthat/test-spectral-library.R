testthat::test_that("packaged panel-builder libraries load by supported cytometer", {
    ids <- spectreasy:::.spectral_panel_libraries()$id

    testthat::expect_equal(ids, c("aurora", "discover", "id7000", "xenith"))

    aurora <- spectreasy:::.load_spectral_library("aurora")
    testthat::expect_true("FITC" %in% rownames(aurora))
    testthat::expect_true("B2-A" %in% colnames(aurora))
    testthat::expect_gt(nrow(aurora), 20)
    testthat::expect_gt(ncol(aurora), 20)
})

testthat::test_that("spectral library detector matching handles display names", {
    discover <- spectreasy:::.load_spectral_library(
        "discover",
        fluorophores = "FITC",
        detectors = c("B2-A", "YG5-A"),
        strict = TRUE
    )

    testthat::expect_equal(rownames(discover), "FITC")
    testthat::expect_equal(colnames(discover), c("B2-A", "YG5-A"))
    testthat::expect_equal(max(discover[1, ]), 1)
})

testthat::test_that("spectral panel data returns spectra, similarity, and complexity", {
    panel <- spectreasy:::.build_spectral_panel_data(
        c("FITC", "PE", "APC"),
        cytometer = "aurora",
        strict = TRUE
    )

    testthat::expect_true(is.matrix(panel$spectra))
    testthat::expect_equal(rownames(panel$spectra), c("FITC", "PE", "APC"))
    testthat::expect_equal(dim(panel$similarity_matrix), c(3L, 3L))
    testthat::expect_equal(names(panel$peak_detectors), c("FITC", "PE", "APC"))
    testthat::expect_true(is.numeric(panel$complexity_index))
    testthat::expect_gte(panel$complexity_index, 1)
})

testthat::test_that("panel payload serializes selected spectra for GUI", {
    payload <- spectreasy:::.spectral_panel_payload(
        cytometer = "id7000",
        fluorophores = c("FITC", "PE")
    )

    testthat::expect_equal(payload$cytometer, "id7000")
    testthat::expect_true(all(c("detector", "laser", "emission", "color") %in% colnames(payload$detectors)))
    testthat::expect_true(all(c("fluorophore", "FITC", "PE") %in% colnames(payload$similarity)))
    testthat::expect_true(all(payload$selected == c("FITC", "PE")))
    testthat::expect_true(is.numeric(payload$complexity_index))
})

testthat::test_that("panel export table and overview PDF are generated", {
    export_table <- spectreasy:::.spectral_panel_export_table(
        fluorophores = c("FITC", "PE"),
        markers = c("CD8", "CD4")
    )

    testthat::expect_equal(colnames(export_table), c("Marker", "Fluorophore"))
    testthat::expect_equal(export_table$Marker, c("CD8", "CD4"))
    testthat::expect_equal(export_table$Fluorophore, c("FITC", "PE"))

    output_file <- tempfile(fileext = ".pdf")
    spectreasy:::.write_spectral_panel_overview_pdf(
        cytometer = "aurora",
        fluorophores = c("FITC", "PE", "APC"),
        markers = c("CD8", "CD4", "CD3"),
        output_file = output_file
    )

    testthat::expect_true(file.exists(output_file))
    testthat::expect_gt(file.info(output_file)$size, 1000)
    sig <- readBin(output_file, what = "raw", n = 4)
    testthat::expect_equal(rawToChar(sig), "%PDF")

    if (requireNamespace("pdftools", quietly = TRUE)) {
        paged_file <- tempfile(fileext = ".pdf")
        spectreasy:::.write_spectral_panel_overview_pdf(
            cytometer = "aurora",
            fluorophores = c("FITC", "PE", "APC", "BV421", "BV510"),
            output_file = paged_file
        )
        testthat::expect_equal(pdftools::pdf_info(paged_file)$pages, 5L)
    }
})

testthat::test_that("panel builder requires the correct ID7000 name", {
    testthat::expect_error(
        spectreasy:::.spectral_panel_payload(cytometer = "id700"),
        regexp = "id7000"
    )
})

testthat::test_that("build_spectral_panel has no preselection parameters", {
    testthat::expect_false("cytometer" %in% names(formals(spectreasy::build_spectral_panel)))
    testthat::expect_false("fluorophores" %in% names(formals(spectreasy::build_spectral_panel)))
    testthat::expect_true(all(c("port", "open_browser", "dev_mode") %in% names(formals(spectreasy::build_spectral_panel))))
})

testthat::test_that("measured spectra compare against packaged library internally", {
    measured <- spectreasy:::.load_spectral_library(
        "aurora",
        fluorophores = c("FITC", "PE"),
        detectors = c("B2-A", "YG1-A", "R1-A"),
        strict = TRUE
    )

    cmp <- spectreasy:::.compare_control_spectra(measured, cytometer = "aurora")
    testthat::expect_equal(cmp$fluorophore, c("FITC", "PE"))
    testthat::expect_true(all(cmp$cosine_similarity > 0.999))
    testthat::expect_true(all(cmp$flag == "ok"))
})
