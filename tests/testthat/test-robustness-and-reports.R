test_that("with_optional_seed creates deterministic local draws", {
    val1 <- local({
        spectreasy:::.with_optional_seed(999)
        runif(3)
    })
    val2 <- local({
        spectreasy:::.with_optional_seed(999)
        runif(3)
    })

    expect_equal(val1, val2)
})

test_that("derive_unmixing_matrix fails for singular OLS matrix", {
    M <- matrix(c(
        1, 2, 3,
        2, 4, 6
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("A", "B")
    colnames(M) <- c("UV1-A", "UV2-A", "UV3-A")

    expect_error(
        spectreasy::derive_unmixing_matrix(M, method = "OLS"),
        regexp = "singular"
    )
})

test_that("calc_residuals fails when detector columns are missing", {
    M <- matrix(c(
        1, 0.2,
        0.1, 1
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B1-A", "YG1-A")

    exprs <- matrix(c(100, 1, 1000, 500), nrow = 1)
    colnames(exprs) <- c("B1-A", "Time", "FSC-A", "SSC-A")
    ff <- flowCore::flowFrame(exprs)

    expect_error(
        spectreasy::calc_residuals(ff, M, method = "OLS"),
        regexp = "Detectors in reference matrix not found"
    )
})

test_that("calculate_nps excludes AF markers by default", {
    df <- data.frame(
        File = rep(c("A", "B"), each = 30),
        FITC = rnorm(60, sd = 0.2),
        AF = rnorm(60, sd = 0.2),
        AF_2 = rnorm(60, sd = 0.2),
        check.names = FALSE
    )

    nps <- spectreasy::calculate_nps(df)
    expect_true("FITC" %in% nps$Marker)
    expect_false(any(grepl("^AF($|_)", nps$Marker, ignore.case = TRUE)))
})

test_that("generate_qc_report has stable pages and no recommendation page", {
    skip_if_not_installed("pdftools")

    set.seed(1)
    n <- 120
    results_df <- data.frame(
        FITC = rnorm(n, 0, 0.3),
        PE = rnorm(n, 0, 0.4),
        AF = rnorm(n, 0, 0.2),
        File = rep(c("SampleA", "SampleB"), each = n / 2),
        check.names = FALSE
    )

    M <- matrix(c(
        1.0, 0.2, 0.1,
        0.1, 1.0, 0.2,
        0.2, 0.2, 1.0
    ), nrow = 3, byrow = TRUE)
    rownames(M) <- c("FITC", "PE", "AF")
    colnames(M) <- c("B1-A", "YG1-A", "R1-A")

    pdf_out <- tempfile(fileext = ".pdf")
    spectreasy::generate_qc_report(results_df = results_df, M = M, output_file = pdf_out)

    expect_true(file.exists(pdf_out))

    info <- pdftools::pdf_info(pdf_out)
    expect_equal(info$pages, 4)

    txt <- paste(pdftools::pdf_text(pdf_out), collapse = "\n")
    expect_false(grepl("Conclusions & Recommendations", txt, fixed = TRUE))
})
