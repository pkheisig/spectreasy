test_that("sample HTML report is self-contained and contains major sections", {
    skip_if_not_installed("base64enc")
    set.seed(11)
    M <- rbind(
        FITC = c(1, 0.2, 0.05),
        PE = c(0.1, 1, 0.2),
        AF_1 = c(0.15, 0.12, 0.1)
    )
    colnames(M) <- c("B1-A", "YG1-A", "R1-A")
    results <- data.frame(
        File = rep(c("sample_a", "sample_b"), each = 24),
        FITC = stats::rnorm(48),
        PE = stats::rnorm(48),
        AF_1 = stats::runif(48)
    )
    path <- tempfile(fileext = ".html")
    report <- qc_samples(
        results,
        M = M,
        output_file = path,
        output_format = "html",
        sample_nxn_rows_per_page = 2,
        sample_nxn_max_points = 20
    )
    html <- paste(readLines(report$output_file, warn = FALSE), collapse = "\n")
    expect_true(file.exists(report$output_file))
    expect_true(report$self_contained)
    expect_match(html, "data:image/png;base64,", fixed = TRUE)
    expect_false(grepl("QC Summary", html, fixed = TRUE))
    expect_match(html, "Input samples", fixed = TRUE)
    expect_false(grepl("Detector residual summary", html, fixed = TRUE))
    expect_match(html, "Per-sample marker plots", fixed = TRUE)
    expect_match(html, "Generated artifacts", fixed = TRUE)
    expect_match(html, "reference spectra", ignore.case = TRUE)
    expect_false(grepl("Warnings and logs", html, fixed = TRUE))
    expect_false(grepl(">Settings<", html, fixed = TRUE))
})

test_that("sample detector residual table includes every available detector", {
    M <- rbind(FITC = c(1, 0.2, 0.05), PE = c(0.1, 1, 0.2))
    colnames(M) <- c("B1-A", "YG1-A", "R1-A")
    sample_data <- data.frame(FITC = stats::rnorm(20), PE = stats::rnorm(20))
    residuals <- matrix(stats::rnorm(60), nrow = 20, dimnames = list(NULL, colnames(M)))
    results <- list(sample = list(data = sample_data, residuals = residuals))
    class(results) <- c("spectreasy_unmixed_results", "list")
    attr(results, "method") <- "OLS"
    report_data <- collect_sample_report_data(results, M, unmixing_method = "OLS", sample_nxn_max_points = 10)
    expect_equal(nrow(report_data$detector_rms), ncol(M))
    expect_setequal(report_data$detector_rms$detector, colnames(M))
})

test_that("metric plots replace redundant raw numeric tables", {
    M <- rbind(FITC = c(1, 0.2, 0.05), PE = c(0.1, 1, 0.2))
    colnames(M) <- c("B1-A", "YG1-A", "R1-A")
    sample_data <- data.frame(FITC = stats::rnorm(20), PE = stats::rnorm(20))
    residuals <- matrix(stats::rnorm(60), nrow = 20, dimnames = list(NULL, colnames(M)))
    results <- list(sample = list(data = sample_data, residuals = residuals))
    class(results) <- c("spectreasy_unmixed_results", "list")
    attr(results, "method") <- "OLS"

    report_data <- collect_sample_report_data(
        results,
        M,
        unmixing_method = "OLS",
        sample_nxn_max_points = 10
    )
    sections <- spectreasy:::.report_sections(report_data)

    expect_false(grepl(">0.2<", sections[["reference"]][2], fixed = TRUE))
    expect_false(grepl("<table", sections[["matrix"]][2], fixed = TRUE))
    expect_false(grepl("<table", sections[["residual"]][2], fixed = TRUE))
    expect_match(sections[["matrix"]][2], "data:image/png;base64,", fixed = TRUE)
    expect_match(sections[["residual"]][2], "data:image/png;base64,", fixed = TRUE)
})

test_that("metric tables remain available when their plots are missing", {
    metric <- matrix(c(1, 0.4, 0.4, 1), nrow = 2)
    colnames(metric) <- rownames(metric) <- c("FITC", "PE")
    html <- spectreasy:::.report_metric_html(metric, tempfile(fileext = ".png"))
    expect_match(html, "<table", fixed = TRUE)
    expect_match(html, "0.4", fixed = TRUE)
})

test_that("control HTML report uses supplied matrix without rebuilding controls", {
    skip_if_not_installed("base64enc")
    M <- rbind(FITC = c(1, 0.1), PE = c(0.2, 1), AF_1 = c(0.1, 0.1))
    colnames(M) <- c("B1-A", "YG1-A")
    missing_scc <- tempfile("missing_scc_")
    path <- tempfile(fileext = ".html")
    report <- qc_controls(
        M = M,
        scc_dir = missing_scc,
        control_file = tempfile(fileext = ".csv"),
        output_file = path,
        output_format = "html",
        unmixing_method = "OLS",
        qc_plot_dir = tempfile("control_html_plots_")
    )
    html <- paste(readLines(report$output_file, warn = FALSE), collapse = "\n")
    expect_true(file.exists(report$output_file))
    expect_false(grepl("Input files and mapping", html, fixed = TRUE))
    expect_false(grepl("Gating summary", html, fixed = TRUE))
    expect_match(html, "AF bank summary", fixed = TRUE)
    expect_false(grepl("SCC unmixing diagnostics", html, fixed = TRUE))
})

test_that("report format and filename extensions are validated", {
    M <- rbind(FITC = c(1, 0.1), PE = c(0.1, 1))
    colnames(M) <- c("B1-A", "YG1-A")
    results <- data.frame(File = rep("sample", 8), FITC = stats::rnorm(8), PE = stats::rnorm(8))
    expect_error(qc_samples(results, M = M, output_file = tempfile(fileext = ".pdf"), output_format = "html"), "conflicts")
    expect_error(qc_samples(results, M = M, output_file = tempfile(fileext = ".html"), output_format = "pdf"), "conflicts")
    expect_error(qc_samples(results, M = M, output_format = "both"), "arg")
    inferred <- qc_samples(results, M = M, output_file = tempfile(fileext = ".html"), sample_nxn_max_points = 8)
    expect_identical(inferred$format, "html")
})

test_that("report output paths must be non-empty scalar strings", {
    invalid_paths <- list(NULL, character(), NA_character_, "", "   ", c("one.html", "two.html"))

    for (path in invalid_paths) {
        expect_error(
            spectreasy:::.report_output_spec(path, output_format = "html"),
            "output_file must be a single non-empty file path",
            fixed = TRUE
        )
    }
})

test_that("HTML report wrappers fail clearly for missing scientific inputs", {
    M <- rbind(FITC = c(1, 0.1), PE = c(0.1, 1))
    colnames(M) <- c("B1-A", "YG1-A")
    results <- data.frame(File = "sample", FITC = 1, PE = 1)
    expect_error(collect_sample_report_data(NULL, M), "No unmixed results")
    expect_error(collect_sample_report_data(results, NULL), "No reference matrix")
    expect_error(
        qc_samples(results, unmixing_matrix_file = tempfile("missing_matrix_"), output_format = "html"),
        "not found"
    )
})

test_that("HTML overwrite policy versions, overwrites, errors, and detects staleness", {
    M <- rbind(FITC = c(1, 0.1), PE = c(0.1, 1))
    colnames(M) <- c("B1-A", "YG1-A")
    results <- data.frame(File = rep("sample", 8), FITC = stats::rnorm(8), PE = stats::rnorm(8))
    target <- tempfile(fileext = ".html")
    first <- qc_samples(results, M = M, output_file = target, output_format = "html", overwrite = "overwrite", sample_nxn_max_points = 8)
    second <- qc_samples(results, M = M, output_file = target, output_format = "html", overwrite = "version", sample_nxn_max_points = 8)
    expect_false(identical(first$output_file, second$output_file))
    expect_match(basename(second$output_file), "_v001\\.html$")
    expect_error(qc_samples(results, M = M, output_file = target, output_format = "html", overwrite = "error", sample_nxn_max_points = 8), "already exists")

    source <- tempfile(fileext = ".csv")
    writeLines("x", source)
    Sys.sleep(1)
    report_file <- tempfile(fileext = ".html")
    writeLines("report", report_file)
    expect_false(spectreasy:::.report_is_stale(report_file, source))
    Sys.sleep(1)
    writeLines("changed", source)
    expect_true(spectreasy:::.report_is_stale(report_file, source))
})
