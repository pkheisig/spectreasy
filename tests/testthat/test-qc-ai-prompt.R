test_that("HTML reports write one compact measurement-only prompt", {
    M <- matrix(c(1, 0.1, 0.2, 1), nrow = 2, byrow = TRUE,
        dimnames = list(c("FITC", "PE"), c("B1-A", "YG1-A")))
    results <- data.frame(
        File = rep("sample.fcs", 20), FITC = rnorm(20), PE = rnorm(20),
        `AF Index` = rep(c(1L, 2L), each = 10), check.names = FALSE
    )
    metrics_dir <- tempfile("qc_metrics_")
    report <- collect_sample_report_data(
        results, M, unmixing_method = "AutoSpectral", qc_metrics_dir = metrics_dir,
        max_events_per_sample = 20, plot_dir = tempfile("plots_")
    )
    path <- tempfile(fileext = ".html")
    rendered <- render_qc_html_report(report, path, overwrite = "overwrite")
    prompt_path <- rendered$ai_qc_prompt_path
    expect_true(file.exists(prompt_path))
    prompt <- paste(readLines(prompt_path, warn = FALSE), collapse = "\n")
    expect_match(prompt, "SPECTREASY NUMERICAL QC INTERPRETATION PROMPT", fixed = TRUE)
    expect_match(prompt, "af_band_usage.csv", fixed = TRUE)
    expect_false(grepl("quality_grade|grade_counts|threshold_profile|not_graded|caution: QC-|Review:", prompt, ignore.case = TRUE))
    expect_lt(file.info(prompt_path)$size, 250000)
    expect_false(dir.exists(file.path(dirname(path), paste0(tools::file_path_sans_ext(basename(path)), "_ai_qc_data"))))
    expect_silent(spectreasy:::.export_report_ai_qc(report, path))
})

test_that("long tables are summarized rather than dumped into the prompt", {
    path <- tempfile(fileext = ".csv")
    utils::write.csv(data.frame(sample = paste0("s", seq_len(500)), value = seq_len(500)), path, row.names = FALSE)
    report_file <- tempfile(fileext = ".html")
    report <- list(
        report_type = "Sample QC", project_path = tempdir(), created_at = Sys.time(),
        version = "test", unmixing_method = "AutoSpectral", cytometer = "test",
        counts = list(samples = 500L), warnings = character(), qc_metric_paths = path
    )
    class(report) <- c("spectreasy_sample_report_data", "spectreasy_report_data", "list")
    prompt <- spectreasy:::.build_report_ai_qc_prompt(report, report_file, path)
    expect_match(prompt, "Numerical distribution summary", fixed = TRUE)
    expect_lt(nchar(prompt), 20000)
})

test_that("prompt export can be disabled without changing the report", {
    M <- matrix(c(1, 0.1, 0.2, 1), nrow = 2, byrow = TRUE,
        dimnames = list(c("FITC", "PE"), c("B1-A", "YG1-A")))
    report <- collect_sample_report_data(
        data.frame(File = "sample.fcs", FITC = 1, PE = 2), M,
        run_settings = list(save_ai_qc = FALSE), plot_dir = tempfile("plots_")
    )
    path <- tempfile(fileext = ".html")
    rendered <- render_qc_html_report(report, path, overwrite = "overwrite")
    expect_true(file.exists(rendered$output_file))
    expect_null(rendered$ai_qc_prompt_path)
    expect_false(file.exists(spectreasy:::.ai_qc_report_companion_path(path)))
})
