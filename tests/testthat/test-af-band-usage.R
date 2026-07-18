test_that("AF usage rows use assignment counts and the eligible-event denominator", {
    usage <- spectreasy:::.af_band_usage_rows(
        "sample_1", c("AF_1", "AF_2", "AF_3"), c(2, 3, 5), 10
    )
    expect_identical(usage$assigned_events, c(2L, 3L, 5L))
    expect_equal(usage$usage_fraction, c(0.2, 0.3, 0.5))
    expect_equal(sum(usage$usage_fraction), 1)
})

test_that("AF usage pages paginate samples with one common color scale", {
    usage <- do.call(rbind, lapply(seq_len(31), function(index) {
        spectreasy:::.af_band_usage_rows(
            paste0("sample_", index), c("AF_1", "AF_2"), c(index, 31 - index), 31
        )
    }))
    pages <- spectreasy:::.build_af_band_usage_pages(usage)
    expect_length(pages, 2L)
    expect_equal(
        attr(pages[[1]], "spectreasy_af_usage_limits"),
        attr(pages[[2]], "spectreasy_af_usage_limits")
    )
})

test_that("sample report places AF spectra and usage directly after reference spectra", {
    skip_if_not_installed("base64enc")
    M <- rbind(FITC = c(1, 0.2), PE = c(0.2, 1), AF_1 = c(0.3, 0.1), AF_2 = c(0.1, 0.3))
    colnames(M) <- c("B1-A", "YG1-A")
    results <- list(sample_a = list(
        data = data.frame(FITC = rnorm(20), PE = rnorm(20), AF_1 = runif(20), AF_2 = runif(20), File = "sample_a"),
        residuals = matrix(rnorm(40), ncol = 2, dimnames = list(NULL, colnames(M)))
    ))
    class(results) <- c("spectreasy_unmixed_results", "list")
    attr(results, "method") <- "AutoSpectral"
    attr(results, "af_band_usage") <- spectreasy:::.af_band_usage_rows("sample_a", c("AF_1", "AF_2"), c(8, 12), 20)
    attr(results, "af_band_usage_status") <- "available"
    metrics <- tempfile("af_usage_metrics_")
    report <- collect_sample_report_data(results, M, qc_metrics_dir = metrics, plot_dir = tempfile("af_usage_plots_"))

    expect_true(file.exists(file.path(metrics, "af_band_usage_by_sample.csv")))
    expect_length(report$plot_manifest$af, 1L)
    expect_length(report$plot_manifest$af_usage, 1L)
    sections <- spectreasy:::.report_sections(spectreasy:::.report_embed_plot_manifest(report))
    expect_identical(names(sections)[1:3], c("spectra", "af", "af_usage"))
})

test_that("sample report remains valid when AF usage is unavailable", {
    M <- rbind(FITC = c(1, 0.2), AF_1 = c(0.1, 0.3))
    colnames(M) <- c("B1-A", "YG1-A")
    results <- data.frame(File = "sample", FITC = rnorm(10), AF_1 = runif(10))
    report <- collect_sample_report_data(results, M, unmixing_method = "OLS", plot_dir = tempfile("af_usage_unavailable_"))
    sections <- spectreasy:::.report_sections(report)
    expect_match(sections$af_usage[2], "not available", fixed = TRUE)
})

test_that("sample PDF writes AF usage plots and metrics after the AF bank", {
    M <- rbind(
        FITC = c(1, 0.2, 0.1),
        PE = c(0.1, 1, 0.2),
        AF = c(0.4, 0.2, 0.1),
        AF_2 = c(0.1, 0.3, 0.5)
    )
    colnames(M) <- c("B1-A", "YG1-A", "R1-A")
    results <- data.frame(
        File = rep(c("Sample_01", "Sample_02"), each = 24),
        FITC = stats::rnorm(48),
        PE = stats::rnorm(48),
        AF = stats::runif(48)
    )
    usage <- rbind(
        spectreasy:::.af_band_usage_rows("Sample_01", c("AF", "AF_2"), c(7, 17), 24),
        spectreasy:::.af_band_usage_rows("Sample_02", c("AF", "AF_2"), c(11, 13), 24)
    )
    attr(results, "af_band_usage") <- usage
    attr(results, "af_band_usage_status") <- "available"
    plot_dir <- tempfile("af_usage_pdf_plots_")
    metrics_dir <- tempfile("af_usage_pdf_metrics_")
    output <- tempfile(fileext = ".pdf")

    report <- qc_samples(
        results,
        M = M,
        output_file = output,
        unmixing_method = "AutoSpectral",
        report_format = "pdf",
        qc_plot_dir = plot_dir,
        save_qc_pngs = TRUE,
        qc_metrics_dir = metrics_dir,
        sample_nxn_max_points = 20
    )

    expect_true(file.exists(report$output_file))
    expect_true(file.exists(file.path(plot_dir, "af_bank.png")))
    expect_true(file.exists(file.path(plot_dir, "af_band_usage_01.png")))
    expect_true(file.exists(file.path(metrics_dir, "af_band_usage_by_sample.csv")))
})
