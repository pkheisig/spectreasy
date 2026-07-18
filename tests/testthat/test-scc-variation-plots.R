test_that("SCC variation uses normalized raw positive events and matrix detector order", {
    M <- rbind(BUV805 = c(2, 1, 0.2), PE = c(0.1, 2, 1), AF_1 = c(0.2, 0.2, 0.2))
    colnames(M) <- c("R1-A", "UV2-A", "B1-A")
    events <- rbind(
        c(4, 2, 0),
        c(2, 2, 1),
        c(8, 1, 2),
        matrix(rep(c(4, 2, 1), 18), ncol = 3, byrow = TRUE),
        c(NA, -3, 0)
    )
    colnames(events) <- colnames(M)
    attr(M, "scc_positive_events") <- list(BUV805 = events)

    variation <- spectreasy:::.collect_scc_variation_data(M)
    expect_identical(variation$metadata$fluorophore, c("BUV805", "PE"))
    expect_identical(variation$metadata$status, c("available", "missing_positive_events"))
    expect_identical(variation$detector_order, c("UV2-A", "B1-A", "R1-A"))
    expect_identical(variation$data$BUV805$detector, variation$detector_order)
    expect_equal(max(variation$data$BUV805$reference), 1)
    expect_equal(variation$data$BUV805$reference, c(0.5, 0.1, 1))
    expect_true(all(variation$data$BUV805$lower <= variation$data$BUV805$upper))
    expect_equal(variation$metadata$event_count_used[1], 21L)
})

test_that("SCC variation records factual omission statuses", {
    M <- rbind(FITC = c(1, 0.2), PE = c(0.2, 1), APC = c(1, 0.5))
    colnames(M) <- c("B1-A", "R1-A")
    mismatch <- matrix(1, nrow = 25, ncol = 1, dimnames = list(NULL, "B1-A"))
    insufficient <- matrix(1, nrow = 19, ncol = 2, dimnames = list(NULL, colnames(M)))
    invalid <- matrix(c(NA, -1), nrow = 25, ncol = 2, byrow = TRUE, dimnames = list(NULL, colnames(M)))
    attr(M, "scc_positive_events") <- list(FITC = mismatch, PE = insufficient, APC = invalid)

    variation <- spectreasy:::.collect_scc_variation_data(M)
    expect_identical(
        variation$metadata$status,
        c("detector_mismatch", "insufficient_events", "invalid_event_spectra")
    )
    expect_length(variation$data, 0)
})

test_that("control report embeds an SCC variation gallery without failing omitted fluorophores", {
    skip_if_not_installed("base64enc")
    M <- rbind(FITC = c(1, 0.2), PE = c(0.2, 1))
    colnames(M) <- c("B1-A", "YG1-A")
    events <- matrix(rep(c(1, 0.2), 25), ncol = 2, byrow = TRUE, dimnames = list(NULL, colnames(M)))
    attr(M, "scc_positive_events") <- list(FITC = events)

    report <- collect_control_report_data(
        M,
        scc_dir = tempfile("missing_scc_"),
        control_file = tempfile(fileext = ".csv"),
        plot_dir = tempfile("scc_variation_report_")
    )
    expect_identical(report$scc_variation_metadata$status, c("available", "missing_positive_events"))
    expect_length(report$plot_manifest$scc_variation, 1L)
    embedded <- spectreasy:::.report_embed_plot_manifest(report)
    expect_match(unname(embedded$plot_manifest$scc_variation), "data:image/png;base64,", fixed = TRUE)
    sections <- spectreasy:::.report_sections(embedded)
    expect_true("scc_variation" %in% names(sections))
    expect_lt(match("scc_variation", names(sections)), match("spectra", names(sections)))
    expect_match(sections$scc_variation[2], "10th–90th percentile", fixed = TRUE)
})

test_that("control PDF renders the complete SCC variation gallery", {
    skip_if_not_installed("png")
    M <- rbind(FITC = c(2, 1), PE = c(0.2, 4), APC = c(1, 3))
    colnames(M) <- c("B1-A", "YG1-A")
    events <- matrix(rep(c(2, 1), 24), ncol = 2, byrow = TRUE)
    colnames(events) <- colnames(M)
    attr(M, "scc_positive_events") <- list(FITC = events, PE = events[, 2:1, drop = FALSE], APC = events)
    plot_dir <- tempfile("scc_pdf_plots_")
    dir.create(plot_dir, recursive = TRUE)
    output <- tempfile(fileext = ".pdf")
    one_control <- function() list(
        data = data.frame(FITC = stats::rnorm(30), PE = stats::rnorm(30), APC = stats::rnorm(30)),
        residuals = matrix(stats::rnorm(60), ncol = 2, dimnames = list(NULL, colnames(M)))
    )
    unmixed <- list(FITC = one_control(), PE = one_control(), APC = one_control())
    class(unmixed) <- c("spectreasy_unmixed_results", "list")

    result <- qc_controls(
        M = M,
        output_file = output,
        qc_summary = data.frame(sample = character(), fluorophore = character()),
        report_plot_dir = plot_dir,
        pd = data.frame(name = colnames(M), desc = colnames(M)),
        unmixed_list = unmixed,
        unmixing_method = "NNLS",
        save_qc_pngs = TRUE,
        report_format = "pdf"
    )

    expect_true(file.exists(result$output_file))
    expect_gt(file.info(result$output_file)$size, 0)
    expect_length(list.files(file.path(plot_dir, "scc_variation"), pattern = "\\.png$"), 3L)
})
