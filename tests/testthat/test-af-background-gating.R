testthat::test_that("AF background extraction prefers manual scatter gates", {
    testthat::skip_if_not_installed("flowCore")
    fcs_file <- tempfile(fileext = ".fcs")
    raw <- cbind(
        `FSC-A` = seq_len(30),
        `SSC-A` = seq_len(30),
        `FSC-H` = seq_len(30),
        `V1-A` = seq_len(30) * 10
    )
    flowCore::write.FCS(flowCore::flowFrame(raw), fcs_file)

    manual_result <- list(
        gated_data = raw[11:20, , drop = FALSE],
        fsc = "FSC-A",
        ssc = "SSC-A"
    )
    testthat::local_mocked_bindings(
        .apply_reference_manual_scatter_gates = function(...) manual_result,
        .compute_reference_scatter_gate = function(...) stop("automatic gate should not run"),
        .package = "spectreasy"
    )

    extracted <- spectreasy:::.extract_reference_af_gated_events(
        fcs_file = fcs_file,
        detector_names = "V1-A",
        config = list(manual_gates = list()),
        sample_type = "cells"
    )

    testthat::expect_equal(nrow(extracted$events), 10L)
    testthat::expect_equal(as.numeric(extracted$events[, "V1-A"]), raw[11:20, "V1-A"])
    testthat::expect_equal(extracted$source$n_scatter_gated, 10L)
})
