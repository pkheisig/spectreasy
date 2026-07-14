test_that("enumerated string arguments are matched case-insensitively", {
    expect_identical(
        spectreasy:::.match_arg_ci("FlOwSeT", c("list", "flowSet", "SingleCellExperiment"), "return_type"),
        "flowSet"
    )
    expect_identical(spectreasy:::.normalize_unmix_method("autospectral"), "AutoSpectral")
    expect_identical(
        spectreasy:::.resolve_cytometer_id("CYTEK AURORA", allow_auto = TRUE),
        "aurora"
    )
    expect_identical(
        spectreasy:::.validate_scc_background_args("SCATTER_KNN", 2L)$method,
        "scatter_knn"
    )

    options <- gating_options(
        histogram_direction_beads = "RIGHT",
        histogram_direction_cells = "BoTh"
    )
    expect_identical(options$histogram_direction_beads, "right")
    expect_identical(options$histogram_direction_cells, "both")

    mat <- cbind(a = seq_len(20), b = rep(1, 20))
    expect_identical(
        gate_positive_cells(mat, histogram_direction = "RIGHT"),
        gate_positive_cells(mat, histogram_direction = "right")
    )
})

test_that("report_format replaces output_format and defaults to HTML", {
    report_functions <- list(unmix_controls, unmix_samples, qc_controls, qc_samples)
    for (fun in report_functions) {
        expect_true("report_format" %in% names(formals(fun)))
        expect_false("output_format" %in% names(formals(fun)))
        expect_identical(formals(fun)$report_format, "html")
    }

    target <- tempfile(fileext = ".html")
    spec <- spectreasy:::.report_output_spec(target, report_format = "HTML")
    expect_identical(spec$format, "html")
    expect_identical(spec$path, target)
})
