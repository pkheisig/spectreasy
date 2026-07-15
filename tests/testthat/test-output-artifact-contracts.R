test_that("saved SCC detector noise contains only persisted columns", {
    detector_noise <- data.frame(
        detector = c("B1-A", "YG1-A"),
        noise_floor = c(125, 400),
        signal_scale = c(1, 1),
        transient_internal_value = c(2, 3)
    )
    path <- tempfile(fileext = ".csv")

    spectreasy:::.save_detector_noise_csv(detector_noise, path)

    saved <- utils::read.csv(path, check.names = FALSE)
    expect_identical(colnames(saved), c("detector", "noise_floor"))
    expect_equal(saved$detector, detector_noise$detector)
    expect_equal(saved$noise_floor, detector_noise$noise_floor)
})

test_that("automatic control report defaults inside the qc_controls directory", {
    output_dir <- tempfile("spectreasy_control_outputs_")
    paths <- spectreasy:::.unmix_output_paths(output_dir)

    expect_identical(
        paths$qc_report_html,
        file.path(output_dir, "qc_controls", "qc_controls_report.html")
    )
})
