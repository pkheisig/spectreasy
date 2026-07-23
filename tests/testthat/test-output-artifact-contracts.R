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
    expect_identical(
        paths$control_mapping_csv,
        file.path(output_dir, "fcs_mapping_used.csv")
    )
})

test_that("control-stage collision policy versions the complete stage directory", {
    output_root <- tempfile("spectreasy_control_stage_")
    canonical <- file.path(output_root, "unmix_controls")

    expect_identical(
        spectreasy:::.resolve_unmix_controls_stage_dir(output_root, "version"),
        canonical
    )
    dir.create(canonical, recursive = TRUE)
    versioned <- spectreasy:::.resolve_unmix_controls_stage_dir(output_root, "version")
    expect_identical(
        versioned,
        paste0(canonical, "_2")
    )
    expect_true(all(startsWith(
        unlist(spectreasy:::.unmix_output_paths(versioned), use.names = FALSE),
        paste0(versioned, .Platform$file.sep)
    )))
    dir.create(paste0(canonical, "_2"))
    expect_identical(
        spectreasy:::.resolve_unmix_controls_stage_dir(output_root, "version"),
        paste0(canonical, "_3")
    )
    expect_identical(
        spectreasy:::.resolve_unmix_controls_stage_dir(output_root, "overwrite"),
        canonical
    )
    expect_error(
        spectreasy:::.resolve_unmix_controls_stage_dir(output_root, "error"),
        "Control-stage output already exists"
    )
})
