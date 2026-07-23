test_that("control metrics use explicit formulas and preserve detector order", {
    set.seed(11)
    negative <- cbind(D1 = rnorm(500, 100, 4), D2 = rnorm(500, 100, 4), D3 = rnorm(500, 100, 4))
    positive <- negative + cbind(D1 = rnorm(500, 900, 10), D2 = rnorm(500, 180, 8), D3 = rnorm(500, 40, 5))
    qc <- calculate_control_qc_metrics(
        positive, negative, expected_peak = "D1",
        detector_ranges = c(D1 = 4095, D2 = 4095, D3 = 4095)
    )
    expect_identical(qc$observed_peak, "D1")
    expect_identical(qc$spectrum_variability$detector, colnames(positive))
    values <- setNames(qc$metrics$value, qc$metrics$metric_id)
    expect_gt(values[["robust_separation"]], 20)
    expect_lt(values[["overlap_misclassification_fraction"]], 0.01)
    expect_equal(values[["peak_agreement"]], 1)
    expect_true(all(c(
        "expected_peak_upper_boundary_fraction",
        "any_detector_upper_boundary_fraction"
    ) %in% names(values)))
    expect_true(all(c("q10", "q25", "median", "q75", "q90") %in% names(qc$spectrum_variability)))
})

test_that("clipping metrics use recorded detector ranges and identify detectors", {
    negative <- cbind(D1 = rep(10, 10), D2 = rep(10, 10))
    positive <- cbind(D1 = c(rep(100, 8), 255, 255), D2 = c(0, rep(40, 9)))
    qc <- calculate_control_qc_metrics(
        positive, negative, expected_peak = "D1",
        detector_ranges = c(D1 = 255, D2 = 255)
    )
    values <- setNames(qc$metrics$value, qc$metrics$metric_id)
    expect_equal(values[["expected_peak_upper_boundary_fraction"]], 0.2)
    expect_equal(values[["any_detector_upper_boundary_fraction"]], 0.2)
    expect_identical(qc$clipping_detectors, "D1")
})

test_that("controlled weak stain and contamination affect intended metrics", {
    set.seed(12)
    negative <- cbind(D1 = rnorm(600, 50, 5), D2 = rnorm(600, 50, 5))
    signal <- cbind(D1 = rnorm(600, 600, 20), D2 = rnorm(600, 80, 8))
    baseline <- calculate_control_qc_metrics(negative + signal, negative, expected_peak = "D1")
    weak <- calculate_control_qc_metrics(negative + signal * 0.25, negative, expected_peak = "D1")
    get <- function(x, id) setNames(x$metrics$value, x$metrics$metric_id)[[id]]
    expect_lt(get(weak, "robust_separation"), get(baseline, "robust_separation"))
    expect_lt(get(weak, "raw_positive_minus_negative"), get(baseline, "raw_positive_minus_negative"))
})

test_that("AF diagnostics expose redundancy and occupancy", {
    af <- rbind(AF_1 = c(1, 0.3, 0.1), AF_2 = c(1.01, 0.29, 0.1), AF_3 = c(0.1, 0.2, 1))
    colnames(af) <- c("D1", "D2", "D3")
    qc <- calculate_af_bank_qc_metrics(af, assignments = c(rep(1, 90), rep(3, 10)), requested_bands = 3)
    values <- setNames(qc$metrics$value, qc$metrics$metric_id)
    expect_gt(values[["maximum_pairwise_similarity"]], 0.99)
    expect_equal(values[["invalid_assignment_count"]], 0)
    expect_equal(qc$occupancy$fraction, c(0.9, 0, 0.1))
})

test_that("AF diagnostics exclude and count invalid assignments", {
    af <- rbind(AF_1 = c(1, 0), AF_2 = c(0, 1))
    colnames(af) <- c("D1", "D2")

    invalid <- calculate_af_bank_qc_metrics(af, assignments = c(0, 3, NA, "unknown"))
    invalid_values <- setNames(invalid$metrics$value, invalid$metrics$metric_id)
    expect_true(all(is.na(invalid$occupancy$fraction)))
    expect_equal(invalid_values[["assignment_count"]], 4)
    expect_equal(invalid_values[["valid_assignment_count"]], 0)
    expect_equal(invalid_values[["invalid_assignment_count"]], 4)
    expect_equal(invalid_values[["invalid_assignment_fraction"]], 1)

    mixed <- calculate_af_bank_qc_metrics(af, assignments = c(1, "AF_2", 0, 3))
    mixed_values <- setNames(mixed$metrics$value, mixed$metrics$metric_id)
    expect_equal(mixed$occupancy$fraction, c(0.5, 0.5))
    expect_equal(mixed_values[["valid_assignment_count"]], 2)
    expect_equal(mixed_values[["invalid_assignment_count"]], 2)
    expect_equal(mixed_values[["invalid_assignment_fraction"]], 0.5)
})

test_that("AF occupancy is unavailable rather than falsely empty without assignments", {
    af <- rbind(AF_1 = c(1, 0.2), AF_2 = c(0.2, 1))
    colnames(af) <- c("D1", "D2")
    qc <- calculate_af_bank_qc_metrics(af)
    expect_true(all(is.na(qc$occupancy$fraction)))
})

test_that("control report adds only consolidated numerical QC tables", {
    M <- matrix(c(1, 0.2, 0.1, 0.4, 0.8, 0.2), nrow = 2, byrow = TRUE,
        dimnames = list(c("FITC", "AF_1"), c("D1", "D2", "D3")))
    set.seed(13)
    negative <- cbind(D1 = rnorm(100, 100, 5), D2 = rnorm(100, 100, 5), D3 = rnorm(100, 100, 5))
    positive <- negative + cbind(D1 = rnorm(100, 900, 10), D2 = rnorm(100, 150, 8), D3 = rnorm(100, 30, 5))
    attr(M, "scc_qc_positive_events") <- list(FITC = positive)
    attr(M, "scc_qc_negative_events") <- list(FITC = negative)
    attr(M, "detector_pd") <- data.frame(name = c("D1", "D2", "D3"), maxRange = rep(4095, 3))
    metrics_dir <- tempfile("control_metrics_")
    report <- collect_control_report_data(
        M, scc_dir = tempfile(), control_file = tempfile(),
        qc_summary = data.frame(fluorophore = "FITC", peak_channel = "D1"),
        qc_metrics_dir = metrics_dir, plot_dir = tempfile("control_plots_")
    )
    expect_setequal(basename(list.files(metrics_dir)), c(
        "reference_spectra.csv", "af_bank_spectra.csv", "control_qc_summary.csv",
        "control_signal_metrics.csv", "control_spectrum_variability.csv", "af_bank_summary.csv"
    ))
    expect_equal(nrow(report$control_signal_metrics), 1L)
    expect_false(any(dir.exists(list.files(metrics_dir, full.names = TRUE))))
})
