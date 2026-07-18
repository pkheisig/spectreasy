benchmark_test_entry <- function(label, method = "AutoSpectral", enabled = TRUE) {
    spectreasy:::.normalize_benchmark_entry(list(
        entry_id = paste0("entry-", tolower(label)),
        label = label,
        enabled = enabled,
        method = method
    ))
}

test_that("benchmark configuration round-trips linked entry IDs in project scope", {
    project <- tempfile("benchmark_project_")
    dir.create(project)
    config <- spectreasy:::.new_benchmark_config("Comparison", "benchmark-test")
    config$entries <- list(benchmark_test_entry("One", "Spectreasy"), benchmark_test_entry("Two", "AutoSpectral"))

    path <- spectreasy:::save_benchmark_config(config, project)
    loaded <- spectreasy:::load_benchmark_config("benchmark-test", project)

    expect_true(file.exists(path))
    expect_identical(vapply(loaded$entries, `[[`, character(1), "entry_id"), c("entry-one", "entry-two"))
    expect_identical(vapply(loaded$entries, `[[`, character(1), "method"), c("Spectreasy", "AutoSpectral"))
    expect_true(startsWith(normalizePath(path), normalizePath(project)))
    expect_error(spectreasy:::.benchmark_root(project, "../escape"), "must contain")
})

test_that("benchmark validation requires two enabled entries only when running", {
    config <- spectreasy:::.new_benchmark_config(benchmark_id = "benchmark-one")
    config$entries <- list(benchmark_test_entry("One"), benchmark_test_entry("Two", enabled = FALSE))
    expect_silent(spectreasy:::validate_benchmark_config(config))
    expect_error(spectreasy:::validate_benchmark_config(config, require_runnable = TRUE), "At least two")
})

test_that("control benchmark continues after one failure and marks linked samples stale", {
    project <- tempfile("benchmark_control_queue_")
    dir.create(project)
    dir.create(file.path(project, "scc"))
    dir.create(file.path(project, "samples"))
    utils::write.csv(data.frame(filename = "x.fcs", fluorophore = "FITC"), file.path(project, "fcs_mapping.csv"), row.names = FALSE)
    config <- spectreasy:::.new_benchmark_config("Queue", "benchmark-queue")
    config$entries <- list(
        benchmark_test_entry("One", "Spectreasy"),
        benchmark_test_entry("Two", "AutoSpectral"),
        benchmark_test_entry("Three", "OLS")
    )
    config$entries[[1]]$sample_status <- spectreasy:::.benchmark_status("complete")
    spectreasy:::save_benchmark_config(config, project)
    calls <- character()
    runner <- function(output_dir, unmixing_method, ...) {
        calls <<- c(calls, unmixing_method)
        if (identical(unmixing_method, "AutoSpectral")) stop("intentional failure")
        stage <- file.path(output_dir, "unmix_controls")
        dir.create(file.path(stage, "qc_controls"), recursive = TRUE)
        paths <- file.path(stage, c("scc_reference_matrix.csv", "scc_detector_noise.csv", "qc_controls/qc_controls_report.html"))
        file.create(paths)
        list(reference_matrix_file = paths[1], detector_noise_file = paths[2], spectral_variant_library_file = NULL, qc_report_file = paths[3])
    }

    result <- spectreasy:::run_benchmark_controls("benchmark-queue", project, runner = runner)
    loaded <- spectreasy:::load_benchmark_config("benchmark-queue", project)

    expect_identical(calls, c("Spectreasy", "AutoSpectral", "OLS"))
    expect_false(result$success)
    expect_identical(vapply(loaded$entries, function(entry) entry$control_status$state, character(1)), c("complete", "failed", "complete"))
    expect_identical(loaded$entries[[1]]$sample_status$state, "stale")
    expect_identical(loaded$entries[[2]]$sample_status$state, "failed")
    expect_match(loaded$entries[[2]]$sample_status$message, "Blocked")
    expect_match(loaded$entries[[1]]$control_outputs$matrix_file, "benchmarks/benchmark-queue/controls/entry-one", fixed = TRUE)
})

test_that("sample benchmark uses the matching control matrix for every entry", {
    project <- tempfile("benchmark_sample_queue_")
    dir.create(project)
    dir.create(file.path(project, "scc"))
    dir.create(file.path(project, "samples"))
    config <- spectreasy:::.new_benchmark_config("Samples", "benchmark-samples")
    config$entries <- list(benchmark_test_entry("One"), benchmark_test_entry("Two"))
    for (index in seq_along(config$entries)) {
        matrix <- file.path(project, "matrices", paste0("matrix-", index, ".csv"))
        dir.create(dirname(matrix), recursive = TRUE, showWarnings = FALSE)
        file.create(matrix)
        config$entries[[index]]$control_status <- spectreasy:::.benchmark_status("complete")
        config$entries[[index]]$control_outputs$matrix_file <- spectreasy:::.benchmark_relative_path(matrix, project)
    }
    spectreasy:::save_benchmark_config(config, project)
    matrices <- character()
    runner <- function(output_dir, unmixing_matrix_file, ...) {
        matrices <<- c(matrices, basename(unmixing_matrix_file))
        report <- file.path(output_dir, "unmix_samples", "qc_samples", "qc_samples_report.html")
        dir.create(dirname(report), recursive = TRUE)
        file.create(report)
        out <- list(sample = data.frame())
        class(out) <- c("spectreasy_unmixed_results", "list")
        attr(out, "qc_report_file") <- report
        out
    }

    result <- spectreasy:::run_benchmark_samples("benchmark-samples", project, runner = runner)
    expect_true(result$success)
    expect_identical(matrices, c("matrix-1.csv", "matrix-2.csv"))
    expect_true(all(vapply(result$benchmark$entries, function(entry) entry$sample_status$state == "complete", logical(1))))
})
