test_that("sample file vectors reject bad inputs and preserve duplicate basenames", {
    dir_a <- tempfile("sample_files_a_")
    dir_b <- tempfile("sample_files_b_")
    dir.create(dir_a)
    dir.create(dir_b)
    file_a <- file.path(dir_a, "sample.fcs")
    file_b <- file.path(dir_b, "sample.fcs")
    file.create(file_a, file_b)

    entries <- spectreasy:::.prepare_unmix_samples_input(c(file_a, file_b))
    expect_equal(vapply(entries, `[[`, character(1), "sample_name"), c("sample", "sample.1"))

    expect_error(
        spectreasy:::.prepare_unmix_samples_input(c(file_a, file.path(dir_b, "missing.fcs"))),
        "missing or invalid FCS file paths"
    )
    expect_error(
        spectreasy:::.prepare_unmix_samples_input(c(file_a, file_a)),
        "duplicate FCS file paths"
    )
})

test_that("control inclusion matching is case-insensitive and rejects ambiguous stems", {
    files <- c("/tmp/Control One.FCS", "/tmp/Control Two.fcs", "/tmp/unused.fcs")
    control_df <- data.frame(
        filename = c("control one.fcs", "CONTROL TWO"),
        stringsAsFactors = FALSE
    )
    selected <- spectreasy:::.control_validation_select_scc_files(control_df, files)
    expect_equal(selected, files[1:2])

    expect_error(
        spectreasy:::.control_validation_select_scc_files(
            data.frame(filename = "duplicate.fcs"),
            c("/tmp/duplicate.fcs", "/tmp/DUPLICATE.FCS")
        ),
        "Multiple FCS files match"
    )
})

test_that("unmixing validates event data and scalar options before numeric kernels", {
    M <- diag(2)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B1-A", "YG1-A")
    ff_bad <- flowCore::flowFrame(matrix(
        c(Inf, 20), nrow = 1,
        dimnames = list(NULL, colnames(M))
    ))

    expect_error(spectreasy::calc_residuals(matrix(1), M), "flowCore::flowFrame")
    expect_error(spectreasy::calc_residuals(ff_bad, M, method = "OLS"), "non-finite values.*B1-A \\(1\\)")
    expect_error(spectreasy:::.normalize_unmix_chunk_size(1.5), "integer >= 1")
    expect_error(spectreasy:::.normalize_unmix_plot_n_events(c(10, 20)), "length 1")
    expect_error(spectreasy:::.normalize_unmix_threads("2"), "integer >= 1")

    M_reserved <- M
    rownames(M_reserved)[1] <- "Time"
    ff_ok <- flowCore::flowFrame(matrix(
        c(10, 20, 1), nrow = 1,
        dimnames = list(NULL, c(colnames(M), "Time"))
    ))
    expect_error(
        spectreasy::calc_residuals(ff_ok, M_reserved, method = "OLS"),
        "reserved output/acquisition columns: Time"
    )
})

test_that("detector sorting supports parameter data without descriptions", {
    pd <- data.frame(name = c("FSC-A", "V2-A", "V1-A", "Time"), stringsAsFactors = FALSE)
    sorted <- spectreasy::get_sorted_detectors(pd)
    expect_equal(sorted$names, c("V1-A", "V2-A"))
    expect_error(spectreasy::get_sorted_detectors(data.frame(desc = "x")), "'name' column")
})

test_that("chunked FCS assembly preserves all events without unsafe output paths", {
    M <- diag(2)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B1-A", "YG1-A")
    expr <- cbind(
        "B1-A" = seq_len(7),
        "YG1-A" = seq_len(7) * 2,
        "Time" = seq_len(7)
    )
    fs <- flowCore::flowSet(list("run/one" = flowCore::flowFrame(expr)))
    output_dir <- tempfile("chunked_output_")

    out <- spectreasy::unmix_samples(
        sample_dir = fs,
        M = M,
        unmixing_method = "OLS",
        output_dir = output_dir,
        chunk_size = 2L,
        plot_n_events = NULL,
        write_fcs = TRUE,
        save_report = FALSE,
        verbose = FALSE
    )

    expect_equal(nrow(out[[1]]$data), 7L)
    output_files <- list.files(output_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
    expect_length(output_files, 1L)
    expect_false(grepl("run/one", output_files, fixed = TRUE))
    written <- flowCore::read.FCS(output_files[[1]], transformation = FALSE, truncate_max_range = FALSE)
    expect_equal(nrow(flowCore::exprs(written)), 7L)
})

test_that("corrupt FCS errors retain the parser reason", {
    bad_fcs <- tempfile(fileext = ".fcs")
    writeLines("not an FCS file", bad_fcs)
    expect_error(spectreasy:::.spectreasy_read_fcs(bad_fcs), "Reason:")
})

test_that("GUI report discovery excludes unrelated project documents", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) {
        api_path <- system.file("api/gui_api.R", package = "spectreasy")
    }
    testthat::skip_if_not(file.exists(api_path))

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)

    project_dir <- tempfile("spectreasy_report_discovery_")
    expected <- c(
        file.path(project_dir, "reports", "panel_overview.pdf"),
        file.path(project_dir, "spectreasy_outputs", "unmix_controls", "qc_controls_report.html")
    )
    unrelated <- c(
        file.path(project_dir, "vignettes", "getting_started.html"),
        file.path(project_dir, "Rplots.pdf")
    )
    for (path in c(expected, unrelated)) {
        dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
        file.create(path)
    }

    discovered <- api_env$gui_project_report_files(project_dir)
    expect_setequal(normalizePath(discovered), normalizePath(expected))
    expect_false(any(normalizePath(unrelated) %in% normalizePath(discovered)))
})
