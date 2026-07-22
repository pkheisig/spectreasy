with_cockpit_af_profiles <- function(names, code) {
    old_dir <- getOption("spectreasy.af_profile_dir")
    profile_dir <- tempfile("cockpit_af_library_")
    options(spectreasy.af_profile_dir = profile_dir)
    on.exit(options(spectreasy.af_profile_dir = old_dir), add = TRUE)
    profile <- matrix(
        c(1, 0.4, 0.1),
        nrow = 1,
        dimnames = list("AF", c("B1-A", "YG1-A", "R1-A"))
    )
    for (name in names) spectreasy::save_af_profile(name, profile)
    force(code)
}

test_that("cockpit derives collision-safe AF profile names from source files", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    skip_if_not(file.exists(api_path))

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)
    expect_identical(api_env$gui_default_af_profile_name("/data/Unstained cells.fcs"), "Unstained_cells")
    expect_identical(
        api_env$gui_default_af_profile_name("/data/Unstained cells.fcs", c("Unstained_cells", "Unstained_cells_2")),
        "Unstained_cells_3"
    )
})

test_that("cockpit active AF profile filters only the obsolete unstained cell control", {
    with_cockpit_af_profiles("PBMC_saved", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    skip_if_not(file.exists(api_path))

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)
    root <- tempfile("cockpit_af_project_")
    dir.create(file.path(root, ".spectreasy"), recursive = TRUE)
    jsonlite::write_json(
        list(profile_name = "PBMC_saved"),
        file.path(root, ".spectreasy", "active_af_profile.json"),
        auto_unbox = TRUE
    )
    mapping <- data.frame(
        filename = c("unstained.fcs", "dead_unstained.fcs", "bead_unstained.fcs", "FITC.fcs"),
        fluorophore = c("AF", "AF_dead", "AF_beads", "FITC"),
        marker = c("Autofluorescence", "Dead cell background", "Bead background", "CD3"),
        channel = c("V1-A", "V2-A", "B1-A", "B2-A"),
        control.type = c("cells", "cells", "beads", "cells"),
        universal.negative = c("", "", "", "unstained.fcs"),
        stringsAsFactors = FALSE
    )

    annotated <- api_env$gui_annotate_active_af_mapping(mapping, root)
    expect_identical(annotated$ignored, c(TRUE, FALSE, FALSE, FALSE))
    expect_match(annotated$ignored_reason[1], "PBMC_saved", fixed = TRUE)

    filtered <- api_env$gui_filter_active_af_mapping(mapping, root)
    expect_equal(filtered$filename, c("dead_unstained.fcs", "bead_unstained.fcs", "FITC.fcs"))
    expect_identical(filtered$universal.negative[filtered$filename == "FITC.fcs"], "")
    })
})

test_that("cockpit active AF profile selection is project specific", {
    with_cockpit_af_profiles("dataset_one", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    skip_if_not(file.exists(api_path))

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)
    first <- tempfile("cockpit_af_first_")
    second <- tempfile("cockpit_af_second_")
    dir.create(file.path(first, ".spectreasy"), recursive = TRUE)
    dir.create(second)
    jsonlite::write_json(
        list(profile_name = "dataset_one"),
        file.path(first, ".spectreasy", "active_af_profile.json"),
        auto_unbox = TRUE
    )

    expect_identical(api_env$gui_read_active_af_profile(first), "dataset_one")
    expect_identical(api_env$gui_read_active_af_profile(second), "")
    })
})

test_that("cockpit can unlink only the AF profile linked to a dataset", {
    with_cockpit_af_profiles("Erlangen_PBMC", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    skip_if_not(file.exists(api_path))

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)
    root <- tempfile("cockpit_af_unlink_")
    dir.create(file.path(root, ".spectreasy"), recursive = TRUE)
    jsonlite::write_json(
        list(profile_name = "Erlangen_PBMC"),
        file.path(root, ".spectreasy", "active_af_profile.json"),
        auto_unbox = TRUE
    )

    expect_error(
        api_env$gui_unlink_active_af_profile("another_profile", root),
        "not linked"
    )
    expect_identical(api_env$gui_read_active_af_profile(root), "Erlangen_PBMC")
    expect_identical(
        api_env$gui_unlink_active_af_profile("Erlangen_PBMC", root),
        "Erlangen_PBMC"
    )
    expect_identical(api_env$gui_read_active_af_profile(root), "")
    expect_identical(
        api_env$gui_unlink_active_af_profile("Erlangen_PBMC", root),
        "Erlangen_PBMC"
    )
    })
})

test_that("cockpit rename keeps the active dataset profile link synchronized", {
    with_cockpit_af_profiles("profile_before", {
        api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
        if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
        skip_if_not(file.exists(api_path))

        router <- plumber::plumb(api_path)
        endpoints <- unlist(router$endpoints, recursive = FALSE, use.names = FALSE)
        endpoint <- endpoints[[which(vapply(endpoints, function(item) identical(item$path, "/af_profiles/rename"), logical(1)))[1]]]
        rename_route <- endpoint$getFunc()
        root <- tempfile("cockpit_af_rename_")
        dir.create(file.path(root, ".spectreasy"), recursive = TRUE)
        on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
        jsonlite::write_json(
            list(profile_name = "profile_before"),
            file.path(root, ".spectreasy", "active_af_profile.json"),
            auto_unbox = TRUE
        )
        request <- list(postBody = jsonlite::toJSON(list(
            profile_name = "profile_before",
            new_name = "profile_after",
            projectPath = root
        ), auto_unbox = TRUE))

        result <- rename_route(request)
        expect_true(result$success)
        expect_identical(result$new_name, "profile_after")
        active <- jsonlite::fromJSON(file.path(root, ".spectreasy", "active_af_profile.json"))
        expect_identical(active$profile_name, "profile_after")
        expect_true("profile_after" %in% spectreasy::list_af_profiles()$name)
        expect_false("profile_before" %in% spectreasy::list_af_profiles()$name)
    })
})

test_that("cockpit removes a stale project link when its global AF profile is gone", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    skip_if_not(file.exists(api_path))

    old_dir <- getOption("spectreasy.af_profile_dir")
    options(spectreasy.af_profile_dir = tempfile("empty_af_library_"))
    on.exit(options(spectreasy.af_profile_dir = old_dir), add = TRUE)
    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)
    root <- tempfile("cockpit_af_stale_")
    dir.create(file.path(root, ".spectreasy"), recursive = TRUE)
    config <- file.path(root, ".spectreasy", "active_af_profile.json")
    jsonlite::write_json(list(profile_name = "deleted_profile"), config, auto_unbox = TRUE)

    expect_identical(api_env$gui_read_active_af_profile(root), "")
    expect_false(file.exists(config))
})

test_that("cockpit workflows resolve project paths without changing the R working directory", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    skip_if_not(file.exists(api_path))

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)
    root <- tempfile("cockpit_workflow_project_")
    dir.create(root)
    launch_dir <- tempfile("cockpit_launch_dir_")
    dir.create(launch_dir)
    old_wd <- setwd(launch_dir)
    on.exit(setwd(old_wd), add = TRUE)
    old_project_option <- getOption("spectreasy.project_dir", NULL)

    expect_identical(
        api_env$gui_workflow_resolve_path(file.path("spectreasy_outputs", "unmix_samples"), root),
        file.path(normalizePath(root, mustWork = TRUE), "spectreasy_outputs", "unmix_samples")
    )
    absolute <- normalizePath(file.path(root, "absolute.csv"), mustWork = FALSE)
    expect_identical(
        api_env$gui_workflow_resolve_path(absolute, root),
        file.path(normalizePath(root, mustWork = TRUE), "absolute.csv")
    )
    outside <- normalizePath(file.path(launch_dir, "outside.csv"), mustWork = FALSE)
    expect_error(api_env$gui_workflow_resolve_path(outside, root), "outside the active project")
    expect_error(api_env$gui_workflow_resolve_path("../outside.csv", root), "outside the active project")
    expect_identical(api_env$gui_workflow_file_or_null("", root), NULL)
    run <- suppressWarnings(api_env$gui_workflow_run(
        list(projectPath = root),
        "test",
        {
            cat("workflow output\n")
            warning("workflow warning")
            getwd()
        }
    ))

    expect_true(run$success)
    expect_identical(run$result, normalizePath(launch_dir, mustWork = TRUE))
    expect_identical(getwd(), normalizePath(launch_dir, mustWork = TRUE))
    expect_identical(run$project_path, normalizePath(root, mustWork = TRUE))
    expect_identical(getOption("spectreasy.project_dir", NULL), old_project_option)
    expect_true(any(grepl("workflow output", run$logs, fixed = TRUE)))
    expect_true(any(grepl("Warning: workflow warning", run$logs, fixed = TRUE)))
})

test_that("cockpit resolves only in-project HTML reports", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    skip_if_not(file.exists(api_path))

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)
    root <- tempfile("cockpit_report_project_")
    dir.create(file.path(root, "spectreasy_outputs", "unmix_controls"), recursive = TRUE)
    report <- file.path(root, "spectreasy_outputs", "unmix_controls", "qc_controls_report.html")
    writeLines("<!doctype html><title>QC</title><h1>Report</h1>", report)

    expect_identical(
        api_env$gui_resolve_project_report("spectreasy_outputs/unmix_controls/qc_controls_report.html", root),
        normalizePath(report)
    )
    expect_error(api_env$gui_resolve_project_report("../outside.html", root), "Invalid project report path")
    expect_error(api_env$gui_resolve_project_report("spectreasy_outputs/unmix_controls/report.pdf", root), "not found")
})

test_that("cockpit discovers and classifies control and sample reports from relative project paths", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    skip_if_not(file.exists(api_path))

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)
    root <- tempfile("sample_named_cockpit_project_")
    control <- file.path(root, "spectreasy_outputs", "unmix_controls", "qc_controls", "qc_controls_report.html")
    control_2 <- file.path(root, "spectreasy_outputs", "unmix_controls", "qc_controls_2", "qc_controls_report.pdf")
    control_stage_2 <- file.path(root, "spectreasy_outputs", "unmix_controls_2", "qc_controls", "qc_controls_report.html")
    sample <- file.path(root, "spectreasy_outputs", "unmix_samples", "qc_samples", "qc_samples_report.html")
    sample_2 <- file.path(root, "spectreasy_outputs", "unmix_samples", "qc_samples_2", "qc_samples_report.pdf")
    custom <- file.path(root, "spectreasy_outputs_065", "unmix_samples", "qc_samples", "qc_samples_report.html")
    dir.create(dirname(control), recursive = TRUE)
    dir.create(dirname(control_2), recursive = TRUE)
    dir.create(dirname(control_stage_2), recursive = TRUE)
    dir.create(dirname(sample), recursive = TRUE)
    dir.create(dirname(sample_2), recursive = TRUE)
    dir.create(dirname(custom), recursive = TRUE)
    writeLines("<!doctype html><title>Controls</title>", control)
    writeLines("PDF", control_2)
    writeLines("<!doctype html><title>Controls rerun</title>", control_stage_2)
    writeLines("<!doctype html><title>Samples</title>", sample)
    writeLines("PDF", sample_2)
    writeLines("<!doctype html><title>Custom output root</title>", custom)

    reports <- api_env$gui_project_report_files(root)
    relative <- vapply(reports, api_env$gui_project_relative_path, character(1), root = root)
    types <- vapply(relative, api_env$gui_project_report_type, character(1))

    expect_setequal(relative, c(
        "spectreasy_outputs/unmix_controls/qc_controls/qc_controls_report.html",
        "spectreasy_outputs/unmix_controls/qc_controls_2/qc_controls_report.pdf",
        "spectreasy_outputs/unmix_controls_2/qc_controls/qc_controls_report.html",
        "spectreasy_outputs/unmix_samples/qc_samples/qc_samples_report.html",
        "spectreasy_outputs/unmix_samples/qc_samples_2/qc_samples_report.pdf",
        "spectreasy_outputs_065/unmix_samples/qc_samples/qc_samples_report.html"
    ))
    expect_identical(unname(types[grepl("unmix_controls", relative)]), c("Control QC", "Control QC", "Control QC"))
    expect_identical(unname(types[grepl("unmix_samples", relative)]), c("Sample QC", "Sample QC", "Sample QC"))
    expect_identical(api_env$gui_project_report_type("spectreasy_outputs/unmix_samples_10/qc_samples/qc_samples_report.html"), "Sample QC")

    exact_control <- api_env$gui_project_report_files(
        root,
        output_root = "spectreasy_outputs",
        report_type = "control"
    )
    exact_sample <- api_env$gui_project_report_files(
        root,
        output_root = "spectreasy_outputs",
        report_type = "sample"
    )
    expect_setequal(normalizePath(exact_control), normalizePath(c(control, control_2, control_stage_2)))
    expect_setequal(normalizePath(exact_sample), normalizePath(c(sample, sample_2)))
    expect_false(normalizePath(custom) %in% normalizePath(exact_sample))
})

test_that("cockpit exports existing HTML through Chromium without rerunning QC", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    skip_if_not(file.exists(api_path))

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)
    skip_if(!nzchar(api_env$gui_find_chromium()), "Chrome or Chromium is not available")
    html <- tempfile(fileext = ".html")
    pdf <- tempfile(fileext = ".pdf")
    writeLines("<!doctype html><style>body{background:#eef7f3}</style><h1>Existing QC report</h1>", html)

    output <- tryCatch(
        api_env$gui_export_html_report_pdf(html, pdf),
        error = function(e) e
    )
    if (inherits(output, "error")) {
        skip(paste("Chromium could not start in this test environment:", conditionMessage(output)))
    }
    expect_true(file.exists(output))
    expect_gt(file.info(output)$size, 1000)
    expect_identical(rawToChar(readBin(output, "raw", n = 4)), "%PDF")
})
