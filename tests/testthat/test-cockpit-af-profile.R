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

test_that("cockpit terminal accepts only the configured cockpit origin and session token", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    skip_if_not(file.exists(api_path))

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)
    old_origins <- getOption("spectreasy.gui_allowed_origins")
    old_token <- getOption("spectreasy.gui_api_token")
    options(
        spectreasy.gui_allowed_origins = "http://127.0.0.1:5173",
        spectreasy.gui_api_token = "test-session-token"
    )
    on.exit(options(
        spectreasy.gui_allowed_origins = old_origins,
        spectreasy.gui_api_token = old_token
    ), add = TRUE)
    request <- new.env(parent = emptyenv())
    request$HTTP_ORIGIN <- "http://127.0.0.1:5173"
    request$HTTP_X_SPECTREASY_TOKEN <- "test-session-token"
    expect_true(api_env$gui_terminal_origin_allowed(request))
    expect_true(api_env$gui_api_token_allowed(request))
    request$HTTP_ORIGIN <- "https://pkheisig.github.io"
    expect_false(api_env$gui_terminal_origin_allowed(request))
    request$HTTP_ORIGIN <- "https://example.com"
    expect_false(api_env$gui_terminal_origin_allowed(request))
    request$HTTP_ORIGIN <- "http://127.0.0.1:5173"
    request$HTTP_X_SPECTREASY_TOKEN <- "wrong-token"
    expect_false(api_env$gui_api_token_allowed(request))
})

test_that("cockpit terminal evaluates persistent R code and intercepts quit", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    skip_if_not(file.exists(api_path))

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)
    root <- tempfile("cockpit_terminal_project_")
    dir.create(root)
    starting_cwd <- getwd()

    assigned <- api_env$gui_terminal_evaluate("terminal_value <- 41", root)
    evaluated <- api_env$gui_terminal_evaluate("terminal_value + 1", root)
    stopped <- api_env$gui_terminal_evaluate("q()", root)

    expect_true(assigned$success)
    expect_true(evaluated$success)
    expect_match(evaluated$output, "42", fixed = TRUE)
    expect_true(evaluated$refresh)
    expect_true(stopped$shutdown_requested)
    expect_false(stopped$refresh)
    expect_identical(getwd(), starting_cwd)
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
