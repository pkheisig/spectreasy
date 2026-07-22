cockpit_api_path <- function() {
    namespace_path <- tryCatch(
        getNamespaceInfo(asNamespace("spectreasy"), "path"),
        error = function(e) ""
    )
    candidates <- c(
        file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R"),
        system.file("api", "gui_api.R", package = "spectreasy"),
        if (nzchar(namespace_path)) file.path(namespace_path, "api", "gui_api.R") else ""
    )
    candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
    if (!length(candidates)) testthat::skip("The installed cockpit API source is unavailable.")
    candidates[[1]]
}

load_cockpit_api_env <- function() {
    api_path <- cockpit_api_path()
    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)
    api_env
}

test_that("cockpit workflow JSON and numeric settings are validated strictly", {
    api_env <- load_cockpit_api_env()
    request <- new.env(parent = emptyenv())
    request$postBody <- "{broken"
    expect_error(api_env$gui_workflow_body(request), "valid JSON")

    expect_equal(api_env$gui_workflow_number(list(events = 12), "events", 5, integer = TRUE, minimum = 1), 12)
    expect_error(api_env$gui_workflow_number(list(events = 1.5), "events", 5, integer = TRUE), "integer")
    expect_error(api_env$gui_workflow_number(list(fraction = 1.2), "fraction", 0.5, maximum = 1), "at most 1")
    expect_error(api_env$gui_workflow_number(list(fraction = Inf), "fraction", 0.5), "Invalid numeric value")
})

test_that("cockpit workflow project roots and session tokens fail closed", {
    api_env <- load_cockpit_api_env()
    project <- tempfile("cockpit_contract_project_")
    dir.create(project)
    on.exit(unlink(project, recursive = TRUE, force = TRUE), add = TRUE)

    expect_equal(api_env$gui_workflow_root(list(projectPath = project)), normalizePath(project))
    expect_error(api_env$gui_workflow_root(list(projectPath = paste0(project, "-missing"))), "not found")

    withr::local_options(list(spectreasy.gui_api_token = "active-token"))
    expect_true(api_env$gui_api_token_value_allowed("active-token"))
    expect_false(api_env$gui_api_token_value_allowed("stale-token"))
    expect_false(api_env$gui_api_token_value_allowed(""))
})

test_that("cockpit shutdown is scoped to its recorded httpuv server", {
    api_path <- cockpit_api_path()
    source_text <- paste(readLines(api_path, warn = FALSE), collapse = "\n")
    expect_false(grepl("stopAllServers", source_text, fixed = TRUE))
    expect_match(source_text, "httpuv::stopServer(server)", fixed = TRUE)
})

test_that("cockpit gating uses shared alternate scatter channel aliases", {
    api_env <- load_cockpit_api_env()
    channels <- c("FSC-A", "BSSC-A", "VSSC-H", "V1-A")
    expect_identical(api_env$gate_pick_channel(channels, "SSC", "A"), "BSSC-A")
    expect_setequal(api_env$gate_scatter_channels(channels), c("FSC-A", "BSSC-A", "VSSC-H"))
})

test_that("cockpit AI-QC readiness and generation remain project-contained", {
    api_env <- load_cockpit_api_env()
    project <- tempfile("cockpit_ai_qc_project_")
    matrix_dir <- file.path(project, "spectreasy_outputs", "unmix_controls")
    dir.create(matrix_dir, recursive = TRUE)
    on.exit(unlink(project, recursive = TRUE, force = TRUE), add = TRUE)
    utils::write.csv(
        data.frame(fluorophore = c("FITC", "PE"), `B1-A` = c(1, 0.1), `YG1-A` = c(0.2, 1), check.names = FALSE),
        file.path(matrix_dir, "scc_reference_matrix.csv"), row.names = FALSE
    )

    readiness <- api_env$gui_ai_qc_readiness(project)
    expect_identical(readiness$status, "not_generated")
    report_dir <- file.path(project, "spectreasy_outputs", "unmix_controls", "qc_controls")
    dir.create(report_dir, recursive = TRUE)
    prompt_path <- file.path(report_dir, "qc_controls_report_ai_qc_prompt.txt")
    writeLines("SPECTREASY NUMERICAL QC INTERPRETATION PROMPT\nNo raw event-level data are included.", prompt_path)

    generated <- api_env$gui_ai_qc_generate(list(
        projectPath = project, output_root = "spectreasy_outputs",
        scope = "control", privacy = "strict", detail = "compact", reference = "none"
    ))
    expect_identical(generated$status, "ready")
    expect_length(generated$artifact_paths, 1L)
    expect_true(all(vapply(generated$artifact_paths, function(path) !startsWith(path, "/"), logical(1))))
    expect_false(grepl(normalizePath(project), generated$prompt, fixed = TRUE))
    expect_match(tolower(generated$prompt), "no raw event-level data", fixed = TRUE)
})

test_that("cockpit report discovery exposes only a precomputed prompt companion", {
    api_env <- load_cockpit_api_env()
    project <- tempfile("cockpit_report_prompt_")
    report_dir <- file.path(project, "spectreasy_outputs", "unmix_samples", "qc_samples")
    dir.create(report_dir, recursive = TRUE)
    on.exit(unlink(project, recursive = TRUE, force = TRUE), add = TRUE)
    report <- file.path(report_dir, "qc_samples_report.html")
    writeLines("<html></html>", report)
    expect_identical(api_env$gui_project_report_prompt_file(report), "")
    prompt <- file.path(report_dir, "qc_samples_report_ai_qc_prompt.txt")
    writeLines("saved prompt", prompt)
    expect_identical(api_env$gui_project_report_prompt_file(report), normalizePath(prompt))
    reports <- api_env$gui_project_report_files(project, output_root = "spectreasy_outputs", report_type = "sample")
    expect_identical(normalizePath(reports), normalizePath(report))
})
