load_cockpit_api_env <- function() {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    testthat::skip_if_not(file.exists(api_path))
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
