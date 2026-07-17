test_that("cockpit project initialization is explicit and file helpers are safe", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    skip_if_not(file.exists(api_path))

    project <- tempfile("cockpit_project_files_")
    dir.create(project)
    on.exit(unlink(project, recursive = TRUE, force = TRUE), add = TRUE)
    old_options <- options()[c(
        "spectreasy.project_dir",
        "spectreasy.matrix_dir",
        "spectreasy.samples_dir",
        "spectreasy.gating_scc_dir",
        "spectreasy.project_selected"
    )]
    on.exit(options(old_options), add = TRUE)

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)

    api_env$gui_set_project_context(project)
    scan <- api_env$gui_project_scan(project)
    expect_setequal(scan$missing_input_dirs, c("scc", "samples"))
    expect_false(dir.exists(file.path(project, "scc")))
    expect_false(dir.exists(file.path(project, "samples")))

    api_env$gui_ensure_project_input_dirs(project)
    expect_true(dir.exists(file.path(project, "scc")))
    expect_true(dir.exists(file.path(project, "samples")))

    files <- file.path(project, "scc", c("Control 10.fcs", "Control 2.fcs", "Control 1.fcs"))
    file.create(files)
    listed <- api_env$gui_project_file_rows("controls")
    expect_identical(listed$name, c("Control 1.fcs", "Control 2.fcs", "Control 10.fcs"))

    valid_fcs <- tempfile(fileext = ".fcs")
    invalid_fcs <- tempfile(fileext = ".fcs")
    writeBin(c(charToRaw("FCS3.1"), charToRaw("payload")), valid_fcs)
    writeBin(charToRaw("not-fcs"), invalid_fcs)
    expect_true(api_env$gui_validate_fcs_upload(valid_fcs))
    expect_error(api_env$gui_validate_fcs_upload(invalid_fcs), "valid FCS header")
})

test_that("project picker routes distinguish opening from creating", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    skip_if_not(file.exists(api_path))

    collect_endpoints <- function(node) {
        if (inherits(node, "PlumberEndpoint")) return(list(node))
        if (!is.list(node)) return(list())
        unlist(lapply(seq_along(node), function(index) collect_endpoints(node[[index]])), recursive = FALSE)
    }
    endpoints <- collect_endpoints(plumber::plumb(api_path)$routes)
    by_path <- function(path) endpoints[[which(vapply(endpoints, function(endpoint) identical(endpoint$path, path), logical(1)))[1]]]
    select_route <- by_path("/project/select")$getFunc()
    create_route <- by_path("/project/create")$getFunc()
    select_body <- paste(deparse(body(select_route)), collapse = "\n")
    create_body <- paste(deparse(body(create_route)), collapse = "\n")
    expect_match(select_body, "allow_create = FALSE", fixed = TRUE)
    expect_match(create_body, "allow_create = TRUE", fixed = TRUE)
})

test_that("cockpit API keeps all critical route families registered", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    skip_if_not(file.exists(api_path))

    router <- plumber::plumb(api_path)
    endpoints <- unlist(router$endpoints, recursive = FALSE, use.names = FALSE)
    paths <- unique(vapply(endpoints, function(endpoint) endpoint$path, character(1)))
    critical_paths <- c(
        "/status",
        "/gui_state",
        "/control_mapping",
        "/gate_files",
        "/spectral_panel",
        "/spectral_panel_metrics",
        "/af_profiles",
        "/af_profiles/delete",
        "/project/context",
        "/project/files",
        "/project/initialize",
        "/project/upload-start",
        "/workflow/control",
        "/workflow/sample"
    )

    expect_true(all(critical_paths %in% paths))
    expect_false("/import_sample_content" %in% paths)
    expect_false("/workflow/synthetic" %in% paths)
})

test_that("cockpit workflows use explicit project paths and never mutate the working directory", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    skip_if_not(file.exists(api_path))

    source_text <- paste(readLines(api_path, warn = FALSE), collapse = "\n")
    expect_false(grepl("setwd(", source_text, fixed = TRUE))
    expect_match(source_text, "gui_workflow_resolve_path(gui_workflow_value(body, \"scc_dir\"", fixed = TRUE)
    expect_match(source_text, "gui_workflow_resolve_path(gui_workflow_value(body, \"sample_dir\"", fixed = TRUE)
    expect_match(source_text, "gui_workflow_resolve_path(gui_workflow_value(body, \"fcs_file\"", fixed = TRUE)
    expect_match(source_text, "gui_workflow_resolve_path(gui_workflow_value(body, \"output_dir\"", fixed = TRUE)
})

test_that("cockpit requests keep project identity scoped to the request", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    skip_if_not(file.exists(api_path))

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)
    project_a <- tempfile("cockpit_project_a_")
    project_b <- tempfile("cockpit_project_b_")
    dir.create(file.path(project_a, "scc"), recursive = TRUE)
    dir.create(file.path(project_b, "scc"), recursive = TRUE)
    file.create(file.path(project_a, "scc", "A.fcs"))
    file.create(file.path(project_b, "scc", "B.fcs"))
    on.exit(unlink(c(project_a, project_b), recursive = TRUE, force = TRUE), add = TRUE)

    old_options <- options(spectreasy.project_dir = project_a, spectreasy.matrix_dir = project_a)
    on.exit(options(old_options), add = TRUE)
    expect_identical(api_env$gui_project_file_rows("controls", project_a)$name, "A.fcs")
    expect_identical(api_env$gui_project_file_rows("controls", project_b)$name, "B.fcs")

    observed <- api_env$gui_with_project_context(project_b, getOption("spectreasy.project_dir"))
    expect_identical(observed, normalizePath(project_b, mustWork = TRUE))
    expect_identical(getOption("spectreasy.project_dir"), project_a)
    api_env$gui_with_project_context(project_a, options(spectreasy.gating_state_cache = list(project = "A")))
    api_env$gui_with_project_context(project_b, options(spectreasy.gating_state_cache = list(project = "B")))
    expect_identical(api_env$gui_with_project_context(project_a, getOption("spectreasy.gating_state_cache"))$project, "A")
    expect_identical(api_env$gui_with_project_context(project_b, getOption("spectreasy.gating_state_cache"))$project, "B")
})

test_that("cockpit project files reject traversal and non-FCS input", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    skip_if_not(file.exists(api_path))

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)
    project <- tempfile("cockpit_project_file_safety_")
    dir.create(project)
    on.exit(unlink(project, recursive = TRUE, force = TRUE), add = TRUE)
    old_options <- options(spectreasy.matrix_dir = project, spectreasy.project_dir = project)
    on.exit(options(old_options), add = TRUE)

    expect_error(api_env$gui_project_file_location("controls", "../escape.fcs"), "Invalid")
    expect_error(api_env$gui_project_file_location("samples", "notes.txt"), "Only FCS")
    expect_error(api_env$gui_project_file_location("other", "test.fcs"), "controls.*samples")
})
