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

test_that("project scans do not initialize or mutate an untouched folder", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    skip_if_not(file.exists(api_path))

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)
    project <- tempfile("cockpit_read_only_scan_")
    dir.create(project)
    on.exit(unlink(project, recursive = TRUE, force = TRUE), add = TRUE)

    api_env$gui_project_scan(project)
    expect_false(file.exists(file.path(project, ".spectreasy", "project.json")))
    expect_false(dir.exists(file.path(project, "scc")))
    expect_false(dir.exists(file.path(project, "samples")))
})

test_that("project data revisions change when applet inputs change", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    skip_if_not(file.exists(api_path))

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)
    project <- tempfile("cockpit_data_revision_")
    dir.create(file.path(project, "scc"), recursive = TRUE)
    dir.create(file.path(project, "samples"), recursive = TRUE)
    on.exit(unlink(project, recursive = TRUE, force = TRUE), add = TRUE)

    control <- file.path(project, "scc", "control.fcs")
    writeBin(charToRaw("FCS3.1-one"), control)
    first <- api_env$gui_project_scan(project)$data_revision
    dir.create(file.path(project, "derived"))
    writeBin(charToRaw("unmixed-output"), file.path(project, "derived", "output.fcs"))
    unchanged <- api_env$gui_project_scan(project)$data_revision
    writeBin(charToRaw("FCS3.1-a-longer-payload"), control)
    second <- api_env$gui_project_scan(project)$data_revision

    expect_type(first, "character")
    expect_lt(nchar(first), 64)
    expect_identical(first, unchanged)
    expect_false(identical(first, second))
})

test_that("project picker routes distinguish opening from creating", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    skip_if_not(file.exists(api_path))

    endpoints <- unname(unlist(plumber::plumb(api_path)$endpoints, recursive = FALSE))
    by_path <- function(path, verb) {
        matches <- which(vapply(endpoints, function(endpoint) {
            identical(endpoint$path, path) && verb %in% endpoint$verbs
        }, logical(1)))
        endpoints[[matches[1]]]
    }
    select_route <- by_path("/project/select", "POST")$getFunc()
    create_route <- by_path("/project/create", "POST")$getFunc()
    select_body <- paste(deparse(body(select_route)), collapse = "\n")
    create_body <- paste(deparse(body(create_route)), collapse = "\n")
    expect_match(select_body, "allow_create = FALSE", fixed = TRUE)
    expect_match(create_body, "allow_create = TRUE", fixed = TRUE)

    project <- tempfile("cockpit_picker_project_")
    dir.create(project)
    on.exit(unlink(project, recursive = TRUE, force = TRUE), add = TRUE)
    select_env <- environment(select_route)
    create_env <- environment(create_route)
    original_select_picker <- get("gui_pick_project_directory", envir = select_env)
    original_create_picker <- get("gui_pick_project_directory", envir = create_env)
    on.exit(assign("gui_pick_project_directory", original_select_picker, envir = select_env), add = TRUE)
    on.exit(assign("gui_pick_project_directory", original_create_picker, envir = create_env), add = TRUE)
    picker_modes <- logical()
    fake_picker <- function(initial_dir, allow_create = FALSE) {
        picker_modes <<- c(picker_modes, allow_create)
        project
    }
    assign("gui_pick_project_directory", fake_picker, envir = select_env)
    assign("gui_pick_project_directory", fake_picker, envir = create_env)

    opened <- select_route(NULL)
    created <- create_route(NULL)

    expect_identical(picker_modes, c(FALSE, TRUE))
    expect_true(opened$success)
    expect_true(created$success)
    expect_identical(as.character(opened$project_path), normalizePath(project))
    expect_identical(as.character(created$project_path), normalizePath(project))
    expect_named(opened$project, "project_path")
    expect_null(opened$project$files)
})

test_that("project context route activates the selected folder", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    skip_if_not(file.exists(api_path))

    endpoints <- unname(unlist(plumber::plumb(api_path)$endpoints, recursive = FALSE))
    matches <- which(vapply(endpoints, function(endpoint) {
        identical(endpoint$path, "/project/context") && "POST" %in% endpoint$verbs
    }, logical(1)))
    context_route <- endpoints[[matches[1]]]$getFunc()
    project <- tempfile("cockpit_context_project_")
    dir.create(project)
    on.exit(unlink(project, recursive = TRUE, force = TRUE), add = TRUE)
    context_options <- c(
        "spectreasy.project_dir", "spectreasy.matrix_dir", "spectreasy.samples_dir",
        "spectreasy.gating_scc_dir", "spectreasy.gating_control_file",
        "spectreasy.gating_gate_file", "spectreasy.project_selected"
    )
    old_options <- options()[context_options]
    on.exit(options(old_options), add = TRUE)
    req <- new.env(parent = emptyenv())
    req$postBody <- jsonlite::toJSON(list(projectPath = project), auto_unbox = TRUE)

    response <- context_route(req)

    expect_true(response$success)
    expect_identical(as.character(response$project$project_path), normalizePath(project))
    expect_identical(getOption("spectreasy.project_dir"), normalizePath(project))
    expect_identical(getOption("spectreasy.matrix_dir"), normalizePath(project))
    expect_true(isTRUE(getOption("spectreasy.project_selected")))
})

test_that("project input layout persists GUI renames and follows manual filesystem renames", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    skip_if_not(file.exists(api_path))

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)
    project <- tempfile("cockpit_project_layout_")
    dir.create(file.path(project, "scc"), recursive = TRUE)
    dir.create(file.path(project, "samples"), recursive = TRUE)
    file.create(file.path(project, "scc", "control.fcs"))
    file.create(file.path(project, "samples", "sample.fcs"))
    on.exit(unlink(project, recursive = TRUE, force = TRUE), add = TRUE)

    layout <- api_env$gui_project_layout(project)
    expect_identical(layout$control_input_dir, "scc")
    expect_identical(layout$sample_input_dir, "samples")
    expect_true(file.exists(file.path(project, "scc", api_env$gui_project_input_marker_name())))
    expect_true(file.exists(file.path(project, "samples", api_env$gui_project_input_marker_name())))

    expect_true(file.rename(file.path(project, "scc"), file.path(project, "single_stains")))
    reconciled <- api_env$gui_project_layout(project)
    expect_identical(reconciled$control_input_dir, "single_stains")
    expect_identical(api_env$gui_project_file_rows("controls", project)$name, "control.fcs")

    updated <- api_env$gui_update_project_input_dir(project, "samples", "data/specimens")
    expect_identical(updated$sample_input_dir, "data/specimens")
    expect_true(file.exists(file.path(project, "data", "specimens", "sample.fcs")))
    expect_identical(api_env$gui_project_file_rows("samples", project)$name, "sample.fcs")

    persisted <- jsonlite::fromJSON(file.path(project, ".spectreasy", "project.json"))
    expect_identical(persisted$control_input_dir, "single_stains")
    expect_identical(persisted$sample_input_dir, "data/specimens")
    scan <- api_env$gui_project_scan(project)
    expect_identical(scan$scan$controls, 1L)
    expect_identical(scan$scan$samples, 1L)
    expect_error(api_env$gui_update_project_input_dir(project, "samples", "../outside"), "invalid path component")
    expect_error(api_env$gui_update_project_input_dir(project, "samples", "single_stains"), "must be different")
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
        "/af_profiles/rename",
        "/project/context",
        "/project/layout",
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
