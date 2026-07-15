test_that("cockpit project file routes create and manage input folders", {
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

    router <- plumber::plumb(api_path)
    context_route <- router$routes$project$context$getFunc()
    file_routes <- router$routes$project$files
    direct_file_routes <- file_routes[names(file_routes) == ""]
    list_route <- direct_file_routes[[which(vapply(direct_file_routes, function(route) identical(route$verbs, "GET"), logical(1)))]]$getFunc()
    upload_route <- direct_file_routes[[which(vapply(direct_file_routes, function(route) identical(route$verbs, "POST"), logical(1)))]]$getFunc()
    delete_route <- direct_file_routes[[which(vapply(direct_file_routes, function(route) identical(route$verbs, "DELETE"), logical(1)))]]$getFunc()
    delete_all_route <- file_routes$all$getFunc()

    request <- new.env(parent = emptyenv())
    request$postBody <- jsonlite::toJSON(list(projectPath = project), auto_unbox = TRUE)
    selected <- context_route(request)
    expect_true(selected$success)
    expect_true(dir.exists(file.path(project, "scc")))
    expect_true(dir.exists(file.path(project, "samples")))

    request$postBody <- jsonlite::toJSON(list(
        kind = "controls",
        filename = "Control 01.fcs",
        content_base64 = jsonlite::base64_enc(charToRaw("FCS-test-payload")),
        overwrite = FALSE
    ), auto_unbox = TRUE)
    uploaded <- upload_route(request)
    expect_true(uploaded$success)
    expect_true(file.exists(file.path(project, "scc", "Control 01.fcs")))

    listed <- list_route("controls")
    expect_true(listed$success)
    expect_identical(listed$files$name, "Control 01.fcs")
    expect_identical(listed$files$kind, "controls")
    expect_gt(listed$files$size, 0)

    duplicate <- upload_route(request)
    expect_false(duplicate$success)
    expect_match(duplicate$error, "already exists")

    removed <- delete_route("controls", "Control 01.fcs")
    expect_true(removed$success)
    expect_false(file.exists(file.path(project, "scc", "Control 01.fcs")))

    file.create(file.path(project, "samples", c("Sample 01.fcs", "Sample 02.FCS")))
    deleted_all <- delete_all_route("samples")
    expect_true(deleted_all$success)
    expect_identical(deleted_all$deleted, 2L)
    expect_length(list.files(file.path(project, "samples"), pattern = "\\.fcs$", ignore.case = TRUE), 0L)
})

test_that("project picker routes distinguish opening from creating", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    skip_if_not(file.exists(api_path))

    router <- plumber::plumb(api_path)
    select_route <- router$routes$project$select$getFunc()
    create_route <- router$routes$project$create$getFunc()
    select_body <- paste(deparse(body(select_route)), collapse = "\n")
    create_body <- paste(deparse(body(create_route)), collapse = "\n")
    expect_match(select_body, "allow_create = FALSE", fixed = TRUE)
    expect_match(create_body, "allow_create = TRUE", fixed = TRUE)
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
