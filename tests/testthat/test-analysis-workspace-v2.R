analysis_v2_api_env <- function() {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api", "gui_api.R", package = "spectreasy")
    if (!nzchar(api_path) || !file.exists(api_path)) testthat::skip("The analysis API source is unavailable.")
    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)
    api_env
}

analysis_v2_fixture <- function() {
    project <- tempfile("analysis_v2_project_")
    sample_dir <- file.path(project, "samples")
    dir.create(sample_dir, recursive = TRUE)
    values <- c(seq(-1, 1, length.out = 60), seq(9, 11, length.out = 60))
    data <- cbind(
        `FSC-A` = seq_len(120),
        `SSC-A` = rep(seq_len(20), 6),
        `CD3-A` = values
    )
    frame <- flowCore::flowFrame(data)
    descriptions <- flowCore::parameters(frame)@data
    descriptions$desc <- c("FSC", "SSC", "CD3")
    flowCore::parameters(frame)@data <- descriptions
    path <- file.path(sample_dir, "sample.fcs")
    flowCore::write.FCS(frame, path)
    list(project = project, relative_file = "samples/sample.fcs", data = data)
}

analysis_v2_lineage_fixture <- function() {
    project <- tempfile("analysis_v2_lineage_")
    sample_dir <- file.path(project, "samples")
    dir.create(sample_dir, recursive = TRUE)
    set.seed(918)
    events <- 180L
    truth <- seq(0, 1, length.out = events)
    noise <- function(sd = 0.04) stats::rnorm(events, 0, sd)
    data <- cbind(
        `FSC-A` = 80 + 35 * truth + noise(0.8),
        `SSC-A` = 45 + 22 * truth + noise(0.7),
        `M1-A` = 9 * truth + noise(),
        `M2-A` = 4 * sin(pi * truth) + noise(),
        `M3-A` = 6 * truth^2 + noise()
    )
    frame <- flowCore::flowFrame(data)
    parameters <- flowCore::parameters(frame)@data
    parameters$desc <- c("FSC", "SSC", "M1", "M2", "M3")
    flowCore::parameters(frame)@data <- parameters
    path <- file.path(sample_dir, "lineage.fcs")
    flowCore::write.FCS(frame, path)
    list(
        project = project,
        relative_file = "samples/lineage.fcs",
        truth = truth,
        markers = colnames(data)
    )
}

analysis_v2_cluster_fixture <- function() {
    project <- tempfile("analysis_v2_clusters_")
    sample_dir <- file.path(project, "samples")
    dir.create(sample_dir, recursive = TRUE)
    set.seed(921)
    group <- rep(seq_len(3L), each = 60L)
    centers <- rbind(
        c(-5, -4, 1, 0, 2),
        c(0, 5, 5, -3, 0),
        c(6, -1, -2, 5, 4)
    )
    data <- centers[group, , drop = FALSE] + matrix(stats::rnorm(180L * 5L, sd = 0.35), 180L, 5L)
    colnames(data) <- c("FSC-A", "SSC-A", "M1-A", "M2-A", "M3-A")
    frame <- flowCore::flowFrame(data)
    parameters <- flowCore::parameters(frame)@data
    parameters$desc <- c("FSC", "SSC", "M1", "M2", "M3")
    flowCore::parameters(frame)@data <- parameters
    path <- file.path(sample_dir, "clusters.fcs")
    flowCore::write.FCS(frame, path)
    list(
        project = project,
        relative_file = "samples/clusters.fcs",
        truth = group,
        markers = colnames(data)
    )
}

analysis_v2_adjusted_rand <- function(truth, predicted) {
    table <- table(truth, predicted)
    choose2 <- function(value) value * (value - 1) / 2
    sum_cells <- sum(choose2(table))
    sum_rows <- sum(choose2(rowSums(table)))
    sum_columns <- sum(choose2(colSums(table)))
    pairs <- choose2(sum(table))
    expected <- sum_rows * sum_columns / pairs
    maximum <- (sum_rows + sum_columns) / 2
    if (maximum == expected) return(1)
    (sum_cells - expected) / (maximum - expected)
}

test_that("analysis v2 samples deterministically without replacement", {
    api <- analysis_v2_api_env()
    indices <- seq_len(1000)
    first <- api$gui_analysis_sample_indices(indices, 75, 42, "sample/population")
    second <- api$gui_analysis_sample_indices(indices, 75, 42, "sample/population")
    changed <- api$gui_analysis_sample_indices(indices, 75, 43, "sample/population")

    expect_identical(first, second)
    expect_length(first, 75L)
    expect_length(unique(first), 75L)
    expect_false(identical(first, changed))
})

test_that("analysis v2 event payload includes marker labels", {
    api <- analysis_v2_api_env()
    fixture <- analysis_v2_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    payload <- api$gui_analysis_event_payload(
        fixture$project, fixture$relative_file,
        x = "FSC-A", y = "SSC-A", max_points = 25L
    )

    expect_equal(nrow(payload$events), 25L)
    expect_identical(payload$channels$channel, c("FSC-A", "SSC-A", "CD3-A"))
    expect_identical(payload$channels$marker, c("FSC", "SSC", "CD3"))
})

test_that("population gate CSV round-trips every geometry and repairs workspace references", {
    api <- analysis_v2_api_env()
    fixture <- analysis_v2_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    template <- api$gui_analysis_default_workspace()
    template$source_path <- "samples"
    template$selected_file <- fixture$relative_file
    template$populations <- c(template$populations, list(
        list(
            id = "rect", name = "Lymphocytes", parent_id = "root",
            type = "rectangle", role = "positive", source_file = NULL,
            x = "FSC-A", y = "SSC-A",
            geometry = list(x_min = 5, x_max = 90, y_min = 1, y_max = 20)
        ),
        list(
            id = "ellipse", name = "Focused", parent_id = "rect",
            type = "ellipse", role = NULL, source_file = fixture$relative_file,
            x = "FSC-A", y = "SSC-A",
            geometry = list(center_x = 40, center_y = 10, radius_x = 12, radius_y = 5)
        ),
        list(
            id = "polygon", name = "Polygon", parent_id = "root",
            type = "polygon", role = "terminal", source_file = NULL,
            x = "FSC-A", y = "SSC-A",
            geometry = list(points = list(
                list(x = 10, y = 2), list(x = 70, y = 3), list(x = 50, y = 17)
            ))
        ),
        list(
            id = "range", name = "CD3 positive", parent_id = "rect",
            type = "range", role = "negative", source_file = NULL,
            x = "CD3-A", y = NULL, geometry = list(min = 0, max = 12)
        )
    ))
    template <- api$gui_analysis_normalize_workspace(template)
    csv <- api$gui_analysis_gate_set_csv(template)
    rows <- api$gui_analysis_gate_set_read_text(csv)
    expect_identical(rows$record_type[[1]], "metadata")
    expect_equal(sum(rows$population_id == "rect"), 2L)
    expect_equal(sum(rows$population_id == "ellipse"), 2L)
    expect_equal(sum(rows$population_id == "polygon"), 3L)
    expect_equal(sum(rows$population_id == "range"), 2L)
    expect_false(any(rows$population_id == "root"))

    current <- api$gui_analysis_default_workspace()
    current$source_path <- "samples"
    current$selected_file <- fixture$relative_file
    current$populations <- c(current$populations, list(list(
        id = "old", name = "Old gate", parent_id = "root",
        type = "range", role = NULL, source_file = NULL,
        x = "FSC-A", y = NULL, geometry = list(min = 1, max = 50)
    )))
    current$active_population_id <- "old"
    current$plots[[1]]$population_id <- "old"
    current$plots[[1]]$overlay_population_id <- "root"
    current$annotations <- list(
        list(population_id = "old", label = "discard"),
        list(population_id = "root", label = "retain")
    )
    current$root_event_id <- 1L
    current$root_population_id <- "old"
    current$root_source_file <- fixture$relative_file
    current <- api$gui_analysis_normalize_workspace(current)

    prepared <- api$gui_analysis_prepare_gate_set_import(fixture$project, csv, current, "template.csv")
    imported <- prepared$workspace
    expect_identical(vapply(imported$populations, `[[`, character(1), "id"), c("root", "rect", "ellipse", "polygon", "range"))
    expect_identical(imported$active_population_id, "root")
    expect_identical(imported$plots[[1]]$population_id, "root")
    expect_null(imported$plots[[1]]$overlay_population_id)
    expect_identical(imported$annotations[[1]]$population_id, "root")
    expect_null(imported$root_event_id)
    expect_equal(prepared$summary$current_gate_count, 1L)
    expect_equal(prepared$summary$imported_gate_count, 4L)
    expect_true(prepared$summary$trajectory_root_cleared)
    expect_gte(length(prepared$warnings), 3L)
    expect_silent(api$gui_analysis_normalize_workspace(imported))

    ellipse <- imported$populations[[3]]
    expect_equal(ellipse$geometry$center_x, 40)
    expect_equal(ellipse$geometry$radius_y, 5)
    polygon <- imported$populations[[4]]
    expect_equal(length(polygon$geometry$points), 3L)
})

test_that("population gate CSV rejects incompatible schemas, hierarchies, channels, files, and coordinates", {
    api <- analysis_v2_api_env()
    fixture <- analysis_v2_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    workspace <- api$gui_analysis_default_workspace()
    workspace$source_path <- "samples"
    workspace$selected_file <- fixture$relative_file
    workspace$populations <- c(workspace$populations, list(
        list(
            id = "a", name = "A", parent_id = "root", type = "rectangle",
            role = NULL, source_file = NULL, x = "FSC-A", y = "SSC-A",
            geometry = list(x_min = 1, x_max = 50, y_min = 1, y_max = 10)
        ),
        list(
            id = "b", name = "B", parent_id = "a", type = "range",
            role = NULL, source_file = NULL, x = "CD3-A", y = NULL,
            geometry = list(min = -1, max = 1)
        )
    ))
    workspace <- api$gui_analysis_normalize_workspace(workspace)
    base_rows <- api$gui_analysis_gate_set_rows(workspace)
    as_csv <- function(rows) paste(capture.output(utils::write.csv(rows, row.names = FALSE, quote = TRUE, na = "")), collapse = "\n")
    prepare <- function(rows) api$gui_analysis_prepare_gate_set_import(fixture$project, as_csv(rows), workspace)

    wrong_version <- base_rows
    wrong_version$schema_version <- "99"
    expect_error(prepare(wrong_version), "Unsupported population gate CSV schema version")
    fractional_version <- base_rows
    fractional_version$schema_version <- "1.5"
    expect_error(prepare(fractional_version), "Unsupported population gate CSV schema version")

    missing_parent <- base_rows
    missing_parent$parent_id[missing_parent$population_id == "b"] <- "missing"
    expect_error(prepare(missing_parent), "existing parent")

    cycle <- base_rows
    cycle$parent_id[cycle$population_id == "a"] <- "b"
    expect_error(prepare(cycle), "cycle")

    duplicate_vertex <- base_rows
    duplicate_vertex$vertex_index[duplicate_vertex$population_id == "a"] <- "1"
    expect_error(prepare(duplicate_vertex), "unique positive integer")

    non_finite <- base_rows
    non_finite$x[non_finite$population_id == "b" & non_finite$vertex_index == "2"] <- "Inf"
    expect_error(prepare(non_finite), "non-finite X")

    missing_channel <- base_rows
    missing_channel$x_channel[missing_channel$population_id == "b"] <- "NOT-A-CHANNEL"
    expect_error(prepare(missing_channel), "missing from the selected source")

    traversal <- base_rows
    traversal$source_file[traversal$population_id == "a"] <- "../outside.fcs"
    expect_error(prepare(traversal), "project-relative path")

    missing_file <- base_rows
    missing_file$source_file[missing_file$population_id == "a"] <- "samples/missing.fcs"
    expect_error(prepare(missing_file), "unavailable source file")

    layout <- api$gui_project_layout(fixture$project, persist = FALSE)
    outside_source <- file.path(fixture$project, layout$control_input_dir, "reference.fcs")
    dir.create(dirname(outside_source), recursive = TRUE, showWarnings = FALSE)
    file.copy(file.path(fixture$project, fixture$relative_file), outside_source)
    file_specific <- base_rows
    file_specific$source_file[file_specific$population_id == "a"] <- api$gui_analysis_relative_path(outside_source, fixture$project)
    prepared <- prepare(file_specific)
    expect_match(paste(prepared$warnings, collapse = " "), "outside the selected source")
})

test_that("analysis v2 pools files while retaining source event identity", {
    api <- analysis_v2_api_env()
    fixture <- analysis_v2_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    second_file <- file.path(fixture$project, "samples", "sample_2.fcs")
    second <- flowCore::flowFrame(fixture$data + matrix(c(0, 0, 2), nrow(fixture$data), 3L, byrow = TRUE))
    parameters <- flowCore::parameters(second)@data
    parameters$desc <- c("FSC", "SSC", "CD3")
    flowCore::parameters(second)@data <- parameters
    flowCore::write.FCS(second, second_file)

    result <- api$gui_analysis_run_method(fixture$project, list(
        files = c(fixture$relative_file, "samples/sample_2.fcs"),
        populationId = "root",
        reductionMethod = "pca",
        clusterMethod = "none",
        markers = c("FSC-A", "SSC-A", "CD3-A"),
        maxEvents = 160L,
        seed = 17L
    ))

    expect_identical(result$metadata$source_files, c(fixture$relative_file, "samples/sample_2.fcs"))
    expect_equal(nrow(result$events), 160L)
    expect_setequal(unique(result$events$source_file), c(fixture$relative_file, "samples/sample_2.fcs"))
    expect_true(all(result$events$source_event_id >= 1L & result$events$source_event_id <= 120L))
    expect_identical(result$events$sample_id, match(result$events$source_file, result$metadata$source_files))

    parameters$desc <- c("FSC", "SSC", "CD19")
    flowCore::parameters(second)@data <- parameters
    flowCore::write.FCS(second, second_file)
    expect_error(
        api$gui_analysis_run_method(fixture$project, list(
            files = c(fixture$relative_file, "samples/sample_2.fcs"),
            populationId = "root", reductionMethod = "pca", clusterMethod = "none",
            markers = c("FSC-A", "SSC-A", "CD3-A"), maxEvents = 160L, seed = 17L
        )),
        "inconsistent marker labels"
    )
})

test_that("analysis v2 discovers cluster markers without inventing identities", {
    api <- analysis_v2_api_env()
    fixture <- analysis_v2_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    directory <- file.path(fixture$project, "spectreasy_outputs", "analysis", "marker-test")
    dir.create(directory, recursive = TRUE)
    events <- data.frame(
        event_id = seq_len(80),
        marker_1 = c(seq(-2, -1, length.out = 40), seq(2, 3, length.out = 40)),
        marker_2 = c(seq(3, 2, length.out = 40), seq(-1, -2, length.out = 40)),
        cluster_id = rep(c(1L, 2L), each = 40L)
    )
    data.table::fwrite(events, file.path(directory, "events.csv"))
    jsonlite::write_json(list(
        analysis_id = "marker-test",
        marker_columns = list(
            list(marker = "CD3", channel = "CD3-A", column = "marker_1"),
            list(marker = "CD19", channel = "CD19-A", column = "marker_2")
        )
    ), file.path(directory, "metadata.json"), auto_unbox = TRUE)

    result <- api$gui_analysis_find_cluster_markers(fixture$project, list(
        analysisId = "marker-test", topN = 1L, minimumAuc = 0.55
    ))

    expect_equal(nrow(result$markers), 2L)
    expect_identical(result$markers$marker, c("CD19", "CD3"))
    expect_true(all(result$markers$auc > 0.99))
    expect_true(file.exists(file.path(fixture$project, result$metadata$full_table$path)))
    expect_identical(result$metadata$method, "cluster-vs-rest-rank-auc")
})

test_that("analysis v2 background jobs expose progress and a completed result", {
    skip_if_not_installed("processx")
    api <- analysis_v2_api_env()
    fixture <- analysis_v2_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    started <- api$gui_analysis_start_job(fixture$project, list(
        file = fixture$relative_file,
        files = fixture$relative_file,
        populationId = "root",
        reductionMethod = "pca",
        clusterMethod = "none",
        markers = c("FSC-A", "SSC-A", "CD3-A"),
        maxEvents = 100L,
        seed = 23L
    ))
    on.exit(try(api$gui_analysis_cancel_job(fixture$project, started$job_id), silent = TRUE), add = TRUE)

    deadline <- Sys.time() + 30
    repeat {
        status <- api$gui_analysis_job_status(fixture$project, started$job_id)
        if (status$state %in% c("completed", "failed", "cancelled")) break
        if (Sys.time() > deadline) fail("Background analysis did not finish in 30 seconds.")
        Sys.sleep(0.1)
    }
    expect_identical(status$state, "completed", info = if (!is.null(status$error)) status$error else status$message)
    expect_equal(nrow(status$result$events), 100L)
    expect_identical(status$result$metadata$reduction_method$id, "pca")
})

test_that("analysis v2 workspace records where a trajectory root was selected", {
    api <- analysis_v2_api_env()
    workspace <- api$gui_analysis_default_workspace()
    expect_identical(workspace$schema_version, 2L)
    expect_null(workspace$root_event_id)
    expect_null(workspace$root_population_id)
    expect_null(workspace$root_source_file)
})

test_that("analysis workspace import validates hierarchy and round-trips complete state", {
    api <- analysis_v2_api_env()
    fixture <- analysis_v2_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    workspace <- api$gui_analysis_default_workspace()
    workspace$seed <- 77L
    workspace$annotations <- list(list(id = "identity-1", label = "T cell"))
    workspace$plots[[1]]$color_palette <- "sunset"
    workspace$populations <- c(workspace$populations, list(list(
        id = "gate-1", name = "Cells", parent_id = "root", type = "ellipse",
        role = "positive", source_file = NULL, x = "FSC-A", y = "SSC-A",
        geometry = list(center_x = 60, center_y = 10, radius_x = 30, radius_y = 6)
    )))

    saved <- api$gui_analysis_write_workspace(fixture$project, workspace)
    restored <- api$gui_analysis_read_workspace(fixture$project)
    expect_identical(restored$seed, 77L)
    expect_identical(restored$plots[[1]]$color_palette, "sunset")
    expect_identical(restored$annotations[[1]]$label, "T cell")
    expect_equal(restored$populations[[2]]$geometry, saved$populations[[2]]$geometry)

    invalid <- workspace
    invalid$populations <- invalid$populations[-1]
    expect_error(api$gui_analysis_write_workspace(fixture$project, invalid), "missing its root")

    missing_parent <- workspace
    missing_parent$populations[[2]]$parent_id <- "missing"
    expect_error(api$gui_analysis_write_workspace(fixture$project, missing_parent), "existing parent")

    cycle <- workspace
    cycle$populations <- c(cycle$populations, list(list(
        id = "gate-2", name = "Child", parent_id = "gate-1", type = "range",
        role = NULL, source_file = NULL, x = "CD3-A", geometry = list(min = 1, max = 2)
    )))
    cycle$populations[[2]]$parent_id <- "gate-2"
    expect_error(api$gui_analysis_write_workspace(fixture$project, cycle), "cycle")

    invalid_geometry <- workspace
    invalid_geometry$populations[[2]]$geometry$radius_x <- 0
    expect_error(api$gui_analysis_write_workspace(fixture$project, invalid_geometry), "radii must be positive")

    dangling_plot <- workspace
    dangling_plot$plots[[1]]$population_id <- "missing"
    expect_error(api$gui_analysis_write_workspace(fixture$project, dangling_plot), "missing population")
})

test_that("analysis requests use the submitted workspace snapshot before autosave completes", {
    api <- analysis_v2_api_env()
    fixture <- analysis_v2_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    workspace <- api$gui_analysis_default_workspace()
    workspace$populations <- c(workspace$populations, list(list(
        id = "unsaved-gate", name = "Unsaved gate", parent_id = "root", type = "range",
        role = NULL, source_file = NULL, x = "FSC-A",
        geometry = list(min = 15, max = 45)
    )))

    statistics <- api$gui_analysis_population_statistics(
        fixture$project, fixture$relative_file, "unsaved-gate", "CD3-A", workspace
    )
    result <- api$gui_analysis_run_method(fixture$project, list(
        file = fixture$relative_file, populationId = "unsaved-gate",
        reductionMethod = "pca", clusterMethod = "none",
        markers = c("FSC-A", "SSC-A", "CD3-A"), maxEvents = 100L, seed = 17L,
        workspace = workspace
    ))

    expect_gt(statistics$count, 0L)
    expect_lt(statistics$count, statistics$total_count)
    expect_identical(result$metadata$population_id, "unsaved-gate")
    expect_equal(nrow(result$events), statistics$count)
    expect_identical(
        vapply(api$gui_analysis_read_workspace(fixture$project)$populations, `[[`, character(1), "id"),
        "root"
    )
})

test_that("analysis export folder route uses the native picker and rejects paths outside the project", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    endpoints <- unname(unlist(plumber::plumb(api_path)$endpoints, recursive = FALSE))
    match <- which(vapply(endpoints, function(endpoint) {
        identical(endpoint$path, "/analysis/export-folder") && "POST" %in% endpoint$verbs
    }, logical(1)))
    route <- endpoints[[match[1]]]$getFunc()
    route_env <- environment(route)
    original_picker <- get("gui_pick_project_directory", envir = route_env)
    on.exit(assign("gui_pick_project_directory", original_picker, envir = route_env), add = TRUE)
    project <- tempfile("analysis_export_picker_")
    output <- file.path(project, "spectreasy_outputs", "analysis", "exports")
    dir.create(output, recursive = TRUE)
    on.exit(unlink(project, recursive = TRUE, force = TRUE), add = TRUE)
    request <- new.env(parent = emptyenv())
    request$postBody <- jsonlite::toJSON(list(
        projectPath = project,
        current = "spectreasy_outputs/analysis/exports"
    ), auto_unbox = TRUE)
    observed <- list()
    assign("gui_pick_project_directory", function(initial_dir, allow_create = FALSE, prompt = NULL) {
        observed <<- list(initial_dir = initial_dir, allow_create = allow_create, prompt = prompt)
        output
    }, envir = route_env)
    response <- route(request)
    expect_true(response$success)
    expect_false(response$cancelled)
    expect_identical(response$path, "spectreasy_outputs/analysis/exports")
    expect_true(observed$allow_create)
    expect_match(observed$prompt, "analysis export")

    assign("gui_pick_project_directory", function(...) tempdir(), envir = route_env)
    rejected <- route(request)
    expect_false(rejected$success)
    expect_match(rejected$error, "outside the active project")
})

test_that("analysis job routes preserve the explicit project for polling and cancellation", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    endpoints <- unname(unlist(plumber::plumb(api_path)$endpoints, recursive = FALSE))
    job_routes <- Filter(function(endpoint) identical(endpoint$path, "/analysis/jobs"), endpoints)
    get_route <- job_routes[[which(vapply(job_routes, function(endpoint) "GET" %in% endpoint$verbs, logical(1)))[1]]]$getFunc()
    delete_route <- job_routes[[which(vapply(job_routes, function(endpoint) "DELETE" %in% endpoint$verbs, logical(1)))[1]]]$getFunc()
    route_env <- environment(get_route)
    project <- tempfile("analysis_job_route_")
    dir.create(project)
    on.exit(unlink(project, recursive = TRUE, force = TRUE), add = TRUE)
    observed <- list()
    original_status <- get("gui_analysis_job_status", envir = route_env)
    original_cancel <- get("gui_analysis_cancel_job", envir = route_env)
    on.exit(assign("gui_analysis_job_status", original_status, envir = route_env), add = TRUE)
    on.exit(assign("gui_analysis_cancel_job", original_cancel, envir = route_env), add = TRUE)
    assign("gui_analysis_job_status", function(root, job_id) {
        observed$status <<- c(root = root, job_id = job_id)
        list(state = "running")
    }, envir = route_env)
    assign("gui_analysis_cancel_job", function(root, job_id) {
        observed$cancel <<- c(root = root, job_id = job_id)
        list(state = "cancelled")
    }, envir = route_env)

    polled <- get_route(job_id = "20260724-deadbeef", project_path = project)
    request <- new.env(parent = emptyenv())
    request$postBody <- jsonlite::toJSON(list(projectPath = project), auto_unbox = TRUE)
    cancelled <- delete_route(request, job_id = "20260724-deadbeef")

    expect_true(polled$success)
    expect_true(cancelled$success)
    expect_identical(unname(observed$status), c(normalizePath(project), "20260724-deadbeef"))
    expect_identical(unname(observed$cancel), c(normalizePath(project), "20260724-deadbeef"))
})

test_that("analysis v2 persists hierarchical gates and calculates robust staining index", {
    api <- analysis_v2_api_env()
    fixture <- analysis_v2_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    workspace <- api$gui_analysis_default_workspace()
    workspace$source_path <- "samples"
    workspace$selected_file <- fixture$relative_file
    workspace$populations <- c(workspace$populations, list(
        list(id = "negative", name = "CD3 negative", parent_id = "root", type = "range", role = "negative", source_file = NULL, x = "CD3-A", geometry = list(min = -2, max = 5)),
        list(id = "positive", name = "CD3 positive", parent_id = "root", type = "range", role = "positive", source_file = NULL, x = "CD3-A", geometry = list(min = 5, max = 12))
    ))
    api$gui_analysis_write_workspace(fixture$project, workspace)

    restored <- api$gui_analysis_read_workspace(fixture$project)
    expect_identical(vapply(restored$populations, `[[`, character(1), "id"), c("root", "negative", "positive"))
    result <- api$gui_analysis_staining_index(fixture$project, fixture$relative_file, "CD3-A")
    expected_sd <- stats::mad(fixture$data[1:60, "CD3-A"], constant = 1.4826)
    expected <- (stats::median(fixture$data[61:120, "CD3-A"]) - stats::median(fixture$data[1:60, "CD3-A"])) / (2 * expected_sd)
    expect_equal(result$negative_count, 60L)
    expect_equal(result$positive_count, 60L)
    expect_equal(result$staining_index, expected, tolerance = 1e-6)
})

test_that("analysis v2 evaluates ellipse gates in raw channel coordinates", {
    api <- analysis_v2_api_env()
    fixture <- analysis_v2_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    workspace <- api$gui_analysis_default_workspace()
    workspace$populations <- c(workspace$populations, list(list(
        id = "ellipse", name = "Ellipse", parent_id = "root", type = "ellipse",
        role = NULL, source_file = NULL, x = "FSC-A", y = "SSC-A",
        geometry = list(center_x = 60, center_y = 10, radius_x = 30, radius_y = 6)
    )))
    api$gui_analysis_write_workspace(fixture$project, workspace)
    frame <- api$gui_analysis_flow_frame(fixture$project, fixture$relative_file)
    data <- flowCore::exprs(frame)
    observed <- api$gui_analysis_gate_mask(data, workspace, "ellipse", fixture$relative_file)
    expected <- ((data[, "FSC-A"] - 60) / 30)^2 + ((data[, "SSC-A"] - 10) / 6)^2 <= 1

    expect_identical(observed, unname(expected))
    expect_gt(sum(observed), 0L)
    expect_lt(sum(observed), nrow(data))
})

test_that("analysis v2 exports a reproducible population with provenance", {
    api <- analysis_v2_api_env()
    fixture <- analysis_v2_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    workspace <- api$gui_analysis_default_workspace()
    workspace$populations <- c(workspace$populations, list(
        list(id = "middle", name = "Middle", parent_id = "root", type = "rectangle", role = NULL, source_file = NULL, x = "FSC-A", y = "SSC-A", geometry = list(x_min = 21, x_max = 100, y_min = 1, y_max = 20))
    ))
    api$gui_analysis_write_workspace(fixture$project, workspace)
    body <- list(files = fixture$relative_file, populationId = "middle", format = "both", maxEvents = 25, seed = 99, outputFolder = "spectreasy_outputs/analysis/exports")
    result <- api$gui_analysis_export(fixture$project, body)
    records <- result$files
    exported <- vapply(records, `[[`, character(1), "path")

    expect_length(exported, 2L)
    expect_true(all(file.exists(file.path(fixture$project, exported))))
    csv_path <- file.path(fixture$project, exported[grepl("[.]csv$", exported)])
    fcs_path <- file.path(fixture$project, exported[grepl("[.]fcs$", exported)])
    csv <- data.table::fread(csv_path)
    expect_equal(nrow(csv), 25L)
    expect_length(unique(csv$event_id), 25L)
    header <- flowCore::read.FCSheader(fcs_path)[[1]]
    expect_identical(unname(header["SpectreasyPopulationPath"]), "All events/Middle")
    expect_identical(unname(header["SpectreasySeed"]), "99")
    expect_true(all(nzchar(vapply(records, `[[`, character(1), "sha256"))))

    statistics <- api$gui_analysis_export_statistics(fixture$project, list(
        files = fixture$relative_file, populationIds = c("root", "middle"),
        markers = c("FSC-A", "CD3-A"), outputFolder = "spectreasy_outputs/analysis/exports"
    ))
    expect_equal(statistics$rows, 4L)
    stats_table <- data.table::fread(file.path(fixture$project, statistics$path))
    expect_setequal(stats_table$population_id, c("root", "middle"))
    expect_setequal(stats_table$marker, c("FSC-A", "CD3-A"))
})

test_that("analysis v2 keeps pooled input stable under file reordering", {
    api <- analysis_v2_api_env()
    fixture <- analysis_v2_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    second_path <- file.path(fixture$project, "samples", "sample_2.fcs")
    second <- flowCore::flowFrame(fixture$data + 0.25)
    second_parameters <- flowCore::parameters(second)@data
    second_parameters$desc <- c("FSC", "SSC", "CD3")
    flowCore::parameters(second)@data <- second_parameters
    flowCore::write.FCS(second, second_path)
    files <- c(fixture$relative_file, "samples/sample_2.fcs")
    request <- list(
        populationId = "root", reductionMethod = "pca", clusterMethod = "none",
        markers = c("FSC-A", "SSC-A", "CD3-A"), maxEvents = 100L, seed = 31L
    )

    forward <- api$gui_analysis_run_method(fixture$project, c(request, list(files = files)))
    reversed <- api$gui_analysis_run_method(fixture$project, c(request, list(files = rev(files))))

    expect_identical(forward$metadata$source_files, sort(files, method = "radix"))
    expect_identical(reversed$metadata$source_files, forward$metadata$source_files)
    expect_identical(
        reversed$events[, c("source_file", "source_event_id")],
        forward$events[, c("source_file", "source_event_id")]
    )
})

test_that("analysis v2 export does not overwrite duplicate source basenames", {
    api <- analysis_v2_api_env()
    fixture <- analysis_v2_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    second_directory <- file.path(fixture$project, "alternate")
    dir.create(second_directory)
    second_path <- file.path(second_directory, "sample.fcs")
    file.copy(file.path(fixture$project, fixture$relative_file), second_path)

    result <- api$gui_analysis_export(fixture$project, list(
        files = c(fixture$relative_file, "alternate/sample.fcs"),
        populationId = "root", format = "csv", maxEvents = 10L, seed = 44L,
        outputFolder = "spectreasy_outputs/analysis/exports"
    ))
    paths <- vapply(result$files, `[[`, character(1), "path")

    expect_length(unique(paths), 2L)
    expect_true(all(file.exists(file.path(fixture$project, paths))))
    expect_true(all(vapply(paths, function(path) nrow(data.table::fread(file.path(fixture$project, path))) == 10L, logical(1))))
})

test_that("analysis v2 export reports files where a scoped population does not apply", {
    api <- analysis_v2_api_env()
    fixture <- analysis_v2_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    second_path <- file.path(fixture$project, "samples", "sample_2.fcs")
    file.copy(file.path(fixture$project, fixture$relative_file), second_path)
    workspace <- api$gui_analysis_default_workspace()
    workspace$populations <- c(workspace$populations, list(list(
        id = "first-only", name = "First only", parent_id = "root", type = "range",
        role = NULL, source_file = fixture$relative_file, x = "FSC-A",
        geometry = list(min = 1, max = 100)
    )))
    api$gui_analysis_write_workspace(fixture$project, workspace)

    result <- api$gui_analysis_export(fixture$project, list(
        files = c(fixture$relative_file, "samples/sample_2.fcs"),
        populationId = "first-only", format = "both", maxEvents = 0L, seed = 45L,
        outputFolder = "spectreasy_outputs/analysis/exports"
    ))

    expect_length(result$files, 2L)
    expect_length(result$skipped_files, 1L)
    expect_identical(result$skipped_files[[1]]$source_file, "samples/sample_2.fcs")
    expect_match(result$skipped_files[[1]]$reason, "no events")
})

test_that("analysis v2 marks a stopped worker as failed instead of polling forever", {
    api <- analysis_v2_api_env()
    project <- tempfile("analysis_dead_worker_")
    dir.create(project)
    on.exit(unlink(project, recursive = TRUE, force = TRUE), add = TRUE)
    job_id <- "20260724-deadbeef"
    directory <- api$gui_analysis_job_directory(project, job_id, create = TRUE)
    dir.create(directory)
    jsonlite::write_json(
        list(state = "queued", message = "Queued", updated_at = as.character(Sys.time()), error = NULL),
        file.path(directory, "status.json"), auto_unbox = TRUE, null = "null"
    )
    writeLines("worker bootstrap failed", file.path(directory, "stderr.log"))
    fake_process <- list(
        is_alive = function() FALSE,
        get_exit_status = function() 42L
    )
    assign(job_id, list(process = fake_process, root = project, directory = directory), envir = api$.gui_analysis_jobs)
    on.exit(if (exists(job_id, envir = api$.gui_analysis_jobs, inherits = FALSE)) {
        rm(list = job_id, envir = api$.gui_analysis_jobs)
    }, add = TRUE)

    status <- api$gui_analysis_job_status(project, job_id)

    expect_identical(status$state, "failed")
    expect_match(status$error, "worker bootstrap failed")
})

test_that("analysis v2 advertises maintained executable adapters and method settings", {
    api <- analysis_v2_api_env()
    methods <- api$gui_analysis_method_registry()
    by_id <- setNames(methods, vapply(methods, `[[`, character(1), "id"))
    expect_true(by_id$pca$available)
    expect_true(by_id$pca$adapter_verified)
    expect_identical(by_id$pca$availability_state, "ready")
    expect_true(by_id$pca$supports_3d)
    expect_true(by_id$slingshot$adapter_verified)
    expect_identical(
        by_id$slingshot$available,
        requireNamespace("slingshot", quietly = TRUE) &&
            requireNamespace("DelayedMatrixStats", quietly = TRUE)
    )
    expect_true(by_id$tscan$adapter_verified)
    expect_identical(by_id$tscan$available, requireNamespace("TSCAN", quietly = TRUE))
    expect_true(by_id$phate$adapter_verified)
    expect_true(by_id$hsne$adapter_verified)
    expect_true(by_id$palantir$adapter_verified)
    expect_true(by_id$`paga-dpt`$adapter_verified)
    builtin_ready <- isTRUE(
        api$gui_analysis_python_packages()[["spectreasy_builtin"]]$available
    )
    expect_true(by_id$wanderlust$visible)
    expect_identical(by_id$wanderlust$available, builtin_ready)
    expect_identical(by_id$wanderlust$package, "spectreasy_builtin")
    expect_true(by_id$wishbone$visible)
    expect_identical(by_id$wishbone$available, builtin_ready)
    expect_identical(by_id$wishbone$package, "spectreasy_builtin")
    expect_true(all(vapply(methods, function(method) is.list(method$parameters), logical(1))))
    expect_setequal(
        vapply(by_id$umap$parameters, `[[`, character(1), "id"),
        c("neighbors", "min_dist", "spread", "metric", "epochs", "learning_rate", "repulsion", "negative_samples")
    )
    expect_identical(by_id$flowsom$outputs, "cluster_labels")
    expect_identical(by_id$flowsom$pipeline, "flowsom")
    expect_false(by_id$flowsom$supports_3d)
    expect_identical(by_id$phenograph$outputs, "cluster_labels")
    expect_false(any(grepl("pca", c(by_id$flowsom$pipeline, by_id$phenograph$pipeline), ignore.case = TRUE)))
    expect_match(by_id$dpt$citation, "Haghverdi")
    expect_identical(by_id$dpt$prerequisites, "diffusion-map")
    expect_identical(by_id$dpt$pipeline, c("diffusion-map", "dpt"))
    expect_setequal(by_id$dpt$outputs, c("embedding", "pseudotime"))
})

test_that("analysis v2 writes cited PCA coordinates and publication plot formats", {
    api <- analysis_v2_api_env()
    fixture <- analysis_v2_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    result <- api$gui_analysis_run_method(fixture$project, list(
        file = fixture$relative_file, populationId = "root", method = "pca",
        markers = c("FSC-A", "SSC-A", "CD3-A"), maxEvents = 100, seed = 123, cofactor = 150
    ))

    expect_equal(nrow(result$events), 100L)
    expect_identical(result$metadata$method$id, "pca")
    expect_match(result$metadata$method$citation, "Pearson")
    formats <- vapply(result$metadata$plot_files, `[[`, character(1), "format")
    expect_true(all(c("png", "pdf", "html") %in% formats))
    if (requireNamespace("svglite", quietly = TRUE)) expect_true("svg" %in% formats)
    paths <- vapply(result$metadata$plot_files, `[[`, character(1), "path")
    expect_true(all(file.exists(file.path(fixture$project, paths))))
    expect_true(file.exists(file.path(fixture$project, result$metadata$events_file)))
    expect_true(file.exists(file.path(fixture$project, result$metadata_file)))
    expect_identical(
        vapply(result$metadata$marker_columns, `[[`, character(1), "marker"),
        c("FSC", "SSC", "CD3")
    )
    expect_identical(
        vapply(result$metadata$marker_columns, `[[`, character(1), "channel"),
        c("FSC-A", "SSC-A", "CD3-A")
    )
    expect_identical(
        vapply(result$metadata$marker_columns, `[[`, character(1), "column"),
        c("marker_1", "marker_2", "marker_3")
    )
    expect_true(all(c("marker_1", "marker_2", "marker_3") %in% names(result$events)))
})

test_that("flow-specific cell identities use signed marker evidence and preserve uncertainty", {
    api <- analysis_v2_api_env()
    fixture <- analysis_v2_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    result <- api$gui_analysis_run_method(fixture$project, list(
        file = fixture$relative_file, populationId = "root", method = "pca",
        markers = c("FSC-A", "SSC-A", "CD3-A"), maxEvents = 120, seed = 123, cofactor = 1
    ))
    signatures <- list(
        list(name = "CD3 low", color = "#197783", positive_markers = character(), negative_markers = "CD3"),
        list(name = "CD3 high", color = "#d06d32", positive_markers = "CD3", negative_markers = character())
    )
    annotated <- api$gui_analysis_annotate_result(fixture$project, list(
        analysisId = result$metadata$analysis_id,
        signatures = signatures,
        minScore = 0.55,
        minMargin = 0.05
    ))

    expect_equal(nrow(annotated$events), 120L)
    expect_setequal(unique(annotated$events$predicted_identity), c("CD3 low", "CD3 high"))
    expect_true(all(annotated$events$identity_score >= 0.55))
    expect_true(all(annotated$events$identity_margin >= 0.05))
    expect_equal(annotated$metadata$assigned_count, 120L)
    expect_equal(annotated$metadata$unassigned_count, 0L)
    expect_true(file.exists(file.path(fixture$project, annotated$metadata$model$path)))
    expect_true(file.exists(file.path(fixture$project, annotated$metadata$scores$path)))
    expect_true(file.exists(file.path(fixture$project, annotated$metadata$metadata$path)))

    conservative <- api$gui_analysis_score_identities(
        cbind(CD3 = c(-1, 0, 9, 10)),
        signatures,
        min_score = 0.99,
        min_margin = 0.9
    )
    expect_true(all(conservative$labels == "Unassigned"))
})

test_that("cell identity signatures reject ambiguous or unavailable marker rules", {
    api <- analysis_v2_api_env()
    matrix <- cbind(`CD3-A` = c(-1, 1), `CD19-A` = c(1, -1))
    expect_error(api$gui_analysis_score_identities(matrix, list(
        list(name = "A", positive_markers = "CD3-A", negative_markers = "CD3-A"),
        list(name = "B", positive_markers = "CD19-A")
    )), "both positive and negative")
    expect_error(api$gui_analysis_score_identities(matrix, list(
        list(name = "A", positive_markers = "CD4-A"),
        list(name = "B", positive_markers = "CD19-A")
    )), "unavailable marker")
})

test_that("analysis v2 only emits a third coordinate when the fitted input supports it", {
    api <- analysis_v2_api_env()
    fixture <- analysis_v2_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    result <- api$gui_analysis_run_method(fixture$project, list(
        file = fixture$relative_file, populationId = "root", method = "pca",
        markers = c("FSC-A", "SSC-A"), maxEvents = 100, seed = 123, cofactor = 150
    ))

    expect_identical(result$metadata$coordinate_count, 2L)
    expect_identical(result$metadata$coordinate_labels, c("PC 1", "PC 2"))
    expect_false("dimension_3" %in% names(result$events))
})

test_that("every executable reduction and trajectory adapter satisfies the small-population contract", {
    api <- analysis_v2_api_env()
    fixture <- analysis_v2_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    methods <- api$gui_analysis_method_registry()
    executable <- Filter(function(method) isTRUE(method$available), methods)
    expect_true(all(vapply(executable, `[[`, logical(1), "adapter_verified")))

    for (method in Filter(function(method) !identical(method$family, "clustering"), executable)) {
        body <- list(
            file = fixture$relative_file,
            populationId = "root",
            method = method$id,
            markers = c("FSC-A", "SSC-A", "CD3-A"),
            maxEvents = 100,
            seed = 731,
            rootEventId = 1L,
            cofactor = 150,
            neighbors = 8L,
            clusters = 4L,
            perplexity = 10
        )
        if (identical(method$family, "trajectory")) {
            if ("clustering" %in% method$prerequisites) body$clusterMethod <- "flowsom"
            if ("reduction" %in% method$prerequisites) body$reductionMethod <- "pca"
        }
        result <- api$gui_analysis_run_method(fixture$project, body)
        expect_identical(result$metadata$method$id, method$id, info = method$id)
        expect_equal(nrow(result$events), 100L, info = method$id)
        expect_true(all(is.finite(result$events$dimension_1)), info = method$id)
        expect_true(all(is.finite(result$events$dimension_2)), info = method$id)
        expected_dimensions <- if (identical(method$id, "hsne")) 2L else 3L
        if (expected_dimensions == 3L) expect_true(all(is.finite(result$events$dimension_3)), info = method$id)
        expect_identical(result$metadata$coordinate_count, expected_dimensions, info = method$id)
        expect_length(result$metadata$coordinate_labels, expected_dimensions)
        expect_true(length(result$metadata$pipeline) >= 2L, info = method$id)
        if (identical(method$family, "trajectory")) {
            expect_true("pseudotime" %in% names(result$events), info = method$id)
            expect_true(all(is.finite(result$events$pseudotime)), info = method$id)
            expect_true(min(result$events$pseudotime) >= 0, info = method$id)
            expect_true(max(result$events$pseudotime) <= 1, info = method$id)
        }
        if (identical(method$id, "dpt")) {
            expect_true("pseudotime" %in% names(result$events))
            expect_true(all(is.finite(result$events$pseudotime)))
            expect_identical(
                vapply(result$metadata$pipeline, `[[`, character(1), "id"),
                c("input", "diffusion-map", "dpt")
            )
            expect_length(result$metadata$intermediate_files, 1L)
            intermediate <- result$metadata$intermediate_files[[1]]
            expect_identical(intermediate$id, "diffusion_map")
            expect_true(file.exists(file.path(fixture$project, intermediate$path)))
            diffusion <- data.table::fread(file.path(fixture$project, intermediate$path))
            expect_identical(names(diffusion), c("event_id", "diffusion_1", "diffusion_2", "diffusion_3"))
            expect_equal(nrow(diffusion), 100L)
        }
        repeated <- api$gui_analysis_run_method(fixture$project, body)
        comparable <- intersect(
            c("event_id", "dimension_1", "dimension_2", "dimension_3", "cluster_id", "pseudotime"),
            names(result$events)
        )
        expect_equal(
            repeated$events[, comparable, drop = FALSE],
            result$events[, comparable, drop = FALSE],
            tolerance = 1e-8,
            info = paste(method$id, "same-seed reproducibility")
        )
        if (identical(method$family, "trajectory") && !identical(method$id, "dpt")) {
            expect_true(repeated$metadata$artifacts$trajectory$reused, info = method$id)
        }
    }
})

test_that("every executable reduction preserves a small synthetic continuum", {
    api <- analysis_v2_api_env()
    fixture <- analysis_v2_lineage_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    reductions <- Filter(
        function(method) identical(method$family, "reduction") && isTRUE(method$available),
        api$gui_analysis_method_registry()
    )
    expect_true(length(reductions) >= 1L)
    expect_true("pca" %in% vapply(reductions, `[[`, character(1), "id"))

    set.seed(919)
    pairs <- cbind(sample.int(180L, 1200L, replace = TRUE), sample.int(180L, 1200L, replace = TRUE))
    pairs <- pairs[pairs[, 1] != pairs[, 2], , drop = FALSE]
    for (method in reductions) {
        result <- api$gui_analysis_run_method(fixture$project, list(
            file = fixture$relative_file,
            populationId = "root",
            reductionMethod = method$id,
            clusterMethod = "none",
            markers = fixture$markers,
            maxEvents = 180L,
            seed = 919L,
            neighbors = 12L,
            perplexity = 18
        ))
        dimensions <- grep("^dimension_", names(result$events), value = TRUE)
        embedding <- as.matrix(result$events[, dimensions, drop = FALSE])
        truth <- fixture$truth[result$events$source_event_id]
        truth_distance <- abs(truth[pairs[, 1]] - truth[pairs[, 2]])
        map_distance <- sqrt(rowSums(
            (embedding[pairs[, 1], , drop = FALSE] - embedding[pairs[, 2], , drop = FALSE])^2
        ))
        preservation <- suppressWarnings(stats::cor(truth_distance, map_distance, method = "spearman"))
        minimum_preservation <- if (identical(method$id, "hsne")) 0.20 else 0.35
        expect_true(
            is.finite(preservation) && preservation > minimum_preservation,
            info = paste(method$id, "continuum distance preservation:", round(preservation, 3))
        )
    }
})

test_that("every executable clustering method recovers separated populations", {
    api <- analysis_v2_api_env()
    fixture <- analysis_v2_cluster_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    clustering <- Filter(
        function(method) identical(method$family, "clustering") && isTRUE(method$available),
        api$gui_analysis_method_registry()
    )
    expect_gte(length(clustering), 2L)

    for (method in clustering) {
        result <- api$gui_analysis_run_method(fixture$project, list(
            file = fixture$relative_file,
            populationId = "root",
            clusterMethod = method$id,
            reductionMethod = "pca",
            markers = fixture$markers,
            maxEvents = 180L,
            seed = 922L,
            neighbors = 15L,
            clusters = 3L
        ))
        truth <- fixture$truth[result$events$source_event_id]
        recovery <- analysis_v2_adjusted_rand(truth, result$events$cluster_id)
        expect_true(
            is.finite(recovery) && recovery > 0.80,
            info = paste(method$id, "adjusted Rand recovery:", round(recovery, 3))
        )
    }
})

test_that("every executable trajectory orders a small seeded lineage from its root", {
    api <- analysis_v2_api_env()
    fixture <- analysis_v2_lineage_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    trajectories <- Filter(
        function(method) identical(method$family, "trajectory") && isTRUE(method$available),
        api$gui_analysis_method_registry()
    )
    expect_true(length(trajectories) >= 1L)

    for (method in trajectories) {
        body <- list(
            file = fixture$relative_file,
            populationId = "root",
            method = method$id,
            markers = fixture$markers,
            maxEvents = 180L,
            seed = 920L,
            rootEventId = 1L,
            neighbors = 12L,
            clusters = 5L,
            perplexity = 18
        )
        if ("clustering" %in% method$prerequisites) body$clusterMethod <- "flowsom"
        if ("reduction" %in% method$prerequisites) body$reductionMethod <- "pca"
        result <- api$gui_analysis_run_method(fixture$project, body)
        truth <- fixture$truth[result$events$source_event_id]
        ordering <- suppressWarnings(stats::cor(truth, result$events$pseudotime, method = "spearman"))
        expect_true(
            is.finite(ordering) && ordering > 0.60,
            info = paste(method$id, "rooted lineage ordering:", round(ordering, 3))
        )
        root_pseudotime <- result$events$pseudotime[[which.min(truth)]]
        expect_true(
            is.finite(root_pseudotime) && root_pseudotime <= 0.12,
            info = paste(method$id, "root pseudotime:", round(root_pseudotime, 3))
        )
    }
})

test_that("clustering and dimensional reduction are separate reusable artifacts", {
    api <- analysis_v2_api_env()
    fixture <- analysis_v2_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    base <- list(
        file = fixture$relative_file,
        populationId = "root",
        clusterMethod = "flowsom",
        markers = c("FSC-A", "SSC-A", "CD3-A"),
        maxEvents = 100,
        seed = 731,
        cofactor = 150,
        neighbors = 8L,
        clusters = 4L,
        perplexity = 10
    )

    pca <- api$gui_analysis_run_method(fixture$project, c(base, list(reductionMethod = "pca")))
    umap <- api$gui_analysis_run_method(fixture$project, c(base, list(reductionMethod = "umap")))

    expect_identical(pca$metadata$display_name, "FlowSOM → PCA")
    expect_identical(umap$metadata$display_name, "FlowSOM → UMAP")
    expect_identical(pca$metadata$cluster_method$id, "flowsom")
    expect_identical(pca$metadata$reduction_method$id, "pca")
    expect_identical(umap$metadata$reduction_method$id, "umap")
    expect_identical(
        pca$metadata$artifacts$clustering$id,
        umap$metadata$artifacts$clustering$id
    )
    expect_false(pca$metadata$artifacts$clustering$reused)
    expect_true(umap$metadata$artifacts$clustering$reused)
    expect_false(identical(
        pca$metadata$artifacts$embedding$id,
        umap$metadata$artifacts$embedding$id
    ))
    expect_identical(pca$events$event_id, umap$events$event_id)
    expect_identical(pca$events$cluster_id, umap$events$cluster_id)
    expect_true(length(unique(pca$events$cluster_id)) >= 2L)
    expect_identical(
        vapply(pca$metadata$pipeline, `[[`, character(1), "id"),
        c("input", "flowsom", "pca")
    )
    expect_false(any(grepl("pca-visualization", vapply(pca$metadata$pipeline, `[[`, character(1), "id"), fixed = TRUE)))

    artifacts <- c(
        pca$metadata$artifacts$clustering[c("object", "values")],
        pca$metadata$artifacts$embedding[c("object", "values")],
        umap$metadata$artifacts$embedding[c("object", "values")]
    )
    artifact_paths <- vapply(artifacts, function(artifact) artifact$path, character(1))
    expect_true(all(file.exists(file.path(fixture$project, artifact_paths))))
    skipped <- api$gui_analysis_run_method(fixture$project, c(
        base[c("file", "populationId", "markers", "maxEvents", "seed", "cofactor")],
        list(clusterMethod = "none", reductionMethod = "pca")
    ))
    expect_null(skipped$metadata$artifacts$clustering)
    expect_false("cluster_id" %in% names(skipped$events))
    expect_error(
        api$gui_analysis_run_method(fixture$project, list(method = "flowsom")),
        "Choose a dimensional-reduction method"
    )
})

test_that("each clustering adapter can be paired with an independent map", {
    api <- analysis_v2_api_env()
    fixture <- analysis_v2_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)

    for (cluster_method in c("flowsom", "phenograph")) {
        result <- api$gui_analysis_run_method(fixture$project, list(
            file = fixture$relative_file,
            populationId = "root",
            clusterMethod = cluster_method,
            reductionMethod = "pca",
            markers = c("FSC-A", "SSC-A", "CD3-A"),
            maxEvents = 100,
            seed = 947,
            cofactor = 150,
            neighbors = 8L,
            clusters = 4L
        ))
        expect_identical(result$metadata$cluster_method$id, cluster_method)
        expect_identical(result$metadata$reduction_method$id, "pca")
        expect_true("cluster_id" %in% names(result$events))
        expect_true(length(unique(result$events$cluster_id)) >= 2L)
        expect_true(all(is.finite(result$events$dimension_1)))
        expect_identical(
            vapply(result$metadata$pipeline, `[[`, character(1), "id"),
            c("input", cluster_method, "pca")
        )
    }
})

test_that("advanced method settings are validated and affect artifact identity", {
    api <- analysis_v2_api_env()
    umap <- api$gui_analysis_method("umap")
    expect_error(
        api$gui_analysis_method_parameters(umap, list(neighbors = 1L)),
        "below its minimum"
    )
    expect_error(
        api$gui_analysis_method_parameters(umap, list(metric = "not-a-distance")),
        "unsupported value"
    )
    expect_error(
        api$gui_analysis_method_parameters(umap, list(unknown = 4)),
        "no advanced setting"
    )
    expect_error(
        api$gui_analysis_method_parameters(
            api$gui_analysis_method("wishbone"),
            list(neighbors = 20L, candidate_neighbors = 20L)
        ),
        "must exceed retained"
    )
    defaults <- api$gui_analysis_method_parameters(umap)
    changed <- defaults
    changed$neighbors <- defaults$neighbors + 1L
    expect_false(identical(
        api$gui_analysis_manifest_id("embedding", list(method = "umap", parameters = defaults)),
        api$gui_analysis_manifest_id("embedding", list(method = "umap", parameters = changed))
    ))
})
