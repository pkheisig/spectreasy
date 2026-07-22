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

test_that("analysis v2 advertises only execution-verified method adapters", {
    api <- analysis_v2_api_env()
    methods <- api$gui_analysis_method_registry()
    by_id <- setNames(methods, vapply(methods, `[[`, character(1), "id"))
    expect_true(by_id$pca$available)
    expect_true(by_id$pca$adapter_verified)
    expect_false(by_id$slingshot$available)
    expect_false(by_id$slingshot$adapter_verified)
    expect_false(by_id$tscan$available)
    expect_false(by_id$wishbone$available)
    expect_match(by_id$dpt$citation, "Haghverdi")
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
})
