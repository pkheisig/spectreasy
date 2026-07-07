testthat::test_that("GUI matrix listing recurses and only shows matrix-named CSVs", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) {
        api_path <- system.file("api/gui_api.R", package = "spectreasy")
    }
    testthat::skip_if_not(file.exists(api_path))

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)

    matrix_dir <- tempfile("spectreasy_gui_matrices_")
    dir.create(file.path(matrix_dir, "nested", "deeper"), recursive = TRUE)

    matrix_df <- data.frame(
        Marker = c("FITC", "PE"),
        `B1-A` = c(1, 0.1),
        `YG1-A` = c(0.1, 1),
        check.names = FALSE
    )
    utils::write.csv(matrix_df, file.path(matrix_dir, "fcs_mapping.csv"), row.names = FALSE)
    utils::write.csv(matrix_df, file.path(matrix_dir, "nested", "deeper", "scc_reference_matrix.csv"), row.names = FALSE)
    utils::write.csv(matrix_df, file.path(matrix_dir, "nested", "deeper", "comp-mat rix.csv"), row.names = FALSE)

    withr::local_options(list(spectreasy.matrix_dir = matrix_dir))

    files <- api_env$list_matrix_csv_files()
    testthat::expect_equal(
        files,
        c("nested/deeper/comp-mat rix.csv", "nested/deeper/scc_reference_matrix.csv")
    )

    paths <- file.path(api_env$get_matrix_dir(), files)
    keep <- vapply(paths, api_env$is_probably_matrix_csv, logical(1))
    testthat::expect_true(all(keep))
    testthat::expect_equal(
        normalizePath(api_env$matrix_path("nested/deeper/scc_reference_matrix.csv"), mustWork = FALSE),
        normalizePath(file.path(matrix_dir, "nested", "deeper", "scc_reference_matrix.csv"), mustWork = FALSE)
    )
    testthat::expect_error(api_env$matrix_path("../escape_matrix.csv"), "Invalid matrix filename")
})
