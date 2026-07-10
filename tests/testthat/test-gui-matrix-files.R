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

testthat::test_that("GUI spectrum gating keeps zero-event gate results empty", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) {
        api_path <- system.file("api/gui_api.R", package = "spectreasy")
    }
    testthat::skip_if_not(file.exists(api_path))

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)

    expr <- data.frame(
        `FSC-A` = c(10, 20, 30),
        `SSC-A` = c(10, 20, 30),
        `Peak-A` = c(100, 200, 300),
        check.names = FALSE
    )
    empty_scatter_gate <- list(
        xChannel = "FSC-A",
        yChannel = "SSC-A",
        vertices = list(
            list(x = 1000, y = 1000),
            list(x = 1100, y = 1000),
            list(x = 1100, y = 1100),
            list(x = 1000, y = 1100)
        )
    )
    empty_positive_gate <- list(
        mode = "positive_1d",
        xChannel = "Peak-A",
        vertices = list(
            list(x = 1000, y = 0),
            list(x = 1100, y = 0)
        )
    )

    testthat::expect_equal(nrow(api_env$gate_apply_polygon_gate(expr, empty_scatter_gate)), 0L)
    testthat::expect_equal(nrow(api_env$gate_apply_positive_gate(expr, empty_positive_gate, "Peak-A")), 0L)
})

testthat::test_that("GUI histogram autogating creates positive and negative intervals with package defaults", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) {
        api_path <- system.file("api/gui_api.R", package = "spectreasy")
    }
    testthat::skip_if_not(file.exists(api_path))

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)

    captured <- NULL
    fake_compute <- function(peak_vals,
                             sample_type,
                             histogram_pct_beads,
                             histogram_direction_beads,
                             histogram_pct_cells,
                             histogram_direction_cells,
                             is_viability) {
        captured <<- list(
            sample_type = sample_type,
            histogram_pct_beads = histogram_pct_beads,
            histogram_direction_beads = histogram_direction_beads,
            histogram_pct_cells = histogram_pct_cells,
            histogram_direction_cells = histogram_direction_cells,
            is_viability = is_viability
        )
        vals_log <- log10(pmax(peak_vals, 1))
        attr(vals_log, "neg_log_min") <- 1
        attr(vals_log, "neg_log_max") <- 2
        attr(vals_log, "negative_gate_present") <- TRUE
        list(vals_log = vals_log, gate_min = 1000, gate_max = 10000)
    }

    ranges <- api_env$gate_histogram_autogate_ranges(
        peak_vals = seq_len(100),
        sample_type = "cells",
        is_viability = TRUE,
        compute_fun = fake_compute
    )

    testthat::expect_equal(ranges$positive, c(1000, 10000))
    testthat::expect_equal(ranges$negative, c(10, 100))
    testthat::expect_equal(api_env$gate_normalize_control_type("bead"), "beads")
    testthat::expect_equal(api_env$gate_normalize_control_type("cell"), "cells")
    testthat::expect_equal(captured$histogram_pct_beads, 0.98)
    testthat::expect_equal(captured$histogram_direction_beads, "right")
    testthat::expect_equal(captured$histogram_pct_cells, 0.35)
    testthat::expect_equal(captured$histogram_direction_cells, "right")
    testthat::expect_true(captured$is_viability)
})

testthat::test_that("GUI spectrum renderer keeps an empty plot for zero events", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) {
        api_path <- system.file("api/gui_api.R", package = "spectreasy")
    }
    testthat::skip_if_not(file.exists(api_path))

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)

    det_info <- list(
        names = c("B1-A", "YG1-A"),
        labels = c("B1-A", "YG1-A")
    )
    expr <- data.frame(`B1-A` = numeric(), `YG1-A` = numeric(), check.names = FALSE)
    image <- api_env$gate_build_spectrum_plot(expr, c("B1-A", "YG1-A"), det_info)

    testthat::expect_true(is.character(image))
    testthat::expect_match(image, "^data:image/png;base64,")
})
