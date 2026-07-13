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

testthat::test_that("GUI gate file listing handles an empty project", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) {
        api_path <- system.file("api/gui_api.R", package = "spectreasy")
    }
    testthat::skip_if_not(file.exists(api_path))

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)

    project_dir <- tempfile("spectreasy_empty_gui_project_")
    scc_dir <- file.path(project_dir, "scc")
    dir.create(scc_dir, recursive = TRUE)
    withr::local_options(list(
        spectreasy.gating_scc_dir = scc_dir,
        spectreasy.gating_control_file = file.path(project_dir, "fcs_mapping.csv")
    ))

    mapping <- api_env$gate_read_mapping()
    testthat::expect_s3_class(mapping, "data.frame")
    testthat::expect_equal(nrow(mapping), 0L)
    testthat::expect_named(
        mapping,
        c(
            "filename", "fluorophore", "marker", "channel", "control.type",
            "universal.negative", "is.viability", "file_exists", "is_af",
            "uses_histogram_gates", "uses_negative_histogram_gate", "is_viability", "id"
        )
    )
})

testthat::test_that("GUI disables bead SCC negative histograms when AF_beads is mapped", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    testthat::skip_if_not(file.exists(api_path))

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)
    project_dir <- tempfile("spectreasy_bead_negative_gui_")
    scc_dir <- file.path(project_dir, "scc")
    dir.create(scc_dir, recursive = TRUE)
    filenames <- c("PE beads.fcs", "AF beads.fcs", "FITC cells.fcs")
    file.create(file.path(scc_dir, filenames))
    mapping <- data.frame(
        filename = filenames,
        fluorophore = c("PE", "AF_beads", "FITC"),
        marker = c("CD4", "Bead background", "CD8"),
        channel = c("YG1-A", "YG1-A", "B1-A"),
        control.type = c("beads", "beads", "cells"),
        universal.negative = c("AF beads.fcs", "", "AF"),
        is.viability = FALSE,
        stringsAsFactors = FALSE
    )
    mapping_path <- file.path(project_dir, "fcs_mapping.csv")
    utils::write.csv(mapping, mapping_path, row.names = FALSE)
    withr::local_options(list(
        spectreasy.gating_scc_dir = scc_dir,
        spectreasy.gating_control_file = mapping_path
    ))

    resolved <- api_env$gate_read_mapping()
    by_file <- split(resolved, resolved$filename)
    testthat::expect_false(by_file[["PE beads.fcs"]]$uses_negative_histogram_gate)
    testthat::expect_true(by_file[["PE beads.fcs"]]$uses_histogram_gates)
    testthat::expect_false(by_file[["AF beads.fcs"]]$uses_histogram_gates)
    testthat::expect_true(by_file[["FITC cells.fcs"]]$uses_negative_histogram_gate)
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
        peak_vals = seq_len(20000),
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

testthat::test_that("GUI histogram autogating repairs zero-event intervals", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) {
        api_path <- system.file("api/gui_api.R", package = "spectreasy")
    }
    testthat::skip_if_not(file.exists(api_path))

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)
    fake_compute <- function(peak_vals, ...) {
        vals_log <- log10(pmax(peak_vals, 1))
        attr(vals_log, "neg_log_min") <- 8
        attr(vals_log, "neg_log_max") <- 9
        attr(vals_log, "negative_gate_present") <- TRUE
        list(vals_log = vals_log, gate_min = 1e8, gate_max = 1e9)
    }

    values <- seq_len(1000)
    ranges <- api_env$gate_histogram_autogate_ranges(
        peak_vals = values,
        sample_type = "cells",
        compute_fun = fake_compute
    )

    testthat::expect_gte(sum(values >= ranges$positive[1] & values <= ranges$positive[2]), 10L)
    testthat::expect_gte(sum(values >= ranges$negative[1] & values <= ranges$negative[2]), 10L)
    testthat::expect_equal(ranges$positive, c(951, 1000))
    testthat::expect_equal(ranges$negative, c(1, 50))
})

testthat::test_that("GUI histogram autogating can generate positive-only bead gates", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    testthat::skip_if_not(file.exists(api_path))
    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)
    fake_compute <- function(peak_vals, ...) {
        vals_log <- log10(pmax(peak_vals, 1))
        attr(vals_log, "negative_gate_present") <- FALSE
        list(vals_log = vals_log, gate_min = 500, gate_max = 1000)
    }

    ranges <- api_env$gate_histogram_autogate_ranges(
        peak_vals = seq_len(1000), sample_type = "beads",
        compute_fun = fake_compute, include_negative = FALSE
    )

    testthat::expect_equal(ranges$positive, c(500, 1000))
    testthat::expect_null(ranges$negative)
})

testthat::test_that("GUI histogram autogates serialize as scalar 1D intervals", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) {
        api_path <- system.file("api/gui_api.R", package = "spectreasy")
    }
    testthat::skip_if_not(file.exists(api_path))

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)

    interval <- api_env$gate_histogram_interval(
        type = "positive",
        filename = "PE-Cy7.fcs",
        peak = "R1-A",
        limits = c(55000, 497000)
    )
    wire <- jsonlite::fromJSON(
        jsonlite::toJSON(interval, auto_unbox = FALSE),
        simplifyVector = FALSE
    )

    testthat::expect_identical(wire$type, "positive")
    testthat::expect_identical(wire$mode, "positive_1d")
    testthat::expect_identical(wire$xChannel, "R1-A")
    testthat::expect_identical(wire$vertices[[1]]$x, 55000L)
    testthat::expect_identical(wire$vertices[[2]]$x, 497000L)
    testthat::expect_identical(wire$vertices[[1]]$y, 0L)
    testthat::expect_identical(wire$vertices[[2]]$y, 0L)
})

testthat::test_that("GUI histogram autogating preserves existing intervals", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) {
        api_path <- system.file("api/gui_api.R", package = "spectreasy")
    }
    testthat::skip_if_not(file.exists(api_path))
    testthat::skip_if_not_installed("flowCore")

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)

    api_env$gate_read_mapping <- function() {
        data.frame(
            filename = c("control.fcs", "unstained.fcs"),
            channel = c("R1-A", "R1-A"),
            control.type = c("cells", "cells"),
            is_af = c(FALSE, TRUE),
            uses_histogram_gates = c(TRUE, FALSE),
            uses_negative_histogram_gate = c(TRUE, FALSE),
            is_viability = c(FALSE, FALSE),
            stringsAsFactors = FALSE
        )
    }
    api_env$gate_read_fcs <- function(path) {
        flowCore::flowFrame(matrix(seq_len(100), ncol = 1, dimnames = list(NULL, "R1-A")))
    }
    api_env$gate_apply_polygon_gate <- function(expr, gate) expr
    api_env$gate_histogram_autogate_ranges <- function(...) {
        list(positive = c(60, 100), negative = c(1, 40))
    }

    polygon <- list(
        mode = "scatter",
        vertices = list(
            list(x = 0, y = 0), list(x = 1, y = 0),
            list(x = 1, y = 1), list(x = 0, y = 1)
        )
    )
    gates <- list()
    gates[["cell:cells"]] <- polygon
    gates[["singlet:cells"]] <- polygon
    gates[["positive:control.fcs"]] <- api_env$gate_histogram_interval(
        "positive", "control.fcs", "R1-A", c(55, 95)
    )

    generated <- api_env$gate_autogenerate_histograms(gates)

    testthat::expect_setequal(names(generated$gates), "negative:control.fcs")
    testthat::expect_equal(generated$gates_generated, 1L)
    testthat::expect_equal(generated$gates_preserved, 1L)

    gates[["negative:control.fcs"]] <- api_env$gate_histogram_interval(
        "negative", "control.fcs", "R1-A", c(1, 40)
    )
    api_env$gate_read_fcs <- function(path) stop("FCS should not be read when both gates exist")

    preserved <- api_env$gate_autogenerate_histograms(gates)

    testthat::expect_length(preserved$gates, 0L)
    testthat::expect_equal(preserved$gates_generated, 0L)
    testthat::expect_equal(preserved$gates_preserved, 2L)
})

testthat::test_that("GUI histogram autogating omits bead negatives supplied by AF_beads", {
    api_path <- file.path(testthat::test_path("../.."), "inst", "api", "gui_api.R")
    if (!file.exists(api_path)) api_path <- system.file("api/gui_api.R", package = "spectreasy")
    testthat::skip_if_not(file.exists(api_path))
    testthat::skip_if_not_installed("flowCore")

    api_env <- new.env(parent = globalenv())
    source(api_path, local = api_env)
    api_env$gate_read_mapping <- function() data.frame(
        filename = "bead-scc.fcs", channel = "R1-A", control.type = "beads",
        is_af = FALSE, uses_histogram_gates = TRUE,
        uses_negative_histogram_gate = FALSE, is_viability = FALSE,
        stringsAsFactors = FALSE
    )
    api_env$gate_read_fcs <- function(path) {
        flowCore::flowFrame(matrix(seq_len(100), ncol = 1, dimnames = list(NULL, "R1-A")))
    }
    api_env$gate_apply_polygon_gate <- function(expr, gate) expr
    api_env$gate_histogram_autogate_ranges <- function(...) {
        list(positive = c(60, 100), negative = c(1, 40))
    }
    polygon <- list(mode = "scatter", vertices = list(
        list(x = 0, y = 0), list(x = 1, y = 0), list(x = 1, y = 1)
    ))
    gates <- list("cell:beads" = polygon, "singlet:beads" = polygon)

    generated <- api_env$gate_autogenerate_histograms(gates)

    testthat::expect_setequal(names(generated$gates), "positive:bead-scc.fcs")
    testthat::expect_equal(generated$gates_generated, 1L)
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
