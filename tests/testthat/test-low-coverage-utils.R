test_that("matrix utility helpers cover key coercion and lookup branches", {
    df_marker <- data.frame(
        Marker = c("FITC", "PE"),
        `B1-A` = c("1", "0.1"),
        `YG1-A` = c("0.2", "1"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    M1 <- spectreasy:::.as_reference_matrix(df_marker, "M")
    expect_equal(rownames(M1), c("FITC", "PE"))
    expect_true(is.numeric(M1))

    df_first_col <- data.frame(
        fluor = c("FITC", "PE"),
        a = c(1, 0.1),
        b = c(0.2, 1),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    M2 <- spectreasy:::.as_reference_matrix(df_first_col, "M")
    expect_equal(rownames(M2), c("FITC", "PE"))

    expect_error(
        spectreasy:::.as_reference_matrix(data.frame(Marker = c("A", "B")), "M"),
        regexp = "no detector columns"
    )
    expect_error(
        spectreasy:::.as_reference_matrix(data.frame(Marker = c("A", "B"), X = c("a", "b")), "M"),
        regexp = "non-numeric detector columns"
    )

    pd <- data.frame(
        name = c("FSC-A", "UV1-A", "V2-A", "B3-A", "YG4-A", "R5-A", "Time"),
        desc = c("FSC-A", "355 nm - 379/28", "405 nm - 450/50", "488 nm - 530/30", "561 nm - 586/15", "640 nm - 670/30", "Time"),
        stringsAsFactors = FALSE
    )
    alias_map <- spectreasy:::.build_channel_alias_map_from_pd(pd)
    expect_equal(unname(alias_map[["UV1"]]), "UV1-A")
    expect_equal(unname(alias_map[["V2"]]), "V2-A")
    expect_equal(unname(alias_map[["B3"]]), "B3-A")
    expect_equal(unname(alias_map[["YG4"]]), "YG4-A")
    expect_equal(unname(alias_map[["R5"]]), "R5-A")

    expect_identical(spectreasy:::.get_passthrough_parameter_names(c("B1-A", "Time", "FSC-A", "SSC-H"), detector_names = "B1-A"), c("Time", "FSC-A", "SSC-H"))
    expect_equal(spectreasy:::.get_result_metadata_columns(c("FITC", "File", "Time", "FSC-A")), c("File", "Time", "FSC-A"))

    primary <- spectreasy:::.get_primary_scatter_channels(c("SSC-W", "FSC-H", "FSC-A", "SSC-A"))
    expect_equal(primary$fsc, "FSC-A")
    expect_equal(primary$ssc, "SSC-A")

    primary_alt <- spectreasy:::.get_primary_scatter_channels(c("FS-H", "Side Scatter-A"))
    expect_equal(primary_alt$fsc, "FS-H")
    expect_equal(primary_alt$ssc, "Side Scatter-A")

    out <- data.frame(FITC = 1:3)
    full_data <- matrix(1:12, nrow = 3)
    colnames(full_data) <- c("B1-A", "Time", "FSC-A", "SSC-A")
    appended <- spectreasy:::.append_passthrough_parameters(out, full_data, detector_names = "B1-A")
    expect_setequal(colnames(appended), c("FITC", "Time", "FSC-A", "SSC-A"))

    expect_error(spectreasy:::.with_optional_seed(NA), regexp = "finite integer")
    draws1 <- local({ spectreasy:::.with_optional_seed(42); rnorm(3) })
    draws2 <- local({ spectreasy:::.with_optional_seed(42); rnorm(3) })
    expect_equal(draws1, draws2)
})

test_that("diagnostic plot helpers handle save and no-residual branches", {
    M <- matrix(c(
        1, 0.2,
        0.1, 1
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B1-A", "YG1-A")

    res_list <- list(
        data = data.frame(FITC = c(1, 2, 3), PE = c(0.2, 0.1, 0.3)),
        residuals = matrix(c(1, -2, 3, 2, -1, 1), ncol = 2, byrow = TRUE, dimnames = list(NULL, colnames(M)))
    )
    pd <- data.frame(
        name = c("B1-A", "YG1-A"),
        desc = c("488 nm - 530/30", "561 nm - 586/15"),
        stringsAsFactors = FALSE
    )

    det_png <- tempfile(fileext = ".png")
    p_det <- spectreasy::plot_detector_residuals(res_list, M, top_n = 2, output_file = det_png, pd = pd)
    expect_s3_class(p_det, "ggplot")
    expect_true(file.exists(det_png))

    expect_warning(
        p_null <- spectreasy::plot_detector_residuals(list(data = data.frame(x = 1), residuals = NULL), M),
        regexp = "No residuals"
    )
    expect_null(p_null)

    nps <- data.frame(File = c("A", "B"), Marker = c("FITC", "PE"), NPS = c(0.1, 0.2))
    nps_png <- tempfile(fileext = ".png")
    p_nps <- spectreasy::plot_nps(nps, output_file = nps_png)
    expect_s3_class(p_nps, "ggplot")
    expect_true(file.exists(nps_png))
})

test_that("launch_gui internal helpers validate packaged assets and dev-mode requirements", {
    paths <- spectreasy:::.prepare_launch_gui_paths()
    expect_true(file.exists(paths$api_path))
    expect_true(dir.exists(paths$gui_path))

    tmp_matrix_dir <- tempfile("spectreasy_gui_dirs_")
    dir.create(tmp_matrix_dir, recursive = TRUE, showWarnings = FALSE)
    dirs <- spectreasy:::.normalize_gui_dirs(tmp_matrix_dir)
    expect_equal(dirs$matrix_dir, normalizePath(tmp_matrix_dir))
    expect_match(dirs$samples_dir, "samples$")

    old_wd <- getwd()
    tmp_wd <- tempfile("spectreasy_gui_default_")
    dir.create(file.path(tmp_wd, "spectreasy_outputs", "unmix_controls"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(tmp_wd, "samples"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(tmp_wd, "spectreasy_outputs", "unmix_samples", "unmixed_fcs"), recursive = TRUE, showWarnings = FALSE)
    on.exit(setwd(old_wd), add = TRUE)
    setwd(tmp_wd)
    expect_equal(
        spectreasy:::.default_launch_gui_matrix_dir(),
        normalizePath(file.path(tmp_wd, "spectreasy_outputs", "unmix_controls"))
    )
    expect_equal(
        spectreasy:::.default_launch_gui_samples_dir(),
        normalizePath(file.path(tmp_wd, "samples"))
    )
    expect_equal(
        spectreasy:::.normalize_gui_dirs(spectreasy:::.default_launch_gui_matrix_dir())$samples_dir,
        normalizePath(file.path(tmp_wd, "samples"))
    )
    expect_equal(
        spectreasy:::.normalize_gui_dirs(file.path(tmp_wd, "spectreasy_outputs", "unmix_controls"))$samples_dir,
        normalizePath(file.path(tmp_wd, "samples"))
    )
    unlink(file.path(tmp_wd, "samples"), recursive = TRUE)
    expect_equal(
        spectreasy:::.default_launch_gui_samples_dir(
            matrix_dir = file.path(tmp_wd, "spectreasy_outputs", "unmix_controls")
        ),
        normalizePath(file.path(tmp_wd, "spectreasy_outputs", "unmix_samples", "unmixed_fcs"))
    )
    unlink(file.path(tmp_wd, "spectreasy_outputs", "unmix_samples", "unmixed_fcs"), recursive = TRUE)
    writeLines("", file.path(tmp_wd, "spectreasy_outputs", "unmix_samples", "sample_unmixed.fcs"))
    expect_equal(
        spectreasy:::.default_launch_gui_samples_dir(),
        normalizePath(file.path(tmp_wd, "spectreasy_outputs", "unmix_samples"))
    )

    frontend_bundled <- spectreasy:::.resolve_launch_gui_frontend(
        gui_path = paths$gui_path,
        dist_path = paths$dist_path,
        port = 8000,
        dev_mode = FALSE
    )
    expect_equal(frontend_bundled$mode, "bundled")
    expect_match(frontend_bundled$frontend_url, "127.0.0.1:8000")

    expect_error(
        spectreasy:::.resolve_launch_gui_frontend(
            gui_path = paths$gui_path,
            dist_path = tempfile("spectreasy_missing_dist_"),
            port = 8000,
            dev_mode = FALSE
        ),
        regexp = "Bundled GUI assets not found"
    )

    expect_error(
        spectreasy:::.resolve_launch_gui_frontend(
            gui_path = paths$gui_path,
            dist_path = paths$dist_path,
            port = 8000,
            dev_mode = TRUE,
            npm_bin = ""
        ),
        regexp = "requires npm"
    )

    tmp_gui_path <- tempfile("spectreasy_gui_dev_")
    dir.create(tmp_gui_path, recursive = TRUE, showWarnings = FALSE)
    expect_error(
        spectreasy:::.resolve_launch_gui_frontend(
            gui_path = tmp_gui_path,
            dist_path = paths$dist_path,
            port = 9000,
            dev_mode = TRUE,
            npm_bin = "npm"
        ),
        regexp = "requires GUI dependencies"
    )

    dir.create(file.path(tmp_gui_path, "node_modules"), recursive = TRUE, showWarnings = FALSE)
    frontend_dev <- spectreasy:::.resolve_launch_gui_frontend(
        gui_path = tmp_gui_path,
        dist_path = paths$dist_path,
        port = 9000,
        dev_mode = TRUE,
        npm_bin = "npm"
    )
    expect_equal(frontend_dev$mode, "dev")
    expect_equal(frontend_dev$frontend_url, "http://127.0.0.1:5174")
    expect_equal(frontend_dev$npm_bin, "npm")
})

test_that("launch_gui dev-server helper starts npm with API base and restores working directory", {
    tmp_gui_path <- tempfile("spectreasy_gui_dev_server_")
    dir.create(tmp_gui_path, recursive = TRUE, showWarnings = FALSE)

    fake_npm <- tempfile("fake_npm_")
    fake_npm_log <- tempfile("fake_npm_log_")
    writeLines(
        c(
            "#!/bin/sh",
            "{",
            "  printf 'pwd=%s\\n' \"$PWD\"",
            "  printf 'args=%s\\n' \"$*\"",
            "  printf 'api=%s\\n' \"$VITE_API_BASE\"",
            paste0("} > ", shQuote(fake_npm_log))
        ),
        fake_npm
    )
    Sys.chmod(fake_npm, "0755")

    old_wd <- getwd()
    on.exit({
        setwd(old_wd)
    }, add = TRUE)

    expect_null(spectreasy:::.start_launch_gui_dev_server(tmp_gui_path, port = 8123, npm_bin = fake_npm))
    expect_equal(getwd(), old_wd)
    expect_true(file.exists(fake_npm_log))

    log_lines <- readLines(fake_npm_log, warn = FALSE)
    expect_equal(normalizePath(sub("^pwd=", "", log_lines[1])), normalizePath(tmp_gui_path))
    expect_equal(log_lines[2], "args=run dev")
    expect_equal(log_lines[3], "api=http://127.0.0.1:8123")
})
