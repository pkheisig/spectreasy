test_that("sample HTML report follows PDF plot order and writes an NxN companion", {
    skip_if_not_installed("base64enc")
    set.seed(11)
    M <- rbind(
        FITC = c(1, 0.2, 0.05),
        PE = c(0.1, 1, 0.2),
        AF_1 = c(0.15, 0.12, 0.1)
    )
    colnames(M) <- c("B1-A", "YG1-A", "R1-A")
    results <- data.frame(
        File = rep(c("sample_a", "sample_b"), each = 24),
        FITC = stats::rnorm(48),
        PE = stats::rnorm(48),
        AF_1 = stats::runif(48)
    )
    path <- tempfile(fileext = ".html")
    png_dir <- tempfile("sample_html_pngs_")
    report <- qc_samples(
        results,
        M = M,
        output_file = path,
        report_format = "html",
        qc_plot_dir = png_dir,
        save_qc_pngs = FALSE,
        sample_nxn_rows_per_page = 2,
        sample_nxn_max_points = 20
    )
    html <- paste(readLines(report$output_file, warn = FALSE), collapse = "\n")
    expect_true(file.exists(report$output_file))
    expect_false(dir.exists(png_dir))
    expect_false(report$self_contained)
    expect_length(report$companion_files, 1)
    expect_true(file.exists(report$companion_files))
    expect_match(html, "data:image/png;base64,", fixed = TRUE)
    expect_false(grepl("<table", html, fixed = TRUE))
    expect_false(grepl("<figcaption", html, fixed = TRUE))
    expect_false(grepl("Input status", html, fixed = TRUE))
    expect_match(html, "<small>NxN rows per page</small><strong>2</strong>", fixed = TRUE)
    expect_match(html, "<small>NxN events per sample</small><strong>20</strong>", fixed = TRUE)
    expect_identical(names(spectreasy:::.report_sections(report$report_data)), c("spectra", "similarity", "nps", "nxn"))
    expect_match(html, "Reference spectra", fixed = TRUE)
    expect_match(html, "Fluorophore similarity", fixed = TRUE)
    expect_match(html, "Negative population spread", fixed = TRUE)
    expect_match(html, "Sample NxN scatter matrices", fixed = TRUE)
    expect_false(grepl("Input samples", html, fixed = TRUE))
    expect_false(grepl("Unmixing settings", html, fixed = TRUE))
    expect_false(grepl("Generated artifacts", html, fixed = TRUE))
    expect_false(grepl("Generated locally by Spectreasy", html, fixed = TRUE))
    expect_false(grepl("Report plots are embedded", html, fixed = TRUE))
    expect_match(html, "sidebar-toggle", fixed = TRUE)
    expect_match(html, "IntersectionObserver", fixed = TRUE)
    expect_match(html, "--accent:#f28b54", fixed = TRUE)
    expect_match(html, "class=\"report-sample\"", fixed = TRUE)
    expect_match(html, "grid-template-columns:156px", fixed = TRUE)
    expect_match(html, "id=\"report-width\"", fixed = TRUE)
    expect_match(html, "value=\"67\"", fixed = TRUE)
    expect_match(html, "--report-width", fixed = TRUE)
    nxn_html <- paste(readLines(report$companion_files, warn = FALSE), collapse = "\n")
    expect_match(nxn_html, "matrix-canvas", fixed = TRUE)
    expect_match(nxn_html, "id=\"sample-select\"", fixed = TRUE)
    expect_match(nxn_html, "id=\"gradient-select\"", fixed = TRUE)
    expect_match(nxn_html, "id=\"plot-type\"", fixed = TRUE)
    expect_match(nxn_html, "id=\"point-size\"", fixed = TRUE)
    expect_match(nxn_html, "id=\"cell-size\"", fixed = TRUE)
    expect_match(nxn_html, '"name":"sample_a"', fixed = TRUE)
    expect_equal(length(strsplit(nxn_html, "id=\"matrix-canvas\"", fixed = TRUE)[[1]]) - 1L, 1L)
    expect_false(grepl("<img", nxn_html, fixed = TRUE))
    expect_false(grepl("<section", nxn_html, fixed = TRUE))
    expect_match(html, "reportUrl.searchParams.set('path',companionPath)", fixed = TRUE)
})

test_that("HTML report metadata uses its explicit project path", {
    project <- tempfile("active_report_project_")
    launch_dir <- tempfile("report_launch_dir_")
    dir.create(project)
    dir.create(launch_dir)
    old_options <- options(spectreasy.project_dir = project)
    old_wd <- setwd(launch_dir)
    on.exit(options(old_options), add = TRUE)
    on.exit(setwd(old_wd), add = TRUE)

    M <- matrix(1, nrow = 1, dimnames = list("FITC" = "FITC", "B1-A" = "B1-A"))
    report <- collect_control_report_data(
        M,
        scc_dir = tempfile("missing_scc_"),
        control_file = tempfile(fileext = ".csv"),
        project_path = project
    )
    expect_identical(report$project_path, normalizePath(project, mustWork = TRUE))
    expect_false(identical(report$project_path, normalizePath(launch_dir, mustWork = TRUE)))
})

test_that("sample detector residual table includes every available detector", {
    M <- rbind(FITC = c(1, 0.2, 0.05), PE = c(0.1, 1, 0.2))
    colnames(M) <- c("B1-A", "YG1-A", "R1-A")
    sample_data <- data.frame(FITC = stats::rnorm(20), PE = stats::rnorm(20))
    residuals <- matrix(stats::rnorm(60), nrow = 20, dimnames = list(NULL, colnames(M)))
    results <- list(sample = list(data = sample_data, residuals = residuals))
    class(results) <- c("spectreasy_unmixed_results", "list")
    attr(results, "method") <- "OLS"
    report_data <- collect_sample_report_data(results, M, unmixing_method = "OLS", sample_nxn_max_points = 10)
    expect_equal(nrow(report_data$detector_rms), ncol(M))
    expect_setequal(report_data$detector_rms$detector, colnames(M))
})

test_that("sample report renders only PDF-equivalent plots and no raw metric tables", {
    M <- rbind(FITC = c(1, 0.2, 0.05), PE = c(0.1, 1, 0.2))
    colnames(M) <- c("B1-A", "YG1-A", "R1-A")
    sample_data <- data.frame(FITC = stats::rnorm(20), PE = stats::rnorm(20))
    residuals <- matrix(stats::rnorm(60), nrow = 20, dimnames = list(NULL, colnames(M)))
    results <- list(sample = list(data = sample_data, residuals = residuals))
    class(results) <- c("spectreasy_unmixed_results", "list")
    attr(results, "method") <- "OLS"

    report_data <- collect_sample_report_data(
        results,
        M,
        unmixing_method = "OLS",
        sample_nxn_max_points = 10,
        run_settings = list(seed = 1, n_threads = 4L, rwls_max_iter = 1L, chunk_size = 50000L)
    )
    sections <- spectreasy:::.report_sections(report_data)
    changed_settings <- spectreasy:::.report_changed_run_settings(report_data)

    expect_false(grepl("<table", paste(vapply(sections, `[[`, character(1), 2), collapse = ""), fixed = TRUE))
    expect_match(sections[["spectra"]][2], "data:image/png;base64,", fixed = TRUE)
    expect_match(sections[["similarity"]][2], "data:image/png;base64,", fixed = TRUE)
    expect_false(grepl("<table", sections[["residual"]][2], fixed = TRUE))
    expect_match(sections[["residual"]][2], "data:image/png;base64,", fixed = TRUE)
    expect_match(sections[["nps"]][2], "data:image/png;base64,", fixed = TRUE)
    expect_match(sections[["reconstruction"]][2], "data:image/png;base64,", fixed = TRUE)
    expect_identical(names(sections), c("spectra", "similarity", "nps", "nxn", "residual", "reconstruction"))
    expect_false(grepl("<figcaption", paste(vapply(sections, `[[`, character(1), 2), collapse = ""), fixed = TRUE))
    expect_identical(names(changed_settings), c("Threads", "Seed"))
    expect_identical(unname(unlist(changed_settings)), c("4", "1"))
})

test_that("missing metric plots do not fall back to raw numeric tables", {
    metric <- matrix(c(1, 0.4, 0.4, 1), nrow = 2)
    colnames(metric) <- rownames(metric) <- c("FITC", "PE")
    html <- spectreasy:::.report_metric_html(metric, tempfile(fileext = ".png"))
    expect_false(grepl("<table", html, fixed = TRUE))
    expect_false(grepl("0.4", html, fixed = TRUE))
    expect_match(html, "plot was not available", fixed = TRUE)
})

test_that("partial NxN pages use their real row count instead of padded blank panels", {
    markers <- paste0("M", seq_len(17))
    results <- as.data.frame(matrix(stats::rnorm(60 * length(markers)), nrow = 60))
    colnames(results) <- markers
    results$File <- "sample"
    pages <- spectreasy:::.build_qc_report_sample_scatter_pages(
        results,
        markers = markers,
        rows_per_page = 8,
        max_points = 50,
        all_samples = TRUE
    )

    expect_length(pages, 5)
    expect_identical(attr(pages[[5]], "spectreasy_nxn_rows"), 1L)
    expect_identical(attr(pages[[5]], "spectreasy_nxn_cols"), 8L)
    expect_equal(unname(spectreasy:::.report_nxn_plot_size(pages[[5]])), c(10, 2.5))
})

test_that("control HTML report uses supplied matrix without rebuilding controls", {
    skip_if_not_installed("base64enc")
    M <- rbind(FITC = c(1, 0.1), PE = c(0.2, 1), AF_1 = c(0.1, 0.1))
    colnames(M) <- c("B1-A", "YG1-A")
    missing_scc <- tempfile("missing_scc_")
    path <- tempfile(fileext = ".html")
    png_dir <- tempfile("control_html_pngs_")
    report <- qc_controls(
        M = M,
        scc_dir = missing_scc,
        control_file = tempfile(fileext = ".csv"),
        output_file = path,
        report_format = "html",
        unmixing_method = "OLS",
        qc_plot_dir = png_dir,
        save_qc_pngs = FALSE
    )
    html <- paste(readLines(report$output_file, warn = FALSE), collapse = "\n")
    expect_true(file.exists(report$output_file))
    expect_false(dir.exists(png_dir))
    expect_false(grepl("<table", html, fixed = TRUE))
    expect_match(html, "Reference spectra", fixed = TRUE)
    expect_match(html, "Fluorophore similarity", fixed = TRUE)
    expect_false(grepl("AF bank spectra", html, fixed = TRUE))
    expect_false(grepl("SCC AF spectra", html, fixed = TRUE))
    expect_false(grepl("SCC unmixing diagnostics", html, fixed = TRUE))
    expect_false(grepl("Generated locally by Spectreasy", html, fixed = TRUE))
    expect_false(grepl("Report plots are embedded", html, fixed = TRUE))
    expect_match(html, "class=\"report-control\"", fixed = TRUE)
    expect_match(html, "--accent:#70afe8", fixed = TRUE)
})

test_that("control HTML groups each PDF control page once without duplicate plot sections", {
    skip_if_not_installed("base64enc")
    plot_root <- tempfile("control_pdf_panels_")
    for (directory in c("fsc_ssc", "histogram", "spectrum")) {
        dir.create(file.path(plot_root, directory), recursive = TRUE, showWarnings = FALSE)
    }
    p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
    ggplot2::ggsave(file.path(plot_root, "fsc_ssc", "ctrl_a_fsc_ssc.png"), p, width = 2, height = 2, dpi = 72)
    ggplot2::ggsave(file.path(plot_root, "histogram", "ctrl_a_histogram.png"), p, width = 2, height = 2, dpi = 72)
    ggplot2::ggsave(file.path(plot_root, "spectrum", "ctrl_a_spectrum.png"), p, width = 2, height = 2, dpi = 72)

    M <- rbind(FITC = c(1, 0.2), PE = c(0.1, 1), AF_1 = c(0.2, 0.1), AF_2 = c(0.1, 0.2))
    colnames(M) <- c("B1-A", "YG1-A")
    report_data <- collect_control_report_data(
        M,
        scc_dir = tempfile("missing_scc_"),
        control_file = tempfile(fileext = ".csv"),
        qc_summary = data.frame(
            sample = "ctrl_a",
            fluorophore = "FITC",
            marker = "CD3",
            type = "beads",
            peak_channel = "B1-A",
            n_total = 100,
            n_scatter_gated = 80,
            scatter_gate_pct = 80,
            n_final = 50,
            histogram_gate_pct = 62.5
        ),
        report_plot_dir = plot_root,
        af_bank_info = list(source_count = 1, derived_bands = 2),
        plots = list(unmixing_scatter = p),
        plot_dir = tempfile("control_report_plots_")
    )
    sections <- spectreasy:::.report_sections(report_data)
    expect_identical(names(sections), c("spectra", "similarity", "nxn", "gating"))
    expect_equal(length(strsplit(sections[["spectra"]][2], "<img", fixed = TRUE)[[1]]) - 1L, 2L)
    expect_length(report_data$plot_manifest$controls, 1)
    expect_length(report_data$plot_manifest$controls[[1]]$paths, 3)
    expect_identical(report_data$plot_manifest$controls[[1]]$image_titles, c("FSC/SSC gate", "Peak-channel histogram gate", "Per-event spectrum distribution"))
    expect_length(unique(report_data$plot_manifest$controls[[1]]$paths), 3)
    expect_false(grepl("<table", paste(vapply(sections, `[[`, character(1), 2), collapse = ""), fixed = TRUE))
    expect_false(grepl("SCC AF spectra", paste(vapply(sections, `[[`, character(1), 2), collapse = ""), fixed = TRUE))

    embedded_report_data <- spectreasy:::.report_embed_plot_manifest(report_data)
    expect_identical(embedded_report_data$plot_manifest$nxn, report_data$plot_manifest$nxn)
    expect_match(embedded_report_data$plot_manifest$reference, "data:image/png;base64,", fixed = TRUE)
    expected_default_nxn <- spectreasy:::.report_nxn_companion_paths("qc_controls_report.html", embedded_report_data)
    expect_identical(basename(expected_default_nxn), "qc_controls_report_nxn.png")
    rendered <- render_qc_html_report(embedded_report_data, tempfile(fileext = ".html"), overwrite = "overwrite")
    expect_length(rendered$companion_files, 1)
    expect_named(rendered$companion_files, "Controls")
    expect_identical(tolower(tools::file_ext(rendered$companion_files)), "png")
    expect_true(file.exists(rendered$companion_files))
    expect_gt(file.info(rendered$companion_files)$size, 0)
    expect_identical(rendered$report_data$plot_manifest$nxn, rendered$companion_files)
    main_html <- paste(readLines(rendered$output_file, warn = FALSE), collapse = "\n")
    expect_match(main_html, "Open high-resolution control NxN matrix", fixed = TRUE)
    expect_match(main_html, basename(rendered$companion_files), fixed = TRUE)
    expect_false(grepl("qc_controls_report_nxn.html", main_html, fixed = TRUE))
    expect_match(main_html, "<section id=\"gating\"><h2><button type=\"button\">Gating</button>", fixed = TRUE)
    nxn_section <- regmatches(main_html, regexpr("<section id=\"nxn\".*?</section>", main_html, perl = TRUE))
    expect_length(nxn_section, 1)
    expect_false(grepl("<img", nxn_section, fixed = TRUE))
})

test_that("control HTML NxN points remain visible in the full high-resolution matrix", {
    set.seed(19)
    unmixed <- list(
        fitc = data.frame(FITC = stats::rnorm(80, 10, 2), PE = stats::rnorm(80)),
        pe = data.frame(FITC = stats::rnorm(80), PE = stats::rnorm(80, 10, 2))
    )
    pages <- spectreasy:::.build_qc_report_control_scatter_pages(
        unmixed,
        sample_to_marker = c(fitc = "FITC", pe = "PE"),
        markers = c("FITC", "PE"),
        rows_per_page = 2,
        point_size = 0.35,
        point_alpha = 0.82
    )
    expect_length(pages, 1)
    point_layers <- Filter(function(layer) inherits(layer$geom, "GeomPoint"), pages[[1]]$layers)
    expect_length(point_layers, 1)
    expect_equal(point_layers[[1]]$aes_params$size, 0.35)
    expect_equal(point_layers[[1]]$aes_params$alpha, 0.82)
})

test_that("HTML writes one interactive viewer with a selector for all selected sample matrices", {
    markers <- paste0("M", seq_len(6))
    M <- matrix(stats::runif(length(markers) * 4), nrow = length(markers), dimnames = list(markers, paste0("D", 1:4)))
    results <- as.data.frame(matrix(stats::rnorm(80 * length(markers)), nrow = 80))
    colnames(results) <- markers
    results$File <- rep(c("sample_a", "sample_b"), each = 40)
    report <- qc_samples(
        results,
        M = M,
        output_file = tempfile(fileext = ".html"),
        report_format = "html",
        nxn_all_samples = TRUE,
        sample_nxn_rows_per_page = 2,
        sample_nxn_max_points = 30
    )
    expect_length(report$report_data$plot_manifest$nxn, 2)
    expect_setequal(names(report$report_data$plot_manifest$nxn), c("sample_a", "sample_b"))
    expect_length(report$companion_files, 1)
    expect_identical(names(report$companion_files), "Samples")
    expect_true(all(file.exists(report$companion_files)))
    companion <- paste(readLines(report$companion_files, warn = FALSE), collapse = "\n")
    expect_equal(length(strsplit(companion, "id=\"matrix-canvas\"", fixed = TRUE)[[1]]) - 1L, 1L)
    expect_match(companion, '"name":"sample_a"', fixed = TRUE)
    expect_match(companion, '"name":"sample_b"', fixed = TRUE)
    expect_match(companion, "Density", fixed = TRUE)
    expect_match(companion, "Hexbin", fixed = TRUE)
    expect_false(grepl("<img", companion, fixed = TRUE))
    expect_false(grepl("<section", companion, fixed = TRUE))
    main_html <- paste(readLines(report$output_file, warn = FALSE), collapse = "\n")
    expect_equal(length(strsplit(main_html, "class=\"companion-link\"", fixed = TRUE)[[1]]) - 1L, 1L)
})

test_that("per-sample HTML reporting uses one report with a sample selector", {
    markers <- c("FITC", "PE")
    M <- matrix(c(1, 0.2, 0.1, 1), nrow = 2, dimnames = list(markers, c("B1-A", "YG1-A")))
    results <- data.frame(
        File = rep(c("sample_a", "sample_b"), each = 12),
        FITC = stats::rnorm(24),
        PE = stats::rnorm(24)
    )
    report <- qc_samples(
        results,
        M = M,
        output_file = tempfile(fileext = ".html"),
        report_per_sample = TRUE,
        sample_nxn_max_points = 10,
        overwrite = "overwrite"
    )
    html <- paste(readLines(report$output_file, warn = FALSE), collapse = "\n")
    expect_match(html, "id=\"report-sample\"", fixed = TRUE)
    expect_match(html, ">sample_a</option>", fixed = TRUE)
    expect_match(html, ">sample_b</option>", fixed = TRUE)
    expect_equal(length(strsplit(html, "data-sample-panel=", fixed = TRUE)[[1]]) - 1L, 2L)
    expect_match(html, "sample-1-spectra", fixed = TRUE)
    expect_match(html, "sample-2-spectra", fixed = TRUE)
    expect_false(formals(unmix_samples)$report_per_sample)
    expect_false(formals(qc_samples)$report_per_sample)
})

test_that("per-sample PDF reporting writes one named file per sample", {
    markers <- c("FITC", "PE")
    M <- matrix(c(1, 0.2, 0.1, 1), nrow = 2, dimnames = list(markers, c("B1-A", "YG1-A")))
    results <- data.frame(
        File = rep(c("sample_a", "sample_b"), each = 8),
        FITC = stats::rnorm(16),
        PE = stats::rnorm(16)
    )
    report <- qc_samples(
        results,
        M = M,
        output_file = file.path(tempdir(), "per_sample_qc.pdf"),
        report_format = "pdf",
        report_per_sample = TRUE,
        sample_nxn_max_points = 8
    )
    expect_identical(names(report$output_file), c("sample_a", "sample_b"))
    expect_true(all(file.exists(report$output_file)))
    expect_match(basename(report$output_file[[1]]), "per_sample_qc_sample_a.pdf", fixed = TRUE)
    expect_match(basename(report$output_file[[2]]), "per_sample_qc_sample_b.pdf", fixed = TRUE)
})

test_that("interactive GUI gating remains the default and its four PDF panels are reused", {
    expect_identical(formals(unmix_controls)$gating_mode, "interactive")
    expect_false("manual_gating" %in% names(formals(unmix_controls)))

    plot_root <- tempfile("manual_gate_panels_")
    for (directory in c("fsc_ssc", "singlet", "histogram", "spectrum")) {
        dir.create(file.path(plot_root, directory), recursive = TRUE, showWarnings = FALSE)
    }
    p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
    ggplot2::ggsave(file.path(plot_root, "fsc_ssc", "ctrl_a_fsc_ssc.png"), p, width = 2, height = 2, dpi = 72)
    ggplot2::ggsave(file.path(plot_root, "singlet", "ctrl_a_singlet.png"), p, width = 2, height = 2, dpi = 72)
    ggplot2::ggsave(file.path(plot_root, "histogram", "ctrl_a_histogram.png"), p, width = 2, height = 2, dpi = 72)
    ggplot2::ggsave(file.path(plot_root, "spectrum", "ctrl_a_spectrum.png"), p, width = 2, height = 2, dpi = 72)

    panels <- spectreasy:::.report_control_panels(
        plot_root,
        qc_summary = data.frame(sample = "ctrl_a", fluorophore = "FITC", stringsAsFactors = FALSE)
    )
    expect_length(panels, 1)
    expect_identical(panels[[1]]$image_titles, c("Cell gate", "Singlet gate", "Histogram", "Per-event spectrum distribution"))
    expect_length(panels[[1]]$paths, 4)
    panel_html <- spectreasy:::.report_control_panels_html(panels)
    title_positions <- vapply(panels[[1]]$image_titles, function(title) regexpr(title, panel_html, fixed = TRUE)[1], integer(1))
    expect_true(all(title_positions > 0L))
    expect_true(all(diff(title_positions) > 0L))
    expect_equal(length(strsplit(panel_html, "class=\"plot-wrap\"", fixed = TRUE)[[1]]) - 1L, 4L)
})

test_that("report format and filename extensions are validated", {
    M <- rbind(FITC = c(1, 0.1), PE = c(0.1, 1))
    colnames(M) <- c("B1-A", "YG1-A")
    results <- data.frame(File = rep("sample", 8), FITC = stats::rnorm(8), PE = stats::rnorm(8))
    expect_error(qc_samples(results, M = M, output_file = tempfile(fileext = ".pdf"), report_format = "html"), "conflicts")
    expect_error(qc_samples(results, M = M, output_file = tempfile(fileext = ".html"), report_format = "pdf"), "conflicts")
    expect_error(qc_samples(results, M = M, report_format = "both"), "report_format")
    inferred <- qc_samples(results, M = M, output_file = tempfile(fileext = ".html"), sample_nxn_max_points = 8)
    expect_identical(inferred$format, "html")
})

test_that("report output paths must be non-empty scalar strings", {
    invalid_paths <- list(NULL, character(), NA_character_, "", "   ", c("one.html", "two.html"))

    for (path in invalid_paths) {
        expect_error(
            spectreasy:::.report_output_spec(path, report_format = "html"),
            "output_file must be a single non-empty file path",
            fixed = TRUE
        )
    }
})

test_that("HTML report wrappers fail clearly for missing scientific inputs", {
    M <- rbind(FITC = c(1, 0.1), PE = c(0.1, 1))
    colnames(M) <- c("B1-A", "YG1-A")
    results <- data.frame(File = "sample", FITC = 1, PE = 1)
    expect_error(collect_sample_report_data(NULL, M), "No unmixed results")
    expect_error(collect_sample_report_data(results, NULL), "No reference matrix")
    expect_error(
        qc_samples(results, unmixing_matrix_file = tempfile("missing_matrix_"), report_format = "html"),
        "not found"
    )
})

test_that("HTML overwrite policy versions, overwrites, errors, and detects staleness", {
    expect_identical(eval(formals(qc_controls)$overwrite), c("version", "overwrite", "error"))
    expect_identical(eval(formals(qc_samples)$overwrite), c("version", "overwrite", "error"))
    existing <- tempfile(fileext = ".html")
    writeLines("existing", existing)
    expect_match(basename(spectreasy:::.report_resolve_overwrite(existing)), "_v001\\.html$")

    M <- rbind(FITC = c(1, 0.1), PE = c(0.1, 1))
    colnames(M) <- c("B1-A", "YG1-A")
    results <- data.frame(File = rep("sample", 8), FITC = stats::rnorm(8), PE = stats::rnorm(8))
    target <- tempfile(fileext = ".html")
    first <- qc_samples(results, M = M, output_file = target, report_format = "html", overwrite = "overwrite", sample_nxn_max_points = 8)
    second <- qc_samples(results, M = M, output_file = target, report_format = "html", overwrite = "version", sample_nxn_max_points = 8)
    expect_false(identical(first$output_file, second$output_file))
    expect_match(basename(second$output_file), "_v001\\.html$")
    expect_error(qc_samples(results, M = M, output_file = target, report_format = "html", overwrite = "error", sample_nxn_max_points = 8), "already exists")

    source <- tempfile(fileext = ".csv")
    writeLines("x", source)
    Sys.sleep(1)
    report_file <- tempfile(fileext = ".html")
    writeLines("report", report_file)
    expect_false(spectreasy:::.report_is_stale(report_file, source))
    Sys.sleep(1)
    writeLines("changed", source)
    expect_true(spectreasy:::.report_is_stale(report_file, source))
})
test_that("AF bank report plots use the cockpit light accent", {
    af <- matrix(
        c(1, 0.5, 0.1, 1, 0.4, 0.08),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(c("AF", "AF_2"), c("UV7-A", "UV8-A", "V1-A"))
    )
    plot <- spectreasy:::.style_af_bank_plot(
        spectreasy::plot_spectra(af, output_file = NULL),
        rownames(af)
    )
    layer <- ggplot2::ggplot_build(plot)$data[[1]]

    expect_identical(unique(toupper(layer$colour)), "#1F7A6D")
    expect_equal(unique(layer$alpha), 0.2)
})
