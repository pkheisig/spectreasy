make_synthetic_ff <- function(primary_detector = c("B1-A" = 1200, "YG1-A" = 100), n = 800) {
    stopifnot(all(c("B1-A", "YG1-A") %in% names(primary_detector)))

    pos_n <- n %/% 2
    neg_n <- n - pos_n
    exprs <- cbind(
        "B1-A" = c(rlnorm(pos_n, log(primary_detector[["B1-A"]]), 0.15), rlnorm(neg_n, log(35), 0.18)),
        "YG1-A" = c(rlnorm(pos_n, log(primary_detector[["YG1-A"]]), 0.15), rlnorm(neg_n, log(25), 0.18)),
        "FSC-A" = rnorm(n, 100000, 4500),
        "SSC-A" = rnorm(n, 50000, 3500)
    )

    ff <- flowCore::flowFrame(exprs)
    pd <- flowCore::pData(flowCore::parameters(ff))
    pd$desc <- pd$name
    flowCore::parameters(ff) <- methods::new("AnnotatedDataFrame", data = pd)
    ff
}

make_synthetic_workflow <- function(include_af = FALSE) {
    set.seed(1)
    scc_dir <- tempfile("spectreasy_covr_scc_")
    dir.create(scc_dir, recursive = TRUE, showWarnings = FALSE)

    fitc_ff <- make_synthetic_ff(c("B1-A" = 1200, "YG1-A" = 120))
    pe_ff <- make_synthetic_ff(c("B1-A" = 150, "YG1-A" = 1300))

    flowCore::write.FCS(fitc_ff, file.path(scc_dir, "FITC (Beads).fcs"))
    flowCore::write.FCS(pe_ff, file.path(scc_dir, "PE (Beads).fcs"))

    control_df <- data.frame(
        filename = c("FITC (Beads).fcs", "PE (Beads).fcs"),
        fluorophore = c("FITC", "PE"),
        marker = c("CD1", "CD2"),
        channel = c("B1-A", "YG1-A"),
        control.type = c("beads", "beads"),
        universal.negative = c("", ""),
        large.gate = c("", ""),
        is.viability = c("", ""),
        stringsAsFactors = FALSE
    )

    if (isTRUE(include_af)) {
        af_ff <- make_synthetic_ff(c("B1-A" = 60, "YG1-A" = 55))
        flowCore::write.FCS(af_ff, file.path(scc_dir, "Unstained (Cells).fcs"))
        control_df <- rbind(
            control_df,
            data.frame(
                filename = "Unstained (Cells).fcs",
                fluorophore = "AF",
                marker = "Autofluorescence",
                channel = "B1-A",
                control.type = "cells",
                universal.negative = "",
                large.gate = "TRUE",
                is.viability = "",
                stringsAsFactors = FALSE
            )
        )
    }

    list(scc_dir = scc_dir, control_df = control_df)
}

test_that("AF helper detects AF-like control rows", {
    expect_true(spectreasy:::.is_af_control_row(fluorophore = "AF"))
    expect_true(spectreasy:::.is_af_control_row(fluorophore = "AF_Internal", marker = "Autofluorescence"))
    expect_true(spectreasy:::.is_af_control_row(marker = "Autofluorescence"))
    expect_true(spectreasy:::.is_af_control_row(filename = "Unstained (Cells).fcs"))
    expect_false(spectreasy:::.is_af_control_row(fluorophore = "FITC", filename = "FITC (Beads).fcs"))

    expect_true(spectreasy:::.is_primary_af_control_row(fluorophore = "AF"))
    expect_true(spectreasy:::.is_primary_af_control_row(fluorophore = "AF_2"))
    expect_false(spectreasy:::.is_primary_af_control_row(
        fluorophore = "AF_Internal",
        marker = "Autofluorescence",
        filename = "scc_cells_AF_UnstainedDead.fcs"
    ))
    expect_true(spectreasy:::.is_dead_af_control_row(
        fluorophore = "AF_dead",
        filename = "scc_cells_AF_UnstainedDead.fcs",
        control_type = "cells"
    ))
    expect_true(spectreasy:::.is_dead_af_control_row(
        fluorophore = "AF_Internal",
        marker = "Autofluorescence",
        filename = "scc_cells_AF_UnstainedDead.fcs",
        control_type = "cells"
    ))
    expect_true(spectreasy:::.is_primary_af_control_row(
        fluorophore = "",
        marker = "Autofluorescence",
        filename = "Unstained (Cells).fcs"
    ))
})

test_that("validate_control_file_mapping validates synthetic SCC setup", {
    wf <- make_synthetic_workflow()

    preflight_ok <- spectreasy::validate_control_file_mapping(
        control_df = wf$control_df,
        scc_dir = wf$scc_dir,
        require_channels = TRUE
    )
    expect_true(preflight_ok$ok)
    expect_length(preflight_ok$errors, 0)

    bad_df <- wf$control_df
    bad_df$channel[1] <- "NOT_A_CHANNEL"
    preflight_bad <- spectreasy::validate_control_file_mapping(
        control_df = bad_df,
        scc_dir = wf$scc_dir,
        require_channels = TRUE
    )
    expect_false(preflight_bad$ok)
    expect_true(any(grepl("Invalid channel", preflight_bad$errors)))
})

test_that("unmix preflight treats the control file as an inclusion list", {
    wf <- make_synthetic_workflow()
    extra_file <- file.path(wf$scc_dir, "Unmapped extra control.fcs")
    expect_true(file.copy(file.path(wf$scc_dir, wf$control_df$filename[1]), extra_file))

    preflight_extra <- spectreasy:::.run_unmix_preflight(wf$control_df, wf$scc_dir)
    expect_true(preflight_extra$ok)
    expect_true(any(grepl("Ignoring SCC files not listed in control file", preflight_extra$warnings)))

    file.remove(file.path(wf$scc_dir, wf$control_df$filename[1]))
    preflight_missing <- spectreasy:::.run_unmix_preflight(wf$control_df, wf$scc_dir)
    expect_false(preflight_missing$ok)
    expect_true(any(grepl("Control files missing from scc_dir", preflight_missing$errors)))
})

test_that("build_reference_matrix works on synthetic SCC files", {
    wf <- make_synthetic_workflow()

    M <- spectreasy::build_reference_matrix(
        input_folder = wf$scc_dir,
        control_df = wf$control_df,
        save_qc_plots = FALSE,
        seed = 1,
        subsample_n = 400
    )

    expect_equal(rownames(M), c("FITC", "PE"))
    expect_setequal(colnames(M), c("B1-A", "YG1-A"))
    expect_gt(M["FITC", "B1-A"], M["FITC", "YG1-A"])
    expect_gt(M["PE", "YG1-A"], M["PE", "B1-A"])

    qc_summary <- attr(M, "qc_summary")
    expect_true(is.data.frame(qc_summary))
    expect_equal(nrow(qc_summary), 2)
})

test_that("build_reference_matrix fails when a mapped SCC cannot produce a spectrum", {
    wf <- make_synthetic_workflow()
    bad_ff <- flowCore::flowFrame(cbind(
        "B1-A" = rnorm(50, mean = 100, sd = 5),
        "YG1-A" = rnorm(50, mean = 20, sd = 2),
        "FSC-A" = rnorm(50, mean = 90000, sd = 3000),
        "SSC-A" = rnorm(50, mean = 45000, sd = 2000)
    ))
    flowCore::write.FCS(bad_ff, file.path(wf$scc_dir, "APC (Beads).fcs"))
    wf$control_df <- rbind(
        wf$control_df,
        data.frame(
            filename = "APC (Beads).fcs",
            fluorophore = "APC",
            marker = "CD3",
            channel = "B1-A",
            control.type = "beads",
            universal.negative = "",
            large.gate = "",
            is.viability = "",
            stringsAsFactors = FALSE
        )
    )

    expect_error(
        spectreasy::build_reference_matrix(
            input_folder = wf$scc_dir,
            control_df = wf$control_df,
            save_qc_plots = FALSE,
            seed = 1,
            subsample_n = 400
        ),
        regexp = "did not produce spectra for all mapped non-AF controls"
    )
})

test_that("build_reference_matrix fails early on SCC detector mismatches", {
    wf <- make_synthetic_workflow()
    mismatch_ff <- flowCore::flowFrame(cbind(
        "B1-A" = rnorm(800, mean = 1200, sd = 50),
        "FSC-A" = rnorm(800, mean = 90000, sd = 3000),
        "SSC-A" = rnorm(800, mean = 45000, sd = 2000)
    ))
    flowCore::write.FCS(mismatch_ff, file.path(wf$scc_dir, "APC (Beads).fcs"))
    wf$control_df <- rbind(
        wf$control_df,
        data.frame(
            filename = "APC (Beads).fcs",
            fluorophore = "APC",
            marker = "CD3",
            channel = "B1-A",
            control.type = "beads",
            universal.negative = "",
            large.gate = "",
            is.viability = "",
            stringsAsFactors = FALSE
        )
    )

    expect_error(
        spectreasy::build_reference_matrix(
            input_folder = wf$scc_dir,
            control_df = wf$control_df,
            save_qc_plots = FALSE,
            seed = 1,
            subsample_n = 400
        ),
        regexp = "Detector set mismatch"
    )
})

test_that("unmix_controls does not overwrite an existing invalid control file", {
    wf <- make_synthetic_workflow()
    output_dir <- tempfile("spectreasy_no_overwrite_")
    control_csv <- tempfile(fileext = ".csv")
    bad_df <- wf$control_df
    bad_df$channel[1] <- "NOT_A_CHANNEL"
    utils::write.csv(bad_df, control_csv, row.names = FALSE, quote = TRUE)
    before <- readLines(control_csv, warn = FALSE)

    expect_error(
        spectreasy::unmix_controls(
            scc_dir = wf$scc_dir,
            control_file = control_csv,
            auto_create_mapping = TRUE,
            output_dir = output_dir,
            seed = 1,
            subsample_n = 400
        ),
        regexp = "preflight failed"
    )
    expect_equal(readLines(control_csv, warn = FALSE), before)
})

test_that("unmix_controls errors early for an explicit missing gating file", {
    wf <- make_synthetic_workflow()
    output_dir <- tempfile("spectreasy_missing_gate_")
    control_csv <- tempfile(fileext = ".csv")
    missing_gate <- tempfile(fileext = ".csv")
    utils::write.csv(wf$control_df, control_csv, row.names = FALSE, quote = TRUE)

    expect_error(
        spectreasy::unmix_controls(
            scc_dir = wf$scc_dir,
            control_file = control_csv,
            output_dir = output_dir,
            unmixing_method = "OLS",
            gating_mode = "REUSE",
            manual_gate_file = missing_gate,
            save_report = FALSE,
            seed = 1,
            subsample_n = 120
        ),
        regexp = "requires an existing manual_gate_file"
    )
})

test_that("automatic gating ignores an explicitly missing gate file", {
    wf <- make_synthetic_workflow()
    output_dir <- tempfile("spectreasy_automatic_gate_")
    control_csv <- tempfile(fileext = ".csv")
    missing_gate <- tempfile(fileext = ".csv")
    utils::write.csv(wf$control_df, control_csv, row.names = FALSE, quote = TRUE)

    expect_no_error(
        spectreasy::unmix_controls(
            scc_dir = wf$scc_dir,
            control_file = control_csv,
            output_dir = output_dir,
            unmixing_method = "OLS",
            gating_mode = "AUTOMATIC",
            manual_gate_file = missing_gate,
            save_report = FALSE,
            seed = 1,
            subsample_n = 120
        )
    )
})

test_that("interactive gating reuses an existing gate file when no GUI session is available", {
    wf <- make_synthetic_workflow()
    output_dir <- tempfile("spectreasy_interactive_reuse_")
    control_csv <- tempfile(fileext = ".csv")
    gate_csv <- tempfile(fileext = ".csv")
    utils::write.csv(wf$control_df, control_csv, row.names = FALSE, quote = TRUE)
    utils::write.csv(
        data.frame(
            gate_type = "setting",
            scope = "global",
            filename = "",
            x_channel = "histogram_transform",
            y_channel = "",
            plot_mode = "setting",
            vertex_index = 0,
            x = "asinh",
            y = "",
            stringsAsFactors = FALSE
        ),
        gate_csv,
        row.names = FALSE,
        quote = TRUE
    )

    result <- NULL
    expect_warning(
        result <- spectreasy::unmix_controls(
            scc_dir = wf$scc_dir,
            control_file = control_csv,
            output_dir = output_dir,
            unmixing_method = "OLS",
            gating_mode = "interactive",
            manual_gate_file = gate_csv,
            save_report = FALSE,
            seed = 1,
            subsample_n = 120
        ),
        "reusing the existing manual_gate_file"
    )
    expect_identical(result$gating_mode, "reuse")
    expect_identical(result$manual_gate_file, normalizePath(gate_csv, mustWork = TRUE))
})

test_that("common missing path errors are user friendly", {
    missing_sample_dir <- tempfile("missing_samples_")
    expect_error(
        spectreasy::unmix_samples(sample_dir = missing_sample_dir),
        regexp = "sample_dir not found:"
    )

    empty_sample_dir <- tempfile("empty_samples_")
    dir.create(empty_sample_dir, recursive = TRUE)
    expect_error(
        spectreasy::unmix_samples(sample_dir = empty_sample_dir, M = matrix(1, nrow = 1, dimnames = list("FITC", "B1-A"))),
        regexp = "No FCS files found in sample_dir:"
    )

    wf <- make_synthetic_workflow()
    missing_matrix <- tempfile(fileext = ".csv")
    expect_error(
        spectreasy::unmix_samples(sample_dir = wf$scc_dir, unmixing_matrix_file = missing_matrix),
        regexp = "unmixing_matrix_file not found:"
    )

    missing_scc_dir <- tempfile("missing_scc_")
    expect_error(
        spectreasy::build_reference_matrix(input_folder = missing_scc_dir),
        regexp = "input_folder not found:"
    )
    expect_error(
        spectreasy::create_control_file(input_folder = missing_scc_dir),
        regexp = "input_folder not found:"
    )
})

test_that("unmix_controls runs end-to-end on synthetic SCC files", {
    wf <- make_synthetic_workflow()
    extra_file <- file.path(wf$scc_dir, "Unmapped extra control.fcs")
    expect_true(file.copy(file.path(wf$scc_dir, wf$control_df$filename[1]), extra_file))
    output_dir <- tempfile("spectreasy_covr_auto_")
    control_csv <- tempfile(fileext = ".csv")
    utils::write.csv(wf$control_df, control_csv, row.names = FALSE, quote = TRUE)

    ctrl <- spectreasy::unmix_controls(
        scc_dir = wf$scc_dir,
        control_file = control_csv,
        output_dir = output_dir,
        unmixing_method = "OLS",
        n_threads = 2,
        save_qc_png = TRUE,
        seed = 1,
        subsample_n = 400
    )

    expect_true(file.exists(ctrl$reference_matrix_file))
    expect_true(file.exists(ctrl$unmixing_matrix_file))
    expect_true(file.exists(ctrl$qc_report_file))
    expect_true(dir.exists(ctrl$qc_controls_dir))
    expect_true(file.exists(ctrl$spectra_file))
    expect_true(file.exists(ctrl$unmixing_scatter_file))
    expect_s3_class(ctrl$unmixing_scatter_plot, "ggplot")
    expect_equal(sort(names(ctrl$unmixed_list)), c("FITC (Beads)", "PE (Beads)"))
    expect_false(file.exists(file.path(output_dir, "unmix_controls", "unmixed_fcs", "Unmapped extra control_OLS-0AF.fcs")))
})

test_that("unmix workflows route a shared output root into stage-specific unmixed_fcs folders", {
    wf <- make_synthetic_workflow()
    project_dir <- tempfile("spectreasy_loop_style_")
    sample_dir <- tempfile("spectreasy_loop_samples_")
    control_csv <- tempfile(fileext = ".csv")
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(wf$control_df, control_csv, row.names = FALSE, quote = TRUE)
    flowCore::write.FCS(make_synthetic_ff(c("B1-A" = 900, "YG1-A" = 150), n = 120), file.path(sample_dir, "sample.fcs"))

    ctrl <- spectreasy::unmix_controls(
        scc_dir = wf$scc_dir,
        control_file = control_csv,
        output_dir = project_dir,
        unmixing_method = "OLS",
        save_report = FALSE,
        seed = 1,
        subsample_n = 120
    )

    expect_true(file.exists(file.path(project_dir, "unmix_controls", "scc_reference_matrix.csv")))
    expect_equal(normalizePath(ctrl$reference_matrix_file), normalizePath(file.path(project_dir, "unmix_controls", "scc_reference_matrix.csv")))
    expect_true(file.exists(file.path(project_dir, "unmix_controls", "unmixed_fcs", "FITC (Beads)_OLS-0AF.fcs")))

    unmixed <- spectreasy::unmix_samples(
        sample_dir = sample_dir,
        unmixing_method = "OLS",
        output_dir = project_dir,
        write_fcs = TRUE,
        save_report = FALSE
    )

    expect_s3_class(unmixed, "spectreasy_unmixed_results")
    expect_equal(names(unmixed), "sample")
    expect_true(file.exists(file.path(project_dir, "unmix_samples", "unmixed_fcs", "sample_OLS-0AF.fcs")))
})

test_that("unmix_controls handles WLS output with AF controls", {
    wf <- make_synthetic_workflow(include_af = TRUE)
    output_dir <- tempfile("spectreasy_covr_auto_wls_")
    control_csv <- tempfile(fileext = ".csv")
    utils::write.csv(wf$control_df, control_csv, row.names = FALSE, quote = TRUE)

    messages <- capture.output(
        ctrl <- spectreasy::unmix_controls(
            scc_dir = wf$scc_dir,
            control_file = control_csv,
            output_dir = output_dir,
            unmixing_method = "WLS",
            af_n_bands = 1,
            seed = 1,
            subsample_n = 400
        ),
        type = "message"
    )

    expect_equal(ctrl$static_unmixing_matrix_method, "WLS")
    expect_true(file.exists(ctrl$detector_noise_file))
    expect_identical(
        colnames(utils::read.csv(ctrl$detector_noise_file, check.names = FALSE)),
        c("detector", "noise_floor")
    )
    expect_true(file.exists(ctrl$qc_report_file))
    expect_equal(normalizePath(dirname(ctrl$qc_report_file)), normalizePath(ctrl$qc_controls_dir))
    expect_true(all(file.exists(ctrl$qc_nxn_files)))
    expect_true(all(normalizePath(dirname(ctrl$qc_nxn_files)) == normalizePath(ctrl$qc_controls_dir)))
    expect_true(dir.exists(ctrl$qc_controls_dir))
    expect_false(is.null(attr(ctrl$M, "detector_noise")))
    expect_true(any(grepl("^AF($|_)", rownames(ctrl$M), ignore.case = TRUE)))
    expect_true(any(grepl("Unstained", names(ctrl$unmixed_list), ignore.case = TRUE)))
    expect_true(file.exists(file.path(output_dir, "unmix_controls", "unmixed_fcs", "Unstained (Cells)_WLS-1AF.fcs")))
    expect_false(any(grepl("scc_report_plots_", list.dirs(output_dir, recursive = TRUE, full.names = FALSE))))
    expect_equal(sum(grepl("^Detectors\\s*: 2 spectral channel", messages)), 1L)
    expect_equal(sum(grepl("^SCC\\s*:", messages)), 2L)
    expect_equal(sum(grepl("^AF bank\\s*: 1 signature", messages)), 1L)
})

test_that("unmix_controls tolerates a missing unstained mapping row", {
    wf <- make_synthetic_workflow(include_af = TRUE)
    file.remove(file.path(wf$scc_dir, "Unstained (Cells).fcs"))

    output_dir <- tempfile("spectreasy_missing_af_")
    control_csv <- tempfile(fileext = ".csv")
    utils::write.csv(wf$control_df, control_csv, row.names = FALSE, quote = TRUE)

    ctrl <- spectreasy::unmix_controls(
        scc_dir = wf$scc_dir,
        control_file = control_csv,
        output_dir = output_dir,
        unmixing_method = "WLS",
        seed = 1,
        subsample_n = 300
    )

    expect_true(file.exists(ctrl$reference_matrix_file))
    expect_true(file.exists(ctrl$qc_report_file))
    expect_true(dir.exists(ctrl$qc_controls_dir))
    expect_false(any(grepl("^AF($|_)", rownames(ctrl$M), ignore.case = TRUE)))
    expect_equal(sort(names(ctrl$unmixed_list)), c("FITC (Beads)", "PE (Beads)"))
})

test_that("unmix_samples runs WLS from a saved reference matrix without variance metadata", {
    wf <- make_synthetic_workflow()
    M <- spectreasy::build_reference_matrix(
        input_folder = wf$scc_dir,
        control_df = wf$control_df,
        save_qc_plots = FALSE,
        seed = 1,
        subsample_n = 400
    )

    ref_file <- tempfile("scc_reference_matrix_", fileext = ".csv")
    ref_df <- as.data.frame(M, check.names = FALSE)
    ref_df$Marker <- rownames(M)
    ref_df <- ref_df[, c("Marker", setdiff(colnames(ref_df), "Marker")), drop = FALSE]
    utils::write.csv(ref_df, ref_file, row.names = FALSE)

    sample_dir <- tempfile("spectreasy_recompute_samples_")
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    flowCore::write.FCS(make_synthetic_ff(c("B1-A" = 900, "YG1-A" = 500), n = 120), file.path(sample_dir, "sample.fcs"))

    output_dir <- tempfile("spectreasy_recompute_unmixed_")
    res <- spectreasy::unmix_samples(
        sample_dir = sample_dir,
        unmixing_matrix_file = ref_file,
        unmixing_method = "WLS",
        output_dir = output_dir,
        write_fcs = TRUE
    )

    expect_s3_class(res, "spectreasy_unmixed_results")
    expect_true(file.exists(file.path(output_dir, "unmix_samples", "unmixed_fcs", "sample_WLS-0AF.fcs")))
})

test_that("unmix_samples writes FCS files by default and returns invisibly", {
    wf <- make_synthetic_workflow()
    sample_dir <- tempfile("spectreasy_covr_samples_")
    output_dir <- tempfile("spectreasy_covr_unmixed_")
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)

    flowCore::write.FCS(make_synthetic_ff(c("B1-A" = 900, "YG1-A" = 150)), file.path(sample_dir, "sample_a.fcs"))
    flowCore::write.FCS(make_synthetic_ff(c("B1-A" = 180, "YG1-A" = 1100)), file.path(sample_dir, "sample_b.fcs"))

    M <- spectreasy::build_reference_matrix(
        input_folder = wf$scc_dir,
        control_df = wf$control_df,
        save_qc_plots = FALSE,
        seed = 1,
        subsample_n = 400
    )

    call_result <- withVisible(
        spectreasy::unmix_samples(
            sample_dir = sample_dir,
            M = M,
            unmixing_method = "OLS",
            output_dir = output_dir
        )
    )

    expect_false(call_result$visible)
    expect_setequal(names(call_result$value), c("sample_a", "sample_b"))
    expect_true(file.exists(file.path(output_dir, "unmix_samples", "unmixed_fcs", "sample_a_OLS-0AF.fcs")))
    expect_true(file.exists(file.path(output_dir, "unmix_samples", "unmixed_fcs", "sample_b_OLS-0AF.fcs")))
})

test_that("unmix_samples accepts samples_dir as an alias for sample_dir", {
    wf <- make_synthetic_workflow()
    sample_dir <- tempfile("spectreasy_alias_samples_")
    output_dir <- tempfile("spectreasy_alias_unmixed_")
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)

    flowCore::write.FCS(make_synthetic_ff(c("B1-A" = 900, "YG1-A" = 150)), file.path(sample_dir, "sample_alias.fcs"))

    M <- spectreasy::build_reference_matrix(
        input_folder = wf$scc_dir,
        control_df = wf$control_df,
        save_qc_plots = FALSE,
        seed = 1,
        subsample_n = 400
    )

    unmixed <- spectreasy::unmix_samples(
        samples_dir = sample_dir,
        M = M,
        unmixing_method = "OLS",
        output_dir = output_dir,
        save_report = FALSE
    )

    expect_s3_class(unmixed, "spectreasy_unmixed_results")
    expect_true(file.exists(file.path(output_dir, "unmix_samples", "unmixed_fcs", "sample_alias_OLS-0AF.fcs")))
})

test_that("qc_controls writes a PDF from synthetic SCC files", {
    wf <- make_synthetic_workflow()
    output_pdf <- tempfile(fileext = ".pdf")
    qc_plot_dir <- tempfile("spectreasy_covr_scc_report_")
    qc_metrics_dir <- tempfile("spectreasy_covr_scc_metrics_")
    control_csv <- tempfile(fileext = ".csv")
    utils::write.csv(wf$control_df, control_csv, row.names = FALSE, quote = TRUE)

    out <- spectreasy::qc_controls(
        scc_dir = wf$scc_dir,
        output_file = output_pdf,
        control_file = control_csv,
        unmixing_method = "NNLS",
        qc_plot_dir = qc_plot_dir,
        save_qc_pngs = TRUE,
        qc_metrics_dir = qc_metrics_dir,
        seed = 1,
        subsample_n = 400
    )

    expect_true(file.exists(output_pdf))
    expect_true(dir.exists(out$qc_plot_dir))
    expect_true(is.matrix(out$M))
    expect_true(is.data.frame(out$qc_summary))
    expect_equal(out$unmixing_method, "NNLS")
    expect_equal(out$qc_metrics_dir, qc_metrics_dir)
    expect_true(file.exists(file.path(qc_metrics_dir, "reference_spectra.csv")))
    expect_true(file.exists(file.path(qc_metrics_dir, "fluorophore_spectral_similarity.csv")))
    expect_true(file.exists(file.path(qc_metrics_dir, "rms_residual_per_detector.csv")))
    expect_true(file.exists(file.path(qc_metrics_dir, "detector_reconstruction_error.csv")))
    expect_false(file.exists(file.path(qc_metrics_dir, "overall_detector_reconstruction_error_per_sample.csv")))
    expect_false(file.exists(file.path(qc_metrics_dir, "negative_population_spread.csv")))
    expect_false(file.exists(file.path(qc_metrics_dir, "sample_qc_summary.csv")))
    expect_false(file.exists(file.path(qc_metrics_dir, "spectral_spread_matrix.csv")))
    expect_false(file.exists(file.path(qc_metrics_dir, "directional_spread_score.csv")))

    # Test default output file behavior
    default_html <- "spectreasy_outputs/unmix_controls/qc_controls_report.html"
    if (file.exists(default_html)) file.remove(default_html)
    
    out_default <- spectreasy::qc_controls(
        scc_dir = wf$scc_dir,
        control_file = control_csv,
        qc_plot_dir = qc_plot_dir,
        save_qc_pngs = FALSE,
        seed = 1,
        subsample_n = 400
    )
    expect_true(file.exists(default_html))
    file.remove(default_html)
})

test_that("qc_controls does not retain QC PNGs unless requested", {
    wf <- make_synthetic_workflow()
    output_pdf <- tempfile(fileext = ".pdf")
    qc_plot_dir <- tempfile("spectreasy_scc_report_no_retain_")
    control_csv <- tempfile(fileext = ".csv")
    utils::write.csv(wf$control_df, control_csv, row.names = FALSE, quote = TRUE)

    out <- spectreasy::qc_controls(
        scc_dir = wf$scc_dir,
        output_file = output_pdf,
        control_file = control_csv,
        qc_plot_dir = qc_plot_dir,
        seed = 1,
        subsample_n = 400
    )

    expect_true(file.exists(output_pdf))
    expect_null(out$qc_plot_dir)
    expect_false(dir.exists(qc_plot_dir))
})

test_that("adjust_matrix starts packaged GUI on localhost", {
    skip_on_os("windows")
    skip_if_not_installed("plumber")

    expect_error(
        spectreasy::adjust_matrix(matrix_dir = tempfile("spectreasy_missing_matrix_"), open_browser = FALSE),
        regexp = "matrix_dir not found:"
    )

    tmp_matrix_dir <- tempfile("spectreasy_gui_matrix_")
    dir.create(tmp_matrix_dir, recursive = TRUE, showWarnings = FALSE)
    port <- sample(18000:18999, 1)

    job <- parallel::mcparallel({
        spectreasy::adjust_matrix(
            matrix_dir = tmp_matrix_dir,
            open_browser = FALSE,
            dev_mode = FALSE,
            port = port,
            unmixing_method = "AutoSpectral"
        )
    })

    on.exit({
        tools::pskill(job$pid)
        Sys.sleep(0.5)
    }, add = TRUE)

    Sys.sleep(2)
    resp <- tryCatch(
        readLines(sprintf("http://127.0.0.1:%s/status", port), warn = FALSE),
        error = function(e) character()
    )

    expect_true(length(resp) > 0)
    expect_true(any(grepl("ok", resp, fixed = TRUE)))
    expect_true(any(grepl("AutoSpectral", resp, fixed = TRUE)))
    expect_true(any(grepl(basename(tmp_matrix_dir), resp, fixed = TRUE)))
    expect_true(any(grepl("project_selected", resp, fixed = TRUE)))
})
