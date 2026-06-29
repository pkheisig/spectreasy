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
                is.viability = "",
                stringsAsFactors = FALSE
            )
        )
    }

    list(scc_dir = scc_dir, control_df = control_df)
}

test_that("AF helper detects AF-like control rows", {
    expect_true(spectreasy:::.is_af_control_row(fluorophore = "AF"))
    expect_true(spectreasy:::.is_af_control_row(marker = "Autofluorescence"))
    expect_true(spectreasy:::.is_af_control_row(filename = "Unstained (Cells).fcs"))
    expect_false(spectreasy:::.is_af_control_row(fluorophore = "FITC", filename = "FITC (Beads).fcs"))
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

test_that("validate_control_file_mapping reports unreadable SCC files", {
    scc_dir <- tempfile("spectreasy_bad_scc_")
    dir.create(scc_dir, recursive = TRUE, showWarnings = FALSE)
    writeBin(as.raw(rep(0, 128)), file.path(scc_dir, "BV510 (Cells).fcs"))

    control_df <- data.frame(
        filename = "BV510 (Cells).fcs",
        fluorophore = "BV510",
        marker = "",
        channel = "",
        control.type = "cells",
        universal.negative = "",
        is.viability = "",
        stringsAsFactors = FALSE
    )

    preflight <- spectreasy::validate_control_file_mapping(
        control_df = control_df,
        scc_dir = scc_dir,
        require_channels = TRUE
    )

    expect_false(preflight$ok)
    expect_true(any(grepl("Missing channel", preflight$errors)))
    expect_true(any(grepl("Could not read FCS file for validation", preflight$errors)))
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

test_that("build_reference_matrix pools duplicate mapped unstained SCC files", {
    wf <- make_synthetic_workflow(include_af = TRUE)
    duplicate_name <- "Unstained duplicate (Cells).fcs"
    duplicate_path <- file.path(wf$scc_dir, duplicate_name)
    control_csv <- file.path(wf$scc_dir, "fcs_mapping.csv")
    testthat::expect_true(file.copy(file.path(wf$scc_dir, "Unstained (Cells).fcs"), duplicate_path))
    on.exit(unlink(c(duplicate_path, control_csv), force = TRUE), add = TRUE)

    control_df <- rbind(
        wf$control_df,
        data.frame(
            filename = duplicate_name,
            fluorophore = "AF_2",
            marker = "Autofluorescence",
            channel = "B1-A",
            control.type = "cells",
            universal.negative = "",
            is.viability = "",
            stringsAsFactors = FALSE
        )
    )
    utils::write.csv(control_df, control_csv, row.names = FALSE, quote = TRUE)

    M <- spectreasy::build_reference_matrix(
        input_folder = wf$scc_dir,
        control_df = control_csv,
        save_qc_plots = FALSE,
        af_n_bands = 2,
        seed = 1,
        subsample_n = 400
    )

    af_info <- attr(M, "af_bank_info")
    expect_equal(af_info$source_count, 2L)
    expect_setequal(af_info$sources$file, c("Unstained (Cells).fcs", duplicate_name))
    expect_true(all(c("AF", "AF_2") %in% rownames(M)))
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
            auto_create_control = TRUE,
            output_dir = output_dir,
            seed = 1,
            subsample_n = 400
        ),
        regexp = "preflight failed"
    )
    expect_equal(readLines(control_csv, warn = FALSE), before)
})

test_that("unmix_controls runs end-to-end on synthetic SCC files", {
    skip_slow_tests("full SCC control workflow with QC report plots")

    wf <- make_synthetic_workflow()
    output_dir <- tempfile("spectreasy_covr_auto_")
    control_csv <- tempfile(fileext = ".csv")
    utils::write.csv(wf$control_df, control_csv, row.names = FALSE, quote = TRUE)

    ctrl <- spectreasy::unmix_controls(
        scc_dir = wf$scc_dir,
        control_file = control_csv,
        output_dir = output_dir,
        unmix_method = "OLS",
        save_qc_plots = TRUE,
        seed = 1,
        subsample_n = 400
    )

    expect_true(file.exists(ctrl$reference_matrix_file))
    expect_true(file.exists(ctrl$unmixing_matrix_file))
    expect_true(file.exists(ctrl$spectra_file))
    expect_true(file.exists(ctrl$unmixing_scatter_file))
    expect_true(file.exists(ctrl$qc_report_file))
    expect_s3_class(ctrl$unmixing_scatter_plot, "ggplot")
    expect_equal(sort(names(ctrl$unmixed_list)), c("FITC (Beads)", "PE (Beads)"))

    regenerated_pdf <- tempfile(fileext = ".pdf")
    expect_no_error(
        suppressMessages(
            spectreasy::qc_controls(
                results = ctrl,
                scc_dir = wf$scc_dir,
                output_file = regenerated_pdf
            )
        )
    )
    expect_true(file.exists(regenerated_pdf))

    output_dir_pdf <- tempfile(fileext = ".pdf")
    expect_no_error(
        suppressMessages(
            spectreasy::qc_controls(
                unmix_controls_dir = output_dir,
                scc_dir = tempfile("missing_scc_"),
                output_file = output_dir_pdf
            )
        )
    )
    expect_true(file.exists(output_dir_pdf))
})

test_that("unmix_controls handles WLS output and exclude_af branch", {
    wf <- make_synthetic_workflow(include_af = TRUE)
    output_dir <- tempfile("spectreasy_covr_auto_wls_")
    control_csv <- tempfile(fileext = ".csv")
    utils::write.csv(wf$control_df, control_csv, row.names = FALSE, quote = TRUE)

    ctrl <- spectreasy::unmix_controls(
        scc_dir = wf$scc_dir,
        control_file = control_csv,
        output_dir = output_dir,
        exclude_af = TRUE,
        unmix_method = "WLS",
        save_report = FALSE,
        seed = 1,
        subsample_n = 400
    )

    expect_equal(ctrl$static_unmixing_matrix_method, "WLS")
    expect_true(file.exists(ctrl$detector_noise_file))
    expect_false(is.null(attr(ctrl$M, "detector_noise")))
    expect_false(any(grepl("^AF($|_)", rownames(ctrl$M), ignore.case = TRUE)))
    expect_false(any(grepl("Unstained", names(ctrl$unmixed_list), ignore.case = TRUE)))
    expect_false(file.exists(file.path(output_dir, "unmixed_fcs", "Unstained (Cells)_unmixed.fcs")))
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
        unmix_method = "WLS",
        save_report = FALSE,
        seed = 1,
        subsample_n = 300
    )

    expect_true(file.exists(ctrl$reference_matrix_file))
    expect_false(any(grepl("^AF($|_)", rownames(ctrl$M), ignore.case = TRUE)))
    expect_equal(sort(names(ctrl$unmixed_list)), c("FITC (Beads)", "PE (Beads)"))
})

test_that("unmix_samples runs WLS without recomputing missing SCC variances", {
    wf <- make_synthetic_workflow()
    M <- spectreasy::build_reference_matrix(
        input_folder = wf$scc_dir,
        control_df = wf$control_df,
        save_qc_plots = FALSE,
        seed = 1,
        subsample_n = 400
    )

    ref_file <- tempfile("scc_reference_matrix_", fileext = ".csv")
    var_file <- tempfile("scc_variances_", fileext = ".csv")
    ref_df <- as.data.frame(M, check.names = FALSE)
    ref_df$Marker <- rownames(M)
    ref_df <- ref_df[, c("Marker", setdiff(colnames(ref_df), "Marker")), drop = FALSE]
    utils::write.csv(ref_df, ref_file, row.names = FALSE)
    if (file.exists(var_file)) file.remove(var_file)

    sample_dir <- tempfile("spectreasy_recompute_samples_")
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    flowCore::write.FCS(make_synthetic_ff(c("B1-A" = 900, "YG1-A" = 500), n = 120), file.path(sample_dir, "sample.fcs"))

    output_dir <- tempfile("spectreasy_recompute_unmixed_")
    res <- spectreasy::unmix_samples(
        sample_dir = sample_dir,
        unmixing_matrix_file = ref_file,
        variances_file = var_file,
        method = "WLS",
        scc_dir = wf$scc_dir,
        control_file = wf$control_df,
        output_dir = output_dir,
        write_fcs = TRUE,
        save_report = FALSE
    )

    expect_s3_class(res, "spectreasy_unmixed_results")
    expect_false(file.exists(var_file))
    expect_true(file.exists(file.path(output_dir, "sample_unmixed.fcs")))
})

test_that("unmix_samples writes by default and suffixes existing outputs", {
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
            method = "OLS",
            output_dir = output_dir
        )
    )

    expect_false(call_result$visible)
    expect_setequal(names(call_result$value), c("sample_a", "sample_b"))
    expect_true(file.exists(file.path(output_dir, "sample_a_unmixed.fcs")))
    expect_true(file.exists(file.path(output_dir, "sample_b_unmixed.fcs")))
    expect_equal(attr(call_result$value, "qc_report_file"), file.path(output_dir, "qc_samples_report.pdf"))
    expect_true(file.exists(file.path(output_dir, "qc_samples_report.pdf")))

    suffixed_result <- spectreasy::unmix_samples(
        sample_dir = sample_dir,
        M = M,
        method = "OLS",
        output_dir = output_dir
    )
    expect_s3_class(suffixed_result, "spectreasy_unmixed_results")
    expect_true(file.exists(file.path(output_dir, "sample_a_unmixed_2.fcs")))
    expect_true(file.exists(file.path(output_dir, "sample_b_unmixed_2.fcs")))
    expect_equal(attr(suffixed_result, "qc_report_file"), file.path(output_dir, "qc_samples_report_2.pdf"))
    expect_true(file.exists(file.path(output_dir, "qc_samples_report_2.pdf")))

    if (run_slow_tests()) {
        qc_png_dir <- tempfile("spectreasy_sample_qc_pngs_")
        res_with_pngs <- spectreasy::unmix_samples(
            sample_dir = sample_dir,
            M = M,
            method = "OLS",
            output_dir = tempfile("spectreasy_covr_unmixed_png_"),
            save_report = TRUE,
            save_qc_plots = TRUE,
            qc_plot_dir = qc_png_dir,
            write_fcs = FALSE
        )
        expect_true(dir.exists(attr(res_with_pngs, "qc_plot_dir")))
        expect_true(length(list.files(qc_png_dir, pattern = "\\.png$", full.names = TRUE)) > 0)
    }
})

test_that("qc_controls writes a PDF from synthetic SCC files", {
    skip_slow_tests("SCC QC PDF generation")

    wf <- make_synthetic_workflow()
    output_pdf <- tempfile(fileext = ".pdf")
    qc_plot_dir <- tempfile("spectreasy_covr_scc_report_")
    control_csv <- tempfile(fileext = ".csv")
    utils::write.csv(wf$control_df, control_csv, row.names = FALSE, quote = TRUE)

    out <- spectreasy::qc_controls(
        scc_dir = wf$scc_dir,
        output_file = output_pdf,
        control_file = control_csv,
        method = "NNLS",
        qc_plot_dir = qc_plot_dir,
        save_qc_pngs = TRUE,
        seed = 1,
        subsample_n = 400
    )

    expect_true(file.exists(output_pdf))
    expect_true(dir.exists(out$qc_plot_dir))
    expect_true(is.matrix(out$M))
    expect_true(is.data.frame(out$qc_summary))
    expect_equal(out$method, "NNLS")

    # Test default output file behavior
    default_pdf <- "spectreasy_outputs/unmix_controls/qc_controls_report.pdf"
    if (file.exists(default_pdf)) file.remove(default_pdf)
    
    out_default <- spectreasy::qc_controls(
        scc_dir = wf$scc_dir,
        control_file = control_csv,
        qc_plot_dir = qc_plot_dir,
        save_qc_pngs = FALSE,
        seed = 1,
        subsample_n = 400
    )
    expect_true(file.exists(default_pdf))
    file.remove(default_pdf)
})

test_that("qc_controls does not retain QC PNGs unless requested", {
    skip_slow_tests("SCC QC PDF generation")

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
    skip_slow_tests("packaged GUI localhost smoke test")
    skip_on_os("windows")
    skip_if_not_installed("plumber")

    expect_error(
        spectreasy::adjust_matrix(matrix_dir = tempfile("spectreasy_missing_matrix_"), open_browser = FALSE),
        regexp = "cannot be found|No such file|mustWork"
    )

    tmp_matrix_dir <- tempfile("spectreasy_gui_matrix_")
    dir.create(tmp_matrix_dir, recursive = TRUE, showWarnings = FALSE)
    port <- sample(18000:18999, 1)

    job <- parallel::mcparallel({
        spectreasy::adjust_matrix(
            matrix_dir = tmp_matrix_dir,
            open_browser = FALSE,
            dev_mode = FALSE,
            port = port
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
})
