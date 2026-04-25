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

test_that("autounmix_controls runs end-to-end on synthetic SCC files", {
    wf <- make_synthetic_workflow()
    output_dir <- tempfile("spectreasy_covr_auto_")

    ctrl <- spectreasy::autounmix_controls(
        scc_dir = wf$scc_dir,
        control_df = wf$control_df,
        output_dir = output_dir,
        unmix_method = "OLS",
        build_qc_plots = FALSE,
        seed = 1,
        subsample_n = 400
    )

    expect_true(file.exists(ctrl$reference_matrix_file))
    expect_true(file.exists(ctrl$unmixing_matrix_file))
    expect_true(file.exists(ctrl$spectra_file))
    expect_true(file.exists(ctrl$unmixing_matrix_plot))
    expect_true(file.exists(ctrl$unmixing_scatter_file))
    expect_equal(sort(names(ctrl$unmixed_list)), c("FITC (Beads)", "PE (Beads)"))
})

test_that("autounmix_controls handles WLS output and exclude_af branch", {
    wf <- make_synthetic_workflow(include_af = TRUE)
    output_dir <- tempfile("spectreasy_covr_auto_wls_")
    control_csv <- tempfile(fileext = ".csv")
    utils::write.csv(wf$control_df, control_csv, row.names = FALSE, quote = TRUE)

    ctrl <- spectreasy::autounmix_controls(
        scc_dir = wf$scc_dir,
        control_df = control_csv,
        output_dir = output_dir,
        exclude_af = TRUE,
        unmix_method = "WLS",
        build_qc_plots = FALSE,
        seed = 1,
        subsample_n = 400
    )

    expect_equal(ctrl$static_unmixing_matrix_method, "WLS")
    expect_false(any(grepl("^AF($|_)", rownames(ctrl$M), ignore.case = TRUE)))
    expect_false(any(grepl("Unstained", names(ctrl$unmixed_list), ignore.case = TRUE)))
    expect_false(file.exists(file.path(output_dir, "scc_unmixed", "Unstained (Cells)_unmixed.fcs")))
})

test_that("generate_scc_report writes a PDF from synthetic SCC files", {
    wf <- make_synthetic_workflow()
    output_pdf <- tempfile(fileext = ".pdf")
    qc_plot_dir <- tempfile("spectreasy_covr_scc_report_")

    out <- spectreasy::generate_scc_report(
        scc_dir = wf$scc_dir,
        output_file = output_pdf,
        control_df = wf$control_df,
        qc_plot_dir = qc_plot_dir,
        save_qc_pngs = TRUE,
        include_ssm = TRUE,
        seed = 1,
        subsample_n = 400
    )

    expect_true(file.exists(output_pdf))
    expect_true(dir.exists(out$qc_plot_dir))
    expect_true(is.matrix(out$M))
    expect_true(is.data.frame(out$qc_summary))
})

test_that("generate_scc_report does not retain QC PNGs unless requested", {
    wf <- make_synthetic_workflow()
    output_pdf <- tempfile(fileext = ".pdf")
    qc_plot_dir <- tempfile("spectreasy_scc_report_no_retain_")

    out <- spectreasy::generate_scc_report(
        scc_dir = wf$scc_dir,
        output_file = output_pdf,
        control_df = wf$control_df,
        qc_plot_dir = qc_plot_dir,
        include_ssm = FALSE,
        seed = 1,
        subsample_n = 400
    )

    expect_true(file.exists(output_pdf))
    expect_null(out$qc_plot_dir)
    expect_false(dir.exists(qc_plot_dir))
})

test_that("launch_gui starts packaged GUI on localhost", {
    skip_on_os("windows")
    skip_if_not_installed("plumber")

    expect_error(
        spectreasy::launch_gui(matrix_dir = tempfile("spectreasy_missing_matrix_"), open_browser = FALSE),
        regexp = "cannot be found|No such file|mustWork"
    )

    tmp_matrix_dir <- tempfile("spectreasy_gui_matrix_")
    dir.create(tmp_matrix_dir, recursive = TRUE, showWarnings = FALSE)
    port <- sample(18000:18999, 1)

    job <- parallel::mcparallel({
        spectreasy::launch_gui(
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
