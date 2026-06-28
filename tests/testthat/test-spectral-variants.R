make_variant_ff <- function(shapes, n_each = 120, bg = 20, noise_sd = 2) {
    shapes <- as.matrix(shapes)
    rows <- lapply(seq_len(nrow(shapes)), function(i) {
        amp <- runif(n_each, 600, 1200)
        sweep(matrix(shapes[i, ], nrow = n_each, ncol = ncol(shapes), byrow = TRUE), 1, amp, "*") +
            bg +
            matrix(rnorm(n_each * ncol(shapes), sd = noise_sd), ncol = ncol(shapes))
    })
    Y <- do.call(rbind, rows)
    colnames(Y) <- colnames(shapes)
    exprs <- cbind(
        Y,
        "FSC-A" = rnorm(nrow(Y), mean = 90000, sd = 4000),
        "SSC-A" = rnorm(nrow(Y), mean = 45000, sd = 2500)
    )
    flowCore::flowFrame(exprs)
}

make_variant_control_set <- function() {
    set.seed(11)
    scc_dir <- tempfile("spectreasy_variant_scc_")
    dir.create(scc_dir, recursive = TRUE, showWarnings = FALSE)
    detectors <- c("B1-A", "YG1-A", "R1-A")
    fitc_shapes <- rbind(
        base = c(1.00, 0.18, 0.02),
        shifted = c(0.78, 0.42, 0.02)
    )
    pe_shapes <- rbind(base = c(0.04, 0.18, 1.00))
    colnames(fitc_shapes) <- colnames(pe_shapes) <- detectors
    flowCore::write.FCS(make_variant_ff(fitc_shapes), file.path(scc_dir, "FITC (Beads).fcs"))
    flowCore::write.FCS(make_variant_ff(pe_shapes), file.path(scc_dir, "PE (Beads).fcs"))
    control_df <- data.frame(
        filename = c("FITC (Beads).fcs", "PE (Beads).fcs"),
        fluorophore = c("FITC", "PE"),
        marker = c("CD4", "CD8"),
        channel = c("B1-A", "R1-A"),
        control.type = c("beads", "beads"),
        universal.negative = c("", ""),
        large.gate = c("", ""),
        is.viability = c("", ""),
        stringsAsFactors = FALSE
    )
    list(scc_dir = scc_dir, control_df = control_df, detectors = detectors)
}

test_that("unmix_controls variant optimization can be disabled", {
    wf <- make_variant_control_set()
    output_dir <- tempfile("spectreasy_variant_disabled_")
    control_csv <- tempfile(fileext = ".csv")
    utils::write.csv(wf$control_df, control_csv, row.names = FALSE, quote = TRUE)

    ctrl <- spectreasy::unmix_controls(
        scc_dir = wf$scc_dir,
        control_file = control_csv,
        output_dir = output_dir,
        unmix_method = "OLS",
        optimize_spectral_variants = FALSE,
        save_qc_plots = FALSE,
        save_report = FALSE,
        seed = 1
    )

    expect_true(file.exists(ctrl$reference_matrix_file))
    expect_null(ctrl$spectral_variant_library)
    expect_null(ctrl$spectral_variant_library_file)
    expect_false(file.exists(file.path(output_dir, "scc_spectral_variants.rds")))
})

test_that("unmix_controls creates a spectral variant library from SCC controls", {
    wf <- make_variant_control_set()
    output_dir <- tempfile("spectreasy_variant_enabled_")
    control_csv <- tempfile(fileext = ".csv")
    utils::write.csv(wf$control_df, control_csv, row.names = FALSE, quote = TRUE)

    ctrl <- spectreasy::unmix_controls(
        scc_dir = wf$scc_dir,
        control_file = control_csv,
        output_dir = output_dir,
        unmix_method = "OLS",
        optimize_spectral_variants = TRUE,
        spectral_variant_som_nodes = 4,
        spectral_variant_cosine_threshold = 0.90,
        spectral_variant_min_events = 20,
        save_report = FALSE,
        seed = 2
    )

    expect_s3_class(ctrl$spectral_variant_library, "spectreasy_spectral_variant_library")
    expect_true(file.exists(ctrl$spectral_variant_library_file))
    expect_true("FITC" %in% names(ctrl$spectral_variant_library$variants))
    expect_gte(nrow(ctrl$spectral_variant_library$variants$FITC), 1)
    expect_true("FITC" %in% ctrl$spectral_variant_info$fluorophore)
})

test_that("spectral variants reuse cleaned positive SCC events from the reference matrix", {
    set.seed(13)
    detectors <- c("B1-A", "YG1-A", "R1-A")
    M <- rbind(
        FITC = c(1.00, 0.18, 0.02),
        PE = c(0.04, 0.18, 1.00)
    )
    colnames(M) <- detectors
    shifted <- c(0.78, 0.42, 0.02)
    events <- rbind(
        sweep(matrix(M["FITC", ], nrow = 80, ncol = length(detectors), byrow = TRUE), 1, runif(80, 700, 1100), "*"),
        sweep(matrix(shifted, nrow = 80, ncol = length(detectors), byrow = TRUE), 1, runif(80, 700, 1100), "*")
    )
    colnames(events) <- detectors
    attr(M, "scc_positive_events") <- list(FITC = events)

    scc_dir <- tempfile("spectreasy_variant_stored_events_")
    dir.create(scc_dir, recursive = TRUE, showWarnings = FALSE)
    control_df <- data.frame(
        filename = "FITC file intentionally missing.fcs",
        fluorophore = "FITC",
        marker = "CD4",
        channel = "B1-A",
        control.type = "cells",
        stringsAsFactors = FALSE
    )

    lib <- spectreasy:::.learn_spectral_variant_library(
        scc_dir = scc_dir,
        control_df = control_df,
        M = M,
        enabled = TRUE,
        som_nodes = 4,
        cosine_threshold = 0.85,
        max_variants = 4,
        min_events = 20,
        warn = FALSE,
        seed = 13
    )

    expect_true("FITC" %in% names(lib$variants))
    expect_gte(nrow(lib$variants$FITC), 1)
    expect_equal(lib$info$event_count[lib$info$fluorophore == "FITC"], nrow(events))
})

test_that("spectral variant learning skips instead of falling back without stored SCC positives", {
    wf <- make_variant_control_set()
    M <- rbind(
        FITC = c(1.00, 0.18, 0.02),
        PE = c(0.04, 0.18, 1.00)
    )
    colnames(M) <- wf$detectors

    expect_warning(
        lib <- spectreasy:::.learn_spectral_variant_library(
            scc_dir = wf$scc_dir,
            control_df = wf$control_df,
            M = M,
            enabled = TRUE,
            som_nodes = 4,
            cosine_threshold = 0.85,
            max_variants = 4,
            min_events = 20,
            seed = 13
        ),
        "No `scc_positive_events` attribute was found"
    )

    expect_s3_class(lib, "spectreasy_spectral_variant_library")
    expect_false(lib$enabled)
    expect_length(lib$variants, 0)
    expect_equal(nrow(lib$info), 0)
})

test_that("per-cell spectral variants improve shifted fluorophore recovery", {
    M <- rbind(
        FITC = c(1.00, 0.15, 0.02),
        PE = c(0.02, 0.20, 1.00)
    )
    colnames(M) <- c("B1-A", "YG1-A", "R1-A")
    shifted <- matrix(c(0.78, 0.42, 0.02), nrow = 1, dimnames = list("FITC_variant_1", colnames(M)))
    lib <- list(
        enabled = TRUE,
        detector_names = colnames(M),
        variants = list(FITC = shifted),
        info = data.frame(fluorophore = "FITC", variants = 1L, reason = "ok"),
        settings = list()
    )
    class(lib) <- c("spectreasy_spectral_variant_library", "list")

    n <- 80
    truth <- cbind(FITC = runif(n, 100, 250), PE = runif(n, 40, 100))
    Y <- sweep(matrix(shifted[1, ], nrow = n, ncol = ncol(M), byrow = TRUE), 1, truth[, "FITC"], "*") +
        sweep(matrix(M["PE", ], nrow = n, ncol = ncol(M), byrow = TRUE), 1, truth[, "PE"], "*") +
        matrix(rnorm(n * ncol(M), sd = 0.5), ncol = ncol(M))
    colnames(Y) <- colnames(M)
    ff <- flowCore::flowFrame(Y)

    base <- spectreasy::calc_residuals(ff, M, method = "OLS")
    opt <- spectreasy::calc_residuals(
        ff,
        M,
        method = "OLS",
        spectral_variant_library = lib,
        optimize_spectral_variants = TRUE,
        spectral_variant_min_abundance = 5,
        return_residuals = TRUE
    )

    base_err <- mean(abs(as.matrix(base[, c("FITC", "PE")]) - truth))
    opt_err <- mean(abs(as.matrix(opt$data[, c("FITC", "PE")]) - truth))
    expect_lt(opt_err, base_err)
    expect_gt(opt$spectral_variant_info$changed_events, 0)
})

test_that("spectral variants do not switch on pure detector noise", {
    M <- rbind(
        FITC = c(1.00, 0.15, 0.02),
        PE = c(0.02, 0.20, 1.00)
    )
    colnames(M) <- c("B1-A", "YG1-A", "R1-A")
    lib <- list(
        enabled = TRUE,
        detector_names = colnames(M),
        variants = list(FITC = matrix(c(0.78, 0.42, 0.02), nrow = 1, dimnames = list("FITC_variant_1", colnames(M)))),
        info = data.frame(fluorophore = "FITC", variants = 1L, reason = "ok"),
        settings = list()
    )
    class(lib) <- c("spectreasy_spectral_variant_library", "list")
    Y <- matrix(rnorm(90, sd = 0.03), ncol = 3)
    colnames(Y) <- colnames(M)

    res <- spectreasy::calc_residuals(
        flowCore::flowFrame(Y),
        M,
        method = "OLS",
        spectral_variant_library = lib,
        optimize_spectral_variants = TRUE,
        spectral_variant_min_abundance = 5,
        return_residuals = TRUE
    )

    expect_equal(res$spectral_variant_info$changed_events, 0)
})

test_that("unmix_samples reuses a saved spectral variant library", {
    M <- rbind(
        FITC = c(1.00, 0.15, 0.02),
        PE = c(0.02, 0.20, 1.00)
    )
    colnames(M) <- c("B1-A", "YG1-A", "R1-A")
    shifted <- matrix(c(0.78, 0.42, 0.02), nrow = 1, dimnames = list("FITC_variant_1", colnames(M)))
    lib <- list(
        enabled = TRUE,
        detector_names = colnames(M),
        variants = list(FITC = shifted),
        info = data.frame(fluorophore = "FITC", variants = 1L, reason = "ok"),
        settings = list()
    )
    class(lib) <- c("spectreasy_spectral_variant_library", "list")

    matrix_dir <- tempfile("spectreasy_variant_matrix_")
    dir.create(matrix_dir, recursive = TRUE)
    ref_file <- file.path(matrix_dir, "scc_reference_matrix.csv")
    spectreasy:::.save_reference_matrix_csv(M, ref_file)
    spectreasy:::.save_spectral_variant_library(lib, file.path(matrix_dir, "scc_spectral_variants.rds"))

    sample_dir <- tempfile("spectreasy_variant_samples_")
    dir.create(sample_dir, recursive = TRUE)
    truth <- c(FITC = 180, PE = 60)
    Y <- matrix(truth[["FITC"]] * shifted[1, ] + truth[["PE"]] * M["PE", ], nrow = 1)
    colnames(Y) <- colnames(M)
    flowCore::write.FCS(flowCore::flowFrame(Y), file.path(sample_dir, "sample.fcs"))

    res <- spectreasy::unmix_samples(
        sample_dir = sample_dir,
        unmixing_matrix_file = ref_file,
        method = "OLS",
        write_fcs = FALSE,
        save_report = FALSE,
        optimize_spectral_variants = TRUE,
        spectral_variant_min_abundance = 5,
        verbose = FALSE
    )

    expect_s3_class(attr(res, "spectral_variant_library"), "spectreasy_spectral_variant_library")
    expect_gt(res$sample$spectral_variant_info$changed_events, 0)
})
