make_autospectral_test_exprs <- function(n = 600, dye = c("B1-A" = 0, "YG1-A" = 0, "R1-A" = 0)) {
    stopifnot(all(c("B1-A", "YG1-A", "R1-A") %in% names(dye)))
    pos <- seq_len(n) <= (n %/% 2)
    fsc <- stats::rnorm(n, mean = 100000, sd = 3500)
    ssc <- stats::rnorm(n, mean = 50000, sd = 2500)
    af_b1 <- 0.0015 * fsc + stats::rnorm(n, sd = 5)
    af_yg1 <- 0.0035 * ssc + stats::rnorm(n, sd = 5)
    af_r1 <- 0.0010 * ssc + stats::rnorm(n, sd = 4)
    amp <- numeric(n)
    amp[pos] <- stats::runif(sum(pos), 650, 950)
    exprs <- cbind(
        "B1-A" = af_b1 + amp * dye[["B1-A"]],
        "YG1-A" = af_yg1 + amp * dye[["YG1-A"]],
        "R1-A" = af_r1 + amp * dye[["R1-A"]],
        "FSC-A" = fsc,
        "SSC-A" = ssc
    )
    exprs[exprs < 1] <- 1
    exprs
}

make_autospectral_test_workflow <- function(n = 600) {
    scc_dir <- tempfile("spectreasy_as_scc_")
    dir.create(scc_dir, recursive = TRUE, showWarnings = FALSE)

    flowCore::write.FCS(
        flowCore::flowFrame(make_autospectral_test_exprs(n = n)),
        file.path(scc_dir, "Unstained (Cells).fcs")
    )
    flowCore::write.FCS(
        flowCore::flowFrame(make_autospectral_test_exprs(n = n, dye = c("B1-A" = 1.0, "YG1-A" = 0.10, "R1-A" = 0.02))),
        file.path(scc_dir, "FITC (Cells).fcs")
    )
    flowCore::write.FCS(
        flowCore::flowFrame(make_autospectral_test_exprs(n = n, dye = c("B1-A" = 0.12, "YG1-A" = 1.0, "R1-A" = 0.03))),
        file.path(scc_dir, "PE (Cells).fcs")
    )

    control_df <- data.frame(
        filename = c("Unstained (Cells).fcs", "FITC (Cells).fcs", "PE (Cells).fcs"),
        fluorophore = c("AF", "FITC", "PE"),
        marker = c("Autofluorescence", "CD4", "CD8"),
        channel = c("B1-A", "B1-A", "YG1-A"),
        control.type = c("cells", "cells", "cells"),
        universal.negative = c("", "", ""),
        is.viability = c("", "", ""),
        stringsAsFactors = FALSE
    )
    list(scc_dir = scc_dir, control_df = control_df)
}

test_that("saved AF profiles provide their bank and background without a mapped unstained file", {
    set.seed(31)
    workflow <- make_autospectral_test_workflow(n = 600)
    unstained <- file.path(workflow$scc_dir, "Unstained (Cells).fcs")
    profile <- spectreasy::extract_af_profile(
        unstained,
        af_n_bands = 1,
        seed = 31,
        show_plot = FALSE,
        verbose = FALSE
    )
    expect_true(unlink(unstained) == 0L)
    marker_controls <- workflow$control_df[workflow$control_df$fluorophore != "AF", , drop = FALSE]

    matrix <- spectreasy::build_reference_matrix(
        input_folder = workflow$scc_dir,
        control_df = marker_controls,
        af_profile = profile,
        af_n_bands = 1,
        unmixing_method = "AutoSpectral",
        save_qc_plots = FALSE,
        seed = 31
    )

    expect_true("AF" %in% rownames(matrix))
    expect_identical(attr(matrix, "af_bank_info")$mode, "saved_profile")
    expect_equal(attr(matrix, "af_bank_info")$source_count, 0L)
})

test_that("AutoSpectral calc_residuals uses OLS when no variants are supplied", {
    M <- matrix(c(
        1, 0.1,
        0.2, 1
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B1-A", "YG1-A")
    exprs <- matrix(c(
        100, 20,
        30, 90,
        40, 45
    ), nrow = 3, byrow = TRUE)
    colnames(exprs) <- colnames(M)
    ff <- flowCore::flowFrame(exprs)

    ols <- spectreasy::calc_residuals(ff, M, method = "OLS")
    autospectral <- spectreasy::calc_residuals(ff, M, method = "AutoSpectral")

    expect_equal(autospectral[, rownames(M)], ols[, rownames(M)], tolerance = 1e-10)
})

test_that("AutoSpectral calc_residuals records per-cell AF assignment for AF banks", {
    M <- rbind(
        FITC = c(1.00, 0.10, 0.05),
        PE = c(0.10, 1.00, 0.05),
        AF = c(0.30, 0.20, 1.00),
        AF_2 = c(0.05, 0.20, 0.35)
    )
    colnames(M) <- c("B1-A", "YG1-A", "R1-A")
    exprs <- rbind(
        c(FITC = 200, PE = 20, AF = 80, AF_2 = 0) %*% M,
        c(FITC = 20, PE = 200, AF = 0, AF_2 = 80) %*% M
    )
    colnames(exprs) <- colnames(M)
    ff <- flowCore::flowFrame(exprs)

    res <- spectreasy::calc_residuals(ff, M, method = "AutoSpectral")

    expect_true("AF Index" %in% colnames(res))
    expect_true(all(res[["AF Index"]] %in% c(1L, 2L)))
    expect_true(all(rowSums(res[, c("AF", "AF_2")] != 0) <= 1))
})

test_that("Spectreasy calc_residuals applies AS-only soft decoder weights", {
    M <- rbind(
        FITC = c(1.00, 0.10, 0.05),
        PE = c(0.12, 1.00, 0.08),
        AF = c(0.80, 0.12, 1.00),
        AF_2 = c(0.10, 0.45, 0.80)
    )
    colnames(M) <- c("B1-A", "YG1-A", "R1-A")
    attr(M, "detector_noise") <- data.frame(
        detector = colnames(M),
        noise_floor = c(20, 20, 20),
        signal_scale = c(0.02, 0.02, 0.02),
        max_weight_ratio = c(100, 100, 100)
    )
    exprs <- rbind(
        c(FITC = 250, PE = 20, AF = 90, AF_2 = 0) %*% M,
        c(FITC = 20, PE = 240, AF = 0, AF_2 = 80) %*% M,
        c(FITC = 120, PE = 90, AF = 50, AF_2 = 60) %*% M
    )
    colnames(exprs) <- colnames(M)
    ff <- flowCore::flowFrame(exprs)

    autospectral <- spectreasy::calc_residuals(ff, M, method = "AutoSpectral", return_residuals = TRUE)
    spectreasy <- spectreasy::calc_residuals(ff, M, method = "Spectreasy", return_residuals = TRUE)
    weights <- spectreasy:::.decoder_projected_af_marker_weights(M, spectreasy_weight_quantile = 0.9)
    baseline <- spectreasy:::.spectreasy_marker_only_baseline_fit(exprs, M)

    expect_equal(spectreasy$spectreasy_decoder_weights, weights, tolerance = 1e-12)
    for (marker in names(weights)) {
        expected <- baseline[, marker] + weights[[marker]] * (autospectral$data[[marker]] - baseline[, marker])
        expect_equal(spectreasy$data[[marker]], expected, tolerance = 1e-5)
    }
    expect_true("AF Index" %in% colnames(spectreasy$data))
    expect_true(all(spectreasy$spectreasy_decoder_weights >= 0 & spectreasy$spectreasy_decoder_weights <= 1))
    expect_equal(attr(spectreasy$spectreasy_decoder_weights, "quantile"), 0.9)
    expect_error(
        spectreasy::calc_residuals(ff, M, method = "Spectreasy", spectreasy_weight_quantile = 1.2),
        "spectreasy_weight_quantile"
    )
    expect_error(
        spectreasy::calc_residuals(ff, M, method = "WLS", spectreasy_weight_quantile = 0.9),
        "method = \"Spectreasy\""
    )
})

test_that("SCC processing follows the selected unmixing method", {
    set.seed(101)
    wf <- make_autospectral_test_workflow(n = 650)

    regular <- spectreasy::build_reference_matrix(
        input_folder = wf$scc_dir,
        control_df = wf$control_df,
        unmixing_method = "WLS",
        save_qc_plots = FALSE,
        seed = 101,
        subsample_n = 400
    )
    autospectral <- spectreasy::build_reference_matrix(
        input_folder = wf$scc_dir,
        control_df = wf$control_df,
        unmixing_method = "AutoSpectral",
        save_qc_plots = FALSE,
        seed = 101,
        subsample_n = 400
    )

    expect_null(attr(regular, "scc_positive_events"))
    positive_events <- attr(autospectral, "scc_positive_events")
    expect_true(is.list(positive_events))
    expect_true(all(c("FITC", "PE") %in% names(positive_events)))
    expect_true(all(vapply(positive_events[c("FITC", "PE")], nrow, integer(1)) > 0))

    qc_summary <- attr(autospectral, "qc_summary")
    expect_true(all(qc_summary$intensity_gate_type[qc_summary$fluorophore %in% c("FITC", "PE")] == "autospectral_external"))
    expect_true(all(qc_summary$scc_background_method[qc_summary$fluorophore %in% c("FITC", "PE")] %in% c("scatter_knn", "external_median")))
})

test_that("unmix_controls learns variants only for AutoSpectral-style methods", {
    set.seed(102)
    wf <- make_autospectral_test_workflow(n = 650)
    control_csv <- tempfile(fileext = ".csv")
    utils::write.csv(wf$control_df, control_csv, row.names = FALSE, quote = TRUE)

    regular_dir <- tempfile("spectreasy_as_regular_")
    regular <- spectreasy::unmix_controls(
        scc_dir = wf$scc_dir,
        control_file = control_csv,
        output_dir = regular_dir,
        unmixing_method = "OLS",
        save_report = FALSE,
        seed = 102,
        subsample_n = 400
    )
    expect_null(regular$spectral_variant_library)
    expect_false(file.exists(file.path(regular_dir, "unmix_controls", "scc_spectral_variants.rds")))

    autospectral_dir <- tempfile("spectreasy_as_enabled_")
    autospectral <- spectreasy::unmix_controls(
        scc_dir = wf$scc_dir,
        control_file = control_csv,
        output_dir = autospectral_dir,
        unmixing_method = "AutoSpectral",
        spectral_variant_min_events = 20,
        save_report = FALSE,
        seed = 102,
        subsample_n = 400
    )
    expect_s3_class(autospectral$spectral_variant_library, "spectreasy_spectral_variant_library")
    expect_true(file.exists(autospectral$spectral_variant_library_file))
    expect_equal(autospectral$static_unmixing_matrix_method, "AutoSpectral")
    expect_equal(attr(autospectral$unmixed_list, "method"), "AutoSpectral")

    spectreasy_dir <- tempfile("spectreasy_enabled_")
    spectreasy <- spectreasy::unmix_controls(
        scc_dir = wf$scc_dir,
        control_file = control_csv,
        output_dir = spectreasy_dir,
        unmixing_method = "Spectreasy",
        spectral_variant_min_events = 20,
        save_report = FALSE,
        seed = 102,
        subsample_n = 400
    )
    expect_s3_class(spectreasy$spectral_variant_library, "spectreasy_spectral_variant_library")
    expect_true(file.exists(spectreasy$spectral_variant_library_file))
    expect_equal(spectreasy$static_unmixing_matrix_method, "Spectreasy")
    expect_equal(attr(spectreasy$unmixed_list, "method"), "Spectreasy")
})

test_that("AF refinement helper returns fixed-size refined k-means bank", {
    set.seed(103)
    M <- rbind(
        FITC = c(1.00, 0.08, 0.02),
        PE = c(0.08, 1.00, 0.02),
        AF = c(0.65, 0.25, 1.00),
        AF_2 = c(0.25, 0.65, 1.00)
    )
    colnames(M) <- c("B1-A", "YG1-A", "R1-A")
    make_events <- function(shape, n) {
        amp <- stats::runif(n, 400, 1600)
        noise <- matrix(stats::rnorm(n * length(shape), sd = 8), ncol = length(shape))
        out <- noise + amp %o% shape
        colnames(out) <- names(shape)
        pmax(out, 1)
    }
    af_events <- rbind(
        make_events(c("B1-A" = 0.95, "YG1-A" = 0.18, "R1-A" = 1.00), 500),
        make_events(c("B1-A" = 0.18, "YG1-A" = 0.95, "R1-A" = 1.00), 500)
    )

    refined <- spectreasy:::.reference_refine_af_bank(
        M = M,
        af_events = af_events,
        af_n_bands = 2,
        af_max_cells = 1000,
        n_threads = 2L,
        seed = 103,
        verbose = FALSE
    )

    expect_false(is.null(refined))
    expect_equal(nrow(refined$signatures), 2)
    expect_equal(rownames(refined$signatures), c("AF", "AF_2"))
    expect_match(refined$selection$method, "refined_fixed")
    expect_equal(refined$selection$requested_bands, 2)
    expect_true(all(refined$signatures >= 0))
    expect_true(all(refined$signatures <= 1))
})

test_that("unmix_controls accepts autospectral_refine only for AutoSpectral", {
    expect_error(
        spectreasy::unmix_controls(
            scc_dir = tempfile("missing_scc_"),
            unmixing_method = "WLS",
            autospectral_refine = TRUE,
            save_report = FALSE
        ),
        "unmixing_method = \"AutoSpectral\""
    )
    expect_error(
        spectreasy::unmix_controls(
            scc_dir = tempfile("missing_scc_"),
            unmixing_method = "Spectreasy",
            autospectral_refine = TRUE,
            save_report = FALSE
        ),
        "unmixing_method = \"AutoSpectral\""
    )
})

test_that("spectreasy_weight_quantile is accepted only for Spectreasy", {
    M <- matrix(c(1, 0.1, 0.2, 1), nrow = 2, byrow = TRUE)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B1-A", "YG1-A")
    exprs <- matrix(c(100, 20, 30, 90), nrow = 2, byrow = TRUE)
    colnames(exprs) <- colnames(M)
    sample_dir <- tempfile("spectreasy_quantile_samples_")
    dir.create(sample_dir, recursive = TRUE)
    flowCore::write.FCS(flowCore::flowFrame(exprs), file.path(sample_dir, "sample.fcs"))

    expect_error(
        spectreasy::unmix_samples(
            sample_dir = sample_dir,
            M = M,
            unmixing_method = "WLS",
            spectreasy_weight_quantile = 0.9,
            write_fcs = FALSE,
            save_report = FALSE
        ),
        "unmixing_method = \"Spectreasy\""
    )
    expect_error(
        spectreasy::unmix_controls(
            scc_dir = tempfile("missing_scc_"),
            unmixing_method = "WLS",
            spectreasy_weight_quantile = 0.9,
            save_report = FALSE
        ),
        "unmixing_method = \"Spectreasy\""
    )
})

test_that("AutoSpectral refine rewrites AF bank metadata while preserving fixed band count", {
    set.seed(104)
    wf <- make_autospectral_test_workflow(n = 650)
    control_csv <- tempfile(fileext = ".csv")
    utils::write.csv(wf$control_df, control_csv, row.names = FALSE, quote = TRUE)

    testthat::local_mocked_bindings(
        .reference_refine_af_bank = function(M, af_events, af_n_bands, ...) {
            signatures <- matrix(0.2, nrow = af_n_bands, ncol = ncol(M))
            signatures[cbind(seq_len(af_n_bands), ((seq_len(af_n_bands) - 1L) %% ncol(M)) + 1L)] <- 1
            rownames(signatures) <- c("AF", if (af_n_bands > 1L) paste0("AF_", seq.int(2L, af_n_bands)) else NULL)
            colnames(signatures) <- colnames(M)
            list(
                signatures = signatures,
                selection = list(
                    method = "kmeans_refined_fixed",
                    n_bands = af_n_bands,
                    requested_bands = af_n_bands,
                    final_bands = af_n_bands
                )
            )
        },
        .env = asNamespace("spectreasy")
    )

    ctrl <- spectreasy::unmix_controls(
        scc_dir = wf$scc_dir,
        control_file = control_csv,
        output_dir = tempfile("spectreasy_as_refine_"),
        unmixing_method = "AutoSpectral",
        autospectral_refine = TRUE,
        af_n_bands = 2,
        spectral_variant_min_events = 20,
        save_report = FALSE,
        seed = 104,
        subsample_n = 400
    )

    af_rows <- grepl("^AF($|_)", rownames(ctrl$M), ignore.case = TRUE)
    af_info <- attr(ctrl$M, "af_bank_info")
    expect_equal(sum(af_rows), 2)
    expect_true(isTRUE(af_info$refine))
    expect_equal(af_info$selection$method, "kmeans_refined_fixed")
    expect_equal(af_info$selection$requested_bands, 2)
    expect_true(isTRUE(ctrl$autospectral_refine))
})
