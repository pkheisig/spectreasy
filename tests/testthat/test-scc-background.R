make_matched_af_exprs <- function(n = 900, positive = FALSE) {
    pop <- if (isTRUE(positive)) {
        c(rep(1, n %/% 2), rep(0, n - n %/% 2))
    } else {
        stats::rbinom(n, size = 1, prob = 0.35)
    }
    fsc <- rnorm(n, mean = ifelse(pop == 1, 108000, 96000), sd = 4200)
    ssc <- rnorm(n, mean = ifelse(pop == 1, 54000, 44000), sd = 3200)
    af_b1 <- 0.0018 * fsc + rnorm(n, sd = 8)
    af_yg1 <- 0.0040 * ssc + rnorm(n, sd = 8)
    dye_amp <- if (isTRUE(positive)) {
        c(runif(n %/% 2, 760, 860), rep(0, n - n %/% 2))
    } else {
        rep(0, n)
    }
    exprs <- cbind(
        "B1-A" = af_b1 + dye_amp,
        "YG1-A" = af_yg1 + 0.10 * dye_amp,
        "FSC-A" = fsc,
        "SSC-A" = ssc
    )
    exprs[exprs < 1] <- 1
    exprs
}

make_matched_af_scc <- function() {
    scc_dir <- tempfile("spectreasy_scc_bg_")
    dir.create(scc_dir, recursive = TRUE, showWarnings = FALSE)
    flowCore::write.FCS(
        flowCore::flowFrame(make_matched_af_exprs(positive = FALSE)),
        file.path(scc_dir, "Unstained (Cells).fcs")
    )
    flowCore::write.FCS(
        flowCore::flowFrame(make_matched_af_exprs(positive = TRUE)),
        file.path(scc_dir, "FITC (Cells).fcs")
    )
    control_df <- data.frame(
        filename = c("Unstained (Cells).fcs", "FITC (Cells).fcs"),
        fluorophore = c("AF", "FITC"),
        marker = c("Autofluorescence", "CD4"),
        channel = c("B1-A", "B1-A"),
        control.type = c("cells", "cells"),
        universal.negative = c("", ""),
        is.viability = c("", ""),
        stringsAsFactors = FALSE
    )
    list(scc_dir = scc_dir, control_df = control_df)
}

test_that("cell SCC spectra can use scatter-matched unstained background", {
    set.seed(44)
    wf <- make_matched_af_scc()

    raw <- spectreasy::build_reference_matrix(
        input_folder = wf$scc_dir,
        control_df = wf$control_df,
        clean_scc_with_unstained = FALSE,
        save_qc_plots = FALSE,
        seed = 44,
        subsample_n = 500
    )
    clean <- spectreasy::build_reference_matrix(
        input_folder = wf$scc_dir,
        control_df = wf$control_df,
        clean_scc_with_unstained = TRUE,
        scc_background_k = 3,
        save_qc_plots = FALSE,
        seed = 44,
        subsample_n = 500
    )

    expect_gt(raw["FITC", "YG1-A"], 0.12)
    expect_lt(clean["FITC", "YG1-A"], raw["FITC", "YG1-A"])
    expect_lt(abs(clean["FITC", "YG1-A"] - 0.10), abs(raw["FITC", "YG1-A"] - 0.10))
    expect_equal(attr(clean, "qc_summary")$scc_background_method[[1]], "scatter_knn")

    positive_events <- attr(clean, "scc_positive_events")
    expect_true("FITC" %in% names(positive_events))
    qc_summary <- attr(clean, "qc_summary")
    expect_equal(nrow(positive_events$FITC), qc_summary$n_final[qc_summary$fluorophore == "FITC"])
    positive_shapes <- spectreasy:::.normalize_spectral_variant_shapes(positive_events$FITC)
    expect_lt(abs(stats::median(positive_shapes[, "YG1-A"]) - 0.10), 0.06)
})

test_that("bead SCCs keep intensity selection even when AF controls exist", {
    set.seed(45)
    scc_dir <- tempfile("spectreasy_bead_gate_")
    dir.create(scc_dir, recursive = TRUE, showWarnings = FALSE)

    af_exprs <- make_matched_af_exprs(positive = FALSE)
    unstained_beads <- cbind(
        "B1-A" = stats::rnorm(600, 200, 12),
        "YG1-A" = stats::rnorm(600, 320, 14),
        "FSC-A" = stats::rnorm(600, 70000, 2500),
        "SSC-A" = stats::rnorm(600, 35000, 1800)
    )
    bead_exprs <- cbind(
        "B1-A" = stats::rnorm(600, 1200, 35),
        "YG1-A" = stats::rnorm(600, 420, 12),
        "FSC-A" = stats::rnorm(600, 70000, 2500),
        "SSC-A" = stats::rnorm(600, 35000, 1800)
    )
    unstained_beads[unstained_beads < 1] <- 1
    bead_exprs[bead_exprs < 1] <- 1

    flowCore::write.FCS(
        flowCore::flowFrame(unstained_beads),
        file.path(scc_dir, "Unstained (Beads).fcs")
    )
    flowCore::write.FCS(
        flowCore::flowFrame(af_exprs),
        file.path(scc_dir, "Unstained (Cells).fcs")
    )
    flowCore::write.FCS(
        flowCore::flowFrame(bead_exprs),
        file.path(scc_dir, "FITC (Beads).fcs")
    )
    control_df <- data.frame(
        filename = c("Unstained (Beads).fcs", "Unstained (Cells).fcs", "FITC (Beads).fcs"),
        fluorophore = c("AF", "AF", "FITC"),
        marker = c("Autofluorescence", "Autofluorescence", "FITC"),
        channel = c("B1-A", "B1-A", "B1-A"),
        control.type = c("beads", "cells", "beads"),
        universal.negative = c("", "", ""),
        is.viability = c("", "", ""),
        stringsAsFactors = FALSE
    )

    M <- spectreasy::build_reference_matrix(
        input_folder = scc_dir,
        control_df = control_df,
        clean_scc_with_unstained = TRUE,
        save_qc_plots = FALSE,
        seed = 45,
        subsample_n = 500
    )

    qc_summary <- attr(M, "qc_summary")
    expect_equal(qc_summary$intensity_gate_type[qc_summary$fluorophore == "FITC"], "scatter")
    expect_equal(qc_summary$scc_background_method[qc_summary$fluorophore == "FITC"], "none")
    expect_gt(M["FITC", "B1-A"], 0.95)
    expect_lt(abs(M["FITC", "YG1-A"] - 0.10), 0.04)
    expect_equal(attr(M, "af_bank_info")$sources$file[1], "Unstained (Cells).fcs")
})

test_that("SCC QC plots include spectral selection spectra", {
    set.seed(46)
    wf <- make_matched_af_scc()
    output_dir <- tempfile("spectreasy_spectral_selection_qc_")

    spectreasy::build_reference_matrix(
        input_folder = wf$scc_dir,
        output_folder = output_dir,
        control_df = wf$control_df,
        clean_scc_with_unstained = TRUE,
        save_qc_plots = TRUE,
        seed = 46,
        subsample_n = 500
    )

    expect_true(file.exists(file.path(output_dir, "spectral_selection", "FITC (Cells)_spectral_selection.png")))
})

test_that("matched background helpers clean local AF from variant events", {
    bg_exprs <- make_matched_af_exprs(n = 300, positive = FALSE)
    pos_exprs <- make_matched_af_exprs(n = 120, positive = TRUE)
    background <- list(
        method = "scatter_knn",
        detector_names = c("B1-A", "YG1-A"),
        scatter_names = c("FSC-A", "SSC-A"),
        scatter = bg_exprs[, c("FSC-A", "SSC-A")],
        spectra = bg_exprs[, c("B1-A", "YG1-A")],
        n_events = nrow(bg_exprs)
    )

    clean <- spectreasy:::.scc_background_clean_events(
        events = pos_exprs[pos_exprs[, "B1-A"] > stats::median(pos_exprs[, "B1-A"]), ],
        detector_names = c("B1-A", "YG1-A"),
        background = background,
        k = 3
    )
    shapes <- spectreasy:::.normalize_spectral_variant_shapes(clean)

    expect_lt(abs(stats::median(shapes[, "YG1-A"]) - 0.10), 0.06)
})
