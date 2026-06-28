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
        large.gate = c("", ""),
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
