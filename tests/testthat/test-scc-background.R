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

test_that("AF cosine SCC selection enforces scatter-matched SCC cleaning", {
    expect_warning(
        args <- spectreasy:::.validate_scc_background_args(
            clean_scc_with_unstained = FALSE,
            scc_background_method = "none",
            scc_background_k = 3L,
            require_for_af_cosine = TRUE
        ),
        "requires scatter-matched unstained SCC cleaning"
    )
    expect_true(args$enabled)
    expect_equal(args$method, "scatter_knn")
    expect_equal(args$k, 3L)

    relaxed <- spectreasy:::.validate_scc_background_args(
        clean_scc_with_unstained = FALSE,
        scc_background_method = "none",
        scc_background_k = 3L
    )
    expect_false(relaxed$enabled)
    expect_equal(relaxed$method, "none")
})

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
    expect_equal(qc_summary$intensity_gate_type[qc_summary$fluorophore == "FITC"], "scatter")
    expect_equal(nrow(positive_events$FITC), qc_summary$n_final[qc_summary$fluorophore == "FITC"])
    positive_shapes <- spectreasy:::.normalize_spectral_variant_shapes(positive_events$FITC)
    expect_lt(abs(stats::median(positive_shapes[, "YG1-A"]) - 0.10), 0.06)
})

test_that("AF cosine SCC selection is opt-in for cell controls", {
    set.seed(47)
    wf <- make_matched_af_scc()

    default <- spectreasy::build_reference_matrix(
        input_folder = wf$scc_dir,
        control_df = wf$control_df,
        clean_scc_with_unstained = TRUE,
        save_qc_plots = FALSE,
        seed = 47,
        subsample_n = 500
    )
    opt_in <- spectreasy::build_reference_matrix(
        input_folder = wf$scc_dir,
        control_df = wf$control_df,
        clean_scc_with_unstained = TRUE,
        use_af_cosine_scc_selection = TRUE,
        save_qc_plots = FALSE,
        seed = 47,
        subsample_n = 500
    )

    default_summary <- attr(default, "qc_summary")
    opt_in_summary <- attr(opt_in, "qc_summary")
    expect_equal(default_summary$intensity_gate_type[default_summary$fluorophore == "FITC"], "scatter")
    expect_equal(opt_in_summary$intensity_gate_type[opt_in_summary$fluorophore == "FITC"], "af_cosine")
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

test_that("SCC report prefers spectral-selection plot only for AF-cosine gates", {
    testthat::skip_if_not_installed("pdftools")

    plot_dir <- tempfile("spectreasy_report_gate_choice_")
    sample_id <- "FITC (Cells)"
    for (subdir in c("fsc_ssc", "intensity_scatter", "spectral_selection", "spectrum")) {
        dir.create(file.path(plot_dir, subdir), recursive = TRUE, showWarnings = FALSE)
    }
    write_blank_png <- function(path) {
        grDevices::png(path, width = 80, height = 80)
        on.exit(grDevices::dev.off(), add = TRUE)
        grid::grid.newpage()
        grid::grid.rect(gp = grid::gpar(col = "grey80", fill = "white"))
    }
    write_blank_png(file.path(plot_dir, "fsc_ssc", paste0(sample_id, "_fsc_ssc.png")))
    write_blank_png(file.path(plot_dir, "intensity_scatter", paste0(sample_id, "_intensity_scatter.png")))
    write_blank_png(file.path(plot_dir, "spectral_selection", paste0(sample_id, "_spectral_selection.png")))
    write_blank_png(file.path(plot_dir, "spectrum", paste0(sample_id, "_spectrum.png")))

    row <- data.table::data.table(
        sample = sample_id,
        fluorophore = "FITC",
        marker = "",
        type = "cells",
        peak_channel = "B1-A",
        n_total = 100L,
        n_scatter_gated = 80L,
        scatter_gate_pct = 80,
        n_final = 40L,
        histogram_gate_pct = 50,
        intensity_gate_type = "scatter"
    )

    scatter_pdf <- tempfile(fileext = ".pdf")
    grDevices::pdf(scatter_pdf, width = 11, height = 8.5)
    spectreasy:::.draw_scc_report_sample_page(row, report_plot_dir = plot_dir, use_scatter_gating = TRUE)
    grDevices::dev.off()
    scatter_text <- paste(pdftools::pdf_text(scatter_pdf), collapse = "\n")
    expect_match(scatter_text, "Scatter/Spectral Event Gate", fixed = TRUE)
    expect_false(grepl("SCC/AF Spectral Selection", scatter_text, fixed = TRUE))

    row$intensity_gate_type <- "af_cosine"
    cosine_pdf <- tempfile(fileext = ".pdf")
    grDevices::pdf(cosine_pdf, width = 11, height = 8.5)
    spectreasy:::.draw_scc_report_sample_page(row, report_plot_dir = plot_dir, use_scatter_gating = TRUE)
    grDevices::dev.off()
    cosine_text <- paste(pdftools::pdf_text(cosine_pdf), collapse = "\n")
    expect_match(cosine_text, "SCC/AF Spectral Selection", fixed = TRUE)
})

test_that("post-unmixing control QC summarizes off-target metrics", {
    markers <- c("FITC", "PE")
    unmixed_list <- list(
        "FITC (Cells)" = list(data = data.frame(
            FITC = seq(100, 200, length.out = 80),
            PE = c(rep(0, 60), rep(15, 20)),
            check.names = FALSE
        )),
        "Unstained (Cells)" = list(data = data.frame(
            FITC = stats::rnorm(80, 0, 1),
            PE = stats::rnorm(80, 0, 1),
            check.names = FALSE
        ))
    )
    qc_summary <- data.frame(
        sample = "FITC (Cells)",
        fluorophore = "FITC",
        type = "cells",
        stringsAsFactors = FALSE
    )

    metrics <- spectreasy:::.compute_scc_post_unmix_qc(
        unmixed_list = unmixed_list,
        qc_summary = qc_summary,
        markers = markers
    )

    expect_s3_class(metrics$overview, "data.frame")
    expect_s3_class(metrics$pairs, "data.frame")
    expect_equal(metrics$overview$target, "FITC")
    expect_equal(metrics$pairs$marker, "PE")
    expect_true(all(c("nps", "bias", "fpr", "slope") %in% colnames(metrics$pairs)))
    expect_gte(metrics$pairs$fpr, 0)
})

test_that("post-unmixing control QC pages render into PDF reports", {
    testthat::skip_if_not_installed("pdftools")

    markers <- c("FITC", "PE")
    unmixed_list <- list(
        "FITC (Beads)" = list(data = data.frame(
            FITC = stats::rnorm(80, 1200, 80),
            PE = stats::rnorm(80, 30, 5),
            check.names = FALSE
        )),
        "PE (Beads)" = list(data = data.frame(
            FITC = stats::rnorm(80, 25, 5),
            PE = stats::rnorm(80, 1100, 75),
            check.names = FALSE
        ))
    )
    qc_summary <- data.frame(
        sample = c("FITC (Beads)", "PE (Beads)"),
        fluorophore = c("FITC", "PE"),
        type = c("beads", "beads"),
        stringsAsFactors = FALSE
    )
    output_pdf <- tempfile(fileext = ".pdf")

    grDevices::pdf(output_pdf, width = 11, height = 8.5)
    on.exit(if (grDevices::dev.cur() > 1) grDevices::dev.off(), add = TRUE)
    spectreasy:::.draw_scc_post_unmix_qc_pages(
        unmixed_list = unmixed_list,
        qc_summary = qc_summary,
        markers = markers
    )
    grDevices::dev.off()

    pdf_text <- paste(pdftools::pdf_text(output_pdf), collapse = "\n")
    plain_text <- gsub("[^[:alnum:] ]+", " ", pdf_text)
    expect_match(plain_text, "Post unmixing control QC")
    expect_match(plain_text, "Pairwise false positive rate")
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
