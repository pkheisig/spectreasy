testthat::test_that("create_control_file recognizes fluor and control type from filename variants", {
    testthat::skip_if_not_installed("spectreasy")

    scc_dir <- tempfile("spectreasy_scc_")
    dir.create(scc_dir, recursive = TRUE, showWarnings = FALSE)

    files <- c(
        "LIVE DEAD NIR (Cells).fcs",
        "Alexa 647 (Beads).fcs",
        "AF594 (Beads).fcs",
        "AF488 CD4 (Beads).fcs",
        "BB515 CD123 (Beads).fcs",
        "BB700 HLA-DR (Beads).fcs",
        "APC-H7 CD279 (Beads).fcs",
        "APC PD-1 (Beads).fcs",
        "PE PD1 (Beads).fcs",
        "FITC CD39 (Beads).fcs",
        "FITC ENTPD1 (Beads).fcs",
        "BV510 CD1a (Beads).fcs",
        "BV421 CD371 (Beads).fcs",
        "FVD780 CD8 (Cells).fcs",
        "Zombie Aqua CCR7 (Cells).fcs",
        "PE-CF594 (Beads).fcs",
        "PECF594 (beads).fcs",
        "pe cy7 (bEaDs).fcs",
        "PE-Fire 700 (BEADS).fcs",
        "PE (Beads).fcs",
        "Unstained (Cells).fcs",
        "Unstained Dead (Cells).fcs",
        "Dead Control (Cells).fcs",
        "Live Dead Control (Cells).fcs",
        "Unstained (Beads).fcs",
        "US Beads.fcs",
        "Negative Beads.fcs",
        "Bead background.fcs",
        "BG CompBeads.fcs"
    )
    created <- file.create(file.path(scc_dir, files))
    testthat::expect_true(all(created))

    out_csv <- tempfile(fileext = ".csv")
    df <- spectreasy::create_control_file(
        input_folder = scc_dir,
        output_file = out_csv
    )

    by_file <- split(df, df$filename)
    testthat::expect_false("large.gate" %in% colnames(df))

    testthat::expect_equal(by_file[["LIVE DEAD NIR (Cells).fcs"]]$fluorophore[[1]], "LIVE/DEAD NIR")
    testthat::expect_equal(by_file[["LIVE DEAD NIR (Cells).fcs"]]$control.type[[1]], "cells")
    testthat::expect_equal(by_file[["LIVE DEAD NIR (Cells).fcs"]]$is.viability[[1]], "TRUE")

    testthat::expect_equal(by_file[["Alexa 647 (Beads).fcs"]]$fluorophore[[1]], "Alexa Fluor 647")
    testthat::expect_equal(by_file[["AF594 (Beads).fcs"]]$fluorophore[[1]], "Alexa Fluor 594")

    testthat::expect_equal(by_file[["AF488 CD4 (Beads).fcs"]]$fluorophore[[1]], "Alexa Fluor 488")
    testthat::expect_equal(by_file[["AF488 CD4 (Beads).fcs"]]$marker[[1]], "CD4")
    testthat::expect_equal(by_file[["AF488 CD4 (Beads).fcs"]]$control.type[[1]], "beads")

    testthat::expect_equal(by_file[["BB515 CD123 (Beads).fcs"]]$fluorophore[[1]], "BB515")
    testthat::expect_equal(by_file[["BB515 CD123 (Beads).fcs"]]$marker[[1]], "CD123")

    testthat::expect_equal(by_file[["BB700 HLA-DR (Beads).fcs"]]$fluorophore[[1]], "BB700")
    testthat::expect_equal(by_file[["BB700 HLA-DR (Beads).fcs"]]$marker[[1]], "HLA-DR")

    testthat::expect_equal(by_file[["APC-H7 CD279 (Beads).fcs"]]$fluorophore[[1]], "APC-H7")
    testthat::expect_equal(by_file[["APC-H7 CD279 (Beads).fcs"]]$marker[[1]], "CD279")

    testthat::expect_equal(by_file[["APC PD-1 (Beads).fcs"]]$marker[[1]], "PD-1")
    testthat::expect_equal(by_file[["PE PD1 (Beads).fcs"]]$marker[[1]], "PD1")
    testthat::expect_equal(by_file[["FITC CD39 (Beads).fcs"]]$marker[[1]], "CD39")
    testthat::expect_equal(by_file[["FITC ENTPD1 (Beads).fcs"]]$marker[[1]], "ENTPD1")
    testthat::expect_equal(by_file[["BV510 CD1a (Beads).fcs"]]$marker[[1]], "CD1a")

    testthat::expect_equal(by_file[["BV421 CD371 (Beads).fcs"]]$marker[[1]], "CD371")

    testthat::expect_equal(by_file[["FVD780 CD8 (Cells).fcs"]]$fluorophore[[1]], "FVD780")
    testthat::expect_equal(by_file[["FVD780 CD8 (Cells).fcs"]]$marker[[1]], "CD8")
    testthat::expect_equal(by_file[["FVD780 CD8 (Cells).fcs"]]$is.viability[[1]], "TRUE")

    testthat::expect_equal(by_file[["Zombie Aqua CCR7 (Cells).fcs"]]$fluorophore[[1]], "Zombie Aqua")
    testthat::expect_equal(by_file[["Zombie Aqua CCR7 (Cells).fcs"]]$marker[[1]], "CCR7")
    testthat::expect_equal(by_file[["Zombie Aqua CCR7 (Cells).fcs"]]$is.viability[[1]], "TRUE")

    testthat::expect_equal(by_file[["PE-CF594 (Beads).fcs"]]$fluorophore[[1]], "PE-CF594")
    testthat::expect_equal(by_file[["PE-CF594 (Beads).fcs"]]$control.type[[1]], "beads")

    testthat::expect_equal(by_file[["PECF594 (beads).fcs"]]$fluorophore[[1]], "PE-CF594")
    testthat::expect_equal(by_file[["PECF594 (beads).fcs"]]$control.type[[1]], "beads")

    testthat::expect_equal(by_file[["pe cy7 (bEaDs).fcs"]]$fluorophore[[1]], "PE-Cy7")
    testthat::expect_equal(by_file[["pe cy7 (bEaDs).fcs"]]$control.type[[1]], "beads")

    testthat::expect_equal(by_file[["PE-Fire 700 (BEADS).fcs"]]$fluorophore[[1]], "PE-Fire 700")
    testthat::expect_equal(by_file[["PE-Fire 700 (BEADS).fcs"]]$control.type[[1]], "beads")

    testthat::expect_equal(by_file[["PE (Beads).fcs"]]$fluorophore[[1]], "PE")
    testthat::expect_equal(by_file[["PE (Beads).fcs"]]$control.type[[1]], "beads")
    testthat::expect_equal(by_file[["PE (Beads).fcs"]]$is.viability[[1]], "")

    testthat::expect_equal(by_file[["Unstained (Cells).fcs"]]$fluorophore[[1]], "AF")
    testthat::expect_equal(by_file[["Unstained (Cells).fcs"]]$marker[[1]], "Autofluorescence")
    testthat::expect_equal(by_file[["Unstained (Cells).fcs"]]$control.type[[1]], "cells")
    testthat::expect_equal(by_file[["Unstained (Cells).fcs"]]$is.dead[[1]], "")

    testthat::expect_equal(by_file[["Unstained Dead (Cells).fcs"]]$fluorophore[[1]], "AF_Dead")
    testthat::expect_equal(by_file[["Unstained Dead (Cells).fcs"]]$marker[[1]], "Viability dead background")
    testthat::expect_equal(by_file[["Unstained Dead (Cells).fcs"]]$control.type[[1]], "cells")
    testthat::expect_equal(by_file[["Unstained Dead (Cells).fcs"]]$is.viability[[1]], "")
    testthat::expect_equal(by_file[["Unstained Dead (Cells).fcs"]]$is.dead[[1]], "TRUE")

    testthat::expect_equal(by_file[["Dead Control (Cells).fcs"]]$fluorophore[[1]], "AF_Dead")
    testthat::expect_equal(by_file[["Dead Control (Cells).fcs"]]$marker[[1]], "Viability dead background")
    testthat::expect_equal(by_file[["Dead Control (Cells).fcs"]]$is.dead[[1]], "TRUE")

    testthat::expect_equal(by_file[["Live Dead Control (Cells).fcs"]]$is.viability[[1]], "TRUE")
    testthat::expect_equal(by_file[["Live Dead Control (Cells).fcs"]]$is.dead[[1]], "")

    bead_negative_files <- c(
        "Unstained (Beads).fcs",
        "US Beads.fcs",
        "Negative Beads.fcs",
        "Bead background.fcs",
        "BG CompBeads.fcs"
    )
    for (fn in bead_negative_files) {
        testthat::expect_equal(by_file[[fn]]$fluorophore[[1]], "AF_Bead", info = fn)
        testthat::expect_equal(by_file[[fn]]$marker[[1]], "Bead background", info = fn)
        testthat::expect_equal(by_file[[fn]]$control.type[[1]], "beads", info = fn)
        testthat::expect_equal(by_file[[fn]]$is.viability[[1]], "", info = fn)
    }
})

testthat::test_that("create_control_file maps multiple unstained SCC files as AF sources", {
    scc_dir <- tempfile("spectreasy_multi_af_scc_")
    dir.create(scc_dir, recursive = TRUE, showWarnings = FALSE)
    files <- c(
        "Unstained lymphocytes (Cells).fcs",
        "Unstained myeloid (Cells).fcs",
        "Unstained tumor (Cells).fcs",
        "FITC CD4 (Beads).fcs"
    )
    testthat::expect_true(all(file.create(file.path(scc_dir, files))))

    df <- spectreasy::create_control_file(
        input_folder = scc_dir,
        output_file = tempfile(fileext = ".csv")
    )

    by_file <- split(df, df$filename)
    testthat::expect_equal(by_file[["Unstained lymphocytes (Cells).fcs"]]$fluorophore[[1]], "AF")
    testthat::expect_equal(by_file[["Unstained myeloid (Cells).fcs"]]$fluorophore[[1]], "AF_2")
    testthat::expect_equal(by_file[["Unstained tumor (Cells).fcs"]]$fluorophore[[1]], "AF_3")
    testthat::expect_true(all(vapply(
        by_file[files[1:3]],
        function(x) identical(x$marker[[1]], "Autofluorescence") && identical(x$control.type[[1]], "cells"),
        logical(1)
    )))
})

testthat::test_that("detector fallback matches whole detector codes", {
    testthat::skip_if_not_installed("spectreasy")

    infer <- spectreasy:::.infer_fluor_from_detector

    testthat::expect_equal(infer("V1-A", fluor_channel_map = character()), "BV421")
    testthat::expect_equal(infer("V10-A", fluor_channel_map = character()), "BV605")
    testthat::expect_equal(infer("V11-A", fluor_channel_map = character()), "BV650")
    testthat::expect_equal(infer("V13-A", fluor_channel_map = character()), "BV711")
    testthat::expect_equal(infer("UV10-A", fluor_channel_map = character()), "BUV615")
    testthat::expect_equal(infer("UV11-A", fluor_channel_map = character()), "BUV661")
    testthat::expect_equal(infer("UV17-A", fluor_channel_map = character()), "BUV395")
    testthat::expect_equal(infer("R10-A", fluor_channel_map = character()), "APC")
})

testthat::test_that("create_control_file warns when peak detection cannot read an FCS", {
    scc_dir <- tempfile("spectreasy_bad_scc_")
    dir.create(scc_dir, recursive = TRUE, showWarnings = FALSE)
    writeBin(as.raw(rep(0, 128)), file.path(scc_dir, "BV510 (Cells).fcs"))

    out_csv <- tempfile(fileext = ".csv")
    testthat::expect_warning(
        df <- spectreasy::create_control_file(
            input_folder = scc_dir,
            output_file = out_csv
        ),
        regexp = "Could not read FCS file while auto-detecting peak channel"
    )

    testthat::expect_equal(df$fluorophore[[1]], "BV510")
    testthat::expect_equal(df$channel[[1]], "V7-A")
})

testthat::test_that("supported cytometer metadata is normalized", {
    ids <- spectreasy::supported_cytometers()

    testthat::expect_true(all(c(
        "aurora", "northern_lights", "id7000", "discover_s8",
        "discover_a8", "a5se", "opteon", "mosaic", "xenith"
    ) %in% ids))
    testthat::expect_true("auto" %in% spectreasy::supported_cytometers(include_auto = TRUE))

    channels <- spectreasy:::.read_fluorophore_channel_dictionary()
    testthat::expect_true(any(
        channels$cytometer == "xenith" &
            channels$fluorophore == "Alexa Fluor 488" &
            channels$channel == "FL37-A"
    ))
    testthat::expect_true(any(
        channels$cytometer == "aurora" &
            channels$fluorophore == "7-AAD" &
            channels$channel == "YG-4"
    ))

    xenith_ref <- spectreasy:::.load_control_file_shipped_reference("Xenith")
    testthat::expect_equal(xenith_ref$channel_map[["FL37-A"]], "FITC")
    aurora_ref <- spectreasy:::.load_control_file_shipped_reference("Aurora")
    testthat::expect_equal(aurora_ref$channel_map[["YG4-A"]], "PE-Fire 640")
})

testthat::test_that("cytometer auto detection recognizes detector naming conventions", {
    xenith_pd <- data.frame(
        name = c("FL07-A", "FL08-A", "FL37-A", "FL36-A", "FSC51-A", "SSC52-A", "Time"),
        desc = c(
            "349nm - 387/11-A", "349nm - 420/10-A", "488nm - 530/30-A",
            "488nm - 515/20-A", "FSC51-A", "SSC52-A", "Time"
        ),
        stringsAsFactors = FALSE
    )
    testthat::expect_equal(spectreasy:::.infer_cytometer_from_pd(xenith_pd), "xenith")

    discover_pd <- data.frame(
        name = c("UV1 (375)-A", "UV2 (390)-A", "B2 (515)-A", "YG5 (655)-A", "FSC-A", "SSC (Violet)-A"),
        desc = c("UV1 (375)-A", "UV2 (390)-A", "B2 (515)-A", "YG5 (655)-A", "FSC-A", "SSC (Violet)-A"),
        stringsAsFactors = FALSE
    )
    testthat::expect_true(spectreasy:::.infer_cytometer_from_pd(discover_pd) %in% c("discover_s8", "discover_a8"))
})

testthat::test_that("custom fluorophore overrides accept common filename forms", {
    testthat::skip_if_not_installed("spectreasy")

    scc_dir <- tempfile("spectreasy_scc_")
    dir.create(scc_dir, recursive = TRUE, showWarnings = FALSE)
    testthat::expect_true(file.create(file.path(scc_dir, "odd-control-name.fcs")))

    out_csv <- tempfile(fileext = ".csv")
    df <- spectreasy::create_control_file(
        input_folder = scc_dir,
        unknown_fluor_policy = "empty",
        output_file = out_csv,
        custom_fluorophores = c("odd-control-name" = "BUV737")
    )

    testthat::expect_equal(df$fluorophore[[1]], "BUV737")
})

testthat::test_that("create_control_file keeps dictionary peak channel for known fluorophores", {
    testthat::skip_if_not_installed("spectreasy")

    scc_dir <- tempfile("spectreasy_scc_")
    dir.create(scc_dir, recursive = TRUE, showWarnings = FALSE)

    n <- 300
    exprs <- cbind(
        "V3-A" = c(stats::rnorm(n - 10, 200, 20), stats::rnorm(10, 80000, 500)),
        "V15-A" = stats::rnorm(n, 6000, 300),
        "FSC-A" = stats::rnorm(n, 90000, 7000),
        "SSC-A" = stats::rnorm(n, 45000, 5000),
        Time = seq_len(n)
    )
    flowCore::write.FCS(
        flowCore::flowFrame(exprs),
        file.path(scc_dir, "CD274 BV786 (Beads).fcs")
    )
    flowCore::write.FCS(
        flowCore::flowFrame(exprs),
        file.path(scc_dir, "Mystery Dye (Beads).fcs")
    )

    df <- spectreasy::create_control_file(
        input_folder = scc_dir,
        cytometer = "aurora",
        unknown_fluor_policy = "empty",
        custom_fluorophores = c("Mystery Dye (Beads)" = "Mystery Dye"),
        output_file = tempfile(fileext = ".csv")
    )
    by_file <- split(df, df$filename)

    testthat::expect_equal(by_file[["CD274 BV786 (Beads).fcs"]]$fluorophore[[1]], "BV786")
    testthat::expect_equal(by_file[["CD274 BV786 (Beads).fcs"]]$channel[[1]], "V15-A")
    testthat::expect_equal(by_file[["Mystery Dye (Beads).fcs"]]$fluorophore[[1]], "Mystery Dye")
    testthat::expect_equal(by_file[["Mystery Dye (Beads).fcs"]]$channel[[1]], "V3-A")
})
