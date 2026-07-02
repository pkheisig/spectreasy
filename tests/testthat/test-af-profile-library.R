with_local_af_profile_dir <- function(code) {
    old <- getOption("spectreasy.af_profile_dir")
    dir_path <- tempfile("spectreasy_af_profiles_")
    options(spectreasy.af_profile_dir = dir_path)
    on.exit({
        options(spectreasy.af_profile_dir = old)
        unlink(dir_path, recursive = TRUE, force = TRUE)
    }, add = TRUE)
    force(code)
}

test_that("save_af_profile stores only AF rows from a full reference matrix", {
    with_local_af_profile_dir({
        M <- matrix(c(
            1.0, 0.1, 0.2,
            0.2, 1.0, 0.3,
            0.9, 0.4, 0.1,
            0.7, 0.2, 0.5
        ), nrow = 4, byrow = TRUE)
        rownames(M) <- c("FITC", "PE", "AF", "AF_2")
        colnames(M) <- c("B1-A", "YG1-A", "R1-A")

        path <- spectreasy::save_af_profile("AF_B16.F10_tumor", M)
        expect_true(file.exists(path))

        afp <- spectreasy::load_af_profile("AF_B16.F10_tumor")
        expect_s3_class(afp, "spectreasy_af_profile")
        expect_equal(rownames(afp$profile), c("AF", "AF_2"))
        expect_false("FITC" %in% rownames(afp$profile))
        expect_s3_class(afp$plot, "ggplot")

        listed <- spectreasy::list_af_profiles()
        expect_equal(listed$name, "AF_B16.F10_tumor")
        expect_equal(listed$bands, 2)
        expect_equal(listed$detectors, 3)
    })
})

test_that("add_af_profile appends loaded AF bands and replaces existing AF by default", {
    with_local_af_profile_dir({
        profile <- matrix(c(
            0.9, 0.4, 0.1,
            0.7, 0.2, 0.5
        ), nrow = 2, byrow = TRUE)
        rownames(profile) <- c("AF", "AF_2")
        colnames(profile) <- c("B1-A", "YG1-A", "R1-A")
        spectreasy::save_af_profile("AF_saved", profile)

        M <- matrix(c(
            1.0, 0.1, 0.2,
            0.2, 1.0, 0.3,
            0.1, 0.1, 0.1
        ), nrow = 3, byrow = TRUE)
        rownames(M) <- c("FITC", "PE", "AF")
        colnames(M) <- c("B1-A", "YG1-A", "R1-A")

        out <- spectreasy::add_af_profile(M, "AF_saved")
        expect_equal(rownames(out), c("FITC", "PE", "AF", "AF_2"))
        expect_equal(out["AF", ], profile["AF", ])
        expect_equal(out["AF_2", ], profile["AF_2", ])

        expect_error(
            spectreasy::add_af_profile(M, "AF_saved", replace_existing = FALSE),
            regexp = "already exist"
        )
    })
})

test_that("extract_af_profile builds the default saved k-means AF bank", {
    detector_names <- c("B1-A", "YG1-A", "V1-A")
    af_events <- rbind(
        matrix(rep(c(100, 15, 5), 120), ncol = 3, byrow = TRUE),
        matrix(rep(c(10, 90, 20), 120), ncol = 3, byrow = TRUE)
    )
    colnames(af_events) <- detector_names
    src_unstained <- tempfile(fileext = ".fcs")
    file.create(src_unstained)

    testthat::local_mocked_bindings(
        .prepare_reference_detector_info = function(fcs_file) {
            list(
                detector_names = detector_names,
                pd_meta = data.frame(name = detector_names, stringsAsFactors = FALSE)
            )
        },
        .extract_reference_af_gated_events = function(fcs_file, detector_names, config) {
            list(
                events = af_events,
                source = data.table::data.table(
                    file = basename(fcs_file),
                    path = normalizePath(fcs_file, mustWork = FALSE),
                    n_total = nrow(af_events),
                    n_scatter_gated = nrow(af_events),
                    scatter_gate_pct = 100,
                    fsc_channel = "FSC-A",
                    ssc_channel = "SSC-A"
                )
            )
        },
        .build_af_profile_plot = function(profile, pd = NULL) {
            ggplot2::ggplot()
        },
        .package = "spectreasy"
    )

    afp <- spectreasy::extract_af_profile(
        src_unstained,
        seed = 1,
        show_plot = FALSE,
        verbose = FALSE
    )

    expect_s3_class(afp, "spectreasy_af_profile")
    expect_equal(nrow(afp$profile), 2)
    expect_equal(rownames(afp$profile), c("AF", "AF_2"))
    expect_equal(afp$extraction$auto_selection$method, "kmeans_distinct")
    expect_equal(afp$extraction$auto_selection$n_bands, 2)
    expect_s3_class(afp$plot, "ggplot")
})

test_that("AF profile library rejects missing AF rows and detector mismatches", {
    with_local_af_profile_dir({
        M_no_af <- matrix(c(1, 0.2, 0.1, 1), nrow = 2, byrow = TRUE)
        rownames(M_no_af) <- c("FITC", "PE")
        colnames(M_no_af) <- c("B1-A", "YG1-A")
        expect_error(
            spectreasy::save_af_profile("AF_none", M_no_af),
            regexp = "does not contain AF rows"
        )

        profile <- matrix(c(1, 0.5), nrow = 1)
        rownames(profile) <- "AF"
        colnames(profile) <- c("B1-A", "R1-A")
        spectreasy::save_af_profile("AF_mismatch", profile)

        expect_error(
            spectreasy::add_af_profile(M_no_af, "AF_mismatch"),
            regexp = "detectors do not match"
        )
    })
})
