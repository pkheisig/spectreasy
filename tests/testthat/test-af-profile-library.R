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

test_that("saved AF profiles retain compatible raw and scatter-matched background data", {
    with_local_af_profile_dir({
        profile <- matrix(
            c(1, 0.4, 0.1, 0.8, 0.3, 0.2),
            nrow = 2,
            byrow = TRUE,
            dimnames = list(c("AF", "AF_2"), c("B1-A", "YG1-A", "R1-A"))
        )
        background <- list(
            method = "scatter_knn",
            detector_names = colnames(profile),
            scatter_names = c("FSC-A", "SSC-A"),
            scatter = matrix(seq_len(20), ncol = 2),
            spectra = matrix(seq_len(30), ncol = 3),
            n_events = 10L
        )
        raw_median <- c("B1-A" = 12, "YG1-A" = 8, "R1-A" = 3)
        object <- spectreasy:::.new_af_profile_object(
            profile = profile,
            raw_median = raw_median,
            scc_background = background
        )

        spectreasy::save_af_profile("AF_with_background", object)
        loaded <- spectreasy::load_af_profile("AF_with_background")

        expect_equal(loaded$raw_median, raw_median)
        expect_equal(loaded$scc_background, background)
        expect_equal(loaded$profile_version, 2L)
    })
})

test_that("extract_af_profile can derive a standalone profile from an unstained FCS file", {
    src_unstained <- testthat::test_path("../../scc/Unstained (Cells).fcs")
    testthat::skip_if_not(file.exists(src_unstained))

    afp <- spectreasy::extract_af_profile(
        src_unstained,
        af_n_bands = 1,
        seed = 1,
        show_plot = FALSE,
        verbose = FALSE
    )

    expect_s3_class(afp, "spectreasy_af_profile")
    expect_equal(rownames(afp$profile), "AF")
    expect_gt(ncol(afp$profile), 1)
    expect_s3_class(afp$plot, "ggplot")
    expect_named(afp$raw_median, colnames(afp$profile))
    expect_identical(afp$scc_background$detector_names, colnames(afp$profile))
    expect_equal(ncol(afp$scc_background$spectra), ncol(afp$profile))
    expect_equal(nrow(afp$scc_background$scatter), nrow(afp$scc_background$spectra))
})

test_that("AF profile extraction drops undersized spectral clusters", {
    common <- matrix(rep(c(1, 0.45, 0.1), 980), ncol = 3, byrow = TRUE)
    rare <- matrix(rep(c(0.1, 0.45, 1), 20), ncol = 3, byrow = TRUE)
    shapes <- rbind(common, rare)
    colnames(shapes) <- c("UV7-A", "UV8-A", "V7-A")

    result <- spectreasy:::.reference_kmeans_af_centers(
        af_shape = shapes,
        n_centers = 2,
        min_cluster_events = 50,
        min_cluster_proportion = 0.05
    )

    expect_equal(nrow(result$centers), 1)
    expect_equal(result$cluster_sizes, 980L)
    expect_equal(result$min_cluster_size, 50L)
    expect_equal(unname(which.max(result$centers[1, ])), 1L)
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
