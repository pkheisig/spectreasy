with_metadata_af_profile_dir <- function(code) {
    old <- getOption("spectreasy.af_profile_dir")
    dir_path <- tempfile("spectreasy_af_metadata_")
    dir.create(dir_path, recursive = TRUE)
    options(spectreasy.af_profile_dir = dir_path)
    on.exit({
        options(spectreasy.af_profile_dir = old)
        unlink(dir_path, recursive = TRUE, force = TRUE)
    }, add = TRUE)
    force(code)
}

metadata_test_profile <- function(values = c(1, 0.3, 0.1, 0.2, 1, 0.4)) {
    matrix(
        values,
        nrow = 2,
        byrow = TRUE,
        dimnames = list(c("AF", "AF_2"), c("B1-A", "YG1-A", "R1-A"))
    )
}

test_that("AF profile metadata survives save, load, rename, and overwrite", {
    with_metadata_af_profile_dir({
        first_metadata <- list(
            cytometer = "aurora",
            acquisition_date = "2026-07-18",
            tissue = "Peripheral blood",
            sample_type = "PBMC",
            preprocessing = "Cryopreserved, thawed and washed"
        )
        first <- spectreasy:::.new_af_profile_object(
            profile = metadata_test_profile(),
            extraction = list(fcs_file = "/data/unstained_cells.fcs"),
            metadata = first_metadata,
            created = as.POSIXct("2026-07-18 08:42:00", tz = "UTC")
        )
        spectreasy::save_af_profile("PBMC_AF", first)

        loaded <- spectreasy::load_af_profile("PBMC_AF")
        expect_identical(loaded$schema_version, 1L)
        expect_equal(loaded$spectra, loaded$profile)
        expect_identical(loaded$detectors, colnames(loaded$profile))
        expect_identical(loaded$metadata$cytometer, "aurora")
        expect_identical(loaded$metadata$source_fcs, "unstained_cells.fcs")
        expect_identical(loaded$metadata$band_count, 2L)

        spectreasy::rename_af_profile("PBMC_AF", "PBMC_AF_renamed")
        renamed <- spectreasy::load_af_profile("PBMC_AF_renamed")
        expect_equal(renamed$metadata, loaded$metadata)

        replacement <- spectreasy:::.new_af_profile_object(
            profile = metadata_test_profile(c(1, 0, 0, 0, 1, 0)),
            metadata = list(cytometer = "northern-lights", tissue = "Tumor")
        )
        spectreasy::save_af_profile("PBMC_AF_renamed", replacement, overwrite = TRUE)
        overwritten <- spectreasy::load_af_profile("PBMC_AF_renamed")
        expect_identical(overwritten$metadata$cytometer, "northern-lights")
        expect_identical(overwritten$metadata$tissue, "Tumor")
        expect_true(is.na(overwritten$metadata$sample_type))
        expect_equal(overwritten$profile, replacement$profile)
    })
})

test_that("legacy AF profiles are wrapped with inferred system fields", {
    with_metadata_af_profile_dir({
        profile <- metadata_test_profile()
        saveRDS(profile, file.path(getOption("spectreasy.af_profile_dir"), "legacy.rds"))

        loaded <- spectreasy::load_af_profile("legacy")
        expect_s3_class(loaded, "spectreasy_af_profile")
        expect_identical(loaded$schema_version, 1L)
        expect_equal(loaded$profile, profile)
        expect_identical(loaded$detectors, colnames(profile))
        expect_identical(loaded$metadata$band_count, 2L)
        expect_true(is.na(loaded$metadata$cytometer))
        expect_true(is.na(loaded$metadata$created_at))
        expect_true(is.na(loaded$metadata$spectreasy_version))

        spectreasy::rename_af_profile("legacy", "legacy_renamed")
        expect_equal(spectreasy::load_af_profile("legacy_renamed")$profile, profile)
    })
})

test_that("AF profile metadata validates dates and scalar strings", {
    profile <- metadata_test_profile()
    expect_error(
        spectreasy:::.new_af_profile_object(profile = profile, metadata = "aurora"),
        "named list"
    )
    expect_error(
        spectreasy:::.new_af_profile_object(profile = profile, metadata = list(cytometer = c("a", "b"))),
        "single character string"
    )
    expect_error(
        spectreasy:::.new_af_profile_object(profile = profile, metadata = list(tissue = 42)),
        "single character string"
    )
    expect_error(
        spectreasy:::.new_af_profile_object(profile = profile, metadata = list(acquisition_date = "18-07-2026")),
        "YYYY-MM-DD"
    )
    expect_no_error(
        spectreasy:::.new_af_profile_object(profile = profile, metadata = list(acquisition_date = "2026-07-18"))
    )
})
