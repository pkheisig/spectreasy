test_that("spectreasy_example_data downloads, unzips, and caches local archives", {
    skip_if_not(nzchar(Sys.which("zip")))

    src_root <- tempfile("spectreasy_example_src_")
    dir.create(src_root)
    dir.create(file.path(src_root, "sample"))
    dir.create(file.path(src_root, "scc"))
    writeBin(charToRaw("sample"), file.path(src_root, "sample", "sample.fcs"))
    writeBin(charToRaw("control"), file.path(src_root, "scc", "FITC (Beads).fcs"))

    zip_file <- tempfile(fileext = ".zip")
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(src_root)
    utils::zip(zipfile = zip_file, files = c("sample", "scc"), flags = "-r9Xq")

    cache_dir <- tempfile("spectreasy_example_cache_")
    paths <- spectreasy::spectreasy_example_data(asset = zip_file, cache_dir = cache_dir)

    expect_true(dir.exists(paths$sample_dir))
    expect_true(dir.exists(paths$scc_dir))
    expect_true(file.exists(file.path(paths$sample_dir, "sample.fcs")))
    expect_true(file.exists(file.path(paths$scc_dir, "FITC (Beads).fcs")))
    expect_true(file.exists(paths$zip_file))
    expect_match(paths$zip_file, "\\.zip$")
    expect_true(length(paths$sample_files) == 1)

    file.remove(paths$zip_file)
    cached_paths <- spectreasy::spectreasy_example_data(asset = zip_file, cache_dir = cache_dir)
    expect_equal(normalizePath(cached_paths$sample_dir), normalizePath(paths$sample_dir))
    expect_true(file.exists(file.path(cached_paths$scc_dir, "FITC (Beads).fcs")))
})

test_that("spectreasy_example_data validates archive contents", {
    skip_if_not(nzchar(Sys.which("zip")))

    src_root <- tempfile("spectreasy_example_bad_src_")
    dir.create(src_root)
    dir.create(file.path(src_root, "other"))
    writeBin(charToRaw("x"), file.path(src_root, "other", "file.txt"))

    zip_file <- tempfile(fileext = ".zip")
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(src_root)
    utils::zip(zipfile = zip_file, files = "other", flags = "-r9Xq")

    expect_error(
        spectreasy::spectreasy_example_data(asset = zip_file, cache_dir = tempfile("spectreasy_example_bad_cache_")),
        regexp = "expected 'sample/' and 'scc/' folders"
    )
})
