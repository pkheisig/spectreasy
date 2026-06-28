test_that("default test suite documents the slow-test skip policy", {
    test_files <- list.files(
        testthat::test_path(),
        pattern = "^test-.*[.]R$",
        full.names = TRUE
    )
    test_files <- test_files[basename(test_files) != "test-slow-test-policy.R"]

    slow_skip_calls <- unlist(lapply(test_files, function(path) {
        grep("skip_slow_tests(", readLines(path, warn = FALSE), fixed = TRUE, value = TRUE)
    }), use.names = FALSE)

    expect_equal(length(slow_skip_calls), SPECTREASY_DEFAULT_SLOW_TEST_COUNT)
})
