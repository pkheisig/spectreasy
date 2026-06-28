SPECTREASY_DEFAULT_SLOW_TEST_COUNT <- 13L

run_slow_tests <- function() {
    tolower(Sys.getenv("SPECTREASY_RUN_SLOW_TESTS", "false")) %in% c("1", "true", "yes")
}

skip_slow_tests <- function(reason = "slow integration/report test") {
    # The default developer/testthat run intentionally skips the slow
    # integration/report tests. Set SPECTREASY_RUN_SLOW_TESTS=true when those
    # checks are needed; test-slow-test-policy.R locks the default skip count.
    testthat::skip_if_not(
        run_slow_tests(),
        paste0(reason, "; set SPECTREASY_RUN_SLOW_TESTS=true to run")
    )
}
