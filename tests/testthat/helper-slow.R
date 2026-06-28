run_slow_tests <- function() {
    tolower(Sys.getenv("SPECTREASY_RUN_SLOW_TESTS", "false")) %in% c("1", "true", "yes")
}

skip_slow_tests <- function(reason = "slow integration/report test") {
    testthat::skip_if_not(
        run_slow_tests(),
        paste0(reason, "; set SPECTREASY_RUN_SLOW_TESTS=true to run")
    )
}
