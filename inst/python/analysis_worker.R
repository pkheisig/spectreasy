args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("usage: analysis_worker.R REQUEST_JSON RESULT_RDS STATUS_JSON", call. = FALSE)
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

request_path <- normalizePath(args[[1]], mustWork = TRUE)
result_path <- normalizePath(args[[2]], mustWork = FALSE)
status_path <- normalizePath(args[[3]], mustWork = FALSE)

write_status <- function(state, message, error = NULL) {
    value <- list(
        state = state,
        message = message,
        updated_at = as.character(Sys.time()),
        error = error
    )
    temporary <- tempfile("analysis-job-", tmpdir = dirname(status_path), fileext = ".json")
    jsonlite::write_json(value, temporary, auto_unbox = TRUE, pretty = TRUE, null = "null")
    if (!file.rename(temporary, status_path)) {
        file.copy(temporary, status_path, overwrite = TRUE)
        unlink(temporary, force = TRUE)
    }
}

request <- jsonlite::fromJSON(request_path, simplifyVector = FALSE)
repository <- as.character(request$repository %||% "")[1]
if (nzchar(repository) && file.exists(file.path(repository, "DESCRIPTION")) && requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(repository, quiet = TRUE)
} else {
    suppressPackageStartupMessages(library(spectreasy))
}

write_status("running", "Preparing gated events")
tryCatch({
    result <- spectreasy:::.spectreasy_analysis_run_request(
        as.character(request$project)[1],
        request$request
    )
    saveRDS(result, result_path)
    write_status("completed", "Analysis completed")
}, error = function(error) {
    write_status("failed", "Analysis failed", conditionMessage(error))
    quit(status = 1L)
})
