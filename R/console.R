# Internal console-formatting helpers.

.spectreasy_console_width <- function() {
    58L
}

.spectreasy_console_line <- function(width = .spectreasy_console_width()) {
    strrep("-", width)
}

.spectreasy_console_header <- function(title, width = .spectreasy_console_width()) {
    line <- .spectreasy_console_line(width)
    message("\n", line)
    message("Spectreasy ", title)
    message(line)
}

.spectreasy_console_footer <- function(width = .spectreasy_console_width(), blank = TRUE) {
    suffix <- if (isTRUE(blank)) "\n" else ""
    message(.spectreasy_console_line(width), suffix)
}

.spectreasy_console_path <- function(path) {
    if (is.null(path) || length(path) == 0 || is.na(path[1]) || !nzchar(path[1])) {
        return("")
    }
    path <- as.character(path[1])
    path <- normalizePath(path, mustWork = FALSE)
    wd <- normalizePath(getwd(), mustWork = FALSE)
    if (identical(dirname(path), wd)) basename(path) else path
}

.spectreasy_console_field <- function(label, value = "", width = 9L) {
    label <- as.character(label)[1]
    value <- paste(as.character(value), collapse = " ")
    message(sprintf(paste0("%-", width, "s: %s"), label, value))
}

.spectreasy_console_step <- function(label, detail = NULL, width = 18L) {
    label <- as.character(label)[1]
    if (is.null(detail) || length(detail) == 0 || is.na(detail[1]) || !nzchar(as.character(detail[1]))) {
        message("  - ", label)
    } else {
        detail <- paste(as.character(detail), collapse = " ")
        message(sprintf(paste0("  - %-", width, "s %s"), label, detail))
    }
}

.spectreasy_stop_missing_directory <- function(path, label = "directory") {
    stop(label, " not found: ", path, call. = FALSE)
}

.spectreasy_stop_empty_fcs_directory <- function(path, label = "directory") {
    stop("No FCS files found in ", label, ": ", path, call. = FALSE)
}

.spectreasy_stop_missing_file <- function(path, label = "file") {
    stop(label, " not found: ", path, call. = FALSE)
}

.spectreasy_read_fcs <- function(path, label = "FCS file") {
    tryCatch(
        flowCore::read.FCS(path, transformation = FALSE, truncate_max_range = FALSE),
        error = function(e) {
            stop(
                "Could not read ", label, ": ", path,
                ". Reason: ", conditionMessage(e),
                call. = FALSE
            )
        }
    )
}
