.ai_qc_alias_map <- function(values, prefix) {
    values <- sort(unique(as.character(values[!is.na(values) & nzchar(as.character(values))])))
    stats::setNames(sprintf("%s_%03d", prefix, seq_along(values)), values)
}

.ai_qc_replace_values <- function(x, aliases) {
    if (!length(aliases)) return(x)
    visit <- function(value) {
        if (is.data.frame(value)) {
            value[] <- lapply(value, visit)
            return(value)
        }
        if (is.list(value)) return(lapply(value, visit))
        if (is.character(value)) {
            hit <- match(value, names(aliases))
            value[!is.na(hit)] <- unname(aliases[hit[!is.na(hit)]])
            for (source in names(aliases)[order(nchar(names(aliases)), decreasing = TRUE)]) {
                prefix <- paste0(source, ":")
                prefixed <- startsWith(value, prefix)
                value[prefixed] <- paste0(aliases[[source]], substring(value[prefixed], nchar(source) + 1L))
            }
        }
        value
    }
    visit(x)
}

.ai_qc_strip_paths <- function(x, strict = FALSE) {
    visit <- function(value, key = "") {
        if (is.data.frame(value)) {
            value[] <- lapply(seq_along(value), function(i) visit(value[[i]], names(value)[i]))
            return(value)
        }
        if (!is.list(value)) {
            if (is.character(value) && grepl("(^|_)(path|paths|file|files|directory|artifact|artifacts)($|_)", key, ignore.case = TRUE)) {
                if (isTRUE(strict)) return(NULL)
                return(vapply(value, function(item) if (grepl("^(/|[A-Za-z]:[/\\\\])", item)) basename(item) else item, character(1)))
            }
            return(value)
        }
        nms <- names(value)
        for (i in seq_along(value)) {
            child_key <- if (!is.null(nms)) nms[i] else key
            if (isTRUE(strict) && grepl("operator|date|comment|note|free.?text", child_key, ignore.case = TRUE)) {
                value[i] <- list(NULL)
            } else {
                value[i] <- list(visit(value[[i]], child_key))
            }
        }
        value
    }
    visit(x)
}

.redact_ai_qc <- function(x, mode = c("standard", "strict", "none")) {
    mode <- .match_arg_ci(mode, c("standard", "strict", "none"), "mode")
    if (identical(mode, "none")) {
        x$privacy$redactions <- list(sample_aliases = 0L, control_aliases = 0L, absolute_paths_removed = FALSE)
        return(x)
    }
    samples <- vapply(x$samples$entities %||% list(), function(entity) as.character(entity$id %||% ""), character(1))
    sample_aliases <- .ai_qc_alias_map(samples, "sample")
    x <- .ai_qc_replace_values(x, sample_aliases)
    control_aliases <- character()
    if (identical(mode, "strict")) {
        control_names <- vapply(x$controls$entities %||% list(), function(entity) as.character(entity$sample %||% ""), character(1))
        control_aliases <- .ai_qc_alias_map(control_names, "control")
        x <- .ai_qc_replace_values(x, control_aliases)
    }
    x <- .ai_qc_strip_paths(x, strict = identical(mode, "strict"))
    x$privacy$mode <- mode
    x$privacy$redactions <- list(
        sample_aliases = length(sample_aliases),
        control_aliases = length(control_aliases),
        absolute_paths_removed = TRUE,
        free_text_removed = identical(mode, "strict")
    )
    class(x) <- c("spectreasy_ai_qc", "list")
    x
}
