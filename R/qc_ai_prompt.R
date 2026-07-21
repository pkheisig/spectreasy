.ai_qc_prompt_template <- function() {
    path <- system.file("extdata", "qc_ai_prompt_template.txt", package = "spectreasy")
    if (!nzchar(path)) path <- file.path("inst", "extdata", "qc_ai_prompt_template.txt")
    paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

#' Build a paste-ready AI QC evaluation prompt
#'
#' @param x A `spectreasy_ai_qc` object or arguments accepted by
#'   [collect_ai_qc()].
#' @param detail Prompt detail level.
#' @param context Optional analyst instructions included outside the canonical
#'   data delimiter.
#' @param ... Arguments forwarded to [collect_ai_qc()] when `x` is not already
#'   a canonical object.
#' @return A character scalar containing the local paste-ready prompt.
#' @export
build_ai_qc_prompt <- function(x = NULL, detail = c("standard", "compact", "full"), context = NULL, ...) {
    detail <- .match_arg_ci(detail, c("compact", "standard", "full"), "detail")
    if (!inherits(x, "spectreasy_ai_qc")) x <- do.call(collect_ai_qc, c(list(controls = x), list(...)))
    data_text <- if (identical(detail, "full")) .ai_qc_to_json(x, pretty = TRUE) else .render_ai_qc_text(x, detail = detail, markdown = FALSE)
    template <- .ai_qc_prompt_template()
    template <- gsub("{{SCOPE}}", x$scope$value, template, fixed = TRUE)
    template <- gsub("{{PRIVACY}}", x$privacy$mode, template, fixed = TRUE)
    template <- gsub("{{CONTEXT}}", context %||% "No additional analyst context supplied.", template, fixed = TRUE)
    template <- gsub("{{QC_DATA}}", data_text, template, fixed = TRUE)
    template
}
