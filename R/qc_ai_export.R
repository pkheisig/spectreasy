.ai_qc_to_json <- function(x, pretty = TRUE) {
    ordered <- .ai_qc_stable_order(unclass(x))
    jsonlite::toJSON(
        ordered, auto_unbox = TRUE, dataframe = "rows", matrix = "rowmajor",
        null = "null", na = "null", digits = NA, pretty = isTRUE(pretty),
        POSIXt = "ISO8601", force = TRUE
    )
}

.ai_qc_hash_payload <- function(x) {
    payload <- unclass(x)
    payload$provenance$generated_at <- NULL
    payload$provenance$content_sha256 <- NULL
    payload$provenance$source_artifacts <- lapply(payload$provenance$source_artifacts %||% list(), function(item) {
        item$path <- NULL
        item
    })
    payload$references$artifact_paths <- NULL
    jsonlite::toJSON(.ai_qc_stable_order(payload), auto_unbox = TRUE, dataframe = "rows", matrix = "rowmajor", null = "null", na = "null", digits = NA, pretty = FALSE, force = TRUE)
}

.ai_qc_sha256_text <- function(text) paste0(openssl::sha256(charToRaw(enc2utf8(paste0(text, collapse = "")))))

.ai_qc_file_hash <- function(path) {
    con <- file(path, "rb")
    on.exit(close(con), add = TRUE)
    bytes <- readBin(con, what = "raw", n = file.info(path)$size)
    paste0(openssl::sha256(bytes))
}

.ai_qc_resolve_path <- function(path, overwrite) {
    if (!file.exists(path) || identical(overwrite, "overwrite")) return(path)
    if (identical(overwrite, "error")) stop("AI-QC artifact already exists: ", path, call. = FALSE)
    stem <- tools::file_path_sans_ext(path)
    ext <- tools::file_ext(path)
    i <- 2L
    repeat {
        candidate <- paste0(stem, "_", i, if (nzchar(ext)) paste0(".", ext) else "")
        if (!file.exists(candidate)) return(candidate)
        i <- i + 1L
    }
}

.ai_qc_write_text <- function(text, path) {
    con <- file(path, "wb")
    on.exit(close(con), add = TRUE)
    writeBin(charToRaw(enc2utf8(text)), con)
    path
}

.ai_qc_bundle_stem <- function(scope) switch(scope, control = "spectreasy_ai_qc_controls", sample = "spectreasy_ai_qc_samples", combined = "spectreasy_ai_qc_combined")

#' Export deterministic local AI-ready QC artifacts
#'
#' @param x Optional canonical object. When omitted, one is built with
#'   [collect_ai_qc()].
#' @param output_dir Artifact root directory.
#' @param scope Export scope.
#' @param detail Human-readable and prompt detail.
#' @param privacy Privacy mode.
#' @param reference Reference-profile selection.
#' @param formats Any of `"json"`, `"txt"`, `"md"`, and `"prompt"`.
#' @param prompt Whether to write the paste-ready prompt.
#' @param context Optional prompt context.
#' @param overwrite Collision policy: `"version"`, `"overwrite"`, or
#'   `"error"`.
#' @param generated_at Generation time; tests may supply a fixed value.
#' @param ... Arguments forwarded to [collect_ai_qc()].
#' @return Invisibly returns the canonical object, artifact paths, grade counts,
#'   warnings, schema/profile details, and missing sections.
#' @export
export_ai_qc <- function(
    x = NULL, output_dir = file.path("spectreasy_outputs", "ai_qc"),
    scope = c("auto", "control", "sample", "combined"),
    detail = c("standard", "compact", "full"),
    privacy = c("standard", "strict", "none"), reference = "auto",
    formats = c("json", "txt", "md"), prompt = TRUE, context = NULL,
    overwrite = c("version", "overwrite", "error"),
    generated_at = Sys.time(), ...
) {
    scope <- .match_arg_ci(scope, c("auto", "control", "sample", "combined"), "scope")
    detail <- .match_arg_ci(detail, c("compact", "standard", "full"), "detail")
    privacy <- .match_arg_ci(privacy, c("standard", "strict", "none"), "privacy")
    overwrite <- .match_arg_ci(overwrite, c("version", "overwrite", "error"), "overwrite")
    formats <- unique(tolower(as.character(formats)))
    if (isTRUE(prompt)) formats <- unique(c(formats, "prompt"))
    unsupported <- setdiff(formats, c("json", "txt", "md", "prompt"))
    if (length(unsupported)) stop("Unsupported AI-QC formats: ", paste(unsupported, collapse = ", "), call. = FALSE)
    if (!inherits(x, "spectreasy_ai_qc")) {
        args <- list(...)
        if (is.null(args$controls)) args$controls <- x
        args$scope <- scope
        args$privacy <- privacy
        args$reference <- reference
        args$generated_at <- generated_at
        x <- do.call(collect_ai_qc, args)
    }
    scope <- x$scope$value
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(output_dir)) stop("Could not create AI-QC output directory: ", output_dir, call. = FALSE)
    stem <- .ai_qc_bundle_stem(scope)
    suffixes <- c(json = ".json", txt = ".txt", md = ".md", prompt = "_prompt.txt")
    paths <- stats::setNames(vapply(formats, function(format) .ai_qc_resolve_path(file.path(output_dir, paste0(stem, suffixes[[format]])), overwrite), character(1)), formats)
    content_hash <- .ai_qc_sha256_text(.ai_qc_hash_payload(x))
    x$provenance$content_sha256 <- content_hash
    class(x) <- c("spectreasy_ai_qc", "list")
    content <- list(
        json = .ai_qc_to_json(x, pretty = TRUE),
        txt = .render_ai_qc_text(x, detail = detail, markdown = FALSE),
        md = .render_ai_qc_text(x, detail = detail, markdown = TRUE),
        prompt = build_ai_qc_prompt(x, detail = detail, context = context)
    )
    for (format in formats) .ai_qc_write_text(content[[format]], paths[[format]])
    manifest_path <- .ai_qc_resolve_path(file.path(output_dir, paste0(stem, "_manifest.json")), overwrite)
    artifact_rows <- lapply(names(paths), function(format) list(
        format = format, path = basename(paths[[format]]), sha256 = .ai_qc_file_hash(paths[[format]]), bytes = unname(file.info(paths[[format]])$size)
    ))
    manifest <- list(
        schema = x$schema, generated_at = x$provenance$generated_at,
        content_sha256 = content_hash, package_version = x$provenance$package_version,
        git_commit = x$provenance$git_commit,
        template_version = "1.0.0", rule_version = .ai_qc_metric_version,
        profile = x$quality_reference$profile,
        privacy = x$privacy, scope = scope, detail = detail,
        source_artifacts = x$provenance$source_artifacts,
        artifacts = artifact_rows
    )
    .ai_qc_write_text(jsonlite::toJSON(.ai_qc_stable_order(manifest), auto_unbox = TRUE, null = "null", na = "null", digits = NA, pretty = TRUE), manifest_path)
    paths <- c(paths, manifest = manifest_path)
    missing <- unique(vapply(Filter(function(m) !isTRUE(m$availability), .ai_qc_walk_metrics(x)), `[[`, character(1), "stage"))
    invisible(list(
        object = x, paths = paths, grade_counts = x$grade_summary$counts,
        warnings = character(), schema = x$schema,
        profile = x$quality_reference$profile,
        missing_sections = missing,
        content_sha256 = content_hash,
        prompt_size = list(characters = nchar(content$prompt), estimated_tokens = ceiling(nchar(content$prompt) / 4))
    ))
}
