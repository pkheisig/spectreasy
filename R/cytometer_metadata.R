#' Supported Cytometers
#'
#' Lists the cytometer IDs accepted by spectreasy metadata helpers.
#'
#' @param include_auto Logical; include `"auto"` in the returned values.
#' @return Character vector of supported cytometer IDs.
#' @examples
#' supported_cytometers()
#' @export
supported_cytometers <- function(include_auto = FALSE) {
    ids <- c(
        "aurora",
        "northern_lights",
        "id7000",
        "discover_s8",
        "discover_a8",
        "a5se",
        "opteon",
        "mosaic",
        "xenith"
    )
    if (isTRUE(include_auto)) c("auto", ids) else ids
}

.spectreasy_extdata_file <- function(filename) {
    p <- system.file("extdata", filename, package = "spectreasy")
    if (nzchar(p) && file.exists(p)) return(p)
    local_p <- file.path("inst", "extdata", filename)
    if (file.exists(local_p)) return(local_p)
    ""
}

.normalize_cytometer_token <- function(x) {
    out <- gsub("[^a-z0-9]+", "", tolower(trimws(as.character(x))))
    out[is.na(out)] <- ""
    out
}

.cytometer_alias_map <- function() {
    aliases <- list(
        auto = c("auto", "automatic", ""),
        aurora = c("aurora", "cytek aurora", "cytekaurora"),
        northern_lights = c(
            "northern_lights", "northern lights", "northernlights",
            "auroranl", "aurora nl", "nl", "cytek northern lights"
        ),
        id7000 = c("id7000", "sony id7000", "sonyid7000"),
        discover_s8 = c(
            "discover_s8", "discover s8", "discover-s8", "s8",
            "facsdiscover s8", "facsdiscovers8", "bd facsdiscover s8"
        ),
        discover_a8 = c(
            "discover_a8", "discover a8", "discover-a8", "a8",
            "facsdiscover a8", "facsdiscovera8", "bd facsdiscover a8"
        ),
        a5se = c(
            "a5se", "a5 se", "symphony", "facssymphony",
            "facssymphony a5 se", "bd facssymphony a5 se"
        ),
        opteon = c("opteon", "novocyte opteon", "agilent opteon"),
        mosaic = c("mosaic", "cytoflex mosaic", "beckman coulter mosaic"),
        xenith = c("xenith", "attune xenith", "thermofisher xenith")
    )
    keys <- unlist(lapply(aliases, .normalize_cytometer_token), use.names = FALSE)
    vals <- rep(names(aliases), lengths(aliases))
    stats::setNames(vals, keys)
}

.resolve_cytometer_id <- function(cytometer, allow_auto = TRUE, unknown_as_auto = TRUE) {
    if (missing(cytometer) || is.null(cytometer) || length(cytometer) == 0) {
        return(if (isTRUE(allow_auto)) "auto" else "")
    }
    key <- .normalize_cytometer_token(cytometer[1])
    aliases <- .cytometer_alias_map()
    if (key %in% names(aliases)) {
        id <- unname(aliases[[key]])
        if (identical(id, "auto") && !isTRUE(allow_auto)) return("")
        return(id)
    }
    if (key %in% supported_cytometers()) return(key)
    if (isTRUE(unknown_as_auto) && isTRUE(allow_auto)) "auto" else key
}

.read_cytometer_dictionary <- function() {
    if (exists("cytometer_dictionary", envir = .spectreasy_cache)) {
        return(get("cytometer_dictionary", envir = .spectreasy_cache))
    }
    path <- .spectreasy_extdata_file("cytometer_dictionary.csv")
    if (!nzchar(path)) return(data.frame())
    out <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    assign("cytometer_dictionary", out, envir = .spectreasy_cache)
    out
}

.read_fluorophore_channel_dictionary <- function() {
    if (exists("fluorophore_channel_dictionary", envir = .spectreasy_cache)) {
        return(get("fluorophore_channel_dictionary", envir = .spectreasy_cache))
    }
    path <- .spectreasy_extdata_file("fluorophore_channel_dictionary.csv")
    if (!nzchar(path)) return(data.frame())
    out <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    assign("fluorophore_channel_dictionary", out, envir = .spectreasy_cache)
    out
}

.normalize_detector_token <- function(x) {
    out <- toupper(gsub("\\s+", "", trimws(as.character(x))))
    out <- gsub("([A-Z]+)-([0-9])", "\\1\\2", out, perl = TRUE)
    out[is.na(out)] <- ""
    out
}

.infer_cytometer_from_pd <- function(pd) {
    if (is.null(pd) || !is.data.frame(pd) || !("name" %in% colnames(pd))) {
        return("auto")
    }
    dict <- .read_cytometer_dictionary()
    if (!is.data.frame(dict) || nrow(dict) == 0 || !all(c("cytometer", "detector") %in% colnames(dict))) {
        return("auto")
    }
    det_info <- tryCatch(get_sorted_detectors(pd), error = function(e) NULL)
    detector_names <- if (!is.null(det_info) && length(det_info$names) > 0) det_info$names else pd$name
    actual <- unique(.normalize_detector_token(detector_names))
    actual <- actual[nzchar(actual)]
    if (length(actual) == 0) return("auto")

    scores <- vapply(split(dict$detector, dict$cytometer), function(detectors) {
        expected <- unique(.normalize_detector_token(detectors))
        expected <- expected[nzchar(expected)]
        common <- length(intersect(actual, expected))
        common / max(1, length(actual))
    }, numeric(1))
    if (length(scores) == 0 || all(!is.finite(scores))) return("auto")
    best <- names(which.max(scores))
    best_score <- unname(max(scores, na.rm = TRUE))
    min_score <- if (length(actual) <= 4) 0.75 else 0.35
    if (!is.finite(best_score) || best_score < min_score) "auto" else best
}

.resolve_cytometer_from_pd <- function(cytometer, pd) {
    id <- .resolve_cytometer_id(cytometer, allow_auto = TRUE, unknown_as_auto = TRUE)
    if (!identical(id, "auto")) return(id)
    .infer_cytometer_from_pd(pd)
}

.resolve_cytometer_from_files <- function(cytometer, files) {
    id <- .resolve_cytometer_id(cytometer, allow_auto = TRUE, unknown_as_auto = TRUE)
    if (!identical(id, "auto")) return(id)
    files <- files[file.exists(files)]
    if (length(files) == 0) return(id)
    for (path in files) {
        inferred <- tryCatch({
            ff <- suppressWarnings(flowCore::read.FCS(path, transformation = FALSE, truncate_max_range = FALSE))
            pd <- flowCore::pData(flowCore::parameters(ff))
            .infer_cytometer_from_pd(pd)
        }, error = function(e) "auto")
        if (!identical(inferred, "auto")) return(inferred)
    }
    id
}
