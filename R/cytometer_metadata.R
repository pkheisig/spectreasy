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
        xenith = c(
            "xenith", "attune xenith", "thermofisher xenith",
            "thermo fisher xenith", "thermo fisher attune xenith",
            "thermofisher attune xenith", "thermofisherattunexenith",
            "thermo scientific xenith", "thermo scientific attune xenith"
        )
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
    path <- .spectreasy_extdata_file("cytometer_dictionary.csv")
    if (!nzchar(path)) return(data.frame())
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

.read_fluorophore_channel_dictionary <- function() {
    path <- .spectreasy_extdata_file("fluorophore_channel_dictionary.csv")
    if (!nzchar(path)) return(data.frame())
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

.normalize_detector_token <- function(x) {
    out <- toupper(gsub("\\s+", "", trimws(as.character(x))))
    out <- gsub("([A-Z]+)-([0-9])", "\\1\\2", out, perl = TRUE)
    out[is.na(out)] <- ""
    out
}

.infer_cytometer_from_detector_names <- function(detector_names) {
    dict <- .read_cytometer_dictionary()
    if (!is.data.frame(dict) || nrow(dict) == 0 || !("detector" %in% colnames(dict))) {
        return("auto")
    }
    actual <- unique(.normalize_detector_token(detector_names))
    actual <- actual[nzchar(actual)]
    if (length(actual) == 0) return("auto")
    dict_cytometer <- if ("cytometer" %in% colnames(dict)) {
        vapply(
            dict$cytometer,
            .resolve_cytometer_id,
            character(1),
            allow_auto = FALSE,
            unknown_as_auto = FALSE
        )
    } else {
        rep("", nrow(dict))
    }
    scores <- vapply(split(dict$detector, dict_cytometer), function(detectors) {
        expected <- unique(.normalize_detector_token(detectors))
        expected <- expected[nzchar(expected)]
        length(intersect(actual, expected)) / max(1, length(actual))
    }, numeric(1))
    if (length(scores) == 0 || all(!is.finite(scores))) return("auto")
    best <- names(which.max(scores))
    best_score <- unname(max(scores, na.rm = TRUE))
    min_score <- if (length(actual) <= 4) 0.75 else 0.35
    if (!is.finite(best_score) || best_score < min_score || !nzchar(best)) "auto" else best
}

.detector_laser_excitation_label <- function(cytometer, laser, detector) {
    cytometer <- .resolve_cytometer_id(cytometer, allow_auto = TRUE, unknown_as_auto = TRUE)
    laser <- trimws(as.character(laser)[1])
    detector <- trimws(as.character(detector)[1])
    if (is.na(laser)) laser <- ""
    if (is.na(detector)) detector <- ""
    prefix_nm <- regmatches(detector, regexpr("^[0-9]{3}", detector, perl = TRUE))
    if (length(prefix_nm) > 0 && nzchar(prefix_nm)) return(paste0(prefix_nm, "nm"))

    key <- gsub("\\s+", "", tolower(laser))
    if (!nzchar(key)) {
        det_key <- toupper(gsub("-A$", "", detector, ignore.case = TRUE))
        if (grepl("^UV|^U", det_key)) key <- "uv"
        else if (grepl("^V", det_key)) key <- "violet"
        else if (grepl("^B", det_key)) key <- "blue"
        else if (grepl("^YG|^Y|^G", det_key)) key <- "yellowgreen"
        else if (grepl("^R", det_key)) key <- "red"
        else if (grepl("^IR", det_key)) key <- "ir"
    }

    nm_map_default <- c(deepuv = "320nm", uv = "355nm", violet = "405nm", blue = "488nm", yellowgreen = "561nm", red = "640nm", ir = "808nm")
    nm_map_xenith <- c(deepuv = "349nm", uv = "349nm", violet = "405nm", blue = "488nm", yellowgreen = "561nm", red = "637nm", ir = "781nm")
    nm_map <- if (identical(cytometer, "xenith")) nm_map_xenith else nm_map_default
    if (key %in% names(nm_map)) return(unname(nm_map[[key]]))
    if (nzchar(laser)) laser else "Detector"
}

.aurora_like_detector_emission <- function(detector) {
    clean <- toupper(gsub("-A$", "", trimws(as.character(detector))))
    prefix <- sub("[0-9]+$", "", clean)
    channel <- suppressWarnings(as.integer(regmatches(clean, regexpr("[0-9]+$", clean))))
    if (!is.finite(channel)) return(NA_integer_)
    maps <- list(
        UV = c(370, 395, 420, 440, 450, 480, 480, 500, 520, 550, 570, 580, 600, 660, 750, 800),
        U = c(370, 395, 420, 440, 450, 480, 480, 500, 520, 550, 570, 580, 600, 660, 750, 800),
        V = c(420, 440, 450, 480, 480, 500, 550, 570, 580, 600, 660, 680, 690, 700, 730, 780),
        B = c(500, 520, 550, 550, 570, 580, 600, 600, 660, 680, 690, 700, 750, 780),
        YG = c(570, 580, 600, 600, 660, 680, 700, 730, 750, 780),
        Y = c(570, 580, 600, 600, 660, 680, 700, 730, 750, 780),
        R = c(660, 680, 700, 730, 730, 750, 780, 800),
        IR = c(810, 825, 840, 855, 870, 885)
    )
    values <- maps[[prefix]]
    if (is.null(values) || channel < 1L || channel > length(values)) return(NA_integer_)
    as.integer(values[[channel]])
}

.id7000_like_detector_emission <- function(detector) {
    clean <- toupper(gsub("-A$", "", trimws(as.character(detector))))
    prefix <- regmatches(clean, regexpr("^[0-9]{3}", clean, perl = TRUE))
    channel <- suppressWarnings(as.integer(regmatches(clean, regexpr("(?<=CH)[0-9]+$", clean, perl = TRUE))))
    if (length(prefix) == 0 || !nzchar(prefix) || !is.finite(channel)) return(NA_integer_)
    start_channel <- c(`320` = 1L, `355` = 1L, `405` = 1L, `488` = 4L, `561` = 10L, `637` = 17L, `808` = 36L)
    start_emission <- c(`320` = 350L, `355` = 370L, `405` = 420L, `488` = 500L, `561` = 570L, `637` = 660L, `808` = 810L)
    if (!(prefix %in% names(start_channel))) return(NA_integer_)
    as.integer(start_emission[[prefix]] + (channel - start_channel[[prefix]]) * 15L)
}

.detector_emission_wavelength <- function(cytometer, detector, description = "") {
    cytometer <- .resolve_cytometer_id(cytometer, allow_auto = TRUE, unknown_as_auto = TRUE)
    detector <- trimws(as.character(detector)[1])
    description <- trimws(as.character(description)[1])
    if (is.na(detector)) detector <- ""
    if (is.na(description)) description <- ""

    if (nzchar(description)) {
        hit <- regmatches(description, regexpr("[0-9]{3}(?=/|/LP|-A)", description, perl = TRUE))
        if (length(hit) > 0 && nzchar(hit)) return(as.integer(hit))
    }

    paren <- regmatches(detector, regexpr("(?<=\\()[0-9]{3}(?=\\))", detector, perl = TRUE))
    if (length(paren) > 0 && nzchar(paren)) return(as.integer(paren))

    if (cytometer %in% c("aurora", "northern_lights")) {
        emission <- .aurora_like_detector_emission(detector)
        if (is.finite(emission) && emission > 0) return(as.integer(emission))
    }
    if (identical(cytometer, "id7000")) {
        emission <- .id7000_like_detector_emission(detector)
        if (is.finite(emission) && emission > 0) return(as.integer(emission))
    }

    embedded <- regmatches(detector, regexpr("(?<!^)[0-9]{3}(?=-A$)", detector, perl = TRUE))
    if (length(embedded) > 0 && nzchar(embedded)) return(as.integer(embedded))

    det_key <- toupper(gsub("-A$", "", detector, ignore.case = TRUE))
    channel <- suppressWarnings(as.integer(regmatches(det_key, regexpr("[0-9]+$", det_key, perl = TRUE))))
    if (!is.finite(channel)) return(NA_integer_)
    laser_label <- .detector_laser_excitation_label(cytometer, "", detector)
    start <- if (grepl("^320", laser_label)) 350L else if (grepl("^355|^UV", laser_label, ignore.case = TRUE)) 370L else if (grepl("^405|Violet", laser_label, ignore.case = TRUE)) 420L else if (grepl("^488|Blue", laser_label, ignore.case = TRUE)) 500L else if (grepl("^561|Yellow", laser_label, ignore.case = TRUE)) 570L else if (grepl("^637|640|Red", laser_label, ignore.case = TRUE)) 660L else if (grepl("^781|808|IR", laser_label, ignore.case = TRUE)) 810L else 400L
    as.integer(start + (channel - 1L) * 15L)
}

.format_detector_display_label <- function(cytometer, detector, laser = "", description = "") {
    cytometer <- .resolve_cytometer_id(cytometer, allow_auto = TRUE, unknown_as_auto = TRUE)
    detector <- trimws(as.character(detector)[1])
    description <- trimws(as.character(description)[1])
    if (is.na(detector)) detector <- ""
    if (is.na(description)) description <- ""
    if (cytometer %in% c("aurora", "northern_lights")) return(detector)
    if (nzchar(description)) return(description)
    emission <- .detector_emission_wavelength(cytometer, detector, description = description)
    laser_label <- .detector_laser_excitation_label(cytometer, laser, detector)
    suffix <- if (grepl("-A$", detector, ignore.case = TRUE)) "-A" else ""
    if (is.finite(emission) && emission > 0) {
        paste0(laser_label, " - ", emission, suffix)
    } else {
        gsub("-A$", "", detector, ignore.case = TRUE)
    }
}

.detector_metadata_from_dictionary <- function(detectors, cytometer = "auto") {
    detectors <- trimws(as.character(detectors))
    detectors <- detectors[nzchar(detectors)]
    if (length(detectors) == 0) return(NULL)
    cytometer <- .resolve_cytometer_id(cytometer, allow_auto = TRUE, unknown_as_auto = TRUE)
    if (identical(cytometer, "auto")) {
        cytometer <- .infer_cytometer_from_detector_names(detectors)
    }
    if (identical(cytometer, "auto")) return(NULL)

    dict <- .read_cytometer_dictionary()
    if (!is.data.frame(dict) || nrow(dict) == 0 || !all(c("cytometer", "detector") %in% colnames(dict))) {
        return(NULL)
    }
    dict_cytometer <- vapply(
        dict$cytometer,
        .resolve_cytometer_id,
        character(1),
        allow_auto = FALSE,
        unknown_as_auto = FALSE
    )
    dict <- dict[dict_cytometer == cytometer, , drop = FALSE]
    if (nrow(dict) == 0) return(NULL)

    dict_key <- .normalize_detector_token(dict$detector)
    key <- .normalize_detector_token(detectors)
    hit <- match(key, dict_key)
    ok <- !is.na(hit)
    if (!any(ok)) return(NULL)

    laser <- rep("", length(detectors))
    description <- rep("", length(detectors))
    if ("laser" %in% colnames(dict)) laser[ok] <- as.character(dict$laser[hit[ok]])
    if ("description" %in% colnames(dict)) description[ok] <- as.character(dict$description[hit[ok]])
    label <- mapply(
        .format_detector_display_label,
        cytometer = cytometer,
        detector = detectors,
        laser = laser,
        description = description,
        USE.NAMES = FALSE
    )
    emission <- mapply(
        .detector_emission_wavelength,
        cytometer = cytometer,
        detector = detectors,
        description = description,
        USE.NAMES = FALSE
    )
    data.frame(
        detector = detectors,
        label = label,
        laser = ifelse(nzchar(laser), laser, "Other"),
        emission = as.integer(emission),
        stringsAsFactors = FALSE
    )
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
