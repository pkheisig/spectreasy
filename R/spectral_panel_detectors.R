.spectral_panel_laser_palette <- function() {
    c(
        DeepUV = "#4c1d95",
        UV = "#6f006f",
        Violet = "#9d00d8",
        Blue = "#0757f2",
        YellowGreen = "#9acd2f",
        Red = "#ff140f",
        IR = "#7f1d1d",
        Other = "#64748b"
    )
}

.spectral_detector_laser <- function(cytometer, detector) {
    dict <- .read_cytometer_dictionary()
    id <- .resolve_spectral_panel_cytometer(cytometer)
    if (nrow(dict) > 0 && all(c("cytometer", "detector", "laser") %in% colnames(dict))) {
        candidates <- dict[dict$cytometer %in% c(id, if (identical(id, "discover")) c("discover_s8", "discover_a8") else id), , drop = FALSE]
        detector_keys <- .spectral_detector_keys(detector)
        idx <- which(vapply(candidates$detector, function(x) any(.spectral_detector_keys(x) %in% detector_keys), logical(1)))
        if (length(idx) > 0) return(as.character(candidates$laser[idx[1]]))
    }
    if (grepl("^320", detector)) return("DeepUV")
    if (grepl("^UV|^355", detector, ignore.case = TRUE)) return("UV")
    if (grepl("^V|^405", detector, ignore.case = TRUE)) return("Violet")
    if (grepl("^B|^488", detector, ignore.case = TRUE)) return("Blue")
    if (grepl("^YG|^Y|^561", detector, ignore.case = TRUE)) return("YellowGreen")
    if (grepl("^R|^637|^640", detector, ignore.case = TRUE)) return("Red")
    if (grepl("^IR|^808|^781", detector, ignore.case = TRUE)) return("IR")
    "Other"
}

.spectral_detector_emission <- function(cytometer, detector) {
    id <- .resolve_spectral_panel_cytometer(cytometer)
    dict <- .read_cytometer_dictionary()
    if (nrow(dict) > 0 && all(c("cytometer", "detector", "description") %in% colnames(dict))) {
        candidates <- dict[dict$cytometer %in% c(id, if (identical(id, "discover")) c("discover_s8", "discover_a8") else id), , drop = FALSE]
        detector_keys <- .spectral_detector_keys(detector)
        idx <- which(vapply(candidates$detector, function(x) any(.spectral_detector_keys(x) %in% detector_keys), logical(1)))
        if (length(idx) > 0) {
            desc <- as.character(candidates$description[idx[1]])
            hit <- regmatches(desc, regexpr("[0-9]{3}(?=/|/LP|-A)", desc, perl = TRUE))
            if (length(hit) > 0 && nzchar(hit)) return(as.integer(hit))
        }
    }
    paren <- regmatches(detector, regexpr("(?<=\\()[0-9]{3}(?=\\))", detector, perl = TRUE))
    if (length(paren) > 0 && nzchar(paren)) return(as.integer(paren))
    embedded <- regmatches(detector, regexpr("(?<!^)[0-9]{3}(?=-A$)", detector, perl = TRUE))
    if (length(embedded) > 0 && nzchar(embedded)) return(as.integer(embedded))
    if (identical(id, "aurora")) {
        emission <- .aurora_detector_emission(detector)
        if (is.finite(emission) && emission > 0) return(as.integer(emission))
    }
    if (identical(id, "id7000")) {
        emission <- .id7000_detector_emission(detector)
        if (is.finite(emission) && emission > 0) return(as.integer(emission))
    }
    laser <- .spectral_detector_laser(id, detector)
    offset <- suppressWarnings(as.integer(regmatches(detector, regexpr("[0-9]+(?=(?:-A)?$)", detector, perl = TRUE))))
    if (!is.finite(offset)) offset <- 1L
    starts <- c(DeepUV = 350L, UV = 370L, Violet = 420L, Blue = 500L, YellowGreen = 570L, Red = 660L, IR = 810L, Other = 400L)
    start <- if (laser %in% names(starts)) starts[[laser]] else starts[["Other"]]
    as.integer(start + (offset - 1L) * 15L)
}

.aurora_detector_emission <- function(detector) {
    clean <- toupper(gsub("-A$", "", trimws(as.character(detector))))
    prefix <- sub("[0-9]+$", "", clean)
    channel <- suppressWarnings(as.integer(regmatches(clean, regexpr("[0-9]+$", clean))))
    if (!is.finite(channel)) return(NA_integer_)
    maps <- list(
        UV = c(370, 395, 420, 440, 450, 480, 480, 500, 520, 550, 570, 580, 600, 660, 750, 800),
        V = c(420, 440, 450, 480, 480, 500, 550, 570, 580, 600, 660, 680, 690, 700, 730, 780),
        B = c(500, 520, 550, 550, 570, 580, 600, 600, 660, 680, 690, 700, 750, 780),
        YG = c(570, 580, 600, 600, 660, 680, 700, 730, 750, 780),
        R = c(660, 680, 700, 730, 730, 750, 780, 800)
    )
    values <- maps[[prefix]]
    if (is.null(values) || channel < 1L || channel > length(values)) return(NA_integer_)
    values[[channel]]
}

.id7000_detector_emission <- function(detector) {
    clean <- toupper(gsub("-A$", "", trimws(as.character(detector))))
    prefix <- regmatches(clean, regexpr("^[0-9]{3}", clean, perl = TRUE))
    channel <- suppressWarnings(as.integer(regmatches(clean, regexpr("(?<=CH)[0-9]+$", clean, perl = TRUE))))
    if (length(prefix) == 0 || !nzchar(prefix) || !is.finite(channel)) return(NA_integer_)
    start_channel <- c(`320` = 1L, `355` = 1L, `405` = 1L, `488` = 4L, `561` = 10L, `637` = 17L, `808` = 36L)
    start_emission <- c(`320` = 350L, `355` = 370L, `405` = 420L, `488` = 500L, `561` = 570L, `637` = 660L, `808` = 810L)
    if (!(prefix %in% names(start_channel))) return(NA_integer_)
    as.integer(start_emission[[prefix]] + (channel - start_channel[[prefix]]) * 15L)
}

.spectral_detector_channel_index <- function(detector) {
    clean <- toupper(trimws(gsub("-A$", "", gsub("\\s*\\([^)]*\\)", "", as.character(detector)))))
    hit <- regmatches(clean, regexpr("(?<=CH)[0-9]+$", clean, perl = TRUE))
    if (length(hit) == 0 || !nzchar(hit)) {
        hit <- regmatches(clean, regexpr("[0-9]+$", clean, perl = TRUE))
    }
    if (length(hit) == 0 || !nzchar(hit)) return(0L)
    out <- suppressWarnings(as.integer(hit))
    if (is.finite(out)) out else 0L
}

.spectral_detector_metadata <- function(cytometer, detectors) {
    id <- .resolve_spectral_panel_cytometer(cytometer)
    cache_key <- paste0("det_", id, "_", paste(.normalize_detector_token(detectors), collapse = "_"))
    if (exists(cache_key, envir = .spectreasy_cache)) {
        return(get(cache_key, envir = .spectreasy_cache))
    }

    palette <- .spectral_panel_laser_palette()
    lasers <- vapply(detectors, function(det) .spectral_detector_laser(id, det), character(1))
    emissions <- vapply(detectors, function(det) .spectral_detector_emission(id, det), integer(1))
    dict <- .read_cytometer_dictionary()
    descriptions <- rep("", length(detectors))
    if (is.data.frame(dict) && nrow(dict) > 0 && all(c("cytometer", "detector", "description") %in% colnames(dict))) {
        dict_cytometers <- vapply(
            dict$cytometer,
            .resolve_cytometer_id,
            character(1),
            allow_auto = FALSE,
            unknown_as_auto = FALSE
        )
        candidate_ids <- c(id, if (identical(id, "discover")) c("discover_s8", "discover_a8") else id)
        candidates <- dict[dict_cytometers %in% candidate_ids, , drop = FALSE]
        if (nrow(candidates) > 0) {
            dict_keys <- vapply(candidates$detector, function(x) .spectral_detector_keys(x)[1], character(1))
            for (i in seq_along(detectors)) {
                hit <- match(.spectral_detector_keys(detectors[i])[1], dict_keys)
                if (!is.na(hit)) descriptions[i] <- trimws(as.character(candidates$description[hit]))
            }
        }
    }
    labels <- mapply(
        .format_detector_display_label,
        cytometer = id,
        detector = detectors,
        laser = lasers,
        description = descriptions,
        USE.NAMES = FALSE
    )
    out <- data.frame(
        detector = detectors,
        label = labels,
        laser = lasers,
        emission = emissions,
        color = unname(palette[ifelse(lasers %in% names(palette), lasers, "Other")]),
        stringsAsFactors = FALSE
    )

    laser_order_ref <- c("DeepUV", "UV", "Violet", "Blue", "YellowGreen", "Red", "IR", "Other")
    laser_rank <- match(out$laser, laser_order_ref)
    laser_rank[is.na(laser_rank)] <- length(laser_order_ref)

    channel_num <- vapply(out$detector, .spectral_detector_channel_index, integer(1))

    out <- out[order(laser_rank, out$emission, channel_num, out$detector), , drop = FALSE]
    rownames(out) <- NULL

    assign(cache_key, out, envir = .spectreasy_cache)
    out
}

.spectral_panel_configuration_detectors <- function(cytometer = "aurora",
                                                    configuration = NULL,
                                                    detectors = NULL) {
    id <- .resolve_spectral_panel_cytometer(cytometer)
    if (is.null(detectors)) {
        detectors <- colnames(.load_spectral_library(id, renormalize = FALSE))
    }
    lasers <- .spectral_panel_configuration_lasers(id, configuration)
    detector_info <- .spectral_detector_metadata(id, detectors)
    if (is.null(lasers)) return(detector_info$detector)
    detector_info$detector[detector_info$laser %in% lasers]
}

.spectral_panel_configuration_spectra <- function(cytometer = "aurora",
                                                  configuration = NULL,
                                                  fluorophores = NULL,
                                                  strict = FALSE,
                                                  min_retained_signal = 0.02) {
    id <- .resolve_spectral_panel_cytometer(cytometer)
    full <- .load_spectral_library(
        cytometer = id,
        fluorophores = fluorophores,
        strict = strict,
        renormalize = FALSE
    )
    config_detectors <- .spectral_panel_configuration_detectors(
        cytometer = id,
        configuration = configuration,
        detectors = colnames(full)
    )
    config_detectors <- intersect(config_detectors, colnames(full))
    if (length(config_detectors) == 0) {
        stop("Selected spectral panel configuration has no matching detectors.", call. = FALSE)
    }

    retained_signal <- apply(abs(full[, config_detectors, drop = FALSE]), 1, max, na.rm = TRUE)
    retained_signal[!is.finite(retained_signal)] <- 0
    keep <- retained_signal >= min_retained_signal
    unavailable <- rownames(full)[!keep]
    if (length(unavailable) > 0 && isTRUE(strict) && !is.null(fluorophores)) {
        stop(
            "Fluorophore(s) have too little signal in this detector configuration: ",
            paste(unavailable, collapse = ", "),
            call. = FALSE
        )
    }

    out <- full[keep, config_detectors, drop = FALSE]
    out <- .normalize_spectral_rows(out)
    attr(out, "cytometer") <- id
    attr(out, "configuration") <- .resolve_spectral_panel_configuration(id, configuration)
    attr(out, "retained_signal") <- retained_signal[keep]
    out
}
