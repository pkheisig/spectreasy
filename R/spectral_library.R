.spectral_library_file_map <- function() {
    c(
        aurora = "aurora_spectra.csv",
        discover = "discover_spectra.csv",
        id7000 = "id7000_spectra.csv",
        xenith = "xenith_spectra.csv"
    )
}

.spectral_panel_library_aliases <- function() {
    c(
        aurora = "aurora",
        cytekaurora = "aurora",
        discover = "discover",
        facsdiscover = "discover",
        bdfacsdiscover = "discover",
        discover_s8 = "discover",
        discovers8 = "discover",
        discover_a8 = "discover",
        discovera8 = "discover",
        id7000 = "id7000",
        sonyid7000 = "id7000",
        xenith = "xenith",
        attunexenith = "xenith"
    )
}

.resolve_spectral_panel_cytometer <- function(cytometer) {
    if (missing(cytometer) || is.null(cytometer) || length(cytometer) == 0) {
        return("aurora")
    }
    key <- .normalize_cytometer_token(cytometer[1])
    aliases <- .spectral_panel_library_aliases()
    if (key %in% names(aliases)) return(unname(aliases[[key]]))
    if (key %in% names(.spectral_library_file_map())) return(key)
    stop(
        "Spectral panel builder supports: aurora, discover, id7000, xenith.",
        call. = FALSE
    )
}

.spectral_panel_libraries <- function() {
    data.frame(
        id = c("aurora", "discover", "id7000", "xenith"),
        label = c("Cytek Aurora", "BD FACSDiscover", "Sony ID7000", "Thermo Fisher Attune Xenith"),
        stringsAsFactors = FALSE
    )
}

.spectral_panel_configurations <- function(cytometer = "aurora") {
    id <- .resolve_spectral_panel_cytometer(cytometer)
    if (identical(id, "aurora")) {
        return(data.frame(
            id = c("5l_uv_v_b_yg_r", "4l_uv_v_b_r", "4l_v_b_yg_r", "3l_v_b_r"),
            label = c(
                "Aurora 5L: UV/V/B/YG/R",
                "Aurora 4L: UV/V/B/R",
                "Aurora 4L: V/B/YG/R",
                "Aurora 3L: V/B/R"
            ),
            description = c(
                "16UV-16V-14B-10YG-8R",
                "16UV-16V-14B-8R",
                "16V-14B-10YG-8R",
                "16V-14B-8R"
            ),
            stringsAsFactors = FALSE
        ))
    }
    if (identical(id, "discover")) {
        return(data.frame(
            id = c("discover_s8", "discover_a8"),
            label = c(
                "FACSDiscover S8: UV/V/B/YG/R",
                "FACSDiscover A8: UV/V/B/YG/R"
            ),
            description = c(
                "22UV-20V-16B-12YG-8R",
                "22UV-20V-16B-12YG-8R"
            ),
            stringsAsFactors = FALSE
        ))
    }
    if (identical(id, "id7000")) {
        return(data.frame(
            id = c("id7000_5l", "id7000_4l", "id7000_3l"),
            label = c(
                "ID7000 5L: UV/V/B/YG/R",
                "ID7000 4L: V/B/YG/R",
                "ID7000 3L: V/B/R"
            ),
            description = c(
                "147 fluorescence detectors",
                "112 fluorescence detectors",
                "86 fluorescence detectors"
            ),
            stringsAsFactors = FALSE
        ))
    }

    libs <- .spectral_panel_libraries()
    label <- libs$label[match(id, libs$id)]
    if (is.na(label) || !nzchar(label)) label <- id
    data.frame(
        id = "full",
        label = paste0(label, " full detector set"),
        description = "All packaged detectors",
        stringsAsFactors = FALSE
    )
}

.resolve_spectral_panel_configuration <- function(cytometer = "aurora", configuration = NULL) {
    id <- .resolve_spectral_panel_cytometer(cytometer)
    configs <- .spectral_panel_configurations(id)
    default <- configs$id[1]
    if (missing(configuration) || is.null(configuration) || length(configuration) == 0) {
        return(default)
    }

    key <- .normalize_cytometer_token(configuration[1])
    aliases <- c(
        full = "full",
        discovers8 = "discover_s8",
        discover_s8 = "discover_s8",
        facsdiscovers8 = "discover_s8",
        bdfacsdiscovers8 = "discover_s8",
        discovera8 = "discover_a8",
        discover_a8 = "discover_a8",
        facsdiscovera8 = "discover_a8",
        bdfacsdiscovera8 = "discover_a8",
        id7000 = "id7000_5l",
        id70005l = "id7000_5l",
        id70005laser = "id7000_5l",
        id70005lcompact = "id7000_5l",
        id70004l = "id7000_4l",
        id70004laser = "id7000_4l",
        id70003l = "id7000_3l",
        id70003laser = "id7000_3l",
        aurora5l = "5l_uv_v_b_yg_r",
        aurora5laser = "5l_uv_v_b_yg_r",
        `5l` = "5l_uv_v_b_yg_r",
        `5laser` = "5l_uv_v_b_yg_r",
        `5luv_v_b_yg_r` = "5l_uv_v_b_yg_r",
        `5luv_v_b_ygr` = "5l_uv_v_b_yg_r",
        aurora4luv = "4l_uv_v_b_r",
        aurora4luvvbr = "4l_uv_v_b_r",
        `4luv` = "4l_uv_v_b_r",
        `4luv_v_b_r` = "4l_uv_v_b_r",
        aurora4lyg = "4l_v_b_yg_r",
        aurora4lvbygr = "4l_v_b_yg_r",
        `4lyg` = "4l_v_b_yg_r",
        `4lv_b_yg_r` = "4l_v_b_yg_r",
        aurora3l = "3l_v_b_r",
        aurora3laser = "3l_v_b_r",
        `3l` = "3l_v_b_r",
        `3laser` = "3l_v_b_r",
        `3lv_b_r` = "3l_v_b_r"
    )
    if (key %in% names(aliases)) {
        hit <- unname(aliases[[key]])
        if (hit %in% configs$id) return(hit)
    }
    if (key %in% .normalize_cytometer_token(configs$id)) {
        return(configs$id[match(key, .normalize_cytometer_token(configs$id))])
    }

    default
}

.spectral_panel_configuration_lasers <- function(cytometer = "aurora", configuration = NULL) {
    id <- .resolve_spectral_panel_cytometer(cytometer)
    config <- .resolve_spectral_panel_configuration(id, configuration)
    switch(config,
        `5l_uv_v_b_yg_r` = c("UV", "Violet", "Blue", "YellowGreen", "Red"),
        `4l_uv_v_b_r` = c("UV", "Violet", "Blue", "Red"),
        `4l_v_b_yg_r` = c("Violet", "Blue", "YellowGreen", "Red"),
        `discover_s8` = c("UV", "Violet", "Blue", "YellowGreen", "Red"),
        `discover_a8` = c("UV", "Violet", "Blue", "YellowGreen", "Red"),
        `id7000_5l` = c("UV", "Violet", "Blue", "YellowGreen", "Red"),
        `id7000_4l` = c("Violet", "Blue", "YellowGreen", "Red"),
        `id7000_3l` = c("Violet", "Blue", "Red"),
        `3l_v_b_r` = c("Violet", "Blue", "Red"),
        NULL
    )
}

.spectral_library_path <- function(cytometer) {
    id <- .resolve_spectral_panel_cytometer(cytometer)
    file_map <- .spectral_library_file_map()
    path <- .spectreasy_extdata_file(file_map[[id]])
    if (!nzchar(path) || !file.exists(path)) {
        stop("Spectral library file is missing for cytometer '", id, "'.", call. = FALSE)
    }
    path
}

.read_spectral_library_matrix <- function(cytometer, normalize = TRUE) {
    id <- .resolve_spectral_panel_cytometer(cytometer)
    cache_key <- paste0("lib_", id, "_", normalize)
    if (exists(cache_key, envir = .spectreasy_cache)) {
        return(get(cache_key, envir = .spectreasy_cache))
    }

    path <- .spectral_library_path(id)
    df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    if (ncol(df) < 2) {
        stop("Spectral library has no detector columns: ", path, call. = FALSE)
    }
    names(df)[1] <- "fluorophore"
    fluor <- trimws(as.character(df$fluorophore))
    keep <- !is.na(fluor) & nzchar(fluor)
    df <- df[keep, , drop = FALSE]
    fluor <- fluor[keep]
    values <- as.matrix(df[, -1, drop = FALSE])
    suppressWarnings(storage.mode(values) <- "numeric")
    values[!is.finite(values)] <- 0
    rownames(values) <- fluor
    values <- values[!duplicated(rownames(values)), , drop = FALSE]
    out <- if (isTRUE(normalize)) .normalize_spectral_rows(values) else values

    assign(cache_key, out, envir = .spectreasy_cache)
    out
}

.normalize_spectral_rows <- function(M) {
    if (nrow(M) == 0 || ncol(M) == 0) return(M)
    denom <- apply(abs(M), 1, max, na.rm = TRUE)
    denom[!is.finite(denom) | denom <= 0] <- 1
    sweep(M, 1, denom, "/")
}

.spectral_detector_keys <- function(x) {
    x <- trimws(as.character(x))
    x[is.na(x)] <- ""
    no_paren <- gsub("\\s*\\([^)]*\\)", "", x)
    base <- gsub("-A$", "", no_paren, ignore.case = TRUE)
    keys <- unique(c(
        .normalize_detector_token(x),
        .normalize_detector_token(no_paren),
        .normalize_detector_token(base),
        .normalize_detector_token(paste0(base, "-A"))
    ))
    keys[nzchar(keys)]
}

.match_spectral_detectors <- function(library_detectors, requested_detectors, strict = FALSE) {
    lib_map <- character()
    for (det in library_detectors) {
        keys <- .spectral_detector_keys(det)
        keys <- keys[!(keys %in% names(lib_map))]
        lib_map[keys] <- det
    }

    matched <- character(length(requested_detectors))
    for (i in seq_along(requested_detectors)) {
        keys <- .spectral_detector_keys(requested_detectors[i])
        hit <- keys[keys %in% names(lib_map)]
        matched[i] <- if (length(hit) > 0) lib_map[[hit[1]]] else ""
    }
    names(matched) <- requested_detectors

    missing <- names(matched)[!nzchar(matched)]
    if (length(missing) > 0 && isTRUE(strict)) {
        stop("Spectral library is missing detector(s): ", paste(missing, collapse = ", "), call. = FALSE)
    }
    matched[nzchar(matched)]
}

.spectral_fluor_key <- function(x) {
    gsub("[^a-z0-9]+", "", tolower(trimws(as.character(x))))
}

.match_spectral_fluorophores <- function(library_fluors, requested_fluors, strict = FALSE) {
    lib_map <- stats::setNames(library_fluors, .spectral_fluor_key(library_fluors))
    dict_file <- .spectreasy_extdata_file("fluorophore_dictionary.csv")
    if (nzchar(dict_file) && file.exists(dict_file)) {
        dict <- tryCatch(utils::read.csv(dict_file, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
        if (!is.null(dict) && all(c("fluorophore", "aliases") %in% colnames(dict))) {
            for (i in seq_len(nrow(dict))) {
                canonical <- trimws(as.character(dict$fluorophore[i]))
                if (!nzchar(canonical)) next
                lib_hit <- unname(lib_map[.spectral_fluor_key(canonical)])
                if (length(lib_hit) == 0 || is.na(lib_hit) || !nzchar(lib_hit)) next
                aliases <- unique(c(canonical, .control_file_split_semicolon(dict$aliases[i])))
                alias_keys <- .spectral_fluor_key(aliases)
                alias_keys <- alias_keys[nzchar(alias_keys) & !(alias_keys %in% names(lib_map))]
                lib_map[alias_keys] <- lib_hit
            }
        }
    }

    requested <- trimws(as.character(requested_fluors))
    matched <- unname(lib_map[.spectral_fluor_key(requested)])
    matched[is.na(matched)] <- ""
    names(matched) <- requested_fluors

    missing <- names(matched)[!nzchar(matched)]
    if (length(missing) > 0 && isTRUE(strict)) {
        stop("Spectral library is missing fluorophore(s): ", paste(missing, collapse = ", "), call. = FALSE)
    }
    matched
}

.load_spectral_library <- function(cytometer = "aurora",
                                   fluorophores = NULL,
                                   detectors = NULL,
                                   strict = FALSE,
                                   renormalize = TRUE) {
    id <- .resolve_spectral_panel_cytometer(cytometer)
    M <- .read_spectral_library_matrix(id)

    if (!is.null(fluorophores)) {
        fluor_hits <- .match_spectral_fluorophores(rownames(M), fluorophores, strict = strict)
        keep <- fluor_hits[nzchar(fluor_hits)]
        M <- M[unname(keep), , drop = FALSE]
        rownames(M) <- if (length(keep) > 0) names(keep) else character()
    }

    if (!is.null(detectors)) {
        detector_hits <- .match_spectral_detectors(colnames(M), detectors, strict = strict)
        M <- M[, unname(detector_hits), drop = FALSE]
        colnames(M) <- names(detector_hits)
    }

    attr(M, "cytometer") <- id
    if (isTRUE(renormalize)) .normalize_spectral_rows(M) else M
}

.spectral_cosine <- function(a, b) {
    denom <- sqrt(sum(a^2, na.rm = TRUE) + 1e-9) * sqrt(sum(b^2, na.rm = TRUE) + 1e-9)
    if (!is.finite(denom) || denom <= 0) return(NA_real_)
    as.numeric(sum(a * b, na.rm = TRUE) / denom)
}

.compare_control_spectra <- function(M,
                                     cytometer = "aurora",
                                     similarity_warn = 0.95) {
    M <- .as_reference_matrix(M, "M")
    measured <- .normalize_spectral_rows(M)
    library <- .load_spectral_library(
        cytometer = cytometer,
        detectors = colnames(measured),
        strict = FALSE
    )
    common_detectors <- intersect(colnames(measured), colnames(library))
    if (length(common_detectors) == 0) {
        stop("No shared detectors between measured spectra and spectral library.", call. = FALSE)
    }
    measured <- measured[, common_detectors, drop = FALSE]
    library <- library[, common_detectors, drop = FALSE]

    expected_hits <- .match_spectral_fluorophores(rownames(library), rownames(measured), strict = FALSE)
    rows <- lapply(seq_len(nrow(measured)), function(i) {
        fluor <- rownames(measured)[i]
        y <- as.numeric(measured[i, ])
        expected <- expected_hits[[fluor]]
        best_scores <- vapply(seq_len(nrow(library)), function(j) .spectral_cosine(y, as.numeric(library[j, ])), numeric(1))
        best_idx <- if (all(is.na(best_scores))) NA_integer_ else which.max(best_scores)
        best_name <- if (is.na(best_idx)) "" else rownames(library)[best_idx]
        best_score <- if (is.na(best_idx)) NA_real_ else best_scores[best_idx]

        if (!nzchar(expected)) {
            return(data.frame(
                fluorophore = fluor,
                library_match = "",
                cosine_similarity = NA_real_,
                correlation = NA_real_,
                max_abs_delta = NA_real_,
                measured_peak = colnames(measured)[which.max(y)],
                library_peak = "",
                peak_match = NA,
                best_library_match = best_name,
                best_library_similarity = best_score,
                flag = "missing_library",
                stringsAsFactors = FALSE
            ))
        }

        x <- as.numeric(library[expected, ])
        measured_peak <- colnames(measured)[which.max(y)]
        library_peak <- colnames(library)[which.max(x)]
        cos <- .spectral_cosine(y, x)
        corr <- suppressWarnings(stats::cor(y, x, use = "pairwise.complete.obs"))
        if (!is.finite(corr)) corr <- NA_real_
        delta <- max(abs(y - x), na.rm = TRUE)
        flag <- if (!is.finite(cos)) "missing_library" else if (cos < similarity_warn) "warn" else "ok"
        data.frame(
            fluorophore = fluor,
            library_match = expected,
            cosine_similarity = cos,
            correlation = corr,
            max_abs_delta = delta,
            measured_peak = measured_peak,
            library_peak = library_peak,
            peak_match = identical(measured_peak, library_peak),
            best_library_match = best_name,
            best_library_similarity = best_score,
            flag = flag,
            stringsAsFactors = FALSE
        )
    })

    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    attr(out, "cytometer") <- attr(library, "cytometer")
    attr(out, "detectors") <- common_detectors
    out
}

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
    out <- data.frame(
        detector = detectors,
        label = gsub("-A$", "", gsub("\\s*\\([^)]*\\)", "", detectors)),
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

.calculate_panel_complexity <- function(spectra) {
    spectra <- .as_reference_matrix(spectra, "spectra")
    if (nrow(spectra) < 2) return(1)
    sv <- svd(spectra)$d
    sv <- sv[is.finite(sv) & sv > 0]
    if (length(sv) == 0) return(NA_real_)
    round(max(sv) / min(sv), 2)
}

.matrix_to_named_rows <- function(M, row_name = "fluorophore") {
    if (length(M) == 0 || nrow(M) == 0) return(data.frame())
    out <- as.data.frame(M, check.names = FALSE)
    out[[row_name]] <- rownames(M)
    out[, c(row_name, setdiff(colnames(out), row_name)), drop = FALSE]
}

.build_spectral_panel_data <- function(fluorophores,
                                       cytometer = "aurora",
                                       configuration = NULL,
                                       strict = TRUE) {
    spectra <- .spectral_panel_configuration_spectra(
        cytometer = cytometer,
        configuration = configuration,
        fluorophores = fluorophores,
        strict = strict
    )
    peak_detectors <- apply(spectra, 1, function(x) colnames(spectra)[which.max(x)])
    list(
        spectra = spectra,
        similarity_matrix = calculate_similarity_matrix(spectra),
        complexity_index = .calculate_panel_complexity(spectra),
        peak_detectors = peak_detectors
    )
}

.spectral_panel_export_table <- function(fluorophores, markers = NULL) {
    fluorophores <- trimws(as.character(fluorophores))
    fluorophores <- fluorophores[nzchar(fluorophores)]
    markers <- if (is.null(markers)) rep("", length(fluorophores)) else trimws(as.character(markers))
    length(markers) <- length(fluorophores)
    markers[is.na(markers)] <- ""
    data.frame(
        Marker = markers,
        Fluorophore = fluorophores,
        stringsAsFactors = FALSE
    )
}

.spectral_panel_label_rows <- function(fluorophores, markers = NULL) {
    export_df <- .spectral_panel_export_table(fluorophores, markers)
    labels <- ifelse(
        nzchar(export_df$Marker),
        paste0(export_df$Marker, " / ", export_df$Fluorophore),
        export_df$Fluorophore
    )
    make.unique(labels, sep = " ")
}

.spectral_panel_signature_bins <- function(signature,
                                           detector_labels,
                                           y_power = 1.5,
                                           max_log_intensity = 6,
                                           n_bins = 37L) {
    signature <- as.numeric(signature)
    signature[!is.finite(signature)] <- 0
    signature <- pmax(0, pmin(1, signature))

    offsets <- seq(-0.42, 0.42, length.out = n_bins)
    center_weight <- 1 - abs(seq(-1, 1, length.out = n_bins))^1.6
    center_weight <- pmax(0.08, center_weight)

    rows <- lapply(seq_along(signature), function(i) {
        center_log <- 0.35 + (signature[i]^0.72) * (max_log_intensity - 0.35)
        y_orig <- pmax(0.05, pmin(max_log_intensity, center_log + offsets))
        data.frame(
            ch_idx = i,
            Detector = detector_labels[i],
            y_orig = y_orig,
            y = y_orig^y_power,
            fill = pmin(1, center_weight * pmax(0.12, sqrt(signature[i])) * 1.35),
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    out$Detector <- factor(out$Detector, levels = detector_labels)
    out
}

.plot_spectral_panel_signature <- function(signature,
                                           detector_labels,
                                           title,
                                           y_power = 1.5,
                                           max_log_intensity = 6) {
    dt_c <- .spectral_panel_signature_bins(
        signature = signature,
        detector_labels = detector_labels,
        y_power = y_power,
        max_log_intensity = max_log_intensity
    )
    y_breaks_orig <- 0:max_log_intensity
    y_breaks_trans <- y_breaks_orig^y_power
    y_labels <- vapply(y_breaks_orig, function(x) paste0("10^", x), character(1))

    ggplot2::ggplot(dt_c, ggplot2::aes(ch_idx, y, fill = fill)) +
        ggplot2::geom_tile(width = 0.76, height = 0.055) +
        ggplot2::scale_fill_gradientn(
            colors = c("#0000FF", "#00FFFF", "#00FF00", "#FFFF00", "#FF0000"),
            limits = c(0, 1),
            oob = scales::squish,
            guide = "none"
        ) +
        ggplot2::scale_x_continuous(
            breaks = seq_along(detector_labels),
            labels = detector_labels
        ) +
        ggplot2::scale_y_continuous(
            limits = c(-0.1, (max_log_intensity + 0.35)^y_power),
            breaks = y_breaks_trans,
            labels = y_labels
        ) +
        ggplot2::coord_cartesian(expand = FALSE) +
        ggplot2::labs(title = title, x = NULL, y = "Intensity") +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5),
            axis.text.y = ggplot2::element_text(size = 8),
            axis.title.y = ggplot2::element_text(size = 9, face = "bold"),
            panel.grid.major = ggplot2::element_line(color = "#eeeeee", linewidth = 0.25),
            panel.grid.minor = ggplot2::element_line(color = "#f4f4f4", linewidth = 0.18),
            panel.background = ggplot2::element_rect(fill = "white", color = NA),
            plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
            plot.margin = ggplot2::margin(4, 8, 4, 4)
        )
}

.plot_spectral_panel_bands <- function(spectra, cytometer = "aurora", markers = NULL) {
    spectra <- .as_reference_matrix(spectra, "spectra")
    labels <- .spectral_panel_label_rows(rownames(spectra), markers)
    detector_info <- .spectral_detector_metadata(cytometer, colnames(spectra))
    detector_labels <- detector_info$label

    pages <- vector("list", nrow(spectra))
    for (i in seq_len(nrow(spectra))) {
        pages[[i]] <- .plot_spectral_panel_signature(
            signature = spectra[i, ],
            detector_labels = detector_labels,
            title = labels[i]
        )
    }
    pages
}

.draw_spectral_panel_signature_pages <- function(signature_pages, plots_per_page = 2L) {
    if (length(signature_pages) == 0) return(invisible(NULL))
    plots_per_page <- max(1L, as.integer(plots_per_page[1]))
    starts <- seq(1L, length(signature_pages), by = plots_per_page)
    for (start_idx in starts) {
        page_plots <- signature_pages[start_idx:min(start_idx + plots_per_page - 1L, length(signature_pages))]
        grid::grid.newpage()
        if (length(page_plots) == 1L) {
            grid::grid.draw(grid::editGrob(
                ggplot2::ggplotGrob(page_plots[[1]]),
                vp = grid::viewport(x = 0.5, y = 0.50, width = 0.94, height = 0.72)
            ))
        } else {
            grid::grid.draw(grid::editGrob(
                ggplot2::ggplotGrob(page_plots[[1]]),
                vp = grid::viewport(x = 0.5, y = 0.72, width = 0.94, height = 0.43)
            ))
            grid::grid.draw(grid::editGrob(
                ggplot2::ggplotGrob(page_plots[[2]]),
                vp = grid::viewport(x = 0.5, y = 0.27, width = 0.94, height = 0.43)
            ))
        }
    }
    invisible(NULL)
}

.draw_spectral_panel_similarity_pages <- function(similarity_matrix,
                                                  complexity_index,
                                                  cytometer_label = "",
                                                  configuration_label = "",
                                                  fluorophore_count = NULL,
                                                  max_markers_per_page = NULL) {
    if (is.null(similarity_matrix) || nrow(similarity_matrix) < 2L) return(invisible(NULL))
    n_markers <- nrow(similarity_matrix)
    if (is.null(max_markers_per_page)) {
        max_markers_per_page <- if (n_markers <= 15L) 15L else 8L
    }
    pages <- .build_qc_report_matrix_pages(
        similarity_matrix,
        plot_fun = plot_similarity_matrix,
        max_markers_per_page = max_markers_per_page,
        item_label = "Fluorophores"
    )
    if (length(pages) == 0) return(invisible(NULL))

    panel_label_parts <- c(
        cytometer_label,
        configuration_label,
        if (!is.null(fluorophore_count)) paste0(fluorophore_count, " fluorophore(s)") else NULL
    )
    panel_label <- paste(panel_label_parts[nzchar(panel_label_parts)], collapse = " | ")

    for (p in pages) {
        p <- p +
            ggplot2::theme(
                plot.title = ggplot2::element_blank(),
                plot.subtitle = ggplot2::element_blank(),
                axis.text.x = ggplot2::element_text(size = 8, angle = 45, hjust = 1),
                axis.text.y = ggplot2::element_text(size = 8),
                legend.position = "right",
                plot.margin = ggplot2::margin(8, 8, 8, 8)
            )
        grid::grid.newpage()
        grid::grid.text(
            "Fluorophore Spectral Similarity",
            x = 0.04,
            y = 0.965,
            just = c("left", "top"),
            gp = grid::gpar(fontsize = 18, fontface = "bold")
        )
        if (nzchar(panel_label)) {
            grid::grid.text(
                panel_label,
                x = 0.04,
                y = 0.925,
                just = c("left", "top"),
                gp = grid::gpar(fontsize = 10, col = "#475569")
            )
        }
        grid::grid.text(
            paste0("Complexity Index: ", format(round(complexity_index, 2), nsmall = 2)),
            x = 0.965,
            y = 0.955,
            just = c("right", "top"),
            gp = grid::gpar(fontsize = 12, fontface = "bold")
        )
        grid::grid.draw(grid::editGrob(
            ggplot2::ggplotGrob(p),
            vp = grid::viewport(x = 0.50, y = 0.46, width = 0.94, height = 0.78)
        ))
    }
    invisible(NULL)
}

.write_spectral_panel_overview_pdf <- function(cytometer = "aurora",
                                               configuration = NULL,
                                               fluorophores,
                                               markers = NULL,
                                               output_file) {
    if (missing(output_file) || is.null(output_file) || !nzchar(trimws(as.character(output_file)[1]))) {
        stop("output_file is required.", call. = FALSE)
    }

    requested_fluorophores <- unique(trimws(as.character(fluorophores)))
    requested_fluorophores <- requested_fluorophores[nzchar(requested_fluorophores)]
    requested_markers <- if (is.null(markers)) rep("", length(requested_fluorophores)) else trimws(as.character(markers))
    length(requested_markers) <- length(requested_fluorophores)
    requested_markers[is.na(requested_markers)] <- ""
    marker_map <- stats::setNames(requested_markers, requested_fluorophores)

    panel <- .build_spectral_panel_data(
        fluorophores = requested_fluorophores,
        cytometer = cytometer,
        configuration = configuration,
        strict = FALSE
    )
    spectra <- panel$spectra
    if (nrow(spectra) == 0) {
        stop("Select at least one fluorophore before exporting a panel overview.", call. = FALSE)
    }

    markers <- unname(marker_map[rownames(spectra)])
    markers[is.na(markers)] <- ""

    id <- .resolve_spectral_panel_cytometer(cytometer)
    config <- .resolve_spectral_panel_configuration(id, configuration)
    libs <- .spectral_panel_libraries()
    cyt_label <- libs$label[match(id, libs$id)]
    if (is.na(cyt_label) || !nzchar(cyt_label)) cyt_label <- id
    configs <- .spectral_panel_configurations(id)
    config_label <- configs$label[match(config, configs$id)]
    if (is.na(config_label) || !nzchar(config_label)) config_label <- config

    signature_pages <- .plot_spectral_panel_bands(spectra, cytometer = id, markers = markers)
    sim <- panel$similarity_matrix

    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    grDevices::pdf(output_file, width = 11, height = 8.5)
    on.exit(grDevices::dev.off(), add = TRUE)

    .draw_spectral_panel_similarity_pages(
        similarity_matrix = sim,
        complexity_index = panel$complexity_index,
        cytometer_label = cyt_label,
        configuration_label = config_label,
        fluorophore_count = nrow(spectra)
    )
    .draw_spectral_panel_signature_pages(signature_pages, plots_per_page = 2L)

    invisible(output_file)
}

.spectral_panel_payload <- function(cytometer = "aurora",
                                    fluorophores = character(),
                                    configuration = NULL) {
    id <- .resolve_spectral_panel_cytometer(cytometer)
    config <- .resolve_spectral_panel_configuration(id, configuration)
    full_library <- .load_spectral_library(id, renormalize = FALSE)
    config_detectors <- .spectral_panel_configuration_detectors(
        cytometer = id,
        configuration = config,
        detectors = colnames(full_library)
    )
    retained_signal <- apply(abs(full_library[, config_detectors, drop = FALSE]), 1, max, na.rm = TRUE)
    retained_signal[!is.finite(retained_signal)] <- 0
    available <- retained_signal >= 0.02
    library <- .normalize_spectral_rows(full_library[available, config_detectors, drop = FALSE])
    detector_info <- .spectral_detector_metadata(id, colnames(library))
    fluorophores <- unique(trimws(as.character(fluorophores)))
    fluorophores <- fluorophores[nzchar(fluorophores)]

    selected <- .spectral_panel_configuration_spectra(
        cytometer = id,
        configuration = config,
        fluorophores = fluorophores,
        strict = FALSE
    )
    selected_names <- rownames(selected)
    if (is.null(selected_names)) selected_names <- character()
    sim <- if (length(selected_names) > 0) calculate_similarity_matrix(selected) else matrix(numeric(0), nrow = 0, ncol = 0)
    complexity <- if (length(selected_names) > 0) .calculate_panel_complexity(selected) else NA_real_
    peaks <- if (length(selected_names) > 0) {
        apply(selected, 1, function(x) colnames(selected)[which.max(x)])
    } else {
        character()
    }

    fluor_table <- data.frame(
        fluorophore = rownames(library),
        peak_detector = apply(library, 1, function(x) colnames(library)[which.max(x)]),
        retained_signal = as.numeric(retained_signal[rownames(library)]),
        stringsAsFactors = FALSE
    )
    fluor_table$peak_laser <- detector_info$laser[match(fluor_table$peak_detector, detector_info$detector)]
    fluor_table$peak_color <- detector_info$color[match(fluor_table$peak_detector, detector_info$detector)]
    fluor_table <- fluor_table[order(fluor_table$peak_laser, fluor_table$fluorophore), , drop = FALSE]

    list(
        cytometer = id,
        configuration = config,
        libraries = .spectral_panel_libraries(),
        configurations = .spectral_panel_configurations(id),
        detectors = detector_info,
        fluorophores = fluor_table,
        selected = selected_names,
        spectra = .matrix_to_named_rows(selected),
        similarity = .matrix_to_named_rows(sim),
        complexity_index = complexity,
        peak_detectors = unname(peaks)
    )
}

#' Launch Spectral Panel Builder
#'
#' Opens an interactive browser-based spectral panel builder for packaged
#' theoretical spectra. The panel builder currently supports Aurora, Discover,
#' ID7000, and Xenith libraries.
#'
#' @param port API port (default: 8000).
#' @param open_browser Logical. Open browser automatically?
#' @param dev_mode Logical. Use the Vite dev server instead of bundled assets.
#' @return Invisibly returns NULL. This function blocks while the API is running.
#' @export
#' @examples
#' if (interactive()) {
#'   build_spectral_panel(open_browser = FALSE)
#' }
build_spectral_panel <- function(port = 8000,
                                 open_browser = TRUE,
                                 dev_mode = FALSE) {
    .launch_spectreasy_gui(
        matrix_dir = NULL,
        samples_dir = NULL,
        port = port,
        open_browser = open_browser,
        dev_mode = dev_mode,
        mode = "panel-builder",
        panel_cytometer = "aurora"
    )
    invisible(NULL)
}
