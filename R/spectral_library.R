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
        attunexenith = "xenith",
        thermofisherxenith = "xenith",
        thermofisherattunexenith = "xenith",
        thermoscientificxenith = "xenith",
        thermoscientificattunexenith = "xenith"
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
