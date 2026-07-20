.merge_control_file_alias_map <- function(existing, incoming) {
    if (length(incoming) == 0) return(existing)
    incoming <- incoming[!is.na(names(incoming)) & nzchar(names(incoming))]
    incoming <- incoming[!(names(incoming) %in% names(existing))]
    c(existing, incoming)
}

.control_file_split_filename_tokens <- function(text) {
    parts <- unlist(strsplit(as.character(text), "[^A-Za-z0-9]+"))
    parts <- .control_file_normalize_token(parts)
    parts[nzchar(parts)]
}

.control_file_build_compound_tokens <- function(tokens, max_n = 4) {
    tokens <- tokens[nzchar(tokens)]
    n <- length(tokens)
    if (n == 0) return(character())
    out <- tokens
    max_n <- min(max_n, n)
    if (max_n >= 2) {
        for (w in 2:max_n) {
            for (i in seq_len(n - w + 1)) {
                out <- c(out, paste0(tokens[i:(i + w - 1)], collapse = ""))
            }
        }
    }
    unique(out[nzchar(out)])
}

.control_file_resolve_custom_fluor <- function(fn, custom_fluorophores = NULL) {
    if (is.null(custom_fluorophores) || length(custom_fluorophores) == 0) return("")
    custom_names <- names(custom_fluorophores)
    if (is.null(custom_names)) return("")

    candidates <- unique(c(
        fn,
        basename(fn),
        tools::file_path_sans_ext(fn),
        tools::file_path_sans_ext(basename(fn))
    ))
    idx <- match(candidates, custom_names)
    idx <- idx[!is.na(idx)]
    if (length(idx) == 0) return("")
    as.character(custom_fluorophores[[idx[1]]])
}

.detect_control_file_alias <- function(stem, alias_map, min_substring_n = 4) {
    if (length(alias_map) == 0) return("")
    stem_norm <- .control_file_normalize_token(stem)
    if (!nzchar(stem_norm)) return("")
    token_norm <- .control_file_split_filename_tokens(stem)
    compounds <- .control_file_build_compound_tokens(token_norm, max_n = 4)
    keys <- names(alias_map)
    keys <- keys[order(-nchar(keys))]

    for (k in keys) {
        if (k == stem_norm) return(unname(alias_map[[k]]))
    }
    for (k in keys) {
        if (k %in% compounds) return(unname(alias_map[[k]]))
    }
    for (k in keys) {
        if (nchar(k) >= min_substring_n && grepl(k, stem_norm, fixed = TRUE)) {
            return(unname(alias_map[[k]]))
        }
    }
    ""
}

.control_file_marker_regex <- function(marker) {
    chars <- strsplit(marker, "", fixed = TRUE)[[1]]
    escaped <- vapply(chars, function(ch) {
        if (grepl("[A-Za-z0-9]", ch)) return(ch)
        if (grepl("\\s", ch)) return("\\s+")
        if (ch %in% c("\\", ".", "^", "$", "|", "?", "*", "+", "(", ")", "[", "]", "{", "}")) return(paste0("\\", ch))
        ch
    }, character(1))
    paste0("(^|[^[:alnum:]])", paste0(escaped, collapse = ""), "([^[:alnum:]]|$)")
}

.detect_control_file_marker <- function(stem, marker_names = character(), marker_name_map = character()) {
    marker_guess <- .detect_control_file_alias(stem, marker_name_map, min_substring_n = 3)
    if (nzchar(marker_guess)) return(marker_guess)

    marker_names <- unique(trimws(as.character(marker_names)))
    marker_names <- marker_names[nzchar(marker_names)]
    if (length(marker_names) > 0) {
        marker_names <- marker_names[order(-nchar(marker_names), marker_names)]
        for (marker in marker_names) {
            if (grepl(.control_file_marker_regex(marker), stem, ignore.case = TRUE, perl = TRUE)) {
                return(marker)
            }
        }
    }
    ""
}

.infer_control_type_from_filename <- function(fn, fluor_guess = "", marker_guess = "") {
    stem <- tools::file_path_sans_ext(basename(fn))
    tok <- .control_file_split_filename_tokens(stem)
    if (any(tok %in% c("bead", "beads", "compbead", "compbeads"))) return("beads")
    if (any(tok %in% c("cell", "cells", "pbmc", "lymphocyte", "lymphocytes", "splenocyte", "splenocytes"))) return("cells")

    fluor_val <- trimws(as.character(fluor_guess))
    marker_val <- trimws(as.character(marker_guess))
    if (grepl("^AF($|_)", fluor_val, ignore.case = TRUE)) return("cells")
    if (grepl("autofluorescence", marker_val, ignore.case = TRUE)) return("cells")
    ""
}

.control_file_is_bead_negative_filename <- function(fn) {
    stem <- tools::file_path_sans_ext(basename(fn))
    stem_norm <- .control_file_normalize_token(stem)
    tok <- .control_file_split_filename_tokens(stem)

    has_bead <- any(tok %in% c("bead", "beads", "compbead", "compbeads")) ||
        grepl("compbeads?", stem_norm) ||
        grepl("beads?", stem_norm)
    has_negative <- any(tok %in% c(
        "unstained", "unstain", "us", "ut", "usut", "usut1", "neg",
        "negative", "background", "bg", "blank", "minus"
    )) ||
        grepl("us[_ -]?ut", stem, ignore.case = TRUE) ||
        grepl("unstained|negative|background|blank", stem, ignore.case = TRUE) ||
        grepl("(^|[^[:alnum:]])(?:us|neg|bg)(?:[^[:alnum:]]|$)", stem, ignore.case = TRUE, perl = TRUE)

    has_bead && has_negative
}

.infer_is_viability_from_filename <- function(fn, fluor_guess = "", marker_guess = "") {
    stem <- tools::file_path_sans_ext(basename(fn))
    stem_norm <- .control_file_normalize_token(stem)
    tok <- .control_file_split_filename_tokens(stem)

    if (any(tok %in% c("live", "dead", "ld", "viability", "viable", "fixviab", "fvd", "fvs", "zombie"))) return("TRUE")
    if (grepl("livedead", stem_norm, fixed = TRUE) ||
        grepl("fixableviability", stem_norm, fixed = TRUE) ||
        grepl("fixviab", stem_norm, fixed = TRUE) ||
        grepl("fvd[0-9]*", stem_norm) ||
        grepl("fvs[0-9]*", stem_norm)) return("TRUE")

    fluor_norm <- .control_file_normalize_token(fluor_guess)
    marker_norm <- .control_file_normalize_token(marker_guess)
    if (any(c(fluor_norm, marker_norm) %in% c(
        "livedeadnir", "livedeadaqua", "livedeadviolet", "livedeaduv",
        "livedeadgreen", "livedeadred", "livedeadyellow", "livedeadblue",
        "livedeadfarred",
        "zombienir", "zombieaqua", "zombieviolet", "zombieuv",
        "zombiegreen", "zombiered", "zombieyellow"
    ))) {
        return("TRUE")
    }
    ""
}

.infer_fluor_from_filename <- function(fn, ref, custom_fluorophores = NULL) {
    stem <- tools::file_path_sans_ext(basename(fn))
    stem_norm <- .control_file_normalize_token(stem)

    custom_fluor <- .control_file_resolve_custom_fluor(fn, custom_fluorophores = custom_fluorophores)
    if (nzchar(custom_fluor)) return(custom_fluor)

    alexa_match <- regexec("(alexa(?:fluor)?|af)([0-9]{3})", stem_norm, perl = TRUE)
    alexa_parts <- regmatches(stem_norm, alexa_match)[[1]]
    if (length(alexa_parts) >= 3) {
        return(paste("Alexa Fluor", alexa_parts[3]))
    }

    if (grepl("livedeadnir", stem_norm, fixed = TRUE) ||
        grepl("fixableviabilitynir", stem_norm, fixed = TRUE) ||
        grepl("livedeadfixablenearir", stem_norm, fixed = TRUE)) {
        return("LIVE/DEAD NIR")
    }
    if (grepl("zombienir", stem_norm, fixed = TRUE)) return("Zombie NIR")
    if (grepl("pecf594", stem_norm, fixed = TRUE)) return("PE-CF594")
    if (grepl("pecy7", stem_norm, fixed = TRUE)) return("PE-Cy7")
    if (grepl("pefire700", stem_norm, fixed = TRUE)) return("PE-Fire 700")

    fluor_guess <- .detect_control_file_alias(stem, ref$fluor_name_map, min_substring_n = 4)
    if (nzchar(fluor_guess)) return(fluor_guess)

    if (grepl("^strep[-_]", stem, ignore.case = TRUE)) {
        strep_map <- c(
            "421" = "BV421", "510" = "BV510", "570" = "BV570", "605" = "BV605",
            "650" = "BV650", "661" = "BV661", "711" = "BV711", "737" = "BV737",
            "750" = "BV750", "785" = "BV785"
        )
        for (k in names(strep_map)) {
            if (grepl(k, stem_norm, fixed = TRUE)) return(strep_map[[k]])
        }
        if (grepl("pecf594", stem_norm, fixed = TRUE)) return("PE-CF594")
        if (grepl("apc", stem_norm, fixed = TRUE)) return("APC")
        if (grepl("pe", stem_norm, fixed = TRUE)) return("PE")
    }

    fn_norm <- .control_file_normalize_token(fn)
    for (p in ref$fluor_patterns) {
        key <- .control_file_normalize_token(p)
        p_norm <- .control_file_normalize_token(p)
        literal_hit <- nzchar(p_norm) && grepl(p_norm, fn_norm, fixed = TRUE)
        normalized_hit <- nzchar(key) && grepl(key, fn_norm, fixed = TRUE)
        if (literal_hit || normalized_hit) {
            if (key %in% names(ref$fluor_name_map)) {
                return(unname(ref$fluor_name_map[[key]]))
            }
            return(p)
        }
    }
    ""
}

.infer_marker_from_filename <- function(fn, fluor_guess = "", marker_name_map = character(), marker_names = character()) {
    stem <- tools::file_path_sans_ext(basename(fn))
    stem_lower <- tolower(stem)
    if (grepl("unstained|us_ut|^af($|[-_])", stem_lower)) return("Autofluorescence")
    marker_guess <- .detect_control_file_marker(stem, marker_names = marker_names, marker_name_map = marker_name_map)
    if (!nzchar(marker_guess)) return("")
    if (nzchar(fluor_guess) && .control_file_normalize_token(marker_guess) == .control_file_normalize_token(fluor_guess)) return("")
    marker_guess
}

.control_file_resolve_detector_channel <- function(channel_name, detector_desc = "", channel_alias_map = character(), fluor_channel_map = character()) {
    keys <- .control_file_canonicalize_channel(c(channel_name, detector_desc))
    if (length(keys) == 0 || all(keys == "")) return("")
    keys <- keys[nzchar(keys)]
    keys <- c(keys, gsub("-A$", "", keys), paste0(keys, "-A"))
    keys <- unique(keys[nzchar(keys)])
    if (length(channel_alias_map) > 0) {
        alias_hits <- channel_alias_map[keys]
        alias_hits <- alias_hits[!is.na(alias_hits) & nzchar(alias_hits)]
        if (length(alias_hits) > 0) keys <- unique(c(keys, alias_hits))
    }
    if (length(fluor_channel_map) > 0) {
        hit <- keys[keys %in% names(fluor_channel_map)]
        if (length(hit) > 0) return(hit[1])
    }
    keys[1]
}

.control_file_resolve_channel_alias <- function(channel_name, channel_alias_map = character()) {
    keys <- .control_file_canonicalize_channel(channel_name)
    keys <- unique(c(keys, gsub("-A$", "", keys), paste0(keys, "-A")))
    keys <- keys[nzchar(keys)]
    if (length(keys) == 0) return("")
    if (length(channel_alias_map) > 0) {
        alias_hits <- channel_alias_map[keys]
        alias_hits <- alias_hits[!is.na(alias_hits) & nzchar(alias_hits)]
        if (length(alias_hits) > 0) keys <- unique(c(keys, alias_hits))
    }
    keys[1]
}

.control_file_get_expected_af_channel <- function(cytometer_name) {
    id <- .resolve_cytometer_id(cytometer_name, allow_auto = TRUE, unknown_as_auto = TRUE)
    switch(id,
        aurora = "B1-A",
        northern_lights = "B1-A",
        id7000 = "405CH7-A",
        discover_s8 = "V6 (515)-A",
        discover_a8 = "V6 (515)-A",
        a5se = "UV515-A",
        opteon = "UV508-A",
        mosaic = "V8-A",
        xenith = "FL13-A",
        ""
    )
}

.control_file_extract_detector_codes <- function(...) {
    text <- toupper(paste(..., collapse = " "))
    if (!nzchar(trimws(text))) return(character())
    matches <- gregexpr("(^|[^A-Z0-9])(UV|YG|[VRB])([0-9]+)(-A)?([^A-Z0-9]|$)", text, perl = TRUE)
    raw <- regmatches(text, matches)[[1]]
    if (length(raw) == 1 && raw[1] == "-1") return(character())
    codes <- sub(".*?(UV|YG|[VRB])([0-9]+)(-A)?.*", "\\1\\2", raw, perl = TRUE)
    unique(codes[nzchar(codes)])
}

.control_file_next_af_tag <- function(existing_vals) {
    ex <- trimws(as.character(existing_vals))
    ex <- ex[!is.na(ex) & nzchar(ex)]
    if (!any(toupper(ex) == "AF")) return("AF")
    k <- 2
    repeat {
        cand <- paste0("AF_", k)
        if (!(cand %in% ex)) return(cand)
        k <- k + 1
    }
}

.control_file_expected_channel_for_fluor <- function(fluorophore, fluor_peak_channel_map = character()) {
    fluorophore <- trimws(as.character(fluorophore)[1])
    if (!nzchar(fluorophore) || length(fluor_peak_channel_map) == 0) {
        return("")
    }
    key <- .control_file_normalize_token(fluorophore)
    hit <- fluor_peak_channel_map[key]
    if (is.na(hit) || !nzchar(hit)) {
        return("")
    }
    as.character(hit)
}

.infer_fluor_from_detector <- function(channel_name,
                                       detector_desc = "",
                                       channel_alias_map = character(),
                                       fluor_channel_map = character()) {
    resolved_channel <- .control_file_resolve_detector_channel(
        channel_name,
        detector_desc = detector_desc,
        channel_alias_map = channel_alias_map,
        fluor_channel_map = fluor_channel_map
    )
    if (length(fluor_channel_map) > 0) {
        direct <- fluor_channel_map[resolved_channel]
        if (!is.na(direct) && nzchar(direct)) return(as.character(direct))
    }

    exact_map <- c(
        "YG1" = "PE",
        "YG3" = "PE-CF594",
        "YG5" = "PE-Cy5",
        "YG7" = "PE-Cy5.5",
        "YG9" = "PE-Cy7",
        "R1" = "APC",
        "R2" = "Alexa Fluor 647",
        "R4" = "Alexa Fluor 700",
        "R7" = "APC-Cy7",
        "V1" = "BV421",
        "V5" = "BV480",
        "V7" = "BV510",
        "V8" = "BV570",
        "V10" = "BV605",
        "V11" = "BV650",
        "V13" = "BV711",
        "V14" = "BV750",
        "V15" = "BV785",
        "B2" = "FITC",
        "B8" = "PerCP",
        "B9" = "PerCP-Cy5.5",
        "UV2" = "BUV395",
        "UV7" = "BUV496",
        "UV9" = "BUV563",
        "UV10" = "BUV615",
        "UV11" = "BUV661",
        "UV14" = "BUV737",
        "UV16" = "BUV805"
    )

    detector_codes <- .control_file_extract_detector_codes(channel_name, detector_desc)
    for (code in detector_codes) {
        if (code %in% names(exact_map)) return(exact_map[[code]])
    }

    code_prefix <- sub("[0-9]+$", "", detector_codes)
    if ("YG" %in% code_prefix) return("PE")
    if ("R" %in% code_prefix) return("APC")
    if ("V" %in% code_prefix) return("BV421")
    if ("B" %in% code_prefix) return("FITC")
    if ("UV" %in% code_prefix) return("BUV395")

    token <- toupper(gsub("\\s+", "", paste(channel_name, detector_desc)))
    if (grepl("YELLOWGREEN", token, fixed = TRUE)) return("PE")
    if (grepl("RED", token, fixed = TRUE)) return("APC")
    if (grepl("VIOLET", token, fixed = TRUE)) return("BV421")
    if (grepl("BLUE", token, fixed = TRUE)) return("FITC")
    if (grepl("ULTRAVIOLET", token, fixed = TRUE)) return("BUV395")
    ""
}
