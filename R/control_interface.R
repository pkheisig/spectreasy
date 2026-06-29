.empty_control_file_df <- function() {
    data.frame(
        filename = character(),
        fluorophore = character(),
        marker = character(),
        channel = character(),
        control.type = character(),
        universal.negative = character(),
        is.viability = character(),
        stringsAsFactors = FALSE
    )
}

.control_file_normalize_token <- function(x) {
    out <- gsub("[^a-z0-9]+", "", tolower(trimws(as.character(x))))
    out[is.na(out)] <- ""
    out
}

.control_file_canonicalize_channel <- function(x) {
    out <- toupper(gsub("\\s+", "", trimws(as.character(x))))
    out <- gsub("([A-Z]+)-([0-9])", "\\1\\2", out, perl = TRUE)
    out[is.na(out)] <- ""
    out
}

.control_file_split_semicolon <- function(x) {
    if (is.null(x) || is.na(x)) return(character())
    vals <- trimws(unlist(strsplit(as.character(x), ";", fixed = TRUE)))
    vals[nzchar(vals)]
}

.control_file_extdata_file <- function(filename) {
    .spectreasy_extdata_file(filename)
}

.control_file_preferred_fluors <- function() {
    c(
        "FITC", "BB515", "Alexa Fluor 488",
        "PE", "PE-CF594", "PE-Dazzle 594", "PE-Fire 640", "PE-Fire 700", "PE-Fire 810",
        "PE-Cy5", "PE-Cy5.5", "PE-Cy7",
        "APC", "Alexa Fluor 647", "Alexa Fluor 700", "APC-R700", "APC-Cy7", "APC-H7", "APC-Fire 750", "APC-Fire 810",
        "PerCP", "PerCP-Cy5.5", "BB700",
        "BV421", "BV480", "BV510", "BV570", "BV605", "BV650", "BV661", "BV711", "BV737", "BV750", "BV785",
        "BUV395", "BUV496", "BUV563", "BUV615", "BUV661", "BUV737", "BUV805",
        "Alexa Fluor 350", "Alexa Fluor 405", "Alexa Fluor 430", "Alexa Fluor 514",
        "Alexa Fluor 532", "Alexa Fluor 546", "Alexa Fluor 555", "Alexa Fluor 568",
        "Alexa Fluor 594", "Alexa Fluor 610", "Alexa Fluor 633", "Alexa Fluor 660",
        "Alexa Fluor 680", "Alexa Fluor 750", "Alexa Fluor 790",
        "RB545", "RB613", "RB667", "RB705", "RB744", "RB780",
        "Spark Blue 550", "Spark NIR 685",
        "Super Bright 436", "Super Bright 600", "Super Bright 645", "Super Bright 702", "Super Bright 780"
    )
}

.control_file_choose_preferred_fluor <- function(candidates, preferred_fluors = .control_file_preferred_fluors()) {
    candidates <- unique(trimws(as.character(candidates)))
    candidates <- candidates[nzchar(candidates)]
    if (length(candidates) == 0) return("")
    pref_idx <- match(tolower(candidates), tolower(preferred_fluors))
    if (any(!is.na(pref_idx))) {
        return(candidates[which.min(ifelse(is.na(pref_idx), Inf, pref_idx))])
    }
    sort(candidates)[1]
}

.control_file_generated_marker_map <- function(max_cd = 371L) {
    marker_names <- c(
        paste0("CD", seq_len(max_cd)),
        paste0("CD1", letters[1:5]),
        "CD45RA", "CD45RO"
    )
    marker_names <- unique(marker_names)
    keys <- .control_file_normalize_token(marker_names)
    keys <- keys[nzchar(keys)]
    stats::setNames(marker_names[seq_along(keys)], keys)
}

.load_control_file_shipped_reference <- function(cytometer_name, preferred_fluors = .control_file_preferred_fluors()) {
    fluor_file <- .control_file_extdata_file("fluorophore_dictionary.csv")
    channel_file <- .control_file_extdata_file("fluorophore_channel_dictionary.csv")
    marker_file <- .control_file_extdata_file("marker_dictionary.csv")
    out <- list(name_map = character(), channel_map = character(), marker_map = character(), marker_names = character())
    cytometer_id <- .resolve_cytometer_id(cytometer_name, allow_auto = TRUE, unknown_as_auto = TRUE)

    if (nzchar(fluor_file)) {
        fluor_df <- tryCatch(
            utils::read.csv(fluor_file, stringsAsFactors = FALSE, check.names = FALSE),
            error = function(e) NULL
        )
        if (!is.null(fluor_df) && nrow(fluor_df) > 0 && "fluorophore" %in% colnames(fluor_df)) {
            for (i in seq_len(nrow(fluor_df))) {
                canonical <- trimws(as.character(fluor_df$fluorophore[i]))
                if (!nzchar(canonical)) next
                aliases_raw <- if ("aliases" %in% colnames(fluor_df)) fluor_df$aliases[i] else ""
                aliases <- unique(c(canonical, .control_file_split_semicolon(aliases_raw)))
                alias_keys <- .control_file_normalize_token(aliases)
                alias_keys <- alias_keys[nzchar(alias_keys)]
                for (k in alias_keys) {
                    if (!k %in% names(out$name_map)) out$name_map[k] <- canonical
                }
            }
        }
    }

    if (nzchar(channel_file)) {
        channel_df <- tryCatch(
            utils::read.csv(channel_file, stringsAsFactors = FALSE, check.names = FALSE),
            error = function(e) NULL
        )
        required <- c("fluorophore", "cytometer", "channel")
        if (!is.null(channel_df) && nrow(channel_df) > 0 && all(required %in% colnames(channel_df))) {
            if (!identical(cytometer_id, "auto")) {
                channel_ids <- vapply(
                    channel_df$cytometer,
                    .resolve_cytometer_id,
                    character(1),
                    allow_auto = FALSE,
                    unknown_as_auto = FALSE
                )
                channel_df <- channel_df[channel_ids == cytometer_id, , drop = FALSE]
            }
            for (i in seq_len(nrow(channel_df))) {
                canonical <- trimws(as.character(channel_df$fluorophore[i]))
                if (!nzchar(canonical)) next
                raw_channels <- .control_file_split_semicolon(channel_df$channel[i])
                ch_vals <- .control_file_canonicalize_channel(raw_channels)
                ch_vals <- unique(c(ch_vals, gsub("-A$", "", ch_vals), paste0(gsub("-A$", "", ch_vals), "-A")))
                ch_vals <- ch_vals[nzchar(ch_vals)]
                for (ch in ch_vals) {
                    if (!ch %in% names(out$channel_map)) {
                        out$channel_map[ch] <- canonical
                    } else {
                        out$channel_map[ch] <- .control_file_choose_preferred_fluor(
                            c(out$channel_map[ch], canonical),
                            preferred_fluors = preferred_fluors
                        )
                    }
                }
            }
        }
    }

    if (nzchar(marker_file)) {
        marker_df <- tryCatch(
            utils::read.csv(marker_file, stringsAsFactors = FALSE, check.names = FALSE),
            error = function(e) NULL
        )
        if (!is.null(marker_df) && nrow(marker_df) > 0 && "marker" %in% colnames(marker_df)) {
            for (i in seq_len(nrow(marker_df))) {
                marker <- trimws(as.character(marker_df$marker[i]))
                if (!nzchar(marker)) next
                out$marker_names <- c(out$marker_names, marker)
                aliases_raw <- if ("aliases" %in% colnames(marker_df)) marker_df$aliases[i] else ""
                aliases <- unique(c(marker, .control_file_split_semicolon(aliases_raw)))
                alias_keys <- .control_file_normalize_token(aliases)
                alias_keys <- alias_keys[nzchar(alias_keys)]
                for (k in alias_keys) {
                    if (!k %in% names(out$marker_map)) out$marker_map[k] <- marker
                }
            }
        }
    }

    generated_marker_map <- .control_file_generated_marker_map()
    out$marker_names <- unique(c(out$marker_names, unname(generated_marker_map)))
    generated_marker_map <- generated_marker_map[!(names(generated_marker_map) %in% names(out$marker_map))]
    out$marker_map <- c(out$marker_map, generated_marker_map)

    out
}

.prepare_control_file_reference <- function(cytometer) {
    patterns_list <- get_fluorophore_patterns()
    specific_dyes <- c(patterns_list$beads, patterns_list$cells)
    generic_terms <- c("LIVE/DEAD", "Fixable Viability", "Viability", "FVD")
    specific_dyes <- setdiff(specific_dyes, generic_terms)
    specific_dyes <- specific_dyes[order(-nchar(specific_dyes))]
    fluor_patterns <- c(specific_dyes, generic_terms)
    preferred_fluors <- .control_file_preferred_fluors()
    shipped_ref <- .load_control_file_shipped_reference(cytometer, preferred_fluors = preferred_fluors)

    list(
        fluor_patterns = fluor_patterns,
        sample_patterns = patterns_list,
        fluor_name_map = shipped_ref$name_map,
        fluor_channel_map = shipped_ref$channel_map,
        marker_name_map = shipped_ref$marker_map,
        marker_names = shipped_ref$marker_names,
        preferred_fluors = preferred_fluors
    )
}

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

    marker_guess <- .detect_control_file_alias(stem, marker_name_map, min_substring_n = 3)
    if (nzchar(marker_guess)) return(marker_guess)

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

.control_file_is_cell_af_filename <- function(fn) {
    !.control_file_is_bead_negative_filename(fn) && .is_af_filename(fn)
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
    if (.control_file_is_bead_negative_filename(fn)) return("Bead background")
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
    k <- 1
    repeat {
        cand <- paste0("AF_", k)
        if (!(cand %in% ex)) return(cand)
        k <- k + 1
    }
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

.control_file_row_from_scc <- function(fn, ref, custom_fluorophores = NULL) {
    fluor <- .infer_fluor_from_filename(fn, ref = ref, custom_fluorophores = custom_fluorophores)
    is_bead_negative <- .control_file_is_bead_negative_filename(fn)
    if (is_bead_negative) {
        fluor <- "AF_Bead"
    } else if (.control_file_is_cell_af_filename(fn)) {
        fluor <- "AF_Internal"
    }
    marker <- .infer_marker_from_filename(
        fn,
        fluor_guess = fluor,
        marker_name_map = ref$marker_name_map,
        marker_names = ref$marker_names
    )
    control_type <- .infer_control_type_from_filename(fn, fluor_guess = fluor, marker_guess = marker)
    viability_flag <- .infer_is_viability_from_filename(fn, fluor_guess = fluor, marker_guess = marker)

    data.frame(
        filename = fn,
        fluorophore = fluor,
        marker = marker,
        channel = "",
        control.type = control_type,
        universal.negative = "",
        is.viability = if (is_bead_negative) "" else viability_flag,
        stringsAsFactors = FALSE
    )
}

.build_control_file_scc_df <- function(scc_files, ref, custom_fluorophores = NULL) {
    if (length(scc_files) == 0) {
        return(.empty_control_file_df())
    }
    rows <- lapply(scc_files, .control_file_row_from_scc, ref = ref, custom_fluorophores = custom_fluorophores)
    do.call(rbind, rows)
}

.build_control_file_af_df <- function(af_files) {
    if (length(af_files) == 0) {
        return(.empty_control_file_df())
    }

    rows <- lapply(seq_along(af_files), function(i) {
        fn <- af_files[i]
        tag <- if (i == 1) "AF" else paste0("AF_", i)
        data.frame(
            filename = fn,
            fluorophore = tag,
            marker = "Autofluorescence",
            channel = "",
            control.type = "cells",
            universal.negative = "",
            is.viability = "",
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, rows)
}

.combine_control_file_seed_rows <- function(scc_df, af_df) {
    df <- do.call(rbind, list(af_df, scc_df))
    df[!duplicated(df$filename), , drop = FALSE]
}

.control_file_row_path <- function(fn, input_folder, af_folder, include_af_folder = TRUE) {
    if (isTRUE(include_af_folder) && file.exists(file.path(af_folder, fn))) {
        return(file.path(af_folder, fn))
    }
    file.path(input_folder, fn)
}

.summarize_control_file_fcs <- function(path, channel_alias_map = character()) {
    pd <- NULL
    peak_channel <- ""
    top_channels <- character()
    top_ratio <- NA_real_
    entropy_hi <- NA_real_
    read_error <- ""

    tryCatch({
        ff <- suppressWarnings(flowCore::read.FCS(path, transformation = FALSE, truncate_max_range = FALSE))
        pd <- flowCore::pData(flowCore::parameters(ff))
        channel_alias_map <- .merge_control_file_alias_map(channel_alias_map, .build_channel_alias_map_from_pd(pd))
        fl_pd <- get_sorted_detectors(pd)
        if (length(fl_pd$names) > 0) {
            expr <- flowCore::exprs(ff)[, fl_pd$names, drop = FALSE]

            use_idx <- seq_len(nrow(expr))
            if ("name" %in% colnames(pd)) {
                fsc <- pd$name[grepl("^FSC", pd$name) & grepl("-A$", pd$name)][1]
                ssc <- pd$name[grepl("^SSC", pd$name) & grepl("-A$", pd$name)][1]
                if (!is.na(fsc) && !is.na(ssc) && fsc %in% colnames(flowCore::exprs(ff)) && ssc %in% colnames(flowCore::exprs(ff))) {
                    scatter <- flowCore::exprs(ff)[, c(fsc, ssc), drop = FALSE]
                    idx <- which(
                        scatter[, 1] > quantile(scatter[, 1], 0.10, na.rm = TRUE) &
                        scatter[, 1] < quantile(scatter[, 1], 0.90, na.rm = TRUE) &
                        scatter[, 2] > quantile(scatter[, 2], 0.10, na.rm = TRUE) &
                        scatter[, 2] < quantile(scatter[, 2], 0.90, na.rm = TRUE)
                    )
                    if (length(idx) >= 200) use_idx <- idx
                }
            }
            expr_use <- expr[use_idx, , drop = FALSE]
            q_hi <- apply(expr_use, 2, function(v) quantile(v, 0.999, na.rm = TRUE))
            q_hi[!is.finite(q_hi)] <- NA_real_
            peak_channel <- if (all(is.na(q_hi))) {
                medians <- apply(expr_use, 2, median, na.rm = TRUE)
                names(medians)[which.max(medians)]
            } else {
                names(q_hi)[which.max(q_hi)]
            }

            q_pos <- q_hi[is.finite(q_hi) & q_hi > 0]
            if (length(q_pos) > 0) {
                ord <- order(q_pos, decreasing = TRUE)
                top_channels <- names(q_pos)[ord[seq_len(min(3, length(ord)))]]
                p <- q_pos / sum(q_pos)
                entropy_hi <- -sum(p * log(p))
                if (length(ord) >= 2) {
                    top_ratio <- as.numeric(q_pos[ord[1]] / pmax(q_pos[ord[2]], 1e-9))
                }
            }
        }
    }, error = function(e) {
        read_error <<- conditionMessage(e)
    })

    list(
        pd = pd,
        peak_channel = peak_channel,
        top_channels = top_channels,
        top_ratio = top_ratio,
        entropy_hi = entropy_hi,
        read_error = read_error,
        channel_alias_map = channel_alias_map
    )
}

.control_file_should_warn_read_error <- function(path) {
    size <- suppressWarnings(file.info(path)$size)
    is.finite(size) && !is.na(size) && size > 0
}

.control_file_is_af_like <- function(current_fluor,
                                     current_marker,
                                     top_channels,
                                     top_ratio,
                                     entropy_hi,
                                     expected_af_channel,
                                     channel_alias_map = character()) {
    unresolved_label <- (is.na(current_fluor) || current_fluor == "" || identical(current_fluor, "Unknown")) &&
        (is.na(current_marker) || current_marker == "")
    if (!unresolved_label || !nzchar(expected_af_channel) || length(top_channels) == 0) {
        return(FALSE)
    }

    resolved_top <- unique(vapply(top_channels, .control_file_resolve_channel_alias, character(1), channel_alias_map = channel_alias_map))
    af_in_top <- expected_af_channel %in% resolved_top
    broad_signal <- is.finite(top_ratio) && top_ratio <= 1.7
    diffuse_signal <- is.finite(entropy_hi) && entropy_hi >= 3.2
    af_in_top && broad_signal && diffuse_signal
}

.annotate_control_file_rows <- function(df,
                                        input_folder,
                                        af_folder,
                                        include_af_folder,
                                        cytometer,
                                        unknown_fluor_policy,
                                        ref) {
    channel_alias_map <- character()
    expected_af_channel <- .control_file_resolve_channel_alias(
        .control_file_get_expected_af_channel(cytometer),
        channel_alias_map = channel_alias_map
    )
    auto_af_files <- character(0)

    for (i in seq_len(nrow(df))) {
        fn <- df$filename[i]
        path <- .control_file_row_path(fn, input_folder = input_folder, af_folder = af_folder, include_af_folder = include_af_folder)
        summary <- .summarize_control_file_fcs(path, channel_alias_map = channel_alias_map)
        channel_alias_map <- summary$channel_alias_map

        if (nzchar(summary$read_error) && .control_file_should_warn_read_error(path)) {
            warning(
                "Could not read FCS file while auto-detecting peak channel: ",
                path,
                "\n",
                summary$read_error,
                call. = FALSE
            )
        }

        if (nzchar(summary$peak_channel)) {
            df$channel[i] <- summary$peak_channel
        }

        current_fluor <- trimws(as.character(df$fluorophore[i]))
        current_marker <- trimws(as.character(df$marker[i]))
        af_like <- .control_file_is_af_like(
            current_fluor = current_fluor,
            current_marker = current_marker,
            top_channels = summary$top_channels,
            top_ratio = summary$top_ratio,
            entropy_hi = summary$entropy_hi,
            expected_af_channel = expected_af_channel,
            channel_alias_map = channel_alias_map
        )
        if (af_like) {
            new_tag <- .control_file_next_af_tag(df$fluorophore)
        df$fluorophore[i] <- new_tag
        df$marker[i] <- "Autofluorescence"
        df$control.type[i] <- "cells"
        auto_af_files <- c(auto_af_files, fn)
        }

        current_fluor <- trimws(as.character(df$fluorophore[i]))
        unresolved <- is.na(current_fluor) || current_fluor == "" || identical(current_fluor, "Unknown")
        if (unresolved) {
            if (unknown_fluor_policy == "filename") {
                df$fluorophore[i] <- tools::file_path_sans_ext(fn)
            } else if (unknown_fluor_policy == "by_channel") {
                detector_desc <- ""
                if (!is.null(summary$pd) && "desc" %in% colnames(summary$pd) && "name" %in% colnames(summary$pd)) {
                    idx <- match(df$channel[i], as.character(summary$pd$name))
                    if (!is.na(idx)) {
                        detector_desc <- as.character(summary$pd$desc[idx])
                    }
                }
                df$fluorophore[i] <- .infer_fluor_from_detector(
                    df$channel[i],
                    detector_desc = detector_desc,
                    channel_alias_map = channel_alias_map,
                    fluor_channel_map = ref$fluor_channel_map
                )
            } else {
                df$fluorophore[i] <- ""
            }
        }

        current_marker <- trimws(as.character(df$marker[i]))
        if (is.na(current_marker) || current_marker == "") {
            df$marker[i] <- .infer_marker_from_filename(
                fn,
                fluor_guess = df$fluorophore[i],
                marker_name_map = ref$marker_name_map,
                marker_names = ref$marker_names
            )
        }

        current_type <- trimws(as.character(df$control.type[i]))
        if (is.na(current_type) || current_type == "") {
            df$control.type[i] <- .infer_control_type_from_filename(
                fn,
                fluor_guess = df$fluorophore[i],
                marker_guess = df$marker[i]
            )
        }

        current_viability <- toupper(trimws(as.character(df$is.viability[i])))
        if (is.na(current_viability) || current_viability == "") {
            df$is.viability[i] <- .infer_is_viability_from_filename(
                fn,
                fluor_guess = df$fluorophore[i],
                marker_guess = df$marker[i]
            )
        }
    }

    list(df = df, auto_af_files = auto_af_files)
}

.finalize_control_file_df <- function(df, scc_files, af_files, auto_af_files = character()) {
    bead_negative_files <- scc_files[vapply(scc_files, .control_file_is_bead_negative_filename, logical(1))]
    primary_af_candidates <- scc_files[vapply(scc_files, .control_file_is_cell_af_filename, logical(1))]
    auto_af_files <- auto_af_files[!vapply(auto_af_files, .control_file_is_bead_negative_filename, logical(1))]
    auto_af_primary <- if (length(auto_af_files) > 0) auto_af_files[1] else ""
    primary_af_file <- if (length(af_files) > 0) {
        af_files[1]
    } else if (length(primary_af_candidates) > 0) {
        primary_af_candidates[1]
    } else if (nzchar(auto_af_primary)) {
        auto_af_primary
    } else {
        ""
    }

    df$universal.negative <- ""

    df$control.type <- tolower(trimws(as.character(df$control.type)))
    df$control.type[is.na(df$control.type)] <- ""
    invalid_type <- !(df$control.type %in% c("", "beads", "cells"))
    df$control.type[invalid_type] <- ""

    df$is.viability <- toupper(trimws(as.character(df$is.viability)))
    df$is.viability[is.na(df$is.viability)] <- ""
    df$is.viability[df$is.viability %in% c("T", "TRUE", "1", "YES", "Y")] <- "TRUE"
    df$is.viability[!(df$is.viability %in% c("", "TRUE"))] <- ""
    viability_missing_marker <- df$is.viability == "TRUE" & !nzchar(trimws(as.character(df$marker)))
    df$marker[viability_missing_marker] <- "Live"

    bead_negative_row <- df$filename %in% bead_negative_files |
        (tolower(df$control.type) == "beads" & grepl("^AF_Bead$", as.character(df$fluorophore), ignore.case = TRUE))
    if (any(bead_negative_row)) {
        df$fluorophore[bead_negative_row] <- "AF_Bead"
        df$marker[bead_negative_row] <- "Bead background"
        df$control.type[bead_negative_row] <- "beads"
        df$is.viability[bead_negative_row] <- ""
    }

    is_af_row <- grepl("^AF($|_)", as.character(df$fluorophore), ignore.case = TRUE) & !bead_negative_row
    df$control.type[is_af_row] <- "cells"

    if (primary_af_file != "") {
        df$fluorophore[df$filename == primary_af_file] <- "AF"
        df$marker[df$filename == primary_af_file] <- "Autofluorescence"
        df$control.type[df$filename == primary_af_file] <- "cells"
    }

    preferred_order <- c(
        "filename",
        "fluorophore",
        "marker",
        "channel",
        "control.type",
        "universal.negative",
        "is.viability"
    )
    keep <- intersect(preferred_order, colnames(df))
    df[, c(keep, setdiff(colnames(df), keep)), drop = FALSE]
}

#' Create spectreasy Control File
#' 
#' Generates a spectreasy-compatible control CSV.
#' 
#' @param input_folder Directory containing single-stained control FCS files.
#' @param af_folder Optional AF folder.
#' @param include_af_folder Logical. Include AF folder files in the control file.
#' @param cytometer Cytometer name used as a channel-mapping hint. The default,
#'   `"auto"`, infers the cytometer from FCS detector names when possible.
#' @param default_control_type Deprecated. Kept for backward compatibility and ignored.
#' @param unknown_fluor_policy How to fill unresolved fluorophores:
#'   `"empty"` (recommended), `"by_channel"` (best-effort guess), `"filename"`.
#' @param output_file Path where the CSV will be saved (default: "fcs_mapping.csv").
#' @param custom_fluorophores Optional named vector to map filenames to fluorophore names.
#' @return A data frame containing the control file information.
#'
#' The generated file intentionally leaves `universal.negative` empty by default.
#' It auto-detects `control.type` from filename tokens (for example `"beads"`
#' or `"cells"`). Unstained/negative bead files are denoted as `AF_Bead`
#' with `control.type = "beads"` so they can be used as bead-background
#' negatives without being mixed into the cellular AF bank.
#' @export
#' @examples
#' make_example_ff <- function(main, n = 250) {
#'   exprs <- cbind(
#'     "B2-A" = pmax(rnorm(n, main[1], 40), 1),
#'     "YG1-A" = pmax(rnorm(n, main[2], 15), 1),
#'     "R1-A" = pmax(rnorm(n, main[3], 10), 1),
#'     "FSC-A" = rnorm(n, 90000, 7000),
#'     "SSC-A" = rnorm(n, 45000, 5000),
#'     Time = seq_len(n)
#'   )
#'   flowCore::flowFrame(exprs)
#' }
#'
#' td <- tempfile("spectreasy-")
#' dir.create(td)
#' scc_dir <- file.path(td, "scc")
#' dir.create(scc_dir)
#' flowCore::write.FCS(make_example_ff(c(800, 80, 50)), file.path(scc_dir, "FITC_cells.fcs"))
#' flowCore::write.FCS(make_example_ff(c(80, 820, 60)), file.path(scc_dir, "PE_cells.fcs"))
#'
#' control_df <- create_control_file(
#'   input_folder = scc_dir,
#'   include_af_folder = FALSE,
#'   cytometer = "auto",
#'   output_file = file.path(td, "fcs_mapping.csv")
#' )
#' head(control_df)
create_control_file <- function(input_folder = "scc",
                                af_folder = "af",
                                include_af_folder = TRUE,
                                cytometer = "auto",
                                default_control_type = "cells",
                                unknown_fluor_policy = c("empty", "by_channel", "filename"),
                                output_file = "fcs_mapping.csv",
                                custom_fluorophores = NULL) {
    unknown_fluor_policy <- match.arg(unknown_fluor_policy)
    scc_files <- list.files(input_folder, pattern = "\\.fcs$", full.names = FALSE, ignore.case = TRUE)
    af_files <- if (include_af_folder && dir.exists(af_folder)) {
        list.files(af_folder, pattern = "\\.fcs$", full.names = FALSE, ignore.case = TRUE)
    } else {
        character(0)
    }

    if (length(scc_files) == 0) stop("No FCS files found in ", input_folder)

    cytometer_resolved <- .resolve_cytometer_from_files(
        cytometer,
        files = c(file.path(input_folder, scc_files), file.path(af_folder, af_files))
    )
    ref <- .prepare_control_file_reference(cytometer_resolved)
    scc_df <- .build_control_file_scc_df(scc_files, ref = ref, custom_fluorophores = custom_fluorophores)
    af_df <- .build_control_file_af_df(af_files)
    df <- .combine_control_file_seed_rows(scc_df = scc_df, af_df = af_df)

    annotated <- .annotate_control_file_rows(
        df = df,
        input_folder = input_folder,
        af_folder = af_folder,
        include_af_folder = include_af_folder,
        cytometer = cytometer_resolved,
        unknown_fluor_policy = unknown_fluor_policy,
        ref = ref
    )
    df <- .finalize_control_file_df(
        df = annotated$df,
        scc_files = scc_files,
        af_files = af_files,
        auto_af_files = annotated$auto_af_files
    )

    utils::write.csv(df, output_file, row.names = FALSE, quote = TRUE)
    df
}

#' Get Spectra via Internal Backend (Robust Multi-AF)
#'
#' Extracts SCC signatures using internal logic and optional AF signatures from
#' the `af/` folder, then combines them for downstream unmixing.
#'
#' @param flow_frame A `flowFrame` used to determine detector ordering.
#' @param control_file Path to spectreasy-compatible control CSV.
#' @param control_dir Directory containing SCC FCS files.
#' @param af_dir Directory containing optional AF FCS files.
#' @param method Reserved for future method selection.
#' @param cytometer Cytometer name used as a channel-mapping hint. The default,
#'   `"auto"`, infers the cytometer from FCS detector names when possible.
#' @return Expanded reference matrix aligned to the detectors in `flow_frame`.
#' @export
#' @examples
#' if (interactive()) {
#'   ff <- flowCore::read.FCS("samples/Sample1.fcs", transformation = FALSE)
#'   M_ctrl <- get_control_spectra(
#'     flow_frame = ff,
#'     control_file = "fcs_mapping.csv",
#'     control_dir = "scc",
#'     cytometer = "auto"
#'   )
#'   dim(M_ctrl)
#' }
get_control_spectra <- function(flow_frame,
                                control_file = "fcs_mapping.csv",
                                control_dir = "scc",
                                af_dir = "af",
                                method = "AutoSpectral",
                                cytometer = "auto") {
    # 1. Get detector info
    pd <- flowCore::pData(flowCore::parameters(flow_frame))
    det_info <- get_sorted_detectors(pd)
    detector_names <- det_info$names
    cytometer_resolved <- .resolve_cytometer_from_pd(cytometer, pd)

    control_file <- .resolve_control_file_path(control_file)

    # 2. Extract SCC signatures
    message("  - Extracting reference signatures from single-color controls...")
    control_df <- utils::read.csv(control_file, stringsAsFactors = FALSE, check.names = FALSE)
    M_scc <- build_reference_matrix(input_folder = control_dir, control_df = control_df, cytometer = cytometer_resolved)
    M_scc <- M_scc[, detector_names, drop = FALSE]

    # 3. Extract Multi-AF signatures
    if (dir.exists(af_dir)) {
        message("  - Extracting additional AF signatures from '", af_dir, "' folder...")
        M_af <- extract_af_signatures(af_dir, detector_names)
        M_expanded <- rbind(M_scc, M_af)
    } else {
        M_expanded <- M_scc
    }
    
    return(M_expanded[, detector_names, drop = FALSE])
}

#' Internal Helper: Extract Gated AF Signatures
#' @noRd
extract_af_signatures <- function(af_dir, detector_names) {
    af_files <- list.files(af_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
    if (length(af_files) == 0) return(matrix(0, nrow = 0, ncol = length(detector_names)))
    
    spectra <- list()
    for (f in af_files) {
        ff <- flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)
        raw_data <- flowCore::exprs(ff)
        pd <- flowCore::pData(flowCore::parameters(ff))
        fsc <- pd$name[grepl("^FSC", pd$name) & grepl("-A$", pd$name)][1]
        ssc <- pd$name[grepl("^SSC", pd$name) & grepl("-A$", pd$name)][1]
        
        # 20-80% quantile gate on FSC/SSC to isolate main population
        idx <- which(raw_data[, fsc] > quantile(raw_data[, fsc], 0.2) &
                     raw_data[, fsc] < quantile(raw_data[, fsc], 0.8) &
                     raw_data[, ssc] > quantile(raw_data[, ssc], 0.2) &
                     raw_data[, ssc] < quantile(raw_data[, ssc], 0.8))
        
        if (length(idx) < 10) idx <- seq_len(nrow(raw_data)) 
        
        sig <- apply(raw_data[idx, detector_names, drop = FALSE], 2, median)
        sig_norm <- sig / max(sig)
        
        tag <- if(length(spectra) == 0) "AF" else paste0("AF_", length(spectra) + 1)
        spectra[[tag]] <- sig_norm
    }
    return(do.call(rbind, spectra))
}
