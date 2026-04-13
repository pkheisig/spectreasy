#' Create spectreasy Control File
#' 
#' Generates a spectreasy-compatible control CSV.
#' 
#' @param input_folder Directory containing single-stained control FCS files.
#' @param af_folder Optional AF folder.
#' @param include_af_folder Logical. Include AF folder files in the control file.
#' @param cytometer Cytometer name used for channel-to-fluorophore mapping
#'   (for example `"Aurora"`).
#' @param default_control_type Deprecated. Kept for backward compatibility and ignored.
#' @param unknown_fluor_policy How to fill unresolved fluorophores:
#'   `"empty"` (recommended), `"by_channel"` (best-effort guess), `"filename"`.
#' @param output_file Path where the CSV will be saved (default: "fcs_mapping.csv").
#' @param custom_fluorophores Optional named vector to map filenames to fluorophore names.
#' @return A data frame containing the control file information.
#'
#' The generated file intentionally leaves `universal.negative` empty by default.
#' It auto-detects `control.type` from filename tokens (for example `"beads"`
#' or `"cells"`), and forces AF rows to `cells`.
#' @export
#' @examples
#' \dontrun{
#' control_df <- create_control_file(
#'   input_folder = "scc",
#'   cytometer = "Aurora",
#'   output_file = "fcs_mapping.csv"
#' )
#' head(control_df)
#' }
create_control_file <- function(input_folder = "scc",
                                af_folder = "af",
                                include_af_folder = TRUE,
                                cytometer = "Aurora",
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
    
    # 1. Identify patterns and lookup dictionaries
    patterns_list <- get_fluorophore_patterns()
    specific_dyes <- c(patterns_list$beads, patterns_list$cells)
    generic_terms <- c("LIVE/DEAD", "Fixable Viability", "Viability", "FVD")
    specific_dyes <- setdiff(specific_dyes, generic_terms)
    specific_dyes <- specific_dyes[order(-nchar(specific_dyes))]
    fluor_patterns <- c(specific_dyes, generic_terms)
    sample_patterns <- get_fluorophore_patterns()
    preferred_fluors <- c(
        "FITC", "PE", "PE-CF594", "PE-Cy5", "PE-Cy5.5", "PE-Cy7",
        "APC", "Alexa Fluor 700", "APC-R700", "APC-Cy7",
        "PerCP", "PerCP-Cy5.5",
        "BV421", "BV480", "BV510", "BV570", "BV605", "BV650", "BV661", "BV711", "BV737", "BV750", "BV785",
        "BUV395", "BUV496", "BUV563", "BUV615", "BUV661", "BUV737", "BUV805",
        "Alexa Fluor 488", "Alexa Fluor 594", "Alexa Fluor 647", "Alexa Fluor 700"
    )

    normalize_token <- function(x) {
        out <- gsub("[^a-z0-9]+", "", tolower(trimws(as.character(x))))
        out[is.na(out)] <- ""
        out
    }

    canonicalize_channel <- function(x) {
        out <- toupper(gsub("\\s+", "", trimws(as.character(x))))
        out[is.na(out)] <- ""
        out
    }

    split_semicolon <- function(x) {
        if (is.null(x) || is.na(x)) return(character())
        vals <- trimws(unlist(strsplit(as.character(x), ";", fixed = TRUE)))
        vals[nzchar(vals)]
    }

    extdata_file <- function(filename) {
        p <- system.file("extdata", filename, package = "spectreasy")
        if (nzchar(p) && file.exists(p)) return(p)
        local_p <- file.path("inst", "extdata", filename)
        if (file.exists(local_p)) return(local_p)
        ""
    }

    choose_preferred_fluor <- function(candidates) {
        candidates <- unique(trimws(as.character(candidates)))
        candidates <- candidates[nzchar(candidates)]
        if (length(candidates) == 0) return("")
        pref_idx <- match(tolower(candidates), tolower(preferred_fluors))
        if (any(!is.na(pref_idx))) {
            return(candidates[which.min(ifelse(is.na(pref_idx), Inf, pref_idx))])
        }
        sort(candidates)[1]
    }

    load_shipped_reference <- function(cytometer_name) {
        fluor_file <- extdata_file("fluorophore_dictionary.csv")
        marker_file <- extdata_file("marker_dictionary.csv")
        out <- list(name_map = character(), channel_map = character(), marker_map = character())

        if (nzchar(fluor_file)) {
            fluor_df <- tryCatch(
                utils::read.csv(fluor_file, stringsAsFactors = FALSE, check.names = FALSE),
                error = function(e) NULL
            )
            if (!is.null(fluor_df) && nrow(fluor_df) > 0 && "fluorophore" %in% colnames(fluor_df)) {
                channel_col <- colnames(fluor_df)[tolower(colnames(fluor_df)) == tolower(paste0("channel_", cytometer_name))]
                if (length(channel_col) == 0) {
                    channel_col <- colnames(fluor_df)[tolower(colnames(fluor_df)) == "channel_aurora"]
                }
                for (i in seq_len(nrow(fluor_df))) {
                    canonical <- trimws(as.character(fluor_df$fluorophore[i]))
                    if (!nzchar(canonical)) next
                    aliases <- unique(c(canonical, split_semicolon(fluor_df$aliases[i])))
                    alias_keys <- normalize_token(aliases)
                    alias_keys <- alias_keys[nzchar(alias_keys)]
                    for (k in alias_keys) {
                        if (!k %in% names(out$name_map)) out$name_map[k] <- canonical
                    }
                    if (length(channel_col) > 0) {
                        ch_vals <- canonicalize_channel(split_semicolon(fluor_df[i, channel_col[1]]))
                        ch_vals <- ch_vals[nzchar(ch_vals)]
                        for (ch in ch_vals) {
                            if (!ch %in% names(out$channel_map)) {
                                out$channel_map[ch] <- canonical
                            } else {
                                out$channel_map[ch] <- choose_preferred_fluor(c(out$channel_map[ch], canonical))
                            }
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
                    aliases <- unique(c(marker, split_semicolon(marker_df$aliases[i])))
                    alias_keys <- normalize_token(aliases)
                    alias_keys <- alias_keys[nzchar(alias_keys)]
                    for (k in alias_keys) {
                        if (!k %in% names(out$marker_map)) out$marker_map[k] <- marker
                    }
                }
            }
        }

        out
    }

    shipped_ref <- load_shipped_reference(cytometer)
    fluor_name_map <- shipped_ref$name_map
    fluor_channel_map <- shipped_ref$channel_map
    marker_name_map <- shipped_ref$marker_map
    channel_alias_map <- character()

    merge_alias_map <- function(existing, incoming) {
        if (length(incoming) == 0) return(existing)
        incoming <- incoming[!is.na(names(incoming)) & nzchar(names(incoming))]
        incoming <- incoming[!(names(incoming) %in% names(existing))]
        c(existing, incoming)
    }

    split_filename_tokens <- function(text) {
        parts <- unlist(strsplit(as.character(text), "[^A-Za-z0-9]+"))
        parts <- normalize_token(parts)
        parts[nzchar(parts)]
    }

    build_compound_tokens <- function(tokens, max_n = 4) {
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

    detect_alias <- function(stem, alias_map, min_substring_n = 4) {
        if (length(alias_map) == 0) return("")
        stem_norm <- normalize_token(stem)
        if (!nzchar(stem_norm)) return("")
        token_norm <- split_filename_tokens(stem)
        compounds <- build_compound_tokens(token_norm, max_n = 4)
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

    contains_literal <- function(text, patterns) {
        text_up <- toupper(as.character(text))
        any(vapply(
            patterns,
            function(p) {
                p_up <- toupper(as.character(p))
                nzchar(p_up) && grepl(p_up, text_up, fixed = TRUE)
            },
            logical(1)
        ))
    }

    infer_control_type_from_filename <- function(fn, fluor_guess = "", marker_guess = "") {
        stem <- tools::file_path_sans_ext(basename(fn))
        tok <- split_filename_tokens(stem)
        if (any(tok %in% c("bead", "beads", "compbead", "compbeads"))) return("beads")
        if (any(tok %in% c("cell", "cells", "pbmc", "lymphocyte", "lymphocytes", "splenocyte", "splenocytes"))) return("cells")

        fluor_val <- trimws(as.character(fluor_guess))
        marker_val <- trimws(as.character(marker_guess))
        if (grepl("^AF($|_)", fluor_val, ignore.case = TRUE)) return("cells")
        if (grepl("autofluorescence", marker_val, ignore.case = TRUE)) return("cells")
        ""
    }

    infer_is_viability_from_filename <- function(fn, fluor_guess = "", marker_guess = "") {
        stem <- tools::file_path_sans_ext(basename(fn))
        stem_norm <- normalize_token(stem)
        tok <- split_filename_tokens(stem)

        # Explicit user-facing rule: if live/dead/viability is in the filename, mark viability.
        if (any(tok %in% c("live", "dead", "viability"))) return("TRUE")
        if (grepl("livedead", stem_norm, fixed = TRUE) ||
            grepl("fixableviability", stem_norm, fixed = TRUE)) return("TRUE")

        # Common viability dye families.
        fluor_norm <- normalize_token(fluor_guess)
        marker_norm <- normalize_token(marker_guess)
        if (any(c(fluor_norm, marker_norm) %in% c("zombienir", "zombieaqua", "zombieviolet", "zombieuv", "zombiegreen", "zombiered", "zombieyellow"))) {
            return("TRUE")
        }
        ""
    }

    infer_fluor_from_filename <- function(fn) {
        stem <- tools::file_path_sans_ext(basename(fn))
        stem_norm <- normalize_token(stem)

        if (!is.null(custom_fluorophores) && fn %in% names(custom_fluorophores)) {
            return(as.character(custom_fluorophores[fn]))
        }

        # Hard-map frequent filename variants before generic matching.
        if (grepl("livedeadnir", stem_norm, fixed = TRUE) ||
            grepl("fixableviabilitynir", stem_norm, fixed = TRUE) ||
            grepl("livedeadfixablenearir", stem_norm, fixed = TRUE)) {
            return("LIVE/DEAD NIR")
        }
        if (grepl("zombienir", stem_norm, fixed = TRUE)) return("Zombie NIR")
        if (grepl("pecf594", stem_norm, fixed = TRUE)) return("PE-CF594")
        if (grepl("pecy7", stem_norm, fixed = TRUE)) return("PE-Cy7")
        if (grepl("pefire700", stem_norm, fixed = TRUE)) return("PE-Fire 700")

        fluor_guess <- detect_alias(stem, fluor_name_map, min_substring_n = 4)
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

        fn_norm <- normalize_token(fn)
        for (p in fluor_patterns) {
            key <- normalize_token(p)
            p_norm <- normalize_token(p)
            literal_hit <- nzchar(p_norm) && grepl(p_norm, fn_norm, fixed = TRUE)
            normalized_hit <- nzchar(key) && grepl(key, fn_norm, fixed = TRUE)
            if (literal_hit || normalized_hit) {
                if (key %in% names(fluor_name_map)) {
                    return(unname(fluor_name_map[[key]]))
                }
                return(p)
            }
        }
        ""
    }

    infer_marker_from_filename <- function(fn, fluor_guess = "") {
        stem <- tools::file_path_sans_ext(basename(fn))
        stem_lower <- tolower(stem)
        if (grepl("unstained|us_ut|^af($|[-_])", stem_lower)) return("Autofluorescence")
        marker_guess <- detect_alias(stem, marker_name_map, min_substring_n = 3)
        if (!nzchar(marker_guess)) return("")
        if (nzchar(fluor_guess) && normalize_token(marker_guess) == normalize_token(fluor_guess)) return("")
        marker_guess
    }

    resolve_detector_channel <- function(channel_name, detector_desc = "") {
        keys <- canonicalize_channel(c(channel_name, detector_desc))
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

    resolve_channel_alias <- function(channel_name) {
        keys <- canonicalize_channel(channel_name)
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

    get_expected_af_channel <- function(cytometer_name) {
        # Internal default AF-channel hints by cytometer family.
        # These hints are only used for optional AF-like auto-tagging heuristics.
        cy <- tolower(trimws(as.character(cytometer_name)))
        if (!nzchar(cy)) return("")
        if (grepl("aurora", cy, fixed = TRUE)) return("B1-A")
        if (grepl("northern lights", cy, fixed = TRUE)) return("B1-A")
        ""
    }

    next_af_tag <- function(existing_vals) {
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

    infer_fluor_from_detector <- function(channel_name, detector_desc = "") {
        resolved_channel <- resolve_detector_channel(channel_name, detector_desc)
        if (length(fluor_channel_map) > 0) {
            direct <- fluor_channel_map[resolved_channel]
            if (!is.na(direct) && nzchar(direct)) return(as.character(direct))
        }

        token <- toupper(gsub("\\s+", "", paste(channel_name, detector_desc)))
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

        for (k in names(exact_map)) {
            if (grepl(k, token, fixed = TRUE)) return(exact_map[[k]])
        }
        if (grepl("YG[0-9]", token) || grepl("YELLOWGREEN", token, fixed = TRUE)) return("PE")
        if (grepl("R[0-9]", token) || grepl("RED", token, fixed = TRUE)) return("APC")
        if (grepl("V[0-9]", token) || grepl("VIOLET", token, fixed = TRUE)) return("BV421")
        if (grepl("B[0-9]", token) || grepl("BLUE", token, fixed = TRUE)) return("FITC")
        if (grepl("UV[0-9]", token) || grepl("ULTRAVIOLET", token, fixed = TRUE)) return("BUV395")
        ""
    }

    # 2. Build entries for SCC files
    scc_rows <- list()
    for (fn in scc_files) {
        fluor <- infer_fluor_from_filename(fn)
        
        # Check for internal AF
        if (grepl("Unstained|US_UT", fn, ignore.case = TRUE)) fluor <- "AF_Internal"
        marker <- infer_marker_from_filename(fn, fluor)
        control_type <- infer_control_type_from_filename(fn, fluor_guess = fluor, marker_guess = marker)
        viability_flag <- infer_is_viability_from_filename(fn, fluor_guess = fluor, marker_guess = marker)
        
        scc_rows[[fn]] <- data.frame(
            filename = fn,
            fluorophore = fluor,
            marker = marker,
            channel = "",
            control.type = control_type,
            universal.negative = "",
            large.gate = "",
            is.viability = viability_flag,
            stringsAsFactors = FALSE
        )
    }
    
    # 3. Build entries for AF folder files
    af_rows <- list()
    for (i in seq_along(af_files)) {
        fn <- af_files[i]
        tag <- if(i == 1) "AF" else paste0("AF_", i)
        af_rows[[fn]] <- data.frame(
            filename = fn,
            fluorophore = tag,
            marker = "Autofluorescence",
            channel = "",
            control.type = "cells",
            universal.negative = "",
            large.gate = "TRUE",
            is.viability = "",
            stringsAsFactors = FALSE
        )
    }
    
    # 4. Combine
    # If we have files in 'af' folder, they take precedence over SCC-based detections
    df <- do.call(rbind, c(af_rows, scc_rows))
    
    # Strict uniqueness check: one row per filename
    # We keep the first occurrence (which would be the AF-folder one if duplicated)
    df <- df[!duplicated(df$filename), ]

    expected_af_channel <- resolve_channel_alias(get_expected_af_channel(cytometer))
    auto_af_files <- character(0)

    # Identify channels
    for (i in seq_len(nrow(df))) {
        fn <- df$filename[i]
        path <- if(include_af_folder && file.exists(file.path(af_folder, fn))) file.path(af_folder, fn) else file.path(input_folder, fn)
        pd <- NULL
        q_hi <- numeric()
        top_channels <- character()
        top_ratio <- NA_real_
        entropy_hi <- NA_real_
        tryCatch({
            ff <- flowCore::read.FCS(path, transformation = FALSE, truncate_max_range = FALSE)
            pd <- flowCore::pData(flowCore::parameters(ff))
            channel_alias_map <<- merge_alias_map(channel_alias_map, .build_channel_alias_map_from_pd(pd))
            fl_pd <- get_sorted_detectors(pd)
            if (length(fl_pd$names) > 0) {
                expr <- flowCore::exprs(ff)[, fl_pd$names, drop = FALSE]

                # Robust peak-channel detection:
                # gate to the central FSC/SSC cloud first, then score channels by 99.9th percentile.
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
                if (!is.na(peak_channel) && nzchar(peak_channel)) df$channel[i] <- peak_channel

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
        }, error = function(e) NULL)

        # Heuristic: unlabeled files that match cytometer AF profile become AF.
        # Conditions: no filename-based marker/fluor assignment, expected AF channel among top channels,
        # broad profile (low top/second ratio), and sufficiently diffuse spectrum.
        current_fluor <- trimws(as.character(df$fluorophore[i]))
        current_marker <- trimws(as.character(df$marker[i]))
        unresolved_label <- (is.na(current_fluor) || current_fluor == "" || identical(current_fluor, "Unknown")) &&
            (is.na(current_marker) || current_marker == "")
        af_like <- FALSE
        if (unresolved_label && nzchar(expected_af_channel) && length(top_channels) > 0) {
            resolved_top <- unique(vapply(top_channels, resolve_channel_alias, character(1)))
            af_in_top <- expected_af_channel %in% resolved_top
            broad_signal <- is.finite(top_ratio) && top_ratio <= 1.7
            diffuse_signal <- is.finite(entropy_hi) && entropy_hi >= 3.2
            af_like <- af_in_top && broad_signal && diffuse_signal
        }
        if (af_like) {
            new_tag <- next_af_tag(df$fluorophore)
            df$fluorophore[i] <- new_tag
            df$marker[i] <- "Autofluorescence"
            df$control.type[i] <- "cells"
            df$large.gate[i] <- "TRUE"
            auto_af_files <- c(auto_af_files, fn)
        }

        # Fill unresolved fluorophore after channel detection if requested.
        current_fluor <- trimws(as.character(df$fluorophore[i]))
        unresolved <- is.na(current_fluor) || current_fluor == "" || identical(current_fluor, "Unknown")
        if (unresolved) {
            if (unknown_fluor_policy == "filename") {
                df$fluorophore[i] <- tools::file_path_sans_ext(fn)
            } else if (unknown_fluor_policy == "by_channel") {
                detector_desc <- ""
                if (!is.null(pd) && "desc" %in% colnames(pd) && "name" %in% colnames(pd)) {
                    idx <- match(df$channel[i], as.character(pd$name))
                    if (!is.na(idx)) {
                        detector_desc <- as.character(pd$desc[idx])
                    }
                }
                guess <- infer_fluor_from_detector(df$channel[i], detector_desc)
                df$fluorophore[i] <- guess
            } else {
                df$fluorophore[i] <- ""
            }
        }

        # Fill marker from filename if still empty.
        current_marker <- trimws(as.character(df$marker[i]))
        if (is.na(current_marker) || current_marker == "") {
            df$marker[i] <- infer_marker_from_filename(fn, df$fluorophore[i])
        }

        # Fill control type from filename if still unresolved.
        current_type <- trimws(as.character(df$control.type[i]))
        if (is.na(current_type) || current_type == "") {
            df$control.type[i] <- infer_control_type_from_filename(
                fn,
                fluor_guess = df$fluorophore[i],
                marker_guess = df$marker[i]
            )
        }

        # Fill viability flag from filename/dye cues if still unresolved.
        current_viability <- toupper(trimws(as.character(df$is.viability[i])))
        if (is.na(current_viability) || current_viability == "") {
            df$is.viability[i] <- infer_is_viability_from_filename(
                fn,
                fluor_guess = df$fluorophore[i],
                marker_guess = df$marker[i]
            )
        }
    }

    # Final AF fixes (after channel + profile-based auto detection)
    primary_af_candidates <- scc_files[grep("Unstained|US_UT", scc_files, ignore.case = TRUE)]
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

    # Leave universal.negative empty by default for all rows.
    df$universal.negative <- ""

    # Normalize control-type values and force AF rows to cells.
    df$control.type <- tolower(trimws(as.character(df$control.type)))
    df$control.type[is.na(df$control.type)] <- ""
    invalid_type <- !(df$control.type %in% c("", "beads", "cells"))
    df$control.type[invalid_type] <- ""

    # Normalize viability values to TRUE/blank.
    df$is.viability <- toupper(trimws(as.character(df$is.viability)))
    df$is.viability[is.na(df$is.viability)] <- ""
    df$is.viability[df$is.viability %in% c("T", "TRUE", "1", "YES", "Y")] <- "TRUE"
    df$is.viability[!(df$is.viability %in% c("", "TRUE"))] <- ""

    # AF-tagged rows are always cells.
    is_af_row <- grepl("^AF($|_)", as.character(df$fluorophore), ignore.case = TRUE)
    df$control.type[is_af_row] <- "cells"

    # Normalize primary AF naming
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
        "large.gate",
        "is.viability"
    )
    keep <- intersect(preferred_order, colnames(df))
    df <- df[, c(keep, setdiff(colnames(df), keep)), drop = FALSE]

    # Use base write.csv for broad control-file parser compatibility.
    utils::write.csv(df, output_file, row.names = FALSE, quote = TRUE)
    return(df)
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
#' @param cytometer Cytometer type (for example `"Aurora"`).
#' @return Expanded reference matrix aligned to the detectors in `flow_frame`.
#' @export
#' @examples
#' \dontrun{
#' ff <- flowCore::read.FCS("samples/Sample1.fcs", transformation = FALSE)
#' M_ctrl <- get_control_spectra(
#'   flow_frame = ff,
#'   control_file = "fcs_mapping.csv",
#'   control_dir = "scc",
#'   cytometer = "Aurora"
#' )
#' }
get_control_spectra <- function(flow_frame,
                                control_file = "fcs_mapping.csv",
                                control_dir = "scc",
                                af_dir = "af",
                                method = "WLS",
                                cytometer = "Aurora") {
    cytometer_resolved <- cytometer

    # 1. Get detector info
    pd <- flowCore::pData(flowCore::parameters(flow_frame))
    det_info <- get_sorted_detectors(pd)
    detector_names <- det_info$names

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
