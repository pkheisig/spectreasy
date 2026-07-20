.reference_negative_key <- function(x) {
    tools::file_path_sans_ext(basename(trimws(as.character(x))))
}

.reference_is_bead_negative_file <- function(x) {
    stem <- tools::file_path_sans_ext(basename(trimws(as.character(x))))
    stem_lower <- tolower(stem)
    stem_norm <- gsub("[^[:alnum:]]+", "", stem_lower)
    tokens <- unlist(strsplit(stem_lower, "[^[:alnum:]]+"))
    tokens <- tokens[nzchar(tokens)]

    has_bead <- any(tokens %in% c("bead", "beads", "compbead", "compbeads")) ||
        grepl("compbeads?", stem_norm) ||
        grepl("beads?", stem_norm)
    has_negative <- any(tokens %in% c(
        "unstained", "unstain", "us", "ut", "usut", "usut1", "neg",
        "negative", "background", "bg", "blank", "minus"
    )) ||
        grepl("us[_ -]?ut", stem, ignore.case = TRUE) ||
        grepl("unstained|negative|background|blank", stem, ignore.case = TRUE) ||
        grepl("(^|[^[:alnum:]])(?:us|neg|bg)(?:[^[:alnum:]]|$)", stem, ignore.case = TRUE, perl = TRUE)

    has_bead && has_negative
}

.collect_reference_unstained_bead_negative <- function(control_df,
                                                       fcs_files,
                                                       detector_names,
                                                       config) {
    if (length(fcs_files) == 0) {
        return(NULL)
    }

    file_keys <- .reference_negative_key(fcs_files)
    names(fcs_files) <- file_keys
    bead_keys <- character()

    if (!is.null(control_df) && is.data.frame(control_df) && nrow(control_df) > 0 && "filename" %in% colnames(control_df)) {
        is_af <- .is_af_control_row(
            fluorophore = if ("fluorophore" %in% colnames(control_df)) control_df$fluorophore else NULL,
            marker = if ("marker" %in% colnames(control_df)) control_df$marker else NULL,
            filename = control_df$filename
        )
        control_type <- if ("control.type" %in% colnames(control_df)) {
            tolower(trimws(as.character(control_df$control.type)))
        } else {
            rep("", nrow(control_df))
        }
        file_is_bead <- grepl("beads?", basename(as.character(control_df$filename)), ignore.case = TRUE)
        file_is_bead_negative <- vapply(control_df$filename, .reference_is_bead_negative_file, logical(1))
        idx <- which((is_af & (control_type == "beads" | file_is_bead)) | file_is_bead_negative)
        if (length(idx) > 0) {
            bead_keys <- .reference_negative_key(control_df$filename[idx])
        }
    }

    if (length(bead_keys) == 0) {
        bead_keys <- file_keys[vapply(fcs_files, .reference_is_bead_negative_file, logical(1))]
    }
    bead_keys <- unique(bead_keys[nzchar(bead_keys) & bead_keys %in% names(fcs_files)])
    if (length(bead_keys) == 0) {
        return(NULL)
    }

    bead_events <- list()
    bead_gated_list <- list()
    for (key in bead_keys) {
        fcs_file <- unname(fcs_files[[key]])

        ff <- tryCatch(
            .spectreasy_read_fcs(fcs_file),
            error = function(e) NULL
        )
        if (is.null(ff)) {
            warning("Could not read unstained bead negative file: ", fcs_file)
            next
        }

        pd <- flowCore::pData(flowCore::parameters(ff))
        raw_data <- flowCore::exprs(ff)
        scatter_info <- .reference_background_scatter_gate(
            raw_data = raw_data,
            pd = pd,
            sample_type = "beads",
            filename = basename(fcs_file),
            config = config
        )

        neg_data <- if (!is.null(scatter_info)) scatter_info$gated_data else raw_data
        if (nrow(neg_data) > 0) {
            bead_events[[length(bead_events) + 1L]] <- neg_data[, detector_names, drop = FALSE]
            if (!is.null(scatter_info)) {
                bead_gated_list[[length(bead_gated_list) + 1L]] <- list(
                    events = neg_data[, detector_names, drop = FALSE],
                    scatter = neg_data[, c(scatter_info$fsc, scatter_info$ssc), drop = FALSE],
                    scatter_names = c(scatter_info$fsc, scatter_info$ssc)
                )
            }
            .spectreasy_console_field("Bead neg", basename(fcs_file))
        }
    }

    if (length(bead_events) == 0) {
        return(NULL)
    }

    bead_mat <- do.call(rbind, bead_events)
    bead_negative <- apply(bead_mat[, detector_names, drop = FALSE], 2, stats::median, na.rm = TRUE)
    bead_background <- .scc_background_from_gated_af_list(
        af_gated_list = bead_gated_list,
        detector_names = detector_names
    )
    if (!is.null(bead_background)) {
        attr(bead_negative, "scc_background") <- bead_background
    }
    bead_negative
}

# Identifies and loads the universal/marker-specific negative control profiles.
# If specified in the control mapping, reads the designated negative control files and computes
# their median spectral profiles to be subtracted from target samples.
# Returns a list of numeric vectors mapping fluorophore names to their negative control profiles.
.collect_reference_universal_negatives <- function(control_df,
                                                   fcs_files,
                                                   detector_names,
                                                   sample_patterns,
                                                   config) {
    out <- list()
    if (is.null(control_df) || !is.data.frame(control_df) || !("universal.negative" %in% colnames(control_df))) {
        return(out)
    }

    uv_vals <- trimws(as.character(control_df$universal.negative))
    uv_vals[is.na(uv_vals)] <- ""
    uv_vals <- unique(uv_vals[nzchar(uv_vals) & !toupper(uv_vals) %in% c("FALSE", "TRUE", "AF")])
    if (length(uv_vals) == 0) {
        return(out)
    }

    file_keys <- .reference_negative_key(fcs_files)
    names(fcs_files) <- file_keys

    for (uv in uv_vals) {
        key <- .reference_negative_key(uv)
        if (!nzchar(key) || !key %in% names(fcs_files)) {
            next
        }

        fcs_file <- unname(fcs_files[[key]])
        sn_ext <- basename(fcs_file)
        sn <- tools::file_path_sans_ext(sn_ext)
        row_info <- .get_control_rows_for_reference(control_df, c(sn_ext, sn))
        sample_info <- .resolve_reference_sample_type(
            filename = sn,
            row_info = row_info,
            patterns = sample_patterns,
            default = config$default_sample_type
        )

        ff <- tryCatch(
            .spectreasy_read_fcs(fcs_file),
            error = function(e) NULL
        )
        if (is.null(ff)) {
            warning("Could not read universal.negative file: ", fcs_file)
            next
        }

        pd <- flowCore::pData(flowCore::parameters(ff))
        raw_data <- flowCore::exprs(ff)
        scatter_info <- .reference_background_scatter_gate(
            raw_data = raw_data,
            pd = pd,
            sample_type = sample_info$type,
            filename = basename(fcs_file),
            config = config
        )

        neg_data <- if (!is.null(scatter_info)) scatter_info$gated_data else raw_data
        negative <- apply(neg_data[, detector_names, drop = FALSE], 2, stats::median, na.rm = TRUE)
        if (!is.null(scatter_info)) {
            background <- .scc_background_from_gated_af_list(
                af_gated_list = list(list(
                    events = neg_data[, detector_names, drop = FALSE],
                    scatter = neg_data[, c(scatter_info$fsc, scatter_info$ssc), drop = FALSE],
                    scatter_names = c(scatter_info$fsc, scatter_info$ssc)
                )),
                detector_names = detector_names
            )
            if (!is.null(background)) attr(negative, "scc_background") <- background
        }
        out[[key]] <- negative
        .spectreasy_console_field("Negative", basename(fcs_file))
    }

    out
}

.resolve_reference_scc_negative_source <- function(row_info,
                                                   sample_type,
                                                   af_data_raw = NULL,
                                                   scc_background = NULL,
                                                   universal_negatives = NULL,
                                                   bead_negative = NULL) {
    if (identical(sample_type, "beads")) {
        return(list(
            negative = bead_negative,
            background = attr(bead_negative, "scc_background", exact = TRUE),
            source = "AF_beads"
        ))
    }

    uv_val <- if (nrow(row_info) > 0L && "universal.negative" %in% colnames(row_info)) {
        trimws(as.character(row_info$universal.negative[1]))
    } else {
        ""
    }
    uv_upper <- toupper(uv_val)
    uv_key <- .reference_negative_key(uv_val)
    named_negative <- if (
        nzchar(uv_key) &&
        !uv_upper %in% c("FALSE", "TRUE", "AF") &&
        !is.null(universal_negatives) &&
        uv_key %in% names(universal_negatives)
    ) {
        universal_negatives[[uv_key]]
    } else {
        NULL
    }

    if (!is.null(named_negative)) {
        return(list(
            negative = named_negative,
            background = attr(named_negative, "scc_background", exact = TRUE),
            source = uv_val
        ))
    }

    list(
        negative = af_data_raw,
        background = scc_background,
        source = "AF"
    )
}

# Helper to extract the description attribute of a channel for plotting.
# Falls back to the channel name if no description column or value is found.
# Returns the description string.
