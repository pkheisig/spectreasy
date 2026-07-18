# Internal helpers for control-file preflight.
.control_validation_as_chr <- function(x) {
    if (is.null(x)) return(character())
    out <- trimws(as.character(x))
    out[is.na(out)] <- ""
    out
}

.control_validation_file_key <- function(x) {
    tolower(tools::file_path_sans_ext(basename(.control_validation_as_chr(x))))
}

.control_validation_normalize_channel <- function(x) {
    out <- toupper(gsub("\\s+", "", trimws(as.character(x))))
    out <- gsub("([A-Z]+)-([0-9])", "\\1\\2", out, perl = TRUE)
    out[is.na(out)] <- ""
    out
}

.control_validation_channel_in_detectors <- function(channel_value, detector_names, alias_map = character()) {
    key <- .control_validation_normalize_channel(channel_value)
    if (!nzchar(key)) return(TRUE)
    det_norm <- .control_validation_normalize_channel(detector_names)
    alias_targets <- character()
    if (length(alias_map) > 0) {
        hits <- alias_map[c(key, gsub("-A$", "", key), paste0(key, "-A"))]
        hits <- hits[!is.na(hits) & nzchar(hits)]
        alias_targets <- .control_validation_normalize_channel(unname(hits))
    }
    candidates <- unique(c(key, gsub("-A$", "", key), paste0(key, "-A"), alias_targets))
    any(candidates %in% det_norm)
}

.control_validation_effective_scc_files <- function(dir_path) {
    files <- list.files(dir_path, pattern = "\\.fcs$", full.names = FALSE, ignore.case = TRUE)
    if (length(files) == 0) return(character())
    bases <- tools::file_path_sans_ext(files)
    aligned_idx <- grep("_aligned$", bases)
    if (length(aligned_idx) > 0) {
        keep <- rep(TRUE, length(files))
        for (idx in aligned_idx) {
            base_original <- sub("_aligned$", "", bases[idx])
            dup_idx <- which(bases == base_original)
            if (length(dup_idx) > 0) keep[dup_idx] <- FALSE
            keep[idx] <- TRUE
        }
        files <- files[keep]
    }
    files <- unique(files)
    files
}

.validate_control_file_structure <- function(control_df, scc_dir, required_cols = c("filename", "fluorophore", "channel")) {
    errors <- character()

    if (!is.data.frame(control_df)) {
        errors <- c(errors, "control_df must be a data.frame.")
    } else {
        missing_cols <- setdiff(required_cols, colnames(control_df))
        if (length(missing_cols) > 0) {
            errors <- c(errors, paste0("Missing required columns: ", paste(missing_cols, collapse = ", ")))
        }
    }

    if (!dir.exists(scc_dir)) {
        errors <- c(errors, paste0("scc_dir does not exist: ", scc_dir))
    }

    errors
}

.normalize_control_validation_df <- function(control_df) {
    df <- control_df
    df$filename <- .control_validation_as_chr(df$filename)
    df$fluorophore <- .control_validation_as_chr(df$fluorophore)
    df$channel <- .control_validation_as_chr(df$channel)
    if (!("control.type" %in% colnames(df))) {
        df$control.type <- ""
    }
    df$control.type <- tolower(.control_validation_as_chr(df$control.type))
    if (!("universal.negative" %in% colnames(df))) {
        df$universal.negative <- ""
    }
    df$universal.negative <- .control_validation_as_chr(df$universal.negative)
    df
}

.validate_control_validation_rows <- function(df, require_channels = TRUE) {
    errors <- character()

    if (any(df$filename == "" | is.na(df$filename))) {
        bad <- which(df$filename == "" | is.na(df$filename))
        errors <- c(errors, paste0("Empty filename in control file rows: ", paste(bad, collapse = ", ")))
    }

    dup_filenames <- unique(df$filename[duplicated(df$filename)])
    if (length(dup_filenames) > 0) {
        errors <- c(errors, paste0("Duplicate filename entries in control file: ", paste(dup_filenames, collapse = ", ")))
    }
    filename_keys <- .control_validation_file_key(df$filename)
    dup_filename_keys <- unique(filename_keys[nzchar(filename_keys) & duplicated(filename_keys)])
    if (length(dup_filename_keys) > 0) {
        duplicate_rows <- df$filename[filename_keys %in% dup_filename_keys]
        errors <- c(
            errors,
            paste0(
                "Duplicate filename stems in control file (case and .fcs extension are ignored): ",
                paste(unique(duplicate_rows), collapse = ", ")
            )
        )
    }

    if (any(df$fluorophore == "" | is.na(df$fluorophore))) {
        bad <- which(df$fluorophore == "" | is.na(df$fluorophore))
        errors <- c(errors, paste0("Empty fluorophore in control file rows: ", paste(bad, collapse = ", ")))
    }

    invalid_control_type <- which(nzchar(df$control.type) & !df$control.type %in% c("beads", "cells"))
    if (length(invalid_control_type) > 0) {
        errors <- c(errors, paste0("Invalid control.type in rows: ", paste(invalid_control_type, collapse = ", "), ". Use 'beads', 'cells', or leave empty."))
    }

    is_af <- .is_af_control_row(
        fluorophore = df$fluorophore,
        marker = if ("marker" %in% colnames(df)) df$marker else NULL,
        filename = df$filename
    )
    if (isTRUE(require_channels)) {
        missing_channel_rows <- which(!is_af & (df$channel == "" | is.na(df$channel)))
        if (length(missing_channel_rows) > 0) {
            bad_files <- df$filename[missing_channel_rows]
            errors <- c(errors, paste0("Missing channel for non-AF rows: ", paste(bad_files, collapse = ", ")))
        }
    }

    list(errors = errors, is_af = is_af)
}

.collect_control_validation_known_files <- function(scc_dir) {
    scc_files <- .control_validation_effective_scc_files(scc_dir)
    list(scc_files = scc_files, known_files = unique(scc_files))
}

.control_validation_select_scc_files <- function(control_df, fcs_files) {
    if (is.null(control_df) || !is.data.frame(control_df) ||
        nrow(control_df) == 0 || !("filename" %in% colnames(control_df))) {
        return(fcs_files)
    }

    mapped <- .control_validation_as_chr(control_df$filename)
    if ("universal.negative" %in% colnames(control_df)) {
        universal_negatives <- .control_validation_as_chr(control_df$universal.negative)
        universal_negatives <- universal_negatives[
            nzchar(universal_negatives) &
                !(toupper(universal_negatives) %in% c("TRUE", "FALSE", "AF"))
        ]
        mapped <- c(mapped, universal_negatives)
    }
    mapped <- unique(mapped[nzchar(mapped)])

    file_keys <- .control_validation_file_key(fcs_files)
    mapped_keys <- unique(.control_validation_file_key(mapped))
    ambiguous <- intersect(unique(file_keys[duplicated(file_keys)]), mapped_keys)
    if (length(ambiguous) > 0L) {
        stop(
            "Multiple FCS files match the same mapped filename stem: ",
            paste(basename(fcs_files[file_keys %in% ambiguous]), collapse = ", "),
            call. = FALSE
        )
    }
    fcs_files[file_keys %in% mapped_keys]
}

.validate_control_file_mapping_coverage <- function(df, known_files, scc_files, require_all_scc_mapped = TRUE) {
    errors <- character()
    warnings <- character()

    control_keys <- .control_validation_file_key(df$filename)
    known_keys <- .control_validation_file_key(known_files)
    scc_keys <- .control_validation_file_key(scc_files)
    unknown_in_control <- df$filename[!(control_keys %in% known_keys)]
    if (length(unknown_in_control) > 0) {
        errors <- c(errors, paste0("Control files missing from scc_dir: ", paste(unknown_in_control, collapse = ", ")))
    }

    unmapped_scc <- scc_files[!(scc_keys %in% control_keys)]
    if (length(unmapped_scc) > 0) {
        if (isTRUE(require_all_scc_mapped)) {
            errors <- c(errors, paste0("SCC files missing from control file: ", paste(unmapped_scc, collapse = ", ")))
        } else {
            warnings <- c(warnings, paste0("Ignoring SCC files not listed in control file: ", paste(unmapped_scc, collapse = ", ")))
        }
    }

    active_rows <- control_keys %in% scc_keys
    list(errors = errors, warnings = warnings, active_rows = active_rows)
}

.validate_control_file_universal_negative <- function(df, active_rows, known_files) {
    errors <- character()
    uv_vals <- df$universal.negative[active_rows]
    keyword_vals <- c("", "TRUE", "FALSE", "AF")
    known_keys <- unique(.control_validation_file_key(known_files))
    uv_keys <- .control_validation_file_key(uv_vals)
    uv_upper <- toupper(uv_vals)
    unresolved_uv <- uv_vals[!(uv_upper %in% keyword_vals) & !(uv_keys %in% known_keys)]
    if (length(unresolved_uv) > 0) {
        errors <- c(errors, paste0("universal.negative points to unknown files for active controls: ", paste(unique(unresolved_uv), collapse = ", ")))
    }
    errors
}

.validate_control_file_channels_exist <- function(df, active_rows, scc_files, scc_dir) {
    errors <- character()
    files_to_check <- unique(df$filename[active_rows])

    for (fn in files_to_check) {
        file_idx <- match(.control_validation_file_key(fn), .control_validation_file_key(scc_files))
        if (is.na(file_idx)) next
        actual_name <- scc_files[[file_idx]]
        path <- file.path(scc_dir, actual_name)
        ff <- tryCatch(
            suppressWarnings(.spectreasy_read_fcs(path)),
            error = function(e) NULL
        )
        if (is.null(ff)) {
            errors <- c(errors, paste0("Could not read FCS file for validation: ", path))
            next
        }
        pd <- flowCore::pData(flowCore::parameters(ff))
        det_names <- .control_validation_as_chr(pd$name)
        alias_map <- .build_channel_alias_map_from_pd(pd)
        ch_vals <- unique(df$channel[
            active_rows & .control_validation_file_key(df$filename) == .control_validation_file_key(fn)
        ])
        ch_vals <- ch_vals[nzchar(ch_vals) & !is.na(ch_vals)]
        if (length(ch_vals) == 0) next
        missing_channels <- ch_vals[!vapply(
            ch_vals,
            .control_validation_channel_in_detectors,
            logical(1),
            detector_names = det_names,
            alias_map = alias_map
        )]
        if (length(missing_channels) > 0) {
            errors <- c(errors, paste0("Invalid channel in ", fn, ": ", paste(missing_channels, collapse = ", ")))
        }
    }

    errors
}

.stop_control_validation_errors <- function(out) {
    stop(
        paste(
            c(
                "Control mapping preflight failed:",
                paste0(" - ", out$errors),
                if (length(out$warnings) > 0) c("Preflight notes:", paste0(" - ", out$warnings)) else NULL
            ),
            collapse = "\n"
        )
    )
}

#' Validate SCC Control Mapping Before Unmixing
#'
#' Performs strict preflight checks on the control mapping used for SCC workflows.
#' This catches common setup issues early (missing files, wrong channel names,
#' malformed universal negatives, incomplete mappings).
#'
#' @param control_df Control mapping data.frame.
#' @param scc_dir SCC directory containing FCS files.
#' @param require_all_scc_mapped Logical. If TRUE, every effective SCC file must
#'   have a row in `control_df`.
#' @param require_channels Logical. If TRUE, non-AF rows must define `channel`.
#' @param stop_on_error Logical. If TRUE, stop with a combined error message.
#' @return A list with `ok`, `errors`, and `warnings`.
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
#' control_df <- data.frame(
#'   filename = c("FITC_cells.fcs", "PE_cells.fcs"),
#'   fluorophore = c("FITC", "PE"),
#'   channel = c("B2-A", "YG1-A"),
#'   stringsAsFactors = FALSE
#' )
#'
#' preflight <- validate_control_file_mapping(control_df, scc_dir = scc_dir)
#' preflight$ok
#' @export
validate_control_file_mapping <- function(
    control_df,
    scc_dir = "scc",
    require_all_scc_mapped = TRUE,
    require_channels = TRUE,
    stop_on_error = FALSE
) {
    errors <- .validate_control_file_structure(control_df, scc_dir = scc_dir)
    warnings <- character()

    if (length(errors) == 0) {
        df <- .normalize_control_validation_df(control_df)
        row_checks <- .validate_control_validation_rows(df, require_channels = require_channels)
        errors <- c(errors, row_checks$errors)

        known <- .collect_control_validation_known_files(scc_dir = scc_dir)
        coverage <- .validate_control_file_mapping_coverage(
            df = df,
            known_files = known$known_files,
            scc_files = known$scc_files,
            require_all_scc_mapped = require_all_scc_mapped
        )
        errors <- c(errors, coverage$errors)
        warnings <- c(warnings, coverage$warnings)
        errors <- c(errors, .validate_control_file_universal_negative(df, active_rows = coverage$active_rows, known_files = known$known_files))
        errors <- c(
            errors,
            .validate_control_file_channels_exist(
                df,
                active_rows = coverage$active_rows,
                scc_files = known$scc_files,
                scc_dir = scc_dir
            )
        )
    }

    out <- list(ok = length(errors) == 0, errors = unique(errors), warnings = unique(warnings))
    if (isTRUE(stop_on_error) && !out$ok) {
        .stop_control_validation_errors(out)
    }

    out
}
