# Internal helpers for control-file preflight.
.control_validation_as_chr <- function(x) {
    if (is.null(x)) return(character())
    out <- trimws(as.character(x))
    out[is.na(out)] <- ""
    out
}

.control_validation_normalize_channel <- function(x) {
    out <- toupper(gsub("\\s+", "", trimws(as.character(x))))
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

.control_validation_effective_scc_files <- function(dir_path, exclude_af = FALSE) {
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
    if (isTRUE(exclude_af)) {
        files <- files[!.is_af_filename(files)]
    }
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

    if (any(df$fluorophore == "" | is.na(df$fluorophore))) {
        bad <- which(df$fluorophore == "" | is.na(df$fluorophore))
        errors <- c(errors, paste0("Empty fluorophore in control file rows: ", paste(bad, collapse = ", ")))
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

.collect_control_validation_known_files <- function(scc_dir, af_dir = "af", include_multi_af = FALSE, exclude_af = FALSE) {
    scc_files <- .control_validation_effective_scc_files(scc_dir, exclude_af = exclude_af)
    af_files <- if (!isTRUE(exclude_af) && isTRUE(include_multi_af) && dir.exists(af_dir)) {
        list.files(af_dir, pattern = "\\.fcs$", full.names = FALSE, ignore.case = TRUE)
    } else {
        character()
    }

    list(scc_files = scc_files, af_files = af_files, known_files = unique(c(scc_files, af_files)))
}

.validate_control_file_mapping_coverage <- function(df, is_af, known_files, scc_files, af_files, require_all_scc_mapped = TRUE, exclude_af = FALSE) {
    errors <- character()
    warnings <- character()

    eligible_rows <- if (isTRUE(exclude_af)) !is_af else rep(TRUE, nrow(df))

    unknown_in_control <- setdiff(df$filename[eligible_rows], known_files)
    if (length(unknown_in_control) > 0) {
        warnings <- c(warnings, paste0("Ignoring control rows for files not found in selected SCC/AF dirs: ", paste(unknown_in_control, collapse = ", ")))
    }

    if (isTRUE(require_all_scc_mapped)) {
        unmapped_scc <- setdiff(scc_files, df$filename[eligible_rows])
        if (length(unmapped_scc) > 0) {
            errors <- c(errors, paste0("SCC files missing from control file: ", paste(unmapped_scc, collapse = ", ")))
        }
    }

    active_rows <- df$filename %in% c(scc_files, af_files) & eligible_rows
    list(errors = errors, warnings = warnings, active_rows = active_rows)
}

.validate_control_file_universal_negative <- function(df, active_rows, known_files) {
    errors <- character()
    uv_vals <- df$universal.negative[active_rows]
    keyword_vals <- c("", "TRUE", "FALSE", "AF")
    unresolved_uv <- uv_vals[!(uv_vals %in% keyword_vals) & !(uv_vals %in% known_files)]
    if (length(unresolved_uv) > 0) {
        errors <- c(errors, paste0("universal.negative points to unknown files for active controls: ", paste(unique(unresolved_uv), collapse = ", ")))
    }
    errors
}

.validate_control_file_channels_exist <- function(df, active_rows, scc_files, af_files, scc_dir, af_dir) {
    errors <- character()
    files_to_check <- unique(df$filename[active_rows & df$channel != "" & !is.na(df$channel)])

    for (fn in files_to_check) {
        path <- if (fn %in% scc_files) file.path(scc_dir, fn) else file.path(af_dir, fn)
        if (!file.exists(path)) next
        ff <- tryCatch(
            flowCore::read.FCS(path, transformation = FALSE, truncate_max_range = FALSE),
            error = function(e) NULL
        )
        if (is.null(ff)) {
            errors <- c(errors, paste0("Could not read FCS file for validation: ", path))
            next
        }
        pd <- flowCore::pData(flowCore::parameters(ff))
        det_names <- .control_validation_as_chr(pd$name)
        alias_map <- .build_channel_alias_map_from_pd(pd)
        ch_vals <- unique(df$channel[df$filename == fn])
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
#' @param include_multi_af Logical. Whether AF files from `af_dir` are expected.
#' @param exclude_af Logical. If `TRUE`, AF/unstained controls are optional and
#'   ignored during mapping completeness checks.
#' @param af_dir Optional AF directory with additional FCS files.
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
    include_multi_af = FALSE,
    exclude_af = FALSE,
    af_dir = "af",
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

        known <- .collect_control_validation_known_files(
            scc_dir = scc_dir,
            af_dir = af_dir,
            include_multi_af = include_multi_af,
            exclude_af = exclude_af
        )
        coverage <- .validate_control_file_mapping_coverage(
            df = df,
            is_af = row_checks$is_af,
            known_files = known$known_files,
            scc_files = known$scc_files,
            af_files = known$af_files,
            require_all_scc_mapped = require_all_scc_mapped,
            exclude_af = exclude_af
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
                af_files = known$af_files,
                scc_dir = scc_dir,
                af_dir = af_dir
            )
        )
    }

    out <- list(ok = length(errors) == 0, errors = unique(errors), warnings = unique(warnings))
    if (isTRUE(stop_on_error) && !out$ok) {
        .stop_control_validation_errors(out)
    }

    out
}
