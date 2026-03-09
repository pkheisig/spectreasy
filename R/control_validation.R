#' Validate SCC Control Mapping Before Unmixing
#'
#' Performs strict preflight checks on the control mapping used for SCC workflows.
#' This catches common setup issues early (missing files, wrong channel names,
#' malformed universal negatives, incomplete mappings).
#'
#' @param control_df Control mapping data.frame.
#' @param scc_dir SCC directory containing FCS files.
#' @param include_multi_af Logical. Whether AF files from `af_dir` are expected.
#' @param af_dir Optional AF directory with additional FCS files.
#' @param require_all_scc_mapped Logical. If TRUE, every effective SCC file must
#'   have a row in `control_df`.
#' @param require_channels Logical. If TRUE, non-AF rows must define `channel`.
#' @param stop_on_error Logical. If TRUE, stop with a combined error message.
#' @return A list with `ok`, `errors`, and `warnings`.
#' @examples
#' \dontrun{
#' control_df <- read.csv("fcs_mapping.csv", stringsAsFactors = FALSE)
#' preflight <- validate_control_file_mapping(control_df, scc_dir = "scc")
#' preflight$ok
#' }
#' @export
validate_control_file_mapping <- function(
    control_df,
    scc_dir = "scc",
    include_multi_af = FALSE,
    af_dir = "af",
    require_all_scc_mapped = TRUE,
    require_channels = TRUE,
    stop_on_error = FALSE
) {
    errors <- character()
    warnings <- character()
    required_cols <- c("filename", "fluorophore", "channel")

    as_chr <- function(x) {
        if (is.null(x)) return(character())
        out <- trimws(as.character(x))
        out[is.na(out)] <- ""
        out
    }

    normalize_channel <- function(x) {
        out <- toupper(gsub("\\s+", "", trimws(as.character(x))))
        out[is.na(out)] <- ""
        out
    }

    channel_in_detectors <- function(channel_value, detector_names, alias_map = character()) {
        key <- normalize_channel(channel_value)
        if (!nzchar(key)) return(TRUE)
        det_norm <- normalize_channel(detector_names)
        alias_targets <- character()
        if (length(alias_map) > 0) {
            hits <- alias_map[c(key, gsub("-A$", "", key), paste0(key, "-A"))]
            hits <- hits[!is.na(hits) & nzchar(hits)]
            alias_targets <- normalize_channel(unname(hits))
        }
        candidates <- unique(c(key, gsub("-A$", "", key), paste0(key, "-A"), alias_targets))
        any(candidates %in% det_norm)
    }

    effective_scc_files <- function(dir_path) {
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
        unique(files)
    }

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

    if (length(errors) == 0) {
        df <- control_df
        df$filename <- as_chr(df$filename)
        df$fluorophore <- as_chr(df$fluorophore)
        df$channel <- as_chr(df$channel)
        if (!("universal.negative" %in% colnames(df))) {
            df$universal.negative <- ""
        }
        df$universal.negative <- as_chr(df$universal.negative)

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

        is_af <- grepl("^AF($|_)", df$fluorophore, ignore.case = TRUE)
        if (require_channels) {
            missing_channel_rows <- which(!is_af & (df$channel == "" | is.na(df$channel)))
            if (length(missing_channel_rows) > 0) {
                bad_files <- df$filename[missing_channel_rows]
                errors <- c(errors, paste0("Missing channel for non-AF rows: ", paste(bad_files, collapse = ", ")))
            }
        }

        scc_files <- effective_scc_files(scc_dir)
        af_files <- if (include_multi_af && dir.exists(af_dir)) {
            list.files(af_dir, pattern = "\\.fcs$", full.names = FALSE, ignore.case = TRUE)
        } else {
            character()
        }
        known_files <- unique(c(scc_files, af_files))

        unknown_in_control <- setdiff(df$filename, known_files)
        if (length(unknown_in_control) > 0) {
            warnings <- c(warnings, paste0("Ignoring control rows for files not found in selected SCC/AF dirs: ", paste(unknown_in_control, collapse = ", ")))
        }

        if (require_all_scc_mapped) {
            unmapped_scc <- setdiff(scc_files, df$filename)
            if (length(unmapped_scc) > 0) {
                errors <- c(errors, paste0("SCC files missing from control file: ", paste(unmapped_scc, collapse = ", ")))
            }
        }

        # Check universal.negative values for files that are actually processed.
        active_rows <- df$filename %in% c(scc_files, af_files)
        uv_vals <- df$universal.negative[active_rows]
        keyword_vals <- c("", "TRUE", "FALSE", "AF")
        unresolved_uv <- uv_vals[!(uv_vals %in% keyword_vals) & !(uv_vals %in% known_files)]
        if (length(unresolved_uv) > 0) {
            errors <- c(errors, paste0("universal.negative points to unknown files for active controls: ", paste(unique(unresolved_uv), collapse = ", ")))
        }

        # Per-file channel existence check
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
            det_names <- as_chr(pd$name)
            alias_map <- .build_channel_alias_map_from_pd(pd)
            ch_vals <- unique(df$channel[df$filename == fn])
            missing_channels <- ch_vals[!vapply(
                ch_vals,
                channel_in_detectors,
                logical(1),
                detector_names = det_names,
                alias_map = alias_map
            )]
            if (length(missing_channels) > 0) {
                errors <- c(errors, paste0("Invalid channel in ", fn, ": ", paste(missing_channels, collapse = ", ")))
            }
        }
    }

    ok <- length(errors) == 0
    out <- list(ok = ok, errors = unique(errors), warnings = unique(warnings))

    if (stop_on_error && !ok) {
        stop(
            paste(
                c(
                    "Control mapping preflight failed:",
                    paste0(" - ", out$errors),
                    if (length(out$warnings) > 0) c("Warnings:", paste0(" - ", out$warnings)) else NULL
                ),
                collapse = "\n"
            ),
            call. = FALSE
        )
    }

    out
}
