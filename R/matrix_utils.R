.as_reference_matrix <- function(M, arg_name = "M") {
    if (is.null(M)) {
        stop(arg_name, " must not be NULL.")
    }
    
    if (is.matrix(M)) {
        if (!is.numeric(M)) {
            suppressWarnings(storage.mode(M) <- "numeric")
        }
        if (!is.numeric(M)) {
            stop(arg_name, " must be numeric.")
        }
        return(M)
    }
    
    if (inherits(M, "data.frame")) {
        df <- as.data.frame(M, stringsAsFactors = FALSE)
        if (ncol(df) == 0) {
            stop(arg_name, " has no columns.")
        }
        
        lower_names <- tolower(colnames(df))
        marker_idx <- match("marker", lower_names)
        if (is.na(marker_idx)) marker_idx <- match("fluorophore", lower_names)
        if (is.na(marker_idx)) marker_idx <- match("file", lower_names)
        if (is.na(marker_idx) && ncol(df) > 0 && !is.numeric(df[[1]])) {
            marker_idx <- 1L
        }
        
        if (!is.na(marker_idx)) {
            marker_names <- as.character(df[[marker_idx]])
            df <- df[, -marker_idx, drop = FALSE]
        } else {
            marker_names <- rownames(df)
            if (is.null(marker_names) || !any(nzchar(marker_names))) {
                marker_names <- as.character(seq_len(nrow(df)))
            }
        }
        
        if (ncol(df) == 0) {
            stop(arg_name, " has no detector columns after removing marker labels.")
        }
        
        numeric_cols <- lapply(df, function(x) suppressWarnings(as.numeric(x)))
        bad_cols <- vapply(seq_along(numeric_cols), function(i) {
            original <- df[[i]]
            converted <- numeric_cols[[i]]
            any(!is.na(original)) && all(is.na(converted))
        }, logical(1))
        
        if (any(bad_cols)) {
            stop(
                arg_name,
                " contains non-numeric detector columns: ",
                paste(names(df)[bad_cols], collapse = ", ")
            )
        }
        
        mat <- as.matrix(as.data.frame(numeric_cols, stringsAsFactors = FALSE))
        rownames(mat) <- marker_names
        colnames(mat) <- colnames(df)
        return(mat)
    }
    
    stop(arg_name, " must be a numeric matrix or a data.frame with detector columns.")
}

.normalize_channel_token <- function(x) {
    out <- toupper(gsub("\\s+", "", trimws(as.character(x))))
    out[is.na(out)] <- ""
    out
}

.build_channel_alias_map_from_pd <- function(pd) {
    if (is.null(pd) || !is.data.frame(pd) || !all(c("name", "desc") %in% colnames(pd))) {
        return(character())
    }

    det_info <- tryCatch(get_sorted_detectors(pd), error = function(e) NULL)
    if (is.null(det_info) || length(det_info$names) == 0) {
        return(character())
    }
    detector_names <- as.character(det_info$names)
    if (length(detector_names) == 0) {
        return(character())
    }

    alias_map <- character()
    add_alias <- function(alias, target) {
        key <- .normalize_channel_token(alias)
        if (!nzchar(key) || is.na(target) || !nzchar(target)) return(invisible(NULL))
        if (!key %in% names(alias_map)) alias_map[[key]] <<- target
        invisible(NULL)
    }

    # Direct aliases from detector names.
    for (det in detector_names) {
        det_key <- .normalize_channel_token(det)
        add_alias(det_key, det)
        add_alias(gsub("-A$", "", det_key), det)
        add_alias(paste0(gsub("-A$", "", det_key), "-A"), det)
    }

    # Build code-style aliases (UV1, V7, B2, YG1, R4) from detector metadata.
    idx <- match(detector_names, as.character(pd$name))
    desc <- trimws(as.character(pd$desc[idx]))
    desc[is.na(desc)] <- ""

    parse_laser_nm <- function(x) {
        if (!nzchar(x)) return(NA_real_)
        m <- regexec("([0-9]{3})\\s*NM", toupper(x))
        g <- regmatches(toupper(x), m)[[1]]
        if (length(g) < 2) return(NA_real_)
        as.numeric(g[2])
    }

    parse_center_nm <- function(x) {
        if (!nzchar(x)) return(NA_real_)
        m <- regexec("-\\s*([0-9]{3})\\s*/", toupper(x))
        g <- regmatches(toupper(x), m)[[1]]
        if (length(g) < 2) return(NA_real_)
        as.numeric(g[2])
    }

    laser_nm <- vapply(desc, parse_laser_nm, numeric(1))
    center_nm <- vapply(desc, parse_center_nm, numeric(1))

    laser_prefix <- vapply(laser_nm, function(v) {
        if (!is.finite(v)) return(NA_character_)
        if (v < 380) return("UV")
        if (v < 460) return("V")
        if (v < 530) return("B")
        if (v < 600) return("YG")
        "R"
    }, character(1))

    code_df <- data.frame(
        detector = detector_names,
        prefix = laser_prefix,
        center = center_nm,
        stringsAsFactors = FALSE
    )
    code_df <- code_df[!is.na(code_df$prefix) & nzchar(code_df$prefix), , drop = FALSE]
    if (nrow(code_df) > 0) {
        split_idx <- split(seq_len(nrow(code_df)), code_df$prefix)
        for (prefix in names(split_idx)) {
            rows <- split_idx[[prefix]]
            block <- code_df[rows, , drop = FALSE]
            ord <- order(
                ifelse(is.finite(block$center), block$center, Inf),
                block$detector
            )
            block <- block[ord, , drop = FALSE]
            for (i in seq_len(nrow(block))) {
                alias_base <- paste0(prefix, i)
                add_alias(alias_base, block$detector[i])
                add_alias(paste0(alias_base, "-A"), block$detector[i])
            }
        }
    }

    alias_map
}

.is_passthrough_parameter <- function(param_names) {
    if (length(param_names) == 0) {
        return(logical(0))
    }

    param_names <- trimws(as.character(param_names))
    is_scatter <- grepl("^(FSC|SSC)", param_names, ignore.case = TRUE)
    is_time <- grepl("^TIME($|[^A-Z0-9])", param_names, ignore.case = TRUE)

    is_scatter | is_time
}

.get_passthrough_parameter_names <- function(param_names, detector_names = character()) {
    param_names <- as.character(param_names)
    detector_names <- as.character(detector_names)

    keep <- param_names[.is_passthrough_parameter(param_names)]
    keep[!(keep %in% detector_names)]
}

.append_passthrough_parameters <- function(out, full_data, detector_names = character()) {
    keep <- .get_passthrough_parameter_names(colnames(full_data), detector_names = detector_names)
    if (length(keep) == 0) {
        return(out)
    }

    out[keep] <- as.data.frame(full_data[, keep, drop = FALSE], check.names = FALSE)
    out
}

.get_result_metadata_columns <- function(col_names) {
    unique(c(
        "File",
        .get_passthrough_parameter_names(col_names)
    ))
}

.get_primary_scatter_channels <- function(col_names) {
    pick_primary <- function(prefix) {
        exact <- col_names[toupper(col_names) == paste0(prefix, "-A")]
        if (length(exact) > 0) {
            return(exact[[1]])
        }

        area <- col_names[grepl(paste0("^", prefix, ".*-A$"), col_names, ignore.case = TRUE)]
        if (length(area) > 0) {
            return(area[[1]])
        }

        any_match <- col_names[grepl(paste0("^", prefix), col_names, ignore.case = TRUE)]
        if (length(any_match) > 0) {
            return(any_match[[1]])
        }

        NA_character_
    }

    list(
        fsc = pick_primary("FSC"),
        ssc = pick_primary("SSC")
    )
}

.resolve_control_file_path <- function(control_file = "fcs_mapping.csv") {
    path <- as.character(control_file)[1]
    if (is.na(path) || !nzchar(trimws(path))) {
        return(path)
    }

    # Backward compatibility for older projects.
    if (identical(path, "fcs_mapping.csv") &&
        !file.exists(path) &&
        file.exists("fcs_control_file.csv")) {
        message("Default control file 'fcs_mapping.csv' not found. Using legacy 'fcs_control_file.csv'.")
        return("fcs_control_file.csv")
    }

    path
}
