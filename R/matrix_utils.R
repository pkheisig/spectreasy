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

.default_wls_background_noise <- function() 125

.default_wls_signal_scale <- function() 1

.default_wls_max_weight_ratio <- function() 100

.default_wls_noise_tail_fraction <- function() 0.20

.coerce_wls_vector <- function(x, n_detectors, default, arg_name) {
    if (length(x) == 0 || all(is.na(x))) {
        x <- default
    }
    x <- suppressWarnings(as.numeric(x))
    if (length(x) == 1L) {
        x <- rep(x, n_detectors)
    }
    if (length(x) != n_detectors) {
        stop(arg_name, " must be a scalar or have one value per detector.")
    }
    x
}

.attach_detector_noise <- function(M, detector_noise, source = "detector_noise") {
    M <- .as_reference_matrix(M, "M")
    if (is.null(detector_noise)) {
        return(M)
    }

    noise <- as.data.frame(detector_noise, stringsAsFactors = FALSE, check.names = FALSE)
    lower_names <- tolower(colnames(noise))
    detector_idx <- match("detector", lower_names)
    noise_idx <- match("noise_floor", lower_names)
    if (is.na(noise_idx)) noise_idx <- match("background_noise", lower_names)
    if (is.na(noise_idx)) noise_idx <- match("sigma0", lower_names)
    signal_idx <- match("signal_scale", lower_names)
    if (is.na(signal_idx)) signal_idx <- match("alpha", lower_names)

    if (is.na(detector_idx) || is.na(noise_idx)) {
        stop(source, " must contain 'detector' and 'noise_floor' columns.")
    }

    detector <- trimws(as.character(noise[[detector_idx]]))
    noise_floor <- suppressWarnings(as.numeric(noise[[noise_idx]]))
    signal_scale <- if (!is.na(signal_idx)) {
        suppressWarnings(as.numeric(noise[[signal_idx]]))
    } else {
        rep(.default_wls_signal_scale(), length(detector))
    }

    keep <- nzchar(detector) & !is.na(detector)
    detector <- detector[keep]
    noise_floor <- noise_floor[keep]
    signal_scale <- signal_scale[keep]
    if (length(detector) == 0) {
        stop(source, " does not contain any detector rows.")
    }
    if (any(duplicated(detector))) {
        stop(source, " contains duplicated detector rows: ", paste(unique(detector[duplicated(detector)]), collapse = ", "))
    }

    missing <- setdiff(colnames(M), detector)
    if (length(missing) > 0) {
        stop(source, " is missing detector rows: ", paste(missing, collapse = ", "))
    }

    idx <- match(colnames(M), detector)
    attr(M, "detector_noise") <- data.frame(
        detector = colnames(M),
        noise_floor = noise_floor[idx],
        signal_scale = signal_scale[idx],
        stringsAsFactors = FALSE
    )
    M
}

.estimate_wls_detector_noise <- function(scc_dir,
                                         detectors,
                                         fallback = .default_wls_background_noise(),
                                         tail_fraction = .default_wls_noise_tail_fraction(),
                                         signal_scale = .default_wls_signal_scale()) {
    detectors <- as.character(detectors)
    if (length(detectors) == 0) {
        stop("detectors must contain at least one detector name.")
    }
    if (!dir.exists(scc_dir)) {
        stop("scc_dir not found: ", scc_dir)
    }

    tail_fraction <- suppressWarnings(as.numeric(tail_fraction)[1])
    if (!is.finite(tail_fraction) || tail_fraction <= 0 || tail_fraction >= 1) {
        tail_fraction <- .default_wls_noise_tail_fraction()
    }

    fallback <- suppressWarnings(as.numeric(fallback)[1])
    if (!is.finite(fallback) || fallback <= 0) {
        fallback <- .default_wls_background_noise()
    }

    fcs_files <- list.files(scc_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
    if (length(fcs_files) == 0) {
        stop("No FCS files found in scc_dir: ", scc_dir)
    }

    per_file <- list()
    for (path in fcs_files) {
        ff <- flowCore::read.FCS(path, transformation = FALSE, truncate_max_range = FALSE)
        expr <- flowCore::exprs(ff)
        missing <- setdiff(detectors, colnames(expr))
        if (length(missing) > 0) {
            next
        }

        per_file[[basename(path)]] <- vapply(detectors, function(det) {
            x <- as.numeric(expr[, det])
            x <- x[is.finite(x)]
            if (length(x) < 10) {
                return(NA_real_)
            }
            cutoff <- stats::quantile(x, tail_fraction, na.rm = TRUE, names = FALSE)
            low_tail <- x[x <= cutoff]
            if (length(low_tail) < 10) {
                return(NA_real_)
            }
            stats::mad(low_tail, constant = 1.4826, na.rm = TRUE)
        }, numeric(1))
    }

    if (length(per_file) == 0) {
        stop("No SCC files contained all requested detectors.")
    }

    per_file_mat <- do.call(rbind, per_file)
    noise_floor <- apply(per_file_mat, 2, stats::median, na.rm = TRUE)
    noise_floor[!is.finite(noise_floor) | noise_floor <= 0] <- fallback

    data.frame(
        detector = detectors,
        noise_floor = as.numeric(noise_floor[detectors]),
        signal_scale = rep(signal_scale, length(detectors)),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

.attach_estimated_wls_detector_noise <- function(M,
                                                scc_dir,
                                                fallback = .default_wls_background_noise(),
                                                tail_fraction = .default_wls_noise_tail_fraction(),
                                                signal_scale = .default_wls_signal_scale(),
                                                warn = TRUE) {
    M <- .as_reference_matrix(M, "M")
    detector_noise <- tryCatch(
        .estimate_wls_detector_noise(
            scc_dir = scc_dir,
            detectors = colnames(M),
            fallback = fallback,
            tail_fraction = tail_fraction,
            signal_scale = signal_scale
        ),
        error = function(e) {
            if (isTRUE(warn)) {
                warning(
                    "Could not estimate WLS detector noise from SCC files; using fallback noise floor ",
                    fallback,
                    ". Reason: ",
                    conditionMessage(e),
                    call. = FALSE
                )
            }
            NULL
        }
    )
    if (is.null(detector_noise)) {
        return(M)
    }
    .attach_detector_noise(M, detector_noise, source = "estimated detector noise")
}

.save_detector_noise_csv <- function(detector_noise, path) {
    if (is.null(detector_noise)) {
        return(invisible(NULL))
    }
    detector_noise <- as.data.frame(detector_noise, stringsAsFactors = FALSE, check.names = FALSE)
    utils::write.csv(detector_noise, path, row.names = FALSE, quote = TRUE)
    invisible(path)
}

.resolve_wls_noise_parameters <- function(M,
                                          background_noise = .default_wls_background_noise(),
                                          signal_scale = .default_wls_signal_scale(),
                                          max_weight_ratio = .default_wls_max_weight_ratio()) {
    M <- .as_reference_matrix(M, "M")
    n_detectors <- ncol(M)

    noise_floor <- .coerce_wls_vector(
        background_noise,
        n_detectors = n_detectors,
        default = .default_wls_background_noise(),
        arg_name = "background_noise"
    )
    signal_scale_vec <- .coerce_wls_vector(
        signal_scale,
        n_detectors = n_detectors,
        default = .default_wls_signal_scale(),
        arg_name = "signal_scale"
    )

    detector_noise <- attr(M, "detector_noise")
    if (!is.null(detector_noise)) {
        M_tmp <- .attach_detector_noise(M, detector_noise)
        detector_noise <- attr(M_tmp, "detector_noise")
        if (!is.null(detector_noise$noise_floor)) {
            noise_floor <- .coerce_wls_vector(
                detector_noise$noise_floor,
                n_detectors = n_detectors,
                default = .default_wls_background_noise(),
                arg_name = "detector_noise$noise_floor"
            )
        }
        if (!is.null(detector_noise$signal_scale)) {
            signal_scale_vec <- .coerce_wls_vector(
                detector_noise$signal_scale,
                n_detectors = n_detectors,
                default = .default_wls_signal_scale(),
                arg_name = "detector_noise$signal_scale"
            )
        }
    }

    noise_floor[!is.finite(noise_floor) | noise_floor <= 0] <- .default_wls_background_noise()
    signal_scale_vec[!is.finite(signal_scale_vec) | signal_scale_vec < 0] <- .default_wls_signal_scale()

    max_weight_ratio <- suppressWarnings(as.numeric(max_weight_ratio)[1])
    if (!is.finite(max_weight_ratio) || max_weight_ratio < 1) {
        max_weight_ratio <- .default_wls_max_weight_ratio()
    }

    list(
        noise_floor = as.numeric(noise_floor),
        signal_scale = as.numeric(signal_scale_vec),
        max_weight_ratio = max_weight_ratio
    )
}

.wls_event_weights <- function(y,
                               noise_floor,
                               signal_scale,
                               max_weight_ratio = .default_wls_max_weight_ratio()) {
    y <- as.numeric(y)
    noise_floor <- .coerce_wls_vector(noise_floor, length(y), .default_wls_background_noise(), "noise_floor")
    signal_scale <- .coerce_wls_vector(signal_scale, length(y), .default_wls_signal_scale(), "signal_scale")

    denom <- noise_floor + signal_scale * pmax(y, 0)
    positive_denom <- denom[is.finite(denom) & denom > 0]
    fallback_denom <- if (length(positive_denom) > 0) {
        stats::median(positive_denom, na.rm = TRUE)
    } else {
        .default_wls_background_noise()
    }
    if (!is.finite(fallback_denom) || fallback_denom <= 0) {
        fallback_denom <- .default_wls_background_noise()
    }
    denom[!is.finite(denom) | denom <= 0] <- fallback_denom
    weights <- 1 / denom

    med <- stats::median(weights[is.finite(weights) & weights > 0], na.rm = TRUE)
    if (!is.finite(med) || med <= 0) med <- 1
    weights <- weights / med

    cap <- suppressWarnings(as.numeric(max_weight_ratio)[1])
    if (is.finite(cap) && cap > 1) {
        half_cap <- sqrt(cap)
        weights <- pmin(pmax(weights, 1 / half_cap), half_cap)
    }
    weights
}

.wls_static_detector_weights <- function(M,
                                         background_noise = .default_wls_background_noise(),
                                         signal_scale = .default_wls_signal_scale(),
                                         max_weight_ratio = .default_wls_max_weight_ratio()) {
    params <- .resolve_wls_noise_parameters(
        M = M,
        background_noise = background_noise,
        signal_scale = signal_scale,
        max_weight_ratio = max_weight_ratio
    )
    .wls_event_weights(
        y = rep(0, ncol(M)),
        noise_floor = params$noise_floor,
        signal_scale = params$signal_scale,
        max_weight_ratio = params$max_weight_ratio
    )
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
    add_alias <- function(alias_map, alias, target) {
        key <- .normalize_channel_token(alias)
        if (!nzchar(key) || is.na(target) || !nzchar(target)) return(alias_map)
        if (!key %in% names(alias_map)) alias_map[[key]] <- target
        alias_map
    }

    # Direct aliases from detector names.
    for (det in detector_names) {
        det_key <- .normalize_channel_token(det)
        alias_map <- add_alias(alias_map, det_key, det)
        alias_map <- add_alias(alias_map, gsub("-A$", "", det_key), det)
        alias_map <- add_alias(alias_map, paste0(gsub("-A$", "", det_key), "-A"), det)
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
                alias_map <- add_alias(alias_map, alias_base, block$detector[i])
                alias_map <- add_alias(alias_map, paste0(alias_base, "-A"), block$detector[i])
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
    path
}

.with_optional_seed <- function(seed = NULL, .local_envir = parent.frame()) {
    if (is.null(seed)) {
        return(invisible(NULL))
    }
    seed <- as.integer(seed[1])
    if (!is.finite(seed) || is.na(seed)) {
        stop("seed must be a finite integer when provided.")
    }

    withr::local_seed(seed, .local_envir = .local_envir)
    invisible(NULL)
}
