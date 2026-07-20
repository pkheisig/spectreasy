.extract_reference_af_profiles <- function(ff_af = NULL,
                                           detector_names,
                                           n_bands = 10,
                                           max_cells = 50000,
                                           af_events = NULL) {
    if (is.null(af_events)) {
        if (is.null(ff_af)) {
            stop("Either ff_af or af_events must be provided.")
        }
        raw <- flowCore::exprs(ff_af)
        af_events <- raw[, detector_names, drop = FALSE]
    }

    af_events <- af_events[stats::complete.cases(af_events), , drop = FALSE]
    if (nrow(af_events) == 0) {
        return(list(raw_median = NULL, signatures = NULL))
    }

    if (nrow(af_events) > max_cells) {
        af_events <- af_events[sample.int(nrow(af_events), max_cells), , drop = FALSE]
    }

    raw_median <- apply(af_events, 2, stats::median, na.rm = TRUE)

    n_bands <- suppressWarnings(as.integer(n_bands[1]))
    if (!is.finite(n_bands) || is.na(n_bands) || n_bands < 1) {
        stop("n_bands must be an integer >= 1.")
    }

    af_pos <- pmax(af_events, 0)
    row_scale <- apply(af_pos, 1, max, na.rm = TRUE)
    keep <- is.finite(row_scale) & row_scale > 0
    if (!any(keep)) {
        if (n_bands != 1L) {
            stop(
                "Cannot derive exactly ", n_bands,
                " AF bands because the AF events contain no positive spectral shapes.",
                call. = FALSE
            )
        }
        sig <- pmax(raw_median, 0)
        sig_max <- max(sig, na.rm = TRUE)
        if (is.finite(sig_max) && sig_max > 0) {
            sig <- sig / sig_max
        } else {
            sig <- rep(0, length(sig))
            names(sig) <- detector_names
        }
        sig_mat <- matrix(sig, nrow = 1)
        rownames(sig_mat) <- "AF"
        colnames(sig_mat) <- detector_names
        return(list(raw_median = raw_median, signatures = sig_mat))
    }

    af_shape <- af_pos[keep, , drop = FALSE] / row_scale[keep]
    km <- .reference_kmeans_af_centers(
        af_shape = af_shape,
        n_centers = n_bands
    )
    centers <- .reference_normalize_af_centers(km$centers)
    if (ncol(centers) != length(detector_names) && nrow(centers) == length(detector_names)) {
        centers <- t(centers)
    }
    if (ncol(centers) != length(detector_names)) {
        stop(
            "AF center detector count mismatch: expected ",
            length(detector_names),
            " detector(s), got ",
            ncol(centers),
            ".",
            call. = FALSE
        )
    }
    if (nrow(centers) != n_bands) {
        stop(
            "AF extraction requested ", n_bands, " bands but clustering returned ",
            nrow(centers), ". No reduced AF bank was accepted.",
            call. = FALSE
        )
    }
    selection <- list(
        method = if (identical(km$center_method, "median")) "median_fixed" else "kmeans_fixed",
        n_bands = nrow(centers),
        requested_bands = n_bands,
        raw_center_count = nrow(km$centers),
        final_bands = nrow(centers),
        cluster_sizes = km$cluster_sizes
    )
    rownames(centers) <- c("AF", if (nrow(centers) > 1) paste0("AF_", seq.int(2, nrow(centers))) else NULL)
    colnames(centers) <- detector_names

    list(raw_median = raw_median, signatures = centers, selection = selection)
}

# Determines the name of the autofluorescence (unstained) control file.
# Looks it up in the control mapping or falls back to identifying files matching AF file naming heuristics.
# Returns the AF file basename (without extension) or NULL.
.resolve_reference_af_name <- function(control_df, fcs_files) {
    af_fn <- NULL
    if (!is.null(control_df)) {
        af_rows <- if ("fluorophore" %in% colnames(control_df)) {
            control_df[.is_primary_af_control_row(
                fluorophore = control_df$fluorophore,
                marker = if ("marker" %in% colnames(control_df)) control_df$marker else NULL,
                filename = control_df$filename
            ), ]
        } else {
            data.frame()
        }
        if (nrow(af_rows) > 0) {
            control_type <- if ("control.type" %in% colnames(af_rows)) {
                tolower(trimws(as.character(af_rows$control.type)))
            } else {
                rep("", nrow(af_rows))
            }
            bead_negative <- vapply(af_rows$filename, .reference_is_bead_negative_file, logical(1))
            af_rows <- af_rows[control_type != "beads" & !bead_negative, , drop = FALSE]
        }
        if (nrow(af_rows) > 0) af_fn <- tools::file_path_sans_ext(basename(af_rows$filename[1]))
    }
    if (is.null(af_fn)) {
        af_idx_tmp <- which(.is_af_filename(fcs_files) & !vapply(fcs_files, .reference_is_bead_negative_file, logical(1)))
        if (length(af_idx_tmp) > 0) af_fn <- tools::file_path_sans_ext(basename(fcs_files[af_idx_tmp[1]]))
    }

    af_fn
}

# Safe retriever for key values from the configuration list.
# Returns the configured value if present, otherwise returns the specified default.
.get_reference_config_value <- function(config, name, default) {
    if (!is.null(config) && name %in% names(config)) {
        return(config[[name]])
    }
    default
}

.reference_external_negative_summary <- function(detector_names,
                                                 af_data_raw = NULL,
                                                 scc_background = NULL) {
    bg <- if (!is.null(scc_background) && !is.null(scc_background$spectra)) {
        as.matrix(scc_background$spectra[, detector_names, drop = FALSE])
    } else {
        NULL
    }
    if (!is.null(bg) && nrow(bg) > 0L) {
        bg <- bg[stats::complete.cases(bg), , drop = FALSE]
        bg <- pmax(bg, 0)
        bg <- bg[rowSums(bg, na.rm = TRUE) > 0, , drop = FALSE]
        if (nrow(bg) > 0L) {
            return(list(
                mean = apply(bg, 2, mean, na.rm = TRUE),
                median = apply(bg, 2, stats::median, na.rm = TRUE),
                spectra = bg
            ))
        }
    }

    if (!is.null(af_data_raw) && all(detector_names %in% names(af_data_raw))) {
        af_vec <- pmax(as.numeric(af_data_raw[detector_names]), 0)
        names(af_vec) <- detector_names
        if (sum(af_vec^2, na.rm = TRUE) > 0) {
            return(list(
                mean = af_vec,
                median = af_vec,
                spectra = matrix(af_vec, nrow = 1L, dimnames = list(NULL, detector_names))
            ))
        }
    }

    NULL
}

.reference_cosine_to_vector <- function(mat, vec) {
    mat <- pmax(as.matrix(mat), 0)
    vec <- pmax(as.numeric(vec), 0)
    denom <- sqrt(rowSums(mat^2, na.rm = TRUE)) * sqrt(sum(vec^2, na.rm = TRUE))
    out <- rep(NA_real_, nrow(mat))
    ok <- is.finite(denom) & denom > 0
    if (any(ok)) {
        out[ok] <- as.numeric(mat[ok, , drop = FALSE] %*% vec) / denom[ok]
    }
    out
}

.reference_background_scatter_gate <- function(raw_data,
                                               pd,
                                               sample_type,
                                               filename,
                                               config,
                                               auto_sample_type = sample_type) {
    scatter_info <- .apply_reference_manual_scatter_gates(
        raw_data = raw_data,
        sample_type = sample_type,
        filename = filename,
        manual_gates = .get_reference_config_value(config, "manual_gates", NULL)
    )
    if (!is.null(scatter_info)) return(scatter_info)

    .compute_reference_scatter_gate(
        raw_data = raw_data,
        pd = pd,
        sample_type = auto_sample_type,
        outlier_percentile = .get_reference_config_value(config, "outlier_percentile", 0.02),
        debris_percentile = .get_reference_config_value(config, "debris_percentile", 0.08),
        subsample_n = .get_reference_config_value(config, "subsample_n", 5000),
        max_clusters = .get_reference_config_value(config, "max_clusters", 10),
        min_cluster_proportion = .get_reference_config_value(config, "min_cluster_proportion", 0.03),
        gate_contour_beads = .get_reference_config_value(config, "gate_contour_beads", 0.95),
        gate_contour_cells = .get_reference_config_value(config, "gate_contour_cells", 0.90),
        bead_gate_scale = .get_reference_config_value(config, "bead_gate_scale", 1.3)
    )
}

# Reads and gates an AF / unstained control FCS file.
# Applies scatter-gating on FSC-SSC to isolate the main cell/bead population from debris.
# Returns a list containing the gated event data across all spectral detectors and gating metadata.
.extract_reference_af_gated_events <- function(fcs_file,
                                               detector_names,
                                               config,
                                               sample_type = "cells") {
    ff_af <- tryCatch(
        .spectreasy_read_fcs(fcs_file),
        error = function(e) NULL
    )
    if (is.null(ff_af)) {
        warning("Could not read AF control file: ", fcs_file)
        return(NULL)
    }

    pd_af <- flowCore::pData(flowCore::parameters(ff_af))
    raw_af <- flowCore::exprs(ff_af)
    sn <- tools::file_path_sans_ext(basename(fcs_file))
    if (!.validate_reference_raw_data(raw_af, sn, detector_names = detector_names)) {
        return(NULL)
    }

    sample_type <- if (identical(tolower(trimws(as.character(sample_type)[1])), "beads")) "beads" else "cells"
    scatter_info <- .reference_background_scatter_gate(
        raw_data = raw_af,
        pd = pd_af,
        sample_type = sample_type,
        filename = basename(fcs_file),
        config = config,
        auto_sample_type = if (identical(sample_type, "beads")) "beads" else "unstained"
    )
    if (is.null(scatter_info)) {
        .spectreasy_console_step("Skip AF", paste0(basename(fcs_file), " has too few scatter-gated cells"))
        return(NULL)
    }

    list(
        events = scatter_info$gated_data[, detector_names, drop = FALSE],
        scatter = scatter_info$gated_data[, c(scatter_info$fsc, scatter_info$ssc), drop = FALSE],
        scatter_names = c(scatter_info$fsc, scatter_info$ssc),
        source = data.table::data.table(
            file = basename(fcs_file),
            path = normalizePath(fcs_file, mustWork = FALSE),
            n_total = nrow(raw_af),
            n_scatter_gated = nrow(scatter_info$gated_data),
            scatter_gate_pct = round(100 * nrow(scatter_info$gated_data) / max(nrow(raw_af), 1), 1),
            fsc_channel = scatter_info$fsc,
            ssc_channel = scatter_info$ssc
        )
    )
}

# Gathers all autofluorescence and unstained control files and extracts their signatures.
# Locates primary and extra AF files, loads and gates them, pools events, and computes
# the specified number of AF basis signatures.
# Returns a list containing raw medians, normalized basis matrices, and banking metadata.
.collect_reference_af_profiles <- function(control_df,
                                           fcs_files,
                                           detector_names,
                                           af_n_bands,
                                           af_max_cells,
                                           fcs_files_all = fcs_files,
                                           config = NULL) {
    af_data_raw <- NULL
    af_signatures_norm <- NULL
    af_bank_info <- NULL
    scc_background <- NULL
    af_events <- NULL
    af_fn <- .resolve_reference_af_name(control_df = control_df, fcs_files = fcs_files)

    af_paths <- character()
    af_source_types <- character()
    fcs_keys <- tools::file_path_sans_ext(basename(fcs_files_all))
    if (!is.null(control_df) && nrow(control_df) > 0) {
        af_rows <- .is_primary_af_control_row(
            fluorophore = if ("fluorophore" %in% colnames(control_df)) control_df$fluorophore else NULL,
            marker = if ("marker" %in% colnames(control_df)) control_df$marker else NULL,
            filename = control_df$filename
        )
        if (any(af_rows)) {
            af_df <- control_df[af_rows, , drop = FALSE]
            control_type <- if ("control.type" %in% colnames(af_df)) {
                tolower(trimws(as.character(af_df$control.type)))
            } else {
                rep("", nrow(af_df))
            }
            bead_negative <- vapply(af_df$filename, .reference_is_bead_negative_file, logical(1))
            af_df <- af_df[control_type != "beads" & !bead_negative, , drop = FALSE]
            if (nrow(af_df) > 0) {
                af_file_keys <- tools::file_path_sans_ext(basename(as.character(af_df$filename)))
                hits <- match(af_file_keys, fcs_keys)
                hits <- hits[!is.na(hits)]
                if (length(hits) > 0) {
                    af_paths <- c(af_paths, fcs_files_all[hits])
                    af_source_types <- c(af_source_types, rep("mapped_unstained", length(hits)))
                }
            }
        }
    }
    if (length(af_paths) == 0 && !is.null(af_fn)) {
        af_path <- fcs_files_all[match(af_fn, fcs_keys)]
        af_path <- af_path[!is.na(af_path)]
        if (length(af_path) > 0) {
            af_paths <- c(af_paths, af_path[1])
            af_source_types <- c(af_source_types, "primary_unstained")
        }
    }
    keep_unique <- !duplicated(normalizePath(af_paths, mustWork = FALSE))
    af_paths <- af_paths[keep_unique]
    af_source_types <- af_source_types[keep_unique]

    if (length(af_paths) > 0) {
        af_gated_list <- lapply(
            af_paths,
            .extract_reference_af_gated_events,
            detector_names = detector_names,
            config = config,
            sample_type = "cells"
        )
        keep_gated <- vapply(af_gated_list, function(x) !is.null(x) && !is.null(x$events) && nrow(x$events) > 0, logical(1))
        af_gated_list <- af_gated_list[keep_gated]
        af_source_types <- af_source_types[keep_gated]
        if (length(af_gated_list) == 0) {
            return(list(
                af_data_raw = af_data_raw,
                af_signatures_norm = af_signatures_norm,
                af_bank_info = af_bank_info,
                scc_background = scc_background,
                af_events = af_events
            ))
        }

        af_events <- do.call(rbind, lapply(af_gated_list, `[[`, "events"))
        scc_background <- .scc_background_from_gated_af_list(
            af_gated_list = af_gated_list,
            detector_names = detector_names
        )
        n_af_sources <- length(af_gated_list)
        requested_bands <- af_n_bands
        af_profiles <- .extract_reference_af_profiles(
            detector_names = detector_names,
            n_bands = requested_bands,
            max_cells = af_max_cells,
            af_events = af_events
        )
        af_data_raw <- af_profiles$raw_median
        af_signatures_norm <- af_profiles$signatures
        af_sources <- data.table::rbindlist(lapply(af_gated_list, `[[`, "source"))
        af_sources$source_type <- af_source_types
        data.table::setcolorder(af_sources, c("file", "source_type", "n_total", "n_scatter_gated", "scatter_gate_pct", "fsc_channel", "ssc_channel", "path"))
        af_bank_info <- list(
            source_count = n_af_sources,
            sources = af_sources,
            pooled_events = nrow(af_events),
            requested_bands = requested_bands,
            derived_bands = if (!is.null(af_signatures_norm)) nrow(af_signatures_norm) else 0L,
            selection = af_profiles$selection,
            mode = if (n_af_sources > 1) "pooled_af_sources" else "single_af"
        )
        if (!is.null(af_signatures_norm)) {
            msg <- if (n_af_sources == 1) {
                "primary unstained control"
            } else {
                paste0(n_af_sources, " pooled AF control files")
            }
            selection_msg <- if (!is.null(af_profiles$selection)) {
                paste0(" (selected by ", af_profiles$selection$method, ")")
            } else {
                ""
            }
            .spectreasy_console_field(
                "AF bank",
                paste0(nrow(af_signatures_norm), " signature(s) from ", msg, selection_msg)
            )
        }
    }

    list(
        af_data_raw = af_data_raw,
        af_signatures_norm = af_signatures_norm,
        af_bank_info = af_bank_info,
        scc_background = scc_background,
        af_events = af_events
    )
}
