.get_control_rows_for_reference <- function(df, filenames) {
    if (is.null(df) || !("filename" %in% colnames(df))) {
        return(data.frame())
    }
    fn <- as.character(df$filename)
    df[fn %in% filenames, ]
}

# Prepares the complete list of FCS files to be processed.
# Locates FCS files in the input folder.
# Returns a list containing 'fcs_files' and 'fcs_files_all'.
.prepare_reference_file_set <- function(input_folder, control_df = NULL) {
    if (!dir.exists(input_folder)) {
        .spectreasy_stop_missing_directory(input_folder, label = "input_folder")
    }
    fcs_files_all <- list.files(input_folder, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
    if (length(fcs_files_all) == 0) {
        .spectreasy_stop_empty_fcs_directory(input_folder, label = "input_folder")
    }
    fcs_files_all <- .control_validation_select_scc_files(control_df, fcs_files_all)
    if (length(fcs_files_all) == 0) {
        stop("None of the FCS files listed in control_df were found in input_folder.", call. = FALSE)
    }

    list(fcs_files = fcs_files_all, fcs_files_all = fcs_files_all)
}

# Prepares and returns the normalized path of the output directory.
# If save_qc_plots is TRUE, recursively creates subdirectories for FSC/SSC plots,
# intensity gates, and spectra.
# Returns the absolute path string.
.prepare_reference_output_path <- function(output_folder, save_qc_plots = FALSE) {
    out_path <- normalizePath(output_folder, mustWork = FALSE)
    if (isTRUE(save_qc_plots)) {
        dir.create(file.path(out_path, "fsc_ssc"), showWarnings = FALSE, recursive = TRUE)
        dir.create(file.path(out_path, "singlet"), showWarnings = FALSE, recursive = TRUE)
        dir.create(file.path(out_path, "histogram"), showWarnings = FALSE, recursive = TRUE)
        dir.create(file.path(out_path, "spectrum"), showWarnings = FALSE, recursive = TRUE)
    }
    out_path
}

# Helper to normalize a channel name by converting to uppercase, removing whitespace,
# and handling NA values.
# Returns the normalized character vector.
.normalize_reference_channel <- function(x) {
    out <- toupper(gsub("\\s+", "", trimws(as.character(x))))
    out <- gsub("([A-Z]+)-([0-9])", "\\1\\2", out, perl = TRUE)
    out[is.na(out)] <- ""
    out
}

# Resolves a fluorophore's target channel name to an actual detector in the FCS file.
# Normalizes string cases/whitespace and tries matching the channel name directly, with/without
# standard suffix/prefix options (e.g., "-A"), or using a cytometer-specific channel alias map.
# Returns the matched detector name string or an empty string if unresolved.
.resolve_reference_control_channel <- function(channel_value, det_names, channel_alias_map = character()) {
    if (is.null(channel_value) || is.na(channel_value) || trimws(channel_value) == "") return("")
    if (channel_value %in% det_names) return(channel_value)

    det_norm <- .normalize_reference_channel(det_names)
    name_by_norm <- stats::setNames(det_names, det_norm)
    key <- .normalize_reference_channel(channel_value)
    candidates <- unique(c(
        key,
        gsub("-A$", "", key),
        paste0(key, "-A"),
        if (length(channel_alias_map) > 0 && key %in% names(channel_alias_map)) channel_alias_map[[key]] else NULL
    ))
    candidates <- candidates[nzchar(candidates)]
    for (cand in candidates) {
        if (cand %in% det_norm) return(name_by_norm[[cand]])
    }
    ""
}

# Extracts detector information from the first FCS file in the set.
# Reads the file metadata to determine detector names, labels, and channel mappings.
# Returns a list containing the parameter data table (pd_meta), sorted detector info, and alias map.
.prepare_reference_detector_info <- function(first_fcs_file) {
    ff_meta <- .spectreasy_read_fcs(first_fcs_file, label = "SCC FCS file")
    pd_meta <- flowCore::pData(flowCore::parameters(ff_meta))
    det_info <- get_sorted_detectors(pd_meta)

    list(
        pd_meta = pd_meta,
        det_info = det_info,
        detector_names = det_info$names,
        detector_labels = det_info$labels,
        channel_alias_map = .build_channel_alias_map_from_pd(pd_meta)
    )
}

# Automatically identifies the primary Forward Scatter (FSC) and Side Scatter (SSC) channels.
# Matches names in the raw data columns using standard patterns for FSC and SSC.
# Returns a list with 'fsc' and 'ssc' channel names, or NULL if they cannot be resolved.
.resolve_reference_scatter_channels <- function(raw_data) {
    primary <- .get_primary_scatter_channels(colnames(raw_data))

    if (is.na(primary$fsc) || is.na(primary$ssc) || !all(c(primary$fsc, primary$ssc) %in% colnames(raw_data))) {
        return(NULL)
    }

    primary
}

.read_reference_manual_gates <- function(manual_gate_file = NULL) {
    if (is.null(manual_gate_file) || length(manual_gate_file) == 0 || is.na(manual_gate_file[1])) return(NULL)
    path <- normalizePath(as.character(manual_gate_file)[1], mustWork = FALSE)
    if (!file.exists(path)) return(NULL)
    df <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    required <- c("gate_type", "scope", "filename", "x_channel", "y_channel", "plot_mode", "vertex_index", "x", "y")
    if (is.null(df) || nrow(df) == 0) {
        stop("Could not read a valid manual gate CSV: ", path, call. = FALSE)
    }
    missing_cols <- setdiff(required, colnames(df))
    if (length(missing_cols) > 0) {
        stop("Manual gate CSV is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    settings <- list()
    setting_rows <- df[df$gate_type == "setting", , drop = FALSE]
    if (nrow(setting_rows) > 0) {
        for (i in seq_len(nrow(setting_rows))) {
            key <- trimws(as.character(setting_rows$x_channel[i]))
            value <- trimws(as.character(setting_rows$x[i]))
            if (!nzchar(key)) next
            if (identical(key, "point_size")) {
                settings$point_size <- suppressWarnings(as.numeric(value))
            } else if (identical(key, "max_points")) {
                settings$max_points <- suppressWarnings(as.integer(value))
            } else if (identical(key, "histogram_bins")) {
                settings$histogram_bins <- suppressWarnings(as.integer(value))
            } else if (identical(key, "histogram_transform")) {
                settings$histogram_transform <- tolower(value)
            }
        }
    }
    keep <- df$gate_type != "setting" &
        (df$plot_mode == "blocked" | (df$plot_mode != "missing" & suppressWarnings(as.integer(df$vertex_index)) > 0))
    df <- df[keep, , drop = FALSE]
    if (nrow(df) == 0) {
        out <- list()
    } else {
        out <- split(df, paste(df$gate_type, df$scope, df$filename, sep = "\r"))
    }
    attr(out, "settings") <- settings
    out
}

.reference_manual_gate_settings <- function(manual_gates) {
    settings <- if (!is.null(manual_gates)) attr(manual_gates, "settings") else NULL
    if (is.null(settings)) settings <- list()
    bins_value <- if (is.null(settings$histogram_bins)) 100L else settings$histogram_bins
    bins <- suppressWarnings(as.integer(bins_value))
    if (!is.finite(bins) || is.na(bins)) bins <- 100L
    bins <- min(max(bins, 5L), 500L)
    transform_value <- if (is.null(settings$histogram_transform)) "auto" else settings$histogram_transform
    transform <- tolower(trimws(as.character(transform_value)))
    if (!transform %in% c("auto", "asinh", "linear", "log10", "biexponential")) transform <- "auto"
    max_points_value <- if (is.null(settings$max_points)) 50000L else settings$max_points
    max_points <- suppressWarnings(as.integer(max_points_value))
    if (!is.finite(max_points) || is.na(max_points) || max_points <= 0L) max_points <- 50000L
    point_size_value <- if (is.null(settings$point_size)) 1.5 else settings$point_size
    point_size <- suppressWarnings(as.numeric(point_size_value))
    if (!is.finite(point_size) || is.na(point_size) || point_size <= 0) point_size <- 1.5
    list(
        histogram_bins = bins,
        histogram_transform = transform,
        max_points = max_points,
        point_size = point_size
    )
}

.reference_manual_gate_key <- function(gate_type, scope, filename = "") {
    paste(gate_type, scope, filename, sep = "\r")
}

.reference_manual_gate <- function(manual_gates, gate_type, filename, sample_type) {
    if (is.null(manual_gates)) return(NULL)
    keys <- c(
        .reference_manual_gate_key(gate_type, "file", filename),
        .reference_manual_gate_key(gate_type, sample_type, ""),
        .reference_manual_gate_key(gate_type, "cells", "")
    )
    hit <- keys[keys %in% names(manual_gates)][1]
    if (is.na(hit) || !nzchar(hit)) return(NULL)
    gate <- manual_gates[[hit]]
    gate[order(suppressWarnings(as.integer(gate$vertex_index))), , drop = FALSE]
}

.reference_gate_vertices <- function(gate) {
    if (is.null(gate) || nrow(gate) == 0) return(NULL)
    out <- data.frame(
        x = suppressWarnings(as.numeric(gate$x)),
        y = suppressWarnings(as.numeric(gate$y))
    )
    out <- out[stats::complete.cases(out), , drop = FALSE]
    if (nrow(out) == 0) return(NULL)
    out
}

.reference_gate_channel <- function(gate, column, fallback, raw_data) {
    channel <- if (!is.null(gate) && column %in% colnames(gate)) as.character(gate[[column]][1]) else fallback
    if (is.na(channel) || !nzchar(channel) || !channel %in% colnames(raw_data)) fallback else channel
}

.apply_reference_manual_scatter_gates <- function(raw_data, sample_type, filename, manual_gates) {
    scatter <- .resolve_reference_scatter_channels(raw_data)
    if (is.null(scatter)) return(NULL)
    fsc <- scatter$fsc
    ssc <- scatter$ssc
    fsc_h <- sub("-A$", "-H", fsc, ignore.case = TRUE)
    if (!fsc_h %in% colnames(raw_data)) {
        fsc_h <- grep("^FSC.*-H$", colnames(raw_data), value = TRUE, ignore.case = TRUE)[1]
    }

    keep <- rep(TRUE, nrow(raw_data))
    cell_keep <- rep(TRUE, nrow(raw_data))
    singlet_keep <- rep(TRUE, nrow(raw_data))
    final_gate <- NULL
    cell_gate <- .reference_manual_gate(manual_gates, "cell", filename, sample_type)
    cell_vertices <- .reference_gate_vertices(cell_gate)
    cell_x <- fsc
    cell_y <- ssc
    if (!is.null(cell_vertices) && nrow(cell_vertices) >= 3) {
        cell_x <- .reference_gate_channel(cell_gate, "x_channel", fsc, raw_data)
        cell_y <- .reference_gate_channel(cell_gate, "y_channel", ssc, raw_data)
        cell_keep <- sp::point.in.polygon(raw_data[, cell_x], raw_data[, cell_y], cell_vertices$x, cell_vertices$y) > 0
        keep <- keep & cell_keep
        if (identical(cell_x, fsc) && identical(cell_y, ssc)) {
            final_gate <- cell_vertices
        }
    }
    singlet_gate <- .reference_manual_gate(manual_gates, "singlet", filename, sample_type)
    singlet_vertices <- .reference_gate_vertices(singlet_gate)
    manual_attempted <- !is.null(cell_gate) || !is.null(singlet_gate)
    singlet_x <- fsc_h
    singlet_y <- fsc
    if (!is.null(singlet_vertices) && nrow(singlet_vertices) >= 3 && !is.na(fsc_h) && fsc_h %in% colnames(raw_data)) {
        singlet_x <- .reference_gate_channel(singlet_gate, "x_channel", fsc_h, raw_data)
        singlet_y <- .reference_gate_channel(singlet_gate, "y_channel", fsc, raw_data)
        singlet_keep <- sp::point.in.polygon(raw_data[, singlet_x], raw_data[, singlet_y], singlet_vertices$x, singlet_vertices$y) > 0
        keep <- keep & singlet_keep
    }
    if (sum(keep, na.rm = TRUE) < 10) return(NULL)
    if (all(keep) && !manual_attempted) return(NULL)
    gated_data <- raw_data[keep, , drop = FALSE]
    if (is.null(final_gate)) {
        x_min <- min(gated_data[, fsc], na.rm = TRUE)
        x_max <- max(gated_data[, fsc], na.rm = TRUE)
        y_min <- min(gated_data[, ssc], na.rm = TRUE)
        y_max <- max(gated_data[, ssc], na.rm = TRUE)
        final_gate <- data.frame(
            x = c(x_min, x_max, x_max, x_min, x_min),
            y = c(y_min, y_min, y_max, y_max, y_min)
        )
    }
    list(
        gated_data = gated_data,
        final_gate = final_gate,
        fsc = fsc,
        ssc = ssc,
        fsc_max = stats::quantile(raw_data[, fsc], 0.98, na.rm = TRUE),
        ssc_max = stats::quantile(raw_data[, ssc], 0.98, na.rm = TRUE),
        manual = TRUE,
        manual_gate_info = list(
            cell = list(vertices = cell_vertices, x_channel = cell_x, y_channel = cell_y, keep_count = sum(cell_keep, na.rm = TRUE)),
            singlet = list(
                vertices = singlet_vertices,
                x_channel = singlet_x,
                y_channel = singlet_y,
                source_count = sum(cell_keep, na.rm = TRUE),
                keep_count = sum(cell_keep & singlet_keep, na.rm = TRUE)
            ),
            settings = .reference_manual_gate_settings(manual_gates)
        )
    )
}

.apply_reference_manual_positive_gate <- function(gated_data, peak_channel, filename, sample_type, manual_gates) {
    gate <- .reference_manual_gate(manual_gates, "positive", filename, sample_type)
    verts <- .reference_gate_vertices(gate)
    if (is.null(gate) || is.null(verts) || nrow(verts) == 0) return(NULL)
    mode <- as.character(gate$plot_mode[1])
    peak_vals <- gated_data[, peak_channel]
    keep <- NULL
    if (identical(mode, "separator")) {
        keep <- peak_vals >= verts$x[1]
    } else if (identical(mode, "positive_1d") && nrow(verts) >= 2) {
        lim <- range(verts$x, na.rm = TRUE)
        keep <- peak_vals >= lim[1] & peak_vals <= lim[2]
    } else if (nrow(verts) >= 3) {
        x_channel <- as.character(gate$x_channel[1])
        y_channel <- as.character(gate$y_channel[1])
        if (!x_channel %in% colnames(gated_data)) x_channel <- peak_channel
        if (!y_channel %in% colnames(gated_data)) {
            scatter <- .resolve_reference_scatter_channels(gated_data)
            y_channel <- if (!is.null(scatter)) scatter$fsc else peak_channel
        }
        keep <- sp::point.in.polygon(gated_data[, x_channel], gated_data[, y_channel], verts$x, verts$y) > 0
    }
    if (is.null(keep) || sum(keep, na.rm = TRUE) < 10) return(NULL)
    vals_log <- log10(pmax(peak_vals, 1))
    attr(vals_log, "gate_type") <- paste0("manual_", mode)
    attr(vals_log, "gate_method") <- paste0("manual_", mode)
    attr(vals_log, "positive_gate_present") <- TRUE
    if (identical(mode, "positive_1d") && nrow(verts) >= 2) {
        positive_lim <- range(verts$x, na.rm = TRUE)
        if (all(is.finite(positive_lim)) && positive_lim[2] > positive_lim[1]) {
            attr(vals_log, "pos_raw_min") <- positive_lim[1]
            attr(vals_log, "pos_raw_max") <- positive_lim[2]
        }
    }
    negative_gate <- .reference_manual_gate(manual_gates, "negative", filename, sample_type)
    negative_verts <- .reference_gate_vertices(negative_gate)
    if (!is.null(negative_gate) && !is.null(negative_verts) && nrow(negative_verts) >= 2) {
        negative_mode <- as.character(negative_gate$plot_mode[1])
        if (identical(negative_mode, "negative_1d")) {
            negative_lim <- range(negative_verts$x, na.rm = TRUE)
            if (all(is.finite(negative_lim)) && negative_lim[2] > negative_lim[1]) {
                attr(vals_log, "negative_gate_present") <- TRUE
                attr(vals_log, "neg_log_min") <- log10(max(negative_lim[1], 1))
                attr(vals_log, "neg_log_max") <- log10(max(negative_lim[2], 1))
                attr(vals_log, "neg_raw_min") <- negative_lim[1]
                attr(vals_log, "neg_raw_max") <- negative_lim[2]
            }
        }
    }
    list(
        final_gated_data = gated_data[keep, , drop = FALSE],
        hist_info = list(
            vals_log = vals_log,
            gate_min = min(peak_vals[keep], na.rm = TRUE),
            gate_max = max(peak_vals[keep], na.rm = TRUE),
            positive_raw_min = attr(vals_log, "pos_raw_min"),
            positive_raw_max = attr(vals_log, "pos_raw_max"),
            negative_raw_min = attr(vals_log, "neg_raw_min"),
            negative_raw_max = attr(vals_log, "neg_raw_max")
        )
    )
}

.resolve_reference_positive_histogram_gate <- function(gated_data,
                                                       peak_channel,
                                                       filename,
                                                       sample_type,
                                                       manual_gates,
                                                       config,
                                                       row_info) {
    manual <- .apply_reference_manual_positive_gate(
        gated_data = gated_data,
        peak_channel = peak_channel,
        filename = filename,
        sample_type = sample_type,
        manual_gates = manual_gates
    )
    if (!is.null(manual)) return(manual)

    peak_vals <- gated_data[, peak_channel]
    hist_info <- .compute_reference_histogram_gate(
        peak_vals = peak_vals,
        sample_type = sample_type,
        histogram_pct_beads = config$histogram_pct_beads,
        histogram_direction_beads = config$histogram_direction_beads,
        histogram_pct_cells = config$histogram_pct_cells,
        histogram_direction_cells = config$histogram_direction_cells,
        is_viability = nrow(row_info) > 0 &&
            "is.viability" %in% colnames(row_info) &&
            toupper(trimws(as.character(row_info$is.viability[1]))) == "TRUE"
    )
    keep <- is.finite(peak_vals) &
        peak_vals >= hist_info$gate_min & peak_vals <= hist_info$gate_max
    keep[is.na(keep)] <- FALSE
    if (sum(keep) < 10L) return(NULL)
    list(
        final_gated_data = gated_data[keep, , drop = FALSE],
        hist_info = hist_info
    )
}

.reference_histogram_negative_source <- function(gated_data,
                                                  peak_channel,
                                                  filename,
                                                  sample_type,
                                                  manual_gates,
                                                  detector_names,
                                                  fsc,
                                                  ssc,
                                                  hist_info = NULL) {
    if (!all(c(peak_channel, detector_names, fsc, ssc) %in% colnames(gated_data))) return(NULL)

    gate <- .reference_manual_gate(manual_gates, "negative", filename, sample_type)
    verts <- .reference_gate_vertices(gate)
    keep <- NULL
    source <- "lower_tail"
    if (!is.null(gate) && !is.null(verts) && nrow(verts) >= 2L &&
        identical(as.character(gate$plot_mode[1]), "negative_1d")) {
        limits <- range(verts$x, na.rm = TRUE)
        if (all(is.finite(limits)) && limits[2] > limits[1]) {
            keep <- gated_data[, peak_channel] >= limits[1] & gated_data[, peak_channel] <= limits[2]
            source <- "manual_negative_gate"
        }
    }
    if (is.null(keep) && !is.null(hist_info) && !is.null(hist_info$vals_log)) {
        neg_min <- as.numeric(attr(hist_info$vals_log, "neg_raw_min", exact = TRUE))[1]
        neg_max <- as.numeric(attr(hist_info$vals_log, "neg_raw_max", exact = TRUE))[1]
        if (!is.finite(neg_min) || !is.finite(neg_max) || neg_max <= neg_min) {
            neg_min_log <- as.numeric(attr(hist_info$vals_log, "neg_log_min", exact = TRUE))[1]
            neg_max_log <- as.numeric(attr(hist_info$vals_log, "neg_log_max", exact = TRUE))[1]
            if (is.finite(neg_min_log) && is.finite(neg_max_log) && neg_max_log > neg_min_log) {
                neg_min <- 10^neg_min_log
                neg_max <- 10^neg_max_log
            }
        }
        if (is.finite(neg_min) && is.finite(neg_max) && neg_max > neg_min) {
            keep <- gated_data[, peak_channel] >= neg_min & gated_data[, peak_channel] <= neg_max
            source <- "automatic_negative_gate"
        }
    }
    if (is.null(keep) || sum(keep, na.rm = TRUE) < 10L) {
        peak_vals <- gated_data[, peak_channel]
        finite <- is.finite(peak_vals)
        if (sum(finite) < 10L) return(NULL)
        cutoff <- as.numeric(stats::quantile(peak_vals[finite], 0.15, na.rm = TRUE))
        keep <- finite & peak_vals <= cutoff
        source <- "lower_tail"
    }
    keep[is.na(keep)] <- FALSE
    negative_events <- gated_data[keep, , drop = FALSE]
    if (nrow(negative_events) < 10L) return(NULL)

    negative <- apply(negative_events[, detector_names, drop = FALSE], 2, stats::median, na.rm = TRUE)
    background <- .scc_background_from_gated_af_list(
        af_gated_list = list(list(
            events = negative_events[, detector_names, drop = FALSE],
            scatter = negative_events[, c(fsc, ssc), drop = FALSE],
            scatter_names = c(fsc, ssc)
        )),
        detector_names = detector_names
    )
    if (!is.null(background)) attr(negative, "scc_background") <- background
    attr(negative, "source") <- source
    negative
}

# Computes the 2D ellipse coordinates at a given confidence level.
# Used for gating populations in FSC/SSC scatter plots based on GMM component mean and variance.
# Returns a data.table containing the ellipse points.
.get_reference_ellipse <- function(mean, sigma, level = 0.95, n = 100, scale = 1.0) {
    if (length(mean) < 2 || any(!is.finite(mean)) ||
        !is.matrix(sigma) || any(dim(sigma) < 2) || any(!is.finite(sigma))) {
        return(NULL)
    }
    chi2_val <- qchisq(level, df = 2)
    eig <- eigen(sigma)
    eig_values <- pmax(Re(eig$values), 0)
    if (length(eig_values) < 2 || any(!is.finite(eig_values)) || sum(eig_values > 0) == 0) {
        return(NULL)
    }
    a <- sqrt(eig_values[1] * chi2_val) * scale
    b <- sqrt(eig_values[2] * chi2_val) * scale
    angle <- atan2(eig$vectors[2, 1], eig$vectors[1, 1])
    theta <- seq(0, 2 * pi, length.out = n)
    ellipse_x <- a * cos(theta)
    ellipse_y <- b * sin(theta)
    rot_x <- ellipse_x * cos(angle) - ellipse_y * sin(angle)
    rot_y <- ellipse_x * sin(angle) + ellipse_y * cos(angle)
    data.table::data.table(x = rot_x + mean[1], y = rot_y + mean[2])
}

# Fits a Gaussian Mixture Model (GMM) on scatter data using mclust.
# Tests components from 1 to max_k to automatically cluster cells or beads based on FSC/SSC density.
# Returns a list with the model fit, mixture proportions, cluster means, covariances, and main components.
.fit_reference_gmm_populations <- function(data, max_k = 5, min_prop = 0.05) {
    mclustBIC <- get("mclustBIC", envir = asNamespace("mclust"))
    fit <- mclust::Mclust(data, G = seq_len(max_k), verbose = FALSE)
    if (is.null(fit)) {
        return(NULL)
    }
    sigmas <- lapply(seq_len(fit$G), function(k) {
        if (fit$modelName %in% c("EII", "VII")) {
            diag(fit$parameters$variance$sigmasq[k], nrow = 2)
        } else {
            fit$parameters$variance$sigma[, , k]
        }
    })
    list(
        fit = fit,
        proportions = fit$parameters$pro,
        means = fit$parameters$mean,
        sigmas = sigmas,
        main_populations = which(fit$parameters$pro >= min_prop)
    )
}

# Resolves sample type (cells/beads/unstained) from the FCS filename based on pre-defined regex patterns.
# Returns a list with the detected type and matching pattern.
.get_reference_sample_type <- function(filename, patterns, default) {
    for (type in names(patterns)) {
        pats <- patterns[[type]][order(-nchar(patterns[[type]]))]
        for (p in pats) {
            if (grepl(p, filename, fixed = FALSE, ignore.case = TRUE)) {
                return(list(type = type, pattern = p))
            }
        }
    }
    list(type = default, pattern = "default")
}

# Resolves control sample type, giving priority to control_df metadata over filename heuristics.
# Returns the final sample type configuration.
.resolve_reference_sample_type <- function(filename, row_info, patterns, default) {
    sample_info <- .get_reference_sample_type(filename, patterns, default)
    if (nrow(row_info) > 0 && "control.type" %in% colnames(row_info)) {
        control_type <- tolower(trimws(as.character(row_info$control.type[1])))
        if (control_type %in% c("beads", "cells")) {
            sample_info$type <- control_type
        }
    }
    sample_info
}

# Selects the most likely bead population from GMM components.
# Uses local density in a normalized FSC/SSC space to identify the actual peak/bead population
# rather than debris or background.
# Returns a list with the index of the selected population.
.select_reference_bead_population <- function(gmm_result) {
    if (length(gmm_result$main_populations) == 0) {
        return(NULL)
    }

    # Calculate local density of raw data around each GMM mean to ensure we select the actual peak
    data_fit <- gmm_result$fit$data
    means <- gmm_result$means

    sd_fsc <- stats::sd(data_fit[, 1])
    sd_ssc <- stats::sd(data_fit[, 2])

    # Avoid division by zero
    if (is.na(sd_fsc) || sd_fsc == 0) sd_fsc <- 1
    if (is.na(sd_ssc) || sd_ssc == 0) sd_ssc <- 1

    # Normalized data
    data_norm <- data_fit
    data_norm[, 1] <- data_fit[, 1] / sd_fsc
    data_norm[, 2] <- data_fit[, 2] / sd_ssc

    # Compute local density within r = 0.10 in normalized space
    local_densities <- sapply(gmm_result$main_populations, function(k) {
        mu_norm <- c(means[1, k] / sd_fsc, means[2, k] / sd_ssc)
        dists <- sqrt((data_norm[, 1] - mu_norm[1])^2 + (data_norm[, 2] - mu_norm[2])^2)
        sum(dists <= 0.10)
    })

    best <- gmm_result$main_populations[which.max(local_densities)]
    list(selected = best)
}

# Selects valid cell populations from GMM components.
# Filters components using minimum and maximum FSC/SSC thresholds and ratio boundaries
# to exclude debris and keep only the intact single-cell population.
# Returns a list with the indices of all valid cell populations.
.select_reference_cell_populations <- function(gmm_result,
                                               fsc_min,
                                               fsc_max = Inf,
                                               ssc_max = Inf,
                                               ratio_max = Inf) {
    valid <- c()
    means <- gmm_result$means
    fsc <- means[1, ]
    ssc <- means[2, ]
    ratio <- ssc / pmax(fsc, 1)
    if (!is.finite(fsc_max)) fsc_max <- Inf
    if (!is.finite(ssc_max)) ssc_max <- Inf
    if (!is.finite(ratio_max)) ratio_max <- Inf

    for (k in gmm_result$main_populations) {
        if (fsc[k] >= fsc_min &&
            fsc[k] <= fsc_max &&
            ssc[k] <= ssc_max &&
            ratio[k] <= ratio_max) {
            valid <- c(valid, k)
        }
    }
    list(selected = valid)
}

# Creates a polygon gate enclosing the selected GMM populations.
# Computes confidence ellipses at the specified level. If multiple populations are selected,
# computes the convex hull around the union of their ellipses. Optionally clips the gate boundaries.
# Returns a data.table of (x, y) coordinates representing the gate polygon.
.create_reference_merged_gate <- function(gmm_result, populations, level, scale = 1.0, clip_x = Inf, clip_y = Inf) {
    if (length(populations) == 0) {
        return(NULL)
    }
    if (length(populations) == 1) {
        ell <- .get_reference_ellipse(gmm_result$means[, populations], gmm_result$sigmas[[populations]], level, scale = scale)
    } else {
        ellipses <- lapply(populations, function(k) {
            .get_reference_ellipse(gmm_result$means[, k], gmm_result$sigmas[[k]], level, scale = scale)
        })
        ellipses <- ellipses[!vapply(ellipses, is.null, logical(1))]
        if (length(ellipses) == 0) {
            return(NULL)
        }
        all_pts <- data.table::rbindlist(
            ellipses
        )
        hull_idx <- grDevices::chull(all_pts$x, all_pts$y)
        ell <- all_pts[hull_idx, ]
    }
    if (is.null(ell) || nrow(ell) == 0) {
        return(NULL)
    }
    ell$x <- pmax(0, pmin(ell$x, clip_x))
    ell$y <- pmax(0, pmin(ell$y, clip_y))
    if (nrow(ell) > 0 && (ell$x[1] != ell$x[nrow(ell)] || ell$y[1] != ell$y[nrow(ell)])) {
        ell <- rbind(ell, ell[1, , drop = FALSE])
    }
    ell
}

# Isolates the peak population from a 1D density distribution of log-transformed values.
# Detects peaks and troughs in the density curve, selects the rightmost peak exceeding
# the height cutoff, and extracts values bounded by the adjacent troughs.
# Returns the vector of values corresponding to this peak population.
.select_reference_hist_peak_population <- function(vals_log, rel_height_cutoff = 0.15, min_points = 20) {
    vals_log <- vals_log[is.finite(vals_log)]
    if (length(vals_log) < min_points) return(vals_log)

    d <- density(vals_log, n = 1024)
    x <- d$x
    y <- d$y

    peak_idx <- which(diff(sign(diff(y))) == -2) + 1
    trough_idx <- which(diff(sign(diff(y))) == 2) + 1
    if (length(peak_idx) == 0) return(vals_log)

    h_cut <- max(y, na.rm = TRUE) * rel_height_cutoff
    keep_peaks <- peak_idx[y[peak_idx] >= h_cut]
    if (length(keep_peaks) == 0) keep_peaks <- peak_idx
    sel_peak <- keep_peaks[which.max(x[keep_peaks])]

    left_trough <- trough_idx[trough_idx < sel_peak]
    right_trough <- trough_idx[trough_idx > sel_peak]
    left_x <- if (length(left_trough) > 0) x[max(left_trough)] else min(x, na.rm = TRUE)
    right_x <- if (length(right_trough) > 0) x[min(right_trough)] else max(x, na.rm = TRUE)

    peak_vals <- vals_log[vals_log >= left_x & vals_log <= right_x]
    if (length(peak_vals) < min_points) return(vals_log)
    peak_vals
}

# Extracts normalized Autofluorescence (AF) profiles from gated unstained events.
# Uses the median normalized shape for one AF band and k-means on scale-free
# spectral shapes for multi-band AF banks.
# Each center is normalized by its peak detector to yield a relative signature from 0 to 1.
# Returns a list with the raw median spectrum and a matrix of normalized AF basis signatures.
