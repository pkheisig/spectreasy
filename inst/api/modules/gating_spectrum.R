gate_detector_labels <- function(det_info, spec_channels) {
    labels <- det_info$labels[match(spec_channels, det_info$names)]
    labels[is.na(labels) | !nzchar(labels)] <- spec_channels[is.na(labels) | !nzchar(labels)]
    labels
}

gate_detector_signature <- function(pd) {
    desc <- if ("desc" %in% colnames(pd)) as.character(pd$desc) else rep("", nrow(pd))
    list(name = as.character(pd$name), desc = desc)
}

gate_detector_info <- function(ff, compute_fun = NULL) {
    pd <- flowCore::pData(flowCore::parameters(ff))
    signature <- gate_detector_signature(pd)
    cache <- getOption("spectreasy.gating_detector_cache", list())
    hit <- which(vapply(cache, function(entry) identical(entry$signature, signature), logical(1)))
    if (length(hit) > 0) return(cache[[hit[1]]]$info)
    if (is.null(compute_fun)) compute_fun <- spectreasy::get_sorted_detectors
    info <- tryCatch(compute_fun(pd), error = function(e) NULL)
    cache[[length(cache) + 1L]] <- list(signature = signature, info = info)
    options(spectreasy.gating_detector_cache = cache)
    info
}

gate_spectrum_counts_payload <- function(counts_mat) {
    connection <- rawConnection(raw(), open = "wb")
    on.exit(close(connection), add = TRUE)
    writeBin(as.integer(counts_mat), connection, size = 4L, endian = "little")
    list(
        format = jsonlite::unbox("uint32-column-major"),
        rows = jsonlite::unbox(nrow(counts_mat)),
        columns = jsonlite::unbox(ncol(counts_mat)),
        data = jsonlite::unbox(jsonlite::base64_enc(rawConnectionValue(connection)))
    )
}

gate_build_spectrum_data <- function(expr_filtered, spec_channels, det_info) {
    labels <- gate_detector_labels(det_info, spec_channels)
    y_power <- 1.5
    event_count <- nrow(expr_filtered)
    if (event_count == 0) {
        counts_mat <- matrix(0L, nrow = 150L, ncol = length(spec_channels))
        return(list(
            format = jsonlite::unbox("spectrum-histogram-v1"),
            channels = as.list(spec_channels),
            labels = as.list(labels),
            bin_mid = as.list(seq(0, 6, length.out = 150L)),
            bin_height = jsonlite::unbox(6 / 150),
            max_y = jsonlite::unbox(6),
            y_power = jsonlite::unbox(y_power),
            min_bin_count = jsonlite::unbox(1L),
            fill_limits = list(0, 1),
            event_count = jsonlite::unbox(0L),
            counts = gate_spectrum_counts_payload(counts_mat)
        ))
    }

    spectral_values <- if (identical(colnames(expr_filtered), spec_channels)) {
        as.matrix(expr_filtered)
    } else {
        as.matrix(expr_filtered[, spec_channels, drop = FALSE])
    }
    log_mat <- log10(pmax(spectral_values, 1e-3))
    finite <- log_mat[is.finite(log_mat)]
    if (length(finite) == 0) {
        return(gate_build_spectrum_data(expr_filtered[0, , drop = FALSE], spec_channels, det_info))
    }
    min_y <- floor(min(finite))
    max_y <- ceiling(max(finite))
    if (!is.finite(min_y)) min_y <- 0
    if (!is.finite(max_y)) max_y <- 6
    if (min_y == max_y) max_y <- min_y + 1
    breaks <- seq(min_y, max_y, length.out = 151L)
    bin_mid <- (breaks[-1] + breaks[-length(breaks)]) / 2
    bin_height <- breaks[2] - breaks[1]
    counts_mat <- vapply(seq_len(ncol(log_mat)), function(j) {
        values <- log_mat[, j]
        values <- values[is.finite(values)]
        as.integer(graphics::hist(values, breaks = breaks, plot = FALSE)$counts)
    }, integer(length(bin_mid)))
    if (is.null(dim(counts_mat))) counts_mat <- matrix(counts_mat, ncol = 1L)
    min_bin_count <- if (event_count <= 3000L) 1L else 3L
    visible_counts <- counts_mat[counts_mat >= min_bin_count]
    if (length(visible_counts) == 0) {
        fill_limits <- c(0, 1)
    } else {
        fills <- log10(visible_counts + 1)
        fill_limits <- c(min(fills), unname(stats::quantile(fills, 0.96, na.rm = TRUE)))
        if (!is.finite(fill_limits[2]) || fill_limits[2] <= fill_limits[1]) {
            fill_limits[2] <- fill_limits[1] + 1
        }
    }
    list(
        format = jsonlite::unbox("spectrum-histogram-v1"),
        channels = as.list(spec_channels),
        labels = as.list(labels),
        bin_mid = as.list(unname(bin_mid)),
        bin_height = jsonlite::unbox(unname(bin_height)),
        max_y = jsonlite::unbox(max(0, max_y)),
        y_power = jsonlite::unbox(y_power),
        min_bin_count = jsonlite::unbox(min_bin_count),
        fill_limits = as.list(unname(fill_limits)),
        event_count = jsonlite::unbox(event_count),
        counts = gate_spectrum_counts_payload(counts_mat)
    )
}

gate_spectrum_data_from_loaded <- function(expr, ff, spectrum_gates, peak, is_af = FALSE, scatter_mask = NULL) {
    det_info <- gate_detector_info(ff)
    if (is.null(det_info) || length(det_info$names) == 0) return(NULL)
    spec_channels <- intersect(det_info$names, colnames(expr))
    if (length(spec_channels) == 0) return(NULL)
    if (is.null(scatter_mask)) {
        scatter_mask <- gate_polygon_mask(expr, spectrum_gates$cell) &
            gate_polygon_mask(expr, spectrum_gates$singlet)
    }
    keep <- scatter_mask
    if (!isTRUE(is_af)) keep <- keep & gate_positive_mask(expr, spectrum_gates$positive, peak)
    spectral_expr <- expr[keep, spec_channels, drop = FALSE]
    gate_build_spectrum_data(spectral_expr, spec_channels, det_info)
}

gate_cache_spectrum <- function(filename, gates, spectrum) {
    cache <- getOption("spectreasy.gating_spectrum_cache", list())
    cache[[gate_safe_basename(filename)]] <- list(gates = gates, spectrum = spectrum)
    options(spectreasy.gating_spectrum_cache = cache)
    spectrum
}

gate_spectrum_for_file <- function(filename, gates = NULL) {
    df <- gate_read_mapping()
    name <- gate_safe_basename(filename)
    row <- df[df$filename == name, , drop = FALSE]
    if (nrow(row) == 0) stop("File is not present in mapping: ", name, call. = FALSE)
    if (is.null(gates)) {
        cache_data <- getOption("spectreasy.gating_state_cache")
        gates <- if (!is.null(cache_data)) cache_data$gates else NULL
    }
    control_type <- gate_normalize_control_type(row$control.type[1])
    spectrum_gates <- gate_spectrum_gates(
        gates,
        filename = name,
        control_type = control_type,
        is_af = isTRUE(row$is_af[1])
    )
    if (!isTRUE(spectrum_gates$ready)) return(NULL)

    spectrum_cache <- getOption("spectreasy.gating_spectrum_cache", list())
    gate_signature <- spectrum_gates[c("cell", "singlet", "positive")]
    cached <- spectrum_cache[[name]]
    if (!is.null(cached) && identical(cached$gates, gate_signature)) return(cached$spectrum)

    path <- file.path(get_gate_scc_dir(), name)
    ff <- gate_read_fcs(path)
    expr <- flowCore::exprs(ff)
    peak <- as.character(row$channel[1])
    if (!peak %in% colnames(expr)) {
        det_info <- gate_detector_info(ff)
        spec_channels <- if (is.null(det_info)) character() else intersect(det_info$names, colnames(expr))
        if (length(spec_channels) == 0) return(NULL)
        peak <- spec_channels[which.max(apply(expr[, spec_channels, drop = FALSE], 2, stats::var, na.rm = TRUE))]
    }
    spectrum <- gate_spectrum_data_from_loaded(
        expr = expr,
        ff = ff,
        spectrum_gates = spectrum_gates,
        peak = peak,
        is_af = isTRUE(row$is_af[1])
    )
    gate_cache_spectrum(name, gate_signature, spectrum)
}

gate_payload_for_file <- function(filename, max_points = 3000L, cache_result = TRUE) {
    name <- gate_safe_basename(filename)
    cache_key <- paste(name, max_points, sep = "::")
    cache <- getOption("spectreasy.gating_payload_cache", list())
    if (isTRUE(cache_result) && !is.null(cache[[cache_key]])) return(cache[[cache_key]])

    df <- gate_read_mapping()
    row <- df[df$filename == name, , drop = FALSE]
    if (nrow(row) == 0) stop("File is not present in mapping: ", name, call. = FALSE)

    path <- file.path(get_gate_scc_dir(), name)
    ff <- gate_read_fcs(path)
    expr <- flowCore::exprs(ff)
    meta <- gate_fcs_metadata(path)
    scatter_channels <- gate_scatter_channels(colnames(expr))
    peak <- as.character(row$channel[1])
    if (!peak %in% colnames(expr)) {
        det_info <- gate_detector_info(ff)
        det <- if (is.null(det_info)) colnames(expr) else det_info$names
        det <- intersect(det, colnames(expr))
        peak <- det[which.max(apply(expr[, det, drop = FALSE], 2, stats::var, na.rm = TRUE))]
    }

    n <- nrow(expr)
    max_points <- suppressWarnings(as.integer(max_points[1]))
    use_all <- is.finite(max_points) && !is.na(max_points) && max_points <= 0
    if (!use_all && (!is.finite(max_points) || is.na(max_points) || max_points < 500)) max_points <- 3000L
    if (!use_all && n > max_points) {
        seed <- sum(utf8ToInt(name)) %% .Machine$integer.max
        set.seed(seed)
        idx <- sort(sample.int(n, max_points))
    } else {
        idx <- seq_len(n)
    }
    out <- data.frame(
        event_index = idx,
        fsc_a = as.numeric(expr[idx, meta$fsc_a]),
        ssc_a = as.numeric(expr[idx, meta$ssc_a]),
        fsc_h = as.numeric(expr[idx, meta$fsc_h]),
        peak = as.numeric(expr[idx, peak]),
        stringsAsFactors = FALSE
    )
    for (channel in setdiff(scatter_channels, colnames(out))) {
        out[[channel]] <- as.numeric(expr[idx, channel])
    }
    out <- out[stats::complete.cases(out), , drop = FALSE]
    scatter_domains <- stats::setNames(lapply(scatter_channels, function(channel) gate_compute_extent(out[[channel]])), scatter_channels)
    payload <- list(
        file = row[1, , drop = FALSE],
        total_events = n,
        returned_events = nrow(out),
        channels = list(fsc_a = meta$fsc_a, fsc_h = meta$fsc_h, ssc_a = meta$ssc_a, peak = peak, scatter = scatter_channels),
        labels = meta$labels,
        domains = list(
            fsc_a = gate_compute_extent(out$fsc_a),
            ssc_a = gate_compute_extent(out$ssc_a),
            fsc_h = gate_compute_extent(out$fsc_h),
            peak = gate_compute_extent(out$peak),
            scatter = scatter_domains
        ),
        events = out
    )
    if (isTRUE(cache_result)) {
        cache[[cache_key]] <- payload
        options(spectreasy.gating_payload_cache = cache)
    }
    payload
}

gate_compact_payload <- function(payload) {
    events <- payload$events
    matrix_values <- as.matrix(events)
    storage.mode(matrix_values) <- "double"
    connection <- rawConnection(raw(), open = "wb")
    on.exit(close(connection), add = TRUE)
    writeBin(as.numeric(matrix_values), connection, size = 4L, endian = "little")
    encoded <- jsonlite::base64_enc(rawConnectionValue(connection))
    payload$events <- NULL
    payload$events_compact <- list(
        format = "float32-column-major",
        fields = colnames(matrix_values),
        rows = nrow(matrix_values),
        data = encoded
    )
    payload
}
