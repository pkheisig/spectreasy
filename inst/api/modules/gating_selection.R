gate_compute_extent <- function(v) {
    finite <- v[is.finite(v)]
    if (length(finite) == 0) return(c(0, 1))
    val_max <- as.numeric(stats::quantile(finite, 0.998, na.rm = TRUE))
    if (!is.finite(val_max) || val_max <= 0) val_max <- max(finite, na.rm = TRUE)
    if (!is.finite(val_max) || val_max <= 0) val_max <- 1
    c(0, val_max + val_max * 0.04)
}

gate_scatter_channels <- function(col_names) {
    spectreasy:::.get_scatter_channels(col_names)
}

gate_point_in_polygon <- function(x, y, poly) {
    n <- nrow(poly)
    inside <- rep(FALSE, length(x))
    if (n < 3) return(inside)
    j <- n
    for (i in seq_len(n)) {
        xi <- poly[i, 1]
        yi <- poly[i, 2]
        xj <- poly[j, 1]
        yj <- poly[j, 2]
        intersect <- ((yi > y) != (yj > y)) & (x < (xj - xi) * (y - yi) / (yj - yi + 1e-12) + xi)
        inside[intersect] <- !inside[intersect]
        j <- i
    }
    inside
}

gate_vertices_matrix <- function(gate) {
    if (is.null(gate$vertices) || length(gate$vertices) == 0) return(NULL)
    verts <- do.call(rbind, lapply(gate$vertices, function(v) c(as.numeric(v$x), as.numeric(v$y))))
    if (!is.matrix(verts) || nrow(verts) == 0 || ncol(verts) < 2) return(NULL)
    verts[stats::complete.cases(verts), , drop = FALSE]
}

gate_value <- function(gate, primary, fallback = NULL) {
    value <- gate[[primary]]
    if (!is.null(value)) return(value)
    alt <- switch(primary,
        xChannel = "x_channel",
        yChannel = "y_channel",
        NULL
    )
    if (!is.null(alt) && !is.null(gate[[alt]])) return(gate[[alt]])
    fallback
}

gate_cached_gate <- function(cache, gate_type, filename, control_type) {
    if (is.null(cache)) return(NULL)
    file_key <- paste0(gate_type, ":", filename)
    global_key <- paste0(gate_type, ":", control_type)
    if (file_key %in% names(cache)) return(cache[[file_key]])
    if (global_key %in% names(cache)) return(cache[[global_key]])
    NULL
}

gate_apply_polygon_gate <- function(expr, gate) {
    expr[gate_polygon_mask(expr, gate), , drop = FALSE]
}

gate_polygon_mask <- function(expr, gate) {
    keep_all <- rep(TRUE, nrow(expr))
    verts <- gate_vertices_matrix(gate)
    if (is.null(verts) || nrow(verts) < 3) return(keep_all)
    x_channel <- as.character(gate_value(gate, "xChannel", ""))[1]
    y_channel <- as.character(gate_value(gate, "yChannel", ""))[1]
    if (!x_channel %in% colnames(expr) || !y_channel %in% colnames(expr)) return(keep_all)
    keep <- gate_point_in_polygon(expr[, x_channel], expr[, y_channel], verts)
    !is.na(keep) & keep
}

gate_positive_mask <- function(expr, gate, peak) {
    keep_all <- rep(TRUE, nrow(expr))
    if (is.null(gate) || !peak %in% colnames(expr)) return(keep_all)
    verts <- gate_vertices_matrix(gate)
    if (is.null(verts) || nrow(verts) == 0) return(keep_all)
    mode <- as.character(gate_value(gate, "mode", ""))[1]
    keep <- NULL
    if (identical(mode, "separator")) {
        keep <- expr[, peak] >= verts[1, 1]
    } else if (identical(mode, "positive_1d") && nrow(verts) >= 2) {
        lim <- range(verts[, 1], na.rm = TRUE)
        keep <- expr[, peak] >= lim[1] & expr[, peak] <= lim[2]
    } else if (nrow(verts) >= 3) {
        x_channel <- as.character(gate_value(gate, "xChannel", peak))[1]
        y_channel <- as.character(gate_value(gate, "yChannel", ""))[1]
        if (!x_channel %in% colnames(expr)) x_channel <- peak
        if (!y_channel %in% colnames(expr)) return(keep_all)
        keep <- gate_point_in_polygon(expr[, x_channel], expr[, y_channel], verts)
    }
    if (is.null(keep)) return(keep_all)
    !is.na(keep) & keep
}

gate_apply_positive_gate <- function(expr, gate, peak) {
    expr[gate_positive_mask(expr, gate, peak), , drop = FALSE]
}

gate_is_finalized_polygon <- function(gate) {
    if (is.null(gate) || identical(as.character(gate_value(gate, "mode", ""))[1], "blocked")) {
        return(FALSE)
    }
    verts <- gate_vertices_matrix(gate)
    !is.null(verts) && nrow(verts) >= 3 && all(is.finite(verts))
}

gate_is_finalized_histogram <- function(gate) {
    if (is.null(gate)) return(FALSE)
    mode <- as.character(gate_value(gate, "mode", ""))[1]
    if (!nzchar(mode) || mode %in% c("missing", "blocked")) return(FALSE)
    verts <- gate_vertices_matrix(gate)
    if (is.null(verts) || !all(is.finite(verts))) return(FALSE)
    if (identical(mode, "separator")) return(nrow(verts) >= 1L)
    if (mode %in% c("positive_1d", "negative_1d")) return(nrow(verts) >= 2L)
    nrow(verts) >= 3L
}

gate_spectrum_gates <- function(cache, filename, control_type, is_af = FALSE) {
    gates <- list(
        cell = gate_cached_gate(cache, "cell", filename, control_type),
        singlet = gate_cached_gate(cache, "singlet", filename, control_type),
        positive = gate_cached_gate(cache, "positive", filename, control_type)
    )
    gates$ready <- gate_is_finalized_polygon(gates$cell) &&
        gate_is_finalized_polygon(gates$singlet) &&
        (isTRUE(is_af) || gate_is_finalized_histogram(gates$positive))
    gates
}

gate_histogram_autogate_ranges <- function(peak_vals,
                                           sample_type,
                                           is_viability = FALSE,
                                           compute_fun = NULL,
                                           include_negative = TRUE) {
    peak_vals <- as.numeric(peak_vals)
    peak_vals <- peak_vals[is.finite(peak_vals)]
    if (length(peak_vals) < 10L) {
        stop("Too few gated events to auto-generate histogram gates.", call. = FALSE)
    }
    if (is.null(compute_fun)) {
        compute_fun <- get(".compute_reference_histogram_gate", envir = asNamespace("spectreasy"))
    }
    opts <- spectreasy::gating_options()
    result <- compute_fun(
        peak_vals = peak_vals,
        sample_type = sample_type,
        histogram_pct_beads = opts$histogram_pct_beads,
        histogram_direction_beads = opts$histogram_direction_beads,
        histogram_pct_cells = opts$histogram_pct_cells,
        histogram_direction_cells = opts$histogram_direction_cells,
        is_viability = isTRUE(is_viability)
    )
    positive <- sort(as.numeric(c(result$gate_min, result$gate_max)))
    negative <- sort(10^as.numeric(c(
        attr(result$vals_log, "neg_log_min"),
        attr(result$vals_log, "neg_log_max")
    )))
    if (length(positive) != 2L || any(!is.finite(positive)) || positive[2] <= positive[1]) {
        stop("The histogram autogater did not produce a valid positive interval.", call. = FALSE)
    }
    if (isTRUE(include_negative) && (!isTRUE(attr(result$vals_log, "negative_gate_present")) ||
        length(negative) != 2L || any(!is.finite(negative)) || negative[2] <= negative[1])) {
        stop("The histogram autogater did not produce a valid negative interval.", call. = FALSE)
    }
    repair_empty_interval <- function(limits, direction) {
        in_gate <- peak_vals >= limits[1] & peak_vals <= limits[2]
        minimum_count <- min(10L, length(peak_vals))
        if (sum(in_gate, na.rm = TRUE) >= minimum_count) return(limits)

        ordered <- sort(peak_vals)
        fallback_count <- min(
            length(ordered),
            max(minimum_count, as.integer(ceiling(length(ordered) * 0.05)))
        )
        fallback <- if (identical(direction, "high")) {
            utils::tail(ordered, fallback_count)
        } else {
            utils::head(ordered, fallback_count)
        }
        repaired <- range(fallback)
        if (repaired[1] == repaired[2]) {
            padding <- max(abs(repaired[1]) * 1e-8, .Machine$double.eps^0.5)
            repaired <- repaired + c(-padding, padding)
        }
        repaired
    }
    positive <- repair_empty_interval(positive, "high")
    if (isTRUE(include_negative)) {
        negative <- repair_empty_interval(negative, "low")
    } else {
        negative <- NULL
    }
    list(positive = positive, negative = negative)
}

gate_histogram_interval <- function(type, filename, peak, limits) {
    list(
        type = jsonlite::unbox(as.character(type)[1]),
        scope = jsonlite::unbox("file"),
        filename = jsonlite::unbox(as.character(filename)[1]),
        xChannel = jsonlite::unbox(as.character(peak)[1]),
        yChannel = jsonlite::unbox(""),
        mode = jsonlite::unbox(paste0(as.character(type)[1], "_1d")),
        vertices = list(
            list(x = jsonlite::unbox(unname(as.numeric(limits[1]))), y = jsonlite::unbox(0)),
            list(x = jsonlite::unbox(unname(as.numeric(limits[2]))), y = jsonlite::unbox(0))
        )
    )
}

gate_autogenerate_histograms <- function(gates) {
    df <- gate_read_mapping()
    uses_histogram <- if ("uses_histogram_gates" %in% colnames(df)) {
        as.logical(df$uses_histogram_gates)
    } else {
        !as.logical(df$is_af)
    }
    uses_histogram[is.na(uses_histogram)] <- FALSE
    targets <- df[uses_histogram, , drop = FALSE]
    if (nrow(targets) == 0) {
        stop("No non-AF single-color controls are available for histogram autogating.", call. = FALSE)
    }

    resolved <- lapply(seq_len(nrow(targets)), function(i) {
        row <- targets[i, , drop = FALSE]
        filename <- as.character(row$filename[1])
        control_type <- gate_normalize_control_type(row$control.type[1])
        cell_gate <- gate_cached_gate(gates, "cell", filename, control_type)
        singlet_gate <- gate_cached_gate(gates, "singlet", filename, control_type)
        list(
            row = row,
            filename = filename,
            control_type = control_type,
            cell_gate = cell_gate,
            singlet_gate = singlet_gate,
            needs_negative = isTRUE(as.logical(row$uses_negative_histogram_gate[1])),
            ready = gate_is_finalized_polygon(cell_gate) && gate_is_finalized_polygon(singlet_gate)
        )
    })
    missing <- vapply(resolved, function(item) if (item$ready) "" else item$filename, character(1))
    missing <- missing[nzchar(missing)]
    if (length(missing) > 0) {
        stop(
            "Complete FSC/SSC and singlet gates before histogram autogating: ",
            paste(missing, collapse = ", "),
            call. = FALSE
        )
    }

    generated <- list()
    spectra <- list()
    preserved_count <- 0L
    for (item in resolved) {
        positive_key <- paste0("positive:", item$filename)
        negative_key <- paste0("negative:", item$filename)
        keep_positive <- gate_is_finalized_histogram(gates[[positive_key]])
        keep_negative <- item$needs_negative && gate_is_finalized_histogram(gates[[negative_key]])
        preserved_count <- preserved_count + as.integer(keep_positive) + as.integer(keep_negative)

        path <- file.path(get_gate_scc_dir(), item$filename)
        ff <- gate_read_fcs(path)
        expr <- flowCore::exprs(ff)
        peak <- as.character(item$row$channel[1])
        if (!peak %in% colnames(expr)) {
            det_info <- gate_detector_info(ff)
            detectors <- if (is.null(det_info)) colnames(expr) else det_info$names
            detectors <- intersect(detectors, colnames(expr))
            if (length(detectors) == 0) {
                stop("Could not resolve a peak channel for ", item$filename, ".", call. = FALSE)
            }
            peak <- detectors[which.max(apply(expr[, detectors, drop = FALSE], 2, stats::var, na.rm = TRUE))]
        }
        scatter_mask <- gate_polygon_mask(expr, item$cell_gate) &
            gate_polygon_mask(expr, item$singlet_gate)
        if (!keep_positive || (item$needs_negative && !keep_negative)) {
            ranges <- gate_histogram_autogate_ranges(
                expr[scatter_mask, peak],
                sample_type = item$control_type,
                is_viability = isTRUE(item$row$is_viability[1]),
                include_negative = item$needs_negative
            )
            if (!keep_positive) {
                generated[[positive_key]] <- gate_histogram_interval(
                    "positive", item$filename, peak, ranges$positive
                )
            }
            if (item$needs_negative && !keep_negative) {
                generated[[negative_key]] <- gate_histogram_interval(
                    "negative", item$filename, peak, ranges$negative
                )
            }
        }
        positive_gate <- if (keep_positive) gates[[positive_key]] else generated[[positive_key]]
        spectrum_gates <- list(
            cell = item$cell_gate,
            singlet = item$singlet_gate,
            positive = positive_gate
        )
        spectrum <- gate_spectrum_data_from_loaded(
            expr = expr,
            ff = ff,
            spectrum_gates = spectrum_gates,
            peak = peak,
            is_af = FALSE,
            scatter_mask = scatter_mask
        )
        spectra[[item$filename]] <- gate_cache_spectrum(
            item$filename,
            spectrum_gates,
            spectrum
        )
    }
    list(
        gates = generated,
        spectra = spectra,
        files_processed = nrow(targets),
        gates_generated = length(generated),
        gates_preserved = preserved_count
    )
}
