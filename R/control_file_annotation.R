.control_file_row_from_scc <- function(fn, ref, custom_fluorophores = NULL) {
    fluor <- .infer_fluor_from_filename(fn, ref = ref, custom_fluorophores = custom_fluorophores)
    is_bead_negative <- .control_file_is_bead_negative_filename(fn)
    if (is_bead_negative) {
        fluor <- "AF_beads"
    } else if (.is_dead_af_filename(fn)) {
        fluor <- "AF_dead"
    } else if (grepl("Unstained|US_UT", fn, ignore.case = TRUE)) {
        fluor <- "AF"
    }
    marker <- .infer_marker_from_filename(
        fn,
        fluor_guess = fluor,
        marker_name_map = ref$marker_name_map,
        marker_names = ref$marker_names
    )
    if (is_bead_negative) {
        marker <- "Bead background"
    } else if (.is_dead_af_filename(fn)) {
        marker <- "Dead cell background"
    }
    control_type <- .infer_control_type_from_filename(fn, fluor_guess = fluor, marker_guess = marker)
    viability_flag <- .infer_is_viability_from_filename(fn, fluor_guess = fluor, marker_guess = marker)

    data.frame(
        filename = fn,
        fluorophore = fluor,
        marker = marker,
        channel = "",
        control.type = control_type,
        universal.negative = "",
        is.viability = viability_flag,
        stringsAsFactors = FALSE
    )
}

.build_control_file_scc_df <- function(scc_files, ref, custom_fluorophores = NULL) {
    if (length(scc_files) == 0) {
        return(.empty_control_file_df())
    }
    rows <- lapply(scc_files, .control_file_row_from_scc, ref = ref, custom_fluorophores = custom_fluorophores)
    do.call(rbind, rows)
}

.control_file_row_path <- function(fn, input_folder) {
    file.path(input_folder, fn)
}

.summarize_control_file_fcs <- function(path, channel_alias_map = character()) {
    pd <- NULL
    peak_channel <- ""
    top_channels <- character()
    top_ratio <- NA_real_
    entropy_hi <- NA_real_
    read_error <- ""

    tryCatch({
        ff <- suppressWarnings(.spectreasy_read_fcs(path))
        pd <- flowCore::pData(flowCore::parameters(ff))
        channel_alias_map <- .merge_control_file_alias_map(channel_alias_map, .build_channel_alias_map_from_pd(pd))
        fl_pd <- get_sorted_detectors(pd)
        if (length(fl_pd$names) > 0) {
            expr <- flowCore::exprs(ff)[, fl_pd$names, drop = FALSE]

            use_idx <- seq_len(nrow(expr))
            if ("name" %in% colnames(pd)) {
                fsc <- pd$name[grepl("^FSC", pd$name) & grepl("-A$", pd$name)][1]
                ssc <- pd$name[grepl("^SSC", pd$name) & grepl("-A$", pd$name)][1]
                if (!is.na(fsc) && !is.na(ssc) && fsc %in% colnames(flowCore::exprs(ff)) && ssc %in% colnames(flowCore::exprs(ff))) {
                    scatter <- flowCore::exprs(ff)[, c(fsc, ssc), drop = FALSE]
                    idx <- which(
                        scatter[, 1] > quantile(scatter[, 1], 0.10, na.rm = TRUE) &
                        scatter[, 1] < quantile(scatter[, 1], 0.90, na.rm = TRUE) &
                        scatter[, 2] > quantile(scatter[, 2], 0.10, na.rm = TRUE) &
                        scatter[, 2] < quantile(scatter[, 2], 0.90, na.rm = TRUE)
                    )
                    if (length(idx) >= 200) use_idx <- idx
                }
            }
            expr_use <- expr[use_idx, , drop = FALSE]
            q_hi <- apply(expr_use, 2, function(v) quantile(v, 0.999, na.rm = TRUE))
            q_hi[!is.finite(q_hi)] <- NA_real_
            peak_channel <- if (all(is.na(q_hi))) {
                medians <- apply(expr_use, 2, median, na.rm = TRUE)
                names(medians)[which.max(medians)]
            } else {
                names(q_hi)[which.max(q_hi)]
            }

            q_pos <- q_hi[is.finite(q_hi) & q_hi > 0]
            if (length(q_pos) > 0) {
                ord <- order(q_pos, decreasing = TRUE)
                top_channels <- names(q_pos)[ord[seq_len(min(3, length(ord)))]]
                p <- q_pos / sum(q_pos)
                entropy_hi <- -sum(p * log(p))
                if (length(ord) >= 2) {
                    top_ratio <- as.numeric(q_pos[ord[1]] / pmax(q_pos[ord[2]], 1e-9))
                }
            }
        }
    }, error = function(e) {
        read_error <<- conditionMessage(e)
    })

    list(
        pd = pd,
        peak_channel = peak_channel,
        top_channels = top_channels,
        top_ratio = top_ratio,
        entropy_hi = entropy_hi,
        read_error = read_error,
        channel_alias_map = channel_alias_map
    )
}

.control_file_should_warn_read_error <- function(path) {
    size <- suppressWarnings(file.info(path)$size)
    is.finite(size) && !is.na(size) && size > 0
}

.control_file_is_af_like <- function(current_fluor,
                                     current_marker,
                                     top_channels,
                                     top_ratio,
                                     entropy_hi,
                                     expected_af_channel,
                                     channel_alias_map = character()) {
    unresolved_label <- (is.na(current_fluor) || current_fluor == "" || identical(current_fluor, "Unknown")) &&
        (is.na(current_marker) || current_marker == "")
    if (!unresolved_label || !nzchar(expected_af_channel) || length(top_channels) == 0) {
        return(FALSE)
    }

    resolved_top <- unique(vapply(top_channels, .control_file_resolve_channel_alias, character(1), channel_alias_map = channel_alias_map))
    af_in_top <- expected_af_channel %in% resolved_top
    broad_signal <- is.finite(top_ratio) && top_ratio <= 1.7
    diffuse_signal <- is.finite(entropy_hi) && entropy_hi >= 3.2
    af_in_top && broad_signal && diffuse_signal
}

.annotate_control_file_rows <- function(df,
                                        input_folder,
                                        cytometer,
                                        unknown_fluor_policy,
                                        ref) {
    channel_alias_map <- character()
    expected_af_channel <- .control_file_resolve_channel_alias(
        .control_file_get_expected_af_channel(cytometer),
        channel_alias_map = channel_alias_map
    )

    for (i in seq_len(nrow(df))) {
        fn <- df$filename[i]
        path <- .control_file_row_path(fn, input_folder = input_folder)
        summary <- .summarize_control_file_fcs(path, channel_alias_map = channel_alias_map)
        channel_alias_map <- summary$channel_alias_map

        if (nzchar(summary$read_error) && .control_file_should_warn_read_error(path)) {
            warning(
                "Could not read FCS file while auto-detecting peak channel: ",
                path,
                "\n",
                summary$read_error,
                call. = FALSE
            )
        }

        current_fluor <- trimws(as.character(df$fluorophore[i]))
        expected_channel <- .control_file_expected_channel_for_fluor(
            current_fluor,
            fluor_peak_channel_map = ref$fluor_peak_channel_map
        )
        if (nzchar(expected_channel)) {
            df$channel[i] <- expected_channel
        } else if (nzchar(summary$peak_channel)) {
            df$channel[i] <- summary$peak_channel
        }

        current_marker <- trimws(as.character(df$marker[i]))
        af_like <- .control_file_is_af_like(
            current_fluor = current_fluor,
            current_marker = current_marker,
            top_channels = summary$top_channels,
            top_ratio = summary$top_ratio,
            entropy_hi = summary$entropy_hi,
            expected_af_channel = expected_af_channel,
            channel_alias_map = channel_alias_map
        )
        if (af_like) {
            new_tag <- .control_file_next_af_tag(df$fluorophore)
            df$fluorophore[i] <- new_tag
            df$marker[i] <- "Autofluorescence"
            df$control.type[i] <- "cells"
        }

        current_fluor <- trimws(as.character(df$fluorophore[i]))
        unresolved <- is.na(current_fluor) || current_fluor == "" || identical(current_fluor, "Unknown")
        if (unresolved) {
            if (unknown_fluor_policy == "filename") {
                df$fluorophore[i] <- tools::file_path_sans_ext(fn)
            } else if (unknown_fluor_policy == "by_channel") {
                detector_desc <- ""
                if (!is.null(summary$pd) && "desc" %in% colnames(summary$pd) && "name" %in% colnames(summary$pd)) {
                    idx <- match(df$channel[i], as.character(summary$pd$name))
                    if (!is.na(idx)) {
                        detector_desc <- as.character(summary$pd$desc[idx])
                    }
                }
                df$fluorophore[i] <- .infer_fluor_from_detector(
                    df$channel[i],
                    detector_desc = detector_desc,
                    channel_alias_map = channel_alias_map,
                    fluor_channel_map = ref$fluor_channel_map
                )
            } else {
                df$fluorophore[i] <- ""
            }
        }

        current_marker <- trimws(as.character(df$marker[i]))
        if (is.na(current_marker) || current_marker == "") {
            df$marker[i] <- .infer_marker_from_filename(
                fn,
                fluor_guess = df$fluorophore[i],
                marker_name_map = ref$marker_name_map,
                marker_names = ref$marker_names
            )
        }

        current_type <- trimws(as.character(df$control.type[i]))
        if (is.na(current_type) || current_type == "") {
            df$control.type[i] <- .infer_control_type_from_filename(
                fn,
                fluor_guess = df$fluorophore[i],
                marker_guess = df$marker[i]
            )
        }

        current_viability <- toupper(trimws(as.character(df$is.viability[i])))
        if (is.na(current_viability) || current_viability == "") {
            df$is.viability[i] <- .infer_is_viability_from_filename(
                fn,
                fluor_guess = df$fluorophore[i],
                marker_guess = df$marker[i]
            )
        }
    }

    df
}

.finalize_control_file_df <- function(df, scc_files) {
    if (!("universal.negative" %in% colnames(df))) {
        df$universal.negative <- ""
    }
    bead_negative_files <- scc_files[vapply(scc_files, .control_file_is_bead_negative_filename, logical(1))]

    df$control.type <- tolower(trimws(as.character(df$control.type)))
    df$control.type[is.na(df$control.type)] <- ""
    invalid_type <- !(df$control.type %in% c("", "beads", "cells"))
    df$control.type[invalid_type] <- ""

    df$is.viability <- toupper(trimws(as.character(df$is.viability)))
    df$is.viability[is.na(df$is.viability)] <- ""
    df$is.viability[df$is.viability %in% c("T", "TRUE", "1", "YES", "Y")] <- "TRUE"
    df$is.viability[!(df$is.viability %in% c("", "TRUE"))] <- ""
    viability_missing_marker <- df$is.viability == "TRUE" & !nzchar(trimws(as.character(df$marker)))
    df$marker[viability_missing_marker] <- "Live"

    bead_negative_row <- df$filename %in% bead_negative_files |
        (tolower(df$control.type) == "beads" & grepl("^AF_beads?$", as.character(df$fluorophore), ignore.case = TRUE))
    if (any(bead_negative_row)) {
        df$fluorophore[bead_negative_row] <- "AF_beads"
        df$marker[bead_negative_row] <- "Bead background"
        df$control.type[bead_negative_row] <- "beads"
        df$universal.negative[bead_negative_row] <- ""
        df$is.viability[bead_negative_row] <- ""
    }

    dead_af_row <- .is_dead_af_control_row(
        fluorophore = df$fluorophore,
        marker = if ("marker" %in% colnames(df)) df$marker else NULL,
        filename = df$filename,
        control_type = df$control.type
    )
    if (any(dead_af_row)) {
        df$fluorophore[dead_af_row] <- "AF_dead"
        df$marker[dead_af_row] <- "Dead cell background"
        df$control.type[dead_af_row] <- "cells"
        df$is.viability[dead_af_row] <- ""
        df$universal.negative[dead_af_row] <- ""
    }

    is_af_row <- grepl("^AF($|_)", as.character(df$fluorophore), ignore.case = TRUE) & !bead_negative_row
    df$control.type[is_af_row] <- "cells"

    df <- .canonicalize_primary_af_labels(df)
    df <- .assign_reference_automatic_negatives(df)

    preferred_order <- c(
        "filename",
        "fluorophore",
        "marker",
        "channel",
        "control.type",
        "universal.negative",
        "is.viability"
    )
    keep <- intersect(preferred_order, colnames(df))
    df[, c(keep, setdiff(colnames(df), keep)), drop = FALSE]
}
