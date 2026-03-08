#' Build a Reference Matrix from Single-Color Controls
#'
#' Reads SCC FCS files, performs FSC/SSC gating and positive-peak histogram gating,
#' then computes one normalized spectrum per fluorophore.
#'
#' This is the core matrix-construction step used before unmixing.
#'
#' @param input_folder Directory containing SCC `.fcs` files.
#' @param output_folder Directory where gating/spectrum plots are written.
#' @param save_qc_plots Logical; if `TRUE`, write FSC/SSC, histogram, and spectrum plots.
#' @param control_df Optional control mapping as a data.frame or CSV path.
#'   Expected columns: `filename`, `fluorophore`, `channel`; `universal.negative` is optional.
#' @param include_multi_af Logical; if `TRUE`, include additional AF controls from `af_dir`.
#' @param af_dir Directory with extra AF controls when `include_multi_af = TRUE`.
#' @param default_sample_type Fallback type when filename heuristics are ambiguous (`"beads"` or `"cells"`).
#' @param cytometer Cytometer name used for channel alias resolution (for example `"Aurora"`).
#' @param histogram_pct_beads Histogram gate width for bead controls.
#' @param histogram_direction_beads Histogram gate direction for bead controls: `"both"`, `"left"`, or `"right"`.
#' @param histogram_pct_cells Histogram gate width for cell controls.
#' @param histogram_direction_cells Histogram gate direction for cell controls: `"both"`, `"left"`, or `"right"`.
#' @param outlier_percentile Upper-tail FSC/SSC filtering percentile.
#' @param debris_percentile Debris filtering percentile for cell controls.
#' @param bead_gate_scale Ellipse scaling factor for bead FSC/SSC gate.
#' @param histogram_min_x_log Reserved histogram lower x-limit parameter.
#' @param max_clusters Maximum number of GMM components tested.
#' @param min_cluster_proportion Minimum population proportion kept from GMM fit.
#' @param gate_contour_beads Contour probability for bead gating ellipse/hull.
#' @param gate_contour_cells Contour probability for cell gating ellipse/hull.
#' @param subsample_n Maximum number of events used for GMM fitting per file.
#'
#' @return Numeric matrix with rows = fluorophores and columns = detectors (normalized spectra).
#' @export
#' @examples
#' \dontrun{
#' # Build matrix from SCC folder and control CSV
#' M <- build_reference_matrix(
#'   input_folder = "scc",
#'   output_folder = "spectreasy_outputs/build_reference_plots",
#'   control_df = "fcs_control_file.csv",
#'   cytometer = "Aurora"
#' )
#'
#' # Build quickly without saving QC plots
#' M_fast <- build_reference_matrix(
#'   input_folder = "scc",
#'   save_qc_plots = FALSE
#' )
#' }
build_reference_matrix <- function(
  input_folder = "scc",
  output_folder = "gating_and_spectrum_plots",
  save_qc_plots = TRUE,
  control_df = NULL,
  include_multi_af = FALSE,
  af_dir = "af",
  default_sample_type = "beads",
  cytometer = "Aurora",
  histogram_pct_beads = 0.98,
  histogram_direction_beads = "both",
  histogram_pct_cells = 0.35,
  histogram_direction_cells = "both",
  outlier_percentile = 0.02,
  debris_percentile = 0.02,
  bead_gate_scale = 1.3,
  histogram_min_x_log = 2,
  max_clusters = 6,
  min_cluster_proportion = 0.03,
  gate_contour_beads = 0.999999999,
  gate_contour_cells = 0.95,
  subsample_n = 5000
) {
    if (is.character(control_df) && length(control_df) == 1 && !is.na(control_df)) {
        if (!file.exists(control_df)) stop("control_df file not found: ", control_df)
        control_df <- utils::read.csv(control_df, stringsAsFactors = FALSE, check.names = FALSE)
    } else if (!is.null(control_df) && !is.data.frame(control_df)) {
        stop("control_df must be either a data.frame or a single CSV path.")
    }
    if (is.data.frame(control_df)) {
        if (!("filename" %in% colnames(control_df))) {
            stop("control_df is missing required column: filename")
        }
        if (!("fluorophore" %in% colnames(control_df))) control_df$fluorophore <- ""
        if (!("channel" %in% colnames(control_df))) control_df$channel <- ""
        if (!("universal.negative" %in% colnames(control_df))) control_df$universal.negative <- ""

        control_df$filename <- trimws(as.character(control_df$filename))
        control_df$fluorophore <- trimws(as.character(control_df$fluorophore))
        control_df$channel <- trimws(as.character(control_df$channel))
        control_df$universal.negative <- trimws(as.character(control_df$universal.negative))
    }

    sample_patterns <- get_fluorophore_patterns()
    fcs_files <- list.files(input_folder, pattern = "\\.fcs$", full.names = TRUE)
    if (length(fcs_files) == 0) stop("No FCS files found in ", input_folder)

    get_control_rows <- function(df, filenames) {
        if (is.null(df) || !("filename" %in% colnames(df))) {
            return(data.frame())
        }
        fn <- as.character(df$filename)
        df[fn %in% filenames, ]
    }

    # Filter to prioritize _aligned versions if duplicates exist
    basenames <- tools::file_path_sans_ext(basename(fcs_files))
    aligned_indices <- grep("_aligned$", basenames)

    if (length(aligned_indices) > 0) {
        kept_files <- fcs_files
        for (idx in aligned_indices) {
            aligned_file <- fcs_files[idx]
            base_original <- sub("_aligned$", "", basenames[idx])
            # pattern match exact basename (e.g. "File.fcs" vs "File_aligned.fcs")
            # We want to remove "File.fcs" if "File_aligned.fcs" is present
            dup_idx <- which(basenames == base_original)
            if (length(dup_idx) > 0) {
                message("  Using aligned file: ", basename(aligned_file), " (Skipping original)")
                kept_files <- setdiff(kept_files, fcs_files[dup_idx])
            }
        }
        fcs_files <- kept_files
    }
    out_path <- normalizePath(output_folder, mustWork = FALSE)
    if (save_qc_plots) {
        dir.create(file.path(out_path, "fsc_ssc"), showWarnings = FALSE, recursive = TRUE)
        dir.create(file.path(out_path, "histogram"), showWarnings = FALSE, recursive = TRUE)
        dir.create(file.path(out_path, "spectrum"), showWarnings = FALSE, recursive = TRUE)
    }

    # 1. Establish Detector Sorting and Labels
    ff_meta <- flowCore::read.FCS(fcs_files[1], transformation = FALSE, truncate_max_range = FALSE)
    det_info <- get_sorted_detectors(flowCore::pData(flowCore::parameters(ff_meta)))
    detector_names <- det_info$names
    detector_labels <- det_info$labels

    normalize_channel <- function(x) {
        out <- toupper(gsub("\\s+", "", trimws(as.character(x))))
        out[is.na(out)] <- ""
        out
    }

    load_cytometer_alias_map <- function(cytometer_name) {
        if (!requireNamespace("AutoSpectral", quietly = TRUE)) return(character())
        cdb <- system.file("extdata/cytometer_database.csv", package = "AutoSpectral")
        if (!file.exists(cdb)) return(character())
        cyto_db <- tryCatch(
            utils::read.csv(cdb, stringsAsFactors = FALSE, check.names = FALSE),
            error = function(e) NULL
        )
        if (is.null(cyto_db) || nrow(cyto_db) == 0) return(character())
        target_col <- colnames(cyto_db)[tolower(colnames(cyto_db)) == tolower(cytometer_name)]
        if (length(target_col) == 0) return(character())

        alias_map <- character()
        targets <- normalize_channel(cyto_db[[target_col[1]]])
        for (i in seq_len(nrow(cyto_db))) {
            target <- targets[i]
            if (!nzchar(target)) next
            aliases <- normalize_channel(unlist(cyto_db[i, , drop = FALSE]))
            aliases <- aliases[nzchar(aliases)]
            for (alias in aliases) {
                if (!alias %in% names(alias_map)) alias_map[alias] <- target
            }
        }
        alias_map
    }

    channel_alias_map <- load_cytometer_alias_map(cytometer)

    resolve_control_channel <- function(channel_value, det_names) {
        if (is.null(channel_value) || is.na(channel_value) || trimws(channel_value) == "") return("")
        if (channel_value %in% det_names) return(channel_value)

        det_norm <- normalize_channel(det_names)
        name_by_norm <- stats::setNames(det_names, det_norm)
        key <- normalize_channel(channel_value)
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

    message("Found ", length(detector_names), " spectral detectors. Sorting by laser...")

    get_ellipse <- function(mean, sigma, level = 0.95, n = 100, scale = 1.0) {
        chi2_val <- qchisq(level, df = 2)
        eig <- eigen(sigma)
        a <- sqrt(eig$values[1] * chi2_val) * scale
        b <- sqrt(eig$values[2] * chi2_val) * scale
        angle <- atan2(eig$vectors[2, 1], eig$vectors[1, 1])
        theta <- seq(0, 2 * pi, length.out = n)
        ellipse_x <- a * cos(theta)
        ellipse_y <- b * sin(theta)
        rot_x <- ellipse_x * cos(angle) - ellipse_y * sin(angle)
        rot_y <- ellipse_x * sin(angle) + ellipse_y * cos(angle)
        data.table::data.table(x = rot_x + mean[1], y = rot_y + mean[2])
    }

    fit_gmm_populations <- function(data, max_k = 5, min_prop = 0.05) {
        # mclust::Mclust evaluates mclustBIC() in parent.frame(), so bind it here.
        mclustBIC <- get("mclustBIC", envir = asNamespace("mclust"))
        fit <- mclust::Mclust(data, G = 1:max_k, verbose = FALSE)
        if (is.null(fit)) {
            return(NULL)
        }
        sigmas <- lapply(1:fit$G, function(k) {
            if (fit$modelName %in% c("EII", "VII")) {
                diag(fit$parameters$variance$sigmasq[k], nrow = 2)
            } else {
                fit$parameters$variance$sigma[, , k]
            }
        })
        list(fit = fit, proportions = fit$parameters$pro, means = fit$parameters$mean, sigmas = sigmas, main_populations = which(fit$parameters$pro >= min_prop))
    }

    get_sample_type <- function(filename, patterns, default) {
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

    select_bead_population <- function(gmm_result) {
        if (length(gmm_result$main_populations) == 0) {
            return(NULL)
        }
        best <- gmm_result$main_populations[which.max(gmm_result$proportions[gmm_result$main_populations])]
        list(selected = best)
    }

    select_cell_populations <- function(gmm_result, debris_threshold, detection_limit) {
        valid <- c()
        for (k in gmm_result$main_populations) {
            if (gmm_result$means[1, k] >= debris_threshold && gmm_result$means[1, k] <= detection_limit * 0.8) valid <- c(valid, k)
        }
        list(selected = valid)
    }

    create_merged_gate <- function(gmm_result, populations, level, scale = 1.0, clip_x = Inf, clip_y = Inf) {
        if (length(populations) == 0) {
            return(NULL)
        }
        if (length(populations) == 1) {
            ell <- get_ellipse(gmm_result$means[, populations], gmm_result$sigmas[[populations]], level, scale = scale)
        } else {
            all_pts <- data.table::rbindlist(lapply(populations, function(k) get_ellipse(gmm_result$means[, k], gmm_result$sigmas[[k]], level, scale = scale)))
            ell <- all_pts[grDevices::chull(x, y), ]
        }
        ell$x <- pmax(0, pmin(ell$x, clip_x))
        ell$y <- pmax(0, pmin(ell$y, clip_y))
        return(ell)
    }

    select_hist_peak_population <- function(vals_log, rel_height_cutoff = 0.15, min_points = 20) {
        vals_log <- vals_log[is.finite(vals_log)]
        if (length(vals_log) < min_points) return(vals_log)

        d <- density(vals_log, n = 1024)
        x <- d$x
        y <- d$y

        # Local maxima/minima on smoothed density
        peak_idx <- which(diff(sign(diff(y))) == -2) + 1
        trough_idx <- which(diff(sign(diff(y))) == 2) + 1
        if (length(peak_idx) == 0) return(vals_log)

        # Keep meaningful peaks and pick the brightest-rightmost one.
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

    results_list <- list()
    qc_summary_list <- list()

    # 1.1 Add files from AF directory if requested
    fcs_files_all <- fcs_files
    if (include_multi_af && dir.exists(af_dir)) {
        af_files <- list.files(af_dir, pattern = "\\.fcs$", full.names = TRUE)
        message("Found ", length(af_files), " extra AF files in '", af_dir, "'")
        # Pre-pend AF files so they are processed
        fcs_files_all <- c(af_files, fcs_files)
    }

    af_data_raw <- NULL
    af_fn <- NULL
    if (!is.null(control_df)) {
        af_rows <- if ("fluorophore" %in% colnames(control_df)) {
            control_df[as.character(control_df$fluorophore) == "AF", ]
        } else {
            data.frame()
        }
        if (nrow(af_rows) > 0) af_fn <- tools::file_path_sans_ext(basename(af_rows$filename[1]))
    }
    if (is.null(af_fn)) {
        af_idx_tmp <- grep("Unstained|US_UT", fcs_files, ignore.case = TRUE)
        if (length(af_idx_tmp) > 0) af_fn <- tools::file_path_sans_ext(basename(fcs_files[af_idx_tmp[1]]))
    }
    if (!is.null(af_fn)) {
        af_path <- fcs_files[grep(af_fn, fcs_files, fixed = TRUE)]
        if (length(af_path) > 0) {
            ff_af <- flowCore::read.FCS(af_path[1], transformation = FALSE, truncate_max_range = FALSE)
            af_data_raw <- apply(flowCore::exprs(ff_af)[, detector_names, drop = FALSE], 2, median)
        }
    }

    for (fcs_file in fcs_files_all) {
        sn_ext <- basename(fcs_file)
        sn <- tools::file_path_sans_ext(sn_ext)

        # Determine if this is an extra AF file
        is_extra_af <- FALSE
        if (include_multi_af && grepl(normalizePath(af_dir), normalizePath(fcs_file), fixed = TRUE)) {
            is_extra_af <- TRUE
        }

        row_info <- get_control_rows(control_df, c(sn_ext, sn))
        sample_info <- get_sample_type(sn, sample_patterns, default_sample_type)

        fluor_name <- if (nrow(row_info) > 0 && !is.na(row_info$fluorophore[1])) row_info$fluorophore[1] else sample_info$pattern

        if (is_extra_af) {
            # Assign a unique AF name if not in control_df
            if (nrow(row_info) == 0) {
                fluor_name <- paste0("AF_", sn)
            }
            sample_info$type <- "cells" # Force cells for AF folder
        }

        if (fluor_name == "AF" && !is_extra_af) next

        message("Processing SCC: ", fluor_name, " (", sn, ")")
        ff <- flowCore::read.FCS(fcs_file, transformation = FALSE, truncate_max_range = FALSE)
        pd <- flowCore::pData(flowCore::parameters(ff))
        raw_data <- flowCore::exprs(ff)

        # Sanity check for data integrity
        if (any(is.infinite(raw_data))) {
            message("  Skipping ", sn, ": Contains infinite values.")
            next
        }
        na_prop <- sum(is.na(raw_data)) / length(raw_data)
        if (na_prop > 0.1) {
            message("  Skipping ", sn, ": Too many NAs (", round(na_prop * 100, 1), "%).")
            next
        }
        max_val <- max(raw_data, na.rm = TRUE)
        if (max_val > 1e9) {
            message("  Skipping ", sn, ": Extreme values detected (max > 1e9). File may be corrupted.")
            next
        }

        fsc <- pd$name[grepl("^FSC", pd$name) & grepl("-A$", pd$name)][1]
        ssc <- pd$name[grepl("^SSC", pd$name) & grepl("-A$", pd$name)][1]

        data_raw_scatter <- raw_data[, c(fsc, ssc)]
        fsc_max <- quantile(data_raw_scatter[, 1], 1 - outlier_percentile, na.rm = TRUE)
        ssc_max <- quantile(data_raw_scatter[, 2], 1 - outlier_percentile, na.rm = TRUE)
        valid_idx <- which(data_raw_scatter[, 1] < fsc_max & data_raw_scatter[, 2] < ssc_max & data_raw_scatter[, 1] > 0 & data_raw_scatter[, 2] > 0)
        data_filtered <- data_raw_scatter[valid_idx, ]
        if (sample_info$type %in% c("cells", "unstained")) {
            debris_threshold <- quantile(data_filtered[, 1], debris_percentile, na.rm = TRUE)
            data_filtered <- data_filtered[data_filtered[, 1] >= debris_threshold, ]
        } else {
            debris_threshold <- 0
        }

        if (!is.null(subsample_n) && nrow(data_filtered) > subsample_n) {
            data_fit <- data_filtered[sample(nrow(data_filtered), subsample_n), ]
        } else {
            data_fit <- data_filtered
        }

        gmm_result <- fit_gmm_populations(data_fit, max_k = max_clusters, min_prop = min_cluster_proportion)
        if (is.null(gmm_result)) next
        if (sample_info$type == "beads") {
            selected_pops <- select_bead_population(gmm_result)$selected
            gate_level <- gate_contour_beads
        } else {
            selected_pops <- select_cell_populations(gmm_result, debris_threshold, min(fsc_max, ssc_max))$selected
            gate_level <- gate_contour_cells
            if (length(selected_pops) == 0) selected_pops <- select_bead_population(gmm_result)$selected
        }
        if (length(selected_pops) == 0) next
        final_gate <- create_merged_gate(gmm_result, selected_pops, gate_level, scale = if (sample_info$type == "beads") bead_gate_scale else 1.0, clip_x = Inf, clip_y = Inf)

        gated_data <- raw_data[sp::point.in.polygon(raw_data[, fsc], raw_data[, ssc], final_gate$x, final_gate$y) > 0, ]
        if (nrow(gated_data) < 100) next

        q999_by_channel <- apply(
            gated_data[, detector_names, drop = FALSE],
            2,
            function(x) stats::quantile(x, 0.999, na.rm = TRUE)
        )
        inferred_peak_channel <- detector_names[which.max(q999_by_channel)]
        peak_channel <- inferred_peak_channel
        if (nrow(row_info) > 0 && !is.na(row_info$channel[1]) && row_info$channel[1] != "") {
            resolved_channel <- resolve_control_channel(row_info$channel[1], detector_names)
            if (nzchar(resolved_channel)) {
                # Guardrail: if control-file channel is clearly weak vs inferred peak,
                # prefer inferred channel to avoid collapsing bimodal histograms.
                inferred_val <- as.numeric(q999_by_channel[inferred_peak_channel])
                resolved_val <- as.numeric(q999_by_channel[resolved_channel])
                ranked <- names(sort(q999_by_channel, decreasing = TRUE))
                resolved_rank <- match(resolved_channel, ranked)
                use_inferred <- is.finite(inferred_val) &&
                    is.finite(resolved_val) &&
                    inferred_val > 0 &&
                    !is.na(resolved_rank) &&
                    resolved_rank > 8 &&
                    (resolved_val / inferred_val) < 0.25

                if (use_inferred) {
                    warning(
                        "Control channel '", row_info$channel[1], "' for ", sn_ext,
                        " appears inconsistent with signal profile (rank ",
                        resolved_rank, ", 99.9% ratio ",
                        round(resolved_val / inferred_val, 3),
                        "). Falling back to inferred channel ",
                        inferred_peak_channel, "."
                    )
                    peak_channel <- inferred_peak_channel
                } else {
                    peak_channel <- resolved_channel
                }
            } else {
                warning("Control channel '", row_info$channel[1], "' for ", sn_ext, " not found in file. Falling back to inferred channel ", inferred_peak_channel, ".")
            }
        }
        message("  Peak channel: ", peak_channel)
        peak_vals <- gated_data[, peak_channel]
        vals_log <- log10(pmax(peak_vals, 1))
        vals_for_gate <- if (sample_info$type %in% c("unstained", "cells")) vals_log else select_hist_peak_population(vals_log)
        pct <- if (sample_info$type %in% c("unstained", "cells")) histogram_pct_cells else histogram_pct_beads
        dir <- if (sample_info$type %in% c("unstained", "cells")) histogram_direction_cells else histogram_direction_beads
        pct <- min(max(pct, 0.01), 0.999)
        if (dir == "right") {
            lq <- 0.5
            uq <- min(1, 0.5 + pct)
        } else if (dir == "left") {
            lq <- max(0, 0.5 - pct)
            uq <- 0.5
        } else {
            lq <- 0.5 - pct / 2
            uq <- 0.5 + pct / 2
        }
        gate_min <- 10^quantile(vals_for_gate, max(0, lq), na.rm = TRUE)
        gate_max <- 10^quantile(vals_for_gate, min(1, uq), na.rm = TRUE)
        if (!is.finite(gate_min) || !is.finite(gate_max) || gate_max <= gate_min) {
            gate_min <- 10^min(vals_for_gate, na.rm = TRUE)
            gate_max <- 10^max(vals_for_gate, na.rm = TRUE)
        }
        final_gated_data <- gated_data[peak_vals >= gate_min & peak_vals <= gate_max, ]
        if (nrow(final_gated_data) < 10) next

        pos_spectrum_raw <- apply(final_gated_data[, detector_names, drop = FALSE], 2, median, na.rm = TRUE)
        neg_spectrum_raw <- apply(gated_data[peak_vals <= 10^quantile(vals_log, 0.15, na.rm = TRUE), detector_names, drop = FALSE], 2, median, na.rm = TRUE)
        uv_val <- if (nrow(row_info) > 0 && "universal.negative" %in% colnames(row_info)) {
            trimws(as.character(row_info$universal.negative[1]))
        } else {
            ""
        }
        use_univ <- toupper(uv_val) %in% c("TRUE", "AF")
        final_neg <- if (use_univ && !is.null(af_data_raw)) af_data_raw else neg_spectrum_raw
        sig_pure <- pmax(pos_spectrum_raw - final_neg, 0)
        if (max(sig_pure, na.rm = TRUE) <= 0) sig_pure <- pmax(pos_spectrum_raw, 0)
        spectrum_norm <- sig_pure / max(sig_pure, na.rm = TRUE)

        if (save_qc_plots) {
            # FSC/SSC plot
            get_axis_label <- function(ch_name, pd_tbl) {
                idx <- match(ch_name, as.character(pd_tbl$name))
                if (!is.na(idx) && "desc" %in% colnames(pd_tbl)) {
                    desc_val <- trimws(as.character(pd_tbl$desc[idx]))
                    if (!is.na(desc_val) && nzchar(desc_val)) {
                        return(desc_val)
                    }
                }
                ch_name
            }
            fsc_desc <- get_axis_label(fsc, pd)
            ssc_desc <- get_axis_label(ssc, pd)
            dt_plot_gating <- data.table::data.table(x = raw_data[, fsc], y = raw_data[, ssc])
            x_breaks <- seq(0, max(fsc_max, ssc_max) * 1.05, length.out = 201)
            y_breaks <- x_breaks
            x_bin <- findInterval(dt_plot_gating$x, x_breaks)
            y_bin <- findInterval(dt_plot_gating$y, y_breaks)
            keep_bin <- x_bin >= 1 & x_bin <= 200 & y_bin >= 1 & y_bin <= 200
            dt2d <- data.table::as.data.table(as.data.frame(table(
                x_bin = x_bin[keep_bin],
                y_bin = y_bin[keep_bin]
            )))
            data.table::setnames(dt2d, c("x_bin", "y_bin", "count"))
            dt2d$x_bin <- as.integer(as.character(dt2d$x_bin))
            dt2d$y_bin <- as.integer(as.character(dt2d$y_bin))
            dt2d$count <- as.integer(as.character(dt2d$count))
            dt2d <- dt2d[dt2d$count > 0, ]
            dt2d$x <- x_breaks[dt2d$x_bin]
            dt2d$y <- y_breaks[dt2d$y_bin]
            dt2d$fill <- log10(dt2d$count + 1)
            p1 <- ggplot2::ggplot(dt2d, ggplot2::aes(x, y, fill = fill)) +
                ggplot2::geom_tile(width = diff(x_breaks)[1], height = diff(y_breaks)[1]) +
                ggplot2::scale_fill_gradientn(colors = c("#0000FF", "#00FFFF", "#00FF00", "#FFFF00", "#FF0000"), guide = "none") +
                ggplot2::geom_path(data = final_gate, ggplot2::aes(x, y), inherit.aes = FALSE, color = "red", linewidth = 1) +
                ggplot2::labs(title = paste0(sn, " - FSC/SSC"), subtitle = paste0(round(100 * nrow(gated_data) / nrow(raw_data), 1), "% gated"), x = fsc_desc, y = ssc_desc) +
                ggplot2::theme_minimal() +
                ggplot2::theme(legend.position = "none", panel.grid = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = "white", color = NA)) +
                ggplot2::coord_cartesian(xlim = c(0, max(fsc_max, ssc_max) * 1.05), ylim = c(0, max(fsc_max, ssc_max) * 1.05))
            ggplot2::ggsave(file.path(out_path, "fsc_ssc", paste0(sn, "_fsc_ssc.png")), p1, width = 5, height = 5, dpi = 300)

            p2 <- ggplot2::ggplot(data.table::data.table(x = vals_log), ggplot2::aes(x)) +
                ggplot2::geom_density(fill = "grey80", color = "grey40") +
                ggplot2::geom_vline(xintercept = log10(gate_min), color = "red", linewidth = 1) +
                ggplot2::geom_vline(xintercept = log10(gate_max), color = "red", linewidth = 1) +
                ggplot2::annotate("rect", xmin = log10(gate_min), xmax = log10(gate_max), ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "red") +
                ggplot2::labs(title = paste0(sn, " - ", peak_channel), subtitle = paste0(round(100 * nrow(final_gated_data) / nrow(gated_data), 1), "% gated"), x = paste0("log10(", peak_channel, ")")) +
                ggplot2::theme_minimal() +
                ggplot2::theme(legend.position = "none")
            ggplot2::ggsave(file.path(out_path, "histogram", paste0(sn, "_histogram.png")), p2, width = 6.5, height = 4, dpi = 300)

            # Spectrum plot logic - based on working code from autogate_contour.R
            log_mat <- log10(pmax(final_gated_data[, detector_names, drop = FALSE], 1e-3))
            min_y <- floor(min(log_mat, na.rm = TRUE))
            max_y <- ceiling(max(log_mat, na.rm = TRUE))
            breaks <- seq(min_y, max_y, length.out = 151)
            bin_mid <- (breaks[-1] + breaks[-length(breaks)]) / 2
            bin_height <- breaks[2] - breaks[1]
            counts_mat <- vapply(seq_len(ncol(log_mat)), function(j) as.numeric(graphics::hist(log_mat[, j], breaks = breaks, plot = FALSE)$counts), numeric(length(bin_mid)))
            rownames(counts_mat) <- as.character(seq_along(bin_mid))
            colnames(counts_mat) <- as.character(seq_len(ncol(log_mat)))
            dt_c <- data.table::as.data.table(as.table(counts_mat))
            data.table::setnames(dt_c, c("bin_idx", "ch_idx", "count"))
            dt_c$bin_idx <- as.integer(as.character(dt_c$bin_idx))
            dt_c$ch_idx <- as.integer(as.character(dt_c$ch_idx))
            dt_c$y_orig <- bin_mid[dt_c$bin_idx]
            dt_c$fill <- log10(dt_c$count + 1)
            min_bin_count <- 3
            dt_c <- dt_c[dt_c$count >= min_bin_count, ]
            y_power <- 1.5
            dt_c$y <- dt_c$y_orig^y_power

            if (nrow(dt_c) == 0) {
                message("  Warning: dt_c is empty for ", sn, ". Max count: ", max(counts_mat, na.rm = TRUE))
            } else {
                message("  Plotting ", nrow(dt_c), " bins for ", sn)
            }

            vlines <- which(diff(det_info$laser_nm) != 0) + 0.5
            fill_lo <- min(dt_c$fill, na.rm = TRUE)
            fill_hi <- quantile(dt_c$fill, 0.96, na.rm = TRUE)
            y_breaks_orig <- 0:ceiling(max_y)
            y_breaks_trans <- y_breaks_orig^y_power
            y_labels <- sapply(y_breaks_orig, function(x) bquote(10^.(x)))

            p3 <- ggplot2::ggplot(dt_c, ggplot2::aes(ch_idx, y, fill = fill)) +
                ggplot2::geom_tile(width = 0.7, height = bin_height * 3) +
                ggplot2::scale_fill_gradientn(
                    colors = c("#0000FF", "#00FFFF", "#00FF00", "#FFFF00", "#FF0000"),
                    limits = c(fill_lo, fill_hi), oob = scales::squish
                ) +
                ggplot2::scale_x_continuous(breaks = seq_along(detector_names), labels = detector_labels) +
                ggplot2::scale_y_continuous(limits = c(0, (max_y + 0.5)^y_power), breaks = y_breaks_trans, labels = y_labels) +
                ggplot2::coord_cartesian(expand = FALSE) +
                ggplot2::labs(title = paste0(sn, " - Spectrum"), x = NULL, y = "Intensity") +
                ggplot2::theme_minimal() +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6), legend.position = "none", panel.background = ggplot2::element_rect(fill = "white", color = NA))
            ggplot2::ggsave(file.path(out_path, "spectrum", paste0(sn, "_spectrum.png")), p3, width = 300, height = 120, units = "mm", dpi = 600)
        }

        results_list[[sn]] <- data.table::data.table(sample = sn, fluorophore = fluor_name, type = sample_info$type, n_total = nrow(raw_data), n_final = nrow(final_gated_data), spectrum = list(spectrum_norm))
        qc_summary_list[[sn]] <- data.table::data.table(
            sample = sn,
            fluorophore = fluor_name,
            type = sample_info$type,
            peak_channel = peak_channel,
            fsc_channel = fsc,
            ssc_channel = ssc,
            n_total = nrow(raw_data),
            n_scatter_gated = nrow(gated_data),
            n_final = nrow(final_gated_data),
            scatter_gate_pct = round(100 * nrow(gated_data) / max(nrow(raw_data), 1), 1),
            histogram_gate_pct = round(100 * nrow(final_gated_data) / max(nrow(gated_data), 1), 1)
        )
    }
    results_dt <- data.table::rbindlist(results_list)
    spectra_list <- results_dt$spectrum
    names(spectra_list) <- results_dt$fluorophore
    if (!is.null(af_data_raw)) spectra_list[["AF"]] <- af_data_raw / max(af_data_raw, na.rm = TRUE)

    if (length(spectra_list) == 0) {
        warning("No valid spectra found!")
        return(NULL)
    }

    M <- do.call(rbind, spectra_list)
    colnames(M) <- detector_names
    M_df <- as.data.frame(M)
    M_df$file <- rownames(M_df)
    output_root <- file.path("spectreasy_outputs")
    dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
    reference_matrix_file <- file.path(output_root, "reference_matrix.csv")
    utils::write.csv(M_df[, c("file", colnames(M)), drop = FALSE], reference_matrix_file, row.names = FALSE, quote = TRUE)
    message("Reference matrix (", nrow(M), " markers) saved to: ", reference_matrix_file)
    if (length(qc_summary_list) > 0) {
        attr(M, "qc_summary") <- data.table::rbindlist(qc_summary_list)
    }
    if (isTRUE(save_qc_plots)) {
        attr(M, "qc_plot_dir") <- out_path
    }
    return(M)
}
