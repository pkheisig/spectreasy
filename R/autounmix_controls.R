# Internal helpers for autounmix_controls().
.autounmix_confirm_created_control_file <- function(path) {
    msg <- paste(
        c(
            "A new control file was created:",
            paste0(" - ", path),
            "Please review it before unmixing.",
            "Check at least these columns:",
            " - fluorophore / marker / channel mappings",
            " - control.type: auto-detected from filename tokens ('beads'/'cells'); fix unknown rows as needed",
            " - universal.negative: leave empty unless you explicitly use it"
        ),
        collapse = "\n"
    )

    if (!interactive()) {
        stop(
            msg, "\n",
            "R session is non-interactive, so confirmation prompt is not possible.\n",
            "After reviewing the file, rerun autounmix_controls().",
            call. = FALSE
        )
    }

    message(msg)
    repeat {
        ans <- tolower(trimws(readline("Proceed with autounmix_controls now? [y/n]: ")))
        if (ans %in% c("y", "yes")) return(invisible(TRUE))
        if (ans %in% c("n", "no", "")) {
            stop("Stopped after control-file creation. Review and rerun autounmix_controls() when ready.", call. = FALSE)
        }
        message("Please answer 'y' or 'n'.")
    }
}

.autounmix_read_control_df <- function(control_file) {
    tryCatch(
        utils::read.csv(control_file, stringsAsFactors = FALSE, check.names = FALSE),
        error = function(e) NULL
    )
}

.autounmix_generate_control_file <- function(scc_dir,
                                             control_file,
                                             cytometer,
                                             auto_default_control_type,
                                             auto_unknown_fluor_policy) {
    message("Control file not found: ", control_file)
    message("Auto-generating control file from SCC filenames and peak channels...")
    create_control_file(
        input_folder = scc_dir,
        include_af_folder = FALSE,
        cytometer = cytometer,
        default_control_type = auto_default_control_type,
        unknown_fluor_policy = auto_unknown_fluor_policy,
        output_file = control_file
    )
    message("Auto-generated control file: ", control_file, " (please review marker/fluorophore/channel columns).")
    invisible(control_file)
}

.prepare_autounmix_control_df <- function(control_df,
                                          control_file,
                                          scc_dir,
                                          auto_create_control,
                                          cytometer,
                                          auto_default_control_type,
                                          auto_unknown_fluor_policy) {
    created_control_file <- FALSE

    if (is.character(control_df) && length(control_df) == 1 && !is.na(control_df)) {
        control_file <- control_df
        control_df <- NULL
    }
    if (!is.null(control_df) && !is.data.frame(control_df)) {
        stop("control_df must be either a data.frame or a single CSV path.")
    }

    if (is.null(control_df)) {
        if (!file.exists(control_file)) {
            if (!isTRUE(auto_create_control)) {
                stop(
                    "Control file not found: ", control_file, "\n",
                    "Set auto_create_control = TRUE to auto-generate it from SCC files."
                )
            }
            .autounmix_generate_control_file(
                scc_dir = scc_dir,
                control_file = control_file,
                cytometer = cytometer,
                auto_default_control_type = auto_default_control_type,
                auto_unknown_fluor_policy = auto_unknown_fluor_policy
            )
            created_control_file <- TRUE
        }
        control_df <- .autounmix_read_control_df(control_file)
        if (is.null(control_df)) {
            stop("Could not read control file: ", control_file)
        }
    }

    list(control_df = control_df, control_file = control_file, created_control_file = created_control_file)
}

.normalize_autounmix_control_df <- function(control_df, required_cols = c("filename", "fluorophore", "channel")) {
    missing_cols <- setdiff(required_cols, colnames(control_df))
    if (length(missing_cols) > 0) {
        stop(
            "Control mapping is missing required columns: ", paste(missing_cols, collapse = ", "), "\n",
            "Required columns: ", paste(required_cols, collapse = ", "), "\n",
            "Provide at least filename, fluorophore, and channel."
        )
    }
    if (!("universal.negative" %in% colnames(control_df))) {
        control_df$universal.negative <- ""
    }
    control_df
}

.filter_autounmix_af_rows <- function(control_df, exclude_af = FALSE, emit_message = TRUE) {
    excluded_af_filenames <- character()
    if (!isTRUE(exclude_af)) {
        return(list(control_df = control_df, excluded_af_filenames = excluded_af_filenames))
    }

    af_rows <- .is_af_control_row(
        fluorophore = control_df$fluorophore,
        marker = if ("marker" %in% colnames(control_df)) control_df$marker else NULL,
        filename = control_df$filename
    )
    if (any(af_rows)) {
        excluded_af_filenames <- unique(as.character(control_df$filename[af_rows]))
        if (isTRUE(emit_message)) {
            message("exclude_af = TRUE: ignoring ", sum(af_rows), " AF/unstained control row(s) from control mapping.")
        }
        control_df <- control_df[!af_rows, , drop = FALSE]
    }

    list(control_df = control_df, excluded_af_filenames = excluded_af_filenames)
}

.run_autounmix_preflight <- function(control_df, scc_dir, exclude_af = FALSE) {
    validate_control_file_mapping(
        control_df = control_df,
        scc_dir = scc_dir,
        include_multi_af = FALSE,
        exclude_af = exclude_af,
        af_dir = "af",
        require_all_scc_mapped = TRUE,
        require_channels = TRUE,
        stop_on_error = FALSE
    )
}

.autounmix_regenerate_control_df <- function(control_file,
                                             scc_dir,
                                             cytometer,
                                             auto_default_control_type,
                                             auto_unknown_fluor_policy,
                                             exclude_af = FALSE) {
    message("Control preflight failed; attempting automatic control-file regeneration...")
    create_control_file(
        input_folder = scc_dir,
        include_af_folder = FALSE,
        cytometer = cytometer,
        default_control_type = auto_default_control_type,
        unknown_fluor_policy = auto_unknown_fluor_policy,
        output_file = control_file
    )
    control_df <- utils::read.csv(control_file, stringsAsFactors = FALSE, check.names = FALSE)
    filtered <- .filter_autounmix_af_rows(control_df, exclude_af = exclude_af, emit_message = FALSE)
    list(control_df = filtered$control_df, excluded_af_filenames = filtered$excluded_af_filenames)
}

.stop_autounmix_preflight <- function(preflight, auto_unknown_fluor_policy) {
    hint <- NULL
    if (any(grepl("^Empty fluorophore", preflight$errors)) && identical(auto_unknown_fluor_policy, "empty")) {
        hint <- "Tip: rerun with auto_unknown_fluor_policy = \"by_channel\" to auto-fill common fluorophore names from detected channels."
    }
    stop(
        paste(
            c(
                "autounmix_controls preflight failed:",
                paste0(" - ", preflight$errors),
                if (length(preflight$warnings) > 0) c("Preflight notes:", paste0(" - ", preflight$warnings)) else NULL,
                hint,
                "Fix the control file and rerun autounmix_controls()."
            ),
            collapse = "\n"
        )
    )
}

.autounmix_output_paths <- function(output_dir) {
    list(
        build_plots_dir = file.path(output_dir, "build_reference_plots"),
        unmixed_dir = file.path(output_dir, "scc_unmixed"),
        spectra_file = file.path(output_dir, "scc_spectra.png"),
        reference_matrix_csv = file.path(output_dir, "scc_reference_matrix.csv"),
        unmixing_matrix_csv = file.path(output_dir, "scc_unmixing_matrix.csv"),
        unmixing_matrix_png = file.path(output_dir, "scc_unmixing_matrix.png"),
        unmixing_scatter_png = file.path(output_dir, "scc_unmixing_scatter_matrix.png")
    )
}

.save_reference_matrix_csv <- function(M, path) {
    M_df <- as.data.frame(M, check.names = FALSE)
    M_df$Marker <- rownames(M)
    M_df <- M_df[, c("Marker", setdiff(colnames(M_df), "Marker")), drop = FALSE]
    utils::write.csv(M_df, path, row.names = FALSE, quote = TRUE)
    invisible(path)
}

.read_autounmix_metadata_pd <- function(scc_dir) {
    fcs_files <- list.files(scc_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
    if (length(fcs_files) == 0) {
        stop("No FCS files found in scc_dir: ", scc_dir)
    }
    ff_meta <- flowCore::read.FCS(fcs_files[1], transformation = FALSE, truncate_max_range = FALSE)
    list(fcs_files = fcs_files, pd = flowCore::pData(flowCore::parameters(ff_meta)))
}

.filter_autounmix_excluded_af_outputs <- function(unmixed_list, scc_dir, excluded_af_filenames, unmixed_dir, exclude_af = FALSE) {
    if (!isTRUE(exclude_af)) {
        return(unmixed_list)
    }

    scc_basenames <- list.files(scc_dir, pattern = "\\.fcs$", full.names = FALSE, ignore.case = TRUE)
    excluded_af_filenames <- unique(c(excluded_af_filenames, scc_basenames[.is_af_filename(scc_basenames)]))
    excluded_af_keys <- tools::file_path_sans_ext(basename(excluded_af_filenames))
    excluded_af_keys <- excluded_af_keys[!is.na(excluded_af_keys) & excluded_af_keys != ""]
    if (length(excluded_af_keys) == 0) {
        return(unmixed_list)
    }

    unmixed_list <- unmixed_list[!(names(unmixed_list) %in% excluded_af_keys)]
    for (sample_key in excluded_af_keys) {
        out_file <- file.path(unmixed_dir, paste0(sample_key, "_unmixed.fcs"))
        if (file.exists(out_file)) file.remove(out_file)
    }

    unmixed_list
}

.estimate_autounmix_wls_global_weights <- function(files, detector_names, background_noise = 25, max_events_per_file = 50000) {
    sums <- rep(0, length(detector_names))
    counts <- rep(0, length(detector_names))

    for (fp in files) {
        ff <- flowCore::read.FCS(fp, transformation = FALSE, truncate_max_range = FALSE)
        raw <- flowCore::exprs(ff)
        common <- intersect(detector_names, colnames(raw))
        if (length(common) == 0) next

        y <- raw[, common, drop = FALSE]
        if (nrow(y) > max_events_per_file) {
            y <- y[sample.int(nrow(y), max_events_per_file), , drop = FALSE]
        }

        sig <- colMeans(pmax(y, 0), na.rm = TRUE)
        idx <- match(common, detector_names)
        keep <- is.finite(sig)
        sums[idx[keep]] <- sums[idx[keep]] + sig[keep]
        counts[idx[keep]] <- counts[idx[keep]] + 1
    }

    mean_sig <- sums / pmax(counts, 1)
    mean_sig[!is.finite(mean_sig)] <- 0
    1 / (pmax(mean_sig, 0) + background_noise)
}

.derive_autounmix_static_matrix <- function(M, fcs_files, unmix_method) {
    static_unmixing_matrix_method <- toupper(unmix_method)
    if (static_unmixing_matrix_method == "WLS") {
        wls_weights <- .estimate_autounmix_wls_global_weights(fcs_files, colnames(M))
        W <- derive_unmixing_matrix(M, method = "WLS", global_weights = wls_weights)
    } else {
        W <- derive_unmixing_matrix(M, method = static_unmixing_matrix_method)
    }

    list(W = W, static_unmixing_matrix_method = static_unmixing_matrix_method)
}

.resolve_autounmix_marker_mappings <- function(control_df) {
    sample_to_marker <- NULL
    marker_display <- NULL

    if (is.data.frame(control_df) && all(c("filename", "fluorophore") %in% colnames(control_df))) {
        sample_keys <- tools::file_path_sans_ext(basename(as.character(control_df$filename)))
        primary_vals <- trimws(as.character(control_df$fluorophore))
        secondary_vals <- if ("marker" %in% colnames(control_df)) {
            trimws(as.character(control_df$marker))
        } else {
            rep("", length(primary_vals))
        }
        secondary_vals[is.na(secondary_vals)] <- ""
        keep <- !is.na(sample_keys) & sample_keys != "" & !is.na(primary_vals) & primary_vals != ""
        if (any(keep)) {
            sample_to_marker <- stats::setNames(primary_vals[keep], sample_keys[keep])
            sample_to_marker <- sample_to_marker[!duplicated(names(sample_to_marker))]

            display_vals <- primary_vals
            show_secondary <- secondary_vals != "" &
                toupper(secondary_vals) != "AUTOFLUORESCENCE" &
                tolower(secondary_vals) != tolower(primary_vals)
            display_vals[show_secondary] <- paste0(primary_vals[show_secondary], " / ", secondary_vals[show_secondary])
            marker_display <- stats::setNames(display_vals[keep], primary_vals[keep])
            marker_display <- marker_display[names(marker_display) != ""]
            marker_display <- marker_display[!duplicated(names(marker_display))]
        }
    }

    list(sample_to_marker = sample_to_marker, marker_display = marker_display)
}

#' Auto-Unmix Single-Color Controls
#'
#' SCC workflow helper that:
#' 1) validates/creates the control file,
#' 2) builds the reference matrix from SCCs,
#' 3) unmixed SCC files,
#' 4) saves reference/unmixing matrices,
#' 5) plots spectra, unmixing matrix, and SCC unmixing scatter matrix.
#'
#' This function is intended as the control-stage step before GUI adjustment
#' and downstream sample unmixing.
#'
#' @param scc_dir Directory containing SCC FCS files.
#' @param control_df Optional control mapping data.frame, or a path to control CSV.
#' @param control_file Path to control mapping CSV used when `control_df` is `NULL`.
#' @param auto_create_control Logical; auto-generate control file when missing.
#' @param cytometer Cytometer name (for example `"Aurora"`).
#' @param auto_default_control_type Deprecated and ignored.
#' @param auto_unknown_fluor_policy Auto-fill policy for unresolved fluorophores
#'   when creating controls (`"by_channel"`, `"empty"`, `"filename"`).
#' @param output_dir Output directory for SCC workflow artifacts.
#' @param exclude_af Logical; if `TRUE`, ignore AF/unstained controls even when
#'   they are present in the SCC folder or control mapping.
#' @param unmix_method SCC unmixing method (`"WLS"`, `"OLS"`, `"NNLS"`).
#' @param build_qc_plots Logical; keep detailed build_reference_matrix plots.
#' @param unmix_scatter_panel_size_mm Panel size for SCC unmixing scatter matrix plot.
#' @param seed Optional integer seed for deterministic subsampling and plotting.
#' @param ... Additional arguments forwarded to [build_reference_matrix()].
#'
#' @return List with `M`, `W`, `unmixed_list`, and key output file paths.
#' @export
#' @examples
#' if (interactive()) {
#'   ctrl <- autounmix_controls(
#'     scc_dir = "scc",
#'     control_file = "fcs_mapping.csv",
#'     auto_create_control = TRUE,
#'     cytometer = "Aurora",
#'     output_dir = "spectreasy_outputs/autounmix_controls"
#'   )
#'   ctrl$unmixing_matrix_file
#' }
autounmix_controls <- function(
    scc_dir = "scc",
    control_df = NULL,
    control_file = "fcs_mapping.csv",
    auto_create_control = TRUE,
    cytometer = "Aurora",
    auto_default_control_type = "beads",
    auto_unknown_fluor_policy = c("by_channel", "empty", "filename"),
    output_dir = "spectreasy_outputs/autounmix_controls",
    exclude_af = FALSE,
    unmix_method = "WLS",
    build_qc_plots = FALSE,
    unmix_scatter_panel_size_mm = 30,
    seed = NULL,
    ...
) {
    auto_unknown_fluor_policy <- match.arg(auto_unknown_fluor_policy)
    .with_optional_seed(seed)

    user_supplied_control_df <- !is.null(control_df)
    control_file <- .resolve_control_file_path(control_file)

    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    if (!dir.exists(scc_dir)) stop("scc_dir not found: ", scc_dir)

    control_info <- .prepare_autounmix_control_df(
        control_df = control_df,
        control_file = control_file,
        scc_dir = scc_dir,
        auto_create_control = auto_create_control,
        cytometer = cytometer,
        auto_default_control_type = auto_default_control_type,
        auto_unknown_fluor_policy = auto_unknown_fluor_policy
    )
    control_df <- .normalize_autounmix_control_df(control_info$control_df)
    control_file <- control_info$control_file

    filtered <- .filter_autounmix_af_rows(control_df, exclude_af = exclude_af)
    control_df <- filtered$control_df
    excluded_af_filenames <- filtered$excluded_af_filenames

    if (isTRUE(control_info$created_control_file)) {
        .autounmix_confirm_created_control_file(control_file)
    }

    preflight <- .run_autounmix_preflight(control_df, scc_dir = scc_dir, exclude_af = exclude_af)
    if (!preflight$ok && isTRUE(auto_create_control) && !user_supplied_control_df) {
        regenerated <- .autounmix_regenerate_control_df(
            control_file = control_file,
            scc_dir = scc_dir,
            cytometer = cytometer,
            auto_default_control_type = auto_default_control_type,
            auto_unknown_fluor_policy = auto_unknown_fluor_policy,
            exclude_af = exclude_af
        )
        control_df <- regenerated$control_df
        excluded_af_filenames <- unique(c(excluded_af_filenames, regenerated$excluded_af_filenames))
        preflight <- .run_autounmix_preflight(control_df, scc_dir = scc_dir, exclude_af = exclude_af)
    }
    if (!preflight$ok) {
        .stop_autounmix_preflight(preflight, auto_unknown_fluor_policy = auto_unknown_fluor_policy)
    }

    output_paths <- .autounmix_output_paths(output_dir)
    M <- build_reference_matrix(
        input_folder = scc_dir,
        output_folder = output_paths$build_plots_dir,
        save_qc_plots = build_qc_plots,
        control_df = control_df,
        cytometer = cytometer,
        exclude_af = exclude_af,
        seed = seed,
        ...
    )
    if (is.null(M) || nrow(M) == 0) stop("No valid spectra found while building reference matrix.")

    .save_reference_matrix_csv(M, output_paths$reference_matrix_csv)
    meta_info <- .read_autounmix_metadata_pd(scc_dir)
    p_spectra <- plot_spectra(M, pd = meta_info$pd, output_file = output_paths$spectra_file)

    unmixed_list <- unmix_samples(
        sample_dir = scc_dir,
        M = M,
        method = unmix_method,
        cytometer = cytometer,
        output_dir = output_paths$unmixed_dir,
        write_fcs = TRUE
    )
    unmixed_list <- .filter_autounmix_excluded_af_outputs(
        unmixed_list = unmixed_list,
        scc_dir = scc_dir,
        excluded_af_filenames = excluded_af_filenames,
        unmixed_dir = output_paths$unmixed_dir,
        exclude_af = exclude_af
    )

    static_info <- .derive_autounmix_static_matrix(M, fcs_files = meta_info$fcs_files, unmix_method = unmix_method)
    W <- static_info$W
    save_unmixing_matrix(W, output_paths$unmixing_matrix_csv)

    p_unmix <- plot_unmixing_matrix(W, pd = meta_info$pd)
    ggplot2::ggsave(output_paths$unmixing_matrix_png, p_unmix, width = 200, height = 150, units = "mm")

    marker_mapping <- .resolve_autounmix_marker_mappings(control_df)
    plot_unmixing_scatter_matrix(
        unmixed_list = unmixed_list,
        sample_to_marker = marker_mapping$sample_to_marker,
        markers = rownames(M),
        marker_display = marker_mapping$marker_display,
        output_file = output_paths$unmixing_scatter_png,
        transform = "none",
        panel_size_mm = unmix_scatter_panel_size_mm,
        seed = seed
    )

    invisible(list(
        M = M,
        W = W,
        unmixed_list = unmixed_list,
        reference_matrix_file = output_paths$reference_matrix_csv,
        unmixing_matrix_file = output_paths$unmixing_matrix_csv,
        spectra_file = output_paths$spectra_file,
        unmixing_matrix_plot = output_paths$unmixing_matrix_png,
        unmixing_scatter_file = output_paths$unmixing_scatter_png,
        static_unmixing_matrix_method = static_info$static_unmixing_matrix_method,
        spectra_plot = p_spectra,
        unmixing_plot = p_unmix
    ))
}
