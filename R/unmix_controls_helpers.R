# Internal helpers for unmix_controls().
.unmix_confirm_created_control_file <- function(path) {
    msg <- paste(
        c(
            "A new control file was created:",
            paste0(" - ", path),
            "Please review it before unmixing.",
            "Check at least these columns:",
            " - fluorophore / marker / channel mappings",
            " - control.type: auto-detected from filename tokens ('beads'/'cells'); fix unknown rows as needed"
        ),
        collapse = "\n"
    )

    if (!interactive()) {
        stop(
            msg, "\n",
            "R session is non-interactive, so confirmation prompt is not possible.\n",
            "After reviewing the file, rerun unmix_controls().",
            call. = FALSE
        )
    }

    .spectreasy_console_header("control file created")
    .spectreasy_console_field("File", .spectreasy_console_path(path))
    .spectreasy_console_field("Review", "marker, fluorophore, channel, and control.type columns")
    .spectreasy_console_footer(blank = FALSE)
    repeat {
        ans <- tolower(trimws(readline("Proceed with unmix_controls now? [y/n]: ")))
        if (ans %in% c("y", "yes")) return(invisible(TRUE))
        if (ans %in% c("n", "no", "")) {
            stop("Stopped after control-file creation. Review and rerun unmix_controls() when ready.", call. = FALSE)
        }
        message("Please answer 'y' or 'n'.")
    }
}

.unmix_read_control_df <- function(control_file) {
    tryCatch(
        utils::read.csv(control_file, stringsAsFactors = FALSE, check.names = FALSE),
        error = function(e) NULL
    )
}

.unmix_generate_control_file <- function(scc_dir,
                                         control_file,
                                         cytometer,
                                         auto_unknown_fluor_policy) {
    .spectreasy_console_field("Controls", paste0(.spectreasy_console_path(control_file), " not found"))
    .spectreasy_console_step("Create mapping", "from SCC filenames and peak channels")
    create_control_file(
        input_folder = scc_dir,
        cytometer = cytometer,
        unknown_fluor_policy = auto_unknown_fluor_policy,
        output_file = control_file
    )
    .spectreasy_console_field("Controls", paste0(.spectreasy_console_path(control_file), " created"))
    invisible(control_file)
}

.prepare_unmix_control_df <- function(control_file,
                                      scc_dir,
                                      auto_create_mapping,
                                      cytometer,
                                      auto_unknown_fluor_policy) {
    created_control_file <- FALSE

    if (!file.exists(control_file)) {
        if (!isTRUE(auto_create_mapping)) {
            stop(
                "Control file not found: ", control_file, "\n",
                "Set auto_create_mapping = TRUE to auto-generate it from SCC files."
            )
        }
        .unmix_generate_control_file(
            scc_dir = scc_dir,
            control_file = control_file,
            cytometer = cytometer,
            auto_unknown_fluor_policy = auto_unknown_fluor_policy
        )
        created_control_file <- TRUE
    }
    control_df <- .unmix_read_control_df(control_file)
    if (is.null(control_df)) {
        stop("Could not read control file: ", control_file)
    }

    list(control_df = control_df, control_file = control_file, created_control_file = created_control_file)
}

.normalize_unmix_control_df <- function(control_df, required_cols = c("filename", "fluorophore", "channel")) {
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

.drop_missing_unmix_af_rows <- function(control_df, scc_dir, emit_message = TRUE) {
    if (is.null(control_df) || nrow(control_df) == 0 || !("filename" %in% colnames(control_df))) {
        return(control_df)
    }

    af_rows <- .is_af_control_row(
        fluorophore = if ("fluorophore" %in% colnames(control_df)) control_df$fluorophore else NULL,
        marker = if ("marker" %in% colnames(control_df)) control_df$marker else NULL,
        filename = control_df$filename
    )
    if (!any(af_rows)) {
        return(control_df)
    }

    af_files <- file.path(scc_dir, as.character(control_df$filename[af_rows]))
    missing_af <- af_rows
    missing_af[af_rows] <- !file.exists(af_files)
    if (!any(missing_af)) {
        return(control_df)
    }

    if (isTRUE(emit_message)) {
        .spectreasy_console_step("AF controls", "missing; building marker-only reference matrix")
        .spectreasy_console_step("Next option", "use unmix_samples(estimate_af = TRUE) for stained-sample AF")
    }
    control_df[!missing_af, , drop = FALSE]
}

.run_unmix_preflight <- function(control_df, scc_dir) {
    validate_control_file_mapping(
        control_df = control_df,
        scc_dir = scc_dir,
        require_all_scc_mapped = FALSE,
        require_channels = TRUE,
        stop_on_error = FALSE
    )
}

.stop_unmix_preflight <- function(preflight, auto_unknown_fluor_policy) {
    hint <- NULL
    if (any(grepl("^Empty fluorophore", preflight$errors)) && identical(auto_unknown_fluor_policy, "empty")) {
        hint <- "Tip: rerun with auto_unknown_fluor_policy = \"by_channel\" to auto-fill common fluorophore names from detected channels."
    }
    stop(
        paste(
            c(
                "unmix_controls preflight failed:",
                paste0(" - ", preflight$errors),
                if (length(preflight$warnings) > 0) c("Preflight notes:", paste0(" - ", preflight$warnings)) else NULL,
                hint,
                "Fix the control file and rerun unmix_controls()."
            ),
            collapse = "\n"
        )
    )
}

.unmix_controls_stage_dir <- function(output_dir) {
    output_dir <- as.character(output_dir)[1]
    file.path(output_dir, "unmix_controls")
}

.unmix_output_paths <- function(output_dir) {
    list(
        unmixed_dir = file.path(output_dir, "unmixed_fcs"),
        qc_controls_dir = file.path(output_dir, "qc_controls"),
        spectra_file = file.path(output_dir, "scc_spectra.png"),
        af_spectra_file = file.path(output_dir, "scc_af_spectra.png"),
        reference_matrix_csv = file.path(output_dir, "scc_reference_matrix.csv"),
        detector_noise_csv = file.path(output_dir, "scc_detector_noise.csv"),
        unmixing_matrix_csv = file.path(output_dir, "scc_unmixing_matrix.csv"),
        control_mapping_csv = file.path(output_dir, "fcs_mapping_used.csv"),
        spectral_variant_library_rds = file.path(output_dir, "scc_spectral_variants.rds"),
        unmixing_scatter_png = file.path(output_dir, "scc_unmixing_scatter_matrix.png"),
        qc_report_html = file.path(output_dir, "qc_controls", "qc_controls_report.html")
    )
}

.run_optional_unmix_artifact <- function(label, expr) {
    tryCatch(
        force(expr),
        error = function(e) {
            warning(
                label, " could not be created: ",
                conditionMessage(e),
                ". Continuing unmixing.",
                call. = FALSE
            )
            NULL
        }
    )
}

.save_reference_matrix_csv <- function(M, path) {
    M_df <- as.data.frame(M, check.names = FALSE)
    M_df$Marker <- rownames(M)
    M_df <- M_df[, c("Marker", setdiff(colnames(M_df), "Marker")), drop = FALSE]
    utils::write.csv(M_df, path, row.names = FALSE, quote = TRUE)
    invisible(path)
}

.read_unmix_metadata_pd <- function(scc_dir, control_df = NULL) {
    if (!dir.exists(scc_dir)) {
        .spectreasy_stop_missing_directory(scc_dir, label = "scc_dir")
    }
    fcs_files <- list.files(scc_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
    if (length(fcs_files) == 0) {
        .spectreasy_stop_empty_fcs_directory(scc_dir, label = "scc_dir")
    }
    fcs_files <- .control_validation_select_scc_files(control_df, fcs_files)
    if (length(fcs_files) == 0) {
        stop("None of the FCS files listed in the control file were found in scc_dir.", call. = FALSE)
    }
    ff_meta <- .spectreasy_read_fcs(fcs_files[1], label = "SCC FCS file")
    list(fcs_files = fcs_files, pd = flowCore::pData(flowCore::parameters(ff_meta)))
}

.derive_unmix_static_matrix <- function(M, fcs_files, unmixing_method) {
    # Static exported matrices cannot represent per-event AF selection reliably.
    # Drop all AF rows and keep AF modeling in dynamic unmixing results.
    marker_names <- rownames(M)
    af_rows <- grepl("^AF($|_)", marker_names, ignore.case = TRUE)
    if (any(af_rows)) {
        M_static <- M[!af_rows, , drop = FALSE]
        detector_noise <- attr(M, "detector_noise")
        if (!is.null(detector_noise)) {
            attr(M_static, "detector_noise") <- detector_noise
        }
    } else {
        M_static <- M
    }

    static_unmixing_matrix_method <- .normalize_unmix_method(unmixing_method)
    static_solver_method <- .solver_method_for_unmix(static_unmixing_matrix_method)
    if (static_solver_method %in% c("WLS", "RWLS")) {
        W <- tryCatch(
            derive_unmixing_matrix(M_static, method = "WLS"),
            error = function(e) {
                if (!grepl("singular", conditionMessage(e), ignore.case = TRUE)) {
                    stop(e)
                }
                .derive_regularized_static_wls_matrix(M_static)
            }
        )
    } else {
        W <- derive_unmixing_matrix(M_static, method = static_solver_method)
    }

    list(W = W, static_unmixing_matrix_method = static_unmixing_matrix_method)
}

.derive_regularized_static_wls_matrix <- function(M) {
    detector_weights <- .wls_static_detector_weights(
        M = M,
        background_noise = .default_wls_background_noise(),
        signal_scale = .default_wls_signal_scale(),
        max_weight_ratio = .default_wls_max_weight_ratio()
    )
    V_inv <- diag(detector_weights)
    Mt <- t(M)
    MVMt <- M %*% V_inv %*% Mt
    ridge <- max(c(diag(MVMt), 1), na.rm = TRUE) * 1e-8
    W <- solve(MVMt + diag(ridge, nrow(MVMt)), M %*% V_inv)
    rownames(W) <- rownames(M)
    colnames(W) <- colnames(M)
    W
}

.resolve_unmix_marker_mappings <- function(control_df) {
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
