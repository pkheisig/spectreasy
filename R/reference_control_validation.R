.active_reference_control_rows <- function(control_df, fcs_files_all) {
    if (is.null(control_df) || !is.data.frame(control_df) || nrow(control_df) == 0 ||
        !all(c("filename", "fluorophore") %in% colnames(control_df))) {
        return(data.frame())
    }

    known_files <- basename(fcs_files_all)
    known_keys <- tools::file_path_sans_ext(known_files)
    control_keys <- tools::file_path_sans_ext(basename(as.character(control_df$filename)))
    active <- as.character(control_df$filename) %in% known_files | control_keys %in% known_keys
    is_af <- .is_af_control_row(
        fluorophore = control_df$fluorophore,
        marker = if ("marker" %in% colnames(control_df)) control_df$marker else NULL,
        filename = control_df$filename
    )
    active <- active & !is_af

    control_df[active, , drop = FALSE]
}

.validate_reference_complete_controls <- function(control_df, fcs_files_all, processed_results) {
    active_rows <- .active_reference_control_rows(
        control_df = control_df,
        fcs_files_all = fcs_files_all
    )
    if (nrow(active_rows) == 0) {
        return(invisible(TRUE))
    }

    processed_fluors <- if (length(processed_results) > 0) {
        trimws(as.character(vapply(processed_results, function(x) x$fluorophore[[1]], character(1))))
    } else {
        character()
    }
    expected_fluors <- trimws(as.character(active_rows$fluorophore))
    expected_fluors <- expected_fluors[nzchar(expected_fluors)]
    missing <- setdiff(unique(expected_fluors), unique(processed_fluors))
    if (length(missing) == 0) {
        return(invisible(TRUE))
    }

    missing_rows <- active_rows[trimws(as.character(active_rows$fluorophore)) %in% missing, , drop = FALSE]
    details <- paste0(
        as.character(missing_rows$filename),
        " -> ",
        as.character(missing_rows$fluorophore)
    )
    stop(
        paste(
            c(
                "Reference matrix construction did not produce spectra for all mapped non-AF controls.",
                "This usually means one or more controls failed reading, scatter gating, or histogram gating.",
                paste0(" - ", unique(details)),
                "Fix the listed controls or remove them from the control file before continuing."
            ),
            collapse = "\n"
        ),
        call. = FALSE
    )
}

# Computes a 2D scatter gate in FSC-SSC space for a sample.
# Filters out extreme outliers and debris, fits a Gaussian Mixture Model (GMM) to find clusters,
# selects the appropriate cell/bead population component(s), and builds a polygon gating boundary.
# Returns a list containing gated event data, gate coordinates, and scatter parameter bounds.
