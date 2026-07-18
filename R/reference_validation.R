# Internal validation for raw reference-control FCS data.
.validate_reference_raw_data <- function(raw_data, sn, detector_names = NULL) {
    if (!is.matrix(raw_data) && !is.data.frame(raw_data)) {
        .spectreasy_console_step("Skip SCC", paste0(sn, " does not contain a rectangular event matrix"))
        return(FALSE)
    }
    if (nrow(raw_data) == 0L || ncol(raw_data) == 0L) {
        .spectreasy_console_step("Skip SCC", paste0(sn, " contains no events or measured channels"))
        return(FALSE)
    }

    checked_data <- raw_data
    if (!is.null(detector_names)) {
        scatter_channels <- unlist(.get_primary_scatter_channels(colnames(raw_data)), use.names = FALSE)
        scatter_channels <- scatter_channels[!is.na(scatter_channels)]
        checked_names <- intersect(unique(c(as.character(detector_names), scatter_channels)), colnames(raw_data))
        if (length(checked_names) > 0L) {
            checked_data <- raw_data[, checked_names, drop = FALSE]
        }
    }
    if (any(is.infinite(checked_data))) {
        .spectreasy_console_step("Skip SCC", paste0(sn, " has infinite values"))
        return(FALSE)
    }
    na_prop <- sum(is.na(checked_data)) / length(checked_data)
    if (na_prop > 0.1) {
        .spectreasy_console_step("Skip SCC", paste0(sn, " has too many NAs (", round(na_prop * 100, 1), "%)"))
        return(FALSE)
    }
    max_val <- max(checked_data, na.rm = TRUE)
    if (max_val > 1e9) {
        .spectreasy_console_step("Skip SCC", paste0(sn, " has extreme values (max > 1e9)"))
        return(FALSE)
    }
    TRUE
}

.validate_reference_detector_consistency <- function(fcs_files, detector_names) {
    if (length(fcs_files) == 0 || length(detector_names) == 0) {
        return(invisible(TRUE))
    }

    mismatches <- character()
    for (fcs_file in fcs_files) {
        ff <- tryCatch(
            .spectreasy_read_fcs(fcs_file),
            error = function(e) {
                stop("Could not read FCS file while checking detector consistency: ", fcs_file, call. = FALSE)
            }
        )
        pd <- flowCore::pData(flowCore::parameters(ff))
        current_detectors <- tryCatch(
            get_sorted_detectors(pd)$names,
            error = function(e) character()
        )
        if (length(current_detectors) == 0) {
            current_detectors <- intersect(detector_names, colnames(flowCore::exprs(ff)))
        }
        missing <- setdiff(detector_names, current_detectors)
        extra <- setdiff(current_detectors, detector_names)
        if (length(missing) > 0 || length(extra) > 0) {
            mismatches <- c(
                mismatches,
                paste0(
                    basename(fcs_file),
                    if (length(missing) > 0) paste0(" is missing: ", paste(missing, collapse = ", ")) else "",
                    if (length(missing) > 0 && length(extra) > 0) "; " else "",
                    if (length(extra) > 0) paste0("has extra detectors: ", paste(extra, collapse = ", ")) else ""
                )
            )
        }
    }

    if (length(mismatches) > 0) {
        stop(
            paste(
                c(
                    "Detector set mismatch across SCC/AF files.",
                    "All files used to build a reference matrix must contain the detectors found in the first SCC file.",
                    paste0(" - ", mismatches)
                ),
                collapse = "\n"
            ),
            call. = FALSE
        )
    }

    invisible(TRUE)
}
