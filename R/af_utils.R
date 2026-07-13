.is_af_fluorophore <- function(x) {
    out <- trimws(as.character(x))
    out[is.na(out)] <- ""
    grepl("^AF($|_)", out, ignore.case = TRUE)
}

.is_primary_af_fluorophore <- function(x) {
    out <- trimws(as.character(x))
    out[is.na(out)] <- ""
    grepl("^AF($|_[0-9]+$)", out, ignore.case = TRUE)
}

.is_dead_af_fluorophore <- function(x) {
    out <- trimws(as.character(x))
    out[is.na(out)] <- ""
    grepl("^AF_DEAD$", out, ignore.case = TRUE)
}

.is_af_marker <- function(x) {
    out <- trimws(as.character(x))
    out[is.na(out)] <- ""
    grepl("autofluorescence", out, ignore.case = TRUE)
}

.is_af_filename <- function(x) {
    out <- trimws(as.character(x))
    out[is.na(out)] <- ""
    out <- tools::file_path_sans_ext(basename(out))

    grepl(
        "(^|[^[:alnum:]])US(_UT)?([^[:alnum:]]|$)|UNSTAINED|AF[ _-]?ONLY|AUTOFLUORESCENCE|(^|[^[:alnum:]])BLANK([^[:alnum:]]|$)",
        out,
        ignore.case = TRUE
    )
}

.is_dead_af_filename <- function(x) {
    stem <- tools::file_path_sans_ext(basename(trimws(as.character(x))))
    stem[is.na(stem)] <- ""
    stem_norm <- gsub("[^[:alnum:]]+", "", tolower(stem))
    has_af_context <- grepl("unstained|autofluorescence|(^|[^[:alnum:]])af([^[:alnum:]]|$)", stem, ignore.case = TRUE) |
        grepl("unstained|autofluorescence", stem_norm, ignore.case = TRUE)
    has_dead_context <- grepl("dead|viability|livedead|liveanddead", stem_norm, ignore.case = TRUE)

    has_af_context & has_dead_context
}

.is_af_control_row <- function(fluorophore = NULL, marker = NULL, filename = NULL) {
    n <- max(length(fluorophore), length(marker), length(filename), 1L)

    fluor_flag <- if (length(fluorophore) > 0) rep_len(.is_af_fluorophore(fluorophore), n) else rep(FALSE, n)
    marker_flag <- if (length(marker) > 0) rep_len(.is_af_marker(marker), n) else rep(FALSE, n)
    file_flag <- if (length(filename) > 0) rep_len(.is_af_filename(filename), n) else rep(FALSE, n)

    fluor_flag | marker_flag | file_flag
}

.is_primary_af_control_row <- function(fluorophore = NULL, marker = NULL, filename = NULL) {
    n <- max(length(fluorophore), length(marker), length(filename), 1L)

    fluor_text <- if (length(fluorophore) > 0) {
        out <- trimws(as.character(rep_len(fluorophore, n)))
        out[is.na(out)] <- ""
        out
    } else {
        rep("", n)
    }
    fluor_known <- nzchar(fluor_text)
    primary_fluor_flag <- .is_primary_af_fluorophore(fluor_text)
    marker_flag <- if (length(marker) > 0) rep_len(.is_af_marker(marker), n) else rep(FALSE, n)
    file_flag <- if (length(filename) > 0) rep_len(.is_af_filename(filename), n) else rep(FALSE, n)

    primary_fluor_flag | (!fluor_known & (marker_flag | file_flag))
}

.is_dead_af_control_row <- function(fluorophore = NULL, marker = NULL, filename = NULL, control_type = NULL) {
    n <- max(length(fluorophore), length(marker), length(filename), length(control_type), 1L)

    fluor_text <- if (length(fluorophore) > 0) {
        out <- trimws(as.character(rep_len(fluorophore, n)))
        out[is.na(out)] <- ""
        out
    } else {
        rep("", n)
    }
    fluor_flag <- .is_dead_af_fluorophore(fluor_text)
    legacy_internal_flag <- grepl("^AF_INTERNAL$", fluor_text, ignore.case = TRUE)
    marker_flag <- if (length(marker) > 0) rep_len(.is_af_marker(marker), n) else rep(FALSE, n)
    file_flag <- if (length(filename) > 0) rep_len(.is_dead_af_filename(filename), n) else rep(FALSE, n)
    type_text <- if (length(control_type) > 0) {
        out <- tolower(trimws(as.character(rep_len(control_type, n))))
        out[is.na(out)] <- ""
        out
    } else {
        rep("", n)
    }
    is_bead <- type_text == "beads"

    (fluor_flag | (legacy_internal_flag & file_flag) | (marker_flag & file_flag) | file_flag) & !is_bead
}

.canonicalize_primary_af_labels <- function(control_df) {
    if (is.null(control_df) || !is.data.frame(control_df) || nrow(control_df) == 0L ||
            !("fluorophore" %in% colnames(control_df))) {
        return(control_df)
    }

    n <- nrow(control_df)
    fluor <- trimws(as.character(control_df$fluorophore))
    fluor[is.na(fluor)] <- ""
    marker <- if ("marker" %in% colnames(control_df)) {
        out <- trimws(as.character(control_df$marker))
        out[is.na(out)] <- ""
        out
    } else {
        rep("", n)
    }
    filename <- if ("filename" %in% colnames(control_df)) {
        out <- trimws(as.character(control_df$filename))
        out[is.na(out)] <- ""
        out
    } else {
        rep("", n)
    }
    control_type <- if ("control.type" %in% colnames(control_df)) {
        out <- tolower(trimws(as.character(control_df$control.type)))
        out[is.na(out)] <- ""
        out
    } else {
        rep("", n)
    }

    af_context <- .is_af_marker(marker) | .is_af_filename(filename)
    legacy_internal <- grepl("^AF_INTERNAL$", fluor, ignore.case = TRUE)
    candidates <- .is_primary_af_fluorophore(fluor) |
        (legacy_internal & af_context) |
        (!nzchar(fluor) & af_context)
    special_rows <- .is_dead_af_control_row(
        fluorophore = fluor,
        marker = marker,
        filename = filename,
        control_type = control_type
    ) | control_type == "beads" |
        grepl("^AF_BEADS$", fluor, ignore.case = TRUE)
    candidates <- candidates & !special_rows

    idx <- which(candidates)
    if (length(idx) > 0L) {
        control_df$fluorophore[idx] <- c(
            "AF",
            if (length(idx) > 1L) paste0("AF_", seq.int(2L, length(idx))) else NULL
        )
        if ("marker" %in% colnames(control_df)) {
            control_df$marker[idx] <- "Autofluorescence"
        }
        if ("control.type" %in% colnames(control_df)) {
            control_df$control.type[idx] <- "cells"
        }
    }

    control_df
}
