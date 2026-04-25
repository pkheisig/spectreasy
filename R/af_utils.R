.is_af_fluorophore <- function(x) {
    out <- trimws(as.character(x))
    out[is.na(out)] <- ""
    grepl("^AF($|_)", out, ignore.case = TRUE)
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
        "(^|[^[:alnum:]])US_UT([^[:alnum:]]|$)|UNSTAINED|AF[ _-]?ONLY|AUTOFLUORESCENCE|(^|[^[:alnum:]])BLANK([^[:alnum:]]|$)",
        out,
        ignore.case = TRUE
    )
}

.is_af_control_row <- function(fluorophore = NULL, marker = NULL, filename = NULL) {
    n <- max(length(fluorophore), length(marker), length(filename), 1L)

    fluor_flag <- if (length(fluorophore) > 0) rep_len(.is_af_fluorophore(fluorophore), n) else rep(FALSE, n)
    marker_flag <- if (length(marker) > 0) rep_len(.is_af_marker(marker), n) else rep(FALSE, n)
    file_flag <- if (length(filename) > 0) rep_len(.is_af_filename(filename), n) else rep(FALSE, n)

    fluor_flag | marker_flag | file_flag
}
