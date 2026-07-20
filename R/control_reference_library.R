.empty_control_file_df <- function() {
    data.frame(
        filename = character(),
        fluorophore = character(),
        marker = character(),
        channel = character(),
        control.type = character(),
        universal.negative = character(),
        is.viability = character(),
        stringsAsFactors = FALSE
    )
}

.control_file_normalize_token <- function(x) {
    out <- gsub("[^a-z0-9]+", "", tolower(trimws(as.character(x))))
    out[is.na(out)] <- ""
    out
}

.control_file_canonicalize_channel <- function(x) {
    out <- toupper(gsub("\\s+", "", trimws(as.character(x))))
    out <- gsub("([A-Z]+)-([0-9])", "\\1\\2", out, perl = TRUE)
    out[is.na(out)] <- ""
    out
}

.control_file_split_semicolon <- function(x) {
    if (is.null(x) || is.na(x)) return(character())
    vals <- trimws(unlist(strsplit(as.character(x), ";", fixed = TRUE)))
    vals[nzchar(vals)]
}

.control_file_extdata_file <- function(filename) {
    .spectreasy_extdata_file(filename)
}

.control_file_preferred_fluors <- function() {
    c(
        "FITC", "BB515", "Alexa Fluor 488",
        "PE", "PE-CF594", "PE-Dazzle 594", "PE-Fire 640", "PE-Fire 700", "PE-Fire 810",
        "PE-Cy5", "PE-Cy5.5", "PE-Cy7",
        "APC", "Alexa Fluor 647", "Alexa Fluor 700", "APC-R700", "APC-Cy7", "APC-H7", "APC-Fire 750", "APC-Fire 810",
        "PerCP", "PerCP-Cy5.5", "BB700",
        "BV421", "BV480", "BV510", "BV570", "BV605", "BV650", "BV661", "BV711", "BV737", "BV750", "BV785",
        "BUV395", "BUV496", "BUV563", "BUV615", "BUV661", "BUV737", "BUV805",
        "Alexa Fluor 350", "Alexa Fluor 405", "Alexa Fluor 430", "Alexa Fluor 514",
        "Alexa Fluor 532", "Alexa Fluor 546", "Alexa Fluor 555", "Alexa Fluor 568",
        "Alexa Fluor 594", "Alexa Fluor 610", "Alexa Fluor 633", "Alexa Fluor 660",
        "Alexa Fluor 680", "Alexa Fluor 750", "Alexa Fluor 790",
        "RB545", "RB613", "RB667", "RB705", "RB744", "RB780",
        "Spark Blue 550", "Spark NIR 685",
        "Super Bright 436", "Super Bright 600", "Super Bright 645", "Super Bright 702", "Super Bright 780"
    )
}

.control_file_choose_preferred_fluor <- function(candidates, preferred_fluors = .control_file_preferred_fluors()) {
    candidates <- unique(trimws(as.character(candidates)))
    candidates <- candidates[nzchar(candidates)]
    if (length(candidates) == 0) return("")
    pref_idx <- match(tolower(candidates), tolower(preferred_fluors))
    if (any(!is.na(pref_idx))) {
        return(candidates[which.min(ifelse(is.na(pref_idx), Inf, pref_idx))])
    }
    sort(candidates)[1]
}

.control_file_generated_marker_map <- function(max_cd = 371L) {
    marker_names <- c(
        paste0("CD", seq_len(max_cd)),
        paste0("CD1", letters[1:5]),
        "CD45RA", "CD45RO"
    )
    marker_names <- unique(marker_names)
    keys <- .control_file_normalize_token(marker_names)
    keys <- keys[nzchar(keys)]
    stats::setNames(marker_names[seq_along(keys)], keys)
}

.control_file_read_dictionary <- function(path) {
    if (!nzchar(path)) return(NULL)
    tryCatch(
        utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
        error = function(e) NULL
    )
}

.control_file_load_fluor_names <- function(path) {
    data <- .control_file_read_dictionary(path)
    result <- character()
    if (is.null(data) || nrow(data) == 0 || !"fluorophore" %in% colnames(data)) return(result)
    for (i in seq_len(nrow(data))) {
        canonical <- trimws(as.character(data$fluorophore[i]))
        if (!nzchar(canonical)) next
        aliases <- if ("aliases" %in% colnames(data)) data$aliases[i] else ""
        keys <- .control_file_normalize_token(unique(c(canonical, .control_file_split_semicolon(aliases))))
        keys <- keys[nzchar(keys) & !keys %in% names(result)]
        result[keys] <- canonical
    }
    result
}

.control_file_add_peak_channel <- function(reference, fluorophore, peak_channel, overwrite = FALSE) {
    fluorophore <- trimws(as.character(fluorophore)[1])
    peak_channel <- .control_file_canonicalize_channel(peak_channel)[1]
    if (!nzchar(fluorophore) || is.na(peak_channel) || !nzchar(peak_channel)) return(reference)
    canonical <- reference$name_map[.control_file_normalize_token(fluorophore)]
    if (is.na(canonical) || !nzchar(canonical)) canonical <- fluorophore
    keys <- unique(c(
        .control_file_normalize_token(fluorophore),
        .control_file_normalize_token(canonical),
        names(reference$name_map)[reference$name_map == canonical]
    ))
    keys <- keys[nzchar(keys)]
    if (!isTRUE(overwrite)) keys <- keys[!keys %in% names(reference$fluor_peak_channel_map)]
    reference$fluor_peak_channel_map[keys] <- peak_channel
    reference
}

.control_file_filter_channels_by_cytometer <- function(data, cytometer_id) {
    if (identical(cytometer_id, "auto")) return(data)
    ids <- vapply(
        data$cytometer,
        .resolve_cytometer_id,
        character(1),
        allow_auto = FALSE,
        unknown_as_auto = FALSE
    )
    data[ids == cytometer_id, , drop = FALSE]
}

.control_file_load_channels <- function(reference, path, cytometer_id, preferred_fluors) {
    data <- .control_file_read_dictionary(path)
    required <- c("fluorophore", "cytometer", "channel")
    if (is.null(data) || nrow(data) == 0 || !all(required %in% colnames(data))) return(reference)
    data <- .control_file_filter_channels_by_cytometer(data, cytometer_id)
    for (i in seq_len(nrow(data))) {
        canonical <- trimws(as.character(data$fluorophore[i]))
        if (!nzchar(canonical)) next
        channels <- .control_file_canonicalize_channel(.control_file_split_semicolon(data$channel[i]))
        channels <- unique(c(channels, gsub("-A$", "", channels), paste0(gsub("-A$", "", channels), "-A")))
        channels <- channels[nzchar(channels)]
        peak <- channels[grepl("-A$", channels)][1]
        if (is.na(peak) || !nzchar(peak)) peak <- channels[1]
        reference <- .control_file_add_peak_channel(reference, canonical, peak)
        for (channel in channels) {
            reference$channel_map[channel] <- if (!channel %in% names(reference$channel_map)) {
                canonical
            } else {
                .control_file_choose_preferred_fluor(
                    c(reference$channel_map[channel], canonical),
                    preferred_fluors = preferred_fluors
                )
            }
        }
    }
    reference
}

.control_file_add_library_peaks <- function(reference, cytometer_id) {
    if (identical(cytometer_id, "auto") || !cytometer_id %in% names(.spectral_library_file_map())) return(reference)
    library <- tryCatch(
        .read_spectral_library_matrix(cytometer_id, normalize = FALSE),
        error = function(e) NULL
    )
    if (is.null(library) || nrow(library) == 0 || ncol(library) == 0) return(reference)
    peak_channels <- colnames(library)[max.col(library, ties.method = "first")]
    for (i in seq_len(nrow(library))) {
        reference <- .control_file_add_peak_channel(reference, rownames(library)[i], peak_channels[i])
    }
    reference
}

.control_file_load_markers <- function(reference, path) {
    data <- .control_file_read_dictionary(path)
    if (!is.null(data) && nrow(data) > 0 && "marker" %in% colnames(data)) {
        for (i in seq_len(nrow(data))) {
            marker <- trimws(as.character(data$marker[i]))
            if (!nzchar(marker)) next
            reference$marker_names <- c(reference$marker_names, marker)
            aliases <- if ("aliases" %in% colnames(data)) data$aliases[i] else ""
            keys <- .control_file_normalize_token(unique(c(marker, .control_file_split_semicolon(aliases))))
            keys <- keys[nzchar(keys) & !keys %in% names(reference$marker_map)]
            reference$marker_map[keys] <- marker
        }
    }
    generated <- .control_file_generated_marker_map()
    reference$marker_names <- unique(c(reference$marker_names, unname(generated)))
    generated <- generated[!names(generated) %in% names(reference$marker_map)]
    reference$marker_map <- c(reference$marker_map, generated)
    reference
}

.load_control_file_shipped_reference <- function(cytometer_name, preferred_fluors = .control_file_preferred_fluors()) {
    reference <- list(
        name_map = .control_file_load_fluor_names(
            .control_file_extdata_file("fluorophore_dictionary.csv")
        ),
        channel_map = character(),
        fluor_peak_channel_map = character(),
        marker_map = character(),
        marker_names = character()
    )
    cytometer_id <- .resolve_cytometer_id(
        cytometer_name,
        allow_auto = TRUE,
        unknown_as_auto = TRUE
    )
    reference <- .control_file_load_channels(
        reference,
        .control_file_extdata_file("fluorophore_channel_dictionary.csv"),
        cytometer_id,
        preferred_fluors
    )
    reference <- .control_file_add_library_peaks(reference, cytometer_id)
    .control_file_load_markers(
        reference,
        .control_file_extdata_file("marker_dictionary.csv")
    )
}

.prepare_control_file_reference <- function(cytometer) {
    patterns_list <- get_fluorophore_patterns()
    specific_dyes <- c(patterns_list$beads, patterns_list$cells)
    generic_terms <- c("LIVE/DEAD", "Fixable Viability", "Viability", "FVD")
    specific_dyes <- setdiff(specific_dyes, generic_terms)
    specific_dyes <- specific_dyes[order(-nchar(specific_dyes))]
    fluor_patterns <- c(specific_dyes, generic_terms)
    preferred_fluors <- .control_file_preferred_fluors()
    shipped_ref <- .load_control_file_shipped_reference(cytometer, preferred_fluors = preferred_fluors)

    list(
        fluor_patterns = fluor_patterns,
        sample_patterns = patterns_list,
        fluor_name_map = shipped_ref$name_map,
        fluor_channel_map = shipped_ref$channel_map,
        fluor_peak_channel_map = shipped_ref$fluor_peak_channel_map,
        marker_name_map = shipped_ref$marker_map,
        marker_names = shipped_ref$marker_names,
        preferred_fluors = preferred_fluors
    )
}
