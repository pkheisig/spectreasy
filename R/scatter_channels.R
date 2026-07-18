# Internal scatter-channel aliases shared by workflow and GUI code.
.scatter_channel_aliases <- function(type = c("fsc", "ssc")) {
    type <- match.arg(type)
    if (identical(type, "fsc")) {
        c("FSC", "FS", "Forward Scatter")
    } else {
        c("SSC", "BSSC", "VSSC", "USSC", "YSSC", "RSSC", "IRSSC", "SS", "Side Scatter")
    }
}

.get_scatter_channel <- function(col_names, type = c("fsc", "ssc"), suffix = "A") {
    type <- match.arg(type)
    col_names <- as.character(col_names)
    normalized <- toupper(gsub("[^A-Z0-9]", "", col_names))
    prefixes <- toupper(gsub("[^A-Z0-9]", "", .scatter_channel_aliases(type)))
    suffix <- toupper(gsub("[^A-Z0-9]", "", as.character(suffix)[1]))
    is_match <- Reduce(`|`, lapply(prefixes, function(prefix) startsWith(normalized, prefix)))

    for (prefix in prefixes) {
        exact <- col_names[normalized == paste0(prefix, suffix)]
        if (length(exact) > 0L) return(exact[[1]])
    }

    suffixed <- col_names[is_match & grepl(paste0("-", suffix, "$"), col_names, ignore.case = TRUE)]
    if (length(suffixed) > 0L) return(suffixed[[1]])

    any_match <- col_names[is_match]
    if (length(any_match) > 0L) return(any_match[[1]])
    NA_character_
}

.get_scatter_channels <- function(col_names, suffixes = c("A", "H", "W")) {
    col_names <- as.character(col_names)
    normalized <- toupper(gsub("[^A-Z0-9]", "", col_names))
    prefixes <- toupper(gsub(
        "[^A-Z0-9]", "", c(.scatter_channel_aliases("fsc"), .scatter_channel_aliases("ssc"))
    ))
    is_scatter <- Reduce(`|`, lapply(prefixes, function(prefix) startsWith(normalized, prefix)))
    suffix_pattern <- paste0("-(", paste(toupper(suffixes), collapse = "|"), ")$")
    unique(col_names[is_scatter & grepl(suffix_pattern, col_names, ignore.case = TRUE)])
}

.get_primary_scatter_channels <- function(col_names) {
    list(
        fsc = .get_scatter_channel(col_names, "fsc", "A"),
        ssc = .get_scatter_channel(col_names, "ssc", "A")
    )
}
