#' Gating options
#'
#' Helper constructor for SCC event-selection settings.
#'
#' @param use_scatter_gating Logical; if `TRUE` (default), keep broad FSC/SSC
#'   cleanup and use the intensity-vs-FSC GMM/EM selector for SCC events. If
#'   `FALSE`, use the legacy one-dimensional histogram gate.
#' @param use_af_cosine_scc_selection Logical; if `TRUE` (default), use the
#'   adaptive AF projection/cosine selector for cell SCCs when AF controls are
#'   available. This mode scores events by low AF similarity, residual
#'   target-channel brightness, and target-channel dominance, then keeps the
#'   high-score component. Bead SCCs continue to use the GMM/EM
#'   intensity-vs-FSC selector after broad FSC/SSC cleanup.
#' @param histogram_pct_beads Quantile width for the bead histogram gate.
#' @param histogram_direction_beads Gate direction for beads: `"right"` starts at the median,
#'   `"both"` centers on the median, and `"left"` ends at the median.
#' @param histogram_pct_cells Quantile width for the default cell histogram
#'   gate. When `use_af_cosine_scc_selection = TRUE`, this acts as a
#'   conservative maximum cap for the experimental adaptive AF-score selector
#'   rather than a fixed fraction to keep.
#' @param histogram_direction_cells Gate direction for cells: `"right"` starts at the median,
#'   `"both"` centers on the median, and `"left"` ends at the median.
#' @return Named list with gating settings.
#' @examples
#' opts <- gating_options(histogram_pct_beads = 0.98, histogram_pct_cells = 0.35)
#' str(opts)
#' @export
gating_options <- function(use_scatter_gating = TRUE,
                           use_af_cosine_scc_selection = TRUE,
                           histogram_pct_beads = 0.98,
                           histogram_direction_beads = "right",
                           histogram_pct_cells = 0.35,
                           histogram_direction_cells = "right") {
    list(
        use_scatter_gating = use_scatter_gating,
        use_af_cosine_scc_selection = use_af_cosine_scc_selection,
        histogram_pct_beads = histogram_pct_beads,
        histogram_direction_beads = histogram_direction_beads,
        histogram_pct_cells = histogram_pct_cells,
        histogram_direction_cells = histogram_direction_cells
    )
}

#' Get Standard Fluorophore Filename Patterns
#'
#' Returns pattern dictionaries used for auto-detection from SCC filenames.
#'
#' @return List with `unstained`, `beads`, and `cells` pattern vectors.
#' @examples
#' pats <- get_fluorophore_patterns()
#' names(pats)
#' @export
get_fluorophore_patterns <- function() {
    list(
        unstained = c("US_UT", "Unstained", "unstained", "blank", "AF only", "AF"),
        beads = c(
            # Alexa dyes
            "Alexa Fluor 350", "Alexa Fluor 405", "Alexa Fluor 430", "Alexa Fluor 488", "Alexa Fluor 514",
            "Alexa Fluor 532", "Alexa Fluor 546", "Alexa Fluor 555", "Alexa Fluor 568", "Alexa Fluor 594",
            "Alexa Fluor 610", "Alexa Fluor 633", "Alexa Fluor 647", "Alexa Fluor 660", "Alexa Fluor 680",
            "Alexa Fluor 700", "Alexa Fluor 750", "Alexa Fluor 790",
            "AlexaFluor 594", "AlexaFluor", "Alexa Fluor", "AF594",
            "Alexa 350", "Alexa 405", "Alexa 430", "Alexa 488", "Alexa 514",
            "Alexa 532", "Alexa 546", "Alexa 555", "Alexa 568", "Alexa 594",
            "Alexa 610", "Alexa 633", "Alexa 647", "Alexa 660", "Alexa 680",
            "Alexa 700", "Alexa 750", "Alexa 790",
            # Proteins
            "FITC", "BB515", "GFP", "YFP", "CFP", "mCherry", "DsRed", "tdTomato",
            "PE-Cy5.5", "PE-Cy5", "PE-Cy7", "PE-CF594", "PE-Texas Red",
            "PE-Dazzle 594", "PE-Fire 640", "PE-Fire 700", "PE-Fire 810", "PE",
            "APC-Cy7", "APC-H7", "APC-R700", "APC-Fire 750", "APC-Fire 810", "APC",
            "PerCP-Cy5.5", "PerCP-eFluor 710", "PerCP", "BB700",
            # Brilliant Violet
            "BV421", "BV480", "BV510", "BV570", "BV605", "BV650", "BV661", "BV711", "BV737", "BV750", "BV785",
            # Brilliant UltraViolet
            "BUV395", "BUV496", "BUV563", "BUV615", "BUV661", "BUV737", "BUV805",
            # Others
            "Pacific Blue", "Pacific Orange", "Cy5.5", "Cy5", "Cy7", "Cy3",
            "RB545", "RB613", "RB667", "RB705", "RB744", "RB780",
            "Spark Blue 550", "Spark NIR 685",
            "Super Bright 436", "Super Bright 600", "Super Bright 645", "Super Bright 702", "Super Bright 780",
            "eFluor 450", "eFluor 520", "eFluor 660", "V450", "V500",
            "SB436", "SB600", "SB645", "SB702", "SB780"
        ),
        cells = c(
            "DAPI", "Hoechst", "PI", "7-AAD", "FVD", "FVD450", "FVD506", "FVD520", "FVD780", "FVS", "FixViab",
            "eFluor 506", "eFluor 780", "eFluor 455UV",
            "Zombie Aqua", "Zombie NIR", "Zombie UV", "Zombie Violet",
            "Zombie Green", "Zombie Red", "Zombie Yellow",
            "LIVE/DEAD", "LIVE DEAD", "LIVE/DEAD NIR", "LIVE DEAD NIR",
            "LIVE/DEAD Aqua", "LIVE DEAD AQUA", "LIVE/DEAD Violet", "LIVE DEAD VIOLET",
            "LIVE/DEAD Yellow", "LIVE DEAD YELLOW", "LIVE/DEAD Green", "LIVE DEAD GREEN",
            "LIVE/DEAD Blue", "LIVE DEAD BLUE", "LIVE/DEAD Red", "LIVE DEAD RED",
            "LIVE/DEAD Far Red", "LIVE DEAD FAR RED",
            "Fixable Viability", "Viability", "Calcein", "SYTOX"
        )
    )
}

#' Get Sorted Detectors by Laser and Wavelength
#'
#' @param pd pData from flowFrame parameters
#' @return A list with [[names]] (FL...) and [[labels]] (405nm - 420/10) sorted by laser.
#' @examples
#' pd <- data.frame(
#'   name = c("FSC-A", "V1-A", "B2-A", "R1-A"),
#'   desc = c("FSC-A", "405 nm 450/50", "488 nm 530/30", "640 nm 670/30"),
#'   stringsAsFactors = FALSE
#' )
#' get_sorted_detectors(pd)$names
#' @export
get_sorted_detectors <- function(pd) {
    # 1. Identify spectral detectors
    # Common patterns: FL[0-9]+-A, [A-Z][0-9]+-A, etc.
    # Exclude FSC, SSC, Time
    exclude_patterns <- "^FSC|^SSC|^Time|^Event|^ID"
    matches <- grep(exclude_patterns, pd$name, ignore.case = TRUE, invert = TRUE)

    # Further filter for Area channels if they exist, otherwise use all
    area_matches <- grep("-A$", pd$name[matches])
    if (length(area_matches) > 0) {
        matches <- matches[area_matches]
    }

    fl_pd <- pd[matches, ]

    # 2. Extract descriptions/labels
    # Use 'desc' if it looks like a filter name, otherwise fallback to 'name'
    desc <- as.character(fl_pd$desc)
    names <- trimws(as.character(fl_pd$name))

    # Final labels for plotting
    labels <- ifelse(!is.na(desc) & desc != "" & desc != names, desc, names)

    # 3. Parse laser and wavelength for sorting
    # Try to find laser nm (e.g., "405nm" or "405-")
    laser_nm <- suppressWarnings(as.integer(ifelse(
        grepl("[0-9]{3}\\s*nm", labels),
        gsub(".*?([0-9]{3})\\s*nm.*", "\\1", labels),
        NA_character_
    )))
    # If failed, try name (e.g., V1, B2)
    if (all(is.na(laser_nm))) {
        laser_code <- substr(names, 1, 1)
        laser_nm <- ifelse(laser_code == "U", 355,
            ifelse(laser_code == "V", 405,
                ifelse(laser_code == "B", 488,
                    ifelse(laser_code == "Y" | laser_code == "G", 561,
                        ifelse(laser_code == "R", 640, 999)
                    )
                )
            )
        )
    }

    # Laser order: UV (~355), V (~405), B (~488), YG (~561), R (~640)
    laser_priority <- ifelse(laser_nm < 360, 1, # UV
        ifelse(laser_nm < 420, 2, # V
            ifelse(laser_nm < 500, 3, # B
                ifelse(laser_nm < 600, 4, # YG
                    5
                )
            )
        )
    ) # R

    # Wavelength: look for numbers after the laser name
    label_wo_laser <- sub("^[0-9]{3}\\s*nm\\s*", "", labels)
    wavelength <- suppressWarnings(as.integer(ifelse(
        grepl("[0-9]{3}", label_wo_laser),
        gsub(".*?([0-9]{3}).*", "\\1", label_wo_laser),
        NA_character_
    )))
    if (all(is.na(wavelength))) wavelength <- seq_along(names) # Fallback to index

    # 4. Sort
    ord <- order(laser_priority, wavelength, na.last = TRUE)

    return(list(
        names = names[ord],
        labels = labels[ord],
        laser_nm = laser_nm[ord]
    ))
}
