.ai_qc_profile_dirs <- function(project_dir = NULL) {
    built_in <- system.file("extdata", "qc_reference_profiles", package = "spectreasy")
    local <- if (!is.null(project_dir)) file.path(project_dir, ".spectreasy", "qc_reference_profiles") else character()
    unique(c(local, built_in[nzchar(built_in)]))
}

.ai_qc_profile_compatibility <- function(profile, context) {
    fields <- c("cytometer", "detector_schema_hash", "method", "af_configuration")
    details <- lapply(fields, function(field) {
        expected <- profile$strata[[field]] %||% NULL
        actual <- context[[field]] %||% NULL
        compatible <- is.null(expected) || is.null(actual) || identical(as.character(expected), as.character(actual))
        list(field = field, expected = expected, actual = actual, compatible = compatible)
    })
    names(details) <- fields
    panel <- suppressWarnings(as.numeric(context$panel_size %||% NA_real_))
    range <- suppressWarnings(as.numeric(profile$strata$panel_range %||% c(NA, NA)))
    panel_ok <- length(range) != 2L || any(!is.finite(range)) || !is.finite(panel) || (panel >= range[1] && panel <= range[2])
    list(
        compatible = all(vapply(details, `[[`, logical(1), "compatible")) && panel_ok,
        details = details,
        panel_range = list(expected = range, actual = panel, compatible = panel_ok)
    )
}

.validate_qc_reference_profile <- function(profile) {
    if (!is.list(profile)) stop("QC reference profile must be a list.", call. = FALSE)
    required <- c("name", "version", "schema_version", "strata", "metrics", "cohort")
    missing <- setdiff(required, names(profile))
    if (length(missing)) stop("QC reference profile is missing: ", paste(missing, collapse = ", "), call. = FALSE)
    class(profile) <- c("spectreasy_qc_reference_profile", "list")
    profile
}

#' List Spectreasy QC reference profiles
#'
#' @param project_dir Optional project directory containing
#'   `.spectreasy/qc_reference_profiles`.
#' @return A data frame describing available profiles.
#' @export
list_qc_reference_profiles <- function(project_dir = NULL) {
    files <- unlist(lapply(.ai_qc_profile_dirs(project_dir), function(path) {
        if (!dir.exists(path)) return(character())
        list.files(path, pattern = "\\.(json|ya?ml|rds)$", full.names = TRUE, ignore.case = TRUE)
    }), use.names = FALSE)
    if (!length(files)) {
        return(data.frame(name = character(), version = character(), path = character(), cohort_n = integer(), stringsAsFactors = FALSE))
    }
    rows <- lapply(files, function(path) {
        profile <- tryCatch(load_qc_reference_profile(path), error = function(e) NULL)
        if (is.null(profile)) return(NULL)
        data.frame(name = profile$name, version = profile$version, path = path, cohort_n = as.integer(profile$cohort$n_clean %||% 0L), stringsAsFactors = FALSE)
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (!length(rows)) return(data.frame(name = character(), version = character(), path = character(), cohort_n = integer(), stringsAsFactors = FALSE))
    do.call(rbind, rows)
}

#' Load a Spectreasy QC reference profile
#'
#' @param reference Profile path, profile name, profile object, `"auto"`, or
#'   `"none"`.
#' @param project_dir Optional project directory used for profile discovery.
#' @param context Optional matching context used when `reference = "auto"`.
#' @return A `spectreasy_qc_reference_profile`, or `NULL` for no profile.
#' @export
load_qc_reference_profile <- function(reference = "auto", project_dir = NULL, context = list()) {
    if (inherits(reference, "spectreasy_qc_reference_profile")) return(reference)
    if (is.list(reference)) return(.validate_qc_reference_profile(reference))
    reference <- as.character(reference)[1]
    if (tolower(reference) %in% c("none", "off", "disabled")) return(NULL)
    candidates <- list_qc_reference_profiles(project_dir)
    if (identical(tolower(reference), "auto")) {
        if (!nrow(candidates)) return(NULL)
        loaded <- lapply(candidates$path, function(path) tryCatch(load_qc_reference_profile(path), error = function(e) NULL))
        loaded <- loaded[!vapply(loaded, is.null, logical(1))]
        matches <- Filter(function(profile) .ai_qc_profile_compatibility(profile, context)$compatible, loaded)
        if (!length(matches)) return(NULL)
        specificity <- vapply(matches, function(profile) sum(vapply(profile$strata, function(x) !is.null(x) && length(x), logical(1))), numeric(1))
        return(matches[[which.max(specificity)]])
    }
    path <- reference
    if (!file.exists(path) && nrow(candidates)) {
        hit <- which(tolower(candidates$name) == tolower(reference))
        if (length(hit)) path <- candidates$path[hit[1]]
    }
    if (!file.exists(path)) stop("QC reference profile not found: ", reference, call. = FALSE)
    ext <- tolower(tools::file_ext(path))
    profile <- switch(ext,
        rds = readRDS(path),
        json = jsonlite::fromJSON(path, simplifyVector = FALSE),
        yml = yaml::read_yaml(path), yaml = yaml::read_yaml(path),
        stop("Unsupported QC reference profile format: ", ext, call. = FALSE)
    )
    profile <- .validate_qc_reference_profile(profile)
    profile$path <- normalizePath(path, mustWork = FALSE)
    profile
}
