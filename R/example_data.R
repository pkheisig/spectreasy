#' Download and cache the spectreasy example FCS files
#'
#' Downloads the example-data zip archive from the latest GitHub release,
#' unzips it automatically, and returns the local paths to the extracted
#' `sample/` and `scc/` folders.
#'
#' The downloaded archive is cached under [tools::R_user_dir()] so it only needs
#' to be fetched once unless `force = TRUE`.
#'
#' @param asset Optional URL to a zip archive containing `sample/` and `scc/`
#'   folders. Defaults to the latest `spectreasy` GitHub release asset. A local
#'   zip file path can also be supplied for testing or offline reuse.
#' @param cache_dir Directory used for the downloaded zip and extracted files.
#'   Defaults to a package-specific user cache directory.
#' @param force Logical; if `TRUE`, redownload and re-extract the archive even if
#'   a cached copy is already available.
#' @param quiet Logical; if `TRUE`, suppress progress messages where possible.
#'
#' @return A named list with elements `root_dir`, `zip_file`, `sample_dir`,
#'   `sample_files`, and `scc_dir`.
#' @export
#'
#' @examples
#' \dontrun{
#' paths <- spectreasy_example_data()
#' list.files(paths$scc_dir)
#' list.files(paths$sample_dir)
#' }
spectreasy_example_data <- function(
    asset = NULL,
    cache_dir = file.path(tools::R_user_dir("spectreasy", which = "cache"), "example-data"),
    force = FALSE,
    quiet = FALSE
) {
    if (is.null(asset)) {
        asset <- .spectreasy_default_example_asset()
    }
    asset <- as.character(asset)[1]
    cache_dir <- as.character(cache_dir)[1]

    if (is.na(asset) || !nzchar(trimws(asset))) {
        stop("asset must be a non-empty URL or local zip file path.", call. = FALSE)
    }
    if (is.na(cache_dir) || !nzchar(trimws(cache_dir))) {
        stop("cache_dir must be a non-empty directory path.", call. = FALSE)
    }

    asset_name <- .spectreasy_example_asset_name(asset)
    extract_root <- file.path(cache_dir, tools::file_path_sans_ext(asset_name))
    zip_file <- file.path(cache_dir, asset_name)

    if (isTRUE(force) && dir.exists(extract_root)) {
        unlink(extract_root, recursive = TRUE, force = TRUE)
    }
    if (isTRUE(force) && file.exists(zip_file)) {
        unlink(zip_file, force = TRUE)
    }

    sample_dir <- .spectreasy_find_example_subdir(extract_root, "sample")
    scc_dir <- .spectreasy_find_example_subdir(extract_root, "scc")

    needs_extract <- is.null(sample_dir) || is.null(scc_dir)
    if (needs_extract) {
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
        if (!file.exists(zip_file)) {
            .spectreasy_download_example_asset(asset = asset, destfile = zip_file, quiet = quiet)
        }

        dir.create(extract_root, recursive = TRUE, showWarnings = FALSE)
        utils::unzip(zip_file, exdir = extract_root)

        sample_dir <- .spectreasy_find_example_subdir(extract_root, "sample")
        scc_dir <- .spectreasy_find_example_subdir(extract_root, "scc")
        if (is.null(sample_dir) || is.null(scc_dir)) {
            stop(
                "Example data archive did not contain the expected 'sample/' and 'scc/' folders: ",
                asset,
                call. = FALSE
            )
        }
    }

    sample_files <- list.files(sample_dir, pattern = "\\.[Ff][Cc][Ss]$", full.names = TRUE)

    list(
        root_dir = extract_root,
        zip_file = zip_file,
        sample_dir = sample_dir,
        sample_files = sample_files,
        scc_dir = scc_dir
    )
}

.spectreasy_default_example_asset <- function() {
    getOption(
        "spectreasy.example_data_asset",
        Sys.getenv(
            "SPECTREASY_EXAMPLE_DATA_ASSET",
            unset = "https://github.com/pkheisig/spectreasy/releases/latest/download/spectreasy-example-data.zip"
        )
    )
}

.spectreasy_example_asset_name <- function(asset) {
    path_part <- strsplit(as.character(asset)[1], "?", fixed = TRUE)[[1]][1]
    basename(path_part)
}

.spectreasy_download_example_asset <- function(asset, destfile, quiet = FALSE) {
    if (file.exists(asset)) {
        ok <- file.copy(asset, destfile, overwrite = TRUE)
        if (!ok) {
            stop("Could not copy local example data archive: ", asset, call. = FALSE)
        }
        return(invisible(destfile))
    }

    status <- tryCatch(
        utils::download.file(asset, destfile = destfile, mode = "wb", quiet = quiet),
        error = function(e) e
    )
    if (inherits(status, "error") || !file.exists(destfile)) {
        stop(
            "Could not download spectreasy example data from: ", asset, "\n",
            "Upload the release asset or supply a local zip file path via asset = ...",
            call. = FALSE
        )
    }

    invisible(destfile)
}

.spectreasy_find_example_subdir <- function(root_dir, target_name) {
    if (!dir.exists(root_dir)) {
        return(NULL)
    }

    direct <- file.path(root_dir, target_name)
    if (dir.exists(direct)) {
        return(normalizePath(direct, mustWork = TRUE))
    }

    all_dirs <- list.dirs(root_dir, recursive = TRUE, full.names = TRUE)
    hits <- all_dirs[basename(all_dirs) == target_name]
    if (length(hits) == 0) {
        return(NULL)
    }

    hits <- hits[order(nchar(hits), hits)]
    normalizePath(hits[[1]], mustWork = TRUE)
}
