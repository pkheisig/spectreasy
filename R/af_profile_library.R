# Autofluorescence profile extraction and global profile library helpers.

.af_profile_name_regex <- "^[A-Za-z0-9._-]+$"

.validate_af_profile_name <- function(name) {
    if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(trimws(name))) {
        stop("name must be a single non-empty character string.", call. = FALSE)
    }
    name <- trimws(name)
    if (!grepl(.af_profile_name_regex, name)) {
        stop(
            "AF profile names may only contain letters, numbers, '.', '_', and '-'.",
            call. = FALSE
        )
    }
    name
}

.af_profile_dir_path <- function(create = TRUE) {
    dir_path <- getOption("spectreasy.af_profile_dir", NULL)
    if (is.null(dir_path) || !is.character(dir_path) || length(dir_path) != 1L || is.na(dir_path) || !nzchar(dir_path)) {
        dir_path <- file.path(tools::R_user_dir("spectreasy", which = "data"), "af_profiles")
    }
    dir_path <- normalizePath(dir_path, mustWork = FALSE)
    if (isTRUE(create)) {
        dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    }
    dir_path
}

.af_profile_file <- function(name, create_dir = TRUE) {
    name <- .validate_af_profile_name(name)
    file.path(.af_profile_dir_path(create = create_dir), paste0(name, ".rds"))
}

.is_af_profile_object <- function(x) {
    inherits(x, "spectreasy_af_profile") || (is.list(x) && "profile" %in% names(x))
}

.coerce_af_profile_matrix <- function(x, arg_name = "x") {
    if (.is_af_profile_object(x)) {
        x <- x$profile
    }
    M <- .as_reference_matrix(x, arg_name = arg_name)
    af_rows <- grepl("^AF($|_)", rownames(M), ignore.case = TRUE)
    if (!any(af_rows)) {
        stop(arg_name, " does not contain AF rows named 'AF', 'AF_2', 'AF_3', ...", call. = FALSE)
    }
    M[af_rows, , drop = FALSE]
}

.build_af_profile_plot <- function(profile, pd = NULL) {
    plot_spectra(profile, pd = pd, output_file = NULL)
}

.new_af_profile_object <- function(name = NULL,
                                   profile,
                                   plot = NULL,
                                   source = NULL,
                                   extraction = NULL,
                                   created = Sys.time()) {
    profile <- .coerce_af_profile_matrix(profile, arg_name = "profile")
    if (is.null(plot)) {
        plot <- .build_af_profile_plot(profile)
    }
    out <- list(
        name = name,
        profile = profile,
        plot = plot,
        source = source,
        extraction = extraction,
        created = created,
        spectreasy_version = as.character(utils::packageVersion("spectreasy"))
    )
    class(out) <- c("spectreasy_af_profile", "list")
    out
}

#' Get the global AF profile directory
#'
#' @param create Logical; create the directory if it does not exist.
#'
#' @return Invisibly returns the AF profile directory path.
#' @export
af_profile_dir <- function(create = TRUE) {
    dir_path <- .af_profile_dir_path(create = create)
    message("spectreasy AF profiles are stored at:\n", dir_path)
    invisible(dir_path)
}

#' Extract an AF profile from one unstained FCS file
#'
#' @param fcs_file Path to an unstained/autofluorescence `.fcs` file.
#' @param af_n_bands Number of AF bands to extract. The default, `"auto"`,
#'   builds a SOM bank up to `af_auto_max_bands` nodes plus a mean AF row.
#' @param af_max_cells Maximum number of scatter-gated AF events used.
#' @param af_auto_max_bands Maximum SOM nodes that `"auto"` may create before
#'   prepending the mean AF row.
#' @param seed Optional integer seed for deterministic subsampling/clustering.
#' @param show_plot Logical; print the AF spectra plot after extraction.
#' @param verbose Logical; print progress updates while extracting.
#'
#' @return A `spectreasy_af_profile` object containing `$profile` and `$plot`.
#' @export
extract_af_profile <- function(fcs_file,
                               af_n_bands = "auto",
                               af_max_cells = 50000,
                               af_auto_max_bands = 100,
                               seed = NULL,
                               show_plot = TRUE,
                               verbose = TRUE) {
    if (!is.character(fcs_file) || length(fcs_file) != 1L || is.na(fcs_file) || !nzchar(fcs_file)) {
        stop("fcs_file must be a single path to an FCS file.", call. = FALSE)
    }
    if (!file.exists(fcs_file)) {
        stop("fcs_file not found: ", fcs_file, call. = FALSE)
    }

    if (isTRUE(verbose)) {
        message("Extracting AF profile from: ", basename(fcs_file))
    }
    af_args <- .validate_build_reference_af_args(
        af_n_bands = af_n_bands,
        af_max_cells = af_max_cells,
        af_bands_per_file = 1,
        af_auto_max_bands = af_auto_max_bands
    )
    .with_optional_seed(seed)

    if (isTRUE(verbose)) {
        message("  Detecting spectral detectors...")
    }
    metadata <- .prepare_reference_detector_info(fcs_file)
    if (isTRUE(verbose)) {
        message("  Found ", length(metadata$detector_names), " spectral detector(s).")
        message("  Scatter-gating AF events...")
    }
    config <- list(
        outlier_percentile = 0.02,
        debris_percentile = 0.08,
        subsample_n = 5000,
        max_clusters = 10,
        min_cluster_proportion = 0.03,
        gate_contour_beads = 0.95,
        gate_contour_cells = 0.90,
        bead_gate_scale = 1.3
    )
    gated <- .extract_reference_af_gated_events(
        fcs_file = fcs_file,
        detector_names = metadata$detector_names,
        config = config
    )
    if (is.null(gated) || is.null(gated$events) || nrow(gated$events) == 0) {
        stop("Could not extract scatter-gated AF events from: ", fcs_file, call. = FALSE)
    }

    if (isTRUE(verbose)) {
        message(
            "  Building SOM AF bank",
            if (identical(af_args$af_n_bands, "auto")) paste0(" with ", af_args$af_auto_max_bands, " nodes plus the mean AF row") else paste0(" with ", af_args$af_n_bands, " node(s) plus the mean AF row"),
            "..."
        )
    }
    af_profiles <- withCallingHandlers(
        .extract_reference_af_profiles(
            detector_names = metadata$detector_names,
            n_bands = af_args$af_n_bands,
            max_cells = af_args$af_max_cells,
            af_events = gated$events,
            auto_max_bands = af_args$af_auto_max_bands
        ),
        warning = function(w) {
            if (grepl("Quick-TRANSfer stage steps exceeded maximum", conditionMessage(w), fixed = TRUE)) {
                invokeRestart("muffleWarning")
            }
        }
    )
    if (is.null(af_profiles$signatures) || nrow(af_profiles$signatures) == 0) {
        stop("AF profile extraction did not produce any AF bands.", call. = FALSE)
    }

    if (isTRUE(verbose)) {
        message("  Building AF spectra plot...")
    }
    p <- .build_af_profile_plot(af_profiles$signatures, pd = metadata$pd_meta)
    out <- .new_af_profile_object(
        profile = af_profiles$signatures,
        plot = p,
        source = gated$source,
        extraction = list(
            fcs_file = normalizePath(fcs_file, mustWork = FALSE),
            af_n_bands = af_args$af_n_bands,
            af_max_cells = af_args$af_max_cells,
            af_auto_max_bands = af_args$af_auto_max_bands,
            auto_selection = af_profiles$selection
        )
    )
    if (isTRUE(show_plot)) {
        print(out$plot)
    }
    if (isTRUE(verbose)) {
        message(
            "Done. Extracted ",
            nrow(out$profile),
            " AF band(s) from ",
            nrow(gated$events),
            " scatter-gated event(s)."
        )
    }
    out
}

#' Save an AF profile globally
#'
#' @param name Profile name. Use letters, numbers, `.`, `_`, and `-`.
#' @param x A full reference matrix, AF profile matrix, or
#'   `spectreasy_af_profile` object. Only AF rows are saved.
#' @param overwrite Logical; replace an existing profile with the same name.
#'
#' @return Invisibly returns the saved file path.
#' @export
save_af_profile <- function(name, x, overwrite = FALSE) {
    name <- .validate_af_profile_name(name)
    profile <- .coerce_af_profile_matrix(x, arg_name = "x")
    file_path <- .af_profile_file(name, create_dir = TRUE)
    if (file.exists(file_path) && !isTRUE(overwrite)) {
        stop(
            "AF profile already exists: ", name,
            ". Use overwrite = TRUE to replace it.",
            call. = FALSE
        )
    }

    plot_obj <- if (.is_af_profile_object(x) && inherits(x$plot, "ggplot")) {
        x$plot
    } else {
        .build_af_profile_plot(profile, pd = attr(x, "detector_pd"))
    }
    source <- if (.is_af_profile_object(x) && "source" %in% names(x)) x$source else NULL
    extraction <- if (.is_af_profile_object(x) && "extraction" %in% names(x)) x$extraction else NULL

    out <- .new_af_profile_object(
        name = name,
        profile = profile,
        plot = plot_obj,
        source = source,
        extraction = extraction
    )
    saveRDS(out, file = file_path, version = 3)
    message(
        "Saved AF profile '", name, "' with ",
        nrow(profile), " band(s) and ", ncol(profile), " detector(s) to:\n",
        file_path
    )
    invisible(file_path)
}

#' Load a saved AF profile
#'
#' @param name Profile name.
#' @param show_plot Logical; print the saved AF spectra plot after loading.
#'
#' @return A `spectreasy_af_profile` object.
#' @export
load_af_profile <- function(name, show_plot = FALSE) {
    file_path <- .af_profile_file(name, create_dir = FALSE)
    if (!file.exists(file_path)) {
        stop("AF profile not found: ", name, call. = FALSE)
    }
    out <- readRDS(file_path)
    out <- .new_af_profile_object(
        name = if (!is.null(out$name)) out$name else .validate_af_profile_name(name),
        profile = out$profile,
        plot = if (inherits(out$plot, "ggplot")) out$plot else NULL,
        source = out$source,
        extraction = out$extraction,
        created = if (!is.null(out$created)) out$created else Sys.time()
    )
    if (isTRUE(show_plot)) {
        print(out$plot)
    }
    out
}

#' List saved AF profiles
#'
#' @return A data frame with saved profile metadata and file paths.
#' @export
list_af_profiles <- function() {
    dir_path <- .af_profile_dir_path(create = FALSE)
    if (!dir.exists(dir_path)) {
        return(data.frame(
            name = character(),
            bands = integer(),
            detectors = integer(),
            created = as.POSIXct(character()),
            path = character(),
            stringsAsFactors = FALSE
        ))
    }
    files <- list.files(dir_path, pattern = "\\.rds$", full.names = TRUE)
    rows <- lapply(files, function(path) {
        obj <- tryCatch(readRDS(path), error = function(e) NULL)
        if (is.null(obj) || is.null(obj$profile)) {
            return(NULL)
        }
        profile <- .coerce_af_profile_matrix(obj$profile, arg_name = basename(path))
        data.frame(
            name = if (!is.null(obj$name) && nzchar(obj$name)) obj$name else tools::file_path_sans_ext(basename(path)),
            bands = nrow(profile),
            detectors = ncol(profile),
            created = if (!is.null(obj$created)) as.POSIXct(obj$created) else as.POSIXct(NA),
            path = path,
            stringsAsFactors = FALSE
        )
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows) == 0) {
        return(data.frame(
            name = character(),
            bands = integer(),
            detectors = integer(),
            created = as.POSIXct(character()),
            path = character(),
            stringsAsFactors = FALSE
        ))
    }
    out <- do.call(rbind, rows)
    out[order(out$name), , drop = FALSE]
}

#' Delete a saved AF profile
#'
#' @param name Profile name.
#'
#' @return Invisibly returns `TRUE` if the profile was deleted.
#' @export
delete_af_profile <- function(name) {
    file_path <- .af_profile_file(name, create_dir = FALSE)
    if (!file.exists(file_path)) {
        stop("AF profile not found: ", name, call. = FALSE)
    }
    ok <- unlink(file_path) == 0
    if (!isTRUE(ok)) {
        stop("Could not delete AF profile: ", name, call. = FALSE)
    }
    message("Deleted AF profile '", name, "'.")
    invisible(TRUE)
}

#' Plot an AF profile
#'
#' @param x Profile name, AF profile object, AF profile matrix, or full
#'   reference matrix. Full matrices are filtered to AF rows before plotting.
#' @param show Logical; print the plot.
#'
#' @return A `ggplot` object.
#' @export
plot_af_profile <- function(x, show = TRUE) {
    if (is.character(x) && length(x) == 1L && grepl(.af_profile_name_regex, x)) {
        obj <- load_af_profile(x, show_plot = FALSE)
        p <- obj$plot
    } else if (.is_af_profile_object(x) && inherits(x$plot, "ggplot")) {
        p <- x$plot
    } else {
        p <- .build_af_profile_plot(.coerce_af_profile_matrix(x, arg_name = "x"))
    }
    if (isTRUE(show)) {
        print(p)
    }
    p
}

#' Add a saved AF profile to a reference matrix
#'
#' @param M Reference matrix.
#' @param profile Profile name, AF profile object, or AF profile matrix.
#' @param replace_existing Logical; remove existing AF rows from `M` first.
#'
#' @return A reference matrix with AF profile rows appended.
#' @export
add_af_profile <- function(M, profile, replace_existing = TRUE) {
    M <- .as_reference_matrix(M, "M")
    profile_mat <- if (is.character(profile) && length(profile) == 1L && grepl(.af_profile_name_regex, profile)) {
        load_af_profile(profile, show_plot = FALSE)$profile
    } else {
        .coerce_af_profile_matrix(profile, arg_name = "profile")
    }

    missing_from_M <- setdiff(colnames(profile_mat), colnames(M))
    missing_from_profile <- setdiff(colnames(M), colnames(profile_mat))
    if (length(missing_from_M) > 0 || length(missing_from_profile) > 0) {
        stop(
            "AF profile detectors do not match M.\n",
            if (length(missing_from_M) > 0) paste0("Only in profile: ", paste(missing_from_M, collapse = ", "), "\n") else "",
            if (length(missing_from_profile) > 0) paste0("Only in M: ", paste(missing_from_profile, collapse = ", ")) else "",
            call. = FALSE
        )
    }
    profile_mat <- profile_mat[, colnames(M), drop = FALSE]

    af_rows <- grepl("^AF($|_)", rownames(M), ignore.case = TRUE)
    if (any(af_rows) && isTRUE(replace_existing)) {
        M <- M[!af_rows, , drop = FALSE]
    }
    if (any(rownames(profile_mat) %in% rownames(M))) {
        stop(
            "AF profile row names already exist in M. Use replace_existing = TRUE or remove existing AF rows first.",
            call. = FALSE
        )
    }

    out <- rbind(M, profile_mat)
    attrs <- attributes(M)
    keep_attrs <- setdiff(names(attrs), c("dim", "dimnames", "names"))
    for (nm in keep_attrs) {
        attr(out, nm) <- attrs[[nm]]
    }
    out
}

#' @export
print.spectreasy_af_profile <- function(x, ...) {
    profile <- .coerce_af_profile_matrix(x, arg_name = "x")
    cat("spectreasy AF profile\n")
    if (!is.null(x$name) && nzchar(x$name)) {
        cat("  name: ", x$name, "\n", sep = "")
    }
    cat("  bands: ", nrow(profile), "\n", sep = "")
    cat("  detectors: ", ncol(profile), "\n", sep = "")
    if (!is.null(x$created)) {
        cat("  created: ", format(x$created), "\n", sep = "")
    }
    cat("  matrix: $profile\n")
    cat("  plot: $plot\n")
    invisible(x)
}
