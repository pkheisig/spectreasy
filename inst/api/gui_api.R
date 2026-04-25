# gui_api_adjust.R
# API for the spectreasy Interactive Tuner (Adjustment/Crosstalk Correction)

get_matrix_dir <- function() {
    normalizePath(getOption("spectreasy.matrix_dir", getwd()), mustWork = FALSE)
}

get_samples_dir <- function() {
    default_samples <- file.path(get_matrix_dir(), "samples")
    normalizePath(getOption("spectreasy.samples_dir", default_samples), mustWork = FALSE)
}

get_config_dir <- function() {
    cfg_dir <- file.path(get_matrix_dir(), "gui_configs")
    if (!dir.exists(cfg_dir)) {
        dir.create(cfg_dir, recursive = TRUE, showWarnings = FALSE)
    }
    cfg_dir
}

normalize_config_filename <- function(filename) {
    if (is.null(filename) || !nzchar(trimws(filename))) {
        return("gui_config.json")
    }
    out <- basename(trimws(filename))
    if (!grepl("\\.json$", out, ignore.case = TRUE)) {
        out <- paste0(out, ".json")
    }
    out
}

is_probably_matrix_csv <- function(path) {
    df <- tryCatch(
        utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, nrows = 200),
        error = function(e) NULL
    )
    if (is.null(df) || !is.data.frame(df) || ncol(df) < 2) {
        return(FALSE)
    }

    first_name <- tolower(trimws(colnames(df)[1]))
    first_col <- df[[1]]
    has_marker_col <- first_name %in% c("marker", "fluorophore", "file") || !is.numeric(first_col)
    mat_df <- if (has_marker_col) df[, -1, drop = FALSE] else df
    if (ncol(mat_df) == 0) {
        return(FALSE)
    }

    numeric_ok <- vapply(mat_df, function(col) {
        converted <- suppressWarnings(as.numeric(col))
        all(is.na(col) | !is.na(converted))
    }, logical(1))

    mean(numeric_ok) >= 0.8
}

#* @filter logger
function(req) {
    cat(as.character(Sys.time()), "-", req$REQUEST_METHOD, req$PATH_INFO, "\n")
    plumber::forward()
}

#* @filter cors
function(req, res) {
    res$setHeader("Access-Control-Allow-Origin", "*")
    res$setHeader("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
    res$setHeader("Access-Control-Allow-Headers", "Content-Type")
    if (req$REQUEST_METHOD == "OPTIONS") {
        res$status <- 200
        return(list())
    }
    plumber::forward()
}

#* Health Check
#* @get /status
function() {
    return(list(status = "ok", time = Sys.time(), wd = getwd(), matrix_dir = get_matrix_dir()))
}

#* List available matrices
#* @get /matrices
function() {
    files <- sort(list.files(get_matrix_dir(), pattern = ".*\\.csv$", ignore.case = TRUE))
    if (length(files) == 0) {
        return(character(0))
    }
    paths <- file.path(get_matrix_dir(), files)
    keep <- vapply(paths, is_probably_matrix_csv, logical(1))
    if (!any(keep)) {
        return(as.character(files))
    }
    return(as.character(files[keep]))
}

#* List available sample files
#* @get /samples
function() {
    samples_dir <- get_samples_dir()
    if (!dir.exists(samples_dir)) return(character(0))
    files <- list.files(samples_dir, pattern = "\\.fcs$", ignore.case = TRUE)
    return(as.character(sort(files)))
}

#* CORS preflight for import_sample_content
#* @options /import_sample_content
function(res) {
    return("")
}

#* Import a sample from uploaded binary content
#* @post /import_sample_content
function(req) {
    body <- jsonlite::fromJSON(req$postBody, simplifyVector = FALSE)
    if (is.null(body)) {
        return(list(error = "Missing request body"))
    }

    filename <- if (is.null(body$filename)) "" else as.character(body$filename)[1]
    content_b64 <- if (is.null(body$content_base64)) "" else as.character(body$content_base64)[1]
    filename <- trimws(filename)

    if (!nzchar(filename)) {
        return(list(error = "Missing filename"))
    }
    if (!nzchar(content_b64)) {
        return(list(error = "Missing content_base64"))
    }

    safe_name <- basename(filename)
    if (!grepl("\\.fcs$", safe_name, ignore.case = TRUE)) {
        safe_name <- paste0(safe_name, ".fcs")
    }

    payload <- tryCatch(
        jsonlite::base64_dec(content_b64),
        error = function(e) NULL
    )
    if (is.null(payload) || length(payload) == 0) {
        return(list(error = "Invalid or empty sample content"))
    }

    samples_dir <- get_samples_dir()
    if (!dir.exists(samples_dir)) {
        dir.create(samples_dir, recursive = TRUE, showWarnings = FALSE)
    }

    dest <- file.path(samples_dir, safe_name)
    con <- file(dest, open = "wb")
    on.exit(close(con), add = TRUE)
    writeBin(payload, con)

    return(list(success = TRUE, filename = safe_name, path = dest))
}

#* List GUI config presets
#* @get /configs
function() {
    cfg_dir <- get_config_dir()
    files <- list.files(cfg_dir, pattern = "\\.json$", ignore.case = TRUE)
    return(as.character(sort(files)))
}

#* Load GUI config preset
#* @get /load_config
#* @param filename
function(filename) {
    cfg_name <- normalize_config_filename(filename)
    path <- file.path(get_config_dir(), cfg_name)
    if (!file.exists(path)) {
        return(list(error = paste("Config file not found:", path)))
    }
    cfg <- jsonlite::fromJSON(path, simplifyVector = TRUE)
    return(cfg)
}

#* Save GUI config preset
#* @post /save_config
function(req) {
    body <- jsonlite::fromJSON(req$postBody)
    cfg_name <- normalize_config_filename(body$filename)
    cfg <- body$config_json
    if (is.null(cfg)) {
        return(list(error = "Missing config_json in request body"))
    }
    path <- file.path(get_config_dir(), cfg_name)
    jsonlite::write_json(cfg, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
    return(list(success = TRUE, filename = cfg_name, path = path))
}

#* Load a specific matrix
#* @get /load_matrix
#* @param filename
function(filename) {
    path <- file.path(get_matrix_dir(), filename)
    if (!file.exists(path)) {
        return(list(error = paste("File not found:", path)))
    }

    df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

    # Check if it's a spillover/reference matrix (M) or unmixing matrix (W)
    # Reference/Spillover usually has Markers as rows, Detectors as cols
    # Unmixing usually has Detectors as rows, Markers as cols (or vice versa depending on convention)
    # spectreasy convention:
    # M (Reference): Rows=Markers, Cols=Detectors
    # W (Unmixing): Rows=Markers, Cols=Detectors (so Unmixed = Raw %*% t(W))

    # Ensure the first column is named 'Marker' for the frontend
    # If the CSV has a header but first col is unnamed or "V1", fix it
    if (colnames(df)[1] %in% c("V1", "")) {
        colnames(df)[1] <- "Marker"
    } else if (colnames(df)[1] != "Marker") {
        # Assume first column is Marker if it contains strings
        if (is.character(df[[1]])) {
            colnames(df)[1] <- "Marker"
        }
    }
    return(df)
}

#* Save the adjusted matrix
#* @post /save_matrix
function(req) {
    body <- jsonlite::fromJSON(req$postBody)
    filename <- body$filename
    matrix_data <- body$matrix_json
    df <- as.data.frame(matrix_data, check.names = FALSE)
    path <- file.path(get_matrix_dir(), filename)
    utils::write.csv(df, path, row.names = FALSE, quote = TRUE)
    return(list(success = TRUE, path = path))
}

#* CORS preflight for import_matrix
#* @options /import_matrix
function(res) {
    return("")
}

#* Import a matrix from an external path
#* @post /import_matrix
#* @param path Absolute path to the CSV file
function(path) {
    if (!file.exists(path)) {
        return(list(error = paste("File not found:", path)))
    }
    filename <- basename(path)
    dest <- file.path(get_matrix_dir(), filename)
    file.copy(path, dest, overwrite = TRUE)
    return(list(success = TRUE, filename = filename))
}

#* CORS preflight for import_matrix_content
#* @options /import_matrix_content
function(res) {
    return("")
}

#* Import a matrix from uploaded content
#* @post /import_matrix_content
#* @param filename The filename
#* @param content The CSV content as text
function(filename, content) {
    dest <- file.path(get_matrix_dir(), filename)
    writeLines(content, dest)
    return(list(success = TRUE, filename = filename))
}

#* Get unmixed data for a subset of cells from a sample
#* @get /data
#* @param sample_name The name of the sample to load
function(sample_name = "") {
    # If no sample provided, take the first one in samples/
    samples_dir <- get_samples_dir()
    files <- sort(list.files(samples_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE))
    if (length(files) == 0) {
        return(list(error = paste0("No FCS files found in samples directory: ", samples_dir)))
    }

    if (is.null(sample_name) || !nzchar(trimws(as.character(sample_name)))) {
        sample_path <- files[1]
    } else {
        sn <- trimws(as.character(sample_name))
        candidates <- c(
            file.path(samples_dir, sn),
            file.path(samples_dir, paste0(sn, ".fcs"))
        )
        sample_path <- candidates[file.exists(candidates)][1]
        if (is.na(sample_path) || !nzchar(sample_path)) {
            target <- tolower(tools::file_path_sans_ext(basename(sn)))
            file_ids <- tolower(tools::file_path_sans_ext(basename(files)))
            idx <- which(file_ids == target)
            if (length(idx) > 0) {
                sample_path <- files[idx[1]]
            } else {
                return(list(error = paste("Sample not found:", sn)))
            }
        }
    }

    ff <- flowCore::read.FCS(sample_path, transformation = FALSE, truncate_max_range = FALSE)
    raw_data <- flowCore::exprs(ff)

    # Subsample for speed - smaller for fast interactive updates
    n_sub <- 2000
    if (nrow(raw_data) > n_sub) {
        set.seed(123)
        raw_data <- raw_data[sample(nrow(raw_data), n_sub), ]
    }

    pd <- flowCore::pData(flowCore::parameters(ff))
    # Helper to get sorted detectors (copying logic from spectreasy if not exported)
    # Assuming spectreasy is loaded or we implement basic logic
    det_info <- tryCatch(
        {
            spectreasy::get_sorted_detectors(pd)
        },
        error = function(e) {
            # Fallback if function not accessible
            fl_cols <- grep("FL", pd$name, value = TRUE)
            list(names = fl_cols, labels = pd$desc[match(fl_cols, pd$name)])
        }
    )

    return(list(
        raw_data = as.data.frame(raw_data),
        sample_name = basename(sample_path),
        detector_names = det_info$names,
        detector_labels = det_info$labels
    ))
}

#* CORS preflight for unmix
#* @options /unmix
function(res) {
    return("")
}

#* Run unmixing (On-demand unmixing endpoint)
#* @post /unmix
#* @param matrix_json The matrix (M or W)
#* @param raw_data_json The raw data
#* @param type "reference" (M) or "unmixing" (W)
function(matrix_json, raw_data_json, type = "reference") {
    # matrix_json format: {MarkerName: {det1: val, det2: val, ...}, ...}
    markers <- names(matrix_json)
    detectors <- names(matrix_json[[1]])

    mat <- matrix(0, nrow = length(markers), ncol = length(detectors))
    for (i in seq_along(markers)) {
        mat[i, ] <- as.numeric(unlist(matrix_json[[markers[i]]][detectors]))
    }
    rownames(mat) <- markers
    colnames(mat) <- detectors

    Y <- as.matrix(as.data.frame(raw_data_json))

    # Matching columns
    common_dets <- intersect(colnames(Y), colnames(mat))
    if (length(common_dets) == 0) {
        return(list(error = "No matching detectors found between data and matrix"))
    }

    Y_sub <- Y[, common_dets, drop = FALSE]
    mat_sub <- mat[, common_dets, drop = FALSE]

    # Perform Unmixing
    if (grepl("unmixing", tolower(type))) {
        # Provided matrix IS the unmixing matrix (W)
        # Unmixed = Raw * t(W)
        unmixed <- Y_sub %*% t(mat_sub)
    } else {
        # Provided matrix IS the reference matrix (M), derive W via OLS
        # W = (M M^t)^-1 M  ? No.
        # Simple OLS: Y ~ U * M  => U = Y * M^T * (M * M^T)^-1 ??
        # Standard OLS unmixing: U = Y * pseudoinverse(M)
        # U = Y * pinv(M) = Y * t(M) * (M * t(M))^-1  (if M is tall? No M is Markers x Detectors, usually Fat)
        # Actually in spectral: Y (Cells x Dets) = U (Cells x Markers) * M (Markers x Dets)
        # U = Y * M_pinv
        # M_pinv = M^T (M M^T)^-1
        # So W^T = M^T (M M^T)^-1
        # W = (M M^T)^-1 M

        W <- tryCatch(
            {
                spectreasy::derive_unmixing_matrix(mat_sub, method = "OLS")
            },
            error = function(e) {
                # Fallback OLS
                t(MASS::ginv(t(mat_sub)))
            }
        )

        unmixed <- Y_sub %*% t(W)
    }

    return(as.data.frame(unmixed))
}
