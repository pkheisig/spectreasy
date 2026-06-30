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

matrix_path <- function(filename) {
    file.path(get_matrix_dir(), basename(trimws(as.character(filename)[1])))
}

read_matrix_csv <- function(path) {
    df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    if (ncol(df) > 0 && colnames(df)[1] %in% c("V1", "")) {
        colnames(df)[1] <- "Marker"
    } else if (ncol(df) > 0 && colnames(df)[1] != "Marker" && is.character(df[[1]])) {
        colnames(df)[1] <- "Marker"
    }
    df
}

is_af_matrix_row <- function(df) {
    if (!is.data.frame(df) || ncol(df) == 0) return(logical(0))
    grepl("^AF($|_)", as.character(df[[1]]), ignore.case = TRUE)
}

merge_hidden_af_rows <- function(df, source_paths = character()) {
    if (!is.data.frame(df) || ncol(df) == 0) {
        return(df)
    }
    if (any(is_af_matrix_row(df))) {
        return(df)
    }

    source_paths <- unique(source_paths[file.exists(source_paths)])
    for (source_path in source_paths) {
        existing_df <- read_matrix_csv(source_path)
        af_rows <- existing_df[is_af_matrix_row(existing_df), , drop = FALSE]
        if (nrow(af_rows) == 0) next

        colnames(af_rows)[1] <- colnames(df)[1]
        af_rows_aligned <- af_rows[, intersect(colnames(df), colnames(af_rows)), drop = FALSE]
        missing_cols <- setdiff(colnames(df), colnames(af_rows_aligned))
        for (m_col in missing_cols) {
            af_rows_aligned[[m_col]] <- 0
        }
        af_rows_aligned <- af_rows_aligned[, colnames(df), drop = FALSE]
        return(rbind(df, af_rows_aligned))
    }

    df
}

raw_data_to_df <- function(raw_data_json) {
    if (is.data.frame(raw_data_json)) {
        return(as.data.frame(raw_data_json, check.names = FALSE))
    }
    if (is.list(raw_data_json) && length(raw_data_json) > 0 && all(vapply(raw_data_json, is.list, logical(1)))) {
        rows <- lapply(raw_data_json, function(row) {
            as.data.frame(row, check.names = FALSE, stringsAsFactors = FALSE)
        })
        return(do.call(rbind, rows))
    }
    as.data.frame(raw_data_json, check.names = FALSE)
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
    res$setHeader("Cache-Control", "no-store, no-cache, must-revalidate, max-age=0")
    res$setHeader("Pragma", "no-cache")
    res$setHeader("Expires", "0")
    if (req$REQUEST_METHOD == "OPTIONS") {
        res$status <- 200
        return(list())
    }
    plumber::forward()
}

#* Health Check
#* @get /status
function() {
    return(list(
        status = "ok",
        time = Sys.time(),
        wd = getwd(),
        matrix_dir = get_matrix_dir(),
        samples_dir = get_samples_dir(),
        gui_mode = getOption("spectreasy.gui_mode", "tuner"),
        panel_cytometer = getOption("spectreasy.panel_cytometer", "aurora")
    ))
}

#* Spectral panel builder metadata and current selection
#* @get /spectral_panel
#* @param cytometer
#* @param configuration
function(cytometer = "", configuration = "") {
    selected_cytometer <- if (is.null(cytometer) || !nzchar(trimws(as.character(cytometer)[1]))) {
        getOption("spectreasy.panel_cytometer", "aurora")
    } else {
        cytometer
    }
    selected_configuration <- if (is.null(configuration) || !nzchar(trimws(as.character(configuration)[1]))) NULL else configuration
    tryCatch(
        spectreasy:::.spectral_panel_payload(
            cytometer = selected_cytometer,
            fluorophores = character(),
            configuration = selected_configuration
        ),
        error = function(e) list(error = conditionMessage(e))
    )
}

#* CORS preflight for spectral_panel_metrics
#* @options /spectral_panel_metrics
function(res) {
    return("")
}

#* Recalculate spectral panel metrics for selected fluorophores
#* @post /spectral_panel_metrics
function(req) {
    body <- jsonlite::fromJSON(req$postBody, simplifyVector = TRUE)
    cytometer <- if (!is.null(body$cytometer)) body$cytometer else getOption("spectreasy.panel_cytometer", "aurora")
    configuration <- if (!is.null(body$configuration)) body$configuration else NULL
    fluorophores <- if (!is.null(body$fluorophores)) body$fluorophores else character()
    tryCatch(
        spectreasy:::.spectral_panel_payload(
            cytometer = cytometer,
            fluorophores = fluorophores,
            configuration = configuration
        ),
        error = function(e) list(error = conditionMessage(e))
    )
}

#* CORS preflight for export_spectral_panel_overview
#* @options /export_spectral_panel_overview
function(res) {
    return("")
}

#* Export spectral panel overview PDF
#* @post /export_spectral_panel_overview
function(req) {
    body <- jsonlite::fromJSON(req$postBody, simplifyVector = TRUE)
    cytometer <- if (!is.null(body$cytometer)) body$cytometer else getOption("spectreasy.panel_cytometer", "aurora")
    configuration <- if (!is.null(body$configuration)) body$configuration else NULL
    fluorophores <- if (!is.null(body$fluorophores)) body$fluorophores else character()
    markers <- if (!is.null(body$markers)) body$markers else character()

    tryCatch({
        output_file <- tempfile("spectreasy_panel_overview_", fileext = ".pdf")
        on.exit(unlink(output_file), add = TRUE)
        spectreasy:::.write_spectral_panel_overview_pdf(
            cytometer = cytometer,
            configuration = configuration,
            fluorophores = fluorophores,
            markers = markers,
            output_file = output_file
        )
        payload <- readBin(output_file, what = "raw", n = file.info(output_file)$size)
        list(
            filename = paste0("spectreasy_", cytometer, "_", ifelse(is.null(configuration), "panel", configuration), "_overview.pdf"),
            content_type = "application/pdf",
            content_base64 = jsonlite::base64_enc(payload)
        )
    }, error = function(e) list(error = conditionMessage(e)))
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
    path <- matrix_path(filename)
    if (!file.exists(path)) {
        return(list(error = paste("File not found:", path)))
    }

    df <- read_matrix_csv(path)

    # Filter out AF rows (matching ^AF($|_) case-insensitively)
    df <- df[!is_af_matrix_row(df), , drop = FALSE]

    return(df)
}

#* Save the adjusted matrix
#* @post /save_matrix
function(req) {
    body <- jsonlite::fromJSON(req$postBody)
    filename <- body$filename
    source_filename <- body$source_filename
    matrix_data <- body$matrix_json
    df <- as.data.frame(matrix_data, check.names = FALSE)
    path <- matrix_path(filename)
    source_paths <- c(
        path,
        if (!is.null(source_filename) && nzchar(trimws(as.character(source_filename)[1]))) {
            matrix_path(source_filename)
        } else {
            character()
        }
    )
    df <- merge_hidden_af_rows(df, source_paths = source_paths)

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
    dest <- matrix_path(filename)
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
function(matrix_json, raw_data_json, type = "reference", matrix_filename = "", method = "WLS") {
    # matrix_json format: {MarkerName: {det1: val, det2: val, ...}, ...}
    markers <- names(matrix_json)
    detectors <- names(matrix_json[[1]])

    mat <- matrix(0, nrow = length(markers), ncol = length(detectors))
    for (i in seq_along(markers)) {
        mat[i, ] <- as.numeric(unlist(matrix_json[[markers[i]]][detectors]))
    }
    rownames(mat) <- markers
    colnames(mat) <- detectors

    Y <- as.matrix(raw_data_to_df(raw_data_json))
    suppressWarnings(storage.mode(Y) <- "numeric")

    if (!grepl("unmixing", tolower(type)) &&
        !is.null(matrix_filename) &&
        nzchar(trimws(as.character(matrix_filename)[1]))) {
        mat_df <- as.data.frame(mat, check.names = FALSE)
        mat_df$Marker <- rownames(mat)
        mat_df <- mat_df[, c("Marker", colnames(mat)), drop = FALSE]
        mat_df <- merge_hidden_af_rows(mat_df, source_paths = matrix_path(matrix_filename))
        markers <- as.character(mat_df[[1]])
        mat <- as.matrix(mat_df[, -1, drop = FALSE])
        storage.mode(mat) <- "numeric"
        rownames(mat) <- markers
    }

    # Matching columns
    common_dets <- intersect(colnames(Y), colnames(mat))
    if (length(common_dets) == 0) {
        marker_cols <- intersect(colnames(Y), markers)
        if (length(marker_cols) > 0) {
            return(as.data.frame(Y[, marker_cols, drop = FALSE]))
        }
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
        method_upper <- toupper(trimws(as.character(method)[1]))
        if (!method_upper %in% c("WLS", "OLS", "NNLS")) {
            method_upper <- "WLS"
        }
        ff <- flowCore::flowFrame(Y_sub)
        res <- tryCatch(
            spectreasy::calc_residuals(ff, mat_sub, method = method_upper),
            error = function(e) {
                return(list(error = conditionMessage(e)))
            }
        )
        if (is.list(res) && !is.data.frame(res) && !is.null(res$error)) {
            return(res)
        }
        unmixed <- as.matrix(res[, rownames(mat_sub), drop = FALSE])
    }

    return(as.data.frame(unmixed))
}
