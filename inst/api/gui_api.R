# gui_api_adjust.R
# API for the spectreasy Interactive Tuner (Adjustment/Crosstalk Correction)

get_gui_default_project_dir <- function() {
    configured <- getOption("spectreasy.project_dir", "")
    if (!is.null(configured) && length(configured) > 0 && nzchar(trimws(as.character(configured[1])))) {
        return(normalizePath(as.character(configured[1]), mustWork = FALSE))
    }

    cwd <- normalizePath(getwd(), mustWork = FALSE)
    package_root <- normalizePath(file.path(cwd, "..", ".."), mustWork = FALSE)
    is_source_package <- file.exists(file.path(package_root, "DESCRIPTION")) &&
        file.exists(file.path(package_root, "inst", "api", "gui_api.R"))
    if (is_source_package) package_root else cwd
}

get_matrix_dir <- function() {
    configured <- getOption("spectreasy.matrix_dir", getOption("spectreasy.project_dir", ""))
    if (is.null(configured) || length(configured) == 0 || !nzchar(trimws(as.character(configured[1])))) {
        configured <- get_gui_default_project_dir()
    }
    normalizePath(as.character(configured[1]), mustWork = FALSE)
}

get_samples_dir <- function() {
    default_samples <- file.path(get_matrix_dir(), "samples")
    normalizePath(getOption("spectreasy.samples_dir", default_samples), mustWork = FALSE)
}

get_unmixing_method <- function() {
    method <- getOption("spectreasy.unmixing_method", "Spectreasy")
    spectreasy:::.normalize_unmix_method(method)
}

get_config_dir <- function() {
    get_user_gui_config_dir()
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

matrix_filename_contains_matrix <- function(filename) {
    if (is.null(filename) || length(filename) == 0 || is.na(filename[1])) {
        return(FALSE)
    }
    normalized <- tolower(gsub("[^[:alnum:]]+", "", basename(as.character(filename)[1])))
    isTRUE(grepl("matrix", normalized, fixed = TRUE))
}

list_matrix_csv_files <- function() {
    matrix_dir <- get_matrix_dir()
    if (!dir.exists(matrix_dir)) {
        return(character(0))
    }

    files <- list.files(
        matrix_dir,
        pattern = "\\.csv$",
        recursive = TRUE,
        full.names = FALSE,
        ignore.case = TRUE
    )
    files <- gsub("\\\\", "/", files)
    files <- files[vapply(files, matrix_filename_contains_matrix, logical(1))]
    sort(files)
}

matrix_path <- function(filename) {
    if (is.null(filename) || length(filename) == 0 || is.na(filename[1])) {
        stop("Invalid matrix filename")
    }
    rel <- trimws(as.character(filename)[1])
    rel <- gsub("\\\\", "/", rel)
    rel <- sub("^/+", "", rel)
    parts <- strsplit(rel, "/", fixed = TRUE)[[1]]
    if (!isTRUE(nzchar(rel)) || any(parts %in% c("", ".", ".."))) {
        stop("Invalid matrix filename")
    }
    file.path(get_matrix_dir(), do.call(file.path, as.list(parts)))
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

get_user_gui_config_dir <- function() {
    cfg_dir <- file.path(tools::R_user_dir("spectreasy", which = "config"), "gui_configs")
    if (!dir.exists(cfg_dir)) {
        dir.create(cfg_dir, recursive = TRUE, showWarnings = FALSE)
    }
    cfg_dir
}

normalize_gui_module <- function(module) {
    module <- trimws(as.character(module)[1])
    if (is.na(module) || !nzchar(module)) module <- "matrix_tuner"
    module <- gsub("[^A-Za-z0-9_-]+", "_", module)
    module
}

user_gui_config_path <- function(module) {
    file.path(get_user_gui_config_dir(), paste0(normalize_gui_module(module), ".json"))
}

get_gate_scc_dir <- function() {
    normalizePath(getOption("spectreasy.gating_scc_dir", file.path(get_gui_default_project_dir(), "scc")), mustWork = FALSE)
}

get_gate_control_file <- function() {
    normalizePath(getOption("spectreasy.gating_control_file", file.path(get_gui_default_project_dir(), "fcs_mapping.csv")), mustWork = FALSE)
}

get_gate_file <- function() {
    normalizePath(getOption("spectreasy.gating_gate_file", file.path(get_gui_default_project_dir(), "ssc_gate_config.csv")), mustWork = FALSE)
}

set_gate_file <- function(path) {
    path <- normalizePath(path, mustWork = FALSE)
    options(spectreasy.gating_gate_file = path)
    path
}

gate_working_dir <- function() {
    get_gui_default_project_dir()
}

gate_applescript_quote <- function(x) {
    x <- gsub("\\\\", "\\\\\\\\", as.character(x)[1])
    x <- gsub('"', '\\\\"', x, fixed = TRUE)
    paste0('"', x, '"')
}

gate_has_system_file_picker <- function() {
    sysname <- Sys.info()[["sysname"]]
    if (identical(sysname, "Darwin")) return(nzchar(Sys.which("osascript")))
    if (identical(sysname, "Windows")) return(nzchar(Sys.which("powershell")))
    if (identical(sysname, "Linux")) return(nzchar(Sys.which("zenity")) || nzchar(Sys.which("kdialog")))
    FALSE
}

gate_powershell_quote <- function(x) {
    paste0("'", gsub("'", "''", as.character(x)[1], fixed = TRUE), "'")
}

gate_pick_csv_file_windows <- function(mode, initial_dir, default_name) {
    title <- if (mode == "save") "Save gate CSV" else "Load gate CSV"
    dialog_class <- if (mode == "save") "SaveFileDialog" else "OpenFileDialog"
    ps_cmd <- paste0(
        "Add-Type -AssemblyName System.Windows.Forms; ",
        "$f = New-Object System.Windows.Forms.", dialog_class, "; ",
        "$f.Title = ", gate_powershell_quote(title), "; ",
        "$f.Filter = 'CSV files (*.csv)|*.csv'; ",
        "$f.InitialDirectory = ", gate_powershell_quote(initial_dir), "; ",
        if (mode == "save") paste0("$f.FileName = ", gate_powershell_quote(default_name), "; "),
        "if ($f.ShowDialog() -eq [System.Windows.Forms.DialogResult]::OK) { $f.FileName } else { write-output 'CANCEL' }"
    )
    res <- tryCatch(
        system2("powershell", c("-NoProfile", "-Command", ps_cmd), stdout = TRUE, stderr = FALSE),
        error = function(e) character()
    )
    if (length(res) > 0) {
        return(trimws(res[length(res)]))
    }
    character()
}

gate_pick_csv_file_linux <- function(mode, initial_dir, default_name) {
    if (nzchar(Sys.which("zenity"))) {
        args <- c("--file-selection", "--title", if (mode == "save") "Save gate CSV" else "Load gate CSV")
        if (mode == "save") {
            args <- c(args, "--save", "--confirm-overwrite", "--filename", file.path(initial_dir, default_name))
        } else {
            args <- c(args, "--filename", initial_dir)
        }
        args <- c(args, "--file-filter=CSV files (*.csv) | *.csv")
        res <- tryCatch(
            system2("zenity", args, stdout = TRUE, stderr = TRUE),
            error = function(e) character()
        )
        status <- attr(res, "status")
        if (identical(status, 1L)) return("CANCEL")
        if (length(res) > 0 && (is.null(status) || identical(status, 0L))) {
            return(trimws(res[1]))
        }
    } else if (nzchar(Sys.which("kdialog"))) {
        args <- if (mode == "save") {
            c("--getsavefilename", file.path(initial_dir, default_name), "*.csv", "--title", "Save gate CSV")
        } else {
            c("--getopenfilename", initial_dir, "*.csv", "--title", "Load gate CSV")
        }
        res <- tryCatch(
            system2("kdialog", args, stdout = TRUE, stderr = TRUE),
            error = function(e) character()
        )
        status <- attr(res, "status")
        if (identical(status, 1L)) return("CANCEL")
        if (length(res) > 0 && (is.null(status) || identical(status, 0L))) {
            return(trimws(res[1]))
        }
    }
    character()
}

gate_pick_csv_file <- function(mode = c("open", "save")) {
    mode <- match.arg(tolower(mode[1]), c("open", "save"))
    initial_dir <- gate_working_dir()
    if (!dir.exists(initial_dir)) {
        initial_dir <- gate_working_dir()
    }
    default_name <- basename(get_gate_file())
    sysname <- Sys.info()[["sysname"]]

    if (identical(sysname, "Darwin")) {
        if (nzchar(Sys.which("osascript"))) {
            prompt <- if (mode == "save") "Save gate CSV" else "Load gate CSV"
            action <- if (mode == "save") {
                paste0(
                    "choose file name with prompt ", gate_applescript_quote(prompt),
                    " default name ", gate_applescript_quote(default_name),
                    " default location defaultFolder"
                )
            } else {
                paste0(
                    "choose file with prompt ", gate_applescript_quote(prompt),
                    " default location defaultFolder"
                )
            }
            fallback_action <- if (mode == "save") {
                paste0(
                    "choose file name with prompt ", gate_applescript_quote(prompt),
                    " default name ", gate_applescript_quote(default_name)
                )
            } else {
                paste0(
                    "choose file with prompt ", gate_applescript_quote(prompt)
                )
            }
            script <- c(
                "try",
                paste0("  set defaultFolder to POSIX file ", gate_applescript_quote(paste0(initial_dir, "/")), " as alias"),
                paste0("  set chosenPath to ", action),
                "on error err number errNum",
                "  if errNum is -128 then",
                "    error number -128",
                "  else",
                paste0("    set chosenPath to ", fallback_action),
                "  end if",
                "end try",
                "POSIX path of chosenPath"
            )
            selected <- tryCatch(
                system2("osascript", as.vector(rbind("-e", script)), stdout = TRUE, stderr = TRUE),
                error = function(e) character()
            )
            status <- attr(selected, "status")
            if (any(grepl("User canceled", selected, ignore.case = TRUE))) {
                return("CANCEL")
            }
            if (length(selected) > 0 && (is.null(status) || identical(status, 0L))) {
                path <- trimws(selected[length(selected)])
                if (nzchar(path)) {
                    return(normalizePath(path, mustWork = mode == "open"))
                }
            }
            return("")
        }
    }

    if (identical(sysname, "Windows")) {
        if (nzchar(Sys.which("powershell"))) {
            path <- gate_pick_csv_file_windows(mode, initial_dir, default_name)
            if (identical(path, "CANCEL")) return("CANCEL")
            if (nzchar(path)) {
                return(normalizePath(path, mustWork = mode == "open"))
            }
            return("")
        }
    }

    if (identical(sysname, "Linux")) {
        path <- gate_pick_csv_file_linux(mode, initial_dir, default_name)
        if (identical(path, "CANCEL")) return("CANCEL")
        if (nzchar(path)) {
            return(normalizePath(path, mustWork = mode == "open"))
        }
        return("")
    }

    ""
}

gate_ensure_csv_extension <- function(path) {
    if (!grepl("\\.csv$", path, ignore.case = TRUE)) {
        path <- paste0(path, ".csv")
    }
    normalizePath(path, mustWork = FALSE)
}

gate_write_config_csv <- function(rows, path) {
    rows <- gate_normalize_config_rows(rows)
    if (nrow(rows) == 0) rows <- gate_empty_config()
    path <- gate_ensure_csv_extension(path)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile("gate_config_", tmpdir = dirname(path), fileext = ".csv")
    utils::write.csv(rows, tmp, row.names = FALSE, quote = TRUE)
    if (!file.rename(tmp, path)) {
        file.copy(tmp, path, overwrite = TRUE)
        unlink(tmp)
    }
    list(path = path, rows = nrow(rows))
}

gate_empty_control_df <- function() {
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

gate_empty_config <- function() {
    data.frame(
        gate_type = character(),
        scope = character(),
        filename = character(),
        x_channel = character(),
        y_channel = character(),
        plot_mode = character(),
        vertex_index = integer(),
        x = numeric(),
        y = numeric(),
        stringsAsFactors = FALSE
    )
}

gate_normalize_config_rows <- function(rows) {
    template <- gate_empty_config()
    if (is.null(rows) || !is.data.frame(rows) || nrow(rows) == 0) return(template)

    rows <- as.data.frame(rows, stringsAsFactors = FALSE, check.names = FALSE)
    for (nm in colnames(rows)) {
        if (is.list(rows[[nm]])) {
            rows[[nm]] <- vapply(rows[[nm]], function(value) {
                value <- unlist(value, recursive = TRUE, use.names = FALSE)
                if (length(value) == 0 || is.null(value[1])) return(NA_character_)
                as.character(value[1])
            }, character(1))
        }
    }
    for (nm in setdiff(colnames(template), colnames(rows))) {
        rows[[nm]] <- NA
    }
    rows <- rows[, colnames(template), drop = FALSE]

    char_cols <- c("gate_type", "scope", "filename", "x_channel", "y_channel", "plot_mode")
    rows[char_cols] <- lapply(rows[char_cols], function(x) {
        x <- as.character(x)
        x[is.na(x)] <- ""
        x
    })
    rows$vertex_index <- suppressWarnings(as.integer(rows$vertex_index))
    rows$x <- suppressWarnings(as.numeric(rows$x))
    rows$y <- suppressWarnings(as.numeric(rows$y))

    keep <- rows$plot_mode %in% c("missing", "blocked") |
        rows$gate_type == "setting" |
        (is.finite(rows$x) & is.finite(rows$y) & !is.na(rows$vertex_index) & rows$vertex_index > 0)
    rows[keep, , drop = FALSE]
}

gate_safe_basename <- function(x) {
    x <- basename(trimws(as.character(x)[1]))
    if (!grepl("\\.fcs$", x, ignore.case = TRUE) && nzchar(x)) x <- paste0(x, ".fcs")
    x
}

gate_read_fcs <- function(path) {
    flowCore::read.FCS(path, transformation = FALSE, truncate_max_range = FALSE)
}

gate_pick_channel <- function(col_names, prefix, suffix) {
    normalized <- toupper(col_names)
    exact <- grep(paste0("^", prefix, "[0-9]*-", suffix, "$"), normalized, value = TRUE)
    if (length(exact) > 0) return(exact[1])
    loose <- grep(paste0("^", prefix, ".*-", suffix, "$"), normalized, value = TRUE)
    if (length(loose) > 0) return(loose[1])
    fallback <- grep(paste0("^", prefix), normalized, value = TRUE)
    if (length(fallback) > 0) return(fallback[1])
    NA_character_
}

gate_paired_height_channel <- function(area_channel, col_names) {
    candidate <- sub("-A$", "-H", area_channel, ignore.case = TRUE)
    if (!is.na(candidate) && candidate %in% col_names) return(candidate)
    gate_pick_channel(col_names, "FSC", "H")
}

gate_fcs_metadata <- function(path) {
    ff <- gate_read_fcs(path)
    expr <- flowCore::exprs(ff)
    pd <- flowCore::pData(flowCore::parameters(ff))
    fsc_a <- gate_pick_channel(colnames(expr), "FSC", "A")
    ssc_a <- gate_pick_channel(colnames(expr), "SSC", "A")
    fsc_h <- gate_paired_height_channel(fsc_a, colnames(expr))
    desc <- setNames(as.character(pd$desc), as.character(pd$name))
    list(
        total_events = nrow(expr),
        columns = colnames(expr),
        labels = as.list(desc),
        fsc_a = fsc_a,
        fsc_h = fsc_h,
        ssc_a = ssc_a
    )
}

gate_guess_fcs_info <- function(filename) {
    path <- file.path(get_gate_scc_dir(), filename)
    ff <- tryCatch(gate_read_fcs(path), error = function(e) NULL)
    if (is.null(ff)) {
        return(list(fluorophore = "Unknown", marker = "Unknown", channel = ""))
    }
    pd <- flowCore::pData(flowCore::parameters(ff))
    expr <- flowCore::exprs(ff)
    det <- tryCatch(spectreasy::get_sorted_detectors(pd)$names, error = function(e) character())
    det <- intersect(det, colnames(expr))
    if (length(det) == 0) {
        det <- grep("-A$", colnames(expr), value = TRUE)
        det <- det[!grepl("^FSC|^SSC|^Time|^Event", det, ignore.case = TRUE)]
    }
    peak_channel <- if (length(det) > 0) {
        vars <- apply(expr[, det, drop = FALSE], 2, stats::var, na.rm = TRUE)
        names(vars)[which.max(vars)]
    } else {
        colnames(expr)[1]
    }
    desc <- pd$desc[match(peak_channel, pd$name)]
    marker <- if (!is.na(desc) && nzchar(desc)) desc else peak_channel
    fluorophore <- marker
    fn_lower <- tolower(filename)
    if (grepl("af_|unstained|autofluor", fn_lower)) {
        fluorophore <- "AF"
        marker <- "Autofluorescence"
    }
    list(fluorophore = fluorophore, marker = marker, channel = peak_channel)
}

gate_normalize_control_type <- function(value) {
    clean <- tolower(trimws(as.character(value)[1]))
    if (startsWith(clean, "bead")) "beads" else "cells"
}

gate_negative_source_key <- function(value) {
    tools::file_path_sans_ext(basename(trimws(as.character(value))))
}

gate_external_negative_available <- function(df) {
    if (is.null(df) || nrow(df) == 0L) return(logical())

    fluor <- trimws(as.character(df$fluorophore))
    marker <- trimws(as.character(df$marker))
    control_type <- as.character(df$control.type)
    viability <- tolower(trimws(as.character(df$is.viability))) %in% c("true", "t", "1", "yes", "y")
    bead_af <- control_type == "beads" & (
        grepl("^AF_beads?$", fluor, ignore.case = TRUE) |
            vapply(df$filename, spectreasy:::.reference_is_bead_negative_file, logical(1)) |
            grepl("bead.*background", marker, ignore.case = TRUE)
    )
    dead_af <- spectreasy:::.is_dead_af_control_row(
        fluorophore = fluor,
        marker = marker,
        filename = df$filename,
        control_type = control_type
    )
    primary_af <- control_type != "beads" & !dead_af &
        spectreasy:::.is_primary_af_control_row(
            fluorophore = fluor,
            marker = marker,
            filename = df$filename
        )

    has_primary_af <- any(primary_af)
    has_dead_af <- any(dead_af)
    has_bead_af <- any(bead_af)
    file_keys <- gate_negative_source_key(df$filename)
    universal_negative <- trimws(as.character(df$universal.negative))
    universal_upper <- toupper(universal_negative)
    universal_keys <- gate_negative_source_key(universal_negative)
    explicitly_mapped <- (universal_upper %in% c("TRUE", "AF") & has_primary_af) |
        (nzchar(universal_keys) &
            !universal_upper %in% c("FALSE", "TRUE", "AF") &
            universal_keys %in% file_keys)

    corresponding_source <- ifelse(
        control_type == "beads",
        has_bead_af,
        ifelse(viability, has_dead_af, has_primary_af)
    )
    !df$is_af & (explicitly_mapped | corresponding_source)
}

gate_read_mapping <- function() {
    path <- get_gate_control_file()
    df <- if (file.exists(path)) {
        utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
        gate_empty_control_df()
    }
    if (!all(c("filename", "fluorophore", "marker", "channel") %in% colnames(df))) {
        df <- gate_empty_control_df()
    }
    df <- df[!startsWith(basename(df$filename), "._"), , drop = FALSE]
    fcs_files <- list.files(get_gate_scc_dir(), pattern = "\\.fcs$", ignore.case = TRUE)
    fcs_files <- fcs_files[!startsWith(fcs_files, "._")]
    missing_files <- setdiff(fcs_files, df$filename)
    if (length(missing_files) > 0) {
        new_rows <- lapply(missing_files, function(f) {
            info <- gate_guess_fcs_info(f)
            data.frame(
                filename = f,
                fluorophore = info$fluorophore,
                marker = info$marker,
                channel = info$channel,
                control.type = "cells",
                universal.negative = "",
                is.viability = "FALSE",
                stringsAsFactors = FALSE
            )
        })
        df <- rbind(df, do.call(rbind, new_rows))
    }
    df <- df[df$filename %in% fcs_files, , drop = FALSE]
    df <- gui_filter_active_af_mapping(df)
    if (!("control.type" %in% colnames(df))) df$control.type <- "cells"
    df$control.type <- vapply(df$control.type, gate_normalize_control_type, character(1))
    if (!("universal.negative" %in% colnames(df))) df$universal.negative <- ""
    if (!("is.viability" %in% colnames(df))) df$is.viability <- "FALSE"
    df$file_exists <- rep(TRUE, nrow(df))
    df$is_af <- grepl("^AF($|_|\\b)", df$fluorophore, ignore.case = TRUE) |
        grepl("autofluorescence|background", df$marker, ignore.case = TRUE)
    df$uses_histogram_gates <- !df$is_af
    external_negative_available <- gate_external_negative_available(df)
    df$uses_negative_histogram_gate <- df$uses_histogram_gates &
        !external_negative_available
    df$is_viability <- tolower(as.character(df$is.viability)) %in% c("true", "t", "1", "yes")
    df$id <- tools::file_path_sans_ext(basename(df$filename))
    df
}

gate_compute_extent <- function(v) {
    finite <- v[is.finite(v)]
    if (length(finite) == 0) return(c(0, 1))
    val_max <- as.numeric(stats::quantile(finite, 0.998, na.rm = TRUE))
    if (!is.finite(val_max) || val_max <= 0) val_max <- max(finite, na.rm = TRUE)
    if (!is.finite(val_max) || val_max <= 0) val_max <- 1
    c(0, val_max + val_max * 0.04)
}

gate_scatter_channels <- function(col_names) {
    hits <- grep("^(FSC|SSC).*-(A|H|W)$", col_names, value = TRUE, ignore.case = TRUE)
    unique(hits)
}

gate_point_in_polygon <- function(x, y, poly) {
    n <- nrow(poly)
    inside <- rep(FALSE, length(x))
    if (n < 3) return(inside)
    j <- n
    for (i in seq_len(n)) {
        xi <- poly[i, 1]
        yi <- poly[i, 2]
        xj <- poly[j, 1]
        yj <- poly[j, 2]
        intersect <- ((yi > y) != (yj > y)) & (x < (xj - xi) * (y - yi) / (yj - yi + 1e-12) + xi)
        inside[intersect] <- !inside[intersect]
        j <- i
    }
    inside
}

gate_vertices_matrix <- function(gate) {
    if (is.null(gate$vertices) || length(gate$vertices) == 0) return(NULL)
    verts <- do.call(rbind, lapply(gate$vertices, function(v) c(as.numeric(v$x), as.numeric(v$y))))
    if (!is.matrix(verts) || nrow(verts) == 0 || ncol(verts) < 2) return(NULL)
    verts[stats::complete.cases(verts), , drop = FALSE]
}

gate_value <- function(gate, primary, fallback = NULL) {
    value <- gate[[primary]]
    if (!is.null(value)) return(value)
    alt <- switch(primary,
        xChannel = "x_channel",
        yChannel = "y_channel",
        NULL
    )
    if (!is.null(alt) && !is.null(gate[[alt]])) return(gate[[alt]])
    fallback
}

gate_cached_gate <- function(cache, gate_type, filename, control_type) {
    if (is.null(cache)) return(NULL)
    file_key <- paste0(gate_type, ":", filename)
    global_key <- paste0(gate_type, ":", control_type)
    if (file_key %in% names(cache)) return(cache[[file_key]])
    if (global_key %in% names(cache)) return(cache[[global_key]])
    NULL
}

gate_apply_polygon_gate <- function(expr, gate) {
    expr[gate_polygon_mask(expr, gate), , drop = FALSE]
}

gate_polygon_mask <- function(expr, gate) {
    keep_all <- rep(TRUE, nrow(expr))
    verts <- gate_vertices_matrix(gate)
    if (is.null(verts) || nrow(verts) < 3) return(keep_all)
    x_channel <- as.character(gate_value(gate, "xChannel", ""))[1]
    y_channel <- as.character(gate_value(gate, "yChannel", ""))[1]
    if (!x_channel %in% colnames(expr) || !y_channel %in% colnames(expr)) return(keep_all)
    keep <- gate_point_in_polygon(expr[, x_channel], expr[, y_channel], verts)
    !is.na(keep) & keep
}

gate_positive_mask <- function(expr, gate, peak) {
    keep_all <- rep(TRUE, nrow(expr))
    if (is.null(gate) || !peak %in% colnames(expr)) return(keep_all)
    verts <- gate_vertices_matrix(gate)
    if (is.null(verts) || nrow(verts) == 0) return(keep_all)
    mode <- as.character(gate_value(gate, "mode", ""))[1]
    keep <- NULL
    if (identical(mode, "separator")) {
        keep <- expr[, peak] >= verts[1, 1]
    } else if (identical(mode, "positive_1d") && nrow(verts) >= 2) {
        lim <- range(verts[, 1], na.rm = TRUE)
        keep <- expr[, peak] >= lim[1] & expr[, peak] <= lim[2]
    } else if (nrow(verts) >= 3) {
        x_channel <- as.character(gate_value(gate, "xChannel", peak))[1]
        y_channel <- as.character(gate_value(gate, "yChannel", ""))[1]
        if (!x_channel %in% colnames(expr)) x_channel <- peak
        if (!y_channel %in% colnames(expr)) return(keep_all)
        keep <- gate_point_in_polygon(expr[, x_channel], expr[, y_channel], verts)
    }
    if (is.null(keep)) return(keep_all)
    !is.na(keep) & keep
}

gate_apply_positive_gate <- function(expr, gate, peak) {
    expr[gate_positive_mask(expr, gate, peak), , drop = FALSE]
}

gate_is_finalized_polygon <- function(gate) {
    if (is.null(gate) || identical(as.character(gate_value(gate, "mode", ""))[1], "blocked")) {
        return(FALSE)
    }
    verts <- gate_vertices_matrix(gate)
    !is.null(verts) && nrow(verts) >= 3 && all(is.finite(verts))
}

gate_is_finalized_histogram <- function(gate) {
    if (is.null(gate)) return(FALSE)
    mode <- as.character(gate_value(gate, "mode", ""))[1]
    if (!nzchar(mode) || mode %in% c("missing", "blocked")) return(FALSE)
    verts <- gate_vertices_matrix(gate)
    if (is.null(verts) || !all(is.finite(verts))) return(FALSE)
    if (identical(mode, "separator")) return(nrow(verts) >= 1L)
    if (mode %in% c("positive_1d", "negative_1d")) return(nrow(verts) >= 2L)
    nrow(verts) >= 3L
}

gate_spectrum_gates <- function(cache, filename, control_type, is_af = FALSE) {
    gates <- list(
        cell = gate_cached_gate(cache, "cell", filename, control_type),
        singlet = gate_cached_gate(cache, "singlet", filename, control_type),
        positive = gate_cached_gate(cache, "positive", filename, control_type)
    )
    gates$ready <- gate_is_finalized_polygon(gates$cell) &&
        gate_is_finalized_polygon(gates$singlet) &&
        (isTRUE(is_af) || gate_is_finalized_histogram(gates$positive))
    gates
}

gate_histogram_autogate_ranges <- function(peak_vals,
                                           sample_type,
                                           is_viability = FALSE,
                                           compute_fun = NULL,
                                           include_negative = TRUE) {
    peak_vals <- as.numeric(peak_vals)
    peak_vals <- peak_vals[is.finite(peak_vals)]
    if (length(peak_vals) < 10L) {
        stop("Too few gated events to auto-generate histogram gates.", call. = FALSE)
    }
    if (is.null(compute_fun)) {
        compute_fun <- get(".compute_reference_histogram_gate", envir = asNamespace("spectreasy"))
    }
    opts <- spectreasy::gating_options()
    result <- compute_fun(
        peak_vals = peak_vals,
        sample_type = sample_type,
        histogram_pct_beads = opts$histogram_pct_beads,
        histogram_direction_beads = opts$histogram_direction_beads,
        histogram_pct_cells = opts$histogram_pct_cells,
        histogram_direction_cells = opts$histogram_direction_cells,
        is_viability = isTRUE(is_viability)
    )
    positive <- sort(as.numeric(c(result$gate_min, result$gate_max)))
    negative <- sort(10^as.numeric(c(
        attr(result$vals_log, "neg_log_min"),
        attr(result$vals_log, "neg_log_max")
    )))
    if (length(positive) != 2L || any(!is.finite(positive)) || positive[2] <= positive[1]) {
        stop("The histogram autogater did not produce a valid positive interval.", call. = FALSE)
    }
    if (isTRUE(include_negative) && (!isTRUE(attr(result$vals_log, "negative_gate_present")) ||
        length(negative) != 2L || any(!is.finite(negative)) || negative[2] <= negative[1])) {
        stop("The histogram autogater did not produce a valid negative interval.", call. = FALSE)
    }
    repair_empty_interval <- function(limits, direction) {
        in_gate <- peak_vals >= limits[1] & peak_vals <= limits[2]
        minimum_count <- min(10L, length(peak_vals))
        if (sum(in_gate, na.rm = TRUE) >= minimum_count) return(limits)

        ordered <- sort(peak_vals)
        fallback_count <- min(
            length(ordered),
            max(minimum_count, as.integer(ceiling(length(ordered) * 0.05)))
        )
        fallback <- if (identical(direction, "high")) {
            utils::tail(ordered, fallback_count)
        } else {
            utils::head(ordered, fallback_count)
        }
        repaired <- range(fallback)
        if (repaired[1] == repaired[2]) {
            padding <- max(abs(repaired[1]) * 1e-8, .Machine$double.eps^0.5)
            repaired <- repaired + c(-padding, padding)
        }
        repaired
    }
    positive <- repair_empty_interval(positive, "high")
    if (isTRUE(include_negative)) {
        negative <- repair_empty_interval(negative, "low")
    } else {
        negative <- NULL
    }
    list(positive = positive, negative = negative)
}

gate_histogram_interval <- function(type, filename, peak, limits) {
    list(
        type = jsonlite::unbox(as.character(type)[1]),
        scope = jsonlite::unbox("file"),
        filename = jsonlite::unbox(as.character(filename)[1]),
        xChannel = jsonlite::unbox(as.character(peak)[1]),
        yChannel = jsonlite::unbox(""),
        mode = jsonlite::unbox(paste0(as.character(type)[1], "_1d")),
        vertices = list(
            list(x = jsonlite::unbox(unname(as.numeric(limits[1]))), y = jsonlite::unbox(0)),
            list(x = jsonlite::unbox(unname(as.numeric(limits[2]))), y = jsonlite::unbox(0))
        )
    )
}

gate_autogenerate_histograms <- function(gates) {
    df <- gate_read_mapping()
    uses_histogram <- if ("uses_histogram_gates" %in% colnames(df)) {
        as.logical(df$uses_histogram_gates)
    } else {
        !as.logical(df$is_af)
    }
    uses_histogram[is.na(uses_histogram)] <- FALSE
    targets <- df[uses_histogram, , drop = FALSE]
    if (nrow(targets) == 0) {
        stop("No non-AF single-color controls are available for histogram autogating.", call. = FALSE)
    }

    resolved <- lapply(seq_len(nrow(targets)), function(i) {
        row <- targets[i, , drop = FALSE]
        filename <- as.character(row$filename[1])
        control_type <- gate_normalize_control_type(row$control.type[1])
        cell_gate <- gate_cached_gate(gates, "cell", filename, control_type)
        singlet_gate <- gate_cached_gate(gates, "singlet", filename, control_type)
        list(
            row = row,
            filename = filename,
            control_type = control_type,
            cell_gate = cell_gate,
            singlet_gate = singlet_gate,
            needs_negative = isTRUE(as.logical(row$uses_negative_histogram_gate[1])),
            ready = gate_is_finalized_polygon(cell_gate) && gate_is_finalized_polygon(singlet_gate)
        )
    })
    missing <- vapply(resolved, function(item) if (item$ready) "" else item$filename, character(1))
    missing <- missing[nzchar(missing)]
    if (length(missing) > 0) {
        stop(
            "Complete FSC/SSC and singlet gates before histogram autogating: ",
            paste(missing, collapse = ", "),
            call. = FALSE
        )
    }

    generated <- list()
    spectra <- list()
    preserved_count <- 0L
    for (item in resolved) {
        positive_key <- paste0("positive:", item$filename)
        negative_key <- paste0("negative:", item$filename)
        keep_positive <- gate_is_finalized_histogram(gates[[positive_key]])
        keep_negative <- item$needs_negative && gate_is_finalized_histogram(gates[[negative_key]])
        preserved_count <- preserved_count + as.integer(keep_positive) + as.integer(keep_negative)

        path <- file.path(get_gate_scc_dir(), item$filename)
        ff <- gate_read_fcs(path)
        expr <- flowCore::exprs(ff)
        peak <- as.character(item$row$channel[1])
        if (!peak %in% colnames(expr)) {
            det_info <- gate_detector_info(ff)
            detectors <- if (is.null(det_info)) colnames(expr) else det_info$names
            detectors <- intersect(detectors, colnames(expr))
            if (length(detectors) == 0) {
                stop("Could not resolve a peak channel for ", item$filename, ".", call. = FALSE)
            }
            peak <- detectors[which.max(apply(expr[, detectors, drop = FALSE], 2, stats::var, na.rm = TRUE))]
        }
        scatter_mask <- gate_polygon_mask(expr, item$cell_gate) &
            gate_polygon_mask(expr, item$singlet_gate)
        if (!keep_positive || (item$needs_negative && !keep_negative)) {
            ranges <- gate_histogram_autogate_ranges(
                expr[scatter_mask, peak],
                sample_type = item$control_type,
                is_viability = isTRUE(item$row$is_viability[1]),
                include_negative = item$needs_negative
            )
            if (!keep_positive) {
                generated[[positive_key]] <- gate_histogram_interval(
                    "positive", item$filename, peak, ranges$positive
                )
            }
            if (item$needs_negative && !keep_negative) {
                generated[[negative_key]] <- gate_histogram_interval(
                    "negative", item$filename, peak, ranges$negative
                )
            }
        }
        positive_gate <- if (keep_positive) gates[[positive_key]] else generated[[positive_key]]
        spectrum_gates <- list(
            cell = item$cell_gate,
            singlet = item$singlet_gate,
            positive = positive_gate
        )
        spectrum <- gate_spectrum_data_from_loaded(
            expr = expr,
            ff = ff,
            spectrum_gates = spectrum_gates,
            peak = peak,
            is_af = FALSE,
            scatter_mask = scatter_mask
        )
        spectra[[item$filename]] <- gate_cache_spectrum(
            item$filename,
            spectrum_gates,
            spectrum
        )
    }
    list(
        gates = generated,
        spectra = spectra,
        files_processed = nrow(targets),
        gates_generated = length(generated),
        gates_preserved = preserved_count
    )
}

gate_detector_labels <- function(det_info, spec_channels) {
    labels <- det_info$labels[match(spec_channels, det_info$names)]
    labels[is.na(labels) | !nzchar(labels)] <- spec_channels[is.na(labels) | !nzchar(labels)]
    labels
}

gate_detector_signature <- function(pd) {
    desc <- if ("desc" %in% colnames(pd)) as.character(pd$desc) else rep("", nrow(pd))
    list(name = as.character(pd$name), desc = desc)
}

gate_detector_info <- function(ff, compute_fun = NULL) {
    pd <- flowCore::pData(flowCore::parameters(ff))
    signature <- gate_detector_signature(pd)
    cache <- getOption("spectreasy.gating_detector_cache", list())
    hit <- which(vapply(cache, function(entry) identical(entry$signature, signature), logical(1)))
    if (length(hit) > 0) return(cache[[hit[1]]]$info)
    if (is.null(compute_fun)) compute_fun <- spectreasy::get_sorted_detectors
    info <- tryCatch(compute_fun(pd), error = function(e) NULL)
    cache[[length(cache) + 1L]] <- list(signature = signature, info = info)
    options(spectreasy.gating_detector_cache = cache)
    info
}

gate_spectrum_counts_payload <- function(counts_mat) {
    connection <- rawConnection(raw(), open = "wb")
    on.exit(close(connection), add = TRUE)
    writeBin(as.integer(counts_mat), connection, size = 4L, endian = "little")
    list(
        format = jsonlite::unbox("uint32-column-major"),
        rows = jsonlite::unbox(nrow(counts_mat)),
        columns = jsonlite::unbox(ncol(counts_mat)),
        data = jsonlite::unbox(jsonlite::base64_enc(rawConnectionValue(connection)))
    )
}

gate_build_spectrum_data <- function(expr_filtered, spec_channels, det_info) {
    labels <- gate_detector_labels(det_info, spec_channels)
    y_power <- 1.5
    event_count <- nrow(expr_filtered)
    if (event_count == 0) {
        counts_mat <- matrix(0L, nrow = 150L, ncol = length(spec_channels))
        return(list(
            format = jsonlite::unbox("spectrum-histogram-v1"),
            channels = as.list(spec_channels),
            labels = as.list(labels),
            bin_mid = as.list(seq(0, 6, length.out = 150L)),
            bin_height = jsonlite::unbox(6 / 150),
            max_y = jsonlite::unbox(6),
            y_power = jsonlite::unbox(y_power),
            min_bin_count = jsonlite::unbox(1L),
            fill_limits = list(0, 1),
            event_count = jsonlite::unbox(0L),
            counts = gate_spectrum_counts_payload(counts_mat)
        ))
    }

    spectral_values <- if (identical(colnames(expr_filtered), spec_channels)) {
        as.matrix(expr_filtered)
    } else {
        as.matrix(expr_filtered[, spec_channels, drop = FALSE])
    }
    log_mat <- log10(pmax(spectral_values, 1e-3))
    finite <- log_mat[is.finite(log_mat)]
    if (length(finite) == 0) {
        return(gate_build_spectrum_data(expr_filtered[0, , drop = FALSE], spec_channels, det_info))
    }
    min_y <- floor(min(finite))
    max_y <- ceiling(max(finite))
    if (!is.finite(min_y)) min_y <- 0
    if (!is.finite(max_y)) max_y <- 6
    if (min_y == max_y) max_y <- min_y + 1
    breaks <- seq(min_y, max_y, length.out = 151L)
    bin_mid <- (breaks[-1] + breaks[-length(breaks)]) / 2
    bin_height <- breaks[2] - breaks[1]
    counts_mat <- vapply(seq_len(ncol(log_mat)), function(j) {
        values <- log_mat[, j]
        values <- values[is.finite(values)]
        as.integer(graphics::hist(values, breaks = breaks, plot = FALSE)$counts)
    }, integer(length(bin_mid)))
    if (is.null(dim(counts_mat))) counts_mat <- matrix(counts_mat, ncol = 1L)
    min_bin_count <- if (event_count <= 3000L) 1L else 3L
    visible_counts <- counts_mat[counts_mat >= min_bin_count]
    if (length(visible_counts) == 0) {
        fill_limits <- c(0, 1)
    } else {
        fills <- log10(visible_counts + 1)
        fill_limits <- c(min(fills), unname(stats::quantile(fills, 0.96, na.rm = TRUE)))
        if (!is.finite(fill_limits[2]) || fill_limits[2] <= fill_limits[1]) {
            fill_limits[2] <- fill_limits[1] + 1
        }
    }
    list(
        format = jsonlite::unbox("spectrum-histogram-v1"),
        channels = as.list(spec_channels),
        labels = as.list(labels),
        bin_mid = as.list(unname(bin_mid)),
        bin_height = jsonlite::unbox(unname(bin_height)),
        max_y = jsonlite::unbox(max(0, max_y)),
        y_power = jsonlite::unbox(y_power),
        min_bin_count = jsonlite::unbox(min_bin_count),
        fill_limits = as.list(unname(fill_limits)),
        event_count = jsonlite::unbox(event_count),
        counts = gate_spectrum_counts_payload(counts_mat)
    )
}

gate_spectrum_data_from_loaded <- function(expr, ff, spectrum_gates, peak, is_af = FALSE, scatter_mask = NULL) {
    det_info <- gate_detector_info(ff)
    if (is.null(det_info) || length(det_info$names) == 0) return(NULL)
    spec_channels <- intersect(det_info$names, colnames(expr))
    if (length(spec_channels) == 0) return(NULL)
    if (is.null(scatter_mask)) {
        scatter_mask <- gate_polygon_mask(expr, spectrum_gates$cell) &
            gate_polygon_mask(expr, spectrum_gates$singlet)
    }
    keep <- scatter_mask
    if (!isTRUE(is_af)) keep <- keep & gate_positive_mask(expr, spectrum_gates$positive, peak)
    spectral_expr <- expr[keep, spec_channels, drop = FALSE]
    gate_build_spectrum_data(spectral_expr, spec_channels, det_info)
}

gate_cache_spectrum <- function(filename, gates, spectrum) {
    cache <- getOption("spectreasy.gating_spectrum_cache", list())
    cache[[gate_safe_basename(filename)]] <- list(gates = gates, spectrum = spectrum)
    options(spectreasy.gating_spectrum_cache = cache)
    spectrum
}

gate_spectrum_for_file <- function(filename, gates = NULL) {
    df <- gate_read_mapping()
    name <- gate_safe_basename(filename)
    row <- df[df$filename == name, , drop = FALSE]
    if (nrow(row) == 0) stop("File is not present in mapping: ", name, call. = FALSE)
    if (is.null(gates)) {
        cache_data <- getOption("spectreasy.gating_state_cache")
        gates <- if (!is.null(cache_data)) cache_data$gates else NULL
    }
    control_type <- gate_normalize_control_type(row$control.type[1])
    spectrum_gates <- gate_spectrum_gates(
        gates,
        filename = name,
        control_type = control_type,
        is_af = isTRUE(row$is_af[1])
    )
    if (!isTRUE(spectrum_gates$ready)) return(NULL)

    spectrum_cache <- getOption("spectreasy.gating_spectrum_cache", list())
    gate_signature <- spectrum_gates[c("cell", "singlet", "positive")]
    cached <- spectrum_cache[[name]]
    if (!is.null(cached) && identical(cached$gates, gate_signature)) return(cached$spectrum)

    path <- file.path(get_gate_scc_dir(), name)
    ff <- gate_read_fcs(path)
    expr <- flowCore::exprs(ff)
    peak <- as.character(row$channel[1])
    if (!peak %in% colnames(expr)) {
        det_info <- gate_detector_info(ff)
        spec_channels <- if (is.null(det_info)) character() else intersect(det_info$names, colnames(expr))
        if (length(spec_channels) == 0) return(NULL)
        peak <- spec_channels[which.max(apply(expr[, spec_channels, drop = FALSE], 2, stats::var, na.rm = TRUE))]
    }
    spectrum <- gate_spectrum_data_from_loaded(
        expr = expr,
        ff = ff,
        spectrum_gates = spectrum_gates,
        peak = peak,
        is_af = isTRUE(row$is_af[1])
    )
    gate_cache_spectrum(name, gate_signature, spectrum)
}

gate_payload_for_file <- function(filename, max_points = 3000L, cache_result = TRUE) {
    name <- gate_safe_basename(filename)
    cache_key <- paste(name, max_points, sep = "::")
    cache <- getOption("spectreasy.gating_payload_cache", list())
    if (isTRUE(cache_result) && !is.null(cache[[cache_key]])) return(cache[[cache_key]])

    df <- gate_read_mapping()
    row <- df[df$filename == name, , drop = FALSE]
    if (nrow(row) == 0) stop("File is not present in mapping: ", name, call. = FALSE)

    path <- file.path(get_gate_scc_dir(), name)
    ff <- gate_read_fcs(path)
    expr <- flowCore::exprs(ff)
    meta <- gate_fcs_metadata(path)
    scatter_channels <- gate_scatter_channels(colnames(expr))
    peak <- as.character(row$channel[1])
    if (!peak %in% colnames(expr)) {
        det_info <- gate_detector_info(ff)
        det <- if (is.null(det_info)) colnames(expr) else det_info$names
        det <- intersect(det, colnames(expr))
        peak <- det[which.max(apply(expr[, det, drop = FALSE], 2, stats::var, na.rm = TRUE))]
    }

    n <- nrow(expr)
    max_points <- suppressWarnings(as.integer(max_points[1]))
    use_all <- is.finite(max_points) && !is.na(max_points) && max_points <= 0
    if (!use_all && (!is.finite(max_points) || is.na(max_points) || max_points < 500)) max_points <- 3000L
    if (!use_all && n > max_points) {
        seed <- sum(utf8ToInt(name)) %% .Machine$integer.max
        set.seed(seed)
        idx <- sort(sample.int(n, max_points))
    } else {
        idx <- seq_len(n)
    }
    out <- data.frame(
        event_index = idx,
        fsc_a = as.numeric(expr[idx, meta$fsc_a]),
        ssc_a = as.numeric(expr[idx, meta$ssc_a]),
        fsc_h = as.numeric(expr[idx, meta$fsc_h]),
        peak = as.numeric(expr[idx, peak]),
        stringsAsFactors = FALSE
    )
    for (channel in setdiff(scatter_channels, colnames(out))) {
        out[[channel]] <- as.numeric(expr[idx, channel])
    }
    out <- out[stats::complete.cases(out), , drop = FALSE]
    scatter_domains <- stats::setNames(lapply(scatter_channels, function(channel) gate_compute_extent(out[[channel]])), scatter_channels)
    payload <- list(
        file = row[1, , drop = FALSE],
        total_events = n,
        returned_events = nrow(out),
        channels = list(fsc_a = meta$fsc_a, fsc_h = meta$fsc_h, ssc_a = meta$ssc_a, peak = peak, scatter = scatter_channels),
        labels = meta$labels,
        domains = list(
            fsc_a = gate_compute_extent(out$fsc_a),
            ssc_a = gate_compute_extent(out$ssc_a),
            fsc_h = gate_compute_extent(out$fsc_h),
            peak = gate_compute_extent(out$peak),
            scatter = scatter_domains
        ),
        events = out
    )
    if (isTRUE(cache_result)) {
        cache[[cache_key]] <- payload
        options(spectreasy.gating_payload_cache = cache)
    }
    payload
}

gate_compact_payload <- function(payload) {
    events <- payload$events
    matrix_values <- as.matrix(events)
    storage.mode(matrix_values) <- "double"
    connection <- rawConnection(raw(), open = "wb")
    on.exit(close(connection), add = TRUE)
    writeBin(as.numeric(matrix_values), connection, size = 4L, endian = "little")
    encoded <- jsonlite::base64_enc(rawConnectionValue(connection))
    payload$events <- NULL
    payload$events_compact <- list(
        format = "float32-column-major",
        fields = colnames(matrix_values),
        rows = nrow(matrix_values),
        data = encoded
    )
    payload
}

#* @filter logger
function(req) {
    if (isTRUE(getOption("spectreasy.gui_request_log", FALSE))) {
        cat(as.character(Sys.time()), "-", req$REQUEST_METHOD, req$PATH_INFO, "\n")
    }
    plumber::forward()
}

gui_request_origin_allowed <- function(req) {
    origin <- req$HTTP_ORIGIN
    if (is.null(origin) || !nzchar(trimws(origin))) return(TRUE)
    allowed <- trimws(as.character(getOption("spectreasy.gui_allowed_origins", character())))
    allowed <- allowed[nzchar(allowed)]
    trimws(origin) %in% allowed
}

gui_api_token_value_allowed <- function(token) {
    expected <- as.character(getOption("spectreasy.gui_api_token", ""))[1]
    if (is.null(token)) token <- ""
    nzchar(expected) && identical(as.character(token)[1], expected)
}

gui_api_token_allowed <- function(req) {
    gui_api_token_value_allowed(req$HTTP_X_SPECTREASY_TOKEN)
}

#* @filter cors
function(req, res) {
    origin <- if (is.null(req$HTTP_ORIGIN)) "" else trimws(req$HTTP_ORIGIN)
    if (!gui_request_origin_allowed(req)) {
        res$status <- 403
        return(list(error = "This GUI origin is not authorized for the local Spectreasy session."))
    }
    mutating_method <- toupper(req$REQUEST_METHOD) %in% c("POST", "PUT", "PATCH", "DELETE")
    supplied_token <- req$HTTP_X_SPECTREASY_TOKEN
    validates_session <- identical(req$PATH_INFO, "/status") &&
        !is.null(supplied_token) && nzchar(trimws(as.character(supplied_token)[1]))
    if ((mutating_method || validates_session) && !gui_api_token_allowed(req)) {
        res$status <- 403
        return(list(error = "This request is not authorized for the active local Spectreasy session."))
    }
    if (nzchar(origin)) {
        res$setHeader("Access-Control-Allow-Origin", origin)
        res$setHeader("Vary", "Origin")
    }
    res$setHeader("Access-Control-Allow-Methods", "GET, POST, DELETE, OPTIONS")
    res$setHeader("Access-Control-Allow-Headers", "Content-Type, X-Spectreasy-Token")
    res$setHeader("Access-Control-Allow-Private-Network", "true")
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
    project_selected <- isTRUE(getOption("spectreasy.project_selected", TRUE))
    project_path <- if (project_selected) get_matrix_dir() else ""
    return(list(
        status = "ok",
        time = Sys.time(),
        wd = getwd(),
        matrix_dir = get_matrix_dir(),
        samples_dir = get_samples_dir(),
        unmixing_method = get_unmixing_method(),
        gui_mode = getOption("spectreasy.gui_mode", "tuner"),
        panel_cytometer = getOption("spectreasy.panel_cytometer", "aurora"),
        project_selected = project_selected,
        project_name = if (nzchar(project_path)) basename(project_path) else ""
    ))
}

#* Load persistent GUI state for one GUI module
#* @get /gui_state
#* @param module
function(module = "matrix_tuner") {
    path <- user_gui_config_path(module)
    if (!file.exists(path)) {
        return(list(module = normalize_gui_module(module), path = path, config = list()))
    }
    cfg <- tryCatch(jsonlite::fromJSON(path, simplifyVector = TRUE), error = function(e) list())
    list(module = normalize_gui_module(module), path = path, config = cfg)
}

#* Save persistent GUI state for one GUI module
#* @post /gui_state
function(req) {
    body <- jsonlite::fromJSON(req$postBody, simplifyVector = FALSE)
    module <- if (!is.null(body$module)) body$module else "matrix_tuner"
    cfg <- if (!is.null(body$config_json)) body$config_json else list()
    path <- user_gui_config_path(module)
    jsonlite::write_json(cfg, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
    list(success = TRUE, module = normalize_gui_module(module), path = path)
}

#* List SCC control files for manual gating
#* @get /gate_files
function() {
    df <- gate_read_mapping()
    first <- df$filename[df$file_exists][1]
    meta <- if (!is.na(first)) gate_fcs_metadata(file.path(get_gate_scc_dir(), first)) else list()
    list(files = df, metadata = meta, gate_file = get_gate_file())
}

#* Read the saved control mapping without synthesizing placeholder rows
#* @get /control_mapping
function() {
    path <- get_gate_control_file()
    if (!file.exists(path)) return(list(rows = data.frame(), exists = FALSE, path = path))
    rows <- tryCatch(
        utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
        error = function(e) NULL
    )
    if (is.null(rows)) return(list(rows = data.frame(), exists = TRUE, path = path, error = "The existing fcs_mapping.csv could not be read."))
    rows <- gui_annotate_active_af_mapping(rows)
    list(rows = rows, exists = TRUE, path = path)
}

#* Create fcs_mapping.csv from the active project's SCC folder
#* @post /control_mapping/create
function(req) {
    tryCatch({
        cytometer <- getOption("spectreasy.panel_cytometer", "auto")
        if (is.null(cytometer) || length(cytometer) == 0L || is.na(cytometer[1]) || !nzchar(trimws(as.character(cytometer[1])))) cytometer <- "auto"
        rows <- spectreasy::create_control_file(
            input_folder = get_gate_scc_dir(),
            cytometer = cytometer,
            unknown_fluor_policy = "by_channel",
            output_file = get_gate_control_file()
        )
        list(success = TRUE, path = get_gate_control_file(), rows = rows)
    }, error = function(e) list(success = FALSE, error = conditionMessage(e)))
}

#* Save the control mapping edited in the cockpit
#* @post /control_mapping
function(req) {
    body <- jsonlite::fromJSON(req$postBody, simplifyVector = FALSE)
    rows <- body$rows
    if (is.null(rows) || length(rows) == 0) {
        return(list(success = FALSE, error = "No control mapping rows were supplied."))
    }
    row_value <- function(row, keys, fallback = "") {
        for (key in keys) {
            value <- row[[key]]
            if (!is.null(value) && length(value) > 0 && !is.na(value[1]) && nzchar(as.character(value[1]))) {
                return(as.character(value[1]))
            }
        }
        fallback
    }
    mapping <- do.call(rbind, lapply(rows, function(row) {
        data.frame(
            filename = row_value(row, c("file", "filename")),
            fluorophore = row_value(row, "fluorophore"),
            marker = row_value(row, "marker"),
            channel = row_value(row, "channel"),
            control.type = if (grepl("bead", row_value(row, c("controlType", "control.type"), "cell"), ignore.case = TRUE)) "beads" else "cells",
            universal.negative = row_value(row, c("universalNegative", "universal.negative")),
            is.viability = row_value(row, c("isViability", "is.viability"), "FALSE"),
            stringsAsFactors = FALSE
        )
    }))
    path <- get_gate_control_file()
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(mapping, path, row.names = FALSE, quote = TRUE)
    list(success = TRUE, path = path, rows = nrow(mapping))
}

#* Load one downsampled SCC control payload
#* @get /gate_events
#* @param filename
#* @param max_points
function(filename, max_points = 3000) {
    gate_payload_for_file(filename = filename, max_points = max_points)
}

#* Preload all downsampled SCC control payloads
#* @get /gate_preload
#* @param max_points
function(max_points = 3000) {
    df <- gate_read_mapping()
    payloads <- lapply(df$filename, function(filename) {
        tryCatch(gate_payload_for_file(filename = filename, max_points = max_points), error = function(e) {
            list(error = conditionMessage(e), filename = filename)
        })
    })
    list(payloads = payloads, max_points = as.integer(max_points))
}

#* Preload all SCC controls in a compact float32 representation
#* @get /gate_preload_compact
#* @param max_points
function(max_points = 3000) {
    df <- gate_read_mapping()
    payloads <- lapply(df$filename, function(filename) {
        tryCatch({
            payload <- gate_payload_for_file(
                filename = filename,
                max_points = max_points,
                cache_result = FALSE
            )
            compact <- gate_compact_payload(payload)
            rm(payload)
            compact
        }, error = function(e) {
            list(error = conditionMessage(e), filename = filename)
        })
    })
    list(payloads = payloads, max_points = as.integer(max_points))
}

#* Auto-generate the required histogram gates for all non-AF controls
#* @post /gate_histogram_autogate
function(req) {
    tryCatch({
        body <- jsonlite::fromJSON(req$postBody, simplifyVector = FALSE)
        result <- gate_autogenerate_histograms(body$gates)
        list(
            success = TRUE,
            gates = result$gates,
            spectra = result$spectra,
            files_processed = result$files_processed,
            gates_generated = result$gates_generated,
            gates_preserved = result$gates_preserved
        )
    }, error = function(e) {
        list(success = FALSE, error = conditionMessage(e))
    })
}

#* Compute selected SCC spectrum bins for the current gate cache
#* @get /gate_spectrum
#* @param filename
#* @param dark
function(filename, dark = "false") {
    tryCatch(
        list(spectrum = gate_spectrum_for_file(filename = filename)),
        error = function(e) list(error = conditionMessage(e), spectrum = NULL)
    )
}

#* Compute and cache spectrum bins for one or more SCC controls using explicit gates
#* @post /gate_spectra
function(req) {
    tryCatch({
        body <- jsonlite::fromJSON(req$postBody, simplifyVector = FALSE)
        filenames <- unlist(body$filenames, use.names = FALSE)
        filenames <- as.character(filenames[nzchar(as.character(filenames))])
        gates <- body$gates
        spectra <- stats::setNames(lapply(filenames, function(filename) {
            gate_spectrum_for_file(
                filename = filename,
                gates = gates
            )
        }), filenames)
        list(success = TRUE, spectra = spectra)
    }, error = function(e) {
        list(success = FALSE, error = conditionMessage(e), spectra = list())
    })
}

#* List gate CSV files
#* @get /gate_configs
function() {
    dir.create(dirname(get_gate_file()), recursive = TRUE, showWarnings = FALSE)
    files <- list.files(dirname(get_gate_file()), pattern = "\\.csv$", full.names = FALSE, ignore.case = TRUE)
    list(configs = sort(files), config_dir = dirname(get_gate_file()), active = basename(get_gate_file()))
}

#* Load gate CSV
#* @get /gate_config
#* @param filename
function(filename = "") {
    path <- if (is.null(filename) || !nzchar(trimws(as.character(filename)[1]))) {
        get_gate_file()
    } else {
        file.path(dirname(get_gate_file()), basename(as.character(filename)[1]))
    }
    if (!file.exists(path)) {
        return(list(path = path, rows = gate_empty_config()))
    }
    rows <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    list(path = path, rows = rows)
}

#* Save gate CSV
#* @post /gate_config
function(req) {
    body <- jsonlite::fromJSON(req$postBody, simplifyVector = TRUE)
    written <- gate_write_config_csv(body$rows, get_gate_file())
    list(success = TRUE, path = written$path, rows = written$rows)
}

#* Save gate CSV with a system file picker
#* @post /gate_config_save_dialog
function(req) {
    body <- jsonlite::fromJSON(req$postBody, simplifyVector = TRUE)
    if (!gate_has_system_file_picker()) {
        return(list(success = FALSE, cancelled = FALSE, fallback = TRUE, message = "No backend system file picker is available."))
    }
    path <- gate_pick_csv_file(mode = "save")
    if (identical(path, "CANCEL")) {
        return(list(success = FALSE, cancelled = TRUE, fallback = FALSE, message = "Save cancelled."))
    }
    if (!nzchar(path)) {
        return(list(success = FALSE, cancelled = FALSE, fallback = TRUE, message = "System file picker failed. Falling back to browser picker."))
    }
    written <- gate_write_config_csv(body$rows, path)
    set_gate_file(written$path)
    list(success = TRUE, cancelled = FALSE, path = written$path, rows = written$rows)
}

#* Load gate CSV with a system file picker
#* @get /gate_config_load_dialog
function() {
    if (!gate_has_system_file_picker()) {
        return(list(success = FALSE, cancelled = FALSE, fallback = TRUE, message = "No backend system file picker is available."))
    }
    path <- gate_pick_csv_file(mode = "open")
    if (identical(path, "CANCEL")) {
        return(list(success = FALSE, cancelled = TRUE, fallback = FALSE, message = "Load cancelled."))
    }
    if (!nzchar(path)) {
        return(list(success = FALSE, cancelled = FALSE, fallback = TRUE, message = "System file picker failed. Falling back to browser picker."))
    }
    if (!file.exists(path)) {
        return(list(success = FALSE, cancelled = FALSE, fallback = FALSE, message = paste("File does not exist:", path)))
    }
    rows <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    set_gate_file(path)
    list(success = TRUE, cancelled = FALSE, path = path, rows = rows)
}

#* Load in-memory gate cache
#* @get /gate_cache
function() {
    cache <- getOption("spectreasy.gating_state_cache")
    if (is.null(cache)) {
        return(list(gates = list(), pointSize = 1.5, maxPoints = 50000, histogramBins = 100, histogramTransform = "asinh", viewSettings = list(), eventCountVersion = 2))
    }
    cache
}

#* Save in-memory gate cache
#* @post /gate_cache
function(req) {
    body <- jsonlite::fromJSON(req$postBody, simplifyVector = FALSE)
    options(spectreasy.gating_state_cache = list(
        gates = body$gates,
        pointSize = body$pointSize,
        maxPoints = body$maxPoints,
        histogramBins = body$histogramBins,
        histogramTransform = body$histogramTransform,
        viewSettings = body$viewSettings,
        eventCountVersion = body$eventCountVersion
    ))
    list(success = TRUE)
}

#* Shut down manual gating GUI
#* @post /gate_shutdown
function(req) {
    options(spectreasy.gui_shutdown_requested = TRUE)
    server <- getOption("spectreasy.gui_server", NULL)
    later::later(function() {
        if (!is.null(server)) {
            try(httpuv::stopServer(server), silent = TRUE)
        }
        try(httpuv::stopAllServers(), silent = TRUE)
    }, delay = 0.1)
    list(success = TRUE, message = "Gate config saved. Manual gating GUI is shutting down.")
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

#* List saved AF profiles
#* @get /af_profiles
function() {
    profiles <- tryCatch(spectreasy::list_af_profiles(), error = function(e) data.frame())
    if (nrow(profiles) > 0L) profiles$active <- profiles$name == gui_read_active_af_profile()
    list(profiles = profiles)
}

#* Return the detector-wise spectra for one saved AF profile
#* @get /af_profiles/data
#* @param name Profile name
function(name = "") {
    name <- trimws(as.character(name)[1])
    if (!nzchar(name)) return(list(error = "A profile name is required."))
    tryCatch({
        profile <- spectreasy::load_af_profile(name, show_plot = FALSE)$profile
        list(
            name = name,
            detectors = colnames(profile),
            spectra = lapply(seq_len(nrow(profile)), function(index) list(
                name = rownames(profile)[index],
                values = as.numeric(profile[index, , drop = TRUE])
            ))
        )
    }, error = function(e) list(error = conditionMessage(e)))
}

gui_pick_af_source_file <- function(initial_dir = file.path(get_matrix_dir(), "scc")) {
    initial_dir <- normalizePath(initial_dir, mustWork = FALSE)
    sysname <- Sys.info()[["sysname"]]
    if (identical(sysname, "Darwin")) {
        script <- paste0(
            "POSIX path of (choose file with prompt \"Select unstained FCS file\" default location POSIX file ",
            gate_applescript_quote(initial_dir), ")"
        )
        result <- suppressWarnings(system2("osascript", c("-e", shQuote(script)), stdout = TRUE, stderr = FALSE))
        if (!is.null(attr(result, "status")) || length(result) == 0L) return(NULL)
        return(normalizePath(trimws(result[[1]]), mustWork = TRUE))
    }
    if (.Platform$OS.type == "windows") {
        selected <- tryCatch(file.choose(new = FALSE), error = function(e) "")
        if (!nzchar(selected)) return(NULL)
        return(normalizePath(selected, mustWork = TRUE))
    }
    picker <- Sys.which("zenity")
    if (!nzchar(picker)) picker <- Sys.which("kdialog")
    if (!nzchar(picker)) stop("No graphical file picker is available. Install zenity or kdialog.", call. = FALSE)
    args <- if (grepl("zenity$", picker)) c("--file-selection", "--title=Select unstained FCS file", paste0("--filename=", initial_dir, "/"), "--file-filter=FCS files | *.fcs *.FCS") else c("--getopenfilename", initial_dir, "FCS files (*.fcs *.FCS)")
    result <- suppressWarnings(system2(picker, args, stdout = TRUE, stderr = FALSE))
    if (!is.null(attr(result, "status")) || length(result) == 0L) return(NULL)
    normalizePath(trimws(result[[1]]), mustWork = TRUE)
}

#* Open the native FCS picker for standalone AF extraction
#* @post /af_profiles/select-source
function() {
    tryCatch({
        selected <- gui_pick_af_source_file()
        if (is.null(selected) || !nzchar(selected)) return(list(success = FALSE, cancelled = TRUE))
        if (!grepl("\\.fcs$", selected, ignore.case = TRUE)) return(list(success = FALSE, cancelled = FALSE, error = "Select an FCS file."))
        list(success = TRUE, cancelled = FALSE, path = selected)
    }, error = function(e) list(success = FALSE, cancelled = FALSE, error = conditionMessage(e)))
}

#* Use a saved AF profile as the active dataset unstained control
#* @post /af_profiles/activate
function(req) {
    body <- jsonlite::fromJSON(req$postBody, simplifyVector = TRUE)
    name <- if (!is.null(body$profile_name)) trimws(as.character(body$profile_name)[1]) else ""
    if (!nzchar(name)) return(list(success = FALSE, error = "A profile name is required."))
    tryCatch({
        active <- gui_write_active_af_profile(name)
        options(
            spectreasy.gating_payload_cache = list(),
            spectreasy.gating_spectrum_cache = list(),
            spectreasy.gating_detector_cache = list()
        )
        list(success = TRUE, profile_name = active)
    }, error = function(e) list(success = FALSE, error = conditionMessage(e)))
}

#* Stop using a saved AF profile for the active dataset
#* @post /af_profiles/deactivate
function(req) {
    body <- jsonlite::fromJSON(req$postBody, simplifyVector = TRUE)
    name <- if (!is.null(body$profile_name)) trimws(as.character(body$profile_name)[1]) else ""
    if (!nzchar(name)) return(list(success = FALSE, error = "A profile name is required."))
    tryCatch({
        removed <- gui_unlink_active_af_profile(name)
        options(
            spectreasy.gating_payload_cache = list(),
            spectreasy.gating_spectrum_cache = list(),
            spectreasy.gating_detector_cache = list()
        )
        list(success = TRUE, profile_name = removed)
    }, error = function(e) list(success = FALSE, error = conditionMessage(e)))
}

#* Delete a saved AF profile
#* @delete /af_profiles/delete
#* @param name Profile name
function(name = "") {
    if (is.null(name) || !nzchar(trimws(as.character(name)[1]))) {
        return(list(success = FALSE, error = "A profile name is required."))
    }
    tryCatch({
        spectreasy::delete_af_profile(as.character(name)[1])
        if (identical(gui_read_active_af_profile(), as.character(name)[1])) {
            unlink(gui_active_af_config_path(), force = TRUE)
        }
        list(success = TRUE, name = as.character(name)[1])
    }, error = function(e) list(success = FALSE, error = conditionMessage(e)))
}

#* Apply a saved AF profile to a matrix
#* @post /af_profiles/apply
function(req) {
    body <- jsonlite::fromJSON(req$postBody, simplifyVector = TRUE)
    body_value <- function(name, fallback = "") {
        value <- body[[name]]
        if (is.null(value) || length(value) == 0 || is.na(value[1])) fallback else value[1]
    }
    matrix_filename <- trimws(as.character(body_value("matrix_filename"))[1])
    profile_name <- trimws(as.character(body_value("profile_name"))[1])
    output_filename <- trimws(as.character(body_value("output_filename", matrix_filename))[1])
    if (!nzchar(matrix_filename) || !nzchar(profile_name) || !nzchar(output_filename)) {
        return(list(success = FALSE, error = "Matrix, profile, and output names are required."))
    }
    tryCatch({
        matrix_file <- matrix_path(matrix_filename)
        output_file <- matrix_path(output_filename)
        if (!file.exists(matrix_file)) stop("Matrix file not found: ", matrix_filename, call. = FALSE)
        matrix_data <- read_matrix_csv(matrix_file)
        profile <- spectreasy::load_af_profile(profile_name, show_plot = FALSE)
        adjusted <- spectreasy::add_af_profile(matrix_data, profile, replace_existing = TRUE)
        dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
        utils::write.csv(adjusted, output_file, row.names = FALSE, quote = TRUE)
        list(success = TRUE, path = output_file, filename = output_filename)
    }, error = function(e) list(success = FALSE, error = conditionMessage(e)))
}

#* List available matrices
#* @get /matrices
function() {
    files <- list_matrix_csv_files()
    if (length(files) == 0) {
        return(character(0))
    }
    paths <- file.path(get_matrix_dir(), files)
    keep <- vapply(paths, is_probably_matrix_csv, logical(1))
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
    source_path <- if (!is.null(source_filename) && nzchar(trimws(as.character(source_filename)[1]))) {
        matrix_path(source_filename)
    } else {
        ""
    }
    if (nzchar(source_path) && identical(
        normalizePath(path, mustWork = FALSE),
        normalizePath(source_path, mustWork = FALSE)
    )) {
        stop("Adjusted matrix filename must differ from the source matrix filename.", call. = FALSE)
    }
    source_paths <- c(
        path,
        if (nzchar(source_path)) source_path else character()
    )
    df <- merge_hidden_af_rows(df, source_paths = source_paths)

    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
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
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
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
function(matrix_json, raw_data_json, type = "reference", matrix_filename = "", method = "") {
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
        method_resolved <- tryCatch(
            spectreasy:::.normalize_unmix_method(if (!nzchar(trimws(as.character(method)[1]))) get_unmixing_method() else method),
            error = function(e) e
        )
        if (inherits(method_resolved, "error")) {
            return(list(error = conditionMessage(method_resolved)))
        }
        ff <- flowCore::flowFrame(Y_sub)
        res <- tryCatch(
            spectreasy::calc_residuals(ff, mat_sub, method = method_resolved),
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

# Workflow wrappers keep the browser orchestration layer small. All scientific
# work remains delegated to the exported Spectreasy functions below.
gui_workflow_body <- function(req) {
    tryCatch(
        jsonlite::fromJSON(req$postBody, simplifyVector = FALSE),
        error = function(e) stop("Invalid JSON request body: ", conditionMessage(e), call. = FALSE)
    )
}

gui_workflow_value <- function(body, key, fallback = NULL) {
    value <- body[[key]]
    if (is.null(value) || length(value) == 0 || !nzchar(trimws(as.character(value[1])))) {
        return(fallback)
    }
    as.character(value[1])
}

gui_workflow_bool <- function(body, key, fallback = FALSE) {
    value <- body[[key]]
    if (is.null(value) || length(value) == 0 || is.na(value[1])) return(isTRUE(fallback))
    if (is.logical(value[1])) return(isTRUE(value[1]))
    normalized <- tolower(trimws(as.character(value[1])))
    if (normalized %in% c("true", "1", "yes", "on")) return(TRUE)
    if (normalized %in% c("false", "0", "no", "off")) return(FALSE)
    isTRUE(fallback)
}

gui_workflow_number <- function(body, key, fallback = 0, integer = FALSE, minimum = NULL, maximum = NULL) {
    supplied <- !is.null(body[[key]]) && length(body[[key]]) > 0L
    value <- suppressWarnings(as.numeric(gui_workflow_value(body, key, fallback)))
    if (length(value) == 0 || !is.finite(value[1])) {
        if (supplied) stop("Invalid numeric value for '", key, "'.", call. = FALSE)
        value <- fallback
    }
    value <- as.numeric(value[1])
    if (!is.null(minimum) && value < minimum) {
        stop("'", key, "' must be at least ", minimum, ".", call. = FALSE)
    }
    if (!is.null(maximum) && value > maximum) {
        stop("'", key, "' must be at most ", maximum, ".", call. = FALSE)
    }
    if (isTRUE(integer)) {
        if (!isTRUE(all.equal(value, round(value)))) stop("'", key, "' must be an integer.", call. = FALSE)
        value <- as.integer(value)
    }
    value
}

gui_workflow_path <- function(body, key, fallback = "", allow_empty = TRUE) {
    value <- gui_workflow_value(body, key, fallback)
    if (is.null(value) || is.na(value)) return(if (isTRUE(allow_empty)) "" else fallback)
    value <- trimws(as.character(value[1]))
    if (!nzchar(value) && isTRUE(allow_empty)) return("")
    value
}

gui_method_optional_args <- function(method, body) {
    resolved <- tryCatch(spectreasy:::.normalize_unmix_method(method), error = function(e) "")
    if (identical(resolved, "Spectreasy")) {
        return(list(spectreasy_weight_quantile = gui_workflow_number(body, "spectreasy_weight_quantile", 0.65, minimum = 0, maximum = 1)))
    }
    list()
}

gui_workflow_root <- function(body) {
    root <- gui_workflow_value(body, "projectPath", get_matrix_dir())
    if (!dir.exists(root)) stop("Project folder not found: ", root, call. = FALSE)
    normalizePath(root, mustWork = TRUE)
}

gui_restore_working_directory <- function(path) {
    if (is.null(path) || length(path) != 1L || is.na(path) || !nzchar(path) || !dir.exists(path)) {
        return(invisible(FALSE))
    }
    tryCatch(
        {
            setwd(path)
            invisible(TRUE)
        },
        error = function(e) invisible(FALSE)
    )
}

gui_workflow_run <- function(body, action, expr) {
    root <- gui_workflow_root(body)
    gui_set_project_context(root)
    old_wd <- getwd()
    on.exit(gui_restore_working_directory(old_wd), add = TRUE)
    setwd(root)
    messages <- character()
    warnings <- character()
    tryCatch(
        {
            result <- NULL
            output <- capture.output(
                result <- withCallingHandlers(
                    force(expr),
                    message = function(condition) {
                        messages <<- c(messages, conditionMessage(condition))
                    },
                    warning = function(condition) {
                        warnings <<- c(warnings, paste0("Warning: ", conditionMessage(condition)))
                    }
                ),
                type = "output"
            )
            if (length(output)) cat(paste0(output, collapse = "\n"), "\n")
            list(
                success = TRUE,
                action = action,
                project_path = root,
                finished_at = as.character(Sys.time()),
                logs = c(output, messages, warnings),
                result = result
            )
        },
        error = function(e) {
            list(
                success = FALSE,
                action = action,
                project_path = root,
                logs = c(messages, warnings),
                error = conditionMessage(e),
                suggested_next_step = "Review the project inputs and the full action log, then retry from the relevant workflow card."
            )
        }
    )
}

gui_workflow_file_or_null <- function(path, root = getwd()) {
    if (is.null(path) || length(path) == 0 || is.na(path[1])) return(NULL)
    path <- trimws(as.character(path[1]))
    if (!nzchar(path)) return(NULL)
    candidate <- if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) path else file.path(root, path)
    if (!file.exists(candidate)) return(NULL)
    normalizePath(candidate, mustWork = TRUE)
}

gui_active_af_config_path <- function(root = get_matrix_dir()) {
    file.path(normalizePath(root, mustWork = FALSE), ".spectreasy", "active_af_profile.json")
}

gui_read_active_af_profile <- function(root = get_matrix_dir()) {
    path <- gui_active_af_config_path(root)
    if (!file.exists(path)) return("")
    value <- tryCatch(jsonlite::fromJSON(path, simplifyVector = TRUE), error = function(e) NULL)
    name <- if (is.list(value)) value$profile_name else NULL
    if (is.null(name) || length(name) == 0L || is.na(name[1])) return("")
    name <- trimws(as.character(name[1]))
    profile_exists <- nzchar(name) && tryCatch(
        file.exists(spectreasy:::.af_profile_file(name, create_dir = FALSE)),
        error = function(e) FALSE
    )
    if (!profile_exists) {
        unlink(path, force = TRUE)
        return("")
    }
    name
}

gui_write_active_af_profile <- function(name, root = get_matrix_dir()) {
    name <- trimws(as.character(name)[1])
    if (!nzchar(name)) stop("A saved AF profile name is required.", call. = FALSE)
    profile <- spectreasy::load_af_profile(name, show_plot = FALSE)
    scc_dir <- file.path(root, "scc")
    fcs_files <- if (dir.exists(scc_dir)) list.files(scc_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE) else character()
    fcs_files <- fcs_files[!startsWith(basename(fcs_files), "._")]
    if (length(fcs_files) > 0L) {
        detector_names <- spectreasy:::.prepare_reference_detector_info(fcs_files[1])$detector_names
        profile_detectors <- colnames(profile$profile)
        if (!setequal(detector_names, profile_detectors)) {
            stop("Saved AF profile detectors do not match this dataset's SCC detector set.", call. = FALSE)
        }
    }
    path <- gui_active_af_config_path(root)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile("active_af_profile_", tmpdir = dirname(path), fileext = ".json")
    jsonlite::write_json(
        list(profile_name = name, updated = as.character(Sys.time())),
        tmp,
        auto_unbox = TRUE,
        pretty = TRUE
    )
    if (!file.rename(tmp, path)) {
        if (!file.copy(tmp, path, overwrite = TRUE)) stop("Could not save the active AF profile selection.", call. = FALSE)
        unlink(tmp)
    }
    name
}

gui_unlink_active_af_profile <- function(name, root = get_matrix_dir()) {
    name <- trimws(as.character(name)[1])
    if (!nzchar(name)) stop("A saved AF profile name is required.", call. = FALSE)
    active <- gui_read_active_af_profile(root)
    if (!nzchar(active)) return(name)
    if (!identical(active, name)) {
        stop(name, " is not linked to this dataset.", call. = FALSE)
    }
    path <- gui_active_af_config_path(root)
    if (file.exists(path)) {
        status <- unlink(path, force = TRUE)
        if (!identical(status, 0L)) stop("Could not unlink the active AF profile.", call. = FALSE)
    }
    name
}

gui_primary_unstained_rows <- function(rows) {
    if (is.null(rows) || !is.data.frame(rows) || nrow(rows) == 0L) return(logical(0))
    fluorophore <- if ("fluorophore" %in% colnames(rows)) as.character(rows$fluorophore) else rep("", nrow(rows))
    marker <- if ("marker" %in% colnames(rows)) as.character(rows$marker) else rep("", nrow(rows))
    filename <- if ("filename" %in% colnames(rows)) as.character(rows$filename) else rep("", nrow(rows))
    control_type <- if ("control.type" %in% colnames(rows)) tolower(trimws(as.character(rows$control.type))) else rep("cells", nrow(rows))
    is_af <- grepl("^AF($|_|\\b)", trimws(fluorophore), ignore.case = TRUE) |
        grepl("autofluorescence|unstained", trimws(marker), ignore.case = TRUE) |
        grepl("unstained", basename(filename), ignore.case = TRUE)
    is_dead <- grepl("dead", paste(fluorophore, marker, filename), ignore.case = TRUE)
    is_bead <- grepl("bead", control_type, ignore.case = TRUE) | grepl("bead", basename(filename), ignore.case = TRUE)
    is_af & !is_dead & !is_bead
}

gui_annotate_active_af_mapping <- function(rows, root = get_matrix_dir()) {
    rows <- as.data.frame(rows, stringsAsFactors = FALSE, check.names = FALSE)
    active <- gui_read_active_af_profile(root)
    ignored <- rep(FALSE, nrow(rows))
    if (nzchar(active)) ignored <- gui_primary_unstained_rows(rows)
    rows$ignored <- ignored
    rows$ignored_reason <- ifelse(
        ignored,
        paste0("Ignored because saved AF profile '", active, "' is used as the unstained cell control."),
        ""
    )
    rows
}

gui_filter_active_af_mapping <- function(rows, root = get_matrix_dir()) {
    if (!nzchar(gui_read_active_af_profile(root)) || is.null(rows) || nrow(rows) == 0L) return(rows)
    ignored <- gui_primary_unstained_rows(rows)
    ignored_files <- if ("filename" %in% colnames(rows)) as.character(rows$filename[ignored]) else character()
    rows <- rows[!ignored, , drop = FALSE]
    if ("universal.negative" %in% colnames(rows)) {
        refs <- trimws(as.character(rows$universal.negative))
        rows$universal.negative[refs %in% c("AF", ignored_files)] <- ""
    }
    rows
}

gui_filtered_control_file <- function(root = get_matrix_dir()) {
    source <- file.path(root, "fcs_mapping.csv")
    if (!file.exists(source)) return(source)
    rows <- utils::read.csv(source, stringsAsFactors = FALSE, check.names = FALSE)
    filtered <- gui_filter_active_af_mapping(rows, root = root)
    if (nrow(filtered) == nrow(rows)) return(source)
    dir.create(file.path(root, ".spectreasy"), recursive = TRUE, showWarnings = FALSE)
    target <- tempfile("fcs_mapping_active_af_", tmpdir = file.path(root, ".spectreasy"), fileext = ".csv")
    utils::write.csv(filtered, target, row.names = FALSE, quote = TRUE)
    target
}

gui_project_scan <- function(root) {
    files <- if (dir.exists(root)) list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE) else character()
    files <- normalizePath(files[file.exists(files)], mustWork = FALSE)
    relative <- if (length(files) > 0) {
        sub(paste0("^", gsub("([.|(){}+*?^$\\[\\]\\\\])", "\\\\\\1", root), "[/\\\\]?"), "", files)
    } else {
        character()
    }
    relative <- gsub("\\\\", "/", relative)
    count_matches <- function(pattern) sum(grepl(pattern, relative, ignore.case = TRUE, perl = TRUE))
    controls <- count_matches("(^|/)(scc|controls?)/.*\\.fcs$")
    samples <- count_matches("(^|/)samples?/.*\\.fcs$")
    matrices <- count_matches("matrix.*\\.csv$|unmixing.*\\.csv$|detector_noise.*\\.csv$")
    reports <- count_matches("\\.(html?|pdf)$")
    gates <- count_matches("gate.*\\.csv$")
    qc_metrics <- count_matches("metric.*\\.csv$")
    spectral_variants <- count_matches("variant.*\\.(rds|csv)$")
    summary <- if (length(files) == 0) {
        "empty project"
    } else if (matrices > 0 && reports > 0 && samples > 0) {
        "mixed/partial project"
    } else if (matrices > 0) {
        "reference built"
    } else if (controls > 0) {
        "controls imported"
    } else {
        "project detected"
    }
    list(
        project_path = root,
        missing_input_dirs = basename(gui_missing_project_input_dirs(root)),
        files = relative,
        scan = list(
            controls = controls,
            samples = samples,
            matrices = matrices,
            reports = reports,
            gates = gates,
            qc_metrics = qc_metrics,
            spectral_variants = spectral_variants
        ),
        summary = summary,
        recommended_next_action = if (matrices == 0) "Review controls and build a reference matrix" else if (samples > 0) "Run sample unmixing" else "Import samples"
    )
}

gui_missing_project_input_dirs <- function(project_path) {
    project_path <- normalizePath(project_path, mustWork = TRUE)
    paths <- file.path(project_path, c("scc", "samples"))
    paths[!dir.exists(paths)]
}

gui_ensure_project_input_dirs <- function(project_path) {
    project_path <- normalizePath(project_path, mustWork = TRUE)
    paths <- file.path(project_path, c("scc", "samples"))
    for (path in paths) {
        if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
            stop("Could not create project input folder: ", path, call. = FALSE)
        }
    }
    invisible(paths)
}

gui_project_file_location <- function(kind, filename = NULL) {
    kind <- tolower(trimws(as.character(kind)[1]))
    folder <- switch(kind, controls = "scc", samples = "samples", NULL)
    if (is.null(folder)) stop("File kind must be 'controls' or 'samples'.", call. = FALSE)
    root <- normalizePath(get_matrix_dir(), mustWork = TRUE)
    directory <- file.path(root, folder)
    if (is.null(filename)) return(list(kind = kind, folder = folder, directory = directory))
    candidate <- gsub("\\\\", "/", trimws(as.character(filename)[1]))
    safe_name <- basename(candidate)
    if (!nzchar(safe_name) || !identical(candidate, safe_name) || safe_name %in% c(".", "..")) {
        stop("Invalid project filename.", call. = FALSE)
    }
    if (!grepl("\\.fcs$", safe_name, ignore.case = TRUE)) {
        stop("Only FCS files can be managed here.", call. = FALSE)
    }
    list(kind = kind, folder = folder, directory = directory, filename = safe_name, path = file.path(directory, safe_name))
}

gui_project_file_rows <- function(kind) {
    location <- gui_project_file_location(kind)
    if (!dir.exists(location$directory)) return(data.frame())
    files <- list.files(location$directory, pattern = "\\.fcs$", ignore.case = TRUE, full.names = TRUE)
    if (!length(files)) return(data.frame())
    info <- file.info(files)
    rows <- data.frame(
        name = basename(files),
        size = as.numeric(info$size),
        modified = format(info$mtime, "%Y-%m-%dT%H:%M:%S%z"),
        modified_epoch = as.numeric(info$mtime),
        kind = location$kind,
        stringsAsFactors = FALSE
    )
    rows[order(gui_natural_sort_key(rows$name), tolower(rows$name)), , drop = FALSE]
}

.gui_project_uploads <- new.env(parent = emptyenv())

gui_project_upload_id <- function() {
    paste(sample(c(letters, LETTERS, 0:9), 32L, replace = TRUE), collapse = "")
}

gui_project_upload_session <- function(upload_id) {
    upload_id <- trimws(as.character(upload_id)[1])
    if (!nzchar(upload_id) || !exists(upload_id, envir = .gui_project_uploads, inherits = FALSE)) {
        stop("Upload session not found or expired.", call. = FALSE)
    }
    get(upload_id, envir = .gui_project_uploads, inherits = FALSE)
}

gui_project_upload_discard <- function(upload_id) {
    upload_id <- trimws(as.character(upload_id)[1])
    if (nzchar(upload_id) && exists(upload_id, envir = .gui_project_uploads, inherits = FALSE)) {
        session <- get(upload_id, envir = .gui_project_uploads, inherits = FALSE)
        unlink(session$temporary, force = TRUE)
        rm(list = upload_id, envir = .gui_project_uploads)
    }
    invisible(NULL)
}

gui_project_upload_discard_stale <- function(max_age_seconds = 3600) {
    upload_ids <- ls(envir = .gui_project_uploads, all.names = TRUE)
    now <- Sys.time()
    for (upload_id in upload_ids) {
        session <- get(upload_id, envir = .gui_project_uploads, inherits = FALSE)
        age <- suppressWarnings(as.numeric(difftime(now, session$created, units = "secs")))
        if (!is.finite(age) || age > max_age_seconds) gui_project_upload_discard(upload_id)
    }
    invisible(NULL)
}

gui_validate_fcs_upload <- function(path) {
    size <- file.info(path)$size
    if (!is.finite(size) || size < 6) stop("Uploaded FCS file is empty or truncated.", call. = FALSE)
    connection <- file(path, open = "rb")
    on.exit(close(connection), add = TRUE)
    header <- rawToChar(readBin(connection, what = "raw", n = 6L))
    if (!grepl("^FCS[0-9]\\.[0-9]$", header)) {
        stop("Uploaded file does not contain a valid FCS header.", call. = FALSE)
    }
    invisible(TRUE)
}

gui_natural_sort_key <- function(x) {
    vapply(as.character(x), function(value) {
        parts <- regmatches(value, gregexpr("[0-9]+|[^0-9]+", value, perl = TRUE))[[1]]
        paste(vapply(parts, function(part) {
            if (grepl("^[0-9]+$", part)) sprintf("%020.0f", as.numeric(part)) else tolower(part)
        }, character(1)), collapse = "")
    }, character(1))
}

gui_set_project_context <- function(project_path) {
    if (!dir.exists(project_path)) stop("Project folder not found: ", project_path, call. = FALSE)
    project_path <- normalizePath(project_path, mustWork = TRUE)
    options(
        spectreasy.project_dir = project_path,
        spectreasy.matrix_dir = project_path,
        spectreasy.samples_dir = file.path(project_path, "samples"),
        spectreasy.gating_scc_dir = file.path(project_path, "scc"),
        spectreasy.gating_control_file = file.path(project_path, "fcs_mapping.csv"),
        spectreasy.gating_gate_file = file.path(project_path, "ssc_gate_config.csv"),
        spectreasy.gating_payload_cache = list(),
        spectreasy.gating_spectrum_cache = list(),
        spectreasy.gating_detector_cache = list()
    )
    options(spectreasy.project_selected = TRUE)
    project_path
}

gui_pick_project_directory <- function(initial_dir = path.expand("~"), allow_create = FALSE) {
    initial_dir <- normalizePath(initial_dir, mustWork = FALSE)
    prompt <- if (isTRUE(allow_create)) "Create a Spectreasy project folder" else "Open a Spectreasy project folder"
    if (identical(Sys.info()[["sysname"]], "Darwin")) {
        script <- paste0(
            "POSIX path of (choose folder with prompt ", gate_applescript_quote(prompt), " default location POSIX file ",
            gate_applescript_quote(initial_dir), ")"
        )
        result <- suppressWarnings(system2("osascript", c("-e", shQuote(script)), stdout = TRUE, stderr = FALSE))
        if (!is.null(attr(result, "status")) || length(result) == 0L) return(NULL)
        return(sub("/$", "", trimws(result[[1]])))
    }
    if (.Platform$OS.type == "windows") {
        result <- utils::choose.dir(default = initial_dir, caption = prompt)
        if (is.na(result)) return(NULL)
        return(result)
    }
    picker <- Sys.which("zenity")
    if (!nzchar(picker)) picker <- Sys.which("kdialog")
    if (!nzchar(picker)) stop("No graphical folder picker is available. Install zenity or kdialog.", call. = FALSE)
    args <- if (grepl("zenity$", picker)) {
        c("--file-selection", "--directory", paste0("--title=", prompt), paste0("--filename=", initial_dir, "/"))
    } else {
        c("--getexistingdirectory", initial_dir, "--title", prompt)
    }
    result <- suppressWarnings(system2(picker, args, stdout = TRUE, stderr = FALSE))
    if (!is.null(attr(result, "status")) || length(result) == 0L) return(NULL)
    trimws(result[[1]])
}

#* Project scan and workflow prerequisites
#* @get /project/status
#* @param project_path Optional project directory
function(project_path = "") {
    if (!isTRUE(getOption("spectreasy.project_selected", TRUE)) &&
        (is.null(project_path) || !nzchar(trimws(as.character(project_path)[1])))) {
        return(list(
            project_path = "",
            files = character(),
            scan = list(controls = 0, samples = 0, matrices = 0, reports = 0, gates = 0, qc_metrics = 0, spectral_variants = 0),
            summary = "no project selected",
            recommended_next_action = "Choose a project folder"
        ))
    }
    root <- if (is.null(project_path) || !nzchar(trimws(as.character(project_path)[1]))) get_matrix_dir() else as.character(project_path)[1]
    root <- normalizePath(root, mustWork = FALSE)
    gui_project_scan(root)
}

gui_project_relative_path <- function(path, root) {
    root_pattern <- gsub("([.|(){}+*?^$\\[\\]\\\\])", "\\\\\\1", normalizePath(root, mustWork = FALSE))
    gsub("\\\\", "/", sub(paste0("^", root_pattern, "[/\\\\]?"), "", path))
}

gui_project_report_files <- function(root, files = NULL, output_root = "", report_type = "") {
    root <- normalizePath(root, mustWork = FALSE)
    output_root <- trimws(as.character(output_root)[1])
    report_type <- tolower(trimws(as.character(report_type)[1]))
    if (nzchar(output_root) && report_type %in% c("control", "sample")) {
        output_base <- if (grepl("^(/|[A-Za-z]:[/\\\\])", output_root)) output_root else file.path(root, output_root)
        output_base <- normalizePath(output_base, mustWork = FALSE)
        root_prefix <- paste0(root, .Platform$file.sep)
        if (!identical(output_base, root) && !startsWith(output_base, root_prefix)) {
            stop("Report output folder is outside the active project.", call. = FALSE)
        }
        stage_dir <- file.path(output_base, if (identical(report_type, "sample")) "unmix_samples" else "unmix_controls")
        report_dir_name <- if (identical(report_type, "sample")) "qc_samples" else "qc_controls"
        report_dirs <- if (dir.exists(stage_dir)) {
            candidates <- list.dirs(stage_dir, recursive = FALSE, full.names = TRUE)
            candidates[grepl(paste0("^", report_dir_name, "(?:_(?:[2-9]|[1-9][0-9]+))?$"), basename(candidates), perl = TRUE)]
        } else character()
        report_name <- if (identical(report_type, "sample")) "qc_samples_report" else "qc_controls_report"
        return(unique(unlist(lapply(report_dirs, function(report_dir) {
            list.files(
                report_dir,
                pattern = paste0("^", report_name, "\\.(html?|pdf)$"),
                full.names = TRUE,
                ignore.case = TRUE
            )
        }), use.names = FALSE)))
    }
    if (is.null(files)) {
        top_level <- if (dir.exists(root)) list.dirs(root, recursive = FALSE, full.names = TRUE) else character()
        report_roots <- top_level[
            grepl("^(reports|spectreasy_outputs(?:_[^/]+)?)$", basename(top_level), ignore.case = TRUE, perl = TRUE)
        ]
        files <- unique(unlist(lapply(report_roots, function(directory) {
            list.files(directory, recursive = TRUE, full.names = TRUE, all.files = FALSE)
        }), use.names = FALSE))
    }
    relative_files <- gui_project_relative_path(files, root)
    files[
        grepl("\\.(html?|pdf)$", relative_files, ignore.case = TRUE) &
            grepl("(^|/)(reports|spectreasy_outputs(?:_[^/]+)?)/", relative_files, ignore.case = TRUE, perl = TRUE)
    ]
}

gui_project_report_source_files <- function(root, sample = FALSE, output_root = "") {
    input_dir <- file.path(root, if (isTRUE(sample)) "samples" else "scc")
    input_files <- if (dir.exists(input_dir)) {
        list.files(input_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
    } else {
        character()
    }
    shared <- file.path(root, c("fcs_mapping.csv", "ssc_gate_config.csv"))
    shared <- shared[file.exists(shared)]
    output_root <- trimws(as.character(output_root)[1])
    output_roots <- if (nzchar(output_root)) {
        candidate <- if (grepl("^(/|[A-Za-z]:[/\\\\])", output_root)) output_root else file.path(root, output_root)
        candidate <- normalizePath(candidate, mustWork = FALSE)
        if (dir.exists(candidate)) candidate else character()
    } else {
        candidates <- list.dirs(root, recursive = FALSE, full.names = TRUE)
        candidates[grepl("^spectreasy_outputs(?:_[^/]+)?$", basename(candidates), ignore.case = TRUE, perl = TRUE)]
    }
    artifact_pattern <- if (isTRUE(sample)) {
        "(reference_matrix|detector_noise|spectral_variant).*\\.(csv|rds)$"
    } else {
        "(reference_matrix|detector_noise).*\\.csv$"
    }
    artifacts <- unique(unlist(lapply(output_roots, function(directory) {
        list.files(directory, pattern = artifact_pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    }), use.names = FALSE))
    unique(c(input_files, shared, artifacts))
}

gui_project_report_type <- function(relative_path) {
    relative_path <- gsub("\\\\", "/", as.character(relative_path)[1])
    if (grepl("panel", relative_path, ignore.case = TRUE)) return("Panel overview")
    if (grepl("(^|/)(unmix_samples|qc_samples)(/|$)|qc_samples_report", relative_path, ignore.case = TRUE, perl = TRUE)) {
        return("Sample QC")
    }
    "Control QC"
}

gui_resolve_project_report <- function(path, root = get_matrix_dir()) {
    root <- normalizePath(root, mustWork = FALSE)
    relative <- gsub("\\\\", "/", trimws(as.character(path)[1]))
    parts <- strsplit(relative, "/", fixed = TRUE)[[1]]
    if (!nzchar(relative) || any(parts %in% c("", ".", ".."))) {
        stop("Invalid project report path.", call. = FALSE)
    }
    target <- normalizePath(file.path(root, do.call(file.path, as.list(parts))), mustWork = FALSE)
    root_prefix <- paste0(root, .Platform$file.sep)
    if (!identical(target, root) && !startsWith(target, root_prefix)) {
        stop("Project report is outside the active project.", call. = FALSE)
    }
    if (!file.exists(target) || dir.exists(target)) {
        stop("Project report not found.", call. = FALSE)
    }
    if (!tolower(tools::file_ext(target)) %in% c("html", "htm")) {
        stop("Only an existing HTML report can be exported to PDF.", call. = FALSE)
    }
    target
}

gui_find_chromium <- function() {
    candidates <- unique(c(
        unname(Sys.which(c("google-chrome", "google-chrome-stable", "chromium", "chromium-browser", "chrome"))),
        "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
        "/Applications/Chromium.app/Contents/MacOS/Chromium",
        file.path(Sys.getenv("PROGRAMFILES"), "Google", "Chrome", "Application", "chrome.exe"),
        file.path(Sys.getenv("PROGRAMFILES(X86)"), "Google", "Chrome", "Application", "chrome.exe")
    ))
    candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
    if (length(candidates)) candidates[[1]] else ""
}

gui_export_html_report_pdf <- function(html_file, output_file = tempfile(fileext = ".pdf")) {
    browser <- gui_find_chromium()
    if (!nzchar(browser)) {
        stop("Chrome or Chromium is required to export the existing HTML report as PDF.", call. = FALSE)
    }
    if (!requireNamespace("chromote", quietly = TRUE)) {
        stop("The optional 'chromote' package is required to export the existing HTML report as PDF.", call. = FALSE)
    }
    html_file <- normalizePath(html_file, mustWork = TRUE)
    output_file <- normalizePath(output_file, mustWork = FALSE)
    previous_browser <- Sys.getenv("CHROMOTE_CHROME", unset = NA_character_)
    Sys.setenv(CHROMOTE_CHROME = browser)
    on.exit({
        if (is.na(previous_browser)) Sys.unsetenv("CHROMOTE_CHROME") else Sys.setenv(CHROMOTE_CHROME = previous_browser)
    }, add = TRUE)
    session <- chromote::ChromoteSession$new()
    on.exit(try(session$close(), silent = TRUE), add = TRUE)
    session$Page$navigate(utils::URLencode(paste0("file://", html_file), reserved = FALSE))
    for (attempt in seq_len(40L)) {
        state <- tryCatch(
            session$Runtime$evaluate("document.readyState", returnByValue = TRUE)$result$value,
            error = function(e) ""
        )
        if (identical(state, "complete")) break
        Sys.sleep(0.1)
    }
    Sys.sleep(0.4)
    pdf <- session$Page$printToPDF(printBackground = TRUE, preferCSSPageSize = TRUE)
    base64enc::base64decode(pdf$data, output = output_file)
    if (!file.exists(output_file) || file.info(output_file)$size < 5L) {
        stop("Chrome could not export the HTML report as PDF.", call. = FALSE)
    }
    output_file
}

#* Discover report artifacts and compare them with upstream project files
#* @get /project/reports
function(project_path = "", output_root = "", report_type = "") {
    root <- if (is.null(project_path) || !nzchar(trimws(as.character(project_path)[1]))) get_matrix_dir() else as.character(project_path)[1]
    root <- normalizePath(root, mustWork = FALSE)
    reports <- gui_project_report_files(root, output_root = output_root, report_type = report_type)
    if (length(reports) == 0L) return(list(reports = data.frame()))
    relative <- function(path) gui_project_relative_path(path, root)
    control_sources <- gui_project_report_source_files(root, sample = FALSE, output_root = output_root)
    sample_sources <- gui_project_report_source_files(root, sample = TRUE, output_root = output_root)
    rows <- lapply(reports, function(report) {
        relative_report <- relative(report)
        report_type <- gui_project_report_type(relative_report)
        is_sample <- identical(report_type, "Sample QC")
        source_pattern <- if (is_sample) {
            "(^|/)(samples?/.*\\.fcs$|.*reference_matrix.*\\.csv$|.*detector_noise.*\\.csv$|.*spectral_variant.*\\.rds$)"
        } else {
            "(^|/)(scc/.*\\.fcs$|fcs_mapping\\.csv$|.*gate.*\\.csv$|.*reference_matrix.*\\.csv$|.*detector_noise.*\\.csv$)"
        }
        source_files <- if (is_sample) sample_sources else control_sources
        source_files <- source_files[grepl(source_pattern, relative(source_files), ignore.case = TRUE, perl = TRUE)]
        source_files <- setdiff(source_files, report)
        stale <- length(source_files) > 0L && any(file.info(source_files)$mtime > file.info(report)$mtime)
        data.frame(
            path = relative_report,
            report_type = report_type,
            format = if (grepl("\\.pdf$", report, ignore.case = TRUE)) "PDF" else "HTML",
            created = format(file.info(report)$mtime, "%Y-%m-%dT%H:%M:%S%z"),
            created_epoch = as.numeric(file.info(report)$mtime),
            status = if (stale) "stale" else "current",
            stringsAsFactors = FALSE
        )
    })
    list(reports = do.call(rbind, rows))
}

#* Stream a project artifact through the local R backend
#* @get /project/file
#* @param path Relative path inside the active project
#* @param project_path Active cockpit project directory
#* @param token Active cockpit session token
function(path = "", project_path = "", token = "", req, res) {
    if (!gui_api_token_value_allowed(token) && !gui_api_token_allowed(req)) {
        res$status <- 403
        return(list(error = "This project artifact belongs to a different or inactive Spectreasy session."))
    }
    root <- if (is.null(project_path) || !nzchar(trimws(as.character(project_path)[1]))) get_matrix_dir() else as.character(project_path)[1]
    root <- normalizePath(root, mustWork = FALSE)
    relative <- gsub("\\\\", "/", trimws(as.character(path)[1]))
    parts <- strsplit(relative, "/", fixed = TRUE)[[1]]
    if (!nzchar(relative) || any(parts %in% c("", ".", ".."))) {
        res$status <- 400
        return(list(error = "Invalid project artifact path."))
    }
    target <- normalizePath(file.path(root, do.call(file.path, as.list(parts))), mustWork = FALSE)
    root_prefix <- paste0(root, .Platform$file.sep)
    if (!identical(target, root) && !startsWith(target, root_prefix)) {
        res$status <- 403
        return(list(error = "Artifact path is outside the active project."))
    }
    if (!file.exists(target) || dir.exists(target)) {
        res$status <- 404
        return(list(error = "Project artifact not found."))
    }
    extension <- tolower(tools::file_ext(target))
    content_type <- switch(
        extension,
        html = "text/html; charset=utf-8",
        htm = "text/html; charset=utf-8",
        pdf = "application/pdf",
        csv = "text/csv; charset=utf-8",
        json = "application/json; charset=utf-8",
        "application/octet-stream"
    )
    res$setHeader("Content-Type", content_type)
    res$body <- readBin(target, what = "raw", n = file.info(target)$size)
    res
}

#* Export an existing HTML QC report to PDF without rerunning QC
#* @post /project/report/export-pdf
function(req) {
    body <- gui_workflow_body(req)
    tryCatch({
        root <- gui_workflow_root(body)
        report_file <- gui_resolve_project_report(gui_workflow_value(body, "path", ""), root = root)
        pdf_file <- gui_export_html_report_pdf(report_file)
        on.exit(unlink(pdf_file, force = TRUE), add = TRUE)
        list(
            success = TRUE,
            filename = paste0(tools::file_path_sans_ext(basename(report_file)), ".pdf"),
            content_type = "application/pdf",
            content_base64 = base64enc::base64encode(pdf_file)
        )
    }, error = function(e) list(success = FALSE, error = conditionMessage(e)))
}

#* Change the active project folder from the cockpit
#* @post /project/context
function(req) {
    body <- gui_workflow_body(req)
    project_path <- gui_workflow_path(body, "projectPath", "", allow_empty = FALSE)
    tryCatch({
        project_path <- gui_set_project_context(project_path)
        list(success = TRUE, project = gui_project_scan(project_path))
    }, error = function(e) list(success = FALSE, error = conditionMessage(e)))
}

#* Create the standard input folders after explicit cockpit confirmation
#* @post /project/initialize
function(req) {
    body <- gui_workflow_body(req)
    tryCatch({
        root <- gui_workflow_root(body)
        created <- basename(gui_missing_project_input_dirs(root))
        gui_ensure_project_input_dirs(root)
        list(success = TRUE, created = created, project = gui_project_scan(root))
    }, error = function(e) list(success = FALSE, error = conditionMessage(e)))
}

#* Open the native folder picker and activate the selected project
#* @post /project/select
function(req) {
    tryCatch({
        selected <- gui_pick_project_directory(get_matrix_dir(), allow_create = FALSE)
        if (is.null(selected) || !nzchar(selected)) return(list(success = FALSE, cancelled = TRUE))
        selected <- gui_set_project_context(selected)
        list(success = TRUE, cancelled = FALSE, project = gui_project_scan(selected))
    }, error = function(e) list(success = FALSE, cancelled = FALSE, error = conditionMessage(e)))
}

#* Create and activate a project folder with the native folder picker
#* @post /project/create
function(req) {
    tryCatch({
        selected <- gui_pick_project_directory(get_matrix_dir(), allow_create = TRUE)
        if (is.null(selected) || !nzchar(selected)) return(list(success = FALSE, cancelled = TRUE))
        selected <- gui_set_project_context(selected)
        list(success = TRUE, cancelled = FALSE, project = gui_project_scan(selected))
    }, error = function(e) list(success = FALSE, cancelled = FALSE, error = conditionMessage(e)))
}

#* List FCS files in the active project's controls or samples folder
#* @get /project/files
#* @param kind Either controls or samples
function(kind = "controls") {
    tryCatch(
        list(success = TRUE, files = gui_project_file_rows(kind)),
        error = function(e) list(success = FALSE, error = conditionMessage(e), files = data.frame())
    )
}

#* Begin a bounded-memory FCS upload
#* @post /project/upload-start
function(req) {
    body <- gui_workflow_body(req)
    tryCatch({
        gui_project_upload_discard_stale()
        location <- gui_project_file_location(
            gui_workflow_value(body, "kind", ""),
            gui_workflow_value(body, "filename", "")
        )
        if (!dir.exists(location$directory)) {
            stop("Create the project's scc and samples folders before adding files.", call. = FALSE)
        }
        size <- gui_workflow_number(body, "size", 0, integer = TRUE, minimum = 1)
        maximum <- as.numeric(getOption("spectreasy.gui_max_upload_bytes", 50 * 1024^3))
        if (is.finite(maximum) && size > maximum) {
            stop("Uploaded file exceeds the configured size limit.", call. = FALSE)
        }
        if (file.exists(location$path)) stop("A file named '", location$filename, "' already exists.", call. = FALSE)
        upload_id <- gui_project_upload_id()
        temporary <- tempfile("spectreasy-upload-", tmpdir = location$directory)
        if (!file.create(temporary)) stop("Could not prepare the upload destination.", call. = FALSE)
        assign(upload_id, list(
            location = location,
            temporary = temporary,
            size = size,
            written = 0,
            created = Sys.time()
        ), envir = .gui_project_uploads)
        list(success = TRUE, upload_id = upload_id)
    }, error = function(e) list(success = FALSE, error = conditionMessage(e)))
}

#* Append one encoded chunk to an FCS upload
#* @post /project/upload-chunk
function(req) {
    body <- gui_workflow_body(req)
    upload_id <- gui_workflow_value(body, "upload_id", "")
    tryCatch({
        session <- gui_project_upload_session(upload_id)
        offset <- gui_workflow_number(body, "offset", 0, integer = TRUE, minimum = 0)
        if (!identical(as.numeric(offset), as.numeric(session$written))) {
            stop("Upload chunk is out of sequence; restart the upload.", call. = FALSE)
        }
        content <- gui_workflow_value(body, "content_base64", "")
        payload <- tryCatch(jsonlite::base64_dec(content), error = function(e) NULL)
        if (is.null(payload) || !length(payload)) stop("Upload chunk is empty or invalid.", call. = FALSE)
        if (session$written + length(payload) > session$size) stop("Upload exceeds the declared file size.", call. = FALSE)
        connection <- file(session$temporary, open = "ab")
        on.exit(try(close(connection), silent = TRUE), add = TRUE)
        writeBin(payload, connection)
        close(connection)
        session$written <- session$written + length(payload)
        assign(upload_id, session, envir = .gui_project_uploads)
        list(success = TRUE, written = session$written)
    }, error = function(e) {
        gui_project_upload_discard(upload_id)
        list(success = FALSE, error = conditionMessage(e))
    })
}

#* Validate and commit an FCS upload
#* @post /project/upload-finish
function(req) {
    body <- gui_workflow_body(req)
    upload_id <- gui_workflow_value(body, "upload_id", "")
    tryCatch({
        session <- gui_project_upload_session(upload_id)
        if (!identical(as.numeric(session$written), as.numeric(session$size))) {
            stop("Upload is incomplete.", call. = FALSE)
        }
        gui_validate_fcs_upload(session$temporary)
        if (file.exists(session$location$path)) stop("The destination file now exists; the upload was not overwritten.", call. = FALSE)
        if (!file.rename(session$temporary, session$location$path)) stop("Could not save uploaded file.", call. = FALSE)
        rm(list = upload_id, envir = .gui_project_uploads)
        rows <- gui_project_file_rows(session$location$kind)
        row <- rows[rows$name == session$location$filename, , drop = FALSE]
        list(success = TRUE, file = row)
    }, error = function(e) {
        gui_project_upload_discard(upload_id)
        list(success = FALSE, error = conditionMessage(e))
    })
}

#* Cancel an incomplete FCS upload
#* @post /project/upload-abort
function(req) {
    body <- gui_workflow_body(req)
    upload_id <- gui_workflow_value(body, "upload_id", "")
    gui_project_upload_discard(upload_id)
    list(success = TRUE)
}

#* Delete one FCS file from the active project's controls or samples folder
#* @delete /project/files
#* @param kind Either controls or samples
#* @param filename FCS filename inside that folder
function(kind = "", filename = "") {
    tryCatch({
        location <- gui_project_file_location(kind, filename)
        if (!file.exists(location$path) || dir.exists(location$path)) stop("Project file not found.", call. = FALSE)
        removed <- unlink(location$path, force = TRUE)
        if (!identical(removed, 0L)) stop("Could not delete project file.", call. = FALSE)
        list(success = TRUE, filename = location$filename)
    }, error = function(e) list(success = FALSE, error = conditionMessage(e)))
}

#* Delete all FCS files from one active project input folder
#* @delete /project/files/all
#* @param kind Either controls or samples
function(kind = "") {
    tryCatch({
        location <- gui_project_file_location(kind)
        if (!dir.exists(location$directory)) {
            return(list(success = TRUE, deleted = 0L))
        }
        files <- list.files(location$directory, pattern = "\\.fcs$", ignore.case = TRUE, full.names = TRUE)
        if (length(files)) {
            removed <- unlink(files, force = TRUE)
            if (!identical(removed, 0L)) stop("Could not delete all project files.", call. = FALSE)
        }
        list(success = TRUE, deleted = length(files))
    }, error = function(e) list(success = FALSE, error = conditionMessage(e)))
}

#* CORS preflight for workflow_control
#* @options /workflow/control
function(res) {
    return("")
}

#* Run the complete control-stage workflow through Spectreasy R
#* @post /workflow/control
function(req) {
    body <- gui_workflow_body(req)
    root <- gui_workflow_root(body)
    scc_dir <- gui_workflow_value(body, "scc_dir", "scc")
    control_file <- gui_workflow_value(body, "control_file", "fcs_mapping.csv")
    active_af_profile <- gui_read_active_af_profile(root)
    if (nzchar(active_af_profile)) {
        control_file <- gui_filtered_control_file(root)
        if (!identical(normalizePath(control_file, mustWork = FALSE), normalizePath(file.path(root, "fcs_mapping.csv"), mustWork = FALSE))) {
            on.exit(unlink(control_file, force = TRUE), add = TRUE)
        }
    }
    output_dir <- gui_workflow_value(body, "output_dir", "spectreasy_outputs")
    method <- gui_workflow_value(body, "method", get_unmixing_method())
    cytometer <- gui_workflow_value(body, "cytometer", "auto")
    gate_file <- gui_workflow_file_or_null(
        gui_workflow_value(body, "manual_gate_file", ""),
        root = gui_workflow_root(body)
    )
    control_args <- c(
        list(
            scc_dir = scc_dir,
            control_file = control_file,
            auto_create_mapping = gui_workflow_bool(body, "auto_create_mapping", TRUE),
            cytometer = cytometer,
            auto_unknown_fluor_policy = gui_workflow_value(body, "auto_unknown_fluor_policy", "by_channel"),
            output_dir = output_dir,
            unmixing_method = method,
            unmix_scatter_panel_size_mm = gui_workflow_number(body, "unmix_scatter_panel_size_mm", 30, minimum = 1),
            seed = gui_workflow_number(body, "seed", 1, integer = TRUE, minimum = 1),
            af_n_bands = gui_workflow_number(body, "af_n_bands", 100, integer = TRUE, minimum = 1),
            af_max_cells = gui_workflow_number(body, "af_max_cells", 50000, integer = TRUE, minimum = 1),
            default_sample_type = gui_workflow_value(body, "default_sample_type", "beads"),
            histogram_pct_beads = gui_workflow_number(body, "histogram_pct_beads", 0.98, minimum = 0, maximum = 1),
            histogram_direction_beads = gui_workflow_value(body, "histogram_direction_beads", "right"),
            histogram_pct_cells = gui_workflow_number(body, "histogram_pct_cells", 0.35, minimum = 0, maximum = 1),
            histogram_direction_cells = gui_workflow_value(body, "histogram_direction_cells", "right"),
            outlier_percentile = gui_workflow_number(body, "outlier_percentile", 0.02, minimum = 0, maximum = 1),
            debris_percentile = gui_workflow_number(body, "debris_percentile", 0.08, minimum = 0, maximum = 1),
            bead_gate_scale = gui_workflow_number(body, "bead_gate_scale", 1.3, minimum = 0),
            max_clusters = gui_workflow_number(body, "max_clusters", 10, integer = TRUE, minimum = 1),
            min_cluster_proportion = gui_workflow_number(body, "min_cluster_proportion", 0.03, minimum = 0, maximum = 1),
            gate_contour_beads = gui_workflow_number(body, "gate_contour_beads", 0.95, minimum = 0, maximum = 1),
            gate_contour_cells = gui_workflow_number(body, "gate_contour_cells", 0.9, minimum = 0, maximum = 1),
            subsample_n = gui_workflow_number(body, "subsample_n", 5000, integer = TRUE, minimum = 1),
            rwls_max_iter = gui_workflow_number(body, "rwls_max_iter", 1, integer = TRUE, minimum = 1),
            n_threads = gui_workflow_number(body, "n_threads", 1, integer = TRUE, minimum = 1),
            save_qc_png = gui_workflow_bool(body, "save_qc_png", TRUE),
            save_report = gui_workflow_bool(body, "save_report", TRUE),
            report_format = gui_workflow_value(body, "report_format", "html"),
            gating_mode = gui_workflow_value(
                body,
                "gating_mode",
                "interactive"
            ),
            manual_gate_file = gate_file,
            scc_background_method = gui_workflow_value(body, "scc_background_method", "scatter_knn"),
            scc_background_k = gui_workflow_number(body, "scc_background_k", 2, integer = TRUE, minimum = 1),
            spectral_variant_som_nodes = gui_workflow_number(body, "spectral_variant_som_nodes", 16, integer = TRUE, minimum = 1),
            spectral_variant_top_k = gui_workflow_number(body, "spectral_variant_top_k", 3, integer = TRUE, minimum = 1),
            spectral_variant_cosine_threshold = gui_workflow_number(body, "spectral_variant_cosine_threshold", 0.98, minimum = 0, maximum = 1),
            spectral_variant_max_variants = gui_workflow_number(body, "spectral_variant_max_variants", 8, integer = TRUE, minimum = 1),
            spectral_variant_min_events = gui_workflow_number(body, "spectral_variant_min_events", 50, integer = TRUE, minimum = 1),
            autospectral_n_candidates = gui_workflow_number(body, "autospectral_n_candidates", 1000, integer = TRUE, minimum = 1),
            autospectral_n_spectral = gui_workflow_number(body, "autospectral_n_spectral", 200, integer = TRUE, minimum = 1),
            autospectral_min_events = gui_workflow_number(body, "autospectral_min_events", 10, integer = TRUE, minimum = 1),
            autospectral_refine = gui_workflow_bool(body, "autospectral_refine", FALSE)
        ),
        if (nzchar(active_af_profile)) list(af_profile = active_af_profile) else list(),
        gui_method_optional_args(method, body)
    )
    run <- gui_workflow_run(
        body,
        "control",
        do.call(spectreasy::unmix_controls, control_args)
    )
    if (!isTRUE(run$success)) return(run)
    result <- run$result
    run$result <- list(
        reference_matrix_file = result$reference_matrix_file,
        detector_noise_file = result$detector_noise_file,
        unmixing_matrix_file = result$unmixing_matrix_file,
        spectral_variant_library_file = result$spectral_variant_library_file,
        qc_report_file = result$qc_report_file,
        spectra_file = result$spectra_file,
        unmixing_scatter_file = result$unmixing_scatter_file
    )
    run
}

#* CORS preflight for workflow_sample
#* @options /workflow/sample
function(res) {
    return("")
}

#* Run sample unmixing and sample QC through Spectreasy R
#* @post /workflow/sample
function(req) {
    body <- gui_workflow_body(req)
    sample_dir <- gui_workflow_value(body, "sample_dir", "samples")
    matrix_file <- gui_workflow_value(body, "matrix_file", file.path("spectreasy_outputs", "unmix_controls", "scc_reference_matrix.csv"))
    noise_file <- gui_workflow_file_or_null(gui_workflow_value(body, "detector_noise_file", ""), root = gui_workflow_root(body))
    output_dir <- gui_workflow_value(body, "output_dir", "spectreasy_outputs")
    method <- gui_workflow_value(body, "method", get_unmixing_method())
    sample_args <- c(
        list(
            sample_dir = sample_dir,
            unmixing_matrix_file = matrix_file,
            detector_noise_file = noise_file,
            unmixing_method = method,
            rwls_max_iter = gui_workflow_number(body, "rwls_max_iter", 1, integer = TRUE, minimum = 1),
            n_threads = gui_workflow_number(body, "n_threads", 1, integer = TRUE, minimum = 1),
            spectral_variant_top_k = gui_workflow_number(body, "spectral_variant_top_k", 3, integer = TRUE, minimum = 1),
            spectral_variant_min_abundance = gui_workflow_number(body, "spectral_variant_min_abundance", 1, minimum = 0),
            spectral_variant_positive_fraction = gui_workflow_number(body, "spectral_variant_positive_fraction", 0.02, minimum = 0, maximum = 1),
            spectral_variant_min_improvement = gui_workflow_number(body, "spectral_variant_min_improvement", 0.01, minimum = 0),
            spectral_variant_library_file = gui_workflow_file_or_null(gui_workflow_value(body, "spectral_variant_library_file", ""), root = gui_workflow_root(body)),
            estimate_af = gui_workflow_bool(body, "estimate_af", FALSE),
            output_dir = output_dir,
            write_fcs = gui_workflow_bool(body, "write_fcs", TRUE),
            save_report = gui_workflow_bool(body, "save_report", TRUE),
            report_format = gui_workflow_value(body, "report_format", "html"),
            report_per_sample = gui_workflow_bool(body, "report_per_sample", FALSE),
            save_qc_plots = gui_workflow_bool(body, "save_qc_plots", TRUE),
            plot_n_events = gui_workflow_number(body, "plot_n_events", 10000, integer = TRUE, minimum = 1),
            chunk_size = gui_workflow_number(body, "chunk_size", 50000, integer = TRUE, minimum = 1),
            seed = gui_workflow_number(body, "seed", 1, integer = TRUE, minimum = 1),
            return_type = gui_workflow_value(body, "return_type", "list"),
            verbose = FALSE
        ),
        gui_method_optional_args(method, body)
    )
    run <- gui_workflow_run(
        body,
        "sample",
        do.call(spectreasy::unmix_samples, sample_args)
    )
    if (!isTRUE(run$success)) return(run)
    result <- run$result
    run$result <- list(
        qc_report_file = attr(result, "qc_report_file"),
        qc_samples_dir = attr(result, "qc_samples_dir"),
        qc_metrics_dir = attr(result, "qc_metrics_dir"),
        output_dir = output_dir
    )
    run
}

#* CORS preflight for workflow_report
#* @options /workflow/report
function(res) {
    return("")
}

#* Compare selected unmixing methods on one active sample
#* @post /workflow/compare
function(req) {
    body <- gui_workflow_body(req)
    raw_methods <- body$methods
    methods <- unique(as.character(unlist(raw_methods, recursive = TRUE, use.names = FALSE)))
    methods <- methods[nzchar(methods)]
    if (length(methods) == 0) methods <- c("Spectreasy", "OLS", "WLS", "NNLS")
    matrix_input <- gui_workflow_value(body, "matrix_file", "")
    sample_dir_input <- gui_workflow_value(body, "sample_dir", "samples")
    root <- gui_workflow_root(body)
    matrix_file <- gui_workflow_file_or_null(matrix_input, root = root)
    sample_dir <- if (grepl("^(/|[A-Za-z]:[/\\\\])", sample_dir_input)) sample_dir_input else file.path(root, sample_dir_input)
    if (is.null(matrix_file) || !file.exists(matrix_file)) {
        return(list(success = FALSE, error = "Select a readable reference matrix before comparing methods."))
    }
    sample_files <- sort(list.files(sample_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE))
    if (length(sample_files) == 0) {
        return(list(success = FALSE, error = paste("No FCS sample files found in", sample_dir)))
    }
    matrix_df <- read_matrix_csv(matrix_file)
    if (ncol(matrix_df) < 2) return(list(success = FALSE, error = "The selected matrix has no detector columns."))
    marker_names <- as.character(matrix_df[[1]])
    matrix_values <- as.matrix(matrix_df[, -1, drop = FALSE])
    rownames(matrix_values) <- marker_names
    sample_frame <- flowCore::read.FCS(sample_files[1], transformation = FALSE, truncate_max_range = FALSE)
    events <- flowCore::exprs(sample_frame)
    if (nrow(events) > 2000) {
        set.seed(gui_workflow_number(body, "seed", 1, integer = TRUE, minimum = 1))
        events <- events[sort(sample.int(nrow(events), 2000)), , drop = FALSE]
    }
    common_detectors <- intersect(colnames(events), colnames(matrix_values))
    if (length(common_detectors) == 0) return(list(success = FALSE, error = "The sample and matrix share no detector channels."))
    comparison <- lapply(methods, function(method) {
        resolved <- tryCatch(spectreasy:::.normalize_unmix_method(method), error = function(e) e)
        if (inherits(resolved, "error")) return(data.frame(method = method, status = "error", residual_rms = NA_real_, message = conditionMessage(resolved), stringsAsFactors = FALSE))
        result <- tryCatch(
            spectreasy::calc_residuals(
                flowCore::flowFrame(events[, common_detectors, drop = FALSE]),
                matrix_values[, common_detectors, drop = FALSE],
                method = resolved
            ),
            error = function(e) e
        )
        if (inherits(result, "error")) return(data.frame(method = method, status = "error", residual_rms = NA_real_, message = conditionMessage(result), stringsAsFactors = FALSE))
        numeric_result <- suppressWarnings(as.matrix(result))
        data.frame(method = method, status = "complete", residual_rms = sqrt(mean(numeric_result^2, na.rm = TRUE)), message = "", stringsAsFactors = FALSE)
    })
    comparison_df <- do.call(rbind, comparison)
    output_dir <- file.path(root, "spectreasy_outputs", "method_comparison")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    output_file <- file.path(output_dir, "method_comparison.csv")
    utils::write.csv(comparison_df, output_file, row.names = FALSE, quote = TRUE)
    list(success = TRUE, result = list(sample = basename(sample_files[1]), matrix = matrix_file, output_file = output_file, rows = comparison_df))
}

#* Render a report from existing Spectreasy outputs
#* @post /workflow/report
function(req) {
    body <- gui_workflow_body(req)
    report_type <- tolower(gui_workflow_value(body, "report_type", "control"))
    report_format <- tolower(gui_workflow_value(body, "report_format", "html"))
    if (!report_format %in% c("html", "pdf")) {
        return(list(success = FALSE, error = "Report format must be 'html' or 'pdf'."))
    }
    overwrite <- tolower(gui_workflow_value(body, "overwrite", "overwrite"))
    if (!overwrite %in% c("version", "overwrite", "error")) {
        return(list(success = FALSE, error = "Overwrite behavior must be 'version', 'overwrite', or 'error'."))
    }
    root <- gui_workflow_root(body)
    run <- gui_workflow_run(
        body,
        "report",
        if (identical(report_type, "sample")) {
            stop("Sample report rendering requires the current sample results object. Run sample unmixing from the Samples workspace first.")
        } else {
            matrix_file <- file.path(root, "spectreasy_outputs", "unmix_controls", "scc_reference_matrix.csv")
            if (!file.exists(matrix_file)) stop("Reference matrix not found: ", matrix_file)
            M <- spectreasy:::.read_unmixing_matrix_csv(matrix_file)
            report_dir <- file.path(root, "spectreasy_outputs", "unmix_controls", "qc_controls")
            dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
            report_file <- file.path(report_dir, paste0("qc_controls_report.", report_format))
            spectreasy::qc_controls(
                M = M,
                unmixing_matrix_file = matrix_file,
                scc_dir = file.path(root, "scc"),
                control_file = file.path(root, "fcs_mapping.csv"),
                output_file = report_file,
                report_format = report_format,
                overwrite = overwrite
            )
        }
    )
    if (!isTRUE(run$success)) return(run)
    report_path <- if (is.list(run$result) && !is.null(run$result$output_file)) run$result$output_file else report_file
    run$result <- list(report_file = report_path, report_format = report_format)
    run
}

#* CORS preflight for workflow_af
#* @options /workflow/af
function(res) {
    return("")
}

gui_default_af_profile_name <- function(fcs_file, existing_names = character()) {
    stem <- tools::file_path_sans_ext(basename(as.character(fcs_file)[1]))
    stem <- gsub("[^A-Za-z0-9._-]+", "_", trimws(stem))
    stem <- gsub("^_+|_+$", "", stem)
    if (!nzchar(stem)) stem <- "af_profile"
    candidate <- stem
    suffix <- 2L
    while (candidate %in% existing_names) {
        candidate <- paste0(stem, "_", suffix)
        suffix <- suffix + 1L
    }
    candidate
}

#* Extract one AF profile through Spectreasy R
#* @post /workflow/af
function(req) {
    body <- gui_workflow_body(req)
    fcs_file <- gui_workflow_value(body, "fcs_file", file.path("scc", "unstained_cells.fcs"))
    bands <- gui_workflow_number(body, "af_n_bands", 100, integer = TRUE, minimum = 1)
    run <- gui_workflow_run(
        body,
        "af",
        spectreasy::extract_af_profile(
            fcs_file = fcs_file,
            af_n_bands = bands,
            af_max_cells = gui_workflow_number(body, "af_max_cells", 50000, integer = TRUE, minimum = 1),
            seed = gui_workflow_number(body, "seed", 1, integer = TRUE, minimum = 1),
            show_plot = FALSE,
            verbose = FALSE
        )
    )
    if (!isTRUE(run$success)) return(run)
    profile <- run$result
    save_name <- trimws(as.character(gui_workflow_value(body, "save_name", ""))[1])
    save_overwrite <- isTRUE(gui_workflow_bool(body, "save_overwrite", FALSE))
    if (!nzchar(save_name)) {
        existing_profiles <- tryCatch(spectreasy::list_af_profiles(), error = function(e) data.frame())
        existing_names <- if (nrow(existing_profiles) > 0L) existing_profiles$name else character()
        save_name <- gui_default_af_profile_name(fcs_file, existing_names)
    }
    saved_path <- tryCatch(
        spectreasy::save_af_profile(save_name, profile, overwrite = save_overwrite),
        error = function(e) e
    )
    if (inherits(saved_path, "error")) {
        run$success <- FALSE
        run$error <- conditionMessage(saved_path)
        return(run)
    }
    saved_path <- as.character(saved_path)
    run$result <- list(
        bands = nrow(profile$profile),
        detectors = ncol(profile$profile),
        source = fcs_file,
        profile_name = save_name,
        profile_path = saved_path
    )
    run
}
