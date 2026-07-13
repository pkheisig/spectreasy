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
    mode <- match.arg(mode)
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

gate_plot_to_base64 <- function(plot, width = 11.5, height = 3.2, dpi = 180) {
    tmp <- tempfile(fileext = ".png")
    on.exit(unlink(tmp), add = TRUE)
    ggplot2::ggsave(tmp, plot = plot, width = width, height = height, dpi = dpi, device = "png")
    raw_bytes <- readBin(tmp, "raw", file.info(tmp)$size)
    paste0("data:image/png;base64,", jsonlite::base64_enc(raw_bytes))
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
    verts <- gate_vertices_matrix(gate)
    if (is.null(verts) || nrow(verts) < 3) return(expr)
    x_channel <- as.character(gate_value(gate, "xChannel", ""))[1]
    y_channel <- as.character(gate_value(gate, "yChannel", ""))[1]
    if (!x_channel %in% colnames(expr) || !y_channel %in% colnames(expr)) return(expr)
    keep <- gate_point_in_polygon(expr[, x_channel], expr[, y_channel], verts)
    keep <- !is.na(keep) & keep
    expr[keep, , drop = FALSE]
}

gate_apply_positive_gate <- function(expr, gate, peak) {
    if (is.null(gate) || !peak %in% colnames(expr)) return(expr)
    verts <- gate_vertices_matrix(gate)
    if (is.null(verts) || nrow(verts) == 0) return(expr)
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
        if (!y_channel %in% colnames(expr)) return(expr)
        keep <- gate_point_in_polygon(expr[, x_channel], expr[, y_channel], verts)
    }
    if (is.null(keep)) return(expr)
    keep <- !is.na(keep) & keep
    expr[keep, , drop = FALSE]
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
    preserved_count <- 0L
    for (item in resolved) {
        positive_key <- paste0("positive:", item$filename)
        negative_key <- paste0("negative:", item$filename)
        keep_positive <- gate_is_finalized_histogram(gates[[positive_key]])
        keep_negative <- item$needs_negative && gate_is_finalized_histogram(gates[[negative_key]])
        preserved_count <- preserved_count + as.integer(keep_positive) + as.integer(keep_negative)
        if (keep_positive && (!item$needs_negative || keep_negative)) next

        path <- file.path(get_gate_scc_dir(), item$filename)
        ff <- gate_read_fcs(path)
        expr <- flowCore::exprs(ff)
        peak <- as.character(item$row$channel[1])
        if (!peak %in% colnames(expr)) {
            detectors <- tryCatch(
                spectreasy::get_sorted_detectors(flowCore::pData(flowCore::parameters(ff)))$names,
                error = function(e) colnames(expr)
            )
            detectors <- intersect(detectors, colnames(expr))
            if (length(detectors) == 0) {
                stop("Could not resolve a peak channel for ", item$filename, ".", call. = FALSE)
            }
            peak <- detectors[which.max(apply(expr[, detectors, drop = FALSE], 2, stats::var, na.rm = TRUE))]
        }
        gated <- gate_apply_polygon_gate(expr, item$cell_gate)
        gated <- gate_apply_polygon_gate(gated, item$singlet_gate)
        ranges <- gate_histogram_autogate_ranges(
            gated[, peak],
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
    list(
        gates = generated,
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

gate_build_spectrum_plot <- function(expr_filtered, spec_channels, det_info) {
    labels <- gate_detector_labels(det_info, spec_channels)
    y_power <- 1.5
    if (nrow(expr_filtered) == 0) {
        max_y <- 6
        y_breaks_orig <- 0:max_y
        y_breaks_trans <- y_breaks_orig^y_power
        y_labels <- vapply(y_breaks_orig, function(x) paste0("10^", x), character(1))
        p_empty <- ggplot2::ggplot(data.frame(ch_idx = seq_along(spec_channels), y = 0), ggplot2::aes(ch_idx, y)) +
            ggplot2::scale_x_continuous(breaks = seq_along(spec_channels), labels = labels) +
            ggplot2::scale_y_continuous(limits = c(0, (max_y + 0.5)^y_power), breaks = y_breaks_trans, labels = y_labels) +
            ggplot2::coord_cartesian(expand = FALSE) +
            ggplot2::labs(title = NULL, x = NULL, y = "Intensity") +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
                legend.position = "none",
                panel.background = ggplot2::element_rect(fill = "white", color = NA),
                panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank()
            )
        return(gate_plot_to_base64(p_empty))
    }

    log_mat <- log10(pmax(expr_filtered[, spec_channels, drop = FALSE], 1e-3))
    min_y <- floor(min(log_mat, na.rm = TRUE))
    max_y <- ceiling(max(log_mat, na.rm = TRUE))
    breaks <- seq(min_y, max_y, length.out = 151)
    bin_mid <- (breaks[-1] + breaks[-length(breaks)]) / 2
    bin_height <- breaks[2] - breaks[1]
    counts_mat <- vapply(seq_len(ncol(log_mat)), function(j) {
        as.numeric(graphics::hist(log_mat[, j], breaks = breaks, plot = FALSE)$counts)
    }, numeric(length(bin_mid)))
    rownames(counts_mat) <- as.character(seq_along(bin_mid))
    colnames(counts_mat) <- as.character(seq_len(ncol(log_mat)))
    dt_c <- as.data.frame(as.table(counts_mat), stringsAsFactors = FALSE)
    names(dt_c) <- c("bin_idx", "ch_idx", "count")
    dt_c$bin_idx <- as.integer(dt_c$bin_idx)
    dt_c$ch_idx <- as.integer(dt_c$ch_idx)
    dt_c$y_orig <- bin_mid[dt_c$bin_idx]
    dt_c$fill <- log10(dt_c$count + 1)
    min_bin_count <- if (nrow(expr_filtered) <= 3000) 1 else 3
    dt_c <- dt_c[dt_c$count >= min_bin_count, , drop = FALSE]
    if (nrow(dt_c) == 0) {
        empty_expr <- expr_filtered[0, , drop = FALSE]
        return(gate_build_spectrum_plot(empty_expr, spec_channels, det_info))
    }
    dt_c$y <- dt_c$y_orig^y_power
    fill_lo <- min(dt_c$fill, na.rm = TRUE)
    fill_hi <- stats::quantile(dt_c$fill, 0.96, na.rm = TRUE)
    y_breaks_orig <- 0:ceiling(max_y)
    y_breaks_trans <- y_breaks_orig^y_power
    y_labels <- vapply(y_breaks_orig, function(x) paste0("10^", x), character(1))

    p_spec <- ggplot2::ggplot(dt_c, ggplot2::aes(ch_idx, y, fill = fill)) +
        ggplot2::geom_tile(width = 0.7, height = bin_height * 3) +
        ggplot2::scale_fill_gradientn(colors = c("#0000FF", "#00FFFF", "#00FF00", "#FFFF00", "#FF0000"), limits = c(fill_lo, fill_hi), oob = scales::squish) +
        ggplot2::scale_x_continuous(breaks = seq_along(spec_channels), labels = labels) +
        ggplot2::scale_y_continuous(limits = c(0, (max_y + 0.5)^y_power), breaks = y_breaks_trans, labels = y_labels) +
        ggplot2::coord_cartesian(expand = FALSE) +
        ggplot2::labs(title = NULL, x = NULL, y = "Intensity") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
            legend.position = "none",
            panel.background = ggplot2::element_rect(fill = "white", color = NA),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank()
        )
    gate_plot_to_base64(p_spec)
}

gate_spectrum_for_file <- function(filename) {
    df <- gate_read_mapping()
    name <- gate_safe_basename(filename)
    row <- df[df$filename == name, , drop = FALSE]
    if (nrow(row) == 0) stop("File is not present in mapping: ", name, call. = FALSE)
    path <- file.path(get_gate_scc_dir(), name)
    ff <- gate_read_fcs(path)
    expr <- flowCore::exprs(ff)
    pd <- flowCore::pData(flowCore::parameters(ff))
    det_info <- tryCatch(spectreasy::get_sorted_detectors(pd), error = function(e) NULL)
    if (is.null(det_info) || length(det_info$names) == 0) return(NULL)
    spec_channels <- intersect(det_info$names, colnames(expr))
    if (length(spec_channels) == 0) return(NULL)

    peak <- as.character(row$channel[1])
    if (!peak %in% colnames(expr)) peak <- spec_channels[which.max(apply(expr[, spec_channels, drop = FALSE], 2, stats::var, na.rm = TRUE))]

    cache_data <- getOption("spectreasy.gating_state_cache")
    cache <- if (!is.null(cache_data)) cache_data$gates else NULL
    control_type <- gate_normalize_control_type(row$control.type[1])

    expr_filtered <- expr
    expr_filtered <- gate_apply_polygon_gate(expr_filtered, gate_cached_gate(cache, "cell", name, control_type))
    expr_filtered <- gate_apply_polygon_gate(expr_filtered, gate_cached_gate(cache, "singlet", name, control_type))
    if (!isTRUE(row$is_af[1])) {
        expr_filtered <- gate_apply_positive_gate(expr_filtered, gate_cached_gate(cache, "positive", name, control_type), peak)
    }
    gate_build_spectrum_plot(expr_filtered, spec_channels, det_info)
}

gate_payload_for_file <- function(filename, max_points = 3000L) {
    name <- gate_safe_basename(filename)
    cache_key <- paste(name, max_points, sep = "::")
    cache <- getOption("spectreasy.gating_payload_cache", list())
    if (!is.null(cache[[cache_key]])) return(cache[[cache_key]])

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
        det <- tryCatch(spectreasy::get_sorted_detectors(flowCore::pData(flowCore::parameters(ff)))$names, error = function(e) colnames(expr))
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
    cache[[cache_key]] <- payload
    options(spectreasy.gating_payload_cache = cache)
    payload
}

#* @filter logger
function(req) {
    if (isTRUE(getOption("spectreasy.gui_request_log", FALSE))) {
        cat(as.character(Sys.time()), "-", req$REQUEST_METHOD, req$PATH_INFO, "\n")
    }
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
        unmixing_method = get_unmixing_method(),
        gui_mode = getOption("spectreasy.gui_mode", "tuner"),
        panel_cytometer = getOption("spectreasy.panel_cytometer", "aurora")
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
            control.type = row_value(row, c("controlType", "control.type"), "cells"),
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

#* Auto-generate positive and negative histogram gates for all non-AF controls
#* @post /gate_histogram_autogate
function(req) {
    tryCatch({
        body <- jsonlite::fromJSON(req$postBody, simplifyVector = FALSE)
        result <- gate_autogenerate_histograms(body$gates)
        list(
            success = TRUE,
            gates = result$gates,
            files_processed = result$files_processed,
            gates_generated = result$gates_generated,
            gates_preserved = result$gates_preserved
        )
    }, error = function(e) {
        list(success = FALSE, error = conditionMessage(e))
    })
}

#* Render selected SCC spectrum for the current gate cache
#* @get /gate_spectrum
#* @param filename
function(filename) {
    tryCatch(
        list(spectrum = gate_spectrum_for_file(filename = filename)),
        error = function(e) list(error = conditionMessage(e), spectrum = NULL)
    )
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
        return(list(gates = list(), pointSize = 1.5, maxPoints = 50000, histogramBins = 100, histogramTransform = "asinh", eventCountVersion = 2))
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
    list(profiles = profiles)
}

#* Delete a saved AF profile
#* @delete /af_profiles
#* @param name Profile name
function(name = "") {
    if (is.null(name) || !nzchar(trimws(as.character(name)[1]))) {
        return(list(success = FALSE, error = "A profile name is required."))
    }
    tryCatch({
        spectreasy::delete_af_profile(as.character(name)[1])
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
        error = function(e) list()
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

gui_workflow_number <- function(body, key, fallback = 0, integer = FALSE, minimum = NULL) {
    value <- suppressWarnings(as.numeric(gui_workflow_value(body, key, fallback)))
    if (length(value) == 0 || !is.finite(value[1])) value <- fallback
    value <- as.numeric(value[1])
    if (!is.null(minimum)) value <- max(value, minimum)
    if (isTRUE(integer)) value <- as.integer(round(value))
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
        return(list(spectreasy_weight_quantile = gui_workflow_number(body, "spectreasy_weight_quantile", 0.9, minimum = 0)))
    }
    list()
}

gui_workflow_root <- function(body) {
    root <- gui_workflow_value(body, "projectPath", get_matrix_dir())
    if (!dir.exists(root)) return(get_matrix_dir())
    normalizePath(root, mustWork = TRUE)
}

gui_workflow_run <- function(body, action, expr) {
    root <- gui_workflow_root(body)
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(root)
    tryCatch(
        {
            result <- force(expr)
            list(
                success = TRUE,
                action = action,
                project_path = root,
                finished_at = as.character(Sys.time()),
                result = result
            )
        },
        error = function(e) {
            list(
                success = FALSE,
                action = action,
                project_path = root,
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

#* Project scan and workflow prerequisites
#* @get /project/status
#* @param project_path Optional project directory
function(project_path = "") {
    root <- if (is.null(project_path) || !nzchar(trimws(as.character(project_path)[1]))) get_matrix_dir() else as.character(project_path)[1]
    root <- normalizePath(root, mustWork = FALSE)
    gui_project_scan(root)
}

gui_project_relative_path <- function(path, root) {
    root_pattern <- gsub("([.|(){}+*?^$\\[\\]\\\\])", "\\\\\\1", normalizePath(root, mustWork = FALSE))
    gsub("\\\\", "/", sub(paste0("^", root_pattern, "[/\\\\]?"), "", path))
}

gui_project_report_files <- function(root, files = NULL) {
    root <- normalizePath(root, mustWork = FALSE)
    if (is.null(files)) {
        files <- if (dir.exists(root)) {
            list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
        } else {
            character()
        }
    }
    relative_files <- gui_project_relative_path(files, root)
    files[
        grepl("\\.(html?|pdf)$", relative_files, ignore.case = TRUE) &
            grepl("(^|/)(reports|spectreasy_outputs)/", relative_files, ignore.case = TRUE, perl = TRUE)
    ]
}

#* Discover report artifacts and compare them with upstream project files
#* @get /project/reports
function(project_path = "") {
    root <- if (is.null(project_path) || !nzchar(trimws(as.character(project_path)[1]))) get_matrix_dir() else as.character(project_path)[1]
    root <- normalizePath(root, mustWork = FALSE)
    files <- if (dir.exists(root)) list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE) else character()
    reports <- gui_project_report_files(root, files = files)
    if (length(reports) == 0L) return(list(reports = data.frame()))
    relative <- function(path) gui_project_relative_path(path, root)
    rows <- lapply(reports, function(report) {
        is_sample <- grepl("sample|qc_samples", report, ignore.case = TRUE)
        is_panel <- grepl("panel", report, ignore.case = TRUE)
        source_pattern <- if (is_sample) {
            "(^|/)(samples?/.*\\.fcs$|.*reference_matrix.*\\.csv$|.*detector_noise.*\\.csv$|.*spectral_variant.*\\.rds$)"
        } else {
            "(^|/)(scc/.*\\.fcs$|fcs_mapping\\.csv$|.*gate.*\\.csv$|.*reference_matrix.*\\.csv$|.*detector_noise.*\\.csv$)"
        }
        source_files <- files[grepl(source_pattern, relative(files), ignore.case = TRUE, perl = TRUE)]
        source_files <- setdiff(source_files, report)
        stale <- length(source_files) > 0L && any(file.info(source_files)$mtime > file.info(report)$mtime)
        data.frame(
            path = relative(report),
            report_type = if (is_panel) "Panel overview" else if (is_sample) "Sample QC" else "Control QC",
            format = if (grepl("\\.pdf$", report, ignore.case = TRUE)) "PDF" else "HTML",
            created = format(file.info(report)$mtime, "%Y-%m-%d %H:%M:%S %Z"),
            status = if (stale) "stale" else "current",
            stringsAsFactors = FALSE
        )
    })
    list(reports = do.call(rbind, rows))
}

#* Stream a project artifact through the local R backend
#* @get /project/file
#* @param path Relative path inside the active project
function(path = "", res) {
    root <- normalizePath(get_matrix_dir(), mustWork = FALSE)
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

#* Change the active project folder from the cockpit
#* @post /project/context
function(req) {
    body <- gui_workflow_body(req)
    project_path <- gui_workflow_path(body, "projectPath", "", allow_empty = FALSE)
    if (!dir.exists(project_path)) {
        return(list(success = FALSE, error = paste("Project folder not found:", project_path)))
    }
    project_path <- normalizePath(project_path, mustWork = TRUE)
    options(
        spectreasy.project_dir = project_path,
        spectreasy.matrix_dir = project_path,
        spectreasy.samples_dir = file.path(project_path, "samples"),
        spectreasy.gating_scc_dir = file.path(project_path, "scc"),
        spectreasy.gating_control_file = file.path(project_path, "fcs_mapping.csv"),
        spectreasy.gating_gate_file = file.path(project_path, "ssc_gate_config.csv")
    )
    list(success = TRUE, project = gui_project_scan(project_path))
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
    scc_dir <- gui_workflow_value(body, "scc_dir", "scc")
    control_file <- gui_workflow_value(body, "control_file", "fcs_mapping.csv")
    output_dir <- gui_workflow_value(body, "output_dir", file.path("spectreasy_outputs", "unmix_controls"))
    method <- gui_workflow_value(body, "method", get_unmixing_method())
    cytometer <- gui_workflow_value(body, "cytometer", "auto")
    gate_file <- gui_workflow_file_or_null(
        gui_workflow_value(body, "gate_file", "ssc_gate_config.csv"),
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
            af_min_cluster_events = gui_workflow_number(body, "af_min_cluster_events", 20, integer = TRUE, minimum = 1),
            af_min_cluster_proportion = gui_workflow_number(body, "af_min_cluster_proportion", 0.005, minimum = 0),
            default_sample_type = gui_workflow_value(body, "default_sample_type", "beads"),
            histogram_pct_beads = gui_workflow_number(body, "histogram_pct_beads", 0.98, minimum = 0),
            histogram_direction_beads = gui_workflow_value(body, "histogram_direction_beads", "right"),
            histogram_pct_cells = gui_workflow_number(body, "histogram_pct_cells", 0.35, minimum = 0),
            histogram_direction_cells = gui_workflow_value(body, "histogram_direction_cells", "right"),
            outlier_percentile = gui_workflow_number(body, "outlier_percentile", 0.02, minimum = 0),
            debris_percentile = gui_workflow_number(body, "debris_percentile", 0.08, minimum = 0),
            bead_gate_scale = gui_workflow_number(body, "bead_gate_scale", 1.3, minimum = 0),
            max_clusters = gui_workflow_number(body, "max_clusters", 10, integer = TRUE, minimum = 1),
            min_cluster_proportion = gui_workflow_number(body, "min_cluster_proportion", 0.03, minimum = 0),
            gate_contour_beads = gui_workflow_number(body, "gate_contour_beads", 0.95, minimum = 0),
            gate_contour_cells = gui_workflow_number(body, "gate_contour_cells", 0.9, minimum = 0),
            subsample_n = gui_workflow_number(body, "subsample_n", 5000, integer = TRUE, minimum = 1),
            rwls_max_iter = gui_workflow_number(body, "rwls_max_iter", 1, integer = TRUE, minimum = 1),
            n_threads = gui_workflow_number(body, "unmix_threads", 1, integer = TRUE, minimum = 1),
            save_qc_plots = gui_workflow_bool(body, "save_qc_plots", TRUE),
            save_report = gui_workflow_bool(body, "save_report", TRUE),
            output_format = tolower(gui_workflow_value(body, "output_format", "pdf")),
            use_scatter_gating = gui_workflow_bool(body, "use_scatter_gating", TRUE),
            manual_gating = FALSE,
            manual_gate_file = gate_file,
            gating_file = gate_file,
            clean_scc_with_unstained = gui_workflow_bool(body, "clean_scc_with_unstained", TRUE),
            scc_background_method = gui_workflow_value(body, "scc_background_method", "scatter_knn"),
            scc_background_k = gui_workflow_number(body, "scc_background_k", 2, integer = TRUE, minimum = 1),
            spectral_variant_som_nodes = gui_workflow_number(body, "spectral_variant_som_nodes", 16, integer = TRUE, minimum = 1),
            spectral_variant_top_k = gui_workflow_number(body, "spectral_variant_top_k", 3, integer = TRUE, minimum = 1),
            spectral_variant_cosine_threshold = gui_workflow_number(body, "spectral_variant_cosine_threshold", 0.98, minimum = 0),
            spectral_variant_max_variants = gui_workflow_number(body, "spectral_variant_max_variants", 8, integer = TRUE, minimum = 1),
            spectral_variant_min_events = gui_workflow_number(body, "spectral_variant_min_events", 50, integer = TRUE, minimum = 1),
            autospectral_n_candidates = gui_workflow_number(body, "autospectral_n_candidates", 1000, integer = TRUE, minimum = 1),
            autospectral_n_spectral = gui_workflow_number(body, "autospectral_n_spectral", 200, integer = TRUE, minimum = 1),
            autospectral_min_events = gui_workflow_number(body, "autospectral_min_events", 10, integer = TRUE, minimum = 1),
            refine = gui_workflow_bool(body, "refine", FALSE)
        ),
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
    output_dir <- gui_workflow_value(body, "output_dir", file.path("spectreasy_outputs", "unmix_samples", "unmixed_fcs"))
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
            spectral_variant_positive_fraction = gui_workflow_number(body, "spectral_variant_positive_fraction", 0.02, minimum = 0),
            spectral_variant_min_improvement = gui_workflow_number(body, "spectral_variant_min_improvement", 0.01, minimum = 0),
            spectral_variant_library_file = gui_workflow_file_or_null(gui_workflow_value(body, "spectral_variant_library_file", ""), root = gui_workflow_root(body)),
            estimate_af = gui_workflow_bool(body, "estimate_af", FALSE),
            output_dir = output_dir,
            write_fcs = gui_workflow_bool(body, "write_fcs", TRUE),
            save_report = gui_workflow_bool(body, "save_report", TRUE),
            output_format = tolower(gui_workflow_value(body, "output_format", "pdf")),
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

#* Generate a synthetic SCC FCS from one matrix signature
#* @post /workflow/synthetic
function(req) {
    body <- gui_workflow_body(req)
    root <- gui_workflow_root(body)
    matrix_input <- gui_workflow_value(body, "matrix_file", "")
    matrix_file <- gui_workflow_file_or_null(matrix_input, root = root)
    if (is.null(matrix_file) || !file.exists(matrix_file)) {
        return(list(success = FALSE, error = "Select a readable reference matrix before generating synthetic SCC."))
    }
    matrix_df <- read_matrix_csv(matrix_file)
    if (ncol(matrix_df) < 2 || nrow(matrix_df) == 0) return(list(success = FALSE, error = "The selected matrix has no usable signatures."))
    marker <- gui_workflow_value(body, "marker", as.character(matrix_df[[1]][1]))
    marker_index <- match(marker, as.character(matrix_df[[1]]))
    if (is.na(marker_index)) return(list(success = FALSE, error = paste("Marker not found in matrix:", marker)))
    detectors <- colnames(matrix_df)[-1]
    signature <- suppressWarnings(as.numeric(matrix_df[marker_index, -1, drop = TRUE]))
    events_n <- gui_workflow_number(body, "events", 1000, integer = TRUE, minimum = 1)
    noise <- gui_workflow_number(body, "noise", 0.02, minimum = 0)
    seed <- gui_workflow_number(body, "seed", 1, integer = TRUE, minimum = 1)
    set.seed(seed)
    synthetic <- matrix(rep(signature, each = events_n), nrow = events_n, ncol = length(signature), byrow = FALSE)
    synthetic <- matrix(pmax(0, synthetic * exp(matrix(stats::rnorm(events_n * length(signature), mean = 0, sd = noise), nrow = events_n))), nrow = events_n, ncol = length(signature))
    colnames(synthetic) <- detectors
    output_input <- gui_workflow_value(body, "output_file", file.path("spectreasy_outputs", "synthetic_scc", paste0("synthetic_", gsub("[^A-Za-z0-9._-]+", "_", marker), ".fcs")))
    output_file <- if (grepl("^(/|[A-Za-z]:[/\\\\])", output_input)) output_input else file.path(root, output_input)
    output_file <- normalizePath(output_file, mustWork = FALSE)
    root_prefix <- paste0(normalizePath(root, mustWork = FALSE), .Platform$file.sep)
    if (!startsWith(output_file, root_prefix)) return(list(success = FALSE, error = "Synthetic output must stay inside the active project."))
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    flowCore::write.FCS(flowCore::flowFrame(synthetic), output_file)
    truth_file <- sub("\\.fcs$", "_truth.csv", output_file, ignore.case = TRUE)
    utils::write.csv(data.frame(marker = marker, detector = detectors, signature = signature), truth_file, row.names = FALSE, quote = TRUE)
    list(success = TRUE, result = list(output_file = output_file, truth_file = truth_file, marker = marker, events = events_n, detectors = length(detectors)))
}

#* Render a report from existing Spectreasy outputs
#* @post /workflow/report
function(req) {
    body <- gui_workflow_body(req)
    report_type <- tolower(gui_workflow_value(body, "report_type", "control"))
    output_format <- tolower(gui_workflow_value(body, "output_format", "html"))
    if (!output_format %in% c("html", "pdf")) {
        return(list(success = FALSE, error = "Report format must be 'html' or 'pdf'."))
    }
    overwrite <- tolower(gui_workflow_value(body, "overwrite", "version"))
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
            report_file <- file.path(report_dir, paste0("qc_controls_report.", output_format))
            spectreasy::qc_controls(
                M = M,
                unmixing_matrix_file = matrix_file,
                scc_dir = file.path(root, "scc"),
                control_file = file.path(root, "fcs_mapping.csv"),
                output_file = report_file,
                output_format = output_format,
                overwrite = overwrite
            )
        }
    )
    if (!isTRUE(run$success)) return(run)
    report_path <- if (is.list(run$result) && !is.null(run$result$output_file)) run$result$output_file else report_file
    run$result <- list(report_file = report_path, output_format = output_format)
    run
}

#* CORS preflight for workflow_af
#* @options /workflow/af
function(res) {
    return("")
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
            af_min_cluster_events = gui_workflow_number(body, "af_min_cluster_events", 20, integer = TRUE, minimum = 1),
            af_min_cluster_proportion = gui_workflow_number(body, "af_min_cluster_proportion", 0.005, minimum = 0),
            seed = gui_workflow_number(body, "seed", 1, integer = TRUE, minimum = 1),
            show_plot = FALSE,
            verbose = FALSE
        )
    )
    if (!isTRUE(run$success)) return(run)
    profile <- run$result
    save_name <- trimws(as.character(gui_workflow_value(body, "save_name", ""))[1])
    save_overwrite <- isTRUE(gui_workflow_bool(body, "save_overwrite", FALSE))
    saved_path <- NULL
    if (nzchar(save_name)) {
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
    }
    run$result <- list(
        bands = nrow(profile$profile),
        detectors = ncol(profile$profile),
        source = fcs_file,
        profile_name = if (nzchar(save_name)) save_name else NULL,
        profile_path = saved_path
    )
    run
}
