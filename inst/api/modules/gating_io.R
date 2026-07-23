get_gate_scc_dir <- function() {
    root <- get_gui_default_project_dir()
    normalizePath(getOption("spectreasy.gating_scc_dir", gui_project_input_path(root, "controls")), mustWork = FALSE)
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

gate_pick_csv_file_windows <- function(mode, initial_dir, default_name, title = NULL) {
    if (is.null(title) || !nzchar(title)) title <- if (mode == "save") "Save gate CSV" else "Load gate CSV"
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

gate_pick_csv_file_linux <- function(mode, initial_dir, default_name, title = NULL) {
    if (is.null(title) || !nzchar(title)) title <- if (mode == "save") "Save gate CSV" else "Load gate CSV"
    if (nzchar(Sys.which("zenity"))) {
        args <- c("--file-selection", "--title", title)
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
            c("--getsavefilename", file.path(initial_dir, default_name), "*.csv", "--title", title)
        } else {
            c("--getopenfilename", initial_dir, "*.csv", "--title", title)
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

gate_pick_csv_file <- function(
    mode = c("open", "save"),
    initial_dir = gate_working_dir(),
    default_name = basename(get_gate_file()),
    title = NULL
) {
    mode <- match.arg(tolower(mode[1]), c("open", "save"))
    if (!dir.exists(initial_dir)) {
        initial_dir <- gate_working_dir()
    }
    if (is.null(title) || !nzchar(title)) title <- if (mode == "save") "Save gate CSV" else "Load gate CSV"
    sysname <- Sys.info()[["sysname"]]

    if (identical(sysname, "Darwin")) {
        if (nzchar(Sys.which("osascript"))) {
            prompt <- title
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
            path <- gate_pick_csv_file_windows(mode, initial_dir, default_name, title)
            if (identical(path, "CANCEL")) return("CANCEL")
            if (nzchar(path)) {
                return(normalizePath(path, mustWork = mode == "open"))
            }
            return("")
        }
    }

    if (identical(sysname, "Linux")) {
        path <- gate_pick_csv_file_linux(mode, initial_dir, default_name, title)
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
    spectreasy:::.spectreasy_read_fcs(path, label = "GUI FCS file")
}

gate_pick_channel <- function(col_names, prefix, suffix) {
    scatter_type <- switch(toupper(as.character(prefix)[1]), FSC = "fsc", SSC = "ssc", NULL)
    if (!is.null(scatter_type)) {
        return(spectreasy:::.get_scatter_channel(col_names, scatter_type, suffix))
    }
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
