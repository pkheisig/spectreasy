args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: Rscript analysis_v2_real_pkh.R DATASET_ROOT [OUTPUT_DIRECTORY]", call. = FALSE)
dataset_root <- normalizePath(args[[1]], mustWork = TRUE)
output_directory <- if (length(args) >= 2L) args[[2]] else tempfile("analysis-v2-real-pkh-")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

script_argument <- grep("^--file=", commandArgs(), value = TRUE)[1]
repository <- normalizePath(file.path(dirname(sub("^--file=", "", script_argument)), "..", ".."))
api <- new.env(parent = globalenv())
source(file.path(repository, "inst", "api", "gui_api.R"), local = api)

preferred <- c(
    file.path(dataset_root, "samples_full", "1.fcs"),
    file.path(dataset_root, "samples", "1.fcs")
)
sample_input <- preferred[file.exists(preferred)][1]
if (is.na(sample_input)) {
    candidates <- list.files(dataset_root, pattern = "[.]fcs$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    sample_input <- candidates[!grepl("unstained", basename(candidates), ignore.case = TRUE)][1]
}
if (is.na(sample_input) || !file.exists(sample_input)) stop("No PKH sample FCS file was found.", call. = FALSE)
unstained_candidates <- list.files(dataset_root, pattern = "unstained.*[.]fcs$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)

project <- tempfile("spectreasy-real-pkh-project-")
on.exit(unlink(project, recursive = TRUE, force = TRUE), add = TRUE)
dir.create(file.path(project, "samples"), recursive = TRUE)
dir.create(file.path(project, "scc"), recursive = TRUE)
file.copy(sample_input, file.path(project, "samples", "PKH_1.fcs"), overwrite = TRUE)
if (length(unstained_candidates)) {
    file.copy(unstained_candidates[[1]], file.path(project, "scc", "PKH_Unstained.fcs"), overwrite = TRUE)
}
relative_file <- "samples/PKH_1.fcs"
frame <- api$gui_analysis_flow_frame(project, relative_file)
values <- flowCore::exprs(frame)
parameters <- flowCore::parameters(frame)@data
channels <- colnames(values)
descriptions <- as.character(parameters$desc)
descriptions[is.na(descriptions)] <- ""
fsc <- channels[match(TRUE, grepl("^FSC.*-A$", channels, ignore.case = TRUE))]
ssc <- channels[match(TRUE, grepl("SSC.*-A$", channels, ignore.case = TRUE))]
if (is.na(fsc)) fsc <- channels[[1]]
if (is.na(ssc)) ssc <- setdiff(channels, fsc)[[1]]
marker_indices <- which(
    nzchar(descriptions) &
        !grepl("^(FSC|SSC|Time)$", descriptions, ignore.case = TRUE) &
        !grepl("^(FSC|SSC|Time)", channels, ignore.case = TRUE)
)
markers <- channels[marker_indices]
if (length(markers) < 3L) markers <- setdiff(channels, c(fsc, ssc, grep("^Time$", channels, value = TRUE, ignore.case = TRUE)))
markers <- head(markers, 8L)
if (length(markers) < 3L) stop("Real PKH file has too few usable marker channels.", call. = FALSE)

qx <- stats::quantile(values[, fsc], c(0.02, 0.25, 0.5, 0.75, 0.98), na.rm = TRUE)
qy <- stats::quantile(values[, ssc], c(0.02, 0.25, 0.5, 0.75, 0.98), na.rm = TRUE)
qm <- stats::quantile(values[, markers[[1]]], c(0.02, 0.35, 0.65, 0.98), na.rm = TRUE)
workspace <- api$gui_analysis_default_workspace()
workspace$source_path <- "samples"
workspace$selected_file <- relative_file
workspace$active_population_id <- "cells"
workspace$root_event_id <- 1L
workspace$root_population_id <- "root"
workspace$root_source_file <- relative_file
workspace$populations <- c(workspace$populations, list(
    list(
        id = "cells", name = "Cells", parent_id = "root", type = "rectangle",
        role = NULL, source_file = NULL, x = fsc, y = ssc,
        geometry = list(x_min = unname(qx[[1]]), x_max = unname(qx[[5]]), y_min = unname(qy[[1]]), y_max = unname(qy[[5]]))
    ),
    list(
        id = "ellipse", name = "Central ellipse", parent_id = "cells", type = "ellipse",
        role = NULL, source_file = NULL, x = fsc, y = ssc,
        geometry = list(
            center_x = unname(qx[[3]]), center_y = unname(qy[[3]]),
            radius_x = unname((qx[[4]] - qx[[2]]) * 1.2), radius_y = unname((qy[[4]] - qy[[2]]) * 1.2)
        )
    ),
    list(
        id = "polygon", name = "Central polygon", parent_id = "cells", type = "polygon",
        role = NULL, source_file = NULL, x = fsc, y = ssc,
        geometry = list(points = list(
            list(x = unname(qx[[2]]), y = unname(qy[[2]])),
            list(x = unname(qx[[4]]), y = unname(qy[[2]])),
            list(x = unname(qx[[4]]), y = unname(qy[[4]])),
            list(x = unname(qx[[2]]), y = unname(qy[[4]]))
        ))
    ),
    list(
        id = "negative", name = "Marker negative", parent_id = "root", type = "range",
        role = "negative", source_file = NULL, x = markers[[1]], y = NULL,
        geometry = list(min = unname(qm[[1]]), max = unname(qm[[2]]))
    ),
    list(
        id = "positive", name = "Marker positive", parent_id = "root", type = "range",
        role = "positive", source_file = NULL, x = markers[[1]], y = NULL,
        geometry = list(min = unname(qm[[3]]), max = unname(qm[[4]]))
    )
))
workspace$plots <- list(
    list(id = "plot-1", type = "scatter", population_id = "root", x = fsc, y = ssc, color_by = "density", x_transform = "linear", y_transform = "linear"),
    list(id = "plot-2", type = "histogram", population_id = "root", x = markers[[1]], y = ssc, color_by = "density", x_transform = "asinh", y_transform = "linear")
)
saved <- api$gui_analysis_write_workspace(project, workspace)
restored <- api$gui_analysis_read_workspace(project)
stopifnot(length(restored$populations) == 6L, length(restored$plots) == 2L)

population_ids <- vapply(restored$populations, `[[`, character(1), "id")
statistics <- lapply(population_ids, function(population_id) {
    api$gui_analysis_population_statistics(project, relative_file, population_id, markers)
})
counts <- vapply(statistics, `[[`, numeric(1), "count")
names(counts) <- population_ids
stopifnot(counts[["root"]] == nrow(values), counts[["cells"]] > 0, counts[["cells"]] < counts[["root"]])
staining <- api$gui_analysis_staining_index(project, relative_file, markers[[1]])
stopifnot(
    is.finite(staining$staining_index),
    staining$positive_count > 0,
    staining$negative_count > 0
)
export <- api$gui_analysis_export(project, list(
    files = relative_file,
    populationId = "cells",
    format = "both",
    maxEvents = 2000L,
    seed = 20260723L,
    outputFolder = "spectreasy_outputs/analysis/acceptance"
))
stopifnot(all(file.exists(file.path(project, vapply(export$files, `[[`, character(1), "path")))))

registry <- api$gui_analysis_method_registry()
available <- Filter(function(method) isTRUE(method$available), registry)
method_results <- lapply(available, function(method) {
    body <- list(
        file = relative_file,
        populationId = "root",
        method = method$id,
        markers = markers,
        maxEvents = 600L,
        seed = 20260723L,
        rootEventId = 1L,
        rootSourceFile = relative_file,
        cofactor = 150,
        neighbors = 12L,
        clusters = 6L,
        perplexity = 20
    )
    if (identical(method$family, "clustering")) {
        body$clusterMethod <- method$id
        body$reductionMethod <- "pca"
    } else if (identical(method$family, "reduction")) {
        body$clusterMethod <- "none"
        body$reductionMethod <- method$id
    } else {
        if ("clustering" %in% method$prerequisites) body$clusterMethod <- "flowsom"
        if ("reduction" %in% method$prerequisites) body$reductionMethod <- "pca"
    }
    started <- proc.time()[["elapsed"]]
    result <- tryCatch(api$gui_analysis_run_method(project, body), error = identity)
    seconds <- proc.time()[["elapsed"]] - started
    if (inherits(result, "error")) {
        return(data.frame(method = method$id, family = method$family, status = "failed", seconds, rows = NA_integer_, error = conditionMessage(result)))
    }
    dimensions <- grep("^dimension_", names(result$events), value = TRUE)
    valid <- nrow(result$events) == 600L &&
        length(dimensions) >= 2L &&
        all(is.finite(as.matrix(result$events[, dimensions, drop = FALSE]))) &&
        (!"pseudotime" %in% names(result$events) || all(is.finite(result$events$pseudotime)))
    data.frame(method = method$id, family = method$family, status = if (valid) "passed" else "invalid", seconds, rows = nrow(result$events), error = "")
})
method_results <- do.call(rbind, method_results)
utils::write.csv(method_results, file.path(output_directory, "real-pkh-methods.csv"), row.names = FALSE)
report <- list(
    generated_at = as.character(Sys.time()),
    dataset = basename(dataset_root),
    source_file = basename(sample_input),
    source_events = nrow(values),
    parameters = ncol(values),
    selected_markers = unname(markers),
    population_counts = as.list(counts),
    staining_index = staining,
    export_files = export$files,
    methods = split(method_results, seq_len(nrow(method_results))),
    all_methods_passed = all(method_results$status == "passed")
)
jsonlite::write_json(report, file.path(output_directory, "real-pkh-acceptance.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
print(method_results, row.names = FALSE)
if (!isTRUE(report$all_methods_passed)) quit(status = 1L)
