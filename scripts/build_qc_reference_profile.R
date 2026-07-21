#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
    index <- match(flag, args)
    if (is.na(index) || index == length(args)) default else args[index + 1L]
}
split_roots <- function(value) trimws(strsplit(value, ";", fixed = TRUE)[[1]])

configured_roots <- Sys.getenv("SPECTREASY_ERLANGEN_ROOTS", unset = "")
default_roots <- unique(c(
    if (nzchar(configured_roots)) split_roots(configured_roots) else character(),
    Sys.glob(path.expand("~/Library/CloudStorage/GoogleDrive-*/My Drive/PhD/Project_Erlangen/Erlangen_Experiments")),
    Sys.glob(path.expand("~/Library/CloudStorage/GoogleDrive-*/My Drive/PhD/spectreasy_datasets"))
))
roots <- split_roots(value_after("--roots", paste(default_roots, collapse = ";")))
roots <- roots[dir.exists(roots)]
output_dir <- normalizePath(value_after("--output-dir", file.path("scratch", "ai_qc_validation", "reference_profiles")), mustWork = FALSE)
manifest_path <- value_after("--manifest", file.path(output_dir, "candidate_manifest.csv"))
profile_name <- value_after("--profile-name", "erlangen-panelflip-candidate")

if (!length(roots)) stop("No current dataset roots exist. Supply --roots root1;root2", call. = FALSE)
if (file.exists("DESCRIPTION") && requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(quiet = TRUE)
} else if (!requireNamespace("spectreasy", quietly = TRUE)) {
    stop("Install spectreasy or run this script from the package source tree.", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

matrix_files <- unique(unlist(lapply(roots, function(root) {
    list.files(root, pattern = "^scc_reference_matrix\\.csv$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
}), use.names = FALSE))
run_dir <- dirname(dirname(matrix_files))
run_id <- sprintf("run_%04d", seq_along(matrix_files))
is_panel_flip <- grepl("PBMC_test_culture_titrations_panelflip", matrix_files, fixed = TRUE)
seed_clean <- is_panel_flip & grepl("spectreasy_outputs_065[/\\\\]unmix_controls[/\\\\]scc_reference_matrix\\.csv$", matrix_files, ignore.case = TRUE)
if (!any(seed_clean) && any(is_panel_flip)) seed_clean[which(is_panel_flip)[1]] <- TRUE

manifest <- data.frame(
    run_id = run_id,
    source_root = vapply(matrix_files, function(path) {
        hit <- roots[startsWith(normalizePath(path, mustWork = FALSE), paste0(normalizePath(roots, mustWork = TRUE), .Platform$file.sep))]
        if (length(hit)) basename(hit[1]) else "unknown"
    }, character(1)),
    source_alias = sprintf("dataset_%04d", seq_along(matrix_files)),
    matrix_file = matrix_files,
    method = ifelse(grepl("065|spectreasy", matrix_files, ignore.case = TRUE), "Spectreasy", "unknown"),
    schema_status = ifelse(file.exists(file.path(dirname(matrix_files), "spectreasy_ai_qc_combined.json")), "ai_qc_1.0.0", "legacy_machine_artifacts"),
    review_label = ifelse(seed_clean, "clean", "unreviewed"),
    rationale = ifelse(seed_clean, "Explicit clean designation for the panel-flip project in the AI-ready QC feature specification.", ifelse(is_panel_flip, "Alternate output from the designated project; review separately to avoid duplicate-run inclusion.", "Human review required before cohort inclusion.")),
    reviewer = ifelse(seed_clean, "feature-specification", ""),
    review_date = ifelse(seed_clean, format(Sys.Date(), "%Y-%m-%d"), ""),
    include = seed_clean,
    exclusion_reason = ifelse(seed_clean, "", "unreviewed"),
    stringsAsFactors = FALSE
)
allowed_labels <- c("clean", "acceptable_with_notes", "exclude", "unreviewed")
if (file.exists(manifest_path)) {
    reviewed <- utils::read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(reviewed$review_label %in% allowed_labels)) stop("Manifest contains an unsupported review_label.", call. = FALSE)
    matching <- match(manifest$run_id, reviewed$run_id)
    for (field in intersect(c("review_label", "rationale", "reviewer", "review_date"), names(reviewed))) {
        use <- !is.na(matching)
        manifest[[field]][use] <- reviewed[[field]][matching[use]]
    }
    manifest$include <- manifest$review_label == "clean"
    manifest$exclusion_reason <- ifelse(manifest$include, "", manifest$review_label)
}
utils::write.csv(manifest, manifest_path, row.names = FALSE, quote = TRUE)

clean <- manifest[manifest$include & manifest$method != "unknown", , drop = FALSE]
metric_rows <- list()
contexts <- list()
for (i in seq_len(nrow(clean))) {
    matrix <- tryCatch(spectreasy:::.read_unmixing_matrix_csv(clean$matrix_file[i]), error = function(e) NULL)
    if (is.null(matrix)) next
    object <- spectreasy::collect_ai_qc(M = matrix, controls = list(M = matrix), scope = "control", privacy = "strict", reference = "none", context = list(method = clean$method[i]))
    metrics <- Filter(function(metric) isTRUE(metric$availability) && is.numeric(metric$value) && length(metric$value) == 1L && is.finite(metric$value), spectreasy:::.ai_qc_walk_metrics(object))
    metric_rows <- c(metric_rows, lapply(metrics, function(metric) data.frame(run_id = clean$run_id[i], metric_id = metric$metric_id, value = metric$value, stringsAsFactors = FALSE)))
    contexts[[length(contexts) + 1L]] <- object$experiment
}
metrics <- if (length(metric_rows)) do.call(rbind, metric_rows) else data.frame(run_id = character(), metric_id = character(), value = numeric())
profile_metrics <- if (nrow(metrics)) lapply(split(metrics$value, metrics$metric_id), function(values) list(
    n = length(values), median = stats::median(values), mad = stats::mad(values),
    interval = unname(stats::quantile(values, c(0.05, 0.95), names = FALSE)),
    aggregation = "median_mad_q05_q95"
)) else list()
profile <- list(
    name = profile_name, version = "0.1.0-candidate", schema_version = "1.0.0",
    strata = list(method = if (length(unique(clean$method)) == 1L) unique(clean$method) else NULL),
    cohort = list(n_discovered = nrow(manifest), n_reviewed = sum(manifest$review_label != "unreviewed"), n_clean = nrow(clean)),
    metrics = profile_metrics,
    source_versions = list(spectreasy = as.character(utils::packageVersion("spectreasy")), metric_definitions = "1.0.0"),
    de_identification = list(source_aliases_only = TRUE, raw_events = FALSE),
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)
profile$profile_hash <- spectreasy:::.ai_qc_sha256_text(jsonlite::toJSON(spectreasy:::.ai_qc_stable_order(profile), auto_unbox = TRUE, null = "null", na = "null", digits = NA))
profile_dir <- file.path(output_dir, ".spectreasy", "qc_reference_profiles")
dir.create(profile_dir, recursive = TRUE, showWarnings = FALSE)
profile_file <- file.path(profile_dir, paste0(profile_name, ".json"))
jsonlite::write_json(profile, profile_file, auto_unbox = TRUE, null = "null", na = "null", pretty = TRUE, digits = NA)

report_lines <- c(
    paste0("# QC reference-profile calibration: ", profile_name), "",
    paste0("- Discovered runs: ", nrow(manifest)),
    paste0("- Reviewed runs: ", sum(manifest$review_label != "unreviewed")),
    paste0("- Included clean runs: ", nrow(clean)),
    paste0("- Excluded or unreviewed runs: ", sum(!manifest$include)),
    paste0("- Profile grading enabled: ", if (nrow(clean) >= 5L) "yes" else "no; fewer than five matched clean runs"),
    paste0("- Profile hash: `", profile$profile_hash, "`"), "",
    "## Inclusion and exclusion", "",
    paste0("| Run | Source alias | Method | Review | Included | Reason |\n|---|---|---|---|---:|---|\n", paste(apply(manifest, 1, function(row) paste0("| ", row[["run_id"]], " | ", row[["source_alias"]], " | ", row[["method"]], " | ", row[["review_label"]], " | ", row[["include"]], " | ", row[["rationale"]], " |")), collapse = "\n")), "",
    "No raw events or source paths are stored in the profile. This candidate remains descriptive until at least five compatible reviewed clean runs are included."
)
writeLines(report_lines, file.path(output_dir, "calibration_report.md"), useBytes = TRUE)
html <- paste0("<!doctype html><meta charset='utf-8'><title>QC calibration</title><style>body{max-width:900px;margin:40px auto;font:15px/1.5 system-ui}table{border-collapse:collapse}td,th{padding:7px;border:1px solid #ccc}</style><pre>", paste(report_lines, collapse = "\n"), "</pre>")
writeLines(html, file.path(output_dir, "calibration_report.html"), useBytes = TRUE)
cat("Discovered:", nrow(manifest), " Included clean:", nrow(clean), "\nProfile:", profile_file, "\n")
