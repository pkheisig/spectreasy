args <- commandArgs(trailingOnly = TRUE)
output_directory <- if (length(args)) args[[1]] else tempfile("analysis-v2-full-benchmark-")
profile <- if (length(args) >= 2L) match.arg(args[[2]], c("small", "full")) else "full"
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

script_argument <- grep("^--file=", commandArgs(), value = TRUE)[1]
repository <- normalizePath(file.path(dirname(sub("^--file=", "", script_argument)), "..", ".."))
api <- new.env(parent = globalenv())
source(file.path(repository, "inst", "api", "gui_api.R"), local = api)

project <- tempfile("spectreasy-analysis-full-performance-")
on.exit(unlink(project, recursive = TRUE, force = TRUE), add = TRUE)
dir.create(file.path(project, "samples"), recursive = TRUE)

maximum_events <- if (identical(profile, "full")) 20000L else 2000L
set.seed(20260723)
truth <- seq(0, 1, length.out = maximum_events)
group <- cut(truth, breaks = c(-Inf, 0.24, 0.5, 0.76, Inf), labels = FALSE)
centers <- matrix(c(
    -4, -3, 0, 1, 2, 0, 1, -1, 0, 2, -2, 1,
    0, 5, 4, -2, 0, 1, -2, 2, 3, 0, 1, -1,
    5, -1, -3, 4, 3, -2, 1, 3, -1, -2, 2, 4,
    1, 2, 2, 0, -4, 5, 3, -2, 4, 1, 0, -3
), nrow = 4L, byrow = TRUE)
matrix <- centers[group, , drop = FALSE] +
    outer(truth, seq(0.2, 1.3, length.out = 12L)) +
    matrix(stats::rnorm(maximum_events * 12L, sd = 0.55), maximum_events, 12L)
colnames(matrix) <- paste0("M", seq_len(12L), "-A")
frame <- flowCore::flowFrame(matrix)
parameters <- flowCore::parameters(frame)@data
parameters$desc <- paste0("M", seq_len(12L))
flowCore::parameters(frame)@data <- parameters
flowCore::write.FCS(frame, file.path(project, "samples", "benchmark.fcs"))

method_registry <- api$gui_analysis_method_registry()
available <- Filter(function(method) isTRUE(method$available), method_registry)
reductions <- vapply(Filter(function(method) identical(method$family, "reduction"), available), `[[`, character(1), "id")
clustering <- vapply(Filter(function(method) identical(method$family, "clustering"), available), `[[`, character(1), "id")
trajectories <- vapply(Filter(function(method) identical(method$family, "trajectory"), available), `[[`, character(1), "id")

sizes <- if (identical(profile, "full")) c(500L, 2000L, 5000L) else c(500L, 2000L)
cases <- do.call(rbind, c(
    lapply(sizes, function(events) data.frame(events, family = "reduction", method = reductions)),
    lapply(sizes, function(events) data.frame(events, family = "clustering", method = clustering)),
    lapply(sizes, function(events) data.frame(events, family = "trajectory", method = trajectories)),
    if (identical(profile, "full")) list(data.frame(
        events = 20000L,
        family = c(rep("reduction", 3L), rep("clustering", 2L)),
        method = c("pca", "umap", "diffusion-map", "flowsom", "phenograph")
    )) else list()
))
cases <- unique(cases)
cases <- cases[cases$method %in% vapply(available, `[[`, character(1), "id"), , drop = FALSE]

memory_rss_mb <- function() {
    if (!requireNamespace("ps", quietly = TRUE)) return(NA_real_)
    unname(ps::ps_memory_info(ps::ps_handle())[["rss"]]) / 1024^2
}

run_once <- function(case, iteration) {
    method <- method_registry[[match(case$method, vapply(method_registry, `[[`, character(1), "id"))]]
    body <- list(
        file = "samples/benchmark.fcs",
        populationId = "root",
        method = case$method,
        markers = colnames(matrix),
        maxEvents = case$events,
        seed = 20260723L,
        rootEventId = 1L,
        cofactor = 150,
        neighbors = min(15L, max(5L, case$events %/% 40L)),
        clusters = 4L,
        perplexity = min(30, max(5, floor((case$events - 1) / 4)))
    )
    if (identical(case$family, "clustering")) {
        body$clusterMethod <- case$method
        body$reductionMethod <- "pca"
    } else if (identical(case$family, "reduction")) {
        body$clusterMethod <- "none"
        body$reductionMethod <- case$method
    } else {
        if ("clustering" %in% method$prerequisites) body$clusterMethod <- "flowsom"
        if ("reduction" %in% method$prerequisites) body$reductionMethod <- "pca"
    }
    gc()
    rss_before <- memory_rss_mb()
    started <- proc.time()
    result <- tryCatch(api$gui_analysis_run_method(project, body), error = identity)
    timing <- proc.time() - started
    rss_after <- memory_rss_mb()
    if (inherits(result, "error")) {
        return(data.frame(
            events = case$events, family = case$family, method = case$method,
            iteration, status = "failed", seconds = unname(timing[["elapsed"]]),
            user_seconds = unname(timing[["user.self"]]), system_seconds = unname(timing[["sys.self"]]),
            events_per_second = NA_real_, rss_delta_mb = rss_after - rss_before,
            returned_events = NA_integer_, finite_coordinates = FALSE,
            valid_pseudotime = FALSE, cache_reused = FALSE,
            error = conditionMessage(result), stringsAsFactors = FALSE
        ))
    }
    coordinates <- grep("^dimension_[123]$", names(result$events), value = TRUE)
    valid_pseudotime <- if ("pseudotime" %in% names(result$events)) {
        all(is.finite(result$events$pseudotime)) &&
            min(result$events$pseudotime) >= 0 &&
            max(result$events$pseudotime) <= 1
    } else {
        NA
    }
    artifact <- if (identical(case$family, "clustering")) {
        result$metadata$artifacts$clustering
    } else if (identical(case$method, "dpt")) {
        result$metadata$artifacts$embedding
    } else if (identical(case$family, "trajectory")) {
        result$metadata$artifacts$trajectory
    } else {
        result$metadata$artifacts$embedding
    }
    data.frame(
        events = case$events, family = case$family, method = case$method,
        iteration, status = "passed", seconds = unname(timing[["elapsed"]]),
        user_seconds = unname(timing[["user.self"]]), system_seconds = unname(timing[["sys.self"]]),
        events_per_second = case$events / unname(timing[["elapsed"]]),
        rss_delta_mb = rss_after - rss_before, returned_events = nrow(result$events),
        finite_coordinates = length(coordinates) >= 2L &&
            all(is.finite(as.matrix(result$events[, coordinates, drop = FALSE]))),
        valid_pseudotime, cache_reused = isTRUE(artifact$reused),
        error = "", stringsAsFactors = FALSE
    )
}

results <- list()
index <- 0L
for (row in seq_len(nrow(cases))) {
    case <- cases[row, , drop = FALSE]
    for (iteration in c("cold", "warm")) {
        index <- index + 1L
        cat(sprintf("[%d/%d] %s %s %s events\n", index, nrow(cases) * 2L, iteration, case$method, case$events))
        results[[index]] <- run_once(case, iteration)
        utils::write.csv(do.call(rbind, results), file.path(output_directory, "benchmark-results.csv"), row.names = FALSE)
    }
}
results <- do.call(rbind, results)
required <- results$status == "passed" &
    results$returned_events == results$events &
    results$finite_coordinates &
    (is.na(results$valid_pseudotime) | results$valid_pseudotime)
summary <- list(
    generated_at = as.character(Sys.time()),
    profile = profile,
    runtime = R.version.string,
    platform = R.version$platform,
    methods_available = vapply(available, `[[`, character(1), "id"),
    cases = nrow(cases),
    executions = nrow(results),
    all_passed = all(required),
    warm_cache_reused = all(results$cache_reused[results$iteration == "warm" & results$status == "passed"]),
    median_cold_seconds = stats::setNames(
        vapply(split(results$seconds[results$iteration == "cold"], results$method[results$iteration == "cold"]), stats::median, numeric(1)),
        unique(results$method[results$iteration == "cold"])
    )
)
jsonlite::write_json(summary, file.path(output_directory, "benchmark-summary.json"), auto_unbox = TRUE, pretty = TRUE)
print(aggregate(seconds ~ family + iteration, results[results$status == "passed", ], stats::median), row.names = FALSE)
cat(sprintf(
    "%d configurations, %d executions, all outputs valid: %s, all warm objects reused: %s\n",
    summary$cases, summary$executions, summary$all_passed, summary$warm_cache_reused
))
if (!isTRUE(summary$all_passed) || !isTRUE(summary$warm_cache_reused)) quit(status = 1L)
