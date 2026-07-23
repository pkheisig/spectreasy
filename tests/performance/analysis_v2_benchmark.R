args <- commandArgs(trailingOnly = TRUE)
events <- if (length(args)) as.integer(args[[1]]) else 20000L
output <- if (length(args) >= 2L) args[[2]] else tempfile("analysis-v2-benchmark-", fileext = ".csv")
if (!is.finite(events) || events < 1000L) stop("EVENTS must be at least 1000.", call. = FALSE)

repository <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "..", ".."))
api <- new.env(parent = globalenv())
source(file.path(repository, "inst", "api", "gui_api.R"), local = api)

project <- tempfile("spectreasy-analysis-performance-")
on.exit(unlink(project, recursive = TRUE, force = TRUE), add = TRUE)
dir.create(file.path(project, "samples"), recursive = TRUE)

set.seed(20260723)
group <- rep(seq_len(4L), length.out = events)
centers <- matrix(c(
    -4, -3, 0, 1, 2, 0, 1, -1, 0, 2, -2, 1,
    0, 5, 4, -2, 0, 1, -2, 2, 3, 0, 1, -1,
    5, -1, -3, 4, 3, -2, 1, 3, -1, -2, 2, 4,
    1, 2, 2, 0, -4, 5, 3, -2, 4, 1, 0, -3
), nrow = 4L, byrow = TRUE)
matrix <- centers[group, , drop = FALSE] + matrix(stats::rnorm(events * 12L, sd = 0.7), events, 12L)
colnames(matrix) <- paste0("M", seq_len(12L), "-A")
frame <- flowCore::flowFrame(matrix)
parameters <- flowCore::parameters(frame)@data
parameters$desc <- paste0("M", seq_len(12L))
flowCore::parameters(frame)@data <- parameters
fcs <- file.path(project, "samples", "benchmark.fcs")
flowCore::write.FCS(frame, fcs)

run <- function(name, cluster, reduction) {
    gc()
    started <- proc.time()[["elapsed"]]
    result <- api$gui_analysis_run_method(project, list(
        file = "samples/benchmark.fcs",
        populationId = "root",
        clusterMethod = cluster,
        reductionMethod = reduction,
        markers = colnames(matrix),
        maxEvents = events,
        seed = 20260723L,
        neighbors = 15L,
        clusters = 4L
    ))
    elapsed <- proc.time()[["elapsed"]] - started
    stopifnot(nrow(result$events) == events)
    data.frame(
        events = events,
        markers = ncol(matrix),
        benchmark = name,
        clustering = cluster,
        reduction = reduction,
        seconds = elapsed,
        events_per_second = events / elapsed,
        artifact_reused = isTRUE(result$metadata$artifacts$embedding$reused),
        stringsAsFactors = FALSE
    )
}

results <- rbind(
    run("pca", "none", "pca"),
    run("flowsom-pca", "flowsom", "pca"),
    if (events <= 20000L) run("umap", "none", "umap") else NULL
)
utils::write.csv(results, output, row.names = FALSE)
print(results, row.names = FALSE)
cat(normalizePath(output), "\n")
