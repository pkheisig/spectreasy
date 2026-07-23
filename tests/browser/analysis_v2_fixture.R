args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: analysis_v2_fixture.R PROJECT", call. = FALSE)

project <- normalizePath(args[[1]], mustWork = FALSE)
unlink(file.path(project, ".spectreasy"), recursive = TRUE, force = TRUE)
unlink(file.path(project, "spectreasy_outputs"), recursive = TRUE, force = TRUE)
samples <- file.path(project, "samples")
dir.create(samples, recursive = TRUE, showWarnings = FALSE)

set.seed(20260723)
events <- 600L
group <- rep(seq_len(3L), each = events / 3L)
data <- cbind(
    `FSC-A` = stats::rnorm(events, c(72, 105, 138)[group], 10),
    `SSC-A` = stats::rnorm(events, c(42, 82, 118)[group], 8),
    `CD3-A` = stats::rnorm(events, c(1, 8, 12)[group], 0.8),
    `CD19-A` = stats::rnorm(events, c(11, 4, 1)[group], 0.8),
    `CD56-A` = stats::rnorm(events, c(2, 3, 11)[group], 0.8)
)
frame <- flowCore::flowFrame(data)
parameters <- flowCore::parameters(frame)@data
parameters$desc <- c("FSC", "SSC", "CD3", "CD19", "CD56")
flowCore::parameters(frame)@data <- parameters
flowCore::write.FCS(frame, file.path(samples, "PKH_browser_sample.fcs"))

set.seed(20260724)
second_data <- data + matrix(stats::rnorm(length(data), 0, 0.12), nrow(data), ncol(data))
second <- flowCore::flowFrame(second_data)
second_parameters <- flowCore::parameters(second)@data
second_parameters$desc <- c("FSC", "SSC", "CD3", "CD19", "CD56")
flowCore::parameters(second)@data <- second_parameters
flowCore::write.FCS(second, file.path(samples, "PKH_browser_sample_repeat.fcs"))
