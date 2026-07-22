args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("usage: analysis_v2_server.R PORT REPOSITORY PROJECT", call. = FALSE)
port <- as.integer(args[[1]])
repository <- normalizePath(args[[2]], mustWork = TRUE)
project <- normalizePath(args[[3]], mustWork = TRUE)

devtools::load_all(repository, quiet = TRUE)
options(
    spectreasy.project_dir = project,
    spectreasy.matrix_dir = project,
    spectreasy.samples_dir = file.path(project, "samples"),
    spectreasy.gating_scc_dir = file.path(project, "scc"),
    spectreasy.project_selected = TRUE,
    spectreasy.gui_mode = "cockpit",
    spectreasy.gui_api_token = "analysis-v2-smoke-token",
    spectreasy.gui_allowed_origins = paste0("http://127.0.0.1:", port)
)
api <- plumber::plumb(file.path(repository, "inst", "api", "gui_api.R"))
api <- plumber::pr_static(api, "/", file.path(repository, "inst", "gui", "dist"))
api$run(port = port, host = "127.0.0.1", docs = FALSE, quiet = TRUE)
