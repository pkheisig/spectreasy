#' Inspect available population-analysis methods
#'
#' Returns the same runtime-aware method registry used by the Spectreasy GUI.
#' Each entry includes availability, prerequisites, outputs, and the complete
#' advanced-parameter schema.
#'
#' @param refresh Re-probe optional Python analysis dependencies.
#'
#' @return A list of method definitions.
#' @export
analysis_methods <- function(refresh = FALSE) {
    backend <- .spectreasy_analysis_backend()
    if (isTRUE(refresh)) {
        rm(list = ls(backend$.gui_analysis_python_cache, all.names = TRUE),
           envir = backend$.gui_analysis_python_cache)
    }
    backend$gui_analysis_method_registry()
}

#' Inspect the managed Python analysis runtime
#'
#' @return A list containing the selected interpreter, lock file, and adapter
#'   availability reported by the same bridge used for analysis.
#' @export
analysis_runtime_status <- function() {
    backend <- .spectreasy_analysis_backend()
    list(
        python = backend$gui_analysis_python_executable(),
        lock_file = .analysis_runtime_lock_file(),
        packages = backend$gui_analysis_python_packages(refresh = TRUE)
    )
}

#' Install the managed Python analysis runtime
#'
#' Creates a package-owned user-cache virtual environment and installs the
#' fully pinned analysis dependency snapshot. Nothing is installed into the
#' system Python environment.
#'
#' @param python Bootstrap Python executable (3.11 or newer).
#' @param path Optional virtual-environment directory.
#' @param rebuild Remove and recreate an existing managed environment.
#'
#' @return The result of [analysis_runtime_status()].
#' @export
install_analysis_runtime <- function(
    python = Sys.which("python3"),
    path = NULL,
    rebuild = FALSE
) {
    python <- .analysis_scalar_character(python, "python")
    if (!file.exists(python)) stop("Bootstrap Python was not found: ", python, call. = FALSE)
    version_output <- suppressWarnings(system2(
        python,
        c("-c", shQuote("import sys; print('.'.join(map(str, sys.version_info[:3])))")),
        stdout = TRUE,
        stderr = TRUE
    ))
    version_status <- attr(version_output, "status") %||% 0L
    python_version <- tryCatch(package_version(trimws(version_output[[1]])), error = function(e) NULL)
    if (!identical(as.integer(version_status), 0L) || is.null(python_version) || python_version < package_version("3.11")) {
        stop("Population analysis requires Python 3.11 or newer; selected interpreter: ", python, call. = FALSE)
    }
    if (is.null(path)) {
        cache_root <- if (identical(Sys.info()[["sysname"]], "Darwin")) {
            file.path(path.expand("~/Library/Caches"), "spectreasy")
        } else {
            file.path(path.expand(Sys.getenv("XDG_CACHE_HOME", "~/.cache")), "spectreasy")
        }
        path <- file.path(cache_root, "analysis-python")
    }
    path <- path.expand(.analysis_scalar_character(path, "path"))
    if (isTRUE(rebuild) && dir.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
    if (!dir.exists(path)) {
        dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
        status <- system2(python, c("-m", "venv", shQuote(path)))
        if (!identical(status, 0L)) stop("Could not create the managed Python environment.", call. = FALSE)
    }
    managed_python <- file.path(path, if (.Platform$OS.type == "windows") "Scripts/python.exe" else "bin/python")
    lock_file <- .analysis_runtime_lock_file()
    status <- system2(managed_python, c("-m", "pip", "install", "--disable-pip-version-check", "-r", shQuote(lock_file)))
    if (!identical(status, 0L)) stop("Managed analysis dependency installation failed.", call. = FALSE)
    options(spectreasy.analysis_python = managed_python)
    rm(list = ls(.spectreasy_analysis_backend_cache, all.names = TRUE), envir = .spectreasy_analysis_backend_cache)
    status <- analysis_runtime_status()
    missing <- names(Filter(function(item) !isTRUE(item$available), status$packages))
    if (length(setdiff(missing, "spectreasy_builtin"))) {
        stop("Managed runtime installed, but adapters failed validation: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    status
}

.analysis_dependency_plan <- function() {
    list(
        cran = c("Rtsne", "uwot"),
        bioconductor = c("DelayedMatrixStats", "destiny", "FlowSOM", "slingshot", "TSCAN"),
        github = c(Rphenograph = "JinmiaoChenLab/Rphenograph@0298487f0ee13aac55eb77d19992f6bd878ba2fc")
    )
}

#' Install optional population-analysis dependencies
#'
#' Installs the complete optional R and managed Python runtime used by the
#' population-analysis GUI and code API. System Python is never modified.
#' Already installed packages are left untouched unless `reinstall` is true.
#'
#' @param include_python Install the pinned managed Python runtime.
#' @param python Bootstrap Python executable (3.11 or newer).
#' @param ask Pass interactive update questions to Bioconductor.
#' @param reinstall Reinstall optional R packages and rebuild Python.
#'
#' @return A list containing the installation plan, installed package versions,
#'   managed Python status, and method availability.
#' @export
install_analysis_dependencies <- function(
    include_python = TRUE,
    python = Sys.which("python3"),
    ask = FALSE,
    reinstall = FALSE
) {
    plan <- .analysis_dependency_plan()
    install_missing <- function(packages) {
        if (isTRUE(reinstall)) packages else packages[
            !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
        ]
    }

    cran <- install_missing(plan$cran)
    if (length(cran)) {
        utils::install.packages(cran, repos = "https://cloud.r-project.org")
    }

    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        utils::install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    bioconductor <- install_missing(plan$bioconductor)
    if (length(bioconductor)) {
        BiocManager::install(bioconductor, ask = ask, update = FALSE, force = isTRUE(reinstall))
    }

    github_names <- names(plan$github)
    github <- github_names[isTRUE(reinstall) | !vapply(github_names, requireNamespace, logical(1), quietly = TRUE)]
    if (length(github)) {
        if (!requireNamespace("remotes", quietly = TRUE)) {
            utils::install.packages("remotes", repos = "https://cloud.r-project.org")
        }
        for (package in github) {
            remotes::install_github(unname(plan$github[[package]]), upgrade = "never", force = isTRUE(reinstall))
        }
    }

    runtime <- if (isTRUE(include_python)) {
        install_analysis_runtime(python = python, rebuild = isTRUE(reinstall))
    } else {
        analysis_runtime_status()
    }
    packages <- unique(c(plan$cran, plan$bioconductor, names(plan$github)))
    versions <- stats::setNames(vapply(packages, function(package) {
        if (requireNamespace(package, quietly = TRUE)) as.character(utils::packageVersion(package)) else NA_character_
    }, character(1)), packages)
    methods <- analysis_methods(refresh = TRUE)
    list(
        plan = plan,
        package_versions = versions,
        python = runtime,
        methods = methods,
        ready = all(vapply(methods, function(method) isTRUE(method$available), logical(1)))
    )
}

.analysis_runtime_lock_file <- function() {
    candidates <- c(
        file.path(getwd(), "inst", "python", "requirements-analysis.lock"),
        file.path(getwd(), "..", "inst", "python", "requirements-analysis.lock"),
        system.file("python", "requirements-analysis.lock", package = "spectreasy")
    )
    hit <- candidates[file.exists(candidates)][1]
    if (is.na(hit) || !nzchar(hit)) stop("The managed analysis runtime lock file is unavailable.", call. = FALSE)
    normalizePath(hit, mustWork = TRUE)
}

#' Read or save the sample-analysis gating workspace
#'
#' @param project_path Spectreasy project directory.
#' @param workspace A workspace list returned by `analysis_workspace()`.
#'
#' @return The normalized workspace.
#' @export
analysis_workspace <- function(project_path) {
    project_path <- .analysis_project_path(project_path)
    .spectreasy_analysis_backend()$gui_analysis_read_workspace(project_path)
}

#' @rdname analysis_workspace
#' @export
save_analysis_workspace <- function(project_path, workspace) {
    project_path <- .analysis_project_path(project_path)
    .spectreasy_analysis_backend()$gui_analysis_write_workspace(project_path, workspace)
}

#' Export or replace reusable population gates as CSV
#'
#' Gate-set CSV files contain only population IDs, names, parent relationships,
#' semantic roles, optional file scope, channels, and gate geometry. They do not
#' contain plots, annotations, analysis settings, seeds, or trajectory state;
#' use [save_analysis_workspace()] for the complete internal workspace instead.
#'
#' `import_population_gates()` validates the CSV and its referenced channels and
#' source files against the active project, then replaces the current gate
#' hierarchy. Plot layout and other analysis settings are retained, while
#' references to populations that no longer exist are repaired or removed.
#'
#' Schema version 1 records `format`, `schema_version`, `record_type`,
#' `population_id`, `parent_id`, `population_name`, `gate_type`,
#' `semantic_role`, `source_file`, `x_channel`, `y_channel`, `vertex_index`,
#' `x`, and `y`. Polygon vertices use one row each; rectangles and ellipses use
#' two bounding-box corners; range gates use two X coordinates. The root
#' population is recreated automatically.
#'
#' @param project_path Spectreasy project directory.
#' @param path CSV file to write or read.
#'
#' @return `export_population_gates()` returns the saved path and gate count.
#'   `import_population_gates()` returns the normalized, saved workspace with
#'   import warnings and a replacement summary in attributes named `warnings`
#'   and `import_summary`.
#' @export
export_population_gates <- function(project_path, path) {
    project_path <- .analysis_project_path(project_path)
    path <- .analysis_scalar_character(path, "path")
    backend <- .spectreasy_analysis_backend()
    result <- backend$gui_analysis_write_gate_set(
        project_path,
        path,
        workspace = backend$gui_analysis_read_workspace(project_path)
    )
    c(result, schema_version = backend$gui_analysis_gate_set_schema_version())
}

#' @rdname export_population_gates
#' @export
import_population_gates <- function(project_path, path) {
    project_path <- .analysis_project_path(project_path)
    path <- .analysis_scalar_character(path, "path")
    backend <- .spectreasy_analysis_backend()
    current <- backend$gui_analysis_read_workspace(project_path)
    prepared <- backend$gui_analysis_read_gate_set_file(project_path, path, workspace = current)
    saved <- backend$gui_analysis_write_workspace(project_path, prepared$workspace)
    attr(saved, "warnings") <- prepared$warnings
    attr(saved, "import_summary") <- prepared$summary
    saved
}

#' Create a population gate in code
#'
#' @param project_path Spectreasy project directory.
#' @param name Population name.
#' @param parent_id Parent population identifier.
#' @param type One of `"rectangle"`, `"ellipse"`, `"polygon"`, or `"range"`.
#' @param x,y FCS channels. `y` is not used for range gates.
#' @param geometry Named gate geometry in raw channel coordinates.
#' @param role Optional semantic role: `"positive"`, `"negative"`, `"root"`,
#'   or `"terminal"`.
#' @param source_file Optional project-relative FCS path for a file-specific gate.
#' @param id Optional stable population identifier.
#'
#' @return The saved workspace.
#' @export
add_population_gate <- function(
    project_path,
    name,
    parent_id = "root",
    type,
    x,
    y = NULL,
    geometry,
    role = NULL,
    source_file = NULL,
    id = NULL
) {
    project_path <- .analysis_project_path(project_path)
    backend <- .spectreasy_analysis_backend()
    workspace <- backend$gui_analysis_read_workspace(project_path)
    parent_id <- .analysis_scalar_character(parent_id, "parent_id")
    backend$gui_analysis_population(workspace, parent_id)
    type <- match.arg(type, c("rectangle", "ellipse", "polygon", "range"))
    name <- .analysis_scalar_character(name, "name")
    x <- .analysis_scalar_character(x, "x")
    if (!identical(type, "range")) y <- .analysis_scalar_character(y, "y")
    if (!is.list(geometry)) stop("`geometry` must be a named list.", call. = FALSE)
    role <- if (is.null(role) || !length(role) || !nzchar(role[[1]])) NULL else match.arg(as.character(role[[1]]), c("positive", "negative", "root", "terminal"))
    id <- if (is.null(id)) {
        paste0("population-", substr(backend$gui_analysis_stable_id("gate", name, parent_id, as.numeric(Sys.time()), stats::runif(1)), 6L, 13L))
    } else {
        .analysis_scalar_character(id, "id")
    }
    if (id %in% vapply(workspace$populations, `[[`, character(1), "id")) {
        stop("Population identifier already exists: ", id, call. = FALSE)
    }
    workspace$populations <- c(workspace$populations, list(list(
        id = id, name = name, parent_id = parent_id, type = type,
        role = role, source_file = source_file, x = x, y = y, geometry = geometry
    )))
    backend$gui_analysis_write_workspace(project_path, workspace)
}

#' Update or delete a saved population gate
#'
#' @param project_path Spectreasy project directory.
#' @param population_id Population identifier.
#' @param ... Named population fields to update.
#'
#' @return The saved workspace.
#' @export
update_population_gate <- function(project_path, population_id, ...) {
    project_path <- .analysis_project_path(project_path)
    backend <- .spectreasy_analysis_backend()
    workspace <- backend$gui_analysis_read_workspace(project_path)
    population_id <- .analysis_scalar_character(population_id, "population_id")
    if (identical(population_id, "root")) stop("The root population cannot be edited.", call. = FALSE)
    index <- match(population_id, vapply(workspace$populations, `[[`, character(1), "id"))
    if (is.na(index)) stop("Population not found: ", population_id, call. = FALSE)
    patch <- list(...)
    allowed <- c("name", "role", "source_file", "x", "y", "geometry")
    unknown <- setdiff(names(patch), allowed)
    if (length(unknown)) stop("Unsupported population field: ", unknown[[1]], call. = FALSE)
    workspace$populations[[index]] <- utils::modifyList(workspace$populations[[index]], patch, keep.null = TRUE)
    backend$gui_analysis_write_workspace(project_path, workspace)
}

#' @rdname update_population_gate
#' @export
delete_population_gate <- function(project_path, population_id) {
    project_path <- .analysis_project_path(project_path)
    backend <- .spectreasy_analysis_backend()
    workspace <- backend$gui_analysis_read_workspace(project_path)
    population_id <- .analysis_scalar_character(population_id, "population_id")
    if (identical(population_id, "root")) stop("The root population cannot be deleted.", call. = FALSE)
    if (!population_id %in% vapply(workspace$populations, `[[`, character(1), "id")) {
        stop("Population not found: ", population_id, call. = FALSE)
    }
    removing <- population_id
    repeat {
        children <- vapply(workspace$populations, function(node) {
            if (as.character(node$parent_id %||% "")[1] %in% removing) as.character(node$id)[1] else ""
        }, character(1))
        next_ids <- setdiff(children[nzchar(children)], removing)
        if (!length(next_ids)) break
        removing <- c(removing, next_ids)
    }
    workspace$populations <- Filter(function(node) !as.character(node$id)[1] %in% removing, workspace$populations)
    workspace$plots <- lapply(workspace$plots, function(plot) {
        if (as.character(plot$population_id %||% "root")[1] %in% removing) plot$population_id <- "root"
        if (as.character(plot$overlay_population_id %||% "")[1] %in% removing) plot$overlay_population_id <- NULL
        plot
    })
    if (workspace$active_population_id %in% removing) workspace$active_population_id <- "root"
    if (as.character(workspace$root_population_id %||% "")[1] %in% removing) {
        workspace$root_event_id <- NULL
        workspace$root_population_id <- NULL
        workspace$root_source_file <- NULL
    }
    backend$gui_analysis_write_workspace(project_path, workspace)
}

#' Summarize or export a gated population in code
#'
#' @param project_path Spectreasy project directory.
#' @param file Project-relative FCS path.
#' @param population_id Population identifier.
#' @param markers Optional channels to summarize.
#' @param marker Marker channel used for the staining index.
#' @param files One or more project-relative FCS paths.
#' @param format `"fcs"`, `"csv"`, or `"both"`.
#' @param max_events Zero for all gated events, or a seeded maximum.
#' @param seed Master sampling seed.
#' @param output_folder Project-relative output folder.
#'
#' @return Population statistics or export metadata.
#' @export
population_statistics <- function(project_path, file, population_id = "root", markers = character()) {
    .spectreasy_analysis_backend()$gui_analysis_population_statistics(
        .analysis_project_path(project_path),
        .analysis_scalar_character(file, "file"),
        .analysis_scalar_character(population_id, "population_id"),
        markers
    )
}

#' @rdname population_statistics
#' @export
staining_index <- function(project_path, file, marker) {
    .spectreasy_analysis_backend()$gui_analysis_staining_index(
        .analysis_project_path(project_path),
        .analysis_scalar_character(file, "file"),
        .analysis_scalar_character(marker, "marker")
    )
}

#' @rdname population_statistics
#' @export
export_gated_population <- function(
    project_path,
    files,
    population_id = "root",
    format = "fcs",
    max_events = 0L,
    seed = 20260723L,
    output_folder = "spectreasy_outputs/analysis/exports"
) {
    .spectreasy_analysis_backend()$gui_analysis_export(
        .analysis_project_path(project_path),
        list(
            files = as.character(files),
            populationId = population_id,
            format = match.arg(format, c("fcs", "csv", "both")),
            maxEvents = max_events,
            seed = seed,
            outputFolder = output_folder
        )
    )
}

#' Run clustering, dimensional reduction, or trajectory analysis
#'
#' This is the code-facing equivalent of the population-analysis window. The
#' GUI calls the same backend function, so validation, seeded sampling, cached
#' fitted objects, provenance, and saved artifacts are identical in both
#' interfaces.
#'
#' @param project_path Spectreasy project directory.
#' @param file One or more project-relative FCS paths. Multiple files are
#'   pooled after applying the same saved population gate, with source file and
#'   original event row retained in the result.
#' @param markers FCS channel names used for analysis.
#' @param population_id Population identifier from the saved gating hierarchy.
#' @param clustering Clustering method, either `"none"`, `"flowsom"`, or
#'   `"phenograph"`.
#' @param reduction Dimensional-reduction method. Required unless `trajectory`
#'   supplies its own coordinates.
#' @param trajectory Optional trajectory method: `"dpt"`, `"slingshot"`,
#'   `"tscan"`, `"palantir"`, `"paga-dpt"`, `"wanderlust"`, or `"wishbone"`.
#' @param max_events Maximum number of deterministically sampled events.
#' @param cofactor Arcsinh cofactor.
#' @param seed Positive master seed.
#' @param root_event_id Original FCS row number used as the trajectory root.
#' @param cluster_settings Named advanced settings for `clustering`.
#' @param reduction_settings Named advanced settings for `reduction`.
#' @param trajectory_settings Named advanced settings for `trajectory`.
#' @param advanced_settings Optional named list keyed by method identifier. This
#'   is useful when constructing pipelines programmatically; method-specific
#'   arguments above take precedence.
#' @param neighbors Compatibility default used by methods without an explicit
#'   advanced neighbor setting.
#' @param clusters Compatibility default used by clustering adapters.
#' @param perplexity Compatibility default used by t-SNE.
#'
#' @return A `spectreasy_analysis` object containing metadata and event-level
#'   coordinates, marker values, cluster labels, branches, or pseudotime.
#' @export
analyze_population <- function(
    project_path,
    file,
    markers,
    population_id = "root",
    clustering = "none",
    reduction = NULL,
    trajectory = NULL,
    max_events = 20000L,
    cofactor = 150,
    seed = 20260723L,
    root_event_id = NULL,
    cluster_settings = list(),
    reduction_settings = list(),
    trajectory_settings = list(),
    advanced_settings = list(),
    neighbors = 15L,
    clusters = 10L,
    perplexity = 30
) {
    project_path <- .analysis_project_path(project_path)
    files <- unique(as.character(file))
    files <- files[nzchar(trimws(files))]
    if (!length(files)) stop("`file` must contain at least one FCS path.", call. = FALSE)
    markers <- unique(as.character(markers))
    markers <- markers[nzchar(trimws(markers))]
    if (length(markers) < 2L) stop("Select at least two marker channels.", call. = FALSE)
    clustering <- .analysis_method_id(clustering %||% "none", "clustering")
    if (!clustering %in% c("none", "flowsom", "phenograph")) {
        stop("Unsupported clustering method: ", clustering, call. = FALSE)
    }
    reduction <- if (is.null(reduction) || !length(reduction) || !nzchar(reduction[[1]])) {
        ""
    } else {
        .analysis_method_id(reduction, "reduction")
    }
    trajectory <- if (is.null(trajectory) || !length(trajectory) || !nzchar(trajectory[[1]])) {
        ""
    } else {
        .analysis_method_id(trajectory, "trajectory")
    }
    if (!nzchar(trajectory) && !nzchar(reduction)) {
        stop("A dimensional-reduction method is required.", call. = FALSE)
    }
    if (nzchar(trajectory) && !trajectory %in% c(
        "dpt", "slingshot", "tscan", "palantir", "paga-dpt", "wanderlust", "wishbone"
    )) {
        stop("Unsupported trajectory method: ", trajectory, call. = FALSE)
    }
    if (nzchar(reduction) && !reduction %in% c(
        "pca", "tsne", "umap", "diffusion-map", "phate", "hsne"
    )) {
        stop("Unsupported dimensional-reduction method: ", reduction, call. = FALSE)
    }
    settings <- .analysis_named_list(advanced_settings, "advanced_settings")
    if (!identical(clustering, "none")) {
        settings[[clustering]] <- utils::modifyList(
            settings[[clustering]] %||% list(),
            .analysis_named_list(cluster_settings, "cluster_settings")
        )
    }
    if (nzchar(reduction)) {
        settings[[reduction]] <- utils::modifyList(
            settings[[reduction]] %||% list(),
            .analysis_named_list(reduction_settings, "reduction_settings")
        )
    }
    if (nzchar(trajectory)) {
        settings[[trajectory]] <- utils::modifyList(
            settings[[trajectory]] %||% list(),
            .analysis_named_list(trajectory_settings, "trajectory_settings")
        )
    }
    request <- list(
        file = files[[1]],
        files = files,
        populationId = .analysis_scalar_character(population_id, "population_id"),
        clusterMethod = clustering,
        reductionMethod = reduction,
        markers = markers,
        maxEvents = max_events,
        cofactor = cofactor,
        seed = seed,
        neighbors = neighbors,
        clusters = clusters,
        perplexity = perplexity,
        advancedSettings = settings
    )
    if (nzchar(trajectory)) request$method <- trajectory
    if (!is.null(root_event_id)) request$rootEventId <- root_event_id
    result <- .spectreasy_analysis_run_request(project_path, request)
    class(result) <- unique(c("spectreasy_analysis", class(result)))
    result
}

#' Load a saved population-analysis result
#'
#' @param project_path Spectreasy project directory.
#' @param analysis_id Analysis identifier returned by [analyze_population()].
#' @param identity_id Optional identity-annotation identifier to merge.
#'
#' @return A `spectreasy_analysis` object.
#' @export
load_population_analysis <- function(project_path, analysis_id, identity_id = NULL) {
    project_path <- .analysis_project_path(project_path)
    analysis_id <- .analysis_artifact_id(analysis_id, "analysis_id")
    directory <- file.path(project_path, "spectreasy_outputs", "analysis", analysis_id)
    metadata_path <- file.path(directory, "metadata.json")
    events_path <- file.path(directory, "events.csv")
    if (!file.exists(metadata_path) || !file.exists(events_path)) {
        stop("Analysis result not found: ", analysis_id, call. = FALSE)
    }
    metadata <- jsonlite::fromJSON(metadata_path, simplifyVector = FALSE)
    events <- data.table::fread(events_path, data.table = FALSE)
    if (!is.null(identity_id)) {
        identity_id <- .analysis_artifact_id(identity_id, "identity_id")
        identity_directory <- file.path(directory, "identities", identity_id)
        identity_metadata_path <- file.path(identity_directory, "metadata.json")
        scores_path <- file.path(identity_directory, "scores.csv")
        if (!file.exists(identity_metadata_path) || !file.exists(scores_path)) {
            stop("Identity annotation not found: ", identity_id, call. = FALSE)
        }
        identity_metadata <- jsonlite::fromJSON(identity_metadata_path, simplifyVector = FALSE)
        scores <- data.table::fread(scores_path, data.table = FALSE)
        events <- merge(events, scores, by = "event_id", all.x = TRUE, sort = FALSE)
        metadata$identity_annotation <- identity_metadata
    }
    result <- list(
        metadata = metadata,
        events = events,
        metadata_file = file.path("spectreasy_outputs", "analysis", analysis_id, "metadata.json")
    )
    class(result) <- c("spectreasy_analysis", "list")
    result
}

#' Discover distinguishing markers for every cluster
#'
#' Performs a cluster-versus-rest rank-AUC comparison on the transformed marker
#' values stored with a completed population analysis. This is the
#' flow-cytometry analogue of marker discovery: it reports what distinguishes
#' each cluster, but does not pretend that marker discovery alone supplies a
#' biological cell-type name.
#'
#' @param analysis A `spectreasy_analysis` object or analysis identifier.
#' @param project_path Project directory; required when `analysis` is an ID or
#'   an in-memory result.
#' @param top_n Maximum number of positive markers returned per cluster.
#' @param minimum_auc Minimum cluster-versus-rest rank AUC.
#'
#' @return Marker-discovery metadata and a data frame of ranked cluster markers.
#' @export
find_population_markers <- function(
    analysis,
    project_path = NULL,
    top_n = 10L,
    minimum_auc = 0.55
) {
    if (inherits(analysis, "spectreasy_analysis")) {
        analysis_id <- .analysis_scalar_character(analysis$metadata$analysis_id, "analysis_id")
    } else {
        analysis_id <- .analysis_artifact_id(analysis, "analysis")
    }
    project_path <- .analysis_project_path(project_path)
    .spectreasy_analysis_marker_request(project_path, list(
        analysisId = analysis_id,
        topN = top_n,
        minimumAuc = minimum_auc
    ))
}

#' Build panel-aware common immune-cell identity patterns
#'
#' Creates editable marker patterns for common CD4 T, CD8 T, B, NK, and
#' monocyte identities. An identity is returned only when every required
#' positive marker group is represented in `markers`. Alternative marker groups
#' such as CD19/CD20 use the first available marker. These patterns are a
#' starting point and must be reviewed for the tissue, panel, and experiment.
#'
#' @param markers Marker display names available in an analysis result. Area,
#'   height, and width suffixes such as `-A` are ignored while matching, but the
#'   returned patterns preserve the supplied names.
#'
#' @return A list of editable identity signatures accepted by
#'   [annotate_population()].
#' @export
population_identity_templates <- function(markers) {
    markers <- unique(as.character(markers))
    markers <- markers[nzchar(markers)]
    if (!length(markers)) stop("`markers` must contain at least one marker name.", call. = FALSE)
    canonical <- function(value) {
        value <- toupper(gsub("\\s+", "", value))
        value <- sub("-(A|H|W)$", "", value)
        gsub("[^A-Z0-9]", "", value)
    }
    available <- stats::setNames(markers, canonical(markers))
    first_available <- function(group) {
        hits <- available[canonical(group)]
        hits <- unname(hits[!is.na(hits) & nzchar(hits)])
        if (length(hits)) hits[[1]] else NULL
    }
    definitions <- list(
        list(name = "CD4 T cell", positive = list("CD3", "CD4"), negative = list("CD8", c("CD19", "CD20"), "CD56", "CD14")),
        list(name = "CD8 T cell", positive = list("CD3", "CD8"), negative = list("CD4", c("CD19", "CD20"), "CD56", "CD14")),
        list(name = "B cell", positive = list(c("CD19", "CD20")), negative = list("CD3", "CD56", "CD14")),
        list(name = "NK cell", positive = list("CD56"), negative = list("CD3", c("CD19", "CD20"), "CD14")),
        list(name = "Monocyte", positive = list(c("CD14", "CD16")), negative = list("CD3", c("CD19", "CD20"), "CD56"))
    )
    colors <- c("#197783", "#d06d32", "#7d58a6", "#4f8f45", "#bc4d65")
    Filter(Negate(is.null), lapply(seq_along(definitions), function(index) {
        definition <- definitions[[index]]
        positive <- lapply(definition$positive, first_available)
        if (any(vapply(positive, is.null, logical(1)))) return(NULL)
        positive <- unique(unlist(positive, use.names = FALSE))
        negative <- unique(Filter(
            function(marker) !is.null(marker) && !marker %in% positive,
            lapply(definition$negative, first_available)
        ))
        list(
            name = definition$name,
            color = colors[[index]],
            positive_markers = positive,
            negative_markers = unlist(negative, use.names = FALSE)
        )
    }))
}

#' Annotate cells with signed flow-cytometry marker signatures
#'
#' @param analysis A `spectreasy_analysis` object or analysis identifier.
#' @param project_path Project directory; required when `analysis` is an ID.
#' @param signatures List of identity definitions. Each definition requires
#'   `name`, `positive_markers`, and/or `negative_markers`; `color` is optional.
#' @param min_score Minimum winning signed-marker score.
#' @param min_margin Minimum separation from the runner-up identity.
#' @param evidence_sensitivity Positive scaling applied to robust marker
#'   evidence.
#'
#' @return Annotation metadata and event-level identity scores.
#' @export
annotate_population <- function(
    analysis,
    project_path = NULL,
    signatures,
    min_score = 0.55,
    min_margin = 0.08,
    evidence_sensitivity = 1
) {
    if (inherits(analysis, "spectreasy_analysis")) {
        analysis_id <- .analysis_scalar_character(analysis$metadata$analysis_id, "analysis_id")
        if (is.null(project_path)) {
            metadata_file <- analysis$metadata_file %||% ""
            stop(
                "`project_path` is required when annotating an in-memory result.",
                if (nzchar(metadata_file)) paste0(" Result metadata: ", metadata_file) else "",
                call. = FALSE
            )
        }
    } else {
        analysis_id <- .analysis_artifact_id(analysis, "analysis")
    }
    project_path <- .analysis_project_path(project_path)
    .spectreasy_analysis_annotate_request(project_path, list(
        analysisId = analysis_id,
        signatures = signatures,
        minScore = min_score,
        minMargin = min_margin,
        evidenceSensitivity = evidence_sensitivity
    ))
}

#' Plot a population-analysis result
#'
#' Produces a square `ggplot2` scatter plot in two dimensions or an interactive
#' Plotly plot in three dimensions. Coordinates, continuous marker coloring,
#' clusters, pseudotime, and predicted identities can all be selected directly.
#'
#' @param analysis A `spectreasy_analysis` object.
#' @param x,y,z Coordinate numbers, stored column names, or displayed coordinate
#'   labels. Supplying `z` requests a 3D Plotly result.
#' @param color_by One of `"density"`, `"cluster"`, `"pseudotime"`,
#'   `"identity"`, a marker display name, or a stored marker column.
#' @param palette Continuous palette name. See
#'   `names(population_analysis_palettes())`.
#' @param point_size Point size.
#' @param alpha Point opacity.
#' @param title Optional plot title.
#' @param background Plot background color.
#'
#' @return A `ggplot` or `plotly` object.
#' @export
plot_population_analysis <- function(
    analysis,
    x = 1L,
    y = 2L,
    z = NULL,
    color_by = "density",
    palette = "control-density",
    point_size = 0.7,
    alpha = 0.82,
    title = NULL,
    background = "#f8f7f3"
) {
    analysis <- .analysis_result(analysis)
    coordinates <- grep("^dimension_[0-9]+$", names(analysis$events), value = TRUE)
    if (length(coordinates) < 2L) stop("Analysis result has fewer than two coordinates.", call. = FALSE)
    x_name <- .analysis_coordinate(analysis, x, coordinates)
    y_name <- .analysis_coordinate(analysis, y, coordinates)
    z_name <- if (is.null(z)) NULL else .analysis_coordinate(analysis, z, coordinates)
    color <- .analysis_color_values(analysis, color_by, x_name, y_name)
    palette_values <- population_analysis_palettes()[[palette]]
    if (is.null(palette_values)) stop("Unknown analysis palette: ", palette, call. = FALSE)
    data <- analysis$events
    data$.analysis_color <- color$values
    data$.analysis_x <- data[[x_name]]
    data$.analysis_y <- data[[y_name]]
    if (!is.null(z_name)) data$.analysis_z <- data[[z_name]]
    if (!is.null(z_name)) {
        if (!requireNamespace("plotly", quietly = TRUE)) {
            stop("Install the suggested 'plotly' package for 3D analysis plots.", call. = FALSE)
        }
        marker <- list(size = point_size * 4, opacity = alpha)
        if (color$continuous) {
            marker$color <- data$.analysis_color
            marker$colorscale <- .analysis_plotly_palette(palette_values)
            marker$showscale <- TRUE
            marker$colorbar <- list(title = color$label)
        } else {
            marker$color <- .analysis_discrete_colors(data$.analysis_color)
        }
        return(plotly::plot_ly(
            x = data$.analysis_x,
            y = data$.analysis_y,
            z = data$.analysis_z,
            type = "scatter3d", mode = "markers", marker = marker,
            text = paste("Event", data$event_id), hoverinfo = "text"
        ) |> plotly::layout(
            title = title %||% analysis$metadata$display_name,
            paper_bgcolor = background,
            scene = list(
                aspectmode = "cube",
                xaxis = list(title = .analysis_coordinate_label(analysis, x_name)),
                yaxis = list(title = .analysis_coordinate_label(analysis, y_name)),
                zaxis = list(title = .analysis_coordinate_label(analysis, z_name))
            )
        ))
    }
    mapping <- ggplot2::aes(
        x = .analysis_x, y = .analysis_y, color = .analysis_color
    )
    plot <- ggplot2::ggplot(data, mapping) +
        ggplot2::geom_point(size = point_size, alpha = alpha, stroke = 0) +
        ggplot2::labs(
            title = title %||% analysis$metadata$display_name,
            x = .analysis_coordinate_label(analysis, x_name),
            y = .analysis_coordinate_label(analysis, y_name),
            color = color$label
        ) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
            aspect.ratio = 1,
            panel.grid = ggplot2::element_blank(),
            panel.background = ggplot2::element_rect(fill = background, color = "#c7c3ba"),
            plot.background = ggplot2::element_rect(fill = background, color = NA)
        )
    if (color$continuous) {
        plot + ggplot2::scale_color_gradientn(colours = palette_values)
    } else {
        plot + ggplot2::scale_color_manual(values = .analysis_discrete_palette(data$.analysis_color))
    }
}

#' Export population-analysis tables, objects, and plots
#'
#' @param analysis A `spectreasy_analysis` object.
#' @param file Output filename.
#' @param format Optional format inferred from `file`: `"csv"`, `"rds"`,
#'   `"png"`, `"svg"`, `"pdf"`, or `"html"`.
#' @param width,height Plot dimensions in inches.
#' @param dpi Raster resolution.
#' @param ... Plot customizations passed to [plot_population_analysis()].
#'
#' @return The normalized output path, invisibly.
#' @export
export_population_analysis <- function(
    analysis,
    file,
    format = tools::file_ext(file),
    width = 6,
    height = 6,
    dpi = 300,
    ...
) {
    analysis <- .analysis_result(analysis)
    file <- normalizePath(path.expand(file), mustWork = FALSE)
    format <- tolower(.analysis_scalar_character(format, "format"))
    if (!format %in% c("csv", "rds", "png", "svg", "pdf", "html")) {
        stop("Unsupported analysis export format: ", format, call. = FALSE)
    }
    directory <- dirname(file)
    if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE, showWarnings = FALSE)) {
        stop("Could not create export directory: ", directory, call. = FALSE)
    }
    if (identical(format, "csv")) {
        data.table::fwrite(.analysis_export_table(analysis), file)
    } else if (identical(format, "rds")) {
        saveRDS(analysis, file)
    } else {
        plot <- plot_population_analysis(analysis, ...)
        if (identical(format, "html")) {
            if (!inherits(plot, "plotly")) {
                if (!requireNamespace("plotly", quietly = TRUE)) {
                    stop("Install the suggested 'plotly' package for HTML exports.", call. = FALSE)
                }
                plot <- plotly::ggplotly(plot)
            }
            if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
                stop("Install the suggested 'htmlwidgets' package for HTML exports.", call. = FALSE)
            }
            htmlwidgets::saveWidget(plot, file, selfcontained = TRUE)
        } else {
            if (inherits(plot, "plotly")) {
                stop("Static 3D export requires Plotly's external image renderer; use HTML instead.", call. = FALSE)
            }
            ggplot2::ggsave(file, plot = plot, width = width, height = height, dpi = dpi)
        }
    }
    invisible(file)
}

#' Continuous palettes for population-analysis plots
#'
#' @return A named list of hexadecimal color vectors.
#' @export
population_analysis_palettes <- function() {
    list(
        `control-density` = c("#293378", "#28579d", "#2781a3", "#42a56f", "#96c93d", "#e8df2e", "#f59f25", "#d94727"),
        viridis = c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725"),
        sunset = c("#30123b", "#6d2f84", "#b73779", "#eb5a5a", "#f89f5b", "#f9e78a"),
        plasma = c("#0d0887", "#7e03a8", "#cc4778", "#f89540", "#f0f921"),
        inferno = c("#000004", "#420a68", "#932667", "#dd513a", "#fca50a", "#fcffa4"),
        magma = c("#000004", "#3b0f70", "#8c2981", "#de4968", "#fe9f6d", "#fcfdbf"),
        cividis = c("#00204c", "#2e4a7d", "#666970", "#a48b58", "#e1c63b", "#ffea46"),
        turbo = c("#30123b", "#4669e8", "#1bcfd4", "#75fe5c", "#f9e721", "#f77b16", "#7a0403"),
        `ice-fire` = c("#00204c", "#2676b8", "#9bd7e4", "#f5f2ea", "#f4a582", "#c83b3b", "#5d0018"),
        spectral = c("#5e4fa2", "#3288bd", "#66c2a5", "#e6f598", "#fee08b", "#f46d43", "#9e0142")
    )
}

.spectreasy_analysis_backend_cache <- new.env(parent = emptyenv())

.spectreasy_analysis_backend <- function() {
    key <- "shared"
    if (exists(key, envir = .spectreasy_analysis_backend_cache, inherits = FALSE)) {
        return(get(key, envir = .spectreasy_analysis_backend_cache, inherits = FALSE))
    }
    candidates <- c(
        file.path(getwd(), "inst", "api", "modules", "analysis_workspace.R"),
        file.path(getwd(), "..", "inst", "api", "modules", "analysis_workspace.R"),
        system.file("api", "modules", "analysis_workspace.R", package = "spectreasy")
    )
    candidates <- unique(candidates[nzchar(candidates)])
    module <- candidates[file.exists(candidates)][1]
    if (is.na(module) || !nzchar(module)) {
        stop("Could not locate the Spectreasy population-analysis backend.", call. = FALSE)
    }
    backend <- new.env(parent = environment())
    backend$gui_request_project_root <- function(project_path = "", fallback = TRUE) {
        value <- trimws(as.character(project_path %||% "")[1])
        if (!nzchar(value) && isTRUE(fallback)) value <- getOption("spectreasy.project_dir", "")
        .analysis_project_path(value)
    }
    backend$gui_project_layout <- function(root, reconcile = TRUE, ensure_markers = TRUE, persist = TRUE) {
        .spectreasy_project_layout(
            backend$gui_request_project_root(root),
            reconcile = reconcile,
            ensure_markers = ensure_markers,
            persist = persist
        )
    }
    backend$gui_workflow_value <- function(body, key, fallback = NULL) {
        value <- body[[key]]
        if (is.null(value) || !length(value) || !nzchar(trimws(as.character(value[[1]])))) fallback else as.character(value[[1]])
    }
    backend$gui_workflow_number <- function(body, key, fallback = 0, integer = FALSE, minimum = NULL, maximum = NULL) {
        supplied <- !is.null(body[[key]]) && length(body[[key]])
        value <- suppressWarnings(as.numeric(backend$gui_workflow_value(body, key, fallback))[[1]])
        if (!is.finite(value)) {
            if (supplied) stop("Invalid numeric value for '", key, "'.", call. = FALSE)
            value <- fallback
        }
        if (!is.null(minimum) && value < minimum) stop("'", key, "' must be at least ", minimum, ".", call. = FALSE)
        if (!is.null(maximum) && value > maximum) stop("'", key, "' must be at most ", maximum, ".", call. = FALSE)
        if (isTRUE(integer)) {
            if (!isTRUE(all.equal(value, round(value)))) stop("'", key, "' must be an integer.", call. = FALSE)
            value <- as.integer(value)
        }
        value
    }
    backend$gui_workflow_path <- function(body, key, fallback = "", allow_empty = TRUE) {
        value <- trimws(as.character(backend$gui_workflow_value(body, key, fallback))[[1]])
        if (!nzchar(value) && !isTRUE(allow_empty)) fallback else value
    }
    backend$gui_workflow_resolve_path <- function(path, root, allow_empty = FALSE) {
        path <- trimws(as.character(path %||% "")[1])
        if (!nzchar(path)) {
            if (isTRUE(allow_empty)) return("")
            stop("Project-relative path is empty.", call. = FALSE)
        }
        root <- .analysis_project_path(root)
        candidate <- normalizePath(if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) path else file.path(root, path), mustWork = FALSE)
        prefix <- paste0(root, .Platform$file.sep)
        if (!identical(candidate, root) && !startsWith(candidate, prefix)) {
            stop("Analysis path is outside the active project.", call. = FALSE)
        }
        candidate
    }
    sys.source(module, envir = backend)
    assign(key, backend, envir = .spectreasy_analysis_backend_cache)
    backend
}

.spectreasy_analysis_run_request <- function(project_path, request) {
    .spectreasy_analysis_backend()$gui_analysis_run_method(
        .analysis_project_path(project_path),
        request
    )
}

.spectreasy_analysis_annotate_request <- function(project_path, request) {
    .spectreasy_analysis_backend()$gui_analysis_annotate_result(
        .analysis_project_path(project_path),
        request
    )
}

.spectreasy_analysis_marker_request <- function(project_path, request) {
    .spectreasy_analysis_backend()$gui_analysis_find_cluster_markers(
        .analysis_project_path(project_path),
        request
    )
}

.analysis_project_path <- function(path) {
    if (is.null(path) || !length(path) || !nzchar(trimws(as.character(path[[1]])))) {
        stop("A Spectreasy project directory is required.", call. = FALSE)
    }
    path <- path.expand(as.character(path[[1]]))
    if (!dir.exists(path)) stop("Project directory not found: ", path, call. = FALSE)
    normalizePath(path, mustWork = TRUE)
}

.analysis_scalar_character <- function(value, name) {
    if (is.null(value) || length(value) != 1L || is.na(value[[1]]) || !nzchar(trimws(as.character(value[[1]])))) {
        stop("`", name, "` must be one non-empty value.", call. = FALSE)
    }
    trimws(as.character(value[[1]]))
}

.analysis_method_id <- function(value, name) {
    tolower(.analysis_scalar_character(value, name))
}

.analysis_named_list <- function(value, name) {
    if (is.null(value)) return(list())
    if (!is.list(value) || (length(value) && (is.null(names(value)) || any(!nzchar(names(value)))))) {
        stop("`", name, "` must be a named list.", call. = FALSE)
    }
    value
}

.analysis_artifact_id <- function(value, name) {
    value <- .analysis_scalar_character(value, name)
    if (!grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", value)) {
        stop("`", name, "` contains unsupported characters.", call. = FALSE)
    }
    value
}

.analysis_result <- function(value) {
    if (!inherits(value, "spectreasy_analysis") || !is.list(value) || !is.data.frame(value$events)) {
        stop("`analysis` must be a result returned by analyze_population() or load_population_analysis().", call. = FALSE)
    }
    value
}

.analysis_coordinate <- function(analysis, value, coordinates) {
    labels <- unlist(analysis$metadata$coordinate_labels %||% character(), use.names = FALSE)
    if (is.numeric(value) && length(value) == 1L && is.finite(value)) {
        index <- as.integer(value)
        if (index < 1L || index > length(coordinates)) stop("Coordinate index is unavailable.", call. = FALSE)
        return(coordinates[[index]])
    }
    value <- .analysis_scalar_character(value, "coordinate")
    if (value %in% coordinates) return(value)
    index <- match(value, labels)
    if (!is.na(index) && index <= length(coordinates)) return(coordinates[[index]])
    stop("Unknown analysis coordinate: ", value, call. = FALSE)
}

.analysis_coordinate_label <- function(analysis, coordinate) {
    coordinates <- grep("^dimension_[0-9]+$", names(analysis$events), value = TRUE)
    labels <- unlist(analysis$metadata$coordinate_labels %||% character(), use.names = FALSE)
    index <- match(coordinate, coordinates)
    if (!is.na(index) && index <= length(labels) && nzchar(labels[[index]])) labels[[index]] else coordinate
}

.analysis_color_values <- function(analysis, color_by, x, y) {
    color_by <- .analysis_scalar_character(color_by, "color_by")
    events <- analysis$events
    if (identical(color_by, "density")) {
        bins <- max(12L, min(60L, as.integer(round(sqrt(nrow(events)) / 2))))
        x_bin <- cut(events[[x]], breaks = bins, labels = FALSE, include.lowest = TRUE)
        y_bin <- cut(events[[y]], breaks = bins, labels = FALSE, include.lowest = TRUE)
        key <- paste(x_bin, y_bin, sep = ":")
        counts <- table(key)
        return(list(values = as.numeric(counts[key]), label = "Density", continuous = TRUE))
    }
    fixed <- list(
        cluster = c("cluster_id", "Cluster", FALSE),
        pseudotime = c("pseudotime", "Pseudotime", TRUE),
        identity = c("predicted_identity", "Predicted identity", FALSE)
    )
    if (!is.null(fixed[[color_by]])) {
        column <- fixed[[color_by]][[1]]
        if (!column %in% names(events)) stop("This result has no ", fixed[[color_by]][[2]], " values.", call. = FALSE)
        return(list(
            values = events[[column]],
            label = fixed[[color_by]][[2]],
            continuous = identical(fixed[[color_by]][[3]], "TRUE")
        ))
    }
    marker_columns <- analysis$metadata$marker_columns %||% list()
    markers <- vapply(marker_columns, function(item) as.character(item$marker %||% "")[[1]], character(1))
    columns <- vapply(marker_columns, function(item) as.character(item$column %||% "")[[1]], character(1))
    marker_index <- match(color_by, c(markers, columns))
    if (is.na(marker_index)) stop("Unknown analysis color value: ", color_by, call. = FALSE)
    if (marker_index > length(markers)) marker_index <- marker_index - length(markers)
    column <- columns[[marker_index]]
    list(values = events[[column]], label = markers[[marker_index]], continuous = TRUE)
}

.analysis_discrete_palette <- function(values) {
    levels <- unique(as.character(values))
    colors <- c(
        "#0c7c86", "#d39a1e", "#cf5148", "#569a43", "#8b55b1", "#ce7131",
        "#416fae", "#bd4f83", "#238b63", "#8c6d31", "#6b6ecf", "#e6550d"
    )
    stats::setNames(rep(colors, length.out = length(levels)), levels)
}

.analysis_discrete_colors <- function(values) {
    unname(.analysis_discrete_palette(values)[as.character(values)])
}

.analysis_plotly_palette <- function(colors) {
    lapply(seq_along(colors), function(index) {
        list((index - 1) / max(1, length(colors) - 1), colors[[index]])
    })
}

.analysis_export_table <- function(analysis) {
    table <- analysis$events
    marker_columns <- analysis$metadata$marker_columns %||% list()
    for (marker in marker_columns) {
        column <- as.character(marker$column %||% "")[[1]]
        label <- as.character(marker$marker %||% "")[[1]]
        if (nzchar(column) && nzchar(label) && column %in% names(table) && !label %in% names(table)) {
            names(table)[names(table) == column] <- label
        }
    }
    coordinate_columns <- grep("^dimension_[0-9]+$", names(table), value = TRUE)
    coordinate_labels <- unlist(analysis$metadata$coordinate_labels %||% character(), use.names = FALSE)
    for (index in seq_along(coordinate_columns)) {
        if (index <= length(coordinate_labels) && nzchar(coordinate_labels[[index]])) {
            names(table)[names(table) == coordinate_columns[[index]]] <- coordinate_labels[[index]]
        }
    }
    table
}
