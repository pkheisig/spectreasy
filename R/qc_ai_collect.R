.ai_qc_read_metric_registry <- function() {
    path <- system.file("extdata", "qc_metric_references.yml", package = "spectreasy")
    if (!nzchar(path)) {
        source_path <- file.path("inst", "extdata", "qc_metric_references.yml")
        if (file.exists(source_path)) path <- source_path
    }
    registry <- if (nzchar(path) && file.exists(path)) yaml::read_yaml(path) else list(version = .ai_qc_metric_version, metrics = list())
    registry
}

.ai_qc_matrix_metrics <- function(M = NULL) {
    if (is.null(M)) {
        return(list(metrics = list(.ai_qc_metric("QC-MATRIX-AVAILABILITY", "reference_matrix", "reference_matrix", missing_reason = "No compatible reference matrix was supplied.")), pairwise_similarity = list(), singular_values = numeric()))
    }
    M <- .as_reference_matrix(M, "M")
    marker_rows <- !grepl("^AF($|_)", rownames(M), ignore.case = TRUE)
    marker_M <- M[marker_rows, , drop = FALSE]
    af_M <- M[!marker_rows, , drop = FALSE]
    sv <- svd(marker_M, nu = 0, nv = 0)$d
    tolerance <- max(dim(marker_M)) * max(sv %||% 0) * .Machine$double.eps
    rank <- sum(sv > tolerance)
    condition <- if (length(sv) && min(sv) > 0) max(sv) / min(sv) else Inf
    rank_metric <- .ai_qc_metric("QC-MATRIX-RANK", "reference_matrix", "reference_matrix", rank, "rank", direction = "higher_better")
    rank_metric <- .ai_qc_hard_grade(rank_metric, failure = rank < min(dim(marker_M)), explanation = if (rank < min(dim(marker_M))) "The marker reference matrix is rank deficient." else "The marker reference matrix has full numerical rank.")
    similarity <- if (nrow(marker_M) > 1L) calculate_similarity_matrix(marker_M) else NULL
    pair_rows <- list()
    if (!is.null(similarity)) {
        idx <- which(upper.tri(similarity), arr.ind = TRUE)
        pair_rows <- lapply(seq_len(nrow(idx)), function(i) list(
            marker_a = rownames(similarity)[idx[i, 1]], marker_b = colnames(similarity)[idx[i, 2]],
            cosine_similarity = as.numeric(similarity[idx[i, 1], idx[i, 2]]),
            interpretation = "panel_complexity_evidence_not_automatic_control_defect"
        ))
        pair_rows <- pair_rows[order(vapply(pair_rows, `[[`, numeric(1), "cosine_similarity"), decreasing = TRUE)]
    }
    decoder <- tryCatch({
        decomp <- svd(marker_M)
        keep <- decomp$d > max(dim(marker_M)) * max(decomp$d) * .Machine$double.eps
        if (!any(keep)) NULL else decomp$v[, keep, drop = FALSE] %*%
            (t(decomp$u[, keep, drop = FALSE]) / decomp$d[keep])
    }, error = function(e) NULL)
    decoder_norms <- if (!is.null(decoder)) sqrt(colSums(decoder^2)) else numeric()
    detector_leverage <- if (nrow(marker_M)) colSums(marker_M^2) else numeric()
    af_similarity <- if (nrow(af_M) && nrow(marker_M)) {
        all_sim <- calculate_similarity_matrix(rbind(marker_M, af_M))
        all_sim[rownames(marker_M), rownames(af_M), drop = FALSE]
    } else NULL
    list(
        dimensions = list(rows = nrow(M), markers = nrow(marker_M), af_bands = nrow(af_M), detectors = ncol(M)),
        singular_values = as.numeric(sv), effective_rank = rank,
        metrics = list(
            rank = rank_metric,
            condition_number = .ai_qc_metric("QC-MATRIX-CONDITION", "reference_matrix", "reference_matrix", condition, "ratio", direction = "lower_better", metadata = list(panel_complexity_evidence = TRUE)),
            maximum_similarity = .ai_qc_metric("QC-MATRIX-MAX-SIMILARITY", "reference_matrix", "reference_matrix", if (!is.null(similarity)) max(similarity[upper.tri(similarity)], na.rm = TRUE) else NULL, "cosine_similarity", direction = "lower_better", missing_reason = "At least two marker spectra are required.", metadata = list(panel_complexity_evidence = TRUE))
        ),
        pairwise_similarity = pair_rows,
        decoder_norms = as.list(decoder_norms),
        detector_leverage = as.list(detector_leverage),
        af_similarity = if (is.null(af_similarity)) NULL else unname(split(as.data.frame(af_similarity), rownames(af_similarity))),
        adjustment = list(adjusted = isTRUE(attr(M, "adjusted")), provenance = attr(M, "adjustment_provenance") %||% NULL)
    )
}

.ai_qc_detector_schema_hash <- function(M) {
    if (is.null(M)) return(NULL)
    raw <- charToRaw(paste(colnames(M), collapse = "\n"))
    paste0(openssl::sha256(raw))
}

.ai_qc_af_metrics <- function(M = NULL, info = list()) {
    if (is.null(M)) {
        return(list(
            available = FALSE, source = NULL, source_count = NULL,
            requested_bands = NULL, derived_bands = 0L, sources = NULL,
            selection = NULL, bank_similarity = NULL, effective_rank = NULL,
            representative_band = NULL, extreme_band = NULL,
            marker_interactions = list(), metrics = list(.ai_qc_metric(
                "QC-AF-AVAILABILITY", "af_bank", "af", missing_reason = "No reference matrix was supplied."
            ))
        ))
    }
    af_rows <- grepl("^AF($|_)", rownames(M), ignore.case = TRUE)
    af_matrix <- M[af_rows, , drop = FALSE]
    marker_matrix <- M[!af_rows, , drop = FALSE]
    if (!nrow(af_matrix)) {
        return(list(
            available = FALSE, source = info$mode %||% NULL,
            source_count = info$source_count %||% NULL,
            requested_bands = info$requested_bands %||% NULL,
            derived_bands = 0L, sources = info$sources %||% NULL,
            selection = info$selection %||% NULL, bank_similarity = NULL,
            effective_rank = 0L, representative_band = NULL,
            extreme_band = NULL, marker_interactions = list(),
            metrics = list(.ai_qc_metric(
                "QC-AF-AVAILABILITY", "af_bank", "af", missing_reason = "The reference matrix contains no AF bands."
            ))
        ))
    }
    norms <- sqrt(rowSums(af_matrix^2))
    normalized <- af_matrix / pmax(norms, .Machine$double.eps)
    centroid <- colMeans(normalized)
    distance <- sqrt(rowSums((normalized - matrix(centroid, nrow(normalized), ncol(normalized), byrow = TRUE))^2))
    representative <- rownames(af_matrix)[which.min(distance)]
    extreme <- rownames(af_matrix)[which.max(distance)]
    similarity_values <- numeric()
    if (nrow(af_matrix) > 1L) {
        similarity <- calculate_similarity_matrix(af_matrix)
        similarity_values <- similarity[upper.tri(similarity)]
    }
    singular_values <- svd(af_matrix, nu = 0, nv = 0)$d
    tolerance <- max(dim(af_matrix)) * max(singular_values) * .Machine$double.eps
    effective_rank <- sum(singular_values > tolerance)
    interactions <- list()
    if (nrow(marker_matrix)) {
        combined <- calculate_similarity_matrix(rbind(marker_matrix, af_matrix))
        cross <- combined[rownames(marker_matrix), rownames(af_matrix), drop = FALSE]
        index <- which(is.finite(cross), arr.ind = TRUE)
        interactions <- lapply(seq_len(nrow(index)), function(i) list(
            marker = rownames(cross)[index[i, 1]],
            af_band = colnames(cross)[index[i, 2]],
            cosine_similarity = as.numeric(cross[index[i, 1], index[i, 2]])
        ))
        interactions <- interactions[order(vapply(interactions, `[[`, numeric(1), "cosine_similarity"), decreasing = TRUE)]
    }
    list(
        available = TRUE, source = info$mode %||% NULL,
        source_count = info$source_count %||% NULL,
        requested_bands = info$requested_bands %||% NULL,
        derived_bands = info$derived_bands %||% nrow(af_matrix),
        sources = info$sources %||% NULL, selection = info$selection %||% NULL,
        bank_similarity = if (length(similarity_values)) list(
            minimum = min(similarity_values), median = stats::median(similarity_values),
            maximum = max(similarity_values)
        ) else NULL,
        singular_values = as.numeric(singular_values), effective_rank = effective_rank,
        representative_band = representative, extreme_band = extreme,
        marker_interactions = interactions,
        metrics = list(
            .ai_qc_metric("QC-AF-EFFECTIVE-RANK", "af_bank", "af", effective_rank, "rank", direction = "higher_better"),
            .ai_qc_metric(
                "QC-AF-MAX-SIMILARITY", "af_bank", "af",
                if (length(similarity_values)) max(similarity_values) else NULL,
                "cosine_similarity", direction = "lower_better",
                missing_reason = "At least two AF bands are required for bank similarity."
            )
        )
    )
}

.ai_qc_discover_existing <- function(existing_output_dir) {
    empty <- list(M = NULL, control_report_data = NULL, sample_report_data = NULL, source_paths = character())
    if (is.null(existing_output_dir) || !dir.exists(existing_output_dir)) return(empty)
    find_latest <- function(pattern) {
        files <- list.files(existing_output_dir, pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
        if (!length(files)) return(NULL)
        files[order(file.info(files)$mtime, decreasing = TRUE)][1]
    }
    matrix_path <- find_latest("^scc_reference_matrix\\.csv$")
    M <- if (!is.null(matrix_path)) tryCatch(.read_unmixing_matrix_csv(matrix_path), error = function(e) NULL) else NULL
    noise_path <- if (!is.null(matrix_path)) file.path(dirname(matrix_path), "scc_detector_noise.csv") else NULL
    if (!is.null(M) && !is.null(noise_path) && file.exists(noise_path)) {
        M <- tryCatch(.attach_detector_noise(M, utils::read.csv(noise_path, check.names = FALSE), source = noise_path), error = function(e) M)
    }
    read_metric <- function(pattern) {
        path <- find_latest(pattern)
        if (is.null(path)) return(NULL)
        tryCatch(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    }
    detector_rms <- read_metric("^rms_residual_per_detector\\.csv$")
    sample_report <- list(
        matrix = M,
        nps = read_metric("^negative_population_spread\\.csv$"),
        detector_rms = detector_rms,
        reconstruction_error = read_metric("^detector_reconstruction_error\\.csv$"),
        unmixing_method = if (!is.null(detector_rms) && "unmixing_method" %in% names(detector_rms)) as.character(detector_rms$unmixing_method[1]) else NULL
    )
    class(sample_report) <- c("spectreasy_sample_report_data", "spectreasy_report_data", "list")
    paths <- unique(Filter(Negate(is.null), c(matrix_path, noise_path,
        find_latest("^negative_population_spread\\.csv$"),
        find_latest("^rms_residual_per_detector\\.csv$"),
        find_latest("^detector_reconstruction_error\\.csv$")
    )))
    list(M = M, control_report_data = NULL, sample_report_data = sample_report, source_paths = unlist(paths, use.names = FALSE))
}

#' Collect a canonical AI-ready QC object
#'
#' Builds a local, deterministic, schema-versioned QC object from in-memory
#' workflow results and compatible machine-readable artifacts. It does not call
#' an AI service and never includes raw event matrices.
#'
#' @param controls Optional result from [unmix_controls()].
#' @param samples Optional result from [unmix_samples()].
#' @param M Optional reference matrix.
#' @param control_report_data Optional control report data object.
#' @param sample_report_data Optional sample report data object.
#' @param project_dir Optional Spectreasy project directory.
#' @param existing_output_dir Optional existing output directory used for
#'   machine-readable artifact discovery.
#' @param scope One of `"auto"`, `"control"`, `"sample"`, or `"combined"`.
#' @param privacy Privacy mode: `"standard"`, `"strict"`, or `"none"`.
#' @param reference QC reference profile selection.
#' @param generated_at Generation time. Tests may pass a fixed time.
#' @param context Optional experiment metadata.
#' @return A `spectreasy_ai_qc` object.
#' @export
collect_ai_qc <- function(
    controls = NULL, samples = NULL, M = NULL,
    control_report_data = NULL, sample_report_data = NULL,
    project_dir = NULL, existing_output_dir = NULL,
    scope = c("auto", "control", "sample", "combined"),
    privacy = c("standard", "strict", "none"), reference = "auto",
    generated_at = Sys.time(), context = list()
) {
    scope <- .match_arg_ci(scope, c("auto", "control", "sample", "combined"), "scope")
    privacy <- .match_arg_ci(privacy, c("standard", "strict", "none"), "privacy")
    discovered <- .ai_qc_discover_existing(existing_output_dir)
    if (is.null(M)) M <- controls$M %||% control_report_data$matrix %||% sample_report_data$matrix %||% attr(samples, "reference_matrix") %||% discovered$M %||% NULL
    if (is.null(control_report_data)) control_report_data <- discovered$control_report_data
    if (is.null(sample_report_data) && !is.null(discovered$M)) sample_report_data <- discovered$sample_report_data
    if (identical(scope, "auto")) scope <- if (!is.null(controls) || !is.null(control_report_data)) {
        if (!is.null(samples) || !is.null(sample_report_data)) "combined" else "control"
    } else "sample"
    method <- sample_report_data$unmixing_method %||% attr(samples, "method") %||% control_report_data$unmixing_method %||% context$method %||% "AutoSpectral"
    matrix_part <- .ai_qc_matrix_metrics(M)
    control_part <- .ai_qc_control_metrics(M, control_report_data, controls)
    sample_part <- .ai_qc_sample_metrics(samples, M, sample_report_data, method)
    pairwise <- .ai_qc_pairwise_metrics(samples)
    af_info <- if (!is.null(M)) attr(M, "af_bank_info") %||% list() else list()
    af_rows <- if (!is.null(M)) grepl("^AF($|_)", rownames(M), ignore.case = TRUE) else logical()
    profile_context <- list(
        cytometer = context$cytometer %||% control_report_data$cytometer %||% NULL,
        detector_schema_hash = .ai_qc_detector_schema_hash(M),
        method = method,
        panel_size = if (!is.null(M)) sum(!af_rows) else NULL,
        af_configuration = if (any(af_rows)) paste0(sum(af_rows), "_bands") else "none"
    )
    profile <- load_qc_reference_profile(reference, project_dir = project_dir, context = profile_context)
    registry <- .ai_qc_read_metric_registry()
    source_paths <- unique(Filter(function(x) is.character(x) && length(x) == 1L && nzchar(x), c(
        project_dir, existing_output_dir, control_report_data$matrix_source,
        sample_report_data$matrix_source, controls$reference_matrix_file, discovered$source_paths
    )))
    structural_metrics <- list(
        schema = .ai_qc_hard_grade(.ai_qc_metric("QC-STRUCT-SCHEMA", "export", "structural_validation", TRUE, "logical"), explanation = "The canonical object matches schema version 1.0.0."),
        matrix = matrix_part$metrics$rank
    )
    components <- list(
        schema = list(name = .ai_qc_schema_name, version = .ai_qc_schema_version),
        provenance = list(
            generated_at = format(as.POSIXct(generated_at, tz = "UTC"), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC"),
            package_version = as.character(utils::packageVersion("spectreasy")),
            git_commit = .ai_qc_git_commit(),
            metric_implementation_version = .ai_qc_metric_version,
            source_artifacts = .ai_qc_source_records(source_paths)
        ),
        privacy = list(mode = privacy, local_only = TRUE, raw_events_included = FALSE, notice = "Generated locally by Spectreasy. Spectreasy sends nothing and contains no embedded AI model, provider integration, API key, or upload step."),
        scope = list(value = scope, stages = intersect(c("controls", "samples"), c(if (scope %in% c("control", "combined")) "controls", if (scope %in% c("sample", "combined")) "samples"))),
        experiment = c(profile_context, context),
        quality_reference = list(
            mode = if (is.null(profile)) "none" else "profile",
            profile = if (is.null(profile)) NULL else list(name = profile$name, version = profile$version, cohort = profile$cohort, compatibility = .ai_qc_profile_compatibility(profile, profile_context))
        ),
        overall_summary = list(status = "partial", top_findings = list()),
        structural_validation = list(metrics = structural_metrics),
        controls = control_part,
        reference_matrix = matrix_part,
        af = .ai_qc_af_metrics(M, af_info),
        samples = sample_part,
        detectors = list(metrics = sample_part$detector_metrics, noise = if (!is.null(M)) attr(M, "detector_noise") else NULL),
        pairwise_findings = list(within_stage = pairwise, cross_stage = .ai_qc_cross_stage_findings(control_part, sample_part, pairwise)),
        grade_summary = list(),
        limitations = c(
            "AI interpretation is advisory and requires scientific review.",
            "Historical outputs lacking machine-readable metadata remain partial and are not inferred from plots, HTML, or PDF.",
            "Within-run associations are descriptive and do not establish causation.",
            if (is.null(profile)) "No reviewed compatible reference cohort was available; most continuous metrics are not graded." else NULL
        ),
        metric_definitions = registry$metrics %||% list(),
        references = registry$references %||% list()
    )
    components <- .ai_qc_apply_profile(components, profile)
    components$grade_summary <- .ai_qc_grade_summary(components)
    poor <- components$grade_summary$counts$poor %||% 0L
    review <- components$grade_summary$counts$review %||% 0L
    components$overall_summary$status <- if (poor > 0L) "poor" else if (review > 0L) "review" else if ((components$grade_summary$counts$not_graded %||% 0L) > 0L) "partial" else "good"
    components$overall_summary$top_findings <- .ai_qc_top_findings(components)
    normalized <- .ai_qc_normalize_nonfinite(components)
    components <- normalized$value
    components$provenance$nonfinite_normalization <- normalized$nonfinite
    object <- .new_spectreasy_ai_qc(components)
    .redact_ai_qc(object, mode = privacy)
}

.ai_qc_top_findings <- function(x, n = 8L) {
    metrics <- .ai_qc_walk_metrics(x)
    priority <- c(poor = 1L, review = 2L, good = 3L, not_graded = 4L)
    metrics <- metrics[order(vapply(metrics, function(m) priority[[m$grade]], integer(1)), vapply(metrics, `[[`, character(1), "metric_id"))]
    lapply(head(metrics, n), function(m) list(metric_id = m$metric_id, entity = m$entity, stage = m$stage, grade = m$grade, explanation = m$grade_provenance$explanation))
}

.ai_qc_git_commit <- function() {
    out <- tryCatch(suppressWarnings(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)), error = function(e) character())
    if (length(out) && grepl("^[0-9a-f]{40}$", out[1])) out[1] else NULL
}

.ai_qc_source_records <- function(paths) {
    lapply(sort(unique(paths)), function(path) list(
        path = normalizePath(path, mustWork = FALSE),
        exists = file.exists(path),
        sha256 = if (file.exists(path) && !dir.exists(path)) .ai_qc_file_hash(path) else NULL
    ))
}
