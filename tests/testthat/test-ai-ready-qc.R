make_ai_qc_fixture <- function() {
    M <- rbind(FITC = c(1, 0.2, 0.05), PE = c(0.1, 1, 0.2), AF_1 = c(0.2, 0.15, 1))
    colnames(M) <- c("B1-A", "YG1-A", "R1-A")
    attr(M, "qc_summary") <- data.frame(
        sample = c("FITC_control.fcs", "PE_control.fcs"),
        fluorophore = c("FITC", "PE"), marker = c("CD3", "CD4"), type = "beads",
        peak_channel = c("B1-A", "YG1-A"), n_total = c(10000, 10000),
        n_scatter_gated = c(8000, 7900), n_final = c(200, 200),
        stain_index = c(12.5, 9.1), saturated = c("OK", "YES"),
        actual_spectral_events = c(200, 200), stringsAsFactors = FALSE
    )
    attr(M, "af_bank_info") <- list(mode = "pooled_af_sources", source_count = 1L, requested_bands = 1L, derived_bands = 1L)
    attr(M, "detector_noise") <- data.frame(detector = colnames(M), noise_floor = c(20, 22, 19))
    set.seed(7)
    samples <- list(
        patient_alpha = list(
            data = data.frame(File = "patient_alpha", FITC = rnorm(80), PE = rnorm(80), `AF Index` = sample(1:2, 80, TRUE), `FSC-A` = rnorm(80), check.names = FALSE),
            residuals = matrix(rnorm(240), ncol = 3, dimnames = list(NULL, colnames(M)))
        )
    )
    class(samples) <- c("spectreasy_unmixed_results", "list")
    attr(samples, "method") <- "WLS"
    attr(samples, "reference_matrix") <- M
    list(M = M, samples = samples)
}

test_that("canonical AI-QC schema has exact components and typed grade provenance", {
    fixture <- make_ai_qc_fixture()
    qc <- collect_ai_qc(
        samples = fixture$samples, controls = list(M = fixture$M), M = fixture$M,
        scope = "combined", reference = "none", privacy = "standard",
        generated_at = as.POSIXct("2026-01-02 03:04:05", tz = "UTC")
    )
    expect_s3_class(qc, "spectreasy_ai_qc")
    expect_identical(names(qc), spectreasy:::.ai_qc_top_level)
    expect_true(validate_ai_qc(qc))
    expect_true(qc$af$available)
    expect_identical(qc$af$effective_rank, 1L)
    expect_true(length(qc$af$marker_interactions) > 0L)
    expect_true(length(qc$samples$negative_bias) > 0L)
    metrics <- spectreasy:::.ai_qc_walk_metrics(qc)
    expect_true(all(vapply(metrics, function(metric) metric$grade %in% c("good", "review", "poor", "not_graded"), logical(1))))
    expect_true(all(vapply(metrics, function(metric) all(c("basis", "threshold", "interval", "profile_name", "profile_version", "reference_n", "matching_strata", "direction", "explanation") %in% names(metric$grade_provenance)), logical(1))))
    expect_false(any(vapply(qc, function(value) is.numeric(value) && any(is.infinite(value)), logical(1))))
})

test_that("AI-QC serialization and hashes are stable with a fixed timestamp", {
    fixture <- make_ai_qc_fixture()
    qc <- collect_ai_qc(samples = fixture$samples, controls = list(M = fixture$M), M = fixture$M, scope = "combined", reference = "none", generated_at = as.POSIXct("2026-01-02", tz = "UTC"))
    first <- export_ai_qc(qc, output_dir = tempfile("ai_qc_a_"), overwrite = "overwrite")
    second <- export_ai_qc(qc, output_dir = tempfile("ai_qc_b_"), overwrite = "overwrite")
    expect_identical(readBin(first$paths[["json"]], "raw", n = file.info(first$paths[["json"]])$size), readBin(second$paths[["json"]], "raw", n = file.info(second$paths[["json"]])$size))
    expect_identical(first$content_sha256, second$content_sha256)
    third <- export_ai_qc(first$object, output_dir = tempfile("ai_qc_c_"), overwrite = "overwrite")
    expect_identical(first$content_sha256, third$content_sha256)
    expect_true(all(file.exists(first$paths)))
    json <- paste(readLines(first$paths[["json"]], warn = FALSE), collapse = "\n")
    expect_false(grepl("NaN|Infinity|-Infinity", json))
})

test_that("control, sample, and combined bundles use exact artifact names", {
    fixture <- make_ai_qc_fixture()
    expected <- c(
        control = "spectreasy_ai_qc_controls",
        sample = "spectreasy_ai_qc_samples",
        combined = "spectreasy_ai_qc_combined"
    )
    for (scope in names(expected)) {
        object <- collect_ai_qc(
            controls = if (scope %in% c("control", "combined")) list(M = fixture$M) else NULL,
            samples = if (scope %in% c("sample", "combined")) fixture$samples else NULL,
            M = fixture$M, scope = scope, reference = "none",
            generated_at = as.POSIXct("2026-01-02", tz = "UTC")
        )
        bundle <- export_ai_qc(object, output_dir = tempfile(paste0("ai_qc_", scope, "_")), overwrite = "overwrite")
        expect_setequal(
            basename(bundle$paths),
            paste0(expected[[scope]], c(".json", ".txt", ".md", "_prompt.txt", "_manifest.json"))
        )
        manifest <- jsonlite::fromJSON(bundle$paths[["manifest"]], simplifyVector = FALSE)
        expect_identical(manifest$scope, scope)
        expect_true(all(vapply(manifest$artifacts, function(item) !startsWith(item$path, "/"), logical(1))))
    }
})

test_that("privacy redaction aliases samples and strict mode aliases control files", {
    fixture <- make_ai_qc_fixture()
    standard <- collect_ai_qc(samples = fixture$samples, controls = list(M = fixture$M), M = fixture$M, scope = "combined", privacy = "standard", reference = "none")
    strict <- collect_ai_qc(samples = fixture$samples, controls = list(M = fixture$M), M = fixture$M, scope = "combined", privacy = "strict", reference = "none")
    expect_identical(standard$samples$entities[[1]]$id, "sample_001")
    expect_identical(strict$controls$entities[[1]]$sample, "control_001")
    expect_false(grepl("patient_alpha", spectreasy:::.ai_qc_to_json(standard), fixed = TRUE))
})

test_that("AF Index and acquisition fields are excluded from marker diagnostics", {
    columns <- c("File", "FITC", "PE", "AF Index", "FSC-A", "SSC-A", "Time", "AF_1")
    excluded <- spectreasy:::.get_result_metadata_columns(columns)
    expect_true(all(c("File", "AF Index", "FSC-A", "SSC-A", "Time", "AF_1") %in% excluded))
    df <- data.frame(File = rep("x", 20), FITC = rnorm(20), PE = rnorm(20), `AF Index` = 1, `FSC-A` = rnorm(20), check.names = FALSE)
    expect_setequal(unique(calculate_nps(df)$Marker), c("FITC", "PE"))
})

test_that("fixed SCC selection is not penalized for percentage retention", {
    fixture <- make_ai_qc_fixture()
    qc <- collect_ai_qc(controls = list(M = fixture$M), M = fixture$M, scope = "control", reference = "none")
    final <- qc$controls$entities[[1]]$metrics$final_events
    retention <- qc$controls$entities[[1]]$metrics$retention
    expect_identical(final$grade, "good")
    expect_true(final$metadata$intentional_fixed_selection)
    expect_identical(retention$grade, "not_graded")
})

test_that("reference profiles require five matched clean runs for empirical grades", {
    metric <- spectreasy:::.ai_qc_metric("QC-DETECTOR-RMS", "B1-A", "samples", 12, "residual_rms", direction = "lower_better")
    profile <- structure(list(name = "test", version = "1", strata = list(method = "WLS"), metrics = list(), cohort = list(n_clean = 4)), class = c("spectreasy_qc_reference_profile", "list"))
    reference <- list(n = 4, median = 2, mad = 1, interval = c(1, 3))
    graded <- spectreasy:::.ai_qc_profile_grade(metric, reference, profile)
    expect_identical(graded$grade, "not_graded")
    expect_match(graded$grade_provenance$explanation, "Fewer than five")
    reference$n <- 6
    graded <- spectreasy:::.ai_qc_profile_grade(metric, reference, profile)
    expect_identical(graded$grade, "poor")
    expect_identical(graded$grade_provenance$basis, "empirical_reference")
})

test_that("grading bases and inclusive threshold boundaries are explicit", {
    metric <- spectreasy:::.ai_qc_metric("QC-TEST", "entity", "samples", 2.5, "ratio", direction = "lower_better")
    literature <- spectreasy:::.ai_qc_threshold_grade(
        metric, review_threshold = 2.5, poor_threshold = 4,
        basis = "literature_rule", assumptions_met = TRUE, event_reliable = TRUE
    )
    expect_identical(literature$grade, "review")
    expect_identical(literature$grade_provenance$basis, "literature_rule")

    metric$value <- 4
    user <- spectreasy:::.ai_qc_threshold_grade(
        metric, review_threshold = 2.5, poor_threshold = 4,
        basis = "user_reference"
    )
    expect_identical(user$grade, "poor")
    expect_identical(user$grade_provenance$basis, "user_reference")

    unreliable <- spectreasy:::.ai_qc_threshold_grade(
        metric, review_threshold = 2.5, poor_threshold = 4,
        basis = "literature_rule", event_reliable = FALSE
    )
    expect_identical(unreliable$grade, "not_graded")
    expect_identical(unreliable$grade_provenance$basis, "not_graded")

    hard <- spectreasy:::.ai_qc_hard_grade(metric, review = TRUE)
    expect_identical(hard$grade_provenance$basis, "hard_validation")
    within <- lapply(c(1, 1.1, 0.9, 1.05, 20), function(value) spectreasy:::.ai_qc_metric("QC-TEST", "entity", "samples", value))
    within <- spectreasy:::.ai_qc_apply_within_run(within)
    expect_identical(within[[5]]$grade_provenance$basis, "relative_within_run")
})

test_that("human renderers have 21 sections and prompt delimiters", {
    fixture <- make_ai_qc_fixture()
    qc <- collect_ai_qc(samples = fixture$samples, M = fixture$M, scope = "sample", reference = "none")
    expect_length(spectreasy:::.ai_qc_sections(qc), 21L)
    prompt <- build_ai_qc_prompt(qc, detail = "compact")
    expect_match(prompt, "<<<SPECTREASY_AI_QC_DATA_BEGIN>>>", fixed = TRUE)
    expect_match(prompt, "<<<SPECTREASY_AI_QC_DATA_END>>>", fixed = TRUE)
    expect_match(prompt, "Spectreasy sends nothing", fixed = TRUE)
    contract <- readLines(testthat::test_path("fixtures", "ai_qc_prompt_contract.txt"), warn = FALSE)
    expect_true(all(contract[nzchar(contract)] %in% strsplit(prompt, "\n", fixed = TRUE)[[1]]))
})

test_that("HTML reports embed the canonical AI-QC summary without external calls", {
    fixture <- make_ai_qc_fixture()
    report <- collect_sample_report_data(fixture$samples, fixture$M, max_events_per_sample = 20, plot_dir = tempfile("plots_"))
    path <- tempfile(fileext = ".html")
    render_qc_html_report(report, path, overwrite = "overwrite")
    html <- paste(readLines(path, warn = FALSE), collapse = "\n")
    expect_match(html, "AI-ready QC", fixed = TRUE)
    expect_match(html, "Copy prompt", fixed = TRUE)
    expect_false(grepl("https://", html, fixed = TRUE))
})
