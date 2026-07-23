population_code_fixture <- function() {
    project <- tempfile("population_code_api_")
    dir.create(file.path(project, "samples"), recursive = TRUE)
    set.seed(20260724)
    group <- rep(seq_len(3L), each = 40L)
    values <- cbind(
        `FSC-A` = stats::rnorm(120L, c(40, 80, 120)[group], 3),
        `SSC-A` = stats::rnorm(120L, c(25, 60, 95)[group], 3),
        `CD3-A` = stats::rnorm(120L, c(1, 7, 12)[group], 0.4),
        `CD19-A` = stats::rnorm(120L, c(12, 5, 1)[group], 0.4)
    )
    frame <- flowCore::flowFrame(values)
    parameters <- flowCore::parameters(frame)@data
    parameters$desc <- c("FSC", "SSC", "CD3", "CD19")
    flowCore::parameters(frame)@data <- parameters
    path <- file.path(project, "samples", "sample.fcs")
    flowCore::write.FCS(frame, path)
    list(
        project = project,
        file = "samples/sample.fcs",
        markers = colnames(values)
    )
}

test_that("exported gating API creates, updates, summarizes, and deletes gates", {
    fixture <- population_code_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    initial <- analysis_workspace(fixture$project)
    expect_identical(initial$populations[[1]]$id, "root")
    saved <- add_population_gate(
        fixture$project,
        name = "Low FSC",
        type = "range",
        x = "FSC-A",
        geometry = list(min = 20, max = 60),
        role = "negative",
        id = "low-fsc"
    )
    expect_true("low-fsc" %in% vapply(saved$populations, `[[`, character(1), "id"))
    updated <- update_population_gate(fixture$project, "low-fsc", name = "FSC low")
    expect_identical(updated$populations[[2]]$name, "FSC low")
    statistics <- population_statistics(fixture$project, fixture$file, "low-fsc", c("FSC-A", "CD3-A"))
    expect_gt(statistics$count, 0L)
    expect_lt(statistics$count, statistics$total_count)
    expect_equal(nrow(statistics$markers), 2L)
    deleted <- delete_population_gate(fixture$project, "low-fsc")
    expect_identical(vapply(deleted$populations, `[[`, character(1), "id"), "root")
})

test_that("the exported method registry exposes every executable parameter schema", {
    methods <- analysis_methods(refresh = TRUE)
    ids <- vapply(methods, `[[`, character(1), "id")
    expect_setequal(ids, c(
        "pca", "tsne", "umap", "diffusion-map", "phate", "hsne",
        "flowsom", "phenograph", "dpt", "slingshot", "tscan",
        "palantir", "paga-dpt", "wanderlust", "wishbone"
    ))
    by_id <- stats::setNames(methods, ids)
    expect_setequal(
        vapply(by_id$umap$parameters, `[[`, character(1), "id"),
        c("neighbors", "min_dist", "spread", "metric", "epochs", "learning_rate", "repulsion", "negative_samples")
    )
    expect_setequal(
        vapply(by_id$wishbone$parameters, `[[`, character(1), "id"),
        c("neighbors", "candidate_neighbors", "graphs", "waypoints", "iterations", "diffusion_components", "branch_confidence")
    )
    expect_true(all(vapply(methods, function(method) {
        all(vapply(method$parameters, function(parameter) {
            all(c("id", "label", "type", "default") %in% names(parameter))
        }, logical(1)))
    }, logical(1))))
})

test_that("the full analysis installer has a complete platform-neutral dependency plan", {
    plan <- spectreasy:::.analysis_dependency_plan()
    expect_setequal(plan$cran, c("Rtsne", "uwot"))
    expect_setequal(plan$bioconductor, c(
        "DelayedMatrixStats", "destiny", "FlowSOM", "slingshot", "TSCAN"
    ))
    expect_identical(
        unname(plan$github[["Rphenograph"]]),
        "JinmiaoChenLab/Rphenograph@0298487f0ee13aac55eb77d19992f6bd878ba2fc"
    )
    expect_true(is.function(install_analysis_dependencies))
})

test_that("exported R pipelines use advanced settings and reuse fitted objects", {
    fixture <- population_code_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    common <- list(
        project_path = fixture$project,
        file = fixture$file,
        markers = fixture$markers,
        clustering = "flowsom",
        max_events = 100L,
        seed = 41L,
        cluster_settings = list(
            clusters = 3L, xdim = 4L, ydim = 4L, rlen = 5L,
            alpha_start = 0.05, alpha_end = 0.01
        )
    )
    pca <- do.call(analyze_population, c(common, list(
        reduction = "pca",
        reduction_settings = list(center = TRUE, scale = TRUE)
    )))
    diffusion <- do.call(analyze_population, c(common, list(
        reduction = "diffusion-map",
        reduction_settings = list(
            neighbors = 8L, eigenvectors = 5L,
            distance = "euclidean", density_normalization = TRUE
        )
    )))

    expect_s3_class(pca, "spectreasy_analysis")
    expect_identical(pca$metadata$display_name, "FlowSOM \u2192 PCA")
    expect_identical(pca$metadata$parameters$clustering$clusters, 3L)
    expect_true(pca$metadata$parameters$method$center)
    expect_true(pca$metadata$parameters$method$scale)
    expect_identical(
        pca$metadata$artifacts$clustering$id,
        diffusion$metadata$artifacts$clustering$id
    )
    expect_false(pca$metadata$artifacts$clustering$reused)
    expect_true(diffusion$metadata$artifacts$clustering$reused)
    expect_identical(pca$events$event_id, diffusion$events$event_id)
    expect_identical(pca$events$cluster_id, diffusion$events$cluster_id)

    loaded <- load_population_analysis(fixture$project, pca$metadata$analysis_id)
    expect_s3_class(loaded, "spectreasy_analysis")
    expect_identical(loaded$events$event_id, pca$events$event_id)
})

test_that("exported trajectory API channels prerequisites and method settings", {
    fixture <- population_code_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    result <- analyze_population(
        project_path = fixture$project,
        file = fixture$file,
        markers = fixture$markers,
        trajectory = "dpt",
        root_event_id = 1L,
        max_events = 100L,
        seed = 52L,
        trajectory_settings = list(
            neighbors = 8L,
            eigenvectors = 5L,
            distance = "euclidean",
            density_normalization = TRUE
        )
    )
    expect_identical(result$metadata$method$id, "dpt")
    expect_identical(result$metadata$root_event_id, 1L)
    expect_identical(result$metadata$parameters$method$neighbors, 8L)
    expect_true(all(is.finite(result$events$pseudotime)))
    expect_equal(range(result$events$pseudotime), c(0, 1))

    expect_error(
        analyze_population(
            fixture$project, fixture$file, fixture$markers,
            trajectory = "dpt", max_events = 100L
        ),
        "requires a root event"
    )
})

test_that("common immune identity templates are panel aware and editable", {
    signatures <- population_identity_templates(c(
        "CD3-A", "CD4-A", "CD8-A", "CD19-A", "CD56-A", "CD14-A"
    ))
    expect_identical(
        vapply(signatures, `[[`, character(1), "name"),
        c("CD4 T cell", "CD8 T cell", "B cell", "NK cell", "Monocyte")
    )
    expect_identical(signatures[[1]]$positive_markers, c("CD3-A", "CD4-A"))
    reduced <- population_identity_templates(c("CD3", "CD8", "CD56"))
    expect_identical(vapply(reduced, `[[`, character(1), "name"), c("CD8 T cell", "NK cell"))
    expect_error(population_identity_templates(character()), "at least one")
})

test_that("R annotation, plotting, palettes, and table exports are first-class", {
    fixture <- population_code_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    result <- analyze_population(
        fixture$project, fixture$file, fixture$markers,
        reduction = "pca", max_events = 100L, seed = 63L,
        reduction_settings = list(center = TRUE, scale = TRUE)
    )
    annotation <- annotate_population(
        result,
        project_path = fixture$project,
        signatures = list(
            list(name = "T-like", positive_markers = "CD3", negative_markers = "CD19"),
            list(name = "B-like", positive_markers = "CD19", negative_markers = "CD3")
        ),
        min_score = 0.5,
        min_margin = 0.05,
        evidence_sensitivity = 1.2
    )
    expect_identical(annotation$metadata$method, "robust-signed-marker-score")
    expect_equal(annotation$metadata$event_count, 100L)

    loaded <- load_population_analysis(
        fixture$project,
        result$metadata$analysis_id,
        annotation$metadata$identity_id
    )
    expect_true(all(c("predicted_identity", "identity_score", "identity_margin") %in% names(loaded$events)))

    marker_plot <- plot_population_analysis(
        loaded, x = "PC 1", y = "PC 2",
        color_by = "CD3", palette = "sunset", point_size = 0.5
    )
    identity_plot <- plot_population_analysis(loaded, color_by = "identity")
    expect_s3_class(marker_plot, "ggplot")
    expect_s3_class(identity_plot, "ggplot")
    expect_setequal(
        names(population_analysis_palettes()),
        c("control-density", "viridis", "sunset", "plasma", "inferno", "magma", "cividis", "turbo", "ice-fire", "spectral")
    )

    csv <- file.path(fixture$project, "exports", "analysis.csv")
    svg <- file.path(fixture$project, "exports", "analysis.svg")
    rds <- file.path(fixture$project, "exports", "analysis.rds")
    export_population_analysis(loaded, csv)
    export_population_analysis(loaded, svg, color_by = "CD19", palette = "viridis")
    export_population_analysis(loaded, rds)
    expect_true(all(file.exists(c(csv, svg, rds))))
    table <- data.table::fread(csv)
    expect_true(all(c(
        "event_id", "PC 1", "PC 2", "CD3", "CD19",
        "predicted_identity", "identity_score", "identity_margin"
    ) %in% names(table)))

    if (requireNamespace("plotly", quietly = TRUE)) {
        plot_3d <- plot_population_analysis(loaded, z = 3L, color_by = "CD3", palette = "viridis")
        expect_s3_class(plot_3d, "plotly")
    }
})

test_that("public settings fail before execution when invalid", {
    fixture <- population_code_fixture()
    on.exit(unlink(fixture$project, recursive = TRUE, force = TRUE), add = TRUE)
    expect_error(
        analyze_population(
            fixture$project, fixture$file, fixture$markers,
            clustering = "flowsom", reduction = "umap", max_events = 100L,
            cluster_settings = list(alpha_start = 0.01, alpha_end = 0.05)
        ),
        "final learning rate"
    )
    wishbone <- analysis_methods()
    wishbone <- wishbone[[match(
        "wishbone",
        vapply(wishbone, `[[`, character(1), "id")
    )]]
    wishbone_error <- if (isTRUE(wishbone$available)) {
        "candidate neighbors"
    } else {
        "not executable"
    }
    expect_error(
        analyze_population(
            fixture$project, fixture$file, fixture$markers,
            trajectory = "wishbone", root_event_id = 1L, max_events = 100L,
            trajectory_settings = list(neighbors = 20L, candidate_neighbors = 15L)
        ),
        wishbone_error
    )
    expect_error(
        analyze_population(
            fixture$project, fixture$file, fixture$markers,
            reduction = "umap", max_events = 100L,
            reduction_settings = list(not_a_setting = 1)
        ),
        "no advanced setting"
    )
})
