test_that("with_optional_seed creates deterministic local draws", {
    val1 <- local({
        spectreasy:::.with_optional_seed(999)
        runif(3)
    })
    val2 <- local({
        spectreasy:::.with_optional_seed(999)
        runif(3)
    })

    expect_equal(val1, val2)
})

test_that("derive_unmixing_matrix fails for singular OLS matrix", {
    M <- matrix(c(
        1, 2, 3,
        2, 4, 6
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("A", "B")
    colnames(M) <- c("UV1-A", "UV2-A", "UV3-A")

    expect_error(
        spectreasy::derive_unmixing_matrix(M, method = "OLS"),
        regexp = "singular"
    )
})

test_that("calc_residuals fails when detector columns are missing", {
    M <- matrix(c(
        1, 0.2,
        0.1, 1
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B1-A", "YG1-A")

    exprs <- matrix(c(100, 1, 1000, 500), nrow = 1)
    colnames(exprs) <- c("B1-A", "Time", "FSC-A", "SSC-A")
    ff <- flowCore::flowFrame(exprs)

    expect_error(
        spectreasy::calc_residuals(ff, M, method = "OLS"),
        regexp = "Detectors in reference matrix not found"
    )
})

test_that("calculate_nps excludes AF markers by default", {
    df <- data.frame(
        File = rep(c("A", "B"), each = 30),
        FITC = rnorm(60, sd = 0.2),
        AF = rnorm(60, sd = 0.2),
        AF_2 = rnorm(60, sd = 0.2),
        check.names = FALSE
    )

    nps <- spectreasy::calculate_nps(df)
    expect_true("FITC" %in% nps$Marker)
    expect_false(any(grepl("^AF($|_)", nps$Marker, ignore.case = TRUE)))
})

test_that("qc_samples accepts unmix_samples results directly", {
    M <- matrix(c(
        1.0, 0.2, 0.1,
        0.1, 1.0, 0.2,
        0.2, 0.2, 1.0
    ), nrow = 3, byrow = TRUE)
    rownames(M) <- c("FITC", "PE", "APC")
    colnames(M) <- c("B1-A", "YG1-A", "R1-A")

    simulate_sample <- function(dominant_marker, M, n_cells = 60) {
        markers <- rownames(M)
        marker_signal <- matrix(rexp(n_cells * length(markers), rate = 8), ncol = length(markers))
        colnames(marker_signal) <- markers
        marker_signal[, dominant_marker] <- rexp(n_cells, rate = 0.6) + 2
        raw_signal <- marker_signal %*% M + matrix(rnorm(n_cells * ncol(M), sd = 0.03), ncol = ncol(M))
        exprs_mat <- cbind(
            raw_signal,
            Time = seq_len(n_cells),
            "FSC-A" = rnorm(n_cells, mean = 90000, sd = 7000),
            "SSC-A" = rnorm(n_cells, mean = 45000, sd = 5000)
        )
        colnames(exprs_mat)[seq_len(ncol(M))] <- colnames(M)
        flowCore::flowFrame(exprs_mat)
    }

    toy_fs <- flowCore::flowSet(list(
        FITC_sample = simulate_sample("FITC", M),
        PE_sample = simulate_sample("PE", M)
    ))

    unmixed <- spectreasy::unmix_samples(
        toy_fs,
        M = M,
        method = "OLS",
        output_dir = tempdir(),
        write_fcs = FALSE
    )

    expect_s3_class(unmixed, "spectreasy_unmixed_results")

    combined <- as.data.frame(unmixed)
    expect_true("File" %in% colnames(combined))
    expect_equal(sort(unique(combined$File)), c("FITC_sample", "PE_sample"))

    pdf_out <- tempfile(fileext = ".pdf")
    expect_no_error(
        suppressMessages(
            spectreasy::qc_samples(
                results = unmixed,
                M = M,
                output_file = pdf_out
            )
        )
    )
    expect_true(file.exists(pdf_out))

    # Test default output file behavior
    default_pdf <- "spectreasy_outputs/unmix_samples/qc_samples_report.pdf"
    if (file.exists(default_pdf)) file.remove(default_pdf)
    expect_no_error(
        suppressMessages(
            spectreasy::qc_samples(
                results = unmixed,
                M = M
            )
        )
    )
    expect_true(file.exists(default_pdf))
    file.remove(default_pdf)
})

test_that("qc_samples has stable pages and no recommendation page", {
    skip_if_not_installed("pdftools")

    set.seed(1)
    n <- 120
    results <- data.frame(
        FITC = rnorm(n, 0, 0.3),
        PE = rnorm(n, 0, 0.4),
        AF = rnorm(n, 0, 0.2),
        File = rep(c("SampleA", "SampleB"), each = n / 2),
        check.names = FALSE
    )

    M <- matrix(c(
        1.0, 0.2, 0.1,
        0.1, 1.0, 0.2,
        0.2, 0.2, 1.0
    ), nrow = 3, byrow = TRUE)
    rownames(M) <- c("FITC", "PE", "AF")
    colnames(M) <- c("B1-A", "YG1-A", "R1-A")

    pdf_out <- tempfile(fileext = ".pdf")
    spectreasy::qc_samples(results = results, M = M, output_file = pdf_out)

    expect_true(file.exists(pdf_out))

    info <- pdftools::pdf_info(pdf_out)
    expect_equal(info$pages, 6)

    txt <- paste(pdftools::pdf_text(pdf_out), collapse = "\n")
    expect_false(grepl("Conclusions & Recommendations", txt, fixed = TRUE))
    expect_true(grepl("Sample NxN Scatter Matrix: SampleA", txt, fixed = TRUE))
    expect_false(grepl("Sample NxN Scatter Matrix: SampleB", txt, fixed = TRUE))
    expect_false(grepl("Good: populations remain compact", txt, fixed = TRUE))
})

test_that("qc_samples skips negative population spread for NNLS", {
    skip_if_not_installed("pdftools")

    set.seed(1)
    n <- 120
    results <- data.frame(
        FITC = abs(rnorm(n, 0, 0.3)),
        PE = abs(rnorm(n, 0, 0.4)),
        AF = abs(rnorm(n, 0, 0.2)),
        File = rep(c("SampleA", "SampleB"), each = n / 2),
        check.names = FALSE
    )

    M <- matrix(c(
        1.0, 0.2, 0.1,
        0.1, 1.0, 0.2,
        0.2, 0.2, 1.0
    ), nrow = 3, byrow = TRUE)
    rownames(M) <- c("FITC", "PE", "AF")
    colnames(M) <- c("B1-A", "YG1-A", "R1-A")

    pdf_out <- tempfile(fileext = ".pdf")
    spectreasy::qc_samples(results = results, M = M, output_file = pdf_out, method = "NNLS")

    info <- pdftools::pdf_info(pdf_out)
    expect_equal(info$pages, 5)

    txt <- paste(pdftools::pdf_text(pdf_out), collapse = "\n")
    expect_false(grepl("Negative Population Spread", txt, fixed = TRUE))
})

test_that("qc_samples can include NxN pages for all samples", {
    skip_if_not_installed("pdftools")

    set.seed(1)
    n <- 120
    results <- data.frame(
        FITC = rnorm(n, 0, 0.3),
        PE = rnorm(n, 0, 0.4),
        AF = rnorm(n, 0, 0.2),
        File = rep(c("SampleA", "SampleB"), each = n / 2),
        check.names = FALSE
    )

    M <- matrix(c(
        1.0, 0.2, 0.1,
        0.1, 1.0, 0.2,
        0.2, 0.2, 1.0
    ), nrow = 3, byrow = TRUE)
    rownames(M) <- c("FITC", "PE", "AF")
    colnames(M) <- c("B1-A", "YG1-A", "R1-A")

    pdf_out <- tempfile(fileext = ".pdf")
    spectreasy::qc_samples(
        results = results,
        M = M,
        output_file = pdf_out,
        nxn_all_samples = TRUE
    )

    info <- pdftools::pdf_info(pdf_out)
    expect_equal(info$pages, 7)

    txt <- paste(pdftools::pdf_text(pdf_out), collapse = "\n")
    expect_true(grepl("Sample NxN Scatter Matrix: SampleA", txt, fixed = TRUE))
    expect_true(grepl("Sample NxN Scatter Matrix: SampleB", txt, fixed = TRUE))
})

test_that("qc_samples overview batches while moderate matrix pages stay complete", {
    marker_names <- paste0("M", seq_len(16))
    detector_names <- paste0("D", seq_len(16), "-A")
    M <- diag(16)
    rownames(M) <- marker_names
    colnames(M) <- detector_names

    res_list <- lapply(seq_len(16), function(i) {
        data <- as.data.frame(matrix(rnorm(20 * 16), nrow = 20))
        colnames(data) <- marker_names
        data$File <- paste0("Sample", i)
        list(
            data = data,
            residuals = matrix(rnorm(20 * 16, sd = 0.1), nrow = 20, dimnames = list(NULL, detector_names))
        )
    })
    names(res_list) <- paste0("Sample", seq_len(16))

    results_df <- do.call(rbind, lapply(res_list, function(x) x$data))
    nps_scores <- spectreasy::calculate_nps(results_df)
    sim_mat <- spectreasy:::calculate_similarity_matrix(M)
    ssm <- spectreasy::calculate_ssm(M)

    rms_pages <- spectreasy:::.build_qc_report_rms_pages(res_list, M = M)
    nps_pages <- spectreasy:::.build_qc_report_nps_pages(nps_scores)
    sim_pages <- spectreasy:::.build_qc_report_matrix_pages(sim_mat, plot_fun = spectreasy:::plot_similarity_matrix)
    ssm_pages <- spectreasy:::.build_qc_report_matrix_pages(ssm, plot_fun = spectreasy::plot_ssm)

    expect_equal(length(rms_pages), 2)
    expect_equal(length(nps_pages), 2)
    expect_equal(length(sim_pages), 1)
    expect_equal(length(ssm_pages), 1)

    sim_chunk_pages <- spectreasy:::.build_qc_report_matrix_pages(
        sim_mat,
        plot_fun = spectreasy:::plot_similarity_matrix,
        max_markers_per_page = 10
    )
    ssm_chunk_pages <- spectreasy:::.build_qc_report_matrix_pages(
        ssm,
        plot_fun = spectreasy::plot_ssm,
        max_markers_per_page = 10
    )
    expect_equal(length(sim_chunk_pages), 2)
    expect_equal(length(ssm_chunk_pages), 2)
})

test_that("qc_samples matrix pages split markers into balanced groups", {
    batches_20 <- spectreasy:::.split_qc_report_matrix_marker_batches(paste0("M", seq_len(20)))
    batches_30 <- spectreasy:::.split_qc_report_matrix_marker_batches(paste0("M", seq_len(30)))
    batches_36 <- spectreasy:::.split_qc_report_matrix_marker_batches(paste0("M", seq_len(36)))
    batches_41 <- spectreasy:::.split_qc_report_matrix_marker_batches(paste0("M", seq_len(41)))

    expect_equal(lengths(batches_20), 20L)
    expect_equal(lengths(batches_30), c(15L, 15L))
    expect_equal(lengths(batches_36), c(18L, 18L))
    expect_equal(lengths(batches_41), c(14L, 14L, 13L))

    make_square <- function(n) {
        M <- diag(n)
        rownames(M) <- paste0("M", seq_len(n))
        colnames(M) <- paste0("M", seq_len(n))
        M
    }

    expect_equal(length(spectreasy:::.build_qc_report_matrix_pages(
        make_square(30),
        plot_fun = spectreasy:::plot_similarity_matrix
    )), 2)
    expect_equal(length(spectreasy:::.build_qc_report_matrix_pages(
        make_square(41),
        plot_fun = spectreasy::plot_ssm
    )), 3)
})

test_that("qc_samples loads M from unmixing_matrix_file", {
    set.seed(1)
    n <- 120
    results <- data.frame(
        FITC = rnorm(n, 0, 0.3),
        PE = rnorm(n, 0, 0.4),
        AF = rnorm(n, 0, 0.2),
        File = rep(c("SampleA", "SampleB"), each = n / 2),
        check.names = FALSE
    )

    M <- matrix(c(
        1.0, 0.2, 0.1,
        0.1, 1.0, 0.2,
        0.2, 0.2, 1.0
    ), nrow = 3, byrow = TRUE)
    rownames(M) <- c("FITC", "PE", "AF")
    colnames(M) <- c("B1-A", "YG1-A", "R1-A")

    csv_file <- tempfile(fileext = ".csv")
    M_df <- as.data.frame(M, check.names = FALSE)
    M_df$Marker <- rownames(M)
    M_df <- M_df[, c("Marker", setdiff(colnames(M_df), "Marker")), drop = FALSE]
    utils::write.csv(M_df, csv_file, row.names = FALSE, quote = TRUE)

    pdf_out <- tempfile(fileext = ".pdf")
    expect_no_error(
        suppressMessages(
            spectreasy::qc_samples(
                results = results,
                unmixing_matrix_file = csv_file,
                output_file = pdf_out
            )
        )
    )
    expect_true(file.exists(pdf_out))
})

test_that("per-cell AF selection unmix selects exactly one AF band per cell", {
    # Create reference matrix with multiple AF bands
    M <- matrix(c(
        1.0, 0.2, 0.1, 0.05,
        0.1, 1.0, 0.2, 0.05,
        0.2, 0.2, 1.0, 0.05,
        0.05, 0.05, 0.05, 1.0
    ), nrow = 4, byrow = TRUE)
    rownames(M) <- c("FITC", "PE", "AF", "AF_2")
    colnames(M) <- c("B1-A", "YG1-A", "R1-A", "V1-A")
    V <- matrix(1, nrow = nrow(M), ncol = ncol(M), dimnames = dimnames(M))
    attr(M, "variances") <- V

    # Construct some cells where we know which AF band they should select
    # Cell 1: high FITC, some PE, high AF (R1-A)
    # Cell 2: high PE, some FITC, high AF_2 (V1-A)
    exprs <- matrix(c(
        100, 10, 50, 0,
        10, 100, 0, 50
    ), nrow = 2, byrow = TRUE)
    colnames(exprs) <- c("B1-A", "YG1-A", "R1-A", "V1-A")
    
    exprs_full <- cbind(
        exprs,
        Time = c(1, 2),
        "FSC-A" = c(90000, 92000),
        "SSC-A" = c(45000, 46000)
    )
    ff <- flowCore::flowFrame(exprs_full)

    # Run calc_residuals with different methods: OLS, WLS, RWLS, NNLS
    for (m in c("OLS", "WLS", "RWLS", "NNLS")) {
        res <- spectreasy::calc_residuals(ff, M, method = m)
        
        # Verify that for each cell, at most one AF band is non-zero
        expect_equal(nrow(res), 2)
        
        # For Cell 1
        expect_true(sum(res[1, c("AF", "AF_2")] != 0) <= 1)
        # For Cell 2
        expect_true(sum(res[2, c("AF", "AF_2")] != 0) <= 1)
    }
})

test_that("unmix_samples can dynamically construct reference matrix and handle multi-AF", {
    # 1. Create a temporary SCC directory
    scc_dir <- tempfile("scc_dir_")
    dir.create(scc_dir)
    
    # 2. Generate a flowFrame with some channels: "B1-A", "YG1-A"
    # FITC control
    fitc_exprs <- cbind(
        "B1-A" = rnorm(200, mean = 1000, sd = 50),
        "YG1-A" = rnorm(200, mean = 10, sd = 2),
        "FSC-A" = rnorm(200, mean = 100000, sd = 5000),
        "SSC-A" = rnorm(200, mean = 50000, sd = 2500)
    )
    fitc_ff <- flowCore::flowFrame(fitc_exprs)
    flowCore::write.FCS(fitc_ff, file.path(scc_dir, "FITC (Beads).fcs"))
    
    # PE control
    pe_exprs <- cbind(
        "B1-A" = rnorm(200, mean = 10, sd = 2),
        "YG1-A" = rnorm(200, mean = 1500, sd = 75),
        "FSC-A" = rnorm(200, mean = 100000, sd = 5000),
        "SSC-A" = rnorm(200, mean = 50000, sd = 2500)
    )
    pe_ff <- flowCore::flowFrame(pe_exprs)
    flowCore::write.FCS(pe_ff, file.path(scc_dir, "PE (Beads).fcs"))
    
    # Unstained / AF control
    af_exprs <- cbind(
        "B1-A" = rnorm(200, mean = 15, sd = 3),
        "YG1-A" = rnorm(200, mean = 15, sd = 3),
        "FSC-A" = rnorm(200, mean = 100000, sd = 5000),
        "SSC-A" = rnorm(200, mean = 50000, sd = 2500)
    )
    af_ff <- flowCore::flowFrame(af_exprs)
    flowCore::write.FCS(af_ff, file.path(scc_dir, "Unstained (Cells).fcs"))
    
    # 3. Create a control_file mapping
    control_df <- data.frame(
        filename = c("FITC (Beads).fcs", "PE (Beads).fcs", "Unstained (Cells).fcs"),
        fluorophore = c("FITC", "PE", "AF"),
        marker = c("CD4", "CD8", "Autofluorescence"),
        channel = c("B1-A", "YG1-A", "B1-A"),
        control.type = c("beads", "beads", "cells"),
        universal.negative = c("", "", ""),
        large.gate = c("", "", ""),
        is.viability = c("", "", ""),
        stringsAsFactors = FALSE
    )
    control_file <- tempfile(fileext = ".csv")
    write.csv(control_df, control_file, row.names = FALSE)
    
    # 4. Create sample flowSet to unmix
    sample_exprs <- matrix(c(
        500, 200, 1, 100000, 50000,
        100, 800, 2, 100000, 50000
    ), nrow = 2, byrow = TRUE)
    colnames(sample_exprs) <- c("B1-A", "YG1-A", "Time", "FSC-A", "SSC-A")
    sample_ff <- flowCore::flowFrame(sample_exprs)
    toy_fs <- flowCore::flowSet(list(Sample1 = sample_ff))
    
    # 5. Call unmix_samples with scc_dir and af_n_bands = 2
    # Ensure unmixing_matrix_file points to a non-existent temp file so it triggers build
    temp_matrix_file <- tempfile(fileext = ".csv")
    
    unmixed <- spectreasy::unmix_samples(
        sample_dir = toy_fs,
        unmixing_matrix_file = temp_matrix_file,
        scc_dir = scc_dir,
        control_file = control_file,
        af_n_bands = 2,
        write_fcs = FALSE,
        return_type = "list"
    )
    
    expect_s3_class(unmixed, "spectreasy_unmixed_results")
    unmixed_df <- as.data.frame(unmixed)
    
    # Check that columns FITC, PE, AF, and AF_2 exist in the output
    expect_true(all(c("FITC", "PE", "AF", "AF_2") %in% colnames(unmixed_df)))
    
    # Verify that selection unmix worked: for each cell, at most one of AF and AF_2 is non-zero
    expect_true(sum(unmixed_df[1, c("AF", "AF_2")] != 0) <= 1)
    expect_true(sum(unmixed_df[2, c("AF", "AF_2")] != 0) <= 1)
})

test_that("unmix_controls supports af_n_bands and include_multi_af", {
    scc_dir <- tempfile("scc_dir_controls_")
    dir.create(scc_dir)
    
    fitc_exprs <- cbind(
        "B1-A" = rnorm(200, mean = 1000, sd = 50),
        "YG1-A" = rnorm(200, mean = 10, sd = 2),
        "V1-A" = rnorm(200, mean = 5, sd = 1),
        "R1-A" = rnorm(200, mean = 5, sd = 1),
        "FSC-A" = rnorm(200, mean = 100000, sd = 5000),
        "SSC-A" = rnorm(200, mean = 50000, sd = 2500)
    )
    flowCore::write.FCS(flowCore::flowFrame(fitc_exprs), file.path(scc_dir, "FITC (Beads).fcs"))
    
    pe_exprs <- cbind(
        "B1-A" = rnorm(200, mean = 10, sd = 2),
        "YG1-A" = rnorm(200, mean = 1500, sd = 75),
        "V1-A" = rnorm(200, mean = 5, sd = 1),
        "R1-A" = rnorm(200, mean = 5, sd = 1),
        "FSC-A" = rnorm(200, mean = 100000, sd = 5000),
        "SSC-A" = rnorm(200, mean = 50000, sd = 2500)
    )
    flowCore::write.FCS(flowCore::flowFrame(pe_exprs), file.path(scc_dir, "PE (Beads).fcs"))
    
    af_exprs <- cbind(
        "B1-A" = rnorm(200, mean = 15, sd = 3),
        "YG1-A" = rnorm(200, mean = 15, sd = 3),
        "V1-A" = rnorm(200, mean = 15, sd = 3),
        "R1-A" = rnorm(200, mean = 15, sd = 3),
        "FSC-A" = rnorm(200, mean = 100000, sd = 5000),
        "SSC-A" = rnorm(200, mean = 50000, sd = 2500)
    )
    flowCore::write.FCS(flowCore::flowFrame(af_exprs), file.path(scc_dir, "Unstained (Cells).fcs"))
    
    control_df <- data.frame(
        filename = c("FITC (Beads).fcs", "PE (Beads).fcs", "Unstained (Cells).fcs"),
        fluorophore = c("FITC", "PE", "AF"),
        marker = c("CD4", "CD8", "Autofluorescence"),
        channel = c("B1-A", "YG1-A", "B1-A"),
        control.type = c("beads", "beads", "cells"),
        universal.negative = c("", "", ""),
        large.gate = c("", "", ""),
        is.viability = c("", "", ""),
        stringsAsFactors = FALSE
    )
    control_file <- tempfile(fileext = ".csv")
    write.csv(control_df, control_file, row.names = FALSE)
    
    output_dir <- tempfile("controls_out_")
    
    res <- spectreasy::unmix_controls(
        scc_dir = scc_dir,
        control_file = control_file,
        auto_create_control = FALSE,
        output_dir = output_dir,
        af_n_bands = 2,
        include_multi_af = FALSE
    )
    
    expect_true(all(c("FITC", "PE", "AF", "AF_2") %in% rownames(res$M)))
})

test_that("Plumber gui_api load_matrix and save_matrix filter and merge AF rows", {
    api_file <- system.file("api/gui_api.R", package = "spectreasy")
    if (api_file == "") {
        api_file <- "../../inst/api/gui_api.R"
    }
    
    expect_true(file.exists(api_file))
    
    pr <- plumber::plumb(api_file)
    load_matrix_fn <- pr$routes$load_matrix$getFunc()
    save_matrix_fn <- pr$routes$save_matrix$getFunc()
    
    tmp_matrix_dir <- tempfile("matrix_dir_")
    dir.create(tmp_matrix_dir)
    
    old_opt <- getOption("spectreasy.matrix_dir")
    options(spectreasy.matrix_dir = tmp_matrix_dir)
    on.exit(options(spectreasy.matrix_dir = old_opt), add = TRUE)
    
    dummy_csv <- data.frame(
        Marker = c("FITC", "PE", "AF", "AF_2"),
        "B1-A" = c(1, 0.1, 0.05, 0.06),
        "YG1-A" = c(0.05, 1, 0.04, 0.03),
        "R1-A" = c(0.02, 0.03, 1, 0.7),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    write.csv(dummy_csv, file.path(tmp_matrix_dir, "ref_matrix.csv"), row.names = FALSE)
    
    loaded_df <- load_matrix_fn("ref_matrix.csv")
    
    expect_equal(nrow(loaded_df), 2)
    expect_equal(loaded_df$Marker, c("FITC", "PE"))
    
    loaded_df[loaded_df$Marker == "PE", "B1-A"] <- 0.15
    
    req <- new.env()
    req$postBody <- jsonlite::toJSON(list(
        filename = "ref_matrix.csv",
        matrix_json = loaded_df
    ), auto_unbox = TRUE)
    
    res <- save_matrix_fn(req)
    expect_true(res$success)
    
    saved_df <- read.csv(file.path(tmp_matrix_dir, "ref_matrix.csv"), stringsAsFactors = FALSE, check.names = FALSE)
    
    expect_equal(nrow(saved_df), 4)
    expect_true(all(c("FITC", "PE", "AF", "AF_2") %in% saved_df$Marker))
    expect_equal(saved_df[saved_df$Marker == "PE", "B1-A"], 0.15)
    expect_equal(saved_df[saved_df$Marker == "AF", "B1-A"], 0.05)
    expect_equal(saved_df[saved_df$Marker == "AF_2", "YG1-A"], 0.03)

    req_adjusted <- new.env()
    req_adjusted$postBody <- jsonlite::toJSON(list(
        filename = "ref_matrix_adjusted.csv",
        source_filename = "ref_matrix.csv",
        matrix_json = loaded_df
    ), auto_unbox = TRUE)

    adjusted_res <- save_matrix_fn(req_adjusted)
    expect_true(adjusted_res$success)
    adjusted_df <- read.csv(file.path(tmp_matrix_dir, "ref_matrix_adjusted.csv"), stringsAsFactors = FALSE, check.names = FALSE)
    expect_true(all(c("FITC", "PE", "AF", "AF_2") %in% adjusted_df$Marker))
    expect_equal(adjusted_df[adjusted_df$Marker == "AF", "R1-A"], 1)

    unmix_fn <- pr$routes$unmix[[2]]$getFunc()
    matrix_json <- list(
        FITC = as.list(loaded_df[loaded_df$Marker == "FITC", setdiff(colnames(loaded_df), "Marker"), drop = FALSE]),
        PE = as.list(loaded_df[loaded_df$Marker == "PE", setdiff(colnames(loaded_df), "Marker"), drop = FALSE])
    )
    raw_data_json <- list(
        list("B1-A" = 100, "YG1-A" = 20, "R1-A" = 50),
        list("B1-A" = 10, "YG1-A" = 120, "R1-A" = 40)
    )
    preview <- unmix_fn(
        matrix_json = matrix_json,
        raw_data_json = raw_data_json,
        type = "reference",
        matrix_filename = "ref_matrix.csv",
        method = "OLS"
    )
    expect_true(is.data.frame(preview))
    expect_true("AF" %in% colnames(preview))
})

test_that("Plumber gui_api serves spectral panel payloads", {
    api_file <- system.file("api/gui_api.R", package = "spectreasy")
    if (api_file == "") {
        api_file <- "../../inst/api/gui_api.R"
    }

    expect_true(file.exists(api_file))

    pr <- plumber::plumb(api_file)
    spectral_panel_fn <- pr$routes$spectral_panel$getFunc()
    spectral_metrics_fn <- pr$routes$spectral_panel_metrics[[2]]$getFunc()

    old_opts <- options(spectreasy.panel_cytometer = "aurora")
    on.exit(options(old_opts), add = TRUE)

    payload <- spectral_panel_fn("")
    expect_equal(payload$cytometer, "aurora")
    expect_equal(payload$selected, character())
    expect_true(nrow(payload$detectors) > 20)

    req <- new.env()
    req$postBody <- jsonlite::toJSON(list(
        cytometer = "id7000",
        fluorophores = c("FITC", "PE", "APC")
    ), auto_unbox = TRUE)

    recalculated <- spectral_metrics_fn(req)
    expect_equal(recalculated$cytometer, "id7000")
    expect_equal(recalculated$selected, c("FITC", "PE", "APC"))
    expect_true(is.numeric(recalculated$complexity_index))
})

test_that("unmix_samples supports in-memory subsampling via subsample_n", {
    M <- matrix(c(
        1.0, 0.2,
        0.1, 1.0
    ), nrow = 2, byrow = TRUE)
    rownames(M) <- c("FITC", "PE")
    colnames(M) <- c("B1-A", "YG1-A")

    exprs <- matrix(rnorm(100, mean = 200, sd = 20), ncol = 2)
    colnames(exprs) <- c("B1-A", "YG1-A")
    exprs_full <- cbind(
        exprs,
        Time = seq_len(50),
        "FSC-A" = rnorm(50, mean = 90000, sd = 7000),
        "SSC-A" = rnorm(50, mean = 45000, sd = 5000)
    )
    ff <- flowCore::flowFrame(exprs_full)
    toy_fs <- flowCore::flowSet(list(Sample1 = ff))

    tmp_dir <- tempfile("unmix_subsample_out_")
    dir.create(tmp_dir)

    unmixed <- spectreasy::unmix_samples(
        toy_fs,
        M = M,
        method = "OLS",
        output_dir = tmp_dir,
        write_fcs = TRUE,
        subsample_n = 10,
        seed = 42
    )

    unmixed_df <- as.data.frame(unmixed)
    expect_equal(nrow(unmixed_df), 10)

    fcs_path <- file.path(tmp_dir, "Sample1_unmixed.fcs")
    expect_true(file.exists(fcs_path))
    ff_written <- flowCore::read.FCS(fcs_path, transformation = FALSE, truncate_max_range = FALSE)
    expect_equal(nrow(flowCore::exprs(ff_written)), 50)
})

test_that("unmix_samples excludes secondary AF bands from written FCS files", {
    M <- matrix(c(
        1.0, 0.2, 0.1, 0.05,
        0.1, 1.0, 0.2, 0.05,
        0.2, 0.2, 1.0, 0.05,
        0.05, 0.05, 0.05, 1.0
    ), nrow = 4, byrow = TRUE)
    rownames(M) <- c("FITC", "PE", "AF", "AF_2")
    colnames(M) <- c("B1-A", "YG1-A", "R1-A", "V1-A")

    exprs <- matrix(rnorm(200, mean = 200, sd = 20), ncol = 4)
    colnames(exprs) <- c("B1-A", "YG1-A", "R1-A", "V1-A")
    exprs_full <- cbind(
        exprs,
        Time = seq_len(50),
        "FSC-A" = rnorm(50, mean = 90000, sd = 7000),
        "SSC-A" = rnorm(50, mean = 45000, sd = 5000)
    )
    ff <- flowCore::flowFrame(exprs_full)
    toy_fs <- flowCore::flowSet(list(Sample1 = ff))

    tmp_dir <- tempfile("unmix_af_out_")
    dir.create(tmp_dir)

    unmixed <- spectreasy::unmix_samples(
        toy_fs,
        M = M,
        method = "OLS",
        output_dir = tmp_dir,
        write_fcs = TRUE
    )

    fcs_path <- file.path(tmp_dir, "Sample1_unmixed.fcs")
    expect_true(file.exists(fcs_path))
    ff_written <- flowCore::read.FCS(fcs_path, transformation = FALSE, truncate_max_range = FALSE)
    written_cols <- colnames(flowCore::exprs(ff_written))
    expect_true("AF" %in% written_cols)
    expect_false("AF_2" %in% written_cols)
})
