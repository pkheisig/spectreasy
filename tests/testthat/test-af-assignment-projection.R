test_that("projection AF assignment recovers structured AF candidates", {
    set.seed(11)
    marker_M <- rbind(
        F1 = c(1, 0.2, 0),
        F2 = c(0, 0.2, 1)
    )
    af_M <- rbind(
        AF = c(0.8, 1, 0.2),
        AF_2 = c(0.1, 1, 0.9)
    )
    colnames(marker_M) <- colnames(af_M) <- c("D1", "D2", "D3")

    n <- 300
    true_af <- rep(1:2, each = n / 2)
    af_amount <- runif(n, 50, 150)
    Y <- af_amount * af_M[true_af, ] + matrix(rnorm(n * 3, sd = 1.5), ncol = 3)
    colnames(Y) <- colnames(marker_M)

    projection <- spectreasy:::.assign_af_by_projection(Y, marker_M, af_M, return_details = TRUE)

    expect_gt(mean(projection$assignments == true_af), 0.85)
    expect_equal(dim(projection$score_matrix), c(n, 2))
    expect_equal(colnames(projection$projected_af_fluor), rownames(marker_M))
})

test_that("projection AF assignment C++ fast path matches R detailed path", {
    set.seed(111)
    marker_M <- rbind(
        F1 = c(1, 0.2, 0),
        F2 = c(0, 0.2, 1)
    )
    af_M <- rbind(
        AF = c(0.8, 1, 0.2),
        AF_2 = c(0.1, 1, 0.9)
    )
    colnames(marker_M) <- colnames(af_M) <- c("D1", "D2", "D3")
    true_af <- rep(1:2, each = 15)
    Y <- runif(length(true_af), 20, 80) * af_M[true_af, ] + matrix(rnorm(length(true_af) * 3, sd = 0.01), ncol = 3)
    colnames(Y) <- colnames(marker_M)

    fast <- spectreasy:::.assign_af_by_projection(Y, marker_M, af_M)
    detailed <- spectreasy:::.assign_af_by_projection(Y, marker_M, af_M, return_details = TRUE)

    expect_equal(fast, detailed$assignments)
})

test_that("projection AF assignment does not create strong fluorophore correction from pure detector noise", {
    set.seed(12)
    marker_M <- rbind(
        F1 = c(1, 0.2, 0),
        F2 = c(0, 0.2, 1)
    )
    af_M <- rbind(
        AF = c(0.8, 1, 0.2),
        AF_2 = c(0.1, 1, 0.9)
    )
    colnames(marker_M) <- colnames(af_M) <- c("D1", "D2", "D3")

    Y <- matrix(rnorm(400 * 3, sd = 0.03), ncol = 3)
    colnames(Y) <- colnames(marker_M)
    projection <- spectreasy:::.assign_af_by_projection(Y, marker_M, af_M, return_details = TRUE)

    expect_lt(stats::median(abs(projection$af_abundance)), 0.1)
    expect_lt(stats::median(rowSums(abs(projection$projected_af_fluor))), 0.1)
    expect_lt(abs(mean(projection$assignments == 1L) - 0.5), 0.35)
})

test_that("projection multi-AF unmixing does not steal real fluorophore signal", {
    set.seed(13)
    M <- rbind(
        F1 = c(1, 0.2, 0),
        F2 = c(0, 0.2, 1),
        AF = c(0.8, 1, 0.2),
        AF_2 = c(0.1, 1, 0.9)
    )
    colnames(M) <- c("D1", "D2", "D3")

    n <- 240
    true_af <- rep(1:2, each = n / 2)
    true_fluor <- cbind(F1 = runif(n, 20, 80), F2 = runif(n, 10, 40))
    true_fluor[true_af == 1, "F2"] <- 0
    true_fluor[true_af == 2, "F1"] <- 0
    true_af_amount <- runif(n, 40, 120)
    Y <- true_fluor %*% M[c("F1", "F2"), ] +
        true_af_amount * M[c("AF", "AF_2")[true_af], ] +
        matrix(rnorm(n * 3, sd = 0.5), ncol = 3)
    colnames(Y) <- colnames(M)
    ff <- flowCore::flowFrame(Y)

    projection <- spectreasy::calc_residuals(ff, M, method = "OLS")

    projection_error <- mean(abs(as.matrix(projection[, c("F1", "F2")]) - true_fluor))
    expect_lt(projection_error, 1)
})
