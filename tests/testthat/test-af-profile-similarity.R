similarity_profile <- function(name, matrix) {
    spectreasy:::.new_af_profile_object(name = name, profile = matrix)
}

test_that("one-profile AF similarity is symmetric with a unit diagonal", {
    profile <- rbind(AF = c(1, 0, 0), AF_2 = c(0.5, 0.5, 0), AF_3 = c(0, 1, 0))
    colnames(profile) <- c("B1-A", "YG1-A", "R1-A")
    result <- spectreasy:::.af_profile_similarity(similarity_profile("PBMC", profile))

    expect_equal(result$similarity, t(result$similarity))
    expect_equal(unname(diag(result$similarity)), rep(1, 3))
    expect_equal(unname(result$similarity[1, 3]), 0)
    expect_equal(result$labels, paste("PBMC", rownames(profile), sep = " · "))
    expect_identical(result$profile_membership, rep("PBMC", 3))
})

test_that("two-profile AF similarity aligns identical detector sets", {
    primary <- rbind(AF = c(1, 0, 0), AF_2 = c(0, 1, 0))
    colnames(primary) <- c("B1-A", "YG1-A", "R1-A")
    comparison <- rbind(AF = c(0, 1, 0), AF_2 = c(0, 0, 1))
    colnames(comparison) <- c("R1-A", "B1-A", "YG1-A")

    result <- spectreasy:::.af_profile_similarity(
        similarity_profile("Primary", primary),
        similarity_profile("Comparison", comparison)
    )
    expect_equal(dim(result$similarity), c(4, 4))
    expect_equal(result$similarity, t(result$similarity))
    expect_equal(unname(diag(result$similarity)), rep(1, 4))
    expect_equal(unname(result$similarity[1, 3]), 1)
    expect_equal(unname(result$similarity[2, 3]), 0)
    expect_identical(result$detectors, colnames(primary))
})

test_that("AF similarity rejects detector mismatches and handles zero vectors", {
    primary <- rbind(AF = c(1, 0), AF_2 = c(0, 0))
    colnames(primary) <- c("B1-A", "R1-A")
    mismatch <- rbind(AF = c(1, 0))
    colnames(mismatch) <- c("B1-A", "YG1-A")

    result <- spectreasy:::.af_profile_similarity(similarity_profile("Primary", primary))
    expect_equal(unname(diag(result$similarity)), c(1, 1))
    expect_equal(unname(result$similarity[1, 2]), 0)
    expect_true(all(is.finite(result$similarity)))
    expect_error(
        spectreasy:::.af_profile_similarity(
            similarity_profile("Primary", primary),
            similarity_profile("Mismatch", mismatch)
        ),
        "different detector sets"
    )
})
