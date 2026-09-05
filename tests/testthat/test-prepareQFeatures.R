test_that("prepareQFeatures only prefixes duplicated feature names", {
    se1 <- SummarizedExperiment::SummarizedExperiment(
        matrix(seq_len(4), nrow = 2,
               dimnames = list(c("a", "b"), c("s1", "s2")))
    )
    se2 <- SummarizedExperiment::SummarizedExperiment(
        matrix(seq_len(4), nrow = 2,
               dimnames = list(c("c", "d"), c("s1", "s2")))
    )

    qf <- QFeatures::QFeatures(list(set1 = se1, set2 = se2))
    expect_no_warning(prepared <- prepareQFeatures(qf))
    expect_identical(rownames(prepared[["set1"]]), c("a", "b"))

    rownames(se2) <- c("a", "c")
    qf <- QFeatures::QFeatures(list(set1 = se1, set2 = se2))
    expect_warning(prepared <- prepareQFeatures(qf),
                   "Prefixed feature row names")
    expect_identical(rownames(prepared[["set1"]]), c("set1:a", "set1:b"))
    expect_identical(rownames(prepared[["set2"]]), c("set2:a", "set2:c"))
})
