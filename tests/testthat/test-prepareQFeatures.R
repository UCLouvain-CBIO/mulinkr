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

test_that("prepareQFeatures convert date column", {
    se <- SummarizedExperiment::SummarizedExperiment(
        matrix(seq_len(4), nrow = 2,
               dimnames = list(c("a", "b"), c("s1", "s2"))),
        rowData = S4Vectors::DataFrame(
            feature_date = as.Date(c("2024-01-01", "2024-01-02")),
            feature_name = c("a", "b")
        ),
        colData = S4Vectors::DataFrame(
            acquisition_time = as.POSIXct(
                c("2024-02-03 04:05:06", "2024-02-04 05:06:07"),
                tz = "UTC"
            )
        )
    )
    qf <- QFeatures::QFeatures(list(set1 = se))
    SummarizedExperiment::colData(qf) <- S4Vectors::DataFrame(
        sample_date = as.Date(c("2024-03-01", "2024-03-02")),
        sample_group = c("A", "B"),
        row.names = c("s1", "s2")
    )

    prepared <- prepareQFeatures(qf)

    expect_identical(
        SummarizedExperiment::colData(prepared)$sample_date,
        c("2024-03-01", "2024-03-02")
    )
    expect_identical(
        SummarizedExperiment::rowData(prepared[["set1"]])$feature_date,
        c("2024-01-01", "2024-01-02")
    )
    expect_identical(
        SummarizedExperiment::colData(prepared[["set1"]])$acquisition_time,
        c("2024-02-03 04:05:06", "2024-02-04 05:06:07")
    )
    expect_identical(
        SummarizedExperiment::rowData(prepared[["set1"]])$feature_name,
        c("a", "b")
    )
    expect_identical(
        SummarizedExperiment::colData(prepared)$sample_group,
        c("A", "B")
    )
})

test_that("prepareQFeatures transform scpModel slot in metadata", {
    data("leduc_minimal", package = "scp")
    qf <- QFeatures::QFeatures(list(scp = leduc_minimal))

    expect_s4_class(S4Vectors::metadata(qf[["scp"]])$model, "ScpModel")

    expect_warning(prepared <- prepareQFeatures(qf), "has a scpModel")
    model <- S4Vectors::metadata(prepared[["scp"]])$model

    expect_type(model, "list")
    expect_named(
        model,
        c("scpModelFormula", "scpModelInputIndex", "scpModelFilterThreshold",
          "scpModelFitList")
    )
    expect_false(inherits(model, "ScpModel"))
})
