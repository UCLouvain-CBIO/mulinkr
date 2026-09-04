library(QFeatures)

data(feat3)

test_that(".h5mu file is created", {
    filePath <- tempfile(fileext = ".h5mu")
    writeLinkH5MU(feat3, filePath)

    testthat::expect_true(file.exists(filePath))
})

test_that("QFeatures after write/read is identical", {
    filePath <- tempfile(fileext = ".h5mu")
    writeLinkH5MU(feat3, filePath)

    newQFeatures <- readLinkH5MU(filePath)
    testthat::expect_equal(names(newQFeatures), names(feat3))
    testthat::expect_equal(colnames(newQFeatures), colnames(feat3))
    testthat::expect_equal(rownames(newQFeatures), rownames(feat3))
    testthat::expect_equal(newQFeatures@assayLinks, feat3@assayLinks)
})
