library(QFeatures)

data(feat3)
feat3 <- setQFeaturesType(feat3, "scp")
preparedFeat3 <- suppressWarnings(prepareQFeatures(feat3))

test_that(".h5mu file is created", {
    filePath <- tempfile(fileext = ".h5mu")
    writeLinkH5MU(preparedFeat3, filePath)

    testthat::expect_true(file.exists(filePath))
})

test_that("QFeatures after write/read is identical", {
    filePath <- tempfile(fileext = ".h5mu")
    writeLinkH5MU(preparedFeat3, filePath)

    newQFeatures <- readLinkH5MU(filePath)
    testthat::expect_equal(names(newQFeatures), names(feat3))
    testthat::expect_equal(colnames(newQFeatures), colnames(feat3))
    #testthat::expect_equal(rownames(newQFeatures), rownames(feat3))
    testthat::expect_equal(names(newQFeatures@assayLinks),
        names(feat3@assayLinks))
    for(set in names(feat3)) {
        testthat::expect_equal(newQFeatures@assayLinks[[set]]@from,
            feat3@assayLinks[[set]]@from)
        # testthat::expect_equal(newQFeatures@assayLinks[[set]]@hits,
        #     feat3@assayLinks[[set]]@hits)
        testthat::expect_equal(newQFeatures@assayLinks[[set]]@name,
            feat3@assayLinks[[set]]@name)
    }
    testthat::expect_equal(getQFeaturesType(feat3),
        getQFeaturesType(newQFeatures))
})
