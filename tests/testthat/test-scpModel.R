data(leduc_minimal, package = "scp")

test_that("Convert scpModel to/from list", {
    mod0 <- metadata(leduc_minimal)$model
    mod0@scpModelFitList <- mod0@scpModelFitList[1:10]
    expect_true(validObject(mod0))
    ll <- scpModelToList(mod0)
    expect_true(validObject(mod <- listToScpModel(ll)))
    expect_equal(mod0, mod, ignore_attr = TRUE)
})

test_that("hasScpModel works", {
    expect_true(hasScpModel(leduc_minimal))
    metadata(leduc_minimal) <- list()
    expect_false(hasScpModel(leduc_minimal))
    metadata(leduc_minimal) <- list(1, letters, ls)
    expect_false(hasScpModel(leduc_minimal))
})

## reset object
data(leduc_minimal)

test_that("scpModelAsList works", {
    expect_s4_class(metadata(leduc_minimal)[[1]], "ScpModel")
    se <- scpModelAsList(leduc_minimal)
    expect_type(metadata(se)[[1]], "list")
    metadata(leduc_minimal) <- list()
    se <- scpModelAsList(leduc_minimal)
    expect_identical(leduc_minimal, se)
    metadata(leduc_minimal) <- list(1, letters, ls)
    se <- scpModelAsList(leduc_minimal)
    expect_identical(leduc_minimal, se)
})