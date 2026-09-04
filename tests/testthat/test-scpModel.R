data(leduc_minimal)

test_that("Convert scpModel to/from list", {
    mod0 <- metadata(leduc_minimal)$model
    mod0@scpModelFitList <- mod0@scpModelFitList[1:10]
    expect_true(validObject(mod0))
    ll <- scpModelToList(mod0)
    expect_true(validObject(mod <- listToScpModel(ll)))
    expect_equivalent(mod0, mod)
})
