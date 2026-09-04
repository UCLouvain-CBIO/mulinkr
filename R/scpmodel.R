scpModelFitListToList <- function(fit) {
    stopifnot(inherits(fit, "ScpModelFit"))
    list(n = fit@n,
         p = fit@p,
         coefficients = fit@coefficients,
         residuals = fit@residuals,
         effects = as.list(fit@effects),
         df = fit@df,
         var = fit@var,
         uvcov = fit@uvcov,
         levels = as.list(fit@levels))
}

scpModelToList <- function(mod) {
    stopifnot(inherits(mod, "ScpModel"))
    list(
        scpModelFormula = paste(as.character(mod@scpModelFormula), collapse = " "),
        scpModelInputIndex = mod@scpModelInputIndex,
        scpModelFilterThreshold = mod@scpModelFilterThreshold,
        scpModelFitList = lapply(mod@scpModelFitList,
                                 scpModelFitListToList))
}

listToScpModelFit <- function(x)
    new("ScpModelFit",
        n = x$n,
        p = x$p,
        coefficients = x$coefficients,
        residuals = x$residuals,
        effects = List(x$effects),
        df = x$df,
        var = x$var,
        uvcov = x$uvcov,
        levels = List(x$levels))


listToScpModel <- function(x)
    new("ScpModel",
        scpModelFormula = as.formula(x$scpModelFormula),
        scpModelInputIndex = x$scpModelInputIndex,
        scpModelFilterThreshold = x$scpModelFilterThreshold,
        scpModelFitList = List(lapply(ll$scpModelFitList, listToScpModelFit)))
