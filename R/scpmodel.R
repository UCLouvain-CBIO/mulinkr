##' @importFrom stats as.formula
##' @importClassesFrom scp ScpModel ScpModelFit
scpModelToList <- function(mod) {
    stopifnot(inherits(mod, "ScpModel"))
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
    list(
        scpModelFormula = paste(as.character(mod@scpModelFormula), collapse = " "),
        scpModelInputIndex = mod@scpModelInputIndex,
        scpModelFilterThreshold = mod@scpModelFilterThreshold,
        scpModelFitList = lapply(mod@scpModelFitList,
                                 scpModelFitListToList))
}

##' @importFrom S4Vectors List
listToScpModel <- function(x) {
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
    new("ScpModel",
        scpModelFormula = as.formula(x$scpModelFormula),
        scpModelInputIndex = x$scpModelInputIndex,
        scpModelFilterThreshold = x$scpModelFilterThreshold,
        scpModelFitList = List(lapply(x$scpModelFitList, listToScpModelFit)))
}

##' @importFrom S4Vectors metadata
hasScpModel <- function(object) {
    stopifnot(inherits(object, "SummarizedExperiment"))
    if (!length(metadata(object)))
        return(FALSE)
    any(sapply(metadata(object),
               function(x) inherits(x, "ScpModel")))
}

##' @importFrom S4Vectors metadata<-
scpModelAsList <- function(object) {
    stopifnot(inherits(object, "SummarizedExperiment"))
    if (hasScpModel(object)) {
        i <- which(sapply(metadata(object),
                          function(x) inherits(x, "ScpModel")))
        for (j in i)
            metadata(object)[[j]] <- scpModelToList(metadata(object)[[j]])
    }
    object
}