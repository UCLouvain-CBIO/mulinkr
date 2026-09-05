#' Prepare a `QFeatures` object for .h5mu writing.
#'
#' `prepareQFeatures()` makes feature row names globally unique across assays
#' by prefixing them with their assay name when needed. It also change any date
#' columns of the global `colData`, individual `rowData` and individual
#' `colData` to a formatted character type.
#'
#' @param object A `QFeatures` object.
#' @param sep Separator between the assay name and feature name.
#'
#' @return A prepared `QFeatures` object.
#' @importFrom MultiAssayExperiment ExperimentList experiments
#' @importFrom SummarizedExperiment "colData<-" "rowData<-" colData rowData
#' @examples
#' data("feat3", package = "QFeatures")
#' rownames(feat3)
#' preparedFeat3 <- prepareQFeatures(feat3)
#' rownames(preparedFeat3)
#' @export
prepareQFeatures <- function(object, sep = ":") {
    object <- makeRownamesUnique(object, sep)
    colData(object) <- convertDateColumns(colData(object))
    cleaned <- lapply(experiments(object), function(set) {
        rowData(set) <- convertDateColumns(rowData(set))
        colData(set) <- convertDateColumns(colData(set))
        if (hasScpModel(set)) {
            warning("Your QFeatures object has a scpModel object in one of ",
                "his set, it will be converted to a simplier list",
                "representation")
            set <- scpModelAsList(set)
        }
        set
    })
    object@ExperimentList <- ExperimentList(cleaned)
    object
}
#' @noRd
convertDateColumns <- function(dfr) {
    for (col in names(dfr)) {
        x <- dfr[[col]]

        if (inherits(x, "POSIXt")) {
            dfr[[col]] <- format(as.POSIXct(x), "%Y-%m-%d %H:%M:%S")
        } else if (inherits(x, "Date")) {
            dfr[[col]] <- format(x, "%Y-%m-%d")
        }
    }

    dfr
}

#' @importClassesFrom QFeatures QFeatures
#' @importFrom methods is validObject
#' @importFrom MultiAssayExperiment ExperimentList experiments
#' @noRd
makeRownamesUnique <- function(object, sep) {
    if (!is(object, "QFeatures"))
        stop("'object' must be a QFeatures object.")

    assay_names <- names(object)
    old_names <- lapply(assay_names, function(name) rownames(object[[name]]))
    names(old_names) <- assay_names

    if (!anyDuplicated(unlist(old_names, use.names = FALSE)))
        return(object)

    new_names <- mapply(function(name, rows) paste0(name, sep, rows),
                        assay_names, old_names, SIMPLIFY = FALSE)

    prepared <- mapply(function(se, rows) {
        rownames(se) <- rows
        se
    }, experiments(object), new_names, SIMPLIFY = FALSE)

    object@ExperimentList <- ExperimentList(prepared)
    object@assayLinks <- rename_assay_link_features(object@assayLinks,
                                                    old_names, new_names)

    validObject(object)
    warning("Prefixed feature row names with assay names because they ",
            "were not globally unique across the QFeatures object.",
            call. = FALSE)
    object
}


#' @importFrom methods is
#' @importClassesFrom S4Vectors List
#' @importFrom stats setNames
#' @noRd
rename_assay_link_features <- function(assay_links, old_names, new_names) {
    name_map <- mapply(setNames, new_names, old_names, SIMPLIFY = FALSE)
    for (name in names(assay_links)) {
        link <- assay_links[[name]]
        if (is(link@hits, "List")) {
            hits <- as.list(link@hits)
            hits <- mapply(function(hit, from) {
                rename_hit_features(hit, from, link@name, name_map)
            }, hits, link@from, SIMPLIFY = FALSE)
            link@hits <- S4Vectors::List(hits)
            names(link@hits) <- names(hits)
        } else {
            link@hits <- rename_hit_features(link@hits, link@from,
                                             link@name, name_map)
        }
        assay_links@listData[[name]] <- link
    }
    assay_links
}


#' @importFrom S4Vectors "mcols<-" mcols
#' @noRd
rename_hit_features <- function(hit, from, to, name_map) {
    md <- mcols(hit)
    if ("names_from" %in% colnames(md) && length(from) == 1L && !is.na(from))
        md$names_from <- unname(name_map[[from]][md$names_from])
    if ("names_to" %in% colnames(md) && length(to) == 1L && !is.na(to))
        md$names_to <- unname(name_map[[to]][md$names_to])
    mcols(hit) <- md
    hit
}
