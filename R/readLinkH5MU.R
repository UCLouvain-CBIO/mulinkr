# Code from:
# https://github.com/lucas-diedrich/ms-scverse-bioconductor/blob/main/R/io.R

## Internals of MuData that a patched readH5MU() would call directly.
.h5autoclose <- MuData:::h5autoclose
.read_matrix <- MuData:::read_matrix
.read_with_index <- MuData:::read_with_index


#' Global feature names of an .h5mu file, in file order.
#' @noRd
.var_names <- function(h5) {
    rownames(.read_with_index(.h5autoclose(h5 & "var")))
}


#' Read the global `.varp` group of an .h5mu file.
#'
#' Mirrors the `varp` branch of `MuData:::read_modality()`, which reads a
#' modality's pairwise feature matrices into `rowPair()`. Here the matrices are
#' returned instead, since they span all modalities and belong to no single one.
#'
#' `read_matrix()` transposes to maintain directionality as observation/feature directionality is inversed between
#' scverse and Bioconductor functions
#'
#' @param file Path to the .h5mu file.
#' @param keys Names to read from `.varp`. All of them if `NULL`.
#'
#' @return A named list of sparse matrices, with the global feature names as
#'     dimnames.
#' @importFrom methods is
#' @importClassesFrom Matrix dsparseMatrix
#' @importMethodsFrom Matrix t
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists h5ls
#' @importFrom stats setNames
#' @noRd
read_varp <- function(file, keys = NULL) {
    h5 <- H5Fopen(file, flags = "H5F_ACC_RDONLY", native = FALSE)
    on.exit(H5Fclose(h5), add = TRUE)

    if (!H5Lexists(h5, "varp")) return(list())

    var_names <- .var_names(h5)
    available <- h5ls(.h5autoclose(h5 & "varp"), recursive = FALSE)$name
    if (is.null(keys)) keys <- available

    missing_keys <- setdiff(keys, available)
    if (length(missing_keys))
        warning("No '", paste(missing_keys, collapse = "', '"), "' in .varp.")

    matrices <- lapply(setNames(nm = intersect(keys, available)), function(key) {
        m <- .read_matrix(.h5autoclose(h5 & paste("varp", key, sep = "/")))
        if (!is(m, "dsparseMatrix")) {
            warning("Pairwise varp matrix ", key, " is not a sparse matrix. ",
                    "Only sparse matrices are currently supported, skipping...")
            return(NULL)
        }
        m <- t(m)
        dimnames(m) <- list(var_names, var_names)
        m
    })

    matrices[!vapply(matrices, is.null, logical(1))]
}


#' Build a `Hits` object from an adjacency submatrix.
#'
#' Rows are parent features, columns are child features. `names_from`/`names_to`
#' are required by QFeatures' validity checks (`.checkLinksInHits`), which verify
#' that every linked feature exists in the corresponding assay.
#'
#' `sort.by.query = TRUE` yields a `SortedByQueryHits`, which is what
#' `findMatches()` and therefore the rest of QFeatures produces.
#' @importFrom methods as
#' @importClassesFrom Matrix TsparseMatrix
#' @importFrom S4Vectors Hits mcols
#' @noRd
hits_from_adjacency <- function(adj) {
    nz <- Matrix::summary(as(adj, "TsparseMatrix"))
    nz <- nz[order(nz$i, nz$j), , drop = FALSE]

    hits <- S4Vectors::Hits(from = nz$i,
                            to = nz$j,
                            nLnode = nrow(adj),
                            nRnode = ncol(adj),
                            sort.by.query = TRUE)

    S4Vectors::mcols(hits)$names_from <- rownames(adj)[nz$i]
    S4Vectors::mcols(hits)$names_to <- colnames(adj)[nz$j]
    hits
}


#' Turn a feature-mapping matrix into one `AssayLink` per assay.
#'
#' An edge u -> v in `.varp` means v's assay was derived from u's, so edges are
#' grouped by their *target* modality: that is the child, and the modalities the
#' edges come from are its parents. Grouping matters — `AssayLinks` is keyed by
#' child name, so adding links pair-by-pair would keep only the last parent of a
#' fan-in.
#' @importFrom QFeatures AssayLink AssayLinks
#' @importFrom stats setNames
#' @noRd
assay_links_from_feature_mapping <- function(experiments,
                                             feature_mapping,
                                             fcol = NA_character_) {
    mod_names <- names(experiments)
    var_names <- rownames(feature_mapping)

    ## The mulink convention requires globally unique feature names, so matching
    ## each modality's rownames against the global index is unambiguous.
    positions <- lapply(experiments, function(se) match(rownames(se), var_names))
    unmatched <- vapply(positions, anyNA, logical(1))
    if (any(unmatched))
        stop("Features of modality/modalities '",
             paste(mod_names[unmatched], collapse = "', '"),
             "' are absent from the global .var index.")

    has_edges <- function(adj) length(adj@x) > 0L && any(adj@x != 0)

    edges <- list()
    for (from in mod_names) {
        for (to in setdiff(mod_names, from)) {
            adj <- feature_mapping[positions[[from]], positions[[to]], drop = FALSE]
            if (!has_edges(adj)) next

            dimnames(adj) <- list(rownames(experiments[[from]]),
                                  rownames(experiments[[to]]))
            edges[[to]] <- c(edges[[to]], setNames(list(adj), from))
        }
    }

    within_modality <- vapply(mod_names, function(m) {
        has_edges(feature_mapping[positions[[m]], positions[[m]], drop = FALSE])
    }, logical(1))
    if (any(within_modality))
        warning("Ignoring edges within modality/modalities '",
                paste(mod_names[within_modality], collapse = "', '"),
                "'; an assay cannot be linked to itself.")

    ## One AssayLink per assay, in assay order: QFeatures' validity checks
    ## compare names(assayLinks) to names(object) with identical().
    AssayLinks(lapply(mod_names, function(to) {
        parents <- edges[[to]]
        if (is.null(parents)) return(AssayLink(name = to))

        hits <- lapply(parents, hits_from_adjacency)
        if (length(hits) > 1L) {
            hits <- S4Vectors::List(hits)
            names(hits) <- names(parents)
        } else {
            hits <- hits[[1L]]
        }

        AssayLink(name = to,
                  from = names(parents),
                  fcol = rep(fcol, length(parents)),
                  hits = hits)
    }))
}


#' Read an .h5mu file and create a `QFeatures` object.
#'
#' @section Limitations:
#' Missing numeric values written as `NA` or `NaN` are read as `NaN`.
#' The original distinction between `NA` and `NaN` is not preserved.
#' See [writeLinkH5MU()] for other limitations of the write/read conversion.
#'
#' @param path Path to the .h5mu file.
#' @param feature_mapping_key Key of the feature graph in the global `.varp`.
#' @param backed Passed to `MuData::readH5MU()`.
#'
#' @return A `QFeatures` object. Assays and their links come from the file; if
#'     the `.varp` key is absent the sets are returned unlinked.
#' @importClassesFrom QFeatures QFeatures
#' @importFrom methods new validObject
#' @importFrom MultiAssayExperiment experiments sampleMap
#' @importFrom MuData readH5MU
#' @importFrom QFeatures AssayLinks
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment colData
#' @examples
#' data("feat3", package = "QFeatures")
#' preparedFeat3 <- prepareQFeatures(feat3)
#' filePath <- tempfile(fileext = ".h5mu")
#' writeLinkH5MU(preparedFeat3, filePath)
#' newQFeatures <- readLinkH5MU(filePath)
#' newQFeatures
#' unlink(filePath)
#' @export
readLinkH5MU <- function(path,
                                feature_mapping_key = "feature_mapping",
                                backed = FALSE) {
    mae <- MuData::readH5MU(path, backed = backed)
    feature_mapping <- read_varp(path, keys = feature_mapping_key)[[feature_mapping_key]]

    assay_links <- if (is.null(feature_mapping)) {
        AssayLinks(names = names(experiments(mae)))
    } else {
        ## fcol normally names the rowData column a link was derived from. The
        ## link comes from .varp here, so record that key instead; the rowData
        ## column that produced it upstream lives in .uns (see ROADMAP.md).
        assay_links_from_feature_mapping(experiments(mae),
                                         feature_mapping,
                                         fcol = feature_mapping_key)
    }

    ## Promotion mirrors the QFeatures() constructor, which builds a
    ## MultiAssayExperiment and copies its slots across.
    qf <- new("QFeatures",
              ExperimentList = experiments(mae),
              colData = colData(mae),
              sampleMap = sampleMap(mae),
              metadata = metadata(mae),
              assayLinks = assay_links)

    validObject(qf)
    qf
}
