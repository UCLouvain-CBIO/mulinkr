.write_matrix <- MuData:::write_matrix
.write_data_frame <- MuData:::write_data_frame

#' Global feature axis of a `QFeatures` object.
#'
#' @param object A `QFeatures` object.
#' @param assays The assays that become modalities, in file order.
#'
#' @return A list of the global feature `key`s and, parallel to them, the `assay`
#'     each came from.
#' @noRd
feature_index <- function(object, assays = names(object)) {
    ids <- lapply(assays, function(nm) rownames(object[[nm]]))
    names(ids) <- assays
    assay <- rep(names(ids), lengths(ids))
    key <- unlist(ids, use.names = FALSE)

    if (anyDuplicated(key))
        stop("Global feature names are not unique. Run prepareQFeatures() ",
             "before writeLinkH5MU(). Duplicates: '",
             paste(unique(key[duplicated(key)]), collapse = "', '"), "'.")

    list(key = key, assay = assay)
}


#' Build the feature-mapping matrix from a `QFeatures` object's `AssayLinks`.
#'
#' Inverse of `assay_links_from_feature_mapping()`. Every `Hits` object in
#' `AssayLinks` is already a sparse edge list, so the only work is relabelling its
#' endpoints from per-assay row names to positions on the global feature axis.
#'
#' Edges are directed parent -> child: `A[u, v] != 0` means the feature `v` was
#' derived from the feature `u`. Only direct edges are stored. Hit metadata is
#' used for feature names because integer hit positions can become stale after
#' subsetting a `QFeatures` object.
#'
#' @param object A `QFeatures` object.
#' @param index The global feature axis, from `feature_index()`.
#' @param assays The assays that become modalities, in file order. 
#'
#' @return A `p x p` sparse matrix over the global feature axis.
#' @importClassesFrom S4Vectors List
#' @importFrom methods is
#' @importFrom Matrix sparseMatrix
#' @importFrom QFeatures assayLink
#' @importFrom S4Vectors mcols
#' @importFrom stats setNames
#' @noRd
feature_mapping_from_assay_links <- function(object,
                                             index = feature_index(object),
                                             assays = names(object)) {
    p <- length(index$key)

    positions <- split(seq_len(p), factor(index$assay, levels = assays))
    positions <- mapply(function(pos, nm) setNames(pos, rownames(object[[nm]])),
                        positions, names(positions), SIMPLIFY = FALSE)

    edges <- lapply(assays, function(to) {
        al <- assayLink(object, to)
        hits <- if (is(al@hits, "List")) as.list(al@hits) else list(al@hits)

        do.call(rbind, lapply(seq_along(hits), function(k) {
            from <- al@from[k]
            h <- hits[[k]]
            if (is.na(from) || !(from %in% assays) || length(h) == 0L)
                return(NULL)
            md <- mcols(h)

            data.frame(i = positions[[from]][md$names_from],
                       j = positions[[to]][md$names_to],
                       x = if ("weigths" %in% colnames(md)) as.numeric(md$weigths)
                           else rep(1, length(h)))
        }))
    })

    edges <- do.call(rbind, edges)
    if (is.null(edges))
        edges <- data.frame(i = integer(0), j = integer(0), x = numeric(0))

    sparseMatrix(i = edges$i, j = edges$j, x = edges$x,
                 dims = c(p, p), dimnames = list(index$key, index$key),
                 repr = "R")
}


#' Cast NA-bearing integer and logical columns to double.
#'
#' Integer and logical columns with missing values are widened so MuData writes
#' them in a Python-readable form. Values and missingness are preserved, but the
#' original integer/logical type is not.
#'
#' @param df A `DataFrame`.
#' @param context Label used to name the cast columns in the warning.
#'
#' @return A list of the possibly modified `df` and the names of the `cast`
#'     columns.
#' @noRd
cast_nullable_columns <- function(df, context) {
    affected <- vapply(df, function(x) (is.integer(x) || is.logical(x)) && anyNA(x),
                       logical(1))
    if (!any(affected))
        return(list(df = df, cast = character(0)))

    df[affected] <- lapply(df[affected], as.double)
    list(df = df, cast = paste0(context, "$", colnames(df)[affected]))
}


#' Rewrite the global `.var` index of an .h5mu file.
#'
#' MuData can leave the global feature index empty when modalities have
#' zero-column `rowData`. This rewrites `/var` after `writeH5MU()` so global
#' `.varp` matrices have a named feature axis.
#'
#' @param file Path to the .h5mu file.
#' @param keys The global feature names, in file order.
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists H5Ldelete
#' @noRd
write_var_index <- function(file, keys) {
    h5 <- H5Fopen(file, flags = "H5F_ACC_RDWR", native = FALSE)
    on.exit(H5Fclose(h5), add = TRUE)

    if (H5Lexists(h5, "var")) H5Ldelete(h5, "var")
    .write_data_frame(h5, "var", data.frame(row.names = keys))
    invisible(NULL)
}


#' Append a matrix to the global `.varp` group of an .h5mu file.
#'
#' `MuData::writeH5MU()` does not write global `.varp`, so this adds or replaces
#' one sparse pairwise feature matrix after the base file has been written.
#'
#' @param file Path to the .h5mu file.
#' @param key Name to store the matrix under in `.varp`.
#' @param mat A sparse matrix over the global feature axis.
#' @importClassesFrom Matrix RsparseMatrix
#' @importFrom methods as
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists H5Ldelete H5Gopen H5Gcreate H5Gclose
#' @noRd
write_varp <- function(file, key, mat) {
    h5 <- H5Fopen(file, flags = "H5F_ACC_RDWR", native = FALSE)
    on.exit(H5Fclose(h5), add = TRUE)

    if (H5Lexists(h5, "varp")) {
        grp <- H5Gopen(h5, "varp")
        if (H5Lexists(grp, key)) H5Ldelete(grp, key)
    } else {
        grp <- H5Gcreate(h5, "varp")
    }
    on.exit(H5Gclose(grp), add = TRUE, after = FALSE)

    .write_matrix(grp, key, as(mat, "RsparseMatrix"), needTranspose = TRUE)
    invisible(NULL)
}


#' Prepare an output path for MuData writing.
#'
#' `MuData::writeH5MU()` currently ignores its `overwrite` argument, so this
#' helper enforces the requested behavior before the file is opened.
#'
#' @param path Path of the .h5mu file to create.
#' @param overwrite Whether to replace `path` if it already exists.
#'
#' @return `NULL`, invisibly.
#' @noRd
prepare_output_path <- function(path, overwrite) {
    if (!file.exists(path))
        return(invisible(NULL))

    if (!overwrite)
        stop("'", path, "' already exists; pass overwrite = TRUE to replace it.")

    unlink(path)
    invisible(NULL)
}


#' Clean one assay before writing with MuData.
#'
#' Matrix-valued `rowData` columns are removed because MuData cannot round-trip
#' them as data-frame columns. NA-bearing integer/logical annotation columns are
#' cast with `cast_nullable_columns()`.
#'
#' @param se A `SummarizedExperiment`.
#' @param name Assay name, used in warning labels.
#'
#' @return A list containing the cleaned assay, dropped columns, and cast columns.
#' @importClassesFrom Matrix Matrix
#' @importFrom methods is
#' @importFrom SummarizedExperiment colData "colData<-" rowData "rowData<-"
#' @noRd
clean_experiment_for_h5mu <- function(se, name) {
    dropped <- character(0)
    cast <- character(0)

    rd <- rowData(se)
    is_matrix <- vapply(rd, function(x) is.matrix(x) || is(x, "Matrix"),
                        logical(1))
    if (any(is_matrix)) {
        dropped <- paste0(name, "$", colnames(rd)[is_matrix])
        rowData(se) <- rd[, !is_matrix, drop = FALSE]
    }

    fixed <- cast_nullable_columns(rowData(se), paste0(name, "/rowData"))
    rowData(se) <- fixed$df
    cast <- c(cast, fixed$cast)

    fixed <- cast_nullable_columns(colData(se), paste0(name, "/colData"))
    colData(se) <- fixed$df
    cast <- c(cast, fixed$cast)

    list(experiment = se, dropped = dropped, cast = cast)
}


#' Clean a `QFeatures` object before writing with MuData.
#'
#' Applies the per-assay cleanup and the same nullable-column cast to global
#' `colData`, which MuData writes as `/obs`.
#'
#' @param object A `QFeatures` object.
#'
#' @return A list containing the cleaned object, dropped columns, and cast columns.
#' @importFrom MultiAssayExperiment ExperimentList experiments
#' @importFrom SummarizedExperiment colData
#' @noRd
clean_qfeatures_for_h5mu <- function(object) {
    assay_names <- names(object)
    cleaned <- mapply(clean_experiment_for_h5mu,
                      experiments(object),
                      assay_names,
                      SIMPLIFY = FALSE)

    prepared <- lapply(cleaned, `[[`, "experiment")
    names(prepared) <- assay_names
    dropped <- unlist(lapply(cleaned, `[[`, "dropped"), use.names = FALSE)
    cast <- unlist(lapply(cleaned, `[[`, "cast"), use.names = FALSE)

    fixed <- cast_nullable_columns(colData(object), "colData")
    object@colData <- fixed$df
    object@ExperimentList <- ExperimentList(prepared)

    list(object = object,
         dropped = dropped,
         cast = c(cast, fixed$cast))
}


#' Warn about H5MU write-time cleanup.
#'
#' Reports type casts and dropped row annotations performed to keep the resulting
#' file readable across the R and Python MuData stacks.
#'
#' @param cast Labels for columns cast to double.
#' @param dropped Labels for matrix-valued `rowData` columns that were removed.
#'
#' @return `NULL`, invisibly.
#' @noRd
warn_h5mu_cleanup <- function(cast, dropped) {
    if (length(cast))
        warning("Cast NA-bearing integer/logical column(s) '",
                paste(cast, collapse = "', '"), "' to double, which MuData ",
                "encodes in a form Python can read. See cast_nullable_columns().")

    if (length(dropped))
        warning("Dropped matrix-valued rowData column(s) '",
                paste(dropped, collapse = "', '"), "', which MuData cannot ",
                "encode. The links they define are written to .varp.")

    invisible(NULL)
}


#' Read and validate the global `.var` index of an .h5mu file.
#'
#' The written feature names are read back so `.varp` can be ordered against the
#' file's actual global feature axis.
#'
#' @param path Path of the .h5mu file.
#' @param expected Expected feature names.
#'
#' @return Feature names from the file, in file order.
#' @importFrom rhdf5 H5Fopen H5Fclose
#' @noRd
read_var_names_checked <- function(path, expected) {
    h5 <- H5Fopen(path, flags = "H5F_ACC_RDONLY", native = FALSE)
    on.exit(H5Fclose(h5), add = TRUE)

    var_names <- .var_names(h5)
    if (!setequal(var_names, expected))
        stop("The global .var index in '", path, "' does not match the ",
             "expected feature axis.")

    var_names
}


#' Write a `QFeatures` object to an .h5mu file.
#'
#' Inverse of `readLinkH5MU()`. `MuData::writeH5MU()` covers the
#' MultiAssayExperiment skeleton; the feature graph held in `AssayLinks` is
#' written separately to the global `.varp` as a `p x p` sparse adjacency
#' matrix over the global `.var` index.
#'
#'
#' @param object A `QFeatures` object prepared with `prepareQFeatures()`.
#' @param path Path of the .h5mu file to create.
#' @param feature_mapping_key Key to store the feature graph under in `.varp`.
#' @param overwrite Whether to replace `path` if it already exists.
#'
#' @return `path`, invisibly.
#' @importClassesFrom Matrix Matrix
#' @importClassesFrom QFeatures QFeatures
#' @importFrom methods is
#' @importFrom MuData writeH5MU
#' @importFrom rhdf5 H5Iis_valid H5Iget_type
#' @export
writeLinkH5MU <- function(object, path,
                               feature_mapping_key = "feature_mapping",
                               overwrite = FALSE) {
    if (!is(object, "QFeatures"))
        stop("'object' must be a QFeatures object.")

    assay_names <- names(object)
    prepare_output_path(path, overwrite)

    index <- feature_index(object, assays = assay_names)
    feature_mapping <- feature_mapping_from_assay_links(object, index,
                                                        assays = assay_names)

    cleaned <- clean_qfeatures_for_h5mu(object)
    object <- cleaned$object
    warn_h5mu_cleanup(cleaned$cast, cleaned$dropped)
    MuData::writeH5MU(object, path)
    write_var_index(path, index$key)

    var_names <- read_var_names_checked(path, index$key)
    write_varp(path, feature_mapping_key,
               feature_mapping[var_names, var_names, drop = FALSE])
    invisible(path)
}
