# Code from:
# https://github.com/lucas-diedrich/ms-scverse-bioconductor/blob/main/R/io.R

## Internals of MuData that a patched readH5MU() would call directly.
.write_matrix <- MuData:::write_matrix
.write_data_frame <- MuData:::write_data_frame


#' Global feature axis of a `QFeatures` object.
#'
#' `MuData::writeH5MU()` builds the global `.var` index by `rbind`-ing one
#' zero-column data.frame per modality in assay order, so the axis is the assays'
#' row names concatenated in that order. This is the axis a `.varp` matrix has to
#' be aligned to.
#'
#' QFeatures enforces row-name uniqueness only *within* an assay
#' (`.unique_row_names`), while mudata keys the global axis by name and mulink
#' queries resolve it with `var_names.get_indexer()`, which fails on duplicates.
#' Colliding names are therefore prefixed with the assay name and the original ID
#' is kept in `.var` (see ROADMAP.md). Prefixing is all-or-nothing rather than
#' per-feature, so that the axis stays predictable.
#'
#' @param object A `QFeatures` object.
#' @param prefix When to prefix feature names with the assay name: `"collision"`
#'     only if the plain names are not globally unique, `"always"`, or `"never"`.
#' @param sep Separator between the assay name and the feature name.
#' @param assays The assays that become modalities, in file order.
#'
#' @return A list of the global feature `key`s and, parallel to them, the `assay`
#'     and original `id` each came from.
#' @noRd
feature_index <- function(object,
                          prefix = c("collision", "always", "never"),
                          sep = ":",
                          assays = names(object)) {
    prefix <- match.arg(prefix)

    ids <- lapply(assays, function(nm) rownames(object[[nm]]))
    names(ids) <- assays
    assay <- rep(names(ids), lengths(ids))
    ids <- unlist(ids, use.names = FALSE)

    key <- switch(prefix,
                  always = paste0(assay, sep, ids),
                  never = ids,
                  collision = if (anyDuplicated(ids)) paste0(assay, sep, ids) else ids)

    if (anyDuplicated(key))
        stop("Global feature names are not unique, even after prefixing with ",
             "the assay name. Duplicates: '",
             paste(unique(key[duplicated(key)]), collapse = "', '"), "'.")

    list(key = key, assay = assay, id = ids)
}


#' Build the feature-mapping matrix from a `QFeatures` object's `AssayLinks`.
#'
#' Inverse of `assay_links_from_feature_mapping()`. Every `Hits` object in
#' `AssayLinks` is already a sparse edge list, so the only work is relabelling its
#' endpoints from per-assay row names to positions on the global feature axis.
#'
#' Edges are directed parent -> child: `A[u, v] != 0` means the feature `v` was
#' derived from the feature `u`. Only direct edges are stored, never the
#' transitive closure — mulink traverses multiple hops at query time, and a
#' closed graph would not survive the round trip, since transitive reduction is
#' not unique (see rfc.md and ROADMAP.md).
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

    ## Assay name -> global positions, named by that assay's row names, so that a
    ## feature name resolves to a row/column of the matrix by lookup.
    positions <- split(seq_len(p), factor(index$assay, levels = assays))
    positions <- mapply(function(pos, nm) setNames(pos, rownames(object[[nm]])),
                        positions, names(positions), SIMPLIFY = FALSE)

    edges <- lapply(assays, function(to) {
        al <- assayLink(object, to)

        ## @hits is either a single Hits or a List of Hits parallel to @from.
        ## Both shapes occur even with a single parent: the adjacency-matrix
        ## branch of QFeatures:::.create_assay_link() wraps its result with
        ## `if (length(hits) > 1)` while `hits` is still a Hits object, so the
        ## test counts hits rather than parents and yields a length-1 HitsList.
        hits <- if (is(al@hits, "List")) as.list(al@hits) else list(al@hits)

        do.call(rbind, lapply(seq_along(hits), function(k) {
            from <- al@from[k]
            h <- hits[[k]]
            if (is.na(from) || !(from %in% assays) || length(h) == 0L)
                return(NULL)  # root assay, or a collapsed block
            md <- mcols(h)

            ## names_from/names_to are the authoritative endpoints.
            ## QFeatures:::.pruneHits() drops hits by row name when the object is
            ## subset but never renumbers the integer from/to indices, nor
            ## nLnode/nRnode, so those are stale on any subset object.
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
#' `MuData:::write_matrix()` routes an integer or logical column that contains NA
#' through the `nullable-integer`/`nullable-boolean` encoding, whose mask rhdf5
#' writes as `int8`. anndata requires a genuine boolean mask and rejects the file
#' with "mask should be boolean numpy array", so a single NA in such a column
#' makes the whole .h5mu unreadable from Python.
#'
#' Doubles avoid that encoding entirely -- `write_matrix()` turns their NAs into
#' NaN, which both ecosystems read as missing -- so the affected columns are cast
#' to double. Values and missingness survive; the integer/logical type does not.
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
#' `MuData::writeH5MU()` builds the global index as
#' `do.call(rbind, vars)` over one *zero-column* data.frame per modality
#' (`write_h5mu.R:184`). `rbind` on zero-column data.frames returns zero rows, so
#' the group it writes is always empty, whatever the object contained. Verified
#' against MuData 1.14.0.
#'
#' The index is not optional here: mudata keys the global axis by name, and both
#' `read_varp()` and mulink resolve `.varp` against it. So it is written again,
#' correctly, after `writeH5MU()` has returned.
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
#' Counterpart of `read_varp()`. `MuData::writeH5MU()` writes no global `.varp`,
#' so the group is added afterwards. By then the file is closed and finalized
#' (`finalize_mudata()`); reopening it read/write preserves the 512-byte user
#' block that holds the MuData signature.
#'
#' `write_matrix()` stores a `dgRMatrix`'s row-compressed arrays under the
#' `csr_matrix` encoding, which is exactly how Python reads them back, so
#' `mdata.varp[key][u, v]` in Python is `mat[u, v]` here. `read_varp()` inverts
#' this: `read_sparse_matrix()` reads a `csr_matrix` as its transpose, and
#' `read_varp()` transposes again.
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
    ## after = FALSE so the group is closed before the file.
    on.exit(H5Gclose(grp), add = TRUE, after = FALSE)

    .write_matrix(grp, key, as(mat, "RsparseMatrix"), needTranspose = TRUE)
    invisible(NULL)
}


#' Write a `QFeatures` object to an .h5mu file.
#'
#' Inverse of `readLinkH5MU()`. `MuData::writeH5MU()` covers the
#' MultiAssayExperiment skeleton; the feature graph held in `AssayLinks` has no
#' place there, so it is written separately to the global `.varp` as a `p x p`
#' sparse adjacency matrix over the global `.var` index.
#'
#'
#' @param object A `QFeatures` object.
#' @param path Path of the .h5mu file to create.
#' @param feature_mapping_key Key to store the feature graph under in `.varp`.
#' @param prefix,sep Passed to `feature_index()`.
#' @param overwrite Whether to replace `path` if it already exists.
#'
#' @return `path`, invisibly.
#' @importClassesFrom Matrix Matrix
#' @importClassesFrom QFeatures QFeatures
#' @importFrom methods is
#' @importFrom MuData writeH5MU
#' @importFrom MultiAssayExperiment ExperimentList experiments
#' @importFrom rhdf5 H5Fopen H5Fclose
#' @importFrom SummarizedExperiment colData "colData<-" rowData "rowData<-"
#' @export
writeLinkH5MU <- function(object, path,
                               feature_mapping_key = "feature_mapping",
                               prefix = c("collision", "always", "never"),
                               sep = ":",
                               overwrite = FALSE) {
    if (!is(object, "QFeatures"))
        stop("'object' must be a QFeatures object.")
    if (file.exists(path)) {
        ## MuData::writeH5MU() takes an `overwrite` argument but never reads it,
        ## and its H5Fcreate() call defaults to H5F_ACC_TRUNC, so it destroys an
        ## existing file whatever is passed. Guard explicitly. Verified against
        ## MuData 1.14.0.
        if (!overwrite)
            stop("'", path, "' already exists; pass overwrite = TRUE to replace it.")
        unlink(path)
    }

    index <- feature_index(object, prefix = prefix, sep = sep,
                           assays = names(object))
    ## Built before any renaming, since the AssayLinks are keyed on the row names
    ## the object currently carries.
    feature_mapping <- feature_mapping_from_assay_links(object, index,
                                                        assays = names(object))

    keys <- split(index$key, factor(index$assay, levels = names(object)))
    dropped <- character(0)
    cast <- character(0)

    prepared <- mapply(function(se, name, key) {
        ## MuData cannot encode a matrix-valued rowData column: writeH5AD()
        ## stores it as a nested sparse-matrix group inside `var`, which
        ## read_dataframe() then cannot coerce back to a DataFrame column, so the
        ## file it produces is unreadable. Such a column is an adjacency matrix
        ## (see QFeatures:::.create_assay_link()) -- the same graph being written
        ## to .varp -- so it is dropped rather than left to break the file.
        rd <- rowData(se)
        is_matrix <- vapply(rd, function(x) is.matrix(x) || is(x, "Matrix"),
                            logical(1))
        if (any(is_matrix)) {
            dropped <<- c(dropped, paste0(name, "$", colnames(rd)[is_matrix]))
            rowData(se) <- rd[, !is_matrix, drop = FALSE]
        }
        if (!identical(rownames(se), key)) {
            rowData(se)$mulink_feature_id <- rownames(se)
            rownames(se) <- key
        }

        fixed <- cast_nullable_columns(rowData(se), paste0(name, "/rowData"))
        rowData(se) <- fixed$df
        cast <<- c(cast, fixed$cast)

        fixed <- cast_nullable_columns(colData(se), paste0(name, "/colData"))
        colData(se) <- fixed$df
        cast <<- c(cast, fixed$cast)

        se
    }, experiments(object), names(object), keys, SIMPLIFY = FALSE)

    ## The global colData is written by writeH5MU() as /obs and is just as
    ## exposed to the nullable encoding as the per-assay annotations.
    fixed <- cast_nullable_columns(colData(object), "colData")
    object@colData <- fixed$df
    cast <- c(cast, fixed$cast)

    if (length(cast))
        warning("Cast NA-bearing integer/logical column(s) '",
                paste(cast, collapse = "', '"), "' to double, which MuData ",
                "encodes in a form Python can read. See cast_nullable_columns().")

    if (length(dropped))
        warning("Dropped matrix-valued rowData column(s) '",
                paste(dropped, collapse = "', '"), "', which MuData cannot ",
                "encode. The links they define are written to .varp.")

    ## Renaming and dropping invalidate the AssayLinks, whose hits still
    ## reference the old row names. writeH5MU() does not read them and the graph
    ## has already been extracted, so the slot is assigned directly to skip the
    ## validity check that would reject this write-only object.
    object@ExperimentList <- ExperimentList(prepared)

    ## would simply never be looked at; they are pruned anyway to keep the
    ## write-only object internally consistent.

    MuData::writeH5MU(object, path)
    write_var_index(path, index$key)

    ## Read the index back rather than assuming it, so that .varp is aligned to
    ## what the file actually contains.
    h5 <- H5Fopen(path, flags = "H5F_ACC_RDONLY", native = FALSE)
    var_names <- .var_names(h5)
    H5Fclose(h5)

    if (!setequal(var_names, index$key))
        stop("The global .var index in '", path, "' does not match the ",
             "expected feature axis.")

    write_varp(path, feature_mapping_key,
               feature_mapping[var_names, var_names, drop = FALSE])
    invisible(path)
}
