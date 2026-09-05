# Write a `QFeatures` object to an .h5mu file.

Inverse of
[`readLinkH5MU()`](https://uclouvain-cbio.github.io/mulinkr/reference/readLinkH5MU.md).
[`MuData::writeH5MU()`](https://rdrr.io/pkg/MuData/man/writeH5MU.html)
covers the MultiAssayExperiment skeleton; the feature graph held in
`AssayLinks` is written separately to the global `.varp` as a `p x p`
sparse adjacency matrix over the global `.var` index.

## Usage

``` r
writeLinkH5MU(
  object,
  path,
  feature_mapping_key = "feature_mapping",
  overwrite = FALSE
)
```

## Arguments

- object:

  A `QFeatures` object prepared with
  [`prepareQFeatures()`](https://uclouvain-cbio.github.io/mulinkr/reference/prepareQFeatures.md).

- path:

  Path of the .h5mu file to create.

- feature_mapping_key:

  Key to store the feature graph under in `.varp`.

- overwrite:

  Whether to replace `path` if it already exists.

## Value

`path`, invisibly.

## Limitations

Writing to .h5mu does not preserve all R data types and structures:

- Date and date-time classes (`Date`, `POSIXct`, and `POSIXlt`) are not
  supported natively and must be converted before writing.

- The distinction between numeric `NA` and `NaN` is lost: both are read
  back as `NaN` by
  [`readLinkH5MU()`](https://uclouvain-cbio.github.io/mulinkr/reference/readLinkH5MU.md).

- Complex data structures in `rowData` or `metadata` are not fully
  supported. Matrix-valued `rowData` columns are dropped with a warning.
  Unsupported metadata objects, such as `scpModel` objects (class
  `ScpModel`), must be simplified before writing.

- Unused factor levels are not preserved.

- Feature row names must be globally unique across all sets. Duplicate
  row names cause an error.

- Integer and logical columns containing `NA` in global `colData`,
  set-specific `colData`, or `rowData` are converted to doubles with a
  warning. Their original types are lost; logical values become `1` and
  `0`, while missing values remain missing.

Use `object <- prepareQFeatures(object)` before writing to address some
of these limitations.
[`prepareQFeatures()`](https://uclouvain-cbio.github.io/mulinkr/reference/prepareQFeatures.md)
converts date and date-time columns in global `colData` and set-specific
`rowData` and `colData` to formatted character strings, converts
`scpModel` objects in set metadata to lists, and prefixes feature row
names with set names when needed to make them globally unique, updating
the set links accordingly. It does not restore the original date
classes, missing-value distinction, unused factor levels, or
integer/logical column types when the file is read back.

## Examples

``` r
data("feat3", package = "QFeatures")
preparedFeat3 <- prepareQFeatures(feat3)
#> Warning: Prefixed feature row names with assay names because they were not globally unique across the QFeatures object.
filePath <- tempfile(fileext = ".h5mu")
writeLinkH5MU(preparedFeat3, filePath)
file.exists(filePath)
#> [1] TRUE
unlink(filePath)
```
