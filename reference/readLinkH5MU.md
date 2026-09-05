# Read an .h5mu file and create a `QFeatures` object.

Read an .h5mu file and create a `QFeatures` object.

## Usage

``` r
readLinkH5MU(path, feature_mapping_key = "feature_mapping", backed = FALSE)
```

## Arguments

- path:

  Path to the .h5mu file.

- feature_mapping_key:

  Key of the feature graph in the global `.varp`.

- backed:

  Passed to
  [`MuData::readH5MU()`](https://rdrr.io/pkg/MuData/man/readH5MU.html).

## Value

A `QFeatures` object. Assays and their links come from the file; if the
`.varp` key is absent the sets are returned unlinked.

## Limitations

Missing numeric values written as `NA` or `NaN` are read as `NaN`. The
original distinction between `NA` and `NaN` is not preserved. See
[`writeLinkH5MU()`](https://uclouvain-cbio.github.io/mulinkr/reference/writeLinkH5MU.md)
for other limitations of the write/read conversion.

## Examples

``` r
data("feat3", package = "QFeatures")
preparedFeat3 <- prepareQFeatures(feat3)
#> Warning: Prefixed feature row names with assay names because they were not globally unique across the QFeatures object.
filePath <- tempfile(fileext = ".h5mu")
writeLinkH5MU(preparedFeat3, filePath)
newQFeatures <- readLinkH5MU(filePath)
#> Warning: sampleMap[['assay']] coerced with as.factor()
newQFeatures
#> An instance of class QFeatures (type: bulk) with 7 sets:
#> 
#>  [1] psms1: SummarizedExperiment with 7 rows and 2 columns 
#>  [2] psms2: SummarizedExperiment with 8 rows and 2 columns 
#>  [3] psmsall: SummarizedExperiment with 10 rows and 4 columns 
#>  [4] peptides: SummarizedExperiment with 3 rows and 4 columns 
#>  [5] proteins: SummarizedExperiment with 2 rows and 4 columns 
#>  [6] normpeptides: SummarizedExperiment with 3 rows and 4 columns 
#>  [7] normproteins: SummarizedExperiment with 2 rows and 4 columns 
unlink(filePath)
```
