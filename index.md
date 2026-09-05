# mulinkr

## Install

``` r

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("UCLouvain-CBIO/mulinkr")
library(mulinkr)
```

## Main functions

- [`prepareQFeatures()`](https://uclouvain-cbio.github.io/mulinkr/reference/prepareQFeatures.md):
  Prepares a `QFeatures` object for export by making feature row names
  unique across assays, converting dates to strings, and simplifying
  `scpModel` metadata to lists.
- [`writeLinkH5MU()`](https://uclouvain-cbio.github.io/mulinkr/reference/writeLinkH5MU.md):
  Writes a prepared `QFeatures` object and its assay links to an `.h5mu`
  file.
- [`readLinkH5MU()`](https://uclouvain-cbio.github.io/mulinkr/reference/readLinkH5MU.md):
  Reads an `.h5mu` file into a `QFeatures` object, restoring assay links
  when present.
