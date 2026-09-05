# mulinkr

## Install

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("UCLouvain-CBIO/mulinkr")
library(mulinkr)
```

## Main functions

- `prepareQFeatures()`: Prepares a `QFeatures` object for export by making
  feature row names unique across assays, converting dates to strings, and
  simplifying `scpModel` metadata to lists.
- `writeLinkH5MU()`: Writes a prepared `QFeatures` object and its assay links
  to an `.h5mu` file.
- `readLinkH5MU()`: Reads an `.h5mu` file into a `QFeatures` object,
  restoring assay links when present.
