# Prepare a `QFeatures` object for .h5mu writing.

`prepareQFeatures()` makes feature row names globally unique across
assays by prefixing them with their assay name when needed. It also
change any date columns of the global `colData`, individual `rowData`
and individual `colData` to a formatted character type.

## Usage

``` r
prepareQFeatures(object, sep = ":")
```

## Arguments

- object:

  A `QFeatures` object.

- sep:

  Separator between the assay name and feature name.

## Value

A prepared `QFeatures` object.

## Examples

``` r
data("feat3", package = "QFeatures")
rownames(feat3)
#> CharacterList of length 7
#> [["psms1"]] PSM1 PSM2 PSM3 PSM4 PSM5 PSM6 PSM7
#> [["psms2"]] PSM3 PSM4 PSM5 PSM6 PSM7 PSM8 PSM9 PSM10
#> [["psmsall"]] PSM3 PSM4 PSM5 PSM6 PSM7 PSM1 PSM2 PSM8 PSM9 PSM10
#> [["peptides"]] ELGNDAYK IAEESNFPFIK SYGFNAAR
#> [["proteins"]] ProtA ProtB
#> [["normpeptides"]] ELGNDAYK IAEESNFPFIK SYGFNAAR
#> [["normproteins"]] ProtA ProtB
preparedFeat3 <- prepareQFeatures(feat3)
#> Warning: Prefixed feature row names with assay names because they were not globally unique across the QFeatures object.
rownames(preparedFeat3)
#> CharacterList of length 7
#> [["psms1"]] psms1:PSM1 psms1:PSM2 psms1:PSM3 ... psms1:PSM6 psms1:PSM7
#> [["psms2"]] psms2:PSM3 psms2:PSM4 psms2:PSM5 ... psms2:PSM9 psms2:PSM10
#> [["psmsall"]] psmsall:PSM3 psmsall:PSM4 ... psmsall:PSM9 psmsall:PSM10
#> [["peptides"]] peptides:ELGNDAYK peptides:IAEESNFPFIK peptides:SYGFNAAR
#> [["proteins"]] proteins:ProtA proteins:ProtB
#> [["normpeptides"]] normpeptides:ELGNDAYK normpeptides:IAEESNFPFIK normpeptides:SYGFNAAR
#> [["normproteins"]] normproteins:ProtA normproteins:ProtB
```
