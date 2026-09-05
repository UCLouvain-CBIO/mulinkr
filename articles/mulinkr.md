# mulinkr

The `mulinkr` package writes `QFeatures` objects and their assay links
to `.h5mu` files and reads them back into R.

Here we will use the example data `feat3` from `QFeatures`.

``` r

library(mulinkr)
data("feat3", package = "QFeatures")
feat3
#> An instance of class QFeatures (type: bulk) with 7 sets:
#> 
#>  [1] psms1: SummarizedExperiment with 7 rows and 2 columns 
#>  [2] psms2: SummarizedExperiment with 8 rows and 2 columns 
#>  [3] psmsall: SummarizedExperiment with 10 rows and 4 columns 
#>  [4] peptides: SummarizedExperiment with 3 rows and 4 columns 
#>  [5] proteins: SummarizedExperiment with 2 rows and 4 columns 
#>  [6] normpeptides: SummarizedExperiment with 3 rows and 4 columns 
#>  [7] normproteins: SummarizedExperiment with 2 rows and 4 columns
```

We first prepare the object for writing. This makes feature names unique
across assays when needed and converts date columns to character.

``` r

preparedFeat3 <- prepareQFeatures(feat3)
#> Warning: Prefixed feature row names with assay names because they were not
#> globally unique across the QFeatures object.
```

Then we write the prepared object to a temporary `.h5mu` file.

``` r

filePath <- tempfile(fileext = ".h5mu")
writeLinkH5MU(preparedFeat3, filePath)
```

The created `.h5mu` file can then be used in python or in R.

In python, the file can be read as a [`mulink`
object](https://mulink.readthedocs.io/).

``` python
import mudata as md
import mulink

mdata = md.read_h5mu("yourFile.h5mu")

print(mdata)
```

In R, the object can be reimported back as a `QFeatures` object.

``` r

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
```

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] mulinkr_0.99.1   BiocStyle_2.40.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] tidyselect_1.2.1            dplyr_1.2.1                
#>  [3] farver_2.1.2                S7_0.2.2                   
#>  [5] fastmap_1.2.0               SingleCellExperiment_1.34.0
#>  [7] lazyeval_0.2.3              nipals_1.0                 
#>  [9] digest_0.6.39               lifecycle_1.0.5            
#> [11] cluster_2.1.8.2             ProtGenerics_1.44.0        
#> [13] magrittr_2.0.5              compiler_4.6.1             
#> [15] rlang_1.3.0                 sass_0.4.10                
#> [17] tools_4.6.1                 igraph_2.3.3               
#> [19] yaml_2.3.12                 knitr_1.51                 
#> [21] S4Arrays_1.12.0             htmlwidgets_1.6.4          
#> [23] DelayedArray_0.38.2         plyr_1.8.9                 
#> [25] RColorBrewer_1.1-3          abind_1.4-8                
#> [27] purrr_1.2.2                 BiocGenerics_0.58.1        
#> [29] desc_1.4.3                  grid_4.6.1                 
#> [31] stats4_4.6.1                Rhdf5lib_2.0.0             
#> [33] ggplot2_4.0.3               scales_1.4.0               
#> [35] MASS_7.3-65                 MultiAssayExperiment_1.38.0
#> [37] SummarizedExperiment_1.42.0 cli_3.6.6                  
#> [39] rmarkdown_2.32              ragg_1.5.2                 
#> [41] generics_0.1.4              metapod_1.20.0             
#> [43] otel_0.2.0                  scp_1.22.0                 
#> [45] reshape2_1.4.5              BiocBaseUtils_1.14.2       
#> [47] cachem_1.1.0                rhdf5_2.56.0               
#> [49] stringr_1.6.0               AnnotationFilter_1.36.0    
#> [51] BiocManager_1.30.27         XVector_0.52.0             
#> [53] matrixStats_1.5.0           vctrs_0.7.3                
#> [55] Matrix_1.7-5                slam_0.1-56                
#> [57] jsonlite_2.0.0              bookdown_0.48              
#> [59] IRanges_2.46.0              S4Vectors_0.50.2           
#> [61] IHW_1.40.0                  ggrepel_0.9.8              
#> [63] clue_0.3-68                 systemfonts_1.3.2          
#> [65] jquerylib_0.1.4             tidyr_1.3.2                
#> [67] glue_1.8.1                  pkgdown_2.2.1              
#> [69] QFeatures_1.22.0            stringi_1.8.9              
#> [71] gtable_0.3.6                GenomicRanges_1.64.0       
#> [73] lpsymphony_1.40.0           tibble_3.3.1               
#> [75] pillar_1.11.1               htmltools_0.5.9            
#> [77] Seqinfo_1.2.0               rhdf5filters_1.24.1        
#> [79] MuData_1.16.0               R6_2.6.1                   
#> [81] textshaping_1.0.5           evaluate_1.0.5             
#> [83] lattice_0.22-9              Biobase_2.72.0             
#> [85] bslib_0.12.0                Rcpp_1.1.2                 
#> [87] fdrtool_1.2.18              SparseArray_1.12.2         
#> [89] xfun_0.60                   MsCoreUtils_1.24.0         
#> [91] fs_2.1.0                    MatrixGenerics_1.24.0      
#> [93] pkgconfig_2.0.3
```
