# Takahide Nejo, MD, PhD
# Okada Lab, Department of Neurological Surgery, UC San Francisco


# note ---------------------------------------------------------------------

# Nejo T, et al., Glioma-neuronal circuit remodeling induces regional immunosuppression
# Analysis of SPATAData human GBM 10X Visium st-RNA-seq data

# Below is the script for one sample ("260_T"). The data from the other samples were processed and analyzed similarly. 

# Using different versions of packages, such as SPATA2, may lead to varying outcomes, which we have not thoroughly validated.

# ref: https://themilolab.github.io/SPATA2/articles/spata-data.html


# rm all ------------------------------------------------------------------

rm(list = ls(all.names = T))


# packages ----------------------------------------------------------------

library(SPATA2) # devtools::install_github("theMILOlab/SPATA2", ref = 'v2.0.4')
library(tidyverse)
library(tidylog)
library(SPATAData) 
# library(EBImage)
# library(hdf5r)


# check -------------------------------------------------------------------

packageVersion("SPATA2")
# [1] ‘2.0.4’

packageVersion("SPATAData")
# [1] ‘0.0.0.9000’


# dir ---------------------------------------------------------------------

dir.1 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/2_10_human_stRNAseq_data_setup/"
if(dir.exists(dir.1) == F){dir.create(dir.1)}


# check the data availability ------------------------------------------------------

SPATAData::validSampleNames() %>% print()


# download the data ---------------------------------------------

START.TIME <- Sys.time() 

spata.obj <- downloadSpataObject(sample_name = "260_T")

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 21.20504 secs


# check -------------------------------------------------------------------

spata.obj 
# An object of class 'spata2' that contains 1 sample named '260_T'.

spata.obj %>% 
  getCoordsDf() %>% 
  colnames() %>% 
  print()
# [1] "barcodes" "sample"   "x"        "y"        "section"  "outline" 

spata.obj %>% 
  getCoordsDf() %>% 
  print()
# # A tibble: 2,997 × 6
# barcodes           sample     x     y section outline
# <chr>              <chr>  <dbl> <dbl> <chr>   <lgl>  
# 1 AAACAATCTACTAGCA-1 260_T  231.   532. 1       TRUE   
# 2 AAACAGAGCGACTCCT-1 260_T  420.   459. 1       FALSE  
# 3 AAACAGCTTTCAGAAG-1 260_T  101.   274. 1       FALSE  
# 4 AAACAGGGTCTATATT-1 260_T  116.   248. 1       FALSE  
# 5 AAACATGGTGAGAGGA-1 260_T   66.7  151. 1       TRUE   
# 6 AAACCGGGTAGGTACC-1 260_T  172.   280. 1       FALSE  
# 7 AAACCGTTCGTCCAGG-1 260_T  224.   214. 1       FALSE  
# 8 AAACCTCATGAAGTTG-1 260_T  139.   312. 1       FALSE  
# 9 AAACGAAGAACATACC-1 260_T  309.   512. 1       FALSE  
# 10 AAACGAAGATGGAGTA-1 260_T   81.8  177. 1       FALSE  
# # ℹ 2,987 more rows
# # ℹ Use `print(n = ...)` to see more rows


# save --------------------------------------------------------------------

START.TIME <- Sys.time() 

setwd(dir.1)
out.f <- "visium_SPATAData_260_T.rds"
saveRDS(spata.obj, out.f)

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 28.38495 secs


# si ----------------------------------------------------------------------

Sys.time()

sessionInfo()
# R version 4.1.3 (2022-03-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Rocky Linux 8.10 (Green Obsidian)
# 
# Matrix products: default
# BLAS:   /software/c4/cbi/software/_rocky8/R-4.1.3/lib64/R/lib/libRblas.so
# LAPACK: /software/c4/cbi/software/_rocky8/R-4.1.3/lib64/R/lib/libRlapack.so
# 
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] SPATAData_0.0.0.9000 tidylog_1.0.2        forcats_0.5.1       
# [4] stringr_1.5.1        dplyr_1.1.4          purrr_1.0.2         
# [7] readr_2.1.5          tidyr_1.3.1          tibble_3.2.1        
# [10] ggplot2_3.5.1        tidyverse_1.3.1      SPATA2_2.0.4        
# 
# loaded via a namespace (and not attached):
# [1] readxl_1.3.1                backports_1.5.0            
# [3] plyr_1.8.7                  igraph_1.3.1               
# [5] lazyeval_0.2.2              sp_2.1-4                   
# [7] splines_4.1.3               listenv_0.8.0              
# [9] scattermore_1.2             GenomeInfoDb_1.28.1        
# [11] digest_0.6.37               yulab.utils_0.1.8          
# [13] htmltools_0.5.8.1           tiff_0.1-11                
# [15] magrittr_2.0.3              tensor_1.5                 
# [17] cluster_2.1.2               ROCR_1.0-11                
# [19] tzdb_0.4.0                  globals_0.15.0             
# [21] modelr_0.1.8                matrixStats_0.62.0         
# [23] timechange_0.3.0            spatstat.sparse_3.0-2      
# [25] jpeg_0.1-9                  colorspace_2.1-1           
# [27] rvest_1.0.1                 ggrepel_0.9.1              
# [29] haven_2.4.3                 crayon_1.5.3               
# [31] RCurl_1.98-1.3              jsonlite_1.8.9             
# [33] progressr_0.10.0            spatstat.data_3.0-1        
# [35] survival_3.2-12             zoo_1.8-10                 
# [37] glue_1.8.0                  polyclip_1.10-7            
# [39] gtable_0.3.6                zlibbioc_1.38.0            
# [41] XVector_0.32.0              leiden_0.3.9               
# [43] V8_6.0.0                    DelayedArray_0.18.0        
# [45] future.apply_1.8.1          SingleCellExperiment_1.14.1
# [47] BiocGenerics_0.38.0         abind_1.4-5                
# [49] scales_1.3.0                DBI_1.2.3                  
# [51] spatstat.random_3.1-5       miniUI_0.1.1.1             
# [53] Rcpp_1.0.13-1               viridisLite_0.4.2          
# [55] xtable_1.8-4                units_0.8-5                
# [57] gridGraphics_0.5-1          reticulate_1.40.0          
# [59] clisymbols_1.2.0            stats4_4.1.3               
# [61] htmlwidgets_1.6.4           httr_1.4.7                 
# [63] anndata_0.7.5.6             RColorBrewer_1.1-3         
# [65] Seurat_4.3.0.1              ica_1.0-2                  
# [67] pkgconfig_2.0.3             farver_2.1.2               
# [69] uwot_0.1.16                 dbplyr_2.1.1               
# [71] deldir_1.0-6                utf8_1.2.4                 
# [73] locfit_1.5-9.4              ggplotify_0.1.2            
# [75] tidyselect_1.2.1            rlang_1.1.4                
# [77] reshape2_1.4.4              later_1.4.1                
# [79] cellranger_1.1.0            munsell_0.5.1              
# [81] tools_4.1.3                 cli_3.6.3                  
# [83] dbscan_1.2-0                generics_0.1.3             
# [85] broom_1.0.7                 ggridges_0.5.6             
# [87] fastmap_1.2.0               fftwtools_0.9-11           
# [89] goftest_1.2-2               fs_1.6.5                   
# [91] fitdistrplus_1.1-5          RANN_2.6.1                 
# [93] pbapply_1.5-0               future_1.25.0              
# [95] nlme_3.1-152                mime_0.12                  
# [97] xml2_1.3.2                  concaveman_1.1.0           
# [99] compiler_4.1.3              rstudioapi_0.17.1          
# [101] curl_6.0.1                  plotly_4.10.4              
# [103] png_0.1-8                   spatstat.utils_3.1-1       
# [105] reprex_2.0.1                confuns_1.0.3              
# [107] stringi_1.8.4               lattice_0.20-44            
# [109] Matrix_1.6-4                vctrs_0.6.5                
# [111] pillar_1.10.0               lifecycle_1.0.4            
# [113] spatstat.geom_3.2-4         lmtest_0.9-40              
# [115] RcppAnnoy_0.0.22            data.table_1.16.4          
# [117] cowplot_1.1.1               bitops_1.0-9               
# [119] irlba_2.3.5.1               httpuv_1.6.15              
# [121] patchwork_1.3.0             GenomicRanges_1.44.0       
# [123] R6_2.5.1                    promises_1.3.2             
# [125] KernSmooth_2.23-20          gridExtra_2.3              
# [127] IRanges_2.26.0              parallelly_1.31.1          
# [129] codetools_0.2-18            MASS_7.3-54                
# [131] assertthat_0.2.1            SummarizedExperiment_1.22.0
# [133] withr_3.0.2                 SeuratObject_4.1.3         
# [135] sctransform_0.3.5           S4Vectors_0.30.2           
# [137] GenomeInfoDbData_1.2.6      parallel_4.1.3             
# [139] hms_1.1.3                   EBImage_4.34.0             
# [141] grid_4.1.3                  MatrixGenerics_1.4.3       
# [143] Rtsne_0.16                  spatstat.explore_3.2-1     
# [145] Biobase_2.52.0              shiny_1.10.0               
# [147] lubridate_1.9.4            


# end ---------------------------------------------------------------------


