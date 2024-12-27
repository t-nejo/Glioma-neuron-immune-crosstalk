# Takahide Nejo, MD, PhD
# Okada Lab, Department of Neurological Surgery, UC San Francisco


# note ---------------------------------------------------------------------

# Nejo T, et al., Glioma-neuronal circuit remodeling induces regional immunosuppression
# Analysis of SB28 tumor-bearing mouse brain 10X Visium st-RNA-seq data

# Below is the script for the first sample. The data from the second sample was processed and analyzed similarly. 

# Using different versions of packages, such as SPATA2, may lead to varying outcomes, which we have not thoroughly validated.


# rm all ------------------------------------------------------------------

rm(list = ls(all.names = T))


# packages ----------------------------------------------------------------

library(SPATA2) # devtools::install_github("theMILOlab/SPATA2", ref = 'v2.0.4')
# devtools::install_github("kueckelj/confuns")
library(tidyverse)
library(tidylog)
library(EBImage)
library(hdf5r)


# dir ---------------------------------------------------------------------

dir.1 <- "/okadalab/data1/tnejo/neuro_immune_proj/10x_visium_01/processed_data/"
dir.2 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/2_00_stRNAseq_data_SB28_setup/"
if(dir.exists(dir.2) == F){dir.create(dir.2)}


# 1. Create a spata-object ------------------------------------------------

# ref: https://themilolab.github.io/SPATA2/articles/spata-v2-object-initiation-and-manipulation.html#from-10x-visium-folders

DIR <- paste0(dir.1, "Visium_FFPE_Mouse_Brain_Okada-TN-3776_A1")
SAMPLE <- "A1_166L_Cas9"

START.TIME <- Sys.time() 

spata.obj <-
  initiateSpataObject_10X(
    directory_10X = DIR, # the directory that contains the data
    sample_name = SAMPLE #, 
    # to.upper = TRUE
  )

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 38.72924 secs


# check -------------------------------------------------------------------

spata.obj 
# An object of class 'spata2' that contains 1 sample named 'A1_166L_Cas9'.

spata.obj %>% 
  getCoordsDf() %>% 
  print()
# # A tibble: 1,975 × 11
# barcodes               x     y tissue   row   col imagerow imagecol sample       section outline
# <chr>              <dbl> <dbl>  <int> <int> <int>    <dbl>    <dbl> <chr>        <chr>   <lgl>  
# 1 AAACAATCTACTAGCA-1  64.5  221.      1   589    43    -8400     1562 A1_166L_Cas9 1       FALSE  
# 2 AAACCGGGTAGGTACC-1 320.   163.      1   550    28    -9792     7742 A1_166L_Cas9 1       FALSE  
# 3 AAACCGTTCGTCCAGG-1 386.   215.      1   540    42    -8522     9333 A1_166L_Cas9 1       FALSE  
# 4 AAACCTCATGAAGTTG-1 287.   129.      1   555    19   -10609     6945 A1_166L_Cas9 1       FALSE  
# 5 AAACGAAGAACATACC-1  84.5  300.      1   586    64    -6488     2045 A1_166L_Cas9 1       TRUE   
# 6 AAACGAGACGGTTGAT-1 275.   355.      1   557    79    -5139     6650 A1_166L_Cas9 1       FALSE  
# 7 AAACTGCTGGCTCCAA-1 340.   310.      1   547    67    -6239     8232 A1_166L_Cas9 1       FALSE  
# 8 AAACTTGCAAACGTAT-1 339.   129.      1   547    19   -10614     8214 A1_166L_Cas9 1       FALSE  
# 9 AAAGACCCAAGTCGCG-1 110.   239.      1   582    48    -7949     2674 A1_166L_Cas9 1       FALSE  
# 10 AAAGACTGGGCGCTTT-1 234.   114.      1   563    15   -10969     5675 A1_166L_Cas9 1       FALSE  
# # ℹ 1,965 more rows
# # ℹ Use `print(n = ...)` to see more rows


# save --------------------------------------------------------------------

setwd(dir.2)
out.f <- paste0("visium_", "SB28_1", ".rds")
saveRDS(spata.obj, out.f)


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
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] hdf5r_1.3.9     EBImage_4.34.0  tidylog_1.0.2   forcats_0.5.1   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2    
# [8] readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_1.3.1 SPATA2_2.0.4   
# 
# loaded via a namespace (and not attached):
# [1] utf8_1.2.4                  spatstat.explore_3.2-1      reticulate_1.40.0           tidyselect_1.2.1           
# [5] htmlwidgets_1.6.4           grid_4.1.3                  Rtsne_0.16                  devtools_2.4.2             
# [9] munsell_0.5.1               codetools_0.2-18            ica_1.0-2                   units_0.8-5                
# [13] future_1.25.0               miniUI_0.1.1.1              withr_3.0.2                 spatstat.random_3.1-5      
# [17] colorspace_2.1-1            progressr_0.10.0            Biobase_2.52.0              rstudioapi_0.17.1          
# [21] Seurat_4.3.0.1              stats4_4.1.3                SingleCellExperiment_1.14.1 ROCR_1.0-11                
# [25] tensor_1.5                  listenv_0.8.0               MatrixGenerics_1.4.3        GenomeInfoDbData_1.2.6     
# [29] polyclip_1.10-7             bit64_4.5.2                 farver_2.1.2                rprojroot_2.0.4            
# [33] parallelly_1.31.1           vctrs_0.6.5                 generics_0.1.3              timechange_0.3.0           
# [37] R6_2.5.1                    GenomeInfoDb_1.28.1         locfit_1.5-9.4              concaveman_1.1.0           
# [41] bitops_1.0-9                spatstat.utils_3.1-1        cachem_1.1.0                DelayedArray_0.18.0        
# [45] assertthat_0.2.1            promises_1.3.2              scales_1.3.0                gtable_0.3.6               
# [49] globals_0.15.0              processx_3.5.2              goftest_1.2-2               clisymbols_1.2.0           
# [53] rlang_1.1.4                 splines_4.1.3               lazyeval_0.2.2              spatstat.geom_3.2-4        
# [57] broom_1.0.7                 reshape2_1.4.4              abind_1.4-5                 modelr_0.1.8               
# [61] backports_1.5.0             httpuv_1.6.15               tools_4.1.3                 usethis_2.0.1              
# [65] ellipsis_0.3.2              RColorBrewer_1.1-3          BiocGenerics_0.38.0         sessioninfo_1.1.1          
# [69] ggridges_0.5.6              Rcpp_1.0.13-1               plyr_1.8.7                  zlibbioc_1.38.0            
# [73] RCurl_1.98-1.3              ps_1.6.0                    prettyunits_1.2.0           dbscan_1.2-0               
# [77] deldir_1.0-6                pbapply_1.5-0               cowplot_1.1.1               S4Vectors_0.30.2           
# [81] zoo_1.8-10                  SeuratObject_4.1.3          SummarizedExperiment_1.22.0 haven_2.4.3                
# [85] ggrepel_0.9.1               cluster_2.1.2               fs_1.6.5                    magrittr_2.0.3             
# [89] data.table_1.16.4           scattermore_1.2             lmtest_0.9-40               reprex_2.0.1               
# [93] RANN_2.6.1                  fitdistrplus_1.1-5          anndata_0.7.5.6             matrixStats_0.62.0         
# [97] pkgload_1.2.1               hms_1.1.3                   patchwork_1.3.0             mime_0.12                  
# [101] fftwtools_0.9-11            xtable_1.8-4                jpeg_0.1-9                  readxl_1.3.1               
# [105] IRanges_2.26.0              gridExtra_2.3               testthat_3.0.4              compiler_4.1.3             
# [109] V8_6.0.0                    KernSmooth_2.23-20          crayon_1.5.3                SPATAData_0.0.0.9000       
# [113] htmltools_0.5.8.1           later_1.4.1                 tzdb_0.4.0                  tiff_0.1-11                
# [117] lubridate_1.9.4             DBI_1.2.3                   dbplyr_2.1.1                MASS_7.3-54                
# [121] Matrix_1.6-4                cli_3.6.3                   parallel_4.1.3              igraph_1.3.1               
# [125] GenomicRanges_1.44.0        pkgconfig_2.0.3             sp_2.1-4                    confuns_1.0.3              
# [129] plotly_4.10.4               spatstat.sparse_3.0-2       xml2_1.3.2                  XVector_0.32.0             
# [133] rvest_1.0.1                 callr_3.7.0                 digest_0.6.37               sctransform_0.3.5          
# [137] RcppAnnoy_0.0.22            spatstat.data_3.0-1         cellranger_1.1.0            leiden_0.3.9               
# [141] uwot_0.1.16                 curl_6.0.1                  shiny_1.10.0                lifecycle_1.0.4            
# [145] nlme_3.1-152                jsonlite_1.8.9              desc_1.3.0                  viridisLite_0.4.2          
# [149] pillar_1.10.0               lattice_0.20-44             fastmap_1.2.0               httr_1.4.7                 
# [153] pkgbuild_1.2.0              survival_3.2-12             glue_1.8.0                  remotes_2.4.2.1            
# [157] png_0.1-8                   bit_4.5.0.1                 stringi_1.8.4               memoise_2.0.1              
# [161] irlba_2.3.5.1               future.apply_1.8.1         


# end ---------------------------------------------------------------------


