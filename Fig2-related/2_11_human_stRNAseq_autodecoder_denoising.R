# Takahide Nejo, MD, PhD
# Okada Lab, Department of Neurological Surgery, UC San Francisco


# note ---------------------------------------------------------------------

# Nejo T, et al., Glioma-neuronal circuit remodeling induces regional immunosuppression
# Analysis of SPATAData human GBM 10X Visium st-RNA-seq data

# Below is the script for one sample ("260_T"). The data from the other samples were processed and analyzed similarly. 

# Using different versions of packages, such as SPATA2, may lead to varying outcomes, which we have not thoroughly validated.

# The RDS object "visium_SPATAData_260_T.rds" is used as input. For more details, please refer to the previous step, "2_10_human_stRNAseq_data_setup.R".

# ref: https://themilolab.github.io/SPATA2/articles/spata-data.html


# rm all ------------------------------------------------------------------

rm(list = ls(all.names = T))


# packages ----------------------------------------------------------------

library(SPATA2)
library(tidyverse)
library(tidylog)
library(tensorflow)
library(keras)


# check -------------------------------------------------------------------

packageVersion("SPATA2")
# [1] ‘2.0.4’

packageVersion("tensorflow")
# [1] ‘2.6.0’

packageVersion("keras")
# [1] ‘2.6.0’


# dir ---------------------------------------------------------------------

dir.1 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/2_10_human_stRNAseq_data_setup/"
dir.2 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/2_11_human_stRNAseq_autodecoder_denoising/"
if(dir.exists(dir.2) == F){dir.create(dir.2)}


# 1. load the spata-object ------------------------------------------------

START.TIME <- Sys.time()

setwd(dir.1)
in.f <- "visium_SPATAData_260_T.rds"
spata.obj <- readRDS(in.f)

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 5.253172 secs


# check -------------------------------------------------------------------

spata.obj 
# An object of class 'spata2' that contains 1 sample named '260_T'.

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


# check -------------------------------------------------------------------

# all expression matrices before denoising

spata.obj %>% 
  getExpressionMatrixNames() %>% 
  print()
# [1] "scaled"

# active expression matrix before denoising
spata.obj %>% 
  getActiveMatrixName() %>% 
  print()
# [1] "scaled"


# find the optimal neural network set up -----------------------------------

# ref: https://themilolab.github.io/SPATA2/reference/runAutoencoderAssessment.html

START.TIME <- Sys.time() 

spata.obj <- spata.obj %>% 
  runAutoencoderAssessment(
    activations = c("relu", "selu", "sigmoid"), 
    bottlenecks = c(32, 40, 48, 56, 64),
    epochs = 20, 
    layers = c(128, 64, 32), 
    dropout = 0.1
  )

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 23.09346 mins


# check -------------------------------------------------------------------

getAutoencoderAssessment(spata.obj, of_sample = NA)[[1]] %>% 
  as_tibble() %>% 
  arrange(desc(total_var)) %>% 
  print(n = Inf)
# # A tibble: 15 × 3
# activation bottleneck total_var
# <chr>      <fct>          <dbl>
# 1 selu       32              477.
# 2 selu       64              476.
# 3 selu       40              435.
# 4 selu       56              426.
# 5 selu       48              393.
# 6 relu       40              357.
# 7 relu       64              344.
# 8 relu       32              338.
# 9 relu       48              309.
# 10 relu       56              283.
# 11 sigmoid    64              258.
# 12 sigmoid    32              249.
# 13 sigmoid    48              235.
# 14 sigmoid    56              232.
# 15 sigmoid    40              228.


# denoising the data --------------------------------------------------

# https://themilolab.github.io/SPATA2/reference/runAutoencoderDenoising.html

START.TIME <- Sys.time() 

setwd(dir.2)
out.f <- "test_01_plot_runAutoencoderDenoising.pdf"
pdf(out.f, w = 8, h = 4.5)

spata.obj <- spata.obj %>%   
  runAutoencoderDenoising(
    activation = "selu", 
    bottleneck = 32, 
    epochs = 20, 
    layers = c(128, 64, 32), 
    dropout = 0.1, 
    display_plot = T, 
    genes = c("AIF1", "PTPRC")
  )

dev.off()

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 1.656749 mins


# check ------------------------------------------------------------------

# all expression matrices after denoising

spata.obj %>% 
  getExpressionMatrixNames() %>% 
  print()
# [1] "scaled"   "denoised"


# active expression matrix after denoising

spata.obj %>% 
  getActiveMatrixName() %>% 
  print()
# [1] "denoised"

# summary

spata.obj %>% 
  printAutoencoderSummary() %>% 
  print()

# The neural network that generated matrix 'denoised' was constructed with the following adjustments: 
#   
# Activation function: selu
# Bottleneck neurons: 32
# Dropout: 0.1
# Epochs: 20
# Layers: 128, 64 and 32


# save --------------------------------------------------------------------

START.TIME <- Sys.time() 

setwd(dir.2)
out.f <- "visium_SPATAData_260_T_denoised.rds"
saveRDS(spata.obj, out.f)

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 1.744315 mins


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
# [1] keras_2.6.0      tensorflow_2.6.0 tidylog_1.0.2    forcats_0.5.1    stringr_1.5.1    dplyr_1.1.4     
# [7] purrr_1.0.2      readr_2.1.5      tidyr_1.3.1      tibble_3.2.1     ggplot2_3.5.1    tidyverse_1.3.1 
# [13] SPATA2_2.0.4    
# 
# loaded via a namespace (and not attached):
# [1] utf8_1.2.4                  spatstat.explore_3.2-1      reticulate_1.40.0           tidyselect_1.2.1           
# [5] htmlwidgets_1.6.4           grid_4.1.3                  Rtsne_0.16                  munsell_0.5.1              
# [9] codetools_0.2-18            ica_1.0-2                   units_0.8-5                 future_1.25.0              
# [13] miniUI_0.1.1.1              withr_3.0.2                 spatstat.random_3.1-5       colorspace_2.1-1           
# [17] progressr_0.10.0            Biobase_2.52.0              rstudioapi_0.17.1           Seurat_4.3.0.1             
# [21] stats4_4.1.3                SingleCellExperiment_1.14.1 ROCR_1.0-11                 ggsignif_0.6.2             
# [25] tensor_1.5                  listenv_0.8.0               labeling_0.4.3              MatrixGenerics_1.4.3       
# [29] GenomeInfoDbData_1.2.6      mnormt_2.1.1                polyclip_1.10-7             farver_2.1.2               
# [33] parallelly_1.31.1           vctrs_0.6.5                 generics_0.1.3              timechange_0.3.0           
# [37] R6_2.5.1                    GenomeInfoDb_1.28.1         locfit_1.5-9.4              bitops_1.0-9               
# [41] spatstat.utils_3.1-1        DelayedArray_0.18.0         assertthat_0.2.1            promises_1.3.2             
# [45] scales_1.3.0                gtable_0.3.6                globals_0.15.0              goftest_1.2-2              
# [49] rlang_1.1.4                 clisymbols_1.2.0            zeallot_0.1.0               splines_4.1.3              
# [53] rstatix_0.7.0               lazyeval_0.2.2              spatstat.geom_3.2-4         broom_1.0.7                
# [57] reshape2_1.4.4              abind_1.4-5                 modelr_0.1.8                backports_1.5.0            
# [61] httpuv_1.6.15               tools_4.1.3                 psych_2.4.12                RColorBrewer_1.1-3         
# [65] BiocGenerics_0.38.0         ggridges_0.5.6              Rcpp_1.0.13-1               plyr_1.8.7                 
# [69] base64enc_0.1-3             zlibbioc_1.38.0             RCurl_1.98-1.3              ggpubr_0.4.0               
# [73] deldir_1.0-6                pbapply_1.5-0               cowplot_1.1.1               S4Vectors_0.30.2           
# [77] zoo_1.8-10                  SeuratObject_4.1.3          SummarizedExperiment_1.22.0 haven_2.4.3                
# [81] ggrepel_0.9.1               cluster_2.1.2               fs_1.6.5                    magrittr_2.0.3             
# [85] data.table_1.16.4           scattermore_1.2             openxlsx_4.2.4              lmtest_0.9-40              
# [89] reprex_2.0.1                RANN_2.6.1                  whisker_0.4                 fitdistrplus_1.1-5         
# [93] anndata_0.7.5.6             matrixStats_0.62.0          hms_1.1.3                   patchwork_1.3.0            
# [97] mime_0.12                   fftwtools_0.9-11            xtable_1.8-4                rio_0.5.27                 
# [101] jpeg_0.1-9                  readxl_1.3.1                IRanges_2.26.0              gridExtra_2.3              
# [105] tfruns_1.5.0                compiler_4.1.3              KernSmooth_2.23-20          crayon_1.5.3               
# [109] SPATAData_0.0.0.9000        htmltools_0.5.8.1           mgcv_1.8-36                 later_1.4.1                
# [113] tzdb_0.4.0                  tiff_0.1-11                 lubridate_1.9.4             DBI_1.2.3                  
# [117] dbplyr_2.1.1                MASS_7.3-54                 Matrix_1.6-4                car_3.0-11                 
# [121] cli_3.6.3                   parallel_4.1.3              igraph_1.3.1                GenomicRanges_1.44.0       
# [125] pkgconfig_2.0.3             foreign_0.8-81              sp_2.1-4                    confuns_1.0.3              
# [129] plotly_4.10.4               spatstat.sparse_3.0-2       xml2_1.3.2                  XVector_0.32.0             
# [133] rvest_1.0.1                 digest_0.6.37               sctransform_0.3.5           RcppAnnoy_0.0.22           
# [137] spatstat.data_3.0-1         cellranger_1.1.0            leiden_0.3.9                uwot_0.1.16                
# [141] curl_6.0.1                  shiny_1.10.0                EBImage_4.34.0              lifecycle_1.0.4            
# [145] nlme_3.1-152                jsonlite_1.8.9              carData_3.0-4               viridisLite_0.4.2          
# [149] pillar_1.10.0               lattice_0.20-44             fastmap_1.2.0               httr_1.4.7                 
# [153] survival_3.2-12             glue_1.8.0                  zip_2.2.0                   png_0.1-8                  
# [157] stringi_1.8.4               irlba_2.3.5.1               future.apply_1.8.1         


# end ---------------------------------------------------------------------


