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

library(SPATA2) 
library(tidyverse)
library(tidylog)
library(EBImage)
library(hdf5r)
library(tensorflow)
library(keras)


# check -------------------------------------------------------------------

packageVersion("tensorflow")
# [1] ‘2.6.2’

packageVersion("keras")
# [1] ‘2.6.0’


# dir ---------------------------------------------------------------------

dir.1 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/2_00_stRNAseq_data_SB28_setup/"
dir.2 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/2_01_stRNAseq_autodecoder_denoising/"
if(dir.exists(dir.2) == F){dir.create(dir.2)}


# 1. load the spata-object ------------------------------------------------

setwd(dir.1)
in.f <- "visium_SB28_1.rds"
spata.obj <- readRDS(in.f)


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
# Time difference of 17.7155 mins


# check -------------------------------------------------------------------

getAutoencoderAssessment(spata.obj, of_sample = NA)[[1]] %>% 
  as_tibble() %>% 
  arrange(desc(total_var)) %>% 
  print(n = Inf)
# # A tibble: 15 × 3
# activation bottleneck total_var
# <chr>      <fct>          <dbl>
# 1 selu       56             1928.
# 2 selu       32             1835.
# 3 selu       40             1829.
# 4 relu       40             1807.
# 5 selu       48             1707.
# 6 selu       64             1619.
# 7 relu       56             1570.
# 8 sigmoid    40             1534.
# 9 relu       64             1471.
# 10 relu       48             1470.
# 11 relu       32             1466.
# 12 sigmoid    56             1366.
# 13 sigmoid    64             1328.
# 14 sigmoid    32             1274.
# 15 sigmoid    48             1171.


# denoising the data --------------------------------------------------

# https://themilolab.github.io/SPATA2/reference/runAutoencoderDenoising.html

START.TIME <- Sys.time() 

setwd(dir.2)
out.f <- "01_plot_runAutoencoderDenoising.pdf"
pdf(out.f, w = 8, h = 4.5)

spata.obj <- spata.obj %>%   
  runAutoencoderDenoising(
    activation = "selu", 
    bottleneck = 56, 
    epochs = 20, 
    layers = c(128, 64, 32), 
    dropout = 0.1, 
    display_plot = T, 
    genes = c("Aif1", "Ptprc")
  )

dev.off()

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 1.340611 mins


#  check ------------------------------------------------------------------

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
# Bottleneck neurons: 56
# Dropout: 0.1
# Epochs: 20
# Layers: 128, 64 and 32


# save --------------------------------------------------------------------

START.TIME <- Sys.time() 

setwd(dir.2)
out.f <- "visium_SB28_1_denoised.rds"
saveRDS(spata.obj, out.f)

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 1.175425 mins


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
# [1] keras_2.6.0      tensorflow_2.6.0 hdf5r_1.3.9      EBImage_4.34.0   tidylog_1.0.2    forcats_0.5.1    stringr_1.5.1   
# [8] dplyr_1.1.4      purrr_1.0.2      readr_2.1.5      tidyr_1.3.1      tibble_3.2.1     ggplot2_3.5.1    tidyverse_1.3.1 
# [15] SPATA2_2.0.4    
# 
# loaded via a namespace (and not attached):
# [1] readxl_1.3.1                backports_1.5.0             plyr_1.8.7                  igraph_1.3.1               
# [5] lazyeval_0.2.2              sp_2.1-4                    splines_4.1.3               listenv_0.8.0              
# [9] tfruns_1.5.0                scattermore_1.2             GenomeInfoDb_1.28.1         digest_0.6.37              
# [13] htmltools_0.5.8.1           tiff_0.1-11                 magrittr_2.0.3              SPATAData_0.0.0.9000       
# [17] tensor_1.5                  cluster_2.1.2               ROCR_1.0-11                 tzdb_0.4.0                 
# [21] globals_0.15.0              modelr_0.1.8                matrixStats_0.62.0          timechange_0.3.0           
# [25] spatstat.sparse_3.0-2       jpeg_0.1-9                  colorspace_2.1-1            rvest_1.0.1                
# [29] ggrepel_0.9.1               haven_2.4.3                 crayon_1.5.3                RCurl_1.98-1.3             
# [33] jsonlite_1.8.9              progressr_0.10.0            spatstat.data_3.0-1         zeallot_0.1.0              
# [37] survival_3.2-12             zoo_1.8-10                  glue_1.8.0                  polyclip_1.10-7            
# [41] gtable_0.3.6                zlibbioc_1.38.0             XVector_0.32.0              leiden_0.3.9               
# [45] DelayedArray_0.18.0         future.apply_1.8.1          SingleCellExperiment_1.14.1 BiocGenerics_0.38.0        
# [49] abind_1.4-5                 scales_1.3.0                DBI_1.2.3                   spatstat.random_3.1-5      
# [53] miniUI_0.1.1.1              Rcpp_1.0.13-1               viridisLite_0.4.2           xtable_1.8-4               
# [57] units_0.8-5                 reticulate_1.40.0           bit_4.5.0.1                 clisymbols_1.2.0           
# [61] stats4_4.1.3                htmlwidgets_1.6.4           httr_1.4.7                  anndata_0.7.5.6            
# [65] RColorBrewer_1.1-3          Seurat_4.3.0.1              ica_1.0-2                   pkgconfig_2.0.3            
# [69] farver_2.1.2                uwot_0.1.16                 dbplyr_2.1.1                deldir_1.0-6               
# [73] utf8_1.2.4                  locfit_1.5-9.4              labeling_0.4.3              tidyselect_1.2.1           
# [77] rlang_1.1.4                 reshape2_1.4.4              later_1.4.1                 cellranger_1.1.0           
# [81] munsell_0.5.1               tools_4.1.3                 cli_3.6.3                   generics_0.1.3             
# [85] broom_1.0.7                 ggridges_0.5.6              fastmap_1.2.0               fftwtools_0.9-11           
# [89] goftest_1.2-2               bit64_4.5.2                 fs_1.6.5                    fitdistrplus_1.1-5         
# [93] RANN_2.6.1                  pbapply_1.5-0               future_1.25.0               nlme_3.1-152               
# [97] whisker_0.4                 mime_0.12                   xml2_1.3.2                  compiler_4.1.3             
# [101] rstudioapi_0.17.1           plotly_4.10.4               png_0.1-8                   spatstat.utils_3.1-1       
# [105] reprex_2.0.1                confuns_1.0.3               stringi_1.8.4               lattice_0.20-44            
# [109] Matrix_1.6-4                psych_2.4.12                vctrs_0.6.5                 pillar_1.10.0              
# [113] lifecycle_1.0.4             spatstat.geom_3.2-4         lmtest_0.9-40               RcppAnnoy_0.0.22           
# [117] data.table_1.16.4           cowplot_1.1.1               bitops_1.0-9                irlba_2.3.5.1              
# [121] httpuv_1.6.15               patchwork_1.3.0             GenomicRanges_1.44.0        R6_2.5.1                   
# [125] promises_1.3.2              KernSmooth_2.23-20          gridExtra_2.3               IRanges_2.26.0             
# [129] parallelly_1.31.1           codetools_0.2-18            MASS_7.3-54                 assertthat_0.2.1           
# [133] SummarizedExperiment_1.22.0 withr_3.0.2                 SeuratObject_4.1.3          mnormt_2.1.1               
# [137] sctransform_0.3.5           S4Vectors_0.30.2            GenomeInfoDbData_1.2.6      mgcv_1.8-36                
# [141] parallel_4.1.3              hms_1.1.3                   grid_4.1.3                  MatrixGenerics_1.4.3       
# [145] Rtsne_0.16                  spatstat.explore_3.2-1      base64enc_0.1-3             Biobase_2.52.0             
# [149] shiny_1.10.0                lubridate_1.9.4            


# end ---------------------------------------------------------------------
