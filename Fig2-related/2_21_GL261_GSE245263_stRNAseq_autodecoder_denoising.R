# Takahide Nejo, MD, PhD
# Okada Lab, Department of Neurological Surgery, UC San Francisco


# note ---------------------------------------------------------------------

# Nejo T, et al., Glioma-neuronal circuit remodeling induces regional immunosuppression

# Analysis of GL261 tumor-bearing mouse brain 10X Visium st-RNA-seq data (GSE245263)
# original article: García-Vicente L, 2023 bioRxiv. 
# https://www.biorxiv.org/content/10.1101/2023.10.26.564166v1.full

# original data analysis info: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7839621

# Using different versions of packages, such as SPATA2, may lead to varying outcomes, which we have not thoroughly validated.

# The RDS object "visium_GSE245263_GL261_coord_adjusted.rds" is used as input. For more details, please refer to the previous step, "2_20_GL261_GSE245263_stRNAseq_data_setup.R".


# rm all ------------------------------------------------------------------

rm(list = ls(all.names = T))


# packages ----------------------------------------------------------------

library(SPATA2) 
library(tidyverse)
library(tidylog)
library(tensorflow)
library(keras)
library(ggpubr)


# check -------------------------------------------------------------------

packageVersion("SPATA2")
# [1] ‘2.0.4’

packageVersion("tensorflow")
# [1] ‘2.6.2’

packageVersion("keras")
# [1] ‘2.6.0’


# dir ---------------------------------------------------------------------

dir.1 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/2_20_GL261_GSE245263_stRNAseq_data_setup/"
dir.2 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/2_21_GL261_GSE245263_stRNAseq_autodecoder_denoising/"
if(dir.exists(dir.2) == F){dir.create(dir.2)}


# 1. load the spata-object ------------------------------------------------

START.TIME <- Sys.time() 

setwd(dir.1)
in.f <- "visium_GSE245263_GL261_coord_adjusted.rds"
spata.obj <- readRDS(in.f)

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 7.351321 secs


# check -------------------------------------------------------------------

spata.obj 
# An object of class 'spata2' that contains 1 sample named 'GSM7839621'.

spata.obj %>% 
  getCoordsDf() %>% 
  print()

# # A tibble: 2,158 × 6
# barcodes           sample     imagerow imagecol     x     y
# <chr>              <chr>         <dbl>    <dbl> <dbl> <dbl>
# 1 AAACAAGTATCTCCCA-1 GSM7839621     9669     7965  653. 1465.
# 2 AAACAATCTACTAGCA-1 GSM7839621    21459    16591 1446.  885.
# 3 AAACCCGAACGAAATC-1 GSM7839621    10942     6090  738. 1592.
# 4 AAACCGTTCGTCCAGG-1 GSM7839621     9107    16653  615.  881.
# 5 AAACGAAGAACATACC-1 GSM7839621    20723    13544 1397. 1090.
# 6 AAACGAGACGGTTGAT-1 GSM7839621    13428    11322  906. 1239.
# 7 AAACGGGCGTACGGGT-1 GSM7839621     5877     9533  397. 1360.
# 8 AAACTGCTGGCTCCAA-1 GSM7839621    10895    13044  735. 1123.
# 9 AAAGACCCAAGTCGCG-1 GSM7839621    19699    15855 1328.  934.
# 10 AAAGGCTCTCGCGCCG-1 GSM7839621     8363    14765  565. 1008.
# # ℹ 2,148 more rows
# # ℹ Use `print(n = ...)` to see more rows


# check -------------------------------------------------------------------

# all expression matrices before denoising

spata.obj %>% 
  getExpressionMatrixNames() %>% 
  print()
# [1] "normalized" "scaled"    

# active expression matrix before denoising
spata.obj %>% 
  getActiveMatrixName() %>% 
  print()
# [1] "scaled"


g <- spata.obj %>% 
  plotSurface(
    display_image = F, 
    color_by = "Aif1", 
    pt_size = 2, 
    pt_alpha = 0.9, 
  )
plot(g)


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
# Time difference of 31.21572 mins


# check -------------------------------------------------------------------

getAutoencoderAssessment(spata.obj)$df %>% 
  as_tibble() %>% 
  arrange(desc(total_var)) %>% 
  print(n = Inf)
# # A tibble: 15 × 3
# activation bottleneck total_var
# <chr>      <fct>          <dbl>
# 1 selu       56             1213.
# 2 selu       48             1208.
# 3 relu       64             1045.
# 4 selu       32             1009.
# 5 relu       32              988.
# 6 selu       64              970.
# 7 sigmoid    40              957.
# 8 relu       56              954.
# 9 selu       40              934.
# 10 relu       40              929.
# 11 relu       48              918.
# 12 sigmoid    56              890.
# 13 sigmoid    64              874.
# 14 sigmoid    32              857.
# 15 sigmoid    48              792.


# denoising the data --------------------------------------------------

# # https://themilolab.github.io/SPATA2/reference/runAutoencoderDenoising.html
# 
# START.TIME <- Sys.time() 
# 
# spata.obj <- spata.obj %>%   
#   runAutoencoderDenoising(
#     activation = "selu", 
#     bottleneck = 64, 
#     epochs = 20, 
#     layers = c(128, 64, 32), 
#     dropout = 0.1
#     )
# 
# FINISH.TIME <- Sys.time() 
# 
# print(FINISH.TIME - START.TIME)
# # Time difference of 1.969001 mins

# did not work, as described here: 
# ref: https://github.com/theMILOlab/SPATA2/issues/109#issuecomment-2040712194


# alternative method ---------------------------------------------------------

START.TIME <- Sys.time() 

ACTIVATION <- "selu"
BOTTLENECK <- 56
EPOCHS <- 20
LAYERS <- c(128, 64, 32)
DROPOUT <- 0.1

x_train <- spata.obj %>% 
  getExpressionMatrix(mtr_name = "scaled")

input_layer <- keras::layer_input(shape = c(ncol(x_train)))

encoder <- input_layer %>% 
  keras::layer_dense(units = LAYERS[1], activation = ACTIVATION) %>% 
  keras::layer_batch_normalization() %>% 
  keras::layer_dropout(rate = DROPOUT) %>% 
  keras::layer_dense(units = LAYERS[2], activation = ACTIVATION) %>% 
  keras::layer_dropout(rate = DROPOUT) %>% 
  keras::layer_dense(units = LAYERS[3], activation = ACTIVATION) %>% 
  keras::layer_dense(units = BOTTLENECK)

decoder <- encoder %>% 
  keras::layer_dense(units = LAYERS[3], activation = ACTIVATION) %>% 
  keras::layer_dropout(rate = DROPOUT) %>% 
  keras::layer_dense(units = LAYERS[2], activation = ACTIVATION) %>% 
  keras::layer_dropout(rate = DROPOUT) %>% 
  keras::layer_dense(units = LAYERS[1], activation = ACTIVATION) %>% 
  keras::layer_dense(units = c(ncol(x_train)))

autoencoder_model <- keras::keras_model(inputs = input_layer, outputs = decoder)

autoencoder_model %>% 
  keras::compile(loss = "mean_squared_error", optimizer = "adam", metrics = c("accuracy"))

history <- autoencoder_model %>% 
  keras::fit(x_train, 
             x_train, 
             epochs = EPOCHS, 
             shuffle = TRUE, 
             validation_data = list(x_train, x_train), 
             verbose = T)

reconstructed_points <- autoencoder_model %>% 
  keras::predict_on_batch(x = x_train)

base::rownames(reconstructed_points) <- base::rownames(x_train)
base::colnames(reconstructed_points) <- base::colnames(x_train)

SET_UP <- list(activation = ACTIVATION, 
               bottleneck = BOTTLENECK, 
               dropout = DROPOUT, 
               epochs = EPOCHS, 
               input_mtr = "scaled", 
               output_mtr = "denoised", 
               layers = LAYERS)

# spata.obj <- spata.obj %>%
#   SPATA2::addAutoencoderSetUp(
#     mtr_name = "denoised",
#     set_up_list = set_up,
#     of_sample = NULL)

# did not work.

# alternative approach:

spata.obj@autoencoder[["GSM7839621"]][["nn_set_ups"]][["denoised"]] <- SET_UP


spata.obj <- spata.obj %>%
  addExpressionMatrix(
    mtr_name = "denoised", 
    expr_mtr = reconstructed_points, 
    of_sample = NULL)


# spata.obj <- spata.obj %>%
#   computeGeneMetaData(
#     mtr_name = "denoised"
#     )

# did not work.

# alternative approach: 

{
  # computeGeneMetaData
  # # function (object, mtr_name = NULL, verbose = TRUE, ...) 
  # # {
  # #   check_object(object)
  # #   deprecated(...)
  # #   expr_mtr <- getExpressionMatrix(object = object, verbose = verbose)
  # #   if (base::is.null(mtr_name)) {
  # #     mtr_name <- getActiveMatrixName(object)
  # #   }
  # #   meta_data <- computeGeneMetaData2(expr_mtr = expr_mtr, verbose = verbose, 
  # #                                     ...)
  # #   object <- addGeneMetaData(object = object, meta_data_list = c(meta_data, 
  # #                                                                 mtr_name = mtr_name))
  # #   return(object)
  # # }
  # # <bytecode: 0x5561ba922600>
  # #   <environment: namespace:SPATA2>
  #   
  # computeGeneMetaData2
  # # function (expr_mtr, verbose = TRUE, ...) 
  # # {
  # #   confuns::give_feedback(msg = glue::glue("Calculating summary statistics for {base::nrow(expr_mtr)} genes."), 
  # #                          verbose = verbose)
  # #   res_df <- psych::describe(x = base::t(expr_mtr)) %>% base::as.data.frame() %>% 
  # #     dplyr::select(-vars) %>% tibble::rownames_to_column(var = "genes")
  # #   res_list <- list(df = res_df, describe_args = list(...))
  # #   return(res_list)
  # # }
  # # <bytecode: 0x5561ba6311b8>
  # #   <environment: namespace:SPATA2>
}


# alternative to "computeGeneMetaData2"

EXPR_MTX <- reconstructed_points

RES_DF <- psych::describe(x = base::t(EXPR_MTX)) %>% 
  base::as.data.frame() %>%
  dplyr::select(-vars) %>% 
  tibble::rownames_to_column(var = "genes")

# res_list <- list(df = RES_DF, describe_args = list(...))
# return(res_list)

# alternative to "computeGeneMetaData"

RES_DF$mtr_name <- "denoised"

# spata.obj <- spata.obj %>%
#   addGeneMetaData(meta_data_list = RES_DF)

spata.obj@gdata[["GSM7839621"]][["denoised"]] <- RES_DF

spata.obj <- spata.obj %>% 
  setActiveMatrix("denoised")

FINISH.TIME <- Sys.time()

print(FINISH.TIME - START.TIME)
# Time difference of 2.324915 mins


# check ------------------------------------------------------------------

# all expression matrices after denoising

spata.obj %>% 
  getExpressionMatrixNames() %>% 
  print()
# [1] "normalized" "scaled"     "denoised"  


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
out.f <- "visium_GSE245263_GL261_denoised.rds"
saveRDS(spata.obj, out.f)

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 1.467956 mins


# visualize ---------------------------------------------------------------

spata.obj <- spata.obj %>% 
  setActiveMatrix("scaled")

g <- spata.obj %>% 
  plotSurface(
    display_image = T, 
    color_by = "Aif1", 
    pt_alpha = 0.75,  
    pt_size = 2, 
    smooth_span = 0.25
  ) 
plot(g)
g1 <- g

spata.obj <- spata.obj %>% 
  setActiveMatrix("denoised")

g <- spata.obj %>% 
  plotSurface(
    display_image = T, 
    color_by = "Aif1", 
    pt_alpha = 0.75,  
    pt_size = 2, 
    smooth_span = 0.25
  ) 
plot(g)
g2 <- g

g <- ggarrange(
  g1, g2, 
  ncol = 2, nrow = 1
)

setwd(dir.2)
out.f <- "02_surfaceplot_pre_n_post_denoising.pdf"
ggsave(out.f, g, w = 12, h = 4.5)


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
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] ggpubr_0.4.0     keras_2.6.0      tensorflow_2.6.0 tidylog_1.0.2    forcats_0.5.1    stringr_1.5.1   
# [7] dplyr_1.1.4      purrr_1.0.2      readr_2.1.5      tidyr_1.3.1      tibble_3.2.1     ggplot2_3.5.1   
# [13] tidyverse_1.3.1  SPATA2_2.0.4    
# 
# loaded via a namespace (and not attached):
# [1] utf8_1.2.4                  spatstat.explore_3.2-1      reticulate_1.40.0          
# [4] tidyselect_1.2.1            htmlwidgets_1.6.4           grid_4.1.3                 
# [7] Rtsne_0.16                  munsell_0.5.1               ragg_1.3.3                 
# [10] codetools_0.2-18            ica_1.0-2                   units_0.8-5                
# [13] future_1.25.0               miniUI_0.1.1.1              withr_3.0.2                
# [16] spatstat.random_3.1-5       colorspace_2.1-1            progressr_0.10.0           
# [19] Biobase_2.52.0              rstudioapi_0.17.1           Seurat_4.3.0.1             
# [22] stats4_4.1.3                SingleCellExperiment_1.14.1 ROCR_1.0-11                
# [25] ggsignif_0.6.2              tensor_1.5                  listenv_0.8.0              
# [28] labeling_0.4.3              MatrixGenerics_1.4.3        GenomeInfoDbData_1.2.6     
# [31] mnormt_2.1.1                polyclip_1.10-7             farver_2.1.2               
# [34] parallelly_1.31.1           vctrs_0.6.5                 generics_0.1.3             
# [37] timechange_0.3.0            R6_2.5.1                    GenomeInfoDb_1.30.1        
# [40] locfit_1.5-9.4              bitops_1.0-9                spatstat.utils_3.1-1       
# [43] DelayedArray_0.18.0         assertthat_0.2.1            promises_1.3.2             
# [46] scales_1.3.0                gtable_0.3.6                globals_0.15.0             
# [49] goftest_1.2-2               spam_2.11-0                 rlang_1.1.4                
# [52] clisymbols_1.2.0            zeallot_0.1.0               systemfonts_1.1.0          
# [55] splines_4.1.3               rstatix_0.7.0               lazyeval_0.2.2             
# [58] spatstat.geom_3.2-4         broom_1.0.7                 reshape2_1.4.4             
# [61] abind_1.4-5                 modelr_0.1.8                backports_1.5.0            
# [64] httpuv_1.6.15               tools_4.1.3                 psych_2.4.12               
# [67] RColorBrewer_1.1-3          BiocGenerics_0.38.0         ggridges_0.5.6             
# [70] Rcpp_1.0.13-1               plyr_1.8.7                  base64enc_0.1-3            
# [73] zlibbioc_1.38.0             RCurl_1.98-1.3              deldir_1.0-6               
# [76] pbapply_1.5-0               cowplot_1.1.1               S4Vectors_0.30.2           
# [79] zoo_1.8-10                  SeuratObject_5.0.2          SummarizedExperiment_1.22.0
# [82] haven_2.4.3                 ggrepel_0.9.1               cluster_2.1.2              
# [85] fs_1.6.5                    magrittr_2.0.3              magick_2.8.5               
# [88] data.table_1.16.4           scattermore_1.2             openxlsx_4.2.4             
# [91] lmtest_0.9-40               reprex_2.0.1                RANN_2.6.1                 
# [94] whisker_0.4                 fitdistrplus_1.1-5          anndata_0.7.5.6            
# [97] matrixStats_0.62.0          hms_1.1.3                   patchwork_1.3.0            
# [100] mime_0.12                   fftwtools_0.9-11            xtable_1.8-4               
# [103] rio_0.5.27                  jpeg_0.1-9                  readxl_1.3.1               
# [106] IRanges_2.26.0              gridExtra_2.3               tfruns_1.5.0               
# [109] compiler_4.1.3              KernSmooth_2.23-20          crayon_1.5.3               
# [112] SPATAData_0.0.0.9000        htmltools_0.5.8.1           later_1.4.1                
# [115] tzdb_0.4.0                  tiff_0.1-11                 lubridate_1.9.4            
# [118] DBI_1.2.3                   dbplyr_2.1.1                MASS_7.3-54                
# [121] Matrix_1.6-4                car_3.0-11                  cli_3.6.3                  
# [124] parallel_4.1.3              dotCall64_1.2               igraph_1.3.1               
# [127] GenomicRanges_1.44.0        pkgconfig_2.0.3             foreign_0.8-81             
# [130] sp_2.1-4                    confuns_1.0.3               plotly_4.10.4              
# [133] spatstat.sparse_3.0-2       xml2_1.3.2                  XVector_0.32.0             
# [136] rvest_1.0.1                 digest_0.6.37               sctransform_0.3.5          
# [139] RcppAnnoy_0.0.22            spatstat.data_3.0-1         cellranger_1.1.0           
# [142] leiden_0.3.9                uwot_0.1.16                 curl_6.0.1                 
# [145] shiny_1.10.0                EBImage_4.34.0              lifecycle_1.0.4            
# [148] nlme_3.1-152                jsonlite_1.8.9              carData_3.0-4              
# [151] viridisLite_0.4.2           pillar_1.10.0               lattice_0.20-44            
# [154] fastmap_1.2.0               httr_1.4.7                  survival_3.2-12            
# [157] glue_1.8.0                  zip_2.2.0                   png_0.1-8                  
# [160] stringi_1.8.4               textshaping_0.3.6           irlba_2.3.5.1              
# [163] future.apply_1.8.1         


# end ---------------------------------------------------------------------
