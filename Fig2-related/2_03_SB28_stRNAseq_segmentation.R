# Takahide Nejo, MD, PhD
# Okada Lab, Department of Neurological Surgery, UC San Francisco


# note ---------------------------------------------------------------------

# Nejo T, et al., Glioma-neuronal circuit remodeling induces regional immunosuppression
# Analysis of SB28 tumor-bearing mouse brain 10X Visium st-RNA-seq data

# Below is the script for the first sample. The data from the second sample was processed and analyzed similarly. 

# Using different versions of packages, such as SPATA2, may lead to varying outcomes, which we have not thoroughly validated.

# The RDS object "visium_SB28_1_ftr_edited.rds" is used as input. For more details, please refer to the previous step, "2_02_stRNAseq_overview.R".


# rm all ------------------------------------------------------------------

rm(list = ls(all.names = T))


# packages ----------------------------------------------------------------

library(SPATA2) 
library(tidyverse)
library(tidylog)
library(EBImage)
library(hdf5r)
library(ggsci)
library(ggpubr)


# check -------------------------------------------------------------------

packageVersion("SPATA2")
# [1] ‘2.0.4’


# dir ---------------------------------------------------------------------

dir.1 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/2_02_stRNAseq_overview/"
dir.2 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/2_03_stRNAseq_segmentation/"
if(dir.exists(dir.2) == F){dir.create(dir.2)}


# 1. load the spata-object ------------------------------------------------

START.TIME <- Sys.time()

setwd(dir.1)
in.f <- "visium_SB28_1_ftr_edited.rds"
spata.obj <- readRDS(in.f)

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 6.167014 secs


# check -------------------------------------------------------------------

spata.obj 
# An object of class 'spata2' that contains 1 sample named 'A1_166L_Cas9'.

spata.obj %>% 
  getFeatureDf() %>% 
  colnames() %>% 
  print()
# [1] "barcodes"            "sample"              "orig.ident"          "nCount_Spatial"     
# [5] "nFeature_Spatial"    "percent.mt"          "percent.RB"          "Spatial_snn_res.0.8"
# [9] "seurat_clusters"     "histology"           "RCTM_NEU"            "VRHK_MES"  

spata.obj %>% 
  getExpressionMatrixNames() %>% 
  print()
# [1] "scaled"   "denoised"

spata.obj %>% 
  getActiveMatrixName() %>% 
  print()
# [1] "denoised"



# add new categorical values ----------------------------------------------

gg.list <- list()
for(i in 1:2){
  print(i)
  FEATURE <- c("VRHK_MES", "RCTM_NEU")[i]
  
  spata.obj.i <- spata.obj
  
  feature.df.i <- spata.obj.i %>% 
    getFeatureDf() %>% 
    dplyr::select(barcodes, marker = FEATURE)
  
  percentile <- feature.df.i %>% 
    pull(marker) %>% 
    quantile(probs = seq(0, 1, 0.1)) %>% 
    as.numeric()
  
  if(i == 1){
    PCTL.MES <- c(percentile[8], percentile[10])
  }else{
    PCTL.NEU <- c(percentile[2], percentile[4])
  }
  
  feature.df.i$marker.cl <- 0
  
  for(j in 1:10){
    if(j <= 9){
      feature.df.i <- feature.df.i %>% 
        mutate(marker.cl = ifelse(percentile[j] <= marker & marker < percentile[(j + 1)], j, marker.cl))      
    }else if(j == 10){
      feature.df.i <- feature.df.i %>% 
        mutate(marker.cl = ifelse(percentile[j] <= marker & marker <= percentile[(j + 1)], j, marker.cl))        
    }
    
  }
  feature.df.i <- feature.df.i %>% 
    mutate(marker.cl = factor(marker.cl))
  
  # check
  feature.df.i %>% 
    pull(marker.cl) %>% 
    table() %>% 
    print()
  # .
  # 1   2   3   4   5   6   7   8   9  10 
  # 198 197 198 197 197 198 197 198 197 198 
  
  
  # vizualize
  
  g <- ggplot(feature.df.i, aes(x = marker.cl, y = marker, color = marker.cl))
  g <- g + geom_point(size = 0.75, position = position_jitter(width =.1))
  g <- g + labs(title = FEATURE, x = NULL)
  g <- g + scale_color_npg()
  g <- g + theme_classic()
  g <- g + theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    axis.text = element_blank(), 
    axis.title = element_blank(), 
    legend.position = "none",
    panel.grid.major = element_line(color = "#f0f0f0"), 
  )
  g <- g + coord_flip()
  
  gg.list[[i]] <- g
  
  FEATURE.cl <- paste0(FEATURE, "_cl")
  feature.df.i <- feature.df.i %>% 
    dplyr::select(barcodes, marker.cl)
  colnames(feature.df.i)[2] <- FEATURE.cl
  
  # add the feature
  spata.obj <- spata.obj %>% 
    addFeatures(
      feature_names = FEATURE.cl,
      feature_df = feature.df.i, 
      overwrite = T
    )
}

# update the featureDf ---------------------------------------------

feature.df <- spata.obj %>%
  getFeatureDf() %>% 
  mutate(label = ifelse( (VRHK_MES_cl == 8 | VRHK_MES_cl == 9) & (RCTM_NEU_cl == 2 | RCTM_NEU_cl == 3), "yes", "no")) %>% 
  dplyr::select(barcodes, label)

# add the feature
spata.obj <- spata.obj %>% 
  addFeatures(
    feature_names = "label",
    feature_df = feature.df, 
    overwrite = T
  )


# visualize (Fig 2g) ---------------------------------------------------------------

feature.df <- spata.obj %>%
  getFeatureDf()

g <- ggplot()
g <- g + geom_point(feature.df, mapping = aes(x = VRHK_MES, y = RCTM_NEU), color = "gray90", alpha = 0.5)
g <- g + geom_point(feature.df %>% dplyr::filter((VRHK_MES_cl == 8 | VRHK_MES_cl == 9) & (RCTM_NEU_cl == 2 | RCTM_NEU_cl == 3) ), mapping = aes(x = VRHK_MES, y = RCTM_NEU), color = "deepskyblue3")
g <- g + geom_vline(xintercept = PCTL.MES, linetype = 2)
g <- g + geom_hline(yintercept = PCTL.NEU, linetype = 2)
g <- g + theme_classic()
g <- g + theme(
  axis.text = element_blank(), 
  axis.title = element_blank(), 
  legend.position = "none",
  panel.grid.major = element_line(color = "#f0f0f0"), 
)
gg.list[[3]] <- g

g <- ggarrange(
  plotlist = gg.list[c(2,1,3)], 
  ncol = 1, nrow = 3, 
  heights = c(0.2, 0.2, 0.6)
)
g

setwd(dir.2)
out.f <- "01_scatterPlot_vrhk_mes_n_rctm_neu.pdf"
ggsave(out.f, g, w = 4.5, h = 6)


# add segmentation info (interactive) ---------------------------------------------------

# draw a rough outline of the "peritumoral region"

spata.obj <- spata.obj %>%
  createSegmentation()


# visualize (check segmentation) ---------------------------------------------------------------

g <- spata.obj %>% 
  plotSurface(
    color_by = "segmentation", 
    display_image = T, 
    pt_size = 1.5, 
  ) + ggplot2::labs(title = "segmentation") + 
  ggplot2::scale_color_manual(values = c("gray90", "goldenrod2")) + 
  ggplot2::scale_alpha_manual(values = c(0, 0.75)) +
  ggplot2::theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  ggplot2::theme(legend.position = "right")

setwd(dir.2)
out.f <- "02_surfaceplot_segmentation.pdf"
ggsave(out.f, g, w = 6, h = 4.5)


# update the featureDf ---------------------------------------------

feature.df <- spata.obj %>%
  getFeatureDf() %>% 
  mutate(label.2 = ifelse(label == "yes" & segmentation == "peritumoral_area", "yes", "no")) %>% 
  dplyr::select(barcodes, label.2) 

# add the feature
spata.obj <- spata.obj %>% 
  addFeatures(
    feature_names = "label.2",
    feature_df = feature.df, 
    overwrite = T
  )


# check -------------------------------------------------------------------

feature.df %>% 
  pull(label.2) %>% 
  table() %>% 
  print()


# "Putative glioma-neuronal infiltration area" Fig 2h ------------------------------------------------------------------

g <- plotSurface(
  object = spata.obj, 
  display_image = F, 
  color_by = "label.2", 
  pt_size = 1.5, 
  pt_alpha = 0.9, 
) + 
  ggplot2::scale_color_manual(values = c("gray90", "deepskyblue3")) + 
  ggplot2::theme(
    legend.position = "none"
  )

setwd(dir.2)
out.f <- "03_surfaceplot_infiltration_area.pdf"
ggsave(out.f, g, w = 6, h = 4.5)


# save --------------------------------------------------------------------

START.TIME <- Sys.time() 

setwd(dir.2)
out.f <- "visium_SB28_1_segmentation.rds"
saveRDS(spata.obj, out.f)

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 1.160291 mins


# subset ------------------------------------------------------------------

barcodes.to.retain <- spata.obj %>%
  getFeatureDf() %>%
  dplyr::filter(label.2 == "yes") %>%
  pull(barcodes)

spata.obj.subset <- spata.obj %>% 
  subsetByBarcodes(
    barcodes = barcodes.to.retain,
    verbose = TRUE
  )


# save --------------------------------------------------------------------

START.TIME <- Sys.time() 

setwd(dir.2)
out.f <- "visium_SB28_1_subset_infiltration_area.rds"
saveRDS(spata.obj.subset, out.f)

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 7.18751 secs


# si ----------------------------------------------------------------------

Sys.time()
# 

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
# [1] shiny_1.10.0    ggpubr_0.4.0    ggsci_3.2.0     hdf5r_1.3.9     EBImage_4.34.0 
# [6] tidylog_1.0.2   forcats_0.5.1   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2    
# [11] readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_1.3.1
# [16] SPATA2_2.0.4   
# 
# loaded via a namespace (and not attached):
# [1] shinydashboard_0.7.2        utf8_1.2.4                  spatstat.explore_3.2-1     
# [4] reticulate_1.40.0           tidyselect_1.2.1            htmlwidgets_1.6.4          
# [7] grid_4.1.3                  Rtsne_0.16                  munsell_0.5.1              
# [10] ragg_1.3.3                  codetools_0.2-18            ica_1.0-2                  
# [13] units_0.8-5                 DT_0.33                     future_1.25.0              
# [16] miniUI_0.1.1.1              withr_3.0.2                 spatstat.random_3.1-5      
# [19] colorspace_2.1-1            progressr_0.10.0            Biobase_2.52.0             
# [22] rstudioapi_0.17.1           Seurat_4.3.0.1              stats4_4.1.3               
# [25] SingleCellExperiment_1.14.1 ROCR_1.0-11                 ggsignif_0.6.2             
# [28] tensor_1.5                  shinyWidgets_0.8.7          listenv_0.8.0              
# [31] MatrixGenerics_1.4.3        labeling_0.4.3              GenomeInfoDbData_1.2.6     
# [34] polyclip_1.10-7             bit64_4.5.2                 farver_2.1.2               
# [37] parallelly_1.31.1           vctrs_0.6.5                 generics_0.1.3             
# [40] timechange_0.3.0            R6_2.5.1                    GenomeInfoDb_1.28.1        
# [43] shinybusy_0.3.3             locfit_1.5-9.4              concaveman_1.1.0           
# [46] cachem_1.1.0                bitops_1.0-9                spatstat.utils_3.1-1       
# [49] DelayedArray_0.18.0         assertthat_0.2.1            promises_1.3.2             
# [52] scales_1.3.0                gtable_0.3.6                keys_0.1.1                 
# [55] globals_0.15.0              goftest_1.2-2               rlang_1.1.4                
# [58] clisymbols_1.2.0            systemfonts_1.1.0           splines_4.1.3              
# [61] rstatix_0.7.0               lazyeval_0.2.2              spatstat.geom_3.2-4        
# [64] broom_1.0.7                 reshape2_1.4.4              abind_1.4-5                
# [67] modelr_0.1.8                backports_1.5.0             httpuv_1.6.15              
# [70] tools_4.1.3                 jquerylib_0.1.4             RColorBrewer_1.1-3         
# [73] BiocGenerics_0.38.0         ggridges_0.5.6              Rcpp_1.0.13-1              
# [76] plyr_1.8.7                  zlibbioc_1.38.0             RCurl_1.98-1.3             
# [79] dbscan_1.2-0                deldir_1.0-6                pbapply_1.5-0              
# [82] cowplot_1.1.1               fontawesome_0.5.3           S4Vectors_0.30.2           
# [85] zoo_1.8-10                  SeuratObject_4.1.3          SummarizedExperiment_1.22.0
# [88] haven_2.4.3                 ggrepel_0.9.1               cluster_2.1.2              
# [91] fs_1.6.5                    magrittr_2.0.3              magick_2.8.5               
# [94] data.table_1.16.4           scattermore_1.2             openxlsx_4.2.4             
# [97] lmtest_0.9-40               reprex_2.0.1                RANN_2.6.1                 
# [100] fitdistrplus_1.1-5          anndata_0.7.5.6             matrixStats_0.62.0         
# [103] hms_1.1.3                   patchwork_1.3.0             mime_0.12                  
# [106] fftwtools_0.9-11            xtable_1.8-4                rio_0.5.27                 
# [109] jpeg_0.1-9                  readxl_1.3.1                IRanges_2.26.0             
# [112] gridExtra_2.3               shinyhelper_0.3.2           compiler_4.1.3             
# [115] V8_6.0.0                    KernSmooth_2.23-20          crayon_1.5.3               
# [118] SPATAData_0.0.0.9000        htmltools_0.5.8.1           later_1.4.1                
# [121] tzdb_0.4.0                  tiff_0.1-11                 lubridate_1.9.4            
# [124] DBI_1.2.3                   dbplyr_2.1.1                MASS_7.3-54                
# [127] babelgene_21.4              Matrix_1.6-4                car_3.0-11                 
# [130] cli_3.6.3                   parallel_4.1.3              igraph_1.3.1               
# [133] GenomicRanges_1.44.0        pkgconfig_2.0.3             foreign_0.8-81             
# [136] sp_2.1-4                    confuns_1.0.3               plotly_4.10.4              
# [139] spatstat.sparse_3.0-2       xml2_1.3.2                  bslib_0.8.0                
# [142] XVector_0.32.0              rvest_1.0.1                 digest_0.6.37              
# [145] sctransform_0.3.5           RcppAnnoy_0.0.22            spatstat.data_3.0-1        
# [148] cellranger_1.1.0            leiden_0.3.9                uwot_0.1.16                
# [151] curl_6.0.1                  lifecycle_1.0.4             nlme_3.1-152               
# [154] jsonlite_1.8.9              carData_3.0-4               viridisLite_0.4.2          
# [157] pillar_1.10.0               lattice_0.20-44             fastmap_1.2.0              
# [160] httr_1.4.7                  survival_3.2-12             glue_1.8.0                 
# [163] zip_2.2.0                   png_0.1-8                   bit_4.5.0.1                
# [166] sass_0.4.9                  stringi_1.8.4               textshaping_0.3.6          
# [169] memoise_2.0.1               irlba_2.3.5.1               future.apply_1.8.1         


# end ---------------------------------------------------------------------
