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
library(ggsci)
library(ggpubr)
library(msigdbr)


# check -------------------------------------------------------------------

packageVersion("SPATA2")
# [1] ‘2.0.4’


# dir ---------------------------------------------------------------------

dir.1 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/2_03_stRNAseq_segmentation/"
dir.2 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/2_04_stRNAseq_correlation/"
if(dir.exists(dir.2) == F){dir.create(dir.2)}


# 1. load the spata-object ------------------------------------------------

START.TIME <- Sys.time()

setwd(dir.1)
in.f <- "visium_SB28_1_subset_infiltration_area.rds"
spata.obj <- readRDS(in.f)

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 1.188813 secs


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
# [13] "segmentation"        "VRHK_MES_cl"         "RCTM_NEU_cl"         "label"              
# [17] "label.2"    

spata.obj %>% 
  getExpressionMatrixNames() %>% 
  print()
# [1] "scaled"   "denoised"

spata.obj %>% 
  getActiveMatrixName() %>% 
  print()
# [1] "denoised"


# prep --------------------------------------------------------------------

names.gs.of.interest <- c(
  "GOMF_POSTSYNAPTIC_NEUROTRANSMITTER_RECEPTOR_ACTIVITY", 
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE"
)


# prep MSigDB ------------------------------------------------------------------

msigdbr_collections() %>% print(n = Inf)

test_gs <- NULL
for(i in 1:2){
  print(i)
  
  test_gs.i <-  msigdbr(
    species = "Mus musculus", # "Homo sapiens"
    category = c("C5", "H")[i],
    subcategory = c("GO:MF", "")[i]
  ) 
  test_gs <- test_gs %>% 
    bind_rows(test_gs.i)
}


# add new geneset data ---------------------------------------------------

spata.obj.orig <- spata.obj # spata.obj <- spata.obj.orig

for(i in 1:length(names.gs.of.interest)){
  print(i)
  name.gs.i <- names.gs.of.interest[i]
  
  gs.i <- test_gs %>% 
    dplyr::filter(gs_name == name.gs.i)
  
  genes.i <- gs.i %>% 
    pull(gene_symbol)
  
  spata.obj <- spata.obj %>% 
    addGeneSet(
      class_name = "NEW2", 
      gs_name = name.gs.i, 
      check_genes = F, 
      genes = genes.i, 
      overwrite = T
    ) 
}


# check -------------------------------------------------------------------

spata.obj@used_genesets %>% 
  dplyr::filter(grepl("NEW2_", ont)) %>% 
  distinct(ont) %>% 
  print()
# distinct: removed 677 rows (99%), 4 rows remaining
# ont
# 1 NEW2_GOMF_POSTSYNAPTIC_NEUROTRANSMITTER_RECEPTOR_ACTIVITY
# 2                     NEW2_HALLMARK_TNFA_SIGNALING_VIA_NFKB
# 3                   NEW2_HALLMARK_INTERFERON_GAMMA_RESPONSE
# 4                       NEW2_HALLMARK_INFLAMMATORY_RESPONSE


# update the featureDf ----------------------------------------------------

coords_df <- getCoordsDf(spata.obj)

# join this data.frame with additional info

joined_df <- spata.obj %>%
  joinWith(
    spata_df = coords_df,
    gene_sets = paste0("NEW2_", names.gs.of.interest), # expression values of the gene set Hallmark-Hypoxia
    verbose = FALSE)


# check -------------------------------------------------------------------

joined_df %>% colnames()
# [1] "barcodes"                                                 
# [2] "x"                                                        
# [3] "y"                                                        
# [4] "tissue"                                                   
# [5] "row"                                                      
# [6] "col"                                                      
# [7] "imagerow"                                                 
# [8] "imagecol"                                                 
# [9] "sample"                                                   
# [10] "section"                                                  
# [11] "outline"                                                  
# [12] "NEW2_GOMF_POSTSYNAPTIC_NEUROTRANSMITTER_RECEPTOR_ACTIVITY"
# [13] "NEW2_HALLMARK_TNFA_SIGNALING_VIA_NFKB"                    
# [14] "NEW2_HALLMARK_INTERFERON_GAMMA_RESPONSE"                  
# [15] "NEW2_HALLMARK_INFLAMMATORY_RESPONSE"  

colnames(joined_df)[12:15] <- c("MF_P_SYN", "HM_TNFA", "HM_IFNG", "HM_INFLAM")

joined_df %>% colnames()
# [1] "barcodes"  "x"         "y"         "tissue"    "row"       "col"       "imagerow"  "imagecol"  "sample"   
# [10] "section"   "outline"   "MF_P_SYN"  "HM_TNFA"   "HM_IFNG"   "HM_INFLAM"

joined_df_simple <- joined_df %>% 
  dplyr::select(barcodes, colnames(joined_df)[12:15])


# add the feature

spata.obj <- spata.obj %>% 
  addFeatures(
    feature_names = colnames(joined_df_simple)[-1],
    feature_df = joined_df_simple, 
    overwrite = T
  )


# check -------------------------------------------------------------------

spata.obj %>% print()
# An object of class 'spata2' that contains 1 sample named 'A1_166L_Cas9'.

spata.obj %>% 
  getFeatureDf() %>% 
  colnames() %>% 
  print()
# [1] "barcodes"            "sample"              "orig.ident"          "nCount_Spatial"     
# [5] "nFeature_Spatial"    "percent.mt"          "percent.RB"          "Spatial_snn_res.0.8"
# [9] "seurat_clusters"     "histology"           "RCTM_NEU"            "VRHK_MES"           
# [13] "segmentation"        "VRHK_MES_cl"         "RCTM_NEU_cl"         "label"              
# [17] "label.2"             "MF_P_SYN"            "HM_TNFA"             "HM_IFNG"            
# [21] "HM_INFLAM"     


# save --------------------------------------------------------------------

START.TIME <- Sys.time() 

setwd(dir.2)
out.f <- "visium_SB28_1_subset_infiltration_area_edited.rds"
saveRDS(spata.obj, out.f)

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 6.917964 secs


# visualize (surface plots) (Fig. 2i) ---------------------------------------------------------------


g <- spata.obj %>% 
  plotSurface(
    display_image = F, 
    color_by = "MF_P_SYN", 
    smooth_span = 1,
    pt_size = 3
  ) + ggplot2::labs(title = "Postsynaptic Neurotransmitter\nReceptor Activity", color = "score") + 
    ggplot2::theme(
      plot.title = element_text(hjust = 0, size = 12), 
      legend.position = "right"
      )

setwd(dir.2)
out.f <- "01_surfacePlot_01_MF_P_SYN.pdf"
ggsave(out.f, g, w = 6, h = 4.5)

g1 <- g

g <- spata.obj %>% 
  plotSurface(
    display_image = F, 
    color_by = "HM_TNFA", 
    smooth_span = 1,
    pt_size = 3
  ) + ggplot2::labs(title = "TNFɑ Signaling via-NFκB", color = "score") + 
  ggplot2::theme(
    plot.title = element_text(hjust = 0, size = 12), 
    legend.position = "right"
  )

setwd(dir.2)
out.f <- "02_surfacePlot_02_HM_TNFA.pdf"
ggsave(out.f, g, w = 6, h = 4.5)

g2 <- g


g <- spata.obj %>% 
  plotSurface(
    display_image = F, 
    color_by = "HM_IFNG", 
    smooth_span = 1,
    pt_size = 3
  ) + ggplot2::labs(title = "IFNg Response", color = "score") + 
  ggplot2::theme(
    plot.title = element_text(hjust = 0, size = 12), 
    legend.position = "right"
  )

setwd(dir.2)
out.f <- "03_surfacePlot_03_HM_IFNG.pdf"
ggsave(out.f, g, w = 6, h = 4.5)

g3 <- g


g <- spata.obj %>% 
  plotSurface(
    display_image = F, 
    color_by = "HM_INFLAM", 
    smooth_span = 1,
    pt_size = 3
  ) + ggplot2::labs(title = "Inflammatory Response", color = "score") + 
  ggplot2::theme(
    plot.title = element_text(hjust = 0, size = 12), 
    legend.position = "right"
  )

setwd(dir.2)
out.f <- "04_surfacePlot_04_HM_INFLAM.pdf"
ggsave(out.f, g, w = 6, h = 4.5)

g4 <- g + theme(legend.position = "none")


# visualize (scatter plots - correlation) (Fig. 2i) ---------------------------------------------------------------

gg.list() <- list()
for(i in 1:3){
  print(i)
  gs.i <- c("HM_TNFA", "HM_IFNG", "HM_INFLAM")[i]
  
  joined_df_simple.i <- joined_df_simple %>% 
    dplyr::select(barcodes, MF_P_SYN, gs.i)
  
  colnames(joined_df_simple.i)[3] <- "marker"
  
  col.res.i <- cor.test(joined_df_simple.i$MF_P_SYN, joined_df_simple.i$marker, method = "pearson")
  TITLE.i <- paste0("Pearson p = ", format(as.numeric(col.res.i["p.value"]), digits = 3), "\n", "Corr = ", round(as.numeric(col.res.i["estimate"]), 2))
  
  g <- ggplot(joined_df_simple.i, aes(x = marker, y = MF_P_SYN) )
  g <- g + geom_smooth(alpha = 0.3, method = "lm", color = "gray20")
  g <- g + geom_point(size = 1.5, alpha = 0.8)
  g <- g + theme_classic()
  g <- g + labs(title = TITLE.i, x = gs.i, y = "MS_P_SYN")
  g <- g + theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    legend.position = "none",
    panel.grid.major = element_line(color = "#f0f0f0"),
  )
  
  I <- formatC((i + 4), width = 2, flag = "0")
  setwd(dir.2)
  out.f <- paste0(I, "_scatterPlot_cor_", I, "_", gs.i, ".pdf")
  ggsave(out.f, g, w = 6, h = 4.5)
  
  gg.list[[i]] <- g
}

g5 <- gg.list[[1]] ; g6 <- gg.list[[2]] ; g7 <- gg.list[[3]]


g <- ggarrange(
  plotlist = list(NULL, g2, g3, g4, g1, g5, g6, g7), 
  ncol = 4, nrow = 2) + theme(legend.position = "none")
  
setwd(dir.2)
out.f <- "08_surface_n_scatter_plots_combined.pdf"
ggsave(out.f, g, w = 16, h = 8)


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
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] msigdbr_7.4.1   ggpubr_0.4.0    ggsci_3.2.0     hdf5r_1.3.9     EBImage_4.34.0  tidylog_1.0.2  
# [7] forcats_0.5.1   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2     readr_2.1.5     tidyr_1.3.1    
# [13] tibble_3.2.1    ggplot2_3.5.1   tidyverse_1.3.1 SPATA2_2.0.4   
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
# [31] polyclip_1.10-7             bit64_4.5.2                 farver_2.1.2               
# [34] parallelly_1.31.1           vctrs_0.6.5                 generics_0.1.3             
# [37] timechange_0.3.0            R6_2.5.1                    GenomeInfoDb_1.28.1        
# [40] locfit_1.5-9.4              bitops_1.0-9                spatstat.utils_3.1-1       
# [43] DelayedArray_0.18.0         assertthat_0.2.1            promises_1.3.2             
# [46] scales_1.3.0                gtable_0.3.6                globals_0.15.0             
# [49] goftest_1.2-2               rlang_1.1.4                 clisymbols_1.2.0           
# [52] systemfonts_1.1.0           splines_4.1.3               rstatix_0.7.0              
# [55] lazyeval_0.2.2              spatstat.geom_3.2-4         broom_1.0.7                
# [58] reshape2_1.4.4              abind_1.4-5                 modelr_0.1.8               
# [61] backports_1.5.0             httpuv_1.6.15               tools_4.1.3                
# [64] RColorBrewer_1.1-3          BiocGenerics_0.38.0         ggridges_0.5.6             
# [67] Rcpp_1.0.13-1               plyr_1.8.7                  zlibbioc_1.38.0            
# [70] RCurl_1.98-1.3              deldir_1.0-6                pbapply_1.5-0              
# [73] cowplot_1.1.1               S4Vectors_0.30.2            zoo_1.8-10                 
# [76] SeuratObject_4.1.3          SummarizedExperiment_1.22.0 haven_2.4.3                
# [79] ggrepel_0.9.1               cluster_2.1.2               fs_1.6.5                   
# [82] magrittr_2.0.3              data.table_1.16.4           scattermore_1.2            
# [85] openxlsx_4.2.4              lmtest_0.9-40               reprex_2.0.1               
# [88] RANN_2.6.1                  fitdistrplus_1.1-5          anndata_0.7.5.6            
# [91] matrixStats_0.62.0          hms_1.1.3                   patchwork_1.3.0            
# [94] mime_0.12                   fftwtools_0.9-11            xtable_1.8-4               
# [97] rio_0.5.27                  jpeg_0.1-9                  readxl_1.3.1               
# [100] IRanges_2.26.0              gridExtra_2.3               compiler_4.1.3             
# [103] KernSmooth_2.23-20          crayon_1.5.3                SPATAData_0.0.0.9000       
# [106] htmltools_0.5.8.1           mgcv_1.8-36                 later_1.4.1                
# [109] tzdb_0.4.0                  tiff_0.1-11                 lubridate_1.9.4            
# [112] DBI_1.2.3                   dbplyr_2.1.1                MASS_7.3-54                
# [115] babelgene_21.4              Matrix_1.6-4                car_3.0-11                 
# [118] cli_3.6.3                   parallel_4.1.3              igraph_1.3.1               
# [121] GenomicRanges_1.44.0        pkgconfig_2.0.3             foreign_0.8-81             
# [124] sp_2.1-4                    confuns_1.0.3               plotly_4.10.4              
# [127] spatstat.sparse_3.0-2       xml2_1.3.2                  XVector_0.32.0             
# [130] rvest_1.0.1                 digest_0.6.37               sctransform_0.3.5          
# [133] RcppAnnoy_0.0.22            spatstat.data_3.0-1         cellranger_1.1.0           
# [136] leiden_0.3.9                uwot_0.1.16                 curl_6.0.1                 
# [139] shiny_1.10.0                lifecycle_1.0.4             nlme_3.1-152               
# [142] jsonlite_1.8.9              carData_3.0-4               viridisLite_0.4.2          
# [145] pillar_1.10.0               lattice_0.20-44             fastmap_1.2.0              
# [148] httr_1.4.7                  survival_3.2-12             glue_1.8.0                 
# [151] zip_2.2.0                   png_0.1-8                   bit_4.5.0.1                
# [154] stringi_1.8.4               textshaping_0.3.6           irlba_2.3.5.1              
# [157] future.apply_1.8.1         


# end ---------------------------------------------------------------------
