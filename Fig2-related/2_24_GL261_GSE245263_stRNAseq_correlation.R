# Takahide Nejo, MD, PhD
# Okada Lab, Department of Neurological Surgery, UC San Francisco


# note ---------------------------------------------------------------------

# Nejo T, et al., Glioma-neuronal circuit remodeling induces regional immunosuppression
# Analysis of SB28 tumor-bearing mouse brain 10X Visium st-RNA-seq data

# Below is the script for the first sample. The data from the second sample was processed and analyzed similarly. 

# Using different versions of packages, such as SPATA2, may lead to varying outcomes, which we have not thoroughly validated.

# The RDS object "visium_GSE245263_GL261_subset_infiltration_area.rds" is used as input. For more details, please refer to the previous step, "2_23_GL261_GSE245263_stRNAseq_segmentation.R".


# rm all ------------------------------------------------------------------

rm(list = ls(all.names = T))


# packages ----------------------------------------------------------------

library(SPATA2) 
library(tidyverse)
library(tidylog)
library(ggsci)
library(ggpubr)
library(msigdbr)


# check -------------------------------------------------------------------

packageVersion("SPATA2")
# [1] ‘2.0.4’

packageVersion("msigdbr")
# [1] ‘7.4.1’


# dir ---------------------------------------------------------------------

dir.1 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/2_23_GL261_GSE245263_stRNAseq_segmentation/"
dir.2 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/2_24_GL261_GSE245263_stRNAseq_correlation/"
if(dir.exists(dir.2) == F){dir.create(dir.2)}


# 1. load the spata-object ------------------------------------------------

START.TIME <- Sys.time()

setwd(dir.1)
in.f <- "visium_GSE245263_GL261_subset_infiltration_area.rds"
spata.obj <- readRDS(in.f)

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 2.410955 secs


# check -------------------------------------------------------------------

spata.obj 
# An object of class 'spata2' that contains 1 sample named 'GSM7839621'.

spata.obj %>% 
  getFeatureDf() %>% 
  colnames() %>% 
  print()
# [1] "barcodes"                    "sample"                      "in_tissue"                  
# [4] "array_row"                   "array_col"                   "n_genes_by_counts"          
# [7] "log1p_n_genes_by_counts"     "total_counts"                "log1p_total_counts"         
# [10] "pct_counts_in_top_50_genes"  "pct_counts_in_top_100_genes" "pct_counts_in_top_200_genes"
# [13] "pct_counts_in_top_500_genes" "total_counts_mt"             "log1p_total_counts_mt"      
# [16] "pct_counts_mt"               "condition"                   "coarse"                     
# [19] "fine"                        "leiden"                      "segmentation"               
# [22] "VRHK_MES"                    "RCTM_NEU"                    "VRHK_MES_cl"                
# [25] "RCTM_NEU_cl"                 "label"                       "label.2"   

spata.obj %>% 
  getExpressionMatrixNames() %>% 
  print()
# [1] "normalized" "scaled"     "denoised"  

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


# check / edit -------------------------------------------------------------------

joined_df %>% colnames()
# [1] "barcodes"                                                 
# [2] "sample"                                                   
# [3] "imagerow"                                                 
# [4] "imagecol"                                                 
# [5] "x"                                                        
# [6] "y"                                                        
# [7] "section"                                                  
# [8] "outline"                                                  
# [9] "NEW2_GOMF_POSTSYNAPTIC_NEUROTRANSMITTER_RECEPTOR_ACTIVITY"
# [10] "NEW2_HALLMARK_TNFA_SIGNALING_VIA_NFKB"                    
# [11] "NEW2_HALLMARK_INTERFERON_GAMMA_RESPONSE"                  
# [12] "NEW2_HALLMARK_INFLAMMATORY_RESPONSE"     

colnames(joined_df)[9:12] <- c("MF_P_SYN", "HM_TNFA", "HM_IFNG", "HM_INFLAM")

joined_df %>% colnames()
# [1] "barcodes"  "sample"    "imagerow"  "imagecol"  "x"         "y"         "section"   "outline"  
# [9] "MF_P_SYN"  "HM_TNFA"   "HM_IFNG"   "HM_INFLAM"

joined_df_simple <- joined_df %>% 
  dplyr::select(barcodes, colnames(joined_df)[9:12])

# add the features

spata.obj <- spata.obj %>% 
  addFeatures(
    feature_names = colnames(joined_df_simple)[-1],
    feature_df = joined_df_simple, 
    overwrite = T
  )


# check -------------------------------------------------------------------

spata.obj %>% print()
# An object of class 'spata2' that contains 1 sample named 'GSM7839621'.

spata.obj %>% 
  getFeatureDf() %>% 
  colnames() %>% 
  print()
# [1] "barcodes"                    "sample"                      "in_tissue"                  
# [4] "array_row"                   "array_col"                   "n_genes_by_counts"          
# [7] "log1p_n_genes_by_counts"     "total_counts"                "log1p_total_counts"         
# [10] "pct_counts_in_top_50_genes"  "pct_counts_in_top_100_genes" "pct_counts_in_top_200_genes"
# [13] "pct_counts_in_top_500_genes" "total_counts_mt"             "log1p_total_counts_mt"      
# [16] "pct_counts_mt"               "condition"                   "coarse"                     
# [19] "fine"                        "leiden"                      "segmentation"               
# [22] "VRHK_MES"                    "RCTM_NEU"                    "VRHK_MES_cl"                
# [25] "RCTM_NEU_cl"                 "label"                       "label.2"                    
# [28] "MF_P_SYN"                    "HM_TNFA"                     "HM_IFNG"                    
# [31] "HM_INFLAM"      


# save --------------------------------------------------------------------

START.TIME <- Sys.time() 

setwd(dir.2)
out.f <- "visium_GSE245263_GL261_subset_infiltration_area_edited.rds"
saveRDS(spata.obj, out.f)

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 19.53429 secs


# visualize (surface plots) (Extended Data Fig. S7n) ---------------------------------------------------------------

gg.list <- list()
for(i in 1:4){
  print(i)
  name.i <- c("MF_P_SYN", "HM_TNFA", "HM_IFNG", "HM_INFLAM")[i]
  title.i <- c("Postsynaptic Neurotransmitter\nReceptor Activity", 
               "TNFɑ Signaling via-NFκB", 
               "IFNg Response",
               "Inflammatory Response"
  )[i]
  g <- spata.obj %>% 
    plotSurface(
      display_image = F, 
      color_by = name.i, 
      smooth_span = 1,
      pt_size = 3
    ) + ggplot2::labs(title = title.i, color = "score") + 
    ggplot2::theme(
      plot.title = element_text(hjust = 0, size = 12), 
      legend.position = "right"
    )
  
  I <- formatC(i, width = 2, flag = "0")
  setwd(dir.2)
  out.f <- paste0(I, "_surfacePlot_", I, "_", name.i, ".pdf")
  ggsave(out.f, g, w = 6, h = 4.5)
  
  gg.list[[i]] <- g + theme(legend.position = "none")
}

g1 <- gg.list[[1]] ; g2 <- gg.list[[2]] ; g3 <- gg.list[[3]]; g4 <- gg.list[[4]]


# visualize (scatter plots - correlation) (Fig. 2i) ---------------------------------------------------------------

gg.list <- list()
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
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] msigdbr_7.4.1   ggpubr_0.4.0    ggsci_3.2.0     tidylog_1.0.2   forcats_0.5.1   stringr_1.5.1  
# [7] dplyr_1.1.4     purrr_1.0.2     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1  
# [13] tidyverse_1.3.1 SPATA2_2.0.4   
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
# [31] polyclip_1.10-7             farver_2.1.2                parallelly_1.31.1          
# [34] vctrs_0.6.5                 generics_0.1.3              timechange_0.3.0           
# [37] R6_2.5.1                    GenomeInfoDb_1.30.1         locfit_1.5-9.4             
# [40] bitops_1.0-9                spatstat.utils_3.1-1        DelayedArray_0.18.0        
# [43] assertthat_0.2.1            promises_1.3.2              scales_1.3.0               
# [46] gtable_0.3.6                globals_0.15.0              goftest_1.2-2              
# [49] spam_2.11-0                 rlang_1.1.4                 clisymbols_1.2.0           
# [52] systemfonts_1.1.0           splines_4.1.3               rstatix_0.7.0              
# [55] lazyeval_0.2.2              spatstat.geom_3.2-4         broom_1.0.7                
# [58] reshape2_1.4.4              abind_1.4-5                 modelr_0.1.8               
# [61] backports_1.5.0             httpuv_1.6.15               tools_4.1.3                
# [64] RColorBrewer_1.1-3          BiocGenerics_0.38.0         ggridges_0.5.6             
# [67] Rcpp_1.0.13-1               plyr_1.8.7                  zlibbioc_1.38.0            
# [70] RCurl_1.98-1.3              deldir_1.0-6                pbapply_1.5-0              
# [73] cowplot_1.1.1               S4Vectors_0.30.2            zoo_1.8-10                 
# [76] SeuratObject_5.0.2          SummarizedExperiment_1.22.0 haven_2.4.3                
# [79] ggrepel_0.9.1               cluster_2.1.2               fs_1.6.5                   
# [82] magrittr_2.0.3              magick_2.8.5                data.table_1.16.4          
# [85] scattermore_1.2             openxlsx_4.2.4              lmtest_0.9-40              
# [88] reprex_2.0.1                RANN_2.6.1                  fitdistrplus_1.1-5         
# [91] anndata_0.7.5.6             matrixStats_0.62.0          hms_1.1.3                  
# [94] patchwork_1.3.0             mime_0.12                   fftwtools_0.9-11           
# [97] xtable_1.8-4                rio_0.5.27                  jpeg_0.1-9                 
# [100] readxl_1.3.1                IRanges_2.26.0              gridExtra_2.3              
# [103] compiler_4.1.3              KernSmooth_2.23-20          crayon_1.5.3               
# [106] SPATAData_0.0.0.9000        htmltools_0.5.8.1           mgcv_1.8-36                
# [109] later_1.4.1                 tzdb_0.4.0                  tiff_0.1-11                
# [112] lubridate_1.9.4             DBI_1.2.3                   dbplyr_2.1.1               
# [115] MASS_7.3-54                 babelgene_21.4              Matrix_1.6-4               
# [118] car_3.0-11                  cli_3.6.3                   parallel_4.1.3             
# [121] dotCall64_1.2               igraph_1.3.1                GenomicRanges_1.44.0       
# [124] pkgconfig_2.0.3             foreign_0.8-81              sp_2.1-4                   
# [127] confuns_1.0.3               plotly_4.10.4               spatstat.sparse_3.0-2      
# [130] xml2_1.3.2                  XVector_0.32.0              rvest_1.0.1                
# [133] digest_0.6.37               sctransform_0.3.5           RcppAnnoy_0.0.22           
# [136] spatstat.data_3.0-1         cellranger_1.1.0            leiden_0.3.9               
# [139] uwot_0.1.16                 curl_6.0.1                  shiny_1.10.0               
# [142] EBImage_4.34.0              lifecycle_1.0.4             nlme_3.1-152               
# [145] jsonlite_1.8.9              carData_3.0-4               viridisLite_0.4.2          
# [148] pillar_1.10.0               lattice_0.20-44             fastmap_1.2.0              
# [151] httr_1.4.7                  survival_3.2-12             glue_1.8.0                 
# [154] zip_2.2.0                   png_0.1-8                   stringi_1.8.4              
# [157] textshaping_0.3.6           irlba_2.3.5.1               future.apply_1.8.1         


# end ---------------------------------------------------------------------
