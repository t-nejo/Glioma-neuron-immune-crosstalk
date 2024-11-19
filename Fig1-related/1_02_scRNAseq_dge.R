# Takahide Nejo, MD, PhD
# Okada Lab, Department of Neurological Surgery, UC San Francisco


# note ---------------------------------------------------------------------

# Nejo T, et al., Glioma-neuronal circuit remodeling induces regional immunosuppression

# The RDS object "004_gbm_neuroimmune_v2.rds" is used as an input Seurat object. For details, please refer to the previous steps "1_00_scRNAseq_data_setup.R" and "1_01_scRNAseq_cell_type_annotation.R". 

# Using different versions of packages, such as Seurat, may lead to varying outcomes, which we have not thoroughly validated.


# rm all ------------------------------------------------------------------

rm(list = ls(all.names = T))


# packages ----------------------------------------------------------------

library(tidyverse)
library(tidylog)
library(Seurat)
library(MAST)
library(ggrepel)


# check -------------------------------------------------------------------

packageVersion("Seurat")
# [1] ‘4.0.3’


# dir ---------------------------------------------------------------------

dir.1 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/1_01_scRNAseq_cell_type_annotation/"
dir.2 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/1_02_scRNAseq_dge/"
if(dir.exists(dir.2) == F){dir.create(dir.2)}


# load files --------------------------------------------------------------

START.TIME <- Sys.time() 

setwd(dir.1) 
in.f <- "004_gbm_neuroimmune_v2.rds"
seurat.obj <- readRDS(in.f)

FINISH.TIME <- Sys.time() 

FINISH.TIME - START.TIME 
# Time difference of 27.52176 secs


# check -------------------------------------------------------------------

seurat.obj %>% print()
# An object of class Seurat 
# 47218 features across 13731 samples within 2 assays 
# Active assay: SCT (23036 features, 3000 variable features)
#  1 other assay present: RNA
#  3 dimensional reductions calculated: pca, harmony, umap

seurat.obj %>% class() %>% print() 
# [1] "Seurat"
# attr(,"package")
# [1] "Seurat"

seurat.obj %>% dim() %>% print()
# [1] 23036 13731

seurat.obj@meta.data %>% colnames() %>% print()
# [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "percent.mt"     
# [5] "nCount_SCT"      "nFeature_SCT"    "batch"           "tissue"         
# [9] "SCT_snn_res.0.2" "seurat_clusters" "clust"      

seurat.obj@meta.data %>% 
  dplyr::select(clust, tissue) %>% 
  table()
{
  # tissue
  # clust                hfc  lfc
  # Astrocytes          20   41
  # Endothelial Cells  117  163
  # Lymphoid Cells      86  309
  # Myeloid Cells      851 2924
  # Neurons            151  327
  # Pericytes          116   89
  # Tumor Cells       5325 3212
}



# subset the data ---------------------------------------------------------

# Tumor Cells; Myeloid Cells; Lymphoid Cells; Astrocytes

START.TIME <- Sys.time() 

list.seurat.obj <- list()
for(i in 1:4){
  print(i)
  
  seurat.obj.i <- seurat.obj %>% 
    subset(clust == c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells", "Astrocytes")[i])

  list.seurat.obj[[i]] <- seurat.obj.i
}

FINISH.TIME <- Sys.time() 
FINISH.TIME - START.TIME 

# Time difference of 8.38921 secs


# check -------------------------------------------------------------------

list.data <- list()
for(i in 1:4){
  print(c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells", "Astrocytes")[i])
  list.seurat.obj[[i]] %>% ncol() %>% print()
}
# [1] "Tumor Cells"
# [1] 8537
# [1] "Myeloid Cells"
# [1] 3775
# [1] "Lymphoid Cells"
# [1] 395
# [1] "Astrocytes"
# [1] 61


# Preprocessing ------------------------------------------------------------------

# remove Hemoglobin-related genes
# normalize and scale the data 

HB.GENES <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ", "MB")

START.TIME <- Sys.time() 

list.seurat.obj.2 <- list()
for(i in 1:4){
  print(i)

  seurat.obj.i <- list.seurat.obj[[i]]
  DefaultAssay(seurat.obj.i) <- "RNA"
  
  # ref: https://github.com/satijalab/seurat/issues/377
  
  counts.i <- seurat.obj.i %>% 
    GetAssayData(assay = "RNA")
  
  counts.i <- counts.i[-(which(rownames(counts.i) %in% HB.GENES)), ] 
  
  
  # remove hemoglobin-related genes
  
  seurat.obj.i <- seurat.obj.i %>% 
    subset(features = rownames(counts.i)) 
  
  
  # scale / normalize -------------------------------------------------------
  
  seurat.obj.i <- seurat.obj.i %>% 
    NormalizeData() %>% 
    ScaleData()
  
  list.seurat.obj.2[[i]] <- seurat.obj.i
}

FINISH.TIME <- Sys.time() 
FINISH.TIME - START.TIME 
# Time difference of 35.11414 secs


# check -------------------------------------------------------------------

list.data <- list()
for(i in 1:4){
  print(c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells", "Astrocytes")[i])
  list.seurat.obj.2[[i]] %>% print()
}

{
  # [1] "Tumor Cells"
  # An object of class Seurat 
  # 47200 features across 8537 samples within 2 assays 
  # Active assay: RNA (24173 features, 0 variable features)
  # 1 other assay present: SCT
  # 3 dimensional reductions calculated: pca, harmony, umap
  # [1] "Myeloid Cells"
  # An object of class Seurat 
  # 47200 features across 3775 samples within 2 assays 
  # Active assay: RNA (24173 features, 0 variable features)
  # 1 other assay present: SCT
  # 3 dimensional reductions calculated: pca, harmony, umap
  # [1] "Lymphoid Cells"
  # An object of class Seurat 
  # 47200 features across 395 samples within 2 assays 
  # Active assay: RNA (24173 features, 0 variable features)
  # 1 other assay present: SCT
  # 3 dimensional reductions calculated: pca, harmony, umap
  # [1] "Astrocytes"
  # An object of class Seurat 
  # 47200 features across 61 samples within 2 assays 
  # Active assay: RNA (24173 features, 0 variable features)
  # 1 other assay present: SCT
  # 3 dimensional reductions calculated: pca, harmony, umap
}


# DGE analysis using MAST algorithm ---------------------------------------

START.TIME <- Sys.time() 

list.res.mast <- list()
for(i in 1:4){
  print(i)
  
  seurat.obj.i <- list.seurat.obj.2[[i]]
  
  Idents(seurat.obj.i) <- "tissue" # hfc vs lfc
  
  # ref: https://github.com/ctlab/fgsea/issues/50
  
  list.res.mast[[i]] <- seurat.obj.i %>% 
    FindMarkers(
      ident.1 = "hfc", 
      ident.2 = "lfc", 
      test.use = "MAST", 
      logfc.threshold	= 0, # default = 0.25
      min.pct = 0, 
      verbose = T
    )
}

FINISH.TIME <- Sys.time() 
FINISH.TIME - START.TIME 
# Time difference of 39.67477 mins

names(list.res.mast) <- c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells", "Astrocytes")


# check -------------------------------------------------------------------

for(i in 1:4){
  print(i)
  print(c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells", "Astrocytes")[i])
  
  list.res.mast[[i]] %>% nrow() %>% print()
  list.res.mast[[i]] %>% 
    rownames_to_column("symbol") %>%
    tibble() %>% 
    print()
}

{
  # [1] 1
  # [1] "Tumor Cells"
  # [1] 24173
  # # A tibble: 24,173 × 6
  # symbol     p_val avg_log2FC pct.1 pct.2 p_val_adj
  # <chr>      <dbl>      <dbl> <dbl> <dbl>     <dbl>
  # 1 CCL3   0             -2.35  0.219 0.64  0        
  # 2 CCL4   0             -2.76  0.229 0.647 0        
  # 3 CCL4L2 3.39e-255     -1.43  0.145 0.462 8.20e-251
  # 4 XIST   6.34e-217     -1.58  0.083 0.347 1.53e-212
  # 5 POSTN  1.08e-207      1.27  0.382 0.099 2.61e-203
  # 6 GAPDH  1.38e-172      0.788 0.988 0.946 3.33e-168
  # 7 ITM2B  3.27e-160      0.616 0.892 0.857 7.91e-156
  # 8 CCL3L1 2.16e-154     -0.965 0.108 0.329 5.21e-150
  # 9 C1QB   5.35e-152     -1.17  0.165 0.404 1.29e-147
  # 10 EGFR   7.40e-149     -0.945 0.289 0.553 1.79e-144
  # # … with 24,163 more rows
  # [1] 2
  # [1] "Myeloid Cells"
  # [1] 24173
  # # A tibble: 24,173 × 6
  # symbol     p_val avg_log2FC pct.1 pct.2 p_val_adj
  # <chr>      <dbl>      <dbl> <dbl> <dbl>     <dbl>
  # 1 IFI44L 2.80e-189     -1.53  0.257 0.733 6.78e-185
  # 2 CCL4   1.29e-127     -1.24  0.716 0.966 3.11e-123
  # 3 CH25H  8.62e-127     -1.57  0.396 0.802 2.08e-122
  # 4 XIST   4.32e-116     -0.848 0.295 0.733 1.04e-111
  # 5 XAF1   3.33e-115     -0.980 0.357 0.696 8.05e-111
  # 6 LY6E   7.33e-108     -0.879 0.209 0.585 1.77e-103
  # 7 IFI6   1.21e-104     -0.974 0.311 0.668 2.93e-100
  # 8 CCL3   2.64e-104     -1.06  0.798 0.977 6.39e-100
  # 9 SAMD9L 7.31e- 92     -0.790 0.368 0.666 1.77e- 87
  # 10 PLTP   1.19e- 85      1.36  0.575 0.286 2.88e- 81
  # # … with 24,163 more rows
  # [1] 3
  # [1] "Lymphoid Cells"
  # [1] 24173
  # # A tibble: 24,173 × 6
  # symbol    p_val avg_log2FC pct.1 pct.2 p_val_adj
  # <chr>     <dbl>      <dbl> <dbl> <dbl>     <dbl>
  # 1 CCL3   1.91e-23     -1.39  0.407 0.909  4.62e-19
  # 2 CCL4   1.65e-12     -0.679 0.698 0.974  3.99e- 8
  # 3 SPP1   1.01e-11     -0.839 0.605 0.874  2.44e- 7
  # 4 APOE   1.18e-10     -1.39  0.453 0.754  2.84e- 6
  # 5 RPL17  1.68e-10      0.831 0.779 0.55   4.05e- 6
  # 6 RPS10  2.89e-10      0.818 0.919 0.929  6.99e- 6
  # 7 FTH1   2.79e- 9     -0.821 0.988 0.997  6.74e- 5
  # 8 NR4A3  2.83e- 9     -1.05  0.128 0.479  6.84e- 5
  # 9 CCL3L1 4.84e- 8     -0.973 0.105 0.411  1.17e- 3
  # 10 NR4A2  5.06e- 8     -0.827 0.651 0.877  1.22e- 3
  # # … with 24,163 more rows
  # [1] 4
  # [1] "Astrocytes"
  # [1] 24173
  # # A tibble: 24,173 × 6
  # symbol       p_val avg_log2FC pct.1 pct.2 p_val_adj
  # <chr>        <dbl>      <dbl> <dbl> <dbl>     <dbl>
  # 1 CCL4    0.00000737     -1.53   0.1  0.732     0.178
  # 2 DDX3Y   0.000110        0.996  0.75 0.195     1    
  # 3 HOPX    0.000226        1.14   0.3  0         1    
  # 4 NELFE   0.000284       -0.671  0    0.366     1    
  # 5 XIST    0.000295       -1.37   0.05 0.512     1    
  # 6 CCL3    0.000485       -1.54   0.2  0.707     1    
  # 7 RNF144B 0.000656       -1.43   0.35 0.829     1    
  # 8 ITGB1   0.000752        1.19   0.9  0.439     1    
  # 9 ERO1A   0.000806        1.62   0.6  0.317     1    
  # 10 PPCS    0.000908       -0.582  0    0.317     1    
  # # … with 24,163 more rows
}


# save ---------------------------------------------------------------------

START.TIME <- Sys.time() 

setwd(dir.2)
out.f <- "001_gbm_neuroimmune_v3_subset_list.rds"
saveRDS(list.seurat.obj.2, out.f)

out.f <- "002_res_findmarkers_mast_list.rds"
saveRDS(list.res.mast, out.f)

FINISH.TIME <- Sys.time() 
FINISH.TIME - START.TIME 
# Time difference of 5.654952 mins


# visualize MAST-DGE results in volcano plots -----------------------------------------

gg.list <- list()
for(i in 1:4){
  print(i)

  res.mast.i <- list.res.mast[[i]] %>% 
    rownames_to_column("symbol") %>%
    tibble() 
  
  NR <- res.mast.i %>% nrow() ;
  
  res.mast.i <- res.mast.i %>% 
    arrange(desc(avg_log2FC)) %>% 
    mutate(lab = c(rep("top10", 10), rep("ns", (NR - 20)), rep("bottom10", 10))) %>% 
    mutate(avg_log2FC = ifelse(avg_log2FC > 1.25, 1.25, avg_log2FC)) %>% # flattening
    mutate(avg_log2FC = ifelse(avg_log2FC < -1.25, -1.25, avg_log2FC)) %>% 
    mutate(minus.log10.padj = - log10(p_val_adj))
  
  res.mast.i$lab = factor(res.mast.i$lab, levels = c("top10", "bottom10", "ns"))
  
  
  # volcano plot ------------------------------------------------------------
  
  TITLE <- c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells", "Astrocytes")[i]
  
  g <- ggplot(res.mast.i, aes(x = avg_log2FC, y = minus.log10.padj))
  g <- g + geom_point(colour = "gray70", alpha = 0.5)
  g <- g + geom_vline(xintercept = 0)
  g <- g + geom_hline(yintercept = 0)
  g <- g + geom_text_repel(subset(res.mast.i, lab != "ns"), mapping = aes(x = avg_log2FC, y = minus.log10.padj, label = symbol), size = 2, fontface = "italic", segment.size = 0.2, show.legend = FALSE)
  g <- g + xlim(-1.25, 1.25)
  g <- g + labs(title = TITLE, x = "log2(FC)", y = "- log10 (P.adj)", colour = NULL)
  g <- g + theme(
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "none"
  )
  # plot(g)
  
  gg.list[[i]] <- g
}

g <- ggarrange(
  plotlist = gg.list, 
  ncol = 2, nrow = 2
)

setwd(dir.2)
out.f <- "003_res_findmarkers_mast_volcano.pdf"
ggsave(out.f, g, w = 12, h = 9)


# si ----------------------------------------------------------------------

Sys.time() %>% print()

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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggrepel_0.9.1      cowplot_1.1.1      glmGamPoi_1.4.0    sctransform_0.3.5 
# [5] ggpubr_0.4.0       ggsci_2.9          RColorBrewer_1.1-3 harmony_0.1.0     
# [9] Rcpp_1.0.8.3       tidylog_1.0.2      forcats_0.5.1      stringr_1.4.0     
# [13] dplyr_1.0.9        purrr_0.3.4        readr_2.1.2        tidyr_1.2.0       
# [17] tibble_3.1.7       ggplot2_3.3.6      tidyverse_1.3.1    SeuratObject_4.1.3
# [21] Seurat_4.3.0.1     SPATA2_2.0.4      
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.2.2                  spatstat.explore_3.2-1     
# [3] reticulate_1.25             tidyselect_1.1.2           
# [5] htmlwidgets_1.6.2           grid_4.1.3                 
# [7] Rtsne_0.16                  munsell_0.5.0              
# [9] codetools_0.2-18            ica_1.0-2                  
# [11] units_0.8-0                 future_1.25.0              
# [13] miniUI_0.1.1.1              withr_2.5.0                
# [15] spatstat.random_3.1-5       colorspace_2.0-3           
# [17] progressr_0.10.0            Biobase_2.52.0             
# [19] rstudioapi_0.13             stats4_4.1.3               
# [21] SingleCellExperiment_1.14.1 ROCR_1.0-11                
# [23] ggsignif_0.6.2              tensor_1.5                 
# [25] listenv_0.8.0               labeling_0.4.2             
# [27] MatrixGenerics_1.4.3        GenomeInfoDbData_1.2.6     
# [29] polyclip_1.10-0             farver_2.1.0               
# [31] parallelly_1.31.1           vctrs_0.6.1                
# [33] generics_0.1.2              timechange_0.2.0           
# [35] R6_2.5.1                    GenomeInfoDb_1.28.1        
# [37] locfit_1.5-9.4              bitops_1.0-7               
# [39] spatstat.utils_3.0-3        DelayedArray_0.18.0        
# [41] assertthat_0.2.1            promises_1.2.0.1           
# [43] scales_1.2.0                gtable_0.3.0               
# [45] globals_0.15.0              goftest_1.2-2              
# [47] rlang_1.1.0                 clisymbols_1.2.0           
# [49] systemfonts_1.0.4           splines_4.1.3              
# [51] rstatix_0.7.0               lazyeval_0.2.2             
# [53] spatstat.geom_3.2-4         broom_0.7.9                
# [55] reshape2_1.4.4              abind_1.4-5                
# [57] modelr_0.1.8                backports_1.2.1            
# [59] httpuv_1.6.5                tools_4.1.3                
# [61] ellipsis_0.3.2              BiocGenerics_0.38.0        
# [63] ggridges_0.5.3              plyr_1.8.7                 
# [65] progress_1.2.2              zlibbioc_1.38.0            
# [67] RCurl_1.98-1.3              prettyunits_1.1.1          
# [69] deldir_1.0-6                pbapply_1.5-0              
# [71] S4Vectors_0.30.2            zoo_1.8-10                 
# [73] SummarizedExperiment_1.22.0 haven_2.4.3                
# [75] cluster_2.1.2               fs_1.5.2                   
# [77] magrittr_2.0.3              data.table_1.14.2          
# [79] scattermore_0.7             openxlsx_4.2.4             
# [81] lmtest_0.9-40               reprex_2.0.1               
# [83] RANN_2.6.1                  fitdistrplus_1.1-5         
# [85] anndata_0.7.5.6             matrixStats_0.62.0         
# [87] hms_1.1.0                   patchwork_1.1.1            
# [89] mime_0.12                   fftwtools_0.9-11           
# [91] xtable_1.8-4                rio_0.5.27                 
# [93] jpeg_0.1-9                  readxl_1.3.1               
# [95] IRanges_2.26.0              gridExtra_2.3              
# [97] compiler_4.1.3              KernSmooth_2.23-20         
# [99] crayon_1.5.1                SPATAData_0.0.0.9000       
# [101] htmltools_0.5.5             later_1.3.0                
# [103] tzdb_0.1.2                  tiff_0.1-11                
# [105] lubridate_1.9.2             DBI_1.1.2                  
# [107] dbplyr_2.1.1                MASS_7.3-54                
# [109] MAST_1.18.0                 Matrix_1.6-0               
# [111] car_3.0-11                  cli_3.6.1                  
# [113] parallel_4.1.3              igraph_1.3.1               
# [115] GenomicRanges_1.44.0        pkgconfig_2.0.3            
# [117] foreign_0.8-81              sp_2.0-0                   
# [119] plotly_4.10.0               spatstat.sparse_3.0-2      
# [121] xml2_1.3.2                  XVector_0.32.0             
# [123] rvest_1.0.1                 digest_0.6.29              
# [125] RcppAnnoy_0.0.19            spatstat.data_3.0-1        
# [127] cellranger_1.1.0            leiden_0.3.9               
# [129] uwot_0.1.16                 curl_4.3.2                 
# [131] shiny_1.7.1                 EBImage_4.34.0             
# [133] lifecycle_1.0.3             nlme_3.1-152               
# [135] jsonlite_1.8.0              carData_3.0-4              
# [137] viridisLite_0.4.0           fansi_1.0.3                
# [139] pillar_1.7.0                lattice_0.20-44            
# [141] fastmap_1.1.0               httr_1.4.3                 
# [143] survival_3.2-12             glue_1.6.2                 
# [145] zip_2.2.0                   png_0.1-8                  
# [147] stringi_1.7.6               textshaping_0.3.6          
# [149] irlba_2.3.5                 future.apply_1.8.1        


# end ---------------------------------------------------------------------


