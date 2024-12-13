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
library(ggpubr)


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

names(list.seurat.obj) <- c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells", "Astrocytes")


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

names(list.seurat.obj.2) <- c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells", "Astrocytes")


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


# save ---------------------------------------------------------------------

START.TIME <- Sys.time() 

setwd(dir.2)
out.f <- "001_gbm_neuroimmune_v3_subset_list.rds"
saveRDS(list.seurat.obj.2, out.f)

FINISH.TIME <- Sys.time() 
FINISH.TIME - START.TIME 
# Time difference of 6.251748 mins


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
      latent.vars = "batch", # ref: https://github.com/satijalab/seurat/discussions/4595
      logfc.threshold	= 0, # default = 0.25
      min.pct = 0, 
      verbose = T
    )
}

FINISH.TIME <- Sys.time() 
FINISH.TIME - START.TIME 
# Time difference of 47.93051 mins

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
  # 1 CCL4   2.03e-188     -2.76  0.229 0.647 4.92e-184
  # 2 CCL3   1.52e-172     -2.35  0.219 0.64  3.67e-168
  # 3 CCL4L2 5.62e-142     -1.43  0.145 0.462 1.36e-137
  # 4 S100A1 4.67e-136     -0.944 0.321 0.415 1.13e-131
  # 5 HIF1A  1.13e-131      0.800 0.788 0.699 2.72e-127
  # 6 PLA2R1 1.46e-125     -0.357 0.11  0.286 3.52e-121
  # 7 FKBP1A 2.64e-118      0.798 0.825 0.712 6.38e-114
  # 8 COTL1  1.51e-117      0.782 0.729 0.609 3.64e-113
  # 9 STK32A 6.22e-111      0.599 0.449 0.366 1.50e-106
  # 10 COX4I1 1.60e-107      0.497 0.868 0.809 3.87e-103
  # # … with 24,163 more rows
  # [1] 2
  # [1] "Myeloid Cells"
  # [1] 24173
  # # A tibble: 24,173 × 6
  # symbol    p_val avg_log2FC pct.1 pct.2 p_val_adj
  # <chr>     <dbl>      <dbl> <dbl> <dbl>     <dbl>
  # 1 IFI44L 4.07e-80     -1.53  0.257 0.733  9.84e-76
  # 2 CCL4   1.02e-65     -1.24  0.716 0.966  2.46e-61
  # 3 RPS10  2.34e-56      0.358 0.853 0.908  5.65e-52
  # 4 CCL3   1.30e-48     -1.06  0.798 0.977  3.15e-44
  # 5 IFI6   6.53e-47     -0.974 0.311 0.668  1.58e-42
  # 6 RPL17  3.12e-43      0.376 0.66  0.657  7.54e-39
  # 7 IL1B   6.49e-41     -1.14  0.441 0.752  1.57e-36
  # 8 LY6E   1.86e-40     -0.879 0.209 0.585  4.50e-36
  # 9 CD83   2.01e-39     -0.804 0.76  0.94   4.87e-35
  # 10 CCL4L2 9.17e-39     -0.590 0.566 0.852  2.22e-34
  # # … with 24,163 more rows
  # [1] 3
  # [1] "Lymphoid Cells"
  # [1] 24173
  # # A tibble: 24,173 × 6
  # symbol    p_val avg_log2FC pct.1 pct.2 p_val_adj
  # <chr>     <dbl>      <dbl> <dbl> <dbl>     <dbl>
  # 1 CCL3   3.87e-20     -1.39  0.407 0.909  9.36e-16
  # 2 RPS10  1.00e-12      0.818 0.919 0.929  2.43e- 8
  # 3 CCL4   2.04e-11     -0.679 0.698 0.974  4.93e- 7
  # 4 RPL17  3.09e-10      0.831 0.779 0.55   7.48e- 6
  # 5 NR4A3  6.25e- 9     -1.05  0.128 0.479  1.51e- 4
  # 6 SPP1   1.01e- 8     -0.839 0.605 0.874  2.44e- 4
  # 7 DUSP1  5.59e- 8      0.435 0.791 0.939  1.35e- 3
  # 8 CCL3L1 6.07e- 8     -0.973 0.105 0.411  1.47e- 3
  # 9 APOE   9.25e- 8     -1.39  0.453 0.754  2.24e- 3
  # 10 FTH1   1.16e- 7     -0.821 0.988 0.997  2.82e- 3
  # # … with 24,163 more rows
  # [1] 4
  # [1] "Astrocytes"
  # [1] 24173
  # # A tibble: 24,173 × 6
  # symbol     p_val avg_log2FC pct.1 pct.2 p_val_adj
  # <chr>      <dbl>      <dbl> <dbl> <dbl>     <dbl>
  # 1 NELFE   0.000277     -0.671  0    0.366         1
  # 2 ITGA5   0.000556      1.13   0.35 0.268         1
  # 3 PPCS    0.000557     -0.582  0    0.317         1
  # 4 SMIM14  0.000716     -1.09   0.05 0.439         1
  # 5 PARVB   0.000888     -0.453  0    0.171         1
  # 6 UBE2E1  0.00100      -0.566  0.05 0.293         1
  # 7 MBTPS1  0.00131      -0.533  0    0.195         1
  # 8 TRIM52  0.00143      -0.597  0    0.171         1
  # 9 PAIP1   0.00149      -0.350  0    0.195         1
  # 10 TMEM173 0.00149      -0.484  0    0.195         1
  # # … with 24,163 more rows
}


# save ---------------------------------------------------------------------

START.TIME <- Sys.time() 

out.f <- "002_res_findmarkers_mast_list.rds"
saveRDS(list.res.mast, out.f)

FINISH.TIME <- Sys.time() 
FINISH.TIME - START.TIME 
# Time difference of 0.2180417 secs


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
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] ggrepel_0.9.1      cowplot_1.1.1      glmGamPoi_1.4.0    sctransform_0.3.5 
# [5] ggpubr_0.4.0       ggsci_2.9          RColorBrewer_1.1-3 harmony_0.1.0     
# [9] Rcpp_1.0.8.3       tidylog_1.0.2      forcats_0.5.1      stringr_1.4.0     
# [13] dplyr_1.0.9        purrr_0.3.4        readr_2.1.2        tidyr_1.2.0       
# [17] tibble_3.1.7       ggplot2_3.3.6      tidyverse_1.3.1    SeuratObject_4.1.3
# [21] Seurat_4.3.0.1     SPATA2_2.0.4      
# 
# loaded via a namespace (and not attached):
# [1] utf8_1.2.2                  spatstat.explore_3.2-1     
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


