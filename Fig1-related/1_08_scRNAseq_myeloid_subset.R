# Takahide Nejo, MD, PhD
# Okada Lab, Department of Neurological Surgery, UC San Francisco


# note ---------------------------------------------------------------------

# Nejo T, et al., Glioma-neuronal circuit remodeling induces regional immunosuppression

# Using different versions of packages, such as Seurat, may lead to varying outcomes, which we have not thoroughly validated.

# The RDS object "001_gbm_neuroimmune_v3_subset_list.rds" is used as inputs. For details, please refer to the previous step "1_02_scRNAseq_dge.R". 


# rm all ------------------------------------------------------------------

rm(list = ls(all.names = T))


# packages ----------------------------------------------------------------

library(tidyverse)
library(tidylog)
library(Seurat)
library(glmGamPoi) ; # BiocManager::install("glmGamPoi")
library(harmony)
library(clustree)


# check -------------------------------------------------------------------

packageVersion("Seurat")
# [1] ‘4.3.0.1’
packageVersion("sctransform")
# [1] ‘0.3.5’


# dir ---------------------------------------------------------------------

dir.1 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/1_02_scRNAseq_dge/"
dir.2 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/1_08_scRNAseq_myeloid_subset/"
if(dir.exists(dir.2) == F){dir.create(dir.2)}


# load file(s) --------------------------------------------------------------

START.TIME <- Sys.time() 

setwd(dir.1) 
in.f <- "001_gbm_neuroimmune_v3_subset_list.rds"
list.seurat.obj <- readRDS(in.f)

FINISH.TIME <- Sys.time() 

FINISH.TIME - START.TIME 
# Time difference of 50.48816 secs


# check -------------------------------------------------------------------

list.seurat.obj %>% length()
# [1] 4

list.seurat.obj %>% names()
# [1] "Tumor Cells"    "Myeloid Cells"  "Lymphoid Cells" "Astrocytes"    

list.seurat.obj[["Myeloid Cells"]]
# An object of class Seurat 
# 47200 features across 3775 samples within 2 assays 
# Active assay: RNA (24173 features, 0 variable features)
# 1 other assay present: SCT
# 3 dimensional reductions calculated: pca, harmony, umap


# retain myeloid data only ------------------------------------------------------------------

seurat.obj <- list.seurat.obj[["Myeloid Cells"]]


# check -------------------------------------------------------------------

seurat.obj %>% print()
# An object of class Seurat 
# 47200 features across 3775 samples within 2 assays 
# Active assay: RNA (24173 features, 0 variable features)
# 1 other assay present: SCT
# 3 dimensional reductions calculated: pca, harmony, umap

seurat.obj %>% dim() %>% print()
# [1] 24173  3775


# retain the RNA assay data only  

DefaultAssay(seurat.obj) <- "RNA"
seurat.obj$SCT <- NULL
seurat.obj@meta.data <- seurat.obj@meta.data[, c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "batch", "tissue")]

rm(list.seurat.obj)

# check -------------------------------------------------------------------

seurat.obj@meta.data %>% head()
# orig.ident nCount_RNA nFeature_RNA percent.mt batch tissue
# HFC1_AAACCCAGTAACCCGC       hfc1       8500         2327 19.5176471   one    hfc
# HFC1_AAACGAACAGGTTCAT       hfc1      22492         3938  7.0558421   one    hfc
# HFC1_AAACGAATCCATTGTT       hfc1       5390         1604  0.3523739   one    hfc
# HFC1_AAACGCTAGTAATTGG       hfc1      21957         4273  7.1640024   one    hfc
# HFC1_AAAGGGCCAGCTCGGT       hfc1       4466         1556  5.7072516   one    hfc
# HFC1_AAAGGTAAGAATCTAG       hfc1      47634         6115  7.8699330   one    hfc

seurat.obj@meta.data$orig.ident %>% table()
# hfc1 hfc2 hfc3 lfc1 lfc2 lfc3 
# 309  294  248 2338   81  505 


# scaling and normalizing using SCTransform ---------------------------------------

# ref: https://satijalab.org/seurat/articles/sctransform_vignette.html
# ref: https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/harmony.html

set.seed(121212)

START.TIME <- Sys.time() 

seurat.obj <- seurat.obj %>% 
  SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE) %>% 
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunHarmony(group.by.vars = "batch", assay.use = "SCT", plot_convergence = T, max.iter.harmony = 1000) %>% 
  RunUMAP(assay.use = "SCT", reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>% 
  FindClusters(resolution = 0.2)


FINISH.TIME <- Sys.time() 

FINISH.TIME - START.TIME
# Time difference of 1.439908 mins


# check -------------------------------------------------------------------

seurat.obj %>% print()
# An object of class Seurat 
# 42324 features across 3775 samples within 2 assays 
# Active assay: SCT (18151 features, 3000 variable features)
# 1 other assay present: RNA
# 3 dimensional reductions calculated: pca, harmony, umap


# dimplot ------------------------------------------------------------------

g <- seurat.obj %>% 
  DimPlot(
    group.by = "SCT_snn_res.0.2", 
    combine = FALSE,
    pt.size = .2
  )
g <- lapply(X = g, FUN = function(x) x + theme(legend.position = "top")
            + guides(color = guide_legend(ncol = 6, byrow = TRUE, override.aes = list(size = 4))))
g <- CombinePlots(g, ncol = 1)

setwd(dir.2)
out.f <- "001_myeloid_subset_dimplot_res0.2.pdf"
ggsave(out.f, g, w = 4, h = 4)


# saveRDS -----------------------------------------------------------------

START.TIME <- Sys.time() 

setwd(dir.2)
out.f <- "002_gbm_neuroimmune_myeloid_subset_res0.2.rds"
saveRDS(seurat.obj, out.f)

FINISH.TIME <- Sys.time() 

FINISH.TIME - START.TIME
# Time difference of 57.20986 secs


# find all markers ----------------------------------------------------------

START.TIME <- Sys.time() 

markers <- FindAllMarkers(seurat.obj, 
                          only.pos = TRUE,
                          min.pct = 0.1, 
                          logfc.threshold = 0.1
)

FINISH.TIME <- Sys.time() 

FINISH.TIME - START.TIME
# Time difference of 27.97654 secs


# check -------------------------------------------------------------------

markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  print(n = Inf)


# save --------------------------------------------------------------------

setwd(dir.2)
out.f <- "003_gbm_neuroimmune_myeloid_subset_res0.2_res_find_all_markers.csv"
write.csv(markers, out.f)
  

# heatmap ---------------------------------------------------------------

markers.to.show <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  pull(gene)
  
g <- seurat.obj %>% 
  DoHeatmap(
    features = markers.to.show, 
    slot = "scale.data",
    assay = "RNA")
# g
  
setwd(dir.2)
out.f <- "004_myeloid_subset_heatmap_findmarkers.pdf"
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
# [1] clustree_0.4.3     ggraph_2.0.5       harmony_0.1.0      Rcpp_1.0.8.3       glmGamPoi_1.4.0   
# [6] SeuratObject_4.1.3 Seurat_4.3.0.1     tidylog_1.0.2      forcats_0.5.1      stringr_1.4.0     
# [11] dplyr_1.0.9        purrr_0.3.4        readr_2.1.2        tidyr_1.2.0        tibble_3.1.7      
# [16] ggplot2_3.3.6      tidyverse_1.3.1   
# 
# loaded via a namespace (and not attached):
# [1] readxl_1.3.1                backports_1.2.1             systemfonts_1.0.4          
# [4] plyr_1.8.7                  igraph_1.3.1                lazyeval_0.2.2             
# [7] sp_2.0-0                    splines_4.1.3               listenv_0.8.0              
# [10] scattermore_0.7             GenomeInfoDb_1.28.1         digest_0.6.29              
# [13] htmltools_0.5.5             viridis_0.6.2               fansi_1.0.3                
# [16] magrittr_2.0.3              tensor_1.5                  cluster_2.1.2              
# [19] ROCR_1.0-11                 limma_3.48.3                tzdb_0.1.2                 
# [22] graphlayouts_0.7.1          globals_0.15.0              modelr_0.1.8               
# [25] matrixStats_0.62.0          timechange_0.2.0            spatstat.sparse_3.0-2      
# [28] colorspace_2.0-3            rvest_1.0.1                 ggrepel_0.9.1              
# [31] textshaping_0.3.6           haven_2.4.3                 crayon_1.5.1               
# [34] RCurl_1.98-1.3              jsonlite_1.8.0              progressr_0.10.0           
# [37] spatstat.data_3.0-1         survival_3.2-12             zoo_1.8-10                 
# [40] glue_1.6.2                  polyclip_1.10-0             gtable_0.3.0               
# [43] zlibbioc_1.38.0             XVector_0.32.0              leiden_0.3.9               
# [46] DelayedArray_0.18.0         future.apply_1.8.1          BiocGenerics_0.38.0        
# [49] abind_1.4-5                 scales_1.2.0                DBI_1.1.2                  
# [52] spatstat.random_3.1-5       miniUI_0.1.1.1              viridisLite_0.4.2          
# [55] xtable_1.8-4                reticulate_1.25             clisymbols_1.2.0           
# [58] stats4_4.1.3                htmlwidgets_1.6.2           httr_1.4.3                 
# [61] RColorBrewer_1.1-3          ellipsis_0.3.2              ica_1.0-2                  
# [64] farver_2.1.0                pkgconfig_2.0.3             uwot_0.1.16                
# [67] dbplyr_2.1.1                deldir_1.0-6                utf8_1.2.2                 
# [70] labeling_0.4.2              tidyselect_1.1.2            rlang_1.1.0                
# [73] reshape2_1.4.4              later_1.3.0                 munsell_0.5.0              
# [76] cellranger_1.1.0            tools_4.1.3                 cli_3.6.1                  
# [79] generics_0.1.2              broom_0.7.9                 ggridges_0.5.3             
# [82] fastmap_1.1.0               goftest_1.2-2               fs_1.5.2                   
# [85] tidygraph_1.2.0             fitdistrplus_1.1-5          RANN_2.6.1                 
# [88] pbapply_1.5-0               future_1.25.0               nlme_3.1-152               
# [91] mime_0.12                   xml2_1.3.2                  compiler_4.1.3             
# [94] rstudioapi_0.13             plotly_4.10.0               png_0.1-8                  
# [97] spatstat.utils_3.1-1        reprex_2.0.1                tweenr_1.0.2               
# [100] stringi_1.7.6               lattice_0.20-44             Matrix_1.6-4               
# [103] vctrs_0.6.1                 pillar_1.7.0                lifecycle_1.0.3            
# [106] spatstat.geom_3.2-4         lmtest_0.9-40               RcppAnnoy_0.0.22           
# [109] data.table_1.14.2           cowplot_1.1.1               bitops_1.0-9               
# [112] irlba_2.3.5.1               httpuv_1.6.5                patchwork_1.1.1            
# [115] GenomicRanges_1.44.0        R6_2.5.1                    promises_1.2.0.1           
# [118] KernSmooth_2.23-20          gridExtra_2.3               IRanges_2.26.0             
# [121] parallelly_1.31.1           codetools_0.2-18            MASS_7.3-54                
# [124] assertthat_0.2.1            SummarizedExperiment_1.22.0 withr_2.5.0                
# [127] sctransform_0.3.5           S4Vectors_0.30.2            GenomeInfoDbData_1.2.6     
# [130] parallel_4.1.3              hms_1.1.0                   grid_4.1.3                 
# [133] MatrixGenerics_1.4.3        Rtsne_0.16                  ggforce_0.3.3              
# [136] spatstat.explore_3.2-1      Biobase_2.52.0              shiny_1.7.1                
# [139] lubridate_1.9.2             


# end ---------------------------------------------------------------------


