# Takahide Nejo, MD, PhD
# Okada Lab, Department of Neurological Surgery, UC San Francisco


# note ---------------------------------------------------------------------

# Nejo T, et al., Glioma-neuronal circuit remodeling induces regional immunosuppression

# Please note that we relied on the original annotation directly provided by the authors of the paper (Krishna S, 2023, Nature). 
# Using different versions of packages, such as Seurat, may lead to varying outcomes, which we have not thoroughly validated.


# rm all ------------------------------------------------------------------

rm(list = ls(all.names = T))


# packages ----------------------------------------------------------------

library(tidyverse)
library(tidylog)
library(Seurat)


# check -------------------------------------------------------------------

packageVersion("Seurat")
# [1] ‘4.0.3’


# dir ---------------------------------------------------------------------

dir.1 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/1_00_scRNAseq_data_setup/" # dir.1 <- "/okadalab/data1/tnejo/neuro_immune_proj/dl/20201222_from_abrar/"
dir.2 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/1_01_scRNAseq_cell_type_annotation/"
if(dir.exists(dir.2) == F){dir.create(dir.2)}


# load files --------------------------------------------------------------

START.TIME <- Sys.time() 

setwd(dir.1) 
in.f <- "gbm_neuroimmune.rds"
seurat.obj <- readRDS(in.f)

FINISH.TIME <- Sys.time() 

FINISH.TIME - START.TIME 
# Time difference of 20.84246 secs


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
# [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "percent.mt"      "nCount_SCT"      "nFeature_SCT"    "batch"          
# [8] "tissue"          "SCT_snn_res.0.2" "seurat_clusters" "clust"  

seurat.obj %>% head() %>% print()
{
#                       orig.ident nCount_RNA nFeature_RNA percent.mt nCount_SCT
# HFC1_AAACCCAGTAACCCGC       hfc1       8500         2327 19.5176471      11587
# HFC1_AAACGAACAGGTTCAT       hfc1      22492         3938  7.0558421      13577
# HFC1_AAACGAAGTAGCTTAC       hfc1      19997         5547  2.1053158      13635
# HFC1_AAACGAATCCATTGTT       hfc1       5392         1606  0.3523739      11307
# HFC1_AAACGCTAGTAATTGG       hfc1      21957         4273  7.1640024      13549
# HFC1_AAACGCTCACGCCAGT       hfc1        880          639  0.1136364       9175
# HFC1_AAACGCTGTAGCGCCT       hfc1       5890         2605  0.9337861      11352
# HFC1_AAAGAACCATGTGGTT       hfc1      11132         3805  2.0661157      11714
# HFC1_AAAGGATCAGAAATTG       hfc1      17480         5450  2.4942792      13730
# HFC1_AAAGGGCAGGTTCATC       hfc1       8779         3326  4.3171204      11339
#                       nFeature_SCT batch tissue SCT_snn_res.0.2 seurat_clusters
# HFC1_AAACCCAGTAACCCGC         2327   one    hfc               0               0
# HFC1_AAACGAACAGGTTCAT         3401   one    hfc               0               0
# HFC1_AAACGAAGTAGCTTAC         5422   one    hfc               4               4
# HFC1_AAACGAATCCATTGTT         1828   one    hfc               0               0
# HFC1_AAACGCTAGTAATTGG         3826   one    hfc               0               0
# HFC1_AAACGCTCACGCCAGT         3001   one    hfc               2               2
# HFC1_AAACGCTGTAGCGCCT         2731   one    hfc               2               2
# HFC1_AAAGAACCATGTGGTT         3805   one    hfc               7               7
# HFC1_AAAGGATCAGAAATTG         5446   one    hfc               4               4
# HFC1_AAAGGGCAGGTTCATC         3326   one    hfc               1               1
#                               clust
# HFC1_AAACCCAGTAACCCGC Myeloid Cells
# HFC1_AAACGAACAGGTTCAT Myeloid Cells
# HFC1_AAACGAAGTAGCTTAC   Tumor Cells
# HFC1_AAACGAATCCATTGTT Myeloid Cells
# HFC1_AAACGCTAGTAATTGG Myeloid Cells
# HFC1_AAACGCTCACGCCAGT   Tumor Cells
# HFC1_AAACGCTGTAGCGCCT   Tumor Cells
# HFC1_AAAGAACCATGTGGTT       Neurons
# HFC1_AAAGGATCAGAAATTG   Tumor Cells
# HFC1_AAAGGGCAGGTTCATC   Tumor Cells
}

seurat.obj@meta.data %>% 
  dplyr::select(clust, tissue) %>% 
  table()
{
#                         tissue
# clust                     hfc  lfc
#   Astrocytes               20   41
#   Endothelial Cells       117  163
#   Lymphoid Cells           86  309
#   Myeloid Cells           851 2924
#   Neurons                 151  327
#   Pericytes               116   89
#   Tumor Cells            5062 2905
#   Tumor Cells (Dividing)  263  307
  }

seurat.obj@meta.data %>% 
  dplyr::select(seurat_clusters, tissue) %>% 
  table()
{
#                tissue
# seurat_clusters  hfc  lfc
#              0   851 2924
#              1  1457 1141
#              2  1427  737
#              3   974  353
#              4   388  363
#              5   467  149
#              6   263  307
#              7   151  327
#              8   274  150
#              9    86  309
#              10  117  163
#              11  116   89
#              12   75   12
#              13   20   41
}

seurat.obj@reductions %>% print()
{
# $pca
# A dimensional reduction object with key PC_ 
#  Number of dimensions: 100 
#  Projected dimensional reduction calculated:  FALSE 
#  Jackstraw run: FALSE 
#  Computed using assay: SCT 
# 
# $harmony
# A dimensional reduction object with key harmony_ 
#  Number of dimensions: 100 
#  Projected dimensional reduction calculated:  TRUE 
#  Jackstraw run: FALSE 
#  Computed using assay: SCT 
# 
# $umap
# A dimensional reduction object with key UMAP_ 
#  Number of dimensions: 2 
#  Projected dimensional reduction calculated:  FALSE 
#  Jackstraw run: FALSE 
#  Computed using assay: SCT 
}


seurat.obj@commands %>% print()
{
# $SCTransform.RNA
# Command: SCTransform(funccon, vars.to.regress = c("nCount_RNA", "percent.mt"))
# Time: 2019-12-10 12:18:22
# assay : RNA 
# new.assay.name : SCT 
# do.correct.umi : TRUE 
# variable.features.n : 3000 
# variable.features.rv.th : 1.3 
# vars.to.regress : nCount_RNA percent.mt 
# do.scale : FALSE 
# do.center : TRUE 
# clip.range : -21.39392 21.39392 
# conserve.memory : FALSE 
# return.only.var.genes : TRUE 
# seed.use : 1448145 
# verbose : TRUE 
# 
# $RunPCA.SCT
# Command: RunPCA(funccon, npcs = 100)
# Time: 2019-12-10 12:24:18
# assay : SCT 
# npcs : 100 
# rev.pca : FALSE 
# weight.by.var : TRUE 
# verbose : TRUE 
# ndims.print : 1 2 3 4 5 
# nfeatures.print : 30 
# reduction.name : pca 
# reduction.key : PC_ 
# seed.use : 42 
# 
# $`Seurat::ProjectDim.SCT.harmony`
# Command: Seurat::ProjectDim(object, reduction = "harmony", overwrite = TRUE,     verbose = FALSE)
# Time: 2019-12-10 12:28:45
# reduction : harmony 
# assay : SCT 
# dims.print : 1 2 3 4 5 
# nfeatures.print : 20 
# overwrite : TRUE 
# do.center : FALSE 
# verbose : FALSE 
# 
# $RunUMAP.SCT.harmony
# Command: RunUMAP(fc, dims = 1:30, reduction = "harmony", min.dist = 0.4)
# Time: 2020-01-11 17:31:45
# dims : 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 
# reduction : harmony 
# assay : SCT 
# umap.method : uwot 
# n.neighbors : 30 
# n.components : 2 
# metric : cosine 
# learning.rate : 1 
# min.dist : 0.4 
# spread : 1 
# set.op.mix.ratio : 1 
# local.connectivity : 1 
# repulsion.strength : 1 
# negative.sample.rate : 5 
# uwot.sgd : FALSE 
# seed.use : 42 
# angular.rp.forest : FALSE 
# verbose : TRUE 
# reduction.name : umap 
# reduction.key : UMAP_ 
# 
# $FindNeighbors.SCT.harmony
# Command: FindNeighbors(fc, dims = 1:30, reduction = "harmony")
# Time: 2020-01-11 17:32:15
# reduction : harmony 
# dims : 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 
# assay : SCT 
# k.param : 20 
# compute.SNN : TRUE 
# prune.SNN : 0.06666667 
# nn.method : rann 
# annoy.metric : euclidean 
# nn.eps : 0 
# verbose : TRUE 
# force.recalc : FALSE 
# do.plot : FALSE 
# graph.name : SCT_nn SCT_snn 
# 
# $FindClusters
# Command: FindClusters(fc, resolution = 0.2)
# Time: 2020-01-11 17:32:27
# graph.name : SCT_snn 
# modularity.fxn : 1 
# resolution : 0.2 
# algorithm : 1 
# n.start : 10 
# n.iter : 10 
# random.seed : 0 
# group.singletons : TRUE 
# verbose : TRUE 
}


# edit the metadata --------------------------------------------------------------------

seurat.obj.orig <- seurat.obj

metadata.new <- seurat.obj@meta.data %>% 
  as_tibble() %>% 
  mutate(clust = ifelse(grepl("Tumor", clust), "Tumor Cells", clust)) %>% 
  as.data.frame()

rownames(metadata.new) <- rownames(seurat.obj@meta.data)

seurat.obj@meta.data <- metadata.new

seurat.obj <- SetIdent(seurat.obj, value = "clust")


# check -------------------------------------------------------------------

seurat.obj@meta.data %>% head()
# orig.ident nCount_RNA nFeature_RNA percent.mt nCount_SCT nFeature_SCT batch tissue SCT_snn_res.0.2
# HFC1_AAACCCAGTAACCCGC       hfc1       8500         2327 19.5176471      11587         2327   one    hfc               0
# HFC1_AAACGAACAGGTTCAT       hfc1      22492         3938  7.0558421      13577         3401   one    hfc               0
# HFC1_AAACGAAGTAGCTTAC       hfc1      19997         5547  2.1053158      13635         5422   one    hfc               4
# HFC1_AAACGAATCCATTGTT       hfc1       5392         1606  0.3523739      11307         1828   one    hfc               0
# HFC1_AAACGCTAGTAATTGG       hfc1      21957         4273  7.1640024      13549         3826   one    hfc               0
# HFC1_AAACGCTCACGCCAGT       hfc1        880          639  0.1136364       9175         3001   one    hfc               2
# seurat_clusters         clust
# HFC1_AAACCCAGTAACCCGC               0 Myeloid Cells
# HFC1_AAACGAACAGGTTCAT               0 Myeloid Cells
# HFC1_AAACGAAGTAGCTTAC               4   Tumor Cells
# HFC1_AAACGAATCCATTGTT               0 Myeloid Cells
# HFC1_AAACGCTAGTAATTGG               0 Myeloid Cells
# HFC1_AAACGCTCACGCCAGT               2   Tumor Cells

seurat.obj@meta.data %>% 
  dplyr::select(clust, tissue) %>% 
  table()
# tissue
# clust                hfc  lfc
# Astrocytes          20   41
# Endothelial Cells  117  163
# Lymphoid Cells      86  309
# Myeloid Cells      851 2924
# Neurons            151  327
# Pericytes          116   89
# Tumor Cells       5325 3212


# DimPlot ------------------------------------------------------------------

TITLE = "13,731 cells"

g <- DimPlot(seurat.obj, label = TRUE)
g <- g + NoLegend()
g <- g + theme(
  plot.title = element_text(hjust = 0.5, face = "bold"), 
  legend.position = "none", 
  panel.grid.major = element_line(color = "#f0f0f0"), 
  axis.text = element_blank()
) 
g <- g + labs(title = TITLE, x = "UMAP-1", y = "UMAP-2")
plot(g)

setwd(dir.2)
out.f <- "001_overview_umap_whole.pdf"
ggsave(out.f, g, w = 8, h = 6)


# Extended Data Fig. S1a --------------------------------------------------

g <- DimPlot(seurat.obj, label = TRUE, split.by = "tissue")
g <- g & NoLegend()
g <- g & theme(
  plot.title = element_text(hjust = 0.5, face = "bold"), 
  legend.position = "none", 
  panel.grid.major = element_line(color = "#f0f0f0"), 
  axis.text = element_blank()
) 
g <- g + labs(x = "UMAP-1", y = "UMAP-2")
plot(g)

setwd(dir.2)
out.f <- "002_overview_umap_hfc_lfc.pdf"
ggsave(out.f, g, w = 12, h = 6)


# Table for Extended Data Fig. S1a --------------------------------------------------

# for each patient, get the numbers of the cells in each cluster.

for(i in 1:3){
  print(i)
  
  table.i <- seurat.obj@meta.data %>% 
    dplyr::filter(batch == c("one", "two", "three")[i]) %>% 
    dplyr::select(clust, tissue) %>% 
    table() %>% 
    as_tibble() %>% 
    spread(key = "tissue", value = "n") %>% 
    mutate(total = hfc + lfc)
  
  colnames(table.i)[-1] <- paste0(colnames(table.i)[-1], ".pt", i)
  
  if(i == 1){
    table.combined <- table.i
  }else{
    table.combined <- table.combined %>% 
      left_join(table.i, by = "clust")
  }
}

table.4 <- seurat.obj@meta.data %>% 
  dplyr::select(clust, tissue) %>% 
  table() %>% 
  as_tibble() %>% 
  spread(key = "tissue", value = "n") %>% 
  mutate(total = hfc + lfc)

colnames(table.4)[-1] <- paste0(colnames(table.4)[-1], ".whole")

table.combined <- table.combined %>% 
  left_join(table.4, by = "clust")

table.combined$clust <- factor(table.combined$clust, 
                               levels = c("Tumor Cells", "Myeloid Cells", "Neurons", "Lymphoid Cells", "Endothelial Cells", "Pericytes", "Astrocytes"))

table.combined <- table.combined %>% 
  arrange(clust) %>% 
  rename(annotation = clust)


# check -------------------------------------------------------------------

table.combined %>% dim()
# [1]  7 13

table.combined %>% colnames()
# [1] "annotation"  "hfc.pt1"     "lfc.pt1"     "total.pt1"   "hfc.pt2"     "lfc.pt2"    
# [7] "total.pt2"   "hfc.pt3"     "lfc.pt3"     "total.pt3"   "hfc.whole"   "lfc.whole"  
# [13] "total.whole"

table.combined %>% print()
{
  # # A tibble: 7 × 13
  # annotation        hfc.pt1 lfc.pt1 total.pt1 hfc.pt2 lfc.pt2 total.pt2 hfc.pt3 lfc.pt3
  # <fct>               <int>   <int>     <int>   <int>   <int>     <int>   <int>   <int>
  # 1 Tumor Cells           525    1250      1775    2143     406      2549    2657    1556
  # 2 Myeloid Cells         309    2338      2647     294      81       375     248     505
  # 3 Neurons                37     276       313      54      49       103      60       2
  # 4 Lymphoid Cells         64     284       348      10       7        17      12      18
  # 5 Endothelial Cells      12      67        79       5       1         6     100      95
  # 6 Pericytes              19      33        52      15       1        16      82      55
  # 7 Astrocytes              1      29        30       8       4        12      11       8
  # # … with 4 more variables: total.pt3 <int>, hfc.whole <int>, lfc.whole <int>,
  # #   total.whole <int>
}


# save --------------------------------------------------------------------

setwd(dir.2)
out.f <- "003_cell_numbers_per_pt_n_cluster.tsv"
write_tsv(table.combined, out.f)


# save --------------------------------------------------------------------

START.TIME <- Sys.time() 

setwd(dir.2)
out.f <- "004_gbm_neuroimmune_v2.rds"
saveRDS(seurat.obj, out.f)


FINISH.TIME <- Sys.time() 

FINISH.TIME - START.TIME 
# Time difference of 3.606859 mins


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
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] RColorBrewer_1.1-3 harmony_0.1.0      Rcpp_1.0.8.3       tidylog_1.0.2      forcats_0.5.1      stringr_1.4.0     
# [7] dplyr_1.0.9        purrr_0.3.4        readr_2.1.2        tidyr_1.2.0        tibble_3.1.7       ggplot2_3.3.6     
# [13] tidyverse_1.3.1    SeuratObject_4.1.3 Seurat_4.3.0.1     SPATA2_2.0.4      
# 
# loaded via a namespace (and not attached):
# [1] readxl_1.3.1                backports_1.2.1             plyr_1.8.7                  igraph_1.3.1               
# [5] lazyeval_0.2.2              sp_2.0-0                    splines_4.1.3               listenv_0.8.0              
# [9] scattermore_0.7             GenomeInfoDb_1.28.1         digest_0.6.29               htmltools_0.5.5            
# [13] tiff_0.1-11                 fansi_1.0.3                 magrittr_2.0.3              SPATAData_0.0.0.9000       
# [17] tensor_1.5                  cluster_2.1.2               ROCR_1.0-11                 tzdb_0.1.2                 
# [21] globals_0.15.0              modelr_0.1.8                matrixStats_0.62.0          timechange_0.2.0           
# [25] spatstat.sparse_3.0-2       jpeg_0.1-9                  colorspace_2.0-3            rvest_1.0.1                
# [29] ggrepel_0.9.1               xfun_0.25                   haven_2.4.3                 crayon_1.5.1               
# [33] RCurl_1.98-1.3              jsonlite_1.8.0              progressr_0.10.0            spatstat.data_3.0-1        
# [37] survival_3.2-12             zoo_1.8-10                  glue_1.6.2                  polyclip_1.10-0            
# [41] gtable_0.3.0                zlibbioc_1.38.0             XVector_0.32.0              leiden_0.3.9               
# [45] DelayedArray_0.18.0         future.apply_1.8.1          SingleCellExperiment_1.14.1 BiocGenerics_0.38.0        
# [49] abind_1.4-5                 scales_1.2.0                DBI_1.1.2                   spatstat.random_3.1-5      
# [53] miniUI_0.1.1.1              viridisLite_0.4.0           xtable_1.8-4                units_0.8-0                
# [57] reticulate_1.25             clisymbols_1.2.0            stats4_4.1.3                htmlwidgets_1.6.2          
# [61] httr_1.4.3                  anndata_0.7.5.6             ellipsis_0.3.2              ica_1.0-2                  
# [65] farver_2.1.0                pkgconfig_2.0.3             uwot_0.1.16                 dbplyr_2.1.1               
# [69] deldir_1.0-6                locfit_1.5-9.4              utf8_1.2.2                  labeling_0.4.2             
# [73] tidyselect_1.1.2            rlang_1.1.0                 reshape2_1.4.4              later_1.3.0                
# [77] cellranger_1.1.0            munsell_0.5.0               tools_4.1.3                 cli_3.6.1                  
# [81] generics_0.1.2              broom_0.7.9                 ggridges_0.5.3              evaluate_0.14              
# [85] fastmap_1.1.0               fftwtools_0.9-11            yaml_2.3.5                  goftest_1.2-2              
# [89] knitr_1.33                  fs_1.5.2                    fitdistrplus_1.1-5          RANN_2.6.1                 
# [93] pbapply_1.5-0               future_1.25.0               nlme_3.1-152                mime_0.12                  
# [97] xml2_1.3.2                  compiler_4.1.3              rstudioapi_0.13             plotly_4.10.0              
# [101] png_0.1-8                   spatstat.utils_3.0-3        reprex_2.0.1                confuns_1.0.2              
# [105] stringi_1.7.6               lattice_0.20-44             Matrix_1.6-0                vctrs_0.6.1                
# [109] pillar_1.7.0                lifecycle_1.0.3             spatstat.geom_3.2-4         lmtest_0.9-40              
# [113] RcppAnnoy_0.0.19            data.table_1.14.2           cowplot_1.1.1               bitops_1.0-7               
# [117] irlba_2.3.5                 httpuv_1.6.5                patchwork_1.1.1             GenomicRanges_1.44.0       
# [121] R6_2.5.1                    promises_1.2.0.1            KernSmooth_2.23-20          gridExtra_2.3              
# [125] IRanges_2.26.0              parallelly_1.31.1           codetools_0.2-18            MASS_7.3-54                
# [129] assertthat_0.2.1            SummarizedExperiment_1.22.0 withr_2.5.0                 sctransform_0.3.5          
# [133] S4Vectors_0.30.2            GenomeInfoDbData_1.2.6      parallel_4.1.3              hms_1.1.0                  
# [137] EBImage_4.34.0              grid_4.1.3                  rmarkdown_2.10              MatrixGenerics_1.4.3       
# [141] Rtsne_0.16                  spatstat.explore_3.2-1      lubridate_1.9.2             Biobase_2.52.0             
# [145] shiny_1.7.1                


# end ---------------------------------------------------------------------


