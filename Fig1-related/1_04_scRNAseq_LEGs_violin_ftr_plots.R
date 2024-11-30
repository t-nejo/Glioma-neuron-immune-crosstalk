# Takahide Nejo, MD, PhD
# Okada Lab, Department of Neurological Surgery, UC San Francisco


# note ---------------------------------------------------------------------

# Nejo T, et al., Glioma-neuronal circuit remodeling induces regional immunosuppression

# Using different versions of packages, such as Seurat, may lead to varying outcomes, which we have not thoroughly validated.

# The RDS object "001_gbm_neuroimmune_v3_subset_list.rds," "002_res_findmarkers_mast_list.rds," and "000_res_fgsea_hallmark_list.rds" are used as inputs. For details, please refer to the previous step "1_02_scRNAseq_dge.R" and "1_03_scRNAseq_fgsea.R". 


# rm all ------------------------------------------------------------------

rm(list = ls(all.names = T))


# packages ----------------------------------------------------------------

library(tidyverse)
library(tidylog)
library(Seurat)
library(ggpubr)
library(ggsci)
library(wesanderson)


# check -------------------------------------------------------------------

packageVersion("Seurat")
# [1] ‘4.3.0.1’


# dir ---------------------------------------------------------------------

dir.1 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/1_02_scRNAseq_dge/"
dir.1b <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/1_03_scRNAseq_fgsea/"
dir.2 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/1_04_scRNAseq_LEGs_violin_ftr_plots/"
if(dir.exists(dir.2) == F){dir.create(dir.2)}


# load files --------------------------------------------------------------

START.TIME <- Sys.time() 

setwd(dir.1) 
in.f <- "001_gbm_neuroimmune_v3_subset_list.rds"
list.seurat.obj <- readRDS(in.f)

setwd(dir.1) 
in.f <- "002_res_findmarkers_mast_list.rds"
list.res.mast <- readRDS(in.f)

setwd(dir.1b)
in.f <- "000_res_fgsea_hallmark_list.rds"
list.res.fgsea <- readRDS(in.f)

FINISH.TIME <- Sys.time() 

FINISH.TIME - START.TIME 
# Time difference of 38.6142 secs


# check -------------------------------------------------------------------

list.seurat.obj %>% length()
# [1] 4

list.seurat.obj %>% names()
# [1] "Tumor Cells"    "Myeloid Cells"  "Lymphoid Cells" "Astrocytes"    

list.seurat.obj[[1]]
# An object of class Seurat 
# 47200 features across 8537 samples within 2 assays 
# Active assay: RNA (24173 features, 0 variable features)
# 1 other assay present: SCT
# 3 dimensional reductions calculated: pca, harmony, umap


list.res.mast %>% length()
# [1] 4

list.res.mast %>% names()
# [1] "Tumor Cells"    "Myeloid Cells"  "Lymphoid Cells" "Astrocytes"    

list.res.mast[[1]] %>% head()
# p_val avg_log2FC pct.1 pct.2     p_val_adj
# CCL4   2.033300e-188 -2.7553522 0.229 0.647 4.915096e-184
# CCL3   1.519361e-172 -2.3478363 0.219 0.640 3.672751e-168
# CCL4L2 5.622437e-142 -1.4319096 0.145 0.462 1.359112e-137
# S100A1 4.665051e-136 -0.9444470 0.321 0.415 1.127683e-131
# HIF1A  1.126022e-131  0.8002181 0.788 0.699 2.721932e-127
# PLA2R1 1.456949e-125 -0.3571337 0.110 0.286 3.521883e-121


list.res.fgsea %>% length()
# [1] 4

list.res.fgsea %>% names()
# [1] "Tumor Cells"    "Myeloid Cells"  "Lymphoid Cells" "Astrocytes"    

list.res.fgsea[[1]] %>% head()
# # A tibble: 6 × 8
# pathway                            pval     padj log2err    ES   NES  size leadingEdge
# <chr>                             <dbl>    <dbl>   <dbl> <dbl> <dbl> <int> <list>     
#   1 HALLMARK_OXIDATIVE_PHOSPHORYL… 7.30e-14 3.65e-12   0.955 0.750  2.10   200 <chr [134]>
#   2 HALLMARK_CHOLESTEROL_HOMEOSTA… 1.09e- 4 1.17e- 3   0.538 0.712  1.74    72 <chr [37]> 
#   3 HALLMARK_MTORC1_SIGNALING      1.40e- 4 1.17e- 3   0.519 0.586  1.64   197 <chr [78]> 
#   4 HALLMARK_GLYCOLYSIS            1.78e- 4 1.27e- 3   0.519 0.583  1.63   194 <chr [77]> 
#   5 HALLMARK_ADIPOGENESIS          2.35e- 4 1.47e- 3   0.519 0.584  1.63   195 <chr [84]> 
#   6 HALLMARK_PROTEIN_SECRETION     3.78e- 4 1.89e- 3   0.498 0.664  1.70    95 <chr [56]> 


# check the leading edge genes --------------------------------------------

imm.pathways <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE", 
  "HALLMARK_INTERFERON_GAMMA_RESPONSE", 
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)

for(i in 1:4){
  print(i)
  
  names(list.res.fgsea)[i] %>% print()
  
  res.fgsea.i <- list.res.fgsea[[i]] %>% 
    dplyr::filter(pathway %in% imm.pathways) 
  
  res.fgsea.i %>% 
    pull(leadingEdge) %>% 
    unlist() %>% 
    table() %>% 
    as_tibble() %>% 
    arrange(desc(n)) %>% 
    dplyr::filter(n >= 2) %>% 
    print(n = Inf)
}
{
  # [1] 1
  # [1] "Tumor Cells"
  # # A tibble: 33 × 2
  # .           n
  # <chr>   <int>
  # 1 CCL2        3
  # 2 CD69        3
  # 3 IRF1        3
  # 4 NFKB1       3
  # 5 TNFAIP6     3
  # 6 BST2        2
  # 7 BTG2        2
  # 8 CCL5        2
  # 9 CD40        2
  # 10 CSF1        2
  # 11 EIF2AK2     2
  # 12 FPR1        2
  # 13 GPR183      2
  # 14 ICAM1       2
  # 15 IFIH1       2
  # 16 IFIT2       2
  # 17 IL10RA      2
  # 18 IL18        2
  # 19 IL1B        2
  # 20 INHBA       2
  # 21 KLF6        2
  # 22 LCP2        2
  # 23 LY6E        2
  # 24 NMI         2
  # 25 OLR1        2
  # 26 PTGER4      2
  # 27 PTGS2       2
  # 28 SOD2        2
  # 29 SRI         2
  # 30 TAPBP       2
  # 31 TLR2        2
  # 32 TNFAIP3     2
  # 33 TNFSF10     2
  # [1] 2
  # [1] "Myeloid Cells"
  # # A tibble: 38 × 2
  # .            n
  # <chr>    <int>
  # 1 CCL2         3
  # 2 CCL5         3
  # 3 CD69         3
  # 4 CDKN1A       3
  # 5 CXCL10       3
  # 6 IL6          3
  # 7 NFKB1        3
  # 8 NFKBIA       3
  # 9 BST2         2
  # 10 BTG2         2
  # 11 CCRL2        2
  # 12 CSF1         2
  # 13 CXCL9        2
  # 14 DDX58        2
  # 15 EDN1         2
  # 16 EIF2AK2      2
  # 17 GPR183       2
  # 18 HBEGF        2
  # 19 IFIH1        2
  # 20 IFIT2        2
  # 21 IL18         2
  # 22 IL1A         2
  # 23 IL1B         2
  # 24 INHBA        2
  # 25 IRF7         2
  # 26 KLF6         2
  # 27 LCP2         2
  # 28 LY6E         2
  # 29 MYC          2
  # 30 NMI          2
  # 31 OLR1         2
  # 32 PTGER4       2
  # 33 PTGS2        2
  # 34 PTPRE        2
  # 35 SERPINE1     2
  # 36 TNFAIP3      2
  # 37 TNFSF10      2
  # 38 TNFSF9       2
  # [1] 3
  # [1] "Lymphoid Cells"
  # # A tibble: 44 × 2
  # .            n
  # <chr>    <int>
  # 1 CCL2         3
  # 2 CDKN1A       3
  # 3 IRF1         3
  # 4 NAMPT        3
  # 5 NFKB1        3
  # 6 PDE4B        3
  # 7 RIPK2        3
  # 8 ATP2B1       2
  # 9 BST2         2
  # 10 BTG1         2
  # 11 BTG2         2
  # 12 CSF1         2
  # 13 CXCL10       2
  # 14 DDX58        2
  # 15 EIF2AK2      2
  # 16 F3           2
  # 17 FPR1         2
  # 18 GPR183       2
  # 19 HBEGF        2
  # 20 HIF1A        2
  # 21 IFIH1        2
  # 22 IFIT2        2
  # 23 IL15         2
  # 24 IL18         2
  # 25 IL1A         2
  # 26 IL1B         2
  # 27 IL4R         2
  # 28 IL7R         2
  # 29 IRF7         2
  # 30 LY6E         2
  # 31 MXD1         2
  # 32 MYC          2
  # 33 NMI          2
  # 34 OLR1         2
  # 35 PLAUR        2
  # 36 PTGER4       2
  # 37 PTGS2        2
  # 38 SERPINE1     2
  # 39 SOCS3        2
  # 40 SOD2         2
  # 41 SRI          2
  # 42 TAPBP        2
  # 43 TNFSF10      2
  # 44 TNFSF9       2
  # [1] 4
  # [1] "Astrocytes"
  # # A tibble: 8 × 2
  # .           n
  # <chr>   <int>
  # 1 CCL2        2
  # 2 CXCL10      2
  # 3 NAMPT       2
  # 4 NFKB1       2
  # 5 PTGS2       2
  # 6 SOD2        2
  # 7 TAP1        2
  # 8 TNFAIP6     2
}


# Extended Data Fig. S1c, e, g --------------------------------------------

for(i in 1:3){
  print(i)
  seurat.obj.i <- list.seurat.obj[[i]]
  
  Idents(seurat.obj.i) <- "tissue"
  
  if(i == 1){
    genes <- c("CCL4", "CD74", "CD83", "FOSB", "HLA-DRB1", "IER2", "ISG15", "STAT1")
  }else if(i == 2){
    genes <- c("CCL4", "CD83", "IFI44L", "IL1B", "IRF7", "ISG15", "NFKBIA", "TNF")
  }else if(i == 3){
    genes <- c("CCL4", "CD74", "CD83", "IER2", "IFI44L", "ISG15", "NFKB1", "NR4A3")
  }
  
  gg.list <- list()
  for(j in 1:length(genes)){
    gene.j <- genes[j]
    
    padj.j <- list.res.mast[[i]] %>% 
      rownames_to_column("symbol") %>% 
      as_tibble() %>% 
      dplyr::filter(symbol == gene.j) %>% 
      pull(p_val_adj)
    
    padj.j <- signif(padj.j, digits=2)
    
    g <- seurat.obj.i %>% 
      VlnPlot(
        features = gene.j, 
        assay = "SCT", 
        cols = alpha(pal_nejm()(2), 0.75), # ref: https://www.biostars.org/p/9571401/
        pt.size = 0
      )
    g <- g & labs(title = paste0(genes[j], "\n", "p.adj = ", padj.j))
    g <- g & theme(
      panel.grid.major = element_line(color = "#f0f0f0"), 
      axis.title = element_blank(), 
      legend.position = "none"
    ) 
    # plot(g)
    
    gg.list[[j]] <- g
  }
  
  g <- ggarrange(
    plotlist = gg.list, 
    nrow = 1, ncol = 8
  )
  
  name.i <- c("tumor", "myeloid", "lymphoid")[i] 
  setwd(dir.2)
  out.f <- paste0("001_curated_LEGs_violin_", i, "_", name.i, ".pdf")
  ggsave(out.f, g, width = 15, height = 2.5)  
}


# Extended Data Fig. S1d, f, h --------------------------------------------

for(i in 1:3){
  print(i)
  seurat.obj.i <- list.seurat.obj[[i]]
  
  Idents(seurat.obj.i) <- "tissue"
  
  if(i == 1){
    genes <- c("CCL4", "CD83", "HLA-DRB1")
  }else if(i == 2){
    genes <- c("CCL4", "IFI44L", "IL1B")
  }else if(i == 3){
    genes <- c("CCL4", "ISG15", "NR4A3")
  }
  
  gg.list <- list()
  for(j in 1:length(genes)){
    gene.j <- genes[j]
    
    g <- seurat.obj.i %>% 
      FeaturePlot(
        split.by = "tissue",
        features = gene.j, 
        cols = c("lightgrey", wes_palette("Zissou1", 100, type = "continuous")),
        pt.size = 0.6,
        keep.scale = "all"
      )
    g <- g & labs(subtitle = gene.j)
    g <- g & theme(
      legend.position = "right",
      panel.grid.major = element_line(color = "#f0f0f0"),
      axis.text = element_blank(),
      axis.title = element_blank()
    )
    if(i == 1){
      g <- g & xlim(-8, 8.5) & ylim(-13.5, 3.5)
    }else if(i == 2){
      g <- g & xlim(-4, 4.5) & ylim(4.5, 13.5)
    }else if(i == 3){
      g <- g & xlim(5.5, 12.5) & ylim(5, 10)
    }

    gg.list[[j]] <- g
  }
  
  g <- ggarrange(
    plotlist = gg.list, 
    nrow = 3, ncol = 1
  )
  # plot(g)
  
  name.i <- c("tumor", "myeloid", "lymphoid")[i] 
  setwd(dir.2)
  out.f <- paste0("002_curated_LEGs_ftrplot_", i, "_", name.i, ".pdf")
  ggsave(out.f, g, w = 9.6, h = 12)  
}


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
# [1] wesanderson_0.3.6  fgsea_1.18.0       msigdbr_7.4.1      ggrepel_0.9.1      cowplot_1.1.1      glmGamPoi_1.4.0   
# [7] sctransform_0.3.5  ggpubr_0.4.0       ggsci_2.9          RColorBrewer_1.1-3 harmony_0.1.0      Rcpp_1.0.8.3      
# [13] tidylog_1.0.2      forcats_0.5.1      stringr_1.4.0      dplyr_1.0.9        purrr_0.3.4        readr_2.1.2       
# [19] tidyr_1.2.0        tibble_3.1.7       ggplot2_3.3.6      tidyverse_1.3.1    SeuratObject_4.1.3 Seurat_4.3.0.1    
# [25] SPATA2_2.0.4      
# 
# loaded via a namespace (and not attached):
# [1] utf8_1.2.2                  spatstat.explore_3.2-1      reticulate_1.25             tidyselect_1.1.2           
# [5] htmlwidgets_1.6.2           BiocParallel_1.26.2         grid_4.1.3                  Rtsne_0.16                 
# [9] munsell_0.5.0               codetools_0.2-18            ica_1.0-2                   units_0.8-0                
# [13] future_1.25.0               miniUI_0.1.1.1              withr_2.5.0                 spatstat.random_3.1-5      
# [17] colorspace_2.0-3            progressr_0.10.0            Biobase_2.52.0              rstudioapi_0.13            
# [21] stats4_4.1.3                SingleCellExperiment_1.14.1 ROCR_1.0-11                 ggsignif_0.6.2             
# [25] tensor_1.5                  listenv_0.8.0               labeling_0.4.2              MatrixGenerics_1.4.3       
# [29] GenomeInfoDbData_1.2.6      polyclip_1.10-0             farver_2.1.0                parallelly_1.31.1          
# [33] vctrs_0.6.1                 generics_0.1.2              timechange_0.2.0            R6_2.5.1                   
# [37] GenomeInfoDb_1.28.1         locfit_1.5-9.4              bitops_1.0-7                spatstat.utils_3.0-3       
# [41] DelayedArray_0.18.0         assertthat_0.2.1            promises_1.2.0.1            scales_1.2.0               
# [45] gtable_0.3.0                globals_0.15.0              goftest_1.2-2               rlang_1.1.0                
# [49] clisymbols_1.2.0            systemfonts_1.0.4           splines_4.1.3               rstatix_0.7.0              
# [53] lazyeval_0.2.2              spatstat.geom_3.2-4         broom_0.7.9                 reshape2_1.4.4             
# [57] abind_1.4-5                 modelr_0.1.8                backports_1.2.1             httpuv_1.6.5               
# [61] tools_4.1.3                 ellipsis_0.3.2              BiocGenerics_0.38.0         ggridges_0.5.3             
# [65] plyr_1.8.7                  zlibbioc_1.38.0             RCurl_1.98-1.3              deldir_1.0-6               
# [69] pbapply_1.5-0               S4Vectors_0.30.2            zoo_1.8-10                  SummarizedExperiment_1.22.0
# [73] haven_2.4.3                 cluster_2.1.2               fs_1.5.2                    magrittr_2.0.3             
# [77] data.table_1.14.2           scattermore_0.7             openxlsx_4.2.4              lmtest_0.9-40              
# [81] reprex_2.0.1                RANN_2.6.1                  fitdistrplus_1.1-5          anndata_0.7.5.6            
# [85] matrixStats_0.62.0          hms_1.1.0                   patchwork_1.1.1             mime_0.12                  
# [89] fftwtools_0.9-11            xtable_1.8-4                rio_0.5.27                  jpeg_0.1-9                 
# [93] readxl_1.3.1                IRanges_2.26.0              gridExtra_2.3               compiler_4.1.3             
# [97] KernSmooth_2.23-20          crayon_1.5.1                SPATAData_0.0.0.9000        htmltools_0.5.5            
# [101] later_1.3.0                 tzdb_0.1.2                  tiff_0.1-11                 lubridate_1.9.2            
# [105] DBI_1.1.2                   dbplyr_2.1.1                MASS_7.3-54                 babelgene_21.4             
# [109] Matrix_1.6-0                car_3.0-11                  cli_3.6.1                   parallel_4.1.3             
# [113] igraph_1.3.1                GenomicRanges_1.44.0        pkgconfig_2.0.3             foreign_0.8-81             
# [117] sp_2.0-0                    confuns_1.0.2               plotly_4.10.0               spatstat.sparse_3.0-2      
# [121] xml2_1.3.2                  XVector_0.32.0              rvest_1.0.1                 digest_0.6.29              
# [125] RcppAnnoy_0.0.19            spatstat.data_3.0-1         fastmatch_1.1-3             cellranger_1.1.0           
# [129] leiden_0.3.9                uwot_0.1.16                 curl_4.3.2                  shiny_1.7.1                
# [133] EBImage_4.34.0              lifecycle_1.0.3             nlme_3.1-152                jsonlite_1.8.0             
# [137] carData_3.0-4               viridisLite_0.4.0           fansi_1.0.3                 pillar_1.7.0               
# [141] lattice_0.20-44             fastmap_1.1.0               httr_1.4.3                  survival_3.2-12            
# [145] glue_1.6.2                  zip_2.2.0                   png_0.1-8                   stringi_1.7.6              
# [149] textshaping_0.3.6           irlba_2.3.5                 future.apply_1.8.1         


# end ---------------------------------------------------------------------


