# Takahide Nejo, MD, PhD
# Okada Lab, Department of Neurological Surgery, UC San Francisco


# note ---------------------------------------------------------------------

# Nejo T, et al., Glioma-neuronal circuit remodeling induces regional immunosuppression

# Using different versions of packages, such as Seurat, may lead to varying outcomes, which we have not thoroughly validated.

# The RDS object "001_gbm_neuroimmune_v3_subset_list.rds" is used as inputs. For details, please refer to the previous step "1_02_scRNAseq_dge.R". 

# Please note that due to the random permutations involved, multiple runs of AddModuleScores() may yield slightly different results each time. 


# rm all ------------------------------------------------------------------

rm(list = ls(all.names = T))


# packages ----------------------------------------------------------------

library(tidyverse)
library(tidylog)
library(Seurat)
library(msigdbr)
library(nortest)
library(ggpubr)
library(ggsci)
library(wesanderson)


# check -------------------------------------------------------------------

packageVersion("Seurat")
# [1] ‘4.3.0.1’


# dir ---------------------------------------------------------------------

dir.1 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/1_02_scRNAseq_dge/"
dir.2 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/1_05_scRNAseq_gene_sig_scores/"
if(dir.exists(dir.2) == F){dir.create(dir.2)}


# load files --------------------------------------------------------------

START.TIME <- Sys.time() 

setwd(dir.1) 
in.f <- "001_gbm_neuroimmune_v3_subset_list.rds"
list.seurat.obj <- readRDS(in.f)

FINISH.TIME <- Sys.time() 

FINISH.TIME - START.TIME 
# Time difference of 41.70154 secs


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


# pathways of interest --------------------------------------------

imm.pathways <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE", 
  "HALLMARK_INTERFERON_GAMMA_RESPONSE", 
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)


# prep the genesets of interest ---------------------------------------------------------
    
test_gs <- msigdbr(
  species = "Homo sapiens", 
  category = "H"
) %>% 
  dplyr::filter(gs_name %in% imm.pathways) %>% 
  dplyr::select(gs_name, human_gene_symbol) %>% 
  group_by(gs_name) %>% 
  summarise(human_gene_symbol = list(human_gene_symbol))


# check -------------------------------------------------------------------

test_gs %>% dim()
# [1] 3 2

test_gs %>% head() %>% print()
# # A tibble: 3 × 2
# gs_name                            human_gene_symbol
# <chr>                              <list>           
#   1 HALLMARK_INFLAMMATORY_RESPONSE     <chr [222]>      
#   2 HALLMARK_INTERFERON_GAMMA_RESPONSE <chr [284]>      
#   3 HALLMARK_TNFA_SIGNALING_VIA_NFKB   <chr [227]>    



# AddModuleScores() --------------------------------------------------------

START.TIME <- Sys.time() 

list.seurat.obj.2 <- list()
for(i in 1:3){
  print(i) 
  name.i <- c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells")[i]
  
  seurat.obj.i <- list.seurat.obj[[i]]
  
  seurat.obj.i2 <- seurat.obj.i %>% 
    AddModuleScore(
    features = test_gs$human_gene_symbol, # need to be a "list"
    name = c("hm_inflam", "hm_ifng", "hm_tnfa"),
    seed = 1212,
    assay = "SCT"
  )
  
  list.seurat.obj.2[[i]] <- seurat.obj.i2
}

FINISH.TIME <- Sys.time() 

FINISH.TIME - START.TIME 
# Time difference of 9.399591 secs


# check/edit --------------------------------------------------------------

list.seurat.obj.2[[1]]@meta.data %>% colnames()
# [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "percent.mt"      "nCount_SCT"      "nFeature_SCT"   
# [7] "batch"           "tissue"          "SCT_snn_res.0.2" "seurat_clusters" "clust"           "hm_inflam1"     
# [13] "hm_ifng2"        "hm_tnfa3"   

for(i in 1:3){
  print(i)
  colnames(list.seurat.obj.2[[i]]@meta.data)[12:14] <- c("hm_inflam", "hm_ifng", "hm_tnfa")
}

list.seurat.obj.2[[1]]@meta.data %>% colnames()
# [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "percent.mt"      "nCount_SCT"      "nFeature_SCT"   
# [7] "batch"           "tissue"          "SCT_snn_res.0.2" "seurat_clusters" "clust"           "hm_inflam"      
# [13] "hm_ifng"         "hm_tnfa" 


# edit the table -------------------------------------------------------------

score.table.merged <- NULL
for(i in 1:3){
  print(i)
  name.i <- c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells")[i]
  
  seurat.obj.i <- list.seurat.obj.2[[i]]
  
  score.table.i <- seurat.obj.i[[]] %>% 
    rownames_to_column("cell.id") %>% 
    as_tibble() %>% 
    dplyr::select(tissue, cell.id, clust, c("hm_inflam", "hm_ifng", "hm_tnfa"))
  
  score.table.merged <- score.table.merged %>% 
    bind_rows(score.table.i)
}


# check -------------------------------------------------------------------

score.table.merged %>% dim()
# [1] 12707     6

score.table.merged %>% head()
# # A tibble: 6 × 6
# tissue cell.id               clust       hm_inflam  hm_ifng hm_tnfa
# <chr>  <chr>                 <chr>           <dbl>    <dbl>   <dbl>
# 1 hfc    HFC1_AAACGAAGTAGCTTAC Tumor Cells  -0.0521  -0.0434   0.0273
# 2 hfc    HFC1_AAACGCTCACGCCAGT Tumor Cells   0.0293   0.0831   0.166 
# 3 hfc    HFC1_AAACGCTGTAGCGCCT Tumor Cells   0.0412   0.0548   0.328 
# 4 hfc    HFC1_AAAGGATCAGAAATTG Tumor Cells  -0.0117  -0.00648 -0.0163
# 5 hfc    HFC1_AAAGGGCAGGTTCATC Tumor Cells   0.00345  0.00184  0.249 
# 6 hfc    HFC1_AAAGGGCGTAGACGGT Tumor Cells   0.0474   0.0519   0.179 


# out ---------------------------------------------------------------------

setwd(dir.2)
out.f <- "001_add_module_score_res_table_merged.tsv"
write_tsv(score.table.merged, out.f)


# stat --------------------------------------------------------------------

# check the normality of data distribution ---------------------------------

res.norm.merged <- NULL
for(i in 1:3){
  print(i)
  clust.i <- c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells")[i]
  
  
  for(j in 1:3){
    print(j)
    param.j <- c("hm_inflam", "hm_ifng", "hm_tnfa")[j]
    
    for(k in 1:2){
      print(k)
      tissue.k <- c("hfc", "lfc")[k]
      
      data <- score.table.merged %>% 
        dplyr::filter(clust == clust.i & tissue == tissue.k) %>% 
        pull(param.j)
      
      # Shapiro-Wilk Test (base R stats)
      if(length(data) <= 5000){
        res.shapiro <- shapiro.test(data)
      }else{
        res.shapiro <- NULL
      }
      
      
      # Anderson-Darling Test (nortest)
      res.ad <- ad.test(data)
      
      # Kolmogorov-Smirnov Test (base R)
      res.ks <- ks.test(data, "pnorm", mean = mean(data), sd = sd(data))
      
      res.ijk <- tibble(
        clust = clust.i, 
        param = param.j, 
        tissue = tissue.k, 
        n = length(data), 
        shapiro.W = res.shapiro$statistic, 
        shapiro.p = res.shapiro$p.value, 
        shapiro.sum = if(is.null(res.shapiro)){"n/a"}else if(res.shapiro$p.value >= 0.05){"ns"}else{"p<.05"}, 
        ad.A = res.ad$statistic,
        ad.p = res.ad$p.value,
        ad.sum = if(res.ad$p.value >= 0.05){"ns"}else{"p<.05"},
        ks.D = res.ks$statistic, 
        ks.p = res.ks$p.value, 
        ks.sum = if(res.ks$p.value >= 0.05){"ns"}else{"p<.05"}
      ) %>% 
        mutate(summary = ifelse( (shapiro.sum == "ns" & ks.sum == "ns") | (ad.sum == "ns" & ks.sum == "ns"), "NS", "check"))
      
      res.norm.merged <- res.norm.merged %>% 
        bind_rows(res.ijk)
    }
  }
}

res.norm.merged <- res.norm.merged %>% 
  dplyr::select(param, clust, tissue, n, 
                shapiro.sum, shapiro.W, shapiro.p, shapiro.sum,  
                ad.A, ad.p, ad.sum, 
                ks.D, ks.p, ks.sum, 
                summary
  )


# check -------------------------------------------------------------------

res.norm.merged %>% dim()
# [1] 18 14

res.norm.merged %>% print(n = Inf)
# # A tibble: 18 × 14
# param     clust          tissue     n shapiro.sum shapiro.W shapiro.p   ad.A     ad.p ad.sum   ks.D     ks.p ks.sum summary
# <chr>     <chr>          <chr>  <int> <chr>           <dbl>     <dbl>  <dbl>    <dbl> <chr>   <dbl>    <dbl> <chr>  <chr>  
# 1 hm_inflam Tumor Cells    hfc     5325 n/a            NA     NA        43.3   3.7 e-24 p<.05  0.0668 0        p<.05  check  
# 2 hm_inflam Tumor Cells    lfc     3212 p<.05           0.939  1.57e-34 36.5   3.7 e-24 p<.05  0.0824 0        p<.05  check  
# 3 hm_ifng   Tumor Cells    hfc     5325 n/a            NA     NA        37.1   3.7 e-24 p<.05  0.0569 2.22e-15 p<.05  check  
# 4 hm_ifng   Tumor Cells    lfc     3212 p<.05           0.930  2.19e-36 29.9   3.7 e-24 p<.05  0.0594 2.79e-10 p<.05  check  
# 5 hm_tnfa   Tumor Cells    hfc     5325 n/a            NA     NA        24.3   3.7 e-24 p<.05  0.0538 7.80e-14 p<.05  check  
# 6 hm_tnfa   Tumor Cells    lfc     3212 p<.05           0.979  1.35e-21 11.3   3.7 e-24 p<.05  0.0333 1.62e- 3 p<.05  check  
# 7 hm_inflam Myeloid Cells  hfc      851 p<.05           0.995  5.21e- 3  0.545 1.61e- 1 ns     0.0191 9.17e- 1 ns     NS     
# 8 hm_inflam Myeloid Cells  lfc     2924 p<.05           0.994  3.11e- 9  5.32  3.89e-13 p<.05  0.0334 2.89e- 3 p<.05  check  
# 9 hm_ifng   Myeloid Cells  hfc      851 p<.05           0.949  1.35e-16  7.74  6.83e-19 p<.05  0.0638 1.98e- 3 p<.05  check  
# 10 hm_ifng   Myeloid Cells  lfc     2924 p<.05           0.975  2.17e-22  6.77  1.36e-16 p<.05  0.0340 2.34e- 3 p<.05  check  
# 11 hm_tnfa   Myeloid Cells  hfc      851 p<.05           0.991  5.54e- 5  1.45  9.64e- 4 p<.05  0.0347 2.57e- 1 ns     check  
# 12 hm_tnfa   Myeloid Cells  lfc     2924 p<.05           0.999  2.68e- 2  0.356 4.58e- 1 ns     0.0134 6.72e- 1 ns     NS     
# 13 hm_inflam Lymphoid Cells hfc       86 ns              0.982  2.80e- 1  0.486 2.21e- 1 ns     0.0739 7.08e- 1 ns     NS     
# 14 hm_inflam Lymphoid Cells lfc      309 p<.05           0.983  8.40e- 4  0.354 4.60e- 1 ns     0.0311 9.26e- 1 ns     NS     
# 15 hm_ifng   Lymphoid Cells hfc       86 p<.05           0.951  2.48e- 3  0.955 1.51e- 2 p<.05  0.0934 4.16e- 1 ns     check  
# 16 hm_ifng   Lymphoid Cells lfc      309 p<.05           0.951  1.38e- 8  3.11  7.97e- 8 p<.05  0.0717 8.35e- 2 ns     check  
# 17 hm_tnfa   Lymphoid Cells hfc       86 ns              0.989  6.70e- 1  0.283 6.26e- 1 ns     0.0560 9.36e- 1 ns     NS     
# 18 hm_tnfa   Lymphoid Cells lfc      309 p<.05           0.978  9.43e- 5  1.20  3.86e- 3 p<.05  0.0489 4.52e- 1 ns     check  


# out ---------------------------------------------------------------------

setwd(dir.2)
out.f <- "002_stat_01_normality_test.tsv"
write_tsv(res.norm.merged, out.f)


# stat test (Wilcoxon rank-sum test - non-parametric) --------------------------------------------------------------------

res.stat.merged <- NULL
for(i in 1:3){
  print(i)
  clust.i <- c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells")[i]
  
  res.stat.i <- NULL
  for(j in 1:3){
    param.j <- c("hm_inflam", "hm_ifng", "hm_tnfa")[j]
    print(j)
    
    
    data <- score.table.merged %>% 
      dplyr::filter(clust == clust.i) %>% 
      dplyr::select(tissue, param.j)
    
    colnames(data)[2] <- "score"
    
    res.ij <- wilcox.test(score ~ tissue, data = data)
    
    res.stat.ij <- data %>% 
      group_by(tissue) %>% 
      summarise(median = median(score)) %>% 
      spread(key = "tissue", value = "median") %>% 
      mutate(clust = clust.i) %>% 
      mutate(pathway = param.j) %>% 
      mutate(stat.method = res.ij$method) %>% 
      mutate(wilcox.p = res.ij$p.value) %>% 
      dplyr::select(clust, pathway, median.HFC = hfc, median.LFC = lfc, stat.method, wilcox.p)
    
    res.stat.i <- res.stat.i %>% 
      bind_rows(res.stat.ij)
  }
  res.stat.i <- res.stat.i %>% 
    mutate(wilcox.p.adj = p.adjust(wilcox.p, method = 'bonferroni')) # NEW: 10/24/2024
  
  res.stat.merged <- res.stat.merged %>% 
    bind_rows(res.stat.i)
}


# check -------------------------------------------------------------------

res.stat.merged %>% dim()
# [1] 9 7

res.stat.merged %>% print(n = Inf)
# # A tibble: 9 × 7
# clust          pathway   median.HFC median.LFC stat.method                                       wilcox.p wilcox.p.adj
# <chr>          <chr>          <dbl>      <dbl> <chr>                                                <dbl>        <dbl>
# 1 Tumor Cells    hm_inflam   -0.0354    -0.0337  Wilcoxon rank sum test with continuity correction 1.08e- 8     3.23e- 8
# 2 Tumor Cells    hm_ifng     -0.00798   -0.00662 Wilcoxon rank sum test with continuity correction 2.04e-11     6.11e-11
# 3 Tumor Cells    hm_tnfa      0.0547     0.0830  Wilcoxon rank sum test with continuity correction 2.08e-38     6.23e-38
# 4 Myeloid Cells  hm_inflam    0.164      0.205   Wilcoxon rank sum test with continuity correction 2.94e-30     8.82e-30
# 5 Myeloid Cells  hm_ifng      0.173      0.249   Wilcoxon rank sum test with continuity correction 1.17e-80     3.51e-80
# 6 Myeloid Cells  hm_tnfa      0.368      0.463   Wilcoxon rank sum test with continuity correction 9.97e-40     2.99e-39
# 7 Lymphoid Cells hm_inflam    0.0494     0.0706  Wilcoxon rank sum test with continuity correction 2.63e- 6     7.88e- 6
# 8 Lymphoid Cells hm_ifng      0.180      0.204   Wilcoxon rank sum test with continuity correction 5.71e- 3     1.71e- 2
# 9 Lymphoid Cells hm_tnfa      0.238      0.313   Wilcoxon rank sum test with continuity correction 8.83e-11     2.65e-10


# out ---------------------------------------------------------------------

setwd(dir.2)
out.f <- "003_stat_02_wilcox_test_results.tsv"
write_tsv(res.stat.merged, out.f)


# Violin plots (Extended Data Fig. S2a, c, e) ---------------------------------------------------------------

for(i in 1:3){
  print(i)

  seurat.obj.i <- list.seurat.obj.2[[i]]
  
  g <- seurat.obj.i %>% 
    VlnPlot(
    features = c("hm_inflam", "hm_ifng", "hm_tnfa"), 
    group.by = "tissue", 
    cols = pal_nejm()(2), 
    pt.size = 0
  )
  g <- g & theme(
    panel.grid.major = element_line(color = "#f0f0f0"), 
    axis.title = element_blank(), 
    axis.text.x = element_text(size = 6), 
    axis.text.y = element_text(size = 8)
  ) 
  # g
  
  name.i <- c("tumor", "myeloid", "lymphoid")[i] 
  setwd(dir.2)
  out.f <- paste0("004_violin_gene_sig_scores_", i, "_", name.i, ".pdf")
  ggsave(out.f, g, width = 12, height = 4)
}


# Feature plots (Extended Data Fig. S2b, d, f) -------------------------------------------------------------

for(i in 1:3){
  print(i)
  
  seurat.obj.i <- list.seurat.obj.2[[i]]
  
  name.i <- c("myeloid", "tumor", "lymphoid")[i]
  name.i2 <- c("Myeloid Cells", "Tumor Cells", "Lymphoid Cells")[i]
  
  
  gg.list <- list()
  for(j in 1:3){
    print(j)
    feature.j <- c("hm_inflam", "hm_ifng", "hm_tnfa")[j]
    
    g <- seurat.obj.i %>% 
      FeaturePlot(
        split.by = "tissue",
        features = feature.j,
        cols = c("lightgrey", wes_palette("Zissou1", 100, type = "continuous")),
        min.cutoff = if(i == 2){'q1'}else{'q0'},
        max.cutoff = if(i == 2){'q97'}else{'q98'},
        pt.size = 0.8,
        keep.scale = "all"
      )
    g <- g & labs(subtitle = feature.j)
    g <- g & theme(
      legend.position = "none",
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
    # g
    
    gg.list[[j]] <- g
  }
  g <- ggarrange(
    plotlist = gg.list, 
    nrow = 3, ncol = 1
  )
  
  name.i <- c("tumor", "myeloid", "lymphoid")[i] 
  setwd(dir.2)
  out.f <- paste0("005_ftrplot_gene_sig_scores_", i, "_", name.i, ".pdf")
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
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
# [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] nortest_1.0-4      msigdbr_7.4.1      wesanderson_0.3.6  ggsci_2.9          ggpubr_0.4.0       SeuratObject_4.1.3 Seurat_4.3.0.1    
# [8] tidylog_1.0.2      forcats_0.5.1      stringr_1.4.0      dplyr_1.0.9        purrr_0.3.4        readr_2.1.2        tidyr_1.2.0       
# [15] tibble_3.1.7       ggplot2_3.3.6      tidyverse_1.3.1   
# 
# loaded via a namespace (and not attached):
# [1] readxl_1.3.1           backports_1.2.1        systemfonts_1.0.4      plyr_1.8.7             igraph_1.3.1          
# [6] lazyeval_0.2.2         sp_2.0-0               splines_4.1.3          listenv_0.8.0          scattermore_0.7       
# [11] digest_0.6.29          htmltools_0.5.5        fansi_1.0.3            magrittr_2.0.3         tensor_1.5            
# [16] cluster_2.1.2          ROCR_1.0-11            tzdb_0.1.2             openxlsx_4.2.4         globals_0.15.0        
# [21] modelr_0.1.8           matrixStats_0.62.0     vroom_1.5.7            timechange_0.2.0       spatstat.sparse_3.0-2 
# [26] colorspace_2.0-3       rvest_1.0.1            ggrepel_0.9.1          textshaping_0.3.6      haven_2.4.3           
# [31] crayon_1.5.1           jsonlite_1.8.0         progressr_0.10.0       spatstat.data_3.0-1    survival_3.2-12       
# [36] zoo_1.8-10             glue_1.6.2             polyclip_1.10-0        gtable_0.3.0           leiden_0.3.9          
# [41] car_3.0-11             future.apply_1.8.1     abind_1.4-5            scales_1.2.0           DBI_1.1.2             
# [46] rstatix_0.7.0          spatstat.random_3.1-5  miniUI_0.1.1.1         Rcpp_1.0.8.3           viridisLite_0.4.0     
# [51] xtable_1.8-4           reticulate_1.25        bit_4.0.4              foreign_0.8-81         clisymbols_1.2.0      
# [56] htmlwidgets_1.6.2      httr_1.4.3             RColorBrewer_1.1-3     ellipsis_0.3.2         ica_1.0-2             
# [61] farver_2.1.0           pkgconfig_2.0.3        uwot_0.1.16            dbplyr_2.1.1           deldir_1.0-6          
# [66] utf8_1.2.2             labeling_0.4.2         tidyselect_1.1.2       rlang_1.1.0            reshape2_1.4.4        
# [71] later_1.3.0            munsell_0.5.0          cellranger_1.1.0       tools_4.1.3            cli_3.6.1             
# [76] generics_0.1.2         broom_0.7.9            ggridges_0.5.3         fastmap_1.1.0          goftest_1.2-2         
# [81] bit64_4.0.5            babelgene_21.4         fs_1.5.2               fitdistrplus_1.1-5     zip_2.2.0             
# [86] RANN_2.6.1             pbapply_1.5-0          future_1.25.0          nlme_3.1-152           mime_0.12             
# [91] ggrastr_1.0.2          xml2_1.3.2             compiler_4.1.3         rstudioapi_0.13        beeswarm_0.4.0        
# [96] plotly_4.10.0          curl_4.3.2             png_0.1-8              ggsignif_0.6.2         spatstat.utils_3.0-3  
# [101] reprex_2.0.1           stringi_1.7.6          lattice_0.20-44        Matrix_1.6-0           vctrs_0.6.1           
# [106] pillar_1.7.0           lifecycle_1.0.3        spatstat.geom_3.2-4    lmtest_0.9-40          RcppAnnoy_0.0.19      
# [111] data.table_1.14.2      cowplot_1.1.1          irlba_2.3.5            httpuv_1.6.5           patchwork_1.1.1       
# [116] R6_2.5.1               promises_1.2.0.1       KernSmooth_2.23-20     gridExtra_2.3          rio_0.5.27            
# [121] vipor_0.4.5            parallelly_1.31.1      codetools_0.2-18       MASS_7.3-54            assertthat_0.2.1      
# [126] withr_2.5.0            sctransform_0.3.5      parallel_4.1.3         hms_1.1.0              grid_4.1.3            
# [131] ggupset_0.3.0          carData_3.0-4          Rtsne_0.16             spatstat.explore_3.2-1 shiny_1.7.1           
# [136] lubridate_1.9.2        ggbeeswarm_0.6.0      


# end ---------------------------------------------------------------------


