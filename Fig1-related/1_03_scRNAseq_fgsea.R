# Takahide Nejo, MD, PhD
# Okada Lab, Department of Neurological Surgery, UC San Francisco


# note ---------------------------------------------------------------------

# Nejo T, et al., Glioma-neuronal circuit remodeling induces regional immunosuppression
# The RDS object "002_res_findmarkers_mast_list.rds" is used as an input Seurat object. For details, please refer to the previous step "1_02_scRNAseq_dge.R". 

# Using different versions of packages, such as Seurat, may lead to varying outcomes, which we have not thoroughly validated.

# Please note that due to the random permutations involved, multiple runs of fgsea may yield slightly different results each time. This behavior is documented in the following discussions: https://github.com/ctlab/fgsea/issues/12 and https://support.bioconductor.org/p/9158125/. 


# rm all ------------------------------------------------------------------

rm(list = ls(all.names = T))


# packages ----------------------------------------------------------------

library(tidyverse)
library(tidylog)
library(msigdbr)
library(fgsea)
library(ggpubr)
library(ggsci)


# check -------------------------------------------------------------------

packageVersion("msigdbr")
# [1] ‘7.4.1’


# dir ---------------------------------------------------------------------

dir.1 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/1_02_scRNAseq_dge/"
dir.2 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/1_03_scRNAseq_fgsea/"
if(dir.exists(dir.2) == F){dir.create(dir.2)}


# load files --------------------------------------------------------------

START.TIME <- Sys.time() 

setwd(dir.1) 
in.f <- "002_res_findmarkers_mast_list.rds"
list.res.mast <- readRDS(in.f)

FINISH.TIME <- Sys.time() 

FINISH.TIME - START.TIME 
# Time difference of 0.09501505 secs


# check -------------------------------------------------------------------

list.res.mast %>% length()
# [1] 4

list.res.mast %>% names()
# [1] "Tumor Cells"    "Myeloid Cells"  "Lymphoid Cells" "Astrocytes"    

list.res.mast[["Tumor Cells"]] %>% head()
# p_val avg_log2FC pct.1 pct.2     p_val_adj
# CCL3    0.000000e+00 -2.3478363 0.219 0.640  0.000000e+00
# CCL4    0.000000e+00 -2.7553522 0.229 0.647  0.000000e+00
# CCL4L2 3.392314e-255 -1.4319096 0.145 0.462 8.200240e-251
# XIST   6.344140e-217 -1.5844814 0.083 0.347 1.533569e-212
# POSTN  1.081059e-207  1.2707088 0.382 0.099 2.613245e-203
# GAPDH  1.376818e-172  0.7879826 0.988 0.946 3.328182e-168


# prep the rank objects -----------------------------------------------------------

list.rank <- list()
for(i in 1:4){
  print(i)
  
  res.mast.i <- list.res.mast[[i]] %>% 
    rownames_to_column("symbol") %>% 
    as_tibble() %>% 
    arrange(desc(avg_log2FC)) # ref: https://github.com/ctlab/fgsea/issues/50
  
  rank.i <- res.mast.i %>% 
    pull(avg_log2FC)
  
  names(rank.i) <- res.mast.i %>% 
    pull(symbol)
  
  
  # # check
  # rank.i %>% summary() %>% print() 
  # res.mast.i %>% dplyr::filter(avg_log2FC == 0) %>% nrow()

  list.rank[[i]] <- rank.i
}

names(list.rank) <- c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells", "Astrocytes")


# MSigDB ------------------------------------------------------------------

# ref: https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
# ref: http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

msigdbr_collections() %>% print(n = Inf)
{
  # # A tibble: 23 × 3
  #    gs_cat gs_subcat         num_genesets
  #    <chr>  <chr>                    <int>
  #  1 C1     ""                         278
  #  2 C2     "CGP"                     3368
  #  3 C2     "CP"                        29
  #  4 C2     "CP:BIOCARTA"              292
  #  5 C2     "CP:KEGG"                  186
  #  6 C2     "CP:PID"                   196
  #  7 C2     "CP:REACTOME"             1604
  #  8 C2     "CP:WIKIPATHWAYS"          615
  #  9 C3     "MIR:MIR_Legacy"           221
  # 10 C3     "MIR:MIRDB"               2377
  # 11 C3     "TFT:GTRD"                 523
  # 12 C3     "TFT:TFT_Legacy"           610
  # 13 C4     "CGN"                      427
  # 14 C4     "CM"                       431
  # 15 C5     "GO:BP"                   7481
  # 16 C5     "GO:CC"                    996
  # 17 C5     "GO:MF"                   1708
  # 18 C5     "HPO"                     4813
  # 19 C6     ""                         189
  # 20 C7     "IMMUNESIGDB"             4872
  # 21 C7     "VAX"                      347
  # 22 C8     ""                         671
  # 23 H      ""                          50
}


# note --------------------------------------------------------------------

# focus on following categories: H (hallmark)



# prep for the analysis ---------------------------------------------------

test.gs <- msigdbr(
  species = "Homo sapiens", 
  category = "H", 
  subcategory = ""
)

# check 
test.gs %>% pull(gs_name) %>% unique() %>% length() %>% print()
# 50
# test.gs %>% head() %>% print() ;

test.gs.list <- split(x = test.gs$gene_symbol, f = test.gs$gs_name)

test.gs.list %>% length()
# [1] 50


# main analysis part ------------------------------------------------------

START.TIME <- Sys.time() 

list.res.fgsea <- list()
for(i in 1:4){
  print(i)
  
  rank.i <- list.rank[[i]]
  
  res.fgsea.i <- fgseaMultilevel(
    pathways = test.gs.list, 
    rank.i, 
    minSize = 10, 
    maxSize = 500, 
    eps = 0
  )
    
  res.fgsea.i <- bind_rows(  
    res.fgsea.i %>% tibble() %>% dplyr::filter(ES > 0) %>% arrange(padj, desc(NES)), 
    res.fgsea.i %>% tibble() %>% dplyr::filter(ES < 0) %>% arrange(desc(padj), desc(NES))
  )
  
  list.res.fgsea[[i]] <- res.fgsea.i
}

names(list.res.fgsea) <- c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells", "Astrocytes")


FINISH.TIME <- Sys.time() 

FINISH.TIME - START.TIME 
# Time difference of 4.092789 secs


# Supplementary Table 1 --------------------------------------------------

for(i in 1:4){
  print(i)
  res.fgsea.i <- list.res.fgsea[[i]]
  
  name.i <- c("tumor", "myeloid", "lymphoid", "astro")[i] 
  setwd(dir.2)
  out.f <- paste0("001_res_fgsea_hallmark_", i, "_", name.i, ".tsv")
  data.table::fwrite(res.fgsea.i, out.f, sep="\t", sep2=c("", "|", ""))
}

# Fig 1 a-c & Extended Data Fig. S1b (bar plot, version 1) ----------------------------------------------------------------

for(i in 1:4){
  print(i)
  
  res.fgsea.i <- list.res.fgsea[[i]]
  
  res.fgsea.i <- res.fgsea.i %>% 
    mutate(NES = round(NES, digits = 2)) %>% 
    mutate(pathway = gsub("HALLMARK_", "", pathway)) %>% 
    mutate(pathway = paste0(substr(pathway, 1, 40), ". (NES) ", NES)) %>% 
    mutate(minus.log10.padj = -log10(padj)) %>% 
    dplyr::select(pathway, NES, padj, minus.log10.padj)
  
  res.pos <- res.fgsea.i %>% 
    dplyr::filter(NES > 0) %>% 
    arrange(padj) %>% 
    head(n = 6) %>%
    mutate(lab = "up")
  res.neg <- res.fgsea.i %>% 
    dplyr::filter(NES < 0) %>% 
    arrange(desc(padj)) %>% 
    tail(n = 6) %>%
    mutate(lab = "down")
  
  res.i <- bind_rows(res.pos, res.neg)
  
  res.i$pathway <- factor(res.i$pathway, levels = rev(as.character(res.i$pathway)))
  res.i$lab <- factor(res.i$lab, levels = c("up", "down"))
  
  TITLE <- c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells", "Astrocytes")[i]
  
  g <- ggplot(res.i, aes(x = pathway, y = minus.log10.padj, fill = lab)) 
  g <- g + geom_bar(stat = "identity", color = "black", alpha = 0.8) 
  g <- g + coord_flip()
  g <- g + scale_fill_nejm()
  g <- g + labs(title = TITLE, x = NULL, y = "-log10(P.adj)")
  g <- g + theme_classic()
  g <- g + theme(
    panel.grid.major = element_line(color = "#f0f0f0"), 
    legend.position = "none"
  )
  # plot(g)
  
  name.i <- c("tumor", "myeloid", "lymphoid", "astro")[i] 
  setwd(dir.2)
  out.f <- paste0("002_barplot_res_fgsea_hm_top_n_bottom_6_", i, "_", name.i, ".pdf")
  ggsave(out.f, g, width = 8, height = 6) ;   
}


# bar plot (version 2. more comprehensive) ----------------------------------------------------------------

for(i in 1:4){
  print(i)
  
  res.fgsea.i <- list.res.fgsea[[i]]
  
  res.fgsea.i <- res.fgsea.i %>% 
    mutate(NES = round(NES, digits = 2)) %>% 
    mutate(pathway = gsub("HALLMARK_", "", pathway)) %>% 
    mutate(pathway = paste0(substr(pathway, 1, 40), ". (NES) ", NES)) %>% 
    mutate(minus.log10.padj = -log10(padj)) %>% 
    dplyr::select(pathway, NES, padj, minus.log10.padj)
  
  res.pos <- res.fgsea.i %>% 
    dplyr::filter(NES > 0) %>% 
    arrange(padj) %>% 
    mutate(lab = "up")
  res.neg <- res.fgsea.i %>% 
    dplyr::filter(NES < 0) %>% 
    arrange(desc(padj)) %>% 
    mutate(lab = "down")
  
  res.i <- bind_rows(res.pos, res.neg)
  
  res.i$pathway <- factor(res.i$pathway, levels = rev(as.character(res.i$pathway)))
  res.i$lab <- factor(res.i$lab, levels = c("up", "down"))
  
  TITLE <- c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells", "Astrocytes")[i]
  
  g <- ggplot(res.i, aes(x = pathway, y = minus.log10.padj, fill = lab)) 
  g <- g + geom_bar(stat = "identity", color = "black", alpha = 0.8) 
  g <- g + coord_flip()
  g <- g + scale_fill_nejm()
  g <- g + labs(title = TITLE, x = NULL, y = "-log10(P.adj)")
  g <- g + theme_classic()
  g <- g + theme(
    panel.grid.major = element_line(color = "#f0f0f0"), 
    legend.position = "none"
  )
  # plot(g)
  
  name.i <- c("tumor", "myeloid", "lymphoid", "astro")[i] 
  setwd(dir.2)
  out.f <- paste0("003_barplot_res_fgsea_hm_whole_", i, "_", name.i, ".pdf")
  ggsave(out.f, g, width = 10, height = 10) ;   
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
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] fgsea_1.18.0       msigdbr_7.4.1      ggrepel_0.9.1      cowplot_1.1.1     
# [5] glmGamPoi_1.4.0    sctransform_0.3.5  ggpubr_0.4.0       ggsci_2.9         
# [9] RColorBrewer_1.1-3 harmony_0.1.0      Rcpp_1.0.8.3       tidylog_1.0.2     
# [13] forcats_0.5.1      stringr_1.4.0      dplyr_1.0.9        purrr_0.3.4       
# [17] readr_2.1.2        tidyr_1.2.0        tibble_3.1.7       ggplot2_3.3.6     
# [21] tidyverse_1.3.1    SeuratObject_4.1.3 Seurat_4.3.0.1     SPATA2_2.0.4      
# 
# loaded via a namespace (and not attached):
# [1] utf8_1.2.2                  spatstat.explore_3.2-1     
# [3] reticulate_1.25             tidyselect_1.1.2           
# [5] htmlwidgets_1.6.2           BiocParallel_1.26.2        
# [7] grid_4.1.3                  Rtsne_0.16                 
# [9] munsell_0.5.0               codetools_0.2-18           
# [11] ica_1.0-2                   units_0.8-0                
# [13] future_1.25.0               miniUI_0.1.1.1             
# [15] withr_2.5.0                 spatstat.random_3.1-5      
# [17] colorspace_2.0-3            progressr_0.10.0           
# [19] Biobase_2.52.0              rstudioapi_0.13            
# [21] stats4_4.1.3                SingleCellExperiment_1.14.1
# [23] ROCR_1.0-11                 ggsignif_0.6.2             
# [25] tensor_1.5                  listenv_0.8.0              
# [27] labeling_0.4.2              MatrixGenerics_1.4.3       
# [29] GenomeInfoDbData_1.2.6      polyclip_1.10-0            
# [31] farver_2.1.0                parallelly_1.31.1          
# [33] vctrs_0.6.1                 generics_0.1.2             
# [35] timechange_0.2.0            R6_2.5.1                   
# [37] GenomeInfoDb_1.28.1         locfit_1.5-9.4             
# [39] bitops_1.0-7                spatstat.utils_3.0-3       
# [41] DelayedArray_0.18.0         assertthat_0.2.1           
# [43] promises_1.2.0.1            scales_1.2.0               
# [45] gtable_0.3.0                globals_0.15.0             
# [47] goftest_1.2-2               rlang_1.1.0                
# [49] clisymbols_1.2.0            systemfonts_1.0.4          
# [51] splines_4.1.3               rstatix_0.7.0              
# [53] lazyeval_0.2.2              spatstat.geom_3.2-4        
# [55] broom_0.7.9                 reshape2_1.4.4             
# [57] abind_1.4-5                 modelr_0.1.8               
# [59] backports_1.2.1             httpuv_1.6.5               
# [61] tools_4.1.3                 ellipsis_0.3.2             
# [63] BiocGenerics_0.38.0         ggridges_0.5.3             
# [65] plyr_1.8.7                  zlibbioc_1.38.0            
# [67] RCurl_1.98-1.3              deldir_1.0-6               
# [69] pbapply_1.5-0               S4Vectors_0.30.2           
# [71] zoo_1.8-10                  SummarizedExperiment_1.22.0
# [73] haven_2.4.3                 cluster_2.1.2              
# [75] fs_1.5.2                    magrittr_2.0.3             
# [77] data.table_1.14.2           scattermore_0.7            
# [79] openxlsx_4.2.4              lmtest_0.9-40              
# [81] reprex_2.0.1                RANN_2.6.1                 
# [83] fitdistrplus_1.1-5          anndata_0.7.5.6            
# [85] matrixStats_0.62.0          hms_1.1.0                  
# [87] patchwork_1.1.1             mime_0.12                  
# [89] fftwtools_0.9-11            xtable_1.8-4               
# [91] rio_0.5.27                  jpeg_0.1-9                 
# [93] readxl_1.3.1                IRanges_2.26.0             
# [95] gridExtra_2.3               compiler_4.1.3             
# [97] KernSmooth_2.23-20          crayon_1.5.1               
# [99] SPATAData_0.0.0.9000        htmltools_0.5.5            
# [101] later_1.3.0                 tzdb_0.1.2                 
# [103] tiff_0.1-11                 lubridate_1.9.2            
# [105] DBI_1.1.2                   dbplyr_2.1.1               
# [107] MASS_7.3-54                 babelgene_21.4             
# [109] Matrix_1.6-0                car_3.0-11                 
# [111] cli_3.6.1                   parallel_4.1.3             
# [113] igraph_1.3.1                GenomicRanges_1.44.0       
# [115] pkgconfig_2.0.3             foreign_0.8-81             
# [117] sp_2.0-0                    confuns_1.0.2              
# [119] plotly_4.10.0               spatstat.sparse_3.0-2      
# [121] xml2_1.3.2                  XVector_0.32.0             
# [123] rvest_1.0.1                 digest_0.6.29              
# [125] RcppAnnoy_0.0.19            spatstat.data_3.0-1        
# [127] fastmatch_1.1-3             cellranger_1.1.0           
# [129] leiden_0.3.9                uwot_0.1.16                
# [131] curl_4.3.2                  shiny_1.7.1                
# [133] EBImage_4.34.0              lifecycle_1.0.3            
# [135] nlme_3.1-152                jsonlite_1.8.0             
# [137] carData_3.0-4               viridisLite_0.4.0          
# [139] fansi_1.0.3                 pillar_1.7.0               
# [141] lattice_0.20-44             fastmap_1.1.0              
# [143] httr_1.4.3                  survival_3.2-12            
# [145] glue_1.6.2                  zip_2.2.0                  
# [147] png_0.1-8                   stringi_1.7.6              
# [149] textshaping_0.3.6           irlba_2.3.5                
# [151] future.apply_1.8.1         


# end ---------------------------------------------------------------------


