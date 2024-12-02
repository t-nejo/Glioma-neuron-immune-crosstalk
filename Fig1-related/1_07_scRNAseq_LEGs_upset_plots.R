# Takahide Nejo, MD, PhD
# Okada Lab, Department of Neurological Surgery, UC San Francisco


# note ---------------------------------------------------------------------

# Nejo T, et al., Glioma-neuronal circuit remodeling induces regional immunosuppression

# Using different versions of packages, such as Seurat, may lead to varying outcomes, which we have not thoroughly validated.

# The RDS object "000_res_fgsea_hallmark_list.rds" is used as input. For details, please refer to the previous step "1_03_scRNAseq_fgsea.R". 


# rm all ------------------------------------------------------------------

rm(list = ls(all.names = T))


# packages ----------------------------------------------------------------

library(tidyverse)
library(tidylog)
library(ggupset) # install.packages("ggupset") 


# check -------------------------------------------------------------------

packageVersion("Seurat")
# [1] ‘4.3.0.1’


# dir ---------------------------------------------------------------------

dir.1 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/1_03_scRNAseq_fgsea/"
dir.2 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/1_07_scRNAseq_LEGs_upset_plots/"
if(dir.exists(dir.2) == F){dir.create(dir.2)}


# load files --------------------------------------------------------------

setwd(dir.1) 
in.f <- "000_res_fgsea_hallmark_list.rds"
list.res.fgsea <- readRDS(in.f)


# check -------------------------------------------------------------------

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


# Extended Data Fig. S3 (upset plots) --------------------------------------------

list.res <- list()
for(i in 1:3){
  print(i) 
  name.i <- c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells")[i]
  
  res.fgsea.i <- list.res.fgsea[[i]] %>% 
    dplyr::filter(pathway %in% imm.pathways)
  
  SYMBOLS <- res.fgsea.i %>% 
    pull(leadingEdge) %>% 
    unlist() %>% 
    unique()
  
  res.i <- tibble(
    symbol = SYMBOLS
  )
  
  for(j in 1:3){
    print(j)
    name.j <- imm.pathways[j]
    # short.name.j <- gsub("HALLMARK", "HM", name.j)
    
    le.genes.j <- res.fgsea.i %>% 
      dplyr::filter(pathway == name.j) %>% 
      pull(leadingEdge) %>% 
      unlist()
    
    res.i <- res.i %>% 
      mutate(new.col = symbol %in% le.genes.j)
      # mutate(new.col = is.element(symbol, le.genes.j))
    
    colnames(res.i)[ncol(res.i)] <- gsub("HALLMARK", "HM", name.j)
  }
  
  # res.i %>% nrow()

  res.i.2 <- res.i %>% 
    gather("pathway", "member", - symbol) %>% 
    dplyr::filter(member) %>% 
    dplyr::select(- member) %>% 
    group_by(symbol) %>%
    summarize(pathway = list(pathway))
  
  g <- ggplot(res.i.2, aes(x = pathway))
  g <- g + geom_bar()
  g <- g + scale_x_upset(order_by = "degree")
  g <- g + theme_classic()
  g <- g + theme(
    panel.grid.major = element_line(color = "#f0f0f0")
  )
  g <- g + theme(plot.margin = margin(0.1, 0.1, 0.1, 2.1, "in"))
  g
  
  name.i <- c("tumor", "myeloid", "lymphoid")[i] 
  setwd(dir.2)
  out.f <- paste0("001_LEG_overlap_in_upset_plot_", i, "_", name.i, ".pdf")
  ggsave(out.f, g, w = 8, h = 4.5)

  list.res[[i]] <- res.i.2
}


# out ---------------------------------------------------------------------

setwd(dir.2)
out.f <- "002_LEG_overlap_summary_table_all.rds"
saveRDS(list.res, out.f)


# edit the table --------------------------------------------------------------------

tmp <- tibble(
  cls = c("INFLAM", "IFNG", "TNFA", "INFLAM|IFNG", "INFLAM|TNFA", "IFNG|TNFA", "INFLAM|IFNG|TNFA"), 
  pathway = list(
    "HM_INFLAMMATORY_RESPONSE", 
    "HM_INTERFERON_GAMMA_RESPONSE",
    "HM_TNFA_SIGNALING_VIA_NFKB", 
    c("HM_INFLAMMATORY_RESPONSE", "HM_INTERFERON_GAMMA_RESPONSE"), 
    c("HM_INFLAMMATORY_RESPONSE", "HM_TNFA_SIGNALING_VIA_NFKB"), 
    c("HM_INTERFERON_GAMMA_RESPONSE", "HM_TNFA_SIGNALING_VIA_NFKB"), 
    c("HM_INFLAMMATORY_RESPONSE", "HM_INTERFERON_GAMMA_RESPONSE", "HM_TNFA_SIGNALING_VIA_NFKB")
  )
)

list.res.2 <- list()
for(i in 1:3){
  print(i) 
  name.i <- c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells")[i]
  
  res.i <- list.res[[i]]
  
  res.i <- res.i %>% 
    left_join(tmp, by = "pathway") %>% 
    dplyr::select(- pathway)
  
  res.i$cls <- factor(res.i$cls, levels = tmp$cls)
  
  res.i <- res.i %>% 
    group_by(cls) %>%
    summarize(symbols = paste0(symbol, collapse = "|")) %>% 
    arrange(cls)
  
  list.res.2[[i]] <- res.i
}


# check --------------------------------------------------------------

for(i in 1:3){
  print(i)
  name.i <- c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells")[i]
  
  print(name.i)
  list.res.2[[i]] %>% 
    dplyr::filter(cls == "INFLAM|IFNG|TNFA") %>% 
    pull(symbols) %>% 
    print()
}

# [1] 1
# [1] "Tumor Cells"
# [1] "CCL2|CD69|IRF1|NFKB1|TNFAIP6"
# [1] 2
# [1] "Myeloid Cells"
# [1] "CCL2|CCL5|CD69|CDKN1A|CXCL10|IL6|NFKB1|NFKBIA"
# [1] 3
# [1] "Lymphoid Cells"
# [1] "CCL2|CDKN1A|IRF1|NAMPT|NFKB1|PDE4B|RIPK2"


# out ---------------------------------------------------------------------

for(i in 1:3){
  name.i <- c("tumor", "myeloid", "lymphoid")[i] 
  setwd(dir.2)
  out.f <- paste0("003_LEG_overlap_summary_table_", i, "_", name.i, ".tsv")
  write_tsv(list.res.2[[i]], out.f)
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
# [1] ggupset_0.3.0   tidylog_1.0.2   forcats_0.5.1   stringr_1.4.0   dplyr_1.0.9     purrr_0.3.4     readr_2.1.2    
# [8] tidyr_1.2.0     tibble_3.1.7    ggplot2_3.3.6   tidyverse_1.3.1
# 
# loaded via a namespace (and not attached):
# [1] Seurat_4.3.0.1         Rtsne_0.16             colorspace_2.0-3       deldir_1.0-6           ellipsis_0.3.2        
# [6] ggridges_0.5.3         fs_1.5.2               rstudioapi_0.13        spatstat.data_3.0-1    farver_2.1.0          
# [11] leiden_0.3.9           listenv_0.8.0          bit64_4.0.5            ggrepel_0.9.1          fansi_1.0.3           
# [16] lubridate_1.9.2        xml2_1.3.2             codetools_0.2-18       splines_4.1.3          polyclip_1.10-0       
# [21] jsonlite_1.8.0         broom_0.7.9            ica_1.0-2              cluster_2.1.2          dbplyr_2.1.1          
# [26] png_0.1-8              uwot_0.1.16            shiny_1.7.1            sctransform_0.3.5      spatstat.sparse_3.0-2 
# [31] compiler_4.1.3         httr_1.4.3             backports_1.2.1        assertthat_0.2.1       SeuratObject_4.1.3    
# [36] Matrix_1.6-0           fastmap_1.1.0          lazyeval_0.2.2         cli_3.6.1              later_1.3.0           
# [41] htmltools_0.5.5        tools_4.1.3            igraph_1.3.1           gtable_0.3.0           glue_1.6.2            
# [46] RANN_2.6.1             reshape2_1.4.4         Rcpp_1.0.8.3           scattermore_0.7        cellranger_1.1.0      
# [51] vctrs_0.6.1            spatstat.explore_3.2-1 nlme_3.1-152           progressr_0.10.0       lmtest_0.9-40         
# [56] spatstat.random_3.1-5  globals_0.15.0         rvest_1.0.1            timechange_0.2.0       mime_0.12             
# [61] miniUI_0.1.1.1         lifecycle_1.0.3        irlba_2.3.5            goftest_1.2-2          future_1.25.0         
# [66] MASS_7.3-54            zoo_1.8-10             scales_1.2.0           vroom_1.5.7            clisymbols_1.2.0      
# [71] hms_1.1.0              promises_1.2.0.1       spatstat.utils_3.0-3   parallel_4.1.3         RColorBrewer_1.1-3    
# [76] reticulate_1.25        pbapply_1.5-0          gridExtra_2.3          stringi_1.7.6          systemfonts_1.0.4     
# [81] rlang_1.1.0            pkgconfig_2.0.3        matrixStats_0.62.0     lattice_0.20-44        ROCR_1.0-11           
# [86] tensor_1.5             labeling_0.4.2         patchwork_1.1.1        htmlwidgets_1.6.2      bit_4.0.4             
# [91] cowplot_1.1.1          tidyselect_1.1.2       parallelly_1.31.1      RcppAnnoy_0.0.19       plyr_1.8.7            
# [96] magrittr_2.0.3         R6_2.5.1               generics_0.1.2         DBI_1.1.2              withr_2.5.0           
# [101] pillar_1.7.0           haven_2.4.3            fitdistrplus_1.1-5     survival_3.2-12        abind_1.4-5           
# [106] sp_2.0-0               future.apply_1.8.1     modelr_0.1.8           crayon_1.5.1           KernSmooth_2.23-20    
# [111] utf8_1.2.2             spatstat.geom_3.2-4    plotly_4.10.0          tzdb_0.1.2             readxl_1.3.1          
# [116] grid_4.1.3             data.table_1.14.2      reprex_2.0.1           digest_0.6.29          xtable_1.8-4          
# [121] httpuv_1.6.5           textshaping_0.3.6      munsell_0.5.0          viridisLite_0.4.0     


# end ---------------------------------------------------------------------


