# Takahide Nejo, MD, PhD
# Okada Lab, Department of Neurological Surgery, UC San Francisco


# note ---------------------------------------------------------------------

# Nejo T, et al., Glioma-neuronal circuit remodeling induces regional immunosuppression

# Using different versions of packages, such as Seurat, may lead to varying outcomes, which we have not thoroughly validated.

# The RDS object "002_gbm_neuroimmune_myeloid_subset_res0.2.rds" is used as inputs. For details, please refer to the previous step "1_08_scRNAseq_myeloid_subset.R". 


# rm all ------------------------------------------------------------------

rm(list = ls(all.names = T))


# packages ----------------------------------------------------------------

library(tidyverse)
library(tidylog)
library(Seurat)


# check -------------------------------------------------------------------

packageVersion("Seurat")
# [1] ‘4.3.0.1’
packageVersion("sctransform")
# [1] ‘0.3.5’


# dir ---------------------------------------------------------------------

dir.1 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/1_08_scRNAseq_myeloid_subset/"
dir.2 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/1_10_scRNAseq_myeloid_characterization_2/"
if(dir.exists(dir.2) == F){dir.create(dir.2)}


# load file(s) --------------------------------------------------------------

START.TIME <- Sys.time() 

setwd(dir.1) 
in.f <- "002_gbm_neuroimmune_myeloid_subset_res0.2.rds"
seurat.obj <- readRDS(in.f)

FINISH.TIME <- Sys.time() 

FINISH.TIME - START.TIME 
# Time difference of 8.720711 secs


# check -------------------------------------------------------------------

seurat.obj %>% print()
# An object of class Seurat 
# 42324 features across 3775 samples within 2 assays 
# Active assay: SCT (18151 features, 3000 variable features)
# 1 other assay present: RNA
# 3 dimensional reductions calculated: pca, harmony, umap

seurat.obj@meta.data %>% colnames()
# [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "percent.mt"      "batch"          
# [6] "tissue"          "nCount_SCT"      "nFeature_SCT"    "seurat_clusters" "SCT_snn_res.0.2"


# assign new cell annotations --------------------------------------------------------

new.name.table <- tibble(
  seurat_clusters = factor(0:5), 
  seurat_clusters.new = c("pro-inflammatory", "undetermined", "anti-inflammatory", "anti-inflammatory", "pro-inflammatory", "anti-inflammatory")
)

metadata.old <- seurat.obj@meta.data

metadata.new <- metadata.old %>% 
  as_tibble() %>% 
  left_join(new.name.table, by = "seurat_clusters") %>% 
  as.data.frame()

# left_join: added one column (seurat_clusters.new)
# > rows only in x       0
# > rows only in y  (    0)
# > matched rows     3,775
# >                 =======
# > rows total       3,775

rownames(metadata.new) <- rownames(metadata.old)

seurat.obj@meta.data <- metadata.new

seurat.obj@meta.data$seurat_clusters.new <- 
  factor(seurat.obj@meta.data$seurat_clusters.new,
         levels = c("pro-inflammatory", "anti-inflammatory", "undetermined")
  )

Idents(seurat.obj) <- "seurat_clusters.new"


# add Muller Mg/Mo-TAM module scores -------------------------------------------------------

Muller_MgTAM_features <- c("P2RY12", "CX3CR1", "NAV3", "SIGLEC8", "SLC1A3")
Muller_MoTAM_features <- c("TGFBI", "ITGA4", "IFITM2", "FPR3", "S100A11", "KYNU")

list.TAM.features <- list(Muller_MgTAM_features, Muller_MoTAM_features)
names(list.TAM.features) <- c("Muller_MgTAM", "Muller_MoTAM")

seurat.obj <- seurat.obj %>% 
  AddModuleScore(
  features = list.TAM.features,
  name = names(list.TAM.features),
  seed = 1212,
  assay = "SCT"
)


# check -------------------------------------------------------------------

seurat.obj@meta.data %>%
  as_tibble() %>%
  dplyr::select(Muller_MgTAM1, Muller_MoTAM2) %>%
  head()
# # A tibble: 6 × 2
# Muller_MgTAM1 Muller_MoTAM2
# <dbl>         <dbl>
# 1        0.500       -0.00642
# 2        0.262       -0.273
# 3       -0.306        0.267
# 4        0.438        0.0688
# 5        0.363       -0.173
# 6       -0.0851       0.159



# feature plot (Fig. 1g) --------------------------------------------------

g <- seurat.obj %>% 
  FeaturePlot(
    features = c(
      "Muller_MgTAM1",
      "Muller_MoTAM2"
    ),
    pt.size = 1.2,
    slot = "scale.data",
    max.cutoff = 0.5,
    min.cutoff = 0,
    cols = c("gray95", "deepskyblue4"),
    ncol = 2
  ) & theme(
    legend.position = "right",
    legend.direction = "vertical",   panel.grid.major = element_line(color = "#f0f0f0"),
    axis.text = element_blank(),
    axis.title = element_blank()
  )
plot(g)

setwd(dir.2)
out.f <- "001_ftrplot_myeloid_MgTAM_n_MoTAM_scores_Muller.pdf"
ggsave(out.f, g, width = 10, height = 5)


# Visualize Mg/Mo-TAM ------------------------------------------------------------

metadata.tmp <- seurat.obj@meta.data %>%
  rownames_to_column("cell.id") %>%
  tibble() %>%
  mutate(Muller_cl = ifelse(Muller_MgTAM1 < 0 & Muller_MoTAM2 < 0, "undetermined",
                            ifelse(Muller_MgTAM1 > Muller_MoTAM2, "Mg-TAM", "Mo-TAM")))

metadata.tmp.2 <- metadata.tmp %>%
  dplyr::select(- cell.id) %>%
  as.data.frame()

rownames(metadata.tmp.2) <- metadata.tmp$cell.id
seurat.obj@meta.data <- metadata.tmp.2

seurat.obj@meta.data$Muller_cl <- factor(seurat.obj@meta.data$Muller_cl,
                                         levels = c("Mo-TAM", "Mg-TAM", "undetermined")
)



# dimplot (Fig. 1h) -------------------------------------------------------

plots <- seurat.obj %>% 
  DimPlot(
    group.by = "Muller_cl",
    split.by = "tissue",
    combine = FALSE,
    cols = c("#44AA99", "#AA4499", "#DDCC77"),
    pt.size = .8
  )
plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "bottom")
                + guides(color = guide_legend(ncol = 3, byrow = FALSE, override.aes = list(size = 4))))
g <- CombinePlots(plots, ncol = 1)
g <- g & theme(
  panel.grid.major = element_line(color = "#f0f0f0"),
  axis.text = element_blank(),
  axis.title = element_blank()
)
g <- g + labs(title = NULL, x = "UMAP-1", y = "UMAP-2")
g

setwd(dir.2)
out.f <- "002_dimplot_myeloid_MgTAM_n_MoTAM_Muller_cl.pdf"
ggsave(out.f, g, width = 10, height = 5)


n.table.2 <- metadata.tmp.2 %>%
  dplyr::select(tissue, Muller_cl) %>%
  table() %>%
  as_tibble()

n.table.2 <- n.table.2 %>%
  arrange(tissue) %>%
  mutate(total = ifelse(tissue == "hfc", 851, 2924)) %>%
  mutate(pct = n/total)

n.table.2$Muller_cl <- factor(n.table.2$Muller_cl, levels = c("Mo-TAM", "Mg-TAM", "undetermined"))

n.table.2


g <- ggplot(n.table.2, aes(x = tissue, y = pct, color = Muller_cl)) ;
g <- g + geom_point(size = 1.5)
g <- g + geom_line(aes(group = Muller_cl), size = 1)
g <- g + theme_classic()
g <- g + scale_color_manual(values = c("#44AA99", "#AA4499", "#DDCC77"))
g <- g + ylim(0, 0.8)
g <- g & theme(
  panel.grid.major = element_line(color = "#f0f0f0"),
  legend.position = "right"
)
g <- g + labs(title = NULL, x = NULL, y = "pct", color = NULL) ;
g

setwd(dir.2)
out.f <- "003_lineplot_myeloid_MgTAM_n_MoTAM_Muller_cl.pdf"
ggsave(out.f, g, w = 6, h = 4)



# Bar plot (Extended Data Fig. S4e)----------------------------------------------------------------

n.table.2$Muller_cl <- factor(n.table.2$Muller_cl, levels = rev(c("Mo-TAM", "Mg-TAM", "undetermined")))

g <- ggplot(n.table.2, aes(x = tissue, y = n, fill = Muller_cl)) ;
g <- g + geom_bar(stat = "identity", position = "fill", width = 0.75, color = "black", alpha = 0.9)
g <- g + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 25, 50, 75, 100))
g <- g + scale_fill_manual(values = rev(c("#44AA99", "#AA4499", "#DDCC77")))
g <- g + labs(title = NULL, x = NULL, y = "fraction (%)", fill = NULL) ;
g <- g + theme_classic()
g <- g & theme(
  panel.grid.major = element_line(color = "#f0f0f0"),
  legend.position = "bottom"
)
g

setwd(dir.2)
out.f <- "004_bar_plot_myeloid_MgTAM_n_MoTAM_Muller_cl_hfc_lfc.pdf"
ggsave(out.f, g, width = 4, height = 3)



# Bar plot (Extended Data Fig. S4f)----------------------------------------------------------------

# per patient

metadata.4 <- seurat.obj@meta.data %>% 
  dplyr::select(orig.ident, Muller_cl) %>% 
  table() %>% 
  as_tibble() %>% 
  mutate(tissue = ifelse(grepl("hfc", orig.ident), "HFC", "LFC")) %>% 
  mutate(pt = ifelse(grepl(1, orig.ident), "pt-1", 
                     ifelse(grepl(2, orig.ident), "pt-2", "pt-3")))

metadata.4$tissue <- factor(metadata.4$tissue, levels = c("HFC", "LFC"))
metadata.4$Muller_cl <- 
  factor(metadata.4$Muller_cl, 
         levels = rev(c("Mo-TAM", "Mg-TAM", "undetermined")))

COLORS <- rev(c("#44AA99", "#AA4499", "#DDCC77"))


g <- ggplot(metadata.4, aes(x = tissue, y = n, fill = Muller_cl))
g <- g + geom_bar(position = "fill", stat = "identity", color = "black")
g <- g + scale_fill_manual(values = COLORS)
g <- g + facet_wrap(~ pt, ncol = 3)
g <- g + labs(title = NULL, x = NULL, y = "cell count", fill = NULL, shape = NULL)
g <- g + theme_classic()
g <- g + theme(
  panel.grid.major = element_line(color = "#f0f0f0"), 
  plot.title = element_text(hjust = 0.5, face = "bold"), 
  strip.background = element_blank(), 
  strip.text = element_text(size = 12), 
  legend.position = "bottom", 
)
plot(g)

setwd(dir.2)
out.f <- "005_bar_plot_myeloid_MgTAM_n_MoTAM_Muller_cl_hfc_lfc_per_patient.pdf"
ggsave(out.f, g, width = 8, height = 4.5) 


# Line plot (Extended Data Fig. S4g) --------------------------------------------------

# per patient

sum.4 <- metadata.4 %>% 
  group_by(orig.ident) %>% 
  summarise(total = sum(n)) %>% 
  ungroup()

metadata.4b <- metadata.4 %>% 
  left_join(sum.4, by = "orig.ident") %>% 
  mutate(pct = n / total * 100) %>% 
  arrange(pt, tissue, desc(Muller_cl)) %>% 
  dplyr::filter(Muller_cl != "undetermined") %>% 
  dplyr::select(pt, tissue, Muller_cl, n, total, pct)


COLORS <- rev(c("#44AA99", "#AA4499"))

g <- ggplot(metadata.4b, aes(x = tissue, y = pct, color = Muller_cl))
g <- g + geom_point(size = 1)
g <- g + geom_line(aes(group = Muller_cl), size = 1)
g <- g + theme_classic()
g <- g + scale_color_manual(values = COLORS)
g <- g + ylim(0, 100)
g <- g + facet_wrap(~ pt)
g <- g & theme(
  panel.grid.major = element_line(color = "#f0f0f0"),
  legend.position = "none",    
  strip.background = element_blank(), 
  strip.text = element_text(size = 12), 
)
g <- g + labs(title = NULL, x = NULL, y = "pct", color = NULL) ;
g

setwd(dir.2)
out.f <- "006_line_plot_myeloid_MgTAM_n_MoTAM_Muller_cl_hfc_lfc_per_patient.pdf"
ggsave(out.f, g, width = 8, height = 4.5) 


# si ----------------------------------------------------------------------

Sys.time() %>% print()
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
# [1] SeuratObject_4.1.3 Seurat_4.3.0.1     tidylog_1.0.2      forcats_0.5.1      stringr_1.4.0     
# [6] dplyr_1.0.9        purrr_0.3.4        readr_2.1.2        tidyr_1.2.0        tibble_3.1.7      
# [11] ggplot2_3.3.6      tidyverse_1.3.1   
# 
# loaded via a namespace (and not attached):
# [1] readxl_1.3.1           backports_1.2.1        systemfonts_1.0.4      plyr_1.8.7            
# [5] igraph_1.3.1           lazyeval_0.2.2         sp_2.0-0               splines_4.1.3         
# [9] listenv_0.8.0          scattermore_0.7        digest_0.6.29          htmltools_0.5.5       
# [13] fansi_1.0.3            magrittr_2.0.3         tensor_1.5             cluster_2.1.2         
# [17] ROCR_1.0-11            tzdb_0.1.2             globals_0.15.0         modelr_0.1.8          
# [21] matrixStats_0.62.0     timechange_0.2.0       spatstat.sparse_3.0-2  colorspace_2.0-3      
# [25] rvest_1.0.1            ggrepel_0.9.1          textshaping_0.3.6      haven_2.4.3           
# [29] crayon_1.5.1           jsonlite_1.8.0         progressr_0.10.0       spatstat.data_3.0-1   
# [33] survival_3.2-12        zoo_1.8-10             glue_1.6.2             polyclip_1.10-0       
# [37] gtable_0.3.0           leiden_0.3.9           future.apply_1.8.1     abind_1.4-5           
# [41] scales_1.2.0           DBI_1.1.2              spatstat.random_3.1-5  miniUI_0.1.1.1        
# [45] Rcpp_1.0.8.3           viridisLite_0.4.2      xtable_1.8-4           reticulate_1.25       
# [49] clisymbols_1.2.0       htmlwidgets_1.6.2      httr_1.4.3             RColorBrewer_1.1-3    
# [53] ellipsis_0.3.2         ica_1.0-2              pkgconfig_2.0.3        farver_2.1.0          
# [57] uwot_0.1.16            dbplyr_2.1.1           deldir_1.0-6           utf8_1.2.2            
# [61] tidyselect_1.1.2       labeling_0.4.2         rlang_1.1.0            reshape2_1.4.4        
# [65] later_1.3.0            munsell_0.5.0          cellranger_1.1.0       tools_4.1.3           
# [69] cli_3.6.1              generics_0.1.2         broom_0.7.9            ggridges_0.5.3        
# [73] fastmap_1.1.0          goftest_1.2-2          fs_1.5.2               fitdistrplus_1.1-5    
# [77] RANN_2.6.1             pbapply_1.5-0          future_1.25.0          nlme_3.1-152          
# [81] mime_0.12              xml2_1.3.2             compiler_4.1.3         rstudioapi_0.13       
# [85] plotly_4.10.0          png_0.1-8              spatstat.utils_3.1-1   reprex_2.0.1          
# [89] stringi_1.7.6          lattice_0.20-44        Matrix_1.6-4           vctrs_0.6.1           
# [93] pillar_1.7.0           lifecycle_1.0.3        spatstat.geom_3.2-4    lmtest_0.9-40         
# [97] RcppAnnoy_0.0.22       data.table_1.14.2      cowplot_1.1.1          irlba_2.3.5.1         
# [101] httpuv_1.6.5           patchwork_1.1.1        R6_2.5.1               promises_1.2.0.1      
# [105] KernSmooth_2.23-20     gridExtra_2.3          parallelly_1.31.1      codetools_0.2-18      
# [109] MASS_7.3-54            assertthat_0.2.1       withr_2.5.0            sctransform_0.3.5     
# [113] parallel_4.1.3         hms_1.1.0              grid_4.1.3             Rtsne_0.16            
# [117] spatstat.explore_3.2-1 shiny_1.7.1            lubridate_1.9.2       


# end ---------------------------------------------------------------------


