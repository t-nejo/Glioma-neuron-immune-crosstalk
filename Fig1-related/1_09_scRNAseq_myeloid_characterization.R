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
library(wesanderson)


# check -------------------------------------------------------------------

packageVersion("Seurat")
# [1] ‘4.3.0.1’
packageVersion("sctransform")
# [1] ‘0.3.5’


# dir ---------------------------------------------------------------------

dir.1 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/1_08_scRNAseq_myeloid_subset/"
dir.2 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/1_09_scRNAseq_myeloid_characterization/"
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


# Feature plot (Fig. 1d) -------------------------------------------------------------------

DefaultAssay(seurat.obj) <- "RNA"

markers.to.show <- c(
  "CX3CR1", "P2RY12", 
  sort(c("IL1B", "CCL3", "CD86", "TNF", "TNFAIP3", "NR4A2", "CD83", "IER2")), 
  "CD68", "CD163", 
  sort(c("RNASE1", "S100A10", "VIM", "LGALS1", "LGALS3", "LDHA", "LYZ", "TGFBI"))
)

CUTOFF <- c(
  3, 3, NA, 6, 3, 
  3, NA, 3, 3, NA, 
  NA, NA, NA, NA, NA, 
  NA, NA, NA, 6, NA
)

g <- seurat.obj %>% 
  FeaturePlot(
    features = markers.to.show,
    pt.size = 0.1,
    slot = "scale.data",
    max.cutoff = CUTOFF,
    cols = c("lightgrey", wes_palette("Zissou1", 100, type = "continuous")),
    ncol = 5
  ) & theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "#f0f0f0"),
    axis.text = element_blank(),
    axis.title = element_blank(), 
    plot.title = element_text(size = 24 / 2.845276, face = "italic"), 
  )
plot(g)

setwd(dir.2)
out.f <- "001_ftrplot_myeloid_markers.pdf"
ggsave(out.f, g, w = 12, h = 7.5)


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


# DimPlot (Fig. 1e) -------------------------------------------------------

plots <- seurat.obj %>% 
  DimPlot(
  group.by = "seurat_clusters.new", 
  combine = FALSE,
  cols = c("#E69F00", "#56B4E9", "#999999"), 
  pt.size = .8
)
plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "bottom")
                + guides(color = guide_legend(ncol = 3, byrow = FALSE, override.aes = list(size = 4))))
g <- CombinePlots(plots, ncol = 1)
g <- g & theme(
  panel.grid.major = element_line(color = "#f0f0f0"),
  axis.text = element_blank()
)
g <- g + labs(title = NULL, x = "UMAP-1", y = "UMAP-2")
plot(g)

setwd(dir.2)
out.f <- "002_dimplot_myeloid_01_inflammatory_status.pdf"
ggsave(out.f, g, w = 6, h = 5)



# DimPlot (Fig. 1f) -------------------------------------------------------

# HFC v. LFC

plots <- seurat.obj %>% 
  DimPlot(
    group.by = "seurat_clusters.new", 
    split.by = "tissue",
    combine = FALSE,
    cols = c("#E69F00", "#56B4E9", "#999999"),
    pt.size = .8
  )
plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "none")
                + guides(color = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(size = 4))))
g <- CombinePlots(plots, ncol = 1)
g <- g & theme(
  panel.grid.major = element_line(color = "#f0f0f0"),
  axis.text = element_blank(),
  axis.title = element_blank()
)
g <- g + labs(title = NULL, x = "UMAP-1", y = "UMAP-2")
g

setwd(dir.2)
out.f <- "003_dimplot_myeloid_02_inflammatory_status_hfc_lfc.pdf"
ggsave(out.f, g, w = 10, h = 5)


# Line plot (Fig. 1f) ---------------------------------------------------------------

n.table.1 <- seurat.obj@meta.data %>% 
  dplyr::select(tissue, seurat_clusters.new) %>%
  table() %>%
  as_tibble() %>% 
  arrange(tissue) %>%
  mutate(total = ifelse(tissue == "hfc", 851, 2924)) %>%
  mutate(fraction = n/total)

n.table.1 %>% print()
# # A tibble: 6 × 5
# tissue seurat_clusters.new     n total fraction
# <chr>  <chr>               <int> <dbl>    <dbl>
# 1 hfc    pro-inflammatory      352   851    0.414
# 2 hfc    anti-inflammatory     321   851    0.377
# 3 hfc    undetermined          178   851    0.209
# 4 lfc    pro-inflammatory     1554  2924    0.531
# 5 lfc    anti-inflammatory     647  2924    0.221
# 6 lfc    undetermined          723  2924    0.247


n.table.1 <- n.table.1 %>%
  arrange(tissue) %>%
  mutate(total = ifelse(tissue == "hfc", 851, 2924)) %>%
  mutate(pct = n/total) %>% 
  dplyr::filter(seurat_clusters.new != "undetermined")

n.table.1$seurat_clusters.new <- factor(n.table.1$seurat_clusters.new, 
                                        levels = c("pro-inflammatory", "anti-inflammatory"))


g <- ggplot(n.table.1, aes(x = tissue, y = pct, color = seurat_clusters.new)) ;
g <- g + geom_point(size = 1.5)
g <- g + geom_line(aes(group = seurat_clusters.new), size = 1)
g <- g + theme_classic()
# g <- g + scale_color_manual(values = pal_npg()(3))
g <- g + scale_color_manual(values = c("#E69F00", "#56B4E9"))
g <- g + ylim(0, 0.6)
g <- g & theme(
  panel.grid.major = element_line(color = "#f0f0f0"),
  legend.position = "right"
)
g <- g + guides(color = element_text(size = 8))
g <- g + labs(title = NULL, x = NULL, y = "pct", color = NULL) ;
g

setwd(dir.2)
out.f <- "004_line_plot_myeloid_inflammatory_status_hfc_lfc.pdf"
ggsave(out.f, g, w = 6, h = 4)



# Vln plot (Extended Data Fig. S4a) ---------------------------------------

# 20230626_revisit_myeloid_054b_vlnPlot_20candidates_size_adjusted.pdf

g <- seurat.obj %>% 
  subset(subset = seurat_clusters.new != "undetermined") %>% 
  VlnPlot(
    features = markers.to.show, 
    slot = "scale.data",
    assay = "RNA", 
    ncol = 5, 
    cols = c("#E69F00", "#56B4E9"), 
    pt.size = 0
  ) 
g <- g & theme(
  axis.text.x = element_text(size = 3), 
  axis.text.y = element_text(size = 6), 
  axis.title = element_blank(), 
  plot.title = element_text(size = 10, face = "italic"), 
)
g

setwd(dir.2)
out.f <- "005_vln_plot_myeloid_inflammatory_status_hfc_lfc.pdf"
ggsave(out.f, g, width = 7, height = 4.8) 



# Bar plot (Extended Data Fig. S4b)----------------------------------------------------------------

metadata.1 <- seurat.obj@meta.data %>% 
  dplyr::select(tissue, seurat_clusters.new) %>% 
  table() %>% 
  as_tibble() %>% 
  mutate(tissue = ifelse(grepl("hfc", tissue), "HFC", "LFC")) 

metadata.1$tissue <- factor(metadata.1$tissue, levels = c("HFC", "LFC"))
metadata.1$seurat_clusters.new <- 
  factor(metadata.1$seurat_clusters.new, 
         levels = rev(c("pro-inflammatory", "anti-inflammatory", "undetermined")))

COLORS <- rev(c("#E69F00", "#56B4E9", "#999999"))


# ggplot ------------------------------------------------------------------

g <- ggplot(metadata.1, aes(x = tissue, y = n, fill = seurat_clusters.new))
g <- g + geom_bar(position = "fill", stat = "identity", color = "black")
g <- g + scale_fill_manual(values = COLORS)
g <- g + labs(title = NULL, x = NULL, y = "cell count", fill = NULL, shape = NULL)
g <- g + theme_classic()
g <- g + theme(
  panel.grid.major = element_line(color = "#f0f0f0"), 
  plot.title = element_text(hjust = 0.5, face = "bold"), 
  strip.background = element_blank(), 
  strip.text = element_text(size = 12), 
  legend.position = "right", 
  # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
)
plot(g)

setwd(dir.2)
out.f <- "006_bar_plot_myeloid_inflammatory_status_hfc_lfc.pdf"
ggsave(out.f, g, width = 6, height = 4.5) 


# Bar plot (Extended Data Fig. S4c)----------------------------------------------------------------

# per patient

metadata.2 <- seurat.obj@meta.data %>% 
  as_tibble() %>% 
  dplyr::select(orig.ident, seurat_clusters.new) %>% 
  table() %>% 
  as_tibble() %>% 
  mutate(tissue = ifelse(grepl("hfc", orig.ident), "HFC", "LFC")) %>% 
  mutate(pt = ifelse(grepl(1, orig.ident), "pt-1", 
                     ifelse(grepl(2, orig.ident), "pt-2", "pt-3")))

metadata.2$tissue <- factor(metadata.2$tissue, levels = c("HFC", "LFC"))
metadata.2$seurat_clusters.new <- 
  factor(metadata.2$seurat_clusters.new, 
         levels = rev(c("pro-inflammatory", "anti-inflammatory", "undetermined")))

COLORS <- rev(c("#E69F00", "#56B4E9", "#999999"))

g <- ggplot(metadata.2, aes(x = tissue, y = n, fill = seurat_clusters.new))
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
  legend.position = "right", 
  # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
)
plot(g)

setwd(dir.2)
out.f <- "007_bar_plot_myeloid_inflammatory_status_hfc_lfc_per_patient.pdf"
ggsave(out.f, g, width = 8, height = 4.5) 



# Line plot (Extended Data Fig. S4d) --------------------------------------------------

# per patient

sum.2 <- metadata.2 %>% 
  group_by(orig.ident) %>% 
  summarise(total = sum(n)) %>% 
  ungroup()

metadata.2b <- metadata.2 %>% 
  left_join(sum.2, by = "orig.ident") %>% 
  mutate(pct = n / total * 100) %>% 
  arrange(pt, tissue, desc(seurat_clusters.new)) %>% 
  dplyr::filter(seurat_clusters.new != "undetermined") %>% 
  dplyr::select(pt, tissue, seurat_clusters.new, n, total, pct)


COLORS <- rev(c("#E69F00", "#56B4E9"))

g <- ggplot(metadata.2b, aes(x = tissue, y = pct, color = seurat_clusters.new)) ;
g <- g + geom_point(size = 1)
g <- g + geom_line(aes(group = seurat_clusters.new), size = 1)
g <- g + theme_classic()
# g <- g + scale_color_manual(values = pal_npg()(3))
g <- g + scale_color_manual(values = COLORS)
g <- g + ylim(0, 60)
g <- g + facet_wrap(~ pt)
g <- g & theme(
  panel.grid.major = element_line(color = "#f0f0f0"),
  legend.position = "none", # "right",   
  strip.background = element_blank(), 
  strip.text = element_text(size = 12), 
)
g <- g + labs(title = NULL, x = NULL, y = "pct", color = NULL) ;
g

setwd(dir.2)
out.f <- "008_line_plot_myeloid_inflammatory_status_hfc_lfc_per_patient.pdf"
ggsave(out.f, g, width = 8, height = 4.5) 


# si ----------------------------------------------------------------------

Sys.time() %>% print()
# [1] "2024-12-19 23:30:46 PST"

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
# [1] wesanderson_0.3.6  SeuratObject_4.1.3 Seurat_4.3.0.1     tidylog_1.0.2      forcats_0.5.1     
# [6] stringr_1.4.0      dplyr_1.0.9        purrr_0.3.4        readr_2.1.2        tidyr_1.2.0       
# [11] tibble_3.1.7       ggplot2_3.3.6      tidyverse_1.3.1   
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
# [81] mime_0.12              ggrastr_1.0.2          xml2_1.3.2             compiler_4.1.3        
# [85] rstudioapi_0.13        beeswarm_0.4.0         plotly_4.10.0          png_0.1-8             
# [89] spatstat.utils_3.1-1   reprex_2.0.1           stringi_1.7.6          lattice_0.20-44       
# [93] Matrix_1.6-4           vctrs_0.6.1            pillar_1.7.0           lifecycle_1.0.3       
# [97] spatstat.geom_3.2-4    lmtest_0.9-40          RcppAnnoy_0.0.22       data.table_1.14.2     
# [101] cowplot_1.1.1          irlba_2.3.5.1          httpuv_1.6.5           patchwork_1.1.1       
# [105] R6_2.5.1               promises_1.2.0.1       KernSmooth_2.23-20     gridExtra_2.3         
# [109] vipor_0.4.5            parallelly_1.31.1      codetools_0.2-18       MASS_7.3-54           
# [113] assertthat_0.2.1       withr_2.5.0            sctransform_0.3.5      parallel_4.1.3        
# [117] hms_1.1.0              grid_4.1.3             Rtsne_0.16             spatstat.explore_3.2-1
# [121] shiny_1.7.1            lubridate_1.9.2        ggbeeswarm_0.6.0      


# end ---------------------------------------------------------------------


