# Takahide Nejo, MD, PhD
# Okada Lab, Department of Neurological Surgery, UC San Francisco


# note ---------------------------------------------------------------------

# Nejo T, et al., Glioma-neuronal circuit remodeling induces regional immunosuppression
# Analysis of SB28 tumor-bearing mouse brain 10X Visium st-RNA-seq data

# Below is the script for the first sample. The data from the second sample was processed and analyzed similarly. 

# Using different versions of packages, such as SPATA2, may lead to varying outcomes, which we have not thoroughly validated.

# The RDS object "visium_SB28_1_subset_infiltration_area_edited.rds" is used as input. For more details, please refer to the previous step, "2_04_stRNAseq_correlation.R".


# rm all ------------------------------------------------------------------

rm(list = ls(all.names = T))


# packages ----------------------------------------------------------------

library(SPATA2) 
library(tidyverse)
library(tidylog)
library(EBImage)
library(hdf5r)
library(ggsci)
library(ggpubr)
library(xCell) # devtools::install_github('dviraran/xCell')
library(gplots) 
library(pheatmap)
library(ggplotify)


# check -------------------------------------------------------------------

packageVersion("SPATA2")
# [1] ‘2.0.4’


# dir ---------------------------------------------------------------------

dir.1 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/2_04_stRNAseq_correlation/"
dir.2 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/2_05_stRNAseq_xCell_deconvolution/"
if(dir.exists(dir.2) == F){dir.create(dir.2)}


# 1. load the spata-object ------------------------------------------------

START.TIME <- Sys.time()

setwd(dir.1)
in.f <- "visium_SB28_1_subset_infiltration_area_edited.rds"
spata.obj <- readRDS(in.f)

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 2.417139 secs


# check -------------------------------------------------------------------

spata.obj 
# An object of class 'spata2' that contains 1 sample named 'A1_166L_Cas9'.

spata.obj %>% 
  getFeatureDf() %>% 
  colnames() %>% 
  print()
# [1] "barcodes"            "sample"              "orig.ident"         
# [4] "nCount_Spatial"      "nFeature_Spatial"    "percent.mt"         
# [7] "percent.RB"          "Spatial_snn_res.0.8" "seurat_clusters"    
# [10] "histology"           "RCTM_NEU"            "VRHK_MES"           
# [13] "segmentation"        "VRHK_MES_cl"         "RCTM_NEU_cl"        
# [16] "label"               "label.2"             "MF_P_SYN"           
# [19] "HM_TNFA"             "HM_IFNG"             "HM_INFLAM" 

spata.obj %>% 
  getExpressionMatrixNames() %>% 
  print()
# [1] "scaled"   "denoised"

spata.obj %>% 
  getActiveMatrixName() %>% 
  print()
# [1] "denoised"

spata.obj %>% 
  getExpressionMatrix() %>% 
  dim() %>% 
  print()
# [1] 19223   148


# prep --------------------------------------------------------------------

mat.1 <- spata.obj %>% 
  getExpressionMatrix() 


# check -------------------------------------------------------------------

mat.1 %>% class() %>% print()
# [1] "matrix" "array" 

mat.1[1:3, 1:3]
# AAATGGTCAATGTGCC-1 AACACGACTGTACTGA-1 AAGGAGCGGTTGGTGC-1
# Xkr4          0.18813181        -0.05182727        -0.21926028
# Rp1           0.03685591        -0.09563567        -0.08453160
# Sox17        -0.16101807        -0.35938329        -0.08148826

rownames(mat.1) <- toupper(rownames(mat.1)) # convert gene symbol (mouse) to uppercase


# run xCell ---------------------------------------------------------------

START.TIME <- Sys.time()

list.xcell.res <- list()

xcell.res <- mat.1 %>% 
  xCellAnalysis()

# [1] "Num. of genes: 9427"
# Setting parallel calculations through a MulticoreParam back-end
# with workers=4 and tasks=100.
# Estimating ssGSEA scores for 489 gene sets.
# |========================================================================| 100%

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 14.31452 secs


# check -------------------------------------------------------------------

xcell.res %>% dim() %>% print()
# [1]  67 148

xcell.res[1:5, 1:3] %>% print()
# AAATGGTCAATGTGCC-1 AACACGACTGTACTGA-1 AAGGAGCGGTTGGTGC-1
# aDC               0.001492873       4.093682e-01         0.27793106
# Adipocytes        0.000000000       5.769047e-19         0.02523801
# Astrocytes        0.000000000       4.834579e-17         0.00000000
# B-cells           0.031447913       1.749415e-01         0.00000000
# Basophils         0.000000000       2.794330e-01         0.33662798

xcell.res %>% 
  rownames()
{
  # [1] "aDC"                           "Adipocytes"                   
  # [3] "Astrocytes"                    "B-cells"                      
  # [5] "Basophils"                     "CD4+ memory T-cells"          
  # [7] "CD4+ naive T-cells"            "CD4+ T-cells"                 
  # [9] "CD4+ Tcm"                      "CD4+ Tem"                     
  # [11] "CD8+ naive T-cells"            "CD8+ T-cells"                 
  # [13] "CD8+ Tcm"                      "CD8+ Tem"                     
  # [15] "cDC"                           "Chondrocytes"                 
  # [17] "Class-switched memory B-cells" "CLP"                          
  # [19] "CMP"                           "DC"                           
  # [21] "Endothelial cells"             "Eosinophils"                  
  # [23] "Epithelial cells"              "Erythrocytes"                 
  # [25] "Fibroblasts"                   "GMP"                          
  # [27] "Hepatocytes"                   "HSC"                          
  # [29] "iDC"                           "Keratinocytes"                
  # [31] "ly Endothelial cells"          "Macrophages"                  
  # [33] "Macrophages M1"                "Macrophages M2"               
  # [35] "Mast cells"                    "Megakaryocytes"               
  # [37] "Melanocytes"                   "Memory B-cells"               
  # [39] "MEP"                           "Mesangial cells"              
  # [41] "Monocytes"                     "MPP"                          
  # [43] "MSC"                           "mv Endothelial cells"         
  # [45] "Myocytes"                      "naive B-cells"                
  # [47] "Neurons"                       "Neutrophils"                  
  # [49] "NK cells"                      "NKT"                          
  # [51] "Osteoblast"                    "pDC"                          
  # [53] "Pericytes"                     "Plasma cells"                 
  # [55] "Platelets"                     "Preadipocytes"                
  # [57] "pro B-cells"                   "Sebocytes"                    
  # [59] "Skeletal muscle"               "Smooth muscle"                
  # [61] "Tgd cells"                     "Th1 cells"                    
  # [63] "Th2 cells"                     "Tregs"                        
  # [65] "ImmuneScore"                   "StromaScore"                  
  # [67] "MicroenvironmentScore"    
  }


# edit --------------------------------------------------------------------

xcell.res.2 <- xcell.res %>% 
  as.data.frame() %>% 
  rownames_to_column("cell.type") %>% 
  as_tibble()



# check -------------------------------------------------------------------

xcell.res.2 %>% dim() %>% print()
# [1]  67 149


# save --------------------------------------------------------------------

setwd(dir.2)
out.f <- "001_xCell_out.tsv"
write_tsv(xcell.res.2, out.f)


# curate the cell types to display --------------------------------

immune.cells <- c("CD4+ T-cells", 
                  "CD8+ T-cells", 
                  "Tregs", 
                  "Th1 cells", "Th2 cells", #"Tgd cells", 
                  "NK cells", "NKT", 
                  "B-cells", "Plasma cells", 
                  "Monocytes", "Macrophages", "Macrophages M1", "Macrophages M2", 
                  "DC", 
                  "Neutrophils", "Eosinophils", "Mast cells", "Basophils")

additional.cell.types <- c("Neurons", "Endothelial cells", "Astrocytes", "Pericytes")

immune.cells %>% length() 
# [1] 18

additional.cell.types %>% length()
# [1] 4



# prep cor for correlation matrix ---------------------------------------------------------------------

# ref: Marquez-Galera A, et al., 2022 STAR Protocols. (https://www.sciencedirect.com/science/article/pii/S2666166722000016)

{
  set.seed(1234)
  
  # a. Define the correlation plot palette. 
  
  # i. Create a matrix of 50x10 random values within the range [−1, +1].
  
  random.matrix <- matrix(runif(500, min = -1, max = 1), nrow = 50)
  
  
  # ii. Produce the sample quantiles corresponding to the given probabilities.
  
  quantile.range <- quantile(random.matrix, probs = seq(0, 1, 0.01))
  
  
  # iii. Define the quantiles where the minimum and maximum correlation values were set to the lowest and highest color. (In our case, it was set empirically to 35% and 83% respectively, as these values maximized the contrast for the distribution of values.)
  
  # palette.breaks <- seq(quantile.range["35%"], quantile.range["83%"], 0.06)
  palette.breaks <- seq(from = quantile.range["0%"], to = quantile.range["100%"], by = 0.02)
  
  
  # iv. Create a color ramp that maps the previous interval.
  
  color.palette <- colorRampPalette(c("#0571b0", "#f7f7f7", "#ca0020"))(length(palette.breaks)-1)
  
  
  # c. Define the hierarchical cluster analysis function, with Pearson correlation coefficient as distance matrix and average as agglomeration method.

  clustFunction <- function(x)
    hclust(as.dist(1-cor(t(as.matrix(x)), method = "pearson")), method = "complete")
  
  # ref: https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/hclust
  
  
  # d. Define the heat map correlation function, with Pearson correlation coefficient as numeric matrix and the color palette, breaks and clustering function set above.
  
  heatmapPearson <- function(correlations)
    heatmap.2(x = correlations,
              col = color.palette,
              breaks = palette.breaks,
              trace = "none", symm = T,
              hclustfun = clustFunction, 
              margins = c(8, 8)
    )
}


# draw the correlation matrix heatmap ------------------------------------------------

fdf <- spata.obj %>%
  getFeatureDf() %>%
  as_tibble() %>%
  dplyr::select(barcodes, RCTM_NEU, MF_P_SYN, HM_TNFA, HM_IFNG, HM_INFLAM)

mat <- xcell.res.2 %>% 
  dplyr::filter(is.element(cell.type, c(immune.cells, additional.cell.types)))

mat <- mat %>%
  gather(barcodes, score, -1) %>%
  spread(key = cell.type, value = score) %>%
  inner_join(fdf, by = "barcodes") %>%
  gather(cell.type, score, -1) %>%
  spread(key = barcodes, value = score)

tmp.rownames <- mat$cell.type
mat <- mat[, -1] %>% 
  as.matrix()
rownames(mat) <- tmp.rownames


# check

mat[1:3, 1:3]
# AAATGGTCAATGTGCC-1 AACACGACTGTACTGA-1 AAGGAGCGGTTGGTGC-1
# Astrocytes         0.00000000       4.834579e-17           0.000000
# B-cells            0.03144791       1.749415e-01           0.000000
# Basophils          0.00000000       2.794330e-01           0.336628


# correlation matrix heatmap

cor.res <- cor(method = "pearson", t(mat))

setwd(dir.2)
out.f <- "002_xCell_res_heatmap_inter_cell_type_correlations.pdf"
pdf(out.f, w = 18, h = 16)
heatmapPearson(cor.res)
dev.off()


# si ----------------------------------------------------------------------

Sys.time()
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
# [1] ggplotify_0.1.2 gplots_3.2.0    pheatmap_1.0.12 xCell_1.1.0    
# [5] ggpubr_0.4.0    ggsci_3.2.0     hdf5r_1.3.9     EBImage_4.34.0 
# [9] tidylog_1.0.2   forcats_0.5.1   stringr_1.5.1   dplyr_1.1.4    
# [13] purrr_1.0.2     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1   
# [17] ggplot2_3.5.1   tidyverse_1.3.1 SPATA2_2.0.4   
# 
# loaded via a namespace (and not attached):
# [1] scattermore_1.2             SPATAData_0.0.0.9000       
# [3] SeuratObject_4.1.3          ragg_1.3.3                 
# [5] bit64_4.5.2                 irlba_2.3.5.1              
# [7] DelayedArray_0.18.0         data.table_1.16.4          
# [9] KEGGREST_1.32.0             RCurl_1.98-1.3             
# [11] generics_0.1.3              BiocGenerics_0.38.0        
# [13] ScaledMatrix_1.0.0          cowplot_1.1.1              
# [15] RSQLite_2.2.8               RANN_2.6.1                 
# [17] future_1.25.0               bit_4.5.0.1                
# [19] tzdb_0.4.0                  spatstat.data_3.0-1        
# [21] xml2_1.3.2                  lubridate_1.9.4            
# [23] httpuv_1.6.15               SummarizedExperiment_1.22.0
# [25] assertthat_0.2.1            confuns_1.0.3              
# [27] hms_1.1.3                   promises_1.3.2             
# [29] caTools_1.18.3              dbplyr_2.1.1               
# [31] readxl_1.3.1                igraph_1.3.1               
# [33] DBI_1.2.3                   htmlwidgets_1.6.4          
# [35] spatstat.geom_3.2-4         stats4_4.1.3               
# [37] backports_1.5.0             annotate_1.70.0            
# [39] deldir_1.0-6                sparseMatrixStats_1.4.2    
# [41] MatrixGenerics_1.4.3        vctrs_0.6.5                
# [43] SingleCellExperiment_1.14.1 Biobase_2.52.0             
# [45] ROCR_1.0-11                 abind_1.4-5                
# [47] cachem_1.1.0                withr_3.0.2                
# [49] progressr_0.10.0            vroom_1.6.5                
# [51] sctransform_0.3.5           goftest_1.2-2              
# [53] cluster_2.1.2               lazyeval_0.2.2             
# [55] crayon_1.5.3                spatstat.explore_3.2-1     
# [57] pkgconfig_2.0.3             labeling_0.4.3             
# [59] units_0.8-5                 GenomeInfoDb_1.28.1        
# [61] nlme_3.1-152                rlang_1.1.4                
# [63] globals_0.15.0              lifecycle_1.0.4            
# [65] miniUI_0.1.1.1              modelr_0.1.8               
# [67] rsvd_1.0.5                  cellranger_1.1.0           
# [69] polyclip_1.10-7             GSVA_1.40.1                
# [71] matrixStats_0.62.0          lmtest_0.9-40              
# [73] graph_1.70.0                tiff_0.1-11                
# [75] Matrix_1.6-4                carData_3.0-4              
# [77] Rhdf5lib_1.14.2             zoo_1.8-10                 
# [79] reprex_2.0.1                ggridges_0.5.6             
# [81] png_0.1-8                   viridisLite_0.4.2          
# [83] clisymbols_1.2.0            bitops_1.0-9               
# [85] KernSmooth_2.23-20          rhdf5filters_1.4.0         
# [87] Biostrings_2.60.2           anndata_0.7.5.6            
# [89] blob_1.2.2                  DelayedMatrixStats_1.14.3  
# [91] parallelly_1.31.1           spatstat.random_3.1-5      
# [93] jpeg_0.1-9                  rstatix_0.7.0              
# [95] gridGraphics_0.5-1          S4Vectors_0.30.2           
# [97] ggsignif_0.6.2              beachmat_2.8.1             
# [99] scales_1.3.0                memoise_2.0.1              
# [101] GSEABase_1.54.0             magrittr_2.0.3             
# [103] plyr_1.8.7                  ica_1.0-2                  
# [105] zlibbioc_1.38.0             compiler_4.1.3             
# [107] RColorBrewer_1.1-3          fitdistrplus_1.1-5         
# [109] cli_3.6.3                   XVector_0.32.0             
# [111] listenv_0.8.0               patchwork_1.3.0            
# [113] pbapply_1.5-0               MASS_7.3-54                
# [115] tidyselect_1.2.1            stringi_1.8.4              
# [117] textshaping_0.3.6           BiocSingular_1.8.1         
# [119] locfit_1.5-9.4              ggrepel_0.9.1              
# [121] grid_4.1.3                  tools_4.1.3                
# [123] timechange_0.3.0            future.apply_1.8.1         
# [125] parallel_4.1.3              rio_0.5.27                 
# [127] rstudioapi_0.17.1           foreign_0.8-81             
# [129] gridExtra_2.3               farver_2.1.2               
# [131] Rtsne_0.16                  digest_0.6.37              
# [133] shiny_1.10.0                pracma_2.4.4               
# [135] quadprog_1.5-8              Rcpp_1.0.13-1              
# [137] GenomicRanges_1.44.0        car_3.0-11                 
# [139] broom_1.0.7                 later_1.4.1                
# [141] RcppAnnoy_0.0.22            httr_1.4.7                 
# [143] AnnotationDbi_1.54.1        colorspace_2.1-1           
# [145] rvest_1.0.1                 XML_3.99-0.7               
# [147] fs_1.6.5                    tensor_1.5                 
# [149] reticulate_1.40.0           IRanges_2.26.0             
# [151] splines_4.1.3               uwot_0.1.16                
# [153] yulab.utils_0.1.8           spatstat.utils_3.1-1       
# [155] sp_2.1-4                    plotly_4.10.4              
# [157] systemfonts_1.1.0           xtable_1.8-4               
# [159] jsonlite_1.8.9              R6_2.5.1                   
# [161] pillar_1.10.0               htmltools_0.5.8.1          
# [163] mime_0.12                   glue_1.8.0                 
# [165] fastmap_1.2.0               BiocParallel_1.26.2        
# [167] fftwtools_0.9-11            codetools_0.2-18           
# [169] utf8_1.2.4                  lattice_0.20-44            
# [171] spatstat.sparse_3.0-2       curl_6.0.1                 
# [173] leiden_0.3.9                gtools_3.9.5               
# [175] zip_2.2.0                   openxlsx_4.2.4             
# [177] survival_3.2-12             munsell_0.5.1              
# [179] rhdf5_2.36.0                GenomeInfoDbData_1.2.6     
# [181] HDF5Array_1.20.0            haven_2.4.3                
# [183] reshape2_1.4.4              gtable_0.3.6               
# [185] Seurat_4.3.0.1             


# end ---------------------------------------------------------------------
