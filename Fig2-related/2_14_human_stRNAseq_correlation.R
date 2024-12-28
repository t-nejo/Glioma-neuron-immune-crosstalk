# Takahide Nejo, MD, PhD
# Okada Lab, Department of Neurological Surgery, UC San Francisco


# note ---------------------------------------------------------------------

# Nejo T, et al., Glioma-neuronal circuit remodeling induces regional immunosuppression
# Analysis of SPATAData human GBM 10X Visium st-RNA-seq data

# Below is the script for one sample ("260_T"). The data from the other samples were processed and analyzed similarly. 

# Using different versions of packages, such as SPATA2, may lead to varying outcomes, which we have not thoroughly validated.

# The RDS object "visium_SPATAData_260_T_cna_index_added.rds" is used as input. For more details, please refer to the previous step, "2_13_human_stRNAseq_cnv_index.R".

# ref: https://themilolab.github.io/SPATA2/articles/spata-data.html


# rm all ------------------------------------------------------------------

rm(list = ls(all.names = T))


# packages ----------------------------------------------------------------

library(SPATA2)
library(tidyverse)
library(tidylog)
library(ggpubr)
library(msigdbr)


# check -------------------------------------------------------------------

packageVersion("SPATA2")
# [1] ‘2.0.4’


# dir ---------------------------------------------------------------------

dir.1 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/2_13_human_stRNAseq_cnv_index/"
dir.2 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/2_14_human_stRNAseq_correlation"

if(dir.exists(dir.2) == F){dir.create(dir.2)}


# load the spata-object ------------------------------------------------

START.TIME <- Sys.time()

setwd(dir.1)
in.f <- "visium_SPATAData_260_T_cna_index_added.rds"
spata.obj <- readRDS(in.f)

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 10.43524 secs


# check -------------------------------------------------------------------

spata.obj 
# An object of class 'spata2' that contains 1 sample named '260_T'.

spata.obj %>% 
  getCoordsDf() %>% 
  print()

# # A tibble: 2,997 × 6
# barcodes           sample     x     y section outline
# <chr>              <chr>  <dbl> <dbl> <chr>   <lgl>  
# 1 AAACAATCTACTAGCA-1 260_T  231.   532. 1       TRUE   
# 2 AAACAGAGCGACTCCT-1 260_T  420.   459. 1       FALSE  
# 3 AAACAGCTTTCAGAAG-1 260_T  101.   274. 1       FALSE  
# 4 AAACAGGGTCTATATT-1 260_T  116.   248. 1       FALSE  
# 5 AAACATGGTGAGAGGA-1 260_T   66.7  151. 1       TRUE   
# 6 AAACCGGGTAGGTACC-1 260_T  172.   280. 1       FALSE  
# 7 AAACCGTTCGTCCAGG-1 260_T  224.   214. 1       FALSE  
# 8 AAACCTCATGAAGTTG-1 260_T  139.   312. 1       FALSE  
# 9 AAACGAAGAACATACC-1 260_T  309.   512. 1       FALSE  
# 10 AAACGAAGATGGAGTA-1 260_T   81.8  177. 1       FALSE  
# # ℹ 2,987 more rows
# # ℹ Use `print(n = ...)` to see more rows

# all expression matrices before denoising

spata.obj %>% 
  getExpressionMatrixNames() %>% 
  print()
# [1] "scaled"   "denoised"

# active expression matrix before denoising
spata.obj %>% 
  getActiveMatrixName() %>% 
  print()
# [1] "denoised"

spata.obj %>% 
  getFeatureDf() %>% 
  colnames() %>% 
  print()
# [1] "barcodes"         "sample"           "segmentation"     "nCount_Spatial"   "nFeature_Spatial"
# [6] "Chr0"             "Chr1"             "Chr10"            "Chr11"            "Chr12"           
# [11] "Chr13"            "Chr14"            "Chr15"            "Chr16"            "Chr17"           
# [16] "Chr18"            "Chr19"            "Chr2"             "Chr20"            "Chr21"           
# [21] "Chr22"            "Chr23"            "Chr3"             "Chr4"             "Chr5"            
# [26] "Chr6"             "Chr7"             "Chr8"             "Chr9"             "Chr10p"          
# [31] "Chr10q"           "Chr11p"           "Chr11q"           "Chr12p"           "Chr12q"          
# [36] "Chr13q"           "Chr14q"           "Chr15q"           "Chr16p"           "Chr16q"          
# [41] "Chr17p"           "Chr17q"           "Chr18p"           "Chr18q"           "Chr19p"          
# [46] "Chr19q"           "Chr1p"            "Chr1q"            "Chr20p"           "Chr20q"          
# [51] "Chr21q"           "Chr22q"           "Chr2p"            "Chr2q"            "Chr3p"           
# [56] "Chr3q"            "Chr4p"            "Chr4q"            "Chr5p"            "Chr5q"           
# [61] "Chr6p"            "Chr6q"            "Chr7p"            "Chr7q"            "Chr8p"           
# [66] "Chr8q"            "Chr9p"            "Chr9q"            "CNA.index"     


# prep --------------------------------------------------------------------

names.gs.of.interest <- c(
  "GOMF_POSTSYNAPTIC_NEUROTRANSMITTER_RECEPTOR_ACTIVITY", 
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE"
)


# prep MSigDB ------------------------------------------------------------------

msigdbr_collections() %>% print(n = Inf)

test_gs <- NULL
for(i in 1:2){
  print(i)
  
  test_gs.i <-  msigdbr(
    species = "Homo sapiens", # "Mus musculus"
    category = c("C5", "H")[i],
    subcategory = c("GO:MF", "")[i]
  ) 
  test_gs <- test_gs %>% 
    bind_rows(test_gs.i)
}


# add new geneset data ---------------------------------------------------

for(i in 1:length(names.gs.of.interest)){
  print(i)
  name.gs.i <- names.gs.of.interest[i]
  
  gs.i <- test_gs %>% 
    dplyr::filter(gs_name == name.gs.i)
  
  genes.i <- gs.i %>% 
    pull(gene_symbol)
  
  spata.obj <- spata.obj %>% 
    addGeneSet(
      class_name = "NEW2", 
      gs_name = name.gs.i, 
      check_genes = F, 
      genes = genes.i, 
      overwrite = T
    ) 
}


# check -------------------------------------------------------------------

spata.obj@used_genesets %>% 
  dplyr::filter(grepl("NEW2_", ont)) %>% 
  distinct(ont) %>% 
  print()
# distinct: removed 802 rows (>99%), 4 rows remaining
# ont
# 1 NEW2_GOMF_POSTSYNAPTIC_NEUROTRANSMITTER_RECEPTOR_ACTIVITY
# 2                     NEW2_HALLMARK_TNFA_SIGNALING_VIA_NFKB
# 3                   NEW2_HALLMARK_INTERFERON_GAMMA_RESPONSE
# 4                       NEW2_HALLMARK_INFLAMMATORY_RESPONSE


# update the featureDf ----------------------------------------------------

coords_df <- getCoordsDf(spata.obj)

# join this data.frame with additional info

joined_df <- spata.obj %>%
  joinWith(
    spata_df = coords_df,
    gene_sets = paste0("NEW2_", names.gs.of.interest), # expression values of the gene set Hallmark-Hypoxia
    verbose = FALSE)


# check / edit -------------------------------------------------------------------

joined_df %>% colnames()
# [1] "barcodes"                                                 
# [2] "sample"                                                   
# [3] "x"                                                        
# [4] "y"                                                        
# [5] "section"                                                  
# [6] "outline"                                                  
# [7] "NEW2_GOMF_POSTSYNAPTIC_NEUROTRANSMITTER_RECEPTOR_ACTIVITY"
# [8] "NEW2_HALLMARK_TNFA_SIGNALING_VIA_NFKB"                    
# [9] "NEW2_HALLMARK_INTERFERON_GAMMA_RESPONSE"                  
# [10] "NEW2_HALLMARK_INFLAMMATORY_RESPONSE"       

colnames(joined_df)[7:10] <- c("MF_P_SYN", "HM_TNFA", "HM_IFNG", "HM_INFLAM")

joined_df %>% colnames()
# [1] "barcodes"  "sample"    "x"         "y"         "section"   "outline"   "MF_P_SYN"  "HM_TNFA"  
# [9] "HM_IFNG"   "HM_INFLAM"

joined_df_simple <- joined_df %>% 
  dplyr::select(barcodes, colnames(joined_df)[7:10])


# add the feature

spata.obj <- spata.obj %>% 
  addFeatures(
    feature_names = colnames(joined_df_simple)[-1],
    feature_df = joined_df_simple, 
    overwrite = T
  )


# check -------------------------------------------------------------------

spata.obj %>% print()
# An object of class 'spata2' that contains 1 sample named 'A1_166L_Cas9'.

spata.obj %>% 
  getFeatureDf() %>% 
  colnames() %>% 
  print()
# [1] "barcodes"         "sample"           "segmentation"     "nCount_Spatial"   "nFeature_Spatial"
# [6] "Chr0"             "Chr1"             "Chr10"            "Chr11"            "Chr12"           
# [11] "Chr13"            "Chr14"            "Chr15"            "Chr16"            "Chr17"           
# [16] "Chr18"            "Chr19"            "Chr2"             "Chr20"            "Chr21"           
# [21] "Chr22"            "Chr23"            "Chr3"             "Chr4"             "Chr5"            
# [26] "Chr6"             "Chr7"             "Chr8"             "Chr9"             "Chr10p"          
# [31] "Chr10q"           "Chr11p"           "Chr11q"           "Chr12p"           "Chr12q"          
# [36] "Chr13q"           "Chr14q"           "Chr15q"           "Chr16p"           "Chr16q"          
# [41] "Chr17p"           "Chr17q"           "Chr18p"           "Chr18q"           "Chr19p"          
# [46] "Chr19q"           "Chr1p"            "Chr1q"            "Chr20p"           "Chr20q"          
# [51] "Chr21q"           "Chr22q"           "Chr2p"            "Chr2q"            "Chr3p"           
# [56] "Chr3q"            "Chr4p"            "Chr4q"            "Chr5p"            "Chr5q"           
# [61] "Chr6p"            "Chr6q"            "Chr7p"            "Chr7q"            "Chr8p"           
# [66] "Chr8q"            "Chr9p"            "Chr9q"            "CNA.index"        "MF_P_SYN"        
# [71] "HM_TNFA"          "HM_IFNG"          "HM_INFLAM"       


# save --------------------------------------------------------------------

START.TIME <- Sys.time() 

setwd(dir.2)
out.f <- "visium_SPATAData_260_T_gs_added.rds"
saveRDS(spata.obj, out.f)

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 1.792389 mins


# visualize (surface plots) (Fig. 2c) ---------------------------------------------------------------

gg.list <- list()
for(i in 1:4){
  print(i)
  name.i <- c("MF_P_SYN", "HM_TNFA", "HM_IFNG", "HM_INFLAM")[i]
  title.i <- c("Postsynaptic Neurotransmitter\nReceptor Activity", 
               "TNFɑ Signaling via-NFκB", 
               "IFNg Response",
               "Inflammatory Response"
               )[i]
  g <- spata.obj %>% 
    plotSurface(
      display_image = T, 
      color_by = name.i, 
      smooth_span = 1,
      pt_size = 2
    ) + ggplot2::labs(title = title.i, color = "score") + 
    ggplot2::theme(
      plot.title = element_text(hjust = 0, size = 12), 
      legend.position = "right"
    )
  
  I <- formatC(i, width = 2, flag = "0")
  setwd(dir.2)
  out.f <- paste0(I, "_surfacePlot_", I, "_", name.i, ".pdf")
  out.f <- "01_surfacePlot_01_MF_P_SYN.pdf"
  ggsave(out.f, g, w = 6, h = 4.5)
  
  gg.list[[i]] <- g + theme(legend.position = "none")
}

g1 <- gg.list[[1]] ; g2 <- gg.list[[2]] ; g3 <- gg.list[[3]]; g4 <- gg.list[[4]]


# visualize (scatter plots - correlation) (Fig. 2i) ---------------------------------------------------------------

gg.list <- list()
for(i in 1:3){
  print(i)
  gs.i <- c("HM_TNFA", "HM_IFNG", "HM_INFLAM")[i]
  
  joined_df_simple.i <- joined_df_simple %>% 
    dplyr::select(barcodes, MF_P_SYN, gs.i)
  
  colnames(joined_df_simple.i)[3] <- "marker"
  
  col.res.i <- cor.test(joined_df_simple.i$MF_P_SYN, joined_df_simple.i$marker, method = "pearson")
  TITLE.i <- paste0("Pearson p = ", format(as.numeric(col.res.i["p.value"]), digits = 3), "\n", "Corr = ", round(as.numeric(col.res.i["estimate"]), 2))
  
  g <- ggplot(joined_df_simple.i, aes(x = marker, y = MF_P_SYN) )
  g <- g + geom_smooth(alpha = 0.3, method = "lm", color = "gray20")
  g <- g + geom_point(size = 1.5, alpha = 0.8)
  g <- g + theme_classic()
  g <- g + labs(title = TITLE.i, x = gs.i, y = "MS_P_SYN")
  g <- g + theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    legend.position = "none",
    panel.grid.major = element_line(color = "#f0f0f0"),
  )
  
  I <- formatC((i + 4), width = 2, flag = "0")
  setwd(dir.2)
  out.f <- paste0(I, "_scatterPlot_cor_", I, "_", gs.i, ".pdf")
  ggsave(out.f, g, w = 6, h = 4.5)
  
  gg.list[[i]] <- g
}

g5 <- gg.list[[1]] ; g6 <- gg.list[[2]] ; g7 <- gg.list[[3]]


g <- ggarrange(
  plotlist = list(NULL, g2, g3, g4, g1, g5, g6, g7), 
  ncol = 4, nrow = 2) + theme(legend.position = "none")

setwd(dir.2)
out.f <- "08_surface_n_scatter_plots_combined.pdf"
ggsave(out.f, g, w = 16, h = 8)


# si ----------------------------------------------------------------------

Sys.time()


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
# [1] msigdbr_7.4.1   ggpubr_0.4.0    tidylog_1.0.2   forcats_0.5.1   stringr_1.5.1   dplyr_1.1.4    
# [7] purrr_1.0.2     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_1.3.1
# [13] SPATA2_2.0.4   
# 
# loaded via a namespace (and not attached):
# [1] utf8_1.2.4                  spatstat.explore_3.2-1      reticulate_1.40.0          
# [4] tidyselect_1.2.1            htmlwidgets_1.6.4           grid_4.1.3                 
# [7] Rtsne_0.16                  munsell_0.5.1               ragg_1.3.3                 
# [10] codetools_0.2-18            ica_1.0-2                   units_0.8-5                
# [13] future_1.25.0               miniUI_0.1.1.1              withr_3.0.2                
# [16] spatstat.random_3.1-5       colorspace_2.1-1            progressr_0.10.0           
# [19] Biobase_2.52.0              rstudioapi_0.17.1           Seurat_4.3.0.1             
# [22] stats4_4.1.3                SingleCellExperiment_1.14.1 ROCR_1.0-11                
# [25] ggsignif_0.6.2              tensor_1.5                  listenv_0.8.0              
# [28] labeling_0.4.3              MatrixGenerics_1.4.3        GenomeInfoDbData_1.2.6     
# [31] polyclip_1.10-7             farver_2.1.2                parallelly_1.31.1          
# [34] vctrs_0.6.5                 generics_0.1.3              timechange_0.3.0           
# [37] R6_2.5.1                    GenomeInfoDb_1.28.1         locfit_1.5-9.4             
# [40] bitops_1.0-9                spatstat.utils_3.1-1        DelayedArray_0.18.0        
# [43] assertthat_0.2.1            promises_1.3.2              scales_1.3.0               
# [46] gtable_0.3.6                globals_0.15.0              goftest_1.2-2              
# [49] rlang_1.1.4                 clisymbols_1.2.0            systemfonts_1.1.0          
# [52] splines_4.1.3               rstatix_0.7.0               lazyeval_0.2.2             
# [55] spatstat.geom_3.2-4         broom_1.0.7                 reshape2_1.4.4             
# [58] abind_1.4-5                 modelr_0.1.8                backports_1.5.0            
# [61] httpuv_1.6.15               tools_4.1.3                 RColorBrewer_1.1-3         
# [64] BiocGenerics_0.38.0         ggridges_0.5.6              Rcpp_1.0.13-1              
# [67] plyr_1.8.7                  zlibbioc_1.38.0             RCurl_1.98-1.3             
# [70] deldir_1.0-6                pbapply_1.5-0               cowplot_1.1.1              
# [73] S4Vectors_0.30.2            zoo_1.8-10                  SeuratObject_4.1.3         
# [76] SummarizedExperiment_1.22.0 haven_2.4.3                 ggrepel_0.9.1              
# [79] cluster_2.1.2               fs_1.6.5                    magrittr_2.0.3             
# [82] magick_2.8.5                data.table_1.16.4           scattermore_1.2            
# [85] openxlsx_4.2.4              lmtest_0.9-40               reprex_2.0.1               
# [88] RANN_2.6.1                  fitdistrplus_1.1-5          anndata_0.7.5.6            
# [91] matrixStats_0.62.0          hms_1.1.3                   patchwork_1.3.0            
# [94] mime_0.12                   fftwtools_0.9-11            xtable_1.8-4               
# [97] rio_0.5.27                  jpeg_0.1-9                  readxl_1.3.1               
# [100] IRanges_2.26.0              gridExtra_2.3               compiler_4.1.3             
# [103] KernSmooth_2.23-20          crayon_1.5.3                SPATAData_0.0.0.9000       
# [106] htmltools_0.5.8.1           mgcv_1.8-36                 later_1.4.1                
# [109] tzdb_0.4.0                  tiff_0.1-11                 lubridate_1.9.4            
# [112] DBI_1.2.3                   dbplyr_2.1.1                MASS_7.3-54                
# [115] babelgene_21.4              Matrix_1.6-4                car_3.0-11                 
# [118] cli_3.6.3                   parallel_4.1.3              igraph_1.3.1               
# [121] GenomicRanges_1.44.0        pkgconfig_2.0.3             foreign_0.8-81             
# [124] sp_2.1-4                    confuns_1.0.3               plotly_4.10.4              
# [127] spatstat.sparse_3.0-2       xml2_1.3.2                  XVector_0.32.0             
# [130] rvest_1.0.1                 digest_0.6.37               sctransform_0.3.5          
# [133] RcppAnnoy_0.0.22            spatstat.data_3.0-1         cellranger_1.1.0           
# [136] leiden_0.3.9                uwot_0.1.16                 curl_6.0.1                 
# [139] shiny_1.10.0                EBImage_4.34.0              lifecycle_1.0.4            
# [142] nlme_3.1-152                jsonlite_1.8.9              carData_3.0-4              
# [145] viridisLite_0.4.2           pillar_1.10.0               lattice_0.20-44            
# [148] fastmap_1.2.0               httr_1.4.7                  survival_3.2-12            
# [151] glue_1.8.0                  zip_2.2.0                   png_0.1-8                  
# [154] stringi_1.8.4               textshaping_0.3.6           irlba_2.3.5.1              
# [157] future.apply_1.8.1         


# end ---------------------------------------------------------------------


