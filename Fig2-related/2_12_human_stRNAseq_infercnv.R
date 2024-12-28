# Takahide Nejo, MD, PhD
# Okada Lab, Department of Neurological Surgery, UC San Francisco


# note ---------------------------------------------------------------------

# Nejo T, et al., Glioma-neuronal circuit remodeling induces regional immunosuppression
# Analysis of SPATAData human GBM 10X Visium st-RNA-seq data

# Below is the script for one sample ("260_T"). The data from the other samples were processed and analyzed similarly. 

# Using different versions of packages, such as SPATA2, may lead to varying outcomes, which we have not thoroughly validated.

# The RDS object "visium_SPATAData_260_T_denoised.rds" is used as input. For more details, please refer to the previous step, "2_11_human_stRNAseq_autodecoder_denoising.R".

# ref: https://themilolab.github.io/SPATA2/articles/spata-data.html


# rm all ------------------------------------------------------------------

rm(list = ls(all.names = T))


# packages ----------------------------------------------------------------

library(SPATA2)
library(tidyverse)
library(tidylog)
library(infercnv)
library(ggpubr)


# check -------------------------------------------------------------------

packageVersion("SPATA2")
# [1] ‘2.0.4’

packageVersion("infercnv")
# [1] ‘1.8.1’


# dir ---------------------------------------------------------------------

dir.1 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/2_11_human_stRNAseq_autodecoder_denoising/"
dir.2 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/2_12_human_stRNAseq_infercnv/"
dir.3 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/2_12_human_stRNAseq_infercnv/cnv_results"

if(dir.exists(dir.2) == F){dir.create(dir.2)}
if(dir.exists(dir.3) == F){dir.create(dir.3)}


# load the spata-object ------------------------------------------------

START.TIME <- Sys.time()

setwd(dir.1) 
in.f <- "visium_SPATAData_260_T_denoised.rds"
spata.obj <- readRDS(in.f)

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 10.31541 secs


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


# run CNV analysis -------------------------------------------------------------------

# https://themilolab.github.io/SPATA2/articles/spata-v2-cnv-analysis.html#running-cnv-analysis

setwd(dir.3)

START.TIME <- Sys.time()

spata.obj <- spata.obj %>% 
  runCnvAnalysis(
    directory_cnv_folder = dir.3, 
    cnv_prefix = "Chr"
  )

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 14.88044 mins


# check -------------------------------------------------------------------

cnv_results <- spata.obj %>% 
  getCnvResults()

names(cnv_results) %>% print()
# [1] "prefix"      "cnv_df"      "cnv_mtr"     "gene_pos_df" "regions_df" 


# save --------------------------------------------------------------------

START.TIME <- Sys.time() 

setwd(dir.2)
out.f <- "visium_SPATAData_260_T_infercnv.rds"
saveRDS(spata.obj, out.f)

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 1.798559 mins


# visualize ---------------------------------------------------------------

# heatmap
g <- spata.obj %>% 
  plotCnvHeatmap() 

setwd(dir.2)
out.f <- "01_heatmap_infercnv.pdf"
ggsave(out.f, g, width = 16, height = 9)


# surface plot
g1 <- spata.obj %>% 
  plotSurface(
    color_by = "Chr7", 
    pt_clrsp = "Reds"
  )

g2 <- spata.obj %>% 
  plotSurface(
    color_by = "Chr10", 
    pt_clrsp = "Oslo"
  )

g <- ggarrange(
  plotlist = list(g1, g2), 
  ncol = 2, nrow = 1
) 

setwd(dir.2)
out.f <- "02_surfaceplot_infercnv.pdf"
ggsave(out.f, g, width = 16, height = 9)


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
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] ggpubr_0.4.0    infercnv_1.8.1  tidylog_1.0.2   forcats_0.5.1   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2    
# [8] readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_1.3.1 SPATA2_2.0.4   
# 
# loaded via a namespace (and not attached):
# [1] utf8_1.2.4                  spatstat.explore_3.2-1      reticulate_1.40.0           tidyselect_1.2.1           
# [5] htmlwidgets_1.6.4           grid_4.1.3                  Rtsne_0.16                  munsell_0.5.1              
# [9] ragg_1.3.3                  codetools_0.2-18            ica_1.0-2                   units_0.8-5                
# [13] future_1.25.0               miniUI_0.1.1.1              withr_3.0.2                 argparse_2.1.5             
# [17] spatstat.random_3.1-5       colorspace_2.1-1            progressr_0.10.0            Biobase_2.52.0             
# [21] rstudioapi_0.17.1           Seurat_4.3.0.1              stats4_4.1.3                SingleCellExperiment_1.14.1
# [25] ROCR_1.0-11                 ggsignif_0.6.2              tensor_1.5                  listenv_0.8.0              
# [29] labeling_0.4.3              MatrixGenerics_1.4.3        GenomeInfoDbData_1.2.6      polyclip_1.10-7            
# [33] farver_2.1.2                TH.data_1.1-1               coda_0.19-4                 parallelly_1.31.1          
# [37] vctrs_0.6.5                 generics_0.1.3              lambda.r_1.2.4              timechange_0.3.0           
# [41] fastcluster_1.2.3           R6_2.5.1                    doParallel_1.0.16           GenomeInfoDb_1.28.1        
# [45] locfit_1.5-9.4              gridGraphics_0.5-1          reshape_0.8.9               bitops_1.0-9               
# [49] spatstat.utils_3.1-1        DelayedArray_0.18.0         assertthat_0.2.1            promises_1.3.2             
# [53] scales_1.3.0                multcomp_1.4-19             gtable_0.3.6                globals_0.15.0             
# [57] goftest_1.2-2               sandwich_3.0-1              rlang_1.1.4                 clisymbols_1.2.0           
# [61] systemfonts_1.1.0           splines_4.1.3               rstatix_0.7.0               lazyeval_0.2.2             
# [65] rjags_4-16                  spatstat.geom_3.2-4         broom_1.0.7                 reshape2_1.4.4             
# [69] abind_1.4-5                 modelr_0.1.8                backports_1.5.0             httpuv_1.6.15              
# [73] tools_4.1.3                 ggplotify_0.1.2             gplots_3.2.0                RColorBrewer_1.1-3         
# [77] BiocGenerics_0.38.0         phyclust_0.1-30             ggridges_0.5.6              Rcpp_1.0.13-1              
# [81] plyr_1.8.7                  zlibbioc_1.38.0             RCurl_1.98-1.3              deldir_1.0-6               
# [85] pbapply_1.5-0               cowplot_1.1.1               S4Vectors_0.30.2            zoo_1.8-10                 
# [89] SeuratObject_4.1.3          SummarizedExperiment_1.22.0 haven_2.4.3                 ggrepel_0.9.1              
# [93] cluster_2.1.2               fs_1.6.5                    magrittr_2.0.3              magick_2.8.5               
# [97] data.table_1.16.4           futile.options_1.0.1        scattermore_1.2             openxlsx_4.2.4             
# [101] lmtest_0.9-40               reprex_2.0.1                RANN_2.6.1                  mvtnorm_1.1-3              
# [105] fitdistrplus_1.1-5          anndata_0.7.5.6             matrixStats_0.62.0          hms_1.1.3                  
# [109] patchwork_1.3.0             mime_0.12                   fftwtools_0.9-11            xtable_1.8-4               
# [113] rio_0.5.27                  jpeg_0.1-9                  readxl_1.3.1                IRanges_2.26.0             
# [117] gridExtra_2.3               compiler_4.1.3              KernSmooth_2.23-20          crayon_1.5.3               
# [121] SPATAData_0.0.0.9000        htmltools_0.5.8.1           ggfun_0.0.4                 later_1.4.1                
# [125] tzdb_0.4.0                  tiff_0.1-11                 aplot_0.1.1                 libcoin_1.0-9              
# [129] lubridate_1.9.4             DBI_1.2.3                   formatR_1.11                dbplyr_2.1.1               
# [133] MASS_7.3-54                 car_3.0-11                  Matrix_1.6-4                cli_3.6.3                  
# [137] parallel_4.1.3              igraph_1.3.1                GenomicRanges_1.44.0        pkgconfig_2.0.3            
# [141] foreign_0.8-81              coin_1.4-2                  sp_2.1-4                    confuns_1.0.3              
# [145] plotly_4.10.4               spatstat.sparse_3.0-2       xml2_1.3.2                  foreach_1.5.1              
# [149] XVector_0.32.0              rvest_1.0.1                 yulab.utils_0.1.8           digest_0.6.37              
# [153] sctransform_0.3.5           RcppAnnoy_0.0.22            spatstat.data_3.0-1         cellranger_1.1.0           
# [157] leiden_0.3.9                edgeR_3.34.1                uwot_0.1.16                 curl_6.0.1                 
# [161] gtools_3.9.5                shiny_1.10.0                EBImage_4.34.0              modeltools_0.2-23          
# [165] lifecycle_1.0.4             nlme_3.1-152                jsonlite_1.8.9              carData_3.0-4              
# [169] futile.logger_1.4.3         limma_3.48.3                viridisLite_0.4.2           pillar_1.10.0              
# [173] lattice_0.20-44             fastmap_1.2.0               httr_1.4.7                  survival_3.2-12            
# [177] glue_1.8.0                  zip_2.2.0                   png_0.1-8                   iterators_1.0.13           
# [181] stringi_1.8.4               textshaping_0.3.6           caTools_1.18.3              irlba_2.3.5.1              
# [185] future.apply_1.8.1          ape_5.6                    


# end ---------------------------------------------------------------------


