# Takahide Nejo, MD, PhD
# Okada Lab, Department of Neurological Surgery, UC San Francisco


# note ---------------------------------------------------------------------

# Nejo T, et al., Glioma-neuronal circuit remodeling induces regional immunosuppression
# Analysis of SB28 tumor-bearing mouse brain 10X Visium st-RNA-seq data

# Below is the script for the first sample. The data from the second sample was processed and analyzed similarly. 

# Using different versions of packages, such as SPATA2, may lead to varying outcomes, which we have not thoroughly validated.


# rm all ------------------------------------------------------------------

rm(list = ls(all.names = T))


# packages ----------------------------------------------------------------

library(SPATA2) 
library(tidyverse)
library(tidylog)
library(EBImage)
library(hdf5r)
library(msigdbr)


# check -------------------------------------------------------------------

packageVersion("SPATA2")
# [1] ‘2.0.4’

packageVersion("msigdbr")
# [1] ‘7.4.1’


# dir ---------------------------------------------------------------------

dir.1 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/2_01_stRNAseq_autodecoder_denoising/"
dir.2 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/2_02_stRNAseq_overview/"
if(dir.exists(dir.2) == F){dir.create(dir.2)}


# 1. load the spata-object ------------------------------------------------

setwd(dir.1)
in.f <- "visium_SB28_1_denoised.rds"
spata.obj <- readRDS(in.f)


# check -------------------------------------------------------------------

spata.obj 
# An object of class 'spata2' that contains 1 sample named 'A1_166L_Cas9'.

spata.obj %>% 
  getCoordsDf() %>% 
  print()

# # A tibble: 1,975 × 11
# barcodes               x     y tissue   row   col imagerow imagecol sample       section outline
# <chr>              <dbl> <dbl>  <int> <int> <int>    <dbl>    <dbl> <chr>        <chr>   <lgl>  
# 1 AAACAATCTACTAGCA-1  64.5  221.      1   589    43    -8400     1562 A1_166L_Cas9 1       FALSE  
# 2 AAACCGGGTAGGTACC-1 320.   163.      1   550    28    -9792     7742 A1_166L_Cas9 1       FALSE  
# 3 AAACCGTTCGTCCAGG-1 386.   215.      1   540    42    -8522     9333 A1_166L_Cas9 1       FALSE  
# 4 AAACCTCATGAAGTTG-1 287.   129.      1   555    19   -10609     6945 A1_166L_Cas9 1       FALSE  
# 5 AAACGAAGAACATACC-1  84.5  300.      1   586    64    -6488     2045 A1_166L_Cas9 1       TRUE   
# 6 AAACGAGACGGTTGAT-1 275.   355.      1   557    79    -5139     6650 A1_166L_Cas9 1       FALSE  
# 7 AAACTGCTGGCTCCAA-1 340.   310.      1   547    67    -6239     8232 A1_166L_Cas9 1       FALSE  
# 8 AAACTTGCAAACGTAT-1 339.   129.      1   547    19   -10614     8214 A1_166L_Cas9 1       FALSE  
# 9 AAAGACCCAAGTCGCG-1 110.   239.      1   582    48    -7949     2674 A1_166L_Cas9 1       FALSE  
# 10 AAAGACTGGGCGCTTT-1 234.   114.      1   563    15   -10969     5675 A1_166L_Cas9 1       FALSE  
# # ℹ 1,965 more rows
# # ℹ Use `print(n = ...)` to see more rows

spata.obj %>% 
  printGeneSetOverview()
# Class Available Gene Sets
# 1    BC                 289
# 2 BP.GO                7269
# 3 CC.GO                 972
# 4    HM                  50
# 5 MF.GO                1581
# 6  RCTM                1493


spata.obj %>% 
  getExpressionMatrixNames() %>% 
  print()
# [1] "scaled"   "denoised"

spata.obj %>% 
  getActiveMatrixName() %>% 
  print()
# [1] "denoised"

# spata.obj %>% 
#   plotSurfaceInteractive()

# prep MSigDB ------------------------------------------------------------------

msigdbr_collections() %>% print(n = Inf)

test_gs <-  msigdbr(
  species = "Mus musculus", # "Homo sapiens"
  category = "C2",
  # subcategory = "CGP"
) 


# check -------------------------------------------------------------------

test_gs %>% 
  pull(gs_name) %>% 
  unique() %>% 
  length() %>% 
  print()
# [1] 6290

test_gs %>% 
  colnames() %>% 
  print()
# [1] "gs_cat"               "gs_subcat"            "gs_name"              "gene_symbol"          "entrez_gene"         
# [6] "ensembl_gene"         "human_gene_symbol"    "human_entrez_gene"    "human_ensembl_gene"   "gs_id"               
# [11] "gs_pmid"              "gs_geoid"             "gs_exact_source"      "gs_url"               "gs_description"      
# [16] "taxon_id"             "ortholog_sources"     "num_ortholog_sources"

test_gs %>% 
  head() %>% 
  print()

test_gs %>% 
  dplyr::filter(grepl("VERHAAK_GLIOBLASTOMA_MESENCHYMAL|REACTOME_NEURONAL_SYSTEM", gs_name)) %>% 
  distinct(gs_name) %>% 
  print()
# # A tibble: 2 × 1
# gs_name                         
# <chr>                           
# 1 REACTOME_NEURONAL_SYSTEM        
# 2 VERHAAK_GLIOBLASTOMA_MESENCHYMAL

test_gs %>% 
  dplyr::filter(grepl("VERHAAK_GLIOBLASTOMA_MESENCHYMAL|REACTOME_NEURONAL_SYSTEM", gs_name)) %>% 
  pull(gs_name) %>% 
  table() %>% 
  print()
# REACTOME_NEURONAL_SYSTEM VERHAAK_GLIOBLASTOMA_MESENCHYMAL 
# 406                              221 


# prep the list of genesets -----------------------------------------------

# list.gs.of.interest <- list()
names.gs.of.interest <- test_gs %>% 
  dplyr::filter(grepl("VERHAAK_GLIOBLASTOMA_MESENCHYMAL|REACTOME_NEURONAL_SYSTEM", gs_name)) %>% 
  distinct(gs_name) %>% 
  pull()


# check -------------------------------------------------------------------

names.gs.of.interest %>% print()
# [1] "REACTOME_NEURONAL_SYSTEM"         "VERHAAK_GLIOBLASTOMA_MESENCHYMAL"


# add new geneset data ---------------------------------------------------

spata.obj.orig <- spata.obj
  
for(i in 1:length(names.gs.of.interest)){
  print(i)
  name.gs.i <- names.gs.of.interest[i]
  
  gs.i <- test_gs %>% 
    dplyr::filter(gs_name == name.gs.i)
  
  genes.i <- gs.i %>% 
    pull(gene_symbol)
  
  spata.obj <- spata.obj %>% 
    addGeneSet(
      class_name = "NEW", 
      gs_name = name.gs.i, 
      check_genes = F, 
      genes = genes.i, 
      overwrite = T
    ) 
}


# check -------------------------------------------------------------------

spata.obj@used_genesets %>% head()
# ont     gene
# 1 BC_RELA_PATHWAY    IKBKG
# 2 BC_RELA_PATHWAY   NFKBIA
# 3 BC_RELA_PATHWAY     FADD
# 4 BC_RELA_PATHWAY     CHUK
# 5 BC_RELA_PATHWAY      TNF
# 6 BC_RELA_PATHWAY TNFRSF1B

spata.obj@used_genesets %>% 
  dplyr::filter(grepl("NEW_", ont)) %>% 
  distinct(ont) %>% 
  print()
# distinct: removed 625 rows (>99%), 2 rows remaining
# ont
# 1         NEW_REACTOME_NEURONAL_SYSTEM
# 2 NEW_VERHAAK_GLIOBLASTOMA_MESENCHYMAL


# update the featureDf ----------------------------------------------------

coords_df <- getCoordsDf(spata.obj)
  
# join this data.frame with additional info

joined_df <- spata.obj %>%
  joinWith(
    spata_df = coords_df,
    gene_sets = paste0("NEW_", names.gs.of.interest), # expression values of the gene set Hallmark-Hypoxia
    verbose = FALSE)


# check -------------------------------------------------------------------

joined_df %>% colnames()
# [1] "barcodes"                             "x"                                   
# [3] "y"                                    "tissue"                              
# [5] "row"                                  "col"                                 
# [7] "imagerow"                             "imagecol"                            
# [9] "sample"                               "section"                             
# [11] "outline"                              "NEW_REACTOME_NEURONAL_SYSTEM"        
# [13] "NEW_VERHAAK_GLIOBLASTOMA_MESENCHYMAL"

colnames(joined_df)[12:13] <- c("RCTM_NEU", "VRHK_MES")

joined_df %>% colnames()
# [1] "barcodes" "x"        "y"        "tissue"   "row"      "col"      "imagerow"
# [8] "imagecol" "sample"   "section"  "outline"  "RCTM_NEU" "VRHK_MES"

joined_df_simple <- joined_df %>% 
  dplyr::select(barcodes, colnames(joined_df)[12:13])
  

# add the feature

spata.obj <- spata.obj %>% 
  addFeatures(
    feature_names = colnames(joined_df)[12:13],
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
# [1] "barcodes"            "sample"              "orig.ident"         
# [4] "nCount_Spatial"      "nFeature_Spatial"    "percent.mt"         
# [7] "percent.RB"          "Spatial_snn_res.0.8" "seurat_clusters"    
# [10] "histology"           "RCTM_NEU"            "VRHK_MES"   


# save --------------------------------------------------------------------

START.TIME <- Sys.time() 

setwd(dir.2)
out.f <- "visium_SB28_1_ftr_edited.rds"
saveRDS(spata.obj, out.f)

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 1.202005 mins


# visualize (Fig 2d-f) ---------------------------------------------------------------

# H/E

g <- spata.obj %>% 
  plotSurface(
    display_image = T, 
    # color_by = "seurat_clusters", 
    pt_size = 4, 
    pt_alpha = 0
  )
plot(g)

setwd(dir.2)
out.f <- "01_surfacePlot_01_HE.pdf"
ggsave(out.f, g, w = 6, h = 4.5)


# seurat_clusters

g <- spata.obj %>% 
  plotSurface(
    display_image = F, 
    color_by = "seurat_clusters", 
    pt_size = 2, 
    pt_alpha = 0.9, 
  ) 
plot(g)

setwd(dir.2)
out.f <- "02_surfacePlot_02_seurat_cl.pdf"
ggsave(out.f, g, w = 6, h = 4.5)


# RCTM_NEU

g <- spata.obj %>% 
  plotSurface(
    display_image = F, 
    color_by = "RCTM_NEU", 
    pt_size = 2, 
    pt_alpha = 0.9, 
  )
plot(g)

setwd(dir.2)
out.f <- "03_surfacePlot_03_rctm_neu.pdf"
ggsave(out.f, g, w = 6, h = 4.5)

# RCTM_NEU

g <- spata.obj %>% 
  plotSurface(
    display_image = F, 
    color_by = "VRHK_MES", 
    pt_size = 2, 
    pt_alpha = 0.9, 
  )
plot(g)

setwd(dir.2)
out.f <- "04_surfacePlot_04_vrhk_mes.pdf"
ggsave(out.f, g, w = 6, h = 4.5)


# si ----------------------------------------------------------------------

Sys.time()
# [1] "2024-12-26 14:25:31 PST"

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
# [1] shiny_1.10.0    msigdbr_7.4.1   hdf5r_1.3.9     EBImage_4.34.0  tidylog_1.0.2  
# [6] forcats_0.5.1   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2     readr_2.1.5    
# [11] tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_1.3.1 SPATA2_2.0.4   
# 
# loaded via a namespace (and not attached):
# [1] utf8_1.2.4                  shinydashboard_0.7.2        spatstat.explore_3.2-1     
# [4] reticulate_1.40.0           tidyselect_1.2.1            htmlwidgets_1.6.4          
# [7] grid_4.1.3                  Rtsne_0.16                  munsell_0.5.1              
# [10] ragg_1.3.3                  codetools_0.2-18            ica_1.0-2                  
# [13] units_0.8-5                 future_1.25.0               miniUI_0.1.1.1             
# [16] withr_3.0.2                 spatstat.random_3.1-5       colorspace_2.1-1           
# [19] progressr_0.10.0            Biobase_2.52.0              rstudioapi_0.17.1          
# [22] Seurat_4.3.0.1              stats4_4.1.3                SingleCellExperiment_1.14.1
# [25] ROCR_1.0-11                 tensor_1.5                  shinyWidgets_0.8.7         
# [28] listenv_0.8.0               labeling_0.4.3              MatrixGenerics_1.4.3       
# [31] GenomeInfoDbData_1.2.6      polyclip_1.10-7             bit64_4.5.2                
# [34] farver_2.1.2                parallelly_1.31.1           vctrs_0.6.5                
# [37] generics_0.1.3              timechange_0.3.0            R6_2.5.1                   
# [40] GenomeInfoDb_1.28.1         shinybusy_0.3.3             locfit_1.5-9.4             
# [43] bitops_1.0-9                spatstat.utils_3.1-1        cachem_1.1.0               
# [46] DelayedArray_0.18.0         assertthat_0.2.1            promises_1.3.2             
# [49] scales_1.3.0                gtable_0.3.6                globals_0.15.0             
# [52] goftest_1.2-2               rlang_1.1.4                 clisymbols_1.2.0           
# [55] systemfonts_1.1.0           splines_4.1.3               lazyeval_0.2.2             
# [58] spatstat.geom_3.2-4         broom_1.0.7                 reshape2_1.4.4             
# [61] abind_1.4-5                 modelr_0.1.8                backports_1.5.0            
# [64] httpuv_1.6.15               tools_4.1.3                 jquerylib_0.1.4            
# [67] RColorBrewer_1.1-3          BiocGenerics_0.38.0         ggridges_0.5.6             
# [70] Rcpp_1.0.13-1               plyr_1.8.7                  zlibbioc_1.38.0            
# [73] RCurl_1.98-1.3              deldir_1.0-6                pbapply_1.5-0              
# [76] cowplot_1.1.1               S4Vectors_0.30.2            fontawesome_0.5.3          
# [79] zoo_1.8-10                  SeuratObject_4.1.3          SummarizedExperiment_1.22.0
# [82] haven_2.4.3                 ggrepel_0.9.1               cluster_2.1.2              
# [85] fs_1.6.5                    magrittr_2.0.3              magick_2.8.5               
# [88] data.table_1.16.4           scattermore_1.2             lmtest_0.9-40              
# [91] reprex_2.0.1                RANN_2.6.1                  fitdistrplus_1.1-5         
# [94] anndata_0.7.5.6             matrixStats_0.62.0          hms_1.1.3                  
# [97] patchwork_1.3.0             mime_0.12                   fftwtools_0.9-11           
# [100] xtable_1.8-4                jpeg_0.1-9                  readxl_1.3.1               
# [103] IRanges_2.26.0              gridExtra_2.3               compiler_4.1.3             
# [106] KernSmooth_2.23-20          crayon_1.5.3                SPATAData_0.0.0.9000       
# [109] htmltools_0.5.8.1           later_1.4.1                 tzdb_0.4.0                 
# [112] tiff_0.1-11                 lubridate_1.9.4             DBI_1.2.3                  
# [115] dbplyr_2.1.1                MASS_7.3-54                 babelgene_21.4             
# [118] Matrix_1.6-4                cli_3.6.3                   parallel_4.1.3             
# [121] igraph_1.3.1                GenomicRanges_1.44.0        pkgconfig_2.0.3            
# [124] sp_2.1-4                    confuns_1.0.3               plotly_4.10.4              
# [127] spatstat.sparse_3.0-2       xml2_1.3.2                  bslib_0.8.0                
# [130] XVector_0.32.0              rvest_1.0.1                 digest_0.6.37              
# [133] sctransform_0.3.5           RcppAnnoy_0.0.22            spatstat.data_3.0-1        
# [136] cellranger_1.1.0            leiden_0.3.9                uwot_0.1.16                
# [139] lifecycle_1.0.4             nlme_3.1-152                jsonlite_1.8.9             
# [142] viridisLite_0.4.2           pillar_1.10.0               lattice_0.20-44            
# [145] fastmap_1.2.0               httr_1.4.7                  survival_3.2-12            
# [148] glue_1.8.0                  png_0.1-8                   bit_4.5.0.1                
# [151] stringi_1.8.4               sass_0.4.9                  textshaping_0.3.6          
# [154] memoise_2.0.1               irlba_2.3.5.1               future.apply_1.8.1      


# end ---------------------------------------------------------------------
