# Takahide Nejo, MD, PhD
# Okada Lab, Department of Neurological Surgery, UC San Francisco


# note ---------------------------------------------------------------------

# Nejo T, et al., Glioma-neuronal circuit remodeling induces regional immunosuppression

# Analysis of GL261 tumor-bearing mouse brain 10X Visium st-RNA-seq data (GSE245263)
# original article: García-Vicente L, 2023 bioRxiv. 
# https://www.biorxiv.org/content/10.1101/2023.10.26.564166v1.full

# original data analysis info: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7839621

# Using different versions of packages, such as SPATA2, may lead to varying outcomes, which we have not thoroughly validated.


# rm all ------------------------------------------------------------------

rm(list = ls(all.names = T))


# packages ----------------------------------------------------------------

library(tidyverse)
library(tidylog)
library(SPATA2) # devtools::install_github("theMILOlab/SPATA2", ref = 'v2.0.4')
  # devtools::install_github("kueckelj/confuns")
library(hdf5r)
library(EBImage)
library(ggpubr)
library(R.utils)
library(Seurat)

Sys.setenv(RETICULATE_PYTHON = "/usr/bin/python3")
library(reticulate)
library(SCP) # remotes::install_github("zhanghao-njmu/SCP",upgrade="never")
  # devtools::install_github("jokergoo/ComplexHeatmap") # BiocManager::install("ComplexHeatmap", force = T)
  # devtools::install_github("jokergoo/circlize")


# check -------------------------------------------------------------------

packageVersion("SPATA2")
# [1] ‘2.0.4’

packageVersion("reticulate")
# [1] ‘1.40.0’

packageVersion("SCP")
# [1] ‘0.5.6’

packageVersion("Seurat")
# [1] ‘4.3.0.1’

py_config() %>% print()
# python:         /usr/bin/python3
# libpython:      /usr/lib64/libpython3.11.so
# pythonhome:     //usr://usr
# version:        3.11.7 (main, Jan 26 2024, 20:24:17) [GCC 8.5.0 20210514 (Red Hat 8.5.0-21)]
# numpy:          /c4/home/tnejo/.local/lib/python3.11/site-packages/numpy
# numpy_version:  2.0.2
# 
# NOTE: Python version was forced by RETICULATE_PYTHON


# dir ---------------------------------------------------------------------

dir.1 <- "/okadalab/data1/tnejo/neuro_immune_proj/10x_visium_03_GSE245263_GL261/dl"
dir.2 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/2_20_GL261_GSE245263_stRNAseq_data_setup/"
if(dir.exists(dir.2) == F){dir.create(dir.2)}


# check -------------------------------------------------------------------

setwd(dir.1)
list.files()
{
  # [1] "GSE245263_RAW.tar"                                           
  # [2] "GSM7839621_10xVisium_annotated.h5ad.gz"                      
  # [3] "GSM7839621_aligned_fiducials.jpg.gz"                         
  # [4] "GSM7839621_CMS75-GBMA7_brain_visium_Run1_1_HE_10x_001.jpg.gz"
  # [5] "GSM7839621_CMS75-GBMA7_brain_visium_Run1_1_HE_10x_001.tif.gz"
  # [6] "GSM7839621_detected_tissue_image.jpg.gz"                     
  # [7] "GSM7839621_scalefactors_json.json.gz"                        
  # [8] "GSM7839621_tissue_hires_image.png.gz"                        
  # [9] "GSM7839621_tissue_lowres_image.png.gz"                       
  # [10] "GSM7839621_tissue_positions_list.csv.gz" 
}


# decompress several necessary files ---------------------------------------------------------------

setwd(dir.1)
gunzip("GSM7839621_10xVisium_annotated.h5ad.gz", remove = F)
gunzip("GSM7839621_tissue_hires_image.png.gz", remove = F)
gunzip("GSM7839621_tissue_lowres_image.png.gz", remove = F)


# check -------------------------------------------------------------------

# in.f <- "GSM7839621_10xVisium_annotated.h5ad"  
# rhdf5::h5dump(in.f, load = FALSE) %>% View() # data structure
# rhdf5::h5ls(in.f) %>% as_tibble() %>% View() # data type

# ref: https://github.com/theislab/zellkonverter/issues/45


# load the original anndata-------------------------------------------------------------------

# ref: https://www.biostars.org/p/9546692/

START.TIME <- Sys.time() 

sc <- import("scanpy")

setwd(dir.1)
in.f <- "GSM7839621_10xVisium_annotated.h5ad"
adata.obj <- sc$read_h5ad(in.f)

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 39.79031 secs


# check -------------------------------------------------------------------

adata.obj %>% print()
{
  # AnnData object with n_obs × n_vars = 2158 × 32272
  # obs: 'in_tissue', 'array_row', 'array_col', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'condition', 'coarse', 'fine', 'leiden'
  # var: 'gene_ids', 'feature_types', 'genome', 'mt', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'mean', 'std'
  # uns: 'condition_colors', 'leiden', 'leiden_colors', 'log1p', 'neighbors', 'pca', 'spatial', 't-test', 'umap', 'wilcoxon'
  # obsm: 'X_pca', 'X_spatial', 'X_umap', 'spatial'
  # varm: 'PCs'
  # layers: 'counts'
  # obsp: 'connectivities', 'distances'
}

adata.obj$layers %>% names()
# [1] "counts"


# conversion to a SPATA2 object --------------------------------------------------------------

START.TIME <- Sys.time() 

spata.obj <- SPATA2::asSPATA2(
  adata.obj,
  sample_name = "GSM7839621", 
  normalized_mtr_name = NULL, 
  scaled_mtr_name = NULL
)

# 11:16:12 Setting up new `spata2` object.
# 11:16:12 Transferring data.
# 11:16:13 No spatial trajectories found. Returning input object.
# 11:16:13 No image annotations found. Returning input object.
# The AnnData object contains the following layers: counts
# Using adata$layers['counts'] as count matrix
# 11:16:16 Done.
# Warning messages:
# 1: In load_adata_matrix(adata = object, count_mtr_name = count_mtr_name,  :
#   No normalized matrix found. Using adata$X as normalized matrix. If you want to use a different matrix,
#               specify a name for the normalized matrix via `normalized_mtr_name`
# 2: In load_adata_matrix(adata = object, count_mtr_name = count_mtr_name,  :
#   No scaled matrix found to import. You can specify the scaled matrix AnnData layer via
#               `scaled_mtr_name` (e.g. scaled_mtr_name='scaled_data')
# 3: In value[[3L]](cond) :
#   Could not find or transfer TSNE data. Did you process the AnnData object correctly?
                                                
                                                    
FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 3.414258 secs


# check -------------------------------------------------------------------

spata.obj %>% print()
# An object of class 'spata2' that contains 1 sample named 'GSM7839621'.

getFeatureDf(spata.obj) %>% dim() %>% print()
# [1] 2158   20

getFeatureDf(spata.obj) %>% colnames() %>% print()
# [1] "barcodes"                    "sample"                      "in_tissue"                  
# [4] "array_row"                   "array_col"                   "n_genes_by_counts"          
# [7] "log1p_n_genes_by_counts"     "total_counts"                "log1p_total_counts"         
# [10] "pct_counts_in_top_50_genes"  "pct_counts_in_top_100_genes" "pct_counts_in_top_200_genes"
# [13] "pct_counts_in_top_500_genes" "total_counts_mt"             "log1p_total_counts_mt"      
# [16] "pct_counts_mt"               "condition"                   "coarse"                     
# [19] "fine"                        "leiden"       

getFeatureDf(spata.obj) %>% head() %>% print()
# # A tibble: 6 × 20
# barcodes        sample in_tissue array_row array_col n_genes_by_counts log1p_n_genes_by_cou…¹
# <chr>           <chr>      <dbl>     <dbl>     <dbl>             <dbl>                  <dbl>
# 1 AAACAAGTATCTCC… GSM78…         1        50       102              2033                   9.91
# 2 AAACAATCTACTAG… GSM78…         1         3        43              3428                   9.91
# 3 AAACCCGAACGAAA… GSM78…         1        45       115              2698                   9.91
# 4 AAACCGTTCGTCCA… GSM78…         1        52        42              3378                   9.91
# 5 AAACGAAGAACATA… GSM78…         1         6        64              3939                   9.91
# 6 AAACGAGACGGTTG… GSM78…         1        35        79              2567                   9.91
# # ℹ abbreviated name: ¹​log1p_n_genes_by_counts
# # ℹ 13 more variables: total_counts <dbl>, log1p_total_counts <dbl>,
# #   pct_counts_in_top_50_genes <dbl>, pct_counts_in_top_100_genes <dbl>,
# #   pct_counts_in_top_200_genes <dbl>, pct_counts_in_top_500_genes <dbl>,
# #   total_counts_mt <dbl>, log1p_total_counts_mt <dbl>, pct_counts_mt <dbl>, condition <fct>,
# #   coarse <fct>, fine <fct>, leiden <fct>

getCoordsDf(spata.obj) %>% dim() %>% print()
# [1] 2158   10

getCoordsDf(spata.obj) %>% colnames() %>% print()
# [1] "barcodes"   "imagerow.x" "imagecol.x" "sample.x"   "x"          "y"          "imagerow.y"
# [8] "imagecol.y" "sample.y"   "sample"   

getCoordsDf(spata.obj) %>% head() %>% print()
# # A tibble: 6 × 10
# barcodes     imagerow.x imagecol.x sample.x     x     y imagerow.y imagecol.y sample.y sample
# <chr>             <dbl>      <dbl> <chr>    <dbl> <dbl>      <dbl>      <dbl> <chr>    <chr> 
# 1 AAACAAGTATC…       9669       7965 GSM7839…  158.  402.       9669       7965 GSM7839… GSM78…
# 2 AAACAATCTAC…      21459      16591 GSM7839…  396.  228.      21459      16591 GSM7839… GSM78…
# 3 AAACCCGAACG…      10942       6090 GSM7839…  184.  440.      10942       6090 GSM7839… GSM78…
# 4 AAACCGTTCGT…       9107      16653 GSM7839…  147.  227.       9107      16653 GSM7839… GSM78…
# 5 AAACGAAGAAC…      20723      13544 GSM7839…  382.  289.      20723      13544 GSM7839… GSM78…
# 6 AAACGAGACGG…      13428      11322 GSM7839…  234.  334.      13428      11322 GSM7839… GSM78…


# check

plotSurface(
  spata.obj, 
  pt_alpha = 1, 
  color_by = "n_genes_by_counts"
)


# add image data --------------------------------------------------------------

# ref: https://themilolab.github.io/SPATA2/articles/spata-v2-image-handling.html#additional-images

START.TIME <- Sys.time()

spata.obj <- spata.obj %>% 
  setImageDirLowres(
    dir = paste0(dir.1, "/", "GSM7839621_tissue_lowres_image.png"), 
    name = "lowres"
  ) %>% 
  setImageDirHighres(
    dir = paste0(dir.1, "/", "GSM7839621_tissue_hires_image.png"), 
    name = "hires"
  ) %>% 
  loadImageLowres() %>% 
  loadImageHighres()

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 4.544086 secs


# check -------------------------------------------------------------------

spata.obj %>% 
  getImageDirectories() %>% 
  print()
# lowres
# "/okadalab/data1/tnejo/neuro_immune_proj/10x_visium_03_GSE245263_GL261/dl/GSM7839621_tissue_lowres_image.png" 
# highres 
# "/okadalab/data1/tnejo/neuro_immune_proj/10x_visium_03_GSE245263_GL261/dl/GSM7839621_tissue_hires_image.png" 


# check -------------------------------------------------------------------

spata.obj %>% 
  plotImageGgplot(
    unit = "mm", 
    xrange = c("0mm", "10mm"),
    yrange = c("0mm", "10mm")
  ) 

g <- spata.obj %>% 
  plotSurface(
    display_image = T, 
    pt_alpha = 0.25, 
    pt_clr = "yellow"
    ) + 
  ggpLayerThemeCoords(unit = "mm") + 
  xlim(0, 2000) + ylim(0, 2000)
plot(g)

setwd(dir.2)
out.f <- "GSE245263_01_surfaceplot_orig_coord.pdf"
ggsave(out.f, g, w = 6, h = 4.5)


# note --------------------------------------------------------------------

# coordinates information needs to be corrected. 


# update coordinates info -------------------------------------------------

spata.obj.orig <- spata.obj

coord_df <- spata.obj %>% 
  getCoordsDf()

# check
coord_df %>% print()


coord_df.clean <- coord_df %>% 
  dplyr::select(barcodes, sample, imagerow = imagerow.x, imagecol = imagecol.x, x, y)

# re-adjustment

coord_df.new <- coord_df.clean %>% 
  mutate(x = x + 125) %>% 
  mutate(y = y + 125)

# These correction coefficients (x = 125; y = 125) were manually identified .


spata.obj <- spata.obj %>% 
  setCoordsDf(coord_df.new)


# check -------------------------------------------------------------------

g <- spata.obj.orig %>% 
  plotSurface(
    display_image = T, 
    pt_alpha = 0.25, 
    pt_clr = "yellow"
  ) + 
  ggpLayerThemeCoords(unit = "mm") + 
  xlim(0, 2000) + ylim(0, 2000)
# plot(g)
g1 <- g

g <- spata.obj %>% 
  plotSurface(
    display_image = T, 
    pt_alpha = 0.25, 
    pt_clr = "orange"
  ) + 
  ggpLayerThemeCoords(unit = "mm") + 
  xlim(0, 2000) + ylim(0, 2000)
# plot(g)
g2 <- g

g <- ggarrange(
  g1, g2, 
  ncol = 2, nrow = 1
)

setwd(dir.2)
out.f <- "GSE245263_02_surfaceplot_readjusted_coord.pdf"
ggsave(out.f, g, w = 12, h = 4.5)


# check -------------------------------------------------------------------

spata.obj.orig <- spata.obj

g <- spata.obj.orig %>% 
  plotSurface(
    display_image = T, 
    color_by = "Aif1", 
    pt_alpha = 0.75,  
    pt_size = 1, 
    smooth_span = 0.25
  ) + 
  ggpLayerThemeCoords(unit = "mm") + 
  xlim(0, 2000) + ylim(0, 2000)
plot(g)
g1 <- g


# add the newly processed "scaled" data --------------------------------------------------------------

# ref: https://themilolab.github.io/SPATA2/articles/adding-data.html

count_mtr <- spata.obj %>% 
  getCountMatrix()

# check
count_mtr %>% class() %>% print()
# [1] "dgRMatrix"
# attr(,"package")
# [1] "Matrix"

count_mtr %>% dim() %>% print()
# [1] 32272  2158

count_mtr[1:5, 1:5] %>% print()
# 5 x 5 sparse Matrix of class "dgRMatrix"
# AAACAAGTATCTCCCA-1 AAACAATCTACTAGCA-1 AAACCCGAACGAAATC-1 AAACCGTTCGTCCAGG-1 AAACGAAGAACATACC-1
# Xkr4                     .                  .                  .                  .                  .
# Gm1992                   .                  .                  .                  .                  .
# Gm19938                  .                  .                  .                  1                  .
# Gm37381                  .                  .                  .                  .                  .
# Rp1                      .                  .                  .                  .                  .

seurat.obj <- count_mtr %>% 
  CreateSeuratObject() %>% 
  NormalizeData() %>% 
  ScaleData()
scaled.mtr <- seurat.obj %>% 
  GetAssayData(layer = "scale.data")

# check: show a small subset of the scaled mtr 
scaled.mtr[1:5, 1:5]
# AAACAAGTATCTCCCA-1 AAACAATCTACTAGCA-1 AAACCCGAACGAAATC-1 AAACCGTTCGTCCAGG-1 AAACGAAGAACATACC-1
# Xkr4           -0.04844399        -0.04844399        -0.04844399        -0.04844399        -0.04844399
# Gm1992         -0.02152654        -0.02152654        -0.02152654        -0.02152654        -0.02152654
# Gm19938        -0.06797521        -0.06797521        -0.06797521         9.40862993        -0.06797521
# Gm37381         0.00000000         0.00000000         0.00000000         0.00000000         0.00000000
# Rp1            -0.03044909        -0.03044909        -0.03044909        -0.03044909        -0.03044909


# add the processed matrix to SPATA2 object under the name 'scaled'

spata.obj <- spata.obj %>% 
  addExpressionMatrix(
    expr_mtr = scaled.mtr,
    mtr_name = "scaled"
  )

spata.obj <- spata.obj %>%
  setActiveMatrix("scaled")


# check -------------------------------------------------------------------

# all expression matrices before denoising

spata.obj %>% 
  getExpressionMatrixNames() %>% 
  print()
# [1] "normalized" "scaled"    

# active expression matrix before denoising
spata.obj %>% 
  getActiveMatrixName() %>% 
  print()
# [1] "scaled"


# visualize ---------------------------------------------------------------

g <- spata.obj %>% 
  plotSurface(
    display_image = T, 
    color_by = "Aif1", 
    pt_alpha = 0.75,  
    pt_size = 1, 
    smooth_span = 0.25
  ) + 
  ggpLayerThemeCoords(unit = "mm") + 
  xlim(0, 2000) + ylim(0, 2000)
plot(g)
g2 <- g


g <- ggarrange(
  g1, g2, 
  ncol = 2, nrow = 1
)

setwd(dir.2)
out.f <- "GSE245263_03_surfaceplot_pre_n_post_scaling.pdf"
ggsave(out.f, g, w = 12, h = 4.5)



# save --------------------------------------------------------------------

START.TIME <- Sys.time() 

setwd(dir.2)
out.f <- "visium_GSE245263_GL261_coord_adjusted.rds"
saveRDS(spata.obj, out.f)

FINISH.TIME <- Sys.time() 

print(FINISH.TIME - START.TIME)
# Time difference of 39.4155 secs


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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] SeuratObject_5.0.2 Seurat_4.3.0.1     SCP_0.5.6          reticulate_1.40.0  R.utils_2.11.0    
# [6] R.oo_1.24.0        R.methodsS3_1.8.1  ggpubr_0.4.0       EBImage_4.34.0     hdf5r_1.3.9       
# [11] SPATA2_2.0.4       tidylog_1.0.2      forcats_0.5.1      stringr_1.5.1      dplyr_1.1.4       
# [16] purrr_1.0.2        readr_2.1.5        tidyr_1.3.1        tibble_3.2.1       ggplot2_3.5.1     
# [21] tidyverse_1.3.1   
# 
# loaded via a namespace (and not attached):
# [1] rappdirs_0.3.3              scattermore_1.2             princurve_2.1.6            
# [4] SPATAData_0.0.0.9000        ragg_1.3.3                  bit64_4.5.2                
# [7] irlba_2.3.5.1               DelayedArray_0.18.0         TrajectoryUtils_1.2.0      
# [10] data.table_1.16.4           KEGGREST_1.32.0             RCurl_1.98-1.3             
# [13] doParallel_1.0.16           generics_0.1.3              BiocGenerics_0.38.0        
# [16] cowplot_1.1.1               RSQLite_2.2.8               shadowtext_0.1.0           
# [19] RANN_2.6.1                  future_1.25.0               enrichplot_1.12.3          
# [22] bit_4.5.0.1                 tzdb_0.4.0                  spatstat.data_3.0-1        
# [25] xml2_1.3.2                  lubridate_1.9.4             httpuv_1.6.15              
# [28] SummarizedExperiment_1.22.0 assertthat_0.2.1            confuns_1.0.3              
# [31] viridis_0.6.5               hms_1.1.3                   promises_1.3.2             
# [34] progress_1.2.3              dbplyr_2.1.1                readxl_1.3.1               
# [37] igraph_1.3.1                DBI_1.2.3                   htmlwidgets_1.6.4          
# [40] spatstat.geom_3.2-4         stats4_4.1.3                ggnewscale_0.5.0           
# [43] backports_1.5.0             biomaRt_2.48.3              deldir_1.0-6               
# [46] MatrixGenerics_1.4.3        vctrs_0.6.5                 SingleCellExperiment_1.14.1
# [49] Biobase_2.52.0              ROCR_1.0-11                 abind_1.4-5                
# [52] cachem_1.1.0                withr_3.0.2                 ggforce_0.4.2              
# [55] progressr_0.10.0            sctransform_0.3.5           treeio_1.16.2              
# [58] prettyunits_1.2.0           goftest_1.2-2               DOSE_3.18.3                
# [61] cluster_2.1.2               ape_5.6                     dotCall64_1.2              
# [64] lazyeval_0.2.2              crayon_1.5.3                spatstat.explore_3.2-1     
# [67] labeling_0.4.3              pkgconfig_2.0.3             units_0.8-5                
# [70] tweenr_2.0.3                GenomeInfoDb_1.30.1         nlme_3.1-152               
# [73] rlang_1.1.4                 globals_0.15.0              lifecycle_1.0.4            
# [76] miniUI_0.1.1.1              downloader_0.4              filelock_1.0.2             
# [79] BiocFileCache_2.0.0         modelr_0.1.8                cellranger_1.1.0           
# [82] polyclip_1.10-7             matrixStats_0.62.0          lmtest_0.9-40              
# [85] tiff_0.1-11                 aplot_0.1.1                 Matrix_1.6-4               
# [88] carData_3.0-4               Rhdf5lib_1.14.2             zoo_1.8-10                 
# [91] reprex_2.0.1                ggridges_0.5.6              GlobalOptions_0.1.2        
# [94] png_0.1-8                   viridisLite_0.4.2           rjson_0.2.20               
# [97] clisymbols_1.2.0            bitops_1.0-9                KernSmooth_2.23-20         
# [100] spam_2.11-0                 rhdf5filters_1.4.0          Biostrings_2.60.2          
# [103] anndata_0.7.5.6             blob_1.2.2                  shape_1.4.6                
# [106] qvalue_2.24.0               slingshot_2.2.1             parallelly_1.31.1          
# [109] spatstat.random_3.1-5       R.cache_0.16.0              gridGraphics_0.5-1         
# [112] jpeg_0.1-9                  rstatix_0.7.0               S4Vectors_0.30.2           
# [115] ggsignif_0.6.2              scales_1.3.0                memoise_2.0.1              
# [118] magrittr_2.0.3              plyr_1.8.7                  ica_1.0-2                  
# [121] zlibbioc_1.38.0             scatterpie_0.1.7            compiler_4.1.3             
# [124] RColorBrewer_1.1-3          clue_0.3-60                 fitdistrplus_1.1-5         
# [127] Rsamtools_2.8.0             cli_3.6.3                   XVector_0.32.0             
# [130] listenv_0.8.0               patchwork_1.3.0             pbapply_1.5-0              
# [133] MASS_7.3-54                 tidyselect_1.2.1            stringi_1.8.4              
# [136] textshaping_0.3.6           GOSemSim_2.18.1             locfit_1.5-9.4             
# [139] ggrepel_0.9.1               grid_4.1.3                  fastmatch_1.1-3            
# [142] tools_4.1.3                 timechange_0.3.0            future.apply_1.8.1         
# [145] parallel_4.1.3              rio_0.5.27                  circlize_0.4.16            
# [148] rstudioapi_0.17.1           foreach_1.5.1               foreign_0.8-81             
# [151] gridExtra_2.3               farver_2.1.2                Rtsne_0.16                 
# [154] ggraph_2.0.5                proxyC_0.4.1                digest_0.6.37              
# [157] shiny_1.10.0                Rcpp_1.0.13-1               GenomicRanges_1.44.0       
# [160] car_3.0-11                  broom_1.0.7                 later_1.4.1                
# [163] RcppAnnoy_0.0.22            httr_1.4.7                  AnnotationDbi_1.54.1       
# [166] ComplexHeatmap_2.21.2       colorspace_2.1-1            rvest_1.0.1                
# [169] XML_3.99-0.7                fs_1.6.5                    tensor_1.5                 
# [172] IRanges_2.26.0              splines_4.1.3               yulab.utils_0.1.8          
# [175] uwot_0.1.16                 RcppRoll_0.3.1              tidytree_0.3.6             
# [178] spatstat.utils_3.1-1        graphlayouts_0.7.1          sp_2.1-4                   
# [181] ggplotify_0.1.2             systemfonts_1.1.0           plotly_4.10.4              
# [184] xtable_1.8-4                ggtree_3.0.4                jsonlite_1.8.9             
# [187] tidygraph_1.2.0             ggfun_0.0.4                 R6_2.5.1                   
# [190] pillar_1.10.0               htmltools_0.5.8.1           mime_0.12                  
# [193] clusterProfiler_4.2.2       glue_1.8.0                  fastmap_1.2.0              
# [196] BiocParallel_1.26.2         fftwtools_0.9-11            codetools_0.2-18           
# [199] fgsea_1.18.0                utf8_1.2.4                  Signac_1.14.0              
# [202] lattice_0.20-44             spatstat.sparse_3.0-2       curl_6.0.1                 
# [205] leiden_0.3.9                magick_2.8.5                GO.db_3.13.0               
# [208] zip_2.2.0                   openxlsx_4.2.4              survival_3.2-12            
# [211] munsell_0.5.1               DO.db_2.9                   GetoptLong_1.0.5           
# [214] rhdf5_2.36.0                GenomeInfoDbData_1.2.6      iterators_1.0.13           
# [217] HDF5Array_1.20.0            haven_2.4.3                 reshape2_1.4.4             
# [220] gtable_0.3.6               


# end ---------------------------------------------------------------------


