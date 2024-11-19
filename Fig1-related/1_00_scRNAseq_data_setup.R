# Takahide Nejo, MD, PhD
# Okada Lab, Department of Neurological Surgery, UC San Francisco


# note ---------------------------------------------------------------------

# Nejo T, et al., Glioma-neuronal circuit remodeling induces regional immunosuppression
# Preparation of the single-cell RNA-seq Seurat object, as described in Krishna S, 2023 Nature (10.1038/s41586-023-06036-1).
# Please note that we relied on the original annotation directly provided by the authors of the paper (Krishna S, 2023, Nature).
# Using different versions of packages, such as Seurat, may lead to varying outcomes, which we have not thoroughly validated.


# relevant description -------------------------------------------------------------

# Krishna S, 2023 Nature (10.1038/s41586-023-06036-1)

# One hundred base pair paired-end reads were sequenced on the Illumina NovaSeq 6000 system at the Center for Advanced Technology at the University of California San Francisco, and the resulting FASTQ files were processed using the CellRanger analysis suite (v.3.0.2; https://github.com/10XGenomics/cellranger) for alignment to the hg38 reference genome, identification of empty droplets, and determination of the count threshold for further analysis. A cell quality filter of greater than 500 features but fewer than 10,000 features per cell, and less than 20% of read counts attributed to mitochondrial genes, was used. Single-cell UMI count data were preprocessed in Seurat (v.3.0.1)74,75 using the sctransform workflow76, with scaling based on the regression of UMI count and the percentage of reads attributed to mitochondrial genes per cell. Dimensionality reduction was performed using principal component analysis and then principal component loadings were corrected for batch effects using Harmony77. Uniform manifold approximation and projection was performed on the reduced data with a minimum distance metric of 0.4 and Louvain clustering was performed using a resolution of 0.2. Marker selection was performed in Seurat using a minimum difference in the fraction of detection of 0.5 and a minimum log-transformed fold change of 0.5. We assessed the single-cell transcriptome from 6,666 HFC-region cells and 7,065 LFC-region cells (Supplementary Table 4).


# rm all ------------------------------------------------------------------

rm(list = ls(all.names = T))


# packages ----------------------------------------------------------------

library(tidyverse)
library(tidylog)
library(Seurat)
library(harmony)


# check -------------------------------------------------------------------

packageVersion("Seurat")

packageVersion("harmony")


# dir ---------------------------------------------------------------------

dir.1 <- "/okadalab/data1/tnejo/neuro_immune_proj/dl/GSE223065/"
dir.2 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/1_00_scRNAseq_data_setup/"
if(dir.exists(dir.2) == F){dir.create(dir.2)}


# load files --------------------------------------------------------------

setwd(dir.1)
list.files() %>% print()
# [1] "GSE223065_RAW.tar" "GSM6939133_HFC1"   "GSM6939134_LFC1"   "GSM6939135_HFC2"  
# [5] "GSM6939136_LFC2"   "GSM6939137_HFC3"   "GSM6939138_LFC3"  

START.TIME <- Sys.time() 

list.seurat.obj <- list()
for(i in 1:6){
  print(i)
  dir.i <- c("GSM6939133_HFC1", "GSM6939134_LFC1", "GSM6939135_HFC2", "GSM6939136_LFC2", "GSM6939137_HFC3", "GSM6939138_LFC3")[i]
  data.i <- Read10X(dir.i)
  list.seurat.obj[[i]] <- CreateSeuratObject(data.i, strip.suffix = T)
  
  list.seurat.obj[[i]][["orig.ident"]] <- c("hfc1", "lfc1", "hfc2", "lfc2", "hfc3", "lfc3")[i]
}

FINISH.TIME <- Sys.time() 

FINISH.TIME - START.TIME 
# Time difference of 31.75395 secs


# check -------------------------------------------------------------------

list.seurat.obj %>% length()
# [1] 6

list.seurat.obj %>% print()
{
# [[1]]
# An object of class Seurat 
# 33538 features across 1314 samples within 1 assay 
# Active assay: RNA (33538 features, 0 variable features)
# 
# [[2]]
# An object of class Seurat 
# 33538 features across 6673 samples within 1 assay 
# Active assay: RNA (33538 features, 0 variable features)
# 
# [[3]]
# An object of class Seurat 
# 33538 features across 3056 samples within 1 assay 
# Active assay: RNA (33538 features, 0 variable features)
# 
# [[4]]
# An object of class Seurat 
# 33538 features across 745 samples within 1 assay 
# Active assay: RNA (33538 features, 0 variable features)
# 
# [[5]]
# An object of class Seurat 
# 33538 features across 4733 samples within 1 assay 
# Active assay: RNA (33538 features, 0 variable features)
# 
# [[6]]
# An object of class Seurat 
# 33538 features across 4046 samples within 1 assay 
# Active assay: RNA (33538 features, 0 variable features)
}

list.seurat.obj[[1]]@meta.data %>% 
  head()


# merge -------------------------------------------------------------------

START.TIME <- Sys.time() 

seurat.obj <- merge(x = list.seurat.obj[[1]], 
                    y = list(list.seurat.obj[[2]], list.seurat.obj[[3]], list.seurat.obj[[4]], list.seurat.obj[[5]], list.seurat.obj[[6]]))

FINISH.TIME <- Sys.time() 

FINISH.TIME - START.TIME 
# Time difference of 6.908758 secs


# check -------------------------------------------------------------------

seurat.obj %>% print()

seurat.obj@meta.data %>% 
  as.data.frame() %>% head()

seurat.obj@meta.data %>% 
  as.data.frame() %>% tail()


# percent.mt ---------------------------------------------------------------

seurat.obj[["percent.mt"]] <- seurat.obj %>% 
  PercentageFeatureSet(pattern = "^MT-")


# check -------------------------------------------------------------------

seurat.obj@meta.data %>% head()

seurat.obj %>% dim()


# filtering ---------------------------------------------------------------

seurat.obj <- seurat.obj %>% 
  subset(subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 20)


# check -------------------------------------------------------------------

seurat.obj %>% dim()

seurat.obj@meta.data %>% 
  head()


seurat.obj@meta.data %>% 
  tail()


# SCTransform -------------------------------------------------------------

seurat.obj.orig <- seurat.obj

START.TIME <- Sys.time() 

seurat.obj <- seurat.obj %>% 
  SCTransform(vars.to.regress = c("nCount_RNA", "percent.mt"))

FINISH.TIME <- Sys.time() 

FINISH.TIME - START.TIME 
# Time difference of 2.325124 mins


# check -------------------------------------------------------------------

seurat.obj %>% print()


# RunPCA -------------------------------------------------------------
 
seurat.obj.orig <- seurat.obj
 
START.TIME <- Sys.time() 
 
seurat.obj <- seurat.obj %>% 
  RunPCA(npcs = 100)

FINISH.TIME <- Sys.time() 
FINISH.TIME - START.TIME 
# Time difference of 27.18047 secs


# check -------------------------------------------------------------------

seurat.obj %>% print()

seurat.obj@meta.data %>% head()


# edit metadata -----------------------------------------------------------

metadata.old <- seurat.obj@meta.data

metadata.new <- metadata.old %>% 
  as_tibble() %>% 
  mutate(batch = ifelse(grepl("1", orig.ident), "one", 
                        ifelse(grepl("2", orig.ident), "two", "three"))) %>% 
  mutate(tissue = ifelse(grepl("hfc", orig.ident), "hfc", "lfc")) %>% 
  as.data.frame()

rownames(metadata.new) <- rownames(metadata.old)

seurat.obj@meta.data <- metadata.new


# check -------------------------------------------------------------------

seurat.obj@meta.data %>% head()


# Harmony -----------------------------------------------------------------

seurat.obj.orig <- seurat.obj

START.TIME <- Sys.time() 

seurat.obj <- seurat.obj %>% 
  RunHarmony(group.by.vars = "batch", assay.use = "SCT")
 
FINISH.TIME <- Sys.time() 
FINISH.TIME - START.TIME 

# Time difference of 1.194733 mins


# check -------------------------------------------------------------------

seurat.obj %>% print()


# RunUMAP / FindNeighbors / FindClusters --------------------------------------------------------------

seurat.obj.orig <- seurat.obj

START.TIME <- Sys.time() 

seurat.obj <- seurat.obj %>% 
  RunUMAP(dims = 1:30, reduction = "harmony", min.dist = 0.4) %>% 
  FindNeighbors(dims = 1:30, reduction = "harmony") %>% 
  FindClusters(resolution = 0.3) #

FINISH.TIME <- Sys.time() 
FINISH.TIME - START.TIME 
# Time difference of 20.66064 secs


# check -------------------------------------------------------------------

seurat.obj %>% print()

seurat.obj@meta.data %>% head()


# cell annotation --------------------------------------------------------------------

annot <- tibble(
  seurat_clusters = 0:13, 
  clust = c("Myeloid Cells", rep("Tumor Cells", 5), "Tumor Cells (Dividing)", 
            "Neurons", "Tumor Cells", "Lymphoid Cells", "Endothelial Cells", "Pericytes", "Tumor Cells", 
            "Astrocytes")
)
annot$seurat_clusters <- factor(annot$seurat_clusters)



# update the metadata -----------------------------------------------------

seurat.obj.orig <- seurat.obj

metadata.old <- seurat.obj@meta.data

metadata.new <- metadata.old %>% 
  as_tibble() %>% 
  left_join(annot, by = "seurat_clusters") %>% 
  as.data.frame()

rownames(metadata.new) <- rownames(metadata.old)

seurat.obj@meta.data <- metadata.new


# check ---------------------------------------------------------------------

seurat.obj %>% print()

seurat.obj@meta.data %>% head()


# out ---------------------------------------------------------------------

setwd(dir.2)
out.f <- "gbm_neuroimmune.rds"
saveRDS(seurat.obj, out.f)


# si ----------------------------------------------------------------------

Sys.time()


sessionInfo()
# R version 4.1.2 (2021-11-01)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# 
# Matrix products: default
# BLAS:   /software/c4/cbi/software/R-4.1.2-gcc8/lib64/R/lib/libRblas.so
# LAPACK: /software/c4/cbi/software/R-4.1.2-gcc8/lib64/R/lib/libRlapack.so
# 
# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] RColorBrewer_1.1-2 SeuratObject_4.0.2 Seurat_4.0.3       forcats_0.5.1     
#  [5] stringr_1.4.0      dplyr_1.0.7        purrr_0.3.4        readr_2.0.1       
#  [9] tidyr_1.1.3        tibble_3.1.3       ggplot2_3.3.5      tidyverse_1.3.1   
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15            colorspace_2.0-2      deldir_0.2-10        
#   [4] ellipsis_0.3.2        ggridges_0.5.3        fs_1.5.0             
#   [7] spatstat.data_2.1-0   rstudioapi_0.13       farver_2.1.0         
#  [10] leiden_0.3.9          listenv_0.8.0         ggrepel_0.9.1        
#  [13] fansi_0.5.0           lubridate_1.7.10      xml2_1.3.2           
#  [16] codetools_0.2-18      splines_4.1.2         polyclip_1.10-0      
#  [19] jsonlite_1.7.2        broom_0.7.9           ica_1.0-2            
#  [22] cluster_2.1.2         dbplyr_2.1.1          png_0.1-7            
#  [25] uwot_0.1.10           spatstat.sparse_2.0-0 shiny_1.6.0          
#  [28] sctransform_0.3.2     compiler_4.1.2        httr_1.4.2           
#  [31] backports_1.2.1       assertthat_0.2.1      Matrix_1.3-4         
#  [34] fastmap_1.1.0         lazyeval_0.2.2        cli_3.0.1            
#  [37] later_1.2.0           htmltools_0.5.1.1     tools_4.1.2          
#  [40] igraph_1.2.6          gtable_0.3.0          glue_1.4.2           
#  [43] RANN_2.6.1            reshape2_1.4.4        Rcpp_1.0.7           
#  [46] scattermore_0.7       cellranger_1.1.0      vctrs_0.3.8          
#  [49] nlme_3.1-152          lmtest_0.9-38         globals_0.14.0       
#  [52] rvest_1.0.1           mime_0.11             miniUI_0.1.1.1       
#  [55] lifecycle_1.0.0       irlba_2.3.3           goftest_1.2-2        
#  [58] future_1.21.0         MASS_7.3-54           zoo_1.8-9            
#  [61] scales_1.1.1          spatstat.core_2.3-0   spatstat.utils_2.2-0 
#  [64] hms_1.1.0             promises_1.2.0.1      parallel_4.1.2       
#  [67] gridExtra_2.3         reticulate_1.20       pbapply_1.4-3        
#  [70] rpart_4.1-15          stringi_1.7.3         rlang_0.4.11         
#  [73] pkgconfig_2.0.3       matrixStats_0.60.0    lattice_0.20-44      
#  [76] tensor_1.5            ROCR_1.0-11           labeling_0.4.2       
#  [79] patchwork_1.1.1       htmlwidgets_1.5.3     cowplot_1.1.1        
#  [82] tidyselect_1.1.1      parallelly_1.27.0     RcppAnnoy_0.0.19     
#  [85] plyr_1.8.6            magrittr_2.0.1        R6_2.5.0             
#  [88] generics_0.1.0        DBI_1.1.1             mgcv_1.8-36          
#  [91] pillar_1.6.2          haven_2.4.3           withr_2.4.2          
#  [94] fitdistrplus_1.1-5    abind_1.4-5           survival_3.2-12      
#  [97] future.apply_1.8.1    modelr_0.1.8          crayon_1.4.1         
# [100] KernSmooth_2.23-20    utf8_1.2.2            spatstat.geom_2.2-2  
# [103] plotly_4.9.4.1        tzdb_0.1.2            grid_4.1.2           
# [106] readxl_1.3.1          data.table_1.14.0     reprex_2.0.1         
# [109] digest_0.6.27         xtable_1.8-4          httpuv_1.6.1         
# [112] munsell_0.5.0         viridisLite_0.4.0    


# end ---------------------------------------------------------------------


