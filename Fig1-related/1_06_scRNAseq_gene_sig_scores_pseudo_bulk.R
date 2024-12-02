# Takahide Nejo, MD, PhD
# Okada Lab, Department of Neurological Surgery, UC San Francisco


# note ---------------------------------------------------------------------

# Nejo T, et al., Glioma-neuronal circuit remodeling induces regional immunosuppression

# Using different versions of packages, such as Seurat, may lead to varying outcomes, which we have not thoroughly validated.

# The RDS object "001_gbm_neuroimmune_v3_subset_list.rds" is used as inputs. For details, please refer to the previous step "1_02_scRNAseq_dge.R". 


# rm all ------------------------------------------------------------------

rm(list = ls(all.names = T))


# packages ----------------------------------------------------------------

library(tidyverse)
library(tidylog)
library(Seurat)
library(msigdbr)
library(DESeq2)
library(SingleCellExperiment)
library(Matrix.utils)
library(data.table)
library(GSVA) 
library(GSA) # install.packages("GSA")
library(AUCell) # BiocManager::install("AUCell")
library(ggsci)
library(ggpubr)


# check -------------------------------------------------------------------

packageVersion("Seurat")
# [1] ‘4.3.0.1’

packageVersion("msigdbr")
# [1] ‘7.4.1’


# dir ---------------------------------------------------------------------

dir.1 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/1_02_scRNAseq_dge/"
dir.2 <- "/okadalab/data1/tnejo/neuro_immune_proj/for_github/out/1_06_scRNAseq_gene_sig_scores_pseudo_bulk/"
if(dir.exists(dir.2) == F){dir.create(dir.2)}


# load files --------------------------------------------------------------

START.TIME <- Sys.time() 

setwd(dir.1) 
in.f <- "001_gbm_neuroimmune_v3_subset_list.rds"
list.seurat.obj <- readRDS(in.f)

FINISH.TIME <- Sys.time() 

FINISH.TIME - START.TIME 
# Time difference of 43.9747 secs


# check -------------------------------------------------------------------

list.seurat.obj %>% length()
# [1] 4

list.seurat.obj %>% names()
# [1] "Tumor Cells"    "Myeloid Cells"  "Lymphoid Cells" "Astrocytes"    

list.seurat.obj[[1]]
# An object of class Seurat 
# 47200 features across 8537 samples within 2 assays 
# Active assay: RNA (24173 features, 0 variable features)
# 1 other assay present: SCT
# 3 dimensional reductions calculated: pca, harmony, umap


# Prep the count matrix, metadata, and sce, for each cell type -----------------------------------------------------

list.sce <- list()
for(i in 1:3){
  print(i)
  
  seurat.obj.i <- list.seurat.obj[[i]]
  
  # extract raw counts and metadata to create SingleCellExperiment object
  
  counts.i <- seurat.obj.i@assays$RNA@counts 
  
  metadata.i.orig <- seurat.obj.i@meta.data
  
  metadata.i <- metadata.i.orig %>% 
    as_tibble() %>% 
    mutate(sample.id = gsub("1", "_1", gsub("2", "_2", gsub("3", "_3", orig.ident)))) %>% 
    dplyr::select(clust, sample.id) %>% 
    as.data.frame()
  
  rownames(metadata.i) <- rownames(metadata.i.orig)
  metadata.i$sample.id <- factor(metadata.i$sample.id)
  
  sce.i <- SingleCellExperiment(
    assays = list(counts = counts.i), 
    colData = metadata.i
  )
  list.sce[[i]] <- sce.i
}

names(list.sce) <- c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells")


# check -------------------------------------------------------------------

for(i in 1:3){
  print(i)
  names(list.sce)[i] %>% print()
  list.sce[[i]] %>% dim() %>% print()
}

# [1] 1
# [1] "Tumor Cells"
# [1] 24173  8537
# [1] 2
# [1] "Myeloid Cells"
# [1] 24173  3775
# [1] 3
# [1] "Lymphoid Cells"
# [1] 24173   395


# aggregate counts per sample -----------------

list.aggr.counts <- list()
for(i in 1:3){
  print(i)
  sce.i <- list.sce[[i]]
  groups.i <- colData(sce.i)[, "sample.id"]
  
  aggr.counts.i <- aggregate.Matrix(
    t(counts(sce.i)), 
    groupings = groups.i, 
    fun = "sum") 
  
  # transpose aggregated matrix to have genes as rows and samples as columns
  aggr.counts.i <- t(aggr.counts.i)
  
  list.aggr.counts[[i]] <- aggr.counts.i
}

names(list.aggr.counts) <- c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells")


# check -------------------------------------------------------------------

for(i in 1:3){
  print(i)
  names(list.aggr.counts)[i] %>% print()
  list.aggr.counts[[i]] %>% class() %>% print()
  list.aggr.counts[[i]] %>% dim() %>% print()
  list.aggr.counts[[i]][1:6, 1:6] %>% print()
  list.aggr.counts[[i]] %>% str() %>% print()
}

{
  # [1] 1
  # [1] "Tumor Cells"
  # [1] "dgCMatrix"
  # attr(,"package")
  # [1] "Matrix"
  # [1] 24173     6
  # 6 x 6 sparse Matrix of class "dgCMatrix"
  # hfc_1 hfc_2 hfc_3 lfc_1 lfc_2 lfc_3
  # AL627309.1     7    15    22    10     5    23
  # AL627309.4     .     .    17     1     .    21
  # AL669831.2     1     .     .     .     .     .
  # AL669831.5    80   628   348   117   186   289
  # FAM87B        10    13    30    16     6    10
  # LINC00115     56   112   224    58    40   134
  # Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
  # ..@ i       : int [1:122222] 0 2 3 4 5 6 7 8 9 10 ...
  # ..@ p       : int [1:7] 0 19351 40890 61684 82284 101144 122222
  # ..@ Dim     : int [1:2] 24173 6
  # ..@ Dimnames:List of 2
  # .. ..$ : chr [1:24173] "AL627309.1" "AL627309.4" "AL669831.2" "AL669831.5" ...
  # .. ..$ : chr [1:6] "hfc_1" "hfc_2" "hfc_3" "lfc_1" ...
  # ..@ x       : num [1:122222] 7 1 80 10 56 43 3 5 386 33 ...
  # ..@ factors : list()
  # NULL
  # [1] 2
  # [1] "Myeloid Cells"
  # [1] "dgCMatrix"
  # attr(,"package")
  # [1] "Matrix"
  # [1] 24173     6
  # 6 x 6 sparse Matrix of class "dgCMatrix"
  # hfc_1 hfc_2 hfc_3 lfc_1 lfc_2 lfc_3
  # AL627309.1     3     2     1    11     4     1
  # AL627309.4     3     .     .    10     .     1
  # AL669831.2     2     .     .     .     .     .
  # AL669831.5    55    25    17   166    13    38
  # FAM87B         5     1     3    58     .     7
  # LINC00115     23     9    10   131     3    12
  # Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
  # ..@ i       : int [1:102663] 0 1 2 3 4 5 6 8 9 10 ...
  # ..@ p       : int [1:7] 0 16847 33924 50120 70407 84958 102663
  # ..@ Dim     : int [1:2] 24173 6
  # ..@ Dimnames:List of 2
  # .. ..$ : chr [1:24173] "AL627309.1" "AL627309.4" "AL669831.2" "AL669831.5" ...
  # .. ..$ : chr [1:6] "hfc_1" "hfc_2" "hfc_3" "lfc_1" ...
  # ..@ x       : num [1:102663] 3 3 2 55 5 23 19 3 95 12 ...
  # ..@ factors : list()
  # NULL
  # [1] 3
  # [1] "Lymphoid Cells"
  # [1] "dgCMatrix"
  # attr(,"package")
  # [1] "Matrix"
  # [1] 24173     6
  # 6 x 6 sparse Matrix of class "dgCMatrix"
  # hfc_1 hfc_2 hfc_3 lfc_1 lfc_2 lfc_3
  # AL627309.1     .     .     .     1     .     .
  # AL627309.4     .     .     .     .     .     .
  # AL669831.2     .     .     .     .     .     .
  # AL669831.5     5     .     .    12     1     1
  # FAM87B         .     .     .     .     .     .
  # LINC00115      4     .     .     4     .     .
  # Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
  # ..@ i       : int [1:59639] 3 5 6 9 10 13 14 17 18 19 ...
  # ..@ p       : int [1:7] 0 11933 22221 29326 44645 50734 59639
  # ..@ Dim     : int [1:2] 24173 6
  # ..@ Dimnames:List of 2
  # .. ..$ : chr [1:24173] "AL627309.1" "AL627309.4" "AL669831.2" "AL669831.5" ...
  # .. ..$ : chr [1:6] "hfc_1" "hfc_2" "hfc_3" "lfc_1" ...
  # ..@ x       : num [1:59639] 5 4 1 17 2 56 1 9 7 27 ...
  # ..@ factors : list()
  # NULL
}



# prepare new metadata for each cell type ---------------------------------

list.metadata <- list()
for (i in 1:length(list.aggr.counts)) {
  print(i)
  
  aggr.counts.i <- list.aggr.counts[[i]]
  
  metadata.i <- data.frame(sample.id = colnames(aggr.counts.i))
  
  # number of cells per sample and cluster
  cell.counts <- table(colData(list.sce[[i]])$sample.id)

  # remove samples with zero cell contributing to the cluster
  cell.counts <- cell.counts[cell.counts > 0]
  
  # append cell_counts to data frame
  metadata.i$cell.count <- cell.counts

  metadata.i$tissue <- rep(c("hfc", "lfc"), each = 3)
  metadata.i$pt.id <- rep(c("pt_1", "pt_2", "pt_3"), 2)
  
  # update rownames of metadata to match colnames of count matrix, as needed later for DE
  rownames(metadata.i) <- metadata.i$sample.id
  
  list.metadata[[i]] <- metadata.i
}

names(list.metadata) <- c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells")


# check -------------------------------------------------------------------

for(i in 1:3){
  print(i)
  names(list.metadata)[i] %>% print()
  list.metadata[[i]] %>% dim() %>% print()
  list.metadata[[i]] %>% head() %>% print()
  list.metadata[[i]] %>% str() %>% print()
}

{
  # [1] 1
  # [1] "Tumor Cells"
  # [1] 6 4
  # sample.id cell.count tissue pt.id
  # hfc_1     hfc_1        525    hfc  pt_1
  # hfc_2     hfc_2       2143    hfc  pt_2
  # hfc_3     hfc_3       2657    hfc  pt_3
  # lfc_1     lfc_1       1250    lfc  pt_1
  # lfc_2     lfc_2        406    lfc  pt_2
  # lfc_3     lfc_3       1556    lfc  pt_3
  # 'data.frame':	6 obs. of  4 variables:
  #   $ sample.id : chr  "hfc_1" "hfc_2" "hfc_3" "lfc_1" ...
  # $ cell.count: 'table' int [1:6(1d)] 525 2143 2657 1250 406 1556
  # $ tissue    : chr  "hfc" "hfc" "hfc" "lfc" ...
  # $ pt.id     : chr  "pt_1" "pt_2" "pt_3" "pt_1" ...
  # NULL
  # [1] 2
  # [1] "Myeloid Cells"
  # [1] 6 4
  # sample.id cell.count tissue pt.id
  # hfc_1     hfc_1        309    hfc  pt_1
  # hfc_2     hfc_2        294    hfc  pt_2
  # hfc_3     hfc_3        248    hfc  pt_3
  # lfc_1     lfc_1       2338    lfc  pt_1
  # lfc_2     lfc_2         81    lfc  pt_2
  # lfc_3     lfc_3        505    lfc  pt_3
  # 'data.frame':	6 obs. of  4 variables:
  #   $ sample.id : chr  "hfc_1" "hfc_2" "hfc_3" "lfc_1" ...
  # $ cell.count: 'table' int [1:6(1d)] 309 294 248 2338 81 505
  # $ tissue    : chr  "hfc" "hfc" "hfc" "lfc" ...
  # $ pt.id     : chr  "pt_1" "pt_2" "pt_3" "pt_1" ...
  # NULL
  # [1] 3
  # [1] "Lymphoid Cells"
  # [1] 6 4
  # sample.id cell.count tissue pt.id
  # hfc_1     hfc_1         64    hfc  pt_1
  # hfc_2     hfc_2         10    hfc  pt_2
  # hfc_3     hfc_3         12    hfc  pt_3
  # lfc_1     lfc_1        284    lfc  pt_1
  # lfc_2     lfc_2          7    lfc  pt_2
  # lfc_3     lfc_3         18    lfc  pt_3
  # 'data.frame':	6 obs. of  4 variables:
  #   $ sample.id : chr  "hfc_1" "hfc_2" "hfc_3" "lfc_1" ...
  # $ cell.count: 'table' int [1:6(1d)] 64 10 12 284 7 18
  # $ tissue    : chr  "hfc" "hfc" "hfc" "lfc" ...
  # $ pt.id     : chr  "pt_1" "pt_2" "pt_3" "pt_1" ...
  # NULL
}


# saveRDS -----------------------------------------------------------------

START.TIME <- Sys.time() 

setwd(dir.2)
out.f <- "001_list_sce_obj.rds"
saveRDS(list.sce, out.f)

setwd(dir.2)
out.f <- "002_list_aggr_counts.rds"
saveRDS(list.aggr.counts, out.f)

setwd(dir.2)
out.f <- "003_list_metadata.rds"
saveRDS(list.metadata, out.f)

FINISH.TIME <- Sys.time() 

FINISH.TIME - START.TIME 
# Time difference of 30.95961 secs


# pathways of interest --------------------------------------------

imm.pathways <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE", 
  "HALLMARK_INTERFERON_GAMMA_RESPONSE", 
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)


# prep the genesets of interest ---------------------------------------------------------

test_gs <- msigdbr(
  species = "Homo sapiens", 
  category = "H"
) %>% 
  dplyr::filter(gs_name %in% imm.pathways) %>% 
  dplyr::select(gs_name, human_gene_symbol) %>% 
  group_by(gs_name) %>% 
  summarise(human_gene_symbol = list(human_gene_symbol))


# check -------------------------------------------------------------------

test_gs %>% dim()
# [1] 3 2

test_gs %>% head() %>% print()
# # A tibble: 3 × 2
# gs_name                            human_gene_symbol
# <chr>                              <list>           
#   1 HALLMARK_INFLAMMATORY_RESPONSE     <chr [222]>      
#   2 HALLMARK_INTERFERON_GAMMA_RESPONSE <chr [284]>      
#   3 HALLMARK_TNFA_SIGNALING_VIA_NFKB   <chr [227]>    



# prep functions ----------------------------------------------------------

# prep for JASMINE --------------------------------------------------------

# ref: https://github.com/NNoureen/JASMINE/blob/main/JASMINE_V1_11October2021.r


RankCalculation <- function(x, genes){
  
  subdata = x[x != 0]                                                                      ### Removing Dropouts from single cell
  DataRanksUpdated=rank(subdata)                                                         ### Calculating ranks of each signature gene per cell
  DataRanksSigGenes = DataRanksUpdated[which(names(DataRanksUpdated) %in% genes)]        ### Shortling rank vector for signature genes 
  CumSum = ifelse(length(DataRanksSigGenes), mean(DataRanksSigGenes, na.rm = TRUE), 0 )     ### Calculating Mean of ranks for signature genes
  FinalRawRank = CumSum/length(subdata)                                                  ### Normalizing Means by total coverage
  return(FinalRawRank)                                                
}			

#### Function2:- Calculating enrichment of signature genes across each cell 	(using odds ratio)

ORCalculation <- function(data, genes){
  GE = data[which(rownames(data) %in% genes), ]                                          ### Subsetting data for signature genes
  NGE = data[-which(rownames(data) %in% genes), ]                                        ### Subsetting data for non-signature genes
  SigGenesExp = apply(GE, 2, function(x) length(x[x != 0]))                                 ### Calculating Number of expressed Signature Genes per cell
  NSigGenesExp = apply(NGE, 2, function(x) length(x[x != 0]))                               ### Calculating Number of expressed Non-Signature Genes per cell
  SigGenesNE = nrow(GE) - SigGenesExp                                                   ### Calculating Number of Not expressed Signature Genes per cell
  SigGenesNE = replace(SigGenesNE, SigGenesNE == 0, 1)									  ### Replacing Zero's with 1
  NSigGenesExp = replace(NSigGenesExp, NSigGenesExp == 0, 1)                                ### Replacing Zero's with 1
  NSigGenesNE = nrow(data) - (NSigGenesExp + SigGenesExp)                               ### Calculating Number of Not expressed Non-Signature Genes per cell
  NSigGenesNE = NSigGenesNE - SigGenesNE
  OR = (SigGenesExp * NSigGenesNE) / (SigGenesNE * NSigGenesExp)                         ### Calculating Enrichment (Odds Ratio)
  return(OR)
}


#### Function3:- Calculating enrichment of signature genes across each cell (using likelihood ratio)

LikelihoodCalculation <- function(data, genes){
  GE = data[which(rownames(data) %in% genes), ]
  NGE = data[-which(rownames(data) %in% genes), ]
  SigGenesExp = apply(GE, 2, function(x) length(x[x != 0]))
  NSigGenesExp = apply(NGE, 2, function(x) length(x[x != 0]))
  SigGenesNE = nrow(GE) - SigGenesExp
  SigGenesNE = replace(SigGenesNE, SigGenesNE == 0, 1)			
  NSigGenesExp = replace(NSigGenesExp, NSigGenesExp == 0, 1)
  NSigGenesNE = nrow(data) - (NSigGenesExp + SigGenesExp)
  NSigGenesNE = NSigGenesNE - SigGenesNE
  LR1 = SigGenesExp*(NSigGenesExp + NSigGenesNE)
  LR2 = NSigGenesExp * (SigGenesExp + SigGenesNE)
  LR = LR1/LR2
  return(LR)
}	


### Function 5:- Signature Scoring via JASMINE merging Means and Enrichment

JASMINE <- function(data, genes, method)
{
  idx = match(genes, rownames(data))                                                 
  idx = idx[!is.na(idx)]
  if(length(idx) > 1){
    RM = apply(data, 2, function(x) RankCalculation(x, genes))                              ### Mean RankCalculation for single cell data matrix
    RM = NormalizationJAS(RM)                                                            ### Normalizing Mean Ranks
    
    if(method == "oddsratio"){
      OR = ORCalculation(data, genes)			                                             ### Signature Enrichment Calculation for single cell data matrix (OR)
      OR = NormalizationJAS(OR)															 ### Normalizing Enrichment Scores (OR)
      JAS_Scores = (RM + OR)/2
    }else if(method == "likelihood"){
      
      LR = LikelihoodCalculation(data, genes)			                                     ### Signature Enrichment Calculation for single cell data matrix  (LR)
      LR = NormalizationJAS(LR)															 ### Normalizing Enrichment Scores (LR)
      JAS_Scores = (RM + LR)/2
    }
    FinalScores = data.frame(names(RM), JAS_Scores)                                       ### JASMINE scores
    colnames(FinalScores)[1] = 'SampleID'
    return(FinalScores)
  }
}



###  Function 4:- Scalar [0, 1] Normalization of Means and Enrichment across set of cells

NormalizationJAS <- function(JAS_Scores)
{
  JAS_Scores = (JAS_Scores - min(JAS_Scores))/(max(JAS_Scores) - min(JAS_Scores))
  return(JAS_Scores)
}


# 01. ssGSEA / 02. GSVA / 03. AUCell / 04. JASMINE (LR) --------------------------------------------------------------

# ref: https://github.com/NNoureen/BenchmarkingProtocol/blob/main/ssGSEA.r

res.merged <- NULL
for(i in 1:3){
  print(i)
  name.i <-  c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells")[i]
  
  aggr.counts.i <- list.aggr.counts[[i]]
  
  # exprMatrix %>% dim()
  # # [1] 24182     6
  
  res.i <- NULL
  for(j in 1:3){
    print(j)
    
    # prep the gene list
    
    pathway.j <- imm.pathways[j]
    
    gene.list.j <- test_gs %>% 
      dplyr::filter(gs_name == pathway.j) %>% 
      pull(human_gene_symbol)
    names(gene.list.j)[1] <- pathway.j
    
    # run the analyses
    
    # ssGSEA
    res.ssgsea <- gsva(aggr.counts.i, gene.list.j, method = "ssgsea")
    
    # GSVA
    res.gsva <- gsva(aggr.counts.i, gene.list.j, method = "gsva") 
    
    # AUCell
    cells_rankings <- AUCell_buildRankings(aggr.counts.i, plotStats = FALSE, verbose = F)
    
    res.aucell <- AUCell_calcAUC(
      gene.list.j, 
      cells_rankings, 
      aucMaxRank = nrow(cells_rankings)*0.05) %>% 
      getAUC()
    
    # JASMINE
    res.jas.lr <- JASMINE(aggr.counts.i, gene.list.j[[1]], method = 'likelihood')
    
    # transform
    res.ssgsea <- res.ssgsea %>% 
      as.data.frame() %>% 
      rownames_to_column("pathway") %>% 
      as_tibble() %>% 
      mutate(method = "ssgsea")
    
    res.gsva <- res.gsva %>% 
      as.data.frame() %>% 
      rownames_to_column("pathway") %>% 
      as_tibble() %>% 
      mutate(method = "gsva")
    
    res.aucell <- res.aucell %>% 
      as.data.frame() %>% 
      rownames_to_column("pathway") %>% 
      as_tibble() %>% 
      mutate(method = "AUCell")
    
    res.jas.lr <- res.jas.lr %>% 
      as_tibble() %>% 
      mutate(pathway = pathway.j) %>% 
      spread(key = "SampleID", value = "JAS_Scores") %>% 
      mutate(method = "JASMINE.LR")
    
    res.j <- bind_rows(res.ssgsea, res.gsva, res.aucell, res.jas.lr) 
    
    colnames(res.j)[2:7] <- paste0(rep(c("hfc_", "lfc_"), each = 3), 1:3)
    
    res.j <- res.j %>% 
      mutate(cell.type = name.i) %>% 
      dplyr::select(method, pathway, cell.type, paste0(rep(c("hfc_", "lfc_"), each = 3), 1:3))
    
    res.i <- res.i %>% 
      bind_rows(res.j)
  }
  res.merged <- res.merged %>% 
    bind_rows(res.i)
}

res.merged <- res.merged %>% 
  dplyr::select(method, cell.type, pathway, hfc_1:lfc_3) %>% 
  arrange(method, cell.type, pathway) 


# check -------------------------------------------------------------------

res.merged %>% print(n = Inf)


# save --------------------------------------------------------------------

setwd(dir.2)
out.f <- "004_gene_sig_scores_pseudo_bulk.tsv"
write_tsv(res.merged, out.f)


# transform the table: z-score -----------------------------------------------------------------

res.merged.2 <- res.merged %>% 
  gather("sample.id", "score", hfc_1:lfc_3) %>% 
  group_by(method, cell.type, pathway) %>% 
  mutate(scaled.score = as.vector(scale(score))) %>% # z-normalization
  ungroup() %>% 
  mutate(tissue = gsub("_.*", "", sample.id)) %>% 
  mutate(pathway = gsub("HALLMARK_", "", pathway)) %>% 
  mutate(patient = ifelse(grepl("1", sample.id), "pt_1", 
                          ifelse(grepl("2", sample.id), "pt_2", "pt_3"))) %>% 
  dplyr::select(method, cell.type, pathway, sample.id, tissue, patient, score, scaled.score) 


res.merged.2$method <- factor(res.merged.2$method, 
                              levels = c("gsva", "ssgsea", "AUCell", "JASMINE.LR"))
res.merged.2$cell.type <- factor(res.merged.2$cell.type, 
                                 levels = c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells"))
res.merged.2$pathway <- factor(res.merged.2$pathway, 
                               levels = c("INFLAMMATORY_RESPONSE", 
                                          "INTERFERON_GAMMA_RESPONSE",
                                          "TNFA_SIGNALING_VIA_NFKB"))
res.merged.2$tissue <- factor(res.merged.2$tissue, 
                              levels = c("hfc", "lfc"))

res.merged.2 <- res.merged.2 %>% 
  arrange(method, cell.type, pathway)


# check -------------------------------------------------------------------

res.merged.2 %>% print(n = Inf)


# save --------------------------------------------------------------------

setwd(dir.2)
out.f <- "005_gene_sig_scores_pseudo_bulk_z_norm.tsv"
write_tsv(res.merged.2, out.f)


# statistical tests 1. normality of data distribution-------------------------------------------------------

res.norm.merged <- NULL
for(i in 1:4){
  print(i)
  method.i <- c("gsva", "ssgsea", "AUCell", "JASMINE.LR")[i]
  
  for(j in 1:3){
    print(j)
    cell.type.j <- c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells")[j]
    
    for(k in 1:3){
      print(k)
      pathway.k <- c("INFLAMMATORY_RESPONSE",
                     "INTERFERON_GAMMA_RESPONSE",
                     "TNFA_SIGNALING_VIA_NFKB")[k]
      
      res.ijk <- res.merged.2 %>% 
        dplyr::filter(method == method.i & cell.type == cell.type.j & pathway == pathway.k) 
      
      for(l in 1:2){
        tissue.l <- c("hfc", "lfc")[l]
        data <- res.ijk %>% 
          dplyr::filter(tissue == tissue.l) %>% 
          pull(scaled.score)
        
        # Shapiro-Wilk Test (base R stats)
        res.shapiro <- shapiro.test(data)
        
        # # Anderson-Darling Test (nortest)
        # res.ad <- ad.test(data)
        
        # Kolmogorov-Smirnov Test (base R)
        res.ks <- ks.test(data, "pnorm", mean = mean(data), sd = sd(data))
        
        res.ijkl <- tibble(
          method = method.i, 
          cell.type = cell.type.j, 
          pathway = pathway.k, 
          tissue = tissue.l, 
          n = length(data), 
          shapiro.W = res.shapiro$statistic, 
          shapiro.p = res.shapiro$p.value, 
          shapiro.sum = if(res.shapiro$p.value >= 0.05){"ns"}else{"p<.05"}, 
          ks.D = res.ks$statistic, 
          ks.p = res.ks$p.value, 
          ks.sum = if(res.ks$p.value >= 0.05){"ns"}else{"p<.05"}
        ) %>% 
          mutate(summary = ifelse(shapiro.sum == "ns" & ks.sum == "ns", "NS", "check"))
        
        res.norm.merged <- res.norm.merged %>% 
          bind_rows(res.ijkl)
      }
    }
  }
}


# check -------------------------------------------------------------------

res.norm.merged %>% dim()
# [1] 72 12

res.norm.merged %>% print(n = Inf)


# out ---------------------------------------------------------------------

setwd(dir.2)
out.f <- "006_gene_sig_scores_pseudo_bulk_stat_01_normality_test_res.tsv"
write_tsv(res.norm.merged, out.f)


# note --------------------------------------------------------------------


# statistical test ---------------------------------------------------------------

res.stat.merged <- NULL
for(i in 1:4){
  print(i)
  method.i <- c("gsva", "ssgsea", "AUCell", "JASMINE.LR")[i]
  
  for(j in 1:3){
    print(j)
    cell.type.j <- c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells")[j]
    
    for(k in 1:3){
      print(k)
      pathway.k <- c("INFLAMMATORY_RESPONSE",
                     "INTERFERON_GAMMA_RESPONSE",
                     "TNFA_SIGNALING_VIA_NFKB")[k]
      
      res.ijk <- res.merged.2 %>% 
        dplyr::filter(method == method.i & cell.type == cell.type.j & pathway == pathway.k) 
      
      data.hfc <- res.ijk %>%
        dplyr::filter(tissue == "hfc") %>%
        pull(scaled.score)
      data.lfc <- res.ijk %>%
        dplyr::filter(tissue == "lfc") %>%
        pull(scaled.score)
      
      # paired t-test
      stat.ijk <- t.test(scaled.score ~ tissue, data = res.ijk, paired = T)
      
      res.stat.ijk <- tibble(
        method = method.i, 
        cell.type = cell.type.j, 
        pathway = pathway.k, 
        value.type = "mean", 
        mean.hfc = mean(data.hfc), 
        mean.lfc = mean(data.lfc), 
        sd.hfc = sd(data.hfc), 
        sd.lfc = sd(data.lfc), 
        test.used = stat.ijk$method, 
        p.value = stat.ijk$p.value
      ) %>% 
        mutate(sig = ifelse(p.value < 0.05, "sig", "ns"))
      res.stat.merged <- res.stat.merged %>% 
        bind_rows(res.stat.ijk)
    }
  }
}


# check -------------------------------------------------------------------

res.stat.merged %>% dim()
# [1] 36 11

res.stat.merged %>% print(n = Inf)


# out ---------------------------------------------------------------------

setwd(dir.2)
out.f <- "007_gene_sig_scores_pseudo_bulk_stat_02_stat_test_res.tsv"
write_tsv(res.stat.merged, out.f)


# Box plots (Extended Data Fig. S2g) -----------------------------------------------------------

res.merged.2$pathway <- factor(res.merged.2$pathway, 
                               levels = c("INFLAMMATORY_RESPONSE", 
                                          "INTERFERON_GAMMA_RESPONSE",
                                          "TNFA_SIGNALING_VIA_NFKB"))

COLOURS <- pal_nejm()(2)

gg.list <- list()
for(i in 1:3){
  print(i)
  cell.type.i <- c("Tumor Cells", "Myeloid Cells", "Lymphoid Cells")[i]
  
  res.merged.i <- res.merged.2 %>% 
    dplyr::filter(cell.type == cell.type.i)
  
  g <- ggplot(data = res.merged.i, aes(x = tissue, y = scaled.score, color = tissue, fill = tissue))
  g <- g + geom_line(aes(group = patient), color = "gray80")
  g <- g + geom_boxplot(width = 0.75, alpha = 0.25)
  g <- g + geom_point(aes(shape = patient), size = 1.5, position = position_jitter(width =.05))
  g <- g + scale_color_manual(values = COLOURS)
  g <- g + scale_fill_manual(values = COLOURS)
  g <- g + labs(title = cell.type.i, x = NULL, y = "signature score")
  g <- g + facet_grid(pathway ~ method)
  g <- g + theme_classic()
  g <- g + theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    strip.background = element_blank(), 
    strip.text = element_text(size = 8), 
    legend.position = "none", 
    legend.direction = "vertical",
    # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
    panel.grid.major = element_line(color = "#f0f0f0")
  )
  # plot(g)
  
  gg.list[[i]] <- g
}

g <- ggarrange(
  plotlist = gg.list, 
  ncol = 3, nrow = 1
)
# plot(g)

setwd(dir.2)
out.f <- "008_gene_sig_scores_pseudo_bulk_zscore_boxplot_combined.pdf"
ggsave(out.f, g, w = 12, h = 7.5)  


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
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] ggpubr_0.4.0                ggsci_2.9                   AUCell_1.14.0              
# [4] GSA_1.03.3                  GSVA_1.40.1                 data.table_1.14.2          
# [7] Matrix.utils_0.9.8          Matrix_1.6-0                SingleCellExperiment_1.14.1
# [10] DESeq2_1.32.0               SummarizedExperiment_1.22.0 Biobase_2.52.0             
# [13] MatrixGenerics_1.4.3        matrixStats_0.62.0          GenomicRanges_1.44.0       
# [16] GenomeInfoDb_1.28.1         IRanges_2.26.0              S4Vectors_0.30.2           
# [19] BiocGenerics_0.38.0         msigdbr_7.4.1               SeuratObject_4.1.3         
# [22] Seurat_4.3.0.1              tidylog_1.0.2               forcats_0.5.1              
# [25] stringr_1.4.0               dplyr_1.0.9                 purrr_0.3.4                
# [28] readr_2.1.2                 tidyr_1.2.0                 tibble_3.1.7               
# [31] ggplot2_3.3.6               tidyverse_1.3.1            
# 
# loaded via a namespace (and not attached):
# [1] utf8_1.2.2                R.utils_2.11.0            spatstat.explore_3.2-1   
# [4] reticulate_1.25           tidyselect_1.1.2          RSQLite_2.2.8            
# [7] AnnotationDbi_1.54.1      htmlwidgets_1.6.2         grid_4.1.3               
# [10] BiocParallel_1.26.2       Rtsne_0.16                ScaledMatrix_1.0.0       
# [13] munsell_0.5.0             codetools_0.2-18          ica_1.0-2                
# [16] future_1.25.0             miniUI_0.1.1.1            withr_2.5.0              
# [19] spatstat.random_3.1-5     colorspace_2.0-3          progressr_0.10.0         
# [22] rstudioapi_0.13           ROCR_1.0-11               ggsignif_0.6.2           
# [25] tensor_1.5                listenv_0.8.0             labeling_0.4.2           
# [28] GenomeInfoDbData_1.2.6    polyclip_1.10-0           farver_2.1.0             
# [31] bit64_4.0.5               rhdf5_2.36.0              parallelly_1.31.1        
# [34] vctrs_0.6.1               generics_0.1.2            timechange_0.2.0         
# [37] R6_2.5.1                  rsvd_1.0.5                locfit_1.5-9.4           
# [40] rhdf5filters_1.4.0        bitops_1.0-7              spatstat.utils_3.0-3     
# [43] cachem_1.0.6              DelayedArray_0.18.0       assertthat_0.2.1         
# [46] vroom_1.5.7               promises_1.2.0.1          scales_1.2.0             
# [49] gtable_0.3.0              beachmat_2.8.1            globals_0.15.0           
# [52] goftest_1.2-2             rlang_1.1.0               clisymbols_1.2.0         
# [55] genefilter_1.74.0         systemfonts_1.0.4         splines_4.1.3            
# [58] rstatix_0.7.0             lazyeval_0.2.2            spatstat.geom_3.2-4      
# [61] broom_0.7.9               reshape2_1.4.4            abind_1.4-5              
# [64] modelr_0.1.8              backports_1.2.1           httpuv_1.6.5             
# [67] tools_4.1.3               ellipsis_0.3.2            RColorBrewer_1.1-3       
# [70] ggridges_0.5.3            Rcpp_1.0.8.3              plyr_1.8.7               
# [73] sparseMatrixStats_1.4.2   zlibbioc_1.38.0           RCurl_1.98-1.3           
# [76] deldir_1.0-6              pbapply_1.5-0             cowplot_1.1.1            
# [79] zoo_1.8-10                grr_0.9.5                 haven_2.4.3              
# [82] ggrepel_0.9.1             cluster_2.1.2             fs_1.5.2                 
# [85] magrittr_2.0.3            scattermore_0.7           openxlsx_4.2.4           
# [88] lmtest_0.9-40             reprex_2.0.1              RANN_2.6.1               
# [91] fitdistrplus_1.1-5        hms_1.1.0                 patchwork_1.1.1          
# [94] mime_0.12                 xtable_1.8-4              XML_3.99-0.7             
# [97] rio_0.5.27                readxl_1.3.1              gridExtra_2.3            
# [100] compiler_4.1.3            KernSmooth_2.23-20        crayon_1.5.1             
# [103] R.oo_1.24.0               htmltools_0.5.5           later_1.3.0              
# [106] tzdb_0.1.2                geneplotter_1.70.0        lubridate_1.9.2          
# [109] DBI_1.1.2                 dbplyr_2.1.1              MASS_7.3-54              
# [112] babelgene_21.4            car_3.0-11                cli_3.6.1                
# [115] R.methodsS3_1.8.1         igraph_1.3.1              pkgconfig_2.0.3          
# [118] foreign_0.8-81            sp_2.0-0                  plotly_4.10.0            
# [121] spatstat.sparse_3.0-2     xml2_1.3.2                annotate_1.70.0          
# [124] XVector_0.32.0            rvest_1.0.1               digest_0.6.29            
# [127] sctransform_0.3.5         RcppAnnoy_0.0.19          graph_1.70.0             
# [130] spatstat.data_3.0-1       Biostrings_2.60.2         cellranger_1.1.0         
# [133] leiden_0.3.9              uwot_0.1.16               DelayedMatrixStats_1.14.3
# [136] GSEABase_1.54.0           curl_4.3.2                shiny_1.7.1              
# [139] lifecycle_1.0.3           nlme_3.1-152              jsonlite_1.8.0           
# [142] Rhdf5lib_1.14.2           carData_3.0-4             viridisLite_0.4.0        
# [145] fansi_1.0.3               pillar_1.7.0              lattice_0.20-44          
# [148] KEGGREST_1.32.0           fastmap_1.1.0             httr_1.4.3               
# [151] survival_3.2-12           glue_1.6.2                zip_2.2.0                
# [154] png_0.1-8                 bit_4.0.4                 HDF5Array_1.20.0         
# [157] stringi_1.7.6             blob_1.2.2                textshaping_0.3.6        
# [160] BiocSingular_1.8.1        memoise_2.0.0             irlba_2.3.5              
# [163] future.apply_1.8.1       


# end ---------------------------------------------------------------------


