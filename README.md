# Glioma-neuron-immune crosstalk paper

## Introduction  

This repository contains the code used to analyze the data in Nejo et al. bioRxiv (2023) "Challenges in the discovery of tumor-specific alternative splicing-derived cell-surface antigens in glioma" (doi: [10.1101/2023.10.26.564156](https://www.biorxiv.org/content/10.1101/2023.10.26.564156v2.full)). In this study, we investigated tumor-specific alternative splicing-derived, cell-surface antigens in glioma, through the analyses of transcriptome data of TCGA-GBM/LGG, GTEx, Mayo-GBM-PDX, as well as the clinical samples of the University of California, San Francisco (UCSF) Brain Tumor Center. Moreover, we conducted nanopore full-length transcript sequencing and proteomics analysis of the CPTAC-GBM project. Our investigation illustrated the diverse characteristics of the tumor-specific AS events and the challenges of antigen exploration due to their notable spatiotemporal heterogeneity and elusive nature at the protein levels. 

  
## Primary contact: 
  
**Takahide Nejo, MD, PhD**  
Postdoctoral Scholar  
University of California, San Francisco, Department of Neurological Surgery  
takahide.nejo@ucsf.edu  
  
  
**Hideho Okada, MD, PhD**  
Professor  
University of California, San Francisco, Department of Neurological Surgery  
hideho.okada@ucsf.edu  
  
  
## Questions about the code:  
  
**Takahide Nejo, MD, PhD**  
takahide.nejo@ucsf.edu  

  
## Key analysis tools:  

[Seurat](https://github.com/alexdobin/STAR): Short-read RNA-seq read alignment and junction detection  
[SPATA2](https://themilolab.github.io/SPATA2/): An analytic platform of spatial transcriptomics data.  
[SPATAData](https://themilolab.github.io/SPATA2/articles/spata-data.html): A database of spatial transcriptomic samples that are made publicly available.  
[fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html): fgsea.  
[GSVA](https://www.bioconductor.org/packages/release/bioc/html/GSVA.html): GSVA.  
[AUCell](https://www.bioconductor.org/packages/release/bioc/html/AUCell.html): AUCell.  
[JASMINE](https://github.com/NNoureen/JASMINE): JASMINE.  


## Data availability:  
Human glioblastoma scRNAseq data (Krishna S, 2023 Nature) is available through the NCBI Gene Expression Omnibus (GEO) (accession, [GSE223065](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223065)). 

The mouse spatial transcriptomics and bulk RNA-seq dataset generated in this study will be made available through the NCBI GEO website upon manuscript acceptance.

## Key data visualization tools:  
[BioRender](https://www.biorender.com/)  
[Affinity Designer](https://affinity.serif.com/en-us/designer/)  
[R package ggplot2](https://ggplot2.tidyverse.org/)  

## Citation

**If you use our data or code for your work, please consider citing the following publication:**  

Nejo T, Krishna S, Jimenez C, Yamamichi A, Young JS, Lakshmanachetty S, Chen T, Phyu SSS, Ogino H, Watchmaker P, Diebold D, Choudhury A, Daniel AGS, Raleigh DR, Hervey-Jumper SL, Okada H. Glioma-neuronal circuit remodeling induces regional immunosuppression. bioRxiv [Preprint](https://www.biorxiv.org/content/10.1101/2023.08.04.548295v1). 2023 Aug 6:2023.08.04.548295. doi: 10.1101/2023.08.04.548295. PMID: 37577659 (https://pubmed.ncbi.nlm.nih.gov/37577659/) ; PMCID: PMC10418167.

  
  
