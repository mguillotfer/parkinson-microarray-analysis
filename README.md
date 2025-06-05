# Analysis for the paper: Whole transcriptome analysis of peripheral blood mononuclear cells from de novo and drug-naïve sporadic Parkinson’s disease patients

This repository contains the code used in the analysis described in the manuscript:  
**Whole transcriptome analysis of peripheral blood mononuclear cells from de novo and drug-naïve sporadic Parkinson’s disease patients**  
Francisco Navarrete1,2,3*, Marina Guillot-Fernández1*, Lorena Martínez-Hostyn1,2,3, Daniela Navarro1,2,3, José A. Molina4, Jose P. López-Atalaya1, and Jorge Manzanares1,2,3

## Overview


The goal of this analysis is to:

- Identify differentially expressed genes (DEGs) between PD patients and controls
- Perform gene set enrichment analysis (GSEA)
- Infer cell-type proportions using CIBERSORTx

The pipeline is implemented in R using raw .CEL files from Affymetrix Clariom™ S Human arrays, publicly available on GEO.

## Repository structure

- scripts/ # R scripts for preprocessing, analysis, and plotting
- data/ # Input data (user must download separately)
- results/ # Output tables and intermediate files
- figures/ # Plots used in the paper
- README.md # This file
  
## Data

The raw microarray data used in this study (Affymetrix `.CEL` files) are publicly available from the Gene Expression Omnibus (GEO) under accession number **[GSE290333](https://urldefense.com/v3/__https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290333__;!!D9dNQwwGXtA!Qd4zkJ-MDO8yfX6hW5pvcIqYDnfychXFnNHriNUe4EQb_R57rfAWJ1VGDUAnICAVyobeKQy8mGcX4WwzXQ$)**.

To reproduce the analysis, download all CEL files from the GEO repository and place them in the `data/` folder.

In addition, gene-level annotation is performed using the Clariom™ S Human Array annotation file. This must be downloaded separately from Thermo Fisher Scientific: https://www.thermofisher.com/order/catalog/product/902927 (Current NetAffx Annotation Files: Clariom™ S Array, human Transcript Cluster Annotations, CSV, Release 36). 

## Requirements

- R version: 4.4.2
- Platform: Ubuntu 22.04.5 LTS

## Packages used

The analysis was performed using the following R packages (with versions):

affy_1.84.0
AnnotationDbi_1.68.0
Biobase_2.66.0
BiocGenerics_0.52.0
Biostrings_2.74.1
clusterProfiler_4.14.6
DBI_1.2.3
dplyr_1.1.4
enrichplot_1.26.6
forcats_1.0.0
ggplot2_3.5.2
gplots_3.2.0
IRanges_2.40.1
limma_3.62.2
lubridate_1.9.4
oligo_1.70.1
oligoClasses_1.68.0
org.Hs.eg.db_3.20.0
pd.clariom.s.human_3.14.1
purrr_1.0.4
RColorBrewer_1.1-3
readr_2.1.5
RSQLite_2.3.11
scales_1.4.0
S4Vectors_0.44.0
stringr_1.5.1
tibble_3.2.1
tidyr_1.3.1
tidyverse_2.0.0
XVector_0.46.0
GenomeInfoDb_1.42.3

