# Reproducible bulk RNA-seq analysis of SMARCB1/BAF47 re-expression in G401 cells

## Overview

This project presents a reproducible bulk RNA-seq workflow for analyzing transcriptional changes associated with SMARCB1/BAF47 re-expression in G401 malignant rhabdoid tumor cells.

Using a public RNA-seq dataset, I focused on four G401 Day 3 samples:

- G401_Empty_Day3_A
- G401_Empty_Day3_B
- G401_BAF47_Day3_A
- G401_BAF47_Day3_B

The workflow includes:

- count matrix construction from per-sample read count files
- sample-level QC and PCA
- differential expression analysis with DESeq2
- publication-style visualization
- Hallmark gene set enrichment analysis with fgsea

This repository is designed to demonstrate practical skills in transcriptomics, reproducible research, statistical analysis, and figure generation in R.

## Biological question

What transcriptional programs change after SMARCB1/BAF47 re-expression in G401 malignant rhabdoid tumor cells?

This project compares BAF47 re-expression versus Empty control in G401 cells and identifies differentially expressed genes and enriched biological pathways associated with this perturbation.

## Repository structure

```text
kadoch-rnaseq-g401/
├── data/
│   ├── raw/
│   ├── metadata/
│   └── processed/
├── scripts/
├── results/
│   ├── figures/
│   └── tables/
├── report/
└── README.md
```

### Folder descriptions

- `data/raw/`  
  Raw per-sample read count files (`*.read_cnt.txt`)

- `data/metadata/`  
  Sample metadata file (`metadata.csv`)

- `data/processed/`  
  Processed count matrix built from raw files

- `scripts/`  
  R scripts for each step of the workflow

- `results/figures/`  
  Generated plots and heatmaps

- `results/tables/`  
  Differential expression tables, PCA coordinates, library sizes, and GSEA results

- `report/`  
  Space for a future Quarto or R Markdown report

## Input data

The analysis uses four sample-level read count files:

- `G401_Empty_Day3_A.read_cnt.txt`
- `G401_Empty_Day3_B.read_cnt.txt`
- `G401_BAF47_Day3_A.read_cnt.txt`
- `G401_BAF47_Day3_B.read_cnt.txt`

Each raw count file is expected to contain five columns:

1. `GeneSymbol`
2. `GeneName`
3. `Read.Count`
4. `Length`
5. `RPKM`

The workflow uses the `GeneSymbol` and `Read.Count` columns to construct a merged count matrix for downstream analysis.

## Metadata

A sample metadata file is stored at:

```text
data/metadata/metadata.csv
```

Example format:

```csv
sample_id,condition,timepoint,replicate
G401_Empty_Day3_A,Empty,Day3,A
G401_Empty_Day3_B,Empty,Day3,B
G401_BAF47_Day3_A,BAF47,Day3,A
G401_BAF47_Day3_B,BAF47,Day3,B
```

## Software and R packages

### Software
- R
- RStudio

### Main R packages
- tidyverse
- DESeq2
- ggplot2
- ggrepel
- pheatmap
- msigdbr
- fgsea

You can install the required packages with:

```r
install.packages(c("tidyverse", "ggrepel", "pheatmap", "msigdbr"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "fgsea"))
```

## Workflow

### Step 1. Build merged count matrix

Script:

```text
scripts/01_build_count_matrix.R
```

This script:

- reads the four raw read count files
- assigns column names
- extracts `GeneSymbol` and `Read.Count`
- merges all samples into one count matrix
- writes output to:

```text
data/processed/count_matrix.csv
```

### Step 2. Quality control and PCA

Script:

```text
scripts/02_qc_pca.R
```

This script:

- loads the processed count matrix and metadata
- constructs a `DESeqDataSet`
- filters low-count genes
- calculates library sizes
- applies variance stabilizing transformation
- generates sample-level PCA
- computes sample-to-sample correlation

Outputs:

```text
results/figures/library_size_barplot.png
results/figures/pca_plot.png
results/figures/sample_correlation_heatmap.png

results/tables/library_size.csv
results/tables/pca_coordinates.csv
results/tables/sample_correlation_matrix.csv
```

### Step 3. Differential expression analysis

Script:

```text
scripts/03_deseq2.R
```

This script:

- runs DESeq2 using the design `~ condition`
- compares `BAF47` versus `Empty`
- exports all genes and significant genes separately
- generates an MA plot

Outputs:

```text
results/tables/deseq2_results_all.csv
results/tables/deseq2_results_sig.csv
results/figures/ma_plot.png
```

Significant genes are defined as:

- `padj < 0.05`
- `|log2FoldChange| > 1`

### Step 4. Visualization of differential expression

Script:

```text
scripts/04_visualization.R
```

This script:

- reads DESeq2 results
- classifies genes as Up, Down, or NS
- generates a volcano plot
- labels the top significant genes
- rebuilds a DESeq2 object
- selects the top 30 differentially expressed genes
- creates a heatmap based on variance-stabilized expression values

Outputs:

```text
results/figures/volcano_plot.png
results/figures/top30_heatmap.png
```

### Step 5. Hallmark GSEA

Script:

```text
scripts/05_gsea.R
```

This script:

- ranks genes using the DESeq2 Wald statistic
- retrieves Hallmark gene sets using `msigdbr`
- runs preranked GSEA with `fgsea`
- exports pathway enrichment results

Output:

```text
results/tables/Hallmark_fgsea_results.csv
```

## Recommended script execution order

Run the scripts in the following order:

```text
01_build_count_matrix.R
02_qc_pca.R
03_deseq2.R
04_visualization.R
05_gsea.R
```

## Main outputs

### Figures
- library size bar plot
- PCA plot
- sample correlation heatmap
- MA plot
- volcano plot
- top 30 DE gene heatmap

### Tables
- merged count matrix
- library size table
- PCA coordinates
- sample correlation matrix
- all DESeq2 results
- significant DESeq2 results
- Hallmark GSEA results

## Notes on reproducibility

This repository is organized as a stepwise and reproducible RNA-seq analysis workflow.

To improve full reproducibility further, future versions may include:

- an R Markdown or Quarto report
- package version locking with `renv`
- automated download instructions for public input files
- expanded analysis including additional cell lines or time points

## Current limitation

This first version focuses only on:

- one cell line: `G401`
- one time point: `Day 3`
- one comparison: `BAF47 vs Empty`

This reduced design is intentional so that the project remains clear, reproducible, and easy to interpret.

## Potential future extensions

Possible next steps include:

- adding G401 Day 7 samples
- adding TTC1240 samples
- comparing Day 3 and Day 7 responses
- integrating additional genomics layers
- packaging the workflow into a more formal report

## Summary

This project demonstrates a complete and reproducible RNA-seq analysis workflow in R, including:

- count matrix construction
- sample QC
- PCA
- differential expression analysis with DESeq2
- visualization of key results
- Hallmark pathway enrichment analysis

It is intended as a compact transcriptomics portfolio project for computational biology and bioinformatics applications.
