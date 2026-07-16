# Bulk RNA-seq analysis of SMARCB1/BAF47 restoration in G401 cells

An end-to-end reanalysis of public bulk RNA-seq data examining how
**SMARCB1/BAF47 restoration** is associated with gene-expression changes in
G401 malignant rhabdoid tumor cells.

This portfolio project demonstrates practical experience with count-matrix
construction, sample-level quality control, differential expression, pathway
enrichment, scientific visualization, and interpretation of a small public
transcriptomics dataset.

## Project at a glance

| Item | Description |
|---|---|
| Biological question | Which transcriptional programs differ between BAF47-restored and empty-vector G401 cells? |
| Public dataset | [GEO GSE90633](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90633) |
| Analysis subset | G401, Day 3, two biological replicates per condition |
| Comparison | BAF47 versus Empty |
| Main methods | DESeq2 differential expression and Hallmark preranked GSEA |
| Key limitation | Exploratory analysis with `n = 2` per condition |

## Selected results

The current workflow tested 17,220 genes after filtering. Using
`padj < 0.05` and `|log2 fold change| > 1`, it identified 2,314 genes in the
selected contrast (1,654 with positive and 660 with negative log2 fold
change). The samples separate by condition in the current PCA, while Hallmark
GSEA highlights positive enrichment of EMT- and differentiation-associated
programs and negative enrichment of MYC-, E2F-, and cell-cycle-associated
programs.

These findings should be treated as **exploratory associations**, not as
independent biological validation, because the analysis includes only one cell
line, one time point, and two replicates per condition.

### PCA

![PCA plot](results/figures/pca_plot.png)

### Differential expression

![Volcano plot](results/figures/volcano_plot.png)

Additional figures and a detailed interpretation are available in the
[technical report](report/g401_smarcb1_rnaseq_technical_report.md).

## Dataset and provenance

The analysis uses the RNA-seq subset of GEO series
[GSE90633](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90633), part
of SuperSeries [GSE90634](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90634).
The original experiment compared empty-vector and constitutive BAF47
expression in malignant rhabdoid tumor cell lines. This repository focuses on
the following four Day 3 G401 samples:

| Local sample ID | Condition | Replicate | GEO sample |
|---|---|---|---|
| `G401_Empty_Day3_A` | Empty | A | [GSM2495977](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2495977) |
| `G401_Empty_Day3_B` | Empty | B | [GSM2495978](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2495978) |
| `G401_BAF47_Day3_A` | BAF47 | A | [GSM2495979](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2495979) |
| `G401_BAF47_Day3_B` | BAF47 | B | [GSM2495980](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2495980) |

The sample-level `read_cnt.txt` files in `data/raw/` are processed count files
distributed with the GEO record; this project begins from those counts rather
than FASTQ files. Consequently, read alignment and transcript quantification
are outside the scope of this repository.

## Workflow

```text
GEO sample-level counts
        |
        v
Merged gene-by-sample count matrix
        |
        +--> Library size, PCA, sample correlation
        |
        +--> DESeq2: BAF47 vs Empty
                  |
                  +--> MA plot, volcano plot, top-gene heatmap
                  |
                  +--> Hallmark preranked GSEA
```

Scripts are organized by execution order:

1. `scripts/01_build_count_matrix.R` — merge sample-level count files
2. `scripts/02_qc_pca.R` — library-size QC, VST, PCA, and correlation
3. `scripts/03_deseq2.R` — differential expression analysis
4. `scripts/04_visualization.R` — volcano, MA, and heatmap outputs
5. `scripts/05_gsea.R` — Hallmark gene set enrichment analysis

## Repository structure

```text
kadoch-rnaseq-g401/
|-- data/
|   |-- metadata/       # analysis metadata and GEO series metadata
|   |-- processed/      # merged count matrix
|   `-- raw/            # GEO sample-level count files
|-- scripts/            # ordered R analysis scripts
|-- results/
|   |-- figures/        # analysis figures
|   `-- tables/         # machine-readable result tables
|-- report/             # technical report
`-- README.md
```

## Running the analysis

Run the scripts from the repository root in numeric order. The workflow uses
the following R packages:

- `tidyverse`
- `DESeq2`
- `pheatmap`
- `ggrepel`
- `msigdbr`
- `fgsea`

Expected inputs are the four count files in `data/raw/` and
`data/metadata/metadata.csv`. Generated tables and figures are written to
`results/tables/` and `results/figures/`.

> Reproducibility note: package versions are not yet locked in this version of
> the project. The committed result files document the current analysis output,
> but exact cross-version reproduction is not guaranteed until an environment
> lockfile and session information are added.

## Interpretation boundaries

- The dataset contains only two biological replicates per condition.
- The analysis covers one cell line and one time point.
- The model estimates the condition contrast only; no additional covariates are
  included.
- Pathway enrichment is hypothesis-generating and does not establish pathway
  activity experimentally.
- No independent dataset or orthogonal assay is used for validation.
- This repository starts from processed gene counts and does not benchmark the
  original alignment or quantification workflow.

## Reference

Nakayama RT, Pulice JL, Valencia AM, et al. SMARCB1 is required for widespread
BAF complex-mediated activation of enhancers and bivalent promoters.
*Nature Genetics*. 2017;49(11):1613-1623.
[doi:10.1038/ng.3958](https://doi.org/10.1038/ng.3958) ·
[PMID: 28945250](https://pubmed.ncbi.nlm.nih.gov/28945250/)

## Author

Xiangyang Hu

This repository was created as a computational biology portfolio project. It
is an independent reanalysis of public data and is not the original study's
analysis repository.
