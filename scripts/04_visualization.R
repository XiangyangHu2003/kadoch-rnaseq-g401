library(tidyverse)
library(ggrepel)
library(DESeq2)
library(pheatmap)

res_df <- read.csv("results/tables/deseq2_results_all.csv")

res_df <- res_df |>
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    ),
    neglog10_padj = -log10(padj)
  )

top_labels <- res_df |>
  filter(significance != "NS") |>
  slice_min(padj, n = 15)

ggplot(res_df, aes(log2FoldChange, neglog10_padj, color = significance)) +
  geom_point(alpha = 0.7, size = 0.6) +
  geom_text_repel(data = top_labels, aes(label = gene), size = 3, show.legend = FALSE) +
  theme_bw()

ggsave("results/figures/volcano_plot.png", width = 7, height = 5)

# Heatmap of top 30 DE genes
counts <- read.csv("data/processed/count_matrix.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("data/metadata/metadata.csv")
rownames(meta) <- meta$sample_id
meta <- meta[, c("condition", "timepoint", "replicate")]

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData = meta,
  design = ~ condition
)

dds <- dds[rowSums(counts(dds)) > 10, ]
dds$condition <- relevel(factor(dds$condition), ref = "Empty")
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "BAF47", "Empty"))
res_df <- as.data.frame(res)
top_genes <- rownames(res_df[order(res_df$padj), ])[1:30]

vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)[top_genes, ]
mat <- t(scale(t(mat)))

png("results/figures/top30_heatmap.png", width = 700, height = 500)
pheatmap(mat, annotation_col = meta)
dev.off()
