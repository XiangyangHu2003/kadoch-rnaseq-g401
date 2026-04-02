library(DESeq2)
library(tidyverse)

counts <- read.csv("data/processed/count_matrix.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("data/metadata/metadata.csv", stringsAsFactors = FALSE)

rownames(meta) <- meta$sample_id
meta <- meta[, c("condition", "timepoint", "replicate")]
meta <- meta[colnames(counts), , drop = FALSE]
stopifnot(all(colnames(counts) == rownames(meta)))

meta$condition <- factor(meta$condition, levels = c("Empty", "BAF47"))

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData = meta,
  design = ~ condition
)

dds <- dds[rowSums(counts(dds)) > 10, ]

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "BAF47", "Empty"))

res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  arrange(padj)

write.csv(res_df, "results/tables/deseq2_results_all.csv", row.names = FALSE)

res_sig <- res_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)

write.csv(res_sig, "results/tables/deseq2_results_sig.csv", row.names = FALSE)

png("results/figures/ma_plot.png", width = 800, height = 600)
plotMA(res, main = "DESeq2 MA plot", ylim = c(-10, 10))
dev.off()

cat("Number of significant genes:\\n")
print(nrow(res_sig))