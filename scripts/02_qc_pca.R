library(tidyverse)
library(DESeq2)
library(pheatmap)

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

libsize <- colSums(counts(dds))
lib_df <- data.frame(
  sample = names(libsize),
  library_size = libsize,
  condition = meta[names(libsize), "condition"]
)

ggplot(lib_df, aes(x = sample, y = library_size, fill = condition)) +
  geom_col() +
  coord_flip() +
  theme_bw()

ggsave("results/figures/library_size_barplot.png", width = 7, height = 4)
write.csv(lib_df, "results/tables/library_size.csv", row.names = FALSE)

vsd <- vst(dds, blind = TRUE)

pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  geom_text(aes(label = name), vjust = -1, show.legend = FALSE) +
  theme_bw() +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15))

ggsave("results/figures/pca_plot.png", width = 6, height = 5)
write.csv(pcaData, "results/tables/pca_coordinates.csv", row.names = FALSE)

sample_cor <- cor(assay(vsd), method = "pearson")

png("results/figures/sample_correlation_heatmap.png", width = 900, height = 800)
pheatmap(sample_cor, main = "Sample-to-sample correlation")
dev.off()

write.csv(sample_cor, "results/tables/sample_correlation_matrix.csv")

