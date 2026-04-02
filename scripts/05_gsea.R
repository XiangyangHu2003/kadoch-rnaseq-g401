library(tidyverse)
library(msigdbr)
library(fgsea)

res_df <- read.csv("results/tables/deseq2_results_all.csv")

ranked_stats <- res_df |>
  filter(!is.na(stat)) |>
  select(gene, stat)

ranks <- ranked_stats$stat
names(ranks) <- ranked_stats$gene
ranks <- sort(ranks, decreasing = TRUE)

msig_df <- msigdbr(species = "Homo sapiens", collection = "H")
msig_h  <- split(x = msig_df$gene_symbol, f = msig_df$gs_name)

fgsea_res <- fgsea(pathways = msig_h, stats = ranks)

fgsea_out <- fgsea_res |>
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))

write.csv(fgsea_out, "results/tables/Hallmark_fgsea_results.csv", row.names = FALSE)
