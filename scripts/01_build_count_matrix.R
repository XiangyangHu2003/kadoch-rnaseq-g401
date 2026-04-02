library(tidyverse)

files <- c(
  "data/raw/G401_Empty_Day3_A.read_cnt.txt",
  "data/raw/G401_Empty_Day3_B.read_cnt.txt",
  "data/raw/G401_BAF47_Day3_A.read_cnt.txt",
  "data/raw/G401_BAF47_Day3_B.read_cnt.txt"
)

sample_names <- c(
  "G401_Empty_Day3_A",
  "G401_Empty_Day3_B",
  "G401_BAF47_Day3_A",
  "G401_BAF47_Day3_B"
)

count_list <- map(files, ~ read.delim(.x, header = FALSE))

count_list <- map(count_list, ~ {
  colnames(.x) <- c("GeneSymbol", "GeneName", "Read.Count", "Length", "RPKM")
  .x
})

names(count_list) <- sample_names

counts_merged <- count_list[[1]] |>
  dplyr::select(GeneSymbol, Read.Count) |>
  dplyr::rename(!!sample_names[1] := Read.Count)

for (i in 2:length(count_list)) {
  tmp <- count_list[[i]] |>
    select(GeneSymbol, Read.Count) |>
    dplyr::rename(!!sample_names[i] := Read.Count)
  
  counts_merged <- full_join(counts_merged, tmp, by = "GeneSymbol")
}

write.csv(counts_merged, "data/processed/count_matrix.csv", row.names = FALSE)
