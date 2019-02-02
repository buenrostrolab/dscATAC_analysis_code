library(dplyr)

z <- readRDS("../output/Zeisel_contribution.rds")
s <- readRDS("../output/saunders_contribution.rds")

total_merge <- left_join(s, z, by = c("atac_samples" = "atac_cluster"))

write.table(total_merge, file = "../output/forSupplementalTable.tsv", sep= "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
