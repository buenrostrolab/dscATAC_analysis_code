library(dplyr)

# Create new synthetic bigwigs
cdf <- readRDS("../data/mousebrain-master_dataframe.rds")
odf <- data.frame(
  tSNE1 = round(cdf$X1,2),
  tSNE2 = round(cdf$X2, 2),
  cluster_id = paste0("V", as.character(cdf$clusters)),
  cell_id = cdf$DropBarcode
)
write.table(odf, file = "../../coords_clusters_annotations/mousebrain_tSNE_coordinates_Fig2a.tsv", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)