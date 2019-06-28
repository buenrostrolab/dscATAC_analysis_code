library(dplyr)

# import data
df <- readRDS("../data/all-umap-exp100.rds")
tsne <- read.table("../data/tsneControls_with_info.txt", header = TRUE)
cluster_id <- rbind(read.table("../data/Baseline-Cells.tsv"), read.table("../data/Stimulation-Cells.tsv"))
old <- c("CD34-like","CMP-like","Early-Ery","Late-Ery","pDC-like","CLP-like-1","CLP-like-2","B-like-2","B-like-1","Mono-like","mDC-like","CD4-like","CD8-like","NK-like","Collision")
new <- c("HSPC","HSPC-ery","Ery-early","Ery-late","pDC","CLP","proB","preB","B","Mono-1","Mono-2","CD4","CD8","NK","Collision")
convert <- new; names(convert) <- old

cluster_id$cluster <- convert[as.character(cluster_id$V2)]
umap_df <- merge(df, cluster_id, by.x = "CellID", by.y = "V1")
odf <- umap_df[,c("CellID", "Condition", "Group", "cluster", "umap1", "umap2")]; colnames(odf) <- c("cell_id", "condition", "group", "cluster", "UMAP1", "UMAP2")

write.table(odf, file = "../../coords_clusters_annotations/bonemarrow_umap_annotations_FigS9a.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

mdf <- left_join(tsne, cluster_id, by = c("DropBarcode" = "V1"))
odf2 <- data.frame(
  tSNE1 = round(mdf$tSNE1, 2), 
  tSNE2 = round(mdf$tSNE2,2),
  cluster = mdf$cluster,
  cell_id = mdf$DropBarcode
)

write.table(odf2, file = "../../coords_clusters_annotations/bonemarrow_unstimulated_tSNE_annotations_Fig4c.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
