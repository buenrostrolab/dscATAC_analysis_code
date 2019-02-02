library(SummarizedExperiment)
library(irlba)
library(BuenColors)
library(chromVAR)
library(igraph)
library(Rtsne)
library(ggrastr)
library(stringr)
library(cowplot)
library(ComplexHeatmap)

# Import and normalize counts data
counts <- data.matrix(data.frame(fread("../data/revision_mouseBrain_Countsmatrix.tsv", header = TRUE)))
cpm <- t(t(counts)/colSums(counts)) * 1000000

# Reannotate with the new cluster definitions
mb_anno <- read.table("../data/revision-cluster_annotations.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t", comment.char = "")
mb_trans <- as.character(mb_anno$new); names(mb_trans) <- gsub("V", "X", as.character(mb_anno$old))
colVec <- as.character(mb_anno$colVec); names(colVec) <- as.character(mb_anno$new)
ordering <- as.character(mb_anno$new)

colnames(cpm) <- mb_trans[colnames(cpm)]
boo0 <- rowMeans(cpm) > 1
cpm2 <- cpm[boo0,]
cpm.Z <- t(t((cpm2-rowMeans(cpm2)))/(rowSds(as.matrix(cpm2)))[row(cpm2)])
cpm.Z <- cpm.Z[,ordering]
boo <- rowMaxs(cpm.Z) > 3 & complete.cases(cpm.Z)
subset <- cpm.Z[boo,]


subset[subset > 3] <- 3
subset[subset < -3] <- -3

png(paste0("../output/mouse-brain-diffEnhancer.png"), width=5, height=7, units = "in", res = 300)
hm <- Heatmap(subset, 
              col=as.character(jdb_palette("solar_extra",type="continuous")),
              show_row_names = FALSE, 
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              row_names_gp = gpar(fontsize = 0),
              column_names_gp = gpar(fontsize = 0),
              split = splitMe, 
              show_column_names = TRUE)
hm
dev.off()


