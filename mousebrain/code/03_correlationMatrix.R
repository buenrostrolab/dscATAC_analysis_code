library(data.table)
library(BuenColors)
library(dplyr)
library(cowplot)
library(matrixStats)
library(ComplexHeatmap)

# Gene scores
atac_gene_scores <- fread("../output/revision-mousebrain_rawGeneScores.tsv") %>% data.frame()
rownames(atac_gene_scores) <- as.character(fread("../output/mousebrain_geneAnnotations.tsv", header = FALSE)[[1]])
mb_anno <- read.table("../data/revision-cluster_annotations.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t", comment.char = "")
mb_trans <- as.character(mb_anno$new); names(mb_trans) <- as.character(mb_anno$old)
colnames(atac_gene_scores) <- mb_trans[as.character(colnames(atac_gene_scores))]
ordering <- as.character(mb_anno$new)

# Make a counts per million matrix
cpm2 <- t(t(atac_gene_scores)/colSums(atac_gene_scores)) * 1000000
colVec <- as.character(mb_anno$colVec); names(colVec) <- as.character(mb_anno$new)

# Import DropViz
anno_dropviz <- readRDS("../data/dropviz/annotation.BrainCellAtlas_Saunders_version_2018.04.01.RDS")
data_dropviz <- readRDS("../data/dropviz/metacells.BrainCellAtlas_Saunders_version_2018.04.01.RDS")

# Determine overlapping genes
gene_set <- unique(intersect(rownames(data_dropviz), rownames(atac_gene_scores)))
gene_set2 <- intersect(c(as.character(anno_dropviz$class_marker), as.character(anno_dropviz$type_marker)), gene_set)


cpm.Z <- t(t((cpm2-rowMeans(cpm2)))/(rowSds(as.matrix(cpm2)))[row(cpm2)])
#cpm.Z <- cpm.Z[rowMeans(cpm2) > 50,ordering]
cormatin <- cpm2[gene_set2,ordering]

cormat <- cor(cormatin, use = "pairwise.complete", method = "pearson")

pdf(paste0("../output/mouse-brain-cormat.pdf"), width=3, height=2.5)
hm <- Heatmap(cormat, 
              col=as.character(jdb_palette("solar_extra",type="continuous")),
              show_row_names = FALSE, 
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              show_column_names = FALSE)
hm
dev.off()
