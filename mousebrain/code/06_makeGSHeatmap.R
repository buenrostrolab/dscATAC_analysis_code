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

cpm.Z <- t(t((cpm2-rowMeans(cpm2)))/(rowSds(as.matrix(cpm2)))[row(cpm2)])
#cpm.Z <- cpm.Z[rowMeans(cpm2) > 50,ordering]
cpm.Z <- cpm.Z[,ordering]
boo <- rowMaxs(cpm.Z) > 3 & complete.cases(cpm.Z)
subset <- cpm.Z[boo,]

splitMe <- max.col(subset, 'first')

data.frame(cluster = ordering[splitMe], gene = rownames(cpm.Z)[boo]) %>%
  arrange(cluster) -> odf

write.table(odf, file = "../output/cluster_specific_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

if(FALSE){
  png(paste0("mouse-brain-diffgene.png"), width=5, height=7, units = "in", res = 300)
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
}
# Determine overlapping genes

anno_dropviz <- readRDS("../data/dropviz/annotation.BrainCellAtlas_Saunders_version_2018.04.01.RDS")
data_dropviz <- readRDS("../data/dropviz/metacells.BrainCellAtlas_Saunders_version_2018.04.01.RDS")

type_markers  <- lapply(as.character(anno_dropviz$type_marker), function(x) strsplit(x, "-|[.]")[[1]]) %>% unlist() %>% unique() %>% as.character()
class_markers  <- lapply(as.character(anno_dropviz$class_marker), function(x) strsplit(x, "-|[.]")[[1]]) %>% unlist() %>% unique() %>% as.character()

gene_set2 <- intersect(c(type_markers, class_markers, "Gad1"), rownames(subset))
splitMe2 <- max.col(subset[gene_set2, ], 'first')

if(FALSE){
  write.table(round(cpm2[gene_set2,ordering], 2) , file = "../output/TableS5_geneScores.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
}

subset[subset > 10] <- 10
subset[subset < -10] <- -10

png(paste0("../output/mouse-brain-diffgene.png"), width=3, height=2.8, units = "in", res = 600)
hm <- Heatmap(subset[gene_set2, ], 
              col=as.character(jdb_palette("solar_extra",type="continuous")),
              show_row_names = FALSE, 
              cluster_columns = FALSE,
              cluster_rows = TRUE,
              row_names_gp = gpar(fontsize = 0),
              column_names_gp = gpar(fontsize = 0),
              split = splitMe2, 
              show_column_names = FALSE)
hm
dev.off()

ha1 <- HeatmapAnnotation(df = data.frame(clusters = names(colVec)), 
                         col = list(clusters = colVec)
)

pdf(paste0("../output/mouse-brain-diffgene.pdf"), width=4.5, height=4)
hm <- Heatmap(subset[gene_set2, ], 
              col=as.character(jdb_palette("solar_extra",type="continuous")),
              show_row_names = FALSE, 
              cluster_columns = FALSE,
              cluster_rows = TRUE,
              row_names_gp = gpar(fontsize = 0),
              column_names_gp = gpar(fontsize = 5),
              split = splitMe2, 
              bottom_annotation = ha1,
              show_column_names = TRUE)
hm
dev.off()
