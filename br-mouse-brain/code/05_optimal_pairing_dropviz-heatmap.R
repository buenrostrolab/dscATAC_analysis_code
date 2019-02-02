library(matchingR)
library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(BuenColors)

# Import gene scores
atac_gene_scores <- fread("../output/revision-mousebrain_rawGeneScores.tsv") %>% data.frame()
rownames(atac_gene_scores) <- as.character(fread("../output/mousebrain_geneAnnotations.tsv", header = FALSE)[[1]])
mb_anno <- read.table("../data/revision-cluster_annotations.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t", comment.char = "")
mb_trans <- as.character(mb_anno$new); names(mb_trans) <- as.character(mb_anno$old)
colnames(atac_gene_scores) <- mb_trans[as.character(colnames(atac_gene_scores))]
ordering <- as.character(mb_anno$new)
class <- as.character(mb_anno$type)
atac_gene_scores <- atac_gene_scores[,ordering]

# Import DropViz
anno_dropviz <- readRDS("../data/dropviz/annotation.BrainCellAtlas_Saunders_version_2018.04.01.RDS")
data_dropviz <- readRDS("../data/dropviz/metacells.BrainCellAtlas_Saunders_version_2018.04.01.RDS")

# Determine overlapping genes
type_markers  <- lapply(as.character(anno_dropviz$type_marker), function(x) strsplit(x, "-|[.]")[[1]]) %>% unlist() %>% unique() %>% as.character()
class_markers  <- lapply(as.character(anno_dropviz$class_marker), function(x) strsplit(x, "-|[.]")[[1]]) %>% unlist() %>% unique() %>% as.character()
gene_set <- unique(intersect(rownames(data_dropviz), rownames(atac_gene_scores)))
gene_set2 <- intersect(c(as.character(type_markers), as.character(class_markers)), gene_set)

atac_gene_scores_s <- atac_gene_scores[gene_set2, ]
dropviz_s <- data_dropviz[gene_set2,]

cormat <- t(cor(atac_gene_scores_s, dropviz_s, method = "spearman"))

stopifnot(all(colnames(data_dropviz) == anno_dropviz$tissue_subcluster))
keep_c <- c("Excite-NEURON", "Inhibit-NEURON", "MICROGLIA", "ENDOTHELIAL", "ASTROCYTE", "OLIGODENDROCYTE", "POLYDENDROCYTE")
str_b <- as.character(anno_dropviz$class)
str_b <- ifelse(str_b == "NEURON", ifelse(grepl("G", anno_dropviz$class_marker), "Inhibit-NEURON", "Excite-NEURON"), str_b)
str_b <-  gsub("_TIP", "", gsub("_STALK", "", str_b))
boo <- str_b %in% keep_c
str_go <- str_b[boo]

cormat.mm <- t(apply(cormat, 1, function(x)(x-min(x))/(max(x)-min(x))))

splitMe2 <-  as.numeric(factor(str_go, levels = keep_c))

pdf(paste0("../output/mouse-brain-global-assign.pdf"), width=4.2, height=3)
hm <- Heatmap(cormat.mm[boo,], 
              col=as.character(jdb_palette("brewer_spectra",type="continuous")),
              show_row_names = FALSE, 
              cluster_columns = FALSE,
              cluster_rows = TRUE,
              split = splitMe2, 
              show_column_names = FALSE)
hm
dev.off()

