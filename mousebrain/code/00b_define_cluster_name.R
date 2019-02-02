library(data.table)
library(dplyr)
library(stringr)
library(BuenColors)

# Import gene scores
atac_gene_scores <- fread("../../br-mouse-brain/output/revision-mousebrain_rawGeneScores.tsv") %>% data.frame()
rownames(atac_gene_scores) <- as.character(fread("../../br-mouse-brain/output/mousebrain_geneAnnotations.tsv", header = FALSE)[[1]])

# Import DropViz
anno_dropviz <- readRDS("../dropviz/annotation.BrainCellAtlas_Saunders_version_2018.04.01.RDS")
data_dropviz <- readRDS("../dropviz/metacells.BrainCellAtlas_Saunders_version_2018.04.01.RDS")

# Determine overlapping genes
gene_set <- unique(intersect(rownames(data_dropviz), rownames(atac_gene_scores)))

type_markers  <- lapply(as.character(anno_dropviz$type_marker), function(x) strsplit(x, "-|[.]")[[1]]) %>% unlist() %>% unique() %>% as.character()
class_markers  <- lapply(as.character(anno_dropviz$class_marker), function(x) strsplit(x, "-|[.]")[[1]]) %>% unlist() %>% unique() %>% as.character()

gene_set2 <- intersect(c(type_markers, class_markers, "Gad1"), gene_set)

atac_gene_scores_s <- atac_gene_scores[gene_set2, ]
dropviz_s <- data_dropviz[gene_set2,]
cormat <- cor(atac_gene_scores_s, dropviz_s, method = "spearman")

atac_cor <- cor(atac_gene_scores_s, method = "spearman")

reshape2::melt(cormat) %>% arrange(Var1,desc(value)) %>%  group_by(Var1) %>% top_n(1) %>%
  left_join(., anno_dropviz[,c("tissue_subcluster", "class", "class_marker", "type_marker", "common_name")], by = c("Var2" = "tissue_subcluster")) %>% 
  arrange(class_marker) %>% data.frame() -> top_assign

# From this have the broad classes above, provide unique names
EN <- paste0("V", as.character(c(2,3,4,7,8,10,13,15,17,18,19,20,21,22,23,25,16)))
IN <- paste0("V", as.character(c(5,9,11,12,14)))

atac_cor["V25", EN] %>% sort(decreasing = TRUE) %>% 
  names() -> old_EN

atac_cor["V12", IN] %>% sort(decreasing = TRUE) %>% 
  names() -> old_IN

new_EN <- paste0("EN", str_pad(1:17, 2, pad = "0"))
new_IN <- paste0("IN", str_pad(1:5, 2, pad = "0"))

colVec <- c(
  colorRampPalette(jdb_palette("purple_baby"))(17), # excitatory
  rev(jdb_palette("brewer_red")[4:8]), # inhibitory
  jdb_palette("BottleRocket2")[c(1)], # microglia
  "#CD8500", # endothelial
  "#006400", # astrocytes
  "#1C86EE", # oligodendrocytes
  "#094078" # pd
)

translate_df <- data.frame(
  old = c(old_EN, old_IN, "V26", "V27", "V6", "V1", "V24"),
  new = c(new_EN, new_IN, "MG1", "E1", "A1", "OG1", "P1"),
  colVec, 
  type = c(rep("Excite-NEURON", 17), rep("Inhibit-NEURON", 5), "MICROGLIA", "ENDOTHELIAL", "ASTROCYTE", "OLIGODENDROCYTE", "POLYDENDROCYTE"),
  stringsAsFactors = FALSE)

vec <- as.character(translate_df$colVec); names(vec) <- as.character(translate_df$new)

ggplot(translate_df, aes(x = new, y = 1, fill = new)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = vec)

write.table(translate_df, file = "../../br-mouse-brain/data/revision-cluster_annotations.tsv",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)




