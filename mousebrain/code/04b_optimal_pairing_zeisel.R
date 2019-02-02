library(matchingR)
library(data.table)
library(dplyr)

# Import gene scores
atac_gene_scores <- fread("../output/revision-mousebrain_rawGeneScores.tsv") %>% data.frame()
rownames(atac_gene_scores) <- as.character(fread("../output/mousebrain_geneAnnotations.tsv", header = FALSE)[[1]])
mb_anno <- read.table("../data/revision-cluster_annotations.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t", comment.char = "")
mb_trans <- as.character(mb_anno$new); names(mb_trans) <- as.character(mb_anno$old)
colnames(atac_gene_scores) <- mb_trans[as.character(colnames(atac_gene_scores))]
ordering <- as.character(mb_anno$new)
class <- as.character(mb_anno$type)
atac_gene_scores <- atac_gene_scores[,ordering]

# Import Zeisel data
anno_zz <- read.table("../data/zeisel-data/additional-annotations.txt", header = TRUE, sep = "\t")
data_zz <- data.matrix(readRDS("../data/zeisel-data/zeisel_data.rds"))
data_zz <- data_zz[,as.character(anno_zz$Cluster.name)]

# Import marker genes
read.table("../zeisel-data/markergenes-zeisel.txt", header = TRUE, sep = "\t") %>%
  filter(About == "Marker") %>% reshape2::melt(id.vars = c("Cluster.ID", "About")) %>%
  pull(value) %>% unique() -> marker_genes

# Determine overlapping genes
gene_set <- unique(intersect(rownames(data_zz), rownames(atac_gene_scores)))
gene_set2 <- intersect(marker_genes, gene_set)

atac_gene_scores_s <- atac_gene_scores[gene_set2, ]
zz_s <- data_zz[gene_set2,]

# Run gale-shipley
cormat <- cor(atac_gene_scores_s,zz_s, method = "spearman")
results <- matchingR::galeShapley.collegeAdmissions(studentUtils =  t(cormat), collegeUtils = cormat, slots = 1)
optimal_cluster <- colnames(cormat)[results$matched.students]
atac_assign <- data.frame(
  atac_cluster = rownames(cormat),
  optimal_cluster
)

optimal_df <- left_join(atac_assign, anno_zz[,c(2,3,6)], by = c( "optimal_cluster" = "Cluster.name" ))

top5df <- cormat %>% reshape2::melt() %>% group_by(Var1) %>% top_n(5, value) %>% 
  arrange(Var1, desc(value)) %>% mutate(char = paste0(as.character(Var2), "(", as.character(round(value, 3)), ")")) %>%
  dplyr::select(Var1, char) %>% 
  summarize(Top5clusters = paste(char, collapse = '; '))

total <- left_join(optimal_df, top5df, by = c("atac_cluster" = "Var1"))
saveRDS(total, file = "../output/Zeisel_contribution.rds")

