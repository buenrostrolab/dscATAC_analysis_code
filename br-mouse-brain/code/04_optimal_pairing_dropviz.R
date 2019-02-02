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

# Import DropViz
anno_dropviz <- readRDS("../data/dropviz/annotation.BrainCellAtlas_Saunders_version_2018.04.01.RDS")
data_dropviz <- readRDS("../data/dropviz/metacells.BrainCellAtlas_Saunders_version_2018.04.01.RDS")

# Determine overlapping genes
gene_set <- unique(intersect(rownames(data_dropviz), rownames(atac_gene_scores)))
gene_set2 <- intersect(c(as.character(anno_dropviz$class_marker), as.character(anno_dropviz$type_marker)), gene_set)

atac_gene_scores_s <- atac_gene_scores[gene_set2, ]
dropviz_s <- data_dropviz[gene_set2,]
cormat_all <- cor(atac_gene_scores_s, dropviz_s, method = "spearman")
stopifnot(all(as.character(anno_dropviz$tissue_subcluster) == colnames(cormat_all)))
colnames(cormat_all) <- as.character(anno_dropviz$full_name)

# Function that maps an optimal ATAC cluster
findOptimalMatch_cluster <- function(atac_class, rna_class, rna_gene_start = ""){
  
  atac_samples <- colnames(atac_gene_scores_s)[class == atac_class]
  
  if(rna_gene_start != ""){
    anno_dropviz[grepl(rna_gene_start, anno_dropviz$class_marker),] %>% filter(class %in% rna_class)  %>%
      pull(tissue_subcluster) -> rna_samples
    
  } else {
    anno_dropviz %>% filter(class %in% rna_class) %>%
      pull(tissue_subcluster) -> rna_samples
  }
  
  if(length(atac_samples) > 1){
    cormat <- cor(atac_gene_scores_s[,atac_samples], dropviz_s[,rna_samples], method = "spearman")
    results <- matchingR::galeShapley.collegeAdmissions(studentUtils =  t(cormat), collegeUtils = cormat, slots = 1)
    optimal_cluster <- colnames(cormat)[results$matched.students]
  } else {
    cormat <- cor(atac_gene_scores_s[,atac_samples], dropviz_s[,rna_samples], method = "spearman")
    optimal_cluster <- colnames(cormat)[which.max(cormat)]
  }
  
  return(data.frame(atac_samples, optimal_cluster))
}

rbind(
  findOptimalMatch_cluster(atac_class = "Excite-NEURON", rna_class = "NEURON", rna_gene_start = "Slc17a"),
  findOptimalMatch_cluster(atac_class = "Inhibit-NEURON", rna_class = "NEURON", rna_gene_start = "Gad"),
  findOptimalMatch_cluster(atac_class = "MICROGLIA", rna_class = "MICROGLIA", rna_gene_start = ""),
  findOptimalMatch_cluster(atac_class = "ENDOTHELIAL", rna_class = c("ENDOTHELIAL_TIP", "ENDOTHELIAL_STALK"), rna_gene_start = ""),
  findOptimalMatch_cluster(atac_class = "ASTROCYTE", rna_class = "ASTROCYTE", rna_gene_start = ""),
  findOptimalMatch_cluster(atac_class = "OLIGODENDROCYTE", rna_class = "OLIGODENDROCYTE", rna_gene_start = ""),
  findOptimalMatch_cluster(atac_class = "POLYDENDROCYTE", rna_class = "POLYDENDROCYTE", rna_gene_start = "")
) -> atac_assign

optimal_df <- left_join(atac_assign, anno_dropviz[,c("tissue_subcluster", "class", "class_marker", "type_marker", "common_name")],
          by = c( "optimal_cluster" = "tissue_subcluster" ))

top5df <- cormat_all %>% reshape2::melt() %>% group_by(Var1) %>% top_n(5, value) %>% 
  arrange(Var1, desc(value)) %>% mutate(char = paste0(as.character(Var2), "(", as.character(round(value, 3)), ")")) %>%
  dplyr::select(Var1, char) %>% 
  summarize(Top5clusters = paste(char, collapse = '; '))

total <- left_join(optimal_df, top5df, by = c("atac_samples" = "Var1"))
saveRDS(total, file = "../output/saunders_contribution.rds")

