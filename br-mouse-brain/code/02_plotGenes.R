library(data.table)
library(BuenColors)
library(dplyr)
library(cowplot)

# Gene scores
atac_gene_scores <- fread("../output/revision-mousebrain_rawGeneScores.tsv") %>% data.frame()
rownames(atac_gene_scores) <- as.character(fread("../output/mousebrain_geneAnnotations.tsv", header = FALSE)[[1]])
mb_anno <- read.table("../data/revision-cluster_annotations.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t", comment.char = "")
mb_trans <- as.character(mb_anno$new); names(mb_trans) <- as.character(mb_anno$old)
colnames(atac_gene_scores) <- mb_trans[as.character(colnames(atac_gene_scores))]

# Make a counts per million matrix
cpm <- t(t(atac_gene_scores)/colSums(atac_gene_scores)) * 1000000
colVec <- as.character(mb_anno$colVec); names(colVec) <- as.character(mb_anno$new)


#-------
# define plot genes here
#-------
genes <- c("Rbfox3", "Zfp36l2", "Slc17a7", "Gad2", "C1qb",
            "Flt1", "Aqp4", "Mog", "Cspg4")

# Subset genes and 
cpm[genes,] %>%
  data.frame() %>% data.matrix() %>% 
  reshape2::melt() %>%
  mutate(cluster = factor(as.character(Var2), levels = as.character(mb_anno$new))) %>% 
  arrange(cluster) -> plot_data_frame

# Function for making one gene plot
make_gene_score_plot <- function(gene) {
  ggplot(plot_data_frame %>% filter(Var1 == gene),aes(x = cluster, y = value, fill = cluster)) +
    geom_bar(color = "black", stat = "identity") +
    scale_fill_manual(values = colVec) +
     pretty_plot(fontsize = 6) + L_border() +
    labs(x = "", y = "Gene score" ) +
    scale_y_continuous(expand = c(0,0)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(gene) + theme(legend.position = "none")
}

# Mass export
cowplot::ggsave(plot_grid(
  plotlist = lapply(genes, make_gene_score_plot), ncol = 3),
  height = 5, width = 7.5, file = "../output/genes_scores.pdf"
)
