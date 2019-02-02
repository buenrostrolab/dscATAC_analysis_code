library(GenomicRanges)
library(SummarizedExperiment)
library(data.table)
library(dplyr)
library(BuenColors)
library(Matrix)

#------------
# This file won't run all the way through due to large data on line 48
#------------

gdf <- read.table("../data/mm10.refGenes.2016.1018.csv", stringsAsFactors = FALSE, sep = c(","), header = TRUE)

tss <- data.frame(chr = gdf$chrom, gene = gdf$gene.name, stringsAsFactors = FALSE)
tss$tss <-  gdf$TSS
tss$start <- ifelse(tss$tss - 100000 > 0, tss$tss - 100000, 0)
tss$stop <- tss$tss + 100000

tss_idx <- makeGRangesFromDataFrame(tss, keep.extra.columns = TRUE)

# Import ATAC peaks
adf <- data.frame(fread("../data/mouseBrain_peaks.bed")); colnames(adf) <- c("chr", "start", "end")
adf$mp <- (adf$start + adf$end)/2
atacgranges <- makeGRangesFromDataFrame(adf, start.field = "mp", end.field = "mp")

#Overlap between ATAC peaks and Ranges linker
ov <- findOverlaps(atacgranges, tss_idx)

# Do distance decay for the weights
dist <- abs(mcols(tss_idx)$tss[subjectHits(ov)] - start(atacgranges)[queryHits(ov)])
exp_dist_model <- exp(-1*dist/5000)

# Prep an outcome matrix
m <- Matrix::sparseMatrix(i = c(queryHits(ov), length(atacgranges)),
                          j = c(subjectHits(ov), length(tss_idx)),
                          x = c(exp_dist_model,0))
colnames(m) <- gdf$gene.name
weights <- m[,which(Matrix::colSums(m) != 0)]

# import counts
get_gene_scores_norm <- function(counts){
  atac_gene_scores <- t(weights) %*% counts
  cpm2 <- t(t(atac_gene_scores)/colSums(atac_gene_scores)) * 1000000
  return(cpm2)
}

# Import summarized experiment
ff <- list.files(".", recursive = TRUE, pattern = ".SE.rds")
ffiles <- ff[grepl("_bap/final/", ff) & grepl("Exp11", ff)]
bigSE <- do.call("cbind", lapply(ffiles, readRDS))

# Subset to cells passing filter
cdf <- readRDS("../data/mousebrain-master_dataframe.rds")
bigSE2 <- bigSE[,colData(bigSE)$DropBarcode %in% as.character(cdf$DropBarcode)]
counts <- assays(bigSE2)[["counts"]]

# Get some statistics
colData(bigSE2) %>% data.frame() %>% pull(uniqueNuclearFrags) %>% mean()
colData(bigSE2) %>% data.frame() %>% pull(FRIP) %>% mean()
colSums(counts > 0) %>% mean()

# Define high / low genes
junb_df <- read.table("../data/EN01_junb_scores.tsv", header = TRUE)
high <- as.character(junb_df[junb_df$Junb > 0,"DropBarcode"])
low <- as.character(junb_df[junb_df$Junb < 0,"DropBarcode"])

# Compute the gene scores
counts_two <- data.matrix(
 data.frame(
   high = rowSums(counts[,colData(bigSE2)$DropBarcode %in% high]),
   low = rowSums(counts[,colData(bigSE2)$DropBarcode %in% low])
 ) 
)

gs <- get_gene_scores_norm(counts_two)

# Import genes from FD
fd_genes <- as.character(read.table("IEG_ICs.txt")[,1])
keep_genes <- intersect(rownames(gs), fd_genes)

gs_df <- data.frame(data.matrix(gs[keep_genes,]))
gs_df$diff <- gs_df$high - gs_df$low
gs_df$gene <- rownames(gs_df)
gs_df$mean <- rowMeans(gs_df[,c("high", "low")])
gs_df$percent_difference <- gs_df$diff / gs_df$mean

gs_df %>% arrange(desc(percent_difference)) %>% 
  mutate(rank = 1:n(), up = diff>= 0, class = "IEG") -> gs_df2

gs_df_all <- data.frame(data.matrix(gs)) %>% mutate(diff = high - low, class = "all")

ggplot(gs_df, aes(x = mean, y = percent_difference * 100, label = gene)) +
  geom_text() +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "Mean Expression", y = "% difference") 

pd <- ggplot(rbind(gs_df2[,c("diff", "class")], gs_df_all[,c("diff", "class")]), aes(x = diff, color = class)) +
  geom_density() + scale_color_manual(values = c("black", "red")) + 
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none") +
  labs(y = "Empirical density", x = "Difference in gene score (Junb high - Junb low)")
cowplot::ggsave(pd, file = "../output/IEG_gene_density.pdf", width = 2, height = 2)


ks.test(gs_df2$diff, gs_df_all$diff)



