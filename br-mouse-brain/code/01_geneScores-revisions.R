library(GenomicRanges)
library(SummarizedExperiment)
library(data.table)
library(dplyr)
library(BuenColors)
library(Matrix)

set.seed(1)

# Import gene bodies; restrict to TSS
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

# Make a plot
p1 <- qplot(table(subjectHits(ov)), binwidth = 1) + theme(plot.subtitle = element_text(vjust = 1), 
                                                          plot.caption = element_text(vjust = 1)) +
  labs(title = "Histogram of peaks per gene",  x = "Peaks / gene") + pretty_plot()
#cowplot::ggsave(p1, file = "../output/histo1-br.png", width = 5, height = 4)

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
counts <- data.matrix(data.frame(fread("../data/revision_mouseBrain_Countsmatrix.tsv")))
geneScores <- t(weights) %*% counts

write.table(data.frame(data.matrix(geneScores)), row.names = FALSE, col.names = FALSE, 
            sep = "\t", quote = FALSE, file = "../output/revision-mousebrain_rawGeneScores.tsv")
write.table(data.frame(rownames(geneScores)), row.names = FALSE, col.names = FALSE, 
            sep = "\t", quote = FALSE, file = "../output/revision-mousebrain_geneAnnotations.tsv")



                      