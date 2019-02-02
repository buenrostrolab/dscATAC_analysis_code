library(dplyr)

ff <- list.files(".", recursive = TRUE, pattern = ".SE.rds")
ffiles <- ff[grepl("_bap/final/", ff) & grepl("Exp11", ff)]

bigSE <- do.call("cbind", lapply(ffiles, readRDS))
cdf <- readRDS("mousebrain-master_dataframe.rds")
clusters <- sort(unique(cdf$cluster))
counts <- assays(bigSE)[["counts"]]

sapply(clusters, function(cluster){
  rowSums(counts[,cdf$DropBarcode[cdf$cluster == cluster]])
}) -> mat


write.table(data.frame(data.matrix(mat)),
            file = "revision_mouseBrain_Countsmatrix.tsv", sep = "\t", quote = FALSE, row.names = FALSE, 
            col.names = TRUE)

if(FALSE){
  write.table(data.frame(rowRanges(bigSE))[,c(1,2,3)],
              file = "mouseBrain_peaks.bed", sep = "\t", quote = FALSE, row.names = FALSE, 
              col.names = TRUE)
}