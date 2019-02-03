library(BuenColors)
library(chromVAR)
library(dplyr)
library(reshape2)
library(mc2d)
library(reshape2)
library(Matrix)
library(GenomicRanges)
library(GenomicAlignments)
library(rtracklayer)

cdf <- readRDS("baseline_louvain.rds")

g <- c(
  readRDS("../Exp100-others-reseq/rds_bam/Exp100-Sample1-bap-bam.rds"),
  readRDS("../Exp100-others-reseq/rds_bam/Exp100-Sample2-bap-bam.rds"),
  readRDS("../Exp100-others-reseq/rds_bam/Exp100-Sample3-bap-bam.rds"),
  readRDS("../Exp100-others-reseq/rds_bam/Exp100-Sample4-bap-bam.rds"),
  readRDS("../Exp100-others-reseq/rds_bam/Exp100-Sample5-bap-bam.rds"),
  readRDS("../Exp100-others-reseq/rds_bam/Exp100-Sample6-bap-bam.rds"),
  readRDS("../Exp100-others-reseq/rds_bam/Exp100-Sample7-bap-bam.rds"),
  readRDS("../Exp100-others-reseq/rds_bam/Exp100-Sample8-bap-bam.rds"),
  readRDS("../Exp100-others-reseq/rds_bam/Exp100-Sample9-bap-bam.rds"),
  readRDS("../Exp100-others-reseq/rds_bam/Exp100-Sample10-bap-bam.rds"),
  readRDS("../Exp100-others-reseq/rds_bam/Exp100-Sample11-bap-bam.rds"),
  readRDS("../Exp100-others-reseq/rds_bam/Exp100-Sample12-bap-bam.rds"),
  readRDS("../Exp100-others-reseq/rds_bam/Exp100-Sample13-bap-bam.rds"),
  readRDS("../Exp100-others-reseq/rds_bam/Exp100-Sample14-bap-bam.rds"),
  readRDS("../Exp100-others-reseq/rds_bam/Exp100-Sample15-bap-bam.rds"),
  readRDS("../Exp100-others-reseq/rds_bam/Exp100-Sample16-bap-bam.rds")
)

clusters <- unique(cdf$new_cluster)
sapply(clusters, function(cluster){
  alignment <- (g[g$DB %in% cdf[cdf$new_cluster == cluster,"DropBarcode"]])
  reads_coverage <- coverage(alignment)/length(alignment)*1000000
  export.bw(reads_coverage, con = paste0("cluster_bigwigs/", as.character(cluster), ".bw"))
  cluster
}) -> bulk2

