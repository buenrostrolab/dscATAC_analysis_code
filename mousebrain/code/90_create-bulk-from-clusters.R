library(BuenColors)
library(chromVAR)
library(dplyr)
library(reshape2)
library(mc2d)
library(reshape2)
library(Matrix)
library(mgatk)
library(GenomicRanges)

source("../all_heme/helper_fns.R")

# Create new synthetic bigwigs
cdf <- readRDS("mousebrain-master_dataframe.rds")
mb_anno <- read.table("revision-reanno.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t", comment.char = "")
mb_trans <- as.character(mb_anno$new); names(mb_trans) <- as.character(mb_anno$old)
cdf$new_cluster <- mb_trans[paste0("V", as.character(cdf$clusters))]

clusters <- sort(unique(cdf$new_cluster))

# Import
g <- c(
  readRDS("rds_bam/N701_Exp119_sample1_S1-bap-bam.rds"),
  readRDS("rds_bam/N702_Exp119_sample2_S1-bap-bam.rds"),
  readRDS("rds_bam/N704_Exp119_sample4_S1-bap-bam.rds"),
  readRDS("rds_bam/N707_Exp119_sample7_S1-bap-bam.rds"),
  readRDS("rds_bam/N721_Exp110_sample1_combined_S1-bap-bam.rds"),
  readRDS("rds_bam/N722_Exp110_sample2_combined_S1-bap-bam.rds"),
  readRDS("rds_bam/N723_Exp110_sample3_combined_S1-bap-bam.rds"),
  readRDS("rds_bam/N724_Exp110_sample4_combined_S1-bap-bam.rds"),
  readRDS("rds_bam/N726_Exp110_sample5_combined_S1-bap-bam.rds"),
  readRDS("rds_bam/N727_Exp110_sample6_combined_S1-bap-bam.rds"),
  readRDS("rds_bam/N728_Exp110_sample7_combined_S1-bap-bam.rds"),
  readRDS("rds_bam/N729_Exp110_sample8_combined_S1-bap-bam.rds")
)

# Bulk new synthetic bulk
sapply(clusters, function(cluster){
  alignment <- (g[g$DB %in% cdf[cdf$new_cluster == cluster,"DropBarcode"]])
  reads_coverage <- coverage(alignment)/length(alignment)*1000000
  export.bw(reads_coverage, con = paste0("cluster_bigwigs/", as.character(cluster), ".bw"))
  cluster
}) -> bulk2


