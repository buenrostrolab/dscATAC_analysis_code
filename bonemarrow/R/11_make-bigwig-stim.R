library(dplyr)
library(data.table)
"%ni%" <- Negate("%in%")

# Do the Principle components centering / scaling
PCs <- data.frame(fread("scScores.txt"))
meta <- data.frame(fread("tsneControls_with_info.txt"))
rownames(PCs) <- PCs$V1
PCs <- data.matrix(PCs[,-1])

PCs_control <- data.frame(PCs[meta$DropBarcode,])
PCs_control$group <- read.table("cluster_IDs_exp100_baseline.tsv")[,2]

PCs_control %>% group_by(group) %>% summarize_all(median) %>% 
  data.frame() -> df
rownames(df) <- df[,1]; df <- df[,-1]

# Figure out stim barcodes
wdf <- readRDS("vinay/all-exp100.rds")
new_pcs <- wdf[wdf$CellID %ni% meta$DropBarcode,4:19]
A <- df; B <- new_pcs

# Classify
euklDist <- t(sqrt(apply(array(apply(B,1,function(x){(x-t(A))^2}),c(ncol(A),nrow(A),nrow(B))),2:3,sum)))
colnames(euklDist) <- rownames(A)
colnames(euklDist)[max.col(-1*euklDist, 'first')] -> vec
assign_df <- data.frame(bc = rownames(new_pcs), annotation = vec)

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

clusters <- unique(assign_df$annotation)
sapply(clusters, function(cluster){
  alignment <- (g[g$DB %in% assign_df[assign_df$annotation == cluster,"bc"]])
  reads_coverage <- coverage(alignment)/length(alignment)*1000000
  export.bw(reads_coverage, con = paste0("cluster_bigwigs/", as.character(cluster), "_STIM.bw"))
  cluster
}) -> bulk2


