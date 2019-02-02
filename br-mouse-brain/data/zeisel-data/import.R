library(Matrix)
library(LoomExperiment)

LD <- import("l5_all.agg.loom", type="SingleCellLoomExperiment")
mat <-  assays(LD)[["matrix"]]
rownames(mat) <- LD@rowRanges@elementMetadata[,"Gene"]
colnames(mat) <- LD@colData[,"ClusterName"]

saveRDS(mat, file = "zeisel_data.rds")
