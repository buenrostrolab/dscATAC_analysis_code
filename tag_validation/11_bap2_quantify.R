library(data.table)
library(precrec)
library(dplyr)

# Thresholds that I used in bap1 for the tag thresholds
thresholds_tag <- c(0.01, 0.01, 0.005, 0.005, 0.005)
names(thresholds_tag) <- c("Sample1", "Sample2", "Sample3", "Sample4", "Sample6")

# ZB: update this path
path_to_csvgz_files <- "../../may24_2019_from_ZB/"

# Simple function that takes the raw file name + sample name and coputed AUROC/AUPRCs
compute_metrics <- function(raw_file, sample){
  dt <- fread(paste0(path_to_csvgz_files, "/", raw_file))
  pass <- dt[["jaccard_tags"]] > thresholds_tag[sample]
  mmpr <- mmdata(dt[["jaccard_frag"]], pass, modnames = c("bap2"))
  mscurves <- evalmod(mmpr)
  dfo <- auc(mscurves)
  dfo$Sample <- sample
  dfo
}

raw_files <- c("N701_Exp69_sample1.implicatedBarcodesWithTags.csv.gz", "N702_Exp69_sample2.implicatedBarcodesWithTags.csv.gz", 
               "N703_Exp69_sample3.implicatedBarcodesWithTags.csv.gz", "N704_Exp69_sample4.implicatedBarcodesWithTags.csv.gz", 
               "N706_Exp69_sample6.implicatedBarcodesWithTags.csv.gz")

# Loop over all 
lapply(1:5, function(i){
  compute_metrics(raw_files[i],names(thresholds_tag)[i])
}) %>% rbindlist()
