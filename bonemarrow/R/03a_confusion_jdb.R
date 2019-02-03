library(dplyr)
library(data.table)
"%ni%" <- Negate("%in%")

# Do the Principle components centering / scaling
PCs <- data.frame(fread("../data/scScores.txt"))
meta <- data.frame(fread("../data/tsneControls_with_info.txt"))
rownames(PCs) <- PCs$V1
PCs <- data.matrix(PCs[,-1])

PCs_control <- data.frame(PCs[meta$DropBarcode,])
PCs_control$group <- read.table("../data/cluster_IDs_exp100_baseline.tsv")[,2]

PCs_control %>% group_by(group) %>% summarize_all(median) %>% 
  data.frame() -> df
rownames(df) <- df[,1]; df <- df[,-1]

# Function to get the proportion of each isolate
class_props <- function(A, B){
  
  # Get the pair-wise euclidean distance msot-alike
  euklDist <- t(sqrt(apply(array(apply(B,1,function(x){(x-t(A))^2}),c(ncol(A),nrow(A),nrow(B))),2:3,sum)))
  colnames(euklDist) <- rownames(A)
  colnames(euklDist)[max.col(-1*euklDist, 'first')] %>% table() -> vec
  
  # Find those not in most alike and append zeros
  niboo <- colnames(euklDist) %ni% names(vec)
  vec2 <- c(as.vector(vec), rep(0, sum(niboo))); names(vec2) <- c(names(vec), colnames(euklDist)[niboo])
  # Return in order to facilitate sapply
  round(vec2[colnames(euklDist)]/sum(vec2), 3)*100
}

jdb <- readRDS("../data/projected-pc-scores/CD34JDB.rds")

stringr::str_split_fixed(gsub("BM0106-|BM0828-|-frozen|PB1022-|-160106|-160107|160809-|160818-|160808-|160105-|Frozen-|160819-|160822-|BM1077-|20160617-|BM1214-|LS-|20160726-|singles-|HYC-|BM1137-|scATAC-|.st.bam", "",
                              rownames(jdb)) , pattern="-", n=2)[,1] %>%
  gsub(pattern = "GMP1low", replacement = "GMP-A") %>% 
  gsub(pattern = "GMP2mid", replacement = "GMP-B") %>% 
  gsub(pattern = "GMP3high", replacement = "GMP-C") %>%
  gsub(pattern = "UNK", replacement = "GMP") %>%
  gsub(pattern = "MCP", replacement = "pDC") -> label


odf <- data.frame(
  CLP=class_props(df, jdb[label == "CLP",]),
  CMP=class_props(df, jdb[label == "CMP",]),
  GMP=class_props(df, jdb[label == "GMP",]),
  GMPA=class_props(df, jdb[label == "GMP-A",]),
  GMPB=class_props(df, jdb[label == "GMP-B",]),
  GMPC=class_props(df, jdb[label == "GMP-C",]),
  HSC=class_props(df, jdb[label == "HSC",]),
  LMPP=class_props(df, jdb[label == "LMPP",]),
  MEP=class_props(df, jdb[label == "MEP",]),
  mono=class_props(df, jdb[label == "mono",]),
  MPP=class_props(df, jdb[label == "MPP",]),
  pDC=class_props(df, jdb[label == "pDC",])
)


