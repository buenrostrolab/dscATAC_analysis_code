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

odf <- data.frame(
  CD19=class_props(df, readRDS("../data/projected-pc-scores/CD19.rds")),
  CD34=class_props(df, readRDS("../data/projected-pc-scores/CD34.rds")),
  Mono=class_props(df, readRDS("../data/projected-pc-scores/Mono.rds")),
  CD4=class_props(df, readRDS("../data/projected-pc-scores/CD4.rds")),
  CD8=class_props(df, readRDS("../data/projected-pc-scores/CD8.rds")),
  NK=class_props(df, readRDS("../data/projected-pc-scores/NK.rds")),
  PBMC=class_props(df, readRDS("../data/projected-pc-scores/PBMC.rds"))
)

meltdf <- reshape2::melt(data.matrix(odf))
meltdf$Var2 <- factor(as.character(meltdf$Var2), rev(unique(as.character(meltdf$Var2))))
meltdf<-meltdf %>% group_by(Var2) %>%
  mutate(valmax = max(value)) %>% ungroup() %>% 
  data.frame() %>%
  mutate(value2 = value/valmax)

p1 <- ggplot(meltdf, aes(x = Var1, y = Var2, fill = value2)) + geom_tile() +
  geom_tile()+ scale_fill_gradientn(colors = jdb_palette("brewer_heat")) +
  labs(x = "Louvain clusters", y = "Sorted Populations", fill = "Most alike\nProjection") +
  pretty_plot(fontsize = 8) +
  theme(legend.position = "bottom") + L_border() +
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid = element_blank(), panel.border = element_blank())
cowplot::ggsave(p1, file = "../plots/confusion-heatmap.pdf", width = 4, height = 4)
