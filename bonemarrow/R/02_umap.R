library(data.table)
library(umap)
library(BuenColors)

meta <- data.frame(fread("../data/tsneControls_with_info.txt"))

# Create the embedding form the reference principle components; only need to do this one
if(FALSE){
  PCs <- data.frame(fread("../data/scScores.txt"))
  rownames(PCs) <- PCs$V1
  PCs <- data.matrix(PCs[,-1])
  
  PCs_control <- PCs[meta$DropBarcode,]
  
  umap.config <- umap.defaults
  umap.config$n_neighbors <- 20
  umap.config$min_dist <- 1
  bm.umap <- umap(PCs_control, config = umap.config)
  
  umap_df <- data.frame(meta, bm.umap$layout)
  ggplot(shuf(umap_df), aes(x = X1 , y = X2 , color =  MostAlikeP)) +
    geom_point(size = 0.1) + labs(color = "Most alike") +
    scale_color_manual(values = c(ejc_color_maps)) +
    pretty_plot() + L_border() + labs(x = "UMAP1", y = "UMAP2") + ggtitle("Exp100") 
  saveRDS(bm.umap, file = "../data/bm.umap.rds")
}

bm.umap <- readRDS("../data/bm.umap.rds")
umap_df <- data.frame(meta, bm.umap$layout)

cd19 <- predict(bm.umap, readRDS("../data/projected-pc-scores/CD19.rds"))
cd4 <- predict(bm.umap, readRDS("../data/projected-pc-scores/CD4.rds"))
cd8 <- predict(bm.umap, readRDS("../data/projected-pc-scores/CD8.rds"))
nk <- predict(bm.umap, readRDS("../data/projected-pc-scores/NK.rds"))
mono <- predict(bm.umap, readRDS("../data/projected-pc-scores/Mono.rds"))
cd34 <- predict(bm.umap, readRDS("../data/projected-pc-scores/CD34.rds"))

pdf <- data.frame(
  umap1 = c(umap_df$X1, cd19[,1], cd4[,1], cd8[,1], nk[,1], mono[,1], cd34[,1]),
  umap2 = c(umap_df$X2, cd19[,2], cd4[,2], cd8[,2], nk[,2], mono[,2], cd34[,2]),
  names = c(rownames(umap_df), rownames(cd19), rownames(cd4), rownames(cd8), rownames(nk), rownames(mono), rownames(cd34)),
  projected = c(rep("base", length(umap_df$X1)),
                rep("B", length(cd19[,1])), rep("CD4", length(cd4[,1])), rep("CD8", length(cd8[,1])), rep("NK", length(nk[,1])), rep("Mono", length(mono[,1])), rep("CD34", length(cd34[,1])))
)

pbr <- ggplot((pdf[c(1:60495, sample(60495:dim(pdf)[1])),]), aes(x = umap1 , y = umap2 , color = projected)) +
  geom_point_rast(size = 0.1) + labs(color = "Projected") +
  scale_color_manual(values = c("base" = "grey", ejc_color_maps, "CD34" = "black")) +
  pretty_plot() + L_border() + labs(x = "UMAP1", y = "UMAP2") + ggtitle("Exp100") 
cowplot::ggsave(pbr, file = "../plots/allExp100-exp47.pdf", width = 4.3, height = 3.3)

# Define clusters
if(FALSE){
  library(igraph)
  library(FNN)
  knn <- get.knn(PCs_control, algo="kd_tree", k = 3)[["nn.index"]]
  igraphObj <- igraph::graph_from_adjacency_matrix(igraph::get.adjacency(igraph::graph.edgelist(data.matrix(reshape2::melt(knn)[,c("Var1", "value")]), directed=FALSE)), mode = "undirected")
  clusters <- cluster_louvain(igraphObj)
  table(membership(clusters))
  
  ldf <- data.frame(meta, clusters = paste0("X", as.character(membership(clusters))))
  
  trans <- c("B-like-1", "B-like-2", "CD4-like", "CMP-like", "pDC-like", "CD4-like",
             "CLP-like-1", "Late-Ery","CD4-like", "mDC-like", "CD4-like", "CD8-like",
             "Mono-like", "CD34-like", "CLP-like-2", "NK-like", "Early-Ery", "Collision")
  names(trans) <- paste0("X", as.character(1:18))
  ldf$new_cluster <- trans[paste0( as.character(ldf$clusters))]
  saveRDS(ldf, file = "baseline_louvain.rds")
  
  colors <- c("purple4", "purple2", "dodgerblue2", "pink", "purple", "navyblue", "firebrick",
              "orange", "dodgerblue4", "orange3", "green4", "dodgerblue3", "purple4", "red", "grey")
  names(colors) <- unique(trans)
  
  
  ggplot(shuf(ldf), aes(x = tSNE1 , y = tSNE2 , color =  new_cluster, label = new_cluster)) +
    geom_point(size = 0.1) + labs(color = "Louvain") +
    scale_color_manual(values = colors) +
    pretty_plot() + L_border() + labs(x = "tSNE-1", y = "tSNE-2") + ggtitle("Exp100") 
  
  bm.umap <- readRDS("../data/bm.umap.rds")
  umap_df <- data.frame(bm.umap$layout, cluster = ldf$new_cluster)
  
  trans <- c("B-like-1", "B-like-2", "CD4-like", "CMP-like", "pDC-like", "CD4-like",
             "CLP-like-1", "Late-Ery","CD4-like", "mDC-like", "CD4-like", "CD8-like",
             "Mono-like", "CD34-like", "CLP-like-2", "NK-like", "Early-Ery", "Collision")
  colors <- c("purple4", "purple2", "dodgerblue2", "pink", "purple", "navyblue", "firebrick",
              "orange", "dodgerblue4", "orange3", "green4", "dodgerblue3", "purple4", "red", "grey")
  names(colors) <- unique(trans)
  
  pLou <- ggplot(shuf(umap_df), aes(x = X1 , y = X2 , color =  cluster)) +
    geom_point_rast(size = 0.1) + labs(color = "Louvain") +
    scale_color_manual(values = colors) +
    pretty_plot() + L_border() + labs(x = "UMAP1", y = "UMAP2") + ggtitle("Louvain") 
  cowplot::ggsave(pLou, file = "plots/allExp100-Louvain.pdf", width = 4.3, height = 3.3)
  
  
}



# Project old data data
jdb <- predict(bm.umap, (readRDS("../data/projected-pc-scores/CD34JDB.rds")))

# Fix cell color annotation
stringr::str_split_fixed(gsub("BM0106-|BM0828-|-frozen|PB1022-|-160106|-160107|160809-|160818-|160808-|160105-|Frozen-|160819-|160822-|BM1077-|20160617-|BM1214-|LS-|20160726-|singles-|HYC-|BM1137-|scATAC-|.st.bam", "",
                              rownames(jdb)) , pattern="-", n=2)[,1] %>%
  gsub(pattern = "GMP1low", replacement = "GMP-A") %>% 
  gsub(pattern = "GMP2mid", replacement = "GMP-B") %>% 
  gsub(pattern = "GMP3high", replacement = "GMP-C") %>%
  gsub(pattern = "UNK", replacement = "GMP") %>%
  gsub(pattern = "MCP", replacement = "pDC") -> label

pdf <- data.frame(
  umap1 = c(umap_df$X1, jdb[,1]),
  umap2 = c(umap_df$X2, jdb[,2]),
  projected = c(rep("base", length(umap_df$X1)),
                label)
)

pCp <- ggplot((pdf), aes(x = umap1 , y = umap2 , color = projected)) +
  geom_point_rast(size = 0.1) + labs(color = "Projected") +
  scale_color_manual(values = c(ejc_color_maps, "base" = "grey", "mono" = "orange")) +
  pretty_plot() + L_border() + labs(x = "UMAP1", y = "UMAP2") + ggtitle("Cell paper") 
cowplot::ggsave(pCp, file = "../plots/allExp100-C1bm.pdf", width = 4.3, height = 3.3)


bm.umap <- readRDS("../data/bm.umap.rds")
umap_df <- data.frame(bm.umap$layout)
pbmc <- predict(bm.umap, (readRDS("../data/projected-pc-scores/PBMC.rds")))

exp <- ifelse(grepl("Exp44",rownames(pbmc),fixed = FALSE),"Exp44",
              ifelse(grepl("Exp45",rownames(pbmc),fixed = FALSE),"Exp45",
                     "Exp92"))

pdf <- data.frame(
  umap1 = c(umap_df$X1, pbmc[,1]),
  umap2 = c(umap_df$X2, pbmc[,2]),
  projected = c(rep("base", length(umap_df$X1)),
                exp)
)

pdf <- pdf[pdf$projected != "Exp92",]

pPBMC <- ggplot((pdf)[c(1:60495, sample(60495:dim(pdf)[1])),], aes(x = umap1 , y = umap2 , color = projected)) +
  geom_point_rast(size = 0.1) + labs(color = "Projected") +
  scale_color_manual(values = c("Exp44" = "dodgerblue", "Exp45" = "firebrick",
                                "Exp92" = "green4",  "base" = "grey")) +
  pretty_plot() + L_border() + labs(x = "UMAP1", y = "UMAP2") + ggtitle("PBMCs - experiment") 
cowplot::ggsave(pPBMC, file = "../plots/allExp100-pbmc.pdf", width = 4.3, height = 3.3)


# Project stimulation 

"%ni%" <- Negate("%in%")
bm.umap <- readRDS("../data/bm.umap.rds")
umap_df <- data.frame(bm.umap$layout)
wdf <- readRDS("../data/projected-pc-scores/all-exp100.rds")
new_pcs <- wdf[wdf$CellID %ni% meta$DropBarcode,4:19]
proj_exp100 <- predict(bm.umap,new_pcs)

umap

pdf <- data.frame(
  umap1 = c(umap_df$X1, proj_exp100[,1]),
  umap2 = c(umap_df$X2, proj_exp100[,2]),
  CellID = c(rownames(umap_df), rownames(proj_exp100))
)

wmdf <- merge(wdf, pdf)

p100 <- ggplot(shuf(wmdf), aes(x = umap1 , y = umap2 , color = Group)) +
  geom_point_rast(size = 0.1) + labs(color = "Projected") +
  scale_color_manual(values = c("LPS" = "dodgerblue", "Serum" = "firebrick",
                                "Control" = "grey")) +
  pretty_plot() + L_border() + labs(x = "UMAP1", y = "UMAP2") + ggtitle("All experiment 100") 

cowplot::ggsave(p100, file = "plots/allExp100-stim.pdf", width = 4.3, height = 3.3)


