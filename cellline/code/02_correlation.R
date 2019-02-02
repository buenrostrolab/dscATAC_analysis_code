library(BuenColors)

counts <- do.call("cbind", lapply(list.files("../data/readsInPeaks/", full.names = TRUE), read.table, header = FALSE))
colnames(counts) <-  gsub("K562_", "k", gsub("GM12878_", "g", gsub("-RIP.txt", "", list.files("../readsInPeaks/"))))
counts <- counts
cpm <- counts/colSums(counts) * 1000000
corMat <- round(cor(log2(counts + 1), method = "spearman"), 3)

ordered_vals <- c("gSureCell", "gBuenrostro", "gCusanovich",  "gPliner",  "gPreissl", "gDNase", "gOMNI", "kSureCell","kBuenrostro", "kDNase", "kOMNI")
melted_cormat <- reshape2::melt(corMat)

# Reorder
melted_cormat$Var1 <- factor(as.character(melted_cormat$Var1), ordered_vals)
melted_cormat$Var2 <- factor(as.character(melted_cormat$Var2), rev(ordered_vals))
p1 <- ggplot(data = melted_cormat[melted_cormat$value < 1.1,], aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+ scale_fill_gradientn(colors = jdb_palette("solar_extra")) +
  labs(x = "", y = "", fill = "Reads in Peaks\n(Spearman Correlation)") + pretty_plot(fontsize = 6) +
  theme(legend.position = "bottom") + L_border() +
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid = element_blank(), panel.border = element_blank())
cowplot::ggsave(p1, file = "../output/correlation_matrix.pdf", height = 2.6, width = 2.2)

#counts$density <- get_density(log2(counts$GM12878_C1+1), log2(counts$GM12878_DNase+1))
#ggplot(counts, aes(x = log2(GM12878_C1 +1), y = log2(GM12878_DNase +1), color = density)) + 
#         geom_point() + scale_color_gradientn(colors = jdb_palette("brewer_red"))


