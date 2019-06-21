library(BuenColors)

df <- readRDS("../rdsin/allROCstuff.rds")

if(FALSE){
  df <- rbind(df, 
              data.frame(
                modnames = c("Spearman", "Spearman", "Pearson", "Pearson", "Jaccard_Peak", "Jaccard_Peak"),
                dsids = 1,
                curvetypes = c("ROC", "PRC", "ROC", "PRC", "ROC", "PRC"), 
                aucs = 0, 
                Sample = "Sample6"
              ))
}
df$stuff <- factor(ifelse(df$curvetypes == "PRC", "P-R Curve", "ROC-Curve"), levels = c("P-R Curve", "ROC-Curve"))
p1 <- ggplot(df, aes(x = Sample, y = aucs, fill = modnames)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.7) + facet_wrap( ~ rev(stuff)) + 
  pretty_plot() + L_border() + scale_fill_manual(values = c("black", "firebrick", "dodgerblue3", "purple3")) +
  theme(legend.position = "bottom") + labs(x = "", y = "AUC", fill = "Method") +
  scale_y_continuous(expand = c(0, 0))

cowplot::ggsave(p1, file = "../figures_out/AUC_results.pdf", width = 4, height = 2.5)
