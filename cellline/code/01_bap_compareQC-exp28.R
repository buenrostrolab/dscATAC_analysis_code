library(data.table)
library(dplyr)
library(BuenColors)

getSumStats <- function(file, what, n = 500){
  print(file)
  df <- data.frame(fread(file, header = TRUE))
  df$Exp <- what
  df$nucRate <- df$totalNuclearFrags/(df$totalNuclearFrags + df$totalMitoFrags)
  return(head(df, n))
  
}

rbindlist(
  list(getSumStats("../data/N701_Exp28_Sample1_S1.QCstats.csv", "1X-Tn5"),
       getSumStats("../data/N705_Exp28_Sample5_S1.QCstats.csv", "3X-Tn5"),
       getSumStats("../data/N707_Exp28_Sample7_S1.QCstats.csv", "4X-Tn5"),
       getSumStats("../data/N711_Exp28_Sample9_S1.QCstats.csv", "Custom-Rep1"),
       getSumStats("../data/N714_Exp28_Sample11_S1.QCstats.csv", "Custom-Rep2"),
       getSumStats("../data/N715_Exp28_Sample12_S1.QCstats.csv", "Custom-Rep3")
  ), fill = TRUE) %>% data.frame() -> sumstatsSM

keep <- c("Exp",  "nucRate", "librarySize", "tssProportion", "medianInsertSize", "uniqueNuclearFrags")
keep %in% colnames(sumstatsSM)

mdf <-  sumstatsSM[,keep]
mdf %>% group_by(Exp) %>% mutate(ncells = n()) %>%summarise_all(mean)

p2 <- ggplot(mdf, aes(x = Exp, y = log10(librarySize))) + geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(2, 6))  + scale_y_continuous(expand = c(0,0)) + 
  pretty_plot(fontsize = 10) + theme(legend.position = "bottom") + labs(x = "", color = "", y = "log10 library size") +
  scale_color_manual(values = "black") + L_border()+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

p3 <- ggplot(mdf, aes(x = Exp, y = nucRate*100)) + geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 100)) + scale_y_continuous(expand = c(0,0)) + 
  pretty_plot(fontsize = 10) + theme(legend.position = "bottom") + labs(x = "", color = "", y = "% nuclear fragments") +
  scale_color_manual(values = "black") + L_border() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

mdf$Distal <- (1-mdf$tssProportion) * mdf$uniqueNuclearFrags
mdf$TSS <- (mdf$tssProportion) * mdf$uniqueNuclearFrags
odf <- reshape2::melt(mdf[,c("Exp", "Distal", "TSS")], id.vars = "Exp")

p4 <- ggplot(mdf, aes(x = Exp, y = log10(TSS))) + geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(2, 5.2))  + scale_y_continuous(expand = c(0,0)) + 
  pretty_plot(fontsize = 10) + theme(legend.position = "bottom") + labs(x = "", color = "", y = "log10 TSS reads") +
  scale_color_manual(values = "black") + L_border()+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

p5 <- ggplot(mdf, aes(x = Exp, y = log10(Distal))) + geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(2, 5.2))  + scale_y_continuous(expand = c(0,0)) + 
  pretty_plot(fontsize = 10) + theme(legend.position = "bottom") + labs(x = "", color = "", y = "log10 Distal reads") +
  scale_color_manual(values = "black") + L_border()+ theme(axis.text.x = element_text(angle = 45, hjust = 1))


cowplot::ggsave(cowplot::plot_grid(p2, p3, p4, p5, nrow =2),
                height = 6, width = 5, filename = "../output/Experiment28-view.pdf")

