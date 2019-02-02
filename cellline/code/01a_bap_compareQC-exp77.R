library(data.table)
library(dplyr)
library(BuenColors)

getSumStats <- function(file, what){
  print(file)
  df <- data.frame(fread(file, header = TRUE))
  df$Exp <- what
  df$nucRate <- df$totalNuclearFrags/(df$totalNuclearFrags + df$totalMitoFrags)
  return(df)
  
}

rbindlist(
  list(getSumStats("../data/N714_Exp77_sample11_S1.QCstats.csv", "Sample-11"),
       getSumStats("../data/N715_Exp77_sample12_S1.QCstats.csv", "Sample-12"),
       getSumStats("../data/N716_Exp77_sample13_S1.QCstats.csv", "Sample-13"),
       getSumStats("../data/N718_Exp77_sample14_S1.QCstats.csv", "Sample-14"),
       getSumStats("../data/N719_Exp77_sample15_S1.QCstats.csv", "Sample-15")
  ), fill = TRUE) %>% data.frame() -> sumstatsSM

keep <- c("Exp",  "nucRate", "librarySize", "tssProportion", "medianInsertSize", "uniqueNuclearFrags")
keep %in% colnames(sumstatsSM)

mdf <-  sumstatsSM[,keep]
mdf %>% group_by(Exp) %>% mutate(ncells = n()) %>%summarise_all(mean)

p1 <- ggplot(mdf, aes(x = Exp, y = tssProportion*100)) + geom_boxplot(outlier.shape = NA) +
  pretty_plot(fontsize = 6) + theme(legend.position = "bottom") + labs(x = "", color = "", y = "% TSS fragments") +
  coord_cartesian(ylim = c(0, 100)) + scale_y_continuous(expand = c(0,0)) + 
  scale_color_manual(values = "black") + L_border()+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(mdf, aes(x = Exp, y = log10(librarySize))) + geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(2, 6))  + scale_y_continuous(expand = c(0,0)) + 
  pretty_plot(fontsize = 6) + theme(legend.position = "bottom") + labs(x = "", color = "", y = "log10 library size") +
  scale_color_manual(values = "black") + L_border()+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

p3 <- ggplot(mdf, aes(x = Exp, y = nucRate*100)) + geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 100)) + scale_y_continuous(expand = c(0,0)) + 
  pretty_plot(fontsize = 6) + theme(legend.position = "bottom") + labs(x = "", color = "", y = "% nuclear fragments") +
  scale_color_manual(values = "black") + L_border() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

p4 <- ggplot(mdf, aes(x = Exp, y = medianInsertSize)) + geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(50, 275)) + scale_y_continuous(expand = c(0,0)) + 
  pretty_plot(fontsize = 6) + theme(legend.position = "bottom") + labs(x = "", color = "", y = "Median cell insert size") +
  scale_color_manual(values = "black") + L_border() + theme(axis.text.x = element_text(angle = 45, hjust = 1))


cowplot::ggsave(cowplot::plot_grid(p1, p2, nrow =1),
                height = 1.75, width = 3.5, filename = "../output/Experiment77-view.pdf")

