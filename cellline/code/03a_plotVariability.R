library(BuenColors)
library(dplyr)

var <- readRDS("../data/computedVariability_GM12878_SureCell.rds")
varDF <- data.frame(name = var$name, variability = var$variability) %>% arrange(desc(variability))
varDF$rank <- 1:dim(varDF)[1]
p1 <- ggplot(varDF, aes(rank, variability)) + geom_point(size = 0.2) + pretty_plot(fontsize = 6) + L_border() +
  labs(x = "transcription factor rank", y = "accessibility variability")
  
cowplot::ggsave(p1, height = 1.8, width = 1.8, filename = "../output/variability.pdf")

