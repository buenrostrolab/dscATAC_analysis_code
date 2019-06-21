library(BuenColors)
library(dplyr)
library(precrec)
library(data.table)
library(ggrastr)
source("00_helper.R")

sample <- "Sample1"

thresholds_frag <- c(0.00975973384731702, 0.0107673626829162, 0.00275636620379365, 0.00280180330297872, 0.0124750065703549)
names(thresholds_frag) <- c("Sample1", "Sample2", "Sample3", "Sample4", "Sample6")

thresholds_tag <- c(0.01, 0.01, 0.005, 0.005, 0.005)
names(thresholds_tag) <- c("Sample1", "Sample2", "Sample3", "Sample4", "Sample6")

lapply(c("Sample1", "Sample2", "Sample3", "Sample4", "Sample6"), function(sample){
  message(sample)
  
  df <- readRDS(paste0("../rdsin/",sample,".all.rds"))
  cor(df[,c("jaccard_tags", "spearman", "pearson", "jaccard_peak", "jaccard_frag")])
  
  # Make a data frame for plotting with ticks
  jaccard_threshold_frag <- thresholds_frag[sample]
  jaccard_threshold_tag <- thresholds_tag[sample]
  
  df <- df %>%
    arrange(desc(jaccard_tags)) %>%
    mutate(rank_tag = 1:n()) %>%
    arrange(desc(jaccard_frag)) %>%
    mutate(rank_frag = 1:n(),
           passKnee = jaccard_frag > jaccard_threshold_frag,
           passTag = jaccard_tags > jaccard_threshold_tag) 
  
  # Tick Marks
  ticks_at_y_frag <- get_ticks(df, 'jaccard_frag')
  ticks_at_x_frag <- get_ticks(df, 'rank_frag')
  x_intercept_frag <- sum(df$passKnee)
  
  ticks_at_y_tag <- ticks_at_y_frag
  ticks_at_x_tag <- get_ticks(df, 'rank_tag')
  x_intercept_tag <- sum(df$passTag)
  
  
  ###
  # Make tag plot
  ###
  
  jaccard_tag_plot <- df %>% head(1000000) %>% 
    ggplot(aes(x=rank_tag, y=(jaccard_tags + 0.00001))) +
    scale_y_log10(breaks = ticks_at_y_tag, labels = as.numeric(ticks_at_y_tag)) +
    scale_x_log10(breaks = ticks_at_x_tag, labels = as.integer(ticks_at_x_tag)) +
    pretty_plot() + L_border() +
    xlab("jaccard tag overlap in rank-descending order") +
    ylab("jaccard tag per barcode pair") +
    theme(axis.text.x = element_text(angle=90)) +
    geom_point_rast(aes(color = passKnee)) +
    scale_color_manual(values = c("black", "dodgerblue3")) +
    theme(legend.position ="none") +
    labs(color = "Pass Knee") 
  
  # Also save as a png for good measure
  cowplot::ggsave(jaccard_tag_plot, file = paste0("../figures_out/", sample, "_tag-log10.png"),
                  width = 6, height = 6)
  
  jaccard_frag_plot <- df %>% head(1000000) %>% 
    ggplot(aes(x=rank_frag, y=jaccard_frag)) +
  #  scale_y_log10(breaks = ticks_at_y_frag, labels = as.numeric(ticks_at_y_frag)) +
    scale_x_log10(breaks = ticks_at_x_frag, labels = as.integer(ticks_at_x_frag)) +
    pretty_plot() + L_border() +
    xlab("jaccard frag overlap in rank-descending order") +
    ylab("jaccard frag per barcode pair") +
    theme(axis.text.x = element_text(angle=90)) +
    geom_point_rast(aes(color = passKnee)) +
    scale_color_manual(values = c("black", "dodgerblue3")) +
    theme(legend.position ="none") +
    labs(color = "Pass Knee") 
  
  # Also save as a png for good measure
  cowplot::ggsave(jaccard_frag_plot, file = paste0("../figures_out/", sample, "_frag-log10.png"),
                  width = 6, height = 6)
  
  
  scores1 <- join_scores(df$spearman, df$pearson, df$jaccard_peak, df$jaccard_frag)
  mmpr <- mmdata(scores1, df$passTag, modnames = c("Spearman", "Pearson", "Jaccard_Peak", "bap"))
  mscurves <- evalmod(mmpr)
  
  dfo <- auc(mscurves)
  dfo$Sample <- sample
  dfo
}) %>% rbindlist() %>% data.frame() -> allOut

#saveRDS(allOut, file = "../rdsin/allROCstuff.rds")

