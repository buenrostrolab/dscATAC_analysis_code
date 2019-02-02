library(BuenColors)
library(dplyr)
library(data.table)

base <- "https://s3.us-east-2.amazonaws.com/jasonbuenrostro/2018_mouse_brain/"

# Here we define the max accessibility of the ACTB promoter for normalization
def_height <- c(2.57, 3.11, 2.58, 2.24, 1.99, # EN1-5
                2.86, 2.29, 1.83, 2.30, 2.30, # EN6-10
                2.13, 2.81, 4.14, 3.15, 2.06, # EN11-15
                4.15, 2.66, 7.65, 3.15, 2.84, # EN16-IN03
                2.48, 2.46, 9.56, 6.73, 7.23, #IN04, 05, MG, E1, A1
                12.00, 7.00) # OG1, P1

map_df <- read.table("../input/revisions-labels-MB.tsv", header = TRUE, comment.char = "", sep = "\t")
ordering <- as.character(map_df$new)
colVec <- as.character(map_df$colVec)
rgbmat <- t(col2rgb(colVec))

vals <- 1:27

lapply(vals, function(i_int){
  anno <- ordering[i_int]
  priority <- as.character(which(anno == ordering))
  data.frame(val1 =  paste0("track ", anno, "\n",
                            "bigDataUrl ", base, as.character(anno), ".bw\n",
                            "shortLabel ", anno, "\n",
                            "longLabel ", anno, "\n",
                            "type bigWig\n",
                            "viewLimits 0:", as.character(def_height[i_int]), "\n",
                            "viewLimitsMax 0:2000\n", 
                            "maxHeightPixels 100:32:8\n",
                            "visibility full\n",
                            "color ", paste(rgbmat[i_int,], collapse = ","), "\n",
                            "priority ",priority, "\n"))
}) %>% rbindlist() %>% data.frame() -> odf

write.table(odf, file = "../output/trackDb.txt", sep = "", quote = FALSE, col.names = FALSE, row.names = FALSE)
