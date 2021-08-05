library(plyr)
library(networkD3)

# load labels
df_labels <- read.csv("Data/chrlabels.csv", header = TRUE, stringsAsFactors = FALSE)
# load all data
df_data <- read.delim("Data/chrfate.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

count_data <- ddply(df_data, c("chromosomeClass","chromosomeFate"), summarise, Value = length(chromosomeClass))

sankey_list <- list(df_labels,count_data)
my_color <- 'd3.scaleOrdinal() .domain(["Free", "Ensheathed", "Aligned", "Rescue", "Micronuclei", "Defects", "Death"]) .range(["blue", "orange" , "gray", "green", "black", "black", "black"])'
sankeyNetwork(Links = sankey_list[[2]], Nodes = sankey_list[[1]], Source = "chromosomeClass",
              Target = "chromosomeFate", Value = "Value", NodeID = "Label",
              colourScale = my_color,
              fontSize = 24, nodeWidth = 30, fontFamily = "Helvetica")
rstudioapi::executeCommand("viewerSaveAsWebPage")
