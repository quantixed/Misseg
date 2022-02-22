library(plyr)
library(networkD3)

# load labels
df_labels <- read.csv("Data/labels.csv", header = TRUE, stringsAsFactors = FALSE)
# load all data
df_data <- read.delim("Data/lbr_baf.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# fate is [0,3] we need [3,6] for sankey diagram
df_data$Fate <- df_data$Fate + 3

# lbr data
df_lbr <- subset(df_data, startsWith(as.character(Experiment),"LBR"))
# baf data
df_baf <- subset(df_data, startsWith(as.character(Experiment),"BAF"))

count_lbr <- ddply(df_lbr, c("Class","Fate"), summarise, Value = length(Cell))
count_baf <- ddply(df_baf, c("Class","Fate"), summarise, Value = length(Cell))

my_color <- 'd3.scaleOrdinal() .domain(["Aligned", "Free", "Ensheathed", "Normal", "Micronuclei", "Defects", "Death"]) .range(["gray", "blue" , "orange", "black", "black", "black", "black"])'

sankey_list <- list(df_labels,count_lbr)
sankeyNetwork(Links = sankey_list[[2]], Nodes = sankey_list[[1]], Source = "Class",
              Target = "Fate", Value = "Value", NodeID = "Label",
              colourScale = my_color,
              fontSize = 24, nodeWidth = 30, fontFamily = "Helvetica")
rstudioapi::executeCommand("viewerSaveAsWebPage")

sankey_list <- list(df_labels,count_baf)
sankeyNetwork(Links = sankey_list[[2]], Nodes = sankey_list[[1]], Source = "Class",
              Target = "Fate", Value = "Value", NodeID = "Label",
              colourScale = my_color,
              fontSize = 24, nodeWidth = 30, fontFamily = "Helvetica")
rstudioapi::executeCommand("viewerSaveAsWebPage")
