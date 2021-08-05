library(plyr)
library(networkD3)

# load labels
df_labels <- read.csv("Data/labels.csv", header = TRUE, stringsAsFactors = FALSE)
# load all data
df_data <- read.delim("Data/fate.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# fate is [0,3] we need [3,6] for sankey diagram
df_data$Fate <- df_data$Fate + 3

# gsk data
df_gsk <- subset(df_data, startsWith(as.character(Experiment),"sec61"))
# parental data
df_parental <- subset(df_data, startsWith(as.character(Experiment),"RPE1Control"))
# control data
df_control <- subset(df_data, startsWith(as.character(Experiment),"RPE1_"))

count_gsk <- ddply(df_gsk, c("Class","Fate"), summarise, Value = length(Cell))
count_parental <- ddply(df_parental, c("Class","Fate"), summarise, Value = length(Cell))
count_control <- ddply(df_control, c("Class","Fate"), summarise, Value = length(Cell))

my_color <- 'd3.scaleOrdinal() .domain(["Aligned", "Free", "Ensheathed", "Normal", "Micronuclei", "Defects", "Death"]) .range(["gray", "blue" , "orange", "black", "black", "black", "black"])'

sankey_list <- list(df_labels,count_gsk)
sankeyNetwork(Links = sankey_list[[2]], Nodes = sankey_list[[1]], Source = "Class",
              Target = "Fate", Value = "Value", NodeID = "Label",
              colourScale = my_color,
              fontSize = 24, nodeWidth = 30, fontFamily = "Helvetica")
rstudioapi::executeCommand("viewerSaveAsWebPage")

sankey_list <- list(df_labels,count_parental)
sankeyNetwork(Links = sankey_list[[2]], Nodes = sankey_list[[1]], Source = "Class",
              Target = "Fate", Value = "Value", NodeID = "Label",
              colourScale = my_color,
              fontSize = 24, nodeWidth = 30, fontFamily = "Helvetica")
rstudioapi::executeCommand("viewerSaveAsWebPage")

sankey_list <- list(df_labels,count_control)
sankeyNetwork(Links = sankey_list[[2]], Nodes = sankey_list[[1]], Source = "Class",
              Target = "Fate", Value = "Value", NodeID = "Label",
              colourScale = my_color,
              fontSize = 24, nodeWidth = 30, fontFamily = "Helvetica")
rstudioapi::executeCommand("viewerSaveAsWebPage")
