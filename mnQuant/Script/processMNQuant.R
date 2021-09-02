# process outputs from NucleiSeg.ijm
library(tidyverse)
library(openxlsx)
library(cowplot)
library(ggforce)

# function definition
get_bg_subtracted_values <- function(datadir, fname, result_matrix, ii, mn, pn){
  stemname = gsub(fname, pattern=".tif$", replacement="")
  # if(!file.exists(paste0(datadir,"/bg_",stemname,".csv"))) {
  #   print(stemname)
  #   return(result_matrix)
  # }
  bg_df <- read.csv(file = paste0(datadir,"/bg_",stemname,".csv"))
  c1_df <- read.csv(file = paste0(datadir,"/Q_quant_C1-",stemname,".csv"))
  c2_df <- read.csv(file = paste0(datadir,"/Q_quant_C2-",stemname,".csv"))
  c3_df <- read.csv(file = paste0(datadir,"/Q_quant_C3-",stemname,".csv"))
  c4_df <- read.csv(file = paste0(datadir,"/Q_quant_C4-",stemname,".csv"))
  
  temp_df <- subset(c1_df, Label == mn)
  result_matrix[ii,1] = temp_df$Mean - bg_df$Mean[1]
  temp_df <- subset(c2_df, Label == mn)
  result_matrix[ii,2] = temp_df$Mean - bg_df$Mean[2]
  temp_df <- subset(c3_df, Label == mn)
  result_matrix[ii,3] = temp_df$Mean - bg_df$Mean[3]
  temp_df <- subset(c4_df, Label == mn)
  result_matrix[ii,4] = temp_df$Mean - bg_df$Mean[4]
  
  temp_df <- subset(c1_df, Label == pn)
  result_matrix[ii,5] = temp_df$Mean - bg_df$Mean[1]
  temp_df <- subset(c2_df, Label == pn)
  result_matrix[ii,6] = temp_df$Mean - bg_df$Mean[2]
  temp_df <- subset(c3_df, Label == pn)
  result_matrix[ii,7] = temp_df$Mean - bg_df$Mean[3]
  temp_df <- subset(c4_df, Label == pn)
  result_matrix[ii,8] = temp_df$Mean - bg_df$Mean[4]
  
  return(result_matrix)
}

## from https://github.com/jdmanton/morgenstemning
ametrine <- function(n=256, mincolor=NULL, maxcolor=NULL, invert=FALSE, alpha=1) {
  controlPoints <- matrix(c(
    30, 60, 150,  # cyan
    180, 90, 155, # purple
    230, 85, 65,  # red-ish
    220, 220, 0), # yellow
    ncol=3, byrow=T)
  controlPoints <- controlPoints / 255
  
  k <- c(1, 17, 32, 64)
  
  # For non-isoluminent points, a normal interpolation gives better results
  cmap <- apply(controlPoints, 2, function(y) approx(x=1:length(y), y, xout=seq(1, 4, length.out=64))$y)
  
  # Linearly interpolate to the required number of points
  cmap <- apply(cmap, 2, function(y) approx(x=1:dim(cmap)[1], y, n=n)$y)
  
  # Flip colormap if required
  if(invert) cmap <- cmap[dim(cmap)[1]:1, , drop=FALSE]
  
  # Convert RGB values to hex codes and replace min/max colors
  if(alpha != 1){
    colors <- rgb(cmap, alpha=alpha)
  } else {
    colors <- rgb(cmap)
  }
  if(!is.null(mincolor)) colors[1] <- mincolor
  if(!is.null(maxcolor)) colors[length(colors)] <- maxcolor
  
  colors
}

process_files <- function(annotationsPath, xLabels){
  groupName <- basename(dirname(annotationsPath))
  datadir <- dirname(annotationsPath)
  # load in the annotated data
  annotations <- read.xlsx(annotationsPath, sheet = 1, startRow = 1, skipEmptyRows = FALSE, skipEmptyCols = FALSE)
  all_data <- subset(annotations, error == 0)
  # 4 channels and for each we need mn, pn and the ratio = 3
  all_meas <- matrix(NA, nrow = nrow(all_data), ncol = 4 * 3)
  
  for (i in 1:nrow(all_data)) {
    imagename <- all_data$filename[i]
    mn <- all_data$mn[i]
    pn <- all_data$pn[i]
    all_meas <- get_bg_subtracted_values(datadir, imagename, all_meas, i, mn, pn)
  }
  
  all_meas[,9:12] <- all_meas[,1:4] / all_meas[,5:8]
  headers <- c(names(all_data),"mnc1","mnc2","mnc3","mnc4","pnc1","pnc2","pnc3","pnc4","ratioc1","ratioc2","ratioc3","ratioc4")
  all_meas <- cbind.data.frame(all_data,all_meas)
  names(all_meas) <- headers
  
  write.csv(all_meas, file = paste0("Output/Data/",groupName,"_processed_data.csv"), row.names = F)
  
  # plotting
  df <- gather(all_meas[,13:16], channel, ratio)
  p1 <- ggplot(df, aes(x = channel, y = ratio)) +
    geom_hline(yintercept = 1, linetype = 'dashed', colour = 'black') +
    geom_sina(alpha = 0.5, stroke = 0) +
    stat_summary(fun.data = mean_se, geom = 'point', size=2, aes(group=channel)) +
    stat_summary(fun.data = mean_sdl, fun.args = list(mult=1), geom = 'errorbar', size=0.8, aes(group=channel), width=0) +
    scale_y_continuous(trans = 'log2', limits = c(0.0078125, 128)) +
    scale_x_discrete(labels = xLabels) +
    xlab("") +
    ylab("MN/PN (Log2)") +
    theme_cowplot(font_size = 9, line_size = 0.4) +
    theme(legend.position = "none")
  plotPath = paste0("Output/Plots/",groupName,"_ratio.pdf")
  ggsave(plotPath, plot = p1, width = 120, height = 80, units = "mm")
  
  logdf <- as.data.frame(log2(all_meas[,13:16]))
  # remove any rows with NAs
  logdf <- logdf[complete.cases(logdf), ]
  
  p2 <- ggplot(logdf, aes(x = ratioc3, y = ratioc2, color = ratioc4)) +
    geom_vline(xintercept = 0, colour = "grey", linetype = "dashed") +
    geom_hline(yintercept = 0, colour = "grey", linetype = "dashed") +
    geom_point(alpha = 1, shape = 21, stroke = 0.25, size = 1) +
    scale_color_gradientn(colours = ametrine(256), limits = c(-7,7) ) +
    xlab("GFP-Sec61b") +
    ylab(xLabels[2]) +
    lims(x = c(-2.5,5), y = c(-2.5,5)) +
    coord_fixed() +
    theme_cowplot(font_size = 9)
  plotPath = paste0("Output/Plots/",groupName,"_scatter.pdf")
  ggsave(plotPath, plot = p2 + theme(legend.position = "none"), width = 50, height = 42, units = "mm")
  plotPath = paste0("Output/Plots/",groupName,"_scatter_lgnd.pdf")
  ggsave(plotPath, plot = p2, width = 50, height = 42, units = "mm")
}

process_files(annotationsPath = "Data/LBR/Annotations.xlsx", xLabels = c("DNA","LBR-mCherry","GFP-Sec61b","H3K27ac"))
process_files(annotationsPath = "Data/BAF/Annotations.xlsx", xLabels = c("DNA","mCherry-BAF","GFP-Sec61b","H3K27ac"))

