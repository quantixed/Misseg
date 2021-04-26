# 201208_LD051_BAFcloneC9_GSK_001.nd2 (series 5)-400box1_chr1
# 201123_LD051_BAFcloneC9_GSK_001.nd2 - 201123_LD051_BAFcloneC9_GSK_001.nd2 (series 1)-400box1
# 201125_LD058_LBRcloneC7_GSK_001.nd2 - 201125_LD058_LBRcloneC7_GSK_001.nd2 (series 4)-400box2
# 201116_LD058_LBRcloneC7_GSK_001.nd2 (series 3)-400box3_chr1

library(tidyverse)
library(cowplot)

# script to generate example traces for the images in the figure

##############################################################################################################

# uses normalise script from combine_raw.R

make_example_plot <- function(subdir, stem, title) {
  # get list of "average_measures.txt" files from the correct output folder
  thePath <- paste0("Output/Data/",subdir)
  all_file_names <- list.files(thePath, pattern = ".txt$")
  # filter list for our cell of interest
  the_file_names <- grep(stem, all_file_names, fixed = TRUE, value = TRUE)
  # get the file we want (relies on only one in directory)
  chr_name <- grep("chr_interior", the_file_names, fixed = TRUE, value = TRUE)
  plt_name <- grep("plt_interior", the_file_names, fixed = TRUE, value = TRUE)
  path_to_file <- paste0(thePath, "/", chr_name)
  df <- read.delim(path_to_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  df_chr <- df[ , c("Time_offset", "Average_mean_Ch1", "Average_mean_Ch2", "Average_mean_Ch3")]
  path_to_file <- paste0(thePath, "/", plt_name)
  df <- read.delim(path_to_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  df_plt <- df[ , c("Time_offset", "Average_mean_Ch1", "Average_mean_Ch2", "Average_mean_Ch3")]
  # normalise to anaphase
  df_chr <- normalise_the_df(df_chr)
  df_plt <- normalise_the_df(df_plt)
  
  # make plot
  p <- ggplot(df_chr, aes(x = time)) +
    geom_hline(yintercept = 1, colour = "grey", linetype = "dotted") +
    geom_vline(xintercept = 0, colour = "grey", linetype = "dotted") +
    geom_line(data = df_plt, aes(x = Time_offset, y = Average_mean_Ch2), colour = "#00592c", size = 0.5, alpha = 0.8) +
    geom_line(data = df_chr, aes(x = Time_offset, y = Average_mean_Ch2), colour = "#00a651", size = 0.5, alpha = 0.8) +
    geom_line(data = df_plt, aes(x = Time_offset, y = Average_mean_Ch1), colour = "#a00309", size = 0.5, alpha = 0.8) +
    geom_line(data = df_chr, aes(x = Time_offset, y = Average_mean_Ch1), colour = "#ed1C24", size = 0.5, alpha = 0.8) +
    xlim(-20, 40) +
    ylim(0.6, 3.2) +
    geom_segment(aes(x = -20, y = 0.6, xend = 0, yend = 0.6)) +
    geom_segment(aes(x = -20, y = 0.6, xend = -20, yend = 1.6)) +
    xlab("Time (min)") +
    ylab("Intensity (norm.)") +
    theme_cowplot(font_size = 9, line_size = 0.4) +
    theme(legend.position = "none")
  
  # check for output path
  output_path <- "Output/Plots/examples"
  ifelse(!dir.exists(output_path), dir.create(output_path), FALSE)
  # file_path <- paste0("Output/Plots/examples/", title,".png")
  # ggsave(file_path, plot = p, dpi = 300)
  file_path <- paste0("Output/Plots/examples/", title,".pdf")
  ggsave(file_path, plot = p, width = 40, height = 28, units = "mm")
}

##############################################################################################################


baf_lag <- "201208_LD051_BAFcloneC9_GSK_001.nd2 (series 5)-400box1_chr1"
baf_polar <- "201123_LD051_BAFcloneC9_GSK_001.nd2 (series 1)-400box1"
lbr_lag <- "201125_LD058_LBRcloneC7_GSK_001.nd2 (series 4)-400box2"
lbr_polar <- "201116_LD058_LBRcloneC7_GSK_001.nd2 (series 3)-400box3_chr1"

make_example_plot("BAF/BAF_lagging", baf_lag, "baf_lagging")
make_example_plot("BAF/BAF_polar", baf_polar, "baf_polar")
make_example_plot("LBR/LBR_lagging", lbr_lag, "lbr_lagging")
make_example_plot("LBR/LBR_polar", lbr_polar, "lbr_polar")

