library(tidyverse)
library(cowplot)

# script to generate combined plots of the data for LD051 and LD058

##############################################################################################################

push_combinations <- function(datadir, label, norm) {
  combine_raw(datadir, label, norm, "chr_border")
  combine_raw(datadir, label, norm, "chr_interior")
  combine_raw(datadir, label, norm, "plt_border")
  combine_raw(datadir, label, norm, "plt_interior")
}

combine_raw <- function(datadir, label, norm, subclass) {
  all_ratio_file_names <- list.files(datadir, pattern = ".txt$")
  # initialise empty data frame
  combined_df <- data.frame(Time_offset = integer(),
                            Average_mean_Ch1 = numeric(),
                            Average_mean_Ch2 = numeric(),
                            Average_mean_Ch3 = numeric(),
                            Movie = character())
  
  for(file_numb in 1:length(all_ratio_file_names)){
    each_file_name <- all_ratio_file_names[file_numb]
    # skip if file is missing
    if(identical(each_file_name, character(0)) == TRUE) {
      next
    } else if(str_detect(each_file_name,"ratio")) {
      next
    } else if(!str_detect(each_file_name, subclass)) {
      next
    }
    each_file_path <- file.path(datadir, each_file_name)
    df <- read.delim(each_file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    this_df <- df[ , c("Time_offset", "Average_mean_Ch1", "Average_mean_Ch2", "Average_mean_Ch3")]
    this_df$Movie <- paste0("file_",toString(file_numb))
    if(norm == 1) {
      this_df <- normalise_the_df(this_df)
    }
    combined_df <- rbind(combined_df, this_df)
  }
  #return(combined_df)
  make_plots_raw(combined_df, label, subclass)
  make_mean_plots_raw(combined_df, label, subclass)
}

# this function will normalise the values to the anaphase frame (time = 0)
# if there is no value at time = 0, the first non-NA value is used for normalisation instead
normalise_the_df <- function(df) {
  rowForNorm <- which(df[,1] == 0)
  normVal <- df[rowForNorm,2]
  if(is.na(normVal) || identical(rowForNorm, integer(0)) == TRUE) {
    nonNAindex <- which(!is.na(df[,2]))
    rowForNorm <- min(nonNAindex)
    normVal <- df[rowForNorm,2]
  }
  df[,2] <- df[,2] / normVal
  rowForNorm <- which(df[,1] == 0)
  normVal <- df[rowForNorm,3]
  if(is.na(normVal) || identical(rowForNorm, integer(0)) == TRUE) {
    nonNAindex <- which(!is.na(df[,3]))
    rowForNorm <- min(nonNAindex)
    normVal <- df[rowForNorm,3]
  }
  df[,3] <- df[,3] / normVal
  rowForNorm <- which(df[,1] == 0)
  normVal <- df[rowForNorm,4]
  if(is.na(normVal) || identical(rowForNorm, integer(0)) == TRUE) {
    nonNAindex <- which(!is.na(df[,4]))
    rowForNorm <- min(nonNAindex)
    normVal <- df[rowForNorm,4]
  }
  df[,4] <- df[,4] / normVal
  
  return(df)
}

make_plots_raw <- function(df, label, subclass) {
  sub_plot_raw(df, label, 1, subclass)
  sub_plot_raw(df, label, 2, subclass)
  sub_plot_raw(df, label, 3, subclass)
}

sub_plot_raw <- function(df,label,channel,subclass) {
  # df has time as first column and then channels as following columns
  datacol <- channel + 1
  # get summary stats for df
  df_summary <- data.frame(n = tapply(df[,datacol], df$Time_offset, length),
                           mean = tapply(df[,datacol], df$Time_offset, mean, na.rm = TRUE),
                           sd  = tapply(df[,datacol], df$Time_offset, sd, na.rm = TRUE))
  df_summary$time <- as.numeric(row.names(df_summary))
  # make plot
  p <- ggplot(df, aes_string(x = "Time_offset", y = paste0("Average_mean_Ch",toString(channel)), col = "Movie")) +
    geom_line() +
    geom_line(data = df_summary, aes(x=time, y=mean), colour = "black", size = 1, alpha = 0.8) +
    xlab("Time") +
    ylab("Intensity") +
    ggtitle(paste(label,"Ch",toString(channel))) +
    theme(legend.position='none')
  
  file_path <- paste0("Output/Plots/combine_raw/all_", label, "_", subclass, "_ch",toString(channel),".png")
  ggsave(file_path, plot = p, dpi = 300)
}

make_mean_plots_raw <- function(df, label, subclass) {
  # get summary stats for df
  df_summary <- data.frame(n_1 = tapply(df[,2], df$Time_offset, length),
                           mean_1 = tapply(df[,2], df$Time_offset, mean, na.rm = TRUE),
                           sd_1  = tapply(df[,2], df$Time_offset, sd, na.rm = TRUE),
                           n_2 = tapply(df[,3], df$Time_offset, length),
                           mean_2 = tapply(df[,3], df$Time_offset, mean, na.rm = TRUE),
                           sd_2  = tapply(df[,3], df$Time_offset, sd, na.rm = TRUE))
  df_summary$time <- as.numeric(row.names(df_summary))
  # calculate error bars - use s.e.m.
  df_summary$lower_1 <- df_summary$mean_1 - (df_summary$sd_1 / sqrt(df_summary$n_1))
  df_summary$upper_1 <- df_summary$mean_1 + (df_summary$sd_1 / sqrt(df_summary$n_1))
  df_summary$lower_2 <- df_summary$mean_2 - (df_summary$sd_2 / sqrt(df_summary$n_2))
  df_summary$upper_2 <- df_summary$mean_2 + (df_summary$sd_2 / sqrt(df_summary$n_2))
  # make plot
  p <- ggplot(df_summary, aes(x = time)) +
    geom_line(data = df_summary, aes(x = time, y = mean_1), colour = "#ed1C24", size = 1, alpha = 0.8) +
    geom_ribbon(aes(ymin = df_summary$lower_1, ymax = df_summary$upper_1), alpha = 0.2, fill = "#ed1C24", colour = NA) +
    geom_line(data = df_summary, aes(x = time, y = mean_2), colour = "#00a651", size = 1, alpha = 0.8) +
    geom_ribbon(aes(ymin = df_summary$lower_2, ymax = df_summary$upper_2), alpha = 0.2, fill = "#00a651", colour = NA) +
    xlim(-40, 40) +
    ylim(0.5, 3) +
    xlab("Time") +
    ylab("Intensity") +
    theme(legend.position = "none")
  
  file_path <- paste0("Output/Plots/combine_raw/mean_", label, "_", subclass,".png")
  ggsave(file_path, plot = p, dpi = 300)
  
  # save these dataframes so that we can combine
  output_path <- paste0("Output/Data/mean")
  ifelse(!dir.exists(output_path), dir.create(output_path), FALSE)
  output_name <- paste0(output_path, "/", label, "_mean_", subclass, "_raw_averages.txt")
  write.table(df_summary, file = output_name, sep="\t", row.names=FALSE, col.name=TRUE)
}

make_figure <- function(title, chrPath, pltPath, theYMin, theYMax) {
  thePath <- paste0("Output/Data/mean/",chrPath)
  df_chr <- read.delim(thePath, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  thePath <- paste0("Output/Data/mean/",pltPath)
  df_plt <- read.delim(thePath, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  # make plot
  p1 <- ggplot(df_chr, aes(x = time)) +
    geom_hline(yintercept = 1, colour = "grey", linetype = "dotted") +
    geom_vline(xintercept = 0, colour = "grey", linetype = "dotted") +
    geom_line(data = df_plt, aes(x = time, y = mean_1), colour = "#a00309", size = 0.5, alpha = 0.8) +
    geom_ribbon(aes(ymin = df_plt$lower_1, ymax = df_plt$upper_1), alpha = 0.2, fill = "#a00309", colour = NA) +
    geom_line(data = df_chr, aes(x = time, y = mean_1), colour = "#ed1C24", size = 0.5, alpha = 0.8) +
    geom_ribbon(aes(ymin = df_chr$lower_1, ymax = df_chr$upper_1), alpha = 0.2, fill = "#ed1C24", colour = NA) +
    xlim(-20, 40) +
    ylim(theYMin, theYMax) +
    xlab("Time (min)") +
    ylab("Intensity (norm.)") +
    theme_cowplot(font_size = 9, line_size = 0.4,) +
    theme(legend.position = "none")
  p2 <- ggplot(df_chr, aes(x = time)) +
    geom_hline(yintercept = 1, colour = "grey", linetype = "dotted") +
    geom_vline(xintercept = 0, colour = "grey", linetype = "dotted") +
    geom_line(data = df_plt, aes(x = time, y = mean_2), colour = "#00592c", size = 0.5, alpha = 0.8) +
    geom_ribbon(aes(ymin = df_plt$lower_2, ymax = df_plt$upper_2), alpha = 0.2, fill = "#00592c", colour = NA) +
    geom_line(data = df_chr, aes(x = time, y = mean_2), colour = "#00a651", size = 0.5, alpha = 0.8) +
    geom_ribbon(aes(ymin = df_chr$lower_2, ymax = df_chr$upper_2), alpha = 0.2, fill = "#00a651", colour = NA) +
    xlim(-20, 40) +
    ylim(0.8, 2) +
    xlab("Time (min)") +
    ylab("Intensity (norm.)") +
    theme_cowplot(font_size = 9, line_size = 0.4) +
    theme(legend.position = "none")
  
  # save these dataframes so that we can combine
  output_path <- paste0("Output/Plots/summary")
  ifelse(!dir.exists(output_path), dir.create(output_path), FALSE)
  # file_path <- paste0("Output/Plots/summary/", title,".png")
  # ggsave(file_path, plot = p, dpi = 300)
  file_path <- paste0("Output/Plots/summary/", title,"1.pdf")
  ggsave(file_path, plot = p1, width = 36, height = 40, units = "mm")
  file_path <- paste0("Output/Plots/summary/", title,"2.pdf")
  ggsave(file_path, plot = p2, width = 36, height = 40, units = "mm")
}

##############################################################################################################

# check for output path
output_path <- "Output/Plots/combine_raw"
ifelse(!dir.exists(output_path), dir.create(output_path), FALSE)
# run the analysis on the outputs generated in the previous script
push_combinations("Output/Data/LBR/LBR_lagging", "LBR_lag", 1)
push_combinations("Output/Data/LBR/LBR_polar", "LBR_polar", 1)
push_combinations("Output/Data/BAF/BAF_lagging", "BAF_lag", 1)
push_combinations("Output/Data/BAF/BAF_polar", "BAF_polar", 1)
# use this function to load the averaged data and plot
make_figure("BAF_lag","BAF_lag_mean_chr_interior_raw_averages.txt","BAF_lag_mean_plt_interior_raw_averages.txt",0.8,3)
make_figure("BAF_polar","BAF_polar_mean_chr_interior_raw_averages.txt","BAF_polar_mean_plt_interior_raw_averages.txt",0.8,3)
make_figure("LBR_lag","LBR_lag_mean_chr_interior_raw_averages.txt","LBR_lag_mean_plt_interior_raw_averages.txt",0.8,2)
make_figure("LBR_polar","LBR_polar_mean_chr_interior_raw_averages.txt","LBR_polar_mean_plt_interior_raw_averages.txt",0.8,2)

