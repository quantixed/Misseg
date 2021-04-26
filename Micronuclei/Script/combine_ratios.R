library(tidyverse)

# script to generate combined plots of the data for LD051 and LD058

##############################################################################################################

combine_ratios <- function(datadir, label) {
  all_ratio_file_names <- list.files(datadir, pattern = ".txt$")
  # initialise empty data frame
  combined_df <- data.frame(Time_offset = integer(),
                            Average_mean_Ch1_ratio = numeric(),
                            Average_mean_Ch2_ratio = numeric(),
                            Average_mean_Ch3_ratio = numeric(),
                            Movie = character())
  
  for(file_numb in 1:length(all_ratio_file_names)){
    each_file_name <- all_ratio_file_names[file_numb]
    # skip if file is missing
    if(identical(each_file_name, character(0)) == TRUE) {
      next
    } else if(!str_detect(each_file_name,"ratio")) {
      next
    } else if(!str_detect(each_file_name,"border")) {
      next
    }
    each_file_path <- file.path(datadir, each_file_name)
    df <- read.delim(each_file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    this_df <- df[ , c("Time_offset", "Average_mean_Ch1_ratio", "Average_mean_Ch2_ratio", "Average_mean_Ch3_ratio")]
    this_df$Movie <- paste0("file_",toString(file_numb))
    combined_df <- rbind(combined_df, this_df)
  }
  #return(combined_df)
  make_plots(combined_df, label)
}

make_plots <- function(df, label) {
  sub_plot(df,label,1)
  sub_plot(df,label,2)
  sub_plot(df,label,3)
}

sub_plot <- function(df,label,channel) {
  # df has time as first column and then channels as following columns
  datacol <- channel + 1
  # get summary stats for df
  df_summary <- data.frame(n = tapply(df[,datacol], df$Time_offset, length),
                           mean = tapply(df[,datacol], df$Time_offset, mean, na.rm = TRUE),
                           sd  = tapply(df[,datacol], df$Time_offset, sd, na.rm = TRUE))
  df_summary$time <- as.numeric(row.names(df_summary))
  # make plot
  p <- ggplot(df, aes_string(x = "Time_offset", y = paste0("Average_mean_Ch",toString(channel),"_ratio"), col = "Movie")) +
    geom_line() +
    geom_line(data = df_summary, aes(x=time, y=mean), colour = "black", size = 1, alpha = 0.8) +
    xlab("Time") +
    ylab("Ratio (Chromosome/Main Mass)") +
    scale_y_continuous(trans = "log2",
                       limits=(c(0.25,4))) +
    ggtitle(paste(label,"Ch",toString(channel))) +
    theme(legend.position='none')
  
  file_path <- paste0("Output/Plots/combine_ratios/all_", label, "_ch",toString(channel),".png")
  ggsave(file_path, plot = p, dpi = 300)
}

##############################################################################################################

# check for output path
output_path <- "Output/Plots/combine_ratios"
ifelse(!dir.exists(output_path), dir.create(output_path), FALSE)
# run the analysis on the outputs generated in the previous script
combine_ratios("Output/Data/LBR/LBR_lagging", "LBR_lag")
combine_ratios("Output/Data/LBR/LBR_polar", "LBR_polar")
combine_ratios("Output/Data/BAF/BAF_lagging", "BAF_lag")
combine_ratios("Output/Data/BAF/BAF_polar", "BAF_polar")
