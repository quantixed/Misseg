# Quantify the data from LD051 and LD058 (outputs from ImageJ script)

##############################################################################################################
library(ggplot2)
library(reshape2)

# Master function
live_function_general <- function(datadir, label){
  
  # loop through sample folders within data folder
  all_sample_folder_names <- list.dirs(datadir, recursive = FALSE)
  
  # output data will be saved here
  output_path <- paste0("Output/",datadir)
  setup_output_path(output_path)
  
  for(folder_number in 1:length(all_sample_folder_names)){
    pathname_sample_folder <- all_sample_folder_names[folder_number]
    each_sample_folder <- basename(pathname_sample_folder)
    pathname_measure_folder <- paste0(pathname_sample_folder, "/", "Measures")
    all_measure_file_names <- list.files(pathname_measure_folder, pattern = "*.csv")
    # skip if there were no files in this folder
    if(identical(all_measure_file_names, character(0)) == TRUE) {
      next
    }
    
    for(file_numb in 1:length(all_measure_file_names)){
      each_file_name <- all_measure_file_names[file_numb]
      # skip if file is missing
      if(identical(each_file_name, character(0)) == TRUE) {
        next
      }
      each_file_path <- paste0(pathname_measure_folder, "/", each_file_name)
      
      #make objects for each of the different measure sets
      if(grepl("chr_enlarged.csv", each_file_name) == TRUE){
        chr_enlarged <- read.csv(file=each_file_path, header=TRUE)
      } else if(grepl("chr_orig.csv", each_file_name) == T){
        chr_orig <- read.csv(file=each_file_path, header=TRUE)
      } else if(grepl("plt_enlarged.csv", each_file_name) == T){
        plt_enlarged <- read.csv(file=each_file_path, header=TRUE)
      } else if(grepl("plt_orig.csv", each_file_name) == T){
        plt_orig <- read.csv(file=each_file_path, header=TRUE)
      }
    }
    # find minimum and maximum frame numbers 
    min_frame <- min(chr_orig$Frame, plt_orig$Frame)
    max_frame <- max(chr_orig$Frame, plt_orig$Frame)
    # find anaphase frame
    anaphase_filepath <- paste0(pathname_sample_folder, "/Frame_anaphase.txt")
    anaphase <- find_anaphase_frame(anaphase_filepath)
    
    # rawintden of enlarged - rawintden of orig = count of pixel values in the boundary region
    # round(rawintden of enlarged / mean of enlarged) - round(rawintden of orig / mean of orig) = count of pixels in the boundary region
    # pixel values / pixels = mean of boundary
    # we'll overwrite the Mean column with the boundary mean (so that we can use a single function to get the averages)
    chr_enlarged$BoundaryIntDen <- chr_enlarged$RawIntDen - chr_orig$RawIntDen
    chr_enlarged$BoundaryPixels <- round(chr_enlarged$RawIntDen / chr_enlarged$Mean) - round(chr_orig$RawIntDen / chr_orig$Mean)
    chr_enlarged$Mean <- chr_enlarged$BoundaryIntDen / chr_enlarged$BoundaryPixels
    plt_enlarged$BoundaryIntDen <- plt_enlarged$RawIntDen - plt_orig$RawIntDen
    plt_enlarged$BoundaryPixels <- round(plt_enlarged$RawIntDen / plt_enlarged$Mean) - round(plt_orig$RawIntDen / plt_orig$Mean)
    plt_enlarged$Mean <- plt_enlarged$BoundaryIntDen / plt_enlarged$BoundaryPixels
    
    # run subfunction to output file of average mean measures of all slices measured for a frame/channel
    # assign as dataframe here too
    output_chr_border <- live_sub(chr_enlarged, paste0(label,"_chr_border"), min_frame, max_frame, anaphase, datadir, each_sample_folder)
    output_plt_border <- live_sub(plt_enlarged, paste0(label,"_plt_border"), min_frame, max_frame, anaphase, datadir, each_sample_folder)
    output_chr_interior <- live_sub(chr_orig, paste0(label,"_chr_interior"), min_frame, max_frame, anaphase, datadir, each_sample_folder)
    output_plt_interior <- live_sub(plt_orig, paste0(label,"_plt_interior"), min_frame, max_frame, anaphase, datadir, each_sample_folder)

    # make a column for the Region (chromosome or Plate)
    output_chr_border$Region <- "Chr"
    output_plt_border$Region <- "Plt"
    output_chr_interior$Region <- "Chr"
    output_plt_interior$Region <- "Plt"
    
    # run subfunction to generate ratios of the chromosome/plate measure for datasets
    # assign variable name to outputs
    output_ratios_border <- ratios_sub(output_chr_border, output_plt_border, "border", min_frame, max_frame, anaphase, datadir, each_sample_folder)
    output_ratios_interior <- ratios_sub(output_chr_interior, output_plt_interior, "interior", min_frame, max_frame, anaphase, datadir, each_sample_folder)
    
    #plot points and save plot
    ratio_plot_sub(output_ratios_border, each_sample_folder, label, "Border")
    ratio_plot_sub(output_ratios_interior, each_sample_folder, label, "Interior")
  }
}


##############################################################################################################
# helper functions
##############################################################################################################

# function to setup the output directory if it doesn't exist
setup_output_path <- function(fpath) {
  # example "Output/Data/foo/bar"
  path_components <- unlist(strsplit(fpath, split = "/"))
  
  for(i in 1:length(path_components)){
    test_path <- paste(path_components[1:i], collapse = "/")
    ifelse(!dir.exists(test_path), dir.create(test_path),FALSE)
  }
}

# function get anaphase frame number from first line of file
# cannot use read.table etc because files do not necessarily have EOL
find_anaphase_frame <- function(fname) {
  ref <- file(fname,"r")
  first_line <- readLines(ref, n=1, warn = F)
  close(ref)
  refvec <- unlist(strsplit(first_line, split = "\t"))
  return(as.numeric(refvec[2]))
}

# subfunction - making an average measures df and output
live_sub <- function(measure_set, measure_name, min, max, anaphase_frame, datadir, sample_name){
  # each frame is three minutes
  timestep <- 3
  #loop from minimum to maximum frame number
  for(frame_numb in min:max){
    #subset all slice measures at that frame for each channel
    #take average of the measure across slices
    all_slices <- subset(measure_set, Frame==frame_numb)
    all_slices_Ch1 <- subset(all_slices, Ch==1)
    all_slices_Ch2 <- subset(all_slices, Ch==2)
    all_slices_Ch3 <- subset(all_slices, Ch==3)
    
    #account for case where frame number outwith range of an individual data set
    if(nrow(all_slices) == 0){
      add <- data.frame(frame_numb, NA, NA, NA)
      names(add) <- c("Frame", "Average_mean_Ch1", "Average_mean_Ch2", "Average_mean_Ch3")
    }
    
    if(nrow(all_slices) > 0){
      average_mean_Ch1 <- mean(all_slices_Ch1$Mean)
      average_mean_Ch2 <- mean(all_slices_Ch2$Mean)
      average_mean_Ch3 <- mean(all_slices_Ch3$Mean)
      
      add <- data.frame(frame_numb, average_mean_Ch1, average_mean_Ch2, average_mean_Ch3)
      names(add) <- c("Frame", "Average_mean_Ch1", "Average_mean_Ch2", "Average_mean_Ch3")
    }
    
    #add to an output dataframe for measure set (frames as rows, channels as columns)
    if(frame_numb == min){
      output_df <- data.frame(add)
    }
    if(frame_numb > min){
      output_df <- rbind(output_df, add)
    }
  }
  
  # add time elapsed column and a column offset to anaphase frame
  output_df$Time <- timestep * (output_df$Frame - 1)
  output_df$Time_offset <- timestep * (output_df$Frame - anaphase_frame)

  #output the data frame for the sample measure set constructed before loop to next measure set
  #name the output file
  output_path <- paste0("Output/", datadir)
  ifelse(!dir.exists(output_path), dir.create(output_path), FALSE)
  output_name <- paste0(output_path, "/", sample_name, "_", measure_name, "_average_measures.txt")
  
  #save file of individual region type measure of repeat
  write.table(output_df, file=output_name, sep="\t", row.names=FALSE, col.name=TRUE)  
  return(output_df)
}

# subfunction - making a ratio df and output
ratios_sub <- function(chr_file, plate_file, borderOrInterior, min, max, anaphase_frame, datadir, sample_name){
  # each frame is three minutes
  timestep <- 3
  
  for(frame_number in min:max){
    chr_values <- subset(chr_file, Frame == frame_number)
    plate_values <- subset(plate_file, Frame == frame_number)
    Average_mean_Ch1_ratio <- chr_values$Average_mean_Ch1 / plate_values$Average_mean_Ch1
    Average_mean_Ch2_ratio <- chr_values$Average_mean_Ch2 / plate_values$Average_mean_Ch2
    Average_mean_Ch3_ratio <- chr_values$Average_mean_Ch3 / plate_values$Average_mean_Ch3
    add_ratios <- data.frame(frame_number,Average_mean_Ch1_ratio, Average_mean_Ch2_ratio, Average_mean_Ch3_ratio)
    names(add_ratios) <- c("Frame", "Average_mean_Ch1_ratio", "Average_mean_Ch2_ratio", "Average_mean_Ch3_ratio")
    
    if(frame_number == min){
      output_df_ratios <- data.frame(add_ratios)
    }
    if(frame_number > min){
      output_df_ratios <- rbind(output_df_ratios, add_ratios)
    }
  }
  
  # add time elapsed column and a column offset to anaphase frame
  output_df_ratios$Time <- timestep * (output_df_ratios$Frame - 1)
  output_df_ratios$Time_offset <- timestep * (output_df_ratios$Frame - anaphase_frame)
  
  #output the data frame for the sample measure set constructed before loop to next measure set
  #name the output file
  output_path <- paste0("Output/", datadir)
  ifelse(!dir.exists(output_path), dir.create(output_path), FALSE)
  output_name <- paste0(output_path, "/", sample_name, "_", borderOrInterior, "_ratio_average_measures.txt")
  
  #save file of individual region type measure of repeat
  write.table(output_df_ratios, file=output_name, sep="\t", row.names=FALSE, col.name=TRUE)  
  
  return(output_df_ratios)
}

# subfunction - making plots
ratio_plot_sub <- function(df, sample_name, type, region) {
  df <- melt(df, id.vars="Time_offset", measure.vars = c("Average_mean_Ch1_ratio","Average_mean_Ch2_ratio","Average_mean_Ch3_ratio"))
  
  #plot points and save plot
  p <- ggplot(df, aes(Time_offset, value, col=variable)) +
    geom_point() +
    geom_line() +
    scale_colour_manual(values = c("#ed1C24","#00a651","#2276b9")) +
    xlab("Time") +
    ylab("Ratio (Chromosome/Main Mass)") +
    scale_y_continuous(limits=(c(0,2))) +
    ggtitle(sample_name) +
    annotate("text", x = min(df$Time_offset) + 3, y = 0.2, label = region) +
    theme(legend.position='none')
  
  file_path <- paste0("Output/Plots/", region, "_", type, "_", sample_name, ".png")
  ggsave(file_path, plot = p, dpi = 300)
}


##############################################################################################################

# Script execution is here

live_function_general("Data/LBR/LBR_lagging", "LBR_lag")
live_function_general("Data/LBR/LBR_polar", "LBR_polar")
live_function_general("Data/BAF/BAF_lagging", "BAF_lag")
live_function_general("Data/BAF/BAF_polar", "BAF_polar")
