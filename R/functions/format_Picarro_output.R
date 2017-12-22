# This function reads a raw text output file from a Picarro GasScouter Analyzer and exports the 
# time, CO2, and CH4 concentration data. 
#
# This function requires the chron library to be loaded for the times() function.
#
# The format_Picarro_output() function requires two inputs:
#   1. raw_files = a vector of raw file names to be processed
#   2. match_files = index of Picarro raw file names that match a single day
# Output of this function is a data frame with measurement times and CO2 & CH4 concentration data

format_Picarro_output <- function(raw_files, match_files){
  
  # Aggregate the files if there is more than one file per day.
  if(length(match_files) > 1){
    for(f in 1:length(match_files)){ 
      
      # Read in raw text Picarro data for each file
      dat  <- read.table(raw_files[match_files][f], header=TRUE)
      
      # Format time column
      Pic_tmp_times <- times(dat$TIME) 
      
      # Aggregate each file into a big daily file
      if(f == 1){
        Pic_data  <- dat
        Pic_times <- Pic_tmp_times
      } else {
        Pic_data  <- rbind(Pic_data,dat)
        Pic_times <- c(Pic_times,Pic_tmp_times)
      }
    }
  } else {
    Pic_data  <- read.table(raw_files[match_files], header=TRUE)
    Pic_times <- times(Pic_data$TIME)
  }
  
  # Format output for times, CO2, and CH4
  Pic_df <- data.frame(times = Pic_times,
                        CO2 = Pic_data$CO2_dry, CH4 = Pic_data$CH4_dry)
  return(Pic_df)
}

