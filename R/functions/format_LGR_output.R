# This function reads a raw text output file from a Los Gatos Research (LGR) Greenhouse Gas Analyzer (GGA)
# and exports the measurement times and CO2 and CH4 concentration data. 
#
# This function requires the chron library to be loaded for the times() function.
#
# The format_LGR_output() function requires two input variables:
#   1. LGR_files = a vector of LGR file names to be processed
#   2. LGR_match_files = index of LGR file names in LGR_files that match single day
# Output of this function is a data frame with measurement times and CO2 & CH4 concentration data

format_LGR_output <- function(LGR_files, match_LGR_files){
  
  # Aggregate the LGR files if there is more than one file per day.
  if(length(match_LGR_files)>1){
    for(day.file in 1:length(match_LGR_files)){ 
      
      # Read LGR data and clean out extra columns and rows
      LGR_tmp  <- read.csv(LGR_files[match_LGR_files[day.file]],header=TRUE,skip=1,stringsAsFactors = FALSE)
      LGR_tmp  <- LGR_tmp[,1:14] #only use first 14 columns
      
      if(sum(LGR_tmp[,1]=="-----BEGIN PGP MESSAGE-----")>0){
        LGR_tmp  <- LGR_tmp[1:(which(LGR_tmp[,1]=="-----BEGIN PGP MESSAGE-----")-1),]
      }
      
      # Format messy column 1 into LGR date & time.
      LGR_date_tmp   <- strsplit(sapply(LGR_tmp$Time, as.character)," ")
      LGR_date <- LGR_time <- vector()
      for(i in 1:length(LGR_date_tmp)){
        LGR_date[i] <- LGR_date_tmp[[i]][3]
        LGR_time[i] <- LGR_date_tmp[[i]][4]
      }
      LGR_tmp_times <- times(LGR_time) 
      
      # Aggregate each file into a big daily file.
      if(day.file == 1){
        LGR_data  = LGR_tmp
        LGR.times = LGR_tmp_times
      } else {
        LGR_data  = rbind(LGR_data,LGR_tmp)
        LGR.times = c(LGR.times,LGR_tmp_times)
      }
    }
  } else { # If there is only one file for this day:
    # Read LGR data and clean out unnecessary columns and text rows.
    LGR_data  <- read.csv(LGR_files[match_LGR_files],
                          header=TRUE,skip=1,stringsAsFactors = FALSE)
    LGR_data  <- LGR_data[,1:14] #only use first 14 columns
    if(sum(LGR_data[,1]=="-----BEGIN PGP MESSAGE-----")>0){
      LGR_data  <- LGR_data[1:(which(LGR_data[,1]=="-----BEGIN PGP MESSAGE-----")-1),]
    }
    
    # Format column 1 into LGR date & time.
    LGR_date_tmp   <- strsplit(sapply(LGR_data$Time, as.character)," ")
    LGR_date <- LGR_time <- vector()
    for(i in 1:length(LGR_date_tmp)){
      LGR_date[i] <- LGR_date_tmp[[i]][3]
      LGR_time[i] <- LGR_date_tmp[[i]][4]
    }
    LGR.times <- times(LGR_time) 
  }
  
  # Format output for times, CO2, and CH4
  LGR_dat <- data.frame(times = LGR.times,
                             CO2 = LGR_data$X.CO2.d_ppm, 
                        CH4 = LGR_data$X.CH4.d_ppm)
  
  return(LGR_dat)
}