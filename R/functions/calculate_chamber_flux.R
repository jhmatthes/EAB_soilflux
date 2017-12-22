# Process the raw output from the Los Gatos Research Greenhouse Gas Analyzer
# or the Picarro GasScouter analyzer into CO2 and CH4 flux
# according to the dates and start times provided.
#
# Jackie Matthes, moved from workflow.R to separate function 6 Nov 2017
#
# INPUT:
#  * raw_files: vector of raw LGR/Picarro file names 
#  * date_time: data frame with chamber start times, col names are sampling dates
#  * rep_data: data frame with chamber replicate names and volumes (L)
#  * init: list of flux processing settings (in workflow.R for transparency/access)
# 
# OUTPUT: 
#  * flux_clean: data frame with chamber ID, date info, CO2 & CH4 fluxes and fit
#
# Requires libs: tidyverse, stringr, lubridate, chron

calculate_chamber_flux <- function(raw_files, date_time, rep_data, init){
  
  # Format dates & times
  dates <- date_time %>% 
    colnames() %>% 
    as.Date(., "%m/%d/%y") %>% 
    as.character()
  colnames(date_time) <- dates
  dates_my <- format(strptime(dates, "%Y-%m-%d"), "%Y-%m")
  dates_y  <- format(strptime(dates, "%Y-%m-%d"), "%Y")
  
  # System volume to total mol in the system (for ppm to mol conversion)
  vol_system   <- 0.315 + 0.001 + rep_data$volume
  init$ntotal <- (739*vol_system)/(62.363577*298.15) # n = (P * vol_system) / (R * T)
  
  ### Process the GHG analyzer files into CO2/CH4 flux by matching the chamber     
  ### measurements dates & times to the raw file information 
  
  # Create flux output storage
  empty_vector <- rep(NA, length = init$totalreps*length(dates))
  flux_data <- data.frame(fdate = empty_vector, fmonth = empty_vector, 
                          fid = empty_vector, fyear = empty_vector,
                          CO2_flux = empty_vector, CO2_r2 = empty_vector, 
                          CO2_SE = empty_vector,
                          CH4_flux = empty_vector, CH4_r2 = empty_vector,
                          CH4_SE = empty_vector)
  
  # Flux calculation: Loop over each sampling date in the start times table
  for(d in 1:length(dates)){
    
    # Find raw data files that match each date 
    match_files <- grep(dates[d], raw_files)
    
    if(length(match_files)>0){ # found an LGR file
      analyzer <- "LGR"
      print(paste0(dates[d], ": Working on an ", analyzer, " file."))
      
      # Read raw LGR analyzer CO2/CH4 concentration data
      conc_data <- format_LGR_output(raw_files, match_files)
      analyzer_timestep <- init$lgr_ts
      flux_start <- init$lgr_fluxstart
      flux_end   <- init$lgr_fluxend
      
    } else if(length(match_files)==0){ # found a Picarro file
      new_date <- str_replace_all(dates[d], "-", "")
      match_files <- grep(new_date, raw_files)
      analyzer <- "Picarro"
      print(paste0(dates[d], ": Working on a ", analyzer, " file."))
      
      # Read raw Picarro analyzer CO2/CH4 concentration data
      conc_data <- format_Picarro_output(raw_files, match_files)
      analyzer_timestep <- init$pic_ts
      flux_start <- init$pic_fluxstart
      flux_end   <- init$pic_fluxend
    } 
    if(length(match_files)==0) { # no file found
      print("Couldn't match a raw data file to sampling date!")
    }
    
    # Set up replicate & start time simplified table
    rep_times <- data.frame(rep = rep_data$rep, 
                            start = date_time[,d])
    
    # Loop over chambers measured on each date and calculate CO2/CH4 flux
    for(c in 1:nrow(rep_times)){
      
      # Index for storage: (day number-1) * total reps + nth rep
      rep_index <- (d-1)*init$totalreps + c 
      
      # Replicate & date bookkeeping
      flux_data$fdate[rep_index]  <- dates[d]
      flux_data$fmonth[rep_index] <- dates_my[d]
      flux_data$fyear[rep_index]  <- dates_y[d]
      flux_data$fid[rep_index]    <- rep_data$rep[c]
      
      # If chamber measurement missing
      if(is.na(rep_times[c,2])){
        flux_data$CO2_r2[rep_index]   <- flux_data$CH4_r2[rep_index] <- NA
        flux_data$CO2_flux[rep_index] <- flux_data$CH4_flux[rep_index] <- NA
        flux_data$CO2_SE[rep_index]   <- flux_data$CH4_SE[rep_index] <- NA
        
      } else {
        
        # Find replicate start time and match to nearest time in raw concentration file
        rep_time    <- times(rep_times[c,2])
        match_times <- which(abs(conc_data$times-rep_time)==min(abs(conc_data$times-rep_time)))
        
        # Define the flux time period, get CO2/CH4 concentration data for that period
        flux_period <- (match_times+flux_start):(match_times+flux_end) 
        CO2_conc    <- conc_data$CO2[flux_period]   #matching CO2 concentrations 
        CH4_conc    <- conc_data$CH4[flux_period]   #matching CH4 concentrations
        flux_time   <- conc_data$time[flux_period]  #matching times
        flux_seconds <- seq(1, ((flux_end - flux_start)+1)*analyzer_timestep,
                            by=analyzer_timestep) #seconds for linear fit
        
        # Fit linear model for x = seconds, y = CH4/CO2 concentration
        lm_CO2 <- lm(CO2_conc ~ flux_seconds)
        lm_CH4 <- lm(CH4_conc ~ flux_seconds)
        CO2_sl <- summary(lm_CO2)$coefficients[2]
        CH4_sl <- summary(lm_CH4)$coefficients[2]
        flux_data$CO2_r2[rep_index] <- summary(lm_CO2)$r.squared
        flux_data$CH4_r2[rep_index] <- summary(lm_CH4)$r.squared
        flux_data$CO2_SE[rep_index] <- (coef(summary(lm_CO2))[, "Std. Error"][2]*init$ntotal[c])/init$surfarea
        flux_data$CH4_SE[rep_index] <- (coef(summary(lm_CH4))[, "Std. Error"][2]*init$ntotal[c])/init$surfarea
        
        # Calculate CO2/CH4 chamber flux in umol/(m^2 * s)
        # V_CO2/V_T = n_CO2/n_total, n_CO2 = (ppm*10^-6) * n_total
        flux_data$CO2_flux[rep_index] <- (CO2_sl*init$ntotal[c])/init$surfarea 
        flux_data$CH4_flux[rep_index] <- (CH4_sl*init$ntotal[c])/init$surfarea 
        
        # Make plots of CO2_conc/CH4_conc vs time to visually inspect
        if(init$plotslope == 1){
          plot(flux_seconds, CO2_conc, 
               main=paste(" Date: ",dates[d]," Rep num: ",rep_data$rep[c],sep=""))
          plot(flux_seconds, CH4_conc, 
               main=paste(" Date: ",dates[d]," Rep num: ",rep_data$rep[c],sep=""))
        }
      } # end else.if start time exists
    } #end chamber loop
  } #end day loop
  
  # Cleaning: Remove flux values with R2 from linear fit < 0.9 (CO2), < 0.8 (CH4)
  flux_data[which(flux_data$CO2_r2 < 0.9), grep("CO2", names(flux_data))] <- NA
  flux_data[which(flux_data$CH4_r2 < 0.8), grep("CH4", names(flux_data))] <- NA
  
  # 1. Concatenate July 2016 & Oct 2016 measurement dates (spread across 2 dates)
  # 2. Add column with month, no year
  # 3. Remove chambers 16b & 21 (extreme heterotroph flux)
  flux_clean <- flux_data %>%
    filter(!is.na(CO2_flux) | !is.na(CH4_flux), CO2_flux < 10) %>%
    mutate(fmonthny = substring(fmonth, 6,7)) 
  
  # Write .csv file?
  if(init$outputfile == 1){
    write_csv(flux_clean, "output/AshFlux_5Nov17.csv")
  }
  
  # OUTPUT data frame back to workflow
  return(flux_clean)
}
