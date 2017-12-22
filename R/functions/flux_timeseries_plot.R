# Plot time series of CO2/CH4 flux measurements

flux_timeseries_plot <- function(data){
  data <- mutate(data, CH4_flux = abs(CH4_flux))
  
  # Perform unpaired t-tests for 2016 vs 2017 differences
  months <- c("06", "07", "08", "09", "10")
  co2_p <- ch4_p <- vector()
  for(m in 1:length(months)){
    flux_2016 <- filter(data, fmonthny == months[m], fyear == "2016") 
    flux_2017 <- filter(data, fmonthny == months[m], fyear == "2017") 
    co2_p[m] <- t.test(flux_2016$CO2_flux, flux_2017$CO2_flux)$p.value
    ch4_p[m] <- t.test(flux_2016$CH4_flux, flux_2017$CH4_flux)$p.value
  }
  
  co2_months_sig <- c("May", paste0("June",ifelse(co2_p[1]<0.05,"*","")),
                      paste0("July",ifelse(co2_p[2]<0.05,"*","")),
                      paste0("Aug",ifelse(co2_p[3]<0.05,"*","")),
                      paste0("Sep",ifelse(co2_p[4]<0.05,"*","")),
                      "Oct","Nov")
  
  ch4_months_sig <- c("May", paste0("June",ifelse(ch4_p[1]<0.05,"*","")),
                      paste0("July",ifelse(ch4_p[2]<0.05,"*","")),
                      paste0("Aug",ifelse(ch4_p[3]<0.05,"*","")),
                      paste0("Sep",ifelse(ch4_p[4]<0.05,"*","")),
                      paste0("Oct",ifelse(ch4_p[5]<0.05,"*","")),
                      "Nov")
  
  # Make simpler dataframes for CO2/CH4 flux with dummy high values for 
  # ggplot boxplots to format as same width (despite missing factors)
  flux_CO2 <- data %>%
    select(fmonthny, fyear, CO2_flux) 
  flux_CO2_2 <- rbind(flux_CO2, 
                      data.frame(fmonthny = c("05", "11"), fyear = c("2016","2017"), 
                                 CO2_flux = rep(100,2)))
  
  flux_CH4 <- data %>%
    select(fmonthny, fyear, CH4_flux)
  flux_CH4_2 <- rbind(flux_CH4, 
                      data.frame(fmonthny = c("05", "11"), fyear = c("2016","2017"), 
                                 CH4_flux = rep(-1,2)))
  
  # FIG 2: Boxplots of monthly CO2 flux in 2016 vs 2017
  fig2a <- ggplot(flux_CO2_2) + 
    geom_boxplot(aes(x = fmonthny, y = CO2_flux, col=as.factor(fyear))) +
    scale_colour_brewer(palette = "Set1") +
    coord_cartesian(ylim = range(flux_CO2$CO2_flux) + c(-.25, .25)) + 
    labs(x = "", 
         y = expression(paste(CO[2]," Flux [", mu, "mol ", m^{-2}, s^{-1}, "]")),
         color = "Year") +
    scale_x_discrete(labels = co2_months_sig) +
    theme_minimal() +
    annotate("text",x = 0.5, y = 10, label = "a)") +
    guides(color=FALSE)
  
  
  # FIG 1B: Boxplots of monthly CH4 flux in 2016 vs 2017
  fig2b <- ggplot(flux_CH4_2) + 
    geom_boxplot(aes(x = fmonthny, y = CH4_flux*1e3, col=as.factor(fyear))) +
    scale_colour_brewer(palette = "Set1") +
    coord_cartesian(ylim = c(0, 9) + c(-.25, .25)) + 
    labs(x = "", 
         y = expression(paste(CH[4]," Uptake [nmol ", m^{-2}, s^{-1}, "]")),
         color = "Year") +
    scale_x_discrete(labels = ch4_months_sig) +
    annotate("text", x = 0.5, y = 8, label = "b)") +
    theme_minimal() +
    theme(legend.position="bottom")
  
  grid.arrange(fig2a, fig2b, ncol = 1)
  
}