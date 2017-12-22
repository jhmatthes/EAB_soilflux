# Plot monthly flux time series difference by tree EAB status class: 
# Impacted = transitioned EAB classes during experiment, versus
# Dead = dead from the start of the experiment

flux_timestatus_plot <- function(data){
  
  # Aggregate condition into Impacted (transitioned during experiment)
  # or Dead (dead from the beginning of the experiment)
  data2 <- mutate(data, condition2 = fct_collapse(condition, 
                                                  "Impacted" = c("H-I", "I-I", "I-D"),
                                                  "Dead" = "D-D"))
  
  data2$condition2 <- factor(data2$condition2, levels = c("Impacted", "Dead"))
  cols <- c("Impacted" = "palegreen3", "Dead" = "navy")
  
  # Find data with replicates in each month-year for difference
  reps_bothyears <- data2 %>%
    filter(fmonthny != "05", fmonthny != "11") %>%
    select(fid, fyear, fmonthny, condition2, CO2_flux, CH4_flux) %>%
    mutate(CH4_flux = abs(CH4_flux*1e3)) %>%
    group_by(fid, fmonthny, condition2) %>%
    summarize(count = n()) %>%
    filter(count > 1)
  
  data_yrreps <- right_join(data2, reps_bothyears, by = c("fid", "fmonthny", "condition2")) %>%
    mutate(CH4_flux = abs(CH4_flux)*1e3) %>%
    group_by(fid, fmonthny, condition2) %>%
    summarize(CO2 = diff(CO2_flux, na.rm=TRUE), 
              CH4 = diff(CH4_flux, na.rm=TRUE)) 
  
  fig3a <- ggplot(data_yrreps) + 
    geom_boxplot(aes(x = fmonthny, y = CO2, col=condition2)) +
    labs(x = "", y = expression(paste(Delta~CO[2]," Flux [", mu, "mol ", m^{-2}, s^{-1}, "]")),
         color = "EAB Transition Type") + 
    scale_color_manual(values = cols) +
    scale_x_discrete(labels = c("June", "July", "Aug", "Sep", "Oct")) +
    guides(color = FALSE) +
    theme_minimal() + 
    geom_hline(yintercept = 0, lty = 2, col = "grey", size = 1) + 
    annotate("text", x = 0.5, y = 5, label = "a)")
  
  fig3b <- data_yrreps %>% filter(CH4 < 5) %>%
    ggplot() + 
    geom_boxplot(aes(x = fmonthny, y = CH4, col=condition2)) +
    labs(x = "", y = expression(paste(Delta~CH[4]," Uptake [nmol ", m^{-2}, s^{-1}, "]")),
         color = "EAB Transition Type") + 
    scale_color_manual(values = cols, name = "EAB condition") +
    scale_x_discrete(labels = c("June", "July", "Aug", "Sep", "Oct")) +
    theme_minimal() + 
    geom_hline(yintercept = 0, lty = 2, col = "grey", size = 1) +
    theme(legend.position="bottom") + 
    annotate("text", x = 0.5, y = 2, label = "b)")

  grid.arrange(fig3a, 
               fig3b, 
               ncol = 1)
  
}
