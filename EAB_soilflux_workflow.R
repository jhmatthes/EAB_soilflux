# Workflow for processing raw data and reproducing all the 
# analysis within the manuscript: 
# Matthes, J.H., A.K. Lang, F.V. Jevon, S.J. Russell. Tree stress and mortality
# 	from emerald ash borer does not systematically alter short-term soil carbon flux
# 	in a mixed northeastern U.S. forest. Submitted to Forests, 21 Dec 2017.
# Code submitted to Github by J.H.M on 21 Dec 2017.

library(chron)
library(tidyverse)
library(stringr) 
library(lubridate)
library(forcats)
library(rstan)
library(ggridges)
library(gridExtra)

# Load local functions
file_sources <- list.files("R/functions/", pattern="*.R", full.names = TRUE)
sapply(file_sources, source, .GlobalEnv)

# File names for raw data from the LGR/Picarro analyzers
raw_files <- list.files("input/raw_CO2CH4_data/", full.names=TRUE)

# Load file with sampling dates and start times for flux measurements
date_time <- read_csv("input/times_key.csv")
rep_data  <- tibble(rep = date_time$rep, volume = date_time$rep_volume_L)
date_time <- date_time %>% select(contains("/"))

# Flux processing settings
init <- list()
init$lgr_ts <- 5 # LGR analyzer measures every 5 seconds
init$pic_ts <- 1 # Picarro analyzer measures every 1 second
init$lgr_fluxstart <- 5 # 10 for LGR
init$pic_fluxstart <- 25 # 50 for Picarro
init$lgr_fluxend   <- 2*60/init$lgr_ts # 24 timepoints
init$pic_fluxend   <- 100 # 50 timepoints
init$surfarea  <- pi*(4*2.54/2)^2 / 100^2 #m^2 
init$totalreps <- nrow(date_time) # total chambers in study
init$plotslope <- 0
init$outputfile <- 1

# Calculate soil CO2 & CH4 flux for each measurement date & replicate
flux_clean <- calculate_chamber_flux(raw_files, date_time, rep_data, init)          

# Merge temperature & moisture dataset with clean flux data
temp_moist  <- read.csv("input/TF_insttempmoisture.csv", 
                        header=TRUE, stringsAsFactors = FALSE) 

temp_moist$date <- as.Date(temp_moist$date, "%m/%d/%y") %>% 
  format(.,"%Y-%m") %>% 
  as.character()

flux_env <- merge(flux_clean, temp_moist, by.x = c("fmonth", "fid"),
                     by.y = c("date","collar"), all=TRUE) 

tree_health <- read.csv("input/TF_TreeHealth.csv", header=TRUE) 

flux_health <- merge(flux_env, tree_health, by.x = c("fid", "fyear"),
                     by.y = c("collar", "year"), all = TRUE) %>%
  arrange(tree.x, fid, fdate) %>%
  filter(!is.na(fdate)) %>%
  mutate(date_form = as.Date(fdate)) 

flux_all <- read_csv("input/TF_SoilCovariates.csv") %>%
  right_join(flux_health, by = "fid") %>%
  select(fid, fyear, fmonth, fdate, fmonthny, CO2_flux, CH4_flux, 
         tree = tree.x, 
         Soil.Temp, Soil.Moisture, status, dbh, root_biomass_g,
         pH, sand_frac, condition)

################################################################
# FIGURE 3: Climate in 2016 & 2017 compared to 30-year climatology
################################################################

# 30-year climatology data, 1981-2010: Concord airport
climatology_dat <- read_csv("input/climate_NOAA/ConcordAirport_19812010.csv") %>%
  mutate(year = year(DATE),
         fmonth = format(DATE, "%Y-%m"),
         month = month(DATE),
         TMEAN_C = ((TMAX + TMIN)/2 - 32)*(5/9),
         PRCP_cm = PRCP*2.54) %>%
  filter(year > 1980 & year <= 2010) %>%
  select(DATE, year, fmonth, month, TMEAN_C, PRCP_cm)

# 30-year mean annual temperature and mean annual precipitation
climatology_ann <- climatology_dat %>%
  group_by(year) %>%
  summarize(t_mean = mean(TMEAN_C, na.rm=TRUE), 
            t_miss = sum(is.na(TMEAN_C)),
            p_sum = sum(PRCP_cm, na.rm=TRUE), 
            p_miss = sum(is.na(PRCP_cm))) %>%
  filter(t_miss < 1, p_miss < 1) %>%
  summarize(t30 = mean(t_mean), tsd = sd(t_mean),
            p30 = mean(p_sum)*10, psd = sd(p_sum)*10)

# 30-year monthly mean temp & total rainfall: May - November
mon_climdat <- climatology_dat %>% 
  group_by(fmonth, month, year) %>%
  summarize(monthly_rain = sum(PRCP_cm, na.rm=TRUE),
            monthly_temp = mean(TMEAN_C, na.rm=TRUE),
            rain_miss = sum(is.na(PRCP_cm)),
            temp_miss = sum(is.na(TMEAN_C))) %>%
  filter(rain_miss < 1, temp_miss < 1, month >= 5 & month <= 11) 

# Weather data that overlaps experiment period
exp_dat <- read_csv("input/climate_NOAA/ConcordAirport_1517.csv") %>%
  mutate(year = year(DATE),
         fmonth = format(DATE, "%Y-%m"),
         month = month(DATE),
         TMEAN_C = ((TMAX + TMIN)/2 - 32)*(5/9),
         PRCP_cm = PRCP*2.54) %>%
  filter(month >= 5 & month <= 11, year > 2015) %>%
  select(DATE, year, month, fmonth, TMEAN_C, PRCP_cm)

month_dat <- exp_dat %>%
  group_by(fmonth, month, year) %>%
  summarize(precip = sum(PRCP_cm, na.rm=TRUE),
            monthly_temp = mean(TMEAN_C, na.rm=TRUE),
            monthly_temp_sd = sd(TMEAN_C, na.rm=TRUE))

# Join monthly precip to flux data
flux_all_2 <- left_join(flux_all, 
                        month_dat, by = "fmonth")

# Plot 2016 & 2017 monthly temperature & monthly precip on 30-year climatology
fig3a <- ggplot(mon_climdat, aes(x = factor(month), y = monthly_temp)) + 
  geom_boxplot() +
  labs(x = "", y = "Monthly Mean Temperature [C]") + 
  scale_x_discrete(labels = c("May", "June", "July", "Aug", "Sep", "Oct", "Nov")) +
  geom_point(data = month_dat, 
             aes(x = as.factor(month), y = monthly_temp, 
                 color = as.factor(year)), 
             shape = 17, size = 3) +
  scale_color_brewer(palette = "Set1") +
  labs(color = "Year") +
  theme_minimal() +
  annotate("text", x = 0.5, y = 30, label = "a)") +
  guides(color = FALSE) 

fig3b <- ggplot(mon_climdat, aes(x = factor(month), y = monthly_rain)) + 
  geom_boxplot() +
  labs(x = "", y = "Monthly Precipitation [cm]") + 
  scale_x_discrete(labels = c("May", "June", "July", "Aug", "Sep", "Oct", "Nov")) +
  geom_point(data = exp_dat, aes(x = as.factor(month), 
                                 y = precip, color = as.factor(year)), 
             shape = 17, size = 3) +
  scale_color_brewer(palette = "Set1") +
  labs(color = "Year")  +
  annotate("text", x = 0.5, y = 35, label = "b)") +
  theme_minimal()+
  theme(legend.position="bottom")

# Write Figure 3: Monthly Precip & Temp 2016/2017 vs 30-year climatology
pdf("output/figures/Fig3.pdf")
grid.arrange(fig3a, fig3b, ncol = 1)
dev.off()

################################################################
# FIGURE S1: Soil Moisture versus monthly total precipitation
################################################################
pdf("output/figures/FigS1.pdf")
par(mar=c(5.1,5.1,4.1,2.1))
plot(flux_all_2$precip, flux_all_2$Soil.Moisture,
     xlab = "Monthly Precipitation [cm]",
     ylab = expression(paste("Instantaneous Soil Moisture [", m^{3}, m^{-3},"]")))
dev.off() 

################################################################
# FIGURE 4: Differences in EAB status, 2016 vs 2017 
################################################################
# Re-order EAB status factors so healthy is 1st, dead is last
tree_health$status <- factor(tree_health$status, 
                             levels = c("healthy", "impacted", "dead"))

# Summarize status numbers by year
tree_health_summ <- tree_health %>%
  group_by(status, year) %>%
  summarize(count = n()) 
tree_health_summ$year <- as.factor(tree_health_summ$year)
dummy_healthy2017 <- data.frame(status = "healthy", 
                                year = "2017", count = 0)
tree_health_summ2 <- bind_rows(tree_health_summ,
                               dummy_healthy2017) 
tree_health_summ2$status <- factor(tree_health_summ2$status, 
                                   levels = c("healthy", "impacted", "dead"))

# Plot numbers of trees in EAB status categories by year
fig4 <- ggplot(tree_health_summ2) + 
  geom_bar(aes(x = status, y = count, fill=as.factor(year)), 
           stat = "identity", position = "dodge", alpha = 0.75) +
  scale_fill_brewer(palette = "Set1") +
  scale_x_discrete(labels = c("Healthy", "Impacted", "Dead")) + 
  labs(x = "EAB Status", y = "Tree Count", fill = "Year") +
  theme_minimal() +
  geom_text(aes(x = status, y = count + 0.25, label = count, group = year),
            position = position_dodge(0.9))

pdf("output/figures/Fig4.pdf", width = 5, height = 3)
fig4
dev.off()

################################################################
# FIGURES 5 & 6: CO2 flux & CH4 uptake data by year and EAB class
################################################################
# Fig 5: Monthly flux patterns in 2016 vs 2017
pdf("output/figures/Fig5.pdf")
flux_timeseries_plot(flux_all)
dev.off()

# Fig 6: CO2/CH4 flux time series by EAB class
pdf("output/figures/Fig6.pdf")
flux_timestatus_plot(flux_all)
dev.off()

################################################################
# Model Set 1: Does EAB status impact soil microclimate?
# Two mixed effects models: Soil Moisture, Soil Temperature
################################################################
stan_iter <- 2000 # number of MCMC iterations within each model

# Model 1a: Soil Moisture Mixed Effects model 
# SM_i = rain_i + sand_i + collar_i + status_i
SM_full_data_Stan <- flux_all_2 %>%
  filter(!is.na(Soil.Moisture))

SM_full_data <- list(collar = as.integer(factor(flux_all_2$fid)),
                     status = as.integer(factor(flux_all_2$status)),
                     month = as.integer(factor(flux_all_2$fmonthny)),
                     sm = flux_all_2$Soil.Moisture,
                     rain = flux_all_2$precip,
                     sand = flux_all_2$sand_frac,
                     N = nrow(flux_all_2),
                     C = length(unique(flux_all_2$fid)),
                     M = length(unique(flux_all_2$fmonthny)),
                     S = length(unique(flux_all_2$status)))

SM_fullmodel <- stan(file = "R/stan/AshSM_CollarStatus.stan", 
                     data = SM_full_data, 
                     iter = stan_iter, chains = 4,
                     control = list(adapt_delta = 0.9))

# Write summary output
SM_table <- stan_outstats(SM_fullmodel)
write_csv(SM_table, "output/data/SM_modout_raw.csv")

# Model 1b: Soil Temperature Mixed Effects model 
# Temp_i = month_i + collar_i + status_i
Temp_full_data_Stan <- flux_all_2 %>%
  filter(!is.na(Soil.Temp))

Temp_full_data <- list(collar = as.integer(factor(flux_all_2$fid)),
                     status = as.integer(factor(flux_all_2$status)),
                     month = as.integer(factor(flux_all_2$fmonthny)),
                     temp = flux_all_2$Soil.Temp,
                     N = nrow(flux_all_2),
                     M = length(unique(flux_all_2$fmonthny)),
                     C = length(unique(flux_all_2$fid)),
                     S = length(unique(flux_all_2$status)))

Temp_fullmodel <- stan(file = "R/stan/AshTemp_CollarStatus.stan", 
                     data = Temp_full_data, 
                     iter = stan_iter, chains = 4,
                     control = list(adapt_delta = 0.9))

Temp_table <- stan_outstats(Temp_fullmodel)
write_csv(Temp_table, "output/data/Temp_modout_raw.csv")

################################################################
# Model Set 2: Does EAB status impact soil CO2 flux/CH4 uptake?
################################################################

# Model 2a: CO2 flux
# Random Effects: collar (replicate), month
# Fixed Effects: soil moisture, temp, month rainfall, DBH, root biomass, pH, sand
CO2_full_data <- list(collar = as.integer(factor(flux_all_2$fid)),
                      month = as.integer(factor(flux_all_2$fmonthny)),
                      co = flux_all_2$CO2_flux,
                      sm = flux_all_2$Soil.Moisture,
                      temp = flux_all_2$Soil.Temp,
                      rain = flux_all_2$precip,
                      dbh = flux_all_2$dbh,
                      root = flux_all_2$root_biomass_g,
                      sand = flux_all_2$sand_frac,
                      pH = flux_all_2$pH,
                      condition = as.integer(factor(flux_all_2$condition)),
                      N = nrow(flux_all_2),
                      C = length(unique(flux_all_2$fid)),
                      M = length(unique(flux_all_2$fmonthny)),
                      S = length(unique(flux_all_2$condition)))

CO2_fullmodel <- stan(file = "R/stan/AshCO2_CollarMonthALLFixedEff.stan", 
                      data = CO2_full_data, 
                      iter = stan_iter, chains = 4)

# Write summary output
CO2_table <- stan_outstats(CO2_fullmodel)
write_csv(CO2_table, "output/data/CO2_modout_raw.csv")

# Model 2b: CH4 uptake
# Random Effects: collar (replicate), month
# Fixed EFfects: soil moisture, temperature, monthly rainfall, pH, sand
CH4_full_data_Stan <- flux_all_2 %>%
  filter(!is.na(CH4_flux))

CH4_full_data <- list(collar = as.integer(factor(CH4_full_data_Stan$fid)),
                      month = as.integer(factor(CH4_full_data_Stan$fmonthny)),
                      ch = abs(CH4_full_data_Stan$CH4_flux)*1e3,
                      sm = CH4_full_data_Stan$Soil.Moisture,
                      temp = CH4_full_data_Stan$Soil.Temp,
                      rain = CH4_full_data_Stan$precip,
                      root = CH4_full_data_Stan$root_biomass_g,
                      sand = CH4_full_data_Stan$sand_frac,
                      pH = CH4_full_data_Stan$pH,
                      condition = as.integer(factor(CH4_full_data_Stan$condition)),
                      N = nrow(CH4_full_data_Stan),
                      C = length(unique(CH4_full_data_Stan$fid)),
                      M = length(unique(CH4_full_data_Stan$fmonthny)))

CH4_fullmodel <- stan(file = "R/stan/AshCH4_CollarMonthALLFixedEff.stan", 
                      data = CH4_full_data, 
                      iter = stan_iter, chains = 4)

# Write summary output
CH4_table <- stan_outstats(CH4_fullmodel)
write_csv(CH4_table, "output/data/CH4_modout_raw.csv")

################################################################
# FIGURE 7: Plot EAB transition class on random collar effects
################################################################
condition_colors <- c("H-I" = "darkgreen", "I-I" = "palegreen3",
                      "I-D" = "skyblue3", "D-D" = "navy")

# Tidy the posterior samples of collar random effects 
collarsamps_CO2 <- format_stan_collar_raneff(flux_all_2, 
                                             CO2_fullmodel, unlog = FALSE)

collarsamps_CH4 <- format_stan_collar_raneff(CH4_full_data_Stan, 
                                             CH4_fullmodel, unlog = FALSE)

# Re-order collar random effects by size and condition effects by level
collar_order_CO2 <- collarsamps_CO2 %>% 
  group_by(collar) %>%
  summarize(median = median(samples, na.rm=TRUE)) %>%
  arrange(median)

collarsamps_CO2$collar <- factor(collarsamps_CO2$collar, 
                                      levels = collar_order_CO2$collar)
collarsamps_CO2$condition <- factor(collarsamps_CO2$condition, 
                                         levels = c("H-I", "I-I","I-D","D-D"))

collar_order_CH4 <- collarsamps_CH4 %>% 
  group_by(collar) %>%
  summarize(median = median(samples, na.rm=TRUE)) %>%
  arrange(median)

collarsamps_CH4$collar <- factor(collarsamps_CH4$collar, 
                                 levels = collar_order_CH4$collar)
collarsamps_CH4$condition <- factor(collarsamps_CH4$condition, 
                                    levels = c("H-I", "I-I","I-D","D-D"))

# Plot collar effects by size and EAB transition class
co2_label <- bquote("Replicate effect ["*mu~"mol"~CO[2]~ m^-2~s^-1*"]")
ch4_label <- bquote("Replicate effect [nmol"~CH[4]~ m^-2~s^-1*"]")

fig7a <- ggplot(collarsamps_CO2, aes(x = samples, y = collar, fill = condition)) +
  geom_density_ridges(alpha = 0.8) +
  geom_vline(xintercept = 0, alpha = 0.5) +
  scale_fill_manual(values = condition_colors) +
  theme_ridges(grid = FALSE, font_size = 11) + 
  labs(x = co2_label, y = "Replicate ID") +
  theme(legend.position="bottom") + 
  annotate("text", x = -1, y = 37, label = "a)") +
  labs(fill = "")

fig7b <- ggplot(collarsamps_CH4, 
                aes(x = samples, y = collar, fill = condition)) +
  geom_density_ridges(alpha = 0.8) +
  geom_vline(xintercept = 0, alpha = 0.5) +
  scale_fill_manual(values = condition_colors) +
  theme_ridges(grid = FALSE, font_size = 11) + 
  labs(x = ch4_label, y = "Replicate ID") +
  theme(legend.position="bottom") + 
  annotate("text", x = -3, y = 37, label = "b)") +
  labs(fill = "")

fig7_legend <- g_legend(fig7a)
pdf("output/figures/Fig7.pdf")
grid.arrange(arrangeGrob(fig7a + theme(legend.position="none"), 
                         fig7b + theme(legend.position="none"), ncol = 2),
             fig7_legend, nrow=2, heights=c(10, 1))
dev.off()

################################################################
# FIGURE S2: Plot random month effects
################################################################

# Tidy the posterior samples of month random effects 
monthsamps_CO2 <- format_stan_month_raneff(data = flux_all_2, 
                                           CO2_fullmodel, unlog = FALSE)
month_order_CO2 <- monthsamps_CO2 %>% 
  group_by(month) %>%
  summarize(month_med = median(prob)) %>%
  arrange(desc(month_med))

monthsamps_CO2$month <- factor(monthsamps_CO2$month, levels = month_order_CO2$month) %>%
  plyr::revalue(., c("05"="May", "06"="June", "07"="July","08"="Aug", "09"="Sep", "10"="Oct", "11"="Nov"))

monthsamps_CH4 <- format_stan_month_raneff(data = CH4_full_data_Stan, 
                                           CH4_fullmodel, unlog = FALSE)
month_order_CH4 <- monthsamps_CH4 %>% 
  group_by(month) %>%
  summarize(month_med = median(prob)) %>%
  arrange(desc(month_med))

monthsamps_CH4$month <- factor(monthsamps_CH4$month, levels = month_order_CH4$month) %>%
  plyr::revalue(., c("05"="May", "06"="June", "07"="July","08"="Aug", "09"="Sep", "10"="Oct", "11"="Nov"))

figS2a <- ggplot(monthsamps_CO2, aes(x = prob, y = month)) +
  geom_density_ridges(alpha = 0.8) +
  geom_vline(xintercept = 0, alpha = 0.5) +
  theme_ridges(grid = FALSE, font_size = 11) + 
  labs(x = co2_label, y = "Month") +
  xlim(-0.75, 0.75) +
  annotate("text", x = -0.75, y = 8, label = "a)")

figS2b <- ggplot(monthsamps_CH4, aes(x = prob, y = month)) +
  geom_density_ridges(alpha = 0.8) +
  geom_vline(xintercept = 0, alpha = 0.5) +
  theme_ridges(grid = FALSE, font_size = 11) + 
  labs(x = ch4_label, y = "Month") +
  xlim(-2.5, 2.5) +
  annotate("text", x = -2.5, y = 8, label = "b)")

pdf("output/figures/FigS2.pdf")
grid.arrange(figS2a, figS2b, ncol = 2)
dev.off()

