# Format MCMC output for Month random effects

format_stan_month_raneff <- function(data, fit, unlog){
  AshDat <- data
  
  ash_dates <- unique(AshDat$fmonthny)
  
  
  for(i in 1:length(unique(AshDat$fmonthny))){
    tmp <- rstan::extract(fit, pars = c(paste0("month_int[",i,"]")))$month_int
    
    if(unlog == TRUE){
      tmp <- exp(tmp)
    } else if(unlog == FALSE){
      tmp <- tmp
    }
    
    name <- paste0("month",i)
    assign(name, data.frame(month = rep(ash_dates[i]), prob = tmp))
  }
  
  date_samples <- rbind(month1, month2, month3, month4, month5, month6, month7)
  
  return(date_samples)
}