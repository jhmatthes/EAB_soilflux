# Function to summarize Stan output into a data frame with median 
# and 95% interval of the posterior distribution

stan_outstats <- function(stan_output){
  
  out_df <- as.data.frame(stan_output)
  med <- vector()
  q95 <- matrix(nrow = ncol(out_df), ncol = 2)
  for(c in 1:ncol(out_df)){
    med[c] <- median(out_df[,c])
    q95[c,] <- quantile(out_df[,c], probs = c(0.025, 0.975))
  }
  
  out_table <- data.frame(var = colnames(out_df), med, q95)
  colnames(out_table)[3:4] <- c("q025", "q975")
  
  return(out_table)
  
}