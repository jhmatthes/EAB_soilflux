# Extract random effects of soil respiration replicate collars

format_stan_collar_raneff <- function(data, fit, unlog){
  AshDat <- data
  ash_trees <- unique(AshDat[c("fid", "tree")])
  ash_status <- AshDat %>%
    filter(fyear == "2017") %>%
    group_by(fid, tree) %>%
    summarize(status = unique(status)[1],
              condition = unique(condition)[1])
  n_chambers <- length(unique(AshDat$fid))
  
  for(i in 1:n_chambers){
    tmp <- rstan::extract(fit, pars = c(paste0("collar_int[",i,"]")))$collar_int
    if(unlog == TRUE){
      tmp <- exp(tmp)
    } else if(unlog == FALSE){
      tmp <- tmp
    }
    name <- paste0(i)

    if(i == 1){
      collar_samples <- data.frame(collar = name, 
                                   tree = rep(ash_trees$tree[i]),
                                   status = rep(ash_status[ash_status$fid == ash_trees$fid[i],]$status),
                                   condition = rep(ash_status[ash_status$fid == ash_trees$fid[i],]$condition),
                                   samples = tmp)
    } else {
      tmp2 <- data.frame(collar = name,
                         tree = rep(ash_trees$tree[i]),
                         status = rep(ash_status[ash_status$fid == ash_trees$fid[i],]$status),
                         condition = rep(ash_status[ash_status$fid == ash_trees$fid[i],]$condition),
                         samples = tmp)
      collar_samples <- bind_rows(collar_samples, tmp2)
    }
  }
  return(collar_samples)
}
