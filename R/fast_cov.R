library('tidyverse')
library('missMethods')

fast_cov <- function(X_miss , cor = FALSE){
  
  # Function for fast covariance estimation from missing data
  
  M <- is.na(X_miss)
  obs_mat <- t(1-M)%*%(1-M)
  
  n <- dim(X_miss)[1]
  
  X_imp_mean <- impute_mean(X_miss)
  
  mean_vec <- apply(X_miss, MARGIN = 2, function(x){return(mean(x , na.rm = T))})
  
  X_imp_demean <- t(t(X_imp_mean) - mean_vec)
  
  est_cov_mat <- t(X_imp_demean)%*%X_imp_demean/(obs_mat-1)
  
  if(cor == T){
    
    sd_mat <- matrix(data = sqrt(diag(est_cov_mat)),
                     nrow = dim(est_cov_mat)[1],
                     ncol = dim(est_cov_mat)[1])
    
    est_cor_mat <- est_cov_mat/sd_mat/t(sd_mat)
    return(est_cor_mat)
  }else{
    return(est_cov_mat)
  }
  
}
