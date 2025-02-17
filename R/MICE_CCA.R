library('tidyverse')
library('mice')

source('./smart_pred_mat.R') #Change '.' to folder where file is kept if necessary


MICE_CCA <- function(D_miss,x_inds,y_inds,n_corrs = 1,n_sets = 10, n_preds = 20,
                     MICE_method = 'norm', maxit = 5, verbose = F){
  
  # Function implementing CCA from missing data using MICE
  
  my_pred_mat <- smart_pred_mat(X = D_miss, k = n_preds)
  
  my_MI <- mice(data = D_miss, m = n_sets, method = MICE_method,
                predictorMatrix = my_pred_mat, maxit = maxit, printFlag = verbose)
  
  Sigma_chol <- matrix(data = 0, nrow = dim(D_miss)[2],ncol = dim(D_miss)[2])
  
  for(i in 1:n_sets){
    Sigma_chol <- Sigma_chol+chol(cov(mice::complete(my_MI, action = i)))
  }
  
  Sigma_chol <- Sigma_chol/n_sets
  
  Sigma <- t(Sigma_chol) %*% Sigma_chol
  
  return(geigen_CCA(Sigma = Sigma, x_inds = x_inds, y_inds = y_inds,n_corrs=n_corrs))
}
