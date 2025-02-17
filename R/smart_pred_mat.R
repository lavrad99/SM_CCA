smart_pred_mat <- function(X , k){
  
  # Function for constructing a "smart" predictor matrix for MICE
  
  S_mat_X <- cor(X , use = 'pairwise.complete.obs')^2*(t(1-is.na(X))%*%is.na(X))
  
  top_kth_val_X <- apply(S_mat_X , MARGIN = 1, FUN = function(x){return(Large(x , min(k,dim(X)[2]))[[1]])})
  
  top_kth_mat_X <- matrix(data = top_kth_val_X , nrow = dim(X)[2] , ncol = dim(X)[2])
  
  pred_mat_X <- S_mat_X >= top_kth_mat_X
  
  return(pred_mat_X)
  
}
