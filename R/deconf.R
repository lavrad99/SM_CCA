deconf <- function(B , X){
  
  # Function used for deconfounding
  
  X_deconf <- X
  
  for(i in 1:(dim(X)[2])){
    obs_indic <- !is.na(X[ , i])
    
    temp_df <- as.data.frame(cbind(B[obs_indic , ] , X[obs_indic , i]))
    
    colnames(temp_df)[dim(B)[2]+1] <- 'y'
    
    my_OLS <- lm(y ~ . , data = temp_df)
    
    X_deconf[obs_indic , i] <- residuals(my_OLS) 
    
  }
  return(X_deconf)
}
