library('softImpute')

source('./elbow_rule.R') #Change '.' to folder where file is kept if necessary

SoftImpute <- function(X_miss , debias = FALSE){
  
  # Function for SoftImpute
  
  data_SVD <- impute_mean(X_miss)
  
  singular_vals <- svd(data_SVD)$d
  
  dim_elbow <- elbow_rule(y = singular_vals)
  
  lambda_elbow <- singular_vals[dim_elbow]
  
  if(debias == F){
    return(softImpute::complete(X_miss , softImpute(x = X_miss, rank.max = dim_elbow,
                                                    lambda = lambda_elbow, maxit = 300)))
  }else{
    
    return(softImpute::complete(X_miss ,
                                deBias(X_miss , softImpute(x = X_miss,
                                                           rank.max = dim_elbow,
                                                           lambda = lambda_elbow, maxit = 300))))
  }
  
}
