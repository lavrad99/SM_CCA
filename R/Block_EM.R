Block_EM <- function(s, X_miss, mu_init = apply(X=X_miss,MARGIN = 2,FUN = function(x){return(mean(x,na.rm=T))}),
                     Sigma_init = diag(apply(X=X_miss,MARGIN = 2,FUN = function(x){return(var(x,na.rm=T))})),
                     maxiter = 100 , tol = 0.00001, verbose = 0){
  
  #Function implementing fast EM-estimation of the covariance matrix from missing data
  
  #s must be an object created by mat_summary
  
  mu <- mu_init
  Sigma <- Sigma_init
  
  C <- length(s)
  
  n <- dim(X_miss)[[1]]
  d <- dim(X_miss)[[2]]
  
  no_missing_block <- paste0(rep('0' , d) , collapse = '')
  
  for(iter in 1:maxiter){
    if(verbose > 0){
      print(paste('Starting Iteration' , as.character(iter)))
    }
    X_tilde <- X_miss
    adj_mat <- matrix(data = 0 , nrow = d , ncol = d)
    for(k in 1:C){
      if(verbose > 1){
        print(paste(as.character(round(k/C*100 , 2)) , '% of iteration completed' , sep =''))
      }
      if(s[[k]]$pattern ==  no_missing_block){
        next
      }
      G <- Sigma[s[[k]]$miss_pattern_cols , s[[k]]$miss_pattern_cols == F] %*% solve(Sigma[s[[k]]$miss_pattern_cols == F, s[[k]]$miss_pattern_cols == F])
      d_miss <- sum(s[[k]]$miss_pattern_cols)
      n_miss <- sum(s[[k]]$miss_pattern_rows)
      X_tilde[s[[k]]$miss_pattern_rows , s[[k]]$miss_pattern_cols] <- t(matrix(data = mu[s[[k]]$miss_pattern_cols] ,
                                                                               nrow = d_miss , ncol = n_miss))+ t((t(X_tilde[s[[k]]$miss_pattern_rows , s[[k]]$miss_pattern_cols == F])-mu[s[[k]]$miss_pattern_cols == F]))%*%t(G)
      Sigma_cond <- Sigma[s[[k]]$miss_pattern_cols , s[[k]]$miss_pattern_cols]-G%*%Sigma[s[[k]]$miss_pattern_cols==F , s[[k]]$miss_pattern_cols]
      adj_mat_k <- matrix(data = 0 , nrow = d , ncol = d)
      adj_mat_k[s[[k]]$miss_pattern_cols , s[[k]]$miss_pattern_cols] <- Sigma_cond
      adj_mat <- adj_mat + s[[k]]$prop*adj_mat_k
    }
    mu_hat <- apply(X_tilde , MARGIN = 2 , mean)
    adj_X_tilde <- t(t(X_tilde)-mu_hat)
    Sigma_hat <- t(adj_X_tilde)%*%adj_X_tilde/n  + adj_mat
    diff_score <- euc_length(as.vector(Sigma_hat[lower.tri(Sigma_hat)]) - as.vector(Sigma[lower.tri(Sigma)]))/euc_length(as.vector(Sigma[lower.tri(Sigma)]))
    if(diff_score < tol){
      ret <- list()
      ret$mu <- mu
      ret$Sigma <- Sigma
      return(ret)
    }else{
      mu <- mu_hat
      Sigma <- Sigma_hat
    }
  }
  ret <- list()
  ret$mu <- mu
  ret$Sigma <- Sigma
  return(ret)
}
