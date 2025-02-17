library('tidyverse')
library('geigen')

geigen_CCA <- function(Sigma , x_inds , y_inds, b_inds = NA, n_corrs = 1){
  
  # Function performing CCA on a covariance matrix using generalised eigenvalue decomposition
  
  # x_inds represents columns belonging to vector x, y_inds to vecor y and b_inds to the confounds
  
  if(sum(is.na(b_inds)) > 0){
    
    my_solution <- list()
    
    d <- dim(Sigma)[1]
    
    if(length(x_inds) < d){
      x_inds <- inds_to_bool(x_inds , d)
    }
    
    if(length(y_inds) < d){
      y_inds <- inds_to_bool(y_inds , d)
    }
    
    A <- matrix(data = 0 , nrow = d, ncol = d)
    
    A[x_inds , y_inds] <- Sigma[x_inds , y_inds]
    
    A[y_inds , x_inds] <- Sigma[y_inds , x_inds]
    
    B <- matrix(data = 0 , nrow = d, ncol = d)
    
    B[x_inds , x_inds] <- Sigma[x_inds , x_inds]
    
    B[y_inds , y_inds] <- Sigma[y_inds , y_inds]
    
    my_geigen <- geigen(A = A , B = B)
    
    d_X <- sum(x_inds)
    
    d_Y <- sum(y_inds)
    
    rho_list <- vector(mode = 'numeric' , length = n_corrs)
    
    beta_hat_X <- matrix(nrow = d_X, ncol = n_corrs)
    
    beta_hat_Y <- matrix(nrow = d_Y, ncol = n_corrs)
    
    for(i in 1:n_corrs){
      x_coeffs <- my_geigen$vectors[ , i][x_inds]/euc_length(my_geigen$vectors[ , i][x_inds])
      
      x_coeffs <- x_coeffs/sign(x_coeffs[1])
      
      y_coeffs <- my_geigen$vectors[ , i][y_inds]/euc_length(my_geigen$vectors[ , i][y_inds])
      
      sd_x <- sqrt(t(x_coeffs) %*% Sigma[x_inds,x_inds]%*% x_coeffs)
      
      sd_y <- sqrt(t(y_coeffs) %*% Sigma[y_inds,y_inds]%*% y_coeffs)
      
      rho <- as.numeric(t(x_coeffs) %*% Sigma[x_inds , y_inds]%*% y_coeffs/sd_x/sd_y)
      
      y_coeffs <- y_coeffs / sign(rho)
      
      rho <- rho / sign(rho)
      
      rho_list[i] <- rho
      
      beta_hat_X[ , i] <- x_coeffs
      
      beta_hat_Y[ , i] <- y_coeffs
      
    }
    
    
    my_solution$beta_hat_X <- beta_hat_X
    
    my_solution$beta_hat_Y<- beta_hat_Y
    
    my_solution$rho_list <- rho_list
    
    return(my_solution)
    
  }else{
    my_solution <- list()
    
    d <- dim(Sigma)[1]
    
    if(length(x_inds) < d){
      x_inds <- inds_to_bool(x_inds , d)
    }
    
    if(length(y_inds) < d){
      y_inds <- inds_to_bool(y_inds , d)
    }
    if(length(b_inds) < d){
      b_inds <- inds_to_bool(b_inds , d)
    }
    
    xy_inds = c(bool_to_inds(x_inds) , bool_to_inds(y_inds))
    
    Sigma <- Sigma[xy_inds,xy_inds] - matrix(data = Sigma[xy_inds,b_inds], nrow = length(xy_inds), ncol = sum(b_inds))%*%solve(Sigma[b_inds,b_inds])%*%matrix(data = Sigma[b_inds,xy_inds], nrow = sum(b_inds), ncol = length(xy_inds))
    
    d <- dim(Sigma)[1]
    n_xs <- sum(x_inds)
    
    x_inds <- rep(F,d)
    y_inds <- rep(F,d)
    
    x_inds[1:n_xs] <- T
    
    y_inds[(1+n_xs):d] <- T
    
    A <- matrix(data = 0 , nrow = d, ncol = d)
    
    A[x_inds , y_inds] <- Sigma[x_inds , y_inds]
    
    A[y_inds , x_inds] <- Sigma[y_inds , x_inds]
    
    B <- matrix(data = 0 , nrow = d, ncol = d)
    
    B[x_inds , x_inds] <- Sigma[x_inds , x_inds]
    
    B[y_inds , y_inds] <- Sigma[y_inds , y_inds]
    
    my_geigen <- geigen(A = A , B = B)
    
    d_X <- sum(x_inds)
    
    d_Y <- sum(y_inds)
    
    rho_list <- vector(mode = 'numeric' , length = n_corrs)
    
    beta_hat_X <- matrix(nrow = d_X, ncol = n_corrs)
    
    beta_hat_Y <- matrix(nrow = d_Y, ncol = n_corrs)
    
    for(i in 1:n_corrs){
      x_coeffs <- my_geigen$vectors[ , i][x_inds]/euc_length(my_geigen$vectors[ , i][x_inds])
      
      x_coeffs <- x_coeffs/sign(x_coeffs[1])
      
      y_coeffs <- my_geigen$vectors[ , i][y_inds]/euc_length(my_geigen$vectors[ , i][y_inds])
      
      sd_x <- sqrt(t(x_coeffs) %*% Sigma[x_inds,x_inds]%*% x_coeffs)
      
      sd_y <- sqrt(t(y_coeffs) %*% Sigma[y_inds,y_inds]%*% y_coeffs)
      
      rho <- as.numeric(t(x_coeffs) %*% Sigma[x_inds , y_inds]%*% y_coeffs/sd_x/sd_y)
      
      y_coeffs <- y_coeffs / sign(rho)
      
      rho <- rho / sign(rho)
      
      rho_list[i] <- rho
      
      beta_hat_X[ , i] <- x_coeffs
      
      beta_hat_Y[ , i] <- y_coeffs
      
    }
    
    
    my_solution$beta_hat_X <- beta_hat_X
    
    my_solution$beta_hat_Y<- beta_hat_Y
    
    my_solution$rho_list <- rho_list
    
    return(my_solution)
    
  }
  
}
