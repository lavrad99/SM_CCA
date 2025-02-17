library('tidyverse')
library('mice')
library('missMethods')

source('./geigen_CCA.R') #Change '.' to folder where file is kept if necessary
source('./Block_EM.R') #Change '.' to folder where file is kept if necessary
source('./SoftImpute.R') #Change '.' to folder where file is kept if necessary

SM_CCA <- function(X_list , Y_list, B, pre_deconf = F , n_corrs = 1,
                   method = 'MICE', n_sets = 10, n_preds = 20, maxit = 5,
                   MICE_method = 'norm', verbose = 0){
  
  # Function implementing our proposed method
  
  struct_miss_all <- rep(T , dim(X_list[[1]])[1])
  
  for(Z in X_list){
    struct_miss_temp <- apply(is.na(Z) , MARGIN = 1 , FUN = mean) == 1
    
    struct_miss_all <- struct_miss_all & struct_miss_temp
  }
  
  for(Z in Y_list){
    struct_miss_temp <- apply(is.na(Z) , MARGIN = 1 , FUN = mean) == 1
    
    struct_miss_all <- struct_miss_all & struct_miss_temp
  }
  
  if(pre_deconf == F){
    
    if(verbose > 0){
      print('Starting Deconfounding')
    }
    
    exist_miss_B <- apply(is.na(B) , MARGIN = 1 , FUN = sum) > 0
    
    remove_points <- exist_miss_B | struct_miss_all
    
    B <- B[!remove_points , ]
    
    for(i in 1:length(X_list)){
      X_list[[i]] <- X_list[[i]][!remove_points , ]
    }
    for(i in 1:length(Y_list)){
      Y_list[[i]] <- Y_list[[i]][!remove_points , ]
    }
    deconf_fun <- function(X){
      return(deconf(B=B,X=X))
    }
    
    X_list <- lapply(X_list, deconf_fun)
    if(verbose > 0){
      print('X completed')
    }
    
    Y_list <- lapply(Y_list, deconf_fun)
    
    if(verbose > 0){
      print('Y completed')
      print('Deconfounding finished')
    }
  }
  
  d_X <- sum(unlist(lapply(X_list , FUN = function(x){return(dim(x)[2])})))
  
  d_Y <- sum(unlist(lapply(Y_list , FUN = function(x){return(dim(x)[2])})))
  
  if(method == 'MICE'){
    if(verbose > 0){
      print('Starting Multiple Imputation for X')
    }
    MI_list_X <- list()
    
    SM_list_X <- list()
    
    for(my_set in 1:length(X_list)){
      X <- X_list[[my_set]]
      
      struct_miss_X <- apply(is.na(X), MARGIN = 1, sum) == dim(X)[2]
      
      SM_list_X[[my_set]] <- struct_miss_X
      
      X_block <- X[!struct_miss_X , ]
      
      if(verbose > 1){
        print(paste('Starting MI X, view' , my_set))
      }
      
      
      
      MI_list_X[[my_set]] <- mice(X_block , m = n_sets,
                                  method = MICE_method,
                                  predictorMatrix = smart_pred_mat(X=X_block,k=n_preds),
                                  maxit = maxit,printFlag = verbose>1)
    }
    
    MI_list_Y <- list()
    
    SM_list_Y <- list()
    
    for(my_set in 1:length(Y_list)){
      Y <- Y_list[[my_set]]
      
      struct_miss_Y <- apply(is.na(Y), MARGIN = 1, sum) == dim(Y)[2]
      
      SM_list_Y[[my_set]] <- struct_miss_Y
      
      Y_block <- Y[!struct_miss_Y , ]
      
      if(verbose > 1){
        print(paste('Starting MI Y, view' , my_set))
      }
      
      
      
      MI_list_Y[[my_set]] <- mice(Y_block , m = n_sets,
                                  method = MICE_method,
                                  predictorMatrix = smart_pred_mat(X=Y_block,k=n_preds),printFlag = verbose>1)
    }
    
    Sigma_chol <- matrix(data = 0, nrow = d_X+d_Y, ncol = d_X+d_Y)
    
    for(i in 1:n_sets){
      D_comb <- matrix(nrow = dim(X_list[[1]])[1], ncol = d_X+d_Y)
      
      start_ind <- 1
      
      for(j in 1:length(X_list)){
        D_comb[!SM_list_X[[j]] , start_ind:(start_ind+dim(X_list[[j]])[2]-1)] <- as.matrix(mice::complete(MI_list_X[[j]],
                                                                                                          action = i))
        start_ind <- start_ind + dim(X_list[[j]])[2]
      }
      
      for(j in 1:length(Y_list)){
        D_comb[!SM_list_Y[[j]] , start_ind:(start_ind+dim(Y_list[[j]])[2]-1)] <- as.matrix(mice::complete(MI_list_Y[[j]],
                                                                                                          action = i))
        start_ind <- start_ind + dim(Y_list[[j]])[2]
      }
      my_summary <- mat_summary(X_miss = D_comb)
      if(verbose > 0){
        print(paste('Starting EM' , i))
      }
      Sigma_chol <- Sigma_chol + chol(Block_EM(s=my_summary, X_miss = D_comb)$Sigma)
      if(verbose > 0){
        print('Finished')
      }
    }
    Sigma_chol <- Sigma_chol/n_sets
    
    Sigma <- t(Sigma_chol)%*%Sigma_chol
    
    return(geigen_CCA(Sigma = Sigma , x_inds = 1:d_X, y_inds = (d_X+1):(d_X+d_Y),n_corrs=n_corrs))
    
  }else if(method == 'SoftImpute'){
    D_comb <- matrix(nrow = dim(X_list[[1]])[1], ncol = d_X+d_Y)
    
    start_ind <- 1
    
    for(j in 1:length(X_list)){
      X <- X_list[[j]]
      
      struct_miss_X <- apply(is.na(X), MARGIN = 1, sum) == dim(X)[2]
      
      D_comb[!struct_miss_X , start_ind:(start_ind+dim(X)[2]-1)] <- SoftImpute(X[!struct_miss_X, ])
      
      start_ind <- start_ind + dim(X)[2]
    }
    
    for(j in 1:length(Y_list)){
      Y <- Y_list[[j]]
      
      struct_miss_Y <- apply(is.na(Y), MARGIN = 1, sum) == dim(Y)[2]
      
      D_comb[!struct_miss_Y , start_ind:(start_ind+dim(Y)[2]-1)] <- SoftImpute(Y[!struct_miss_Y, ])
      
      start_ind <- start_ind + dim(Y)[2]
    }
    
    my_summary <- mat_summary(X_miss = D_comb)
    
    return(geigen_CCA(Sigma = Block_EM(s = my_summary , X_miss = D_comb)$Sigma
                      , x_inds = 1:d_X, y_inds = (d_X+1):(d_X+d_Y),n_corrs=n_corrs))
  }else if(method == 'mean'){
    
    D_comb <- matrix(nrow = dim(X_list[[1]])[1], ncol = d_X+d_Y)
    
    start_ind <- 1
    
    for(j in 1:length(X_list)){
      X <- X_list[[j]]
      
      struct_miss_X <- apply(is.na(X), MARGIN = 1, sum) == dim(X)[2]
      
      D_comb[!struct_miss_X , start_ind:(start_ind+dim(X)[2]-1)] <- impute_mean(X[!struct_miss_X, ])
      
      start_ind <- start_ind + dim(X)[2]
    }
    
    for(j in 1:length(Y_list)){
      Y <- Y_list[[j]]
      
      struct_miss_Y <- apply(is.na(Y), MARGIN = 1, sum) == dim(Y)[2]
      
      D_comb[!struct_miss_Y , start_ind:(start_ind+dim(Y)[2]-1)] <- impute_mean(Y[!struct_miss_Y, ])
      
      start_ind <- start_ind + dim(Y)[2]
    }
    
    my_summary <- mat_summary(X_miss = D_comb)
    
    return(geigen_CCA(Sigma = Block_EM(s = my_summary , X_miss = D_comb)$Sigma
                      , x_inds = 1:d_X, y_inds = (d_X+1):(d_X+d_Y),n_corrs=n_corrs))
  }
  
}
