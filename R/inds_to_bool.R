inds_to_bool <- function(inds , d){
  
  # Function transforming integer to Boolean indexing
  
  inds_bool <- rep(FALSE , d)
  inds_bool[inds] <- TRUE
  return(inds_bool)
}
