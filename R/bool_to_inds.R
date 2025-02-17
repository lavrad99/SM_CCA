bool_to_inds <- function(inds_bool){
  
  # Function transforming Boolean to integer indexing
  
  inds <- 1:length(inds_bool)
  return(inds[inds_bool])
}
