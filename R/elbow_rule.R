elbow_rule <- function(y){
  
  # Function performing the elbow-rule on a curve
  
  n <- length(y)
  
  k <- (y[n] - y[1])/(n-1)
  
  m <- y[1]-k
  
  dist_vec <- abs(k*(1:n)-y+m)/sqrt(k^2+1)
  
  return(which.max(dist_vec))
}
