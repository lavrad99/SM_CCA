library('tidyverse')

mat_summary <- function(X_miss){
  
  # Function summarising the missingness patterns of data matrix
  
  my_summary_list <- vector(mode = 'list')
  
  n <- dim(X_miss)[[1]]
  d <- dim(X_miss)[[2]]
  
  patterns_by_row <- is.na(X_miss) %>% 
    apply(MARGIN = 1 , function(x){paste0(as.integer(x),collapse = '')})
  
  my_patterns <- unique(patterns_by_row)
  
  for(c in 1:length(my_patterns)){
    my_summary_list[[c]] <- vector(mode = 'list')
    my_summary_list[[c]]$pattern <- my_patterns[c]
    my_summary_list[[c]]$miss_pattern_cols <- strsplit(my_patterns[c] ,
                                                       split = '')[[1]] == '1'
    my_summary_list[[c]]$miss_pattern_rows <- patterns_by_row == my_patterns[c]
    my_summary_list[[c]]$prop <- mean(my_summary_list[[c]]$miss_pattern_rows)
    
  }
  return(my_summary_list)
}
