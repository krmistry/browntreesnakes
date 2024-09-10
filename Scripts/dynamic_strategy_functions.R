#### Functions for running dynamic management strategies

quarter_name_fun <- function(quarter_time_step,
                             t){
  y <- quarter_time_step-1
  quarter_name <- vector()
  for(Y in 1:y) {
    quarter_name[Y] <- t*quarter_time_step - Y
  }
  quarter_name <- rev(quarter_name)
  quarter_name[quarter_time_step] <- t*quarter_time_step
  return(quarter_name)
}

# # Test
# quarters <- quarter_name_fun(quarter_time_step = 2,
#                  t = 3)

  