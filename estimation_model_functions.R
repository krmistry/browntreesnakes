### Functions related to estimation model

## Function to create data inputs for estimation model
estimation_inputs_fun <- function(observed_list,
                                  effort_list,
                                  #erad_days,
                                  parameters = c("N","p", "N.sum", "r1", "r2", "r3", "r4"),
                                  N_prior = K/2, # for version of model with uniform N_prior
                                  size_class = size_class_names) {
  
  # Estimation model set up
  # Isolating methods used in this version & the quarters:
  all_effort <- bind_rows(effort_list)
  observed_effort <- all_effort[all_effort$method %in% erad_methods[c(2:3)],]
  obs_methods <- unique(observed_effort$method) 
  obs_quarters <- unique(observed_effort$quarter) # double check that this works when there are empty dataframes from quarters that didn't have observations
  obs_days <- unique(observed_effort$day)
  
  # Values needed for array dimensions & loops
  S <- length(size_class) # number of size classes 
  Q <- length(obs_quarters) # Primary sampling periods (quarters)
  J <- length(obs_methods)
  I <- vector()
  for(i in 1:Q) {
    #I[i] <- length(unique(unlist(erad_days[[i]][obs_methods])))*length(obs_methods)
    I[i] <- length(obs_days)*length(obs_methods)
  }  
  # Secondary sampling periods (days within each quarter) - redo this later, once I re-do how erad_days is formatted to allow it to vary between quarters
  
  # Set up array of method days (if more than 1 method occurs, then each method gets its own day in the model, even if they occur on the same day in real life)
  method_days <- list()
  for(j in 1:J) {
    method_days[[j]] <- seq(j, I[1], J)
  }
  
  
  # Reformatting observations into effort and removals arrays for model
  erad_reformatted_v2 <- all_observations_fun(all_observed_list = observed_list,
                                              all_effort_list = effort_list,
                                              all_days = obs_days,
                                              all_quarters = obs_quarters,
                                              methods = obs_methods,
                                              num_methods = J)
  
  removals_array_v2 <- array(dim = c(S, I[1], Q))
  effort_array_v2 <- array(dim = c(I[1], Q))
  for(j in 1:J) { # check that the code below works for 2 methods (only used it with one so far)
    for(day in 1:length(method_days[[j]])) {
      removals_array_v2[,method_days[[j]][day],] <- erad_reformatted_v2$observation[j,,day,]
      effort_array_v2[method_days[[j]][day],] <- erad_reformatted_v2$effort[j,day,]
    }
  }
  
  ## Days between primary sampling periods when methods are used 
  # With current set up, assumed that the days in each quarter are the same 
  date_diff <- vector()
  for(t in 1:(Q-1)) {
    date_diff[t] <- ((obs_quarters[t+1]-1)*91 + obs_days[1]) - 
      ((obs_quarters[t]-1)*91 + max(max(obs_days)))
  }
  
  # Initial values for N (i.e. N[k, 1, 1])
  Y <- removals_array_v2
  N.base <- array(NA_real_, dim = c(S, I[1], Q))
  for(k in 1:S) {
    N.base[k,1,1] <- sum(Y[k,,1:Q])
  }
  
  # initialize D to > than the number that will be removed in the following year
  Y.remove1 <- vector()
  Y.remove2 <- vector()
  Y.remove3 <- vector()
  Y.remove4 <- vector()
  for(t in 1:Q){
    Y.remove1[t] <- sum(Y[1,,t])
    Y.remove2[t] <- sum(Y[2,,t])
    Y.remove3[t] <- sum(Y[3,,t])
    Y.remove4[t] <- sum(Y[4,,t])
  }
  
  D.init <- array(NA,dim = c(S,(Q-1)))
  for(t in 1:(Q-1)){
    D.init[1,t] <- Y.remove1[t+1] + 1
    D.init[2,t] <- Y.remove2[t+1] + 1
    D.init[3,t] <- Y.remove3[t+1] + 1
    D.init[4,t] <- Y.remove4[t+1] + 1
  }
  
  # Initial values for select parameters
  inits <- function (){
    list(D = D.init)
  }
  
  # Convert method_days list into an array (take off the first day, which is set differently in the model)
  j_days <- array(NA, dim = c(length(method_days), (length(method_days[[1]])-1)))
  for(j in 1:J) {
    j_days[j,] <- method_days[[j]][-1]
  }
  
  
  # Bundle data together
  data <- list(Y = removals_array_v2,
               S = S, I = I[1], Q = Q, J = J,
               xi = effort_array_v2, 
               days_btwn = date_diff, 
               N.base = N.base,
               method_days = j_days,
               N_prior = N_prior)
  
  # Parameters monitored
  parameters <- parameters
  
  return(list(data = data,
              inits = inits,
              parameters = parameters,
              obs_quarters = obs_quarters,
              obs_days = obs_days))
}

  
  
  
  