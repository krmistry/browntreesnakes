########################### Model evaluation metrics ########################### 

obs_quarters <- unique(unlist(erad_quarters[erad_methods[c(2,3)]]))

# Function to summarize simulation data into size class N for each quarter
summed_sim_data_fun <- function(simulation_quarter_data,
                                obs_quarters) {
  real_data <- simulation_quarter_data$all_quarters
  real_data_summed <- as.data.frame(matrix(NA, ncol = (length(size_class_names)+2), 
                                           nrow = erad_quarter_time_step))
  colnames(real_data_summed) <- c("Quarter", size_class_names, "total")
  real_data_summed$Quarter <- c(1:erad_quarter_time_step)
  
  for(quarter in 1:erad_quarter_time_step) {
    for(size in 1:length(size_class_names)) {
      real_data_summed[quarter, size+1] <- nrow(real_data[real_data$Quarter == quarter & real_data$size_category == size_class_names[size], ])
    }
    real_data_summed$total[quarter] <- sum(real_data_summed[quarter,c(2:5)])
  }
  # Extracting the quarters with observations
  real_data_summed <- real_data_summed[obs_quarters,]
  
  return(real_data_summed)
}

# Test
summed_data <- summed_sim_data_fun(erad_quarter_results)


# Function to separate (and sum when appropriate) estimated values - mean N, SD of N, and credible intervals for 95th percentile
summed_est_results_fun <- function(output_jags,
                                   obs_quarters) {
  # Separating relevant raw results - mean, sd and 95% CI of N
  N_jags_output <- output_jags$mean$N
  sd_jags_output <- output_jags$sd$N
  estimated_N_95_CI_lower <- output_jags$q2.5$N
  estimated_N_95_CI_upper <- output_jags$q97.5$N
  # Setting up dataframe to hold estimated mean N values for each quarter and size
  estimated_N <- as.data.frame(matrix(NA, ncol = (length(size_class_names)+2), 
                                    nrow = length(obs_quarters)))
  colnames(estimated_N) <- c("Quarter", size_class_names, "total")
  estimated_N$Quarter <- c(obs_quarters)
  # Setting up dataframe to hold SD for estimated N values for each quarter and size
  estimated_N_sd <- as.data.frame(matrix(NA, ncol = (length(size_class_names)+1), 
                                         nrow = length(obs_quarters)))
  colnames(estimated_N_sd) <- c("Quarter", size_class_names)
  estimated_N_sd$Quarter <- c(obs_quarters)
  
  for(quarter in 1:length(obs_quarters)) {
    for(size in 1:length(size_class_names)) {
      estimated_N[quarter, size+1] <-  N_jags_output[size, 1, quarter]
      estimated_N_sd[quarter, size+1] <-  sd_jags_output[size, 1, quarter]
    }
    estimated_N$total[quarter] <- sum(estimated_N[quarter,c(2:5)])
  }

  return(list(estimated_N = estimated_N,
         estimated_N_sd = estimated_N_sd,
         estimated_N_95_CI_lower = estimated_N_95_CI_lower,
         estimated_N_95_CI_upper = estimated_N_95_CI_upper))
}

# Test
summed_results <- summed_est_results_fun(output_jags, obs_quarters)


### Function to calculate:
### - RMSE - comparing average predicted value to "real" value 
### - Percent relative difference - 2*((estimated - real)/(abs(estimated) + abs(real)))
# Calculate for each size class N and for the total N, for each replicate and then summarize results

accuracy_metrics_fun <- function(real_data,
                                 estimated_results,
                                 obs_quarters) {
  y <- list()
  y2 <- list()
  RMSE <- list()
  PRD <- list()
  for(j in 1:(length(size_class_names)+1)) {
    y[[j]] <- vector()
    y2[[j]] <- vector()
    for(i in 1:length(obs_quarters)) {
     y[[j]][i] <- (real_data[i,j+1] - estimated_results[i,j+1])^2
     y2[[j]][i] <- 2*((estimated_results[i,j+1] - real_data[i,j+1])/(abs(estimated_results[i,j+1]) + abs(real_data[i,j+1])))
    }
    RMSE[[j]] <- sqrt(sum(y[[j]])/length(obs_quarters))
    PRD[[j]] <- mean(y2[[j]])
  }
  names(RMSE) <- c(size_class_names, "total")
  names(PRD) <- c(size_class_names, "total")
  
  return(list(RMSE = RMSE, PRD = PRD))
}

# Test
accuracy_metrics <- accuracy_metrics_fun(summed_data,
                                         summed_results$estimated_N,
                                         obs_quarters)


### Function to calculate:
### percent CV -  the mean of each simulation’s standard deviation of the posterior distribution of abundance divided by the estimated mean abundance 
# Calculate for each size class N and for the total N, for each replicate and then summarize results

# Extracting estimated sd for N from jags output
percent_CV_fun <- function(estimated_N_sd,
                           obs_quarters) {
  
  y3 <- list()
  percent_CV <- list()
  for(size in 1:length(size_class_names)) {
    y3[[size]] <- vector()
    for(i in 1:length(obs_quarters)) {
      y3[[size]][i] <- estimated_N_sd[i,size+1]/estimated_N[i,size+1]
    }
    percent_CV[[size]] <- mean(y3[[size]])*100
  }
  names(percent_CV) <- size_class_names
  
  return(percent_CV)
}

# Test
percent_CV <- percent_CV_fun(summed_results$estimated_N_sd, obs_quarters)

### Function to calculate coverage (how often the value falls between the 95% credible intervals)
coverage_fun <- function(lower_CI,
                         upper_CI,
                         obs_quarters,
                         real_data) {
  y4 <- list()
  coverage <- list()
  for(size in 1:length(size_class_names)) {
    y4[[size]] <- vector()
    for(quarter in 1:length(obs_quarters)) {
      if(real_data[quarter, size+1] > estimated_N_95_CI_lower[size, 1, quarter] & 
         real_data[quarter, size+1] < estimated_N_95_CI_upper[size, 1, quarter]) {
        y4[[size]][quarter] <- 1
      } else {
        y4[[size]][quarter] <- 0
      }
    }
    # if(sum(y4[[size]]) == 0) { # If any of the quarters aren't covered, this will be 0
    #   coverage[[size]] <- 0
    # } else {
    #   coverage[[size]] <- 1
    # }
    # Percentage of quarters in which the real value is within the estimated 95% percentile for each size
    coverage[[size]] <- mean(y4[[size]])
  }
  names(coverage) <- size_class_names
  
  return(coverage)
}

# Test
coverage <- coverage_fun(summed_results$estimated_N_95_CI_lower, 
                         summed_results$estimated_N_95_CI_upper,
                         obs_quarters, real_data_summed)


### Putting all of the above functions together to produce a single list of accuracy metrics for
### a single model run

eval_metrics_fun <- function(simulation_quarter_data,
                             output_jags,
                             obs_quarters) {
  # Format simulation data
  summed_data <- summed_sim_data_fun(erad_quarter_results)
  # Format estimation results
  summed_results <- summed_est_results_fun(output_jags, obs_quarters)
  # Calculate accuracy metrics, RMSE and percent relative difference
  accuracy_metrics <- accuracy_metrics_fun(summed_data,
                                           summed_results$estimated_N,
                                           obs_quarters)
  # Calculate percent CV 
  percent_CV <- percent_CV_fun(summed_results$estimated_N_sd, obs_quarters)
  # Calculate coverage
  coverage <- coverage_fun(summed_results$estimated_N_95_CI_lower, 
                           summed_results$estimated_N_95_CI_upper,
                           obs_quarters, real_data_summed)
  
  return(list(accuracy_metrics = accuracy_metrics,
              percent_CV = percent_CV,
              coverage = coverage,
              summed_results = summed_results))
}

# Test
model_metrics <- eval_metrics_fun(erad_quarter_results,
                                  output_jags,
                                  unique(unlist(erad_quarters[erad_methods[c(2,3)]])))

model_cost <- cost_function(methods = names(erad_quarters), 
                            erad_days, 
                            erad_quarters, 
                            num_transects, 
                            num_teams, 
                            area_size)


model_metrics <- list()
for(variant in 1:num_variants) {
  model_metrics[[variant]] <- eval_metrics_fun(IBM_quarter_results,
                                               jags_output_list[[variant]],
                                               unique(unlist(erad_quarters[erad_methods[c(2,3)]])))
}



