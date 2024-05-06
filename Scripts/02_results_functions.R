################ Functions to format & validate model results ###################

### Function to classify SVL into size categories

# Adding size distribution to combined quarters dataframe
size_class_fun <- function(snake_SVL,
                           size_limits = size_class_limits){
  if(snake_SVL <= size_limits[1,2]) {
    return("small")
  } else if(snake_SVL > size_limits[2,1] & snake_SVL <= size_limits[2,2]) {
    return("medium")
  } else if(snake_SVL > size_limits[3,1] & snake_SVL <= size_limits[3,2]) {
    return("large")
  } else if (snake_SVL > size_limits[4, 1]) {
    return("xlarge")
  } 
}


### Function to format quarter time series results for plotting by quarter
## Inputs:
# - quarter_results = list of dataframes with the time series results for quarters
results_format_fun <- function(quarter_results,
                               size_limits) {
  # Melt quarter timeseries into a single dataframe
  all_quarters <- melt(quarter_results, id.vars = c("ID", "sex", "repro_prob", "growth_quant"))
  colnames(all_quarters)[6:7] <- c("SVL", "Quarter")
  
  # Adding size distribution to combined quarters dataframe
  for(snake in 1:nrow(all_quarters)) {
    all_quarters$size_category[snake] <- size_class_fun(all_quarters$SVL[snake],
                                                        size_limits = size_limits)
  }
  # Factoring size category column for graphing
  all_quarters$size_category <- factor(all_quarters$size_category, 
                                       levels = c("small", "medium", "large", "xlarge"))
  
  # Return results 
  return(all_quarters)
}


### Function to create model validation summary stats


model_validation_fun <- function(quarter_timeseries,
                                 all_quarters,
                                 repro_females,
                                 total_quarters,
                                 size_classes) {
  # Create empty list to store results
  model_validation <- list()
  # Proportion of females reproduced each quarter
  model_validation$mom_rate <- vector()
  # Population growth rate each quarter
  model_validation$pop_growth_rate <- vector()
  # Mortality rate per quarter
  model_validation$mortality_rate <- vector()
  # Calculating the above rates for each quarter
  for(quarter in 1:total_quarters) {
    model_validation$mom_rate[quarter] <- nrow(repro_females[[quarter]])/nrow(quarter_timeseries[[quarter + 1]][quarter_timeseries[[quarter + 1]]$sex == "F", ])
    model_validation$pop_growth_rate[quarter] <- nrow(quarter_timeseries[[quarter+1]][!quarter_timeseries[[quarter+1]]$ID %in% quarter_timeseries[[quarter]]$ID,])/nrow(quarter_timeseries[[quarter+1]])
    model_validation$mortality_rate[quarter] <- nrow(quarter_timeseries[[quarter]][!quarter_timeseries[[quarter]]$ID %in% quarter_timeseries[[quarter + 1]]$ID,])/nrow(quarter_timeseries[[quarter]])
  }
  # Lifespan of individuals born within the simulation
  model_validation$lifespan <- vector()
  # Maximum size each snake reached
  model_validation$max_size <- vector()
  # Isolating snakes born in the simulation
  sim_born_snakes <- all_quarters[-grep("A", all_quarters$ID),]
  # Number of quarters that snakes stay in each size class
  age_by_size <- as.data.frame(matrix(NA, nrow = length(unique(sim_born_snakes$ID)), 
                                      ncol = 5))
  colnames(age_by_size) <- c("ID", "num_Q_small", "num_Q_medium", "num_Q_large",
                             "num_Q_xlarge")
  # Calculating above values for all simulation born snakes
  for(snake in 1:length(unique(sim_born_snakes$ID))) {
    model_validation$lifespan[snake] <- nrow(sim_born_snakes[sim_born_snakes$ID == sim_born_snakes$ID[snake],])
    model_validation$max_size[snake] <- max(sim_born_snakes$SVL[sim_born_snakes$ID == sim_born_snakes$ID[snake]])
    snake_data <- sim_born_snakes[sim_born_snakes$ID == unique(sim_born_snakes$ID)[snake], ] 
    age_by_size[snake, 1] <- snake_data$ID[1]
    for(size in 1:length(size_classes)) {
      age_by_size[snake, (size + 1)] <- nrow(snake_data[snake_data$size_category == size_classes[size],])
    }
  }
  # Average number of quarters that simulation born snakes spend in each size class
  model_validation$average_age_by_size <- vector()
  for(size in 1:length(size_classes)) {
    model_validation$average_age_by_size[size] <- mean(age_by_size[, (size + 1)])
  }
  # Proportion of snakes that lived longer than 15 years
  model_validation$prop_over_60_quart <- length(model_validation$lifespan[model_validation$lifespan > 60])/length(model_validation$lifespan)
  # Isolating snakes who live longer than 10 years to look at their growth trajectories
  model_validation$LL_snakes <- all_quarters[0,]
  for(snake in 1:length(unique(all_quarters$ID))) {
    snake_ts <- all_quarters[all_quarters$ID == unique(all_quarters$ID)[snake], ]
    if(nrow(snake_ts) > 40) {
      model_validation$LL_snakes <- rbind(model_validation$LL_snakes, snake_ts)
    } 
  }
  # Proportion of entire population that lives more than 10 years
  model_validation$LL_snakes_prop <- length(unique(model_validation$LL_snakes$ID))/length(unique(all_quarters$ID))
  # Sampling 25 long lived snakes and graphing their growth
  if(length(unique(model_validation$LL_snakes$ID)) >= 25) {
    plot_snake_IDs <- sample(unique(model_validation$LL_snakes$ID), 25)
    plot_snakes <- all_quarters[all_quarters$ID %in% plot_snake_IDs,]
    model_validation$LL_plot <- ggplot(plot_snakes, aes(x = Quarter, y = SVL, color = ID)) +
      geom_path()
  }
  # Summarizing many of the above states to erport mean, min and max values
  mv_table <- as.data.frame(matrix(NA, nrow = 3, ncol = 5))
  colnames(mv_table) <- c("mom_rate_per_quarter", "pop_growth_rate_per_quarter", 
                          "mortality_rate_per_quarter", "lifespan_in_quarters", "max_SVL")
  rownames(mv_table) <- c("min", "mean", "max")
  mv_table$mom_rate_per_quarter <- summary(model_validation$mom_rate)[c(1,4,6)]
  mv_table$pop_growth_rate_per_quarter <- summary(model_validation$pop_growth_rate)[c(1,4,6)]
  mv_table$mortality_rate_per_quarter <- summary(model_validation$mortality_rate)[c(1,4,6)]
  mv_table$lifespan_in_quarters <- summary(model_validation$lifespan)[c(1,4,6)]
  mv_table$max_SVL <- summary(model_validation$max_size)[c(1,4,6)]
  # Returning both full results and summarized results
  return(list(model_validation = model_validation, model_val_summary = mv_table))
}


                          
### Transforming observed dead snakes into summary dataframes (# of snakes) for inputting
### into the estimation model

all_observations_fun <- function(all_observed_list,
                                 all_effort_list,
                                 all_days,
                                 all_quarters,
                                 methods,
                                 num_methods) {
  # Create vectors with quarters and days when eradication occurred
  # all_quarters <- sort(unique(unlist(erad_quarters[methods])))
  # all_days <- sort(unique(unlist(erad_days[methods])))
  # Combine all observed (captured) snakes into one dataframe
  observed_list <- all_observed_list[unlist(lapply(all_observed_list, length) != 0)]
  # Remove empty lists from observed_list (when no snakes were found in a quarter)
  observed_list <- list_drop_empty(observed_list)
  # Combine all observations into one dataframe
  all_observed <- bind_rows(observed_list)
  for(snake in 1:nrow(all_observed)){
    all_observed$size_category[snake] <- size_class_fun(all_observed$SVL[snake], 
                                                        size_class_limits)
  }
  # Creating empty array to summarize results (for input into jags model)
  observations_array <- array(0, dim = c(num_methods, 4, length(all_days), length(all_quarters)))
  # Effort
  effort_list <- all_effort_list[unlist(lapply(all_effort_list, length) != 0)]
  all_effort <- bind_rows(effort_list)
  effort_array <- array(0, dim = c(num_methods,length(all_days), length(all_quarters)))
  # Count the number of observed snake for each quarter & day that effort occurred
  # for(quarter in all_quarters) {
  #   for(day in all_days) {
  for(quarter in 1:length(all_quarters)) {
    for(day in 1:length(all_days)) {
      for(method in 1:length(methods)) {
        for(size in 1:length(size_class_names)) {
          # Summarizing observations
          observations_array[method,size,day,quarter] <-  nrow(all_observed[all_observed$size_category == size_class_names[size] & all_observed$day == all_days[day] & all_observed$quarter == all_quarters[quarter] & all_observed$method == methods[method],])
        }
        # Summarizing effort
        if(length(all_effort$effort[all_effort$method == methods[method] & all_effort$day == all_days[day] & all_effort$quarter == all_quarters[quarter]]) > 0) {
          effort_array[method, day, quarter] <- all_effort$effort[all_effort$method == methods[method] & all_effort$day == all_days[day] & all_effort$quarter == all_quarters[quarter]]
        } else {
          effort_array[method, day, quarter] <- 0
        }
      }
    }
  }
  # If more than 1 method used, combine methods to make array one dimension smaller
  
  observations_array <- 
  
  return(list(observation = observations_array,
              effort = effort_array))
}


# # Test
# q <- all_observations_fun(erad_results_ts = erad_quarter_results,
#                           methods = erad_methods[c(2:3)])


## Create vector with the difference between all eradication effort days 


# 
# effort_days <- list()
# for(method in 1:length(erad_methods)) {
#   effort_days[[method]] <- list()
#   for(quarter in 1:length(erad_quarters[[method]])) {
#     effort_days[[method]][[quarter]] <- vector()
#     for(day in 1:length(erad_days[[method]])) {
#       effort_days[[method]][[quarter]][day] <- (erad_quarters[[method]][quarter]*91)-91 + erad_days[[method]][day]
#     }
#   }
#   names(effort_days[[method]]) <- paste0("quarter_", erad_quarters[[method]])
# }
# names(effort_days) <- erad_methods
# all_effort_days <- sort(unique(unlist(effort_days)))
# # Difference between last visual survey day in quarter 2 week 5 to first visual survey day 
# # in week 9
# effort_days$visual$quarter_2[14] - effort_days$visual$quarter_2[7]
# effort_days$visual$quarter_7[7] - effort_days$visual$quarter_2[14]
# effort_days$visual$quarter_7[14] - effort_days$visual$quarter_7[7]
# 
# all_days_btwn <- vector()
# for(t in 1:(length(all_effort_days)-1)) {
#   all_days_btwn[t] <- all_effort_days[t+1]-all_effort_days[t]
# }



###### Function to calculate true vital rates (survival, fecundity, size transition) between eradication 
#      quarters
# In current version, there are 3 intervals; between quarters 2 & 3, 3 & 6 and 6 & 7

true_vital_rates_v1_fun <- function(all_erad_quarters,
                                    erad_results_df,
                                    sizes = size_class_names) {
  # Eradication quarters when visual survey or trapping occurred
  observed_quarters <- sort(unique(unlist(all_erad_quarters)))
  n.inter_primary <- length(observed_quarters)-1
  # Beginning and end of reproductive quarters between eradication quarters
  repro_quarters <- list()
  for(quarter in 1:n.inter_primary) {
    repro_quarters[[quarter]] <- c(observed_quarters[quarter], observed_quarters[quarter+1]-1)
  }
  # Inter-primary quarter column names (used by multiple data frames)
  inter_primary_colnames <- c("size_class", paste0("inter-primary_", c(1:n.inter_primary)))
  
  ## Create data frames to hold vital rates 
  # Survival (size class X inter-primary periods)
  survival <- as.data.frame(matrix(NA, nrow = length(sizes), ncol = length(observed_quarters)))
  colnames(survival) <- inter_primary_colnames
  survival$size_class <- sizes
  # Size transition ((size class - 1) X inter-primary periods)
  size_transition <- as.data.frame(matrix(NA, nrow = (length(sizes)-1), ncol = length(observed_quarters)))
  colnames(size_transition) <- inter_primary_colnames
  size_transition$size_class <- sizes[-1]
  ## Fecundity components:
  # Part 1: separating out moms from each quarter between primary periods (list of individual snakes in each inter-primary period
  pre_erad_moms <- list()
  # Part 2: Proportion of total moms in each size class ((size class - 1) X inter-primary period)
  mom_prop <- as.data.frame(matrix(NA, nrow = (length(sizes)-1), ncol = length(observed_quarters)))
  colnames(mom_prop) <- inter_primary_colnames
  mom_prop$size_class <- sizes[-1]
  # Part 3: Fecundity ((size class - 1) X inter-primary periods)
  fecundity <- as.data.frame(matrix(NA, nrow = (length(sizes)-1), ncol = length(observed_quarters)))
  colnames(fecundity) <- inter_primary_colnames
  fecundity$size_class <- sizes[-1]
  
  for(quarter in 1:n.inter_primary) {
    # Separating out population at the beginning of the inter-primary period (after an eradication primary period)
    after_erad_pop <- erad_results_df$all_quarters_before_after_erad$pop_after_erad[observed_quarters[quarter]][[1]]
    # Separating out population at the end of the inter-primary period (before the next eradication primary period)
    before_erad_pop <- erad_results_df$all_quarters_before_after_erad$pop_before_erad[observed_quarters[quarter+1]][[1]]
    # Adding size class column to both of the above data frames
    for(snake in 1:nrow(after_erad_pop)) {
      after_erad_pop$size_class[snake] <- size_class_fun(after_erad_pop$SVL[snake])
    }
    for(snake in 1:nrow(before_erad_pop)) {
      before_erad_pop$size_class[snake] <- size_class_fun(before_erad_pop$SVL[snake])
    }
    # Fecundity part 1: adding together all moms from eradication quarter and any other quarters before next eradication quarter
    all_quarters <- c(repro_quarters[[quarter]][1]:repro_quarters[[quarter]][2])
    pre_erad_moms[[quarter]] <- erad_results_df$all_repro_females[[1]][0,]
    for(inter_quarter in all_quarters) {
      pre_erad_moms[[quarter]] <- rbind(pre_erad_moms[[quarter]],
                                        erad_results_df$all_repro_females[inter_quarter][[1]])
    }
    # Adding size class for each mom
    for(snake in 1:nrow(pre_erad_moms[[quarter]])) {
      pre_erad_moms[[quarter]]$size_class[snake] <- size_class_fun(pre_erad_moms[[quarter]]$SVL[snake])
    }
    # Calculating following rates by size class
    for(size in 1:length(sizes)) {
      size_after_pop <- after_erad_pop[after_erad_pop$size_class == sizes[size],]
      #same_size_before_pop <- before_erad_pop[before_erad_pop$size_class == sizes[size],] # for survival
      # Survival
      survival[size, quarter+1] <- sum(size_after_pop$ID %in% before_erad_pop$ID)/nrow(size_after_pop)
      if(size <= 3) {
        next_size_before_pop <- before_erad_pop[before_erad_pop$size_class == sizes[size + 1],] # for size transition
        # Size transition
        size_transition[size, quarter+1] <- sum(size_after_pop$ID %in% next_size_before_pop$ID)/nrow(size_after_pop)
        # Fecundity part 2
        mom_prop[size, quarter+1] <- nrow(pre_erad_moms[[quarter]][pre_erad_moms[[quarter]]$size_class == sizes[size + 1],])/nrow(pre_erad_moms[[quarter]])
      }
    }
    # Fecundity part 3: Separating out offspring per inter-primary period, multiplying by mom_prop per size class and dividing by the total in each size class
    offspring <- erad_results_df$all_quarters[0,]
    inter_q_pop <- erad_results_df$all_quarters[0,]
    for(inter_q in 1:length(all_quarters)) {
      first_quarter <- erad_results_df$quarter_timeseries[[all_quarters[inter_q]]]
      next_quarter <- erad_results_df$quarter_timeseries[[all_quarters[inter_q]+1]]
      # Adding offspring individuals from each inter-quarter period
      offspring_IDs <- which(!(next_quarter$ID %in% first_quarter$ID))
      offspring <- rbind(offspring, next_quarter[offspring_IDs,])
      
      # Creating inter-quarter population
      inter_q_pop <- rbind(inter_q_pop, first_quarter)
    }
    # Adding size class categories to inter-quarter population
    for(snake in 1:nrow(inter_q_pop)) {
      inter_q_pop$size_class[snake] <- size_class_fun(inter_q_pop$SVL[snake])
    }
    # Removing duplicate individuals in inter-quarter population
    inter_q_pop <- inter_q_pop %>% distinct(ID, .keep_all = TRUE)
    for(size in 1:(length(sizes) - 1)) {
      # Separating out size class at beginning of reproduction in quarter
      size_inter_q_pop_N <- nrow(inter_q_pop[inter_q_pop$size_class == sizes[size+1],])
      # For each size, estimate the (proportion of moms per size*offspring)/size_pop
      fecundity[size, quarter+1] <- (mom_prop[size,quarter+1]*nrow(offspring))/size_inter_q_pop_N
    }
    # Checking to see how many snakes transitioned through 2 size classes in the inter-primary periods 
    # (to check to see if we're missing some by not including those transition parameters in the estimation model)
    dups_w_size_change <- vector(length = n.inter_primary)
    dups <- 0
    for(snake in 1:length(unique(inter_q_pop$ID))) {
      snake_data <- inter_q_pop[inter_q_pop$ID == unique(inter_q_pop$ID)[snake],]
      if(length(unique(snake_data$size_class)) > 1) {
        dups <- dups + 1
      } else {
        dups <- dups
      }
    }
    dups_w_size_change[quarter] <- dups/nrow(inter_q_pop)
  }
  names(dups_w_size_change) <- inter_primary_colnames[-1]
  
  return(list(observed_quarters = observed_quarters,
              survival = survival,
              size_transition = size_transition,
              fecundity = fecundity,
              dups_w_size_change = dups_w_size_change))
}
 
# # Test
# x <- true_vital_rates_v1_fun(all_erad_quarters = erad_quarters[c("visual","trap")],
#                           erad_results_df = erad_quarter_results)


## Function to create plots of mean N estimates from removal-based estimation model

estimated_N_plots <- function(jags_output,
                              all_quarters,
                              observed_erad_quarters) {
  # Separating out quarters with observable removals
  #observed_erad_quarters <- sort(unique(unlist(erad_quarters[c(2,3)])))
  
  mean_N <- as.data.frame(matrix(0, nrow = Q, ncol = S))
  colnames(mean_N) <- size_class_names
  mean_N$Quarter <- observed_erad_quarters
  for(size in 1:length(size_class_names)) {
    for(quarter in 1:length(observed_erad_quarters)) {
      mean_N[quarter, size] <- jags_output$mean$N[size,I[quarter],quarter]
    }
  }
  
  # Making data long for plotting
  mean_N_long <- melt(mean_N, id.vars = "Quarter")
  colnames(mean_N_long)[2:3] <- c("size_class", "N")
  # Adding lower and upper bounds for 95% confidence interval
  for(size in 1:length(size_class_names)) {
    for(quarter in 1:length(observed_erad_quarters)) {
      row_counter <- quarter + (size-1)*length(observed_erad_quarters)
      mean_N_long$lower_bound[row_counter] <- jags_output$q2.5$N[size,I[quarter],quarter]
      mean_N_long$upper_bound[row_counter] <- jags_output$q97.5$N[size,I[quarter],quarter]
    }
  }

  ## Adding total # of captured snakes and total # of alive simulated snakes in each size
  ## class in each quarter to graph against the predicted abundance
  mean_N_long$source <- "estimated"
  real_N <- mean_N
  for(size in size_class_names) {
    for(quarter in 1:length(observed_erad_quarters)) {
      real_N[quarter, size] <- nrow(all_quarters[all_quarters$Quarter == (observed_erad_quarters[quarter]+1) & 
                                                                        all_quarters$size_category == size,]) 
    }
  }
  
  
  real_N_long <- melt(real_N, id.vars = "Quarter")
  colnames(real_N_long)[2:3] <- c("size_class", "N")
  real_N_long$source <- "simulated"
  
  # Adding NAs for upper and lower confidence intervals, so it has the same columns as mean_N_long
  real_N_long$lower_bound <- NA
  real_N_long$upper_bound <- NA
  
  # Combining simulated, removed and estimated for each size class in each quarter
  est_vs_sim_plot_1 <- ggplot(data = rbind(mean_N_long, real_N_long), aes(x = Quarter, y = N, color = source)) +
    geom_line(size = 2) +
    geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), alpha = 0.3, linetype = 0, fill = hue_pal()(1)) +
    facet_wrap("size_class", scales = "free_y") +
    guides(fill = "none") +
    theme_bw()
  
  
  ## Plotting total estimated N vs total simulation N
  est_vs_sim_plot_2 <- ggplot(rbind(mean_N_long, real_N_long), aes(y = N, x = Quarter, fill = source)) +
    geom_col(position = "dodge") +
    theme_bw()
  
  return(list(est_vs_sim_plot_1 = est_vs_sim_plot_1,
              est_vs_sim_plot_2 = est_vs_sim_plot_2,
              mean_N_df = mean_N_long))
}

# # Test
# y <- estimated_N_plots(jags_output = output_jags,
#                        erad_quarter_results = erad_quarter_results,
#                        erad_quarters = erad_quarters)


