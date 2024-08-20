

# FUnction to create all quarters from the original IBM outputs (shouldn't be necessary in future runs)
recreate_IBM_all_quarters <- function(strategy_name,
                                      P,
                                      D,
                                      permutation_name,
                                      variant) {
  IBM_list <- list()
  for(y in 1:total_time_steps[variant]) {
    IBM_output <- readRDS(paste0(results_folders[[strategy_name]][[P]][[D]]["IBM"], "/variant_",
                                 variant, "/IBM_results_start_pop_", permutation_name, "_variant-",
                                 variant,"_IBM_set-",
                                 y, ".rds"))
    IBM_list <- c(IBM_list, IBM_output[[y]]$quarter_timeseries[-erad_quarter_time_step+1])
    IBM_all_quarters <- melt(IBM_list, id.vars = colnames(IBM_output[[1]]$quarter_timeseries[[1]]))
    colnames(IBM_all_quarters)[6] <- "Quarter"
    for(snake in 1:nrow(IBM_all_quarters)) {
      IBM_all_quarters$size_category[snake] <- size_class_fun(IBM_all_quarters$SVL[snake])
    }
    IBM_all_quarters$size_category <- factor(IBM_all_quarters$size_category,
                                             levels = size_class_names)
  }
  return(IBM_all_quarters)
}

# Put effort records back together
recreate_IBM_effort_list <- function(strategy_name,
                                     P,
                                     D,
                                     permutation_name,
                                     variant) {
  #all_obs_quarters) {
  IBM_effort_list <- list()
  for(y in 1:total_time_steps[variant]) {
    IBM_output <- readRDS(paste0(results_folders[[strategy_name]][[P]][[D]]["IBM"], "/variant_",
                                 variant, "/IBM_results_start_pop_", permutation_name, "_variant-",
                                 variant,"_IBM_set-",
                                 y, ".rds"))
    IBM_effort_list[[y]] <- IBM_output[[y]]$all_effort
    # Check if there are fewer quarters in the last set than usual 
    if(y == 10) {
      IBM_effort_list[[y]] <- Filter(length, IBM_effort_list[[y]])
      set_quarters <- length(IBM_effort_list[[y]])
    } else {
      set_quarters <- erad_quarter_time_step
    }
    # Correcting quarters before they can be combined
    for(quarter in 1:set_quarters) {
      IBM_effort_list[[y]][[quarter]]$quarter <- IBM_effort_list[[y]][[quarter]]$quarter + 4*(y-1)
    }
  }
  # Simplify the list
  IBM_effort <- unlist(IBM_effort_list, recursive = FALSE)
  return(IBM_effort)
}

# Test


# Function to recreate observation quarters based on condition and set
obs_quarters_fun <- function(variant,
                             variant_effort_records) {
  obs_quarters_by_set <- list()
  for(set in 1:total_time_steps[variant]) {
    raw_obs_quarters <- unlist(method_options[[variant_effort_records[[variant]][set]]]$erad_quarters[erad_methods[c(2,3)]])
    set_quarters <- quarter_name_fun(erad_quarter_time_step,
                                     set_num = set)
    obs_quarters_by_set[[set]] <- set_quarters[raw_obs_quarters]
  }
  return(obs_quarters_by_set)
}


## Function to format each variant's effort data for plotting
format_effort_fun <- function(effort_results){
  # Checking for quarters when no effort occurred, and selecting those out
  effort <- effort_results[unlist(lapply(effort_results, length) != 0)]
  # Binding quarter list into one dataframe
  all_effort <- bind_rows(effort)
  # Adding years
  all_effort$year <- ceiling(all_effort$quarter/4)
  # Calculate quarters in each year
  all_effort$quarter_per_year <- 4 - all_effort$year*4 + all_effort$quarter
  # Renaming quarter & year values
  # all_effort$quarter <- paste0("Quarter ", all_effort$quarter)
  # all_effort$quarter <- factor(all_effort$quarter, 
  #                              levels = paste0("Quarter ", c(1:40)))
  all_effort$year <- paste0("Year ", all_effort$year)
  all_effort$year <- factor(all_effort$year, 
                            levels = paste0("Year ", c(1:10)))
  
  # Adding weeks
  all_effort$week <- ceiling(all_effort$day/7)
  
  return(all_effort)
}

# test_effort <- format_effort_fun(IBM_effort_results)

# Function to rename quarters based on how frequently estimation happens
quarter_name_fun <- function(quarter_time_step,
                             set_num){
  y <- quarter_time_step-1
  quarter_name <- vector()
  for(Y in 1:y) {
    quarter_name[Y] <- set_num*quarter_time_step - Y
  }
  quarter_name <- rev(quarter_name)
  quarter_name[quarter_time_step] <- set_num*quarter_time_step
  return(quarter_name)
}

# # Test
# quarters <- quarter_name_fun(quarter_time_step = 4,set_num = 3)



# 
# ggplot(test_effort, aes(fill = method, x = week)) +
#   geom_bar(stat = "count", position = "dodge") +
#   facet_grid(quarter_per_year ~ year) +
#   scale_fill_hue(labels = plot_labels$method) +
#   labs(y = "Number of days", fill = "Method", x = "Week") +
#   theme_bw()



## Function to calculate the probability of eradication before 10 years are up for a single variant
erad_prob_fun <- function(total_quarters,
                          final_set,
                          n = num_variants) {
  erad_counter <- 0
  for(i in 1:length(total_quarters)) {
    if(total_quarters[i] < (final_set*erad_quarter_time_step)) {
      erad_counter <- erad_counter + 1
    } 
  }
  erad_prob <- erad_counter/n
  erad_quarter <- total_quarters[total_quarters < 40]
  return(list(erad_prob = erad_prob,
              erad_quarter = erad_quarter))
}

# # Test
# erad_prob <- erad_prob_fun(total_quarters = total_quarters,
#                            final_time_step = final_time_step)



# Function to isolate incomplete effort records (if the model run was interrupted for whatever reason)
find_incomplete_effort <- function(variant_effort_records) {
  # Check if there any NAs in the effort condition record, or if the element is empty
  if(TRUE %in% is.na(variant_effort_records)) {
    record_status <- "incomplete"
  } else if (length(variant_effort_records) == 0) {
    record_status <- "incomplete"
  } else {
    record_status <- "complete"
  }
  # }
  # # Deleting null list elements from all of the above 
  # contains_na <- contains_na[!is.na(contains_na)]
  # empty_list <- empty_list[!is.na(empty_list)]
  # incomplete_effort_ind <- sort(c(incomplete_effort_ind, contains_na, empty_list))
  
  return(record_status)
}

# # Test
# effort_record_status <- find_incomplete_effort(variant_effort_records)

# Check if an effort condition record is complete if it does exist


# Function to recreate the condition effort record for instances when it was interrupted and didn't save properly
recreate_effort_conditions <- function(variant) {
  final_set <- total_time_steps[variant]
  #incomplete_variant_conditions <- list()
  #y <- total_time_steps[variant]
  IBM_output <- readRDS(paste0(results_folders[[strategy_name]][[P]][[D]]["IBM"], "/variant_",
                               variant, "/IBM_results_start_pop_", permutation_name, "_variant-",
                               variant,"_IBM_set-", final_set, ".rds"))
  
  all_IBM_effort <- list()
  condition <- vector()
  #t <- list()
  for(set in 1:final_set) {
    all_IBM_effort[[set]] <- IBM_output[[set]]$all_effort
    first_quarter_methods <- list()
    # Extract methods used in each quarter:
    methods <- unique(all_IBM_effort[[set]][[1]]$method)
    # Extract days for each method used in the first quarter:
    for(method in methods) {
      first_quarter_methods[[method]] <- all_IBM_effort[[set]][[1]]$day[all_IBM_effort[[set]][[1]]$method == method]
    } 
    # Check how many possible conditions have the same number of methods as 
    possible_method_options <- condition_num_methods[which(condition_num_methods == length(first_quarter_methods))]
    #print(possible_method_options)
    # If there's only one possible condition, then no need to check if the method days match
    if(length(possible_method_options) == 1) {
      condition[set] <- names(possible_method_options)
    } # For strategies where more than one condition has the same number of methods, I'll need to match method days
    #t[[set]] <- first_quarter_methods
  }
  #incomplete_variant_conditions[[variant]] <- condition 
  
  # incomplete_variant_conditions <- Filter(Negate(is.null), incomplete_variant_conditions)
  #names(incomplete_variant_conditions) <- paste0("variant_", variant)
  
  return(condition)
}

# # Test
# test_condition <- recreate_effort_conditions(variant = incomplete_effort_ind[1])


## Function to plot the method conditions across variants

condition_record_plot_fun <- function(variant_effort_records,
                                      n = num_variants,
                                      strategy = strategy_name) {
  # Start the dataframe by combining the condition vector with the set number for variant 1
  all_effort_df <- cbind(melt(variant_effort_records[[1]]), c(1:total_time_steps[1]))
  colnames(all_effort_df) <- c("condition", "set")
  all_effort_df$variant <- "variant_1"
  # Add each variant's condition list and corresponding set to the dataframe
  for(variant in 2:n) {
    effort_df <- cbind(melt(variant_effort_records[[variant]]), c(1:total_time_steps[variant]))
    colnames(effort_df) <- c("condition", "set")
    effort_df$variant <- paste0("variant_", variant)
    all_effort_df <- bind_rows(all_effort_df, effort_df)
  }
  # Set up y axis based on strategy
  if(strategy == "Strategy_two") {
    y_breaks <- seq(2, 20, 2)
  } else {
    y_breaks <- c(1:10)
  }
  # Plot condition record for all variants
  condition_record_plot <- ggplot(all_effort_df) +
    geom_bar(aes(y = set, fill = condition)) +
    scale_y_continuous(name = "Year", breaks = y_breaks, labels = c(1:10)) +
    scale_x_continuous(breaks = seq(0, 50, 10), labels = seq(0, 50, 10)/50) +
    labs(x = "Proportion of variants", fill = "Condition") +
    theme_bw() 
  return(condition_record_plot)
}

# # Test
# condition_plot <- condition_record_plot_fun(variant_effort_records)



## Function to plot data for all variants in a strategy, for both N and density
variant_data_plot_fun <- function(variant_data,
                                  type_of_y = c("N", "density"),
                                  n = num_variants) {
  # Reformat variant data to plot together
  variant_data_melted <- list()
  for(variant in 1:n) {
    # Melt each variant's data into a long format
    variant_data_melted[[variant]] <- melt(variant_data[[variant]][[type_of_y]], id.vars = "Quarter")
    colnames(variant_data_melted[[variant]])[c(2:3)] <- c("size_class", type_of_y)
  }
  # Melt variants into 1 df
  real_data_all_variants <- melt(variant_data_melted, 
                                 id.vars = colnames(variant_data_melted[[1]]))
  colnames(real_data_all_variants)[4] <- "variant"
  # Y axis label depends on the type_of_y
  if(type_of_y == "N") {
    y_axis_label <- "Total Population"
  } else if (type_of_y == "density") {
    y_axis_label <- "Density per ha"
  }
  
  real_data_plot <- ggplot(real_data_all_variants, aes(x = Quarter, y = .data[[type_of_y]])) +
    geom_path(aes(group = variant)) +
    geom_point(show.legend = FALSE) +
    facet_grid(cols = vars(size_class), scales = "free_y",
               labeller = labeller(size_class = plot_labels$type_of_N)) +
    # scale_color_manual(values = colors$alt,
    #                    labels = plot_labels$alternative, 
    #                    guide = guide_legend(position = "bottom")) +
    labs(y = y_axis_label, color = "") +
    theme_bw() +
    theme(legend.position = "top")
  
  return(list(data_plot = real_data_plot,
              data_all_variants = real_data_all_variants))
}

# # Test
# test_data_plot <- variant_data_plot_fun(variant_data, type_of_y = "density")



## Function to plot estimates vs data for all variants in a strategy, using the final estimation only
results_vs_data_plot_fun <- function(variant_estimates,
                                     melted_variant_data,
                                     n = num_variants) {
  variant_N_estimates <- list()
  for(variant in 1:n) {
    variant_N_estimates[[variant]] <- variant_estimates[[variant]]$N[total_time_steps[variant]][[1]][,-6]
  }
  results_all_variants <- melt(variant_N_estimates, id.vars = colnames(variant_N_estimates[[1]]))
  colnames(results_all_variants)[6] <- "variant"
  
  
  # Adding source to differentiate between data and results
  melted_variant_data$source <- "simulated"
  melted_variant_data$source <- "estimated"
  
  
  all_N_est_vs_data_plot <- ggplot() +
    geom_path(data = results_all_variants, aes(group = variant, x = Quarter, y = N), color = "red") +
    geom_ribbon(data = results_all_variants, aes(y = N, x = Quarter, ymin = lower_CI, ymax = upper_CI, group = variant), fill = "red", alpha = 0.02, linetype = 0) +
    geom_path(data = melted_variant_data, aes(x = Quarter, y = N, group = variant, color = "simulated \n data"), color = "black") +
    facet_grid(vars(size_class), scales = "free_y",
               labeller = labeller(size_class = plot_labels$type_of_N)) +
    theme_bw() +
    # scale_color_manual(name = "", values = c("simulated \n data" = "black",
    #                                          "alt_1" = hue_pal()(2)[1],
    #                                          "alt_2" = hue_pal()(2)[2]),
    #                    labels = c("simulated \n data" = "Simulated \n data",
    #                               "alt_1" = "Alternative 1",
    #                               "alt_2" = "Alternative 2")) +
    guides(color = "none", fill = "none") +
    labs(y = "Population")
  
  return(list(melted_variant_results = results_all_variants,
              results_v_data_plot = all_N_est_vs_data_plot))
}

# # Test
# test_results_v_data <- results_vs_data_plot_fun(variant_estimates,
#                                                 test_data_plot$data_all_variants)

## Function to convert estimation results into a simple dataframe, with mean and
## lower and upper CI estimates for N
estimate_N_summary <- function(jags_output) {
  # # Final primary period and secondary sampling instance in each primary period
  # final_day <- length(observed_erad_days)
  final_quarter <- dim(jags_output$mean$N)[3]
  final_day <- dim(jags_output$mean$N)[2]
  
  mean_N <- as.data.frame(matrix(0, nrow = final_quarter, ncol = 4))
  colnames(mean_N) <- size_class_names
  for(quarter in 1:final_quarter) {
    for(size in 1:length(size_class_names)) {
      mean_N[quarter, size] <- jags_output$mean$N[size,final_day,quarter]
    }
    mean_N$total[quarter] <- jags_output$mean$N.sum[quarter]
  }
  mean_N$Quarter <- c(1:final_quarter)
  
  # Making data long for plotting
  mean_N_long <- melt(mean_N, id.vars = "Quarter")
  colnames(mean_N_long)[2:3] <- c("size_class", "N")
  # Adding lower and upper bounds for 95% confidence interval
  for(quarter in 1:final_quarter) {
    for(size in 1:length(size_class_names)) {
      row_counter <- quarter + (size-1)*final_quarter
      mean_N_long$lower_bound[row_counter] <- jags_output$q2.5$N[size,final_day,quarter]
      mean_N_long$upper_bound[row_counter] <- jags_output$q97.5$N[size,final_day,quarter]
    }
    mean_N_long$lower_bound[quarter+row_counter] <- jags_output$q2.5$N.sum[quarter]
    mean_N_long$upper_bound[quarter+row_counter] <- jags_output$q97.5$N.sum[quarter]
  }
  
  return(mean_N_df = mean_N_long)
}


## Function to check if a variant meets the total suppression goal (reaching and 
## maintaining less than 1 snake/ha)

total_suppression_obj_fun <- function(variant_data,
                                      total_pop_obj = 1) {
  total_suppressed <- vector()
  suppression_reached <- vector()
  for(variant in 1:length(variant_data)) {
    total_density <- variant_data[[variant]]$density$total
    below_total_obj <- vector()
    for(quarter in 1:length(total_density)) {
      if(total_density[quarter] <= total_pop_obj) {
        below_total_obj[quarter] <- 1
      } else {
        below_total_obj[quarter] <- 0
      }
    }
    if(sum(below_total_obj) > 0) {
      # Record the quarter when the total population first goes below the objective
      suppression_reached[variant] <- min(which(below_total_obj == 1))
      # Check to see if it is below for at least the final year, and count it as 
      # successfully suppressed if so
      last_year <- tail(below_total_obj, 4)
      if(sum(last_year) == 4) {
        total_suppressed[variant] <- 1
      } else {
        total_suppressed[variant] <- 0
      }
    } else {
      total_suppressed[variant] <- 0
      suppression_reached[variant] <- NA
    }
  }
  variants_suppressed <- which(total_suppressed == 1)
  suppression_prob <- sum(total_suppressed)/length(variant_data)
  
  return(list(suppession_prob = suppression_prob,
              suppression_maintained = variants_suppressed,
              suppression_reached = suppression_reached))
}

# # Test
# total_suppression_obj_fun(variant_data)


## Function to evaluate the probability of reaching the upper 3 size classes 
## suppression objective (< 1 snake per ha)

upper_3_suppression_obj_fun <- function(variant_data,
                                        upper_3_obj = 1) {
  upper_3_suppressed <- vector()
  suppression_reached <- vector()
  for(variant in 1:length(variant_data)) {
    upper_3_combined <- variant_data[[variant]]$density$medium + variant_data[[variant]]$density$large + variant_data[[variant]]$density$xlarge
    below_total_obj <- vector()
    for(quarter in 1:length(upper_3_combined)) {
      if(upper_3_combined[quarter] <= upper_3_obj) {
        below_total_obj[quarter] <- 1
      } else {
        below_total_obj[quarter] <- 0
      }
    }
    if(sum(below_total_obj) > 0) {
      # Record the quarter when the total population first goes below the objective
      suppression_reached[variant] <- min(which(below_total_obj == 1))
      # Check to see if it is below for at least the final year, and count it as 
      # successfully suppressed if so
      last_year <- tail(below_total_obj, 4)
      if(sum(last_year) == 4) {
        upper_3_suppressed[variant] <- 1
      } else {
        upper_3_suppressed[variant] <- 0
      }
    } else {
      upper_3_suppressed[variant] <- 0
      suppression_reached[variant] <- NA
    }
  }
  variants_suppressed <- which(upper_3_suppressed == 1)
  suppression_prob <- sum(upper_3_suppressed)/length(variant_data)
  
  return(list(suppession_prob = suppression_prob,
              suppression_maintained = variants_suppressed,
              suppression_reached = suppression_reached))
}

# Test
# upper_3_suppression_obj_fun(variant_data)


# Function to separate (and sum when appropriate) estimated encounter probability values - mean and credible intervals for 95th percentile
encounter_prob_results_fun <- function(output_jags,
                                       obs_quarters) {
  # # Separating relevant raw results - mean, sd and 95% CI of N
  # N_jags_output <- output_jags$mean$N
  # sd_jags_output <- output_jags$sd$N
  # estimated_N_95_CI_lower <- output_jags$q2.5$N
  # estimated_N_95_CI_upper <- output_jags$q97.5$N
  
  # Setting up dataframe to hold estimated mean N values for each quarter and size
  estimated_N <- as.data.frame(matrix(NA, ncol = (length(size_class_names)+2), 
                                      nrow = length(obs_quarters)))
  colnames(estimated_N) <- c("Quarter", size_class_names, "total")
  estimated_N$Quarter <- c(obs_quarters)
  # Setting up dataframe to hold SD, lower and upper CI values for estimated N 
  # for each quarter and size
  estimated_N_sd <- estimated_N
  estimated_N_lower_CI <- estimated_N
  estimated_N_upper_CI <- estimated_N
  
  for(quarter in 1:length(obs_quarters)) {
    for(size in 1:length(size_class_names)) {
      estimated_N[quarter, size+1] <-  output_jags$mean$N[size, 1, quarter]
      estimated_N_sd[quarter, size+1] <-  output_jags$sd$N[size, 1, quarter]
      estimated_N_lower_CI[quarter, size+1] <- output_jags$q2.5$N[size, 1, quarter]
      estimated_N_upper_CI[quarter, size+1] <- output_jags$q97.5$N[size, 1, quarter]
      
    }
    estimated_N$total[quarter] <- output_jags$mean$N.sum[quarter]
    estimated_N_sd$total[quarter] <- output_jags$sd$N.sum[quarter]
    estimated_N_lower_CI$total[quarter] <- output_jags$q2.5$N.sum[quarter]
    estimated_N_upper_CI$total[quarter] <- output_jags$q97.5$N.sum[quarter]
    
  }
  
  # Combining the above dataframes to produce a single dataframe of the estimated
  # Ns, with lower and upper CIs
  all_estimated_N <- melt(estimated_N, id.vars = "Quarter")
  colnames(all_estimated_N)[2:3] <- c("size_class", "N")
  lower_CI <- melt(estimated_N_lower_CI, id.vars = "Quarter")
  upper_CI <- melt(estimated_N_upper_CI, id.vars = "Quarter")
  sd <- melt(estimated_N_sd, id.vars = "Quarter")
  all_estimated_N$lower_CI <- lower_CI$value
  all_estimated_N$upper_CI <- upper_CI$value
  all_estimated_N$sd <- sd$value
  
  
  return(list(estimated_N = estimated_N,
              estimated_N_sd = estimated_N_sd,
              estimated_N_95_CI_lower = estimated_N_lower_CI,
              estimated_N_95_CI_upper = estimated_N_upper_CI,
              all_estimated_N = all_estimated_N))
}


# # Test
# summed_results <- summed_est_results_fun(output_jags, obs_quarters)


## Function to separate (and sum when appropriate) estimated encounter probability values
## - mean and credible intervals for 95th percentile, based on the 
encounter_prob_results_fun <- function(output_jags,
                                       IBM_effort) {
  # Construct quarter vector
  final_quarter <- dim(output_jags$mean$p)[3]
  obs_quarters <- c(1:final_quarter)
  # Truncate IBM_effort based on the quarters in jags output 
  IBM_effort <- IBM_effort[obs_quarters]
  # Figure out how many methods were used in each quarter, and create list of method days
  quarter_methods <- list()
  first_method_day <- list()
  for(quarter in 1:final_quarter) {
    quarter_effort <- IBM_effort[[quarter]][IBM_effort[[quarter]]$method %in% erad_methods[c(2,3)],]
    quarter_methods[[quarter]] <-  unique(quarter_effort$method)
    num_days <- length(unique(quarter_effort$day))
    num_methods <- length(quarter_methods[[quarter]])
    first_method_day[[quarter]] <- list() 
    for(method in 1:num_methods) {
      first_method_day[[quarter]][[method]] <- method
    }
    names(first_method_day[[quarter]]) <- quarter_methods[[quarter]]
  }
  
  # Create lists to hold dfs for mean, lower CI and upper CI encounter values
  estimated_encounter <- list()
  encounter_lower_CI <- list()
  encounter_upper_CI <- list()
  
  for(quarter in 1:final_quarter) {
    estimated_encounter[[quarter]] <- as.data.frame(matrix(NA, ncol = (length(size_class_names)+1), 
                                                           nrow = length(first_method_day[[quarter]])))
    colnames(estimated_encounter[[quarter]]) <- c("method", size_class_names)
    encounter_lower_CI[[quarter]] <- estimated_encounter[[quarter]]
    encounter_upper_CI[[quarter]] <- estimated_encounter[[quarter]]
    for(method in 1:length(first_method_day[[quarter]])) {
      estimated_encounter[[quarter]][method, 1] <- names(first_method_day[[quarter]])[method]
      encounter_lower_CI[[quarter]][method, 1] <- names(first_method_day[[quarter]])[method]
      encounter_upper_CI[[quarter]][method, 1] <- names(first_method_day[[quarter]])[method]
      for(size in 1:length(size_class_names)) {
        estimated_encounter[[quarter]][method, size+1] <-  output_jags$mean$p[size, first_method_day[[quarter]][[method]], quarter]
        encounter_lower_CI[[quarter]][method, size+1] <- output_jags$q2.5$p[size, first_method_day[[quarter]][[method]], quarter]
        encounter_upper_CI[[quarter]][method, size+1] <- output_jags$q97.5$p[size, first_method_day[[quarter]][[method]], quarter]
        
      }
    }
  }
  
  # Combining the above dataframes to produce a single dataframe of the estimated
  # encounter probabilities, with associated lower and upper CIs
  all_estimated_encounter <- melt(estimated_encounter, id.vars = "method")
  colnames(all_estimated_encounter)[2:4] <- c("size_class", "encounter_prob", "Quarter")
  lower_CI <- melt(encounter_lower_CI, id.vars = "method")
  upper_CI <- melt(encounter_upper_CI, id.vars = "method")
  all_estimated_encounter$lower_CI <- lower_CI$value
  all_estimated_encounter$upper_CI <- upper_CI$value
  
  
  return(all_estimated_encounter = all_estimated_encounter)
}

# # Test
# test_encounter <- encounter_prob_results_fun(output_jags,
#                                              IBM_effort_list)



## Function to create plot of encounter probabilities for all variants
encounter_prob_plot <- function(variant_estimates,
                                final_set = total_time_steps,
                                n = num_variants) {
  variant_encounter_estimates <- list()
  for(variant in 1:n) {
    variant_encounter_estimates[[variant]] <- variant_estimates[[variant]]$encounter[final_set[variant]][[1]]
  }
  results_all_variants <- melt(variant_encounter_estimates, id.vars = colnames(variant_encounter_estimates[[1]]))
  colnames(results_all_variants)[7] <- "variant"
  
  
  all_encounter_prob_plot <- ggplot(results_all_variants) +
    geom_path(aes(group = variant, x = Quarter, y = encounter_prob), color = "red") +
    geom_ribbon(aes(y = encounter_prob, x = Quarter, ymin = lower_CI, ymax = upper_CI, group = variant), fill = "red", alpha = 0.02, linetype = 0) +
    facet_grid(size_class ~ method, scales = "free_y",
               labeller = labeller(size_class = plot_labels$size_class)) +
    theme_bw() +
    guides(color = "none", fill = "none") +
    labs(y = "Encounter probability")
  
}

