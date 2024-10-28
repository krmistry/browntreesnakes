### Results processing
library(reshape2)
library(here)


source(here("Scripts/all_strategy_set_up.R"))
source(here("Scripts/model_evaluation.R"))
source(here("Scripts/results_processing_functions.R"))

# Plot labels 
plot_labels <- list()
plot_labels$size_class <- c("small" = "Small",
                            "medium" = "Medium",
                            "large" = "Large",
                            "xlarge" = "X-large")
plot_labels$method <- c("ADS" = "ADS",
                        "visual" = "Visual survey",
                        "trap" = "Live traps",
                        "bait_tube" = "Bait tubes")
plot_labels$type_of_N <- c("small" = "Small", 
                           "medium" = "Medium", 
                           "large" = "Large", 
                           "xlarge" = "X-large", 
                           "total" = "Total")


# Checking for and creating if necessary all results folders
strategies <- paste0("Strategy_", c("one", "two", "three", "four"))
types_of_results <- c("IBM", "Estimation", "Processed_results", "")
save_folder <- paste0(here("Results", "alt_strategies"), "/")

# Creating directory of folders (actually create the folders if they don't exist, 
# otherwise just creating an indexed list of the names for file reading and saving)
results_folders <- list()
for(strategy in strategies) {
  results_folders[[strategy]] <- list()
  for(P in names(starting_pop)) {
    results_folders[[strategy]][[P]] <- list()
    for(D in names(starting_size_dist)) {
      results_folders[[strategy]][[P]][[D]] <- vector()
      permutation_name <- paste0(P, "_", D)
      for(result in types_of_results) {
        folder_name <- paste0(save_folder, strategy, "/", permutation_name, "/", result)
        if(!dir.exists(folder_name)) {
          dir.create(folder_name, recursive = TRUE)
        }
        results_folders[[strategy]][[P]][[D]][result] <- folder_name
      }
    }
  }
}


# Process results for a strategy's single permutation 
P <- names(starting_pop)[2]
D <- names(starting_size_dist)[1]
strategy_name <- strategies[4]
permutation_name <- paste0(P, "_", D)
final_time_step <- 10
num_variants <- 50
source(here("strategy_4_set_up.R"))

# Isolating the methods & method days in the first quarter of each condition for the strategy
first_quarter_method_options <- list()
condition_num_methods <- vector()
for(condition in names(method_options)) {
  first_quarter_method_options[[condition]] <- method_options[[condition]]$erad_days$quarter_1
  condition_num_methods[condition] <- length(first_quarter_method_options[[condition]])
}


issues_with_processing <- c()
runs_to_process <- setdiff(c(1:50), issues_with_processing)


total_time_steps <- vector()
total_quarters <- vector()

variant_effort_records <- list()
variant_data <- list()
variant_estimates <- list()
for(variant in runs_to_process) {
  variant_data[[variant]] <- list()
  variant_estimates[[variant]] <- vector(mode = "list", length = 2)
  variant_estimates[[variant]]$N <- list()
  variant_estimates[[variant]]$density <- list()
  variant_estimates[[variant]]$encounter <- list()
  model_results <- vector(mode = "list", length = 2)
  model_results$model_metrics <- list()
  
  # Check how many time steps the variant went through (to see if population was eradicated)
  all_IBM_files <- list.files(paste0(results_folders[[strategy_name]][[P]][[D]]["IBM"], 
                                     "/variant_", variant))
  set_IBM_files <- all_IBM_files[grep(paste0("IBM_results_"), all_IBM_files)]
  total_time_steps[variant] <- length(set_IBM_files)
  # Create all quarters from the original IBM outputs (since the all_quarters that 
  # I saved originally almost all had the wrong quarters assigned - future version should just 
  # be able to read in the IBM_all_quarters file that is saved for each variant)
  IBM_all_quarters <- recreate_IBM_all_quarters(strategy_name,
                                                P,
                                                D,
                                                permutation_name,
                                                variant)
  
  # If the effort_record exists, record that, otherwise indicate that it needs to be reconstructs
  effort_record_file <- paste0(results_folders[[strategy_name]][[P]][[D]]["IBM"],
                               "/variant_", variant, "/effort_record_start_pop_", 
                               permutation_name, "_variant-",
                               variant, ".rds")
  if(file.exists(effort_record_file) == TRUE) {
    effort_record <- readRDS(effort_record_file)
    variant_effort_records[[variant]] <- effort_record$condition
  } else {
    variant_effort_records[[variant]] <- recreate_effort_conditions(variant)
  }
  # Check if the effort_record is incomplete, and if so then recreate it
  effort_record_status <- find_incomplete_effort(variant_effort_records[[variant]])
  if(effort_record_status == "incomplete") {
    variant_effort_records[[variant]] <- recreate_effort_conditions(variant)
  }
  
  # Create vector of all observed quarters
  all_obs_quarters <- unique(IBM_all_quarters$Quarter)
  # Record the maximum quarter reached, for eradication probability later
  total_quarters[variant] <- max(all_obs_quarters)
  
  # Construct effort quarterly records
  IBM_effort_results <- recreate_IBM_effort_list(strategy_name,
                                                 P,
                                                 D,
                                                 permutation_name,
                                                 variant)
  
  for(set in c(1:5, total_time_steps[variant])) {
    
    # Isolate observed quarters based on the current set
    obs_quarters <- all_obs_quarters[c(1:(set*erad_quarter_time_step))]
    obs_quarters <- obs_quarters[!is.na(obs_quarters)]
    
    
    
    # Read in jags results
    jags_output <- readRDS(paste0(results_folders[[strategy_name]][[P]][[D]]["Estimation"],
                                  "/output_jags_start_pop_", permutation_name,
                                  "_variant-",variant, "_est_", set, ".RDS"))
    # Summarize encounter probability for this set
    variant_estimates[[variant]]$encounter[[set]] <- encounter_prob_results_fun(output_jags = jags_output,
                                                                                IBM_effort = IBM_effort_results)
    
    # Calculate model summary & metrics
    model_results$model_metrics[[set]] <- eval_metrics_fun(simulation_quarter_data = IBM_all_quarters,
                                                           output_jags = jags_output,
                                                           obs_quarters = obs_quarters)
    
    
    if(set == total_time_steps[variant]) {
      # Plotting simulated data only (at first)
      final_erad_plot <- ggplot(IBM_all_quarters,
                                aes(x = Quarter, fill = size_category)) +
        geom_bar() +
        #geom_hline(yintercept = K) +
        theme_bw() +
        scale_x_continuous(breaks = unique(IBM_all_quarters$Quarter),
                           labels = unique(IBM_all_quarters$Quarter))
      # Plotting effort through time
      all_effort <- format_effort_fun(IBM_effort_results)
      effort_plot <- ggplot(all_effort, aes(fill = method, x = week)) +
        geom_bar(stat = "count", position = "dodge") +
        facet_grid(quarter_per_year ~ year) +
        scale_fill_hue(labels = plot_labels$method) +
        labs(y = "Number of days", fill = "Method", x = "Week") +
        theme_bw()
      
      # Saving all time steps of IBM data (only in the last step, as all quarters will be there)
      variant_data[[variant]]$N <- model_results$model_metrics[[set]]$summed_data
      variant_data[[variant]]$density <- model_results$model_metrics[[set]]$summed_data/area_size
      variant_data[[variant]]$density$Quarter <- variant_data[[variant]]$N$Quarter
      
      
    }
    # Saving estimates from each saved estimation sets (raw N, and then density)
    variant_estimates[[variant]]$N[[set]] <- model_results$model_metrics[[set]]$summed_results$all_estimated_N
    variant_estimates[[variant]]$density[[set]] <- model_results$model_metrics[[set]]$summed_results$all_estimated_N
    variant_estimates[[variant]]$density[[set]][, c(2:6)] <- variant_estimates[[variant]]$density[[set]][, c(2:6)]/area_size
    
    
  }
  variant_results <- list(IBM_all_quarters = IBM_all_quarters,
                          effort_data = all_effort, 
                          variant_data_plot = final_erad_plot,
                          variant_effort_plot = effort_plot,
                          variant_metrics = model_results$model_metrics,
                          variant_data = variant_data[[variant]],
                          variant_estimates = variant_estimates[[variant]])
  
  # Save metrics
  saveRDS(variant_results, paste0(results_folders[[strategy_name]][[P]][[D]]["Processed_results"],
                                  "/variant-", variant, "_results.RDS"))
  print(paste0("variant ", variant, " processed"))
}


permutation_results <- list()
# Eradication probability of this permutation
# Calculate the probability if eradication for this variant
permutation_results$erad_prob <- erad_prob_fun(total_quarters = total_quarters,
                                               final_set = final_time_step)
# Calculate the probability of reaching and maintaining total pop suppression goal (1 snake/ha)
permutation_results$total_suppress_prob <- total_suppression_obj_fun(variant_data)
# Calculate the probability of reaching and maintaining upper 3 size class pop suppression goal (1 snake/ha)
permutation_results$upper_3_suppress_prob <- upper_3_suppression_obj_fun(variant_data)
# Plotting data for all variants, for both N and density
permutation_results$N_data_plot <- variant_data_plot_fun(variant_data,
                                                         type_of_y = "N")
permutation_results$density_data_plot <- variant_data_plot_fun(variant_data, 
                                                               type_of_y = "density")
# Plotting results vs data for all variants
permutation_results$results_v_data_plot <- results_vs_data_plot_fun(variant_estimates,
                                                                    permutation_results$N_data_plot$data_all_variants)
# Plotting encounter probability from the final estimate for all variants
permutation_results$encounter_prob_plot <- encounter_prob_results_fun(variant_estimates)
# Plotting method conditions across variants
permutation_results$condition_plot <- condition_record_plot_fun(variant_effort_records)

# Save plots & summarized data sets
saveRDS(permutation_results, paste0(results_folders[[strategy_name]][[P]][[D]]["Processed_results"],
                                    "/permutation-", permutation_name, "_results.RDS"))




variant_plots <- list()

for(variant in 1:10) {
  results <- readRDS(paste0(results_folders[[strategy_name]][[P]][[D]]["Processed_results"], 
                            "/variant-", variant, "_results.RDS"))
  variant_plots[[variant]] <- results$variant_plots$final_erad_plot
  
  
}



# Checking if the estimation triggering a threshold was accurate
estimate_wrong_ind <- list()
num_estimate_wrong <- vector()
condition_per_data <- list()
for(variant in 1:50) {
  condition_per_estimate <- variant_effort_records[[variant]]
  var_data <- variant_data_melt[[variant]]
  colnames(var_data)[c(2:3)] <- c("size_class", "N")
  condition_per_data[[variant]] <- vector()
  for(set in 1:total_time_steps[variant]) {
    condition_per_data[[variant]][set] <- strat_3_threshold_fun(mean_N_df = var_data[var_data$Quarter %in% c(1:erad_quarter_time_step*set),])
  }
  estimate_wrong_ind[[variant]] <- which(condition_per_estimate != condition_per_data[[variant]])
  num_estimate_wrong[variant] <- length(estimate_wrong_ind[[variant]])
}





# Reconstructing incomplete condition records for strategy two
effort_record_status <- vector()
for(variant in 1:num_variants) {
  effort_record_status[variant] <- find_incomplete_effort(variant_effort_records[[variant]])
}
incomplete_ind <- which(effort_record_status == "incomplete")

#incomplete_ind <- which(variant_effort_records == "incomplete effort record")

partial_effort_record <- list()
for(variant in incomplete_ind) {
  partial_effort_record[[variant]] <- vector()
  effort_record_file <- paste0(results_folders[[strategy_name]][[P]][[D]]["IBM"],
                               "/variant_", variant, "/effort_record_start_pop_", 
                               permutation_name, "_variant-", variant, ".rds")
  effort_record <- readRDS(effort_record_file)
  if(length(effort_record$condition) > 0) {
    # Isolate the missing sets
    missing_ind <- which(is.na(effort_record$condition))
  } else {
    missing_ind <- c(1:total_time_steps[variant])
  }
  for(set in missing_ind) {
    # the first condition is always initial, if that's missing
    if(set == 1) {
      partial_effort_record[[variant]][1] <- "initial"
    } else if(set %in% (c(1:5, 10, 15)+1)){ 
      # If the estimation results were saved for the set before the missing one, 
      # then the condition can be recalculated from the estimation results
      jags_output <- readRDS(paste0(results_folders[[strategy_name]][[P]][[D]]["Estimation"],
                                    "/output_jags_start_pop_", permutation_name,
                                    "_variant-",variant, "_est_", (set-1), ".RDS"))
      final_quarter <- dim(jags_output$mean$N)[3]
      estimate_summary <- estimate_N_summary(jags_output = jags_output)
      partial_effort_record[[variant]][set] <- strat_2_threshold_fun(estimate_summary)
    } else {
      # If the estimation wasn't saved, then try recreating with the number of methods (works for some instances)
      recreated_conditions <- recreate_effort_conditions(variant)
      partial_effort_record[[variant]][set] <- recreated_conditions[set]
    }
  }
  effort_record$condition[missing_ind] <- partial_effort_record[[variant]][missing_ind]
  variant_effort_records[[variant]] <- effort_record$condition
  print(paste0("variant ", variant, " completed"))
}



# Separate mean, lower and upper CI values for size and method based encounter probability



