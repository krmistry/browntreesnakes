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
types_of_results <- c("IBM", "Estimation", "Processed_results")
save_folder <- paste0(here("Results/alt_strategies"), "/")


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


strategy_name <- strategies[1]
source(here("Scripts/strategy_1_set_up.R"))
num_variants <- 50
# Permutations that are on this computer (boiga)
permutations <- list(c(P = 1, D = 1),
                     c(P = 1, D = 2))

for(permutation in permutations) {
  # Process results for a strategy's single permutation 
  P <- names(starting_pop)[permutation["P"]]
  D <- names(starting_size_dist)[permutation["D"]]
  permutation_name <- paste0(P, "_", D)
  
  total_quarters <- vector()
  variant_data <- list()
  for(variant in 1:num_variants) {
    variant_data[[variant]] <- list()
    IBM_output <- readRDS(paste0(results_folders[[strategy_name]][[P]][[D]]["IBM"],  
                                 "/IBM_start_pop_", permutation_name, "-var_",
                                 variant, ".rds"))
    if(permutation["P"] == 2 & permutation["D"] == 1) {
      IBM_output <- IBM_output
    } else {
      IBM_output <- IBM_output[[1]]
    }
    # Calculate how many quarters the variant ran (capping it at 10 years, i.e. 40 quarters)
    IBM_all_quarters <- IBM_output$all_quarters
    total_quarters[variant] <- min(max(IBM_all_quarters$Quarter), 40)
    # Plot size class pops through time series
    final_erad_plot <- ggplot(IBM_all_quarters,
                              aes(x = Quarter, fill = size_category)) +
      geom_bar() +
      #geom_hline(yintercept = K) +
      theme_bw() +
      scale_x_continuous(breaks = unique(IBM_all_quarters$Quarter),
                         labels = unique(IBM_all_quarters$Quarter))
    # Plot effort through time series
    IBM_effort <- format_effort_fun(IBM_output$all_effort)
    effort_plot <- ggplot(IBM_effort, aes(fill = method, x = week)) +
      geom_bar(stat = "count", position = "dodge") +
      facet_grid(quarter_per_year ~ year) +
      scale_fill_hue(labels = plot_labels$method) +
      labs(y = "Number of days", fill = "Method", x = "Week") +
      theme_bw()
    # Summarize data 
    summed_data <- summed_sim_data_fun(IBM_all_quarters)
    # Saving all time steps of IBM data (only in the last step, as all quarters will be there)
    variant_data[[variant]]$N <- summed_data
    variant_data[[variant]]$density <- summed_data/area_size
    variant_data[[variant]]$density$Quarter <- variant_data[[variant]]$N$Quarter
    # Save variant data & plots
    variant_results <- list(IBM_all_quarters = IBM_all_quarters,
                            effort_data = IBM_effort, 
                            variant_data_plot = final_erad_plot,
                            variant_effort_plot = effort_plot,
                            summed_data = summed_data)
    saveRDS(variant_results, paste0(results_folders[[strategy_name]][[P]][[D]]["Processed_results"],
                                    "/variant-", variant, "_results.RDS"))
    print(paste0("variant ", variant, " processed"))
  }
  
  permutation_results <- list()
  # Eradication probability of this permutation
  # Calculate the probability if eradication for this variant
  permutation_results$erad_prob <- erad_prob_fun(total_quarters = total_quarters,
                                                 final_set = 1)
  # Calculate the probability of reaching and maintaining total pop suppression goal (1 snake/ha)
  permutation_results$total_suppress_prob <- total_suppression_obj_fun(variant_data)
  # Calculate the probability of reaching and maintaining upper 3 size class pop suppression goal (1 snake/ha)
  permutation_results$upper_3_suppress_prob <- upper_3_suppression_obj_fun(variant_data)
  # Plotting data for all variants, for both N and density
  permutation_results$N_data_plot <- variant_data_plot_fun(variant_data,
                                                           type_of_y = "N")
  permutation_results$density_data_plot <- variant_data_plot_fun(variant_data, 
                                                                 type_of_y = "density")
  # Save plots & summarized data sets
  saveRDS(permutation_results, paste0(results_folders[[strategy_name]][[P]][[D]]["Processed_results"],
                                      "/permutation-", permutation_name, "_results.RDS"))
  
  print(paste0("permutation ", permutation_name, " processed"))
  
}




