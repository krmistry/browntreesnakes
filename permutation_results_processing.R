### Processing & plotting permutation results for each strategy

library(reshape2)
library(here)
library(tidyr)

source(here("Scripts/all_strategy_set_up.R"))
source(here("Scripts/model_evaluation.R"))
source(here("Scripts/results_processing_functions.R"))
source(here("Scripts/alt_strategies_cost_calc.R"))

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
save_folder <- paste0(here("Results", "alt_strategies"), "/")

#strategy_name <- strategies[3]
for(strategy_name in strategies) {
  permutation_N_data_list <- list()
  permutation_estimates_list <- list()
  permutation_conditions_list <- list()
  permutation_objectives_probs <- as.data.frame(matrix(NA, 
                                                          nrow = (length(starting_pop)*length(starting_size_dist)),
                                                          ncol = 8))
  colnames(permutation_objectives_probs) <- c("init_pop", "init_size_dist", "erad_prob", "mean_erad_quarter",
                                              "total_suppress_prob", "mean_total_quarter",  "upper_3_suppress_prob",
                                              "mean_upper_3_quarter")
  row_counter <- 1
  for(P in names(starting_pop)) {
    permutation_N_data_list[[P]] <- list()
    permutation_estimates_list[[P]] <- list()
    permutation_conditions_list[[P]] <- list()
    #permutation_erad_prob_list[[P]] <- list()
    for(D in names(starting_size_dist)) {
      permutation_results <- readRDS(paste0(here("Results", "alt_strategies", strategy_name, "permutation_results"),
                                            "/permutation-", P, "_", D, "_results.RDS"))
      # Assign N data 
      permutation_N_data_list[[P]][[D]] <- permutation_results$N_data_plot$data_all_variants
      permutation_N_data_list[[P]][[D]]$init_size_dist <- D
      permutation_N_data_list[[P]][[D]]$init_pop <- P
      # Assign estimated N data
      permutation_estimates_list[[P]][[D]] <- permutation_results$results_v_data_plot$melted_variant_results
      permutation_estimates_list[[P]][[D]]$init_size_dist <- D
      permutation_estimates_list[[P]][[D]]$init_pop <- P
      # Assign eradication probability 
      #permutation_erad_prob_list[[P]][[D]] <- permutation_results$erad_prob
      
      ##### COME BACK TO THIS - HAVE TO RECALCULATE (quarter when eradication occurs) FROM ORIGINAL DATA #####
      # # Save the probability of meeting each of the 3 objectives
      # permutation_objectives_probs$init_pop[row_counter] <- P
      # permutation_objectives_probs$init_size_dist[row_counter] <- D
      # permutation_objectives_probs$erad_prob[row_counter] <- permutation_results$erad_prob
      # permutation_objectives_probs$total_suppress_prob[row_counter] <- permutation_results$total_suppress_prob$suppession_prob
      # permutation_objectives_probs$upper_3_suppress_prob[row_counter] <- permutation_results$upper_3_suppress_prob$suppession_prob
      # # Calculate the mean quarter when the permutation reached an objective, if it did
      # if(permutation_objectives_probs$erad_prob[row_counter] > 0) {
      #   # First, calculate the number of quarters each variant reached
      #   total_quarters <- vector()
      #   for(variant in 1:num_variants) {
      #     variant_data <- permutation_N_data_list[[P]][[D]][permutation_N_data_list[[P]][[D]]$variant == variant,]
      #     total_quarters[variant] <- max(variant_data$Quarter)
      #   }
      #   # Separate out 
      #   erad_quarter <- total_quarters[total_quarters < 40]
      #   permutation_objectives_probs$mean_erad_quarter[row_counter] <- mean(erad_quarter)
      # } else { # if the permutation had 0 eradication probability, this is NA
      #   permutation_objectives_probs$mean_erad_quarter[row_counter] <- NA
      # }
      # # If the total suppression objective probability is higher than 0, calculate the mean quarter when the population went below the goal
      # if(permutation_objectives_probs$total_suppress_prob[row_counter] > 0) {
      #   maintained_ind <- permutation_results$total_suppress_prob$suppression_maintained
      #   permutation_objectives_probs$mean_total_quarter[row_counter] <- mean(permutation_results$total_suppress_prob$suppression_reached[maintained_ind])
      # } else {
      #   permutation_objectives_probs$mean_total_quarter[row_counter] <- NA
      # }
      # # Same as above for the upper 3 size class objective
      # if(permutation_objectives_probs$upper_3_suppress_prob[row_counter] > 0) {
      #   maintained_ind <- permutation_results$upper_3_suppress_prob$suppression_maintained
      #   permutation_objectives_probs$mean_upper_3_quarter[row_counter] <- mean(permutation_results$upper_3_suppress_prob$suppression_reached[maintained_ind])
      # } else {
      #   permutation_objectives_probs$mean_upper_3_quarter[row_counter] <- NA
      # }
      
      # Processing and saving permutation costs
      permutation_costs <- dynamic_strategy_cost_calc(permutation_results,
                                                      strategy_condition_costs,
                                                      strategy = strategy_name)
      
      # Saving the permtuation's condition list
      permutation_conditions_list[[P]][[D]] <- permutation_results$condition_plot$condition_record
    }
  }
  # Melting all eradication probabilities for the permutations
  # permutation_erad_prob <- melt(permutation_erad_prob_list, id.vars = NULL)
  # colnames(permutation_erad_prob) <- c("erad_prob", "init_size_dist", "init_pop")
  
  # Melting condition lists into a single dataframe
  conditions_list <- melt(permutation_conditions_list, id.vars = colnames(permutation_conditions_list[[1]][[1]]))
  colnames(conditions_list)[c(4:5)] <- c("size_dist", "starting_pop")
  # Combining all N data with appropriate permutation labeling
  data_list <- list()
  for(P in 1:length(starting_pop)) {
    data_list[[P]] <- bind_rows(permutation_N_data_list[[P]])
  }
  all_data <- bind_rows(data_list)
  
  # Creating another column so geom_path will isolate each variant in each permutation
  all_data <- mutate(all_data, variant_pop_size_dist = paste0(init_pop, "_", 
                                                              init_size_dist, "-",
                                                              variant))
  # Manually set levels for initial pop for legend order
  all_data$init_pop <- factor(all_data$init_pop, levels = c("low", "medium", "high"))
  
  strategy_plot <- ggplot(all_data, aes(x = Quarter, y = N, color = init_pop)) +
    geom_path(aes(group = variant_pop_size_dist)) +
    #geom_point(show.legend = FALSE) +
    facet_grid(init_size_dist ~ size_class, scales = "free_y",
               labeller = labeller(size_class = plot_labels$type_of_N)) +
    # scale_color_manual(values = colors$alt,
    #                    labels = plot_labels$alternative, 
    #                    guide = guide_legend(position = "bottom")) +
    labs(y = "Population", color = "", title = strategy_name) +
    theme_bw() +
    theme(legend.position = "top")
  
  ggsave(filename = paste0(save_folder, strategy_name, "/permutation_results/", strategy_name, "_permutations_data_plot.png"), 
         strategy_plot, device = 'png', width = 6, height = 4)
}




### By permutation
permutations <- vector()
counter <- 1
for(p in 1:length(starting_pop)) {
  for(d in 1:length(starting_size_dist)) {
    permutations[counter] <- paste0(names(starting_pop)[p], "_", names(starting_size_dist)[d])
    counter <- counter + 1
  }
}

N_data_list <- list()
estimates_list <- list()
conditions_list <- list()
objectives_probs <- list()
permutation_costs <- list()
for(permutation_name in permutations) {
  N_data_list[[permutation_name]] <- list()
  estimates_list[[permutation_name]] <- list()
  conditions_list[[permutation_name]] <- list()
  permutation_costs[[permutation_name]] <- list()
  objectives_probs[[permutation_name]] <- as.data.frame(matrix(NA,
                                                          nrow = length(strategies),
                                                          ncol = 7))
  colnames(objectives_probs[[permutation_name]]) <- c("strategy", "erad_prob", "mean_erad_quarter",
                                                 "total_suppress_prob", "mean_total_quarter",  "upper_3_suppress_prob",
                                                 "mean_upper_3_quarter")
  row_counter <- 1
  for(strategy in strategies) {
    N_data_list[[permutation_name]][[strategy]] <- list()
    estimates_list[[permutation_name]][[strategy]] <- list()
    conditions_list[[permutation_name]][[strategy]] <- list()
    # Read in results for this permutation and strategy
    permutation_results <- readRDS(paste0(here("Results", "alt_strategies", strategy, "permutation_results"),
                                            "/permutation-", permutation_name, "_results.RDS"))
    # Assign N data 
    N_data_list[[permutation_name]][[strategy]] <- permutation_results$N_data_plot$data_all_variants
    # Assign estimated N data
    estimates_list[[permutation_name]][[strategy]] <- permutation_results$N_data_plot$data_all_variants
    # Assign conditions data (except for strategy 1)
    if(strategy == strategies[1]) {
      conditions_list[[permutation_name]][[strategy]] <- NA
    } else {
      conditions_list[[permutation_name]][[strategy]] <- permutation_results$condition_plot$condition_record
    }
    # Assign eradication & threshold probabilities
    objectives_probs[[permutation_name]]$strategy[row_counter] <- strategy
    objectives_probs[[permutation_name]]$erad_prob[row_counter] <- permutation_results$erad_prob$erad_prob
    objectives_probs[[permutation_name]]$mean_erad_quarter[row_counter] <- mean(permutation_results$erad_prob$erad_quarter)
    objectives_probs[[permutation_name]]$total_suppress_prob[row_counter] <- permutation_results$total_suppress_prob$suppession_prob
    objectives_probs[[permutation_name]]$mean_total_quarter[row_counter] <- mean(permutation_results$total_suppress_prob$suppression_reached[(!is.na(permutation_results$total_suppress_prob$suppression_reached))])
    objectives_probs[[permutation_name]]$upper_3_suppress_prob[row_counter] <- permutation_results$upper_3_suppress_prob$suppession_prob
    objectives_probs[[permutation_name]]$mean_upper_3_quarter[row_counter] <- mean(permutation_results$upper_3_suppress_prob$suppression_reached[(!is.na(permutation_results$upper_3_suppress_prob$suppression_reached))])

    # Processing and saving permutation costs (for dynamic strategies, so all but strategy one)
    if(strategy %in% strategies[c(3:4)]) {
      permutation_costs[[permutation_name]][[strategy]] <- dynamic_strategy_cost_calc(permutation_results,
                                                    strategy_condition_costs,
                                                    strategy = strategy)
    } else if (strategy == strategies[1]) {
      permutation_costs[[permutation_name]][[strategy]] <- strategy_condition_costs[[strategy]]$initial
    } else {
      permutation_costs[[permutation_name]][[strategy]] <- strat_2_dynamic_strategy_cost_calc(permutation_results,
                                                                                              strategy_condition_costs)
    }
      
    row_counter <- row_counter + 1
  }
}

 

# Combining all N data with appropriate labeling
all_data <- melt(N_data_list, id.vars = colnames(N_data_list[[1]][[1]]))
colnames(all_data)[c(5:6)] <- c("Strategy", "permutation")
all_data <- separate(all_data, permutation, into = c("init_pop", "init_size_dist"), 
                     sep = "_", extra = "merge")

# Creating another column so geom_path will isolate each variant in each permutation
all_data <- mutate(all_data, variant_strategy = paste0(Strategy, "_",
                                                            variant))
# Manually set levels for initial pop & strategies for legend order
all_data$init_pop <- factor(all_data$init_pop, levels = c("low", "medium", "high"))
all_data$Strategy <- factor(all_data$Strategy, levels = strategies)

permutations_plot <- ggplot(all_data[all_data$size_class == "total",], aes(x = Quarter, y = N, color = Strategy)) +
  geom_path(aes(group = variant_strategy)) +
  #geom_point(show.legend = FALSE) +
  facet_grid(init_size_dist ~ init_pop, scales = "free_y") +
  # scale_color_manual(values = colors$alt,
  #                    labels = plot_labels$alternative,
  #                    guide = guide_legend(position = "bottom")) +
  labs(y = "Population", color = "") +
  theme_bw() +
  theme(legend.position = "top")
# 
ggsave(filename = paste0(save_folder, "_all_strategies_by_permutation-total_pop.png"),
       permutations_plot, device = 'png', width = 6, height = 4)

# Objectives performance
all_objs_probs <- melt(objectives_probs, id.vars = colnames(objectives_probs[[1]]))
colnames(all_objs_probs)[8] <- "permutation"

all_objs_probs$permutation <- factor(all_objs_probs$permutation, levels = c("low_more_small",
                                                                            "medium_more_small",
                                                                            "high_more_small",
                                                                            "low_more_xlarge",
                                                                            "medium_more_xlarge",
                                                                            "high_more_xlarge"))
all_objs_probs$strategy <- factor(all_objs_probs$strategy, levels = strategies)

objectives_plot_by_permutation <- ggplot(all_objs_probs) +
  geom_point(aes(y = erad_prob, x = total_suppress_prob, color = strategy)) +
  facet_wrap(vars(permutation)) +
  theme_bw()+
  labs(x = "Probability of upper 3 size class suppression", y = "Probability of full eradication")

ggsave(filename = paste0(save_folder, "_objs_comparison_by_permutation_plot.png"),
       objectives_plot_by_permutation, device = 'png', width = 6, height = 4)

objectives_plot_by_strategy <- ggplot(all_objs_probs) +
  geom_point(aes(y = erad_prob, x = total_suppress_prob, color = permutation)) +
  facet_wrap(vars(strategy)) +
  theme_bw()+
  labs(x = "Probability of upper 3 size class suppression", y = "Probability of full eradication")

ggsave(filename = paste0(save_folder, "_objs_comparison_by_strategy_plot.png"),
       objectives_plot_by_strategy, device = 'png', width = 6, height = 4)

## Average costs

# for strategies 3 and 4

summed_costs <- as.data.frame(matrix(NA, nrow = num_variants, ncol = 2))
colnames(summed_costs) <- c("summed_cost", "variant")
summed_costs$variant <- c(1:num_variants)
for(variant in 1:num_variants) {
  variant_costs <- permutation_costs$low_more_xlarge$Strategy_three[permutation_costs$low_more_xlarge$Strategy_three$variant == variant,]
  summed_costs$summed_cost[variant] <- sum(variant_costs$dollars)
}
mean_summed_cost <- mean(summed_costs$summed_cost)



