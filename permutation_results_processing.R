### Processing & plotting permutation results for each strategy

library(reshape2)
library(here)
library(tidyr)

source(here("Scripts/all_strategy_set_up.R"))
source(here("Scripts/model_evaluation.R"))
source(here("Scripts/results_processing_functions.R"))
source(here("Scripts/alt_strategies_cost_calc.R"))

# Strategy  & permutation names
strategies <- paste0("Strategy_", c("one", "two", "three", "four"))
permutations <- vector()
counter <- 1
for(p in 1:length(starting_pop)) {
  for(d in 1:length(starting_size_dist)) {
    permutations[counter] <- paste0(names(starting_pop)[p], "_", names(starting_size_dist)[d])
    counter <- counter + 1
  }
}

# Plot labels & colors
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
plot_labels$strategy <- c("ADS only",
                          "all methods",
                          "ground-based methods",
                          "ADS plus monitoring")
names(plot_labels$strategy) <- strategies
plot_labels$permutation <- c("low density,\n more small snakes",
                            "low density,\n more x-large snakes",
                            "medium density,\n more small snakes",
                            "medium density,\n more x-large snakes",
                            "high density,\n more small snakes",
                            "high density,\n more x-large snakes")
names(plot_labels$permutation) <- permutations

plot_colors <- list()
plot_colors$strategy <- hue_pal()(length(strategies))
names(plot_colors$strategy) <- strategies
plot_colors$permutation <- hue_pal()(length(permutations))
names(plot_colors$permutation) <- permutations

# Checking for and creating if necessary all results folders
save_folder <- paste0(here("Results", "alt_strategies"), "/")

#### Analyzing results by permutation & strategy ####

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

##### Objectives performance #####
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

##### Cost calculations for strategies (dynamic strategies use means across variants) #####
indiv_cost_categories <- c("transects", erad_methods)

mean_permutation_costs <- list()
# For each permutation
for(permutation_name in permutations) {
  strategy_costs <- list()
  # For strategy one (static strategy)
  strategy_costs$Strategy_one <- as.data.frame(matrix(NA, nrow = length(indiv_cost_categories), 
                                                      ncol = 2))
  colnames(strategy_costs$Strategy_one) <- c("method", "mean_cost")
  strategy_costs$Strategy_one$method <- indiv_cost_categories
  strategy_costs$Strategy_one$mean_cost <- 0
  strategy_costs$Strategy_one$mean_cost
  strategy_costs$Strategy_one$mean_cost[strategy_costs$Strategy_one$method == "ADS"] <- permutation_costs[[permutation_name]]$Strategy_one$total_cost$ADS
  
  ### For strategy 2 - using max cost, because its really not a huge difference between min and max. 
  # If the totals are close enough to another strategy that it might affect the optimization, 
  # then I'll reconsider how to do this
  strategy <- strategies[2]
  summed_costs <- as.data.frame(matrix(NA, nrow = num_variants, ncol = 6))
  colnames(summed_costs) <- c("variant",indiv_cost_categories)
  summed_costs$variant <- c(1:num_variants)
  for(variant in 1:num_variants) {
    # Combine all costs for each variant
    variant_costs <- permutation_costs[[permutation_name]][[strategy]][permutation_costs[[permutation_name]][[strategy]]$variant == variant,]
    # Combine costs for each method for each variant
    for(method in indiv_cost_categories) {
      summed_costs[[method]][variant] <- sum(variant_costs$max_dollars[variant_costs$method == method])
    }
  }
  # Calculate mean costs for each strategy (across variants)
  strategy_costs[[strategy]] <- as.data.frame(matrix(NA, nrow = length(indiv_cost_categories), 
                                                                ncol = 2))
  colnames(strategy_costs[[strategy]]) <- c("method", "mean_cost")
  strategy_costs[[strategy]]$method <- indiv_cost_categories
  for(method in 1:length(indiv_cost_categories)) {
    strategy_costs[[strategy]]$mean_cost[method] <- mean(summed_costs[[indiv_cost_categories[method]]])
  }
  
  # For strategies 3 and 4
  for(strategy in strategies[c(3:4)]) {
    summed_costs <- as.data.frame(matrix(NA, nrow = num_variants, ncol = 6))
    colnames(summed_costs) <- c("variant", indiv_cost_categories)
    summed_costs$variant <- c(1:num_variants)
    for(variant in 1:num_variants) {
      # Combine all costs for each variant
      variant_costs <- permutation_costs[[permutation_name]][[strategy]][permutation_costs[[permutation_name]][[strategy]]$variant == variant,]
      # summed_costs$summed_cost[variant] <- sum(variant_costs$dollars)
      # Combine costs for each method for each variant
      for(method in indiv_cost_categories) {
        summed_costs[[method]][variant] <- sum(variant_costs$dollars[variant_costs$method == method])
      }
    }
    # Calculate mean costs for each strategy (across variants)
    strategy_costs[[strategy]] <- as.data.frame(matrix(NA, nrow = length(indiv_cost_categories), 
                                                                  ncol = 2))
    colnames(strategy_costs[[strategy]]) <- c("method", "mean_cost")
    strategy_costs[[strategy]]$method <- indiv_cost_categories
    for(method in 1:length(indiv_cost_categories)) {
      strategy_costs[[strategy]]$mean_cost[method] <- mean(summed_costs[[indiv_cost_categories[method]]])
    }
    #strategy_costs[[strategy]]$summed_cost <- mean(summed_costs$summed_cost)
  }
  
  mean_permutation_costs[[permutation_name]] <- strategy_costs
}
# Combining above into one dataframe with all mean costs
all_mean_costs <- melt(mean_permutation_costs, id.vars = colnames(mean_permutation_costs[[1]][[1]]))
colnames(all_mean_costs)[c(3:4)] <- c("strategy", "permutation")

# Combining summed costs and objectives for plotting
costs_vs_obj_probs <- all_objs_probs
for(permutation_name in permutations) {
  for(strategy in strategies) {
    perm_strat_costs <- all_mean_costs[all_mean_costs$strategy == strategy & 
                                         all_mean_costs$permutation == permutation_name,]
    costs_vs_obj_probs$mean_cost[costs_vs_obj_probs$strategy == strategy & costs_vs_obj_probs$permutation == permutation_name] <- sum(perm_strat_costs$mean_cost)
  }
}


##### Calculate RMSE & PRD #####
# - there could be an appendix of how the estimation improves as data accumulates 
#   (so at each estimation step that I saved, 1 - 5, and then 10, and for strategy 2, 15 and 20)



#### Figures & Tables ####

##### Objectives vs cost plots  #####

# Eradication probabilty vs cost plots
cost_vs_erad_objs_plot_by_permutation <- ggplot(costs_vs_obj_probs) +
  geom_point(aes(y = erad_prob, x = mean_cost/1000000, color = strategy)) +
  scale_y_continuous(limits = c(0,1)) +
  facet_wrap(vars(permutation), 
             labeller = labeller(permutation = plot_labels$permutation)) +
  theme_bw()+
  scale_color_manual(values = plot_colors$strategy,
                     labels = plot_labels$strategy) +
  theme(legend.position="bottom") +
  #theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.75)) +
  labs(x = "Projected mean cost (dollars in millions)", 
       y = "Probability of full eradication",
       color = "Strategy")

ggsave(filename = paste0(save_folder, "costs_vs_erad_obj_comparison_by_permutation_plot.png"),
       cost_vs_erad_objs_plot_by_permutation, device = 'png', width = 7, height = 4)

cost_vs_erad_objs_plot_by_strategy <- ggplot(costs_vs_obj_probs) +
  geom_point(aes(y = erad_prob, x = mean_cost/1000000, color = permutation)) +
  scale_y_continuous(limits = c(-0.01,1)) +
  
  facet_wrap(vars(strategy)) +
  theme_bw()+
  labs(x = "Projected mean cost (dollars in millions)", y = "Probability of full eradication")

ggsave(filename = paste0(save_folder, "costs_vs_erad_obj_comparison_by_strategy_plot.png"),
       cost_vs_erad_objs_plot_by_strategy, device = 'png', width = 6, height = 4)

# Upper 3 size class suppression probabilty vs cost plots
cost_vs_upper_3_objs_plot_by_permutation <- ggplot(costs_vs_obj_probs) +
  geom_point(aes(y = upper_3_suppress_prob, x = mean_cost/1000000, color = strategy)) +
  scale_y_continuous(limits = c(0,1)) +
  facet_wrap(vars(permutation), 
             labeller = labeller(permutation = plot_labels$permutation)) +
  scale_color_manual(values = plot_colors$strategy,
                     labels = plot_labels$strategy) +
  theme_bw()+
  labs(x = "Projected mean cost (dollars in millions)", 
       y = "Probability of suppression upper 3 size classes",
       color = "Strategy")

ggsave(filename = paste0(save_folder, "costs_vs_upper_3_obj_comparison_by_permutation_plot.png"),
       cost_vs_upper_3_objs_plot_by_permutation, device = 'png', width = 6, height = 4)

cost_vs_upper_3_objs_plot_by_strategy <- ggplot(costs_vs_obj_probs) +
  geom_point(aes(y = upper_3_suppress_prob, x = mean_cost/1000000, color = permutation)) +
  scale_y_continuous(limits = c(-0.01,1)) +
  facet_wrap(vars(strategy)) +
  theme_bw()+
  labs(x = "Projected mean cost (dollars in millions)", y = "Probability of suppression upper 3 size classes")

ggsave(filename = paste0(save_folder, "costs_vs_upper_3_obj_comparison_by_strategy_plot.png"),
       cost_vs_upper_3_objs_plot_by_strategy, device = 'png', width = 6, height = 4)

# Total suppression probabilty vs cost plots
cost_vs_total_supp_objs_plot_by_permutation <- ggplot(costs_vs_obj_probs) +
  geom_point(aes(y = total_suppress_prob, x = mean_cost/1000000, color = strategy)) +
  scale_y_continuous(limits = c(0,1)) +
  facet_wrap(vars(permutation), 
             labeller = labeller(permutation = plot_labels$permutation)) +
  scale_color_manual(values = plot_colors$strategy,
                     labels = plot_labels$strategy) +
  theme_bw()+
  labs(x = "Projected mean cost (dollars in millions)", 
       y = "Probability of total population suppression",
       color = "Strategy")

ggsave(filename = paste0(save_folder, "costs_vs_total_supp_obj_comparison_by_permutation_plot.png"),
       cost_vs_total_supp_objs_plot_by_permutation, device = 'png', width = 6, height = 4)

cost_vs_total_supp_objs_plot_by_strategy <- ggplot(costs_vs_obj_probs) +
  geom_point(aes(y = total_suppress_prob, x = mean_cost/1000000, color = permutation)) +
  scale_y_continuous(limits = c(-0.01,1)) +
  facet_wrap(vars(strategy)) +
  theme_bw()+
  labs(x = "Projected mean cost (dollars in millions)", y = "Probability of total population suppression")

ggsave(filename = paste0(save_folder, "costs_vs_total_supp_obj_comparison_by_strategy_plot.png"),
       cost_vs_total_supp_objs_plot_by_strategy, device = 'png', width = 6, height = 4)



#### Plots of methods & costs ####


##### Method combinations in each threshold condition for each strategy #####


##### Actual methods used across replicates for each dynamic strategy #####


##### Format dataframe for a table of total costs for each strategy #####



#### Plots of accuracy & precision metrics (RMSE, PRD and coverage) & dataframe for table ####



