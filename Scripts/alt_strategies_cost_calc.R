### Cost calculation for each strategy
# Read in all strategy set up objects & cost set up
#source(here("Scripts", "all_strategy_set_up.R"))
source(here("Scripts", "cost_setup.R"))

# Labels & objects needed for all strategies
strategies <- paste0("Strategy_", c("one", "two", "three", "four"))
strat_set_up_file <- paste0("strategy_", c(1:4), "_set_up.R")
names(strat_set_up_file) <- strategies
num_transects <- c("total" = 113)

# List to hold condition costs for all strategies
strategy_condition_costs <- list()

#### Strategy one (non-dynamic, so all replicates have the same cost)
# Read in strategy one set up objects, to get method quarters & days
source(paste0(here("Scripts"), "/",strat_set_up_file[1]))

strategy_condition_costs[[strategies[1]]] <- list()
strategy_condition_costs[[strategies[1]]]$initial <- cost_function(methods = names(method_options$initial$erad_quarters), 
                                                                   erad_days = method_options$initial$erad_days, 
                                                                   erad_quarters = method_options$initial$erad_quarters, 
                                                                   area_size = area_size)


### Strategies 2 - 4 (dynamic, based condition in each set of quarters)

for(strategy in strategies[c(2:4)]) {
  variant_cost <- list()
  source(paste0(here("Scripts"), "/",strat_set_up_file[strategy]))
  condition_costs <- list()
  # Calculate costs for each of the method conditions 
  for(condition in names(method_options)) {
    condition_methods <- method_options[[condition]]
    condition_costs[[condition]] <- cost_function(methods = condition_methods$methods, 
                                                  erad_days = condition_methods$erad_days, 
                                                  erad_quarters = condition_methods$erad_quarters, 
                                                  area_size = area_size,
                                                  num_transects = num_transects,
                                                  num_teams <- condition_methods$cost_num_teams)
  }
  strategy_condition_costs[[strategy]] <- condition_costs
}





dynamic_strategy_cost_calc <- function(permutation_results,
                                       strategy_condition_costs,
                                       strategy,
                                       n = num_variants) {
  variant_cost <- list()
  condition_costs <- strategy_condition_costs[[strategy]]
  for(variant in 1:n) {
    variant_cost[[variant]] <- list()
    var_name <- paste0("variant_", variant)
    condition_df <- permutation_results$condition_plot$condition_record
    # Isolate the variants' conditions for each set
    variant_conditions <- condition_df$condition[condition_df$variant == var_name]
    set_costs <- list()
    summed_cost <- vector()
    for(set in 1:length(variant_conditions)) {
      if(!is.na(variant_conditions[set])) {
        summed_cost[set] <-  condition_costs[[variant_conditions[set]]]$summed_cost
        set_costs[[set]] <- condition_costs[[variant_conditions[set]]]$total_cost
      }
    }
    #variant_cost[[variant]]$set_summed_cost <- summed_cost
    variant_cost[[variant]] <- melt(set_costs)
    colnames(variant_cost[[variant]]) <- c("dollars", "method", "set")
    #variant_cost[[variant]]$total_cost <- sum(variant_cost[[variant]]$set_summed_cost)
  }
  permutation_cost <- melt(variant_cost, id.vars = colnames(variant_cost[[1]]))
  colnames(permutation_cost)[4] <- "variant"
  
  return(permutation_cost)
}

# # Test
# t <- dynamic_strategy_cost_calc(permutation_results,
#                                 strategy_condition_costs,
#                                 strategy = strategies[2])

# Dealing with strategy two when condition wasn't saved and sometimes couldn't be reconstructed

strat_2_dynamic_strategy_cost_calc <- function(permutation_results,
                                               strategy_condition_costs,
                                               strategy = "Strategy_two",
                                               n = num_variants) {
  min_variant_cost <- list()
  max_variant_cost <- list()
  condition_costs <- strategy_condition_costs[[strategy]]
  for(variant in 1:n) {
    min_variant_cost[[variant]] <- list()
    max_variant_cost[[variant]] <- list()
    var_name <- paste0("variant_", variant)
    condition_df <- permutation_results$condition_plot$condition_record
    # Isolate the variants' conditions for each set
    variant_conditions <- condition_df$condition[condition_df$variant == var_name]
    min_set_costs <- list()
    max_set_costs <- list()
    min_summed_cost <- vector()
    max_summed_cost <- vector()
    for(set in 1:length(variant_conditions)) {
      if(!is.na(variant_conditions[set])) {
        if(variant_conditions[set] == "threshold_3") {
          variant_conditions[set] <- "threshold_1"
        } 
        min_summed_cost[set] <-  condition_costs[[variant_conditions[set]]]$summed_cost
        max_summed_cost[set] <- min_summed_cost[set]
        min_set_costs[[set]] <- condition_costs[[variant_conditions[set]]]$total_cost
        max_set_costs[[set]] <- min_set_costs[[set]]
      } else {
        min_summed_cost[set] <- condition_costs$initial$summed_cost
        max_summed_cost[set] <- condition_costs$threshold_2$summed_cost
        min_set_costs[[set]] <- condition_costs$initial$total_cost
        max_set_costs[[set]] <- condition_costs$threshold_2$total_cost
      }
    }
    #variant_cost[[variant]]$set_summed_cost <- summed_cost
    min_variant_cost[[variant]] <- melt(min_set_costs)
    max_variant_cost[[variant]] <- melt(max_set_costs)
    colnames(min_variant_cost[[variant]]) <- c("dollars", "method", "set")
    colnames(max_variant_cost[[variant]]) <- c("dollars", "method", "set")
    #variant_cost[[variant]]$total_cost <- sum(variant_cost[[variant]]$set_summed_cost)
  }
  min_permutation_cost <- melt(min_variant_cost, id.vars = colnames(min_variant_cost[[1]]))
  colnames(min_permutation_cost)[4] <- "variant"
  max_permutation_cost <- melt(max_variant_cost, id.vars = colnames(max_variant_cost[[1]]))
  colnames(max_permutation_cost)[4] <- "variant"
  # Format min and max values into one dataframe
  permutation_cost <- cbind(min_permutation_cost$dollars, max_permutation_cost)
  colnames(permutation_cost)[c(1:2)] <- c("min_dollars", "max_dollars")
  
  return(permutation_cost)
}

# Test
#t <- strat_2_dynamic_strategy_cost_calc(permutation_results, strategy_condition_costs)

