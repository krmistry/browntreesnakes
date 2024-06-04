# Sourcing files with all set up needed to run the function below
library(here)
source(here("all_strategy_set_up.R"))
#source(here("strategy_2_set_up.R"))
#source(here("strategy_3_setup.r"))

parallel_fun <- function(P, 
                         D, 
                         final_time_step,
                         variant,
                         threshold_fun,
                         strategy_name,
                         quarter_time_step) {
  # Source the set up script for the strategy 
  if(strategy_name == "Strategy_two") {
    source(here("strategy_2_set_up.R"))
  } 
  if(strategy_name == "Strategy_three") {
    source(here("strategy_3_setup.r"))
  }
  
  # Using just one scenario to try running in parallel:
  N_reps <- starting_pop[[P]]
  size_dist <- starting_size_dist[[D]]
  # Setting starting N for this variant
  N <- N_reps[variant]
  # List to save IBM results outside of loop
  IBM_results <- list()
  IBM_observed_results <- list()
  IBM_effort_results <- list()
  IBM_all_quarters <- as.data.frame(matrix(nrow = 0, ncol = 8))
  colnames(IBM_all_quarters) <- c("ID", "sex", "repro_prob", "growth_quant", "variable", "SVL", "Quarter", "size_category")
  # Record of effort details as condition changes
  estimation_effort_record <- list()
  estimation_effort_record$condition <- vector()
  estimation_effort_record$primary_sampling_period <- list()
  estimation_effort_record$erad_quarters <- list()
  estimation_effort_record$erad_days <- list()
  # Starting condition
  condition <- "initial"
  # Starting counter
  t <- 1
  while(t <= final_time_step) {
    if(t == 1) { # First time step 
      # Separating out methods for initial condition
      new_methods <- method_options$initial
      # Setting parameter to run IBM initially
      run_type <- "initial"
      # Creating empty variable so function will run 
      last_quarter_df <- NULL
    } else {
      # Identifying methods condition 
      new_methods <- method_options[[condition]]
      # Isolating population from final time step to make the new initial population
      last_quarter_df <- IBM_results[[t-1]]$quarter_timeseries[[erad_quarter_time_step+1]]
      # Setting parameter to run IBM with existing population
      run_type <- "continuing"
    } 
    # Separating out parameters based on condition 
    methods <- new_methods$methods
    coverage <- new_methods$erad_coverage
    primary_sampling_period <- new_methods$primary_sampling_period
    erad_quarters <- new_methods$erad_quarters
    erad_days <- new_methods$erad_days
    ADS_overlap <- new_methods$ADS_overlap_on_transect
    num_teams <- new_methods$num_teams
    
    # Save method options 
    # Recording condition effort for next estimation round
    estimation_effort_record$condition[t] <- condition
    estimation_effort_record$primary_sampling_period[[t]] <- primary_sampling_period
    estimation_effort_record$erad_quarters[[t]] <- erad_quarters
    estimation_effort_record$erad_days[[t]] <- erad_days
    
    # Running IBM model with condition methods
    IBM_results[[t]] <- quarter_operations(initial_N = N, 
                                           initial_size_dist = size_dist, 
                                           p_g = g_density_prob,
                                           lambda = lambda,
                                           total_quarters = erad_quarter_time_step,
                                           total_days = day_time_step,
                                           erad = "on",
                                           erad_method =  methods,
                                           erad_coverage = coverage,
                                           primary_sampling_period = primary_sampling_period,
                                           erad_quarters = erad_quarters,
                                           erad_days = erad_days,
                                           ADS_overlap_on_transect = ADS_overlap,
                                           type_of_run = run_type,
                                           num_teams = num_teams,
                                           last_quarter_df = last_quarter_df)
    
    # Saving IBM observed and effort from this run
    #if(t == 1) {
    IBM_observed_results[[t]] <- IBM_results[[t]]$all_observed
    IBM_effort_results[[t]] <- IBM_results[[t]]$all_effort
    # Correcting the quarters in the IBM that just ran if t > 1, after removing the final quarter (which is just the starting values)
    IBM_results[[t]]$all_quarters <- IBM_results[[t]]$all_quarters[IBM_results[[t]]$all_quarters$Quarter %in% c(1:erad_quarter_time_step),]
    if(quarter_time_step == 2) {
      IBM_results[[t]]$all_quarters$Quarter[IBM_results[[t]]$all_quarters$Quarter == 1] <- t*2 - 1
      IBM_results[[t]]$all_quarters$Quarter[IBM_results[[t]]$all_quarters$Quarter == 2] <- t*2 
    } 
    if(quarter_time_step == 4) {
      IBM_results[[t]]$all_quarters$Quarter[IBM_results[[t]]$all_quarters$Quarter == 1] <- t*2 - 3
      IBM_results[[t]]$all_quarters$Quarter[IBM_results[[t]]$all_quarters$Quarter == 2] <- t*2 - 2
      IBM_results[[t]]$all_quarters$Quarter[IBM_results[[t]]$all_quarters$Quarter == 3] <- t*2 - 1
      IBM_results[[t]]$all_quarters$Quarter[IBM_results[[t]]$all_quarters$Quarter == 4] <- t*2 
    }
    
    # Save IBM results from this run
    saveRDS(IBM_results, file = here("Results", "alt_strategies", strategy_name, "IBM",
                                     paste0("IBM_results_start_pop_", names(starting_pop)[P], "_", 
                                            names(starting_size_dist)[D], "_variant-", variant, "_IBM_set-", t, ".rds")))
    
    # Saving all IBM quarter results
    IBM_all_quarters <- rbind(IBM_all_quarters, IBM_results[[t]]$all_quarters)
    
    # Combining observations and efforts from all previous runs, and adjusting the quarters 
    combined_observed <- unlist(IBM_observed_results, recursive = FALSE)
    combined_effort <- unlist(IBM_effort_results, recursive = FALSE)
    for(quarter in 1:length(combined_observed)) {
      if(nrow(combined_observed[[quarter]]) > 0 && length(combined_observed[[quarter]]) > 0) {
        combined_observed[[quarter]]$quarter <- quarter
      } 
      combined_effort[[quarter]]$quarter <- quarter
    }

    # Estimation model set up
    estimation_inputs <- estimation_inputs_fun(observed_list = combined_observed,
                                               effort_list = combined_effort)
    # Call JAGS Function
    output_jags <- jags(data = estimation_inputs$data,
                        inits = estimation_inputs$inits,
                        parameters.to.save = estimation_inputs$parameters,
                        "removal_model_alt_strategies.jags",
                        n.chains = nc,
                        n.thin = nt,
                        n.iter = ni,
                        n.burnin = nb)
    if(length(IBM_results[[t]]$all_quarters) < quarter_time_step+1) {
      saveRDS(output_jags, file = here("Results", "alt_strategies", strategy_name, "Estimation",
                                       paste0("output_jags_sum_start_pop_", names(starting_pop)[P], "_", 
                                              names(starting_size_dist)[D], "_variant-",variant, "_est_", t, ".RDS")))
      saveRDS(estimation_effort_record, file = here("Results", "alt_strategies", strategy_name, "IBM", 
                                                    paste0("effort_record_start_pop_", names(starting_pop)[P], "_", 
                                                           names(starting_size_dist)[D], "_variant-", variant, ".rds")))
      saveRDS(IBM_all_quarters, file = here("Results", "alt_strategies", strategy_name, "IBM", 
                                            paste0("IBM_all_quarters_start_pop_", names(starting_pop)[P], "_", 
                                                   names(starting_size_dist)[D], "_variant-", variant, ".rds")))
      break
    }
    
    if(t %in% c(1:5, 10, 15, 20)) {
      saveRDS(output_jags, file = here("Results", "alt_strategies", strategy_name, "Estimation",
                                       paste0("output_jags_sum_start_pop_", names(starting_pop)[P], "_", 
                                              names(starting_size_dist)[D], "_variant-",variant, "_est_", t, ".RDS")))
    }
    
    # Mean N estimates vs simulated real N
    est_v_sim_N_plots <- estimated_N_plots(jags_output = output_jags,
                                           all_quarters = IBM_all_quarters,
                                           observed_erad_quarters = estimation_inputs$obs_quarters, 
                                           observed_erad_days = estimation_inputs$obs_days)
    
    
    ## Testing to see if any thresholds are met
    mean_N_df <- est_v_sim_N_plots$mean_N_df
    condition <- threshold_fun(mean_N_df)
    print(condition)
    # Advancing counter by the number of quarters that have been estimated
    t <- t+1
  }
  # Saving effort record and the IBM all_quarters data frame from the entire time period
  saveRDS(estimation_effort_record, file = here("Results", "alt_strategies", strategy_name, "IBM", 
                                                paste0("effort_record_start_pop_", names(starting_pop)[P], "_", 
                                                       names(starting_size_dist)[D], "_variant-", variant, ".rds")))
  saveRDS(IBM_all_quarters, file = here("Results", "alt_strategies", strategy_name, "IBM", 
                                        paste0("IBM_all_quarters_start_pop_", names(starting_pop)[P], "_", 
                                               names(starting_size_dist)[D], "_variant-", variant, ".rds")))
  
  return(estimation_effort_record$condition)
}

# test <- parallel_fun(P = 2,
#                      D = 1,
#                      final_time_step = 2,
#                      variant = 1,
#                      threshold_fun = strat_3_threshold_fun,
#                      strategy_name = "Strategy_three",
#                      quarter_time_step = 4)
# 
# test_output_IBM <- readRDS(file = here("Results", "alt_strategies", "Strategy_two", "IBM",
#                                         paste0("IBM_results_parallel_test_variant-", 1, "_IBM_set-", 2, ".rds")))
# test_output_jags <- readRDS(file = here("Results", "alt_strategies", "Strategy_two", "Estimation",
#                                          paste0("output_jags_paralell_test_variant-",1, "_est_", 2, ".RDS")))
# 




