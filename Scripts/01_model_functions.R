############## Functions ###################


### Gompertz growth function 
## Inputs
# - snake = individual snake SVL
# - p_g = probability (0 - 1) for the Bernoulli probability for whether or not growth will 
#         occur
daily_gompertz_growth_fun <- function(snake, 
                                      p_g) {
  if (rbinom(1, 1, prob = p_g) == 1) {
    # Selecting correct growth coefficients based on individual's sex and growth quantile
    growth_coefs <- gompertz_coefficients[[snake$growth_quant]][snake$sex]
    # Snake grows
    snake$SVL <- snake$SVL + growth_coefs[1,] + growth_coefs[2,]*log(snake$SVL/700)
  } else {
    snake$SVL <- snake$SVL
  }
  return(snake)
}


########## Reproduction functions ##############################

### Function to determine an individuals' probability of being reproductively mature based on 
### SVL 
## Inputs:
# - snake_SVL - individual snake's SVL
maturity_fun <- function(snake_SVL) {
  if(snake_SVL < 910) {
    maturity <- 0
  } else if (snake_SVL >= 910 & snake_SVL < 1025) {
    maturity <- pnorm(snake_SVL, maturity_mean, maturity_sd)
  } else {
    maturity <- 1
  }
  return(maturity)
}


### Function to select females that will reproduce in this quarter
## Inputs
# - start_pop = dataframe of population at the start of the quarter 
#               (columns are ID, sex, SVL, growth quantile and reproductive probability)
# - density_prob = probability (0 - 1) for the Bernoulli probability for whether or not 
#                  snake reproduces in this quarter; it's density dependent

repro_females_fun <- function(start_pop,
                              density_prob) {
  # Extracting the females from this population
  females <- start_pop[start_pop$sex == "F",]
  # Coin flip for whether or not reproduction occurs
  for(snake in 1:nrow(females)) {
    # Coin flip of whether reproduction occurs is based on both density dependence and the     probablity that the individual snake will be able to reproduce based on its size
    if (rbinom(1, 1, prob = density_prob*females$repro_prob[snake]) == 1) {
      females$mom_status[snake] <- 1
    } else {
      females$mom_status[snake] <- 0
    }
  }
  moms <- females[females$mom_status == 1,]
  return(select(moms, -c("mom_status")))
}


### Function to generate offspring
## Inputs
# - mom_pop = dataframe of reproducing female population at the start of the quarter 
#             (columns are ID, sex, SVL, growth quantile and reproductive probability)
# - time_step = the (numeric) quarter the loop is running, used to help randomly assign IDs

gen_offspring_fun <- function(mom_pop,
                              time_step,
                              offspring_lambda) {
  # Version one - using the same value for lambda for all sizes (could separate this out      into size buckets, or maybe even could make it continuous based on Nafus data)
  total_offspring <- sum(rpois(nrow(mom_pop), offspring_lambda))
  # Creating empty dataframe to put offspring data into
  offspring <- as.data.frame(matrix(NA, nrow = total_offspring, ncol = 3))
  colnames(offspring) <- c("ID", "SVL", "sex")
  # Filling in IDs, hatchling size (with some variation), sex, reproductive probability (0), and randomly selecting a growth quantile
  if(time_step <= 25) {
    offspring$ID <- paste0(letters[time_step + 1], 
                           sample(1000:9000, total_offspring, replace = F))
  } else {
    x <- floor(time_step/26)
    offspring$ID <- paste0(letters[x], 
                           letters[time_step + 1 - (26*x)], 
                           sample(1000:9000, total_offspring, replace = F))
  }
  
  offspring$SVL <- runif(total_offspring, 250, 350) # Lower limit from Rodda and Savidge 2007 
  offspring$sex <- sample(c("M", "F"), total_offspring, replace = T)
  offspring$repro_prob <- 0
  offspring$growth_quant <- sample(growth_quantiles, total_offspring, replace = T)
  
  ## Combining the offspring with mothers and return as one dataframe
  return(offspring)
}  



### Function for natural mortality
## Inputs
# - pop = dataframe of population (columns are ID, sex, SVL, growth quantile and 
#       reproductive probability)
# - size_limits = dataframe with the lower and upper limits of SVL for the 4 size classes
daily_mortality_fun <- function(pop, 
                                size_limits) {
  # Creating empty vector to hold IDs fo dead snakes
  dead_snakes <- vector()
  # Loop to apply Bernoulli function to each snake to see how many survive
  for(snake in 1:nrow(pop)) {
    if(pop$SVL[snake] <= size_limits$upper[1]) {
      if (rbinom(1, 1, prob = N_mortality[1]) == 1) { 
        dead_snakes[snake] <- pop$ID[snake]
      }
    } else if (pop$SVL[snake] > size_limits$lower[2] & 
               pop$SVL[snake] <= size_limits$upper[2]) {
      if (rbinom(1, 1, prob = N_mortality[2]) == 1) { 
        dead_snakes[snake] <- pop$ID[snake]
      }
    } else if (pop$SVL[snake] > size_limits$lower[3] & 
               pop$SVL[snake] <= size_limits$upper[3]) {
      if (rbinom(1, 1, prob = N_mortality[3]) == 1) { 
        dead_snakes[snake] <- pop$ID[snake]
      } 
    } else if (pop$SVL[snake] > size_limits$lower[4]) {
      if (rbinom(1, 1, prob = N_mortality[4]) == 1) { 
        dead_snakes[snake] <- pop$ID[snake]
      }
    }
  }
  # Get rid of NAs in dead snake ID vector
  dead_snakes <- dead_snakes[!is.na(dead_snakes)]
  
  # Extract surviving snakes 
  pop_2 <- pop[!pop$ID %in% dead_snakes,]
  return(pop_2)
}


### Function to create density dependent parameter for reproduction
## Inputs
# - K = carrying capacity-type value (based on very high density/ha), which controls the 
#       variance of the half-normal distribution
# - current_density = the current total population
DD_param_fun <- function(K,
                         current_density) {
  # Calculate sigma for scaling standard deviation using carrying capacity
  sd <- (K*0.975)/1.96
  # Calculate PDF when density is 0 (highest density probability value)
  max_prob_value <- dhalfnorm(0, theta = ((sqrt(pi/2))/sd)) 
  # Calculate PDF value for the current density
  PDf_value <- dhalfnorm(current_density, theta = ((sqrt(pi/2))/sd))
  # Calculate PDF value for the current density
  PDf_value <- dhalfnorm(current_density, theta = ((sqrt(pi/2))/sd))
  # Use maximum value to scale current density's PDF value and produce A parameter
  p <- (PDf_value/max_prob_value)
  return(p)
}


### Function to generate an initial population to start the timeseries 
## Inputs:
# - init_N = size of total population in first time step
# - init_size_dist = distribution of the SVL sizes in the initial population
init_pop_fun <- function(init_N, 
                         init_size_dist,
                         size_limits) {
  # Checking that the init_N multiplied by the init_size_dist values is the same as init_N,   and if it's not than init_N is slightly adjusted (runif will round these values below and   it won't match up if I don't do this)
  if(init_N != sum(round(init_N*init_size_dist))) {
    init_N <- sum(round(init_N*init_size_dist))
  } else {
    init_N <- init_N
  }
  
  # Setting up initial population dataframe:
  init_pop <- as.data.frame(matrix(NA, nrow = init_N, ncol = 5))
  colnames(init_pop) <- c("ID", "SVL", "sex", "repro_prob", "growth_quant")
  # Creating initial population characteristics
  init_pop$ID <- paste0("A", sample(1000:9000, init_N))
  init_pop$SVL <- c(runif(init_N*init_size_dist[1], size_limits$lower[1], size_limits$upper[1]), 
                    runif(init_N*init_size_dist[2], size_limits$lower[2]+1, size_limits$upper[2]),
                    runif(init_N*init_size_dist[3], size_limits$lower[3]+1, size_limits$upper[3]), 
                    runif(init_N*init_size_dist[4], size_limits$lower[4]+1, size_limits$upper[4]))
  init_pop$sex <- sample(c("M", "F"), init_N, replace = T)
  # Using SVL to determine probability of sexual maturity (only relevant for females, but easier to calculate for all - don't keep this long term though, it might be confusing)
  # At the same time, randomly assign each individual to a growth quantile - initially, each of the 5 possible quantiles will be equally probable (0.25, 0.5, 0.75, 0.9, 0.95), depending on if their starting SVL is below the asymptote of the quantile
  for(snake in 1:nrow(init_pop)) {
    # Set up reproductive probability value for each individual
    init_pop$repro_prob[snake] <- maturity_fun(init_pop$SVL[snake])
    # Set up growth quantiles for each individual, based on starting SVL
    if (init_pop$SVL[snake] <= min(gompertz_asymptotes$tau_0.25)) {
      init_pop$growth_quant[snake] <- sample(tau, 1)
    } else if (init_pop$SVL[snake] <= min(gompertz_asymptotes$tau_0.5)) {
      init_pop$growth_quant[snake] <- sample(tau[-1], 1)
    } else if (init_pop$SVL[snake] <= min(gompertz_asymptotes$tau_0.75)) {
      init_pop$growth_quant[snake] <- sample(tau[-c(1:2)], 1)
    } else if (init_pop$SVL[snake] <= min(gompertz_asymptotes$tau_0.9)) {
      init_pop$growth_quant[snake] <- sample(tau[-c(1:3)], 1)
    } else {
      init_pop$growth_quant[snake] <- sample(tau[-c(1:4)], 1)
    }
  }
  return(init_pop)
}

### Function to perform all eradication methods 
## Inputs:
# - day_pop = starting population for that day 
# - coverage = the percentage of the total area that the eradication method will reach
# - size_affected = what size of snake is susceptible to the method (lower and upper 
#                   limit for now - might change into a probability based on size later
# - mortality_prob = mortality rate of method if encountered 

eradication_fun <- function(day_pop, 
                            coverage, 
                            size_affected = c(0, Inf), 
                            mortality_prob) {
  # Separating out the proportion of the population that will be affected by eradication (coverage) from the non-encountering population
  encounter_pop <- slice_sample(day_pop, prop = coverage)
  no_encounter_pop <- day_pop[!day_pop$ID %in% encounter_pop$ID, ]
  # If snake is in the size range to be affected, perform Bernoulli draw to determine if it encounters the method (encounter = mortality for now)
  for(snake in 1:nrow(encounter_pop)) {
    if(encounter_pop$SVL[snake] > size_affected[1] & encounter_pop$SVL[snake] < size_affected[2]) {
      if (rbinom(1, 1, prob = mortality_prob) == 1) { 
        encounter_pop$encounter[snake] <- 1
      } else {
        encounter_pop$encounter[snake] <- 0
      }
    } else {
      encounter_pop$encounter[snake] <- 0
    }
  }
  # Separating dead snakes from alive snakes
  dead_snakes <- encounter_pop[encounter_pop$encounter == 1,]
  alive_snakes <- encounter_pop[encounter_pop$encounter == 0,]
  
  # Recombining surviving snakes with  non-encounter snakes
  day_pop <- rbind(select(alive_snakes, -c("encounter")), no_encounter_pop)
  # Return both dead and live snakes (eventually, only return live snakes)
  return(list(dead_snakes = dead_snakes, alive_snakes = day_pop))
}



erad_timing_fun <- function(day, 
                            quarter,
                            erad_method,
                            growing_pop){
  # Creating empty lists to hold dead snakes, alive snakes and effort lists for each method
  observed <- list()
  unobserved <- list()
  survivors <- list()
  effort <- list()
  # Perform each eradication that occurs on this day & quarter sequentially
  survivors[[1]] <- growing_pop
  for(method in 1:length(erad_method)) {
    method_name <- erad_method[method]
    #pop <- survivors[[method]]
    if(quarter %in% erad_quarters[[method_name]] & day %in% erad_days[[method_name]]){
      erad_pop <- eradication_fun(survivors[[method]], 
                                  coverage = erad_coverage[[method_name]], 
                                  erad_methods_size_affected[[method_name]],
                                  mortality_prob_erad_methods[[method_name]])
      survivors[[method + 1]] <- erad_pop$alive_snakes
      # Eradication with bodies occurs:
      if(nrow(erad_pop$dead_snakes) > 0) {
        if (method_name %in% erad_methods[c(2:3)]) {
          observed[[method_name]] <- erad_pop$dead_snakes
          observed[[method_name]]$day <- day
          observed[[method_name]]$quarter <- quarter
        } else { # Eradication with no bodies
          unobserved[[method_name]] <- erad_pop$dead_snakes
          unobserved[[method_name]]$day <- day
          unobserved[[method_name]]$quarter <- quarter
        }
        effort[[method_name]] <- as.data.frame(cbind("effort" = effort_erad_methods[[method_name]],
                                                     "day" = day,
                                                     "quarter" = quarter))
      }
    }
    else {
      survivors[[method + 1]] <- survivors[[method]]
    }
  }
  final_survivors <- survivors[[length(erad_method)+1]]
  return(list(observed_snakes = observed, 
              unobserved_dead_snakes = unobserved,
              surviving_snakes = final_survivors,
              all_effort = effort))
}

# # Test
# day <- 14
# quarter <- 3
# growing_pop <- init_pop_fun(N, size_dist, size_class_limits)
# x <- erad_timing_fun(day, quarter, erad_methods, growing_pop)
#
# remove(day)
# remove(growing_pop)
# remove(x)


##################### Functions to run daily & quarterly operations #################

### Function to run daily model operations - natural mortality, individual growth and 
### eradication if it occurs
## Inputs:
# - first_day_pop = dataframe of the population at the beginning of the quarter
# - p_g = probability (0-1) of growth occuring using Bernoulli draw
# - moms = the females that are selected to reproduce in this quarter
# - total_days = the number of days to loop over (in case we change from quarter time steps)

daily_operations <- function(first_day_pop,
                             p_g,
                             moms,
                             total_days, 
                             methods,  
                             erad = c("on", "off"), 
                             quarter) {
  daily_pop <- list()
  erad_results <- list()
  ## Growth and mortality on a daily scale within the quarter 
  daily_pop[[1]] <- first_day_pop
  # Loop for each day in 90 days
  for(day in 1:(total_days-1)) {
    #  All snakes are exposed to natural mortality
    pop <- daily_pop[[day]]
    surviving_pop <- daily_mortality_fun(pop, size_class_limits)
    # Excluding reproductive females to create population who has the potential to grow in this quarter
    growing_pop <- surviving_pop[!surviving_pop$ID %in% moms$ID, ]
    # Separating out surviving moms to add back in after growth
    surviving_moms <- surviving_pop[surviving_pop$ID %in% moms$ID, ]
    
    ## Each snake grows (or not) based on it's individual growth quantile
    for(snake in 1:nrow(growing_pop)) {
      growing_pop[snake,] <- daily_gompertz_growth_fun(growing_pop[snake,],
                                                       p_g)
    }
    ## If eradication occurs in this day & quarter
    if(erad == "on") {
      if(quarter %in% unique(unlist(erad_quarters)) & day %in% unique(unlist(erad_days))) {
        erad_results[[day]] <- erad_timing_fun(day, 
                                               quarter, 
                                               methods, 
                                               growing_pop)
        # Surviving snakes are combined with surviving reproducers to populate the next day
        daily_pop[[day + 1]] <- rbind(erad_results[[day]]$surviving_snakes, surviving_moms)
      } else {
        daily_pop[[day + 1]] <- rbind(growing_pop, surviving_moms)
      }
    } else {
      daily_pop[[day + 1]] <- rbind(growing_pop, surviving_moms)
    }
  }
  names(daily_pop) <- paste0("day_", c(1:total_days))
  if(length(erad_results) > 0) {
    names(erad_results) <- paste0("day_", c(1:length(erad_results)))
  }
  
  return(list(daily_timeseries = daily_pop,
              all_erad_results = erad_results))
}

# # Test
# total_days <- 91
# first_day_pop <- init_pop_fun(N, size_dist, size_class_limits)
# moms <- repro_females_fun(start_pop = first_day_pop,
#                           density_prob = DD_param_fun(K, nrow(first_day_pop)))
# z <- daily_operations(first_day_pop,
#                       g_density_prob,
#                       moms,
#                       total_days,
#                       erad_methods,
#                       erad = "on",
#                       quarter)
# remove(first_day_pop)
# remove(moms)
# remove(z)


# Reformatting eradication results
erad_results_format <- function(erad_results) {
  erad_melt <- erad_results[unlist(lapply(erad_results, length) != 0)]
  quarter_unobserved <- as.data.frame(matrix(nrow = 0, ncol = 9))
  colnames(quarter_unobserved) <- c("ID", "SVL", "sex", "repro_prob", "growth_quant",
                                    "encounter", "day", "quarter", "L1") 
  quarter_observed <- as.data.frame(matrix(nrow = 0, ncol = 9))
  colnames(quarter_observed) <- c("ID", "SVL", "sex", "repro_prob", "growth_quant",
                                    "encounter", "day", "quarter", "L1") 
  quarter_effort <- as.data.frame(matrix(nrow = 0, ncol = 4))
  colnames(quarter_effort) <- c("effort", "day", "quarter", "L1")
  for(day in 1:length(erad_melt)){
    # Melting unobserved (dead) snakes into one dataframe
    if (length(erad_melt[[day]]$unobserved_dead_snakes) > 0) {
      # colnames(quarter_unobserved) <- c(colnames(erad_melt[[day]]$unobserved_dead_snakes[[1]]), "L1")
      quarter_unobserved <- rbind(quarter_unobserved, 
                                  melt(erad_melt[[day]]$unobserved_dead_snakes, 
                                       id.vars = colnames(erad_melt[[day]]$unobserved_dead_snakes[[1]])))
    } else {
      quarter_unobserved <- quarter_unobserved
    }
    # Melting observed snakes into one dataframe
    if (length(erad_melt[[day]]$observed_snakes) > 0) {
      colnames(quarter_observed) <- c(colnames(erad_melt[[day]]$observed_snakes[[1]]), "L1")
      quarter_observed <- rbind(quarter_observed, 
                                  melt(erad_melt[[day]]$observed_snakes, 
                                       id.vars = colnames(erad_melt[[day]]$observed_snakes[[1]])))
    } else {
      quarter_observed <- quarter_observed
    }
    # Melting effort into one dataframe
    if (length(erad_melt[[day]]$all_effort) > 0) {
      colnames(quarter_effort) <- c(colnames(erad_melt[[day]]$all_effort[[1]]), "L1")
      quarter_effort <- rbind(quarter_effort, 
                                melt(erad_melt[[day]]$all_effort, 
                                     id.vars = colnames(erad_melt[[day]]$all_effort[[1]])))
    } else {
      quarter_effort <- quarter_effort
    }
  }
  colnames(quarter_unobserved)[ncol(quarter_unobserved)] <- "method"
  colnames(quarter_observed)[ncol(quarter_observed)] <- "method"
  colnames(quarter_effort)[ncol(quarter_effort)] <- "method"
  
  return(list(unobserved = quarter_unobserved,
              observed = quarter_observed,
              effort = quarter_effort))
}

## Test
# y <- erad_results_format(z$all_erad_results)
#
# remove(y)



### Function to run quarterly operations - initialize first quarter, then reproduction 
### and daily operations
## Inputs:
# - initial_N = size of total population in quarter 1
# - initial_size_dist = distribution of size classes in first quarter
# - p_g = probability (0-1) of growth occuring using Bernoulli draw
# - total_quarters = the number of quarters to loop over
# - total_days = the number of days to loop over (in case we change from quarter time steps)
# - erad = if eradication occurs, use "on" to allow eradication methods to run, the 
#         default is "off", which prevents them from running

quarter_operations <- function(initial_N, 
                               initial_size_dist, 
                               p_g,
                               lambda,
                               total_quarters,
                               total_days, 
                               erad = c("on", "off"),
                               erad_method) {
  # Empty lists to put population from each time step in (records the population at the 
  # beginning of each quarter, and for all the days in the final quarter)
  quarter_timeseries_pop <- list()
  daily_results <- list()
  # Setting up empty dataframe to record IDs for reproducing females in each quarter
  repro_females_timeseries <- list()
  # Empty lists to hold the effort that will be expended if eradication occurs (summed in 
  # each quarter):
  effort_record <- list()
  # for(method in erad_methods){
  #   effort_record[[method]] <- list()
  # }
  # Creating a list of lists to keep track of "observed" snakes from visual surveys & traps
  # and "unobserved" but still dead snakes from ADS and bait tubes  
  observed_snakes <- list()
  unobserved_dead_snakes <- list()
  
  # Setting up the initial population in the first time step
  quarter_timeseries_pop[[1]] <- init_pop_fun(init_N = initial_N,
                                              init_size_dist = initial_size_dist,
                                              size_limits = size_class_limits)
  
  for(quarter in 1:total_quarters) {
    tic(paste0("Quarter ", quarter))
    # Setting up empty list to hold any observed and unobserved dead snakes in this quarter, if observation occurs, as well as the effort expended for each method
    observed_snakes[[quarter]] <- list()
    unobserved_dead_snakes[[quarter]] <- list()
    # for(method in erad_methods){
    #   effort_record[[method]][[quarter]] <- list()
    # }
    effort_record[[quarter]] <- list()
    # Isolate reproducing females for this quarter 
    # After first quarter, check for the females who reproduced last quarter and 
    # skip them this quarter
    if (quarter == 1) {
      moms <- repro_females_fun(start_pop = quarter_timeseries_pop[[quarter]],
                                density_prob = DD_param_fun(K, initial_N))
    } else {
      recent_mom_IDs <- repro_females_timeseries[[quarter - 1]]$ID
      # Update density dependent parameter & r_density
      current_N <- nrow(quarter_timeseries_pop[[quarter]])
      r_density_prob <- DD_param_fun(K, current_N)
      # Choose which females will reproduce
      moms <- repro_females_fun(start_pop = quarter_timeseries_pop[[quarter]][!quarter_timeseries_pop[[quarter]]$ID %in% recent_mom_IDs,],
                                density_prob = r_density_prob)
    }
    # If no reproductive females available, then stop loop and go to next starting value
    if(nrow(moms) < 1) {
      print("0 reproducing females")
      break
    }
    # Producing new offspring to join population at the end   of the quarter
    offspring <- gen_offspring_fun(mom_pop = moms,
                                   time_step = quarter,
                                   offspring_lambda = lambda)
    
    # Keeping track of which females reproduced in which quarters (next step will be to use     this to exclude females who've reproduced in the last 2 quarters)
    repro_females_timeseries[[quarter]] <- moms
    
    # Running daily operations - natural mortality and individual growth
    daily_results <- daily_operations(first_day_pop = quarter_timeseries_pop[[quarter]],
                                      p_g = p_g,
                                      moms = moms,
                                      total_days = total_days,
                                      methods = erad_method,  
                                      erad = erad, 
                                      quarter = quarter)
    daily_timeseries_pop <- daily_results$daily_timeseries
    if(length(daily_results$all_erad_results) > 0) {
    # Filling the observed, unobserved dead snakes and effort lists:
    erad_reformatted <- erad_results_format(erad_results = daily_results$all_erad_results)
    observed_snakes[[quarter]] <- erad_reformatted$observed
    unobserved_dead_snakes[[quarter]] <- erad_reformatted$unobserved
    effort_record[[quarter]] <- erad_reformatted$effort
    } 
    # Adding offspring and the reproducing females back into the surviving and grown general   population to be the start of the next quarter
    quarter_timeseries_pop[[quarter + 1]] <- rbind(daily_timeseries_pop[[total_days]], offspring)
    # Updating sexual maturity status after the quarter's worth of growth
    for(snake in 1:nrow(quarter_timeseries_pop[[quarter + 1]])) {
      quarter_timeseries_pop[[quarter + 1]]$repro_prob[snake] <- maturity_fun(quarter_timeseries_pop[[quarter + 1]]$SVL[snake])
    }
    toc()
  }
  # Formating quarter results for plotting
  all_quarters <- results_format_fun(quarter_results = quarter_timeseries_pop,
                                     size_limits = size_class_limits)
  
  # Return the quarter timeseries, the last quarter's daily timeseries, 
  # all reproductive females, and the quarter timeseries formatted for plotting
  return(list(quarter_timeseries = quarter_timeseries_pop,
              last_daily_timeseries = daily_timeseries_pop,
              all_repro_females = repro_females_timeseries,
              all_quarters = all_quarters,
              all_observed = observed_snakes,
              all_unobserved = unobserved_dead_snakes,
              all_effort = effort_record))
}  

# # Test
# test_total_quarters <- 5
# test_N <- 1000
# test_size_dist <- c(0.5, 0.1, 0.2, 0.2)
# p_g <- 0.75
# a <- quarter_operations(initial_N = test_N,
#                         initial_size_dist = test_size_dist,
#                         p_g = p_g,
#                         lambda = lambda,
#                         total_quarters =  test_total_quarters,
#                         total_days = total_days,
#                         erad = "on",
#                         erad_method = erad_methods)
# remove(test_total_quarter)
# remove(test_N)
# remove(test_size_dist)
# remove(p_g)
# remove(a)
# remove(total_days)


