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
    snake$SVL <- snake$SVL*exp(growth_coefs[1,] + growth_coefs[2,]*log(snake$SVL/700))
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
  if(nrow(females) > 0) {
    for(snake in 1:nrow(females)) {
      # Coin flip of whether reproduction occurs is based on both density dependence and the probability that the individual snake will be able to reproduce based on its size
      if (rbinom(1, 1, prob = density_prob*females$repro_prob[snake]) == 1) {
        females$mom_status[snake] <- 1
      } else {
        females$mom_status[snake] <- 0
      }
    }
    moms <- females[females$mom_status == 1,] %>%
      select(., -c("mom_status"))
  } else { # If there are no females in the population, then return an empty dataframe with the right columns 
    moms <- start_pop[0,]
  }
  
  return(moms)
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
  if(total_offspring > 0) {
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
  offspring$growth_quant <- sample(growth_quantiles, total_offspring, replace = T, 
                                   prob = c(0.25, 0.5, 0.25, 0.1, 0.05))
  } else { # If no offspring are produced, return an empty data frame with the right columns
    offspring <- mom_pop[0,]
  }
  
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
  # if(nrow(pop) == 0) {
  #   print("Population eradicated before natural mortality")
  #   break
  # }
  
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
  
  # 
  # if(nrow(pop_2) == 0) {
  #   stop("Population eradicated after natural mortality")
  # }
  
  return(list(surviving_pop = pop_2,
              num_dead_snakes = length(dead_snakes)))
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
  max_prob_value <- fdrtool::dhalfnorm(0, theta = ((sqrt(pi/2))/sd)) 
  # Calculate PDF value for the current density
  PDf_value <- fdrtool::dhalfnorm(current_density, theta = ((sqrt(pi/2))/sd))
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
  # Calculating size class population N and rounding to make sure its a whole number
  size_pop <- vector()
  for(k in 1:nrow(size_limits)) {
    size_pop[k] <- round(init_N*init_size_dist[k])
  }
  total_pop <- sum(size_pop)
  
  # Setting up initial population dataframe:
  init_pop <- as.data.frame(matrix(NA, nrow = total_pop, ncol = 5))
  colnames(init_pop) <- c("ID", "SVL", "sex", "repro_prob", "growth_quant")
  # Creating initial population characteristics
  init_pop$ID <- paste0("A", sample(1000:9000, total_pop))
  init_pop$SVL <- c(runif(size_pop[1], size_limits$lower[1], size_limits$upper[1]), 
                    runif(size_pop[2], size_limits$lower[2]+1, size_limits$upper[2]),
                    runif(size_pop[3], size_limits$lower[3]+1, size_limits$upper[3]), 
                    runif(size_pop[4], size_limits$lower[4]+1, size_limits$upper[4]))
  init_pop$sex <- sample(c("M", "F"), total_pop, replace = T)
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


### Function for size at transition from lizards to mammals, when snakes become 
### susceptible to the toxic baits - current mean and SD values are approximated from a 
### graph in 2019 SLITHER report (ask Shane for actual values)
## Inputs:
# - snake_SVL - individual snake's SVL
diet_trans_fun <- function(snake_SVL,
                           lower_limit = 700,
                           upper_limit = 950) {
  SD <- (upper_limit-lower_limit)/3.29
  mean <- lower_limit + (upper_limit-lower_limit)/2
  if(snake_SVL < lower_limit) {
    diet_trans <- 0
  } else if (snake_SVL >= lower_limit & snake_SVL < upper_limit) {
    diet_trans <- pnorm(snake_SVL, mean, SD)
  } else {
    diet_trans <- 1
  }
  return(diet_trans)
}

### Function for decay when large snakes are less interested or less likely to be killed by
### the dose in toxic bait methods (ADS and bait tubes) - for now, using a guassian decay
### form (probability = exp(-((SVL-1400)^2/d)), this may be changed to an exponential 
### decay form (probability = b^(SVL-1400)) after consultation 
mortality_decay_fun <- function(snake_SVL,
                                min_SVL = 1400,
                                d = 21000) {
  if(snake_SVL > 1400) {
    y <- exp(-((snake_SVL-min_SVL)^2/d))
  } else {
    y <- 1
  }
  return(y)
}

# # Test
# mortality_decay_fun(1500)
  

# # Plot bait susceptibility vs SVL 
# SVL <- c(350:1800)
# mortality_probability <- vector()
# for(i in 1:length(SVL)){
#   mortality_probability[i] <- diet_trans_fun(SVL[i])*mortality_decay_fun(SVL[i])
# }
# ggplot(data = as.data.frame(cbind(SVL, mortality_probability))) +
#   geom_line(aes(x = SVL, y = mortality_probability)) +
#   geom_vline(xintercept = 850, color = "#7CAE00") +
#   geom_vline(xintercept = 950, color = "#00BFC4") +
#   geom_vline(xintercept = 1150, color = "#C77CFF") +
#   annotate(geom = "text",
#            x = size_class_limits[-1,1]+20,
#            y = 0.1,
#            label = size_class_names[-1],
#            angle = 90) +
#   theme_bw()
  

### Function to perform all eradication methods 
## Inputs:
# - day_pop = starting population for that day 
# - coverage = the percentage of the total area that the eradication method will reach
# - size_affected = what size of snake is susceptible to the method (lower and upper 
#                   limit for now - might change into a probability based on size later
# - mortality_prob = mortality rate of method if encountered 

eradication_fun <- function(encounter_pop, 
                            #day_pop,
                            #coverage, 
                            #size_affected = c(0, Inf), 
                            mortality_prob,
                            method) {
  # # Separating out the proportion of the population that will be affected by eradication (coverage) from the non-encountering population
  # encounter_pop <- slice_sample(day_pop, prop = coverage)
  # no_encounter_pop <- day_pop[!day_pop$ID %in% encounter_pop$ID, ]
  # If snake is in the size range to be affected, perform Bernoulli draw to determine if it encounters the method (encounter = mortality for now)
  for(snake in 1:nrow(encounter_pop)) {
    if(method %in% erad_methods[c(1,4)]) {
      if (rbinom(1, 1, 
                 prob = mortality_prob[1]*diet_trans_fun(encounter_pop$SVL[snake])*mortality_decay_fun(encounter_pop$SVL[snake])) == 1) { 
        encounter_pop$encounter[snake] <- 1
      } else {
        encounter_pop$encounter[snake] <- 0
      }
    } else if(method %in% erad_methods[c(2,3)]) {
      size_class <- size_class_fun(encounter_pop$SVL[snake])
      if (rbinom(1, 1, prob = mortality_prob[size_class]) == 1) { 
        encounter_pop$encounter[snake] <- 1
      } else {
        encounter_pop$encounter[snake] <- 0
      }
    }  
    # if(encounter_pop$SVL[snake] > size_affected[1] & encounter_pop$SVL[snake] < size_affected[2]) {
    #   if (rbinom(1, 1, prob = mortality_prob) == 1) { 
    #     encounter_pop$encounter[snake] <- 1
    #   } else {
    #     encounter_pop$encounter[snake] <- 0
    #   }
    # } else {
    #   encounter_pop$encounter[snake] <- 0
    # }
  }
  # Separating dead snakes from alive snakes, and getting rid of unnecessary encounter column
  dead_snakes <- encounter_pop[encounter_pop$encounter == 1,] %>%
    select(, -c("encounter"))
  alive_snakes <- encounter_pop[encounter_pop$encounter == 0,] %>%
    select(, -c("encounter"))
  
  # # Recombining surviving snakes with  non-encounter snakes
  # day_pop <- rbind(select(alive_snakes, -c("encounter")), no_encounter_pop)
  # Return both dead and live snakes (eventually, only return live snakes)
  return(list(dead_snakes = dead_snakes, 
              alive_snakes = alive_snakes))
}

# # Testing
# first_encounter_pop <- slice_sample(first_day_pop, prop = coverage)
# e <- eradication_fun(encounter_pop = first_encounter_pop,
#                      # day_pop = first_day_pop, 
#                      # coverage = 0.5, 
#                      mortality_prob = mortality_prob_erad_methods[2],
#                      method = erad_methods[2])


# erad_timing_fun <- function(day, 
#                             quarter,
#                             erad_method,
#                             pop,
#                             ADS_IDs,
#                             transect_pop_IDs){
#   # Creating empty lists to hold dead snakes, alive snakes and effort lists for each method
#   observed <- list()
#   unobserved <- list()
#   transect_encounter_pop <- list()
#   transect_encounter_pop[[1]] <- transect_pop
#   #ADS_encounter_pop <- list()
#   effort <- list()
#   #################### below operations should be part of the daily_operations_fun ###
#   # # # Perform each eradication that occurs on this day & quarter sequentially
#   # # survivors[[1]] <- growing_pop
#   # # Selecting proportion of population near transects, where some eradication will 
#   # # occur
#   # transect_encounter_pop <- slice_sample(growing_pop, prop = erad_coverage$visual)
#   # # Selecting proportion of population available to ADS, based on how much overlap there
#   # # is with the transects (take the overlap proportion from transect_encounter_pop,
#   # # then the remaining from the rest of the population)
#   # num_ADS_encounter <- nrow(growing_pop)*erad_coverage$ADS
#   # overlap_encounter_pop <- slice_sample(transect_encounter_pop, 
#   #                                       prop = ADS_overlap_w_transect)
#   # remaining_encounter_pop <- slice_sample(growing_pop[!growing_pop$ID %in% transect_encounter_pop$ID, ],
#   #                                 n = (num_ADS_encounter-nrow(overlap_encounter_pop)))
#   # ADS_encounter_pop <- rbind(overlap_encounter_pop, remaining_encounter_pop)
#   # # Separating individuals who are not subject to any eradication in this quarter
#   # no_encounter_pop <- growing_pop[!growing_pop$ID %in% c(transect_encounter_pop$ID, 
#   #                                                        ADS_encounter_pop$ID), ]
#   ############################################################ 
#   for(method in 1:length(erad_method)) {
#     method_name <- erad_method[method]
#     #pop <- survivors[[method]]
#     if(quarter %in% erad_quarters[[method_name]] & day %in% erad_days[[method_name]]){
#       if(method_name == erad_methods[1]) {
#         erad_pop <- eradication_fun(encounter_pop = ADS_pop, 
#                                   #coverage = erad_coverage[[method_name]], 
#                                   mortality_prob = mortality_prob_erad_methods[[method_name]],
#                                   method = method_name)
#         ADS_surviving_pop <- erad_pop$alive_snakes
#       } else {
#         # if ("ADS" %in% erad_method) {
#         # transect_encounter_pop[[method]] <- transect_pop[!transect_pop$ID %in% erad_pop$dead_snakes$ID, ]
#       #   erad_pop <- eradication_fun(encounter_pop = transect_encounter_pop[[method]], 
#       #                               #coverage = erad_coverage[[method_name]], 
#       #                               mortality_prob = mortality_prob_erad_methods[[method_name]],
#       #                               method = method_name)
#       #   transect_encounter_pop[[method+1]] <- erad_pop$alive_snakes
#       # } else {
#         erad_pop <- eradication_fun(encounter_pop = transect_encounter_pop[[method]], 
#                                     #coverage = erad_coverage[[method_name]], 
#                                     mortality_prob = mortality_prob_erad_methods[[method_name]],
#                                     method = method_name)
#         transect_encounter_pop[[method+1]] <- erad_pop$alive_snakes
#         #}
#       }
#       #survivors[[method + 1]] <- erad_pop$alive_snakes
#       effort[[method_name]] <- as.data.frame(cbind("effort" = effort_erad_methods[[method_name]],
#                                                    "day" = day,
#                                                    "quarter" = quarter))
#   
#       # Eradication with bodies occurs:
#       if(nrow(erad_pop$dead_snakes) > 0) {
#         if (method_name %in% erad_methods[c(2:3)]) {
#           observed[[method_name]] <- erad_pop$dead_snakes
#           observed[[method_name]]$day <- day
#           observed[[method_name]]$quarter <- quarter
#         } else { # Eradication with no bodies
#           unobserved[[method_name]] <- erad_pop$dead_snakes
#           unobserved[[method_name]]$day <- day
#           unobserved[[method_name]]$quarter <- quarter
#         }
#       }
#     # else {
#     #   survivors[[method + 1]] <- survivors[[method]]
#     # }
#       }
#   }
#   
#   #final_transect_pop <- transect_encounter_pop[[length(transect_encounter_pop)]]
#   
#   if("ADS" %in% erad_method) {
#     final_survivors <-  rbind(ADS_surviving_pop, 
#                               transect_encounter_pop[[length(transect_encounter_pop)]]) %>%
#       .[!duplicated(.$ID),]
#   } else {
#     final_survivors <- transect_encounter_pop[[length(transect_encounter_pop)]]
#   }
#   #final_survivors <- rbind(survivors, no_encounter_pop)
#   return(list(observed_snakes = observed, 
#               unobserved_dead_snakes = unobserved,
#               surviving_snakes = final_survivors,
#               all_effort = effort))
# }

erad_timing_fun <- function(day, 
                            quarter,
                            erad_method,
                            pop,
                            ADS_IDs,
                            transect_IDs){
  # Creating empty lists to hold dead snakes, alive snakes and effort lists for each method
  observed <- list()
  unobserved <- list()
  transect_encounter_pop <- list()
  # ###### Subsetting by day by method
  # day_transect_IDs <- list()
  # if(day_method_subset <- "on") {
  #   for(method in 1:length(erad_method)
  #  day_transect_IDs[[method]] <-  sample(transect_IDs, 
  #                                        (erad_day_coverage[[method]]/erad_coverage$transects_per_quarter)*length(transect_IDs))
  # }
  ####### Subsetting by day, but with the same subset for all methods used
  # # Subset the susceptible populations by how much area can be covered in any given day by each method
  # day_transect_IDs <- sample(transect_IDs, 
  #                            (erad_coverage$bait_tube/erad_coverage$transects_per_quarter)*length(transect_IDs))
  # Separate out populations susceptible to transect & ADS methods 
  # transect_encounter_pop[[1]] <- pop[pop$ID %in% day_transect_IDs, ]
  
  ###### Original version, with quarter subset only
  transect_encounter_pop[[1]] <- pop[pop$ID %in% transect_IDs, ]
  
  ADS_encounter_pop <- pop[pop$ID %in% ADS_IDs, ]
  # print(paste0("ADS encounters ", nrow(ADS_encounter_pop)))
  effort <- list()
  for(method in 1:length(erad_method)) {
    method_name <- erad_method[method]
    #print(method_name)
    #pop <- survivors[[method]]
    #if(quarter %in% erad_quarters[[method_name]] & day %in% erad_days[[method_name]]){
      if(method_name == "ADS") {
        if(nrow(ADS_encounter_pop) > 0) {
        erad_pop <- eradication_fun(encounter_pop = ADS_encounter_pop, 
                                    mortality_prob = mortality_prob_erad_methods[[method_name]],
                                    method = method_name)
        ADS_surviving_pop <- erad_pop$alive_snakes
        overlap_dead_snake_IDs <- which(transect_encounter_pop[[method]]$ID %in% erad_pop$dead_snakes$ID)
        } else {
          overlap_dead_snake_IDs <- 0
        } 
        if(length(overlap_dead_snake_IDs) > 0) { 
          transect_encounter_pop[[method+1]] <- transect_encounter_pop[[method]][-overlap_dead_snake_IDs,]
        } else {
          transect_encounter_pop[[method+1]] <- transect_encounter_pop[[method]]
        }
      } else {
        if(nrow(transect_encounter_pop[[method]]) > 0){
        erad_pop <- eradication_fun(encounter_pop = transect_encounter_pop[[method]], 
                                    mortality_prob = mortality_prob_erad_methods[[method_name]],
                                    method = method_name)
        transect_encounter_pop[[method+1]] <- erad_pop$alive_snakes
        } else {
          transect_encounter_pop[[method+1]] <- transect_encounter_pop[[method]]
        }
      }
    # Tracking effort for each method  
    effort[[method_name]] <- as.data.frame(cbind("effort" = effort_erad_methods[[method_name]],
                                                   "day" = day,
                                                   "quarter" = quarter))
      
      # Eradication with bodies occurs:
      if(exists("erad_pop") == TRUE) {
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
        }
      }
  }
  
  # Combine the surviving snakes from the different methods (ignoring duplicates for now)
  if(exists("ADS_surviving_pop") == TRUE) {
    final_survivors <- rbind(ADS_surviving_pop, 
                              transect_encounter_pop[[length(transect_encounter_pop)]])
     # print(paste0("ADS survivors = ", nrow(ADS_surviving_pop)))
  } else {
    final_survivors <- rbind(ADS_encounter_pop,
                             transect_encounter_pop[[length(transect_encounter_pop)]])
    # print(paste0("ADS encounter pop = ", nrow(ADS_encounter_pop)))
  }
  # Delete duplicates from the final survivor dataframe
  final_survivors <- final_survivors[!duplicated(final_survivors$ID),]
  
  # Making sure that no dead snakes make it into the final survivors dataframe (if they were in ADS pop but died by another method)
  all_dead_snake_IDs <- c(observed$visual$ID, observed$trap$ID, 
                          unobserved$ADS$ID, unobserved$bait_tube$ID)
  relic_IDs <- which(final_survivors$ID %in% all_dead_snake_IDs)
  if(length(relic_IDs) > 0) {
    final_survivors <- final_survivors[-relic_IDs, ]
  }
  
  return(list(observed_snakes = observed, 
              unobserved_dead_snakes = unobserved,
              surviving_snakes = final_survivors,
              all_effort = effort,
              all_dead_IDs = all_dead_snake_IDs))
}


# Test
# day <- 7
# quarter <- 1
# growing_pop <- init_pop_fun(N, size_dist, size_class_limits)
# transect_pop <- slice_sample(growing_pop, prop = erad_coverage$transects_per_quarter)
# transect_pop_IDs <- transect_pop$ID
# ADS_overlap_w_transect <- 1
# total_ADS_prop <- nrow(growing_pop)*erad_coverage$ADS
# non_transect_ADS <- total_ADS_prop - nrow(slice_sample(transect_pop, prop = ADS_overlap_w_transect))
# ADS_pop <- rbind(slice_sample(growing_pop, n = non_transect_ADS),
#                  slice_sample(transect_pop, prop = ADS_overlap_w_transect))
# ADS_pop_IDs <- ADS_pop$ID
# 
# 
# x <- erad_timing_fun(day = day,
#                      quarter = quarter,
#                      erad_method = erad_methods[2:4],
#                      pop = growing_pop,
#                      ADS_IDs = ADS_pop_IDs,
#                      transect_IDs = transect_pop_IDs)
# 
# removals <- as.data.frame(matrix(NA, nrow = 50, ncol = 3))
# colnames(removals) <- erad_methods[2:4]
# 
# for(i in 1:50) {
#   x <- erad_timing_fun(day = day,
#                        quarter = quarter,
#                        erad_method = erad_methods[2:4],
#                        pop = growing_pop,
#                        ADS_IDs = ADS_pop_IDs,
#                        transect_IDs = transect_pop_IDs)
#   if(length(x$observed_snakes) > 0) {
#     if(length(x$observed_snakes$visual) > 0) {
#     removals$visual[i] <- nrow(x$observed_snakes$visual)
#     } else {
#       removals$visual[i] <- 0
#     }
#     if(length(x$observed_snakes$trap) > 0) {
#     removals$trap[i] <- nrow(x$observed_snakes$trap)
#     } else {
#       removals$trap[i] <- 0
#     }
#   } else {
#     removals$visual[i] <- 0
#     removals$trap[i] <- 0
#   }
#   if(length(x$unobserved_dead_snakes) > 0) {
#     removals$bait_tube[i] <- nrow(x$unobserved_dead_snakes$bait_tube)
#   } else {
#     removals$bait_tube[i] <- 0
#   }
# }
# summary(removals)
## Without daily subsetting, just using a quarter subset (23% of 2750 total pop)
# visual          trap         bait_tube   
# Min.   :0.00   Min.   :0.0   Min.   :0.00  
# 1st Qu.:1.00   1st Qu.:1.0   1st Qu.:0.25  
# Median :2.00   Median :1.5   Median :1.00  
# Mean   :2.02   Mean   :2.0   Mean   :1.30  
# 3rd Qu.:3.00   3rd Qu.:3.0   3rd Qu.:2.00  
# Max.   :5.00   Max.   :6.0   Max.   :4.00 

## With daily subsetting for visual survey (4 transects per day):
# visual          trap        bait_tube   
# Min.   :0.00   Min.   :0.00   Min.   :0.00  
# 1st Qu.:0.00   1st Qu.:0.00   1st Qu.:0.00  
# Median :0.00   Median :0.00   Median :0.00  
# Mean   :0.44   Mean   :0.48   Mean   :0.42  
# 3rd Qu.:1.00   3rd Qu.:1.00   3rd Qu.:1.00  
# Max.   :2.00   Max.   :2.00   Max.   :2.00


## With daily subsetting for traps (6 transects per day):
# visual          trap        bait_tube   
# Min.   :0.00   Min.   :0.00   Min.   :0.00  
# 1st Qu.:0.00   1st Qu.:0.00   1st Qu.:0.00  
# Median :1.00   Median :0.00   Median :0.00  
# Mean   :1.22   Mean   :0.58   Mean   :0.52  
# 3rd Qu.:2.00   3rd Qu.:1.00   3rd Qu.:1.00  
# Max.   :5.00   Max.   :2.00   Max.   :2.00

## With daily subsetting for bait tubes (12 transects per day)
# visual          trap        bait_tube   
# Min.   :0.00   Min.   :0.00   Min.   :0.00  
# 1st Qu.:0.00   1st Qu.:1.00   1st Qu.:0.00  
# Median :1.00   Median :1.00   Median :1.00  
# Mean   :1.28   Mean   :1.32   Mean   :0.92  
# 3rd Qu.:2.00   3rd Qu.:2.00   3rd Qu.:1.00  
# Max.   :4.00   Max.   :4.00   Max.   :5.00  

# #
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
                             quarter,
                             ADS_overlap_w_transect,
                             erad_coverage,
                             primary_sampling_period,
                             erad_quarters,
                             erad_days) {
  daily_pop <- list()
  erad_results <- list()
  ## Growth and mortality on a daily scale within the quarter 
  daily_pop[[1]] <- first_day_pop
  ## Segment IDs for individuals that will be subject to eradication if eradication is 
  # occurring this quarter
  if(erad == "on" & quarter %in% unique(unlist(erad_quarters))) {
    # Proportion vulnerable to transects
    transect_pop_IDs <- sample(daily_pop[[1]]$ID, 
                               nrow(daily_pop[[1]])*erad_coverage$transects_per_quarter)
    if(ADS_overlap_w_transect == 1) { # if transects are completed covered by ADS
    # Proportion vulnerable to ADS
    ADS_pop_IDs <- sample(daily_pop[[1]]$ID,
                          nrow(daily_pop[[1]])*erad_coverage$ADS)
    
    } else if(ADS_overlap_w_transect < 1 & ADS_overlap_w_transect > 0) { # if there is some ADS coverage of transects
    # Total proportion subject to ADS & how much of that overlaps with the transects
    overlap_pop_IDs <- sample(transect_pop_IDs, length(transect_pop_IDs)*ADS_overlap_w_transect)
    total_num_ADS <- nrow(daily_pop[[1]])*erad_coverage$ADS
    non_transect_ADS <- total_num_ADS - length(overlap_pop_IDs)
    # Proportion vulnerable to ADS
    ADS_pop_IDs <- c(sample(daily_pop[[1]]$ID[!daily_pop[[1]]$ID %in% transect_pop_IDs], non_transect_ADS),
                     overlap_pop_IDs)
    # print(paste0("ADS pop ID day 1 = ", length(ADS_pop_IDs)))
    } else if (ADS_overlap_w_transect == 0) { # if there is no overlap between transects & ADS
      # print(paste("transect", length(transect_pop_IDs)))
      # Proportion vulnerable to ADS
      ADS_pop_IDs <- sample(daily_pop[[1]]$ID[!daily_pop[[1]]$ID %in% transect_pop_IDs],
                            nrow(daily_pop[[1]])*erad_coverage$ADS)
    }
    # Distinguishing the population that won't encounter any eradication method
    no_encounter_pop_IDs <- daily_pop[[1]]$ID[!daily_pop[[1]]$ID %in% unique(c(transect_pop_IDs, ADS_pop_IDs))]
    # print(paste("no_encounter", length(no_encounter_pop_IDs)))
  }
  # print(paste0("transect_pop_IDs ", length(transect_pop_IDs)))
  # print(paste0("ADS_pop_IDs ", length(ADS_pop_IDs)))
  # print(paste0("no_encounter_pop_IDs ", length(no_encounter_pop_IDs)))
  # Loop for each day in 90 days
  # Sequence of events: growth, eradication if it occurs, natural mortality if eradication 
  # doesn't occur
  for(day in 1:(total_days-1)) {
    #  All snakes are exposed to natural mortality
    # print(paste0("daily pop: ", nrow(daily_pop[[day]])))
    pop <- daily_pop[[day]]
    if(nrow(pop) == 0) {
      print("Population eradicated at start of day")
      break
    }
    #print(paste0(day, "before erad"))
    ## If eradication occurring during this model run,
    if(erad == "on") {
      # and if eradication is occurring in the current quarter,
      if(quarter %in% unique(unlist(erad_quarters))) {
        # and if the current day is within the primary sampling period, then no natural 
        # mortality occurs, only individual size growth (not including reproducing females)
        if(day > primary_sampling_period[1] & day < primary_sampling_period[2]) {
          # Excluding reproductive females to create population who has the potential to grow in this quarter
          growing_pop <- pop[!pop$ID %in% moms$ID, ]
          ## Each snake grows (or not) based on it's individual growth quantile
          for(snake in 1:nrow(growing_pop)) {
            growing_pop[snake,] <- daily_gompertz_growth_fun(growing_pop[snake,],
                                                             p_g)
          }
          surviving_moms <- moms
        } else { # if the day is outside the primary sampling period, then
          # Natural mortality occurs on entire population
          mort_results <- daily_mortality_fun(pop, size_class_limits)
          #print(paste0("natural mort occurred - current pop ", nrow(mort_results$surviving_pop)))
          if(nrow(mort_results$surviving_pop) == 0) {
            print("Population eradicated from natural mortality")
            break
          }
          surviving_pop <- mort_results$surviving_pop
          # Checking to see what the mortality is from natural causes
          dead_snakes <- mort_results$num_dead_snakes
          #print(paste0("natural mortality is ", dead_snakes))
          # Excluding reproductive females to create population who has the potential to grow in this quarter
          growing_pop <- surviving_pop[!surviving_pop$ID %in% moms$ID, ]
          # Separating out surviving moms to add back in after growth
          surviving_moms <- surviving_pop[surviving_pop$ID %in% moms$ID, ]
          ## Each snake grows (or not) based on it's individual growth quantile
          #print(paste0("growing pop is ", nrow(growing_pop)))
          for(snake in 1:nrow(growing_pop)) {
            growing_pop[snake,] <- daily_gompertz_growth_fun(growing_pop[snake,],
                                                             p_g)
          }
        }
      } else { # If eradication isn't occurring in this quarter, then both natural 
        # mortality and growth occur
        # Natural mortality occurs on entire population
        mort_results <- daily_mortality_fun(pop, size_class_limits)
        surviving_pop <- mort_results$surviving_pop
        # Checking to see what the mortality is from natural causes
        dead_snakes <- mort_results$num_dead_snakes
        #print(paste0("natural mortality is ", dead_snakes))
        # Excluding reproductive females to create population who has the potential to grow in this quarter
        growing_pop <- surviving_pop[!surviving_pop$ID %in% moms$ID, ]
        # Separating out surviving moms to add back in after growth
        surviving_moms <- surviving_pop[surviving_pop$ID %in% moms$ID, ]
        ## Each snake grows (or not) based on it's individual growth quantile
        for(snake in 1:nrow(growing_pop)) {
          growing_pop[snake,] <- daily_gompertz_growth_fun(growing_pop[snake,],
                                                           p_g)
        }
      }
      # Separate out the methods that are occurring in this quarter and on this day
      if(quarter %in% unique(unlist(erad_quarters)) & day %in% unique(unlist(erad_days[[quarter]]))) {
        # Figure out which eradication methods occur on this day in this quarter
        quarter_methods <- vector()
        day_methods <- vector()
        for(method in methods) {
          if(quarter %in% erad_quarters[[method]]) {
            quarter_methods <- c(quarter_methods, method)
          }
          if(day %in% erad_days[[quarter]][[method]]) {
            day_methods <- c(day_methods, method)
          }
        }
        today_methods <- intersect(quarter_methods, day_methods)
      } else {
        today_methods <- vector()
      }
      # print(today_methods)
      if(length(today_methods) > 0) {
        # If any snakes died before the start of the primary sampling period of natural 
        # causes, remove from the ID lists (performed every day, but should only make 
        # any difference outside of the primary sampling period)
        if(nrow(daily_pop[[1]]) > nrow(growing_pop)) {
          dead_snake_IDs <- setdiff(daily_pop[[1]]$ID, growing_pop$ID)
          ADS_remove <- which(ADS_pop_IDs %in% dead_snake_IDs)
          transect_remove <- which(transect_pop_IDs %in% dead_snake_IDs)
          no_encounter_remove <- which(no_encounter_pop_IDs %in% dead_snake_IDs)
          if(length(ADS_remove) > 0) {
            ADS_pop_IDs <- ADS_pop_IDs[-ADS_remove]
          }
          if(length(transect_remove) > 0) {
            transect_pop_IDs <- transect_pop_IDs[-transect_remove]
          }
          if(length(no_encounter_remove) > 0) {
            no_encounter_pop_IDs <- no_encounter_pop_IDs[-no_encounter_remove]
          }
        }
        # If there are still susceptible snakes left for both transects & ADS
        if(length(transect_pop_IDs) > 0 & length(ADS_pop_IDs) > 0) {
        # With updated ID lists, perform all eradications occurring on this day 
        erad_results[[day]] <- erad_timing_fun(day = day, 
                                               quarter = quarter, 
                                               erad_method = today_methods, 
                                               pop = growing_pop,
                                               ADS_IDs = ADS_pop_IDs,
                                               transect_IDs = transect_pop_IDs)
        # print(paste0("ADS pop IDs before reconcile = ", length(ADS_pop_IDs)))
        # Update eradication population ID lists by removing IDs that were removed today
        if(length(erad_results[[day]]$all_dead_IDs) > 0) {
          ADS_remove_2 <- which(ADS_pop_IDs %in% erad_results[[day]]$all_dead_IDs)
          transect_remove_2 <- which(transect_pop_IDs %in% erad_results[[day]]$all_dead_IDs)
          #print(paste0("transect removals = ", length(transect_remove_2)))
          if(length(ADS_remove_2) > 0) {
            ADS_pop_IDs <- ADS_pop_IDs[-ADS_remove_2]
          }
          if(length(transect_remove_2) > 0) {
          # print(erad_results[[day]]$all_dead_IDs)
            transect_pop_IDs <- transect_pop_IDs[-transect_remove_2]
          }
        }
        # print(paste0("transect_pop_IDs after erad = ", length(transect_pop_IDs)))
        # print(paste0("ADS_pop_IDs after erad = ", length(ADS_pop_IDs)))
        
        # Surviving snakes are combined with surviving reproducers and those not subjected 
        # to eradication to populate the next day
        daily_pop[[day + 1]] <- rbind(erad_results[[day]]$surviving_snakes, surviving_moms,
                                      growing_pop[growing_pop$ID %in% no_encounter_pop_IDs, ])
        # print(paste0("No encounter = ",nrow(growing_pop[growing_pop$ID %in% no_encounter_pop_IDs, ])))
        # print(paste0("Survivors = ", nrow(erad_results[[day]]$surviving_snakes)))
        
        # If only transect susceptible snakes are still alive,
      } else if(length(transect_pop_IDs) > 0 & length(ADS_pop_IDs) == 0) {
        erad_results[[day]] <- erad_timing_fun(day = day, 
                                               quarter = quarter, 
                                               erad_method = today_methods, 
                                               pop = growing_pop,
                                               ADS_IDs = ADS_pop_IDs,
                                               transect_IDs = transect_pop_IDs)
        # Update eradication population ID lists by removing IDs that were removed today
        if(length(erad_results[[day]]$all_dead_IDs) > 0) {
          transect_remove_2 <- which(transect_pop_IDs %in% erad_results[[day]]$all_dead_IDs)
          # print(paste0("transect removals = ", length(transect_remove_2)))
          if(length(transect_remove_2) > 0) {
            # print(erad_results[[day]]$all_dead_IDs)
            transect_pop_IDs <- transect_pop_IDs[-transect_remove_2]
          }
        }
        # Surviving snakes are combined with surviving reproducers to populate the next day
        daily_pop[[day + 1]] <- rbind(erad_results[[day]]$surviving_snakes, surviving_moms,
                                      growing_pop[growing_pop$ID %in% no_encounter_pop_IDs, ])
        # print("daily_pop nrow", nrow(daily_pop[[day+1]]))
        
        # If only ADS susceptible snakes are still alive,
      } else if(length(transect_pop_IDs) == 0 & length(ADS_pop_IDs) > 0) {
        erad_results[[day]] <- erad_timing_fun(day = day, 
                                               quarter = quarter, 
                                               erad_method = today_methods, 
                                               pop = growing_pop,
                                               ADS_IDs = ADS_pop_IDs,
                                               transect_IDs = transect_pop_IDs)
        # Update ADS IDs list
        if(length(erad_results[[day]]$all_dead_IDs) > 0) {
          ADS_remove_2 <- which(ADS_pop_IDs %in% erad_results[[day]]$all_dead_IDs)
          if(length(ADS_remove_2) > 0) {
            ADS_pop_IDs <- ADS_pop_IDs[-ADS_remove_2]
          }
        }
        # Surviving snakes are combined with surviving reproducers and non-susceptible 
        # snakes to populate the next day
        daily_pop[[day + 1]] <- rbind(erad_results[[day]]$surviving_snakes, surviving_moms,
                                      growing_pop[growing_pop$ID %in% no_encounter_pop_IDs, ])
      } else {
        # If no susceptible snakes are still alive, combine non-susceptible snakes and 
        # surviving reproducers to populate the next day
        if(nrow(rbind(growing_pop, surviving_moms)) == 0) {
          print("Population eradicated before next day")
          break
        }
        daily_pop[[day + 1]] <- rbind(growing_pop, surviving_moms)
      }
        }  else { # If no eradication is occurring today, combine non-susceptible snakes and 
          # surviving reproducers to populate the next day
          if(nrow(rbind(growing_pop, surviving_moms)) == 0) {
            print("Population eradicated before next day")
            break
          }
        daily_pop[[day + 1]] <- rbind(growing_pop, surviving_moms)
      }
    } else { # If no eradication is occurring at all, then natural mortality occurs on 
      # entire population every day
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
      daily_pop[[day + 1]] <- rbind(growing_pop, surviving_moms)
    }
    
    # # print(paste0("Daily pop = ", nrow(daily_pop[[day+1]])))
    # if(nrow(rbind(growing_pop, surviving_moms)) == 0) {
    #   print("Population eradicated")
    #   break
    # }
  }
   
  names(daily_pop) <- paste0("day_", c(1:length(daily_pop)))
  if(length(erad_results) > 0) {
    names(erad_results) <- paste0("day_", c(1:length(erad_results)))
  }
  
  ## Extracting the population from the day before any eradication occurs and after all eradication is done
  if(erad == "on") {
    observed_quarter_methods <- names(which(sapply(erad_quarters[c(2:3)], function(x) quarter %in% x)))
    if(length(observed_quarter_methods) > 0) { 
      if(min(unlist(erad_days[[quarter]][observed_quarter_methods])) > 1) {
        day_before_erad <- min(unlist(erad_days[[quarter]][observed_quarter_methods])) - 1
      } else {
        day_before_erad <- min(unlist(erad_days[[quarter]][observed_quarter_methods]))
      }
      if((max(unlist(erad_days[[quarter]][observed_quarter_methods])) + 1) < 91) {
        day_after_erad <- max(unlist(erad_days[[quarter]][observed_quarter_methods])) + 1
      } else {
        day_after_erad <- max(unlist(erad_days[[quarter]][observed_quarter_methods]))
      }
      day_before_erad_pop <- daily_pop[[day_before_erad]]
      day_after_erad_pop <- daily_pop[[day_after_erad]]
    } else {
      day_before_erad_pop <- NA
      day_after_erad_pop <- NA
    }
  } else {
    day_before_erad_pop <- NA
    day_after_erad_pop <- NA
  } 

  
  return(list(daily_timeseries = daily_pop,
              all_erad_results = erad_results,
              day_before_erad_pop = day_before_erad_pop,
              day_after_erad_pop = day_after_erad_pop))
}

# Test
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
#                       quarter = 1,
#                       ADS_overlap_w_transect = 1)
# remove(first_day_pop)
# remove(moms)
# remove(z)


# Reformatting eradication results
erad_results_format <- function(erad_results) {
  # Extracting dataframes from list that have removals
  erad_melt <- erad_results[unlist(lapply(erad_results, length) != 0)]
  # Setting up dataframes to hold unobserved & observed dead snakes, as well as the eradication effort put in
  quarter_unobserved <- as.data.frame(matrix(nrow = 0, ncol = 8))
  colnames(quarter_unobserved) <- c("ID", "SVL", "sex", "repro_prob", "growth_quant",
                                    "day", "quarter", "L1") 
  quarter_observed <- as.data.frame(matrix(nrow = 0, ncol = 8))
  colnames(quarter_observed) <- c("ID", "SVL", "sex", "repro_prob", "growth_quant",
                                    "day", "quarter", "L1") 
  quarter_effort <- as.data.frame(matrix(nrow = 0, ncol = 4))
  colnames(quarter_effort) <- c("effort", "day", "quarter", "L1")
  for(day in 1:length(erad_melt)){
    # print(paste0("day ", day))
    # print(paste0("unobserved dead ", length(erad_melt[[day]]$unobserved_dead_snakes)))
    # print(paste0("observed dead ", length(erad_melt[[day]]$observed_dead_snakes)))
    # Melting unobserved (dead) snakes into one dataframe
    if (length(erad_melt[[day]]$unobserved_dead_snakes) > 0) {
      # colnames(quarter_unobserved) <- c(colnames(erad_melt[[day]]$unobserved_dead_snakes[[1]]), "L1")
      x <- melt(erad_melt[[day]]$unobserved_dead_snakes, 
                id.vars = colnames(erad_melt[[day]]$unobserved_dead_snakes[[1]]))
      # print(str(x))
      quarter_unobserved <- rbind(quarter_unobserved, 
                                  melt(erad_melt[[day]]$unobserved_dead_snakes, 
                                       id.vars = colnames(erad_melt[[day]]$unobserved_dead_snakes[[1]])))
    } else {
      quarter_unobserved <- quarter_unobserved
    }
    # Melting observed snakes into one dataframe
    if (length(erad_melt[[day]]$observed_snakes) > 0) {
      #colnames(quarter_observed) <- c(colnames(erad_melt[[day]]$observed_snakes[[1]]), "L1")
      x2 <- melt(erad_melt[[day]]$observed_snakes, 
                id.vars = colnames(erad_melt[[day]]$observed_snakes[[1]]))
      # print(str(x2))
      quarter_observed <- rbind(quarter_observed, 
                                  melt(erad_melt[[day]]$observed_snakes, 
                                       id.vars = colnames(erad_melt[[day]]$observed_snakes[[1]])))
    } else {
      quarter_observed <- quarter_observed
    }
    # Melting effort into one dataframe
    if (length(erad_melt[[day]]$all_effort) > 0) {
      #colnames(quarter_effort) <- c(colnames(erad_melt[[day]]$all_effort[[1]]), "L1")
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

# # Test
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
                               erad_method,
                               erad_coverage,
                               primary_sampling_period,
                               erad_quarters,
                               erad_days) {
  # Empty lists to put population from each time step in (records the population at the 
  # beginning of each quarter, and for all the days in the final quarter)
  quarter_timeseries_pop <- list()
  daily_results <- list()
  # Setting up empty dataframe to record IDs for reproducing females in each quarter
  repro_females_timeseries <- list()
  # Empty list to hold the effort that will be expended if eradication occurs (summed in 
  # each quarter):
  effort_record <- list()
  # Creating a list of lists to keep track of "observed" snakes from visual surveys & traps
  # and "unobserved" but still dead snakes from ADS and bait tubes  
  observed_snakes <- list()
  unobserved_dead_snakes <- list()
  
  # Empty lists to save the pop from the day before and after all eradications in each quarter
  quarter_pops_before_erad <- list()
  quarter_pops_after_erad <- list()
  
  # Setting up the initial population in the first time step
  quarter_timeseries_pop[[1]] <- init_pop_fun(init_N = initial_N,
                                              init_size_dist = initial_size_dist,
                                              size_limits = size_class_limits)
  # print("initial_pop created")
  for(quarter in 1:total_quarters) {
    tic(paste0("Quarter ", quarter))
    # Check if the population has gone to 0
    if(nrow(quarter_timeseries_pop[[quarter]]) == 0) {
      print("Population eradicated")
      break
    } 
    # Setting up empty list to hold any observed and unobserved dead snakes in this quarter, 
    # if observation occurs, as well as the effort expended for each method
    observed_snakes[[quarter]] <- list()
    unobserved_dead_snakes[[quarter]] <- list()
    effort_record[[quarter]] <- list()
    # Update density dependent parameter & r_density for reproduction
    current_N <- nrow(quarter_timeseries_pop[[quarter]])
    r_density_prob <- DD_param_fun(K, current_N)
    
    # Isolate reproducing females for this quarter, if any exist 
    if(length(repro_females_timeseries) > 0) { # This condition should only be relevant for the 1st quarter
      if(nrow(repro_females_timeseries[[quarter - 1]]) > 0) {
        # Isolate females that reproduced last quarter, if there are any
        recent_mom_IDs <- repro_females_timeseries[[quarter - 1]]$ID
        # Choose which females will reproduce
        moms <- repro_females_fun(start_pop = quarter_timeseries_pop[[quarter]][!quarter_timeseries_pop[[quarter]]$ID %in% recent_mom_IDs,],
                                  density_prob = r_density_prob)
      } else {
        moms <- repro_females_fun(start_pop = quarter_timeseries_pop[[quarter]],
                                  density_prob = r_density_prob)
      }
    } else {
        moms <- repro_females_fun(start_pop = quarter_timeseries_pop[[quarter]],
                                  density_prob = r_density_prob)
    }
     #print(paste0("moms ", nrow(moms)))
    # If there are females reproducing this quarter, then create offspring (if not, don't)
    if(nrow(moms) > 0) {
      # Producing new offspring to join population at the end of the quarter
      offspring <- gen_offspring_fun(mom_pop = moms,
                                     time_step = quarter,
                                     offspring_lambda = lambda)
      #print(str(offspring))
      # Keeping track of which females reproduced in which quarters (next step will be to use     this to exclude females who've reproduced in the last 2 quarters)
      repro_females_timeseries[[quarter]] <- moms
    } else { # Create empty dataframes with the right column names for offspring and reproducing females lists if no reproducing females in this quarter
      offspring <- quarter_timeseries_pop[[1]][0,]
      repro_females_timeseries[[quarter]] <- quarter_timeseries_pop[[1]][0,]
    }
     #print(paste0("Offspring = ", nrow(offspring)))
    # print(paste0("Moms = ", nrow(repro_females_timeseries[[quarter]])))
    # Running daily operations - natural mortality and individual growth
    daily_results <- daily_operations(first_day_pop = quarter_timeseries_pop[[quarter]],
                                      p_g = p_g,
                                      moms = moms,
                                      total_days = total_days,
                                      methods = erad_method,  
                                      erad = erad, 
                                      quarter = quarter,
                                      ADS_overlap_w_transect = ADS_overlap_on_transect,
                                      erad_coverage = erad_coverage,
                                      primary_sampling_period = primary_sampling_period,
                                      erad_quarters = erad_quarters,
                                      erad_days = erad_days)
    daily_timeseries_pop <- daily_results$daily_timeseries
    print(paste0("daily timeseries indiv ", nrow(daily_timeseries_pop[[1]])))
    if(length(daily_timeseries_pop) < total_days) {
      print("Population eradication")
      break
    }
    # print(str(daily_timeseries_pop[[total_days]]))
    if(length(daily_results$all_erad_results) > 0) {
    # Filling the observed, unobserved dead snakes and effort lists:
    erad_reformatted <- erad_results_format(erad_results = daily_results$all_erad_results)
    observed_snakes[[quarter]] <- erad_reformatted$observed
    unobserved_dead_snakes[[quarter]] <- erad_reformatted$unobserved
    effort_record[[quarter]] <- erad_reformatted$effort
    } 
    # Adding offspring and the reproducing females back into the surviving and grown general   population to be the start of the next quarter
    quarter_timeseries_pop[[quarter + 1]] <- rbind(daily_timeseries_pop[[total_days]], offspring)
    # print(str(quarter_timeseries_pop))
    # Updating sexual maturity status after the quarter's worth of growth
    for(snake in 1:nrow(quarter_timeseries_pop[[quarter + 1]])) {
      quarter_timeseries_pop[[quarter + 1]]$repro_prob[snake] <- maturity_fun(quarter_timeseries_pop[[quarter + 1]]$SVL[snake])
    }
    # Saving before and after eradication populations
    quarter_pops_before_erad[[quarter]] <- daily_results$day_before_erad_pop
    quarter_pops_after_erad[[quarter]] <- daily_results$day_after_erad_pop
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
              all_effort = effort_record,
              all_quarters_before_after_erad = list(pop_before_erad = quarter_pops_before_erad,
                                                    pop_after_erad = quarter_pops_after_erad)))

}  

# # Test
# test_total_quarters <- 3
# test_N <- 1000
# test_size_dist <- c(0.5, 0.1, 0.2, 0.2)
# p_g <- 0.75
# a <- quarter_operations(initial_N = test_N,
#                         initial_size_dist = test_size_dist,
#                         p_g = p_g,
#                         lambda = lambda,
#                         total_quarters =  test_total_quarters,
#                         total_days = 91,
#                         erad = "on",
#                         erad_method = erad_methods)
# remove(test_total_quarters)
# remove(test_N)
# remove(test_size_dist)
# remove(p_g)
# remove(a)
# remove(total_days)


