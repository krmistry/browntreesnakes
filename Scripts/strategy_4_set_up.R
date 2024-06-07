

### Manually set IBM parameters
# Number of quarters to generate between estimation rounds - 1 year at a time
erad_quarter_time_step <- 4

# Methods for all evenatualities; starting, threshold 1, threshold 2 and threshold 3
method_option_names <- c("initial", "threshold_1")
method_options <- list()
for(option in 1:length(method_option_names)) {
  method_options[[option]] <- list()
}
names(method_options) <- method_option_names

# Identifying the methods that will be used under each of these conditions
method_options$initial$methods <- erad_methods[c(2:4)]
method_options$threshold_1$methods <- erad_methods[2]

## Quarters where eradication methods are used for each condition
# All methods used in all quarters for initial condition
method_options$initial$erad_quarters <- list()
for(method in method_options$initial$methods) {
  method_options$initial$erad_quarters[[method]] <- c(1:erad_quarter_time_step)
}

# Threshold 1: no methods used in the first 2 quarters, visual in later 2 quarters
method_options$threshold_1$erad_quarters <- list()
method_options$threshold_1$erad_quarters$visual <- c(3:4)

## Days where eradication methods are used for each condition
# Initial condition days
method_options$initial$erad_days <- list()
for(quarter in 1:erad_quarter_time_step) {
  method_options$initial$erad_days[[quarter]] <- list()
  # Visual survey: medium effort (9 weeks, 2 teams) per quarter - surveying every other day
  method_options$initial$erad_days[[quarter]]$visual <- seq(2, (7*10 - 1), 2)
  # Trap: medium effort (6 weeks) per quarter - checking traps every 3 days
  method_options$initial$erad_days[[quarter]]$trap <- seq(2, (7*10 - 1), 3)
  # Bait tubes: medium effort (6 weeks) per quarter
  method_options$initial$erad_days[[quarter]]$bait_tube <- seq(7*3, (7*12 - 1), 3)
}
names(method_options$initial$erad_days) <- paste0("quarter_", c(1:erad_quarter_time_step))

# Threshold 1 days (visual surveys only)
method_options$threshold_1$erad_days <- list()
for(quarter in 1:erad_quarter_time_step) {
  method_options$threshold_1$erad_days[[quarter]] <- list()
}
for(quarter in method_options$threshold_1$erad_quarters$visual) {
  # Visual survey: medium effort (3 weeks, 4 teams) per quarter - surveying every other day
  method_options$threshold_1$erad_days[[quarter]]$visual <- seq(2,(7*3 - 1), 2)
}
names(method_options$threshold_1$erad_days) <- paste0("quarter_", c(1:erad_quarter_time_step))


# Bounds of primary sampling periods for each condition (different for each)
for (condition in 1:length(method_option_names)){
  method_options[[condition]]$primary_sampling_period <- c(2,(7*3 - 1))
}

#method_options$threshold_1$primary_sampling_period <- c(2, 7*2)


## Coverage for each method in a quarter for each condition
for(option in 1:length(method_options)) {
  method_options[[option]]$erad_coverage <- list()
  # even though there is no ADS in this strategy, this is a required parameter, and it doesn't actually impact anything
  method_options[[option]]$ADS_overlap_on_transect <- 1 
}
# Initial condition coverage
method_options$initial$erad_coverage$ADS <- 0 # necessary for a function, can fix this later
# Transect coverage (100%)
method_options$initial$erad_coverage$transects_per_quarter <- 1

# Threshold 1 condition coverage
method_options$threshold_1$erad_coverage$ADS <- 0 # necessary for a function, can fix this later
# Transect coverage (100% because effort is doubled) 
method_options$threshold_1$erad_coverage$transects_per_quarter <- 1


# Number of visual survey teams for each method
for(option in 1:length(method_options)) {
  method_options[[option]]$num_teams <- list()
}

method_options$initial$num_teams$visual <- 1 # 2 teams when calculating cost
method_options$threshold_1$num_teams$visual <- 1 # 3 teams when calculating cost





## Function setting threshold to turn ADS on and off (if x-large snakes reach 0, turn off, 
## otherwise keep on)
strat_4_threshold_fun <- function(mean_N_df,
                                  xlarge_density = 0,
                                  area = area_size,
                                  size_class = size_class_names) {
  # Condition (default is initial)
  condition <- "initial"
  # Adding a column with density
  mean_N_df$mean_density <- mean_N_df$N/area
  # Identifying final time step 
  last_quarter <- tail(sort(unique(mean_N_df$Quarter)),1)
  # Separating final time step
  last_quarter_means <- mean_N_df[mean_N_df$Quarter == last_quarter, ]
  
  ## Calculating mean densities 
  # X-large size class
  xlarge_mean_density <- last_quarter_means$N[last_quarter_means$size_class == size_class[4]]/area
  
  
  ## Threshold 1: Estimated mean density is of x-large size class is <=0.01 snake/ha 
  if(xlarge_mean_density <= xlarge_density) {
    condition <- "threshold_1"
  } 
  
  return(condition)
}