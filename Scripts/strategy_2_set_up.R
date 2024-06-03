## Strategy 2 set up

### Manually set IBM parameters
# Number of quarters to generate - 10 years
erad_quarter_time_step <- 2

# Methods for all evenatualities; starting, threshold 1, threshold 2 and threshold 3
method_option_names <- c("initial", paste0("threshold_", c(1:3)))
method_options <- list()
for(option in 1:length(method_option_names)) {
  method_options[[option]] <- list()
}
names(method_options) <- method_option_names

# Identifying the methods that will be used under each of these conditions
method_options$initial$methods <- erad_methods
method_options$threshold_1$methods <- erad_methods[2]
method_options$threshold_2$methods <- erad_methods
method_options$threshold_3$methods <- erad_methods[2]

## Quarters where eradication methods are used for each condition
# All methods used in all quarters for first 3 conditions
for(option in 1:3) {
  method_options[[option]]$erad_quarters <- list()
  for(method in method_options[[option]]$methods) {
    method_options[[option]]$erad_quarters[[method]] <- c(1:erad_quarter_time_step)
  }
}
# Condition 4 only has methods occurring in 2nd quarter
method_options$threshold_3$erad_quarters <- 2

## Days where eradication methods are used for each condition
# Initial condition days
method_options$initial$erad_days <- list()
for(quarter in 1:erad_quarter_time_step) {
  method_options$initial$erad_days[[quarter]] <- list()
  # ADS: low effort (1 treatment) per quarter
  method_options$initial$erad_days[[quarter]]$ADS <- c(45, 48)
  # Visual survey: medium effort (6 weeks, 2 teams) per quarter - surveying every other day
  method_options$initial$erad_days[[quarter]]$visual <- seq(2, (7*6 - 1), 2)
  # Trap: medium effort (6 weeks) per quarter - checking traps every 3 days
  method_options$initial$erad_days[[quarter]]$trap <- seq(2, (7*6 - 1), 3)
  # Bait tubes: medium effort (6 weeks) per quarter
  method_options$initial$erad_days[[quarter]]$bait_tube <- seq(7*6, (7*12 - 1), 3)
}
names(method_options$initial$erad_days) <- paste0("quarter_", c(1:erad_quarter_time_step))

# Threshold 1 days (visual surveys only)
method_options$threshold_1$erad_days <- list()
for(quarter in 1:erad_quarter_time_step) {
  method_options$threshold_1$erad_days[[quarter]] <- list()
  # Visual survey: medium effort (6 weeks, 2 teams) per quarter - surveying every other day
  method_options$threshold_1$erad_days[[quarter]]$visual <- seq(2, (7*6 - 1), 2)
}
names(method_options$threshold_1$erad_days) <- paste0("quarter_", c(1:erad_quarter_time_step))

# Threshold 2 days (same as initial - coverage is what expands)
method_options$threshold_2$erad_days <- method_options$initial$erad_days

# Threshold 3 days (no effort in first quarter, with visual effort in the 2nd quarter)
method_options$threshold_3$erad_days <- list()
for(quarter in 1:erad_quarter_time_step) {
  method_options$threshold_3$erad_days[[quarter]] <- list()
  method_options$threshold_1$erad_days[[quarter]]$visual <- seq(2, (7*6 - 1), 2)
}
names(method_options$threshold_3$erad_days) <- paste0("quarter_", c(1:erad_quarter_time_step))



# Bounds of primary sampling periods for each condition (the same for the first 3)
for(option in 1:length(method_options)) {
  method_options[[option]]$primary_sampling_period <- c(2,(7*6 - 1))
}



## Coverage for each method in a quarter for each condition
for(option in 1:length(method_options)) {
  method_options[[option]]$erad_coverage <- list()
  # overlap of ADS over transects (100%) is the same for all conditions
  method_options[[option]]$ADS_overlap_on_transect <- 1 
}
# Initial condition coverage
method_options$initial$erad_coverage$ADS <- 1
# Transect coverage (~50%)
method_options$initial$erad_coverage$transects_per_quarter <- 0.5

# Threshold 1 condition coverage
method_options$threshold_1$erad_coverage$ADS <- 0 # necessary for a function, can fix this later
# Transect coverage (100% because effort is doubled) 
method_options$threshold_1$erad_coverage$transects_per_quarter <- 1

# Threshold 2 condition coverage
method_options$threshold_2$erad_coverage$ADS <- 1
# Transect coverage (100%)
method_options$threshold_2$erad_coverage$transects_per_quarter <- 1

# Threshold 3 condition coverage
method_options$threshold_3$erad_coverage$ADS <- 0 # necessary for a function,
# Transect coverage (~50%)
method_options$threshold_3$erad_coverage$transects_per_quarter <- 0.5

# Number of visual survey teams for each method
for(option in 1:length(method_options)) {
  method_options[[option]]$num_teams <- list()
}

method_options$initial$num_teams$visual <- 1
method_options$threshold_1$num_teams$visual <- 1
method_options$threshold_2$num_teams$visual <- 1
method_options$threshold_3$num_teams$visual <- 1

## Function to evaluate the thresholds for strategy 2
strat_2_threshold_fun <- function(mean_N_df,
                                  upper_3_threshold = 3,
                                  total_threshold = 0.01,
                                  small_increase_threshold = 0.1,
                                  area = area_size,
                                  size_class = size_class_names) {
  # Condition (default is initial)
  condition <- "initial"
  # Adding a column with density
  mean_N_df$mean_density <- mean_N_df$N/area
  # Identifying final time step (and the one right before)
  last_quarters <- tail(sort(unique(mean_N_df$Quarter)),2)
  # Separating final time step
  last_quarter_means <- mean_N_df[mean_N_df$Quarter == last_quarters[2], ]
  
  ## Calculating mean densities 
  # Upper 3 size classes combined
  upper_3_mean_density <- sum(last_quarter_means$N[last_quarter_means$size_class != size_class[1]])/area
  # Total population
  total_mean_density <- sum(last_quarter_means$N)
  
  # Calculating density increase between second to last and last time step for small snakes
  small_increase <- mean_N_df$mean_density[mean_N_df$size_class == size_class[1] & mean_N_df$Quarter == last_quarters[2]]/mean_N_df$mean_density[mean_N_df$size_class == size_class[1] & mean_N_df$Quarter == last_quarters[1]] - 1
  
  
  ## Threshold 1: Estimated mean density is of upper size classes is <= 3 snake/ha combined but total population is greater than > 0.01
  if(upper_3_mean_density <= upper_3_threshold & total_mean_density > total_threshold) {
    condition <- "threshold_1"
  }  
  
  ## Threshold 2: Estimated mean small size class density is increasing by >=10% compared to previous time step
  if (upper_3_mean_density >= upper_3_threshold & total_mean_density > total_threshold & small_increase >= small_increase_threshold) {
    condition <- "threshold_2"
  }
  
  ## Threshold 3: Estimated mean density of entire population is =< 0.01 snakes/ha for at least 2 quarters
  if(total_mean_density <= total_threshold) {
    condition <- "threshold_3"
  }
  
  return(condition)
}