### Strategy 3 set up 

# Number of quarters to generate - 10 years
erad_quarter_time_step <- 4

# Methods for all evenatualities; starting, threshold 1, threshold 2 and threshold 3
method_option_names <- c("initial", "threshold_1")
method_options <- list()
for(option in 1:length(method_option_names)) {
  method_options[[option]] <- list()
}
names(method_options) <- method_option_names

# Identifying the methods that will be used under each of these conditions
method_options$initial$methods <- erad_methods[c(1:2)]
method_options$threshold_1$methods <- erad_methods[2]


## Quarters where eradication methods are used for each condition
# All methods used in all quarters 
for(option in 1:length(method_options)) {
  method_options[[option]]$erad_quarters <- list()
  for(method in method_options[[option]]$methods) {
    method_options[[option]]$erad_quarters[[method]] <- c(1:erad_quarter_time_step)
  }
}


## Days where eradication methods are used for each condition
# Initial condition days
method_options$initial$erad_days <- list()
for(quarter in 1:erad_quarter_time_step) {
  method_options$initial$erad_days[[quarter]] <- list()
  # ADS: low effort (1 treatment) per quarter
  method_options$initial$erad_days[[quarter]]$ADS <- c(45, 48)
  # Visual survey: medium effort (6 weeks, 2 teams) per quarter - surveying every other day
  method_options$initial$erad_days[[quarter]]$visual <- c(7:21)
}
names(method_options$initial$erad_days) <- paste0("quarter_", c(1:erad_quarter_time_step))

# Threshold 1 days (visual surveys only)
method_options$threshold_1$erad_days <- list()
for(quarter in 1:erad_quarter_time_step) {
  method_options$threshold_1$erad_days[[quarter]] <- list()
  # Visual survey: medium effort (6 weeks, 2 teams) per quarter - surveying every other day
  method_options$threshold_1$erad_days[[quarter]]$visual <- c(7:21)
}
names(method_options$threshold_1$erad_days) <- paste0("quarter_", c(1:erad_quarter_time_step))


# Bounds of primary sampling periods for each condition (the same for the first 3)
for(option in 1:length(method_options)) {
  method_options[[option]]$primary_sampling_period <- c(7, 21)
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
method_options$initial$erad_coverage$transects_per_quarter <- 1

# Threshold 1 condition coverage
method_options$threshold_1$erad_coverage$ADS <- 0 # necessary for a function, can fix this later
# Transect coverage (100% because effort is doubled) 
method_options$threshold_1$erad_coverage$transects_per_quarter <- 1

# Number of visual survey teams for each method
for(option in 1:length(method_options)) {
  method_options[[option]]$num_teams <- list()
  method_options[[option]]$cost_num_teams <- list()
}
# This one is used to calculate the spatial coverage 
# (there might be more than one team, but if they don't overlap spatially, then the encounter
# probability will be the same as if its one team)
method_options$initial$num_teams$visual <- 1
method_options$threshold_1$num_teams$visual <- 1
# This one is used to calculate the cost, so its the actual number of teams per quarter
method_options$initial$cost_num_teams$visual <- 2
method_options$threshold_1$cost_num_teams$visual <- 2



# Function to evaluate if threshold conditions are met
strat_3_threshold_fun <- function(mean_N_df,
                                  upper_3_threshold = 1,
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
  
  ## Threshold 1: Estimated mean density is of upper size classes is <= 1 snake/ha combined 
  if(upper_3_mean_density <= upper_3_threshold) {
    condition <- "threshold_1"
  }  
  
  return(condition)
}