## Strategy 1 set up

# Number of quarters to generate - 10 years
erad_quarter_time_step <- 4*10

# Quarters where eradication methods are used, and which days in that quarter:
erad_quarters <- list()
erad_quarters$ADS <- c(1:erad_quarter_time_step)
erad_days <- list()
for(quarter in 1:erad_quarter_time_step) {
  erad_days[[quarter]] <- list()
  erad_days[[quarter]]$ADS <- c(45, 48)
}
names(erad_days) <- paste0("quarter_", erad_quarters$ADS)

# Bounds of primary sampling period (since there is no monitoring, this doesn't matter but it's a required parameter for a function)
primary_sampling_period <- c(0,0)

# Coverage for each method in a quarter 
erad_coverage <- list()
erad_coverage$ADS <- 1
# Transect coverage and overlap of ADS over transects - there aren't actually any transect methods so neither of these is relevant, but they are required parameters at the moment and its not a high priority to fix that, so keeping them in for now
erad_coverage$transects_per_quarter <- 0
ADS_overlap_on_transect <- 1 

# Another required parameter that doesn't matter for this strategy, because there is no visual survey occurring
num_teams <- list()
num_teams$visual <- 0

