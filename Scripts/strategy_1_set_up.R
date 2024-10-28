## Strategy 1 set up


# The number of quarters in the strategy (methods are the same for every quarter, so just going to who one)
erad_quarter_time_step <- 40

method_options <- list()
method_options$initial <- list()

# Methods used
method_options$initial$methods <- erad_methods[1]

# Quarters where eradication methods are used, and which days in that quarter:
method_options$initial$erad_quarters <- list()
method_options$initial$erad_quarters$ADS <- c(1:erad_quarter_time_step)
method_options$initial$erad_days <- list()
for(quarter in 1:erad_quarter_time_step) {
  method_options$initial$erad_days[[quarter]] <- list()
  method_options$initial$erad_days[[quarter]]$ADS <- c(45, 48)
}
names(method_options$initial$erad_days) <- paste0("quarter_", method_options$initial$erad_quarters$ADS)

# Bounds of primary sampling period (since there is no monitoring, this doesn't matter but it's a required parameter for a function)
method_options$initial$primary_sampling_period <- c(0,0)

# Coverage for each method in a quarter 
method_options$initial$erad_coverage <- list()
method_options$initial$erad_coverage$ADS <- 1
# Transect coverage and overlap of ADS over transects - there aren't actually any transect methods so neither of these is relevant, but they are required parameters at the moment and its not a high priority to fix that, so keeping them in for now
method_options$initial$erad_coverage$transects_per_quarter <- 0
method_options$initial$ADS_overlap_on_transect <- 1 

# Another required parameter that doesn't matter for this strategy, because there is no visual survey occurring
method_options$initial$num_teams <- list()
method_options$initial$num_teams$visual <- 0

