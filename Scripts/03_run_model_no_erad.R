###### Running full model without eradication ################

library(reshape2)
library(ggplot2)
library(dplyr)
library(here)
library(tictoc)
library(fdrtool)

## Loading objects and functions
source(here("Scripts/00_user_inputs.R"))
source(here("Scripts/01_model_functions.R"))
source(here("Scripts/02_results_functions.R"))


# Parameters that may be changed or subject to sensitivity analysis or to simulate over
# Set study/eradication area size, which dictates both K and the initial population
area_size <- 10

# Set carrying capacity for this population:
K <- 119*area_size

# Set up initial population size:
init_N <- 100*area_size

# Initial size distribution:
init_size_dist <- c("small" = 0.1,
                    "medium" = 0.05,
                    "large" = 0.05,
                    "xlarge" = 0.8)

# Two time steps will be used, one for the number of quarters, the other for number of days in the quarter
quarter_time_step <- 80
day_time_step <- 91
# No eradication occurring, so:
erad_quarters <- 0
# Growth probability (p_g)
g_density_prob <- 0.75 

 
# Running model 
quarter_results <- quarter_operations(initial_N = init_N, 
                                      initial_size_dist = init_size_dist, 
                                      p_g = g_density_prob,
                                      lambda = lambda,
                                      total_quarters = quarter_time_step,
                                      total_days = day_time_step,
                                      erad = "off",
                                      erad_method = erad_methods)



# Plotting the total population in each quarter
plot_1 <- ggplot(quarter_results$all_quarters, aes(x = Quarter, fill = size_category)) +
  geom_bar() +
  geom_hline(yintercept = K) +
  theme_bw()

# Model validation stats and plots
model_val <- model_validation_fun(quarter_timeseries = quarter_results$quarter_timeseries,
                          all_quarters = quarter_results$all_quarters,
                          repro_females = quarter_results$all_repro_females,
                          total_quarters = quarter_time_step,
                          size_classes = rownames(size_class_limits))
model_val$model_val_summary
model_val$model_validation$LL_snakes_prop
model_val$model_validation$LL_plot

# Final size class distributions:
all_pop <- nrow(quarter_results$quarter_timeseries[[quarter_time_step+1]])
small_prop <-  nrow(quarter_results$quarter_timeseries[[quarter_time_step+1]][quarter_results$quarter_timeseries[[quarter_time_step+1]]$SVL < 850, ])/all_pop
medium_prop <- nrow(quarter_results$quarter_timeseries[[quarter_time_step+1]][quarter_results$quarter_timeseries[[quarter_time_step+1]]$SVL > 850 & quarter_results$quarter_timeseries[[quarter_time_step+1]]$SVL < 950, ])/all_pop
large_prop <- nrow(quarter_results$quarter_timeseries[[quarter_time_step+1]][quarter_results$quarter_timeseries[[quarter_time_step+1]]$SVL > 950 & quarter_results$quarter_timeseries[[quarter_time_step+1]]$SVL < 1150, ])/all_pop
xlarge_prop <-  nrow(quarter_results$quarter_timeseries[[quarter_time_step+1]][quarter_results$quarter_timeseries[[quarter_time_step+1]]$SVL > 1150, ])/all_pop

### When starting pop is N = 50*area_size & init_size_dist = 0.25 for all sizes, after 20 years:
# Small prop: 0.4126
# Medium prop: 0.0567
# Large prop: 0.1431
# Xlarge prop: 0.3877


### When starting pop is 50*area_size, init_size_dist = above dists, after 20 years:
# Small prop: 0.4121
# Medium prop: 0.0475
# Large prop: 0.1502
# Xlarge prop: 0.3902

## starting pop = 100*area_size, init_size_dist = 0.25 for all, after 20 years:
# Small prop: 0.3872
# Medium prop: 0.0465
# Large prop: 0.1473
# Xlarge prop: 0.4190


## starting pop = 100*area_size, init_size_dist = 0.8 small, 0.05 medium and large, 0.1 xlarge, after 20 years:
# Small prop: 0.4002
# Medium prop: 0.0526
# Large prop: 0.1530
# Xlarge prop: 0.3942

## starting pop = 100*area_size, init_size_dist = 0.1 small, 0.05 medium and large, 0.8 xlarge, after 20 years:
# Small prop: 0.3906
# Medium prop: 0.0522
# Large prop: 0.1504
# Xlarge prop: 0.4069


## so this seems like a pretty stable distribution

