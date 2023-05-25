############## Run full model with eradication ################

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


# # Set up initial population size based on previous model runs (may change this later):
# last_quarter <- quarter_results$quarter_timeseries[[quarter_time_step]]
# N <- nrow(last_quarter)
# 
# # Initial size distribution, based on previous model runs
# size_dist <- c("small" = nrow(last_quarter[last_quarter$SVL <= size_class_limits[1, 2],])/N,
#                 "medium" = nrow(last_quarter[last_quarter$SVL > size_class_limits[2, 1] & last_quarter$SVL <= size_class_limits[2, 2],])/N,
#                 "large" = nrow(last_quarter[last_quarter$SVL > size_class_limits[3, 1] & last_quarter$SVL <= size_class_limits[3, 2],])/N,
#                 "xlarge" = nrow(last_quarter[last_quarter$SVL > size_class_limits[4, 1],])/N)

# For initial population roughly based on 10 ha runs in the past:
N <- 1000
size_dist <- c(0.57, 0.09, 0.15, 0.19)


# Growth probability (p_g)
g_density_prob <- 0.75 

# Number of quarters to generate
erad_quarter_time_step <- 8
day_time_step <- 91

# Quarters where eradication methods are used, and which days in that quarter:
erad_quarters <- list()
erad_quarters$ADS <- c(3, 6)
erad_quarters$visual <- c(2, 3, 6, 7)
erad_quarters$trap <- c(3, 6)
erad_quarters$bait_tube <- c(2, 3, 6)
erad_days <- list()
erad_days$ADS <- c(7,14,21)
erad_days$visual <- c(28:34, 56:62)
erad_days$trap <- c(56:62)
erad_days$bait_tube <- c(7, 28:34)


# Running model 
erad_quarter_results <- quarter_operations(initial_N = N, 
                                           initial_size_dist = size_dist, 
                                           p_g = g_density_prob,
                                           lambda = lambda,
                                           total_quarters = erad_quarter_time_step,
                                           total_days = day_time_step,
                                           erad = "on",
                                           erad_method = erad_methods)


# Plotting the total population in each quarter
erad_plot_1 <- ggplot(erad_quarter_results$all_quarters, 
                      aes(x = Quarter, fill = size_category)) +
  geom_bar() +
  geom_hline(yintercept = K)

# Model validation stats and plots
erad_model_val <- model_validation_fun(quarter_timeseries = erad_quarter_results$quarter_timeseries,
                                  all_quarters = erad_quarter_results$all_quarters,
                                  repro_females = erad_quarter_results$all_repro_females,
                                  total_quarters = erad_quarter_time_step,
                                  size_classes = rownames(size_class_limits))
erad_model_val$model_val_summary

# Plotting the observed, unobserved dead snakes and effort
effort_list <- erad_quarter_results$all_effort[unlist(lapply(erad_quarter_results$all_effort, length) != 0)]
all_effort <- bind_rows(effort_list)

effort_plot_1 <- ggplot(all_effort, aes(x = quarter, y = effort, fill = method)) +
  geom_col() +
  theme_bw()


observed_list <- erad_quarter_results$all_observed[unlist(lapply(erad_quarter_results$all_observed, length) != 0)]
all_observed <- bind_rows(observed_list)
for(snake in 1:nrow(all_observed)){
  all_observed$size_category[snake] <- size_class_fun(all_observed$SVL[snake], 
                                                      size_class_limits)
}
# Factoring size category column for graphing
all_observed$size_category <- factor(all_observed$size_category, 
                                     levels = c("small", "medium", "large", "xlarge"))

observed_plot_1 <- ggplot(all_observed , aes(x = quarter, y = SVL, color = method)) +
  geom_point() +
  theme_bw()


unobserved_list <- erad_quarter_results$all_unobserved[unlist(lapply(erad_quarter_results$all_unobserved, length) != 0)]
all_unobserved <- bind_rows(unobserved_list[c(1:3)])
for(snake in 1:nrow(all_unobserved)){
  all_unobserved$size_category[snake] <- size_class_fun(all_unobserved$SVL[snake], 
                                                      size_class_limits)
}
# Factoring size category column for graphing
all_unobserved$size_category <- factor(all_unobserved$size_category, 
                                     levels = c("small", "medium", "large", "xlarge"))

unobserved_plot_1 <- ggplot(all_unobserved , aes(x = quarter, y = SVL, color = method)) +
  geom_point() +
  theme_bw()


all_dead_snakes <- bind_rows(all_observed, all_unobserved)
all_dead_plot_1 <- ggplot(all_dead_snakes , aes(x = quarter, y = SVL, color = method)) +
  geom_jitter() +
  theme_bw()
