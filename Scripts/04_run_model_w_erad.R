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
area_size <- 20

# Set carrying capacity for this population:
K <- 119*area_size


# # # Set up initial population size based on previous model runs (may change this later):
# qrt_ts_5.11.23 <- readRDS(here("Data/quarter_timeseries_5.11.23_DD.rds"))
# last_quarter <- qrt_ts_5.11.23[[length(qrt_ts_5.11.23)]]
# N <- nrow(last_quarter)
# 
# # Initial size distribution, based on previous model runs
# size_dist <- c("small" = nrow(last_quarter[last_quarter$SVL <= size_class_limits[1, 2],])/N,
#                 "medium" = nrow(last_quarter[last_quarter$SVL > size_class_limits[2, 1] & last_quarter$SVL <= size_class_limits[2, 2],])/N,
#                 "large" = nrow(last_quarter[last_quarter$SVL > size_class_limits[3, 1] & last_quarter$SVL <= size_class_limits[3, 2],])/N,
#                 "xlarge" = nrow(last_quarter[last_quarter$SVL > size_class_limits[4, 1],])/N)

## For initial population roughly based on 10 ha runs in the past:
N <- 100*area_size
size_dist <- c(0.6, 0.1, 0.1, 0.2)


# Growth probability (p_g)
g_density_prob <- 0.75 

# Number of quarters to generate
erad_quarter_time_step <- 8
day_time_step <- 91

# Quarters where eradication methods are used, and which days in that quarter:
erad_quarters <- list()
erad_quarters$ADS <- c(2, 3, 4, 5, 6, 7)
erad_quarters$visual <- c(2, 3, 4, 5, 6, 7)
erad_quarters$trap <- c(3, 6)
erad_quarters$bait_tube <- c(1,8)
erad_days <- list()
erad_days$ADS <- c(7,14,21)
erad_days$visual <- c(28:62)
erad_days$trap <- c(28:34)
erad_days$bait_tube <- c(28:34)

# Coverage for each method (totally arbitrary, adjust later)
erad_coverage <- list()
erad_coverage$ADS <- 0.5
erad_coverage$visual <- 0.5
erad_coverage$trap <- 0.5
erad_coverage$bait_tube <- 0.5
# Overlap of ADS over transects
ADS_overlap_on_transect <- 0.75 # i.e. total area coverage, by at least one method 



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
  geom_hline(yintercept = K) +
  theme_bw() +
  scale_x_continuous(breaks = unique(erad_quarter_results$all_quarters$Quarter), labels = unique(erad_quarter_results$all_quarters$Quarter))

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

# Plotting observed removed snakes (captured through trapping and visual survey)
observed_list <- erad_quarter_results$all_observed[unlist(lapply(erad_quarter_results$all_observed, length) != 0)]
all_observed <- bind_rows(observed_list)
for(snake in 1:nrow(all_observed)){
  all_observed$size_category[snake] <- size_class_fun(all_observed$SVL[snake], 
                                                      size_class_limits)
}
# Factoring size category column for graphing
all_observed$size_category <- factor(all_observed$size_category, 
                                     levels = c("small", "medium", "large", "xlarge"))

# Scatter plot, SVL size on the y-axis, separated by quarter
observed_plot_1 <- ggplot(all_observed, aes(x = day, y = SVL, color = method)) +
  geom_jitter() +
  theme_bw() +
  facet_wrap("quarter")
# Bar plot
observed_plot_2 <- ggplot(all_observed, aes(fill = size_category, x = quarter)) +
  geom_bar(position = "dodge") +
  facet_wrap("method") +
  theme_bw()

# Plotting unobserved removed snakes (killed through ADS and bait tubes)
unobserved_list <- erad_quarter_results$all_unobserved[unlist(lapply(erad_quarter_results$all_unobserved, length) != 0)]
all_unobserved <- bind_rows(unobserved_list[c(1:3)])
for(snake in 1:nrow(all_unobserved)){
  all_unobserved$size_category[snake] <- size_class_fun(all_unobserved$SVL[snake], 
                                                      size_class_limits)
}
# Factoring size category column for graphing
all_unobserved$size_category <- factor(all_unobserved$size_category, 
                                     levels = c("small", "medium", "large", "xlarge"))

unobserved_plot_1 <- ggplot(all_unobserved , aes(x = day, y = SVL, color = method)) +
  geom_jitter() +
  theme_bw() +
  facet_wrap("quarter")


all_dead_snakes <- bind_rows(all_observed, all_unobserved)
all_dead_plot_1 <- ggplot(all_dead_snakes , aes(x = day, y = SVL, color = method)) +
  geom_jitter() +
  theme_bw() +
  facet_wrap("quarter")


### Plot of each method through time in a representative quarter (3) so we can look at
### how the method performs

ggplot(all_dead_snakes[all_dead_snakes$quarter == 3,]) +
  geom_bar(aes(x = day, fill = size_category)) +
  facet_wrap("method", scales = "free_y") +
  theme_bw() +
  labs(title = "Quarter 3")


## Proportion of population eradicated in each quarter when eradication occurred
erad_prop_per_quarter <- vector()
for(quarter in sort(unique(unlist(erad_quarters)))) {
  erad_prop_per_quarter[quarter] <- sum(nrow(erad_quarter_results$all_observed[[quarter]]), nrow(erad_quarter_results$all_unobserved[[quarter]]))/nrow(erad_quarter_results$quarter_timeseries[[quarter]])
}
erad_prop <- as.data.frame(cbind(sort(unique(unlist(erad_quarters))), 
                   erad_prop_per_quarter[ !is.na(erad_prop_per_quarter)]))
colnames(erad_prop) <- c("Quarter", "Erad_proportion")

ggplot(erad_prop) +
  geom_col(aes(y = Erad_proportion, x = Quarter)) +
  theme_bw()




### Encounter probability 
# Because the encounter probability has multiple pieces - the encounter probability per size, and then
# coverage (if coverage isn't 100% - if it is, then it's only the encounter probability per size

true_encounter_prob <- erad_coverage$visual*mortality_prob_erad_methods$visual







