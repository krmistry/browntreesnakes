### Exploring the influence of stochasticity on BTS generating model; this script runs 
### through scenarios with different starting population and starting size distributions
### without adding stochasticity to each time step's population calculations


library(fdrtool)
library(reshape2)
library(ggplot2)
library(tidyr)
library(here)


##### Vital rates function - use for both growth and reproductive rates

rates_fun <- function(carrying_capacity, 
                      peak_vital_rate,
                      current_density) {
  # Calculate sigma for scaling standard deviation using carrying capacity
  sd <- (carrying_capacity*0.975)/1.96
  # Calculate PDF when density is 0 (highest density probability value)
  max_prob_value <- dhalfnorm(0, theta = ((sqrt(pi/2))/sd))
  # Calculate PDF value for the current density
  PDf_value <- dhalfnorm(current_density, theta = ((sqrt(pi/2))/sd))
  # Use maximum value to scale current density's PDF value, and multiply by the peak vital rate
  rate <- (PDf_value/max_prob_value)*peak_vital_rate
  return(rate)
}

### Function to generate population for all size classes in next time step 
generate_fun <- function(t_minus_1,# The previous time step data (1 row, 4 columns for each size class)
                         growth, # growth rate vector for 3 smaller size classes
                         mortality, # mortality rate vector for all size classes
                         reproduction) { # reproductive rate vector for 3 larger size classes
  # Create vector to hold the current time step's population numbers
  t <- vector() 
  # Calculate small class population (if population is less than 1, set to 0)
  t[1] <- (1-(growth[1] + mortality[1]))*t_minus_1[, 1] + 
    reproduction[1]*t_minus_1[, 2] + 
    reproduction[2]*t_minus_1[, 3] +
    reproduction[3]*t_minus_1[, 4]
  if(t[1] >= 1) {
    t[1] <- t[1]
  } else { 
    t[1] <- 0
  }
  # Calculate medium class population (if population is less than 1, set to 0)
  t[2] <- (1-(growth[2] + mortality[2]))*t_minus_1[, 2] +
    growth[1]*t_minus_1[, 1]
  if(t[2] >= 1) {
    t[2] <- t[2]
  } else {
    t[2] <- 0
  }
  # Calculate large class population (if population is less than 1, set to 0)
  t[3] <- (1-(growth[3] + mortality[3]))*t_minus_1[, 3] +
    growth[2]*t_minus_1[, 2]
  if(t[3] >= 1) {
    t[3] <- t[3]
  } else {
    t[3] <- 0
  }
  # Calculate x-large class population (if population is less than 1, set to 0)
  t[4] <- (1 - mortality[4])*t_minus_1[, 4] +
    growth[3]*t_minus_1[, 3]
  if(t[4] >= 1) {
    t[4] <- t[4]
  } else {
    t[4] <- 0
  }
  
  return(t)
}



#########################################################################################################
############################## Starting Population Versions
##########################################################################################

# Starting population list to iterate over 
N_options <- unique(c(seq(100, 1100, 200), seq(500, 7500, 750))) # 15 scenarios
# Vector of names for lists of results by starting population
SP_list_names <- paste0("SP_", N_options)

##### Parameters that won't change
# Number of years (with 4 time steps per year)
Year <- 100
# Loading vital rates parameters
rates_list <- readRDS("Data/vital_rates_lit.rds")
# Separating out vital rates into separate vectors:
growth_peak_rates <- rates_list$growth_rates
repro_vector <- rates_list$repro_rates
mortality_vector <- rates_list$mortality_rates
# Size classes
bins <- c("small",
          "medium",
          "large",
          "xlarge")
n_bins <- length(bins)
# # Standard deviation of stochastic process - not needed for these versions
# sd <- 5

### Lists to store vital rates and population by size as they change through time
all_vital_rates <- list()
all_timeseries <- list()


### Loops
## Within each loop: set starting population, create dataframes to hold growth and reproductive rates, and the population numbers through time. Then, run model for Year*4 time series with each startig population value, and save the population results and the growth & reproduction values through time

for(sp in 1:length(N_options)) {
  # Set starting population
  N <- N_options[sp]
  # Create list of empty dataframes to hold vital rates
  vital_rates <- list()
  vital_rates$growth <- as.data.frame(matrix(NA, nrow = Year*4, ncol = 3))
  colnames(vital_rates$growth) <- bins[-4]
  vital_rates$repro <- as.data.frame(matrix(NA, nrow = Year*4, ncol = 3))
  colnames(vital_rates$repro) <- bins[-1]
  # Create empty dataframe for population timeseries
  timeseries <- as.data.frame(matrix(NA, nrow = Year*4, ncol = n_bins+1))
  colnames(timeseries) <- c(bins, "total")
  # Time 0 in first row
  timeseries[1, ] <- c(rep(N/4, 4), N)
  
  for(ts in 2:(Year*4)) {
    # Generate growth rate for this year
    vital_rates$growth[ts-1, ] <- rates_fun(carrying_capacity = 6000,
                                            peak_vital_rate = growth_peak_rates,
                                            current_density = timeseries$total[ts-1])
    vital_rates$repro[ts-1, ] <- rates_fun(carrying_capacity = 6000,
                                           peak_vital_rate = repro_vector,
                                           current_density = timeseries$total[ts-1])
    # Generate next time step populations
    timeseries[ts, 1:4] <- unlist(generate_fun(t_minus_1 = timeseries[ts-1,],
                                               growth = vital_rates$growth[ts-1, ],
                                               mortality = mortality_vector,
                                               reproduction = vital_rates$repro[ts-1, ])) 
    #+ rnorm(4, 0, sd) - excluding stochasticity for these versions
    timeseries$total[ts] <- sum(timeseries[ts, 1:4])
  }
  # Adding year and quarter columns
  timeseries$year <- rep(1:Year, each = 4)
  timeseries$quarter <- rep(paste0("Q", 1:4), Year)
  
  # Saving timeseries to list
  all_timeseries[[sp]] <- timeseries
  # Saving vital rate dataframes
  all_vital_rates[[sp]] <- vital_rates
  
}
names(all_timeseries) <- SP_list_names
names(all_vital_rates) <- SP_list_names

#### Creating and saving plots of growth rate, reproductive rate and populations through time

for(sp in SP_list_names) {
  # Restructuring growth rate results for plotting
  all_vital_rates[[sp]]$growth$year <- rep(1:Year, each = 4)
  growth_rates_plots <- melt(all_vital_rates[[sp]]$growth[-(Year*4),], id.vars = c("year"))
  colnames(growth_rates_plots)[2:3] <- c("size", "growth_rate")
  # Creating growth rate plot
  growth_plot <- ggplot(growth_rates_plots, aes(x = year, y = growth_rate)) +
    geom_line()+
    facet_wrap(vars(size)) +
    theme_bw() +
    labs(title = sp)
  #### Saving plots to file
  ggsave(filename = paste0(sp, "_growth.png"), plot = growth_plot,
         path = here("Results/sensitivity_analyses/no_stochasticity/starting_pop/Plots/Growth"))
  
  # Restructuring reproductive rate results for plotting
  all_vital_rates[[sp]]$repro$year <- rep(1:Year, each = 4)
  repro_rates_plots <- melt(all_vital_rates[[sp]]$repro[-(Year*4),], id.vars = c("year"))
  colnames(repro_rates_plots)[2:3] <- c("size", "repro_rate")
  
  reproductive_plot <- ggplot(repro_rates_plots, aes(x = year, y = repro_rate)) +
    geom_line()+
    facet_wrap(vars(size)) +
    theme_bw() +
    labs(title = sp)
  #### Saving plots to file
  ggsave(filename = paste0(sp, "_repro.png"), plot = reproductive_plot,
         path = here("Results/sensitivity_analyses/no_stochasticity/starting_pop/Plots/Reproduction"))
  
  ## Plotting proportion of size classes in population through timeseries
  # restructuring timeseries data for plotting, and adding timestep
  tm_plot <- melt(all_timeseries[[sp]][, -5], id.vars = c("year", "quarter"))
  colnames(tm_plot)[3:4] <- c("size", "pop")
  tm_plot$time_step <- rep(1:(Year*4), 4)
  # Creating population timeseries plot
  population_plot <- ggplot(tm_plot, aes(fill = size, y = pop, x = time_step)) +
    geom_bar(position="stack", stat="identity") +
    theme_bw() +
    labs(title = sp)
  #### Saving plots to file
  ggsave(filename = paste0(sp, "_pop.png"), plot = population_plot,
         path = here("Results/sensitivity_analyses/no_stochasticity/starting_pop/Plots/Population"))
}

### Plots comparing the final population and vital rate values from each starting population scenario
# Pulling out final population numbers
final_pop <- as.data.frame(matrix(NA, nrow = length(N_options), ncol = ncol(all_timeseries[[1]])))
colnames(final_pop) <- colnames(all_timeseries[[1]])

for(sp in 1:length(SP_list_names)) {
  final_pop[sp, ] <- all_timeseries[[sp]][Year*4, ]
}
final_pop$scenario <- SP_list_names
# Melting data for plotting
final_pop_melt <- melt(final_pop[, -c(5:7)], id.vars = "scenario")
colnames(final_pop_melt)[2:3] <- c("size", "pop")
# Making scenarios into ordered factors for clarity
final_pop_melt$scenario <- factor(final_pop_melt$scenario, levels = SP_list_names)

final_ts_plot <- ggplot(final_pop_melt, aes(y = pop, x = scenario, fill = size)) +
  geom_bar(position="stack", stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))

# Saving plot
ggsave(filename = "final_ts_pop_comparison.png", plot = final_ts_plot, 
       path = here("Results/sensitivity_analyses/no_stochasticity/starting_pop/Plots"))

# Pulling out final growth rate numbers
final_growth <- as.data.frame(matrix(NA, nrow = length(N_options), ncol = ncol(all_vital_rates[[1]]$growth)))
colnames(final_growth) <- colnames(all_vital_rates[[1]]$growth)

for(sp in 1:length(SP_list_names)) {
  final_growth[sp, ] <- all_vital_rates[[sp]]$growth[(Year*4-1), ]
}
final_growth$scenario <- SP_list_names
# Melting data for plotting
final_growth_melt <- melt(final_growth, id.vars = c("scenario"))
colnames(final_growth_melt)[2:3] <- c("size", "growth_rate")
# Making scenarios into ordered factors for clarity
final_growth_melt$scenario <- factor(final_growth_melt$scenario, levels = SP_list_names)
# Creating growth rate plot
growth_plot <- ggplot(final_growth_melt, aes(x = size, y = growth_rate, fill = scenario)) +
  geom_bar(position="dodge", stat="identity") +
  theme_bw()

# Saving plot
ggsave(filename = "final_ts_growth_comparison.png", plot = growth_plot, 
       path = here("Results/sensitivity_analyses/no_stochasticity/starting_pop/Plots"))

# Pulling out final reproductive rate numbers
final_repro <- as.data.frame(matrix(NA, nrow = length(N_options), ncol = ncol(all_vital_rates[[1]]$repro)))
colnames(final_repro) <- colnames(all_vital_rates[[1]]$repro)

for(sp in 1:length(SP_list_names)) {
  final_repro[sp, ] <- all_vital_rates[[sp]]$repro[(Year*4-1), ]
}
final_repro$scenario <- SP_list_names

# Melting data for plotting
final_repro_melt <- melt(final_repro, id.vars = c("scenario"))
colnames(final_repro_melt)[2:3] <- c("size", "growth_rate")
# Making scenarios into ordered factors for clarity
final_repro_melt$scenario <- factor(final_repro_melt$scenario, levels = SP_list_names)
# Creating growth rate plot
repro_plot <- ggplot(final_repro_melt, aes(x = size, y = growth_rate, fill = scenario)) +
  geom_bar(position="dodge", stat="identity") +
  theme_bw()

# Saving plot
ggsave(filename = "final_ts_repro_comparison.png", plot = repro_plot, 
       path = here("Results/sensitivity_analyses/no_stochasticity/starting_pop/Plots"))

#########################################################################################################
############################## Size Distribution Versions
##########################################################################################

# Starting size distribution
init_size_dist_options <- list()
# One size class starts with 0
init_size_dist_options[[1]] <- c(0, rep(1/3, 3))
init_size_dist_options[[2]] <- c(1/3, 0, rep(1/3, 2))
init_size_dist_options[[3]] <- c(rep(1/3, 2), 0, 1/3)
init_size_dist_options[[4]] <- c(rep(1/3, 3), 0)
# One small size class, the other 3 equally sized
init_size_dist_options[[5]] <- c(0.2, rep((1-0.2)/3, 3))
init_size_dist_options[[6]] <- c((1-0.2)/3, 0.2, rep((1-0.2)/3, 2))
init_size_dist_options[[7]] <- c(rep((1-0.2)/3, 2), 0.2, (1-0.2)/3)
init_size_dist_options[[8]] <- c(rep((1-0.2)/3, 3), 0.2)
# One large size class, the other 3 equally sized
init_size_dist_options[[9]] <- c(0.8, rep((1-0.8)/3, 3))
init_size_dist_options[[10]] <- c((1-0.8)/3, 0.8, rep((1-0.8)/3, 2))
init_size_dist_options[[11]] <- c(rep((1-0.8)/3, 2), 0.8, (1-0.8)/3)
init_size_dist_options[[12]] <- c(rep((1-0.8)/3, 3), 0.8)
# Two small size classes, two large size classes
init_size_dist_options[[13]] <- c(rep(0.8/2, 2), rep((1-0.8)/2, 2))
init_size_dist_options[[14]] <- c(0.8/2, (1-0.8)/2, 0.8/2, (1-0.8)/2)
init_size_dist_options[[15]] <- c(0.8/2, (1-0.8)/2, (1-0.8)/2, 0.8/2)
init_size_dist_options[[16]] <- c((1-0.8)/2, 0.8/2, (1-0.8)/2, 0.8/2)
init_size_dist_options[[17]] <- c((1-0.8)/2, 0.8/2, 0.8/2, (1-0.8)/2)
init_size_dist_options[[18]] <- c(rep((1-0.8)/2, 2), rep(0.8/2, 2))

# Vector of names for lists of results by starting population
SSD_list_names <- vector()
for(i in 1:length(init_size_dist_options)) {
  SSD_list_names[i] <- paste0("SSD_", paste0(round(init_size_dist_options[[i]], 2), sep = "_", collapse = "" ))
}


##### Parameters that won't change
# Starting population 
N <- 1000
# Number of years (with 4 time steps per year)
Year <- 100
# Loading vital rates parameters
rates_list <- readRDS("Data/vital_rates_lit.rds")
# Separating out vital rates into separate vectors:
growth_peak_rates <- rates_list$growth_rates
repro_vector <- rates_list$repro_rates
mortality_vector <- rates_list$mortality_rates
# Size classes
bins <- c("small",
          "medium",
          "large",
          "xlarge")
n_bins <- length(bins)
# # Standard deviation of stochastic process - not needed for these versions
# sd <- 5

### Lists to store vital rates and population by size as they change through time
all_vital_rates <- list()
all_timeseries <- list()


### Loops
## Within each loop: set starting population, create dataframes to hold growth and reproductive rates, and the population numbers through time. 

for(ssd in 1:length(init_size_dist_options)) {
  # Set starting size distribution
  init_size_dist <- init_size_dist_options[[ssd]]
  # Create list of empty dataframes to hold vital rates
  vital_rates <- list()
  vital_rates$growth <- as.data.frame(matrix(NA, nrow = Year*4, ncol = 3))
  colnames(vital_rates$growth) <- bins[-4]
  vital_rates$repro <- as.data.frame(matrix(NA, nrow = Year*4, ncol = 3))
  colnames(vital_rates$repro) <- bins[-1]
  # Create empty dataframe for population timeseries
  timeseries <- as.data.frame(matrix(NA, nrow = Year*4, ncol = n_bins+1))
  colnames(timeseries) <- c(bins, "total")
  # Time 0 in first row
  timeseries[1, ] <- c(N*init_size_dist, N)
  
  for(ts in 2:(Year*4)) {
    # Generate growth rate for this year
    vital_rates$growth[ts-1, ] <- rates_fun(carrying_capacity = 6000,
                                            peak_vital_rate = growth_peak_rates,
                                            current_density = timeseries$total[ts-1])
    vital_rates$repro[ts-1, ] <- rates_fun(carrying_capacity = 6000,
                                           peak_vital_rate = repro_vector,
                                           current_density = timeseries$total[ts-1])
    # Generate next time step populations
    timeseries[ts, 1:4] <- unlist(generate_fun(t_minus_1 = timeseries[ts-1,],
                                               growth = vital_rates$growth[ts-1, ],
                                               mortality = mortality_vector,
                                               reproduction = vital_rates$repro[ts-1, ])) 
    # + rnorm(4, 0, sd) - excluding stochasticity for these versions
    timeseries$total[ts] <- sum(timeseries[ts, 1:4])
  }
  # Adding year and quarter columns
  timeseries$year <- rep(1:Year, each = 4)
  timeseries$quarter <- rep(paste0("Q", 1:4), Year)
  
  # Saving timeseries to list
  all_timeseries[[ssd]] <- timeseries
  # Saving vital rate dataframes
  all_vital_rates[[ssd]] <- vital_rates
  
}

names(all_timeseries) <- SSD_list_names
names(all_vital_rates) <- SSD_list_names

#### Creating and saving plots of growth rate, reproductive rate and populations through time

for(ssd in SSD_list_names) {
  # Restructuring growth rate results for plotting
  all_vital_rates[[ssd]]$growth$year <- rep(1:Year, each = 4)
  growth_rates_plots <- melt(all_vital_rates[[ssd]]$growth[-(Year*4),], id.vars = c("year"))
  colnames(growth_rates_plots)[2:3] <- c("size", "growth_rate")
  # Creating growth rate plot
  growth_plot <- ggplot(growth_rates_plots, aes(x = year, y = growth_rate)) +
    geom_line()+
    facet_wrap(vars(size)) +
    theme_bw() +
    labs(title = ssd)
  #### Saving plots to file
  ggsave(filename = paste0(ssd, "_growth.png"), plot = growth_plot,
         path = here("Results/sensitivity_analyses/no_stochasticity/size_dist/Plots/Growth"))
  
  # Restructuring reproductive rate results for plotting
  all_vital_rates[[ssd]]$repro$year <- rep(1:Year, each = 4)
  repro_rates_plots <- melt(all_vital_rates[[ssd]]$repro[-(Year*4),], id.vars = c("year"))
  colnames(repro_rates_plots)[2:3] <- c("size", "repro_rate")
  
  reproductive_plot <- ggplot(repro_rates_plots, aes(x = year, y = repro_rate)) +
    geom_line()+
    facet_wrap(vars(size)) +
    theme_bw() +
    labs(title = ssd)
  #### Saving plots to file
  ggsave(filename = paste0(ssd, "_repro.png"), plot = reproductive_plot,
         path = here("Results/sensitivity_analyses/no_stochasticity/size_dist/Plots/Reproduction"))
  
  ## Plotting proportion of size classes in population through timeseries
  # restructuring timeseries data for plotting, and adding timestep
  tm_plot <- melt(all_timeseries[[ssd]][, -5], id.vars = c("year", "quarter"))
  colnames(tm_plot)[3:4] <- c("size", "pop")
  tm_plot$time_step <- rep(1:(Year*4), 4)
  # Creating population timeseries plot
  population_plot <- ggplot(tm_plot, aes(fill = size, y = pop, x = time_step)) +
    geom_bar(position="stack", stat="identity") +
    theme_bw() +
    labs(title = ssd)
  #### Saving plots to file
  ggsave(filename = paste0(ssd, "_pop.png"), plot = population_plot,
         path = here("Results/sensitivity_analyses/no_stochasticity/size_dist/Plots/Population"))
}

### Plots comparing the final population and vital rate values from each starting population scenario
# Pulling out final population numbers
final_pop <- as.data.frame(matrix(NA, nrow = length(init_size_dist_options), ncol = ncol(all_timeseries[[1]])))
colnames(final_pop) <- colnames(all_timeseries[[1]])

for(ssd in 1:length(SSD_list_names)) {
  final_pop[ssd, ] <- all_timeseries[[ssd]][Year*4, ]
}
final_pop$scenario <- SSD_list_names
# Melting data for plotting
final_pop_melt <- melt(final_pop[, -c(5:7)], id.vars = "scenario")
colnames(final_pop_melt)[2:3] <- c("size", "pop")
# Making scenarios into ordered factors for clarity
final_pop_melt$scenario <- factor(final_pop_melt$scenario, levels = SSD_list_names)

final_ts_plot <- ggplot(final_pop_melt, aes(y = pop, x = scenario, fill = size)) +
  geom_bar(position="stack", stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))

# Saving plot
ggsave(filename = "final_ts_pop_comparison.png", plot = final_ts_plot, 
       path = here("Results/sensitivity_analyses/no_stochasticity/size_dist/Plots"))

# Pulling out final growth rate numbers
final_growth <- as.data.frame(matrix(NA, nrow = length(init_size_dist_options), ncol = ncol(all_vital_rates[[1]]$growth)))
colnames(final_growth) <- colnames(all_vital_rates[[1]]$growth)

for(ssd in 1:length(SSD_list_names)) {
  final_growth[ssd, ] <- all_vital_rates[[ssd]]$growth[(Year*4-1), ]
}
final_growth$scenario <- SSD_list_names
# Melting data for plotting
final_growth_melt <- melt(final_growth[, -4], id.vars = c("scenario"))
colnames(final_growth_melt)[2:3] <- c("size", "growth_rate")
# Making scenarios into ordered factors for clarity
final_growth_melt$scenario <- factor(final_growth_melt$scenario, levels = SSD_list_names)
# Creating growth rate plot
growth_plot <- ggplot(final_growth_melt, aes(x = size, y = growth_rate, fill = scenario)) +
  geom_bar(position="dodge", stat="identity") +
  theme_bw()

# Saving plot
ggsave(filename = "final_ts_growth_comparison.png", plot = growth_plot, 
       path = here("Results/sensitivity_analyses/no_stochasticity/size_dist/Plots"))

# Pulling out final reproductive rate numbers
final_repro <- as.data.frame(matrix(NA, nrow = length(init_size_dist_options), ncol = ncol(all_vital_rates[[1]]$repro)))
colnames(final_repro) <- colnames(all_vital_rates[[1]]$repro)

for(ssd in 1:length(SSD_list_names)) {
  final_repro[ssd, ] <- all_vital_rates[[ssd]]$repro[(Year*4-1), ]
}
final_repro$scenario <- SSD_list_names

# Melting data for plotting
final_repro_melt <- melt(final_repro[, -4], id.vars = c("scenario"))
colnames(final_repro_melt)[2:3] <- c("size", "growth_rate")
# Making scenarios into ordered factors for clarity
final_repro_melt$scenario <- factor(final_repro_melt$scenario, levels = SSD_list_names)
# Creating growth rate plot
repro_plot <- ggplot(final_repro_melt, aes(x = size, y = growth_rate, fill = scenario)) +
  geom_bar(position="dodge", stat="identity") +
  theme_bw()

# Saving plot
ggsave(filename = "final_ts_repro_comparison.png", plot = repro_plot, 
       path = here("Results/sensitivity_analyses/no_stochasticity/size_dist/Plots"))

