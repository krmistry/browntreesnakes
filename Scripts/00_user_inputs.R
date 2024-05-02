############### User & Modeler Inputs to BTS Simulation Model #################

library(jagsUI)

####################### IBM objects #######################################
size_class_names <- c("small", "medium", "large", "xlarge")

# Size limits of the 4 size classes
size_class_limits <- as.data.frame(matrix(NA, nrow = 4, ncol = 2))
colnames(size_class_limits) <- c("lower", "upper")
rownames(size_class_limits) <- size_class_names
size_class_limits$lower <- c(250, 850, 950, 1150)
size_class_limits$upper <- c(850, 950, 1150, 1500) # upper limit of largest size only used 
# for initializing sizes in first quarter population, there isn't an upper limit on the 
# largest size class otherwise

## Gompertz growth coefficients
gompertz_coefficients <- list()
tau <- paste0("tau_", c(0.25, 0.5, 0.75, 0.9, 0.95))
for(i in 1:length(tau)) {
  gompertz_coefficients[[i]] <- as.data.frame(matrix(NA, nrow = 4, ncol = 2))
  colnames(gompertz_coefficients[[i]]) <- c("F", "M")
  rownames(gompertz_coefficients[[i]]) <- c("intercept", "slope", "intercept_SE", "slope_SE"  )
}
names(gompertz_coefficients) <- tau
# Manually populating each matrix using Cade's coefficient results (intercepts are averages)
# tau = 0.25
gompertz_coefficients[[1]]$F <- c(0.20213130, # intercept
                                  -0.52878937, # slope (same for both sexes)
                                  0.01089347, # intercept SE
                                  0.01922998) # slope SE (same for both sexes)
gompertz_coefficients[[1]]$M <- c(0.20213130 + 0.05100731, 
                                  -0.52878937, 
                                  0.01345507, 
                                  0.01922998)
# tau = 0.5
gompertz_coefficients[[2]]$F <- c(0.28056394, 
                                  -0.55530229, 
                                  0.01284185, 
                                  0.01876979)
gompertz_coefficients[[2]]$M <- c(0.28056394 + 0.07215225, 
                                  -0.55530229, 
                                  0.01593324, 
                                  0.01876979)
# tau = 0.75
gompertz_coefficients[[3]]$F <- c(0.3817238, 
                                  -0.6409609, 
                                  0.01292154, 
                                  0.02458514)
gompertz_coefficients[[3]]$M <- c(0.3817238 + 0.05100731, 
                                  -0.6409609, 
                                  0.01703584, 
                                  0.02458514)
# tau = 0.9
gompertz_coefficients[[4]]$F <- c(0.5498842, 
                                  -0.7398035, 
                                  0.02166973, 
                                  0.03856088)
gompertz_coefficients[[4]]$M <- c(0.5498842+0.05100731, 
                                  -0.7398035, 
                                  0.02581472, 
                                  0.03856088)
# tau = 0.95
gompertz_coefficients[[5]]$F <- c(0.6224156, 
                                  -0.8061105, 
                                  0.04408919, 
                                  0.05613037)
gompertz_coefficients[[5]]$M <- c(0.6224156 + 0.1634791, 
                                  -0.8061105, 
                                  0.05512590, 
                                  0.05613037)

# Turning the above annual values into daily values
for(t in tau) {
  gompertz_coefficients[[t]] <- gompertz_coefficients[[t]]/365
}

# Calculate the asymptotic values for each of these growth curves 
# (in order to make sure that an inappropriate quantile isn't assigned if a snake has a 
# large SVL in the initial population)
gompertz_asymptotes <- list()
for(i in 1:length(tau)) {
  gompertz_asymptotes[[i]] <- vector()
  gompertz_asymptotes[[i]][1] <- 700*exp(-gompertz_coefficients[[i]]$F[1]/gompertz_coefficients[[i]]$F[2])
  gompertz_asymptotes[[i]][2] <- 700*exp(-gompertz_coefficients[[i]]$M[1]/gompertz_coefficients[[i]]$M[2])
  names(gompertz_asymptotes[[i]]) <- c("F", "M")
}
names(gompertz_asymptotes) <- tau

# Possible growth quantiles that can be assigned to each snake
growth_quantiles <- paste0("tau_", c(0.25, 0.5, 0.75, 0.9, 0.95))


# Parameters for reproduction maturity probability distribution
maturity_mean <- 910 + (1025-910)/2
maturity_sd <- (1025-910)/3.29

# Parameter creating offspring - value based on Nafus et al. 2022 data (genetically tracing offspring to parents)
lambda <- 2.2

# Mortality rates for each size class, based on data from Nafus et al. 2022
# (not including life span cutoff for now, we'll see how it goes and add it if snakes are 
# living too long)
N_mortality <- c("small" = 0.001,
                 "medium" = 0.0004,
                 "large" = 0.0004,
                 "xlarge" = 0.0004)


####################### Eradication objects #######################################

# Trap & visual survey encounter/mortality rates from Staci's results files:
load(here("NWFNTRAPALL_SCRpstar.Rdata"))
# Point estimates (mean) from posterior distributions for encounter rates for traps 
# - rates are by size class
trap_rates_all <- out$mean$p0 
mean_trap_rate <- mean(trap_rates_all)
load(here("NWFNVISALL_SCRpstar.Rdata"))
# Point estimates (mean) from posterior distributions for encounter rates for visual survey
# - rates are by size class
visual_rates_all <- out$mean$p0 
mean_vis_rate <- mean(visual_rates_all)

# Names of the eradication methods possible:
erad_methods <- c("ADS", 
                  "visual",
                  "trap",
                  "bait_tube")

# Size limits for snakes affected for each eradication method
erad_methods_size_affected <- list()
erad_methods_size_affected$ADS <- c(850, 1100) # Clark et al 2018 for lower limit
erad_methods_size_affected$visual <- c(0, Inf) # Clark et al 2018, Rodda et al 2007
erad_methods_size_affected$trap <- c(800, 1100) # Clark et al 2018; some captures of 700-900 mm SVL snakes, but only partially effective
erad_methods_size_affected$bait_tube <- c(850, 1100) # Clark et al 2018 for lower limit

# Encounter/mortality rates for each eradication method - for now, encounter = mortality, could change this in the future
mortality_prob_erad_methods <- list()
mortality_prob_erad_methods$ADS <- rep(0.43, 4) # Nafus 2022 - need to think about this further
mortality_prob_erad_methods$visual <- visual_rates_all # overall average from Staci's model results
mortality_prob_erad_methods$trap <- trap_rates_all # overall average from Staci's model results, but there's more variation between size classes than for visual (0.0016-0.0044)
mortality_prob_erad_methods$bait_tube <- rep(mean_trap_rate, 4) # placeholder value, same as trap since they are similar methods
for(method in 1:length(mortality_prob_erad_methods)) {
  names(mortality_prob_erad_methods[[method]]) <- size_class_names
}

# Artificially increasing detection rates while I'm working on making the estimation model work
# mortality_prob_erad_methods[[2]] <- mortality_prob_erad_methods[[2]]*100
# mortality_prob_erad_methods[[3]] <- mortality_prob_erad_methods[[3]]*100
# mortality_prob_erad_methods[[4]] <- mortality_prob_erad_methods[[4]]*10

# Standard amount of effort per day for each eradication method, in hours (placeholders, figure out what the actual number of hours is later)
effort_erad_methods <- list()
effort_erad_methods$ADS <- 4
effort_erad_methods$visual <- 4
effort_erad_methods$trap <- 12 # assumes that the trap is only working at night - check to see what the convention is for BTS
effort_erad_methods$bait_tube <- 12 # Same as trap



