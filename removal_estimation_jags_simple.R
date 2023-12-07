library(jagsUI)
library(here)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tictoc)
library(fdrtool)
library(vctrs)
library(scales)

## Loading objects and functions
source(here("Scripts/00_user_inputs.R"))
source(here("Scripts/01_model_functions.R"))
source(here("Scripts/02_results_functions.R"))

############################ Estimation model version with merged I and J and with survival & size transition as logit linear equations ##########################
# p is size-dependent (not method, yet - do that after this part is working)
# I is quarter-dependent
# saveRDS(erad_quarter_results, file = "consist_test_erad_quarter_results.rds")
# consist_erad_quarter_results <- readRDS("consist_test_erad_quarter_results.rds")

# Isolating methods used in this version & the quarters:
obs_methods <- erad_methods[c(2:3)]
obs_quarters <- unique(unlist(erad_quarters[obs_methods]))

# Using only visual data
erad_reformatted_v2 <- all_observations_fun(erad_results_ts = erad_quarter_results,
                                            methods = obs_methods)

# removals_array_v2 <- erad_reformatted_v2$observation[1,,,] # array(dim = c(4,14,4)), for visual only
# effort_array_v2 <- erad_reformatted_v2$effort[1,,] # array(dim = c(14,4)), for visual only

##### Reformatting to merge I and J into one dimension, with odds being visual and evens being trap

# Values needed for array dimensions & loops
S <- 4 # number of size classes 
Q <- length(unique(unlist(erad_quarters[obs_methods]))) # Primary sampling periods (quarters)
I <- rep(max(length(erad_days[[obs_methods[1]]]), length(erad_days[[obs_methods[2]]]))*2, Q) # Secondary sampling periods (days within each quarter) - redo this later, once I re-do how erad_days is formatted to allow it to vary between quarters


removals_array_v2 <- array(dim = c(S, I[1], Q))
effort_array_v2 <- array(dim = c(I[1], Q))
for(q in 1:Q) {
  for(i in 1:(I[q]/2)) {
    removals_array_v2[, (i*2 - 1):(i*2), q] <- cbind(erad_reformatted_v2$observation[1,,i,q], erad_reformatted_v2$observation[2,,i,q])
    effort_array_v2[(i*2 - 1):(i*2), q] <- cbind(erad_reformatted_v2$effort[1,i,q], erad_reformatted_v2$effort[2, i, q])
  }
}

# Days between primary sampling periods when method is used
# For visual data
date_diff <- vector()
for(t in 1:(Q-1)) {
  date_diff[t] <- ((erad_quarters[[2]][t+1]-1)*91 + erad_days[[2]][1]) - 
    ((erad_quarters[[2]][t]-1)*91 + max(unlist(erad_days[2])))
}


# JAGS Removal Estimation Model - simple growth version 
sink("removal_model_simple.jags")
cat("
model {

# Set up first row of N
for(k in 1:S) {
  p.miss[1,k,1:1000] <- rep(1/1000,1000)
  miss[1,k] ~ dcat(p.miss[1,k,1:1000])
  N[k,1,1] <-  N.base[k,1,1] + miss[1,k]
}

# Parameter priors
# beta.p[1] ~ dunif(0, 10) # encounter slope small
# beta.p[2] ~ dunif(0, 10) # encounter slope medium
# beta.p[3] ~ dunif(0, 10) # encounter slope for large
# beta.p[4] ~ dunif(0, 10) # encounter slope for x-large
# alpha.p[1] ~ dunif(-10,10) # encounter intercept small
# alpha.p[2] ~ dunif(-10,10) # encounter intercept medium
# alpha.p[3] ~ dunif(-10,10) # encounter intercept large
# alpha.p[4] ~ dunif(-10, 10) # encounter slope x-large

for(k in 1:S) {
  beta.p[k, 1] ~ dunif(0, 10)
  beta.p[k, 2] ~ dunif(0, 10)
  alpha.p[k, 1] ~ dunif(-10, 10)
  alpha.p[k, 2] ~ dunif(-10,10)
  for(v in vis_days) {
    beta.p[k, v] <- beta.p[k,1] 
    alpha.p[k, v] <- alpha.p[k,1]
  }
  for(r in trap_days) {
    beta.p[k, r] <- beta.p[k,2]
    alpha.p[k, r] <- alpha.p[k,2]
  }
}

# Growth per size class priors
r1 ~ dgamma(1,0.3) # growth for small size class
r2 ~ dgamma(1,0.3) # growth for medium size class
r3 ~ dgamma(1,0.3) # growth for large size class
r4 ~ dgamma(1,0.3) # growth for xlarge size class



# Transition matrix (used for population size class growth & reproduction between primary sampling periods)
for(t in 1:(Q-1)){
  P[1,1,t] <- r1^days_btwn[t]
  P[1,2,t] <- 0
  P[1,3,t] <- 0
  P[1,4,t] <- 0
  P[2,1,t] <- 0
  P[2,2,t] <- r2^days_btwn[t]
  P[2,3,t] <- 0
  P[2,4,t] <- 0
  P[3,1,t] <- 0
  P[3,2,t] <- 0
  P[3,3,t] <- r3^days_btwn[t]
  P[3,4,t] <- 0
  P[4,1,t] <- 0
  P[4,2,t] <- 0
  P[4,3,t] <- 0
  P[4,4,t] <- r4^days_btwn[t]
}


for(t in 1:Q) { # start primary sampling period loop  
  for(k in 1:S) { # start size class loop
      for(i in 1:I) { # start secondary sampling instances loop
      # Calculate encounter probability for each method, secondary sampling instance and primary sampling period
      # Odd columns are visual, even columns are trap
        logit(p[k,i,t]) <- alpha.p[k,i] + beta.p[k,i] * log(xi[i,t])
      #***** version below runs, but is not method specific *****
      #  logit(p[k,i,t]) <- alpha.p[k] + beta.p[k] * log(xi[i,t]) # effort is not size-dependent
      # Calculate removals based on encounter probability and M (population in within i instance)
        Y[k,i,t] ~ dbin(p[k,i,t],N[k,i,t])
      # Calculate N using last time step N minus summed removals
        N[k,i+1,t] <- N[k,i,t] - Y[k,i,t]
      } # end secondary sampling instances loop
  # Calculate remaining population at the end of the primary sampling period
    R[k,t] <- N[k,I,t] - Y[k,I,t]
  } # end size class loop
} # end primary sampling period loop

for (t in 1:(Q-1)) { # start operations between primary sampling period loop
  # Calculate population at beginning of primary sampling period using remaining population from the end of previous sampling period X transition matrix
    D[1,t] ~ dpois(R[1,t]*P[1,1,t])
    D[2,t] ~ dpois(R[2,t]*P[2,2,t])
    D[3,t] ~ dpois(R[3,t]*P[3,3,t])
    D[4,t] ~ dpois(R[4,t]*P[4,4,t])
    # Set up first sampling instance of next primary sampling period
    for(k in 1:S){ # start size class loop
      N[k,1,t+1] <- D[k,t]
    } # end size class loop
} # end between primary sampling period loop

for(t in 1:Q) { #start primary sampling period loop
  # Summing all size classes into a single N for each primary sampling period
  N.sum[t] <- sum(N[,1,t])
} # end primary sampling period loop

} # end model
", fill= TRUE)
sink()

##################### Creating all inputs for jags model ###########
# Initial values for N (i.e. N[k, 1, 1])
Y <- removals_array_v2
N.base <- array(NA_real_, dim = c(S, I[1], Q))
for(k in 1:S) {
  N.base[k,1,1] <- sum(Y[k,,1:Q])
}

# initialize D to > than the number that will be removed in the following year
Y.remove1 <- vector()
Y.remove2 <- vector()
Y.remove3 <- vector()
Y.remove4 <- vector()
for(t in 1:Q){
  Y.remove1[t] <- sum(Y[1,,t])
  Y.remove2[t] <- sum(Y[2,,t])
  Y.remove3[t] <- sum(Y[3,,t])
  Y.remove4[t] <- sum(Y[4,,t])
}

D.init <- array(NA,dim = c(S,(Q-1)))
for(t in 1:(Q-1)){
  D.init[1,t] <- Y.remove1[t+1] + 1
  D.init[2,t] <- Y.remove2[t+1] + 1
  D.init[3,t] <- Y.remove3[t+1] + 1
  D.init[4,t] <- Y.remove4[t+1] + 1
}

# Initial values for select parameters
inits <- function (){
  list(D = D.init)
  # beta.p = runif(4,0,10),
  # alpha.p = runif(4,-10,10))
  
}

# Creating vectors for trap and vis days
vis_days <- seq(3, I[1], 2)
trap_days <- seq(4, I[1], 2)

# Bundle data together
data <- list(Y = removals_array_v2,
             S = S, I = I[1], Q = Q, 
             xi = effort_array_v2, 
             days_btwn = date_diff, 
             N.base = N.base,
             vis_days = vis_days,
             trap_days = trap_days)

# Parameters monitored
parameters <- c("N", "r1", "r2", "r3", "r4", 
                "p")

# MCMC settings
nt <- 1
nb <- 50000
ni <- 100000 + nb
nc <- 3

# Call JAGS Function
output_jags <- jags(data, 
                    inits, 
                    parameters, 
                    "removal_model_simple.jags", 
                    n.chains = nc, 
                    n.thin = nt, 
                    n.iter = ni,
                    n.burnin = nb)

##### Testing model consistency with the same input results - run 50 times 
consist_test_output_jags <- list()
for(iter in 1:50) {
  output_jags_temp <- jags(data, 
                           inits, 
                           parameters, 
                           "removal_model_simple.jags", 
                           n.chains = nc, 
                           n.thin = nt, 
                           n.iter = ni,
                           n.burnin = nb)
  consist_test_output_jags[[iter]] <- output_jags_temp$mean
  print( paste0("iteration ", iter, " complete"))
}

## Saving and restoring the 50 run consistency check data
saveRDS(consist_test_output_jags, file = "consistency_test_outputs_simple_12.5.23.rds")
# consist_test_output_jags <- readRDS("consistency_test_outputs_8.16.23.rds")

mean_outputs <- list()
mean_outputs$N_sum <- as.data.frame(matrix(NA, nrow = 50, ncol = length(obs_quarters)))
colnames(mean_outputs$N_sum) <- paste0("Quarter_", obs_quarters)
mean_outputs$p_vis <- as.data.frame(matrix(NA, nrow = 50, ncol = 4))
colnames(mean_outputs$p_vis) <- size_class_names
mean_outputs$p_trap <- as.data.frame(matrix(NA, nrow = 50, ncol = 4))
colnames(mean_outputs$p_trap) <- size_class_names
for(iter in 1:50) {
  # Separating out encounter probability for visual survey vs trapping; since effort is the same throughout, the encounter probability for each size class is the same through time
  mean_outputs$p_vis[iter, ] <- consist_test_output_jags[[iter]]$p[,1,1] 
  mean_outputs$p_trap[iter, ] <- consist_test_output_jags[[iter]]$p[,2,2] # depends on trap days
  for(t in 1:Q) {
    mean_outputs$N_sum[iter, t] <- consist_test_output_jags[[iter]]$N.sum[t] # didn't actually sum all 4 years, just the first 3 (fixed for next set of runs)
  }
}

par(mfrow = c(1,3))
hist(mean_outputs$N_sum$Quarter_2)
abline(v = nrow(erad_quarter_results$quarter_timeseries[[3]]), col = "red")
hist(mean_outputs$N_sum$Quarter_3)
abline(v = nrow(erad_quarter_results$quarter_timeseries[[4]]), col = "red")
hist(mean_outputs$N_sum$Quarter_4)
abline(v = nrow(erad_quarter_results$quarter_timeseries[[5]]), col = "red")
hist(mean_outputs$N_sum$Quarter_5)
abline(v = nrow(erad_quarter_results$quarter_timeseries[[6]]), col = "red")
hist(mean_outputs$N_sum$Quarter_6)
abline(v = nrow(erad_quarter_results$quarter_timeseries[[7]]), col = "red")
# Visual survey mortality probability, estimated histogram with real value in red
par(mfrow = c(2, 2)) # limits will likely need to be adjusted to show the real values
hist(mean_outputs$p_vis$small, xlim = c(0.0015, 0.0025)) 
abline(v = mortality_prob_erad_methods$visual[1]*erad_coverage$visual, col = "red")
hist(mean_outputs$p_vis$medium)
abline(v = mortality_prob_erad_methods$visual[2]*erad_coverage$visual, col = "red")
hist(mean_outputs$p_vis$large, xlim = c(0.0008, 0.0018))
abline(v = mortality_prob_erad_methods$visual[3]*erad_coverage$visual, col = "red")
hist(mean_outputs$p_vis$xlarge, xlim = c(0.0014, 0.0028))
abline(v = mortality_prob_erad_methods$visual[4]*erad_coverage$visual, col = "red")
# Trap mortality probability, estimated histogram with real value in red
par(mfrow = c(2, 2)) # limits will likely need to be adjusted to show the real values
hist(mean_outputs$p_trap$small, xlim = c(0.0008, 0.0013)) 
abline(v = mortality_prob_erad_methods$trap[1]*erad_coverage$trap, col = "red")
hist(mean_outputs$p_trap$medium)
abline(v = mortality_prob_erad_methods$trap[2]*erad_coverage$trap, col = "red")
hist(mean_outputs$p_trap$large, xlim = c(0.0022, 0.0065))
abline(v = mortality_prob_erad_methods$trap[3]*erad_coverage$trap, col = "red")
hist(mean_outputs$p_trap$xlarge, xlim = c(0.002, 0.011))
abline(v = mortality_prob_erad_methods$trap[4]*erad_coverage$trap, col = "red")



traceplot(output_jags, parameters = "alpha.p")
traceplot(output_jags, parameters = "beta.p")
traceplot(output_jags, parameters = c("r1", "r2", "r3", "r4"))
traceplot(output_jags, parameters = "p")

# Mean N estimates vs simulated real N
est_v_sim_N_plots <- estimated_N_plots(jags_output = output_jags,
                                       erad_quarter_results = erad_quarter_results,
                                       erad_quarters = erad_quarters)

# Separating out encounter probability by method - odd columns are visual survey, even columns are trap
p_1 <- output_jags$sims.list$p[,,seq(1, 28, 2),]
p_2 <- output_jags$sims.list$p[,,seq(2, 29, 2),] # there are 0s in here from the days when no trapping occurred - I guess I should take out the 0s
p_2_v2 <- p_2[,,8:14,2:3] # trapping only occurred in quarters 2 & 3, in second 2 weeks of quarters


## Model estimated p values:
# visual: small = 0.07383756, medium = 0.08769911, large = 0.12376324, xlarge = 0.14988058
# trap: small = 0.07436109, medium = 0.08983350, large = 0.13077347, xlarge = 0.18043694

## Real p values:
# visual: 0.3032118 0.3409591 0.3534641 0.2989624 
# trap: 0.1630739 0.2573346 0.4453922 0.4016056 


# Calculating true size class growth rates
size_class_growth <- as.data.frame(matrix(NA, nrow = 4, ncol = 4))
colnames(size_class_growth) <- c("size_class", paste0("inter-primary_", c(1:3)))
size_class_growth$size_class <- size_class_names

observed_quarters <- c(2,3,6,7)
for(quarter in 1:3) {
  # Separating out population at the beginning of the inter-primary period (after an eradication primary period)
  after_erad_pop <- erad_quarter_results$all_quarters_before_after_erad$pop_after_erad[observed_quarters[quarter]][[1]]
  # Separating out population at the end of the inter-primary period (before the next eradication primary period)
  before_erad_pop <- erad_quarter_results$all_quarters_before_after_erad$pop_before_erad[observed_quarters[quarter+1]][[1]]
  # Adding size class column to both of the above data frames
  for(snake in 1:nrow(after_erad_pop)) {
    after_erad_pop$size_class[snake] <- size_class_fun(after_erad_pop$SVL[snake])
  }
  for(snake in 1:nrow(before_erad_pop)) {
    before_erad_pop$size_class[snake] <- size_class_fun(before_erad_pop$SVL[snake])
  }
  for(size in 1:length(size_class_names)) {
    size_after_pop <- after_erad_pop[after_erad_pop$size_class == size_class_names[size],]
    size_before_pop <- before_erad_pop[before_erad_pop$size_class == size_class_names[size],]
    size_class_growth[size, quarter+1] <- nrow(size_before_pop)/nrow(size_after_pop)
  }
}


# > size_class_growth
# size_class inter-primary_1 inter-primary_2 inter-primary_3
# 1      small       1.0044188       0.9814815       1.0508876
# 2     medium       0.4240150       0.5806452       1.0708661
# 3      large       0.4386792       0.5328467       1.0204082
# 4     xlarge       0.4341232       0.4797297       0.9807692

# Except, r is constant across time steps (it's exponentiated by # days between eradication periods), so
# this isn't quite right... how do I calculate that? Average it?


