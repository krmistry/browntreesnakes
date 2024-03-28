library(jagsUI)
library(here)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tictoc)
library(fdrtool)
library(scales)
library(vctrs)

## Loading objects and functions
source(here("Scripts/00_user_inputs.R"))
source(here("Scripts/01_model_functions.R"))
source(here("Scripts/02_results_functions.R"))

############################ Estimation model version with merged I and J and with survival & size transition as logit linear equations ##########################
# p is size-dependent (not method, yet - do that after this part is working)
# I is quarter-dependent
# saveRDS(erad_quarter_results, file = "consist_test_erad_quarter_results.rds")
# consist_erad_quarter_results <- readRDS("consist_test_erad_quarter_results.rds")
# Isolating methods used in this version:
obs_methods <- erad_methods[c(2:3)]

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

## Quarters when ADS occurrs 
ADS_quarters <- erad_quarters$ADS
# Sorting out which intermediary periods that ADS occurred (1) and when it didn't (2), so 2 different survival rates can be calculated
obs_quarters <- unique(unlist(erad_quarters[obs_methods]))
ADS_quarter_counter <- c(1,1,2,1,2,1,1,2,1,1,2,1,1,1,2)


# JAGS Removal Estimation Model - version 2
sink("removal_model_v2.jags")
cat("
model {

# Set up first row of N
for(k in 1:S) {
  p.miss[1,k,1:2000] <- rep(1/2000,2000)
  miss[1,k] ~ dcat(p.miss[1,k,1:2000])
  N[k,1,1] <-  N.base[k,1,1] + miss[1,k]
}

# Parameter priors

for(k in 1:S) {
  beta.p[k, 1] ~ dunif(0, 10)
  beta.p[k, 2] ~ dunif(0, 10)
  alpha.p[k, 1] ~ dunif(-20, 10)
  alpha.p[k, 2] ~ dunif(-20,10)
  for(v in vis_days) {
  beta.p[k, v] <- beta.p[k,1] 
  alpha.p[k, v] <- alpha.p[k,1]
  }
  for(r in trap_days) {
  beta.p[k, r] <- beta.p[k,2]
  alpha.p[k, r] <- alpha.p[k,2]
  }
}

## Survival hyperparameter - one set for when ADS has occurred, one for when it hasn't
# 
# for(m in 1:2) {
#   s1.day[m] ~ dbeta(1,1)
#   s2.day[m] ~ dbeta(1,1)
#   s3.day[m] ~ dbeta(1,1)
#   s4.day[m] ~ dbeta(1,1)
# }
s1.day ~ dbeta(1,1)
s2.day ~ dbeta(1,1)
s3.day ~ dbeta(1,1)
s4.day ~ dbeta(1,1)


# # Size transition hyperparameter - currently using a constant for diagnostic purposes
# t1.day <- mean_t_day[1] # small -> medium
# t2.day <- mean_t_day[2] # medium -> large
# t3.day <- mean_t_day[3] # large -> xlarge
t1.day ~ dbeta(1,1)
t2.day ~ dbeta(1,1)
t3.day ~ dbeta(1,1)

# # Fecundity hyperparameter
# f2.day ~ dgamma(1,0.3) # medium
# f3.day ~ dgamma(1,0.3) # for both large and xlarge
# # Intercepts (a) and slopes (b) -  currently using a constant for diagnostic purposes
# a1 <- true_F_intercepts[1]
# a2 <- true_F_intercepts[2]
# b1 <- 1
# b2 <- 1
a1 ~ dunif(0, 5)
a2 ~ dunif(0, 5)
b1 ~ dunif(0, 5)
b2 ~ dunif(0, 5)


# Setting up parameters to be sensitive to days between primary sampling periods
for(t in 1:(Q-1)) {
  # Survival
  # s1[t] <- pow(s1.day[ADS_quarter_counter[t]], days_btwn[t])
  # s2[t] <- pow(s2.day[ADS_quarter_counter[t]], days_btwn[t])
  # s3[t] <- pow(s3.day[ADS_quarter_counter[t]], days_btwn[t])
  # s4[t] <- pow(s4.day[ADS_quarter_counter[t]], days_btwn[t])
  s1[t] <- pow(s1.day, days_btwn[t])
  s2[t] <- pow(s2.day, days_btwn[t])
  s3[t] <- pow(s3.day, days_btwn[t])
  s4[t] <- pow(s4.day, days_btwn[t])
  # Size transition
  t1[t] <- 1 - pow((1 - t1.day), days_btwn[t])
  t2[t] <- 1 - pow((1 - t2.day), days_btwn[t])
  t3[t] <- 1 - pow((1 - t3.day), days_btwn[t])
  # Fecundity
  f2[t] <- a1 + b1*(days_btwn[t]/91)
  f3[t] <- a2 + b2*(days_btwn[t]/91)
  f4[t] <- f3[t]
}

# Transition matrix (used for population size class growth & reproduction between primary sampling periods)
for(t in 1:(Q-1)){
  P[1,1,t] <- s1[t]*(1-t1[t])
  P[1,2,t] <- f2[t]
  P[1,3,t] <- f3[t]
  P[1,4,t] <- f4[t]
  P[2,1,t] <- t1[t]*s1[t]
  P[2,2,t] <- s2[t]*(1-t2[t])
  P[2,3,t] <- 0
  P[2,4,t] <- 0
  P[3,1,t] <- 0
  P[3,2,t] <- t2[t]*s2[t]
  P[3,3,t] <- s3[t]*(1-t3[t])
  P[3,4,t] <- 0
  P[4,1,t] <- 0
  P[4,2,t] <- 0
  P[4,3,t] <- t3[t]*s3[t]
  P[4,4,t] <- s4[t]
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
    D[1,t] ~ dpois(R[1,t]*P[1,1,t] + R[2,t]*P[1,2,t] + R[3,t]*P[1,3,t] + R[4,t]*P[1,4,t])
    D[2,t] ~ dpois(R[1,t]*P[2,1,t] + R[2,t]*P[2,2,t])
    D[3,t] ~ dpois(R[2,t]*P[3,2,t] + R[3,t]*P[3,3,t])
    D[4,t] ~ dpois(R[3,t]*P[4,3,t] + R[4,t]*P[4,4,t])
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
             # mean_t_day = mean_t_day,
             # true_F_intercepts = true_F_intercepts)
             # ADS_quarters = ADS_quarters,
             # ADS_quarter_counter = ADS_quarter_counter)

# Parameters monitored
parameters <- c("N", "N.sum", "alpha.p", "beta.p", "s1", "s2", "s3", "s4", 
                "f2", "f3", "f4", "t1", "t2", "t3", "p", "s1.day", "s2.day", "s3.day", "s4.day", "t1.day", "t2.day", "t3.day",
                "a1", "a2", "b1", "b2")

# MCMC settings
nb <- 50000
nt <- 10
ni <- 100000 + nb
nc <- 3

# Call JAGS Function
output_jags <- jags(data, 
                    inits, 
                    parameters, 
                    "removal_model_v2.jags", 
                    n.chains = nc, 
                    n.thin = nt, 
                    n.iter = ni,
                    n.burnin = nb)

# Traceplots of key parameters
traceplot(output_jags, parameters = c("s1","s2","s3","s4"))
traceplot(output_jags, parameters = c("s1.day", "s2.day", "s3.day", "s4.day"))
traceplot(output_jags, parameters = c("t1","t2","t3"))
traceplot(output_jags, parameters = c("t1.day", "t2.day", "t3.day"))
traceplot(output_jags, parameters = c("f2","f3","f4"))
#traceplot(output_jags, parameters = c("f2.day", "f3.day"))
traceplot(output_jags, parameters = c("a1", "a2"))
traceplot(output_jags, parameters = c("b1", "b2"))
traceplot(output_jags, parameters = c("p"))
traceplot(output_jags, parameters = c("alpha.p"))
traceplot(output_jags, parameters = c("beta.p"))


# Mean N estimates vs simulated real N
est_v_sim_N_plots <- estimated_N_plots(jags_output = output_jags,
                       erad_quarter_results = erad_quarter_results,
                       erad_quarters = erad_quarters)



##### Testing model consistency with the same input results - run 50 times
consist_test_output_jags <- list()
for(iter in 1:50) {
  output_jags_temp <- jags(data, 
                      inits, 
                      parameters, 
                      "removal_model_v2.jags", 
                      n.chains = nc, 
                      n.thin = nt, 
                      n.iter = ni,
                      n.burnin = nb)
  consist_test_output_jags[[iter]] <- output_jags_temp$mean
  print( paste0("iteration ", iter, " complete"))
}

## Saving and restoring the 50 run consistency check data
# saveRDS(consist_test_output_jags, file = "consistency_test_outputs_9.27.23.rds")
# consist_test_output_jags <- readRDS("consistency_test_outputs_9.27.23.rds")

mean_outputs <- list()
mean_outputs$N_sum <- as.data.frame(matrix(NA, nrow = 50, ncol = 3))
colnames(mean_outputs$N_sum) <- paste0("Quarter_", c(1:3))
mean_outputs$p_vis <- as.data.frame(matrix(NA, nrow = 50, ncol = 4))
colnames(mean_outputs$p_vis) <- size_class_names
mean_outputs$p_trap <- as.data.frame(matrix(NA, nrow = 50, ncol = 4))
colnames(mean_outputs$p_trap) <- size_class_names
for(iter in 1:50) {
  # Separating out encounter probability for visual survey vs trapping; since effort is the same throughout, the encounter probability for each size class is the same through time
  mean_outputs$p_vis[iter, ] <- consist_test_output_jags[[iter]]$p[,1,1]
  mean_outputs$p_trap[iter, ] <- consist_test_output_jags[[iter]]$p[,16,2]
  for(t in 1:(Q-1)) {
    mean_outputs$N_sum[iter, t] <- consist_test_output_jags[[iter]]$N.sum[t] # didn't actually sum all 4 years, just the first 3 (fixed for next set of runs)
  }
}

par(mfrow = c(1,3))
hist(mean_outputs$N_sum$Quarter_1)
abline(v = nrow(erad_quarter_results$quarter_timeseries[[2]]), col = "red")
hist(mean_outputs$N_sum$Quarter_2)
abline(v = nrow(erad_quarter_results$quarter_timeseries[[3]]), col = "red")
hist(mean_outputs$N_sum$Quarter_3)
abline(v = nrow(erad_quarter_results$quarter_timeseries[[6]]), col = "red")

par(mfrow = c(2, 2))
hist(mean_outputs$p_vis$small, xlim = c(0, 0.02))
abline(v = mortality_prob_erad_methods$visual[1]*erad_coverage$visual, col = "red")
hist(mean_outputs$p_vis$medium)
abline(v = mortality_prob_erad_methods$visual[2]*erad_coverage$visual, col = "red")
hist(mean_outputs$p_vis$large, xlim = c(0.0002, 0.002))
abline(v = mortality_prob_erad_methods$visual[3]*erad_coverage$visual, col = "red")
hist(mean_outputs$p_vis$xlarge, xlim = c(0.001, 0.020))
abline(v = mortality_prob_erad_methods$visual[4]*erad_coverage$visual, col = "red")

hist(mean_outputs$p_trap$small)
hist(mean_outputs$p_trap$medium)
hist(mean_outputs$p_trap$large)
hist(mean_outputs$p_trap$xlarge)

traceplot(output_jags, parameters = c("nu",
                                      "epsilon",
                                      "rho",
                                      "gamma"))


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

### Calculate true vital rates from input data
true_vital_rates <- true_vital_rates_v1_fun(all_erad_quarters = erad_quarters[c("visual","trap")],
                          erad_results_df = erad_quarter_results)

## Calculate daily size transition rates & linear equations for fecundity in order to restrict the 
## hyperparameters for those vital rates
library(pracma)
# Size transition, decomposed to a daily rate for each size class
t_list <- list()
for(k in 1:(S-1)) {
t_list[[k]] <- vector()
  for(q in 1:(Q-1)) {
    t_list[[k]][q] <- 1 - nthroot((1-true_vital_rates$size_transition[k,q+1]), date_diff[q])
  }
}
mean_t_day <- vector() 
for(k in 1:(S-1)) {
  mean_t_day[k] <- mean(t_list[[k]])
}

# Fecundity, creating linear models using days between quarters as the covariate
fecundity_lm_df <- as.data.frame(cbind(c(true_vital_rates$fecundity[1,2],
                                         true_vital_rates$fecundity[1,3],
                                         true_vital_rates$fecundity[1,4],
                                         true_vital_rates$fecundity[1,5],
                                         true_vital_rates$fecundity[1,6]),
                                         c(true_vital_rates$fecundity[2,2],
                                           true_vital_rates$fecundity[2,3],
                                           true_vital_rates$fecundity[2,4],
                                           true_vital_rates$fecundity[2,5],
                                           true_vital_rates$fecundity[2,6]),
                                         c(true_vital_rates$fecundity[3,2],
                                           true_vital_rates$fecundity[3,3],
                                            true_vital_rates$fecundity[3,4],
                                           true_vital_rates$fecundity[3,5],
                                           true_vital_rates$fecundity[3,6]),
                        date_diff/91))
colnames(fecundity_lm_df) <- c(size_class_names[-1], "date_diff")
true_F <- list()
true_F$medium <- lm(medium ~ date_diff, fecundity_lm_df)
true_F$large <- lm(large ~ date_diff, fecundity_lm_df)
true_F$xlarge <- lm(xlarge ~ date_diff, fecundity_lm_df)

true_F_intercepts <- c(true_F$medium$coefficients[1],
                       true_F$large$coefficients[1],
                       true_F$xlarge$coefficients[1])



# ### For hyperparameter priors for survival and size transition, calculating linear models to get 
# ### slopes and intercepts (days between eradication periods is the only covariate)
# 
# # Survival
# survival_lm_df <- as.data.frame(cbind(c(true_vital_rates$survival[1, 2],
#                                         true_vital_rates$survival[1, 3],
#                                         true_vital_rates$survival[1, 4]),
#                                       c(true_vital_rates$survival[2, 2],
#                                           true_vital_rates$survival[2, 3],
#                                           true_vital_rates$survival[3, 4]),
#                                       c(true_vital_rates$survival[3, 2],
#                                         true_vital_rates$survival[3, 3],
#                                         true_vital_rates$survival[3, 4]),
#                                       c(true_vital_rates$survival[4, 2],
#                                         true_vital_rates$survival[4, 3],
#                                         true_vital_rates$survival[4, 4]),
#                         date_diff))
# colnames(survival_lm_df)[1:4] <- size_class_names
# # Linear models for each size class
# Su <- list()
# Su$small <- lm(small ~ date_diff, survival_lm_df)
# Su$medium <- lm(medium ~ date_diff, survival_lm_df)
# Su$large <- lm(large ~ date_diff, survival_lm_df)
# Su$xlarge <- lm(xlarge ~ date_diff, survival_lm_df)
# 
# # Collecting intercepts and slopes from above sizes to see if they are different enough to warrant different priors
# S_intercepts <- vector()
# S_slopes <- vector()
# for(size in 1:4) {
#   S_intercepts[size] <- Su[[size]]$coefficients[1]
#   S_slopes[size] <- Su[[size]]$coefficients[2]
# }
# ## With test case (55 ha pop)
# # > S_slopes
# # [1] -0.001256043 -0.002589571 -0.001520559 -0.001290379
# # > S_intercepts
# # [1] 0.9777542 0.7986311 0.7575742 0.7723645
# 
# 
# # Size transition
# size_transition_lm_df <- as.data.frame(cbind(c(true_vital_rates$size_transition[1, 2],
#                                         true_vital_rates$size_transition[1, 3], 
#                                         true_vital_rates$size_transition[1, 4]),
#                                       c(true_vital_rates$size_transition[2, 2],
#                                         true_vital_rates$size_transition[2, 3], 
#                                         true_vital_rates$size_transition[3, 4]),
#                                       c(true_vital_rates$size_transition[3, 2],
#                                         true_vital_rates$size_transition[3, 3], 
#                                         true_vital_rates$size_transition[3, 4]),
#                                       date_diff))
# colnames(size_transition_lm_df)[1:3] <- size_class_names[-1]
# # Linear models for each size class
# ST <- list()
# ST$medium <- lm(medium ~ date_diff, size_transition_lm_df)
# ST$large <- lm(large ~ date_diff, size_transition_lm_df)
# ST$xlarge <- lm(xlarge ~ date_diff, size_transition_lm_df)
# 
# # Collecting intercepts and slopes from above sizes to see if they are different enough to warrant different priors
# ST_intercepts <- vector()
# ST_slopes <- vector()
# for(size in 1:3) {
#   ST_intercepts[size] <- ST[[size]]$coefficients[1]
#   ST_slopes[size] <- ST[[size]]$coefficients[2]
# }
# ## With test case (55 ha pop)
# # > ST_intercepts
# # [1]  0.003838999 -0.008183079  0.008495881
# # > ST_slopes
# # [1] 1.161147e-04 7.669384e-04 7.134546e-05



# Fecundity true value (used to evaluate prior choice)
true_vital_rates$fecundity
## Current prior is gamma(alpha = 1, beta = 0.3), and that would encompass these values, 
## so the prior should be fine, I think



###### Reformatted for single method 


# Isolating methods used in this version:
obs_methods <- erad_methods[2]

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

## Quarters when ADS occurrs 
ADS_quarters <- erad_quarters$ADS
# Sorting out which intermediary periods that ADS occurred (1) and when it didn't (2), so 2 different survival rates can be calculated
obs_quarters <- unique(unlist(erad_quarters[obs_methods]))
ADS_quarter_counter <- c(1,1,2,1,2,1,1,2,1,1,2,1,1,1,2)


# JAGS Removal Estimation Model - version 2
sink("removal_model_v2.jags")
cat("
model {

# Set up first row of N
for(k in 1:S) {
  p.miss[1,k,1:2000] <- rep(1/2000,2000)
  miss[1,k] ~ dcat(p.miss[1,k,1:2000])
  N[k,1,1] <-  N.base[k,1,1] + miss[1,k]
}

# Parameter priors

for(k in 1:S) {
  beta.p[k, 1] ~ dunif(0, 10)
  beta.p[k, 2] ~ dunif(0, 10)
  alpha.p[k, 1] ~ dunif(-20, 10)
  alpha.p[k, 2] ~ dunif(-20,10)
  for(v in vis_days) {
  beta.p[k, v] <- beta.p[k,1] 
  alpha.p[k, v] <- alpha.p[k,1]
  }
  for(r in trap_days) {
  beta.p[k, r] <- beta.p[k,2]
  alpha.p[k, r] <- alpha.p[k,2]
  }
}

## Survival hyperparameter - one set for when ADS has occurred, one for when it hasn't
# 
# for(m in 1:2) {
#   s1.day[m] ~ dbeta(1,1)
#   s2.day[m] ~ dbeta(1,1)
#   s3.day[m] ~ dbeta(1,1)
#   s4.day[m] ~ dbeta(1,1)
# }
s1.day ~ dbeta(1,1)
s2.day ~ dbeta(1,1)
s3.day ~ dbeta(1,1)
s4.day ~ dbeta(1,1)


# Size transition hyperparameter
t1.day ~ dbeta(1,3) # small -> medium
t2.day ~ dbeta(1,3) # medium -> large
t3.day ~ dbeta(1,3) # large -> xlarge

# # Fecundity hyperparameter
# f2.day ~ dgamma(1,0.3) # medium
# f3.day ~ dgamma(1,0.3) # for both large and xlarge
# Intercepts (a) and slopes (b)
a1 ~ dunif(0, 5)
a2 ~ dunif(0, 5)
b1 ~ dunif(0, 2)
b2 ~ dunif(0, 2)


# Setting up parameters to be sensitive to days between primary sampling periods
for(t in 1:(Q-1)) {
  # Survival
  # s1[t] <- pow(s1.day[ADS_quarter_counter[t]], days_btwn[t])
  # s2[t] <- pow(s2.day[ADS_quarter_counter[t]], days_btwn[t])
  # s3[t] <- pow(s3.day[ADS_quarter_counter[t]], days_btwn[t])
  # s4[t] <- pow(s4.day[ADS_quarter_counter[t]], days_btwn[t])
  s1[t] <- pow(s1.day, days_btwn[t])
  s2[t] <- pow(s2.day, days_btwn[t])
  s3[t] <- pow(s3.day, days_btwn[t])
  s4[t] <- pow(s4.day, days_btwn[t])
  # Size transition
  t1[t] <- 1 - pow((1 - t1.day), days_btwn[t])
  t2[t] <- 1 - pow((1 - t2.day), days_btwn[t])
  t3[t] <- 1 - pow((1 - t3.day), days_btwn[t])
  # Fecundity
  f2[t] <- a1 + b1*(days_btwn[t]/91)
  f3[t] <- a2 + b2*(days_btwn[t]/91)
  f4[t] <- f3[t]
}

# Transition matrix (used for population size class growth & reproduction between primary sampling periods)
for(t in 1:(Q-1)){
  P[1,1,t] <- s1[t]*(1-t1[t])
  P[1,2,t] <- f2[t]
  P[1,3,t] <- f3[t]
  P[1,4,t] <- f4[t]
  P[2,1,t] <- t1[t]*s1[t]
  P[2,2,t] <- s2[t]*(1-t2[t])
  P[2,3,t] <- 0
  P[2,4,t] <- 0
  P[3,1,t] <- 0
  P[3,2,t] <- t2[t]*s2[t]
  P[3,3,t] <- s3[t]*(1-t3[t])
  P[3,4,t] <- 0
  P[4,1,t] <- 0
  P[4,2,t] <- 0
  P[4,3,t] <- t3[t]*s3[t]
  P[4,4,t] <- s4[t]
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
    D[1,t] ~ dpois(R[1,t]*P[1,1,t] + R[2,t]*P[1,2,t] + R[3,t]*P[1,3,t] + R[4,t]*P[1,4,t])
    D[2,t] ~ dpois(R[1,t]*P[2,1,t] + R[2,t]*P[2,2,t])
    D[3,t] ~ dpois(R[2,t]*P[3,2,t] + R[3,t]*P[3,3,t])
    D[4,t] ~ dpois(R[3,t]*P[4,3,t] + R[4,t]*P[4,4,t])
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
# ADS_quarters = ADS_quarters,
# ADS_quarter_counter = ADS_quarter_counter)

# Parameters monitored
parameters <- c("N", "N.sum", "alpha.p", "beta.p", "s1", "s2", "s3", "s4", 
                "f2", "f3", "f4", "t1", "t2", "t3", "p",  "R", "D", "s1.day", "s2.day", "s3.day", "s4.day", "t1.day", "t2.day", "t3.day",
                "a1", "a2", "b1", "b2")

# MCMC settings
ni <- 20000
nt <- 1
nb <- 10000
nc <- 3

# Call JAGS Function
output_jags <- jags(data, 
                    inits, 
                    parameters, 
                    "removal_model_v2.jags", 
                    n.chains = nc, 
                    n.thin = nt, 
                    n.iter = ni,
                    n.burnin = nb)

# Traceplots of key parameters
traceplot(output_jags, parameters = c("s1","s2","s3","s4"))
traceplot(output_jags, parameters = c("s1.day", "s2.day", "s3.day", "s4.day"))
traceplot(output_jags, parameters = c("t1","t2","t3"))
traceplot(output_jags, parameters = c("t1.day", "t2.day", "t3.day"))
traceplot(output_jags, parameters = c("f2","f3","f4"))
#traceplot(output_jags, parameters = c("f2.day", "f3.day"))
traceplot(output_jags, parameters = c("a1", "a2"))
traceplot(output_jags, parameters = c("b1", "b2"))
traceplot(output_jags, parameters = c("p"))
traceplot(output_jags, parameters = c("alpha.p"))
traceplot(output_jags, parameters = c("beta.p"))


# Mean N estimates vs simulated real N
est_v_sim_N_plots <- estimated_N_plots(jags_output = output_jags,
                                       erad_quarter_results = erad_quarter_results,
                                       erad_quarters = erad_quarters)
est_v_sim_N_plots$est_vs_sim_plot_1




