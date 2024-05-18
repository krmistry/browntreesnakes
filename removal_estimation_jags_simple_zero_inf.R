library(jagsUI)
library(here)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tictoc)
library(fdrtool)
library(vctrs)
library(scales)
library(vctrs)

## Loading objects and functions
source(here("Scripts/00_user_inputs.R"))
source(here("Scripts/01_model_functions.R"))
source(here("Scripts/02_results_functions.R"))

############################ Estimation model version with merged I and J and with survival & size transition as logit linear equations ##########################
# p is size and method-dependent 
# I is quarter-dependent (but not really utilizing this yet)

## Using one of the runs from ADS first scenario (because low effort, lots of 0s) to test
IBM_alt_2_outputs_list <- readRDS(here("Results", "DoD_final_report", "IBM_outputs_alt_2_50.RDS"))
#IBM_alt_1_outputs_list <- readRDS(here("Results", "DoD_final_report", "IBM_outputs_alt_1_50.RDS"))
erad_quarter_results <- IBM_alt_2_outputs_list[[1]]
#erad_quarter_results <- IBM_alt_1_outputs_list[[1]]

# Isolating methods used in this version & the quarters:
all_effort <- bind_rows(erad_quarter_results$all_effort)
observed_effort <- all_effort[all_effort$method %in% erad_methods[c(2:3)],]
obs_methods <- unique(observed_effort$method) 
obs_quarters <- unique(observed_effort$quarter) # double check that this works when there are empty dataframes from quarters that didn't have observations
obs_days <- unique(observed_effort$day)

# Values needed for array dimensions & loops
S <- length(size_class_names) # number of size classes 
Q <- length(obs_quarters) # Primary sampling periods (quarters)
J <- length(obs_methods)
I <- rep(length(unique(unlist(erad_days[obs_methods])))*length(obs_methods), Q) # Secondary sampling periods (days within each quarter) - redo this later, once I re-do how erad_days is formatted to allow it to vary between quarters

# Set up array of method days (if more than 1 method occurs, then each method gets its own day in the model, even if they occur on the same day in real life)
method_days <- list()
for(j in 1:J) {
  method_days[[j]] <- seq(j, I[1], J)
}


# Reformatting observations into effort and removals arrays for model
erad_reformatted_v2 <- all_observations_fun(erad_results_ts = erad_quarter_results,
                                            methods = obs_methods,
                                            num_methods = J)

removals_array_v2 <- array(dim = c(S, I[1], Q))
effort_array_v2 <- array(dim = c(I[1], Q))
for(j in 1:J) { # check that the code below works for 2 methods (only used it with one so far)
  for(day in 1:length(method_days[[j]])) {
    removals_array_v2[,method_days[[j]][day],] <- erad_reformatted_v2$observation[j,,day,]
    effort_array_v2[method_days[[j]][day],] <- erad_reformatted_v2$effort[j,day,]
  }
}

## Days between primary sampling periods when methods are used 
# With current set up, assumed that the days in each quarter are the same 
date_diff <- vector()
for(t in 1:(Q-1)) {
  date_diff[t] <- ((obs_quarters[t+1]-1)*91 + obs_days[1]) - 
    ((obs_quarters[t]-1)*91 + max(max(obs_days)))
}


# JAGS Removal Estimation Model - simple growth with zero inflation
sink("removal_model_simple_zero-inf.jags")
cat("
model {

# Set up first row of N
for(k in 1:S) {
  p.miss[1,k,1:N_prior] <- rep(1/N_prior,N_prior)
  miss[1,k] ~ dcat(p.miss[1,k,1:N_prior])
  N[k,1,1] <-  N.base[k,1,1] + miss[1,k]
}

# Parameter priors

# Encounter probability intercept and slope for each method
for(k in 1:S) {
  for(j in 1:J) {
      beta.p[k, j] ~ dunif(0, 10)
      alpha.p[k, j] ~ dunif(-10, 10)
    for(day in method_days[j,]) {
      beta.p[k, day] <- beta.p[k, j]
      alpha.p[k, day] <- alpha.p[k, j]
    }
  }
}


# Growth per size class priors
r1 ~ dgamma(1,0.3) # growth for small size class
r2 ~ dgamma(1,0.3) # growth for medium size class
r3 ~ dgamma(1,0.3) # growth for large size class
r4 ~ dgamma(1,0.3) # growth for xlarge size class

# Zero-inflation prior
# with different variables for each size class
# for(k in 1:S) {
#   psi[k] ~ dunif(0, 1)
# }
# with two variables, one for small snakes and the other for the 3 upper size classes
psi[1] ~ dunif(0,1)
psi[2] ~ dunif(0,1)
psi[3] <- psi[2]
psi[4] <- psi[2]

# Transition matrix (used for population size class growth between primary sampling periods)
for(t in 1:(Q-1)){
  P[1,t] <- r1^days_btwn[t]
  P[2,t] <- r2^days_btwn[t]
  P[3,t] <- r3^days_btwn[t]
  P[4,t] <- r4^days_btwn[t]
}


for(t in 1:Q) { # start primary sampling period loop  
  for(k in 1:S) { # start size class loop
      for(i in 1:I) { # start secondary sampling instances loop
      # Calculate encounter probability for each method, secondary sampling instance and primary sampling period
      # Odd columns are visual, even columns are trap
        logit(p[k,i,t]) <- alpha.p[k,i] + beta.p[k,i] * log(xi[i,t])
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
  for(k in 1:S) {
    # Zero-inflation step
    z[k, t] ~ dbern(psi[k])
    # Calculate population at beginning of primary sampling period using remaining population from the end of previous sampling period X transition matrix
    D[k,t] ~ dpois((R[k,t]*P[k,t])*z[k, t] + 0.00001)
  # }
  # # Set up first sampling instance of next primary sampling period
  # for(k in 1:S){ # start size class loop
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
}

# Convert method_days list into an array (take off the first day, which is set differently in the model)
j_days <- array(NA, dim = c(length(method_days), (length(method_days[[1]])-1)))
for(j in 1:J) {
  j_days[j,] <- method_days[[j]][-1]
}

# Value for a somewhat informative prior for N
N_prior <- 30*area_size

# Bundle data together
data <- list(Y = removals_array_v2,
             S = S, I = I[1], Q = Q, J = J,
             xi = effort_array_v2, 
             days_btwn = date_diff, 
             N.base = N.base,
             method_days = j_days,
             N_prior = N_prior)

# Parameters monitored
parameters <- c("N", "r1", "r2", "r3", "r4", "alpha.p", "beta.p",
                "p", "N.sum", "psi")



# MCMC settings
nt <- 1
nb <- 10000
ni <- 30000 + nb
nc <- 3

# Call JAGS Function
output_jags <- jags(data, 
                    inits, 
                    parameters, 
                    "removal_model_simple_zero-inf.jags", 
                    n.chains = nc, 
                    n.thin = nt, 
                    n.iter = ni,
                    n.burnin = nb)


# Parameter trace plots
traceplot(output_jags, parameters = "alpha.p")
traceplot(output_jags, parameters = "beta.p")
traceplot(output_jags, parameters = c("r1", "r2", "r3", "r4"))
traceplot(output_jags, parameters = "p")
traceplot(output_jags, parameters = "psi")

# Mean N estimates vs simulated real N
est_v_sim_N_plots <- estimated_N_plots(jags_output = output_jags,
                                       erad_quarter_results = erad_quarter_results,
                                       erad_quarters = erad_quarters)



#saveRDS(output_jags, file = here("Results", "zero-inf", "All_in_var1.rds"))
#saveRDS(output_jags, file = here("Results", "zero-inf", "ADS_first_var1.rds"))
#saveRDS(output_jags, file = here("Results", "zero-inf", "ADS_first_2_psi_var1.rds"))
#saveRDS(output_jags, file = here("Results", "zero-inf", "All_in_2_psi_var1.rds"))

# Reading in results with and without zero inflation and the input data for one run of each alternative strategy
zero_inf_jags_output_alt_1 <- readRDS(here("Results", "zero-inf", "All_in_var1.rds"))
zero_inf_2_psi_jags_output_alt_1 <- readRDS(here("Results", "zero-inf", "All_in_2_psi_var1.rds"))
jags_output_alt_1 <- alt_model_metrics[[1]][[1]]$summed_results
IBM_output_alt_1 <- alt_model_metrics[[1]][[1]]$summed_data

zero_inf_jags_output_alt_2 <- readRDS(here("Results", "zero-inf", "ADS_first_var1.rds"))
zero_inf_2_psi_jags_output_alt_2 <- readRDS(here("Results", "zero-inf", "ADS_first_2_psi_var1.rds"))
jags_output_alt_2 <- alt_model_metrics[[2]][[1]]$summed_results
IBM_output_alt_2 <- alt_model_metrics[[2]][[1]]$summed_data


########## Alternative 1 - all in scenario

## Formatting data into data frames with the same structures

# Zero inflation with size-class based psi (4, one for each size class)
zero_inf_1 <- estimated_N_plots(jags_output = zero_inf_jags_output_alt_1,
                                       erad_quarter_results = IBM_alt_1_outputs_list[[1]],
                                       erad_quarters = erad_quarters)
zero_inf_alt_1 <- zero_inf_1$mean_N_df
# Zero inflation with 2 psi (small, upper 3)
zero_inf_2_psi_1 <- estimated_N_plots(jags_output = zero_inf_2_psi_jags_output_alt_1 ,
                                    erad_quarter_results = IBM_alt_1_outputs_list[[1]],
                                    erad_quarters = erad_quarters)
zero_inf_2_psi_alt_1 <- zero_inf_2_psi_1$mean_N_df
# From IBM output
IBM_alt_1 <- melt(IBM_output_alt_1[-6], id.vars = "Quarter")
colnames(IBM_alt_1)[2:3] <- c("size_class", "N")
# From already formatted jags output (no zero inflation)
no_zero_inf_alt_1 <- melt(jags_output_alt_1$estimated_N[-6], id.vars = "Quarter")
colnames(no_zero_inf_alt_1)[2:3] <- c("size_class", "N")
CIs_list <- list()
CIs_list$lower <- as.data.frame(matrix(0, nrow = length(size_class_names), ncol = length(obs_quarters)))
colnames(CIs_list$lower) <- size_class_names
CIs_list$upper <- as.data.frame(matrix(0, nrow = length(size_class_names), ncol = length(obs_quarters)))
colnames(CIs_list$upper) <- size_class_names
for(size in 1:length(size_class_names)) {
  for(quarter in 1:length(obs_quarters)) {
    CIs_list$lower[quarter, size] <- jags_output_alt_1$estimated_N_95_CI_lower[size, 1, quarter]
    CIs_list$upper[quarter, size] <- jags_output_alt_1$estimated_N_95_CI_upper[size, 1, quarter]
  }
}
no_zero_inf_alt_1$lower_bound <- melt(CIs_list$lower)$value
no_zero_inf_alt_1$upper_bound <- melt(CIs_list$upper)$value


### Plotting results (with and without zero inflation) vs data for both strategies
ggplot() +
  # Zero-inflated with psi for each size class
  geom_line(data = zero_inf_alt_1, aes(y = N, x = Quarter), size = 2, color = hue_pal()(4)[3]) +
  geom_ribbon(data = zero_inf_alt_1, aes(y = N, x = Quarter, ymin = lower_bound, ymax = upper_bound), alpha = 0.3, linetype = 0, fill = hue_pal()(4)[3]) +
  # Zero-inflated with 2 psi (small, upper 3)
  geom_line(data = zero_inf_2_psi_alt_1, aes(y = N, x = Quarter), size = 2, color = hue_pal()(4)[4]) +
  geom_ribbon(data = zero_inf_2_psi_alt_1, aes(y = N, x = Quarter, ymin = lower_bound, ymax = upper_bound), alpha = 0.3, linetype = 0, fill = hue_pal()(4)[4]) +
  # No zero-inflation
  geom_line(data = no_zero_inf_alt_1, aes(y = N, x = Quarter), size = 2, color = hue_pal()(4)[1]) +
  geom_ribbon(data = no_zero_inf_alt_1, aes(y = N, x = Quarter, ymin = lower_bound, ymax = upper_bound), alpha = 0.3, linetype = 0, fill = hue_pal()(4)[1]) +
  # IBM data
  geom_line(data = IBM_alt_1, aes(y = N, x = Quarter), color = "black", size = 2) +
  facet_wrap("size_class", scales = "free_y") +
  guides(fill = "none") +
  theme_bw()



########## Alternative 2 - ADS first scenario

## Formatting data into data frames with the same structures

# Zero inflation with size-class based psi (4, one for each size class)
zero_inf_2 <- estimated_N_plots(jags_output = zero_inf_jags_output_alt_2,
                              erad_quarter_results = IBM_alt_2_outputs_list[[1]],
                              erad_quarters = erad_quarters)
zero_inf_alt_2 <- zero_inf_2$mean_N_df
# Zero inflation with 2 psi (small, upper 3)
zero_inf_2_psi_2 <- estimated_N_plots(jags_output = zero_inf_2_psi_jags_output_alt_2 ,
                                    erad_quarter_results = IBM_alt_2_outputs_list[[1]],
                                    erad_quarters = erad_quarters)
zero_inf_2_psi_alt_2 <- zero_inf_2_psi_2$mean_N_df
# From IBM output
IBM_alt_2 <- melt(IBM_output_alt_2[-6], id.vars = "Quarter")
colnames(IBM_alt_2)[2:3] <- c("size_class", "N")
# From already formatted jags output (no zero inflation)
no_zero_inf_alt_2 <- melt(jags_output_alt_2$estimated_N[-6], id.vars = "Quarter")
colnames(no_zero_inf_alt_2)[2:3] <- c("size_class", "N")
CIs_list <- list()
CIs_list$lower <- as.data.frame(matrix(0, nrow = length(size_class_names), ncol = length(obs_quarters)))
colnames(CIs_list$lower) <- size_class_names
CIs_list$upper <- as.data.frame(matrix(0, nrow = length(size_class_names), ncol = length(obs_quarters)))
colnames(CIs_list$upper) <- size_class_names
for(size in 1:length(size_class_names)) {
  for(quarter in 1:length(obs_quarters)) {
    CIs_list$lower[quarter, size] <- jags_output_alt_2$estimated_N_95_CI_lower[size, 1, quarter]
    CIs_list$upper[quarter, size] <- jags_output_alt_2$estimated_N_95_CI_upper[size, 1, quarter]
  }
}
no_zero_inf_alt_2$lower_bound <- melt(CIs_list$lower)$value
no_zero_inf_alt_2$upper_bound <- melt(CIs_list$upper)$value


### Plotting results (with and without zero inflation) vs data for both strategies
ggplot() +
  # Zero-inflated with psi for each size class
  geom_line(data = zero_inf_alt_2, aes(y = N, x = Quarter), size = 2, color = hue_pal()(4)[3]) +
  geom_ribbon(data = zero_inf_alt_2, aes(y = N, x = Quarter, ymin = lower_bound, ymax = upper_bound), alpha = 0.3, linetype = 0, fill = hue_pal()(4)[3]) +
  # Zero-inflated with 2 psi (small, upper 3)
  geom_line(data = zero_inf_2_psi_alt_2, aes(y = N, x = Quarter), size = 2, color = hue_pal()(4)[4]) +
  geom_ribbon(data = zero_inf_2_psi_alt_2, aes(y = N, x = Quarter, ymin = lower_bound, ymax = upper_bound), alpha = 0.3, linetype = 0, fill = hue_pal()(4)[4]) +
  # No zero-inflation
  geom_line(data = no_zero_inf_alt_2, aes(y = N, x = Quarter), size = 2, color = hue_pal()(4)[1]) +
  geom_ribbon(data = no_zero_inf_alt_2, aes(y = N, x = Quarter, ymin = lower_bound, ymax = upper_bound), alpha = 0.3, linetype = 0, fill = hue_pal()(4)[1]) +
  # IBM data
  geom_line(data = IBM_alt_2, aes(y = N, x = Quarter), color = "black", size = 2) +
  facet_wrap("size_class", scales = "free_y") +
  guides(fill = "none") +
  theme_bw()


### Running zero-inflation (with 2 parameters) for all 50 variants of alternative 2 (low effort 
### scenario) to see if the improvement in estimation is consistent

# List to capture jags output (not saving parameters because of space)
jags_output_list <- list()

for(variant in 1:1) {
  # Reformat simulation results
  erad_reformatted_v2 <- all_observations_fun(erad_results_ts = erad_results_list[[variant]],
                                              methods = obs_methods)
  # Isolate observations and effort per day
  removals_array_v2 <- array(dim = c(S, I[1], Q))
  effort_array_v2 <- array(dim = c(I[1], Q))
  for(q in 1:Q) {
    for(i in 1:(I[q]/2)) {
      removals_array_v2[, (i*2 - 1):(i*2), q] <- cbind(erad_reformatted_v2$observation[1,,i,q], erad_reformatted_v2$observation[2,,i,q])
      effort_array_v2[(i*2 - 1):(i*2), q] <- cbind(erad_reformatted_v2$effort[1,i,q], erad_reformatted_v2$effort[2, i, q])
    }
  }
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
  # Bundle data together
  data <- list(Y = removals_array_v2,
               S = S, I = I[1], Q = Q, 
               xi = effort_array_v2, 
               days_btwn = date_diff, 
               N.base = N.base,
               # vis_days = vis_days,
               # trap_days = trap_days,
               N_prior = N_prior)
  # ADS_quarters = ADS_quarters,
  # ADS_quarter_counter = ADS_quarter_counter)
  
  # Call JAGS Function
  jags_output_list[[variant]] <- jags(data, 
                                      inits, 
                                      parameters, 
                                      "removal_model_simple.jags", 
                                      n.chains = nc, 
                                      n.thin = nt, 
                                      n.iter = ni,
                                      n.burnin = nb)
  
  print(paste0("iteration ", variant))
}

saveRDS(jags_output_list, file = here("Results" , "scenario_1", "jags_output_list_11-15.rds"))




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

## Saving and restoring the 50 run consistency check results and the input data
saveRDS(consist_test_output_jags, file = "consistency_test_outputs_simple_12.21.23.rds")
# saveRDS(erad_quarter_results, file = "consistency_test_input_data_12.21.23.rds")
# consist_test_output_jags <- readRDS("consistency_test_outputs_8.16.23.rds")
# output_jags <- readRDS("consistency_test_outputs_simple_12.21.23.rds")

mean_outputs <- list()
mean_outputs$N_sum <- as.data.frame(matrix(NA, nrow = 50, ncol = length(obs_quarters)))
colnames(mean_outputs$N_sum) <- paste0("Quarter_", obs_quarters)
mean_outputs$p_vis <- as.data.frame(matrix(NA, nrow = 50, ncol = 4))
colnames(mean_outputs$p_vis) <- size_class_names
mean_outputs$p_trap <- as.data.frame(matrix(NA, nrow = 50, ncol = 4))
colnames(mean_outputs$p_trap) <- size_class_names
mean_outputs$r_values <- as.data.frame(matrix(NA, nrow = 50, ncol = 4))
colnames(mean_outputs$r_values) <- size_class_names
for(iter in 1:50) {
  # Separating out encounter probability for visual survey vs trapping; since effort is the same throughout, the encounter probability for each size class is the same through time
  mean_outputs$p_vis[iter, ] <- consist_test_output_jags[[iter]]$p[,1,1] 
  mean_outputs$p_trap[iter, ] <- consist_test_output_jags[[iter]]$p[,2,2] # depends on trap days
  # Saving growth rates
  mean_outputs$r_values[iter, 1] <- consist_test_output_jags[[iter]]$r1
  mean_outputs$r_values[iter, 2] <- consist_test_output_jags[[iter]]$r2
  mean_outputs$r_values[iter, 3] <- consist_test_output_jags[[iter]]$r3
  mean_outputs$r_values[iter, 4] <- consist_test_output_jags[[iter]]$r4
  for(t in 1:Q) {
    mean_outputs$N_sum[iter, t] <- consist_test_output_jags[[iter]]$N.sum[t] # didn't actually sum all 4 years, just the first 3 (fixed for next set of runs)
  }
}

par(mfrow = c(2,3))
hist(mean_outputs$N_sum$Quarter_2)
abline(v = nrow(erad_quarter_results$quarter_timeseries[[3]]), col = "red")
hist(mean_outputs$N_sum$Quarter_3)
abline(v = nrow(erad_quarter_results$quarter_timeseries[[4]]), col = "red")
hist(mean_outputs$N_sum$Quarter_4, xlim = c(900, 1200))
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
size_class_growth <- as.data.frame(matrix(NA, nrow = length(obs_quarters), ncol = 4))
colnames(size_class_growth) <- size_class_names


for(quarter in 1:(length(obs_quarters)-1)) {
  # Separating out population at the beginning of the inter-primary period (after an eradication primary period)
  after_erad_pop <- erad_quarter_results$all_quarters_before_after_erad$pop_after_erad[obs_quarters[quarter]][[1]]
  # Separating out population at the end of the inter-primary period (before the next eradication primary period)
  before_erad_pop <- erad_quarter_results$all_quarters_before_after_erad$pop_before_erad[obs_quarters[quarter+1]][[1]]
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
    size_class_growth[quarter, size] <- nrow(size_before_pop)/nrow(size_after_pop)
  }
}



# > size_class_growth
# size_class inter-primary_1 inter-primary_2 inter-primary_3
# 1      small       1.0044188       0.9814815       1.0508876
# 2     medium       0.4240150       0.5806452       1.0708661
# 3      large       0.4386792       0.5328467       1.0204082
# 4     xlarge       0.4341232       0.4797297       0.9807692

# Since the number of days between quarters is the same for all inter-quarter periods, I can compare two
# histograms - the mean r values from all 50 interations with the 5 real inter-quarter growth periods
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")

p1 <- hist(mean_outputs$r_values$xlarge^57, plot = FALSE)
p2 <- hist(size_class_growth[,4], plot = FALSE)
plot(p2, col = c1)
plot(p1, col = c2, add = TRUE)

