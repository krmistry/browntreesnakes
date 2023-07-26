library(jagsUI)
library(here)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tictoc)
library(fdrtool)

## Loading objects and functions
source(here("Scripts/00_user_inputs.R"))
source(here("Scripts/01_model_functions.R"))
source(here("Scripts/02_results_functions.R"))

# Reading in data generated on 7/7/23
# load(here("est_model_input_data.RData"))
# erad_quarter_results <- input_data$erad_quarter_results
# erad_methods <- input_data$erad_methods
# erad_days <- input_data$erad_days
# erad_quarters <- input_data$erad_quarters
# size_class_names <- input_data$size_class_names
# all_observed <- input_data$all_observed

############################ Estimation model version with multiple data sources (both visual and trap) ##########################

# Reformatting removals and effort data to the format needed by estimation model
erad_reformatted <- all_observations_fun(erad_results_ts = erad_quarter_results,
                                         methods = erad_methods[c(2:3)])
removals_array <- erad_reformatted$observation
effort_array <- erad_reformatted$effort
# Standardized version of effort_array
effort_array_2 <- array(NA, dim = dim(effort_array))
effort_array_2[1,,] <- effort_array[1,,]/12
effort_array_2[2,,] <- effort_array[2,,]/12
#effort_array[2,,] <- effort_array[2,,] + 1

# Values needed for array dimensions & loops
J <- 2 # number of eradication methods (vision & trap)
K <- 4 # number of size classes
I <- 14 # Secondary sampling periods (days within quarters)
T <- 4 # Primary sampling periods (quarters)

# Days between primary sampling periods (quarters 2, 3, 6 and 7)
date_diff <- vector()
for(t in 1:(T-1)) {
  date_diff[t] <- ((erad_quarters[[2]][t+1]-1)*91 + erad_days[[2]][1]) - 
    ((erad_quarters[[2]][t]-1)*91 + max(unlist(erad_days[c(2:3)])))
}
  



# JAGS Removal Estimation Model 
sink("removal_model.jags")
cat("
model {

# Set up first row of N
for(k in 1:K) {
  p.miss[1,k,1:500] <- rep(1/500,500)
  miss[1,k] ~ dcat(p.miss[1,k,1:500])
  N[k,1,1] <-  N.base[k,1,1] + miss[1,k]
}

# Parameter priors
beta.p[1] ~ dunif(0, 100) # encounter slope method 1
beta.p[2] ~ dunif(0, 100) # encounter slope method 2
alpha.p[1] ~ dunif(-100,100) # encounter intercept method 1
alpha.p[2] ~ dunif(-100,100) # encounter intercept method 2
s1 ~ dunif(0.9, 0.9999) # small survival 
s2 ~ dunif(0.9, 0.9999) # medium survival 
s3 ~ dunif(0.9, 0.99999) # large survival
s4 ~ dunif(0.9, 0.99999) # x-large survival 
f2 ~ dgamma(1,0.3) # medium fecundity
f3 ~ dgamma(1,0.3) # large fecundity
f4 ~ dgamma(1,0.3) # x-large fecundity
t1 ~ dunif(0.001, 0.999) # transition small -> medium
t2 ~ dunif(0.001, 0.999) # transition medium -> large
t3 ~ dunif(0.001, 0.999) # transition large -> x-large


# Transition matrix (used for population size class growth & reproduction between primary sampling periods)
for(t in 1:(T-1)){
  P[1,1,t] <- s1^days_btwn[t]*(1-t1^days_btwn[t])
  P[1,2,t] <- f2
  P[1,3,t] <- f3
  P[1,4,t] <- f4
  P[2,1,t] <- s1^days_btwn[t]*t1^days_btwn[t]
  P[2,2,t] <- s2^days_btwn[t]*(1-t2^days_btwn[t])
  P[2,3,t] <- 0
  P[2,4,t] <- 0
  P[3,1,t] <- 0
  P[3,2,t] <- s2^days_btwn[t]*t2^days_btwn[t]
  P[3,3,t] <- s3^days_btwn[t]*(1-t3^days_btwn[t])
  P[3,4,t] <- 0
  P[4,1,t] <- 0
  P[4,2,t] <- 0
  P[4,3,t] <- s3^days_btwn[t]*t3^days_btwn[t]
  P[4,4,t] <- s4^days_btwn[t]
}


for(t in 1:T) { # start primary sampling period loop
  for(i in 2:I+1) { # start secondary sampling instances loop
    for(j in 1:J) { # start eradication method loop
      # Calculate encounter probability for each method, secondary sampling instance and primary sampling period
      logit(p[j,i-1,t]) <- alpha.p[j] + beta.p[j] * log(xi[j,i-1,t])
    } # end eradication method loop
    
    for(k in 1:K) { # start size class loop
      # Create array to hold within i population, to create Y for each j method
      M[1,k] <- N[k,i-1,t]
    } # end size class loop
    
    for(j in 1:J) { # start eradication method loop
       for(k in 1:K) { # start size class loop
          # Calculate removals based on encounter probability and M (population in within i instance)
          Y[j,k,i-1,t] ~ dbin(p[j,i-1,t],M[j,k])
          # Update M population in the next row to do the next removal method
           M[j+1,k] <- M[j,k] - Y[j,k,i-1,t]
       } # end size class loop
      } # end eradication method loop
      
    for(k in 1:K) { # start size class loop
        # Sum removals in each i for each size class
        Y.sum[k,i-1,t] <- sum(Y[,k,i-1,t])
        # Calculate N using last time step N minus summed removals
        N[k,i,t] <- N[k,i-1,t] - Y.sum[k,i-1,t]
    } # end size class loop
  } # end secondary sampling instance loop
  
  for(k in 1:K) { # start size class loop
    # Calculate remaining population at the end of the primary sampling period
    R[k,t] <- N[k,I,t] - Y.sum[k,I-1,t]
  } # end size class loop
} # end primary sampling period loop

for(t in 1:(T-1)){ # start primary sampling period loop
  # Calculate population at beginning of primary sampling period using remaining population from the end of previous sampling period X transition matrix
    D[1,t] ~ dpois(R[1,t]*P[1,1,t] + R[2,t]*P[1,2,t] + R[3,t]*P[1,3,t] + R[4,t]*P[1,4,t])
    D[2,t] ~ dpois(R[1,t]*P[2,1,t] + R[2,t]*P[2,2,t])
    D[3,t] ~ dpois(R[2,t]*P[3,2,t] + R[3,t]*P[3,3,t])
    D[4,t] ~ dpois(R[3,t]*P[4,3,t] + R[4,t]*P[4,4,t])
    for(k in 1:K){ # start size class loop
      N[k,1,t+1] <- D[k,t]
    } # end size class loop
} # end primary sampling period loop

for(t in 1:T){ # start primary sampling period loop
  N.sum[t] <- sum(N[,1,t])
} # end primary sampling period loop


} # end model
", fill= TRUE)
sink()


# #### Version of jags model above in R language, to help with troubleshooting
# # Set up first row of N
# N <- array(NA, dim = c(K,I+1,T))
# miss <- round(runif(K,1,500))
# N[,1,1] <-  N.base[1,,1,1] + miss
# # Parameter priors (1 draw from each)
# beta.p <- vector()
# beta.p[1] <- runif(1, 0, 10) # encounter slope method 1
# beta.p[2] <- runif(1, 0, 10) # encounter slope method 2
# alpha.p <- vector()
# alpha.p[1] <- runif(1, -10,10) # encounter intercept method 1
# alpha.p[2] <- runif(1, -10,10) # encounter intercept method 2
# s1 <- runif(1, 0.9, 0.9999) # small survival 
# s2 <- runif(1, 0.9, 0.9999) # medium survival 
# s3 <- runif(1, 0.9, 0.99999) # large survival
# s4 <- runif(1, 0.9, 0.99999) # x-large survival 
# f2 <- rgamma(1, 1,0.3) # medium fecundity
# f3 <- rgamma(1, 1,0.3) # large fecundity
# f4 <- rgamma(1, 1,0.3) # x-large fecundity
# t1 <- runif(1, 0.001, 0.999) # transition small -> medium
# t2 <- runif(1, 0.001, 0.999) # transition medium -> large
# t3 <- runif(1, 0.001, 0.999) # transition large -> x-large
# # Transition matrix (used for population size class growth & reproduction between primary sampling periods)
# P <- array(NA, dim = c(K,K,T-1))
# for(t in 1:(T-1)){
#   P[1,1,t] <- s1^days_btwn[t]*(1-t1^days_btwn[t])
#   P[1,2,t] <- f2
#   P[1,3,t] <- f3
#   P[1,4,t] <- f4
#   P[2,1,t] <- s1^days_btwn[t]*t1^days_btwn[t]
#   P[2,2,t] <- s2^days_btwn[t]*(1-t2^days_btwn[t])
#   P[2,3,t] <- 0
#   P[2,4,t] <- 0
#   P[3,1,t] <- 0
#   P[3,2,t] <- s2^days_btwn[t]*t2^days_btwn[t]
#   P[3,3,t] <- s3^days_btwn[t]*(1-t3^days_btwn[t])
#   P[3,4,t] <- 0
#   P[4,1,t] <- 0
#   P[4,2,t] <- 0
#   P[4,3,t] <- s3^days_btwn[t]*t3^days_btwn[t]
#   P[4,4,t] <- s4^days_btwn[t]
# }
# # Generate N, Y and R values
# library(boot)
# p <- array(NA, dim = c(J,I,T))
# for(t in 1:T) { # start primary sampling period loop
#   for(i in 1:I) { # start secondary sampling instances loop
#     for(j in 1:J) { # start eradication method loop
#       # Calculate encounter probability for each method, secondary sampling instance and primary sampling period
#       p[j,i,t] <- inv.logit(alpha.p[j] + beta.p[j] * log(effort_array_2[j,i,t]))
#     } # end eradication method loop
#   }
# }
# M <- array(NA, dim = c(J+1,K))
# Y.sum <- array(NA, dim = c(K, I, T))
# for(i in 1:I) {
#   M[1,] <- N[,i,1]
# for (j in 1:J) { # start eradication method loop
#   for (k in 1:K) { # start size class loop
#     # Calculate removals based on encounter probability and M (population in within i instance)
#     Y[j, k, i, 1] <- rbinom(1, M[j, k], p[j, i, 1])
#     M[j + 1, k] <- M[j, k] - Y[j, k, i, 1]
#   } # end eradication method loop
# }
# 
# for (k in 1:K) {
#     # Sum removals in each i for each size class
#     Y.sum[k, i, 1] <- sum(Y[, k, i, 1])
#     # Calculate N using last time step N minus summed removals
#     N[k, i+1, 1] <- N[k, i, 1] - Y.sum[k, i, 1]
#     #N[j,k,i,t] <- 2000 - Y[j-1,k,i,t] # for trouble-shooting
#   } # end size class loop
# } # end secondary sampling instance loop
# R <- array(dim = c(K,T))
#   #for(k in 1:K) { # start size class loop
#     R[,1] <- N[,I+1,1] - Y.sum[,I,1]
#   #} # end size class loop
# #} # end primary sampling period loop

##################### Creating all inputs for jags model ###########
# Initial values for N (i.e. N[1, k, 1, 1])
Y <- removals_array
N.base <- array(NA_real_, dim = c(K, I, T))
N.base[1,1,1] <- sum(Y[,1,,1:T])
N.base[2,1,1] <- sum(Y[,2,,1:T])
N.base[3,1,1] <- sum(Y[,3,,1:T])
N.base[4,1,1] <- sum(Y[,4,,1:T])


# initialize D to > than the number that will be removed in the following year
Y.remove1 <- vector()
Y.remove2 <- vector()
Y.remove3 <- vector()
Y.remove4 <- vector()
for(t in 1:T){
  Y.remove1[t] <- sum(Y[,1,,t])
  Y.remove2[t] <- sum(Y[,2,,t])
  Y.remove3[t] <- sum(Y[,3,,t])
  Y.remove4[t] <- sum(Y[,4,,t])
}

D.init <- array(NA,dim = c(1,K,(T-1)))
for(t in 1:(T-1)){
  D.init[1,1,t] <- Y.remove1[t+1] + 1
  D.init[1,2,t] <- Y.remove2[t+1] + 1
  D.init[1,3,t] <- Y.remove3[t+1] + 1
  D.init[1,4,t] <- Y.remove4[t+1] + 1
}

# Initial values for select parameters
inits <- function (){
  list(D = D.init,
       beta.p = runif(2,0,10),
       alpha.p = runif(2,-10,10))
  
}

# Bundle data together
data <- list(Y = removals_array,
             J = 2, K = 4, I = 14, T = 1, 
             xi = effort_array, 
             days_btwn = date_diff, 
             N.base = N.base)

# Parameters monitored
parameters <- c("N", "N.sum", "alpha.p", "beta.p", "s1", "s2", "s3", "s4", 
                "f2", "f3", "f4", "t1", "t2", "t3", "p",  "R", "D", "M", "Y.sum")

# MCMC settings
ni <- 20000
nt <- 1
nb <- 10000
nc <- 3

# Call JAGS Function
output_jags <- jags(data, 
                    inits, 
                    parameters, 
                    "removal_model.jags", 
                    n.chains = nc, 
                    n.thin = nt, 
                    n.iter = ni,
                    n.burnin = nb)


############################ Estimation model version with one data source (either visual or trap) ##########################
# In this version, j is not required, as there is only one method, so it also requires modifying the input data sets
# Also, p is size-dependent (not method, of course)

# Using only visual data
erad_reformatted_v2 <- all_observations_fun(erad_results_ts = erad_quarter_results,
                                         methods = erad_methods[2])

removals_array_v2 <- erad_reformatted_v2$observation[1,,,] # array(dim = c(4,14,4))
effort_array_v2 <- erad_reformatted_v2$effort[1,,] # array(dim = c(14,4))

# Values needed for array dimensions & loops
K <- 4 # number of size classes
I <- 14 # Secondary sampling periods (days within quarters)
Q <- 4 # Primary sampling periods (quarters)

# Days between primary sampling periods when method is used
# For visual data
date_diff <- vector()
for(t in 1:(Q-1)) {
  date_diff[t] <- ((erad_quarters[[2]][t+1]-1)*91 + erad_days[[2]][1]) - 
    ((erad_quarters[[2]][t]-1)*91 + max(unlist(erad_days[2])))
}


# JAGS Removal Estimation Model - version 2
sink("removal_model_v2.jags")
cat("
model {

# Set up first row of N
for(k in 1:K) {
  p.miss[1,k,1:1000] <- rep(1/1000,1000)
  miss[1,k] ~ dcat(p.miss[1,k,1:1000])
  N[k,1,1] <-  N.base[k,1,1] + miss[1,k]
}

# Parameter priors
beta.p[1] ~ dunif(0, 10) # encounter slope small
beta.p[2] ~ dunif(0, 10) # encounter slope medium
beta.p[3] ~ dunif(0, 10) # encounter slope for large
beta.p[4] ~ dunif(0, 10) # encounter slope for x-large
alpha.p[1] ~ dunif(-10,10) # encounter intercept small
alpha.p[2] ~ dunif(-10,10) # encounter intercept medium
alpha.p[3] ~ dunif(-10,10) # encounter intercept large
alpha.p[4] ~ dunif(-10, 10) # encounter slope x-large
s1 ~ dunif(0.9, 0.9999) # small survival 
s2 ~ dunif(0.9, 0.9999) # medium survival 
s3 ~ dunif(0.9, 0.99999) # large survival
s4 ~ dunif(0.9, 0.99999) # x-large survival 
f2 ~ dgamma(1,0.3) # medium fecundity
f3 ~ dgamma(1,0.3) # large fecundity
f4 ~ dgamma(1,0.3) # x-large fecundity
t1 ~ dunif(0.001, 0.999) # transition small -> medium
t2 ~ dunif(0.001, 0.999) # transition medium -> large
t3 ~ dunif(0.001, 0.999) # transition large -> x-large


# Transition matrix (used for population size class growth & reproduction between primary sampling periods)
for(t in 1:(Q-1)){
  P[1,1,t] <- s1^days_btwn[t]*(1-t1^days_btwn[t])
  P[1,2,t] <- f2
  P[1,3,t] <- f3
  P[1,4,t] <- f4
  P[2,1,t] <- s1^days_btwn[t]*t1^days_btwn[t]
  P[2,2,t] <- s2^days_btwn[t]*(1-t2^days_btwn[t])
  P[2,3,t] <- 0
  P[2,4,t] <- 0
  P[3,1,t] <- 0
  P[3,2,t] <- s2^days_btwn[t]*t2^days_btwn[t]
  P[3,3,t] <- s3^days_btwn[t]*(1-t3^days_btwn[t])
  P[3,4,t] <- 0
  P[4,1,t] <- 0
  P[4,2,t] <- 0
  P[4,3,t] <- s3^days_btwn[t]*t3^days_btwn[t]
  P[4,4,t] <- s4^days_btwn[t]
}


for(t in 1:Q) { # start primary sampling period loop  
  for(k in 1:K) { # start size class loop
      for(i in 1:I) { # start secondary sampling instances loop
      # Calculate encounter probability for each method, secondary sampling instance and primary sampling period
        logit(p[k,i,t]) <- alpha.p[k] + beta.p[k] * log(xi[i,t]) # effort is not size-dependent
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
    for(k in 1:K){ # start size class loop
      N[k,1,t+1] <- D[k,t]
    } # end size class loop
  # Summing all size classes into a single N for each primary sampling period
  N.sum[t] <- sum(N[,1,t])
} # end primary sampling period loop

} # end model
", fill= TRUE)
sink()

##################### Creating all inputs for jags model ###########
# Initial values for N (i.e. N[k, 1, 1])
Y <- removals_array_v2
N.base <- array(NA_real_, dim = c(K, I, Q))
for(k in 1:K) {
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

D.init <- array(NA,dim = c(K,(Q-1)))
for(t in 1:(Q-1)){
  D.init[1,t] <- Y.remove1[t+1] + 1
  D.init[2,t] <- Y.remove2[t+1] + 1
  D.init[3,t] <- Y.remove3[t+1] + 1
  D.init[4,t] <- Y.remove4[t+1] + 1
}

# Initial values for select parameters
inits <- function (){
  list(D = D.init,
       beta.p = runif(4,0,10),
       alpha.p = runif(4,-10,10))
  
}

# Bundle data together
data <- list(Y = removals_array_v2,
             K = 4, I = 14, Q = 4, 
             xi = effort_array_v2, 
             days_btwn = date_diff, 
             N.base = N.base)

# Parameters monitored
parameters <- c("N", "N.sum", "alpha.p", "beta.p", "s1", "s2", "s3", "s4", 
                "f2", "f3", "f4", "t1", "t2", "t3", "p",  "R", "D")

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

######################## Graphing jags outputs #############################

# Checking traceplots for Rhat > 1.05 (per jagsUI manual)
traceplot(output_jags, Rhat_min = 1.05)

## Traceplots of estimated parameters
traceplot(output_jags, parameters = c("beta.p[1]", 
                                      "beta.p[2]",
                                      "beta.p[3]",
                                      "beta.p[4]",
                                      "alpha.p[1]",
                                      "alpha.p[2]",
                                      "alpha.p[3]",
                                      "alpha.p[4]",
                                      "s1",
                                      "s2",
                                      "s3",
                                      "s4",
                                      "f2",
                                      "f3",
                                      "f4",
                                      "t1",
                                      "t2",
                                      "t3"))

## Plotting mean values for N by size category to get quick loop at model outputs
mean_N <- as.data.frame(matrix(0, nrow = Q, ncol = K))
colnames(mean_N) <- size_class_names
mean_N$Quarter <- sort(unique(unlist(erad_quarters)))
for(size in 1:K) {
  for(quarter in 1:Q) {
    mean_N[quarter, size] <- output_jags$mean$N[size,I,quarter]
  }
}

# Making data long for plotting
mean_N_long <- melt(mean_N, id.vars = "Quarter")
colnames(mean_N_long)[2:3] <- c("size_class", "N")
# Plotting size class (mean) estimated abundance through time
mean_N_plot_1 <- ggplot(mean_N_long, aes(x = Quarter, y = N, color = size_class)) +
  geom_line() +
  geom_hline(yintercept = K) +
  theme_bw()
  

## Adding total # of captured snakes and total # of alive simulated snakes in each size
## class in each quarter to graph against the predicted abundance
mean_N_long$source <- "estimated"
observed_erad_quarters <- sort(unique(unlist(erad_quarters[c(2,3)])))
real_N <- mean_N
observed_N <-  mean_N
for(size in size_class_names) {
  for(quarter in 1:length(observed_erad_quarters)) {
    observed_N[quarter,size] <- nrow(all_observed[all_observed$quarter == observed_erad_quarters[quarter] & all_observed$size_category == size,])
    real_N[quarter, size] <- nrow(erad_quarter_results$all_quarters[erad_quarter_results$all_quarters$Quarter == (observed_erad_quarters[quarter]+1) & erad_quarter_results$all_quarters$size_category == size,]) 
  }
}

real_N_long <- melt(real_N, id.vars = "Quarter")
colnames(real_N_long)[2:3] <- c("size_class", "N")
observed_N_long <- melt(observed_N, id.vars = "Quarter")
colnames(observed_N_long)[2:3] <- c("size_class", "N")
real_N_long$source <- "simulated"
observed_N_long$source <- "removed"

# Combining simulated, removed and estimated for each size class in each quarter
est_vs_sim_plot_1 <- ggplot(rbind(mean_N_long, real_N_long), aes(x = Quarter, y = N, color = source)) +
  geom_line(size = 2) +
  facet_wrap("size_class", scales = "free_y") +
  theme_bw()
## Plotting total estimated N vs total simulation N
est_vs_sim_plot_2 <- ggplot(rbind(mean_N_long, real_N_long), aes(y = N, x = Quarter, fill = source)) +
  geom_col(position = "dodge") +
  theme_bw()



