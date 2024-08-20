## Set up common to all strategies

library(reshape2)
library(ggplot2)
library(dplyr)
library(here)
library(tictoc)
library(fdrtool)
library(jagsUI)
library(scales)
library(vctrs)

source(here("Scripts/00_user_inputs.R"))
source(here("Scripts/01_model_functions.R"))
source(here("Scripts/02_results_functions.R"))
source(here("estimation_model_functions.R"))

##### Shared parameters 
# Target area
area_size <- 50
# Carrying capacity
K <- 60*area_size
# Growth probability (p_g)
g_density_prob <- 0.75 
# Number of days in a quarter
day_time_step <- 91

# MCMC settings for estimation model (use the same for all)
nt <- 1
nb <- 10000 # 20000 
ni <- 30000 + nb # 50000
nc <- 3

# Three different mean starting densities
starting_density <- c("low" = 11, "medium" = 22, "high" = 44)
starting_pop <- list()

set.seed(818) # To keep the starting values consistent
for(i in names(starting_density)) {
  starting_pop[[i]] <- round(rnorm(100, starting_density[i]*area_size, 
                                   starting_density[i]*area_size*0.1))
}

# Two starting size distributions (could add more if there's time to run them)
starting_size_dist <- list()
starting_size_dist$more_small <- c(0.42, rep((1-0.42)/3,3))
starting_size_dist$more_xlarge <- c(rep((1-0.42)/3,3), 0.42)

# Number of variants for each permutation
num_variants <- 50

# JAGS Removal Estimation Model - simple growth with zero inflation
sink("removal_model_alt_strategies.jags")
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