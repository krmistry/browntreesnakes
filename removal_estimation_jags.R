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
load(here("est_model_input_data.RData"))
# erad_quarter_results <- input_data$erad_quarter_results
# erad_methods <- input_data$erad_methods
# erad_days <- input_data$erad_days
# erad_quarters <- input_data$erad_quarters
# size_class_names <- input_data$size_class_names
# all_observed <- input_data$all_observed

# Reformatting removals and effort data to the format needed by estimation model
erad_reformatted <- all_observations_fun(erad_results_ts = erad_quarter_results,
                                         methods = erad_methods[c(2:3)])
removals_array <- erad_reformatted$observation
effort_array <- erad_reformatted$effort
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
  N[1,k,1,1] <-  N.base[1,k,1,1] + miss[1,k]
}

# Parameter priors
beta.p ~ dunif(0, 100)
alpha.p[1] ~ dunif(-100,100)
alpha.p[2] ~ dunif(-100,100)
s1 ~ dunif(0.9980, 0.9999)
s2 ~ dunif(0.9988, 0.9999)
s3 ~ dunif(0.9989, 0.99999)
s4 ~ dunif(0.9989, 0.99999)
f2 ~ dgamma(1,0.3)
f3 ~ dgamma(1,0.3)
f4 ~ dgamma(1,0.3)

# Transition matrix (used for population size class growth & reproduction between primary sampling periods)
for(t in 1:(T-1)){
  P[1,1,t] <- s1^days_btwn[t]
  P[2,1,t] <- s1^days_btwn[t]
  P[2,2,t] <- s2^days_btwn[t]
  P[3,2,t] <- s2^days_btwn[t]
  P[3,3,t] <- s3^days_btwn[t]
  P[4,3,t] <- s3^days_btwn[t]
  P[4,4,t] <- s4^days_btwn[t]
  P[1,2,t] <- f2
  P[1,3,t] <- f3
  P[1,4,t] <- f4
  P[2,3,t] <- 0
  P[2,4,t] <- 0
  P[3,1,t] <- 0
  P[3,4,t] <- 0
  P[4,1,t] <- 0
  P[4,2,t] <- 0
}


for(t in 1:T) { # start primary sampling period loop
  for(i in 1:I) { # start secondary sampling instances loop
    for(j in 1:J) { # start eradication method loop
      # Calculate encounter probability for each method, secondary sampling instance and primary sampling period
      logit(p[j,i,t]) <- alpha.p[j] + beta.p * log(xi[j,i,t])
    } # end eradication method loop
    
    for(k in 1:K) { # start size class loop
       for(j in 1:J) { # start eradication method loop
          # Calculate removals based on encounter probability and N
          Y[j,k,i,t] ~ dbin(p[j,i,t], N[j,k,i,t])
          #Y[j,k,i,t] ~ dbin(0.8, 1000) # for trouble-shooting
       } # end eradicaiton method loop
      
        # Calculate first row of N for j = 1 when i > 1 using i = 1 and j = 2 (J) 
        N[1,k,i+1,t] <- N[J,k,i,t] - Y[J,k,i,t]
      
        for(j in 2:J) { # start eradication method loop for j > 1
            # Calculate N using past method step N minus removals
            N[j,k,i,t] <- N[j-1,k,i,t] - Y[j-1,k,i,t]
            #N[j,k,i,t] <- 2000 - Y[j-1,k,i,t] # for trouble-shooting
        } # end eradicaiton loop for j > 1
      } # end size class loop
  } # end secondary sampling instance loop
  for(k in 1:K) { # start size class loop
    R[1,k,t] <- N[J,k,I,t] - Y[J,k,I,t]
  } # end size class loop
} # end primary sampling period loop

for(t in 1:(T-1)){ # start primary sampling period loop
  # Calculate population at beginning of primary sampling period using remaining population from the end of previous sampling period X transition matrix
    D[1,1,t] ~ dpois(R[1,1,t]*P[1,1,t] + R[1,2,t]*P[1,2,t] + R[1,3,t]*P[1,3,t] + R[1,4,t]*P[1,4,t])
    D[1,2,t] ~ dpois(R[1,1,t]*P[2,1,t] + R[1,2,t]*P[2,2,t])
    D[1,3,t] ~ dpois(R[1,2,t]*P[3,2,t] + R[1,3,t]*P[3,3,t])
    D[1,4,t] ~ dpois(R[1,3,t]*P[4,3,t] + R[1,4,t]*P[4,4,t])
    for(k in 1:K){ # start size class loop
      N[1,k,1,t+1] <- D[1,k,t]
      #N[1,k,1,t+1] <- 1000 # for trouble-shooting
    } # end size class loop
} # end primary sampling period loop

for(t in 1:T){ # start primary sampling period loop
  N.sum[t] <- sum(N[1,,1,t])
} # end primary sampling period loop


} # end model
", fill= TRUE)
sink()



# Initial values for N (i.e. N[1, k, 1, 1])
Y <- removals_array
N.base <- array(NA_real_, dim = c(J, K, I, T))
N.base[1,1,1,1] <- sum(Y[,1,,1:3])
N.base[1,2,1,1] <- sum(Y[,2,,1:3])
N.base[1,3,1,1] <- sum(Y[,3,,1:3])
N.base[1,4,1,1] <- sum(Y[,4,,1:3])


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
       beta.p = runif(1,0,10),
       alpha.p = runif(2,-10,10))
  
}

# Bundle data together
data <- list(Y = removals_array,
             J = 2, K = 4, I = 14, T = 4, 
             xi = effort_array, 
             days_btwn = date_diff, 
             N.base = N.base)

# Parameters monitored
parameters <- c("N", "N.sum", "alpha.p", "beta.p", "s1", "s2", "s3", "s4", 
                "f2", "f3", "f4", "p",  "R", "D")

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

# Checking traceplots for Rhat > 1.05 (per jagsUI manual)
traceplot(output_jags, Rhat_min = 1.05)

## Plotting mean values for N by size category to get quick loop at model outputs
mean_N <- as.data.frame(matrix(0, nrow = T, ncol = K))
colnames(mean_N) <- size_class_names
mean_N$Quarter <- sort(unique(unlist(erad_quarters)))
for(size in 1:K) {
  for(quarter in 1:T) {
    mean_N[quarter, size] <- output_jags$mean$N[2,size,1,quarter]
  }
}

# Making data long for plotting
mean_N_long <- melt(mean_N, id.vars = "Quarter")
colnames(mean_N_long)[2:3] <- c("size_class", "N")
# Plotting size class (mean) estimated abundance through time
ggplot(mean_N_long, aes(x = Quarter, y = N, color = size_class)) +
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
ggplot(rbind(mean_N_long, real_N_long, observed_N_long), aes(x = Quarter, y = N, color = source)) +
  geom_line(size = 2) +
  facet_wrap("size_class", scales = "free_y") +
  theme_bw()


  
