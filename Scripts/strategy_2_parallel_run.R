### Running strategy 2 in parallel

# Setting up parallel clusters
library(doParallel)
# Detect the number of clusters available
n_cores <- detectCores()
# Select half of them - broke, so trying fewer cores
cl <- makeCluster(n_cores/2, outfile = "")
registerDoParallel(cl)


# Using parallel function:
# Parallel test
results <- foreach(variant = 1:50)  %dopar% {
  library(here)
  source(here("strategy_2_parallel_function.R"))
  strat_2_parallel_fun(P = 2, 
                       D = 1, 
                       final_time_step = 20,
                       variant = variant)
}

# Stop the cluster
stopCluster(cl = cl)