################ Functions to format & validate model results ###################

### Function to classify SVL into size categories

# Adding size distribution to combined quarters dataframe
size_class_fun <- function(snake_SVL,
                           size_limits){
  if(snake_SVL <= size_limits[1,2]) {
    return("small")
  } else if(snake_SVL > size_limits[2,1] & snake_SVL <= size_limits[2,2]) {
    return("medium")
  } else if(snake_SVL > size_limits[3,1] & snake_SVL <= size_limits[3,2]) {
    return("large")
  } else if (snake_SVL > size_limits[4, 1]) {
    return("xlarge")
  } 
}


### Function to format quarter time series results for plotting by quarter
## Inputs:
# - quarter_results = list of dataframes with the time series results for quarters
results_format_fun <- function(quarter_results,
                               size_limits) {
  # Melt quarter timeseries into a single dataframe
  all_quarters <- melt(quarter_results, id.vars = c("ID", "sex", "repro_prob", "growth_quant"))
  colnames(all_quarters)[6:7] <- c("SVL", "Quarter")
  
  # Adding size distribution to combined quarters dataframe
  for(snake in 1:nrow(all_quarters)) {
    all_quarters$size_category[snake] <- size_class_fun(all_quarters$SVL[snake],
                                                        size_limits = size_limits)
  }
  # Factoring size category column for graphing
  all_quarters$size_category <- factor(all_quarters$size_category, 
                                       levels = c("small", "medium", "large", "xlarge"))
  
  # Return results 
  return(all_quarters)
}


### Function to create model validation summary stats


model_validation_fun <- function(quarter_timeseries,
                                 all_quarters,
                                 repro_females,
                                 total_quarters,
                                 size_classes) {
  # Create empty list to store results
  model_validation <- list()
  # Proportion of females reproduced each quarter
  model_validation$mom_rate <- vector()
  # Population growth rate each quarter
  model_validation$pop_growth_rate <- vector()
  # Mortality rate per quarter
  model_validation$mortality_rate <- vector()
  # Calculating the above rates for each quarter
  for(quarter in 1:total_quarters) {
    model_validation$mom_rate[quarter] <- nrow(repro_females[[quarter]])/nrow(quarter_timeseries[[quarter + 1]][quarter_timeseries[[quarter + 1]]$sex == "F", ])
    model_validation$pop_growth_rate[quarter] <- nrow(quarter_timeseries[[quarter+1]][!quarter_timeseries[[quarter+1]]$ID %in% quarter_timeseries[[quarter]]$ID,])/nrow(quarter_timeseries[[quarter+1]])
    model_validation$mortality_rate[quarter] <- nrow(quarter_timeseries[[quarter]][!quarter_timeseries[[quarter]]$ID %in% quarter_timeseries[[quarter + 1]]$ID,])/nrow(quarter_timeseries[[quarter]])
  }
  # Lifespan of individuals born within the simulation
  model_validation$lifespan <- vector()
  # Maximum size each snake reached
  model_validation$max_size <- vector()
  # Isolating snakes born in the simulation
  sim_born_snakes <- all_quarters[-grep("A", all_quarters$ID),]
  # Number of quarters that snakes stay in each size class
  age_by_size <- as.data.frame(matrix(NA, nrow = length(unique(sim_born_snakes$ID)), 
                                      ncol = 5))
  colnames(age_by_size) <- c("ID", "num_Q_small", "num_Q_medium", "num_Q_large",
                             "num_Q_xlarge")
  # Calculating above values for all simulation born snakes
  for(snake in 1:length(unique(sim_born_snakes$ID))) {
    model_validation$lifespan[snake] <- nrow(sim_born_snakes[sim_born_snakes$ID == sim_born_snakes$ID[snake],])
    model_validation$max_size[snake] <- max(sim_born_snakes$SVL[sim_born_snakes$ID == sim_born_snakes$ID[snake]])
    snake_data <- sim_born_snakes[sim_born_snakes$ID == unique(sim_born_snakes$ID)[snake], ] 
    age_by_size[snake, 1] <- snake_data$ID[1]
    for(size in 1:length(size_classes)) {
      age_by_size[snake, (size + 1)] <- nrow(snake_data[snake_data$size_category == size_classes[size],])
    }
  }
  # Average number of quarters that simulation born snakes spend in each size class
  model_validation$average_age_by_size <- vector()
  for(size in 1:length(size_classes)) {
    model_validation$average_age_by_size[size] <- mean(age_by_size[, (size + 1)])
  }
  # Proportion of snakes that lived longer than 15 years
  model_validation$prop_over_60_quart <- length(model_validation$lifespan[model_validation$lifespan > 60])/length(model_validation$lifespan)
  # Isolating snakes who live longer than 10 years to look at their growth trajectories
  model_validation$LL_snakes <- all_quarters[0,]
  for(snake in 1:length(unique(all_quarters$ID))) {
    snake_ts <- all_quarters[all_quarters$ID == unique(all_quarters$ID)[snake], ]
    if(nrow(snake_ts) > 40) {
      model_validation$LL_snakes <- rbind(model_validation$LL_snakes, snake_ts)
    } 
  }
  # Proportion of entire population that lives more than 10 years
  model_validation$LL_snakes_prop <- length(unique(model_validation$LL_snakes$ID))/length(unique(all_quarters$ID))
  # Sampling 25 long lived snakes and graphing their growth
  if(length(unique(model_validation$LL_snakes$ID)) >= 25) {
    plot_snake_IDs <- sample(unique(model_validation$LL_snakes$ID), 25)
    plot_snakes <- all_quarters[all_quarters$ID %in% plot_snake_IDs,]
    model_validation$LL_plot <- ggplot(plot_snakes, aes(x = Quarter, y = SVL, color = ID)) +
      geom_path()
  }
  # Summarizing many of the above states to erport mean, min and max values
  mv_table <- as.data.frame(matrix(NA, nrow = 3, ncol = 5))
  colnames(mv_table) <- c("mom_rate_per_quarter", "pop_growth_rate_per_quarter", 
                          "mortality_rate_per_quarter", "lifespan_in_quarters", "max_SVL")
  rownames(mv_table) <- c("min", "mean", "max")
  mv_table$mom_rate_per_quarter <- summary(model_validation$mom_rate)[c(1,4,6)]
  mv_table$pop_growth_rate_per_quarter <- summary(model_validation$pop_growth_rate)[c(1,4,6)]
  mv_table$mortality_rate_per_quarter <- summary(model_validation$mortality_rate)[c(1,4,6)]
  mv_table$lifespan_in_quarters <- summary(model_validation$lifespan)[c(1,4,6)]
  mv_table$max_SVL <- summary(model_validation$max_size)[c(1,4,6)]
  # Returning both full results and summarized results
  return(list(model_validation = model_validation, model_val_summary = mv_table))
}


                          
#### Transforming observed dead snakes into summary dataframes (# of snakes) for inputting
#### into the estimation model

all_observions_fun <- function(all_observed,
                             erad_days, # if more than one method, use c() in argument
                             erad_quarters, # if more than one method, use c() in argument
                             methods) {
  all_observations <- list()
  all_quarters <- sort(unique(erad_quarters))
  all_days <- sort(unique(erad_days))
  # Count the number of observed snake for each quarter & day that effort occurred
  for(quarter in all_quarters) {
    all_observations[[quarter]] <- as.data.frame(matrix(NA, nrow = 0, ncol = 2))
    colnames(all_observations[[quarter]] <- c("num_snakes", "day")
    for(day in all_days) {
      all_observations[[quarter]]$num_snakes <- 
  
  
  names(all_observations) <- paste0("Quarter_", )
  
}







