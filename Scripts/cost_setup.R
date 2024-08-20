### Cost set up

#### Itemized costs based on USDA consultation - these values were generated in 1st quarter 
#### of 2024, found in 2022_03_04 Large Landscape BTS ADS Cost Estimator Tool_INTERNAL USE.xlsx 
#### (in Cost_data folder - will want to update periodically)

## Creating lists to store parameters used for multiple methods
num_teams <- list()
num_devices <- list()
transects_checked_per_hour <- list()
transect_hours_per_day <- list()
overhead_hours <- list()
per_device_cost <- list()
init_baits <- list()
maintain_equip_monthly <- list()
device_spacing <- list()
per_bait_cost <- list()
misc_init_equip <- list()

# Values of standard transect design - length is based on CP, and between transect width of 20 meters is 
# the current best practice. If actual transects differ from this, this can be used to scale costs based on 
# these transect values
standard_transect_length <- 220 # meters 
standard_between_transects <- 20 # meters

#### Costs shared by all transect-based methods
# - one time transect cutting for 1 transect X number of transects
# - transect maintenance at regular intervals (say 4 weeks), X number of transects
# - hourly rate for employees
init_transect_cost <- 300 # per transect cost to set up
maint_transect_cost <- 16.65 # every 4 weeks, per transect cost
person_cost_per_hour <- 33.31

##### Visual survey costs 
# - number of searchers per night X search nights
# - one time purchase of equipment (air rifles, headlamps, etc)
# - person hours for any overhead activities, outside of searching, in a timescale that's easy to multiply - daily or weekly
misc_init_equip$visual <- 1350 # headlamps and other one time purchases
transect_hours_per_day$visual <- 4
overhead_hours$visual <- 15 # cumulative hours, not per person

#### Trap costs
# - one time purchase of traps, perhaps with a percentage of replacement based on time
# - person hours to check traps & feed bait, performed every 2-3 days while traps are out
# - bait population maintenance, a flat rate per time unit (weekly would be good) - may not be entirely necessary, leaving this out for now
device_spacing$trap <- 20
per_device_cost$trap <- 100 # for USGS clamshell traps after modification; there is another option for $300 that doesn't require modification
per_bait_cost$trap <- 9.1
maintain_per_week_per_mouse <- 1.5
transect_hours_per_day$trap <- 6 # it takes ~0.1 hours to check a trap, so if less than 60 traps are being checked, this might need to be decreased, but otherwise a standard per day rate should be okay
misc_init_equip$trap <- 500 # Summed together all other one-time purchase equipment costs, like air rifles
overhead_hours$trap <- 15
maintain_equip_monthly$trap <- 100
transects_checked_per_hour$trap <- 1 # Average number of transects checked per hour

### Bait costs
# - one time purchase of bait tubes, with a percentage of replacement based on time
# - baits, to be replaced every 2-3 days 
# - person hours to replace baits, every 2-3 days
device_spacing$bait_tube <- 20
per_device_cost$bait_tube <- 10 
per_bait_cost$bait_tube <- 0.75 # Variable based on the number ordered 
overhead_hours$bait_tube <- 10
transect_hours_per_day$bait_tube <- 6 # how many hours per day it takes to visit traps 
maintain_equip_monthly$bait_tube <- 50
misc_init_equip$bait_tube <- 250 # Summed together all other one-time purchase equipment costs
transects_checked_per_hour$bait_tube <- 2 # Average number of transects checked per hour

##### ADS fixed parameters
# Standard parameters
baits_per_ha <- 120 # maximum rate, but also the recommended one to get recommended saturation, 1 bait per 9 meters
max_ha_per_day <- 167
max_bait_application_rate_hourly <- 6600
num_applications_per_treatment <- 2


# Crew costs
ground_crew_hourly_rate <- 34*3 + 65 # 3 technicians + supervisor
pilot_hourly_rate <- 66
leadership_hourly <- 188
ground_prep_days <- 2
# Maintenance costs & # of days
mechanic_service_daily_rate <- 900
heli_maintenance_hourly_rate <- 800
# 3 prep days, 2 contingency days, 2 required off days for the pilot to rest 
prep_days_per_treatment <- 3
contingency_days <- 2 
pilot_off_days <- 2 
# Travel costs
num_travelers <- 2 # pilot and mechanic 
travel_days <- 2
flight_cost <- 1750
per_day_per_traveler <- 310 # hotel, meal per diem, and rental car - minimum of 2 days
# Bait costs
# Data frame with the per bait cost based on the total number of baits needed
bait_cost_by_amount <- as.data.frame(matrix(NA, nrow = 10, ncol = 3))
colnames(bait_cost_by_amount) <- c("lower_limit", "upper_limit", "bait_cost")
bait_cost_by_amount$lower_limit <- c(0, 216000, 324000, 360000, 432000, 540000, 
                                     648000, 864000, 1080000, 1620000)
bait_cost_by_amount$upper_limit <- c((bait_cost_by_amount[-1,1]-1), Inf)
bait_cost_by_amount$bait_cost <- c(5.280, 4.279, 3.674, 3.509, 3.179, 2.805, 2.486, 
                                   2.178, 1.980, 1.760)

# Misc. costs
freezer_rental_annual <- 32000 # Annual cost from FY22



### Basic cost function for any combination of methods

cost_function <- function(methods,
                          erad_days,
                          erad_quarters,
                          area_size,
                          num_transects,
                          num_teams) {
  # Creating empty lists 
  method_set_up <- list()
  daily_cost <- list()
  weekly_cost <- list()
  monthly_cost <- list()
  days_per_quarter <- list()
  weeks_per_quarter <- list()
  total_quarters <- list()
  total_cost <- list()
  
  ## Calculate total number of days and quarters for each method
  for(method in methods) {
    days_per_quarter[[method]] <- list()
    weeks_per_quarter[[method]] <- list()
    total_quarters[[method]] <- list()
    for(quarter in c(1:length(erad_quarters[[method]]))) {
      days_per_quarter[[method]][[quarter]] <- length(erad_days[[quarter]][[method]])
      weeks_per_quarter[[method]][[quarter]] <- ceiling((erad_days[[quarter]][[method]][length(erad_days[[quarter]][[method]])]-erad_days[[quarter]][[method]][1])/7) 
      total_quarters[[method]]<- length(erad_quarters[[method]])
    }
  }
  
  # If any transect methods included, calculate transect-related costs
  if(erad_methods[2] %in% methods | erad_methods[3] %in% methods | erad_methods[4] %in% methods) {
    # Set up
    transect_set_up <- num_transects[["total"]]*init_transect_cost
    # Per month
    transect_maintenance <- maint_transect_cost*num_transects[["total"]]
    # Total months 
    transect_months <- (max(unlist(erad_quarters[methods])) -  min(unlist(erad_quarters[methods])))*3
    # Total cost
    total_cost$transects <- transect_set_up + transect_maintenance*transect_months
  }
  
  # Visual survey costs
  if(erad_methods[2] %in% methods) {
    # Set up
    method_set_up$visual <- num_teams[["visual"]]*2*misc_init_equip$visual
    # Per day
    daily_cost$visual <- num_teams[["visual"]]*2*transect_hours_per_day$visual*person_cost_per_hour
    # Per week
    weekly_cost$visual <- overhead_hours$visual*person_cost_per_hour
    # Total days and weeks per quarter
    visual_weeks <- 0
    visual_days <- 0
    for(quarter in erad_quarters$visual) {
      quarter_weeks <- round((max(erad_days[[quarter]]$visual) - min(erad_days[[quarter]]$visual))/7)
      quarter_days <- length(erad_days[[quarter]]$visual)
      # Adding each quarters together
      visual_weeks <- visual_weeks + quarter_weeks
      visual_days <- visual_days + quarter_days
    }
    # Total cost
    total_cost$visual <- method_set_up$visual + visual_days*daily_cost$visual + visual_weeks*weekly_cost$visual
  }
  
  # Trap costs
  if(erad_methods[3] %in% methods) {
    # Number of devices
    num_devices$trap <- num_transects[["total"]]*(standard_transect_length/device_spacing$trap)
    # Set up
    method_set_up$trap <- num_devices$trap*(per_device_cost$trap + per_bait_cost$trap) + misc_init_equip$trap
    # Per day
    daily_cost$trap <- num_teams[["trap"]]*2*transect_hours_per_day$trap*person_cost_per_hour
    # Per week
    weekly_cost$trap <- maintain_per_week_per_mouse*num_devices$trap + overhead_hours$trap*person_cost_per_hour
    # Total days and weeks per quarter
    trap_weeks <- 0
    trap_days <- 0
    for(quarter in erad_quarters$trap) {
      quarter_weeks <- round((max(erad_days[[quarter]]$trap) - min(erad_days[[quarter]]$trap))/7)
      quarter_days <- length(erad_days[[quarter]]$trap)
      # Adding each quarters together
      trap_weeks <- trap_weeks + quarter_weeks
      trap_days <- trap_days + quarter_days
    }
    # Total cost
    total_cost$trap <- method_set_up$trap + trap_days*daily_cost$trap + trap_weeks*weekly_cost$trap
  }
  
  # Bait tube costs
  if(erad_methods[4] %in% methods) {
    # Number of devices
    num_devices$bait_tube <- num_transects[["total"]]*(standard_transect_length/device_spacing$bait_tube)
    # Set up
    method_set_up$bait_tube <- num_devices$bait_tube*(per_device_cost$bait_tube + per_bait_cost$bait_tube) + 
      misc_init_equip$bait_tube
    # Per day
    daily_cost$bait_tube <- num_teams[["bait_tube"]]*2*transect_hours_per_day$bait_tube*person_cost_per_hour + 
      per_bait_cost$bait_tube*(transects_checked_per_hour$bait_tube*transect_hours_per_day$bait_tube*(standard_transect_length/standard_between_transects))
    # Per week
    weekly_cost$bait_tube <- overhead_hours$bait_tube*person_cost_per_hour
    # Total days and weeks
    bait_tube_weeks <- 0
    bait_tube_days <- 0
    for(quarter in erad_quarters$bait_tube) {
      quarter_weeks <- round((max(erad_days[[quarter]]$bait_tube) - min(erad_days[[quarter]]$bait_tube))/7)
      quarter_days <- length(erad_days[[quarter]]$bait_tube)
      # Adding each quarters together
      bait_tube_weeks <- bait_tube_weeks + quarter_weeks
      bait_tube_days <- bait_tube_days + quarter_days
    }
    # Total cost
    total_cost$bait_tube <- method_set_up$bait_tube + bait_tube_days*daily_cost$bait_tube + 
      bait_tube_weeks*weekly_cost$bait_tube
  }
  
  # ADS costs
  if(erad_methods[1] %in% methods) {
    # Calculate the number of ADS treatments (each with 3 applications per treatment) & how many years
    num_treatments <- length(erad_quarters$ADS) # need to change this later (not currently separating years)
    num_ADS_years <- ceiling(max(erad_quarters$ADS/4)) # minimum of 1 year - may need to adjust this...
    misc_supplies_annual <- (10000/500)*area_size 
    # Calculate flight days in each treatment
    flight_hrs_per_app <- ((baits_per_ha*area_size)/max_bait_application_rate_hourly) # Per application
    flight_days_per_application <- max(round(area_size/max_ha_per_day, 1), 1) # minimum of 1 day
    flight_days_per_treatment <- flight_days_per_application*num_applications_per_treatment
    # Calculate travel costs for pilot and mechanic
    traveler_costs <- (flight_days_per_treatment + 
                         prep_days_per_treatment + 
                         contingency_days +
                         pilot_off_days + 
                         travel_days)*per_day_per_traveler*num_travelers
    # Calculating total bait cost for all treatments (ordered all at once)
    total_baits <- baits_per_ha*area_size*num_applications_per_treatment*num_treatments
    bait_bulk_category <- max(which(bait_cost_by_amount$lower_limit < total_baits))
    total_bait_cost <- total_baits*bait_cost_by_amount$bait_cost[bait_bulk_category]
    # Helicopter maintenance & crew costs per treatment
    heli_maintenance_per_treatment <- heli_maintenance_hourly_rate*flight_hrs_per_app*num_applications_per_treatment
    mechanic_service_per_treatment <- mechanic_service_daily_rate*(flight_days_per_treatment + 
                                                                      prep_days_per_treatment +
                                                                      contingency_days +
                                                                      pilot_off_days +
                                                                      travel_days)
    pilot_cost_per_treatment <- pilot_hourly_rate*8*(flight_days_per_treatment + 
                                                     prep_days_per_treatment + 
                                                     contingency_days +
                                                     pilot_off_days +
                                                     travel_days*2) # travel days is counted multiple times to be sure that the pilot is fully compensated
    ground_crew_cost_per_treatment <- ground_crew_hourly_rate*(flight_days_per_treatment + 
                                                                 ground_prep_days +
                                                                 contingency_days) + leadership_hourly*flight_days_per_treatment
    # Total cost
    total_cost$ADS <- total_bait_cost + (traveler_costs + heli_maintenance_per_treatment +
                                           mechanic_service_per_treatment + pilot_cost_per_treatment + 
                                           ground_crew_cost_per_treatment)*num_treatments + (freezer_rental_annual + 
                                                                                               misc_supplies_annual)*num_ADS_years
  }
  
  summed_cost <- sum(unlist(total_cost))
  
  return(list(total_cost = total_cost,
              summed_cost = summed_cost))
}


