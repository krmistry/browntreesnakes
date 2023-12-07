#### Cost calculations per method

## Cost values are place holders, waiting for real costs from USDA

#### Costs shared by all transect-based methods
# - one time transect cutting for 1 transect X number of transects
# - transect maintenance at regular intervals (say 4 weeks), X number of transects
# - hourly rate for employees
init_transect <- 1000
maint_transect <- 200 # every 4 weeks
person_cost_per_hour <- 100

##### Visual survey costs 
# - number of searchers per night X search nights
# - one time purchase of equipment (air rifles, headlamps, etc)
# - person hours for any overhead activities, outside of searching, in a timescale that's easy to multiply - daily or weekly
headlamp <- 500
num_teams <- 1 # number of teams of 2 people searching
num_visual_searchers <- 4*2*num_teams*person_cost_per_hour
init_equip <- headlamp*2*num_teams + 500 # in addition to headlamps, other costs can be lumped together unless they're big purchases that might change over time
overhead_hours <- 10*person_cost_per_hour # weekly hours spent by any personnel on overhead activities related to visual surveys (data entry, processing, etc)

#### Trap costs
# - one time purchase of traps, perhaps with a percentage of replacement based on time
# - bait population maintenance, a flat rate per time unit (weekly would be good) - may not be entirely necessary
# - person hours to check traps & feed bait, performed every 2-3 days while traps are out


### Bait costs
# - one time purchase of bait tubes, with a percentage of replacement based on time
# - baits, to be replaced every 2-3 days 
# - person hours to replace baits, every 2-3 days


# Calculating number of transects in area (assuming 220 m long, 20 m between transects)
transect_length <- 220
between_transects <- 20
area_length <- area_size*erad_coverage$visual*10000/transect_length # converting the area covered by transect from hectares into meters squared
num_transects <- round(area_length/20)

erad_days$visual
erad_quarters$visual


