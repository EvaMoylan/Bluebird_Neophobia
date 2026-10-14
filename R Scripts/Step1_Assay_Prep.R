################################################################################
# Step 1: Neophobia Assay Data Preparation
# Date edited: 10/10/26
# Objective: Add missing variables and calculate latencies
################################################################################

# Load packages
library(tidyverse)
library(here)
library(geosphere)

# Load raw neophobia assay data
neo_trials <- read_csv(here("Data","EM_Bluebird_Neophobia_261010.csv")) |>
  slice(1:130)
# ^ this is preliminary data for now, still doing minor excel prep

################################################################################
############################### DATA PREPARATION ###############################
################################################################################

######################### Add missing blind distances ##########################
# Calculate blind distance from the nest for trials that don't have a rangefinder distance

# Import nestbox coordinates
box_coords <- read_csv(here("Data", "Nestbox_Site_Coords.csv")) |>
  rename(nest_box = name,
         box_lat = lat,
         box_lon = lon) # Rename for easy joining

# Import hide coordinates that are missing for 2024 and one set of trials in 2025
hide_coords <- read_csv(here("Data", "Missing_Hide_Dist.csv")) |>
  filter(!is.na(hide_lat)) # Two nests didn't have any coordinates recorded so remove

# Join hide and box coordinates and calculate distances
hide_dist <- left_join(hide_coords, box_coords, by = "nest_box") |>
  # Shortest distance between the two coordinates
  mutate(hide_distance_m = distHaversine(cbind(hide_lon, hide_lat), cbind(box_lon, box_lat))) |>
  select(nest_box, hide_distance_m, date)

# Fill in the missing distances in the main dataset
neo_trials_hide <- left_join(neo_trials, hide_dist, by = c("nest_box", "date")) |>
  mutate(hide_distance_m = coalesce(hide_distance_m.x, hide_distance_m.y),
         hide_distance_m = case_when(nest_id == "CP02B1_2024" ~ 23.85523, .default = hide_distance_m),
         .after = chick_number) |>
  select(-c(hide_distance_m.x, hide_distance_m.y))

########################### Calculate latency times ############################

# Check data structure
str(neo_trials_clean) # Times need to be converted to hms

unique(neo_trials_hide$bird1_first_approach_hms)
unique(neo_trials_clean$bird1_first_approach_hms)

# Prep data
neo_trials_latencies <- neo_trials_hide |>
  mutate(bird1_first_approach_hms = na_if(bird1_first_approach_hms,"???"),
         bird1_first_approach_hms = as_hms(parse_date_time(bird1_first_approach_hms, orders = c("HM", "HMS"))), # do this for now, won't need it later
         first_feed_time = pmin(bird1_return_feed_hms, bird2_return_feed_hms, na.rm = TRUE), # Pick earliest time or single time if other is missing
         first_approach_time = pmin(bird1_first_approach_hms, bird2_first_approach_hms, na.rm = TRUE),
         first_feed_latency_s = as.numeric(difftime(first_feed_time, trial_start_time_hms, units = "secs")),
         first_approach_latency_s = as.numeric(difftime(first_approach_time, trial_start_time_hms, units = "secs"))) |> # earliest time minus trial start time) # earliest time minus trial start time
  relocate(c(first_feed_time, first_approach_time), .after = bird2_feed_latency_s)

#first_approach_sex = if_else(bird1_first_approach_hms < bird2_first_approach_hms,
                             #bird1_sex,
                             #bird2_sex),
#first_feed_sex = 

################################# Censor data ##################################
  
neo_trials_censor <- neo_trials_latencies |>
    mutate(first_feed_latency_s = if_else(!is.na(trial_start_time_hms) & is.na(first_feed_latency_s),
                                             1200,
                                             first_feed_latency_s),
           feed_censor = if_else(first_feed_latency_s == 1200,
                                 1,
                                 0),
           first_approach_latency_s = if_else(!is.na(trial_start_time_hms) & is.na(first_approach_latency_s),
                                          1200,
                                          first_approach_latency_s),
           approach_censor = if_else(first_approach_latency_s == 1200,
                                 1,
                                 0),
           diff_approach_feed_s = first_feed_latency_s-first_approach_latency_s,
           diff_censor = if_else(feed_censor|approach_censor == 1,
                                  1,
                                  0)) |>
  relocate(c(diff_approach_feed_s, feed_censor, approach_censor, diff_censor), .after = first_approach_latency_s)

neo_trials_censor |> count(feed_censor) # 13 censored trials of 112 total completed
neo_trials_censor |> count(approach_censor) # 4 censored trials of 112 total completed  

write_csv(neo_trials_censor, here("Data","Main_BBNeophobia_261013.csv"))

## To Do: ##
# How to control for trials where we were there but they never began?
# Number of trials attempted?
# Trial order?

