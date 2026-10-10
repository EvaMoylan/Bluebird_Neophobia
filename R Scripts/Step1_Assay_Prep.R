################################################################################
# Step 1: Neophobia Assay Data Preparation
# Date edited: 9/5/26
# Objective: Add missing variables and calculate latencies
################################################################################

# Load packages
library(tidyverse)
library(here)
library(glmmTMB)
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

write_csv(neo_trials_hide, here("Data","Main_BBNeophobia_261010.csv"))

########################### Calculate latency times ############################

# Check data structure
str(neo_trials) # Times already in hms

# Prep data
neo_trials_clean <- neo_trials |>
  add_count(nest_id, name = "nest_total_attempts") |>
  relocate(nest_total_attempts, .after = nest_id) |>
  mutate(first_feed_latency_s = if_else(!is.na(trial_start_time_hms) & is.na(first_feed_latency_s),
                                        1200,
                                        first_feed_latency_s),
         diff_approach_feed_s = first_feed_latency_s - first_approach_latency_s)


    #First.Feed.HMS = pmin(Female.Return.Nest.Entry.Time..hr.min.sec., Male.Return.Nest.Entry.Time..hr.min.sec., na.rm = TRUE), # Whichever feed is sooner
        # First.Feed.Latency = First.Feed.HMS - Trial.Start.Time...hr.min.sec., # First feed minus trial start
         #First.Feed.Latency = as.numeric(First.Feed.Latency),

## To Do: ##
# Calculate the latencies for feed / entry / approach / diff between
# Censor trials with no feed
# How to control for trials where we were there but they never began?
# Number of trials attempted? - done

################################# Censor Data ##################################

