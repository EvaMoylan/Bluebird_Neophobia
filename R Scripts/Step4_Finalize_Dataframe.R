################################################################################
# Step 3: Add Additional Variables of Interest
# Updated: 10/25/26
# Objective: Add percent canopy cover, missing hide distances and AudioMoth spl data to trial dataset
# Authors: Eva Moylan
################################################################################

#################################### Set-up #################################### 

# Load packages
library(tidyverse)
library(here)
library(geosphere)

# Load main neophobia assay dataset from Step 1 - probably should convert all these to .rds
neo_trials_censor <- read_csv(here("Data","Censor_BBNeophobia_261025.csv"))
# Maybe re-organize these steps to avoid intermediate files

######################### Add missing blind distances ##########################
# Calculate blind distance from the nest for trials that don't have a rangefinder distance

# Import nestbox coordinates
box_coords <- read_csv(here("Data", "Nestbox_Site_Coords.csv")) |>
  rename(nest_box = name,
         box_lat = lat,
         box_lon = lon) # Rename for easy joining

# Import hide coordinates that are missing for 2024 and one set of trials in 2025
hide_coords <- read_csv(here("Data", "Missing_Hide_Dist.csv"))

# Join hide and box coordinates and calculate distances
hide_dist <- left_join(hide_coords, box_coords, by = "nest_box") |>
  # Shortest distance between the two coordinates
  mutate(hide_distance_m = distHaversine(cbind(hide_lon, hide_lat), cbind(box_lon, box_lat))) |>
  select(nest_box, hide_distance_m, date)

# Fill in the missing distances in the main dataset
neo_trials_main <- left_join(neo_trials_censor, hide_dist, by = c("nest_box", "date")) |>
  mutate(hide_distance_m = coalesce(hide_distance_m.x, hide_distance_m.y),
         hide_distance_m = case_when(nest_id == "CP02B1_2024" ~ 23.85523, .default = hide_distance_m), # One trial for this nest didn't assign the distance so manually assign
         .after = chick_number) |>
  select(-c(hide_distance_m.x, hide_distance_m.y))
# TP04B6_2024 and CP06B5_2024 have no recorded hide coordinates, so distances are unknown

######################## Add percent canopy cover data #########################

# Load 40m radius percent canopy cover data from Step 2
final_lidar <- read_csv(here("Data", "lidar_metrics_40m_canopy.csv"))

# Join with the main dataset
neo_trials_main <- neo_trials_main |>
  left_join(final_lidar, by = "nest_box") |>
  select(-(point:rumple)) |>
  rename(percent_canopy = perc_can) |>
  relocate(percent_canopy, .after = hide_distance_m)

################### Add AudioMoth sound pressure level data ####################

# Load all AM data

# Save main cleaned dataframe
write_csv(neo_trials_main, here("Data","Main_BBNeophobia_261025.csv"))

