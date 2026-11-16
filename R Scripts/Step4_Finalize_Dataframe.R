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

# CP12B13_2026 4/26 - 0
# PC6B8_2026   5/2/26 - 0

# Load main neophobia assay dataset from Step 1
neo_trials_censor <- read_rds(here("Outputs","Censor_BBNeophobia.rds"))

######################### Add missing blind distances ##########################
# Calculate blind distance from the nest for trials that don't have a rangefinder distance

# Import nestbox coordinates
box_coords <- read_csv(here("Data", "Nestbox_Site_Coords.csv")) |>
  rename(nest_box = name,
         box_lat = lat,
         box_lon = lon) # Rename for easy joining

# Import hide coordinates that are missing for 2024 and one set of trials in 2025
hide_coords <- read_csv(here("Data", "Missing_Hide_Dist.csv")) |>
  mutate(date = mdy(date))

# Join hide and box coordinates and calculate distances
hide_dist <- left_join(hide_coords, box_coords, by = "nest_box") |>
  # Shortest distance between the two coordinates
  mutate(hide_distance_m = distHaversine(cbind(hide_lon, hide_lat), cbind(box_lon, box_lat))) |>
  select(nest_box, hide_distance_m, date)

# Fill in the missing distances in the main dataset
neo_trials_hide <- left_join(neo_trials_censor, hide_dist, by = c("nest_box", "date")) |>
  mutate(hide_distance_m = coalesce(hide_distance_m.x, hide_distance_m.y),
         hide_distance_m = case_when(nest_id == "CP02B1_2024" ~ 23.85523, .default = hide_distance_m), # One trial for this nest didn't assign the distance so manually assign
         .after = chick_number) |>
  select(-c(hide_distance_m.x, hide_distance_m.y))
# TP04B6_2024 and CP06B5_2024 have no recorded hide coordinates, so distances are unknown

######################## Add percent canopy cover data #########################

# Load 40m radius percent canopy cover data from Step 2
final_lidar <- read_csv(here("Outputs", "lidar_metrics_40m_canopy.csv"))

# Join with the main dataset
neo_trials_lidar <- neo_trials_hide |>
  left_join(final_lidar, by = "nest_box") |>
  select(-(point:rumple)) |>
  rename(percent_canopy = perc_can) |>
  relocate(percent_canopy, .after = hide_distance_m)

################### Add AudioMoth sound pressure level data ####################

# Load all AM data - or load .rds below instead
AM_24 <- read_csv(here("Data", "2024_AM_AverageSPL.csv"))
AM_25 <- read_csv(here("Data", "2025_AM_AverageSPL.csv"))
AM_26 <- read_csv(here("Data", "2026_AM_AverageSPL.csv"))

# Bind into one dataframe
AM_spl <- bind_rows(AM_24, AM_25, AM_26)

# Save as an RDS or read in previously saved RDS
# write_rds(AM_spl, here("Outputs", "All_Audiomoth_SPL.rds"))
AM_spl <- readRDS(here("Outputs", "All_Audiomoth_SPL.rds"))

# Create a key with just nest_id and trial date and number, from the neophobia dataframe
neo_key <- neo_trials_censor |>
  select(c(nest_id, date, trial_number))

# Filter and inspect AM data
neo_trials_AM <- AM_spl |>
  mutate(nest_id = str_c(nest_id, year(deployment_start_date), sep = "_")) |> # nest_id is nestbox in this dataframe, add year to the end to match neo data
  semi_join(neo_key, by = "nest_id") |> # Keep only nest_ids present in neophobia dataframe
  # 2026 AM data include two deployments, so select the deployment closest to the date of the trial
  left_join(neo_key, by = "nest_id", relationship = "many-to-many") |>
  filter(trial_number == 1) |> # keep only first trials to have one date per nest
  mutate(time_diff = abs(as.numeric(deployment_start_date - date)), # New column with difference in days between trial date and deployment start date
         .after = deployment) |>
  filter(time_diff == min(time_diff), .by = nest_id) |> # Keep shortest difference
  ungroup() |>
  select(-c(time_diff, date, trial_number)) # remove left over columns
  
# Are any nests missing deployments?
anti_join(neo_key, neo_trials_AM, by = "nest_id") |>
  distinct(nest_id) # Yes... 5 nests in 2026 unfortunately didn't get AM noise

# Did adjacent nests get noise?
# TP06B9 is a random box, so no
# CP01 site no
# TP04B11 - paired with B5 which has data!
# TP02B5 - paired with B11 which has data!

# Dataframe of replacement boxes to join with all audiomoth data
AM_replacement <- AM_spl |>
  mutate(nest_id = str_c(nest_id, year(deployment_start_date), sep = "_")) |>
  filter(nest_id %in% c("TP04B5_2026", "TP02B11_2026"),
         deployment == 1) # Only TP02B11 has 2 deployments, 1st is closer to trial date

# Join this into the AM trial data and keep relevant columns
neo_amrepl <- neo_trials_AM |>
  full_join(AM_replacement) |>
  select(c(nest_id, total_rec_time_hrs, L50_A_weighted)) |> # No comments about recordings on these files so don't keep comment column
  rename(am_rec_hrs = total_rec_time_hrs,
         LA50_am_dB = L50_A_weighted) |>
  mutate(nest_id = str_replace(nest_id, "TP04B5_2026", "TP04B11_2026"),
         nest_id = str_replace(nest_id, "TP02B11_2026", "TP02B5_2026"))

# Join audiomoth SPL with neophobia dataframe
neo_trials_main <- neo_trials_lidar |> 
  left_join(neo_amrepl, by = "nest_id") |>
  # Add comments about replacement noise from adjacent boxes
  mutate(comments = if_else(nest_id == "TP04B11_2026",
                       "Audiomoth noise from paired TP04B5 box.",
                       comments),
    comments = if_else(nest_id == "TP02B5_2026",
                       str_c(comments, "Audiomoth noise from paired TP02B11 box.", sep = " "),
                       comments)) |>
  relocate(c(LA50_am_dB, am_rec_hrs), .after = slm_noise_at_trial)
  

### Save main cleaned neophobia dataframe for analysis ###
write_rds(neo_trials_main, here("Outputs","Main_BBNeophobia.rds"))

