################################################################################
# Step 1: Neophobia Assay Data Preparation
# Date edited: 10/10/26
# Objective: Calculate and censor latencies
# Author: Eva Moylan
################################################################################

#################################### Set-up ####################################

# Load packages
library(tidyverse)
library(here)

# Load raw trial data without latencies or missing hide distances calculated
neo_trials <- read_csv(here("Data","EM_Bluebird_Neophobia_261020.csv")) |>
  slice(1:130) # Re-save this sliced version for final repository

################################ Minor cleaning ################################

neo_trials <- neo_trials |>
  mutate(year = year(mdy(date)), # Add a year column
         stimulus_object = case_when(stimulus_object == "Green Ball" ~ "Green ball", # Fix case differences in object names
                                     stimulus_object == "Pink Beads" ~ "Pink beads",
                                     .default = stimulus_object),
         nest_box = str_remove(nest_box, "_2026")) |> # Remove formatting mistake
  relocate(year, .after = nest_box)

########################### Calculate latency times ############################

# Check data structure
str(neo_trials) # Times are in hms

# Calculate latencies for first approach, first feed, and first nest entry relative to trial start time
# Bird specific latency columns already existed in the dataset, so re-calculate the latencies and rename them later
neo_trials_latencies <- neo_trials |>
  mutate(bird1_approach_latency_s = as.numeric(difftime(bird1_first_approach_hms, trial_start_time_hms, units = "min")), # Rename column as minutes later
         bird1_feed_latency_s = as.numeric(difftime(bird1_return_feed_hms, trial_start_time_hms, units = "min")),
         bird1_entry_latency_min = as.numeric(difftime(bird1_return_entry_hms, trial_start_time_hms, units = "min")),
         bird2_approach_latency_s = as.numeric(difftime(bird2_first_approach_hms, trial_start_time_hms, units = "min")), # Rename column as minutes later
         bird2_feed_latency_s = as.numeric(difftime(bird2_return_feed_hms, trial_start_time_hms, units = "min")),
         bird2_entry_latency_min = as.numeric(difftime(bird2_return_entry_hms, trial_start_time_hms, units = "min")),
         first_feed_time = pmin(bird1_return_feed_hms, bird2_return_feed_hms, na.rm = TRUE), # Pick earliest time or only available time in cases where one bird fed
         first_approach_time = pmin(bird1_first_approach_hms, bird2_first_approach_hms, na.rm = TRUE),
         first_entry_time = pmin(bird1_return_entry_hms, bird2_return_entry_hms, na.rm = TRUE),
         first_feed_latency_min = as.numeric(difftime(first_feed_time, trial_start_time_hms, units = "min")), # Earliest time minus trial start time
         first_approach_latency_min = as.numeric(difftime(first_approach_time, trial_start_time_hms, units = "min")),
         first_entry_latency_min = as.numeric(difftime(first_entry_time, trial_start_time_hms, units = "min"))) |> 
  rename(bird1_approach_latency_min = bird1_approach_latency_s,
         bird1_feed_latency_min = bird1_feed_latency_s,
         bird2_approach_latency_min = bird2_approach_latency_s,
         bird2_feed_latency_min = bird2_feed_latency_s,) |>
  relocate(c(bird2_entry_latency_min, first_approach_time, first_feed_time, first_entry_time, first_approach_latency_min, first_feed_latency_min, first_entry_latency_min), 
           .after = bird2_feed_latency_min) |>
  relocate(bird1_entry_latency_min, 
           .after = bird1_feed_latency_min) |>
  select(-c(first_approach_latency_s, first_feed_latency_s, first_entry_latency_s)) |> # Remove some blank columns leftover from initial dataset
  # Add first feed sex
  mutate(first_feed_sex = case_when(is.na(bird2_feed_latency_min) & !is.na(bird1_feed_latency_min) ~ bird1_sex, # Add first feed sex
                                    is.na(bird1_feed_latency_min) & !is.na(bird2_feed_latency_min) ~ bird2_sex,
                                    bird1_feed_latency_min < bird2_feed_latency_min ~ bird1_sex,
                                    bird1_feed_latency_min > bird2_feed_latency_min ~ bird2_sex),
         first_feed_ID = case_when(is.na(bird2_feed_latency_min) & !is.na(bird1_feed_latency_min) ~ bird1_bands, # Add first feed sex
                                   is.na(bird1_feed_latency_min) & !is.na(bird2_feed_latency_min) ~ bird2_bands,
                                   bird1_feed_latency_min < bird2_feed_latency_min ~ bird1_bands,
                                   bird1_feed_latency_min > bird2_feed_latency_min ~ bird2_bands),
         .after = first_feed_latency_min) |>
  # Add first approach sex
  mutate(first_approach_sex = case_when(is.na(bird2_approach_latency_min) & !is.na(bird1_approach_latency_min) ~ bird1_sex,
                                    is.na(bird1_approach_latency_min) & !is.na(bird2_approach_latency_min) ~ bird2_sex,
                                    bird1_approach_latency_min < bird2_approach_latency_min ~ bird1_sex,
                                    bird1_approach_latency_min > bird2_approach_latency_min ~ bird2_sex),
         first_approach_ID = case_when(is.na(bird2_approach_latency_min) & !is.na(bird1_approach_latency_min) ~ bird1_bands,
                                        is.na(bird1_approach_latency_min) & !is.na(bird2_approach_latency_min) ~ bird2_bands,
                                        bird1_approach_latency_min < bird2_approach_latency_min ~ bird1_bands,
                                        bird1_approach_latency_min > bird2_approach_latency_min ~ bird2_bands),
         .after = first_approach_latency_min)

# TP02B5 typo, bird 2 approach before trial start - need to go and figure out what happened there

################################# Censor data ##################################

# Prevent automatic scientific notation during following pipeline
options(scipen = 999)

# Right-censor trials that began but birds never performed certain behaviors before trial ended
neo_trials_censor <- neo_trials_latencies |>
    mutate(first_approach_latency_min = if_else(!is.na(trial_start_time_hms) & is.na(first_approach_latency_min), # If the trial started but there is no first approach...
                                          20, # Assign full trial time of 20 minutes
                                          first_approach_latency_min), # If there is, keep original latency
           approach_censor = if_else(first_approach_latency_min == 20, # If the approach was right-censored...
                                 0, # Assign 0
                                 1), # If uncensored, assign 0
           first_feed_latency_min = if_else(!is.na(trial_start_time_hms) & is.na(first_feed_latency_min),
                                          20,
                                          first_feed_latency_min),
           feed_censor = if_else(first_feed_latency_min == 20,
                                 0,
                                 1),
           first_entry_latency_min = if_else(!is.na(trial_start_time_hms) & is.na(first_entry_latency_min),
                                          20,
                                          first_feed_latency_min),
           entry_censor = if_else(first_entry_latency_min == 20,
                                 0,
                                 1),
           diff_approach_feed_min = as.numeric(first_feed_latency_min - first_approach_latency_min),
           diff_censor = if_else(feed_censor == 0 | approach_censor == 0,
                                  0,
                                  1)) |>
  relocate(c(diff_approach_feed_min, approach_censor, feed_censor, entry_censor, diff_censor), .after = first_entry_latency_min)

# Explore censored trials
neo_trials_censor |> count(feed_censor) # 13 censored trials of 112 total completed
neo_trials_censor |> count(approach_censor) # 4 censored trials of 112 total completed  
neo_trials_censor |> count(entry_censor) # 13 censored trials of 112 total completed  

# Save censored data
write_csv(neo_trials_censor, here("Data","Censor_BBNeophobia_261025.csv"))
