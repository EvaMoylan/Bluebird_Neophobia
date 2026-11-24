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
neo_trials <- read_csv(here("Data","Raw_Bluebird_Neophobia.csv"))

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
         first_feed_bands = case_when(is.na(bird2_feed_latency_min) & !is.na(bird1_feed_latency_min) ~ bird1_bands, # Add first feed sex
                                   is.na(bird1_feed_latency_min) & !is.na(bird2_feed_latency_min) ~ bird2_bands,
                                   bird1_feed_latency_min < bird2_feed_latency_min ~ bird1_bands,
                                   bird1_feed_latency_min > bird2_feed_latency_min ~ bird2_bands),
         .after = first_feed_latency_min) |>
  # Add first approach sex
  mutate(first_approach_sex = case_when(is.na(bird2_approach_latency_min) & !is.na(bird1_approach_latency_min) ~ bird1_sex,
                                    is.na(bird1_approach_latency_min) & !is.na(bird2_approach_latency_min) ~ bird2_sex,
                                    bird1_approach_latency_min < bird2_approach_latency_min ~ bird1_sex,
                                    bird1_approach_latency_min > bird2_approach_latency_min ~ bird2_sex),
         first_approach_bands = case_when(is.na(bird2_approach_latency_min) & !is.na(bird1_approach_latency_min) ~ bird1_bands,
                                        is.na(bird1_approach_latency_min) & !is.na(bird2_approach_latency_min) ~ bird2_bands,
                                        bird1_approach_latency_min < bird2_approach_latency_min ~ bird1_bands,
                                        bird1_approach_latency_min > bird2_approach_latency_min ~ bird2_bands),
         .after = first_approach_latency_min)

################################# Censor data ##################################

# Prevent automatic scientific notation during following pipeline
options(scipen = 999)

# First, test correlation between first provision latencies and first nest entry latencies to see if they're worth looking at separately
cor.test(neo_trials_latencies$first_feed_latency_min, neo_trials_latencies$first_entry_latency_min, na.action = na.omit)
# 99.8% correlated - use feed latency going forward

# Right-censor trials that began but birds never performed certain behaviors before trial ended
neo_trials_censor <- neo_trials_latencies |>
  select(-c(bird1_entry_latency_min, bird2_entry_latency_min, first_entry_latency_min)) |> # Don't need entry columns now 
  # Censor latency between first approach and first feed
  # We cannot take the difference after trials are censored, but trials where a bird approached but never provisioned can be censored. But, latency for trials where a bird never approached and therefore never provisioned need to be NA for this variable
    mutate(diff_approach_prov_min = if_else(!is.na(first_approach_latency_min) & !is.na(first_feed_latency_min), # If there is a first approach and first provision time...
                                            as.numeric(first_feed_latency_min - first_approach_latency_min), # Take the difference between the first approach and provision
                                            NA), # Otherwise, leave NA
           diff_status = if_else(!is.na(diff_approach_prov_min), # If there is an observed latency in this column, assign 1
                                 1,
                                 NA), # All others will remain NA for now
           diff_approach_prov_min = if_else(!is.na(first_approach_latency_min) & is.na(first_feed_latency_min), # If a bird approached but there was no provision...
                                            20 - first_approach_latency_min, # Assign difference until the trial ended
                                            diff_approach_prov_min), # Otherwise, keep the same
           diff_status = if_else(!is.na(diff_approach_prov_min) & is.na(diff_status), # Censored latencies are still NA in the censor column, so if there's now a value in the difference column, but censor is NA...
                                 0, # Then assign a 0 status
                                 diff_status), # Otherwise, keep status
           
           # Now censor approach latencies
           first_approach_latency_min = if_else(!is.na(trial_start_time_hms) & is.na(first_approach_latency_min), # If the trial started but there is no first approach...
                                          20, # Assign full trial time of 20 minutes
                                          first_approach_latency_min), # If there is, keep original latency
           approach_status = if_else(first_approach_latency_min == 20, # If the approach was right-censored...
                                 0, # Assign 0
                                 1), # If uncensored, assign 0
           
           # Now censor provision latencies
           first_feed_latency_min = if_else(!is.na(trial_start_time_hms) & is.na(first_feed_latency_min),
                                          20,
                                          first_feed_latency_min),
           feed_status = if_else(first_feed_latency_min == 20,
                                 0,
                                 1)) |>
  relocate(c(diff_approach_prov_min, approach_status, feed_status, diff_status), .after = first_feed_latency_min)

# Explore censored trials (entire data)
neo_trials_censor |> count(feed_status) # 13 censored trials, 99 observed, 18 NA
neo_trials_censor |> count(approach_status) # 4 censored trials, 108 observed, 18 NA
neo_trials_censor |> count(diff_status) # 9 censored trials of 99 observed, 22 NA

# Save censored data
write_rds(neo_trials_censor, here("Outputs","Censor_BBNeophobia.rds"))
