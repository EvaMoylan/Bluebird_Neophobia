################################################################################
# Step 5: Data Analysis
# Date edited: 11/24/26
# Objective: Data exploration and modeling
# Author: Eva Moylan
################################################################################

#################################### Set-up ####################################
# Load packages
library(tidyverse)
library(here)
library(glmmTMB)
library(performance)
library(see)
library(survival)
library(coxme)
library(modelr)
library(MuMIn)
library(survminer)
library(lme4)
library(ggeffects)
library(DHARMa)

# Load in cleaned and censored dataframe from previous steps
final_neo <- read_rds(here("Outputs","Main_BBNeophobia.rds"))

# Inspect structure
str(final_neo)

# Create WEBL only dataframe
final_webl <- final_neo |>
  filter(species == "WEBL")

final_webl |> count(feed_status) # 13 censored trials, 99 observed, 18 failed

# Create MOBL only dataframe
final_mobl <- final_neo |>
  filter(species == "MOBL")

final_mobl |> count(feed_status) # 8 censored trials, 30 observed, 5 failed

################################################################################
##################### Variable Exploration and Selection #######################
################################################################################

### Entire data distributions ###

# Latency to feed (LTF)
hist(final_neo$first_feed_latency_min) # Right-skew
hist(log(final_neo$first_feed_latency_min + 1)) # Normal

# Latency to approach (LTA)
hist(final_neo$first_approach_latency_min) # Right-skew
hist(log(final_neo$first_approach_latency_min + 1)) # Normal

# Approach to feed latency (AFL)
hist(final_neo$diff_approach_prov_min) # Right-skew, primarily 0s - coxme models can handle this better than glmm
hist(log(final_neo$diff_approach_prov_min + 1)) # Still right-skew but slightly less dramatic

### WEBL data distributions ###

# Latency to feed (LTF)
hist(final_webl$first_feed_latency_min) # Right-skew
hist(sqrt(final_webl$first_feed_latency_min)) # Normal
hist(log(final_webl$first_feed_latency_min + 1)) # Normal

unique(final_webl$first_feed_latency_min)

# Latency to approach (LTA)
hist(final_webl$first_approach_latency_min) # Right-skew
hist(log(final_webl$first_approach_latency_min + 1)) # Normal

# Approach to feed latency (AFL)
hist(final_webl$diff_approach_feed_min) # Right-skew, primarily 0
hist(log(final_webl$diff_approach_feed_min + 1)) # Right-skew but less dramatic

# Hide distance
hist(final_webl$hide_distance_m) # Normal

# Canopy cover
hist(final_webl$percent_canopy) # Normal

# Cluch size
hist(final_webl$chick_number)

# Nest age
hist(final_webl$nest_age) # Right skew

# Visualize rough samples across treatments
ggplot(final_webl, aes(x = LAeq_slm_dB, y = light_level_lux, color = treatment)) +
  geom_jitter(width = 0.12, height = 0.12, size = 2.5, alpha = 0.75) +
  labs(
    x = "Noise",
    y = "Light",
    color = "Treatment type estimate"
  ) +
  theme_classic()

### MOBL data distributions ###

# Latency to feed (LTF)
hist(final_mobl$first_feed_latency_min) # Even

# Latency to approach (LTA)
hist(final_mobl$first_approach_latency_min) # Right-skew
hist(log(final_mobl$first_approach_latency_min + 1)) # Normal-ish, 0 still skewing

# Approach to feed latency (AFL)
hist(final_mobl$diff_approach_feed_min) # Right-skew, primarily 0
hist(log(final_mobl$diff_approach_feed_min + 1)) # Slight right-skew but much more even

# Hide distance
hist(final_mobl$hide_distance_m) # Normal

# Canopy cover
hist(final_mobl$percent_canopy) # Normal

# Visualize rough samples across treatments
ggplot(final_mobl, aes(x = LAeq_slm_dB, y = light_level_lux, color = treatment)) +
  geom_jitter(width = 0.12, height = 0.12, size = 2.5, alpha = 0.75) +
  labs(
    x = "Noise",
    y = "Light",
    color = "Treatment type estimate"
  ) +
  theme_classic()

################################################################################
########################## Western Bluebird Modeling ###########################
################################################################################

# How many total WEBL trials of each are there in the whole dataset (2024-2026)?
final_webl |> 
  filter(!is.na(trial_start_time_hms)) |> # Remove trials that never started
  count(stimulus_object, # Sample good overall: 33 Control (2 cens.), 18 Green Ball (2 cens.), 23 Pink Beads (1 cens.)
        appr_censored = approach_status == 0, 
        feed_censored = feed_status == 0)

# Filter a non-censored data set (only 5 censored trials)
final_webl_nc <- final_webl |>
  mutate(trial_type = as.factor(trial_type)) |>
  filter(feed_status != 0)

####################### Trial categorization exploration #######################

# Primary hypotheses are structured around comparison between trial types (control v. stimulus)
# Confirm that this is still sound with all three years of data
# Nested random effects structure simplified for all models due to consistently low variance across sites and section groups. Year retained to include annual fluctuation in resources and habitat conditions.

# Fit a simple model with purely hypothesis variables using trial type (control vs. stimulus)
type_mod <- coxme(Surv(first_feed_latency_min, feed_status) ~ scale(LAeq_slm_dB)*scale(light_level_lux) + trial_type*scale(LAeq_slm_dB) + trial_type*scale(light_level_lux) + (1 | year/nest_id), data = final_webl)
summary(type_mod)

# Fit a simple model with purely hypothesis variables using object type (control vs. pink beads vs green ball)
obj_mod <- coxme(Surv(first_feed_latency_min, feed_status) ~ scale(LAeq_slm_dB)*scale(light_level_lux) + stimulus_object*scale(LAeq_slm_dB) + stimulus_object*scale(light_level_lux) + (1 | year/nest_id), data = final_webl)
summary(obj_mod) 
  
# Rank trial models
model.sel(type_mod, obj_mod, rank = "AICc") # Trial type strongly out performs object, support for decision to combine them 

################ Exploration of ecologically relevant covariates ###############
# With sample size of 69 events, three two-way interactions and nested random effects, can't safely include top ideal variables that likely influence provisioning rate. 
# Can likely including 1 without over-fitting, so which variable accounts for most variance in this dataset

# Sex of adult that provisioned (unknown for censored trials so fit using uncensored dataset)
sex_mod <- coxme(Surv(first_feed_latency_min, feed_status) ~ scale(LAeq_slm_dB)*scale(light_level_lux) + trial_type*scale(LAeq_slm_dB) + trial_type*scale(light_level_lux) + first_feed_sex + (1 | year/nest_id), data = final_webl_nc)
summary(sex_mod) # Sex has extremely weak effect on provisioning latency

# Clutch size
chx_mod <- coxme(Surv(first_feed_latency_min, feed_status) ~ scale(LAeq_slm_dB)*scale(light_level_lux) + trial_type*scale(LAeq_slm_dB) + trial_type*scale(light_level_lux) + scale(chick_number) + (1 | year/nest_id), data = final_webl)
summary(chx_mod) # Strong

# Nest age
age_mod <- coxme(Surv(first_feed_latency_min, feed_status) ~ scale(LAeq_slm_dB)*scale(light_level_lux) + trial_type*scale(LAeq_slm_dB) + trial_type*scale(light_level_lux) + scale(nest_age) + (1 | year/nest_id), data = final_webl)
summary(age_mod) # Very weak

# Our hide distance from the nest
hide_mod <- coxme(Surv(first_feed_latency_min, feed_status) ~ scale(LAeq_slm_dB)*scale(light_level_lux) + trial_type*scale(LAeq_slm_dB) + trial_type*scale(light_level_lux) + scale(hide_distance_m) + scale(chick_number) + (1 | year/nest_id), data = final_webl)
summary(hide_mod) # Moderately weak

####################### Primary Cox mixed effects models #######################

# Fit primary models: 
# 1 - Two-way interactions between slm noise, light, trial type (entire data)
# 2 - Three-way interaction between slm noise, light, trial type (entire data)

# FIT MODEL #1 - with slm noise, light, trial type two-way interactions model 
webl_LTP_1 <- coxme(Surv(first_feed_latency_min, feed_status) ~ scale(LAeq_slm_dB)*scale(light_level_lux) + trial_type*scale(LAeq_slm_dB) + trial_type*scale(light_level_lux) + scale(chick_number) + (1 | year/nest_id), data = final_webl)
summary(webl_LTP_1)

# Okay so for every 1 db increase in noise, birds feed ~twice as fast in stimulus trials re... 
# In noisier environments, there's a shorter latency to provision in the presence of a novel object relative to control trials

# Support for random effects structure from a previous model including all four nested re:
# Site and section explain near zero variance - remove
#Random effects
#Group                     Variable    Std Dev      Variance    
#year/section/site/nest_id (Intercept) 3.445044e-01 1.186833e-01
#year/section/site         (Intercept) 1.217354e-02 1.481951e-04
#year/section              (Intercept) 9.996748e-03 9.993497e-05
#year                      (Intercept) 3.208670e-01 1.029557e-01

# CHECK MODEL #1 - fit semi-equivalent coxph models to simulate residuals

# Run cox proportional hazard model with just nest_id re to visualize residuals/ proportional hazard assumptions
# datanovia.com/learn/biostatistics/survival-analysis/testing-ph-assumptions
webl_LTP_1ph <- coxph(Surv(first_feed_latency_min, feed_status) ~ scale(LAeq_slm_dB)*scale(light_level_lux) + trial_type*scale(LAeq_slm_dB) + trial_type*scale(light_level_lux) + scale(chick_number), data = final_webl)
summary(webl_LTP_1ph) # Results don't appear to change much

# Is variance different between coxme model and coxph model? (coxme package, Therneau, 2026)
anova(webl_LTP_1, webl_LTP_1ph) # Good, the two models aren't very different

# Simulate scaled Schoenfeld residuals to check hazard assumptions - variable effect should not change over time
ggcoxzph(cox.zph(webl_LTP_1ph)) # Proportional hazard assumptions are not violated
# This should work with coxme model?? Aborts session

# Check for influential trials
ggcoxdiagnostics(webl_LTP_1ph, 
                 type = "dfbeta", # Coefficient change with each observation removed
                 ox.scale = "observation.id", # Shows residuals by rows (i.e. trials)
                 hline = TRUE,
                 ggtheme = theme_bw()) # Looks good overall

# FIT MODEL #2 - Fit slm noise, light, trial type three-way interaction model 
webl_LTP_2 <- coxme(Surv(first_feed_latency_min, feed_status) ~ scale(LAeq_slm_dB)*scale(light_level_lux)*trial_type + scale(chick_number) + (1 | year/nest_id), data = final_webl)
webl_LTP_2

# CHECK MODEL #2 - fit semi-equivalent coxph models to simulate residuals

# Run cox proportional hazard model with just nest_id re to visualize residuals/ proportional hazard assumptions
# datanovia.com/learn/biostatistics/survival-analysis/testing-ph-assumptions
webl_LTP_2ph <- coxph(Surv(first_feed_latency_min, feed_status) ~ scale(LAeq_slm_dB)*scale(light_level_lux)*trial_type + scale(chick_number), data = final_webl)
summary(webl_LTP_2ph) # Results don't appear to change much

# Is variance different between coxme model and coxph model? (coxme package, Therneau, 2026)
anova(webl_LTP_2, webl_LTP_2ph) # Good, the two models aren't very different

# Simulate scaled Schoenfeld residuals to check hazard assumptions - variable effect should not change over time
ggcoxzph(cox.zph(webl_LTP_2ph)) # Proportional hazard assumptions are not violated

# Check for influential trials
ggcoxdiagnostics(webl_LTP_2ph, 
                 type = "dfbeta", # Coefficient change with each observation removed
                 ox.scale = "observation.id", # Shows residuals by rows (i.e. trials)
                 hline = TRUE,
                 ggtheme = theme_bw()) # Looks good overall

# RANK MODEL 1 & MODEL 2
model.sel(webl_LTP_1, webl_LTP_2, rank = "AICc") # Model 1 with two way interactions fits the data better, but moderate support for the three way interaction
#              logLik  AICc  delta weight
# webl_LTP_1 -223.081 488.0  0.00  0.788
# webl_LTP_2 -223.135 490.6  2.63  0.212

# For main model, does removing censored data seem to impact model output
webl_nc_1 <- coxme(Surv(first_feed_latency_min, feed_status) ~ scale(LAeq_slm_dB)*scale(light_level_lux) + trial_type*scale(LAeq_slm_dB) + trial_type*scale(light_level_lux) + scale(chick_number) + (1 | year/nest_id), data = final_webl_nc)
summary(webl_nc_1) # Yes actually it slightly does, however, this is removing a lot of light trials

# Look more closely at the censored trials
final_webl |>
  filter(trial_status != "failed") |>
  group_by(trial_type) |>
  summarise(light_count = sum(light_level_lux != 0, na.rm = TRUE)) 

final_webl_nc |>
  group_by(trial_type) |>
  summarise(light_count = sum(light_level_lux != 0, na.rm = TRUE)) # removes some light trials but doesn't change sample sizes much. Some high light levels are censored though, so coule be skewing relationship

####################### Visualize primary model results ########################

# Noise*trial plot #

# Removing censored data didn't impact magnitude of hazard ratios or effects for noise*trial_type
# So, fit a glmmTMB on uncensored data just for visualizing linear predictions for this interaction
glmm_webl_LTP <- glmmTMB(first_feed_latency_min ~ scale(LAeq_slm_dB)*scale(light_level_lux) + trial_type*scale(LAeq_slm_dB) + trial_type*scale(light_level_lux) + scale(chick_number) + (1 | year/nest_id), REML = FALSE, data = final_webl_nc)
summary(glmm_webl_LTP)
plot(simulateResiduals(glmm_webl_LTP))

# Predict response
noisetrial_predict <- predict_response(glmm_webl_LTP, terms = c("LAeq_slm_dB","trial_type"))
plot(noisetrial_predict, show_data = TRUE)

# Clutch size plot # 

# Look at raw clutch size data 
ggplot(final_webl_nc, aes(x = chick_number, y = first_feed_latency_min, color = trial_type)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE) # Raw control and stimulus slopes are nearly identical

# Predict response
clutch_predict <- predict_response(glmm_webl_LTP, terms =  c("chick_number","trial_type"))
plot(clutch_predict, show_data = TRUE)

### Forest plot using coxme ###

# Survival plotting functions don't support coxme so extract values from coxme output
webl1_output <- data.frame(
  variable = names(fixef(webl_LTP_1)),
  estimate = fixef(webl_LTP_1),
  se = sqrt(diag(as.matrix(vcov(webl_LTP_1))))) |>
  mutate(HR = exp(estimate),
         conf.low95 = exp(estimate - 1.96 * se),
         conf.high95 = exp(estimate + 1.96 * se))

# Reformat names
webl1_output <- webl1_output |>
  mutate(variable = as.character(variable),
         variable = case_when(variable == "scale(LAeq_slm_dB)" ~ "Noise Level (LAeq)",
                       variable == "scale(light_level_lux)" ~ "Light Level (lux)",
                       variable == "trial_typeStimulus" ~ "Novel Object Trial",
                       variable == "scale(chick_number)" ~ "Clutch Size",
                       variable == "scale(LAeq_slm_dB):scale(light_level_lux)" ~ "Noise Level (LAeq) : Light Level (lux)",
                       variable == "scale(LAeq_slm_dB):trial_typeStimulus" ~ "Noise Level (LAeq) : Novel Object Trial",
                       variable == "scale(light_level_lux):trial_typeStimulus" ~ "Light Level (lux) : Novel Object Trial",
                       .default = variable))

# Plot hazard ratios
ggplot(webl1_output, aes(x = HR, y = reorder(variable, HR), color = reorder(variable, HR))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkgray") +
  geom_errorbar(aes(xmin = conf.low95, xmax = conf.high95), width = 0.3) +
  geom_point(size = 3) +
  scale_color_viridis_d() +
  scale_x_log10(breaks = seq(0.5,4.5, by = 0.5)) + # HR are not additive, so plot on log10 scale
  labs(x = "Hazard Ratio (Log-scale)",
       y = "Predictor") +
  guides(color = "none") +
  theme_classic()
