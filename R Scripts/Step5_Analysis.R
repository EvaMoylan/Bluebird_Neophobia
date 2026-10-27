################################################################################
# Step 5: Data Analysis
# Date edited: 10/25/26
# Objective: Data exploration and modeling
# Author: Eva Moylan
################################################################################

#################################### Set-up ####################################
# Load packages
library(tidyverse)
library(here)
library(glmmTMB)
library(lme4)
library(performance)
library(see)
library(rptR)
library(survival)
library(coxme)
library(modelr)

# Load in cleaned and censored dataframe
final_neo <- read_csv(here("Data","Main_BBNeophobia_261025.csv")) |>
  mutate(first_approach_latency_min = ifelse(row_number() == 113, 0, first_approach_latency_min)) # Replace this trial with 0 for now, need to listen to recording 

# Create WEBL only dataframe
final_webl <- final_neo |>
  filter(species == "WEBL")

# Create MOBL only dataframe
final_mobl <- final_neo |>
  filter(species == "MOBL")

################################################################################
##################### Variable Exploration and Selection #######################
################################################################################

### Entire data distributions ###

# Latency to feed (LTF)
hist(final_neo$first_feed_latency_min) # Right-skew
# Can't take log(0), do more research on if adding a constant is appropriate
hist(log(final_neo$first_feed_latency_min + 1)) # Normal

# Latency to approach (LTA)
hist(final_neo$first_approach_latency_min) # Right-skew
hist(log(final_neo$first_approach_latency_min + 1)) # Normal

# Approach to feed latency (AFL)
hist(final_neo$diff_approach_feed_min) # Right-skew, primarily 0s - coxme models can handle this better than glmm
hist(log(final_neo$diff_approach_feed_min + 1)) # Still right-skew but slightly less dramatic

# Test correlation between first nest entry and first feed latencies
cor.test(final_neo$first_feed_latency_min, final_neo$first_entry_latency_min, na.action = na.omit)
# 100% correlated - use feed latency going forward

### WEBL data distributions ###

# Latency to feed (LTF)
hist(final_webl$first_feed_latency_min) # Right-skew
hist(log(final_webl$first_feed_latency_min + 1)) # Normal

# Latency to approach (LTA)
hist(final_webl$first_approach_latency_min) # Right-skew
hist(log(final_webl$first_approach_latency_min + 1)) # Normal

# Approach to feed latency (AFL)
hist(final_webl$diff_approach_feed_min) # Right-skew, primarily 0
hist(log(final_webl$diff_approach_feed_min + 1)) # Right-skew but less dramatic

### MOBL data distributions ###

# Latency to feed (LTF)
hist(final_mobl$first_feed_latency_min) # Even

# Latency to approach (LTA)
hist(final_mobl$first_approach_latency_min) # Right-skew
hist(log(final_mobl$first_approach_latency_min + 1)) # Normal-ish, 0 still skewing

# Approach to feed latency (AFL)
hist(final_mobl$diff_approach_feed_min) # Right-skew, primarily 0
hist(log(final_mobl$diff_approach_feed_min + 1)) # Slight right-skew but much more even

########################## Variance of random effects ##########################

# Look at entire dataset first, if there's justification to reduce nested effects look at individual species

# [LTF] - What is the variance of the ideal random effects structure?
random_LTF <- glmmTMB(log(first_feed_latency_min + 1) ~ 1 + (1 | nest_id/site/section/year), data = final_neo, REML = T)
summary(random_LTF)

# [LTA] - What is the variance of the ideal random effects structure?
random_LTA <- glmmTMB(log(first_approach_latency_min + 1) ~ 1 + (1 | nest_id/site/section/year), data = final_neo, REML = T)
summary(random_LTA)

# [AFL] - What is the variance of the ideal random effects structure?
random_AFL <- glmmTMB(log(diff_approach_feed_min + 1) ~ 1 + (1 | nest_id/site/section/year), data = final_neo, REML = T)
summary(random_AFL)

# Variances are different for each response variable, so use the most complex structure for all models going forward

############################# Trial categorization  ############################

# Should stimulus object (lvls: control, pink beads, green ball) or trial type (lvls: stimulus, control) be used?
# * Re-run these with adjusted latencies
# * Come back to this and check performances

# How many total WEBL trials of each are there in the whole dataset (2024-2026)?
final_webl |>
  filter(!is.na(trial_start_time_hms)) |> # Remove trials that never started
  count(stimulus_object, # Sample good overall: 33 Control (2 cens.), 18 Green Ball (2 cens.), 23 Pink Beads (1 cens.)
        appr_censored = approach_censor == 0, 
        feed_censored = feed_censor == 0)

# [LTF] - Linear correlation between green ball and pink beads
obj_LTFwebl <- lmer(log(first_feed_latency_min + 1) ~ stimulus_object + (1 | section/year), data = final_webl, REML = F) # REML F for AIC ranking later
summary(obj_LTFwebl) # Corr. between objects = 0.365, fails to converge with most complex random effects, removed site and nest ID
# Full random effects structure converges without transformed latency, glmmTMB converges with transformed but lmer prints correlation

# [LTA] - Linear correlation between green ball and pink beads
obj_LTAwebl <- lmer(log(first_approach_latency_min + 1) ~ stimulus_object + (1 | section/year), data = final_webl)
summary(obj_LTAwebl) # Corr. between objects = 0.35, same convergence note

# [AFL] - Linear correlation between green ball and pink beads
obj_AFLwebl <- lmer(log(diff_approach_feed_min + 1) ~ stimulus_object + (1 | section/year), data = final_webl)
summary(obj_AFLwebl) # Corr. between objects = 0.381, same convergence note

# How many total MOBL trials of each are there in the whole dataset (2024-2026)?
final_mobl |>
  filter(!is.na(trial_start_time_hms)) |> # Remove trials that never started
  count(stimulus_object, # Sample maybe okay: 17 Control (2 cens.), 12 Green ball (2 cens.), 9 Pink beads (4 cens.)
        appr_censored = approach_censor == 0, # Green ball: 2 cens., Pink beads: 1 cens.
        feed_censored = feed_censor == 0) # Green ball: 0 cens., Pink beads: 3 cens.

# [LTF] - Linear correlation between green ball and pink beads
obj_LTFmobl <- lmer(first_feed_latency_min ~ stimulus_object + (1 | nest_id/site/section/year), data = final_mobl, REML = F) # REML F for AIC ranking later
summary(obj_LTFmobl) # Corr. between objects = 0.331 but fully failing to converge with any random effect structure

# [LTA] - Linear correlation between green ball and pink beads
obj_LTAmobl <- lmer(log(first_approach_latency_min + 1) ~ stimulus_object + (1 | year), data = final_mobl, REML = F)
summary(obj_LTAmobl) # Corr. between objects = 0.314, REML = F needs super simplified random effects...

# [AFL] - Linear correlation between green ball and pink beads
obj_AFLmobl <- lmer(log(diff_approach_feed_min + 1) ~ stimulus_object + (1 | nest_id/site/section), data = final_mobl)
summary(obj_AFLmobl) # Corr. between objects = 0.297, simplified random effects

# For both species, responses to objects have low correlations.
# Use stimulus object going forward at least for WEBLs, but reduces sample sizes a lot for MOBLs

### How do stimulus object versus trial type models rank against each other? ###
# Look at feed for now

# [LTF] - Fit rest of the WEBL models
trial_LTFwebl <- glmmTMB(log(first_feed_latency_min + 1) ~ trial_type + (1 | nest_id/site/year), data = final_webl, REML = F)
stim_intLTFwebl <- glmmTMB(log(first_feed_latency_min + 1) ~ 1 + (1 | nest_id/site/year), data = final_webl, REML = F)

# [LTF] - Rank WEBL models
AIC(trial_LTFwebl, obj_LTFwebl, stim_intLTFwebl)
#                 df      AIC
# trial_LTFwebl    6 176.6671
# obj_LTFwebl      6 175.1705 - Object performs slightly better 
# stim_intLTFwebl  5 174.6735

# [LTA] - Fit rest of the MOBL models, look at LTA instead, convergence issues with LTF
trial_LTAmobl <- glmmTMB(log(first_approach_latency_min + 1) ~ trial_type + (1 | year), data = final_mobl, REML = F)
stim_intLTAmobl <- glmmTMB(log(first_approach_latency_min + 1) ~ 1 + (1 | year), data = final_mobl, REML = F)

# [LTA] - Rank WEBL models
AIC(trial_LTAmobl, obj_LTAmobl, stim_intLTAmobl)
#                 df      AIC
# trial_LTAmobl    4 109.3967 # Perform the same, trial type converges with full random effect structure...
# obj_LTAmobl      5 109.8686
# stim_intLTAmobl  3 107.7827

################################################################################
########################## Western Bluebird Modeling ###########################
################################################################################

########################## Adjust response latencies  ##########################

# Blind distance, nest age, number of chicks and canopy cover conceptually influence timings - might be too many, think more about this
# Including each as fixed effects over-complicates models so adjust latencies by using residuals from a glmm 
# nts: look into if this is appropriate for censored data. I think it's okay since the residuals from censored variables will still be considered censored 

# [LTF] - Adjust for blind distance, nest age and number of chicks, and canopy cover (might feel safer with more canopy cover and feed sooner?)
adjLTF_webl <- glmmTMB(log(first_feed_latency_min + 1) ~ scale(nest_age) + scale(chick_number) + scale(hide_distance_m) + scale(percent_canopy) + (1 | nest_id/site/section/year), data = final_webl)
summary(adjLTF_webl) 

# [LTF] - Look at residuals separately
weblLTF_resid <- resid(adjLTF_webl)
qqnorm(weblLTF_resid)
qqline(weblLTF_resid) # Looks really good

# [LTF] - Check other model diagnostics
plot(check_model(adjLTF_webl, panel = FALSE)) # Looks good

# [LTA] - Adjust for blind distance, nest age, number of chicks, and canopy cover
adjLTA_webl <- glmmTMB(log(first_approach_latency_min + 1) ~ scale(nest_age) + scale(chick_number) + scale(hide_distance_m) + scale(percent_canopy) + (1 | nest_id/site/section/year), data = final_webl)
summary(adjLTA_webl)

# [LTA] - Look at residuals separately
weblLTA_resid <- resid(adjLTA_webl)
qqnorm(weblLTA_resid)
qqline(weblLTA_resid) # Looks good 

# [LTA] - Check other model diagnostics
plot(check_model(adjLTA_webl, panel = FALSE)) # Looks good

# [AFL] - Adjust for blind distance, nest age, number of chicks, and canopy cover
adjAFL_webl <- glmmTMB(log(diff_approach_feed_min + 1) ~ scale(nest_age) + scale(chick_number) + scale(hide_distance_m) + scale(percent_canopy) + (1 | nest_id/site/section/year), data = final_webl)
summary(adjAFL_webl)

# [LTA] - Look at residuals separately
weblAFL_resid <- resid(adjAFL_webl)
qqnorm(weblAFL_resid)
qqline(weblAFL_resid) # Doesn't look great, S shaped
hist(weblAFL_resid) # Very right skewed

# [LTA] - Check other model diagnostics
plot(check_model(adjAFL_webl, panel = FALSE)) # Residuals look bad

# NOTES: LTF and LTA models seem to perform really well with log() transformed responses. 
# But, all of these fixed effects have very small z values. Canopy cover is almost negligible, with z's < |0.1|
# Can discount canopy cover if sticking with this log(x+1) transformation

# Add residuals for LTF and LTA into the main WEBL dataset (maybe re-write this using tidyverse)
final_webl$adj_LTF_min[as.numeric(names(residuals(adjLTF_webl)))] <- residuals(adjLTF_webl) # Got a warning but looks okay
final_webl$adj_LTA_min[as.numeric(names(residuals(adjLTA_webl)))] <- residuals(adjLTA_webl) # Got a warning but looks okay

hist(final_webl$adj_LTF_min) # Perfect
hist(final_webl$adj_LTA_min) # Perfect

############################# Fit Survival Models ##############################
# Other variables to consider: trial order, ordinal date, sex, individual ID (random)...
# Not confident, figure out how to check performance of these models

### SPL Meter Noise ###

# (1) LTF Cox Proportional-Hazard Model - using adjusted latencies
webl_LTF_coxme <- coxme(Surv(adj_LTF_min, feed_censor) ~ scale(LAeq_slm_dB)*scale(light_level_lux) + stimulus_object*scale(LAeq_slm_dB) +  stimulus_object*scale(light_level_lux) + (1 | nest_id/site/section/year), data = final_webl)
# ! - I think this is showing that for noisier nests with an object present (both pink beads and green ball), WEBLs are twice as likely to feed later relative to controls. 

webl_LTF_coxme2 <- coxme(Surv(adj_LTF_min, feed_censor) ~ scale(LAeq_slm_dB)*scale(light_level_lux) + trial_type*scale(LAeq_slm_dB) +  trial_type*scale(light_level_lux) + (1 | nest_id/site/section/year), data = final_webl)
# Same interaction

# (2) LTA Cox Proportional-Hazard Model - using adjusted latencies
webl_LTA_coxme <- coxme(Surv(adj_LTA_min, approach_censor) ~ scale(LAeq_slm_dB)*scale(light_level_lux) + stimulus_object*scale(LAeq_slm_dB) +  stimulus_object*scale(light_level_lux) + (1 | nest_id/site/section/year), data = final_webl)
# ! Same interaction between noise and pink beads as feed.

webl_LTA_coxme2 <- coxme(Surv(adj_LTA_min, approach_censor) ~ scale(LAeq_slm_dB)*scale(light_level_lux) + trial_type*scale(LAeq_slm_dB) +  trial_type*scale(light_level_lux) + (1 | nest_id/site/section/year), data = final_webl)
# Same interaction

# No difference between approach and feed yet until fixing adjustment model performance

### AudioMoth Noise ###
# ...

################################################################################
########################## Mountain Bluebird Modeling ##########################
################################################################################

# [LTF] - Adjust for blind distance, nest age and number of chicks. Adding canopy cover here (might feel safer with more canopy cover and feed sooner?)
adjLTF_mobl <- glmmTMB(first_feed_latency_min ~ scale(nest_age) + scale(chick_number) + scale(hide_distance_m) + scale(percent_canopy) + (1 | nest_id/site/section/year), data = final_mobl)
summary(adjLTF_mobl) 

# [LTF] - Look at residuals separately
moblLTF_resid <- resid(adjLTF_mobl)
qqnorm(moblLTF_resid)
qqline(moblLTF_resid) # Maybe looks okay? Tail goes off a bit

# [LTF] - Check other model diagnostics
plot(check_model(adjLTF_mobl, panel = FALSE)) # Okay

# [LTA] - Adjust for blind distance, nest age, number of chicks, and canopy cover
adjLTA_mobl <- glmmTMB(log(first_approach_latency_min + 1) ~ scale(nest_age) + scale(chick_number) + scale(hide_distance_m) + scale(percent_canopy) + (1 | nest_id/site/section/year), data = final_mobl)
summary(adjLTA_mobl)

# [LTA] - Look at residuals separately
moblLTA_resid <- resid(adjLTA_mobl)
qqnorm(moblLTA_resid)
qqline(moblLTA_resid) # Looks good 

# [LTA] - Check other model diagnostics
plot(check_model(adjLTA_mobl, panel = FALSE)) # Looks okay

# [AFL] - Adjust for blind distance, nest age, number of chicks, and canopy cover
adjAFL_mobl <- glmmTMB(log(diff_approach_feed_min + 1) ~ scale(nest_age) + scale(chick_number) + scale(hide_distance_m) + scale(percent_canopy) + (1 | nest_id/site/section/year), data = final_mobl)
summary(adjAFL_mobl)

# [LTA] - Look at residuals separately
moblAFL_resid <- resid(adjAFL_mobl)
qqnorm(moblAFL_resid)
qqline(moblAFL_resid) # Looks okay
hist(moblAFL_resid) # Looks okay

# [LTA] - Check other model diagnostics
plot(check_model(adjAFL_mobl, panel = FALSE)) # Okay

# NOTES: Fixed effects seem to have a stronger influence on MOBLs

# Add residuals for LTF and LTA into the main WEBL dataset (maybe re-write this using tidyverse)
final_mobl$adj_LTF_min[as.numeric(names(residuals(adjLTF_mobl)))] <- residuals(adjLTF_mobl) # Got a warning but looks okay
final_mobl$adj_LTA_min[as.numeric(names(residuals(adjLTA_mobl)))] <- residuals(adjLTA_mobl) # Got a warning but looks okay
final_mobl$adj_AFL_min[as.numeric(names(residuals(adjAFL_mobl)))] <- residuals(adjAFL_mobl) # Got a warning but looks okay

############################# Fit Cox mixed effect models ##############################
# Other variables to consider: trial order, ordinal date, sex, individual ID (random)...
# Not confident, figure out how to check performance of these models

### SPL Meter Noise ###

# (1) LTF Cox Proportional-Hazard Model - using adjusted latencies
mobl_LTF_coxme <- coxme(Surv(adj_LTF_min, feed_censor) ~ scale(LAeq_slm_dB)*scale(light_level_lux) + stimulus_object*scale(LAeq_slm_dB) +  stimulus_object*scale(light_level_lux) + (1 | nest_id/site/section/year), data = final_mobl)
# ! - Nothing jumps out as super strong but very weak interaction between noise and pink beads, where they approach sooner in noisy pink bead trials. Don't trust this too much given the sample size

# With trial type instead
mobl_LTF_coxme_2 <- coxme(Surv(adj_LTF_min, feed_censor) ~ scale(LAeq_slm_dB)*scale(light_level_lux) + trial_type*scale(LAeq_slm_dB) +  trial_type*scale(light_level_lux) + (1 | nest_id/site/section/year), data = final_mobl)
# Nothing - I think they two objects are cancelling each other out

# (2) LTA Cox Proportional-Hazard Model - using adjusted latencies
mobl_LTA_coxme <- coxme(Surv(adj_LTA_min, approach_censor) ~ scale(LAeq_slm_dB)*scale(light_level_lux) + stimulus_object*scale(LAeq_slm_dB) +  stimulus_object*scale(light_level_lux) + (1 | nest_id/site/section/year), data = final_mobl)
# Negative interaction between noise and light, approach sooner with green ball

mobl_LTA_coxme2 <- coxme(Surv(adj_LTA_min, approach_censor) ~ scale(LAeq_slm_dB)*scale(light_level_lux) + trial_type*scale(LAeq_slm_dB) +  trial_type*scale(light_level_lux) + (1 | nest_id/site/section/year), data = final_mobl)
# Nothing

# (3) AFL Cox Proportional-Hazard Model - using adjusted latencies
mobl_AFL_coxme <- coxme(Surv(adj_AFL_min, diff_censor) ~ scale(LAeq_slm_dB)*scale(light_level_lux) + stimulus_object*scale(LAeq_slm_dB) +  stimulus_object*scale(light_level_lux) + (1 | nest_id/site/section/year), data = final_mobl)
# Nothing strong, potentially twice as likely to take longer to feed in lit, green ball trials but t = 1.66

mobl_AFL_coxme2 <- coxme(Surv(adj_AFL_min, diff_censor) ~ scale(LAeq_slm_dB)*scale(light_level_lux) + trial_type*scale(LAeq_slm_dB) +  trial_type*scale(light_level_lux) + (1 | nest_id/site/section/year), data = final_mobl)
# Nothing

### AudioMoth Noise ###
# ...

################################################################################
############################### Between Species ################################
################################################################################

# difference between species models

################################################################################
########################### Secondary Analyses #################################
################################################################################

### Repeatability (need to double check how valid the current assumptions are) ###
# could also look across treatments

# [LTF] - What is the repeatability of latency within individual nests?
nest_rpt <- rpt(first_feed_latency_min ~ 1 + (1 | nest_id), grname = "nest_id", data = final_neo, datatype = "Gaussian")
nest_rpt
# Repeatability for nest_id
# R  = 0.256
# SE = 0.117
# CI = [0.008, 0.472]
# P  = 0.0129 [LRT]
# NA [Permutation]

# [LTF] - What is repeatability of latency across the same nest boxes (i.e. potentially same individuals across years)?
box_rpt <- rpt(first_feed_latency_min ~ 1 + (1 | nest_box), grname = "nest_box", data = final_neo, datatype = "Gaussian")
box_rpt
# Repeatability for nest_box
# R  = 0.169
# SE = 0.11
# CI = [0, 0.391]
# P  = 0.0503 [LRT]
# NA [Permutation]

# Very out of date glmm dump - survival is the right way to model these data instead
webl_feed <- lmer(first_feed_latency_min ~ LAeq_slm_dB*light_level_lux*trial_type + scale(hide_distance_m) + nest_age + chick_number + (1 | year/section/nest_id), data = final_webl)
summary(webl_feed)

plot(check_model(webl_feed, panel = FALSE))

webl_app <- lmer(first_approach_latency_s ~ LAeq_slm_dB*light_level_lux*trial_type + hide_distance_m + nest_age + chick_number + (1 | nest_id/site/section/year), data = final_webl)
summary(webl_app)

plot(check_model(webl_app, panel = FALSE))

webl_diff <- lmer(log(diff_approach_feed_s +1) ~ LAeq_slm_dB*light_level_lux*trial_type + hide_distance_m + nest_age + chick_number + (1 | nest_id/section/year), data = final_webl)
summary(webl_diff) # Failed to converge

mobl_feed <- lmer(first_feed_latency_min ~ LAeq_slm_dB*light_level_lux*trial_type + hide_distance_m + nest_age + chick_number + (1 | nest_id/section/year), data = final_mobl)
summary(mobl_feed)
resid(mobl_feed)

plot(check_model(mobl_feed, panel = FALSE))

mobl_app <- lmer(first_approach_latency_s ~ LAeq_slm_dB*light_level_lux*stimulus_object + hide_distance_m + nest_age + chick_number + (1 | nest_id/section/year), data = final_mobl)
summary(mobl_app)

plot(check_model(mobl_app, panel = FALSE))

# Distribution for this contains too many 0s, model fails to converge. 
# But this is interesting, 
mobl_diff <- lmer(diff_approach_feed_s ~ LAeq_slm_dB*light_level_lux*stimulus_object + hide_distance_m + nest_age + chick_number + (1 | nest_id/section/year), data = final_mobl)
summary(mobl_diff) # Failed to converge
