#####################################################################################################

######## Analysis - Main
# (only with final_data_VIF)

# 1) Multiple Imputation 
# 2) Multi-level mixed effects models (metafor)
#     - Full model (baseline + indicators; including variance component pooling, I^2 calculation)
#     - THEORY-LEVEL models
#     - STATE AND REGION-LEVEL MODELS FOR FOREST PLOTS
#     - INTERACTION MODELS (EXPLORATORY)
#     - Sanity checks: LMER and LM models
# 3) Plotting - Preprocessing
#     - Extract effect sizes for forest plots
#     - Calculate number of treatments, studies and participants per state and region
# 4) Plotting (state- and region- level)
#     - Density plots of cooperation rates
#     - Forest plots
#     - Scatter plots: state/region-level indicators against state/region-level P(C)s controlled for baseline
# 5) Other analyses
#     - Flip coding of discussion
#     - FOREST PLOTS WITH AND WITHOUT DISCUSSION
#     - Kruskall-Wallis Test á la Gächter et al.
#     - Plotting of effect sizes


### Load packages
library("tidyverse")
library("mice")
library("metafor")
library("lme4")
library("faraway")
library("reshape2")
library("ggplot2")
library("ggrepel")

eval(metafor:::.mice)

# Set seed
set.seed(11)

### Set WD
setwd("/home/caroline/Desktop/cross_cultural_metaanalysis/US/analysis")

### Read final_data_VIF
final_data_VIF = read.csv("data_frames/final_data_VIF.csv")

### Convert dichotomous variables to factors
final_data_VIF$matching = factor(final_data_VIF$matching)
final_data_VIF$discussion = factor(final_data_VIF$discussion)
final_data_VIF$simultaneous = factor(final_data_VIF$simultaneous)
final_data_VIF$sanction = factor(final_data_VIF$sanction)
final_data_VIF$symmetry = factor(final_data_VIF$symmetry)
final_data_VIF$one_shot = factor(final_data_VIF$one_shot)

#####################################################################################################

#### Imputation with MICE package 

imp <- mice(final_data_VIF, m=2, maxit=0, seed=11)
pred <- imp$predictorMatrix #get the predictor matrix

#Don't use some variables for imputation:
pred[,"year"] <- 0 # don't use year_data_collection for imputing
pred[,"N"] <- 0 # don't use N for imputing
pred[,"P_C"] <- 0 # don't use P_C for imputing
pred[,"state"] <- 0 # don't use state for imputing
pred[,"study_ID"] <- 0 # don't use study_ID for imputing
pred[,"region"] <- 0 # don't use region for imputing
pred[,"yi"] <- 0 # don't use yi for imputing
pred[,"vi"] <- 0 # don't use vi for imputing

pred["study_ID",] <- 0 # don't impute study_ID
pred["state",] <- 0 # don't impute state
pred["region",] <- 0 # don't impute region
pred["P_C",] <- 0 # don't impute P_C
pred["N",] <- 0 # don't impute N
pred["yi",] <- 0 # don't impute yi
pred["vi",] <- 0 # don't impute vi

#generate the multiple imputations
final_data_VIF_imputed = mice(final_data_VIF, m = 5, maxit = 5, predictorMatrix=pred, seed=11)


#############################################################################################################################

### Models

##############################################################################################################################

# Multi-level mixed effects models (metafor)

### 1) ALL DATA, baseline + indicators - WITH VIF

# no multicollinear variables: PC1_disease_mortality_and_violent_crime, residential_mob, tightness_zscore, PC1_government_effectiveness
fit_full_VIF <- with(final_data_VIF_imputed, rma.mv(yi, vi, mods= ~ symmetry+one_shot+matching+group_size+if_PD_K_index+discussion+simultaneous+sanction+trust_reversed_zscore+subsistence_index+self_expression+church_exp_West+church_exp_East+percent_farms_operated_by_owner_zscore+collectivism+historical_subsistence_index+PC_religion+PC1_cousin_marriage+PC2_voter_turnout_and_corruption_control_perception+PC3_institutional_trust_and_perception_of_electoral_integrity+PC4_state_pension_funding_and_corruption_convictions+PC5_economic_wealth_and_inequality_and_health_insurance+PC2_property_crime_and_caloric_suitability+PC3_disaster_freq_cost_and_demanding_geoclimate, method="REML",
                                                    random = list(~ 1 | region,~ 1 | state,~ 1 | study_ID)))
round(summary(pool(fit_full_VIF)), 4)

# No state-level indicators significant (but religion marginally significant: p=0.0856)
# sign. structural chars: group_size, if_PD_K_index, discussion, simultaneous

### Variance component pooling for rma.mv
# - In Package "mitml-Tools for Multiple Imputation in Multilevel Modeling" (see Example 3: manual extraction of variance estimates in documentation of with.mitml.list and testEstimates), indices (e.g., within-group variation, but I^2 would also be example) are computed for each data set separately, and then simply averaged.
# - 1) Pool sigma2: simply average sigma2 values (exactly same result as calulating variance and using pool.scalar(sigma2_state, sigma2_state_vars_vec))
# - 2) Calculate I2 for each model and pool I2 (i.e., simply average)

# - 1) Pool sigma2: simply average sigma2 values

sigma2_region = mean(c(fit_full_VIF$analyses[[1]]$sigma2[1], 
                       fit_full_VIF$analyses[[2]]$sigma2[1], 
                       fit_full_VIF$analyses[[3]]$sigma2[1], 
                       fit_full_VIF$analyses[[4]]$sigma2[1], 
                       fit_full_VIF$analyses[[5]]$sigma2[1]))
sigma2_region
sigma2_state = mean(c(fit_full_VIF$analyses[[1]]$sigma2[2], 
                      fit_full_VIF$analyses[[2]]$sigma2[2], 
                      fit_full_VIF$analyses[[3]]$sigma2[2], 
                      fit_full_VIF$analyses[[4]]$sigma2[2], 
                      fit_full_VIF$analyses[[5]]$sigma2[2]))
sigma2_state
sigma2_study_ID = mean(c(fit_full_VIF$analyses[[1]]$sigma2[3], 
                         fit_full_VIF$analyses[[2]]$sigma2[3], 
                         fit_full_VIF$analyses[[3]]$sigma2[3], 
                         fit_full_VIF$analyses[[4]]$sigma2[3], 
                         fit_full_VIF$analyses[[5]]$sigma2[3]))
sigma2_study_ID

# - 2) Calculate I2 for each model and pool I2 (i.e., simply average)

W <- diag(1/final_data_VIF$vi)

# Model 1: 
X <- model.matrix(fit_full_VIF$analyses[[1]])
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_1 = 100 * fit_full_VIF$analyses[[1]]$sigma2 / (sum(fit_full_VIF$analyses[[1]]$sigma2) + (fit_full_VIF$analyses[[1]]$k-fit_full_VIF$analyses[[1]]$p)/sum(diag(P)))

# Model 2: 
X <- model.matrix(fit_full_VIF$analyses[[2]])
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_2 = 100 * fit_full_VIF$analyses[[2]]$sigma2 / (sum(fit_full_VIF$analyses[[2]]$sigma2) + (fit_full_VIF$analyses[[2]]$k-fit_full_VIF$analyses[[2]]$p)/sum(diag(P)))

# Model 3: 
X <- model.matrix(fit_full_VIF$analyses[[3]])
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_3 = 100 * fit_full_VIF$analyses[[3]]$sigma2 / (sum(fit_full_VIF$analyses[[3]]$sigma2) + (fit_full_VIF$analyses[[3]]$k-fit_full_VIF$analyses[[3]]$p)/sum(diag(P)))

# Model 4: 
X <- model.matrix(fit_full_VIF$analyses[[4]])
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_4 = 100 * fit_full_VIF$analyses[[4]]$sigma2 / (sum(fit_full_VIF$analyses[[4]]$sigma2) + (fit_full_VIF$analyses[[4]]$k-fit_full_VIF$analyses[[4]]$p)/sum(diag(P)))

# Model 5: 
X <- model.matrix(fit_full_VIF$analyses[[5]])
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_5 = 100 * fit_full_VIF$analyses[[5]]$sigma2 / (sum(fit_full_VIF$analyses[[5]]$sigma2) + (fit_full_VIF$analyses[[5]]$k-fit_full_VIF$analyses[[5]]$p)/sum(diag(P)))

I2_between_regions = mean(c(I2_1[1], I2_2[1], I2_3[1], I2_4[1], I2_5[1]))
I2_between_states = mean(c(I2_1[2], I2_2[2], I2_3[2], I2_4[2], I2_5[2]))
I2_within_states = mean(c(I2_1[3], I2_2[3], I2_3[3], I2_4[3], I2_5[3]))
I2_sampling_var = mean(c( (100 - sum(I2_1)), (100 - sum(I2_2)), (100 - sum(I2_3)), (100 - sum(I2_4)), (100 - sum(I2_5)) ))
paste(I2_between_regions, I2_between_states, I2_within_states, I2_sampling_var)
# About 0% of the total variance is estimated to be due to between-region heterogeneity, 
# about 6.05% of the total variance is estimated to be due to between-state heterogeneity, 
# about 74.78% are due to heterogeneity of studies within states. 
# The still remaining 19.17% are sampling variance.

### !!! BIG differences in explained variance compared to non-VIF model (which had 0%, 3.6%, 76.9%, 19.5%)
### Much more (twice as much!) explained variance by state!

### Plot it

percent_variance = c(I2_between_regions, I2_between_states, I2_within_states, I2_sampling_var)
percent_variance = as.data.frame(percent_variance)
percent_variance$Level = c("Between regions: 0%", "Between states: 6.05%", "Within states: 74.78%", "Sampling variance: 19.17%")
percent_variance$temp = c(" ", " ", " ", " ")
percent_variance$percent_variance = percent_variance$percent_variance / 100
percent_variance$Level = factor(percent_variance$Level, levels = c("Between regions: 0%", "Between states: 6.05%", "Within states: 74.78%", "Sampling variance: 19.17%"))

# Stacked barplot with multiple groups
ggplot(data=percent_variance, aes(x=temp, y=percent_variance, fill=Level)) +
  geom_bar(stat="identity")+
  labs(x=" ", y = "Percent of total variance")+
  ggsave("plots/I2_effects.png", width=4, height=5)
  
#######################################################################################################################


##### THEORY-LEVEL models

### 1a) baseline + religion

fit_VIF_R <- with(final_data_VIF_imputed, rma.mv(yi, vi, mods= ~ symmetry+one_shot+matching+group_size+if_PD_K_index+discussion+simultaneous+sanction+PC_religion, method="REML",
                                                 random = list(~ 1 | region,~ 1 | state, ~ 1 | study_ID)))
round(summary(pool(fit_VIF_R)), 4) #, conf.int = TRUE

# Religion significant! - More R, more C

### 1b) baseline + formal instis

fit_VIF_FI <- with(final_data_VIF_imputed, rma.mv(yi, vi, mods= ~ symmetry+one_shot+matching+group_size+if_PD_K_index+discussion+simultaneous+sanction+PC2_voter_turnout_and_corruption_control_perception+PC3_institutional_trust_and_perception_of_electoral_integrity+PC4_state_pension_funding_and_corruption_convictions+PC5_economic_wealth_and_inequality_and_health_insurance, method="REML",
                                                  random = list(~ 1 | region,~ 1 | state, ~ 1 | study_ID)))
round(summary(pool(fit_VIF_FI)), 4) #, conf.int = TRUE

### 1c) baseline + informal instis

fit_VIF_II <- with(final_data_VIF_imputed, rma.mv(yi, vi, mods= ~ symmetry+one_shot+matching+group_size+if_PD_K_index+discussion+simultaneous+sanction+PC1_cousin_marriage, method="REML",
                                                  random = list(~ 1 | region,~ 1 | state, ~ 1 | study_ID)))
round(summary(pool(fit_VIF_II)), 4) #, conf.int = TRUE

### 1d) baseline + subsistence

fit_VIF_subs <- with(final_data_VIF_imputed, rma.mv(yi, vi, mods= ~ symmetry+one_shot+matching+group_size+if_PD_K_index+discussion+simultaneous+sanction+percent_farms_operated_by_owner_zscore+subsistence_index+historical_subsistence_index, method="REML",
                                                    random = list(~ 1 | region,~ 1 | state, ~ 1 | study_ID)))
round(summary(pool(fit_VIF_subs)), 4) #, conf.int = TRUE

### 1e) baseline + threat

fit_VIF_threat <- with(final_data_VIF_imputed, rma.mv(yi, vi, mods= ~ symmetry+one_shot+matching+group_size+if_PD_K_index+discussion+simultaneous+sanction+PC2_property_crime_and_caloric_suitability+PC3_disaster_freq_cost_and_demanding_geoclimate, method="REML",
                                                      random = list(~ 1 | region,~ 1 | state, ~ 1 | study_ID)))
round(summary(pool(fit_VIF_threat)), 4) #, conf.int = TRUE
# PC3_disaster_freq_cost_and_demanding_geoclimate marginally significant! higher threat, more C

### 1f) baseline + church exposure

fit_VIF_church <- with(final_data_VIF_imputed, rma.mv(yi, vi, mods= ~ symmetry+one_shot+matching+group_size+if_PD_K_index+discussion+simultaneous+sanction+church_exp_West+church_exp_East, method="REML",
                                                      random = list(~ 1 | region,~ 1 | state, ~ 1 | study_ID)))
round(summary(pool(fit_VIF_church)), 4) #, conf.int = TRUE

# church_exp_West significant! 
# Less church exposure West - more C

### 1g) baseline + collectivism

fit_VIF_collectivism <- with(final_data_VIF_imputed, rma.mv(yi, vi, mods= ~ symmetry+one_shot+matching+group_size+if_PD_K_index+discussion+simultaneous+sanction+collectivism, method="REML",
                                                            random = list(~ 1 | region,~ 1 | state, ~ 1 | study_ID)))
round(summary(pool(fit_VIF_collectivism)), 4) #, conf.int = TRUE

# collectivism significant! 
# More collectivism - more C

### 1h) baseline + self_expression

fit_VIF_self_expression <- with(final_data_VIF_imputed, rma.mv(yi, vi, mods= ~ symmetry+one_shot+matching+group_size+if_PD_K_index+discussion+simultaneous+sanction+self_expression, method="REML",
                                                               random = list(~ 1 | region,~ 1 | state, ~ 1 | study_ID)))
round(summary(pool(fit_VIF_self_expression)), 4) #, conf.int = TRUE

# self_expression significant! 
# Less self_expression - more C

### 1i) baseline + trust

fit_VIF_trust <- with(final_data_VIF_imputed, rma.mv(yi, vi, mods= ~ symmetry+one_shot+matching+group_size+if_PD_K_index+discussion+simultaneous+sanction+trust_reversed_zscore, method="REML",
                                                     random = list(~ 1 | region,~ 1 | state, ~ 1 | study_ID)))
round(summary(pool(fit_VIF_trust)), 4) #, conf.int = TRUE

# Trust marginally significant: Less trust - more C

#### Adjust p-values for multiple comparisons
p_values_sign = c(0.0158, 0.0464, 0.0245, 0.0432, 0.0206, 0.0492)
p_values_all = c(0.0158, 0.0464, 0.0245, 0.0432, 0.0206, 0.0492, 0.1206, 0.5059, 0.3083, 0.9317, 0.1244, 0.2074, 0.1150, 0.7185, 0.0808, 0.1983)
p.adjust(p_values_sign, method="fdr")
p.adjust(p_values_all, method="fdr")
p.adjust(p_values_sign, method="bonferroni")
p.adjust(p_values_all, method="bonferroni")

# Nothing significant after correction for multiple hypothesis testing (both with false discovery rate, and with Bonferroni correction)



### 1j) ALL DATA, baseline + indicators + year INTERACTIONS - WITH VIF

# Year interactions: 
fit_full_VIF_year_interactions <- with(final_data_VIF_imputed, rma.mv(yi, vi, mods= ~ symmetry+one_shot+matching+group_size+if_PD_K_index+discussion+simultaneous+sanction+year*church_exp_West+year*church_exp_East+year*percent_farms_operated_by_owner_zscore+year*collectivism+year*self_expression+year*subsistence_index+year*historical_subsistence_index+year*PC_religion+year*PC2_voter_turnout_and_corruption_control_perception+year*PC3_institutional_trust_and_perception_of_electoral_integrity+year*PC4_state_pension_funding_and_corruption_convictions+year*PC5_economic_wealth_and_inequality_and_health_insurance+year*PC1_cousin_marriage+year*PC2_property_crime_and_caloric_suitability+year*PC3_disaster_freq_cost_and_demanding_geoclimate, method="REML",
                                                                      random = list(~ 1 | region,~ 1 | state, ~ 1 | study_ID)))
round(summary(pool(fit_full_VIF_year_interactions)), 4) #, conf.int = TRUE

# Marginally significant effects: PC_religion, year, year:PC_religion

# Other interaction models non-sensical because do not converge
# maximum interactions that can be included with convergence is 1 indicator x 27 (all other indicators) (e.g., year x all other indicators)


##############################################################################################################################

####### STATE AND REGION-LEVEL MODELS FOR FOREST PLOTS

### 2) state as moderator + baseline

fit_state <- with(final_data_VIF_imputed, rma.mv(yi, vi, mods= ~ state+symmetry+one_shot+matching+group_size+if_PD_K_index+discussion+simultaneous+sanction, method="REML",
                                       random = ~ 1 | study_ID))
pool_state <- pool(fit_state)
round(summary(pool_state, conf.int=TRUE), 4)

# No state-level indicators sign
# sign. structural chars: matching, group_size, if_PD_K_index, discussion, simultaneous


### 3) Region as moderator + baseline

fit_region <- with(final_data_VIF_imputed, rma.mv(yi, vi, mods= ~ region+symmetry+one_shot+matching+group_size+if_PD_K_index+discussion+simultaneous+sanction, method="REML",
                                       random = ~ 1 | study_ID))
pool_region <- pool(fit_region)
round(summary(pool_region, conf.int=TRUE), 4)


### 4) Only baseline

fit_baseline <- with(final_data_VIF_imputed, rma.mv(yi, vi, mods= ~ symmetry+one_shot+matching+group_size+if_PD_K_index+discussion+simultaneous+sanction, method="REML",
                                               random = ~ 1 | study_ID))
pool_baseline <- pool(fit_baseline)
round(summary(pool_baseline, conf.int=TRUE), 4)


### 5) No moderators (to get observed yi)

fit_null <- with(final_data_VIF_imputed, rma.mv(yi, vi, method="REML", random = ~ 1 | study_ID))
pool_null <- pool(fit_null)
round(summary(pool_null, conf.int=TRUE), 4)


##############################################################################################################

### Other models ###

### 1) Sanity check: LMER (yi non-meta-analysis)

fit_yi <- with(final_data_VIF_imputed, lmer(yi ~ symmetry+one_shot+matching+group_size+if_PD_K_index+discussion+simultaneous+sanction+church_exp_West+church_exp_East+percent_farms_operated_by_owner_zscore+collectivism+self_expression+subsistence_index+historical_subsistence_index+trust_reversed_zscore+PC_religion+PC1_cousin_marriage+PC2_voter_turnout_and_corruption_control_perception+PC3_institutional_trust_and_perception_of_electoral_integrity+PC4_state_pension_funding_and_corruption_convictions+PC5_economic_wealth_and_inequality_and_health_insurance+PC2_property_crime_and_caloric_suitability+PC3_disaster_freq_cost_and_demanding_geoclimate + (1 | region) + (1 | state) + (1 | study_ID), REML = TRUE))
pool_yi <- pool(fit_yi)
round(summary(pool_yi), 4) #, conf.int = TRUE

# Still nothing significant
# (From structural vars: everything significant)

################################################################################################################

#### Plotting - Preprocessing

################################################################################################################

### Extract effect sizes

# state-level
state_intercept = summary(pool_state)$estimate[1]

states <- data.frame(state_level_effects_controlled_for_baseline=rep.int(NA, 38),
                        state_level_P_C_controlled_for_baseline = rep.int(NA, 38),
                        state=rep.int(NA, 38),
                        estimate_conf_lower = rep.int(NA, 38),
                        estimate_conf_upper = rep.int(NA, 38)) 

states$state_level_effects_controlled_for_baseline_uncorrected_for_intercept = summary(pool_state)$estimate[1:38]
states$state_level_effects_controlled_for_baseline[2:38] = state_intercept + states$state_level_effects_controlled_for_baseline_uncorrected_for_intercept[2:38]
states$state_level_effects_controlled_for_baseline[1] = state_intercept # first state IS baseline, so = INTERCEPT
states$state_level_P_C_controlled_for_baseline = transf.ilogit.int(states$state_level_effects_controlled_for_baseline)
states$state = c("Alabama", "Arizona", "Arkansas", "California", "Connecticut", "Delaware", "District of Columbia", "Florida", "Georgia", "Hawaii", "Illinois", "Indiana", "Iowa", "Kansas", "Maryland", "Massachusetts", "Michigan", "Minnesota", "Mississippi", "Missouri", "Nevada", "New Jersey", "New York", "North Carolina", "Ohio", "Oklahoma", "Oregon", "Pennsylvania", "Rhode Island", "South Carolina", "Tennessee", "Texas", "Utah", "Vermont", "Virginia", "Washington", "Wisconsin", "Wyoming")
states$std_error = summary(pool_state)$std.error[1:38]
# Add std.error estimates for P(C)
states$estimate_conf_lower = states$state_level_effects_controlled_for_baseline - states$std_error
states$estimate_conf_upper = states$state_level_effects_controlled_for_baseline + states$std_error
states$estimate_conf_lower_P_C = transf.ilogit.int(states$estimate_conf_lower)
states$estimate_conf_upper_P_C = transf.ilogit.int(states$estimate_conf_upper)

# region-level
region_intercept = summary(pool_region)$estimate[1]

regions <- data.frame(region_level_effects_controlled_for_baseline=rep.int(NA, 9),
                       region_level_P_C_controlled_for_baseline = rep.int(NA, 9),
                       region=rep.int(NA, 9)) 

regions$region_level_effects_controlled_for_baseline_uncorrected_for_intercept = summary(pool_region)$estimate[1:9]
regions$region_level_effects_controlled_for_baseline[2:9] = region_intercept + regions$region_level_effects_controlled_for_baseline_uncorrected_for_intercept[2:9]
regions$region_level_effects_controlled_for_baseline[1] = region_intercept # first region IS baseline, so = INTERCEPT
regions$region_level_P_C_controlled_for_baseline = transf.ilogit.int(regions$region_level_effects_controlled_for_baseline)
regions$region = c("East North Central", "East South Central", "Middle Atlantic", "Mountain", "New England", "Pacific", "South Atlantic", "West North Central", "West South Central")
regions$std_error = summary(pool_region)$std.error[1:9]
# Add std.error estimates for P(C)
regions$estimate_conf_lower = regions$region_level_effects_controlled_for_baseline - regions$std_error
regions$estimate_conf_upper = regions$region_level_effects_controlled_for_baseline + regions$std_error
regions$estimate_conf_lower_P_C = transf.ilogit.int(regions$estimate_conf_lower)
regions$estimate_conf_upper_P_C = transf.ilogit.int(regions$estimate_conf_upper)

# Get overall cooperation rate controlled for baseline
baseline_intercept_yi = summary(pool_baseline)$estimate[1]
baseline_intercept_P_C = transf.ilogit.int(baseline_intercept_yi)

# Get overall cooperation rate NOT controlled for baseline ("observed")
observed_intercept_yi = summary(pool_null)$estimate[1]
observed_intercept_P_C = transf.ilogit.int(observed_intercept_yi)

# Add k and N to region label
regions$region_name = as.factor(c("East North Central; k=88; N=6305", "East South Central; k=10; N=1312", "Middle Atlantic; k=73; N=7074", "Mountain; k=15; N=1217", "New England; k=27; N=1906", "Pacific; k=71; N=7953", "South Atlantic; k=83; N=7528", "West North Central; k=13; N=1219", "West South Central; k=14; N=385"))
states$state_name = as.factor(c("Alabama; k=6; N=1102", "Arizona; k=8; N=587", "Arkansas; k=1; N=31", "California; k=50; N=5772", "Connecticut; k=3; N=230", "Delaware; k=6; N=1022", "District of Columbia; k=1; N=70", "Florida; k=20; N=1535", "Georgia; k=6; N=865", "Hawaii; k=4; N=462", "Illinois; k=40; N=3001", "Indiana; k=1; N=100", "Iowa; k=4; N=180", "Kansas; k=3; N=184", "Maryland; k=1; N=48", "Massachusetts; k=11; N=987", "Michigan; k=24; N=1867", "Minnesota; k=1; N=200", "Mississippi; k=3; N=120", "Missouri; k=5; N=655", "Nevada; k=1; N=96", "New Jersey; k=12; N=482", "New York; k=34; N=4351", "North Carolina; k=35; N=2579", "Ohio; k=22; N=1209", "Oklahoma; k=1; N=18", "Oregon; k=5; N=536", "Pennsylvania; k=27; N=2241", "Rhode Island; k=11; N=633", "South Carolina; k=12; N=1150", "Tennessee; k=1; N=90", "Texas; k=12; N=336", "Utah; k=3; N=271", "Vermont; k=2; N=56", "Virginia; k=2; N=259", "Washington; k=12; N=1183", "Wisconsin; k=1; N=128", "Wyoming; k=3; N=263"))


# Save df
# write.csv(regions, "data_frames/regions.csv", row.names = FALSE)

###########################################################################################################

## Number of treatments, studies and participants per state and region

# Load final coda data
final_coda_data = read.csv("../data/final_coda_data_US.csv")

# 1) number of treatments

num_treatments_region = final_data_VIF%>%
  filter(!is.na(region))%>%
  group_by(region)%>%
  summarise(num_treatments = n())

num_treatments_state = final_data_VIF%>%
  group_by(state)%>%
  summarise(num_treatments = n())

# 2) Number of studies and participants

# Attach BS/WS info to calculate num subjects correctly
final_data_VIF_BS_WS = final_data_VIF%>%
  mutate(BS_WS = final_coda_data$BS_WS)

# Delete duplicates in WS
WS = filter(final_data_VIF_BS_WS, BS_WS==2)
WS = WS[!duplicated(WS$study_ID), ]

BS = filter(final_data_VIF_BS_WS, !BS_WS==2 | is.na(BS_WS))
BS_WS_combined = bind_rows(BS, WS)

num_studies_region = BS_WS_combined%>%
  filter(!is.na(region))%>%
  group_by(region)%>%
  summarise(num_studies = n(),
            total_n = sum(N))

sum(num_studies_region$total_n)
# 34899

num_studies_state = BS_WS_combined%>%
  group_by(state)%>%
  summarise(num_studies = n(),
            total_n = sum(N))

sum(num_studies_state$total_n)
# 34899

num_studies_region = left_join(num_studies_region, num_treatments_region, by="region")
num_studies_state = left_join(num_studies_state, num_treatments_state, by="state")

# write.csv(num_studies_region, "data_frames/num_studies_region.csv", row.names = FALSE)
# write.csv(num_studies_state, "data_frames/num_studies_state.csv", row.names = FALSE)



################################################################################################################

#### PLOTTING

################################################################################################################

# A) Density plots of cooperation rates

final_data_VIF_with_means = left_join(final_data_VIF, states, by="state")
final_data_VIF_with_means = left_join(final_data_VIF_with_means, regions, by="region")

final_data_VIF_with_means = rename(final_data_VIF_with_means, Region=region_name, State=state_name)

# write.csv(final_data_VIF_with_means, "data_frames/final_data_VIF_with_means.csv", row.names = FALSE)


ggplot(final_data_VIF_with_means, aes(x=P_C, color=Region, fill=Region)) + 
  geom_density(alpha=0.3)+
  #geom_point(aes(x=region_level_P_C_controlled_for_baseline, y=0, size=2))+ doesnt make sense to compare baseline-controlled PC with observed PC
  labs(x = "Observed P(C)", y = "Density") +
  ggsave("plots/density_plot_observed_P_C_region.png", width=9, height=5.5)
ggplot(final_data_VIF_with_means, aes(x=P_C, color=State, fill=State)) + 
  geom_density(alpha=0.3)+
  #geom_point(aes(x=state_level_P_C_controlled_for_baseline, y=0, size=2))+ doesnt make sense to compare baseline-controlled PC with observed PC
  labs(x = "Observed P(C)", y = "Density") +
  ggsave("plots/density_plot_observed_P_C_state.png", width=9, height=5.5)



#### B) Forest plots

# 1) yis: state-level

# Calculate +- std.error
states$std_error_yi_upper = states$state_level_effects_controlled_for_baseline+ states$std_error
states$std_error_yi_lower = states$state_level_effects_controlled_for_baseline- states$std_error

# Set factor ordering to preserve alphabetical state order despite flipping
states$state_name = as.factor(c("Alabama", "Arizona", "Arkansas", "California", "Connecticut", "Delaware", "District of Columbia", "Florida", "Georgia", "Hawaii", "Illinois", "Indiana", "Iowa", "Kansas", "Maryland", "Massachusetts", "Michigan", "Minnesota", "Mississippi", "Missouri", "Nevada", "New Jersey", "New York", "North Carolina", "Ohio", "Oklahoma", "Oregon", "Pennsylvania", "Rhode Island", "South Carolina", "Tennessee", "Texas", "Utah", "Vermont", "Virginia", "Washington", "Wisconsin", "Wyoming"))
states$state_name <- factor(states$state_name, levels = states$state_name[order(states$state_name, decreasing = TRUE)])
states$state_name_with_k = as.factor(c("Alabama; k=6", "Arizona; k=8", "Arkansas; k=1", "California; k=50", "Connecticut; k=3", "Delaware; k=6", "District of Columbia; k=1", "Florida; k=20", "Georgia; k=6", "Hawaii; k=4", "Illinois; k=40", "Indiana; k=1", "Iowa; k=4", "Kansas; k=3", "Maryland; k=1", "Massachusetts; k=11", "Michigan; k=24", "Minnesota; k=1", "Mississippi; k=3", "Missouri; k=5", "Nevada; k=1", "New Jersey; k=12", "New York; k=34", "North Carolina; k=35", "Ohio; k=22", "Oklahoma; k=1", "Oregon; k=5", "Pennsylvania; k=27", "Rhode Island; k=11", "South Carolina; k=12", "Tennessee; k=1", "Texas; k=12", "Utah; k=3", "Vermont; k=2", "Virginia; k=2", "Washington; k=12", "Wisconsin; k=1", "Wyoming; k=3"))
states$state_name_with_k <- factor(states$state_name_with_k, levels = states$state_name_with_k[order(states$state_name_with_k, decreasing = TRUE)])
states$state_name_with_k_N = as.factor(c("Alabama; k=6; N=1102", "Arizona; k=8; N=587", "Arkansas; k=1; N=31", "California; k=50; N=5772", "Connecticut; k=3; N=230", "Delaware; k=6; N=1022", "District of Columbia; k=1; N=70", "Florida; k=20; N=1535", "Georgia; k=6; N=865", "Hawaii; k=4; N=462", "Illinois; k=40; N=3001", "Indiana; k=1; N=100", "Iowa; k=4; N=180", "Kansas; k=3; N=184", "Maryland; k=1; N=48", "Massachusetts; k=11; N=987", "Michigan; k=24; N=1867", "Minnesota; k=1; N=200", "Mississippi; k=3; N=120", "Missouri; k=5; N=655", "Nevada; k=1; N=96", "New Jersey; k=12; N=482", "New York; k=34; N=4351", "North Carolina; k=35; N=2579", "Ohio; k=22; N=1209", "Oklahoma; k=1; N=18", "Oregon; k=5; N=536", "Pennsylvania; k=27; N=2241", "Rhode Island; k=11; N=633", "South Carolina; k=12; N=1150", "Tennessee; k=1; N=90", "Texas; k=12; N=336", "Utah; k=3; N=271", "Vermont; k=2; N=56", "Virginia; k=2; N=259", "Washington; k=12; N=1183", "Wisconsin; k=1; N=128", "Wyoming; k=3; N=263"))
states$state_name_with_k_N <- factor(states$state_name_with_k_N, levels = states$state_name_with_k_N[order(states$state_name_with_k_N, decreasing = TRUE)])

ggplot(data=states, aes(color=state_name)) + 
  geom_pointrange(mapping=aes(x=states$state_name, y=states$state_level_effects_controlled_for_baseline, ymin=states$std_error_yi_lower , ymax=states$std_error_yi_upper), size=0.5)+
  coord_flip() + 
  geom_hline(yintercept = baseline_intercept_yi, linetype="twodash", color = "blue", size=0.5) +
  theme(legend.position = "none")+
  xlab("State") +
  ylab("Yi controlled for baseline")+
  ggsave("plots/forest_plot_state_yi_controlled_for_baseline.png", width = 8, height = 6)

ggplot(data=states, aes(color=state_name_with_k)) + 
  geom_pointrange(mapping=aes(x=states$state_name_with_k, y=states$state_level_effects_controlled_for_baseline, ymin=states$std_error_yi_lower , ymax=states$std_error_yi_upper), size=0.5)+
  coord_flip() + 
  geom_hline(yintercept = baseline_intercept_yi, linetype="twodash", color = "blue", size=0.5) +
  theme(legend.position = "none")+
  xlab("State") +
  ylab("Yi controlled for baseline")+
  ggsave("plots/forest_plot_state_yi_controlled_for_baseline_with_k.png", width = 9, height = 6)

ggplot(data=states, aes(color=state_name_with_k_N)) + 
  geom_pointrange(mapping=aes(x=states$state_name_with_k_N, y=states$state_level_effects_controlled_for_baseline, ymin=states$std_error_yi_lower , ymax=states$std_error_yi_upper), size=0.5)+
  coord_flip() + 
  geom_hline(yintercept = baseline_intercept_yi, linetype="twodash", color = "blue", size=0.5) +
  theme(legend.position = "none")+
  xlab("State") +
  ylab("Yi controlled for baseline")+
  ggsave("plots/forest_plot_state_yi_controlled_for_baseline_with_k_N.png", width = 10, height = 6)

# 2) yis: region-level

regions$std_error_yi_upper = regions$region_level_effects_controlled_for_baseline+ regions$std_error
regions$std_error_yi_lower = regions$region_level_effects_controlled_for_baseline- regions$std_error

# Use better names, convert to factor and make factor ordering preserve alphabetical state order despite flipping
regions$region_names = as.factor(c("East North Central", "East South Central", "Middle Atlantic", "Mountain", "New England", "Pacific", "South Atlantic", "West North Central", "West South Central"))
regions$region_names <- factor(regions$region_names, levels = regions$region_names[order(regions$region_names, decreasing = TRUE)])
regions$region_names_with_k = as.factor(c("East North Central; k=88", "East South Central; k=10", "Middle Atlantic; k=73", "Mountain; k=15", "New England; k=27", "Pacific; k=71", "South Atlantic; k=83", "West North Central; k=13", "West South Central; k=14"))
regions$region_names_with_k <- factor(regions$region_names_with_k, levels = regions$region_names_with_k[order(regions$region_names_with_k, decreasing = TRUE)])
regions$region_names_with_k_N = as.factor(c("East North Central; k=88; N=6305", "East South Central; k=10; N=1312", "Middle Atlantic; k=73; N=7074", "Mountain; k=15; N=1217", "New England; k=27; N=1906", "Pacific; k=71; N=7953", "South Atlantic; k=83; N=7528", "West North Central; k=13; N=1219", "West South Central; k=14; N=385"))
regions$region_names_with_k_N <- factor(regions$region_names_with_k_N, levels = regions$region_names_with_k_N[order(regions$region_names_with_k_N, decreasing = TRUE)])

ggplot(data=regions, aes(color=region_names)) + 
  geom_pointrange(mapping=aes(x=regions$region_names, y=regions$region_level_effects_controlled_for_baseline, ymin=regions$std_error_yi_lower , ymax=regions$std_error_yi_upper), size=0.5)+
  coord_flip() + 
  geom_hline(yintercept = baseline_intercept_yi, linetype="twodash", color = "blue", size=0.5) +
  theme(legend.position = "none")+
  xlab("Region") +
  ylab("Yi controlled for baseline")+
  ggsave("plots/forest_plot_region_yi_controlled_for_baseline.png", width = 7, height = 4.5)

ggplot(data=regions, aes(color=region_names_with_k)) + 
  geom_pointrange(mapping=aes(x=regions$region_names_with_k, y=regions$region_level_effects_controlled_for_baseline, ymin=regions$std_error_yi_lower , ymax=regions$std_error_yi_upper), size=0.5)+
  coord_flip() + 
  geom_hline(yintercept = baseline_intercept_yi, linetype="twodash", color = "blue", size=0.5) +
  theme(legend.position = "none")+
  xlab("Region") +
  ylab("Yi controlled for baseline")+
  ggsave("plots/forest_plot_region_yi_controlled_for_baseline_with_k.png", width = 8, height = 4.5)

ggplot(data=regions, aes(color=region_names_with_k_N)) + 
  geom_pointrange(mapping=aes(x=regions$region_names_with_k_N, y=regions$region_level_effects_controlled_for_baseline, ymin=regions$std_error_yi_lower , ymax=regions$std_error_yi_upper), size=0.5)+
  coord_flip() + 
  geom_hline(yintercept = baseline_intercept_yi, linetype="twodash", color = "blue", size=0.5) +
  theme(legend.position = "none")+
  xlab("Region") +
  ylab("Yi controlled for baseline")+
  ggsave("plots/forest_plot_region_yi_controlled_for_baseline_with_k_N.png", width = 9, height = 4.5)



### C) Scatter plots: state/region-level indicators against state/region-level P(C)s controlled for baseline

# 1) Retrieve state-level and region level means

state_level_indicator_means = final_data_VIF %>%
  group_by(state)%>%
  summarise(mean_church_exp_West = mean(church_exp_West, na.rm=TRUE),
            mean_church_exp_East = mean(church_exp_East, na.rm=TRUE),
            mean_trust_reversed_zscore = mean(trust_reversed_zscore, na.rm=TRUE),
            mean_percent_farms_operated_by_owner_zscore = mean(percent_farms_operated_by_owner_zscore, na.rm=TRUE),
            mean_historical_subsistence_index = mean(historical_subsistence_index, na.rm=TRUE),
            mean_subsistence_index = mean(subsistence_index, na.rm=TRUE),
            #mean_residential_mob = mean(residential_mob, na.rm=TRUE),
            #mean_tightness_zscore = mean(tightness_zscore, na.rm=TRUE),
            mean_collectivism = mean(collectivism, na.rm=TRUE),
            mean_self_expression = mean(self_expression, na.rm=TRUE),
            mean_PC_religion = mean(PC_religion, na.rm=TRUE),
            #mean_PC1_government_effectiveness = mean(PC1_government_effectiveness, na.rm=TRUE),
            mean_PC2_voter_turnout_and_corruption_control_perception = mean(PC2_voter_turnout_and_corruption_control_perception, na.rm=TRUE),
            mean_PC3_institutional_trust_and_perception_of_electoral_integrity = mean(PC3_institutional_trust_and_perception_of_electoral_integrity, na.rm=TRUE),
            mean_PC4_state_pension_funding_and_corruption_convictions = mean(PC4_state_pension_funding_and_corruption_convictions, na.rm=TRUE),
            mean_PC5_economic_wealth_and_inequality_and_health_insurance = mean(PC5_economic_wealth_and_inequality_and_health_insurance, na.rm=TRUE),
            #mean_PC1_disease_mortality_and_violent_crime = mean(PC1_disease_mortality_and_violent_crime, na.rm=TRUE),
            mean_PC2_property_crime_and_caloric_suitability = mean(PC2_property_crime_and_caloric_suitability, na.rm=TRUE),
            mean_PC3_disaster_freq_cost_and_demanding_geoclimate = mean(PC3_disaster_freq_cost_and_demanding_geoclimate, na.rm=TRUE))

region_level_indicator_means = final_data_VIF %>%
  filter(!is.na(region))%>%
  group_by(region)%>%
  summarise(mean_church_exp_West = mean(church_exp_West, na.rm=TRUE),
            mean_church_exp_East = mean(church_exp_East, na.rm=TRUE),
            mean_trust_reversed_zscore = mean(trust_reversed_zscore, na.rm=TRUE),
            mean_percent_farms_operated_by_owner_zscore = mean(percent_farms_operated_by_owner_zscore, na.rm=TRUE),
            mean_historical_subsistence_index = mean(historical_subsistence_index, na.rm=TRUE),
            mean_subsistence_index = mean(subsistence_index, na.rm=TRUE),
            #mean_residential_mob = mean(residential_mob, na.rm=TRUE),
            #mean_tightness_zscore = mean(tightness_zscore, na.rm=TRUE),
            mean_collectivism = mean(collectivism, na.rm=TRUE),
            mean_self_expression = mean(self_expression, na.rm=TRUE),
            mean_PC_religion = mean(PC_religion, na.rm=TRUE),
            #mean_PC1_government_effectiveness = mean(PC1_government_effectiveness, na.rm=TRUE),
            mean_PC2_voter_turnout_and_corruption_control_perception = mean(PC2_voter_turnout_and_corruption_control_perception, na.rm=TRUE),
            mean_PC3_institutional_trust_and_perception_of_electoral_integrity = mean(PC3_institutional_trust_and_perception_of_electoral_integrity, na.rm=TRUE),
            mean_PC4_state_pension_funding_and_corruption_convictions = mean(PC4_state_pension_funding_and_corruption_convictions, na.rm=TRUE),
            mean_PC5_economic_wealth_and_inequality_and_health_insurance = mean(PC5_economic_wealth_and_inequality_and_health_insurance, na.rm=TRUE),
            #mean_PC1_disease_mortality_and_violent_crime = mean(PC1_disease_mortality_and_violent_crime, na.rm=TRUE),
            mean_PC2_property_crime_and_caloric_suitability = mean(PC2_property_crime_and_caloric_suitability, na.rm=TRUE),
            mean_PC3_disaster_freq_cost_and_demanding_geoclimate = mean(PC3_disaster_freq_cost_and_demanding_geoclimate, na.rm=TRUE))

## Add yis
state_level_indicator_means = left_join(state_level_indicator_means, states, by="state")
region_level_indicator_means = left_join(region_level_indicator_means, regions, by="region")

state_level_indicator_means[state_level_indicator_means=="NaN"] = NA
region_level_indicator_means[region_level_indicator_means=="NaN"] = NA

# STATE: Scatter plot for every indicator: state-level P(C)s controlled for baseline

ggplot(state_level_indicator_means, aes(x=mean_church_exp_West, y=state_level_effects_controlled_for_baseline, color=state)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(state_level_indicator_means, aes(x=mean_church_exp_East, y=state_level_effects_controlled_for_baseline, color=state)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(state_level_indicator_means, aes(x=mean_percent_farms_operated_by_owner_zscore, y=state_level_effects_controlled_for_baseline, color=state)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")
# WEIRD: more percent_farms --> lower C
# In overall model other way around, but in model including only percent_farms: negative effect like in scatter
# ---> So must be some additive effect with the other indicators

ggplot(state_level_indicator_means, aes(x=mean_historical_subsistence_index, y=state_level_effects_controlled_for_baseline, color=state)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(state_level_indicator_means, aes(x=mean_subsistence_index, y=state_level_effects_controlled_for_baseline, color=state)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

# ggplot(state_level_indicator_means, aes(x=mean_residential_mob, y=state_level_effects_controlled_for_baseline, color=state)) + 
#   geom_point() + 
#   geom_smooth(method=lm, aes(colour=TRUE)) + 
#   labs(y = "yi controlled for baseline")
# 
# ggplot(state_level_indicator_means, aes(x=mean_tightness_zscore, y=state_level_effects_controlled_for_baseline, color=state)) + 
#   geom_point() + 
#   geom_smooth(method=lm, aes(colour=TRUE)) + 
#   labs(y = "yi controlled for baseline")

ggplot(state_level_indicator_means, aes(x=mean_collectivism, y=state_level_effects_controlled_for_baseline, color=state)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(state_level_indicator_means, aes(x=mean_self_expression, y=state_level_effects_controlled_for_baseline, color=state)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(state_level_indicator_means, aes(x=mean_PC_religion, y=state_level_effects_controlled_for_baseline, color=state)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(state_level_indicator_means, aes(x=mean_PC1_government_effectiveness, y=state_level_effects_controlled_for_baseline, color=state)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(state_level_indicator_means, aes(x=mean_PC2_voter_turnout_and_corruption_control_perception, y=state_level_effects_controlled_for_baseline, color=state)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(state_level_indicator_means, aes(x=mean_PC3_institutional_trust_and_perception_of_electoral_integrity, y=state_level_effects_controlled_for_baseline, color=state)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(state_level_indicator_means, aes(x=mean_PC4_state_pension_funding_and_corruption_convictions, y=state_level_effects_controlled_for_baseline, color=state)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(state_level_indicator_means, aes(x=mean_PC5_economic_wealth_and_inequality_and_health_insurance, y=state_level_effects_controlled_for_baseline, color=state)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(state_level_indicator_means, aes(x=mean_PC1_disease_mortality_and_violent_crime, y=state_level_effects_controlled_for_baseline, color=state)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(state_level_indicator_means, aes(x=mean_PC2_property_crime_and_caloric_suitability, y=state_level_effects_controlled_for_baseline, color=state)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(state_level_indicator_means, aes(x=mean_PC3_disaster_freq_cost_and_demanding_geoclimate, y=state_level_effects_controlled_for_baseline, color=state)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")



# REGION: Scatter plot for every indicator: region-level P(C)s controlled for baseline

ggplot(region_level_indicator_means, aes(x=mean_church_exp_West, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(region_level_indicator_means, aes(x=mean_church_exp_East, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(region_level_indicator_means, aes(x=mean_percent_farms_operated_by_owner_zscore, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(region_level_indicator_means, aes(x=mean_historical_subsistence_index, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(region_level_indicator_means, aes(x=mean_subsistence_index, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(region_level_indicator_means, aes(x=mean_residential_mob, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(region_level_indicator_means, aes(x=mean_tightness_zscore, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(region_level_indicator_means, aes(x=mean_collectivism, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(region_level_indicator_means, aes(x=mean_self_expression, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(region_level_indicator_means, aes(x=mean_return_change_reversed_zscore, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(region_level_indicator_means, aes(x=mean_PC_religion, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(region_level_indicator_means, aes(x=mean_PC1_government_effectiveness, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(region_level_indicator_means, aes(x=mean_PC2_voter_turnout_and_corruption_control_perception, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(region_level_indicator_means, aes(x=mean_PC3_institutional_trust_and_perception_of_electoral_integrity, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(region_level_indicator_means, aes(x=mean_PC4_state_pension_funding_and_corruption_convictions, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(region_level_indicator_means, aes(x=mean_PC5_economic_wealth_and_inequality_and_health_insurance, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(region_level_indicator_means, aes(x=mean_PC1_disease_mortality_and_violent_crime, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(region_level_indicator_means, aes(x=mean_PC2_property_crime_and_caloric_suitability, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")

ggplot(region_level_indicator_means, aes(x=mean_PC3_disaster_freq_cost_and_demanding_geoclimate, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(colour=TRUE)) + 
  labs(y = "yi controlled for baseline")



####################################################################################################

### FOREST PLOTS WITH AND WITHOUT DISCUSSION
### Flip discussion coding, keeping other baseline moderators constant, on state and region level

#####################################################################################################

## Flip coding
final_data_VIF$discussion_flipped = ifelse(final_data_VIF$discussion==1, 0, ifelse(final_data_VIF$discussion==0, 1, 0.5))
final_data_VIF$discussion_flipped = as.factor(final_data_VIF$discussion_flipped)

#### Imputation with MICE package 

### 1) Data with PCA: 

imp <- mice(final_data_VIF, m=2, maxit=0, seed=11)
pred <- imp$predictorMatrix #get the predictor matrix

#Don't use some variables for imputation:
pred[,"study_ID"] <- 0 # don't use N for imputing
pred[,"N"] <- 0 # don't use N for imputing
pred[,"P_C"] <- 0 # don't use P_C for imputing
pred[,"state"] <- 0 # don't use state for imputing
pred[,"region"] <- 0 # don't use region for imputing
pred[,"yi"] <- 0 # don't use yi for imputing
pred[,"vi"] <- 0 # don't use vi for imputing
pred[,"discussion"] <- 0 # don't use vi for imputing

#generate the multiple imputations
final_data_VIF_imputed = mice(final_data_VIF, m = 5, maxit = 5, predictorMatrix=pred, seed=11)

##############################################################################################################################

####### STATE AND REGION-LEVEL MODELS FOR FOREST PLOTS - FLIPPED DISCUSSION

# BUT keep in mind: minority of studies with discussion
table(final_data_VIF$discussion)

### 1) state as moderator + baseline

fit_state <- with(final_data_VIF_imputed, rma.mv(yi, vi, mods= ~ state+symmetry+one_shot+matching+group_size+if_PD_K_index+sanction+simultaneous+discussion_flipped, method="REML",
                                               random = ~ 1 | study_ID))
pool_state <- pool(fit_state)
round(summary(pool_state, conf.int=TRUE), 4)


### 2) region as moderator + baseline

fit_region <- with(final_data_VIF_imputed, rma.mv(yi, vi, mods= ~ region+symmetry+one_shot+matching+group_size+if_PD_K_index+sanction+simultaneous+discussion_flipped, method="REML",
                                               random = ~ 1 | study_ID))
pool_region <- pool(fit_region)
round(summary(pool_region, conf.int=TRUE), 4)


### 3) Only baseline

fit_baseline <- with(final_data_VIF_imputed, rma.mv(yi, vi, mods= ~ symmetry+one_shot+matching+group_size+if_PD_K_index+sanction+simultaneous+discussion_flipped, method="REML",
                                                random = ~ 1 | study_ID))
pool_baseline <- pool(fit_baseline)
round(summary(pool_baseline, conf.int=TRUE), 4)

#########################################################################################

### Extract effect sizes - FLIPPED DISCUSSION

# state-level
state_intercept = summary(pool_state)$estimate[1]

states_discussion_baseline1 <- data.frame(state_level_effects_controlled_for_baseline=rep.int(NA, 38),
                        state_level_P_C_controlled_for_baseline = rep.int(NA, 38),
                        state=rep.int(NA, 38),
                        estimate_conf_lower = rep.int(NA, 38),
                        estimate_conf_upper = rep.int(NA, 38)) 

states_discussion_baseline1$state_level_effects_controlled_for_baseline_uncorrected_for_intercept = summary(pool_state)$estimate[1:38]
states_discussion_baseline1$state_level_effects_controlled_for_baseline[2:38] = state_intercept + states_discussion_baseline1$state_level_effects_controlled_for_baseline_uncorrected_for_intercept[2:38]
states_discussion_baseline1$state_level_effects_controlled_for_baseline[1] = state_intercept # first state IS baseline, so = INTERCEPT
states_discussion_baseline1$state_level_P_C_controlled_for_baseline = transf.ilogit.int(states_discussion_baseline1$state_level_effects_controlled_for_baseline)
states_discussion_baseline1$state = c("Alabama", "Arizona", "Arkansas", "California", "Connecticut", "Delaware", "District of Columbia", "Florida", "Georgia", "Hawaii", "Illinois", "Indiana", "Iowa", "Kansas", "Maryland", "Massachusetts", "Michigan", "Minnesota", "Mississippi", "Missouri", "Nevada", "New Jersey", "New York", "North Carolina", "Ohio", "Oklahoma", "Oregon", "Pennsylvania", "Rhode Island", "South Carolina", "Tennessee", "Texas", "Utah", "Vermont", "Virginia", "Washington", "Wisconsin", "Wyoming")
states_discussion_baseline1$std_error = summary(pool_state)$std.error[1:38]
# Add std.error estimates for P(C)
states_discussion_baseline1$estimate_conf_lower = states_discussion_baseline1$state_level_effects_controlled_for_baseline - states_discussion_baseline1$std_error
states_discussion_baseline1$estimate_conf_upper = states_discussion_baseline1$state_level_effects_controlled_for_baseline + states_discussion_baseline1$std_error
states_discussion_baseline1$estimate_conf_lower_P_C = transf.ilogit.int(states_discussion_baseline1$estimate_conf_lower)
states_discussion_baseline1$estimate_conf_upper_P_C = transf.ilogit.int(states_discussion_baseline1$estimate_conf_upper)

# region-level
region_intercept = summary(pool_region)$estimate[1]

regions_discussion_baseline1 <- data.frame(region_level_effects_controlled_for_baseline=rep.int(NA, 9),
                       region_level_P_C_controlled_for_baseline = rep.int(NA, 9),
                       region=rep.int(NA, 9)) 

regions_discussion_baseline1$region_level_effects_controlled_for_baseline_uncorrected_for_intercept = summary(pool_region)$estimate[1:9]
regions_discussion_baseline1$region_level_effects_controlled_for_baseline[2:9] = region_intercept + regions_discussion_baseline1$region_level_effects_controlled_for_baseline_uncorrected_for_intercept[2:9]
regions_discussion_baseline1$region_level_effects_controlled_for_baseline[1] = region_intercept # first region IS baseline, so = INTERCEPT
regions_discussion_baseline1$region_level_P_C_controlled_for_baseline = transf.ilogit.int(regions_discussion_baseline1$region_level_effects_controlled_for_baseline)
regions_discussion_baseline1$region = c("East North Central", "East South Central", "Middle Atlantic", "Mountain", "New England", "Pacific", "South Atlantic", "West North Central", "West South Central")
regions_discussion_baseline1$std_error = summary(pool_region)$std.error[1:9]
# Add std.error estimates for P(C)
regions_discussion_baseline1$estimate_conf_lower = regions_discussion_baseline1$region_level_effects_controlled_for_baseline - regions_discussion_baseline1$std_error
regions_discussion_baseline1$estimate_conf_upper = regions_discussion_baseline1$region_level_effects_controlled_for_baseline + regions_discussion_baseline1$std_error
regions_discussion_baseline1$estimate_conf_lower_P_C = transf.ilogit.int(regions_discussion_baseline1$estimate_conf_lower)
regions_discussion_baseline1$estimate_conf_upper_P_C = transf.ilogit.int(regions_discussion_baseline1$estimate_conf_upper)

# Get overall cooperation rate controlled for baseline
baseline_intercept_yi_discussion1 = summary(pool_baseline)$estimate[1]
baseline_intercept_P_C = transf.ilogit.int(baseline_intercept_yi)


###########################################################################################################

##### PLOTTING: Forest plots - FLIPPED DISCUSSION

# 1) yis: state-level

# Calculate +- std.error
states_discussion_baseline1$std_error_yi_upper = states_discussion_baseline1$state_level_effects_controlled_for_baseline+ states_discussion_baseline1$std_error
states_discussion_baseline1$std_error_yi_lower = states_discussion_baseline1$state_level_effects_controlled_for_baseline- states_discussion_baseline1$std_error

# Set factor ordering to preserve alphabetical state order despite flipping
states_discussion_baseline1$state_name = as.factor(c("Alabama", "Arizona", "Arkansas", "California", "Connecticut", "Delaware", "District of Columbia", "Florida", "Georgia", "Hawaii", "Illinois", "Indiana", "Iowa", "Kansas", "Maryland", "Massachusetts", "Michigan", "Minnesota", "Mississippi", "Missouri", "Nevada", "New Jersey", "New York", "North Carolina", "Ohio", "Oklahoma", "Oregon", "Pennsylvania", "Rhode Island", "South Carolina", "Tennessee", "Texas", "Utah", "Vermont", "Virginia", "Washington", "Wisconsin", "Wyoming"))
states_discussion_baseline1$state_name <- factor(states_discussion_baseline1$state_name, levels = states_discussion_baseline1$state_name[order(states_discussion_baseline1$state_name, decreasing = TRUE)])
states_discussion_baseline1$state_name_with_k = as.factor(c("Alabama; k=6", "Arizona; k=8", "Arkansas; k=1", "California; k=50", "Connecticut; k=3", "Delaware; k=6", "District of Columbia; k=1", "Florida; k=20", "Georgia; k=6", "Hawaii; k=4", "Illinois; k=40", "Indiana; k=1", "Iowa; k=4", "Kansas; k=3", "Maryland; k=1", "Massachusetts; k=11", "Michigan; k=24", "Minnesota; k=1", "Mississippi; k=3", "Missouri; k=5", "Nevada; k=1", "New Jersey; k=12", "New York; k=34", "North Carolina; k=35", "Ohio; k=22", "Oklahoma; k=1", "Oregon; k=5", "Pennsylvania; k=27", "Rhode Island; k=11", "South Carolina; k=12", "Tennessee; k=1", "Texas; k=12", "Utah; k=3", "Vermont; k=2", "Virginia; k=2", "Washington; k=12", "Wisconsin; k=1", "Wyoming; k=3"))
states_discussion_baseline1$state_name_with_k <- factor(states_discussion_baseline1$state_name_with_k, levels = states_discussion_baseline1$state_name_with_k[order(states_discussion_baseline1$state_name_with_k, decreasing = TRUE)])
states_discussion_baseline1$state_name_with_k_N = as.factor(c("Alabama; k=6; N=1102", "Arizona; k=8; N=587", "Arkansas; k=1; N=31", "California; k=50; N=5772", "Connecticut; k=3; N=230", "Delaware; k=6; N=1022", "District of Columbia; k=1; N=70", "Florida; k=20; N=1535", "Georgia; k=6; N=865", "Hawaii; k=4; N=462", "Illinois; k=40; N=3001", "Indiana; k=1; N=100", "Iowa; k=4; N=180", "Kansas; k=3; N=184", "Maryland; k=1; N=48", "Massachusetts; k=11; N=987", "Michigan; k=24; N=1867", "Minnesota; k=1; N=200", "Mississippi; k=3; N=120", "Missouri; k=5; N=655", "Nevada; k=1; N=96", "New Jersey; k=12; N=482", "New York; k=34; N=4351", "North Carolina; k=35; N=2579", "Ohio; k=22; N=1209", "Oklahoma; k=1; N=18", "Oregon; k=5; N=536", "Pennsylvania; k=27; N=2241", "Rhode Island; k=11; N=633", "South Carolina; k=12; N=1150", "Tennessee; k=1; N=90", "Texas; k=12; N=336", "Utah; k=3; N=271", "Vermont; k=2; N=56", "Virginia; k=2; N=259", "Washington; k=12; N=1183", "Wisconsin; k=1; N=128", "Wyoming; k=3; N=263"))
states_discussion_baseline1$state_name_with_k_N <- factor(states_discussion_baseline1$state_name_with_k_N, levels = states_discussion_baseline1$state_name_with_k_N[order(states_discussion_baseline1$state_name_with_k_N, decreasing = TRUE)])

ggplot(data=states_discussion_baseline1, aes(color=state_name)) + 
  geom_pointrange(mapping=aes(x=states_discussion_baseline1$state_name, y=states_discussion_baseline1$state_level_effects_controlled_for_baseline, ymin=states_discussion_baseline1$std_error_yi_lower , ymax=states_discussion_baseline1$std_error_yi_upper), size=0.5)+
  coord_flip() + 
  geom_hline(yintercept = baseline_intercept_yi_discussion1, linetype="twodash", color = "blue", size=0.5) +
  theme(legend.position = "none")+
  xlab("State") +
  ylab("Yi controlled for baseline (discussion = 1)")+
  ggsave("plots/forest_plot_state_yi_controlled_for_baseline_discussion1.png", width = 8, height = 6)

ggplot(data=states_discussion_baseline1, aes(color=state_name_with_k)) + 
  geom_pointrange(mapping=aes(x=states_discussion_baseline1$state_name_with_k, y=states_discussion_baseline1$state_level_effects_controlled_for_baseline, ymin=states_discussion_baseline1$std_error_yi_lower , ymax=states_discussion_baseline1$std_error_yi_upper), size=0.5)+
  coord_flip() + 
  geom_hline(yintercept = baseline_intercept_yi_discussion1, linetype="twodash", color = "blue", size=0.5) +
  theme(legend.position = "none")+
  xlab("State") +
  ylab("Yi controlled for baseline (discussion = 1)")+
  ggsave("plots/forest_plot_state_yi_controlled_for_baseline_discussion1_with_k.png", width = 9, height = 6)

ggplot(data=states_discussion_baseline1, aes(color=state_name_with_k_N)) + 
  geom_pointrange(mapping=aes(x=states_discussion_baseline1$state_name_with_k_N, y=states_discussion_baseline1$state_level_effects_controlled_for_baseline, ymin=states_discussion_baseline1$std_error_yi_lower , ymax=states_discussion_baseline1$std_error_yi_upper), size=0.5)+
  coord_flip() + 
  geom_hline(yintercept = baseline_intercept_yi_discussion1, linetype="twodash", color = "blue", size=0.5) +
  theme(legend.position = "none")+
  xlab("State") +
  ylab("Yi controlled for baseline (discussion = 1)")+
  ggsave("plots/forest_plot_state_yi_controlled_for_baseline_discussion1_with_k_N.png", width = 10, height = 6)


# 2) yis: region-level

regions_discussion_baseline1$std_error_yi_upper = regions_discussion_baseline1$region_level_effects_controlled_for_baseline+ regions_discussion_baseline1$std_error
regions_discussion_baseline1$std_error_yi_lower = regions_discussion_baseline1$region_level_effects_controlled_for_baseline- regions_discussion_baseline1$std_error

# Use better names, convert to factor and make factor ordering preserve alphabetical state order despite flipping
regions_discussion_baseline1$region_names = as.factor(c("East North Central", "East South Central", "Middle Atlantic", "Mountain", "New England", "Pacific", "South Atlantic", "West North Central", "West South Central"))
regions_discussion_baseline1$region_names <- factor(regions_discussion_baseline1$region_names, levels = regions_discussion_baseline1$region_names[order(regions_discussion_baseline1$region_names, decreasing = TRUE)])
regions_discussion_baseline1$region_names_with_k = as.factor(c("East North Central; k=88", "East South Central; k=10", "Middle Atlantic; k=73", "Mountain; k=15", "New England; k=27", "Pacific; k=71", "South Atlantic; k=83", "West North Central; k=13", "West South Central; k=14"))
regions_discussion_baseline1$region_names_with_k <- factor(regions_discussion_baseline1$region_names_with_k, levels = regions_discussion_baseline1$region_names_with_k[order(regions_discussion_baseline1$region_names_with_k, decreasing = TRUE)])
regions_discussion_baseline1$region_names_with_k_N = as.factor(c("East North Central; k=88; N=6305", "East South Central; k=10; N=1312", "Middle Atlantic; k=73; N=7074", "Mountain; k=15; N=1217", "New England; k=27; N=1906", "Pacific; k=71; N=7953", "South Atlantic; k=83; N=7528", "West North Central; k=13; N=1219", "West South Central; k=14; N=385"))
regions_discussion_baseline1$region_names_with_k_N <- factor(regions_discussion_baseline1$region_names_with_k_N, levels = regions_discussion_baseline1$region_names_with_k_N[order(regions_discussion_baseline1$region_names_with_k_N, decreasing = TRUE)])

ggplot(data=regions_discussion_baseline1, aes(color=region_names)) + 
  geom_pointrange(mapping=aes(x=regions_discussion_baseline1$region_names, y=regions_discussion_baseline1$region_level_effects_controlled_for_baseline, ymin=regions_discussion_baseline1$std_error_yi_lower , ymax=regions_discussion_baseline1$std_error_yi_upper), size=0.5)+
  coord_flip() + 
  geom_hline(yintercept = baseline_intercept_yi_discussion1, linetype="twodash", color = "blue", size=0.5) +
  theme(legend.position = "none")+
  xlab("Region") +
  ylab("Yi controlled for baseline (discussion = 1)")+
  ggsave("plots/forest_plot_region_yi_controlled_for_baseline_discussion1.png", width = 7, height = 4.5)

ggplot(data=regions_discussion_baseline1, aes(color=region_names_with_k)) + 
  geom_pointrange(mapping=aes(x=regions_discussion_baseline1$region_names_with_k, y=regions_discussion_baseline1$region_level_effects_controlled_for_baseline, ymin=regions_discussion_baseline1$std_error_yi_lower , ymax=regions_discussion_baseline1$std_error_yi_upper), size=0.5)+
  coord_flip() + 
  geom_hline(yintercept = baseline_intercept_yi_discussion1, linetype="twodash", color = "blue", size=0.5) +
  theme(legend.position = "none")+
  xlab("Region") +
  ylab("Yi controlled for baseline (discussion = 1)")+
  ggsave("plots/forest_plot_region_yi_controlled_for_baseline_discussion1_with_k.png", width = 8, height = 4.5)

ggplot(data=regions_discussion_baseline1, aes(color=region_names_with_k_N)) + 
  geom_pointrange(mapping=aes(x=regions_discussion_baseline1$region_names_with_k_N, y=regions_discussion_baseline1$region_level_effects_controlled_for_baseline, ymin=regions_discussion_baseline1$std_error_yi_lower , ymax=regions_discussion_baseline1$std_error_yi_upper), size=0.5)+
  coord_flip() + 
  geom_hline(yintercept = baseline_intercept_yi_discussion1, linetype="twodash", color = "blue", size=0.5) +
  theme(legend.position = "none")+
  xlab("Region") +
  ylab("Yi controlled for baseline (discussion = 1)")+
  ggsave("plots/forest_plot_region_yi_controlled_for_baseline_discussion1_with_k_N.png", width = 9, height = 4.5)


#### Combine in 1 Forest plot (with k)

# Merge data from discussion=0 and santion=1
states_discussion_merged = bind_rows(states, states_discussion_baseline1)
states_discussion_merged$discussion = c(rep("no discussion",38), rep("discussion",38))
regions_discussion_merged = bind_rows(regions, regions_discussion_baseline1)
regions_discussion_merged$discussion = c(rep("no discussion",9), rep("discussion",9))

# 1) state
ggplot(data=states_discussion_merged, aes(x = state_name, y = state_level_effects_controlled_for_baseline, ymin = std_error_yi_lower, ymax = std_error_yi_upper ))+
  geom_pointrange(aes(col=discussion))+
  geom_hline(yintercept = baseline_intercept_yi_discussion1, linetype="twodash", color = "#F8766D", size=0.5) +
  geom_hline(yintercept = baseline_intercept_yi, linetype="twodash", color = "#00BFC4", size=0.5) +
  xlab('Discussion')+ ylab("Yi controlled for baseline")+
  geom_errorbar(aes(ymin=std_error_yi_lower, ymax=std_error_yi_upper,col=discussion),width=0.5,cex=1)+ 
  facet_wrap(~discussion,strip.position="left",nrow=9,scales = "free_y") +
  theme(legend.position = "none")+
  coord_flip()+
  ggsave("plots/forest_plot_state_yi_controlled_for_baseline_discussion_facet.png", width = 9, height = 7)

# 2) region

ggplot(data=regions_discussion_merged, aes(x = region_names, y = region_level_effects_controlled_for_baseline, ymin = std_error_yi_lower, ymax = std_error_yi_upper ))+
  geom_pointrange(aes(col=discussion))+
  xlab('Discussion')+ ylab("Yi controlled for baseline")+
  geom_errorbar(aes(ymin=std_error_yi_lower, ymax=std_error_yi_upper,col=discussion),width=0.5,cex=1)+ 
  facet_wrap(~discussion,strip.position="left",nrow=9,scales = "free_y") +
  theme(legend.position = "none")+
  geom_hline(yintercept = baseline_intercept_yi_discussion1, linetype="twodash", color = "#F8766D", size=0.5) +
  geom_hline(yintercept = baseline_intercept_yi, linetype="twodash", color = "#00BFC4", size=0.5) +
  coord_flip()+
  ggsave("plots/forest_plot_region_yi_controlled_for_baseline_discussion_facet.png", width = 9, height = 5)


############################################################################################################################################################

##### Fixed effects summary

### RMA.MV (from THEORY-LEVEL models)
# More collectivism - more C
# More Religiosity, more C
# Less self_expression - more C
# Less church exposure West - more C
# Lower institutional trust and lower perception of electoral integrity (but higher trust in army) - more C
# Less trust - more C

### RMA - marg sign. (from THEORY-LEVEL models)
# Higher physical threat, more C

# ### LM:
# # More collectivism, more C
# # higher physical threat, more C
# # less church west, more C
# # less historical herding, more C
# # Less self-expression, more C
# # Better roads but worse schools, more C
# # Lower voter turnout and corruption control perception - more C
# 
# ### LM - marg. sign.
# # more tightness, less C
# # More religiosity, more C
# # higher disease mortality and violent crime rates - higher C

### Plot most prominent trends: yi (corrected for baseline) and ...
# collectivism, religiosity, self_expression, church exposure West, institutional trust , trust, physical threat

### 1) State-level

# yi and collectivism
cor(state_level_indicator_means$state_level_effects_controlled_for_baseline, state_level_indicator_means$mean_collectivism, use="pairwise.complete.obs")
# R = 0.3625831
cor.test(state_level_indicator_means$state_level_effects_controlled_for_baseline, state_level_indicator_means$mean_collectivism, na.action=na.exclude)
# p-value = 0.02742
ggplot(state_level_indicator_means, aes(x=mean_collectivism, y=state_level_effects_controlled_for_baseline, color=state)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(group=1)) + 
  geom_text_repel(aes(label=state,color=state), segment.size = 0.1)+
  annotate("text", x = -1, y = 0.35, label = "r = 0.36; p < 0.05")+
  labs(x="Collectivism",y = "yi controlled for baseline")+
  theme(legend.position = "none") +
  ggsave("plots/corr_yi_collectivism.png", width=8, height=6)

# yi and collectivism (without outlier Hawaii)
temp = filter(state_level_indicator_means, !state=="Hawaii")
cor(temp$state_level_effects_controlled_for_baseline, temp$mean_collectivism, use="pairwise.complete.obs")
# R = 0.3993859
cor.test(temp$state_level_effects_controlled_for_baseline, temp$mean_collectivism, na.action=na.exclude)
# p-value = 0.01582
ggplot(temp, aes(x=mean_collectivism, y=state_level_effects_controlled_for_baseline, color=state)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(group=1)) + 
  geom_text_repel(aes(label=state,color=state), segment.size = 0.1)+
  annotate("text", x = -1.15, y = 0.2, label = "r = 0.40; p < 0.05")+
  labs(x="Collectivism",y = "yi controlled for baseline")+
  theme(legend.position = "none") +
  ggsave("plots/corr_yi_collectivism_no_Hawaii.png", width=8, height=6)

# yi and religiosity
cor(state_level_indicator_means$state_level_effects_controlled_for_baseline, state_level_indicator_means$mean_PC_religion, use="pairwise.complete.obs")
# R = 0.4056147
cor.test(state_level_indicator_means$state_level_effects_controlled_for_baseline, state_level_indicator_means$mean_PC_religion, na.action=na.exclude)
# p-value =  0.01275
ggplot(state_level_indicator_means, aes(x=mean_PC_religion, y=state_level_effects_controlled_for_baseline, color=state)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(group=1)) + 
  geom_text_repel(aes(label=state,color=state), segment.size = 0.1)+
  annotate("text", x = -2.2, y = 0.35, label = "r = 0.41; p < 0.05")+
  labs(x="Religiosity",y = "yi controlled for baseline")+
  theme(legend.position = "none") +
  ggsave("plots/corr_yi_religion.png", width=8, height=6)

# yi and self_expression
cor(state_level_indicator_means$state_level_effects_controlled_for_baseline, state_level_indicator_means$mean_self_expression, use="pairwise.complete.obs")
# R = -0.3802613
cor.test(state_level_indicator_means$state_level_effects_controlled_for_baseline, state_level_indicator_means$mean_self_expression, na.action=na.exclude)
# p-value =  0.02025
ggplot(state_level_indicator_means, aes(x=mean_self_expression, y=state_level_effects_controlled_for_baseline, color=state)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(group=1)) + 
  geom_text_repel(aes(label=state,color=state), segment.size = 0.1)+
  annotate("text", x = -0.65, y = 0.55, label = "r = -0.38; p < 0.05")+
  labs(x="Self expression",y = "yi controlled for baseline")+
  theme(legend.position = "none") +
  ggsave("plots/corr_yi_self_expression.png", width=8, height=6)

# yi and mean_church_exp_West
cor(state_level_indicator_means$state_level_effects_controlled_for_baseline, state_level_indicator_means$mean_church_exp_West, use="pairwise.complete.obs")
# R = -0.3576735
cor.test(state_level_indicator_means$state_level_effects_controlled_for_baseline, state_level_indicator_means$mean_church_exp_West, na.action=na.exclude)
# p-value =  0.02975
ggplot(state_level_indicator_means, aes(x=mean_church_exp_West, y=state_level_effects_controlled_for_baseline, color=state)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(group=1)) + 
  geom_text_repel(aes(label=state,color=state), segment.size = 0.1)+
  annotate("text", x = -2.6, y = 0.45, label = "r = -0.36; p < 0.05")+
  labs(x="Exposure to the Western Church",y = "yi controlled for baseline")+
  theme(legend.position = "none") +
  ggsave("plots/corr_yi_church_exp_West.png", width=8, height=6)

# yi and mean_church_exp_West (without outlier Hawaii)
temp = filter(state_level_indicator_means, !state=="Hawaii")
cor(temp$state_level_effects_controlled_for_baseline, temp$mean_church_exp_West, use="pairwise.complete.obs")
# R = -0.3907198
cor.test(temp$state_level_effects_controlled_for_baseline, temp$mean_church_exp_West, na.action=na.exclude)
# p-value =  0.01847
ggplot(temp, aes(x=mean_church_exp_West, y=state_level_effects_controlled_for_baseline, color=state)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(group=1)) + 
  geom_text_repel(aes(label=state,color=state), segment.size = 0.1)+
  annotate("text", x = -1.5, y = 0.4, label = "r = -0.39; p < 0.05")+
  labs(x="Exposure to the Western Church",y = "yi controlled for baseline")+
  theme(legend.position = "none") +
  ggsave("plots/corr_yi_church_exp_West_no_Hawaii.png", width=8, height=6)

# yi and mean_PC3_institutional_trust_and_perception_of_electoral_integrity
cor(state_level_indicator_means$state_level_effects_controlled_for_baseline, state_level_indicator_means$mean_PC3_institutional_trust_and_perception_of_electoral_integrity, use="pairwise.complete.obs")
# R = -0.2978435
cor.test(state_level_indicator_means$state_level_effects_controlled_for_baseline, state_level_indicator_means$mean_PC3_institutional_trust_and_perception_of_electoral_integrity, na.action=na.exclude)
# p-value =  0.08711
ggplot(state_level_indicator_means, aes(x=mean_PC3_institutional_trust_and_perception_of_electoral_integrity, y=state_level_effects_controlled_for_baseline, color=state)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(group=1)) + 
  geom_text_repel(aes(label=state,color=state), segment.size = 0.1)+
  annotate("text", x = -0.7, y = 0.4, label = "r = -0.30; p = 0.09")+
  labs(x="Institutional trust",y = "yi controlled for baseline")+
  theme(legend.position = "none") +
  ggsave("plots/corr_yi_PC3_institutional_trust.png", width=8, height=6)

# yi and mean_trust_reversed_zscore
cor(state_level_indicator_means$state_level_effects_controlled_for_baseline, state_level_indicator_means$mean_trust_reversed_zscore, use="pairwise.complete.obs")
# R = -0.3331883
cor.test(state_level_indicator_means$state_level_effects_controlled_for_baseline, state_level_indicator_means$mean_trust_reversed_zscore, na.action=na.exclude)
# p-value =  0.04389
ggplot(state_level_indicator_means, aes(x=mean_trust_reversed_zscore, y=state_level_effects_controlled_for_baseline, color=state)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(group=1)) + 
  geom_text_repel(aes(label=state,color=state), segment.size = 0.1)+
  annotate("text", x = -1, y = 0.5, label = "r = -0.33; p < 0.05")+
  labs(x="Trust",y = "yi controlled for baseline")+
  theme(legend.position = "none") +
  ggsave("plots/corr_yi_trust.png", width=8, height=6)


### 2) Region-level

# yi and collectivism
cor(region_level_indicator_means$region_level_effects_controlled_for_baseline, region_level_indicator_means$mean_collectivism, use="pairwise.complete.obs")
# R = 0.5878222
cor.test(region_level_indicator_means$region_level_effects_controlled_for_baseline, region_level_indicator_means$mean_collectivism, na.action=na.exclude)
# p-value = 0.09599
ggplot(region_level_indicator_means, aes(x=mean_collectivism, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(group=1)) + 
  geom_text_repel(aes(label=region,color=region), segment.size = 0.1)+
  annotate("text", x = -0.6, y = -0.35, label = "r = 0.59; p = 0.1")+
  labs(x="Collectivism",y = "yi controlled for baseline")+
  theme(legend.position = "none") +
  ggsave("plots/corr_yi_collectivism_region.png", width=8, height=6)

# yi and religiosity
cor(region_level_indicator_means$region_level_effects_controlled_for_baseline, region_level_indicator_means$mean_PC_religion, use="pairwise.complete.obs")
# R = 0.6488197
cor.test(region_level_indicator_means$region_level_effects_controlled_for_baseline, region_level_indicator_means$mean_PC_religion, na.action=na.exclude)
# p-value =  0.05868
ggplot(region_level_indicator_means, aes(x=mean_PC_religion, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(group=1)) + 
  geom_text_repel(aes(label=region,color=region), segment.size = 0.1)+
  annotate("text", x = -2.1, y = -0.3, label = "r = 0.65; p = 0.06")+
  labs(x="Religiosity",y = "yi controlled for baseline")+
  theme(legend.position = "none") +
  ggsave("plots/corr_yi_religion_region.png", width=8, height=6)

# yi and self_expression
cor(region_level_indicator_means$region_level_effects_controlled_for_baseline, region_level_indicator_means$mean_self_expression, use="pairwise.complete.obs")
# R = -0.5987178
cor.test(region_level_indicator_means$region_level_effects_controlled_for_baseline, region_level_indicator_means$mean_self_expression, na.action=na.exclude)
# p-value =  0.08848
ggplot(region_level_indicator_means, aes(x=mean_self_expression, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(group=1)) + 
  geom_text_repel(aes(label=region,color=region), segment.size = 0.1)+
  annotate("text", x = -0.55, y = -0.3, label = "r = -0.60; p = 0.09")+
  labs(x="Self expression",y = "yi controlled for baseline")+
  theme(legend.position = "none") +
  ggsave("plots/corr_yi_self_expression_region.png", width=8, height=6)

# yi and mean_church_exp_West
cor(region_level_indicator_means$region_level_effects_controlled_for_baseline, region_level_indicator_means$mean_church_exp_West, use="pairwise.complete.obs")
# R = -0.410111
cor.test(region_level_indicator_means$region_level_effects_controlled_for_baseline, region_level_indicator_means$mean_church_exp_West, na.action=na.exclude)
# p-value =  0.2729
ggplot(region_level_indicator_means, aes(x=mean_church_exp_West, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(group=1)) + 
  geom_text_repel(aes(label=region,color=region), segment.size = 0.1)+
  annotate("text", x = -1.2, y = -0.3, label = "r = -0.41; p = 0.27")+
  labs(x="Exposure to the Western Church",y = "yi controlled for baseline")+
  theme(legend.position = "none") +
  ggsave("plots/corr_yi_church_exp_West_region.png", width=8, height=6)

# yi and mean_PC3_institutional_trust_and_perception_of_electoral_integrity
cor(region_level_indicator_means$region_level_effects_controlled_for_baseline, region_level_indicator_means$mean_PC3_institutional_trust_and_perception_of_electoral_integrity, use="pairwise.complete.obs")
# R = -0.8764996
cor.test(region_level_indicator_means$region_level_effects_controlled_for_baseline, region_level_indicator_means$mean_PC3_institutional_trust_and_perception_of_electoral_integrity, na.action=na.exclude)
# p-value =  0.001928
ggplot(region_level_indicator_means, aes(x=mean_PC3_institutional_trust_and_perception_of_electoral_integrity, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(group=1)) + 
  geom_text_repel(aes(label=region,color=region), segment.size = 0.1)+
  annotate("text", x = 0.5, y = -0.25, label = "r = -0.88; p < 0.01")+
  labs(x="Institutional trust",y = "yi controlled for baseline")+
  theme(legend.position = "none") +
  ggsave("plots/corr_yi_PC3_institutional_trust_region.png", width=8, height=6)

# yi and mean_trust_reversed_zscore
cor(region_level_indicator_means$region_level_effects_controlled_for_baseline, region_level_indicator_means$mean_trust_reversed_zscore, use="pairwise.complete.obs")
# R = -0.5748478
cor.test(region_level_indicator_means$region_level_effects_controlled_for_baseline, region_level_indicator_means$mean_trust_reversed_zscore, na.action=na.exclude)
# p-value =  0.1054
ggplot(region_level_indicator_means, aes(x=mean_trust_reversed_zscore, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(group=1)) + 
  geom_text_repel(aes(label=region,color=region), segment.size = 0.1)+
  annotate("text", x = -0.7, y = -0.3, label = "r = -0.57; p = 0.11")+
  labs(x="Trust",y = "yi controlled for baseline")+
  theme(legend.position = "none") +
  ggsave("plots/corr_yi_trust_region.png", width=8, height=6)


###########################################################################################################

##### Kruskall-Wallis Test á la Gächter

# Across-regions kruskal:

kruskal.test(yi ~ region, data = final_data_VIF) 
# ---> significant differences between regions!
# Kruskal-Wallis chi-squared = 22.041, df = 8, p-value = 0.004841

# Within-regions Kruskall:

# Make region-data_frames
East_North_Central = filter(final_data_VIF, region == "East North Central")
East_South_Central = filter(final_data_VIF, region == "East South Central")
Middle_Atlantic = filter(final_data_VIF, region == "Middle Atlantic")
Mountain = filter(final_data_VIF, region == "Mountain")
New_England = filter(final_data_VIF, region == "New England")
Pacific = filter(final_data_VIF, region == "Pacific")
South_Atlantic = filter(final_data_VIF, region == "South Atlantic")
West_North_Central = filter(final_data_VIF, region == "West North Central")
West_South_Central = filter(final_data_VIF, region == "West South Central")

# Within-region Kruskall
kruskal.test(yi ~ state, data = East_North_Central)
kruskal.test(yi ~ state, data = East_South_Central)
kruskal.test(yi ~ state, data = Middle_Atlantic)
kruskal.test(yi ~ state, data = Mountain)
kruskal.test(yi ~ state, data = New_England)
kruskal.test(yi ~ state, data = Pacific)
kruskal.test(yi ~ state, data = South_Atlantic)
kruskal.test(yi ~ state, data = West_North_Central)
kruskal.test(yi ~ state, data = West_South_Central)

New_England_no_RI = filter(New_England, !state=="Rhode Island")
kruskal.test(yi ~ state, data = New_England_no_RI)
# chi-squared = 1.1086, df = 2, p-value = 0.5745

# Only within-region of New England significant (South_Atlantic marginally significant: p=0.105)
# ---> mostly non-significant differences within regions!

# Check which pairs of regions are different: 
# pairwise.wilcox.test(final_data_VIF$yi, final_data_VIF$region, p.adjust.method = "BH")


###########################################################################################################################

### Effect size plotting: Full model

# Extract effects, stds, and labels (and reverse, s.t. coord_flip() works)
labels_fit_full_VIF_rev = c("PC Natural disaster frequency and cost", "PC Property crime and caloric suitability", "PC Economic wealth and inequality; health insurance", "PC Ratio of state pension funded and corruption convictions", "PC Institutional trust and perception of electoral integrity", "PC Voter turnout and corruption control perception", "PC Cousin marriage", "PC Religiosity", "Historical subsistence index", "Collectivism", "Percent farms operated by owner", "Eastern church exposure", "Western church exposure", "Self expression", "Subsistence index", "Trust", "Sanction = 1", "Sanction = 0.5", "Simultaneous = 1", "Simultaneous = 0.5", "Discussion = 1", "Discussion = 0.5", "K index", "Group size", "Matching = 1", "Matching = 0.5", "One shot = 1", "One shot = 0.5", "Symmetry = 1", "Symmetry = 0.5", "Intercept")
effects_fit_full_VIF_rev = rev(round(summary(pool(fit_full_VIF))$estimate, 4))
std_error_fit_full_VIF_rev = rev(round(summary(pool(fit_full_VIF))$std.error, 4))

# Make data frame
effects = data.frame(labels_fit_full_VIF_rev, effects_fit_full_VIF_rev, std_error_fit_full_VIF_rev)

### Add factor levels excplicity for plotting in this order
effects$labels_fit_full_VIF_rev <- factor(effects$labels_fit_full_VIF_rev, levels = c("PC Natural disaster frequency and cost", "PC Property crime and caloric suitability", "PC Economic wealth and inequality; health insurance", "PC Ratio of state pension funded and corruption convictions", "PC Institutional trust and perception of electoral integrity", "PC Voter turnout and corruption control perception", "PC Cousin marriage", "PC Religiosity", "Eastern church exposure", "Western church exposure", "Historical subsistence index", "Percent farms operated by owner", "Subsistence index", "Trust", "Self expression", "Collectivism", "Sanction = 1", "Sanction = 0.5", "Simultaneous = 1", "Simultaneous = 0.5", "Discussion = 1", "Discussion = 0.5", "K index", "Group size", "Matching = 1", "Matching = 0.5", "One shot = 1", "One shot = 0.5", "Symmetry = 1", "Symmetry = 0.5", "Intercept"))

ggplot(effects, aes(labels_fit_full_VIF_rev, effects_fit_full_VIF_rev)) +
  geom_bar(stat = "identity", aes(fill = labels_fit_full_VIF_rev)) + 
  geom_errorbar( aes(ymin=effects_fit_full_VIF_rev-std_error_fit_full_VIF_rev, ymax=effects_fit_full_VIF_rev+std_error_fit_full_VIF_rev), alpha=0.5, size=0.3) +
  guides(fill=FALSE) +
  labs(x = "", y = "beta estimates") +
  coord_flip()+
  ggsave("plots/effects_fit_full_with_error_bars_reordered.png", width=8, height=6)

#### Without intercept

# Extract effects, stds, and labels (and reverse, s.t. coord_flip() works)
labels_fit_full_VIF_rev = c("PC Natural disaster frequency and cost", "PC Property crime and caloric suitability", "PC Economic wealth and inequality; health insurance", "PC Ratio of state pension funded and corruption convictions", "PC Institutional trust and perception of electoral integrity", "PC Voter turnout and corruption control perception", "PC Cousin marriage", "PC Religiosity", "Historical subsistence index", "Collectivism", "Percent farms operated by owner", "Eastern church exposure", "Western church exposure", "Self expression", "Subsistence index", "Trust", "Sanction = 1", "Sanction = 0.5", "Simultaneous = 1", "Simultaneous = 0.5", "Discussion = 1", "Discussion = 0.5", "K index", "Group size", "Matching = 1", "Matching = 0.5", "One shot = 1", "One shot = 0.5", "Symmetry = 1", "Symmetry = 0.5")
effects_fit_full_VIF_rev = rev(round(summary(pool(fit_full_VIF))$estimate[-1], 4))
std_error_fit_full_VIF_rev = rev(round(summary(pool(fit_full_VIF))$std.error[-1], 4))

# Make data frame
effects = data.frame(labels_fit_full_VIF_rev, effects_fit_full_VIF_rev, std_error_fit_full_VIF_rev)

### Add factor levels excplicity for plotting in this order
#effects$labels_fit_full_VIF_rev <- factor(effects$labels_fit_full_VIF_rev, levels = c("PC Natural disaster frequency and cost", "PC Property crime and caloric suitability", "PC Economic wealth and inequality; health insurance", "PC Ratio of state pension funded and corruption convictions", "PC Institutional trust and perception of electoral integrity", "PC Voter turnout and corruption control perception", "PC Cousin marriage", "PC Religiosity", "Historical subsistence index", "Collectivism", "Percent farms operated by owner", "Eastern church exposure", "Western church exposure", "Self expression", "Subsistence index", "Trust", "Sanction = 1", "Sanction = 0.5", "Simultaneous = 1", "Simultaneous = 0.5", "Discussion = 1", "Discussion = 0.5", "K index", "Group size", "Matching = 1", "Matching = 0.5", "One shot = 1", "One shot = 0.5", "Symmetry = 1", "Symmetry = 0.5"))
effects$labels_fit_full_VIF_rev <- factor(effects$labels_fit_full_VIF_rev, levels = c("PC Natural disaster frequency and cost", "PC Property crime and caloric suitability", "PC Economic wealth and inequality; health insurance", "PC Ratio of state pension funded and corruption convictions", "PC Institutional trust and perception of electoral integrity", "PC Voter turnout and corruption control perception", "PC Cousin marriage", "PC Religiosity", "Eastern church exposure", "Western church exposure", "Historical subsistence index", "Percent farms operated by owner", "Subsistence index", "Trust", "Self expression", "Collectivism", "Sanction = 1", "Sanction = 0.5", "Simultaneous = 1", "Simultaneous = 0.5", "Discussion = 1", "Discussion = 0.5", "K index", "Group size", "Matching = 1", "Matching = 0.5", "One shot = 1", "One shot = 0.5", "Symmetry = 1", "Symmetry = 0.5"))

ggplot(effects, aes(labels_fit_full_VIF_rev, effects_fit_full_VIF_rev)) +
  geom_bar(stat = "identity", aes(fill = labels_fit_full_VIF_rev)) + 
  geom_errorbar( aes(ymin=effects_fit_full_VIF_rev-std_error_fit_full_VIF_rev, ymax=effects_fit_full_VIF_rev+std_error_fit_full_VIF_rev), alpha=0.5, size=0.3) +
  guides(fill=FALSE) +
  labs(x = "", y = "beta estimates") +
  coord_flip()+
  ggsave("plots/effects_fit_full_with_error_bars_without_intercept_reordered.png", width=8, height=6)


###########################################################################################################

### Cooperation rates states/regions facet plot

# extract region data
region_data = select(final_data_VIF, state, region)%>%
  filter(!duplicated(state))

region_data = left_join(states, region_data, by="state")
region_data$region = as.character(region_data$region)
region_data$region = ifelse(is.na(region_data$region), "South Atlantic", region_data$region)
region_data$state2 = c("Alabama; k=6; N=1102", "Arizona; k=8; N=587", "Arkansas; k=1; N=31", "California; k=50; N=5772", "Connecticut; k=3; N=230", "Delaware; k=6; N=1022", "District of Columbia; k=1; N=70", "Florida; k=20; N=1535", "Georgia; k=6; N=865", "Hawaii; k=4; N=462", "Illinois; k=40; N=3001", "Indiana; k=1; N=100", "Iowa; k=4; N=180", "Kansas; k=3; N=184", "Maryland; k=1; N=48", "Massachusetts; k=11; N=987", "Michigan; k=24; N=1867", "Minnesota; k=1; N=200", "Mississippi; k=3; N=120", "Missouri; k=5; N=655", "Nevada; k=1; N=96", "New Jersey; k=12; N=482", "New York; k=34; N=4351", "North Carolina; k=35; N=2579", "Ohio; k=22; N=1209", "Oklahoma; k=1; N=18", "Oregon; k=5; N=536", "Pennsylvania; k=27; N=2241", "Rhode Island; k=11; N=633", "South Carolina; k=12; N=1150", "Tennessee; k=1; N=90", "Texas; k=12; N=336", "Utah; k=3; N=271", "Vermont; k=2; N=56", "Virginia; k=2; N=259", "Washington; k=12; N=1183", "Wisconsin; k=1; N=128", "Wyoming; k=3; N=263")
region_data = rename(region_data, Region=region)

# yi
ggplot(region_data, aes(x=state, y=state_level_effects_controlled_for_baseline, fill=Region))+
  geom_bar(stat='identity') +
  geom_errorbar( aes(ymin=state_level_effects_controlled_for_baseline-std_error, ymax=state_level_effects_controlled_for_baseline+std_error), alpha=0.5, size=0.3) +
  facet_grid(.~Region, scales="free", space="free_x") +
  #geom_hline(yintercept = baseline_intercept_yi, linetype="twodash", color = "blue", size=0.5) +
  labs(x = "State", y = "State-level yi controlled for baseline") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggsave("plots/effects_region_state_facet_yi.png", width=15, height=5)

# P_C
ggplot(region_data, aes(x=state, y=state_level_P_C_controlled_for_baseline, fill=Region))+
  geom_bar(stat='identity') +
  facet_grid(.~Region, scales="free", space="free_x") +
  geom_hline(yintercept = baseline_intercept_P_C, linetype="twodash", color = "blue", size=0.5) +
  labs(x = "State", y = "State-level P(C)\ncontrolled for baseline") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=14))+
  ggsave("plots/effects_region_state_facet_P_C_bar.png", width=9, height=4)

# P_C, with num info
ggplot(region_data, aes(x=state2, y=state_level_P_C_controlled_for_baseline, fill=Region))+
  geom_bar(stat='identity') +
  facet_grid(.~Region, scales="free", space="free_x") +
  geom_hline(yintercept = baseline_intercept_P_C, linetype="twodash", color = "blue", size=0.5) +
  labs(x = "State", y = "State-level P(C) controlled for baseline") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=14))+
  ggsave("plots/effects_region_state_facet_P_C_bar_k_N.png", width=9, height=5)

# # P_C point
# ggplot(region_data, aes(x=state, y=state_level_effects_controlled_for_baseline, fill=Region))+
#   geom_point() +
#   geom_errorbar( aes(ymin=state_level_effects_controlled_for_baseline-std_error, ymax=state_level_effects_controlled_for_baseline+std_error), alpha=0.5, size=0.3) +
#   facet_grid(.~Region, scales="free", space="free_x") +
#   geom_hline(yintercept = baseline_intercept_yi, linetype="twodash", color = "blue", size=0.5) +
#   labs(x = "State", y = "State-level P(C) controlled for baseline") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   ggsave("plots/effects_region_state_facet_P_C_point.png", width=15, height=8)


### Plot P(C) against effect size yi

ggplot(final_data_VIF, aes(x=P_C, y=yi))+
  geom_point() + 
  labs(x = "proportion of cooperation P(C)", y = "Effect size yi") +
  ggsave("plots/effects_against_proportions.png", width=5, height=5)
