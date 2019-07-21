#####################################################################################################

######## Analysis - External validity
# (only with final_data_VIF)

### Load packages
library("tidyverse")
library("mice")
library("metafor")
library("lme4")
library("faraway")
library("MuMIn")
library("psych")
library("ggrepel")

# Set seed
set.seed(11)

### Set WD
setwd("/home/caroline/Desktop/cross_cultural_metaanalysis/US/analysis")

### Read final_data_VIF
final_data_VIF = read.csv("data_frames/final_data_VIF.csv")


############# EXTERNAL VALIDITY

# External validity measures: 
# Donate blood
# Give money to homeless
# Give to charity
# Return money
# Tax evasion
# Government benefit fraud
# Lie to police

# only available at region-level, thus only region level analysis
# using individual-level data from GSS

###############################################################################################################

### Load data
external_validity_indiv = read.csv("../indicators/GSS/external_validity_indiv.csv")

### Indicator data without study char, but with external validity measures

final_data_VIF_without_study_char = final_data_VIF%>%
  select(-study_ID, -year:-N, -yi, -vi)%>%
  filter(!is.na(region))

### Reverse indicators
external_validity_indiv_reversed = external_validity_indiv %>%
  mutate(donate_blood_reversed = reverse.code(-1, donate_blood, 1, 6),
         give_money_homeless_reversed = reverse.code(-1, give_money_homeless, 1, 6),
         give_money_charity_reversed = reverse.code(-1, give_money_charity, 1, 6),
         return_change_reversed = reverse.code(-1, return_change, 1, 6),
         lie_police_moral_reversed = reverse.code(-1, lie_police_moral, 1, 4))%>%
  select(-donate_blood, -give_money_homeless, -give_money_charity, -return_change, -lie_police_moral)%>%
  rename(donate_blood=donate_blood_reversed, give_money_homeless=give_money_homeless_reversed, give_money_charity=give_money_charity_reversed,
         return_change=return_change_reversed, lie_police_moral=lie_police_moral_reversed)

# Aggregation: First aggregate by state, then aggregate by region (if immediately aggregating by state, 
# there would be bias toward those states with more studies, eg NY) 

external_validity_indicators_state = final_data_VIF_without_study_char %>%
  group_by(state)%>%
  summarise(mean_church_exp_West = mean(church_exp_West, na.rm=TRUE),
            mean_church_exp_East = mean(church_exp_East, na.rm=TRUE),
            mean_percent_farms_operated_by_owner_zscore = mean(percent_farms_operated_by_owner_zscore, na.rm=TRUE),
            mean_historical_subsistence_index = mean(historical_subsistence_index, na.rm=TRUE),
            mean_subsistence_index = mean(subsistence_index, na.rm=TRUE),
            mean_collectivism = mean(collectivism, na.rm=TRUE),
            mean_self_expression = mean(self_expression, na.rm=TRUE),
            mean_PC_religion = mean(PC_religion, na.rm=TRUE),
            mean_PC2_voter_turnout_and_corruption_control_perception = mean(PC2_voter_turnout_and_corruption_control_perception, na.rm=TRUE),
            mean_PC3_institutional_trust_and_perception_of_electoral_integrity = mean(PC3_institutional_trust_and_perception_of_electoral_integrity, na.rm=TRUE),
            mean_PC4_state_pension_funding_and_corruption_convictions = mean(PC4_state_pension_funding_and_corruption_convictions, na.rm=TRUE),
            mean_PC5_economic_wealth_and_inequality_and_health_insurance = mean(PC5_economic_wealth_and_inequality_and_health_insurance, na.rm=TRUE),
            mean_PC1_cousin_marriage = mean(PC1_cousin_marriage, na.rm=TRUE),
            mean_PC2_property_crime_and_caloric_suitability = mean(PC2_property_crime_and_caloric_suitability, na.rm=TRUE),
            mean_PC3_disaster_freq_cost_and_demanding_geoclimate = mean(PC3_disaster_freq_cost_and_demanding_geoclimate, na.rm=TRUE))

region_data = read.csv("../indicators/state_region_codes.csv")
region_data = select(region_data, state, region)

external_validity_indicators_state = left_join(external_validity_indicators_state, region_data, by="state")
external_validity_indicators_state[external_validity_indicators_state=="NaN"] = NA

external_validity_indicators_region = external_validity_indicators_state%>%
  group_by(region)%>%
  summarise(mean_church_exp_West = mean(mean_church_exp_West, na.rm=TRUE),
            mean_church_exp_East = mean(mean_church_exp_East, na.rm=TRUE),
            mean_percent_farms_operated_by_owner_zscore = mean(mean_percent_farms_operated_by_owner_zscore, na.rm=TRUE),
            mean_historical_subsistence_index = mean(mean_historical_subsistence_index, na.rm=TRUE),
            mean_subsistence_index = mean(mean_subsistence_index, na.rm=TRUE),
            mean_collectivism = mean(mean_collectivism, na.rm=TRUE),
            mean_self_expression = mean(mean_self_expression, na.rm=TRUE),
            mean_PC_religion = mean(mean_PC_religion, na.rm=TRUE),
            mean_PC2_voter_turnout_and_corruption_control_perception = mean(mean_PC2_voter_turnout_and_corruption_control_perception, na.rm=TRUE),
            mean_PC3_institutional_trust_and_perception_of_electoral_integrity = mean(mean_PC3_institutional_trust_and_perception_of_electoral_integrity, na.rm=TRUE),
            mean_PC4_state_pension_funding_and_corruption_convictions = mean(mean_PC4_state_pension_funding_and_corruption_convictions, na.rm=TRUE),
            mean_PC5_economic_wealth_and_inequality_and_health_insurance = mean(mean_PC5_economic_wealth_and_inequality_and_health_insurance, na.rm=TRUE),
            mean_PC1_cousin_marriage = mean(mean_PC1_cousin_marriage, na.rm=TRUE),
            mean_PC2_property_crime_and_caloric_suitability = mean(mean_PC2_property_crime_and_caloric_suitability, na.rm=TRUE),
            mean_PC3_disaster_freq_cost_and_demanding_geoclimate = mean(mean_PC3_disaster_freq_cost_and_demanding_geoclimate, na.rm=TRUE))


#### Since indicators are heavily correlated, use VIF to exclude multicollinear vars
x = select(external_validity_indicators_region, -region)
vif(x)
iter1 = select(x, -mean_PC1_cousin_marriage)
vif(iter1)
iter2 = select(iter1, -mean_PC2_property_crime_and_caloric_suitability)
vif(iter2)
iter3 = select(iter2, -mean_historical_subsistence_index)
vif(iter3)
iter4 = select(iter3, -mean_church_exp_East)
vif(iter4)
iter5 = select(iter4, -mean_PC2_voter_turnout_and_corruption_control_perception)
vif(iter5)
iter6 = select(iter5, -mean_PC3_institutional_trust_and_perception_of_electoral_integrity)
vif(iter6)
iter7 = select(iter6, -mean_PC4_state_pension_funding_and_corruption_convictions)
vif(iter7)
iter8 = select(iter7, -mean_self_expression)
vif(iter8)
iter9 = select(iter8, -mean_percent_farms_operated_by_owner_zscore)
vif(iter9)
# Now all VIF < 10

external_validity_indicators_region_VIF = select(external_validity_indicators_region, -mean_PC1_cousin_marriage, -mean_PC2_property_crime_and_caloric_suitability, -mean_historical_subsistence_index, -mean_church_exp_East, -mean_PC2_voter_turnout_and_corruption_control_perception, -mean_PC3_institutional_trust_and_perception_of_electoral_integrity, -mean_PC4_state_pension_funding_and_corruption_convictions, -mean_self_expression, -mean_percent_farms_operated_by_owner_zscore)

# Merge indicators with individual-level data
external_validity = left_join(external_validity_indiv_reversed, external_validity_indicators_region, by="region")
external_validity_VIF = left_join(external_validity_indiv_reversed, external_validity_indicators_region_VIF, by="region")

# Make one df per EV indicator (without VIF)
donate_blood = select(external_validity, -tax_cheating_moral, -gov_benefits_cheating_moral, -give_money_homeless, -give_money_charity, -return_change, -lie_police_moral)
donate_blood = donate_blood[!is.na(donate_blood$donate_blood), ]
give_money_homeless = select(external_validity, -tax_cheating_moral, -gov_benefits_cheating_moral, -donate_blood, -give_money_charity, -return_change, -lie_police_moral)
give_money_homeless = give_money_homeless[!is.na(give_money_homeless$give_money_homeless), ]
give_money_charity = select(external_validity, -tax_cheating_moral, -gov_benefits_cheating_moral, -donate_blood, -give_money_homeless, -return_change, -lie_police_moral)
give_money_charity = give_money_charity[!is.na(give_money_charity$give_money_charity), ]
return_change = select(external_validity, -tax_cheating_moral, -gov_benefits_cheating_moral, -donate_blood, -give_money_homeless, -give_money_charity, -lie_police_moral)
return_change = return_change[!is.na(return_change$return_change), ]
tax_cheating_moral = select(external_validity, -gov_benefits_cheating_moral, -donate_blood, -give_money_homeless, -give_money_charity, -return_change, -lie_police_moral)
tax_cheating_moral = tax_cheating_moral[!is.na(tax_cheating_moral$tax_cheating_moral), ]
gov_benefits_cheating_moral = select(external_validity, -tax_cheating_moral, -donate_blood, -give_money_homeless, -give_money_charity, -return_change, -lie_police_moral)
gov_benefits_cheating_moral = gov_benefits_cheating_moral[!is.na(gov_benefits_cheating_moral$gov_benefits_cheating_moral), ]
lie_police_moral = select(external_validity, -tax_cheating_moral, -gov_benefits_cheating_moral, -donate_blood, -give_money_homeless, -give_money_charity, -return_change)
lie_police_moral = lie_police_moral[!is.na(lie_police_moral$lie_police_moral), ]

# Make one df per EV indicator (with VIF)
donate_blood_VIF = select(external_validity_VIF, -tax_cheating_moral, -gov_benefits_cheating_moral, -give_money_homeless, -give_money_charity, -return_change, -lie_police_moral)
donate_blood_VIF = donate_blood_VIF[!is.na(donate_blood_VIF$donate_blood), ]
give_money_homeless_VIF = select(external_validity_VIF, -tax_cheating_moral, -gov_benefits_cheating_moral, -donate_blood, -give_money_charity, -return_change, -lie_police_moral)
give_money_homeless_VIF = give_money_homeless_VIF[!is.na(give_money_homeless_VIF$give_money_homeless), ]
give_money_charity_VIF = select(external_validity_VIF, -tax_cheating_moral, -gov_benefits_cheating_moral, -donate_blood, -give_money_homeless, -return_change, -lie_police_moral)
give_money_charity_VIF = give_money_charity_VIF[!is.na(give_money_charity_VIF$give_money_charity), ]
return_change_VIF = select(external_validity_VIF, -tax_cheating_moral, -gov_benefits_cheating_moral, -donate_blood, -give_money_homeless, -give_money_charity, -lie_police_moral)
return_change_VIF = return_change_VIF[!is.na(return_change_VIF$return_change), ]
tax_cheating_moral_VIF = select(external_validity_VIF, -gov_benefits_cheating_moral, -donate_blood, -give_money_homeless, -give_money_charity, -return_change, -lie_police_moral)
tax_cheating_moral_VIF = tax_cheating_moral_VIF[!is.na(tax_cheating_moral_VIF$tax_cheating_moral), ]
gov_benefits_cheating_moral_VIF = select(external_validity_VIF, -tax_cheating_moral, -donate_blood, -give_money_homeless, -give_money_charity, -return_change, -lie_police_moral)
gov_benefits_cheating_moral_VIF = gov_benefits_cheating_moral_VIF[!is.na(gov_benefits_cheating_moral_VIF$gov_benefits_cheating_moral), ]
lie_police_moral_VIF = select(external_validity_VIF, -tax_cheating_moral, -gov_benefits_cheating_moral, -donate_blood, -give_money_homeless, -give_money_charity, -return_change)
lie_police_moral_VIF = lie_police_moral_VIF[!is.na(lie_police_moral_VIF$lie_police_moral), ]

### Making new imputed data set not necessary, because no missing data at region-level

### 1) Tax cheating moral

# Full model (with VIF)
fit_tax_cheating_moral_VIF = lm(tax_cheating_moral ~ mean_church_exp_West+mean_subsistence_index+mean_collectivism+mean_PC_religion + mean_PC5_economic_wealth_and_inequality_and_health_insurance + mean_PC3_disaster_freq_cost_and_demanding_geoclimate, data=tax_cheating_moral_VIF)
summary(fit_tax_cheating_moral_VIF)
# Nothing significant

# Theory-level models (without VIF)
fit_tax_cheating_moral = lm(tax_cheating_moral ~ mean_church_exp_West + mean_church_exp_East, data=tax_cheating_moral)
summary(fit_tax_cheating_moral)
# significant: More church West, more C; more church east, less C

fit_tax_cheating_moral = lm(tax_cheating_moral ~ mean_percent_farms_operated_by_owner_zscore + mean_historical_subsistence_index + mean_subsistence_index, data=tax_cheating_moral)
summary(fit_tax_cheating_moral)

fit_tax_cheating_moral = lm(tax_cheating_moral ~ mean_collectivism, data=tax_cheating_moral)
summary(fit_tax_cheating_moral)
# significant: More collectivism, less C

fit_tax_cheating_moral = lm(tax_cheating_moral ~ mean_self_expression, data=tax_cheating_moral)
summary(fit_tax_cheating_moral)
# significant: More self_expression, more C

fit_tax_cheating_moral = lm(tax_cheating_moral ~ mean_PC_religion, data=tax_cheating_moral)
summary(fit_tax_cheating_moral)

fit_tax_cheating_moral = lm(tax_cheating_moral ~ mean_PC2_voter_turnout_and_corruption_control_perception+mean_PC3_institutional_trust_and_perception_of_electoral_integrity+mean_PC4_state_pension_funding_and_corruption_convictions+mean_PC5_economic_wealth_and_inequality_and_health_insurance, data=tax_cheating_moral)
summary(fit_tax_cheating_moral)
# Significant: Lower economic wealth and inequality, more without healthcare (but better fiscal condition) - more C

fit_tax_cheating_moral = lm(tax_cheating_moral ~ mean_PC1_cousin_marriage, data=tax_cheating_moral)
summary(fit_tax_cheating_moral)

fit_tax_cheating_moral = lm(tax_cheating_moral ~ mean_PC2_property_crime_and_caloric_suitability+mean_PC3_disaster_freq_cost_and_demanding_geoclimate, data=tax_cheating_moral)
summary(fit_tax_cheating_moral)
# Marginally significant: More physical threat, more C

### 2) Donating blood

# Full model
fit_donate_blood_VIF = lm(donate_blood ~ mean_church_exp_West+mean_subsistence_index+mean_collectivism+mean_PC_religion + mean_PC5_economic_wealth_and_inequality_and_health_insurance + mean_PC3_disaster_freq_cost_and_demanding_geoclimate, data=donate_blood_VIF)
summary(fit_donate_blood_VIF)
# Nothing significant

# Theory-level models
fit_donate_blood = lm(donate_blood ~ mean_church_exp_West + mean_church_exp_East, data=donate_blood)
summary(fit_donate_blood)

fit_donate_blood = lm(donate_blood ~ mean_percent_farms_operated_by_owner_zscore + mean_historical_subsistence_index + mean_subsistence_index, data=donate_blood)
summary(fit_donate_blood)

fit_donate_blood = lm(donate_blood ~ mean_collectivism, data=donate_blood)
summary(fit_donate_blood)

fit_donate_blood = lm(donate_blood ~ mean_self_expression, data=donate_blood)
summary(fit_donate_blood)

fit_donate_blood = lm(donate_blood ~ mean_PC_religion, data=donate_blood)
summary(fit_donate_blood)

fit_donate_blood = lm(donate_blood ~ mean_PC2_voter_turnout_and_corruption_control_perception+mean_PC3_institutional_trust_and_perception_of_electoral_integrity+mean_PC4_state_pension_funding_and_corruption_convictions+mean_PC5_economic_wealth_and_inequality_and_health_insurance, data=donate_blood)
summary(fit_donate_blood)
# Significant: More economic wealth and inequality, more with healthcare (but worse fiscal condition) - less C

fit_donate_blood = lm(donate_blood ~ mean_PC1_cousin_marriage, data=donate_blood)
summary(fit_donate_blood)

fit_donate_blood = lm(donate_blood ~ mean_PC2_property_crime_and_caloric_suitability+mean_PC3_disaster_freq_cost_and_demanding_geoclimate, data=donate_blood)
summary(fit_donate_blood)

### 3) Give money to homeless

# Full model
fit_give_money_homeless_VIF = lm(give_money_homeless ~ mean_church_exp_West+mean_subsistence_index+mean_collectivism+mean_PC_religion + mean_PC5_economic_wealth_and_inequality_and_health_insurance + mean_PC3_disaster_freq_cost_and_demanding_geoclimate, data=give_money_homeless_VIF)
summary(fit_give_money_homeless_VIF)
# Significant:
# More collectivism - more C
# More religiosity - less C
# More physical threat - more C
# Marginally sign:
# More economic wealth and inequality, more with healthcare (but worse fiscal condition) - less C

# Theory-level models
fit_give_money_homeless = lm(give_money_homeless ~ mean_church_exp_West + mean_church_exp_East, data=give_money_homeless)
summary(fit_give_money_homeless)
# significant: More church West, less C

fit_give_money_homeless = lm(give_money_homeless ~ mean_percent_farms_operated_by_owner_zscore + mean_historical_subsistence_index + mean_subsistence_index, data=give_money_homeless)
summary(fit_give_money_homeless)
# Significant: More current herding - more C; more percent farms operated by owner - less C

fit_give_money_homeless = lm(give_money_homeless ~ mean_collectivism, data=give_money_homeless)
summary(fit_give_money_homeless)
# Significant: More collectivism - more C

fit_give_money_homeless = lm(give_money_homeless ~ mean_self_expression, data=give_money_homeless)
summary(fit_give_money_homeless)

fit_give_money_homeless = lm(give_money_homeless ~ mean_PC_religion, data=give_money_homeless)
summary(fit_give_money_homeless)

fit_give_money_homeless = lm(give_money_homeless ~ mean_PC2_voter_turnout_and_corruption_control_perception+mean_PC3_institutional_trust_and_perception_of_electoral_integrity+mean_PC4_state_pension_funding_and_corruption_convictions+mean_PC5_economic_wealth_and_inequality_and_health_insurance, data=give_money_homeless)
summary(fit_give_money_homeless)
# Significant: 
# More economic wealth and inequality, more with healthcare (but worse fiscal condition) - less C
# More institutional_trust_and_perception_of_electoral_integrity - less C
# Higher voter turnout and corruption control perception - less C
# Less state pension funding and more corruption convictions - less C

fit_give_money_homeless = lm(give_money_homeless ~ mean_PC1_cousin_marriage, data=give_money_homeless)
summary(fit_give_money_homeless)
# Higher cousin marriage rate - less C

fit_give_money_homeless = lm(give_money_homeless ~ mean_PC2_property_crime_and_caloric_suitability+mean_PC3_disaster_freq_cost_and_demanding_geoclimate, data=give_money_homeless)
summary(fit_give_money_homeless)
# Significant: Higher property crime rate and caloric suitability - more C
# Significant: More threat from natural disaster - more C


### 4) Give money to charity

# Full model
fit_give_money_charity_VIF = lm(give_money_charity ~ mean_church_exp_West+mean_subsistence_index+mean_collectivism+mean_PC_religion + mean_PC5_economic_wealth_and_inequality_and_health_insurance + mean_PC3_disaster_freq_cost_and_demanding_geoclimate, data=give_money_charity_VIF)
summary(fit_give_money_charity_VIF)
# Marginally significant:
# More mean_church_exp_West, more C
# More mean_collectivism, more C
# More religiosity - less C

# Theory-level models
fit_give_money_charity = lm(give_money_charity ~ mean_church_exp_West + mean_church_exp_East, data=give_money_charity)
summary(fit_give_money_charity)
# significant: More church West, more C

fit_give_money_charity = lm(give_money_charity ~ mean_percent_farms_operated_by_owner_zscore + mean_historical_subsistence_index + mean_subsistence_index, data=give_money_charity)
summary(fit_give_money_charity)

fit_give_money_charity = lm(give_money_charity ~ mean_collectivism, data=give_money_charity)
summary(fit_give_money_charity)

fit_give_money_charity = lm(give_money_charity ~ mean_self_expression, data=give_money_charity)
summary(fit_give_money_charity)

fit_give_money_charity = lm(give_money_charity ~ mean_PC_religion, data=give_money_charity)
summary(fit_give_money_charity)

fit_give_money_charity = lm(give_money_charity ~ mean_PC2_voter_turnout_and_corruption_control_perception+mean_PC3_institutional_trust_and_perception_of_electoral_integrity+mean_PC4_state_pension_funding_and_corruption_convictions+mean_PC5_economic_wealth_and_inequality_and_health_insurance, data=give_money_charity)
summary(fit_give_money_charity)
# Significant: Higher voter turnout and corruption control perception - more C
# Marginally sign: Less institutional trust - more C

fit_give_money_charity = lm(give_money_charity ~ mean_PC1_cousin_marriage, data=give_money_charity)
summary(fit_give_money_charity)

fit_give_money_charity = lm(give_money_charity ~ mean_PC2_property_crime_and_caloric_suitability+mean_PC3_disaster_freq_cost_and_demanding_geoclimate, data=give_money_charity)
summary(fit_give_money_charity)

### 5) Return change

# Full model
fit_return_change_VIF = lm(return_change ~ mean_church_exp_West+mean_subsistence_index+mean_collectivism+mean_PC_religion + mean_PC5_economic_wealth_and_inequality_and_health_insurance + mean_PC3_disaster_freq_cost_and_demanding_geoclimate, data=return_change_VIF)
summary(fit_return_change_VIF)
# Nothing significant

# Theory-level models
fit_return_change = lm(return_change ~ mean_church_exp_West + mean_church_exp_East, data=return_change)
summary(fit_return_change)

fit_return_change = lm(return_change ~ mean_percent_farms_operated_by_owner_zscore + mean_historical_subsistence_index + mean_subsistence_index, data=return_change)
summary(fit_return_change)

fit_return_change = lm(return_change ~ mean_collectivism, data=return_change)
summary(fit_return_change)

fit_return_change = lm(return_change ~ mean_self_expression, data=return_change)
summary(fit_return_change)

fit_return_change = lm(return_change ~ mean_PC_religion, data=return_change)
summary(fit_return_change)

fit_return_change = lm(return_change ~ mean_PC2_voter_turnout_and_corruption_control_perception+mean_PC3_institutional_trust_and_perception_of_electoral_integrity+mean_PC4_state_pension_funding_and_corruption_convictions+mean_PC5_economic_wealth_and_inequality_and_health_insurance, data=return_change)
summary(fit_return_change)

fit_return_change = lm(return_change ~ mean_PC1_cousin_marriage, data=return_change)
summary(fit_return_change)
# Significant: Higher cousin marriage rates - less C

fit_return_change = lm(return_change ~ mean_PC2_property_crime_and_caloric_suitability+mean_PC3_disaster_freq_cost_and_demanding_geoclimate, data=return_change)
summary(fit_return_change)

### 6) gov benefits cheating moral

# Full model
fit_gov_benefits_cheating_moral_VIF = lm(gov_benefits_cheating_moral ~ mean_church_exp_West+mean_subsistence_index+mean_collectivism+mean_PC_religion + mean_PC5_economic_wealth_and_inequality_and_health_insurance + mean_PC3_disaster_freq_cost_and_demanding_geoclimate, data=gov_benefits_cheating_moral_VIF)
summary(fit_gov_benefits_cheating_moral_VIF)
# Nothing significant

# Theory-level models
fit_gov_benefits_cheating_moral = lm(gov_benefits_cheating_moral ~ mean_church_exp_West + mean_church_exp_East, data=gov_benefits_cheating_moral)
summary(fit_gov_benefits_cheating_moral)

fit_gov_benefits_cheating_moral = lm(gov_benefits_cheating_moral ~ mean_percent_farms_operated_by_owner_zscore + mean_historical_subsistence_index + mean_subsistence_index, data=gov_benefits_cheating_moral)
summary(fit_gov_benefits_cheating_moral)

fit_gov_benefits_cheating_moral = lm(gov_benefits_cheating_moral ~ mean_collectivism, data=gov_benefits_cheating_moral)
summary(fit_gov_benefits_cheating_moral)

fit_gov_benefits_cheating_moral = lm(gov_benefits_cheating_moral ~ mean_self_expression, data=gov_benefits_cheating_moral)
summary(fit_gov_benefits_cheating_moral)
# significant: More self_expression, more C

fit_gov_benefits_cheating_moral = lm(gov_benefits_cheating_moral ~ mean_PC_religion, data=gov_benefits_cheating_moral)
summary(fit_gov_benefits_cheating_moral)

fit_gov_benefits_cheating_moral = lm(gov_benefits_cheating_moral ~ mean_PC2_voter_turnout_and_corruption_control_perception+mean_PC3_institutional_trust_and_perception_of_electoral_integrity+mean_PC4_state_pension_funding_and_corruption_convictions+mean_PC5_economic_wealth_and_inequality_and_health_insurance, data=gov_benefits_cheating_moral)
summary(fit_gov_benefits_cheating_moral)

fit_gov_benefits_cheating_moral = lm(gov_benefits_cheating_moral ~ mean_PC1_cousin_marriage, data=gov_benefits_cheating_moral)
summary(fit_gov_benefits_cheating_moral)

fit_gov_benefits_cheating_moral = lm(gov_benefits_cheating_moral ~ mean_PC2_property_crime_and_caloric_suitability+mean_PC3_disaster_freq_cost_and_demanding_geoclimate, data=gov_benefits_cheating_moral)
summary(fit_gov_benefits_cheating_moral)


### 7) lie police moral

# Full model
fit_lie_police_moral_VIF = lm(lie_police_moral ~ mean_church_exp_West+mean_subsistence_index+mean_collectivism+mean_PC_religion + mean_PC5_economic_wealth_and_inequality_and_health_insurance + mean_PC3_disaster_freq_cost_and_demanding_geoclimate, data=lie_police_moral_VIF)
summary(fit_lie_police_moral_VIF)
# Nothing significant

# Theory-level models
fit_lie_police_moral = lm(lie_police_moral ~ mean_church_exp_West + mean_church_exp_East, data=lie_police_moral)
summary(fit_lie_police_moral)

fit_lie_police_moral = lm(lie_police_moral ~ mean_percent_farms_operated_by_owner_zscore + mean_historical_subsistence_index + mean_subsistence_index, data=lie_police_moral)
summary(fit_lie_police_moral)

fit_lie_police_moral = lm(lie_police_moral ~ mean_collectivism, data=lie_police_moral)
summary(fit_lie_police_moral)

fit_lie_police_moral = lm(lie_police_moral ~ mean_self_expression, data=lie_police_moral)
summary(fit_lie_police_moral)

fit_lie_police_moral = lm(lie_police_moral ~ mean_PC_religion, data=lie_police_moral)
summary(fit_lie_police_moral)

fit_lie_police_moral = lm(lie_police_moral ~ mean_PC2_voter_turnout_and_corruption_control_perception+mean_PC3_institutional_trust_and_perception_of_electoral_integrity+mean_PC4_state_pension_funding_and_corruption_convictions+mean_PC5_economic_wealth_and_inequality_and_health_insurance, data=lie_police_moral)
summary(fit_lie_police_moral)

fit_lie_police_moral = lm(lie_police_moral ~ mean_PC1_cousin_marriage, data=lie_police_moral)
summary(fit_lie_police_moral)

fit_lie_police_moral = lm(lie_police_moral ~ mean_PC2_property_crime_and_caloric_suitability+mean_PC3_disaster_freq_cost_and_demanding_geoclimate, data=lie_police_moral)
summary(fit_lie_police_moral)


#####################################################################################################################

### Correlations between baseline-corrected yi and external validity measures

# Load baseline-corrected yis
yis_region = read.csv("data_frames/regions.csv")

# Calculate region-level external validity measures
external_validity_agg = external_validity_indiv_reversed%>%
  group_by(region)%>%
  summarise(mean_donate_blood = mean(donate_blood, na.rm=TRUE),
            mean_give_money_homeless = mean(give_money_homeless, na.rm=TRUE),
            mean_give_money_charity = mean(give_money_charity, na.rm=TRUE),
            mean_return_change = mean(return_change, na.rm=TRUE),
            mean_tax_cheating_moral = mean(tax_cheating_moral, na.rm=TRUE),
            mean_gov_benefits_cheating_moral = mean(gov_benefits_cheating_moral, na.rm=TRUE),
            mean_lie_police_moral = mean(lie_police_moral, na.rm=TRUE))

# Join dfs
correlations = left_join(yis_region, external_validity_agg, by="region")

# yi and Donate blood
cor(correlations$region_level_effects_controlled_for_baseline, correlations$mean_donate_blood)
# R = 0.7001744
cor.test(correlations$region_level_effects_controlled_for_baseline, correlations$mean_donate_blood)
# p-value = 0.0357
ggplot(correlations, aes(x=mean_donate_blood, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(group=1)) + 
  geom_text_repel(aes(label=region,color=region), segment.size = 0.5)+
  annotate("text", x = 1.18, y = -0.35, label = "r = 0.70; p < 0.05")+
  labs(x="Donated blood", y = "yi controlled for baseline")+
  theme(legend.position = "none") +
  ggsave("plots/external_validity/corr_yi_donate_blood.png", width=6, height=5)

# yi and Give money to homeless
cor(correlations$region_level_effects_controlled_for_baseline, correlations$mean_give_money_homeless)
# R = 0.313042
cor.test(correlations$region_level_effects_controlled_for_baseline, correlations$mean_give_money_homeless)
# p-value = 0.4121
ggplot(correlations, aes(x=mean_give_money_homeless, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(group=1)) + 
  geom_text_repel(aes(label=region,color=region), segment.size = 0.1)+
  annotate("text", x = 2.27, y = -0.35, label = "r = 0.31; p = 0.41")+
  labs(x="Gave money to homeless",y = "yi controlled for baseline")+
  theme(legend.position = "none") +
  ggsave("plots/external_validity/corr_yi_give_money_homeless.png", width=6, height=5)

# yi and give to charity
cor(correlations$region_level_effects_controlled_for_baseline, correlations$mean_give_money_charity)
# R = -0.6149299
cor.test(correlations$region_level_effects_controlled_for_baseline, correlations$mean_give_money_charity)
# p-value = 0.078
ggplot(correlations, aes(x=mean_give_money_charity, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(group=1)) + 
  geom_text_repel(aes(label=region,color=region), segment.size = 0.8)+
  annotate("text", x = 2.73, y = -0.3, label = "r = -0.61; p = 0.08")+
  labs(x="Gave money to charity",y = "yi controlled for baseline")+
  theme(legend.position = "none") +
  ggsave("plots/external_validity/corr_yi_give_money_charity.png", width=6, height=5)

# yi and Return change
cor(correlations$region_level_effects_controlled_for_baseline, correlations$mean_return_change)
# R = 0.6037173
cor.test(correlations$region_level_effects_controlled_for_baseline, correlations$mean_return_change)
# p-value = 0.08516
ggplot(correlations, aes(x=mean_return_change, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(group=1)) + 
  geom_text_repel(aes(label=region,color=region), segment.size = 0.1)+
  annotate("text", x = 1.64, y = -0.3, label = "r = 0.60; p = 0.09")+
  labs(x="Returned extra change to cashier",y = "yi controlled for baseline")+
  theme(legend.position = "none") +
  ggsave("plots/external_validity/corr_yi_return_change.png", width=6, height=5)

# yi and Tax evasion
cor(correlations$region_level_effects_controlled_for_baseline, correlations$mean_tax_cheating_moral)
# R = -0.5749673
cor.test(correlations$region_level_effects_controlled_for_baseline, correlations$mean_tax_cheating_moral)
# p-value = 0.1053
ggplot(correlations, aes(x=mean_tax_cheating_moral, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(group=1)) + 
  geom_text_repel(aes(label=region,color=region), segment.size = 0.3)+
  annotate("text", x = 3.05, y = -0.3, label = "r = -0.57; p = 0.11")+
  labs(x="Disapproval of tax cheating",y = "yi controlled for baseline")+
  theme(legend.position = "none") +
  ggsave("plots/external_validity/corr_yi_tax_cheating.png", width=6, height=5)

# yi and Government benefit fraud
cor(correlations$region_level_effects_controlled_for_baseline, correlations$mean_gov_benefits_cheating_moral)
# R = -0.5114251
cor.test(correlations$region_level_effects_controlled_for_baseline, correlations$mean_gov_benefits_cheating_moral)
# p-value = 0.1593
ggplot(correlations, aes(x=mean_gov_benefits_cheating_moral, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(group=1)) + 
  geom_text_repel(aes(label=region,color=region), segment.size = 0.3)+
  annotate("text", x = 3.42, y = -0.35, label = "r = -0.51; p = 0.16")+
  labs(x="Disapproval of government benefit cheating", y = "yi controlled for baseline")+
  theme(legend.position = "none") +
  ggsave("plots/external_validity/corr_yi_gov_benefits_cheating.png", width=6, height=5)

# yi and Lie to police
cor(correlations$region_level_effects_controlled_for_baseline, correlations$mean_lie_police_moral)
# R = 0.1359385
cor.test(correlations$region_level_effects_controlled_for_baseline, correlations$mean_lie_police_moral)
# p-value = 0.7273
ggplot(correlations, aes(x=mean_lie_police_moral, y=region_level_effects_controlled_for_baseline, color=region)) + 
  geom_point() + 
  geom_smooth(method=lm, aes(group=1)) + 
  geom_text_repel(aes(label=region,color=region), segment.size = 0.5)+
  annotate("text", x = 3.25, y = -0.45, label = "r = 0.14; p = 0.73")+
  labs(x="Disapproval of lying to police", y = "yi controlled for baseline")+
  theme(legend.position = "none") +
  ggsave("plots/external_validity/corr_yi_lie_police.png", width=6, height=5)



