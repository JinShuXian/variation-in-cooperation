#####################################################################################################

######## Analysis: Preprocessing

### 1) Make final_data df: Year merging
### 2) PCA (+ direction of PCA components)
### 3) Moderator correlations
### 4) VIF test: Explude multicollinear variables
### 5) Check distribution of cooperation rates


### Output:
# final_data.csv
# final_data_VIF.csv
# density_plot_observed_P_C.png
# moderator_correlation.csv


### Load packages
library("tidyverse")
library("lme4")
library("faraway")
library("MuMIn")
library("psych")

# Set seed
set.seed(11)

### Set WD
setwd("/home/caroline/Desktop/cross_cultural_metaanalysis/US/analysis")


### Read indicator_data
indicators = read.csv("../indicators/0_indicator_preprocessing_merging/indicator_descriptives_plotting/indicators_reversed_zscore_no_outliers_with_composite.csv")
indicators = filter(indicators, !is.na(indicators$state))


#########################################################################################################

#### Merge indicator data with CoDa data and for each CoDa study choose the nearest year

## Add indicators with year information based on year of data collection.
# If year not available for indicator, choose year-1,
# if year-1 is not available, year+1, o.w. year-2, o.w. year+2 etc.

indicators_with_year_info = filter(indicators, !is.na(year))
indicators_with_year_info = indicators_with_year_info[colSums(!is.na(indicators_with_year_info)) > 0]

#indicators_with_year_info1 = select(indicators_with_year_info, state, year, percent_cropland_zscore)
#indicators_with_year_info2 = select(indicators_with_year_info, state, year, percent_herdingland_zscore)
indicators_with_year_info3 = select(indicators_with_year_info, state, year, YLL_communicable_diseases_zscore)
indicators_with_year_info4 = select(indicators_with_year_info, state, year, violent_crime_rate_zscore)
indicators_with_year_info5 = select(indicators_with_year_info, state, year, property_crime_rate_zscore)
indicators_with_year_info6 = select(indicators_with_year_info, state, year, infant_mortality_rate_zscore)
indicators_with_year_info7 = select(indicators_with_year_info, state, year, GDP_per_capita_zscore)
#indicators_with_year_info8 = select(indicators_with_year_info, state, year, cash_solveny_zscore)
#indicators_with_year_info9 = select(indicators_with_year_info, state, year, budget_solvency_zscore)
#indicators_with_year_info10 = select(indicators_with_year_info, state, year, long_run_solvency_zscore)
indicators_with_year_info11 = select(indicators_with_year_info, state, year, voter_turnout_zscore)
indicators_with_year_info12 = select(indicators_with_year_info, state, year, GINI_index_zscore)
indicators_with_year_info25 = select(indicators_with_year_info, state, year, state_pension_funded_ratio_zscore)
indicators_with_year_info30 = select(indicators_with_year_info, state, year, religious_attend_percent_zscore)
#indicators_with_year_info31 = select(indicators_with_year_info, state, year, homosexuality_moral_zscore)
indicators_with_year_info32 = select(indicators_with_year_info, state, year, pray_percent_zscore)
indicators_with_year_info35 = select(indicators_with_year_info, state, year, percent_uninsured_reversed_zscore)
indicators_with_year_info36 = select(indicators_with_year_info, state, year, life_expectancy_at_birth_reversed_zscore)
indicators_with_year_info41 = select(indicators_with_year_info, state, year, trust_reversed_zscore)
indicators_with_year_info42 = select(indicators_with_year_info, state, year, conf_court_reversed_zscore)
indicators_with_year_info43 = select(indicators_with_year_info, state, year, conf_government_reversed_zscore)
indicators_with_year_info44 = select(indicators_with_year_info, state, year, conf_parliament_reversed_zscore)
indicators_with_year_info45 = select(indicators_with_year_info, state, year, conf_army_reversed_zscore)
indicators_with_year_info46 = select(indicators_with_year_info, state, year, belief_life_after_death_reversed_zscore)
#indicators_with_year_info50 = select(indicators_with_year_info, state, year, happy_reversed_zscore)
indicators_with_year_info51 = select(indicators_with_year_info, state, year, fiscal_condition_combined_zscore)
indicators_with_year_info52 = select(indicators_with_year_info, state, year, subsistence_index)
indicators_with_year_info53 = select(indicators_with_year_info, state, year, self_expression_combined_zscore)

#indicators_with_year_info1 = filter(indicators_with_year_info1, !(is.na(percent_cropland_zscore)))
#indicators_with_year_info2 = filter(indicators_with_year_info2, !(is.na(percent_herdingland_zscore)))
indicators_with_year_info3 = filter(indicators_with_year_info3, !(is.na(YLL_communicable_diseases_zscore)))
indicators_with_year_info4 = filter(indicators_with_year_info4, !(is.na(violent_crime_rate_zscore)))
indicators_with_year_info5 = filter(indicators_with_year_info5, !(is.na(property_crime_rate_zscore)))
indicators_with_year_info6 = filter(indicators_with_year_info6, !(is.na(infant_mortality_rate_zscore)))
indicators_with_year_info7 = filter(indicators_with_year_info7, !(is.na(GDP_per_capita_zscore)))
#indicators_with_year_info8 = filter(indicators_with_year_info8, !(is.na(cash_solveny_zscore)))
#indicators_with_year_info9 = filter(indicators_with_year_info9, !(is.na(budget_solvency_zscore)))
#indicators_with_year_info10 = filter(indicators_with_year_info10, !(is.na(long_run_solvency_zscore)))
indicators_with_year_info11 = filter(indicators_with_year_info11, !(is.na(voter_turnout_zscore)))
indicators_with_year_info12 = filter(indicators_with_year_info12, !(is.na(GINI_index_zscore)))
indicators_with_year_info25 = filter(indicators_with_year_info25, !(is.na(state_pension_funded_ratio_zscore)))
indicators_with_year_info30 = filter(indicators_with_year_info30, !(is.na(religious_attend_percent_zscore)))
#indicators_with_year_info31 = filter(indicators_with_year_info31, !(is.na(homosexuality_moral_zscore)))
indicators_with_year_info32 = filter(indicators_with_year_info32, !(is.na(pray_percent_zscore)))
indicators_with_year_info35 = filter(indicators_with_year_info35, !(is.na(percent_uninsured_reversed_zscore)))
indicators_with_year_info36 = filter(indicators_with_year_info36, !(is.na(life_expectancy_at_birth_reversed_zscore)))
indicators_with_year_info41 = filter(indicators_with_year_info41, !(is.na(trust_reversed_zscore)))
indicators_with_year_info42 = filter(indicators_with_year_info42, !(is.na(conf_court_reversed_zscore)))
indicators_with_year_info43 = filter(indicators_with_year_info43, !(is.na(conf_government_reversed_zscore)))
indicators_with_year_info44 = filter(indicators_with_year_info44, !(is.na(conf_parliament_reversed_zscore)))
indicators_with_year_info45 = filter(indicators_with_year_info45, !(is.na(conf_army_reversed_zscore)))
indicators_with_year_info46 = filter(indicators_with_year_info46, !(is.na(belief_life_after_death_reversed_zscore)))
#indicators_with_year_info50 = filter(indicators_with_year_info50, !(is.na(happy_reversed_zscore)))
indicators_with_year_info51 = filter(indicators_with_year_info51, !(is.na(fiscal_condition_combined_zscore)))
indicators_with_year_info52 = filter(indicators_with_year_info52, !(is.na(subsistence_index)))
indicators_with_year_info53 = filter(indicators_with_year_info53, !(is.na(self_expression_combined_zscore)))

# Load US CoDa data
final_coda_data_US = read.csv("../data/final_coda_data_US.csv")

# Select and rename variables
final_coda_data_US_for_indicators = final_coda_data_US%>%
  select(study_ID, year_data_collection, state)%>%
  rename(year=year_data_collection)

# Find closest year and add variable to final_coda_data_US_for_indicators
for(i in 1:nrow(final_coda_data_US_for_indicators)) {
  current_year = final_coda_data_US_for_indicators[i, "year"]
  current_state = final_coda_data_US_for_indicators[i, "state"]
  temp_filter = filter(indicators_with_year_info3, state==as.character(current_state))
  current_index = which.min(abs(temp_filter$year - current_year))
  final_coda_data_US_for_indicators[i, "YLL_communicable_diseases_zscore"] = ifelse(is.null(dim(temp_filter)), NA, temp_filter[current_index, "YLL_communicable_diseases_zscore"])
  temp_filter = filter(indicators_with_year_info4, state==as.character(current_state))
  current_index = which.min(abs(temp_filter$year - current_year))
  final_coda_data_US_for_indicators[i, "violent_crime_rate_zscore"] = ifelse(is.null(dim(temp_filter)), NA, temp_filter[current_index, "violent_crime_rate_zscore"])
  temp_filter = filter(indicators_with_year_info5, state==as.character(current_state))
  current_index = which.min(abs(temp_filter$year - current_year))
  final_coda_data_US_for_indicators[i, "property_crime_rate_zscore"] = ifelse(is.null(dim(temp_filter)), NA, temp_filter[current_index, "property_crime_rate_zscore"])
  temp_filter = filter(indicators_with_year_info6, state==as.character(current_state))
  current_index = which.min(abs(temp_filter$year - current_year))
  final_coda_data_US_for_indicators[i, "infant_mortality_rate_zscore"] = ifelse(is.null(dim(temp_filter)), NA, temp_filter[current_index, "infant_mortality_rate_zscore"])
  temp_filter = filter(indicators_with_year_info7, state==as.character(current_state))
  current_index = which.min(abs(temp_filter$year - current_year))
  final_coda_data_US_for_indicators[i, "GDP_per_capita_zscore"] = ifelse(is.null(dim(temp_filter)), NA, temp_filter[current_index, "GDP_per_capita_zscore"])
  temp_filter = filter(indicators_with_year_info11, state==as.character(current_state))
  current_index = which.min(abs(temp_filter$year - current_year))
  final_coda_data_US_for_indicators[i, "voter_turnout_zscore"] = ifelse(is.null(dim(temp_filter)), NA, temp_filter[current_index, "voter_turnout_zscore"])
  temp_filter = filter(indicators_with_year_info12, state==as.character(current_state))
  current_index = which.min(abs(temp_filter$year - current_year))
  final_coda_data_US_for_indicators[i, "GINI_index_zscore"] = ifelse(is.null(dim(temp_filter)), NA, temp_filter[current_index, "GINI_index_zscore"])
  temp_filter = filter(indicators_with_year_info25, state==as.character(current_state))
  current_index = which.min(abs(temp_filter$year - current_year))
  final_coda_data_US_for_indicators[i, "state_pension_funded_ratio_zscore"] = ifelse(is.null(dim(temp_filter)), NA, temp_filter[current_index, "state_pension_funded_ratio_zscore"])
  temp_filter = filter(indicators_with_year_info30, state==as.character(current_state))
  current_index = which.min(abs(temp_filter$year - current_year))
  final_coda_data_US_for_indicators[i, "religious_attend_percent_zscore"] = ifelse(is.null(dim(temp_filter)), NA, temp_filter[current_index, "religious_attend_percent_zscore"])
  temp_filter = filter(indicators_with_year_info32, state==as.character(current_state))
  current_index = which.min(abs(temp_filter$year - current_year))
  final_coda_data_US_for_indicators[i, "pray_percent_zscore"] = ifelse(is.null(dim(temp_filter)), NA, temp_filter[current_index, "pray_percent_zscore"])
  temp_filter = filter(indicators_with_year_info35, state==as.character(current_state))
  current_index = which.min(abs(temp_filter$year - current_year))
  final_coda_data_US_for_indicators[i, "percent_uninsured_reversed_zscore"] = ifelse(is.null(dim(temp_filter)), NA, temp_filter[current_index, "percent_uninsured_reversed_zscore"])
  temp_filter = filter(indicators_with_year_info36, state==as.character(current_state))
  current_index = which.min(abs(temp_filter$year - current_year))
  final_coda_data_US_for_indicators[i, "life_expectancy_at_birth_reversed_zscore"] = ifelse(is.null(dim(temp_filter)), NA, temp_filter[current_index, "life_expectancy_at_birth_reversed_zscore"])
  temp_filter = filter(indicators_with_year_info41, state==as.character(current_state))
  current_index = which.min(abs(temp_filter$year - current_year))
  final_coda_data_US_for_indicators[i, "trust_reversed_zscore"] = ifelse(is.null(dim(temp_filter)), NA, temp_filter[current_index, "trust_reversed_zscore"])
  temp_filter = filter(indicators_with_year_info42, state==as.character(current_state))
  current_index = which.min(abs(temp_filter$year - current_year))
  final_coda_data_US_for_indicators[i, "conf_court_reversed_zscore"] = ifelse(is.null(dim(temp_filter)), NA, temp_filter[current_index, "conf_court_reversed_zscore"])
  temp_filter = filter(indicators_with_year_info43, state==as.character(current_state))
  current_index = which.min(abs(temp_filter$year - current_year))
  final_coda_data_US_for_indicators[i, "conf_government_reversed_zscore"] = ifelse(is.null(dim(temp_filter)), NA, temp_filter[current_index, "conf_government_reversed_zscore"])
  temp_filter = filter(indicators_with_year_info44, state==as.character(current_state))
  current_index = which.min(abs(temp_filter$year - current_year))
  final_coda_data_US_for_indicators[i, "conf_parliament_reversed_zscore"] = ifelse(is.null(dim(temp_filter)), NA, temp_filter[current_index, "conf_parliament_reversed_zscore"])
  temp_filter = filter(indicators_with_year_info45, state==as.character(current_state))
  current_index = which.min(abs(temp_filter$year - current_year))
  final_coda_data_US_for_indicators[i, "conf_army_reversed_zscore"] = ifelse(is.null(dim(temp_filter)), NA, temp_filter[current_index, "conf_army_reversed_zscore"])
  temp_filter = filter(indicators_with_year_info46, state==as.character(current_state))
  current_index = which.min(abs(temp_filter$year - current_year))
  final_coda_data_US_for_indicators[i, "belief_life_after_death_reversed_zscore"] = ifelse(is.null(dim(temp_filter)), NA, temp_filter[current_index, "belief_life_after_death_reversed_zscore"])
  temp_filter = filter(indicators_with_year_info51, state==as.character(current_state))
  current_index = which.min(abs(temp_filter$year - current_year))
  final_coda_data_US_for_indicators[i, "fiscal_condition_combined_zscore"] = ifelse(is.null(dim(temp_filter)), NA, temp_filter[current_index, "fiscal_condition_combined_zscore"])
  temp_filter = filter(indicators_with_year_info52, state==as.character(current_state))
  current_index = which.min(abs(temp_filter$year - current_year))
  final_coda_data_US_for_indicators[i, "subsistence_index"] = ifelse(is.null(dim(temp_filter)), NA, temp_filter[current_index, "subsistence_index"])
  temp_filter = filter(indicators_with_year_info53, state==as.character(current_state))
  current_index = which.min(abs(temp_filter$year - current_year))
  final_coda_data_US_for_indicators[i, "self_expression_combined_zscore"] = ifelse(is.null(dim(temp_filter)), NA, temp_filter[current_index, "self_expression_combined_zscore"])
}

## Add indicators which have no year information
indicators_no_year_info = filter(indicators, is.na(year))
indicators_no_year_info = indicators_no_year_info[colSums(!is.na(indicators_no_year_info)) > 0]
indicators_no_year_info = select(indicators_no_year_info, -signed_petition_reversed_zscore)

final_indicator_data = left_join(final_coda_data_US_for_indicators, indicators_no_year_info, by="state")

################################################################################################################

#### Check INDICATOR correlations

# Pearson's correlation
only_indicators = final_indicator_data%>%
  select(-study_ID:-state, -state_code, -region, -donate_blood_reversed_zscore:-lie_police_moral_reversed_zscore, -tax_cheating_moral_zscore, -gov_benefits_cheating_moral_zscore)%>%
  mutate_if(is.factor, as.character)%>%
  mutate_if(is.character, as.numeric)
indicator_correlation = cor(only_indicators, method = "pearson", use = "pairwise.complete.obs")
write.csv(indicator_correlation, "data_frames/indicator_correlation.csv", row.names = FALSE)



#################################################################################################################

## PCA FOR EACH HYPOTHESIS CLUSTER

### !!! Cut-off is: > 75% explained (cumulative) variance

###H1 Religion
components<-prcomp(~religious_attend_percent_zscore+pray_percent_zscore+belief_life_after_death_reversed_zscore+belief_heaven_reversed_zscore+belief_hell_reversed_zscore, data = final_indicator_data, na.action = na.exclude)
summary(components)
head(unclass(components$rotation)[, 1:4])
H1_religion_new<-components$x[,1:2]
H1_religion_new <- as.data.frame(H1_religion_new)
H1_religion_new <- select(H1_religion_new, PC1)
H1_religion_new <- rename(H1_religion_new, PC_religion=PC1) 

###H2 Formal institutions
components<-prcomp(~conf_army_reversed_zscore+conf_court_reversed_zscore+conf_government_reversed_zscore+conf_parliament_reversed_zscore+GDP_per_capita_zscore+GINI_index_zscore+voter_turnout_zscore+state_pension_funded_ratio_zscore+percent_uninsured_reversed_zscore+perception_electoral_integrity_zscore+percent_invested_school_maintenance_zscore+corruption_convictions_reversed_zscore+corruption_perception_reversed_zscore+percent_roads_in_poor_condition_reversed_zscore+fiscal_condition_combined_zscore, data = final_indicator_data, na.action = na.exclude)
summary(components)
H2_formal_institutions_new<-components$x[,1:5]
H2_formal_institutions_new <- as.data.frame(H2_formal_institutions_new)
H2_formal_institutions_new <- rename(H2_formal_institutions_new, PC1_government_effectiveness = PC1, PC2_voter_turnout_and_corruption_control_perception = PC2, PC3_institutional_trust_and_perception_of_electoral_integrity = PC3, PC4_state_pension_funding_and_corruption_convictions = PC4, PC5_economic_wealth_and_inequality_and_health_insurance = PC5)

###H3 Informal Institutions
components<-prcomp(~KII_zscore+percent_cousin_marriage_log_zscore, data = final_indicator_data, na.action = na.exclude)
summary(components)
H3_informal_institutions_new<-components$x[,1:2]
H3_informal_institutions_new <- as.data.frame(H3_informal_institutions_new)
H3_informal_institutions_new <- select(H3_informal_institutions_new, PC1)
H3_informal_institutions_new <- rename(H3_informal_institutions_new, PC1_cousin_marriage=PC1)

###H4 Threat
components<-prcomp(~YLL_communicable_diseases_zscore+violent_crime_rate_zscore+property_crime_rate_zscore+infant_mortality_rate_zscore+life_expectancy_at_birth_reversed_zscore+demanding_geoclimate_zscore+natural_disaster_cost_zscore+natural_disaster_freq_zscore+caloric_suitability_post1500_reversed_zscore, data = final_indicator_data, na.action = na.exclude)
summary(components)
H4_threat_new<-components$x[,1:3]
H4_threat_new <- as.data.frame(H4_threat_new)
#H4_threat_new <- select(H6_threat_new, PC1)
#H4_threat_new <- rename(H4_threat_new, PC1_mortality_and_natural_disaster_freq=PC1, PC2_caloric_suitability_natural_disaster_cost=PC2, PC3_crime_diseases_geoclimate=PC3)
H4_threat_new <- rename(H4_threat_new, PC1_disease_mortality_and_violent_crime=PC1, PC2_property_crime_and_caloric_suitability=PC2, PC3_disaster_freq_cost_and_demanding_geoclimate=PC3)

indicators_PCA = bind_cols(H1_religion_new, H2_formal_institutions_new, H3_informal_institutions_new, H4_threat_new)

## 3) Merge all indicators for final indicator df:

non_pca_indicators = select(final_indicator_data, - religious_attend_percent_zscore, -pray_percent_zscore, -belief_life_after_death_reversed_zscore, -belief_heaven_reversed_zscore, -belief_hell_reversed_zscore, -conf_army_reversed_zscore, -conf_court_reversed_zscore, -conf_government_reversed_zscore, -conf_parliament_reversed_zscore, -GDP_per_capita_zscore, -GINI_index_zscore, -voter_turnout_zscore, -state_pension_funded_ratio_zscore, -percent_uninsured_reversed_zscore, -perception_electoral_integrity_zscore, -percent_invested_school_maintenance_zscore, -corruption_convictions_reversed_zscore, -corruption_perception_reversed_zscore, -percent_roads_in_poor_condition_reversed_zscore, -fiscal_condition_combined_zscore, -YLL_communicable_diseases_zscore, -violent_crime_rate_zscore, -property_crime_rate_zscore, -infant_mortality_rate_zscore, -life_expectancy_at_birth_reversed_zscore, -demanding_geoclimate_zscore, -natural_disaster_cost_zscore, -natural_disaster_freq_zscore, -caloric_suitability_post1500_reversed_zscore, -KII_zscore, -percent_cousin_marriage_log_zscore)

final_indicators_after_PCA = bind_cols(non_pca_indicators, indicators_PCA)

##############################################################################################

### Check direction of PCA components for documentation and interpretation of model results

pca_indicators = select(final_indicator_data, religious_attend_percent_zscore, pray_percent_zscore, belief_life_after_death_reversed_zscore, belief_heaven_reversed_zscore, belief_hell_reversed_zscore, conf_army_reversed_zscore, conf_court_reversed_zscore, conf_government_reversed_zscore, conf_parliament_reversed_zscore, GDP_per_capita_zscore, GINI_index_zscore, voter_turnout_zscore, state_pension_funded_ratio_zscore, percent_uninsured_reversed_zscore, perception_electoral_integrity_zscore, percent_invested_school_maintenance_zscore, corruption_convictions_reversed_zscore, corruption_perception_reversed_zscore, percent_roads_in_poor_condition_reversed_zscore, fiscal_condition_combined_zscore, YLL_communicable_diseases_zscore, violent_crime_rate_zscore, property_crime_rate_zscore, infant_mortality_rate_zscore, life_expectancy_at_birth_reversed_zscore, demanding_geoclimate_zscore, natural_disaster_cost_zscore, natural_disaster_freq_zscore, caloric_suitability_post1500_reversed_zscore, KII_zscore, percent_cousin_marriage_log_zscore)
all_indicators_after_PCA = bind_cols(final_indicators_after_PCA, pca_indicators)

# Religion
cor(all_indicators_after_PCA$PC_religion, all_indicators_after_PCA$religious_attend_percent_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC_religion, all_indicators_after_PCA$pray_percent_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC_religion, all_indicators_after_PCA$belief_life_after_death_reversed_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC_religion, all_indicators_after_PCA$belief_heaven_reversed_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC_religion, all_indicators_after_PCA$belief_hell_reversed_zscore, use = "pairwise.complete.obs")
# Higher PC_religion --> lower religiosity

# Formal institutions
cor(all_indicators_after_PCA$PC1_government_effectiveness, all_indicators_after_PCA$percent_roads_in_poor_condition_reversed_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC1_government_effectiveness, all_indicators_after_PCA$percent_invested_school_maintenance_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC2_voter_turnout_and_corruption_control_perception, all_indicators_after_PCA$corruption_perception_reversed_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC2_voter_turnout_and_corruption_control_perception, all_indicators_after_PCA$voter_turnout_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC3_institutional_trust_and_perception_of_electoral_integrity, all_indicators_after_PCA$conf_army_reversed_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC3_institutional_trust_and_perception_of_electoral_integrity, all_indicators_after_PCA$conf_court_reversed_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC3_institutional_trust_and_perception_of_electoral_integrity, all_indicators_after_PCA$conf_government_reversed_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC3_institutional_trust_and_perception_of_electoral_integrity, all_indicators_after_PCA$conf_parliament_reversed_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC3_institutional_trust_and_perception_of_electoral_integrity, all_indicators_after_PCA$perception_electoral_integrity_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC4_state_pension_funding_and_corruption_convictions, all_indicators_after_PCA$state_pension_funded_ratio_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC4_state_pension_funding_and_corruption_convictions, all_indicators_after_PCA$corruption_convictions_reversed_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC5_economic_wealth_and_inequality_and_health_insurance, all_indicators_after_PCA$GINI_index_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC5_economic_wealth_and_inequality_and_health_insurance, all_indicators_after_PCA$GDP_per_capita_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC5_economic_wealth_and_inequality_and_health_insurance, all_indicators_after_PCA$percent_uninsured_reversed_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC5_economic_wealth_and_inequality_and_health_insurance, all_indicators_after_PCA$fiscal_condition_combined_zscore, use = "pairwise.complete.obs")
# Higher PC1_government_effectiveness --> better roads but worse schools
# Higher PC2_voter_turnout_and_corruption_control_perception --> higher voter turnout and corruption control perception
# Higher PC3_institutional_trust_and_perception_of_electoral_integrity --> lower institutional trust and perception of electoral integrity (but higher trust in army)
# Higher PC4_state_pension_funding_and_corruption_convictions --> less state pension funding and more corruption convictions
# Higher PC5_economic_wealth_and_inequality_and_health_insurance --> lower economic wealth and inequality, more without healthcare (but better fiscal condition)

# Informal institutions
cor(all_indicators_after_PCA$PC1_cousin_marriage, all_indicators_after_PCA$KII_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC1_cousin_marriage, all_indicators_after_PCA$percent_cousin_marriage_log_zscore, use = "pairwise.complete.obs")
# Higher PC1_cousin_marriage --> lower cousin marriage rates (but higher kinship intensity index)

# Threat
cor(all_indicators_after_PCA$PC1_disease_mortality_and_violent_crime, all_indicators_after_PCA$YLL_communicable_diseases_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC1_disease_mortality_and_violent_crime, all_indicators_after_PCA$violent_crime_rate_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC1_disease_mortality_and_violent_crime, all_indicators_after_PCA$infant_mortality_rate_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC1_disease_mortality_and_violent_crime, all_indicators_after_PCA$life_expectancy_at_birth_reversed_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC2_property_crime_and_caloric_suitability, all_indicators_after_PCA$property_crime_rate_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC2_property_crime_and_caloric_suitability, all_indicators_after_PCA$caloric_suitability_post1500_reversed_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC3_disaster_freq_cost_and_demanding_geoclimate, all_indicators_after_PCA$demanding_geoclimate_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC3_disaster_freq_cost_and_demanding_geoclimate, all_indicators_after_PCA$natural_disaster_freq_zscore, use = "pairwise.complete.obs")
cor(all_indicators_after_PCA$PC3_disaster_freq_cost_and_demanding_geoclimate, all_indicators_after_PCA$natural_disaster_cost_zscore, use = "pairwise.complete.obs")
# Higher PC1_disease_mortality_and_violent_crime --> higher disease mortality and violent crime rates
# Higher PC2_property_crime_and_caloric_suitability --> lower property crime rate and caloric suitability
# Higher PC3_disaster_freq_cost_and_demanding_geoclimate --> higher physical environmental pressure

#########################################################################################################################

#### Reverse those PCA indicators which were flipped by PCA 
# (e.g. religiosity: All religion indicators are negatively correlated with PC_religion)

all_indicators_after_PCA_reversed = all_indicators_after_PCA %>%
  mutate(PC_religion = reverse.code(-1, PC_religion),
         PC3_institutional_trust_and_perception_of_electoral_integrity = reverse.code(-1, PC3_institutional_trust_and_perception_of_electoral_integrity),
         PC4_state_pension_funding_and_corruption_convictions = reverse.code(-1, PC4_state_pension_funding_and_corruption_convictions),
         PC5_economic_wealth_and_inequality_and_health_insurance = reverse.code(-1, PC5_economic_wealth_and_inequality_and_health_insurance),
         PC1_cousin_marriage = reverse.code(-1, PC1_cousin_marriage),
         PC2_property_crime_and_caloric_suitability = reverse.code(-1, PC2_property_crime_and_caloric_suitability))

#########################################################################################################################

#### Make full CoDa and indicator data frame

## Load CODA data
final_coda_data = read.csv("../data/final_coda_data_US.csv")

#### Merge indicators and data
final_data = bind_cols(final_coda_data, all_indicators_after_PCA_reversed)

# Select columns
final_data = select(final_data, study_ID, state, region, year, symmetry:sanction, P_C, N, trust_reversed_zscore:self_expression_combined_zscore, ChurchExpWest_zscore:PC3_disaster_freq_cost_and_demanding_geoclimate)

################################################################################################

###### Calculate effect size and variance
final_data$yi = log(final_data$P_C / (1 - final_data$P_C))
final_data$vi = sqrt( (1 / final_data$N * final_data$P_C) + (1 / final_data$N * (1- final_data$P_C)))

# Take log of group size
final_data$group_size = log(final_data$group_size)

### Convert dichotomous variables, study_ID and country_ID to factors
final_data$matching = factor(final_data$matching)
final_data$discussion = factor(final_data$discussion)
final_data$simultaneous = factor(final_data$simultaneous)
final_data$sanction = factor(final_data$sanction)
final_data$symmetry = factor(final_data$symmetry)
final_data$one_shot = factor(final_data$one_shot)
final_data$study_ID = factor(final_data$study_ID)
final_data$state = factor(final_data$state)
final_data$region = factor(final_data$region)

## Rename to be consistent with Guiliana Global data

final_data = rename(final_data, residential_mob=total_percent_moved_zscore,
                    church_exp_West=ChurchExpWest_zscore,
                    church_exp_East=ChurchExpEast_zscore,
                    collectivism=collectivism_zscore,
                    self_expression=self_expression_combined_zscore)%>%
  select(-return_change_reversed_zscore, -lie_police_moral_reversed_zscore, -tax_cheating_moral_zscore, -gov_benefits_cheating_moral_zscore, -donate_blood_reversed_zscore, -give_money_homeless_reversed_zscore, -give_money_charity_reversed_zscore)
#-N, -P_C, 

############################################################################################

# Check cooperation distribution
ggplot(final_data, aes(x=P_C)) +
  geom_density(alpha=0.3) + 
  labs(x = "Observed P(C)", y = "Density") +
  ggsave("plots/density_plot_observed_P_C.png", width=8, height=5.5)


#################################################################################################################


#### Check MODERATOR correlations

library(corrplot)

# Pearson's correlation
only_moderators = final_data%>%
  select(-study_ID:-year, -N)%>%
  mutate_if(is.factor, as.character)%>%
  mutate_if(is.character, as.numeric)
moderator_correlation = cor(only_moderators, method = "pearson", use = "pairwise.complete.obs")
write.csv(moderator_correlation, "data_frames/moderator_correlation.csv", row.names = FALSE)

# possibly problematic: (correlation higher than .8/-.8)
# - matching and one-shot (.88)
# - tightness_zscore and PC_Religion

# possibly problematic: (correlation higher than .7/-.7)
# PC1_disease_mortality_and_violent_crime and PC_Religion
# PC3_institutional_trust_and_perception_of_electoral_integrity and trust

# Plot correlations (no labels to save space)
corrplot(moderator_correlation, method="color", cl.pos = "n", tl.pos = "n")

#############################################################################################

### VIF

# "Some papers argue that a VIF<10 is acceptable, but others says that the limit value is 5.
# - "10" as the maximum level of VIF (Hair et al., 1995)
# - "5" as the maximum level of VIF (Ringle et al., 2015)"

only_moderators = select(only_moderators, -P_C, -yi, -vi)

# Check VIF
vif(only_moderators)
# Delete multicollinear variables in a step-wise fashion (in each step, delete 1 var with highest VIF)
iter1 = select(only_moderators, -tightness_zscore)
vif(iter1)
iter2 = select(iter1, -residential_mob)
vif(iter2)
iter3 = select(iter2, -PC1_disease_mortality_and_violent_crime)
vif(iter3)
iter4 = select(iter3, -PC1_government_effectiveness)
vif(iter4)
# Now all VIF < 7.5

# Update final data accordingly
final_data_VIF = select(final_data, -tightness_zscore, -residential_mob, -PC1_disease_mortality_and_violent_crime, -PC1_government_effectiveness)

final_data_VIF$region = as.character(final_data_VIF$region)
final_data_VIF$region = ifelse(final_data_VIF$state=="District of Columbia", "South Atlantic", final_data_VIF$region)

write.csv(final_data_VIF, "data_frames/final_data_VIF.csv", row.names = FALSE)


##########################################################################################################

### Final data descriptives for thesis

# 394 treatments/observations
unique(final_data_VIF$study_ID)
# 284 studies
unique(final_data_VIF$state)
# 38 states

# 34899 participants (double treatments from WS designs excluded)

# Structural vars
table(final_data_VIF$symmetry)
table(final_data_VIF$one_shot)
table(final_data_VIF$matching)
table(final_data_VIF$simultaneous)
table(final_data_VIF$sanction)
table(final_data_VIF$discussion)
describe(exp(final_data_VIF$group_size))
describe(final_data_VIF$group_size)
describe(final_data_VIF$if_PD_K_index)

# Crosscultural vars
describe(final_data_VIF$trust_reversed_zscore)
describe(final_data_VIF$subsistence_index)
describe(final_data_VIF$self_expression)
describe(final_data_VIF$church_exp_West)
describe(final_data_VIF$church_exp_East)
describe(final_data_VIF$percent_farms_operated_by_owner_zscore)
describe(final_data_VIF$collectivism)
describe(final_data_VIF$historical_subsistence_index)
describe(final_data_VIF$PC_religion)
describe(final_data_VIF$PC2_voter_turnout_and_corruption_control_perception)
describe(final_data_VIF$PC3_institutional_trust_and_perception_of_electoral_integrity)
describe(final_data_VIF$PC4_state_pension_funding_and_corruption_convictions)
describe(final_data_VIF$PC5_economic_wealth_and_inequality_and_health_insurance)
describe(final_data_VIF$PC1_cousin_marriage)
describe(final_data_VIF$PC2_property_crime_and_caloric_suitability)
describe(final_data_VIF$PC3_disaster_freq_cost_and_demanding_geoclimate)

# NAs
nrow(filter(final_data_VIF, is.na(trust_reversed_zscore)))
nrow(filter(final_data_VIF, is.na(subsistence_index)))
nrow(filter(final_data_VIF, is.na(self_expression)))
nrow(filter(final_data_VIF, is.na(church_exp_West)))
nrow(filter(final_data_VIF, is.na(church_exp_East)))
nrow(filter(final_data_VIF, is.na(percent_farms_operated_by_owner_zscore)))
nrow(filter(final_data_VIF, is.na(collectivism)))
nrow(filter(final_data_VIF, is.na(historical_subsistence_index)))
nrow(filter(final_data_VIF, is.na(PC_religion)))
nrow(filter(final_data_VIF, is.na(PC2_voter_turnout_and_corruption_control_perception)))
nrow(filter(final_data_VIF, is.na(PC3_institutional_trust_and_perception_of_electoral_integrity)))
nrow(filter(final_data_VIF, is.na(PC4_state_pension_funding_and_corruption_convictions)))
nrow(filter(final_data_VIF, is.na(PC5_economic_wealth_and_inequality_and_health_insurance)))
nrow(filter(final_data_VIF, is.na(PC1_cousin_marriage)))
nrow(filter(final_data_VIF, is.na(PC2_property_crime_and_caloric_suitability)))
nrow(filter(final_data_VIF, is.na(PC3_disaster_freq_cost_and_demanding_geoclimate)))


# Other vars
describe(final_data_VIF$N)$kurtosis
describe(final_data_VIF$P_C)
describe(final_data_VIF$yi)
describe(final_data_VIF$vi)

# skew
describe(exp(final_data_VIF$group_size))$skew
describe(final_data_VIF$group_size)$skew
describe(final_data_VIF$if_PD_K_index)$skew
describe(final_data_VIF$trust_reversed_zscore)$skew
describe(final_data_VIF$subsistence_index)$skew
describe(final_data_VIF$self_expression)$skew
describe(final_data_VIF$church_exp_West)$skew
describe(final_data_VIF$church_exp_East)$skew
describe(final_data_VIF$percent_farms_operated_by_owner_zscore)$skew
describe(final_data_VIF$collectivism)$skew
describe(final_data_VIF$historical_subsistence_index)$skew
describe(final_data_VIF$PC_religion)$skew
describe(final_data_VIF$PC2_voter_turnout_and_corruption_control_perception)$skew
describe(final_data_VIF$PC3_institutional_trust_and_perception_of_electoral_integrity)$skew
describe(final_data_VIF$PC4_state_pension_funding_and_corruption_convictions)$skew
describe(final_data_VIF$PC5_economic_wealth_and_inequality_and_health_insurance)$skew
describe(final_data_VIF$PC1_cousin_marriage)$skew
describe(final_data_VIF$PC2_property_crime_and_caloric_suitability)$skew
describe(final_data_VIF$PC3_disaster_freq_cost_and_demanding_geoclimate)$skew
describe(final_data_VIF$N)$skew
describe(final_data_VIF$P_C)$skew
describe(final_data_VIF$yi)$skew
describe(final_data_VIF$vi)$skew

# kurtosis
describe(exp(final_data_VIF$group_size))$kurtosis
describe(final_data_VIF$group_size)$kurtosis
describe(final_data_VIF$if_PD_K_index)$kurtosis
describe(final_data_VIF$trust_reversed_zscore)$kurtosis
describe(final_data_VIF$subsistence_index)$kurtosis
describe(final_data_VIF$self_expression)$kurtosis
describe(final_data_VIF$church_exp_West)$kurtosis
describe(final_data_VIF$church_exp_East)$kurtosis
describe(final_data_VIF$percent_farms_operated_by_owner_zscore)$kurtosis
describe(final_data_VIF$collectivism)$kurtosis
describe(final_data_VIF$historical_subsistence_index)$kurtosis
describe(final_data_VIF$PC_religion)$kurtosis
describe(final_data_VIF$PC2_voter_turnout_and_corruption_control_perception)$kurtosis
describe(final_data_VIF$PC3_institutional_trust_and_perception_of_electoral_integrity)$kurtosis
describe(final_data_VIF$PC4_state_pension_funding_and_corruption_convictions)$kurtosis
describe(final_data_VIF$PC5_economic_wealth_and_inequality_and_health_insurance)$kurtosis
describe(final_data_VIF$PC1_cousin_marriage)$kurtosis
describe(final_data_VIF$PC2_property_crime_and_caloric_suitability)$kurtosis
describe(final_data_VIF$PC3_disaster_freq_cost_and_demanding_geoclimate)$kurtosis
describe(final_data_VIF$N)$kurtosis
describe(final_data_VIF$P_C)$kurtosis
describe(final_data_VIF$yi)$kurtosis
describe(final_data_VIF$vi)$kurtosis