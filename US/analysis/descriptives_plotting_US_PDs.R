### R script for descriptives and plotting of dichotomous US PDs

### Load packages
library("readxl")
library("tidyverse")
library("dplyr")
library("ggmap")
library("maps")
library("mapdata")
library("devtools")
library("ggplot2")

#### Set current working directory

setwd("/home/caroline/Desktop/cross_cultural_metaanalysis/US/analysis")

# #### OLD CODE IS COMMENTED OUT!!!
# 
# #### Get descriptive data of PD US data set
# 
# ## Load full data set
# study_char_full = read.csv("../data/old_data/study_characteristics_06.03.19.csv")
# #state_data_full = read.csv("state_data_full.csv")
# 
# ## Filter out studies with no overall cooperation rate
# study_char_full_valid = filter(study_char_full, overall_P_C != 999 & overall_P_C != "NA")
# 
# ## Only studies from US:
# US_full_valid = filter(study_char_full_valid, country == "USA")
# 
# ## Filter out invalid studies with no location information (e.g., MTurkers)
# US_full_valid = filter(US_full_valid, location.source != 4, location.source != 5)
# 
# ## Only dich PDG
# 
# US_dich_PDG = filter(US_full_valid, dilemma_types == 1 & num_choice_options == 2)
# 
# # Number of studies
# nrow(US_dich_PDG) # 328 studies
# 
# # Number of participants
# US_dich_PDG$N = as.numeric(as.character(US_dich_PDG$N))
# sum(US_dich_PDG$N) # 41261 subjects
# 
# # Number of states
# unique(US_dich_PDG$region.state) # 39 states
# 
# # Range of years of data collection
# US_dich_PDG$year_data_collection = as.numeric(as.character(US_dich_PDG$year_data_collection))
# table(US_dich_PDG$year_data_collection) # 1958 - 2017 (59 years)
# 
# # Proportion of one-shot PD studies
# temp = filter(US_dich_PDG, one_shot != 999 & one_shot != "0; 1" & one_shot != "1; 0")
# temp$one_shot = as.numeric(as.character(temp$one_shot))
# (sum(temp$one_shot))/nrow(temp) # 19,5% one-shot games
# 
# 
# #### Plotting dichotomous PDs US using US_dich_PDG data frame
# 
# ## Get country levels and their frequency
# US_dich_PDG = rename(US_dich_PDG, region_state=region.state)
# temp = group_by(US_dich_PDG, region_state)
# state_freq = summarise(temp, num_studies = length(region_state))
# 
# ## Maps package requires state names in lower case letters: Modify country_freq accordingly:
# state_freq$region_state = tolower(state_freq$region_state)
# 
# ## Rename region_state to region in order to use state_freq with maps package
# state_freq = rename(state_freq, region=region_state)
# 
# ## Plot states for which we have valid PD data
# 
# # Retrieve map data for US from maps package
# usa <- map_data("state")
# 
# # Select states for which there is data
# PD_US <- subset(usa, region %in% c(state_freq$region))
# 
# # Plot data
# ggplot(data = PD_US) + 
#   geom_polygon(aes(x = long, y = lat, group = group), fill = "blue", color = "black") + 
#   coord_fixed(1.3)
# 
# us_base <- ggplot(data = usa, mapping = aes(x = long, y = lat, group = group)) + 
#   coord_fixed(1.3) + 
#   geom_polygon(color = "black", fill = "gray")
# 
# us_base + #theme_nothing() + 
#   geom_polygon(data = PD_US, fill = "blue", color = "white") +
#   geom_polygon(color = "black", fill = NA)  # get the state border back on top
# 
# 
# ### Plot states in different color according to how many studies per location
# 
# # Add num_regions variable to PD_US
# PD_US_num_studies = left_join(PD_US, state_freq, by="region")
# 
# ## Plot with legend, not overlayed on US base
# ggplot(data = PD_US_num_studies) + 
#   geom_polygon(aes(x = long, y = lat, fill = num_studies, group = group), color = "white") + 
#   coord_fixed(1.3)
# 
# ## Plot with legend, overlayed on US base
# us_base <- ggplot(data = usa, mapping = aes(x = long, y = lat, group = group)) + 
#   coord_fixed(1.3) + 
#   geom_polygon(color = "white", fill = "seashell2")
# 
# us_base + theme_light() + 
#   geom_polygon(data = PD_US_num_studies, aes(x = long, y = lat, fill = num_studies, group = group), color = "white")
#   # + ggsave("plots/PD_US_num_studies.png")
# 
# 
# ### Plot states using log scale of number of studies per location:
# 
# PD_US_num_studies$num_studies_log = log(PD_US_num_studies$num_studies)
# 
# ## Plot with legend, not overlayed on US base
# ggplot(data = PD_US_num_studies) + 
#   geom_polygon(aes(x = long, y = lat, fill = num_studies_log, group = group), color = "white") + 
#   coord_fixed(1.3)
# 
# ## Plot with legend, overlayed on US base
# us_base <- ggplot(data = usa, mapping = aes(x = long, y = lat, group = group)) + 
#   coord_fixed(1.3) + 
#   geom_polygon(color = "white", fill = "seashell2")
# 
# us_base + theme_light() + 
#   geom_polygon(data = PD_US_num_studies, aes(x = long, y = lat, fill = num_studies_log, group = group), color = "white") +
#   ggsave("plots/PD_US_num_studies_log.png")
# 
# 
# 
# #### Check how many studies have no location information (e.g., MTurkers)
# 
# ## Only dichotomous PDG studies from US:
# US_studies = filter(study_char_full_valid, country == "USA")
# US_dich_PDG_studies = filter(US_studies, dilemma_types == 1 & num_choice_options == 2)
# 
# # Number of studies
# nrow(US_dich_PDG_studies) # 413 studies
# 
# ## studies with no location information (e.g., MTurkers)
# US_valid_studies =  filter(US_dich_PDG_studies, location.source != 4, location.source != 5)
# US_invalid_studies = filter(US_dich_PDG_studies, location.source != 1, location.source != 2, location.source != 3)
# 
# ## Plot histogram of year frequencies
# ggplot(data = US_invalid_studies, aes(year_data_collection)) + 
#   geom_bar()
# 
# US_invalid_location = select(US_invalid_studies, study_ID:study_num, year_data_collection, country, overall_P_C, IVs, University:Comments.1)
# 
# write.csv(US_invalid_location, "US_invalid_location.csv")
# 



###############################################################################################################
# Load effect sizes corrected for baseline from analysis
final_data_with_means = read.csv("data_frames/final_data_VIF_with_means.csv")

final_data_with_means = rename(final_data_with_means, state=State, region=Region)

state_effect_sizes = select(final_data_with_means, state, state_level_effects_controlled_for_baseline)
region_effect_sizes = select(final_data_with_means, state, region_level_effects_controlled_for_baseline)
state_effect_sizes = state_effect_sizes[!duplicated(state_effect_sizes$state_level_effects_controlled_for_baseline), ]
region_effect_sizes = region_effect_sizes[!duplicated(region_effect_sizes$state), ]

## Maps package requires state names in lower case letters: Modify accordingly:
state_effect_sizes$state = tolower(state_effect_sizes$state)
region_effect_sizes$state = tolower(region_effect_sizes$state)

## Rename state to region in order to use state_freq with maps package
state_effect_sizes = rename(state_effect_sizes, region=state, yi_controlled_for_baseline=state_level_effects_controlled_for_baseline)
region_effect_sizes = rename(region_effect_sizes, region=state, yi_controlled_for_baseline=region_level_effects_controlled_for_baseline)

# Retrieve map data for US from maps package
usa <- map_data("state")

# Select states for which we have effect sizes
states_map <- subset(usa, region %in% c(state_effect_sizes$region))
regions_map <- subset(usa, region %in% c(region_effect_sizes$region))

# Plot data: 1) States for which we have data
ggplot(data = states_map) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "blue", color = "black") + 
  coord_fixed(1.3)

us_base <- ggplot(data = usa, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color = "black", fill = "gray")

us_base + #theme_nothing() + 
  geom_polygon(data = states_map, fill = "blue", color = "white") +
  geom_polygon(color = "black", fill = NA)  # get the state border back on top

##########################################################################

### Plot states in different color according to effect size

# 1) State-level

# Add effect_size variable
states_map_effects = left_join(states_map, state_effect_sizes, by="region")

## Plot with legend, not overlayed on US base
ggplot(data = states_map_effects) + 
  geom_polygon(aes(x = long, y = lat, fill = yi_controlled_for_baseline, group = group), color = "white") + 
  coord_fixed(1.3)

## Plot with legend, overlayed on US base
us_base <- ggplot(data = usa, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color = "white", fill = "seashell2")

us_base + theme_light() + 
  geom_polygon(data = states_map_effects, aes(x = long, y = lat, fill = yi_controlled_for_baseline, group = group), color = "white") + 
  labs(fill="State-level yi\n(controlled for baseline)")+
  ggsave("plots/states_map_effects.png", width=8, height=5.5)

# 2) Region-level

# Add effect_size variable
regions_map_effects = left_join(regions_map, region_effect_sizes, by="region")

## Plot with legend, not overlayed on US base
ggplot(data = regions_map_effects) + 
  geom_polygon(aes(x = long, y = lat, fill = yi_controlled_for_baseline, group = group), color = "white") + 
  coord_fixed(1.3)

## Plot with legend, overlayed on US base
us_base <- ggplot(data = usa, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color = "white", fill = "seashell2")

us_base + theme_light() + 
  geom_polygon(data = regions_map_effects, aes(x = long, y = lat, fill = yi_controlled_for_baseline, group = group), color = "white") + 
  labs(fill="Region-level yi\n(controlled for baseline)")+
  ggsave("plots/regions_map_effects.png", width=8, height=5.5)


#########################################################################################################

### Plot number of studies per state and region

# Calculate num studies per state

num_studies_per_state = final_data_with_means%>%
  group_by(State)%>%
  summarise(num_studies_per_state=n())%>%
  mutate(num_studies_per_state_log = log(num_studies_per_state))%>%
  mutate(State = tolower(State))%>%
  rename(region=State)

usa <- map_data("state")

num_studies_per_state_plot <- subset(usa, region %in% c(num_studies_per_state$region))

num_studies_per_state_plot = left_join(num_studies_per_state_plot, num_studies_per_state, by="region")

### Plot states in different color according to num studies

# 1) Num studies (no log)

## Plot with legend, not overlayed on US base
ggplot(data = num_studies_per_state_plot) + 
  geom_polygon(aes(x = long, y = lat, fill = num_studies_per_state, group = group), color = "white") + 
  coord_fixed(1.3)

## Plot with legend, overlayed on US base
us_base <- ggplot(data = usa, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color = "white", fill = "seashell2")

us_base + theme_light() + 
  geom_polygon(data = num_studies_per_state_plot, aes(x = long, y = lat, fill = num_studies_per_state, group = group), color = "white") + 
  labs(fill="Number of studies")+
  ggsave("plots/states_num_studies.png", width=8, height=5.5)

# 2) Log num studies

## Plot with legend, not overlayed on US base
ggplot(data = num_studies_per_state_plot) + 
  geom_polygon(aes(x = long, y = lat, fill = num_studies_per_state_log, group = group), color = "white") + 
  coord_fixed(1.3)

## Plot with legend, overlayed on US base
us_base <- ggplot(data = usa, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color = "white", fill = "seashell2")

us_base + theme_light() + 
  geom_polygon(data = num_studies_per_state_plot, aes(x = long, y = lat, fill = num_studies_per_state_log, group = group), color = "white") + 
  labs(fill="Log number of studies")+
  ggsave("plots/states_num_studies_log.png", width=8, height=5.5)

