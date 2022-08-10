###
###
#' 
#' Script for:
#' Mothers front-load their investment to the egg stage when helped in a wild cooperative bird
#' Capilla-Lasheras et al. 
#' Preprint: https://doi.org/10.1101/2021.11.11.468195
#' 
#' Latest update: 2022/08/10
#' 
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list=ls())

##
##
##### Script aim: #####
#' Analysis of egg volume, investigating effects of time since last breeding attempt and number
#' of helpers in previous breeding attempt (Supplementary Materials D).
#' 
##
##

##
##### libraries #####
##
pacman::p_load(dplyr, 
               tidyr, 
               lme4)

##
##
##### egg volume dataset #####
##
##
data <- read.csv("./data/egg_volume_dataset.csv") 
head(data)
nrow(data)


##
##
##### Population-level egg volume model - including time since last breeding event and number of helpers in previous breeding even #####
##
##
data <- data %>%
  filter(!is.na(time_since_prev_ba)) %>%
  filter(fate_prev_ba   == "SUCCESS")   # dataset includes only successful previous breeding events

length(unique(data$clutch_ID)) # n clutches
length(unique(data$mother_ID)) # n mothers
length(unique(data$Group)) # n groups 

summary(data$time_since_prev_ba) # range time since
mean(data$time_since_prev_ba)    # mean time since

##
##
##### Model for effects of time since last breeding attempt on egg volume #####
##
##
model_carryover_effects <- lmer(egg_volume ~ 
                                  
                                  # predictors testing for carry over effects
                                  time_since_prev_ba : female_helpers +
                                  time_since_prev_ba : male_helpers +
                                  
                                  # environmental indexes
                                  poly(rainfall,2)[,1] +
                                  poly(rainfall,2)[,2] +
                                  temp_above_35 +
                                  
                                  
                                  # carryover effect predictors main effect
                                  time_since_prev_ba +
                                  
                                  
                                  # other main effects
                                  female_helpers +
                                  male_helpers +
                                  scale(clutch_size) +
                                  scale(egg_position) +
                                  
                                  # random effects
                                  (1|Group) +
                                  (1|mother_ID) +
                                  (1|Season) + 
                                  (1|clutch_ID),
                                data = data,
                                na.action = "na.fail",
                                REML = F)
summary(model_carryover_effects)
drop1(model_carryover_effects, test = "Chisq")

##
## effect of 'time since last breeding attempt' in a model without interactions
summary(update(model_carryover_effects, . ~ . - 
                 time_since_prev_ba : female_helpers -
                 time_since_prev_ba : male_helpers))
drop1(update(model_carryover_effects, . ~ . - 
               time_since_prev_ba : female_helpers -
               time_since_prev_ba : male_helpers),
      test = "Chisq")
