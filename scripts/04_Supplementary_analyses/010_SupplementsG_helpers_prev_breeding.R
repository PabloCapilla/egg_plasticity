###
###
#' 
#' Script for:
#' Mothers in a cooperatively breeding bird increase investment per offspring at the pre-natal stage when they will have more help with post-natal care
#' Capilla-Lasheras et al. 
#' 
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list=ls())

##
##
##### Script aim: #####
#' Analysis of egg volume, investigating effects of the number of helpers in the previous breeding attempt and
#' the time since last breeding attempt (Supplementary Materials E).
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
## dataset for analysis
data <- data %>%
  filter(!is.na(female_helper_prev_ba)) %>% 
  filter(!is.na(time_since_prev_ba)) %>%
  filter(fate_prev_ba   == "SUCCESS")   # dataset includes only successful previous breeding events

##
## sample sizes
length(unique(data$clutch_ID)) # n clutches
length(unique(data$mother_ID)) # n mothers
length(unique(data$Group)) # n groups 

##
##
##### Correlation betweeen number of helpers in current and previous breeding attempt #####
##
##
cor.test(data$female_helpers, data$female_helper_prev_ba) # for female helpers
cor.test(data$male_helpers, data$male_helper_prev_ba) # for female helpers


##
##
##### Population-level egg volume model - including time since last breeding event and number of helpers in previous breeding even #####
##
##

#' We use AIC here to assess the level of support for two models including different sets of predictors that represent two 
#' alternative hypothesis.
#' Specifically, we test whether effects of past help (model 'model_carryover_effects_prev_help') explain the data better 
#' than effects of current help (model 'model_carryover_effects_current_help').

##
## model with current help
model_carryover_effects_current_help <- lmer(egg_volume ~ 
                                               scale(time_since_prev_ba) : female_helpers +
                                               scale(time_since_prev_ba) : male_helpers +
                                               
                                               # environmental indexes
                                               poly(rainfall,2)[,1] +
                                               poly(rainfall,2)[,2] +
                                               temp_above_35 +
                                               
                                               
                                               # carryover effect predictors main effect
                                               scale(time_since_prev_ba) +
                                               female_helpers +
                                               male_helpers +
                                               
                                               # other main effects
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
AIC(model_carryover_effects_current_help)
summary(model_carryover_effects_current_help)

##
## model with help in previous breeding attempt
model_carryover_effects_prev_help <- lmer(egg_volume ~ 
                                            scale(time_since_prev_ba) : female_helper_prev_ba +
                                            scale(time_since_prev_ba) : male_helper_prev_ba +
                                            # environmental indexes
                                            poly(rainfall,2)[,1] +
                                            poly(rainfall,2)[,2] +
                                            temp_above_35 +
                                            
                                            
                                            # carryover effect predictors main effect
                                            scale(time_since_prev_ba) +
                                            female_helper_prev_ba +
                                            male_helper_prev_ba +
                                            
                                            # other main effects
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
AIC(model_carryover_effects_prev_help)
summary(model_carryover_effects_prev_help)

AIC(model_carryover_effects_current_help) - AIC(model_carryover_effects_prev_help)
