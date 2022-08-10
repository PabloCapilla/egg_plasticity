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
#' Analysis of effects of past help on maternal body condition.
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
##### data ######
##
##
data <- read.csv("./data/SuppC_maternal_condition.csv")
head(data)


##
##
##### Summary of sample sizes #####
##
##
nrow(data) # clutches
length(unique(data$mother_ID)) # mothers
length(unique(data$group_ID))    # groups
length(unique(data$clutch_ID)) 

##
##
##### Model Supplementary Materials C #####
##
##
model_maternal_conditions <- lmer(body_index ~ 
                                    prev_female_number +
                                    prev_male_number +
                                    clutch_order +
                                    (1|season) +
                                    (1|mother_ID) +
                                    (1|group_ID), 
                                  data =data, 
                                  REML = F,
                                  na.action = "na.fail")
summary(model_maternal_conditions)
drop1(model_maternal_conditions, test = "Chisq")

