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
##### libraries #####
##
pacman::p_load(dplyr, 
               tidyr, 
               ggplot2, 
               ggpubr,
               MCMCglmm, 
               lme4)

#####

##
##
##### egg volume dataset #####
##
##
data <- read.csv("./data/egg_volume_dataset.csv") 
data <- data %>% 
  mutate(egg_volume_cm = egg_volume/1000)
head(data)
nrow(data)

####

##
##
##### Bivariate model egg volume - female helper number #####
##
##

## prior
prior_egg_female_helper_bivariate <- list(R=list(V=diag(2),nu=0.002),
                                        G=list(G1=list(V=diag(2), nu=1.002),
                                               G2=list(V=diag(2), nu=1.002),
                                               G3=list(V=diag(2), nu=1.002))) 

## model
set.seed(7)
eggs_female_helper_bivariate <- MCMCglmm(cbind(egg_volume_cm, female_helpers) ~ 
                                           trait - 1 +
                                           
                                           ## egg model
                                           at.level(trait, 1):poly(rainfall,2)[,1] +
                                           at.level(trait, 1):poly(rainfall,2)[,2] +
                                           at.level(trait, 1):egg_position + 
                                           at.level(trait, 1):clutch_size +
                                           at.level(trait, 1):temp_above_35 +
                                           at.level(trait, 1):male_helpers,
                                           
                                           ## female helper model
                                           
                                           random = 
                                           ~us(trait):Group +  
                                           us(trait):Season +
                                           us(trait):mother_ID,
                                         rcov =~ us(trait):units,
                                         family = c("gaussian","gaussian"),
                                         prior = prior_egg_female_helper_bivariate,
                                         nitt=50000,
                                         burnin=1000,
                                         thin=10,
                                         verbose = TRUE,
                                         data = data)

summary(eggs_female_helper_bivariate)
plot(eggs_female_helper_bivariate)

#####

##
## calculating within and among female effects

## among effect
among_effect <- eggs_female_helper_bivariate$VCV[,"traitfemale_helpers:traitegg_volume_cm.mother_ID"]/ 
     eggs_female_helper_bivariate$VCV[,"traitfemale_helpers:traitfemale_helpers.mother_ID"]

mean(among_effect) 
HPDinterval(among_effect)

## within effect
wihtin_effect <- eggs_female_helper_bivariate$VCV[,"traitfemale_helpers:traitegg_volume_cm.units"]/ 
     eggs_female_helper_bivariate$VCV[,"traitfemale_helpers:traitfemale_helpers.units"]

mean(wihtin_effect) 
HPDinterval(wihtin_effect)

## difference
diff_effect <- wihtin_effect - among_effect
mean(diff_effect) 
HPDinterval(diff_effect)

#####
