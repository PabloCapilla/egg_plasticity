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

##
##
##### egg volume dataset #####
##
##
data <- read.csv("./data/maternal_provisioning_dataset.csv")
head(data)
nrow(data)

####

##
##
##### Bivariate model egg volume - female helper number #####
##
##

## prior
prior_prov_female_helper_bivariate <- list(R=list(V=diag(2),nu=0.002),
                                          G=list(G1=list(V=diag(2), nu=1.002),
                                                 G2=list(V=diag(2), nu=1.002),
                                                 G3=list(V=diag(2), nu=1.002))) 

## model
set.seed(7)
prov_female_helper_bivariate <- MCMCglmm(cbind(maternal_prov_rate, female_helpers) ~ 
                                           trait - 1 +
                                           
                                           ## egg model
                                           at.level(trait, 1):poly(rainfall,2)[,1] +
                                           at.level(trait, 1):poly(rainfall,2)[,2] +
                                           at.level(trait, 1):temp_above_35 +
                                           at.level(trait, 1):brood_size +
                                           at.level(trait, 1):male_helpers,
                                         
                                         ## female helper model
                                         
                                         random = 
                                           ~us(trait):group_ID +  
                                           us(trait):season +
                                           us(trait):mother_ID,
                                         rcov =~ us(trait):units,
                                         family = c("gaussian","gaussian"),
                                         prior = prior_prov_female_helper_bivariate,
                                         nitt=50000,
                                         burnin=1000,
                                         thin=10,
                                         verbose = TRUE,
                                         data = data)

summary(prov_female_helper_bivariate)
plot(prov_female_helper_bivariate)

#####

##
## calculating within and among female effects

## among effect
among_effect <- prov_female_helper_bivariate$VCV[,"traitfemale_helpers:traitmaternal_prov_rate.mother_ID"]/ 
     prov_female_helper_bivariate$VCV[,"traitfemale_helpers:traitfemale_helpers.mother_ID"]

mean(among_effect) 
HPDinterval(among_effect)

## within effect
wihtin_effect <- prov_female_helper_bivariate$VCV[,"traitfemale_helpers:traitmaternal_prov_rate.units"]/ 
     prov_female_helper_bivariate$VCV[,"traitfemale_helpers:traitfemale_helpers.units"]

mean(wihtin_effect) 
HPDinterval(wihtin_effect)

## difference
diff_effect <- wihtin_effect - among_effect
mean(diff_effect) 
HPDinterval(diff_effect)

#####

