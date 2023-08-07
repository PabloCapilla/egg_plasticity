###
###
#' 
#' Script for:
#' Mothers front-load their investment to the egg stage when helped in a wild cooperative bird
#' Capilla-Lasheras et al. 
#' Preprint: https://doi.org/10.1101/2021.11.11.468195
#' 
#' Latest update: 2023/07/08
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
               data.table,
               ggplot2, 
               extrafont,
               MCMCglmm,
               lme4)
loadfonts()
# font_import() may be needed first to use font types in ggplot
# check with 'fonts()' that font are available

#####

##
##
##### egg volume dataset #####
##
##
egg <- read.csv("./data/egg_volume_dataset.csv") 
egg <- egg %>% 
  mutate(egg_volume_cm = egg_volume/1000,
         Season = as.character(substr(x = Season, start = 1, stop = 4))) %>% 
  dplyr::select(Season, 
                 Group, 
                 clutch_ID, 
                 clutch_size,
                 mother_ID,
                 egg_position,
                 egg_volume_cm,
                 female_helpers_egg = female_helpers,
                 male_helpers_egg = male_helpers,
                 temp_above_35_egg = temp_above_35,
                 rainfall_egg = rainfall)
head(egg)
nrow(egg)

####

##
##
##### Maternal provisioning rate dataset #####
##
##
frate <- read.csv("./data/maternal_provisioning_dataset.csv")
frate <- frate %>% 
  mutate(Season = as.character(season),
         Group = group_ID) %>% 
  dplyr::select(clutch_ID, 
                mother_ID,
                Group,
                Season,
                brood_size,
                maternal_prov_rate,
                female_helpers_frate = female_helpers,
                male_helpers_frate = male_helpers,
                temp_above_35_frate = temp_above_35,
                rainfall_frate = rainfall)
head(frate)
#####

##
##
##### Combining data sets #####
##
##
data <- full_join(x = egg,
                  y = frate,
                  by = c("clutch_ID", "mother_ID", "Group", "Season"))

head(data, n= 20)

table(is.na(data$Group))
table(is.na(data$Season))
table(is.na(data$mother_ID))
table(is.na(data$clutch_ID))

table(is.na(data$rainfall_egg[!is.na(data$egg_volume_cm)]))
table(is.na(data$egg_position[!is.na(data$egg_volume_cm)]))
table(is.na(data$clutch_size[!is.na(data$egg_volume_cm)]))
table(is.na(data$temp_above_35_egg[!is.na(data$egg_volume_cm)]))
table(is.na(data$female_helpers_egg[!is.na(data$egg_volume_cm)]))
table(is.na(data$clutch_ID[!is.na(data$egg_volume_cm)]))
table(is.na(data$male_helpers_egg[!is.na(data$egg_volume_cm)]))

table(is.na(data$rainfall_frate[!is.na(data$maternal_prov_rate)]))
table(is.na(data$temp_above_35_frate[!is.na(data$maternal_prov_rate)]))
table(is.na(data$female_helpers_frate[!is.na(data$maternal_prov_rate)]))
table(is.na(data$male_helpers_frate[!is.na(data$maternal_prov_rate)]))
table(is.na(data$brood_size[!is.na(data$maternal_prov_rate)]))


nrow(data[!is.na(data$maternal_prov_rate),])
nrow(data[!is.na(data$egg_volume_cm),])

nrow(data[is.na(data$maternal_prov_rate),])
nrow(data[is.na(data$egg_volume_cm),])

#####

##
##
##### bivariate model #####
##
##

## priors
prior_egg_frate_bivariate <- list(R=list(V=matrix(c(1,0,0,0.001), nrow = 2), nu=1.002, fix = 2),
                                  G=list(G1=list(V=diag(2), nu=1.002),
                                         G2=list(V=diag(2), nu=1.002),
                                         G3=list(V=diag(2), nu=1.002),
                                         G4=list(V=diag(2), nu=1.002))) 

##
## model no helpers
set.seed(7)
egg_frate_bivariate_nohelpers <- MCMCglmm(cbind(egg_volume_cm, maternal_prov_rate) ~ 
                                            trait - 1 +
                                            
                                            ## egg model
                                            at.level(trait, 1):I(rainfall_egg^2) +
                                            at.level(trait, 1):rainfall_egg,
                                            at.level(trait, 1):egg_position + 
                                            at.level(trait, 1):clutch_size +
                                            at.level(trait, 1):temp_above_35_egg +
                                            #at.level(trait, 1):female_helpers_egg +
                                            #at.level(trait, 1):male_helpers_egg +
                                            
                                            ## prov model
                                            at.level(trait, 2):I(rainfall_frate^2) +
                                            at.level(trait, 2):rainfall_frate +
                                            at.level(trait, 2):temp_above_35_frate +
                                            #at.level(trait, 2):female_helpers_frate +
                                            #at.level(trait, 1):male_helpers_frate +
                                            at.level(trait, 2):brood_size,
                                          
                                          random = 
                                            ~idh(trait):Group +  
                                            idh(trait):Season +
                                            us(trait):mother_ID + 
                                            us(trait):clutch_ID,
                                          rcov =~ idh(trait):units,
                                          family = c("gaussian","gaussian"),
                                          prior = prior_egg_frate_bivariate,
                                          nitt=5100,
                                          burnin=100,
                                          thin=1,
                                          verbose = TRUE,
                                          data = data)

saveRDS(object = egg_frate_bivariate_nohelpers, 
        file = "./models/bivariate_egg_volume_feeding_rate_no_helpers.RDS")

summary(egg_frate_bivariate_nohelpers)
plot(egg_frate_bivariate_nohelpers)

##
## model with helpers
set.seed(7)
egg_frate_bivariate_helpers <- MCMCglmm(cbind(scale(egg_volume_cm), scale(maternal_prov_rate)) ~ 
                                          trait - 1 +
                                          
                                          ## egg model
                                          at.level(trait, 1):poly(rainfall_egg,2)[,1] +
                                          at.level(trait, 1):poly(rainfall_egg,2)[,2] +
                                          at.level(trait, 1):egg_position + 
                                          at.level(trait, 1):clutch_size +
                                          at.level(trait, 1):temp_above_35_egg +
                                          at.level(trait, 1):female_helpers_egg +
                                          at.level(trait, 1):male_helpers_egg +
                                          
                                          ## prov model
                                          at.level(trait, 2):poly(rainfall_frate,2)[,1] +
                                          at.level(trait, 2):poly(rainfall_frate,2)[,2] +
                                          at.level(trait, 2):temp_above_35_frate +
                                          at.level(trait, 2):female_helpers_frate +
                                          at.level(trait, 2):male_helpers_frate +
                                          at.level(trait, 2):brood_size,
                                        
                                        random = 
                                          ~idh(trait):Group +  
                                          idh(trait):Season +
                                          us(trait):mother_ID + 
                                          us(trait):clutch_ID,
                                        rcov =~ idh(trait):units,
                                        family = c("gaussian","gaussian"),
                                        prior = prior_egg_frate_bivariate,
                                        nitt=510000,
                                        burnin=10000,
                                        thin=100,
                                        verbose = TRUE,
                                        data = data_model)
saveRDS(object = egg_frate_bivariate_helpers, 
        file = "./models/bivariate_egg_volume_feeding_rate_with_helpers.RDS")
summary(egg_frate_bivariate_helpers)
plot(egg_frate_bivariate_helpers)



