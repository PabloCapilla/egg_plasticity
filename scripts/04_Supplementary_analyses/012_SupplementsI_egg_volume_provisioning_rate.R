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
               sjPlot,
               extrafont,
               MCMCglmm,
               brms,
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

#####

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
##### data formatting #####
##
##

# variable format management
data$resp1_sub <- ifelse(is.na(data$egg_volume_cm), 0, 1)
table(data$resp1_sub)
data$resp2_sub <- ifelse(is.na(data$maternal_prov_rate), 0, 1)
table(data$resp2_sub)

# standarised quadratic rainfall effects 
data$rainfall_egg1[!is.na(data$rainfall_egg)] <- poly(data$rainfall_egg[!is.na(data$rainfall_egg)],2)[,1]
data$rainfall_egg2[!is.na(data$rainfall_egg)] <- poly(data$rainfall_egg[!is.na(data$rainfall_egg)],2)[,2]
data$rainfall_frate1[!is.na(data$rainfall_frate)] <- poly(data$rainfall_frate[!is.na(data$rainfall_frate)],2)[,1]
data$rainfall_frate2[!is.na(data$rainfall_frate)] <- poly(data$rainfall_frate[!is.na(data$rainfall_frate)],2)[,2]

#####

##
##
##### bivariate model no helpers #####
##
##
bf_egg <- bf(scale(egg_volume_cm)|subset(resp1_sub) ~ 
               rainfall_egg2 +
               rainfall_egg1 +
               scale(temp_above_35_egg) + 
               scale(egg_position) + 
               scale(clutch_size) +
               (1|a|Group) +
               (1|b|Season) +
               (1|c|mother_ID) + 
               (1|d|clutch_ID), 
             decomp = "QR")

bf_prov <- bf(scale(maternal_prov_rate)|subset(resp2_sub) ~ 
                rainfall_frate2 +
                rainfall_frate1 +
                scale(temp_above_35_frate) + 
                scale(brood_size) +
                (1|z|Group) +
                (1|x|Season) +
                (1|c|mother_ID) + 
                (1|d|clutch_ID),
              decomp = "QR")

model_prior <- prior(constant(0.01),class = "sigma",resp = "scalematernalprovrate") +
  prior(normal(0,100), class = "b", resp = "scalematernalprovrate") +
  prior(normal(0,100), class = "b", resp = "scaleeggvolumecm")

bivar_no_helpers <- brm(bf_egg + bf_prov + set_rescor(F), 
                        data = data, 
                        chains = 4, 
                        cores = 4, 
                        iter = 50000, 
                        warmup = 25000, 
                        thin = 10,
                        seed = 7,
                        prior = model_prior,
                        normalize = FALSE,
                        backend = 'cmdstanr')
##
## save model results
#saveRDS(object = bivar_no_helpers, file = './models/SupplementsI_egg_volume_prov_rate_brms_no_helpers.RDS')
bivar_no_helpers <- readRDS(file = './models/SupplementsI_egg_volume_prov_rate_brms_no_helpers.RDS')

## results
summary(bivar_no_helpers)
plot(bivar_no_helpers)
varcor <- VarCorr(bivar_no_helpers)
(varcor$clutch_ID$cov)



tab_model(bivar_no_helpers, 
          show.ci50 = TRUE, 
          show.se = TRUE, 
          show.re.var = TRUE,
          file="./tables/Table S25.doc",
          string.ci = "CI (95%)",
          string.se = "SE",
          digits = 3,
          use.viewer = T)
#####

##
##
##### bivariate model with helpers #####
##
##
bf_egg_helpers <- bf(scale(egg_volume_cm)|subset(resp1_sub) ~ 
                       rainfall_egg2 +
                       rainfall_egg1 +
                       scale(temp_above_35_egg) + 
                       scale(egg_position) + 
                       scale(clutch_size) +
                       male_helpers_egg +
                       female_helpers_egg +
                       (1|a|Group) +
                       (1|b|Season) +
                       (1|c|mother_ID) + 
                       (1|d|clutch_ID), 
                     decomp = "QR")

bf_prov_helpers <- bf(scale(maternal_prov_rate)|subset(resp2_sub) ~ 
                        rainfall_frate2 +
                        rainfall_frate1 +
                        scale(temp_above_35_frate) + 
                        scale(brood_size) +
                        male_helpers_frate +
                        female_helpers_frate +
                        (1|z|Group) +
                        (1|x|Season) +
                        (1|c|mother_ID) + 
                        (1|d|clutch_ID),
                      decomp = "QR")

model_prior <- prior(constant(0.01),class = "sigma",resp = "scalematernalprovrate") +
  prior(normal(0,100), class = "b", resp = "scalematernalprovrate") +
  prior(normal(0,100), class = "b", resp = "scaleeggvolumecm")

bivar_with_helpers <- brm(bf_egg_helpers + bf_prov_helpers + set_rescor(F), 
                          data = data, 
                          chains = 4, 
                          cores = 4, 
                          iter = 50000, 
                          warmup = 25000, 
                          thin = 10,
                          seed = 7,
                          prior = model_prior,
                          normalize = FALSE,
                          backend = 'cmdstanr')
##
## save model results
#saveRDS(object = bivar_with_helpers, file = './models/SupplementsI_egg_volume_prov_rate_brms_with_helpers.RDS')
bivar_with_helpers <- readRDS(file = './models/SupplementsI_egg_volume_prov_rate_brms_with_helpers.RDS')

## results
summary(bivar_with_helpers)
plot(bivar_with_helpers)
varcor <- VarCorr(bivar_with_helpers)
(varcor$clutch_ID$cov)

tab_model(bivar_with_helpers, 
          show.ci50 = TRUE, 
          show.se = TRUE, 
          show.re.var = TRUE,
          file="./tables/Table S24.doc",
          string.ci = "CI (95%)",
          string.se = "SE",
          digits = 3,
          use.viewer = T)
#####

s