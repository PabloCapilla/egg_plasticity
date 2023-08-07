###
###
#' 
#' Script for:
#' Mothers front-load their investment to the egg stage when helped in a wild cooperative bird
#' Capilla-Lasheras et al. 
#' Preprint: https://doi.org/10.1101/2021.11.11.468195
#' 
#' Latest update: 2023/08/02
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
               lubridate,
               sjPlot, 
               gt,
               ggplot2, 
               ggExtra,
               ggpubr,
               extrafont,
               MuMIn, 
               lme4)
loadfonts()
source('./scripts/00_Functions/FUNCTION_d_centering.R')

# font_import() may be needed first to use font types in ggplot
# check with 'fonts()' that font are available

#####

##
##
##### egg volume dataset #####
##
##
data <- read.csv("./data/egg_volume_dataset.csv") 
head(data)

#####

##
##
##### Summary of initial sample sizes #####
##
##

## clutch_IDes
length(unique(data$clutch_ID))
median(data$clutch_size)

## groups
length(unique(data$Group))

## females
length(unique(data$mother_ID))

## clutch_ID per mother
median(
  {data %>%
      group_by(mother_ID, clutch_ID) %>%
      group_by(mother_ID) %>%
      summarise(n_clutch_IDes_female = n())}$n_clutch_IDes_female
)

range(
  {data %>%
      group_by(mother_ID, clutch_ID) %>%
      group_by(mother_ID) %>%
      summarise(n_clutch_IDes_female = n())}$n_clutch_IDes_female
)


## groups
length(unique(data$Group))

## mean, min, max egg volume per mother
data %>% 
  group_by(mother_ID) %>% 
  summarise(mean_per_mother = mean(egg_volume/1000)) %>% 
  summarise(mean_egg_vol = mean(mean_per_mother),
            min_egg_vol = min(mean_per_mother),
            max_egg_vol = max(mean_per_mother))

#####

##
##
##### 1 - Does egg volume predict egg mass? #####
##
##

##
## DATA 

# calculate egg age to do the analysis using eggs weighed the first day after laying
data$egg_age <- dmy(data$egg_weight_date) - dmy(data$laid_date)

# data for this first analysis
data_mass <- data %>% 
  mutate(egg_volume_cm = egg_volume/1000) %>% 
  filter(laid_date_error == 0) %>% # keep eggs for which lay date was accurately known
  filter(egg_age == 0)             # keep eggs that were measured the same say they were laid

# remove NA values from egg volume and egg mass columns
table(is.na(data_mass$egg_volume))
table(is.na(data_mass$egg_weight))

data_mass <- data_mass %>% 
  filter(!is.na(egg_weight))

##
## MODEL

# linear model to calculate whether egg volume predicts egg mass on the day of laying
m1a_egg_volume_egg_mass <- lmer(egg_weight ~ 
                                  egg_volume_cm +
                                (1|mother_ID), 
                             data = data_mass,
                             na.action = "na.fail")
summary(m1a_egg_volume_egg_mass)
drop1(m1a_egg_volume_egg_mass, test = "Chisq") # LR Test of egg volume effects on egg weight

##
## TABLE OF RESULTS
#
## Supplementary table

# table with model coefficients
tab_model(m1a_egg_volume_egg_mass,
          file="./tables/Table S3.doc",
          pred.labels = c("Intercept", 
                          html("Egg volume (cm3)")),
          string.ci = "CI (95%)",
          string.se = "SE",
          digits = 3,
          digits.p = 3,
          show.se = TRUE, 
          show.stat = FALSE,
          show.p = FALSE,
          show.est = TRUE,
          show.intercept = TRUE,
          rm.terms = NULL,
          show.re.var = FALSE,
          show.ngroups = FALSE,
          show.r2 = FALSE,
          show.obs = FALSE,
          ci.hyphen = ",",
          use.viewer = T)

#####

##
##
##### 2 - Does egg volume predict hatchling mass? #####
##
##

##
## DATA
data_weights <- read.csv("./data/egg_volume_hatchling_mass_dataset.csv")
data_weights$X <- NULL
head(data_weights)
data_weights <- data_weights %>% 
  mutate(egg_volume_cm = egg_volume/1000) 

# sample size no NAs
table(is.na(data_weights$egg_volume))
table(is.na(data_weights$hatchling_weight))

##
## MODEL

##
## cross-sectional
##
m2a_egg_volume_hatch_weight <- lmer(hatchling_weight ~ 
                                      egg_volume_cm +
                                      (1|mother_ID),
                                    data = data_weights, 
                                    na.action = "na.fail")
summary(m2a_egg_volume_hatch_weight)
drop1(m2a_egg_volume_hatch_weight, test = "Chisq") # LR Test of egg volume effects on hatchling weight

##
## TABLE

##
## TABLE OF RESULTS

# Supplementary table

# table with model coefficients
tab_model(m2a_egg_volume_hatch_weight,
          file="./tables/Table S2.doc",
          pred.labels = c("Intercept", 
                          html("Egg volume (cm3)")),
          string.ci = "CI (95%)",
          string.se = "SE",
          digits = 3,
          digits.p = 3,
          show.se = TRUE, 
          show.stat = FALSE,
          show.p = FALSE,
          show.est = TRUE,
          show.intercept = TRUE,
          rm.terms = NULL,
          show.re.var = FALSE,
          show.ngroups = FALSE,
          show.r2 = FALSE,
          show.obs = FALSE,
          ci.hyphen = ",",
          use.viewer = T)

##
## within mothers 
##

# partitioning egg volume in its within and among mother components (function in '00_FUNCTIONS')
data_weights <- data_weights %>% 
  mutate(egg_volume_cm = egg_volume/1000)
length(unique(data_weights$mother_ID)) # number of mothers in the analysis

data_weights <- d_centering(data = data_weights, 
                            centering_by = "mother_ID", 
                            cluster = "group_ID",
                            variable = c("egg_volume_cm"))
head(data_weights)

m2b_egg_volume_hatch_weight <- lmer(hatchling_weight ~ 
                                      cent_egg_volume_cm +
                                      mean_egg_volume_cm +
                                      (1|mother_ID),
                                    data = data_weights, 
                                    na.action = "na.fail")
summary(m2b_egg_volume_hatch_weight)
drop1(m2b_egg_volume_hatch_weight, test = "Chisq") # LR Test of within-mother egg volume effect on hatchling weight



##
## Supplementary table

# table with model coefficients
tab_model(m2b_egg_volume_hatch_weight,
          file="./tables/Table S3.doc",
          pred.labels = c("Intercept", 
                          html("&Delta;Egg volume"),
                          html("&mu;Egg volume")),
          string.ci = "CI (95%)",
          string.se = "SE",
          digits = 3,
          digits.p = 3,
          show.se = TRUE, 
          show.stat = FALSE,
          show.p = FALSE,
          show.est = TRUE,
          show.intercept = TRUE,
          rm.terms = NULL,
          show.re.var = FALSE,
          show.ngroups = FALSE,
          show.r2 = FALSE,
          show.obs = FALSE,
          ci.hyphen = ",",
          use.viewer = T)

#####

##
##
##### Plots for Figure 1 #####
##
##

##
## Plot Figure 1a
##
egg_variation_within_among_females <- ggplot(data =data, 
                                             aes(y = mean_volume/1000, 
                                                 x = cent_volume/1000)) +
  geom_point(color = "grey0", size = 2, alpha = 0.075) +
  #geom_density2d(color = "black", size = 1) +
  scale_x_continuous(limits = c(-0.75,0.75)) +
  scale_y_continuous(limits = c(2.5,4.5)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(family = "Arial", size = 12, vjust = -2),
    axis.title.y = element_text(family = "Arial", size = 12, vjust = +3),
    axis.text = element_text(family = "Arial", size = 10)) +
  labs(y = expression('Mean egg volume per mother (cm'^3*')'),
       x = expression(atop(paste(Delta, " Egg volume (cm"^"3", ")", sep = ""),
                           "(within-mother variation)")))

##
## Plot Figure 1b
##
egg_volume_egg_mass_plot <- ggplot(data = data_mass, 
                                   aes(x = egg_volume/1000, y = egg_weight)) +
  geom_point(size = 2, alpha = 0.15) +
  theme_bw() +
  labs(y = "Egg mass on laying date (g)", 
       x = expression(atop('Egg volume (cm'^3*')', " "))) + 
  theme(axis.title.x = element_text(family = "Arial", size = 12, vjust = -2),
        axis.title.y = element_text(family = "Arial", size = 12, vjust = +3),
        axis.text.x = element_text(family = "Arial", size = 10),
        axis.text.y = element_text(family = "Arial", size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = seq(2.0, 5, 0.5), 
                     limits = c(2,5)) +
  scale_x_continuous(limits = c(2.5, 5)) +
  stat_smooth(method = "lm", 
              color = "#a50f15", 
              fill = "#a50f15", 
              size = 1.5)

##
## Plot Figure 1c
##
egg_volume_hatchling_mass_plot <- ggplot(data = data_weights, 
                                         aes(x = egg_volume/1000, y = hatchling_weight)) +
  geom_point(size = 2, alpha = 0.15) +
  theme_bw() +
  labs(y = "Hatchling mass (g)", 
       x = expression(atop('Egg volume (cm'^3*')', " "))) + 
  theme(axis.title.x = element_text(family = "Arial", size = 12, vjust = -2),
        axis.title.y = element_text(family = "Arial", size = 12, vjust = +3),
        axis.text.x = element_text(family = "Arial", size = 10),
        axis.text.y = element_text(family = "Arial", size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = seq(0.5, 5.5, 0.5)) +
  scale_x_continuous(limits = c(2.5, 5)) +
  stat_smooth(method = "lm", 
              color = "#a50f15", 
              fill = "#a50f15", 
              size = 1.5)


##
## Create Panel for Figure 1 and save
##
Figure1 <- ggarrange(egg_variation_within_among_females,
                     egg_volume_egg_mass_plot,
                     egg_volume_hatchling_mass_plot,
                     align = "hv",
                     vjust = 9.5, 
                     hjust = c(-8.5, -7.7, -8.5),
                     labels = c("a", "b", "c"),
                     font.label = list(size = 30, color = "black", face = "bold", family = "Arial"),
                     ncol = 3, nrow = 1)

ggsave(filename = "./plots/Figure 1.png", 
       plot = Figure1,
       device = "png",
       units = "mm",
       width = 180, 
       height = 95)

#####
