###
###
#' 
#' Script for:
#' Mothers front-load their investment to the egg stage when helped in a wild cooperative bird
#' Capilla-Lasheras et al. 
#' Preprint: https://doi.org/10.1101/2021.11.11.468195
#' 
#' Latest update: 2022/08/09
#' 
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list=ls())

##
##
##### Script aim: #####
#' Introductory analysis of egg volume: within and among mother variation in egg volume and 
#' correlations between egg volume, egg mass and hatchling weight. Results presented in Results section 1 and Figure 1
#' 
##
##

##
##### libraries #####
##
pacman::p_load(dplyr, 
               WeaverTools,
               tidyr, 
               ggplot2,
               ggpubr,
               lubridate,
               extrafont,
               MuMIn, 
               lme4)
loadfonts()
# font_import() may be needed first to use font types in ggplot
# check with 'fonts()' that font are available

##
##
##### egg volume dataset #####
##
##
data <- read.csv("./data/egg_volume_dataset.csv") 
head(data)

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

##
##
##### 1 - Does egg volume predict egg mass? #####
##
##

# calculate egg age to do the analysis using eggs weighed the first day after laying
data$egg_age <- dmy(data$egg_weight_date) - dmy(data$laid_date)

# data for this first analysis
data_mass <- data %>% 
  filter(laid_date_error == 0) %>% # keep eggs for which lay date was accurately known
  filter(egg_age == 0)             # keep eggs that were measured the same say they were laid

# remove NA values from egg volume and egg mass columns
table(is.na(data_mass$egg_volume))
table(is.na(data_mass$egg_weight))
data_mass <- data_mass %>% 
  filter(!is.na(egg_weight))

# linear model to calculate whether egg volume predicts egg mass on the day of laying
m1a_egg_volume_egg_mass <- lmer(egg_weight ~ 
                                egg_volume +
                                (1|mother_ID), 
                             data = data_mass,
                             na.action = "na.fail")
summary(m1a_egg_volume_egg_mass)
drop1(m1a_egg_volume_egg_mass, test = "Chisq") # LR Test of egg volume effects on egg weight

##
## within mothers 
##

# partitioning egg volume in its within and among mother components
# (d_centering in WeaverTools R package available on GitHub)
data_mass <- d_centering(data = data_mass, 
                            centering_by = "mother_ID", 
                            cluster = "Group",
                            variable = c("egg_volume"))
head(data_mass)

m1b_egg_volume_egg_mass <- lmer(egg_weight ~ 
                                cent_egg_volume +
                                mean_egg_volume  +
                                (1|mother_ID),
                              data = data_mass, 
                              na.action = "na.fail")
summary(m1b_egg_volume_egg_mass)
drop1(m1b_egg_volume_egg_mass, test = "Chisq") # LR Test of within-mother egg volume effect on egg weight

length(unique(data_mass$mother_ID)) # number of mothers in the analysis

##
## Test of differences in slopes between and within mothers (applying Equation 3 in van de Pol & Wright 2009)
m1c_test_slopes <- lmer(egg_weight ~ 
                        egg_volume + 
                        mean_egg_volume+  # this is explicitly testing differences between within and between-mother slopes
                        (1|mother_ID),       
                      data = data_mass, 
                      na.action = "na.fail")
summary(m1c_test_slopes)
drop1(m1c_test_slopes, test = "Chisq") # LR Test of differences between within and among-mother slopes


##
##
##### 2 - Does egg volume predict hatchling mass? #####
##
##
data_weights <- read.csv("./data/egg_volume_hatchling_mass_dataset.csv")
data_weights$X <- NULL
head(data_weights)

# sample size no NAs
table(is.na(data_weights$egg_volume))
table(is.na(data_weights$hatchling_weight))

##
## cross-sectional
##
m2a_egg_volume_hatch_weight <- lmer(hatchling_weight ~ 
                                egg_volume +
                                (1|mother_ID),
         data = data_weights, 
         na.action = "na.fail")
summary(m2a_egg_volume_hatch_weight)
drop1(m2a_egg_volume_hatch_weight, test = "Chisq") # LR Test of egg volume effects on hatchling weight

##
## within mothers 
##

# partitioning egg volume in its within and among mother components
# (d_centering in WeaverTools R package available on GitHub)
data_weights <- d_centering(data = data_weights, 
                            centering_by = "mother_ID", 
                            cluster = "group_ID",
                            variable = c("egg_volume"))
head(data_weights)

m2b_egg_volume_hatch_weight <- lmer(hatchling_weight ~ 
                                      cent_egg_volume +
                                      mean_egg_volume +
                                      (1|mother_ID),
                                    data = data_weights, 
                                    na.action = "na.fail")
summary(m2b_egg_volume_hatch_weight)
drop1(m2b_egg_volume_hatch_weight, test = "Chisq") # LR Test of within-mother egg volume effect on hatchling weight

length(unique(data_weights$mother_ID)) # number of mothers in the analysis


##
## Test of differences in slopes between and within mothers (applying Equation 3 in van de Pol & Wright 2009)
m2c_test_slopes <- lmer(hatchling_weight ~ 
                        egg_volume + 
                        mean_egg_volume +  # this is explicitly testing differences between within and between-mother slopes
                      (1|mother_ID),
                      data = data_weights, 
                      na.action = "na.fail")
summary(m2c_test_slopes)
drop1(m2c_test_slopes, test = "Chisq")

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

