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
               extrafont)
loadfonts()


##
##
####  data ####
##
##
data <- readRDS("./data/helper_number_predictability.RDS")
head(data)

##
##
##### Correlation laying and feeding #####
##
##

## females
cor.test(data$female_helpers_feeding - 1, data$female_helpers_laying - 1)

m1 <- lm(female_helpers_feeding ~ female_helpers_laying, 
         data = data, na.action = "na.fail")
summary(m1)

## males
cor.test(data$male_helpers_feeding - 1, data$male_helpers_laying - 1)

m2 <- lm(male_helpers_feeding ~ male_helpers_laying, 
         data = data,
         na.action = "na.fail")
summary(m2)


##
##
##### Plot #####
##
##
data_plot00 <- data %>% 
  group_by(female_helpers_laying) %>%
  summarise(mean_feeding = mean(female_helpers_feeding),
            sd_feeding = sd(female_helpers_feeding)) %>%
  mutate(sex = "females")
names(data_plot00) <- c("laying", names(data_plot00)[2:4] )

data_plot01 <- data %>% 
  group_by(male_helpers_laying) %>%
  summarise(mean_feeding = mean(male_helpers_feeding),
            sd_feeding = sd(male_helpers_feeding)) %>%
  mutate(sex = "males")
names(data_plot01) <- c("laying", names(data_plot01)[2:4] )

# data for plot
data_plot <- rbind(data_plot00, data_plot01)

group_cor_plot <- ggplot(data = data_plot, 
                         aes(y = mean_feeding, x = laying, color = sex)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = mean_feeding - sd_feeding,
                    ymax = mean_feeding + sd_feeding),
                width = 0,
                position = position_dodge(width = 0.5),
                size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", size =0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.margin = unit(c(5,5,30,20), "points"),
        legend.position = c(0.2,0.8),
        legend.title = element_blank(),
        legend.text = element_text(family = "Arial", size = 10),
        axis.title.y =  element_text(family = "Arial", size = 10, vjust = 3),
        axis.title.x =  element_text(family = "Arial", size = 10, vjust = -2),
        axis.text = element_text(family = "Arial", size = 10)) +
  scale_x_continuous(breaks = 0:5,
                     labels = 0:5) +
  scale_y_continuous(breaks = 0:5,
                     labels = 0:5) +
  labs(x = expression(atop("Helpers on laying date", " ")),
       y = "Helpers during rearing period") +
  scale_color_manual(values = c("black", "grey"))

  
ggsave(group_cor_plot, 
       filename = "./plots/Figure S1.png", 
       units = "mm",
       height = 89, 
       width = 89, 
       device = "png")
