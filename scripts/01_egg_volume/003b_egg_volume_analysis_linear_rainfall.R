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
               sjPlot, 
               gt,
               gtsummary,
               ggplot2, 
               MetBrewer,
               ggpubr,
               extrafont,
               MuMIn, 
               lme4)
loadfonts()
# font_import() may be needed first to use font types in ggplot
# check with 'fonts()' that font are available

source('./scripts/00_Functions/FUNCTION_drop1_output.R')

#####

##
##
##### egg volume dataset #####
##
##
data <- read.csv("./data/egg_volume_dataset.csv") 
head(data)
nrow(data)

#####

##
##
##### population-level egg volume model rainfall below peak #####
##
##
data <- data %>% 
  mutate(egg_volume_cm = egg_volume/1000) %>% 
  filter(rainfall < 60)

##
## full model - without interactions
model1_main_effects_model <- lmer(egg_volume_cm ~ 
                                    
                                    # rainfall and heat waves
                                    scale(rainfall) +
                                    scale(temp_above_35) +
                                    
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
summary(model1_main_effects_model) # summary of model
drop1(update(object = model1_main_effects_model, REML = F), 
      test = "Chisq") # Likelihood-ratio test for each predictor

##
## table

## base table
base_table <- model1_main_effects_model %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `scale(rainfall)` = 'Rainfall',
                   `scale(temp_above_35)` = 'Heat waves',
                   `female_helpers` = "Number of female helpers",
                   `male_helpers` = "Number of male helpers",
                   `scale(clutch_size)` = "Clutch size", 
                   `scale(egg_position)` = "Egg position"), 
                 estimate_fun = ~ style_number(.x, digits = 3))


## add features
final_table <- base_table %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=model1_main_effects_model) %>% 
                          dplyr::select(variable = term, Chisq=statistic, df),
                        by = "variable")
    output$df <- ifelse(output$row_type == "label",  output$df, NA)
    output$Chisq <- ifelse(output$row_type == "label",  output$Chisq, NA)
    return(output)
  })

final_table_formatted <- final_table %>% 
  modify_fmt_fun(c(Chisq) ~ function(x) style_number(x, digits = 2)) %>%
  modify_fmt_fun(c(std.error) ~ function(x) style_number(x, digits = 3)) %>%
  modify_fmt_fun(c(p.value) ~ function(x) style_number(x, digits = 3)) %>%
  modify_table_body(~.x %>% dplyr::relocate(p.value, .after = df)) %>% 
  modify_header(label ~ "**Fixed effect**") %>% 
  modify_header(std.error ~ "**SE**") %>%
  modify_header(estimate ~ "**Estimate**") %>%
  modify_header(df ~ "**df**") %>% 
  modify_header(Chisq ~ html("<b>&chi;<sup>2</sup></b>")) %>% 
  as_gt() %>% 
  opt_footnote_marks(marks = "LETTERS")

gtsave(data = final_table_formatted, "./tables/Table S12.html")

#####

##
##
##### partitioning egg volume model rainfall below peak #####
##
##
model2_partitioning_below_peak <- lmer(egg_volume_cm ~ 
                                         scale(rainfall) +
                                         scale(temp_above_35) +
                                         
                                         # main effects
                                         cent_male_helpers +
                                         mean_male_helpers +
                                         cent_female_helpers +
                                         mean_female_helpers +
                                         scale(egg_position) +
                                         scale(clutch_size) +
                                         
                                         #random effects
                                         (1|Group) +
                                         (1|mother_ID) +
                                         (1|Season) + 
                                         (1|clutch_ID),
                                       data = data,
                                       na.action = "na.fail",
                                       REML = F)

summary(model2_partitioning_below_peak) # summary of model
drop1(update(object = model2_partitioning_below_peak, REML = F), 
      test = "Chisq") # Likelihood-ratio test for each predictor

##
## table

## base table
base_table <- model2_partitioning_below_peak %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `scale(rainfall)` = 'Rainfall',
                   `scale(temp_above_35)` = 'Heat waves',
                   `cent_female_helpers` = html("&Delta; Number of helping females"),
                   `mean_female_helpers` = html("&mu; Number of helping females"),
                   `cent_male_helpers` = html("&Delta; Number of helping males"),
                   `mean_male_helpers` = html("&mu; Number of helping males"),
                   `scale(clutch_size)` = "Clutch size", 
                   `scale(egg_position)` = "Egg position"), 
                 estimate_fun = ~ style_number(.x, digits = 3))


## add features
final_table <- base_table %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=model2_partitioning_below_peak) %>% 
                          dplyr::select(variable = term, Chisq=statistic, df),
                        by = "variable")
    output$df <- ifelse(output$row_type == "label",  output$df, NA)
    output$Chisq <- ifelse(output$row_type == "label",  output$Chisq, NA)
    return(output)
  })

final_table_formatted <- final_table %>% 
  modify_fmt_fun(c(Chisq) ~ function(x) style_number(x, digits = 2)) %>%
  modify_fmt_fun(c(std.error) ~ function(x) style_number(x, digits = 3)) %>%
  modify_fmt_fun(c(p.value) ~ function(x) style_number(x, digits = 3)) %>%
  modify_table_body(~.x %>% dplyr::relocate(p.value, .after = df)) %>% 
  modify_header(label ~ "**Fixed effect**") %>% 
  modify_header(std.error ~ "**SE**") %>%
  modify_header(estimate ~ "**Estimate**") %>%
  modify_header(df ~ "**df**") %>% 
  modify_header(Chisq ~ html("<b>&chi;<sup>2</sup></b>")) %>% 
  as_gt() %>% 
  opt_footnote_marks(marks = "LETTERS")

gtsave(data = final_table_formatted, "./tables/Table S13.html")
