###
###
#' 
#' Script for:
#' Mothers front-load their investment to the egg stage when helped in a wild cooperative bird
#' Capilla-Lasheras et al. 
#' Preprint: https://doi.org/10.1101/2021.11.11.468195
#' 
#' Latest update: 2023/08/07
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

##
##
##### egg volume dataset #####
##
##
data <- read.csv("./data/egg_volume_dataset906.csv") 


## filter reformat
data <- data %>% 
  mutate(egg_volume_cm = egg_volume/1000) 

#####

##
##
##### Full model presented in main text with extended dataset #####
##
##
model1_full_model <- lmer(egg_volume_cm ~ 
                            poly(rainfall,2)[,1] +
                            poly(rainfall,2)[,2] +
                            scale(temp_above_35) +

                            # interactions
                            cent_female_helpers : scale(clutch_size) +
                            cent_male_helpers : scale(clutch_size) +
                            mean_female_helpers : scale(clutch_size) +
                            mean_male_helpers : scale(clutch_size) +
                            
                            
                            # other main effects
                            cent_female_helpers +
                            cent_male_helpers +
                            mean_female_helpers +
                            mean_male_helpers +
                            scale(clutch_size) +
                            
                            # random effects
                            (1|Group) +
                            (1|clutch_ID) +
                            (1|mother_ID) +
                            (1|Season),
                          data = data,
                          na.action = "na.fail",
                          REML = F)
summary(model1_full_model) # summary of model

##
## table of results

## base table
base_table <- model1_full_model %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `poly(rainfall, 2)[, 1]` = html("Rainfall<sup>1</sup>"),
                   `poly(rainfall, 2)[, 2]` = html("Rainfall<sup>2</sup>"),
                   `scale(temp_above_35)` = "Heat waves",
                   `cent_female_helpers` = "delta Number of female helpers",
                   `mean_female_helpers` = "mean Number of female helpers",
                   `cent_male_helpers` = "delta Number of male helpers",
                   `mean_male_helpers` = "mean Number of male helpers",
                   `scale(clutch_size)` = "Clutch size",
                   `cent_female_helpers:scale(clutch_size)` = "delta Number of female helpers x Clutch size",
                   `scale(clutch_size):mean_female_helpers` = "Mean Number of female helpers x Clutch size",
                   `scale(clutch_size):cent_male_helpers` = "delta Number of male helpers x Clutch size",
                   `scale(clutch_size):mean_male_helpers` = "Mean Number of female helpers x Clutch size"),
                 estimate_fun = ~ style_number(.x, digits = 3))

## add features
final_table <- base_table %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=model1_full_model) %>% 
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
  modify_header(label ~ "**Predictors**") %>% 
  modify_header(std.error ~ "**SE**") %>%
  modify_header(estimate ~ "**Estimates**") %>%
  modify_header(df ~ "**df**") %>% 
  modify_header(Chisq ~ html("<b>&chi;<sup>2</sup></b>")) %>% 
  as_gt() %>% 
  opt_footnote_marks(marks = "LETTERS")

gtsave(data = final_table_formatted, "./tables/Table S26.html")

#####

##
##
##### 1-egg clutches ####
##
##
onegg_clutch <- data %>% 
  filter(clutch_size == 1)

mode_1egg <- lmer(egg_volume_cm ~ 
                            poly(rainfall,2)[,1] +
                            poly(rainfall,2)[,2] +
                            scale(temp_above_35) +

                            # other main effects
                            cent_female_helpers +
                            cent_male_helpers +
                            mean_female_helpers +
                            mean_male_helpers +

                            # random effects
                            (1|Group) +
                            (1|mother_ID) +
                            (1|Season),
                          data = onegg_clutch,
                          na.action = "na.fail",
                          REML = F)
summary(mode_1egg) # summary of model

##
## table of results

## base table
base_table <- mode_1egg %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `poly(rainfall, 2)[, 1]` = html("Rainfall<sup>1</sup>"),
                   `poly(rainfall, 2)[, 2]` = html("Rainfall<sup>2</sup>"),
                   `scale(temp_above_35)` = "Heat waves",
                   `cent_female_helpers` = "delta Number of female helpers",
                   `mean_female_helpers` = "mean Number of female helpers",
                   `cent_male_helpers` = "delta Number of male helpers",
                   `mean_male_helpers` = "mean Number of male helpers"),
                 estimate_fun = ~ style_number(.x, digits = 3))

## add features
final_table <- base_table %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=mode_1egg) %>% 
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
  modify_header(label ~ "**Predictors**") %>% 
  modify_header(std.error ~ "**SE**") %>%
  modify_header(estimate ~ "**Estimates**") %>%
  modify_header(df ~ "**df**") %>% 
  modify_header(Chisq ~ html("<b>&chi;<sup>2</sup></b>")) %>% 
  as_gt() %>% 
  opt_footnote_marks(marks = "LETTERS")

gtsave(data = final_table_formatted, "./tables/Table S27.html")

#####

##
##
##### 3-egg clutches ####
##
##
threegg_clutch <- data %>% 
  filter(clutch_size == 3)

model_3egg <- lmer(egg_volume_cm ~ 
                    poly(rainfall,2)[,1] +
                    poly(rainfall,2)[,2] +
                    scale(temp_above_35) +
                    
                    # other main effects
                    cent_female_helpers +
                    cent_male_helpers +
                    mean_female_helpers +
                    mean_male_helpers +
                    
                    # random effects
                    (1|Group) +
                    (1|clutch_ID) +
                    (1|mother_ID) +
                    (1|Season),
                  data = threegg_clutch,
                  na.action = "na.fail",
                  REML = F)
summary(model_3egg) # summary of model

##
## table of results

## base table
base_table <- model_3egg %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `poly(rainfall, 2)[, 1]` = html("Rainfall<sup>1</sup>"),
                   `poly(rainfall, 2)[, 2]` = html("Rainfall<sup>2</sup>"),
                   `scale(temp_above_35)` = "Heat waves",
                   `cent_female_helpers` = "delta Number of female helpers",
                   `mean_female_helpers` = "mean Number of female helpers",
                   `cent_male_helpers` = "delta Number of male helpers",
                   `mean_male_helpers` = "mean Number of male helpers"),
                 estimate_fun = ~ style_number(.x, digits = 3))

## add features
final_table <- base_table %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=model_3egg) %>% 
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
  modify_header(label ~ "**Predictors**") %>% 
  modify_header(std.error ~ "**SE**") %>%
  modify_header(estimate ~ "**Estimates**") %>%
  modify_header(df ~ "**df**") %>% 
  modify_header(Chisq ~ html("<b>&chi;<sup>2</sup></b>")) %>% 
  as_gt() %>% 
  opt_footnote_marks(marks = "LETTERS")

gtsave(data = final_table_formatted, "./tables/Table S28.html")

#####

##
##
##### 2-egg clutches ####
##
##
twoegg_clutch <- data %>% 
  filter(clutch_size == 2)

model_2egg <- lmer(egg_volume_cm ~ 
                     poly(rainfall,2)[,1] +
                     poly(rainfall,2)[,2] +
                     scale(temp_above_35) +
                     
                     # other main effects
                     cent_female_helpers +
                     cent_male_helpers +
                     mean_female_helpers +
                     mean_male_helpers +
                     
                     # random effects
                     (1|Group) +
                     (1|clutch_ID) +
                     (1|mother_ID) +
                     (1|Season),
                   data = twoegg_clutch,
                   na.action = "na.fail",
                   REML = F)
summary(model_2egg) # summary of model

##
## table of results

## base table
base_table <- model_2egg %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `poly(rainfall, 2)[, 1]` = html("Rainfall<sup>1</sup>"),
                   `poly(rainfall, 2)[, 2]` = html("Rainfall<sup>2</sup>"),
                   `scale(temp_above_35)` = "Heat waves",
                   `cent_female_helpers` = "delta Number of female helpers",
                   `mean_female_helpers` = "mean Number of female helpers",
                   `cent_male_helpers` = "delta Number of male helpers",
                   `mean_male_helpers` = "mean Number of male helpers"),
                 estimate_fun = ~ style_number(.x, digits = 3))

## add features
final_table <- base_table %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=model_2egg) %>% 
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
  modify_header(label ~ "**Predictors**") %>% 
  modify_header(std.error ~ "**SE**") %>%
  modify_header(estimate ~ "**Estimates**") %>%
  modify_header(df ~ "**df**") %>% 
  modify_header(Chisq ~ html("<b>&chi;<sup>2</sup></b>")) %>% 
  as_gt() %>% 
  opt_footnote_marks(marks = "LETTERS")

gtsave(data = final_table_formatted, "./tables/Table S29.html")

#####

##
##
##### Full model with whole dataset, including clutch size as a factor #####
##
##
model_int_full <- lmer(egg_volume_cm ~ 
                            
                            ##
                            ## no environmental effects as sliding window not optimised for whole dataset
                            # rainfall and heat waves
                            poly(rainfall,2)[,1] +
                            poly(rainfall,2)[,2] +
                            scale(temp_above_35) +
                            
                            # interactions
                            cent_female_helpers : as.factor(clutch_size) +
                            cent_male_helpers : as.factor(clutch_size) +
                            mean_female_helpers : as.factor(clutch_size) +
                            mean_male_helpers : as.factor(clutch_size) +
                            
                            
                            # other main effects
                            cent_female_helpers +
                            cent_male_helpers +
                            mean_female_helpers +
                            mean_male_helpers +
                            as.factor(clutch_size) +
                            
                            # random effects
                            (1|Group) +
                            (1|clutch_ID) +
                            (1|mother_ID) +
                            (1|Season),
                          data = data,
                          na.action = "na.fail",
                          REML = F)
summary(model_int_full) # summary of model
drop1(model_int_full, test = "Chisq")

##
## table of results

## base table
base_table <- model_int_full %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `poly(rainfall, 2)[, 1]` = html("Rainfall<sup>1</sup>"),
                   `poly(rainfall, 2)[, 2]` = html("Rainfall<sup>2</sup>"),
                   `scale(temp_above_35)` = "Heat waves",
                   `cent_female_helpers` = "delta Number of female helpers",
                   `mean_female_helpers` = "mean Number of female helpers",
                   `cent_male_helpers` = "delta Number of male helpers",
                   `mean_male_helpers` = "mean Number of male helpers",
                   `as.factor(clutch_size)` = "Clutch size",
                   `cent_female_helpers:as.factor(clutch_size)` = "delta Number of female helpers x Clutch size",
                   `as.factor(clutch_size):mean_female_helpers` = "Mean Number of female helpers x Clutch size",
                   `as.factor(clutch_size):cent_male_helpers` = "delta Number of male helpers x Clutch size",
                   `as.factor(clutch_size):mean_male_helpers` = "Mean Number of male helpers x Clutch size"),
                 estimate_fun = ~ style_number(.x, digits = 3))

## add features
final_table <- base_table %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=model_int_full) %>% 
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
  modify_header(label ~ "**Predictors**") %>% 
  modify_header(std.error ~ "**SE**") %>%
  modify_header(estimate ~ "**Estimates**") %>%
  modify_header(df ~ "**df**") %>% 
  modify_header(Chisq ~ html("<b>&chi;<sup>2</sup></b>")) %>% 
  as_gt() %>% 
  opt_footnote_marks(marks = "LETTERS")

gtsave(data = final_table_formatted, "./tables/Table S30.html")

#####

