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
               sjPlot, 
               gt,
               ggplot2, 
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
##### Maternal provisioning rate dataset #####
##
##
data <- read.csv("./data/maternal_provisioning_dataset.csv")
head(data)
table(data$season)

#####

##
##
##### population-level maternal provisioning rates #####
##
##

##
## full model without rainfall
model3_full_model <- lmer (maternal_prov_rate ~ 
                             scale(temp_above_35) +
                             female_helpers +
                             male_helpers +
                             scale(brood_size) +
                             (1|group_ID) +
                             (1|mother_ID) +
                             (1|season),
                           data = data,
                           na.action = "na.fail",
                           REML = FALSE) # boundary (singular) warning due to group_ID variance = 0. Removing this random effect removes the warning and produces the same estimates
summary(model3_full_model) # summary of model
drop1(model3_full_model, test = "Chisq")


##
## table

## base table
base_table <- model3_full_model %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `scale(temp_above_35)` = 'Heat waves',
                   `female_helpers` = "Number of female helpers",
                   `male_helpers` = "Number of male helpers",
                   `scale(brood_size)` = "Brood size"), 
                 estimate_fun = ~ style_number(.x, digits = 3))

## add features
final_table <- base_table %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=model3_full_model) %>% 
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

gtsave(data = final_table_formatted, "./tables/Table S16.html")

#####

##
##
##### partitioning maternal provisioning rates #####
##
##


# full model partition within and among-mother variation helper number
model4_main_effects_model <- lmer (maternal_prov_rate ~ 
                                     # rainfall and heat waves
                                     scale(temp_above_35) +
                                     
                                     # other main effects
                                     cent_female_helpers +
                                     mean_female_helpers +
                                     cent_male_helpers +
                                     mean_male_helpers +
                                     scale(brood_size) +
                                     (1|group_ID) +
                                     (1|mother_ID) +
                                     (1|season),
                                   data = data,
                                   na.action = "na.fail",
                                   REML = FALSE) 

##
## table

## base table
base_table <- model4_main_effects_model %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `scale(temp_above_35)` = 'Heat waves',
                   `cent_female_helpers` = html("&Delta; Number of helping females"),
                   `mean_female_helpers` = html("&mu; Number of helping females"),
                   `cent_male_helpers` = html("&Delta; Number of helping males"),
                   `mean_male_helpers` = html("&mu; Number of helping males"),
                   `scale(brood_size)` = "Brood size"), 
                 estimate_fun = ~ style_number(.x, digits = 3))

## add features
final_table <- base_table %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=model4_main_effects_model) %>% 
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

gtsave(data = final_table_formatted, "./tables/Table S17.html")
