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
               data.table,
               sjPlot, 
               gt,
               gtsummary,
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
##### total clutches per season and mother dataset #####
##
##
data <- read.csv("./data/total_clutches_per_year.csv") 

data %>% 
  group_by(season) %>% 
  summarise(sd_rainfall = sd(rainfall)) # no variation in rainfall within breeding seasons

#####

## 
##
####clutch number model - population-level #####
##
##
model6_clutches <- glmer(clutches ~ 
                           female_helpers +
                           male_helpers +
                           rainfall +
                           (1|season) +
                           (1|group_ID) +
                           (1|mother_ID),
                         family = "poisson",
                         na.action="na.fail",
                         data = data)
summary(model6_clutches)
drop1(model6_clutches, test = "Chisq")
anova(model6_clutches, update(model6_clutches, .~.-(1|season)), test = 'Chisq')

##
## TABLE

## base table
base_table <- model6_clutches %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `rainfall` = "Rainfall", 
                   `female_helpers` = "Number of female helpers",
                   `male_helpers` = "Number of male helpers"), 
                 estimate_fun = ~ style_number(.x, digits = 3))


## add features
final_table <- base_table %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=model6_clutches) %>% 
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

gtsave(data = final_table_formatted, "./tables/Table S18.html")

#####

## 
##
#### clutch number model - partitioning #####
##
##
model6_clutches_partitioning <- glmer(clutches ~ 
                                        cent_female_helpers +
                                        cent_male_helpers +
                                        mean_female_helpers +
                                        mean_male_helpers +
                                        rainfall +
                                        (1|season) +
                                        (1|group_ID) +
                                        (1|mother_ID),
                                      family = "poisson",
                                      na.action="na.fail",
                                      data = data)
summary(model6_clutches_partitioning)
drop1(model6_clutches_partitioning, test = "Chisq")
anova(model6_clutches_partitioning, update(model6_clutches_partitioning, .~.-(1|season)), test = 'Chisq')




##
## TABLE

## base table
base_table <- model6_clutches_partitioning %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `rainfall` = "Rainfall", 
                   `cent_female_helpers` = html("&Delta; Number of helping females"),
                   `mean_female_helpers` = html("&mu; Number of helping females"),
                   `cent_male_helpers` = html("&Delta; Number of helping males"),
                   `mean_male_helpers` = html("&mu; Number of helping males")), 
                 estimate_fun = ~ style_number(.x, digits = 3))


## add features
final_table <- base_table %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=model6_clutches_partitioning) %>% 
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

gtsave(data = final_table_formatted, "./tables/Table S19.html")

#####
