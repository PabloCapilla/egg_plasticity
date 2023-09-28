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
               gt,
               gtsummary,
               lme4)

source('./scripts/00_Functions/FUNCTION_drop1_output.R')


##
##
##### data ######
##
##
data <- read.csv("./data/maternal_condition.csv")
head(data)

#####

##
##
##### Summary of sample sizes #####
##
##
nrow(data) # clutches
length(unique(data$mother_ID)) # mothers
length(unique(data$group_ID))    # groups
length(unique(data$clutch_ID)) 

##
##
##### Model Supplementary Materials C #####
##
##
model_maternal_conditions <- lmer(body_index ~ 
                                    prev_female_number +
                                    prev_male_number +
                                    clutch_order +
                                    (1|season) +
                                    (1|mother_ID) +
                                    (1|group_ID), 
                                  data =data, 
                                  REML = F,
                                  na.action = "na.fail")
summary(model_maternal_conditions)
drop1(model_maternal_conditions, test = "Chisq")

##
## table

## base table
base_table <- model_maternal_conditions %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `prev_female_number` = 'Number of female helpers in previous breeding attempt',
                   `prev_male_number` = 'Number of male helpers in previous breeding attempt',
                   `clutch_order` = "Clutch order"), 
                 estimate_fun = ~ style_number(.x, digits = 3))

## add features
final_table <- base_table %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=model_maternal_conditions) %>% 
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

gtsave(data = final_table_formatted, "./tables/Table S20.html")

