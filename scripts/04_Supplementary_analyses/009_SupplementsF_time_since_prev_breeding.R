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
               lme4)

source('./scripts/00_Functions/FUNCTION_drop1_output.R')

##
##
##### egg volume dataset #####
##
##
data <- read.csv("./data/egg_volume_dataset.csv") 
data <- data %>% 
  mutate(egg_volume_cm = egg_volume/1000)
head(data)
nrow(data)


##
##
##### Population-level egg volume model - including time since last breeding event and number of helpers in previous breeding even #####
##
##
data <- data %>%
  filter(!is.na(time_since_prev_ba)) %>%
  filter(fate_prev_ba   == "SUCCESS")   # dataset includes only successful previous breeding events

length(unique(data$clutch_ID)) # n clutches
length(unique(data$mother_ID)) # n mothers
length(unique(data$Group)) # n groups 

summary(data$time_since_prev_ba) # range time since
mean(data$time_since_prev_ba)    # mean time since

##
##
##### Model for effects of time since last breeding attempt on egg volume #####
##
##
model_carryover_effects <- lmer(egg_volume_cm ~ 
                                  
                                  # predictors testing for carry over effects
                                  time_since_prev_ba : female_helpers +
                                  time_since_prev_ba : male_helpers +
                                  
                                  # environmental indexes
                                  poly(rainfall,2)[,1] +
                                  poly(rainfall,2)[,2] +
                                  temp_above_35 +
                                  
                                  
                                  # carryover effect predictors main effect
                                  time_since_prev_ba +
                                  
                                  
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
summary(model_carryover_effects)
drop1(model_carryover_effects, test = "Chisq")

##
## effect of 'time since last breeding attempt' in a model without interactions
summary(update(model_carryover_effects, . ~ . - 
                 time_since_prev_ba : female_helpers -
                 time_since_prev_ba : male_helpers))

drop1(update(model_carryover_effects, . ~ . - 
               time_since_prev_ba : female_helpers -
               time_since_prev_ba : male_helpers),
      test = "Chisq")

single_eff_model <- update(model_carryover_effects, . ~ . - 
                             time_since_prev_ba : female_helpers -
                             time_since_prev_ba : male_helpers)
summary(single_eff_model)
##
## table

## base table
base_table <- single_eff_model %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `poly(rainfall, 2)[, 1]` = html("Rainfall<sup>1</sup>"),
                   `poly(rainfall, 2)[, 2]` = html("Rainfall<sup>2</sup>"),
                   `temp_above_35` = "Heat waves",
                   `female_helpers` = "Number of female helpers",
                   `male_helpers` = "Number of male helpers",
                   `scale(clutch_size)` = "Clutch size",
                   `scale(egg_position)` = "Egg position",
                   `time_since_prev_ba` = 'Time since previous breeding attempt'), 
                 estimate_fun = ~ style_number(.x, digits = 3))

## add features
final_table <- base_table %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=single_eff_model) %>% 
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

gtsave(data = final_table_formatted, "./tables/Table S21.html")

