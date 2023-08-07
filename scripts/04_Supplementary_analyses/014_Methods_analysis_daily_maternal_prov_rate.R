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
##### Maternal provisioning rate dataset #####
##
##
data <- read.csv("./data/maternal_provisioning_dataset_daily.csv")
head(data)
table(data$season)

# sample size in this model
nrow(data)            
length(unique(data$clutch_ID)) # 124 clutches
length(unique(data$mother_ID)) # 50 mothers
length(unique(data$group_ID))  # 34 social groups

##
##
##### 0 - correlation consecutive provisioning rates #####
##
##
corre_frate <- data %>% 
  select(date, clutch_ID, maternal_prov_rate) %>% 
  group_by(clutch_ID) %>% 
  arrange((date)) %>% 
  mutate(prev_frate = lag(maternal_prov_rate, 1)) %>% 
  filter(!is.na(prev_frate))

cor.test(corre_frate$maternal_prov_rate, corre_frate$prev_frate)
summary(lm(scale(maternal_prov_rate) ~ scale(prev_frate), 
           data = corre_frate))
plot(corre_frate$maternal_prov_rate, corre_frate$prev_frate)

#####


##
##
##### 1 - population-level maternal provisioning rates #####
##
##

data$obs_ID <- 1:nrow(data)

##
## full model - without interactions
model3_main_effects_model <- glmer (maternal_n_feeds ~ 
                                      # rainfall and heat waves
                                      poly(rainfall,2)[,1] +
                                      poly(rainfall,2)[,2] +
                                      scale(temp_above_35) +
                                      
                                      # other main effects
                                      female_helpers +
                                      male_helpers +
                                      scale(brood_size) +
                                      (1|obs_ID) +
                                      (1|clutch_ID) +
                                      (1|group_ID) +
                                      (1|mother_ID) +
                                      (1|season),
                                    offset = log(session_duration),
                                    family = "poisson",
                                    data = data,
                                    control = glmerControl(optimizer = 'bobyqa'),
                                    na.action = "na.fail")

summary(model3_main_effects_model) # summary of model
drop1(update(object = model3_main_effects_model), 
      test = "Chisq") # Likelihood-ratio test for each predictor

confint(model3_main_effects_model)

# formal test of differences between female helper and male helper effects within-mother from the best model containing both
car::linearHypothesis(model3_main_effects_model, "female_helpers - male_helpers = 0")


##
## table of results

## base table
base_table <- model3_main_effects_model %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `poly(rainfall, 2)[, 1]` = html("Rainfall<sup>1</sup>"),
                   `poly(rainfall, 2)[, 2]` = html("Rainfall<sup>2</sup>"),
                   `scale(temp_above_35)` = "Heat waves",
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
                        y = drop1_output(x=model3_main_effects_model) %>% 
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

gtsave(data = final_table_formatted, "./tables/Table S31.html")
