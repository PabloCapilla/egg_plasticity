###
###
#' 
#' Script for:
#' Mothers front-load their investment to the egg stage when helped in a wild cooperative bird
#' Capilla-Lasheras et al. 
#' Preprint: https://doi.org/10.1101/2021.11.11.468195
#' 
#' Latest update: 2022/03/31
#' 
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list=ls())

##
##
##### Script aim: #####
#' Analysis of total clutches per year. Results presented in Table S10 and S11
#' 
##
##

##
##### libraries #####
##
pacman::p_load(dplyr, 
               tidyr, 
               sjPlot, 
               gt,
               ggplot2, 
               extrafont,
               MuMIn, 
               lme4)
loadfonts()
# font_import() may be needed first to use font types in ggplot
# check with 'fonts()' that font are available

##
##
##### total clutches per season and mother dataset #####
##
##
data <- read.csv("./data/total_clutches_per_year.csv") 

## 
##
#### First model - population-level clutch number model #####
##
##
clutches_model1 <- glmer(clutches ~ 
                               female_helpers +
                               male_helpers +
                               rainfall +
                               (1|group_ID) +
                               (1|mother_ID),
                             family = "poisson",
                             na.action="na.fail",
                             data = data)
summary(clutches_model1)

# model selection based on AIC
clutches_model1_full_aic_table <- dredge(clutches_model1, rank = "AIC", trace = 3)
clutches_model1_d6_subset <- subset(clutches_model1_full_aic_table, delta < 6) 



##
## 
## create table with AIC results
names(clutches_model1_d6_subset)
number_variables <- 4
col_names <- names(clutches_model1_d6_subset)[1:number_variables]

##
##
##### Code to generate Table S10 #####
##
##

## add each model to the table
# list to store data
list_models_table_S10 <- as.list(NA)
for(m in 1:nrow(clutches_model1_d6_subset)){
  
  # template to store data
  table_S10_template <- data.frame(coefficient = names(clutches_model1_d6_subset)[1:number_variables]) 
  
  # add model coefficients
  model_coef <- data.frame(coefTable(get.models(clutches_model1_d6_subset, m)[[1]])) %>% 
    mutate(coefficient = rownames(coefTable(get.models(clutches_model1_d6_subset, m)[[1]]))) %>% 
    rename(estimate = Estimate, 
           SE = Std..Error)
  table_S10_00 <- left_join(x = table_S10_template, 
                           y = model_coef %>% select(-df), 
                           by = "coefficient")
  
  ## put table data in right format
  table_S10_01 <- table_S10_00 %>% 
    pivot_wider(names_from = coefficient, values_from = c(estimate, SE)) %>% 
    mutate(k = clutches_model1_d6_subset$df[m],
           AIC = clutches_model1_d6_subset$AIC[m],
           delta = clutches_model1_d6_subset$delta[m]) %>% 
    relocate(`Intercept` = `estimate_(Intercept)`,
             `Intercept SE` = `SE_(Intercept)`,
             `Number of helping females` = `estimate_female_helpers`,
             `Number of helping females SE` = `SE_female_helpers`,
             `Number of helping males` = `estimate_male_helpers`,
             `Number of helping males SE` = `SE_male_helpers`,
             `Rainfall` = `estimate_rainfall`,            
             `Rainfall SE` = `SE_rainfall`,
             k = k, 
             AIC = AIC, 
             delta = delta)
  
  # save data per model
  list_models_table_S10[[m]] <- table_S10_01           
  
}

# combine data from each model in one table
table_S10_data <- rbindlist(list_models_table_S10)

# remove columns with all NA (remove variables that don't appear in any model in the table)
table_S10_data <- table_S10_data %>%
  select_if(~ !all(is.na(.)))

# form table
table_S10_clean <- table_S10_data %>% 
  gt() %>% 
  fmt_number(columns = c(1:(ncol(table_S10_data)-3), (ncol(table_S10_data)-1):(ncol(table_S10_data))),
             decimals = 3) %>% 
  fmt_number(columns = ncol(table_S10_data)-2,
             decimals = 0) %>% 
  cols_merge_uncert(col_val = `Intercept`, col_uncert = `Intercept SE`) %>% 
  cols_merge_uncert(col_val = `Number of helping females`, col_uncert =`Number of helping females SE`) %>% 
  cols_merge_uncert(col_val = `Number of helping males`, col_uncert =`Number of helping males SE`) %>% 
  cols_merge_uncert(col_val = `Rainfall`, col_uncert =`Rainfall SE`) %>% 
  fmt_missing(columns = c(1:(ncol(table_S10_data)-3), (ncol(table_S10_data)-1):(ncol(table_S10_data))), 
              missing_text = " ") %>% 
  cols_label(`delta` = html("&Delta;AIC")) %>% 
  cols_align(align = c("center")) %>% 
  tab_options(table.font.names = "Arial",
              table.font.size = 11,
              table.font.color = "black",
              column_labels.font.weight = "bold",
              table_body.hlines.width = 1,
              table_body.hlines.color = "black",
              table_body.border.bottom.width = 2,
              table_body.border.bottom.color = "black",
              column_labels.border.top.width = 2,
              column_labels.border.top.color = "black",
              column_labels.border.bottom.width = 2,
              column_labels.border.bottom.color = "black")

# TABLE S10
table_S10_clean

# save Table S10 (saved in html, then imported in docx to include in manuscript)
table_S10_clean %>%
  gtsave(filename = "./tables/Table S10.html")



## 
##
#### Second model - partitioning variation in female and male helper number #####
##
##
clutches_model2 <- glmer(clutches ~ 
                           cent_female_helpers +
                           cent_male_helpers +
                           mean_female_helpers +
                           mean_male_helpers +
                           rainfall +
                           (1|group_ID) +
                           (1|mother_ID),
                         family = "poisson",
                         na.action="na.fail",
                         data = data)
summary(clutches_model2)

# model selection based on AIC
clutches_model2_full_aic_table <- dredge(clutches_model2, rank = "AIC", trace = 3)
clutches_model2_d6_subset <- subset(clutches_model2_full_aic_table, delta < 6) 



##
## 
## create table with AIC results
names(clutches_model2_d6_subset)
number_variables <- 6
col_names <- names(clutches_model2_d6_subset)[1:number_variables]

##
##
##### Code to generate Table S11 #####
##
##

## add each model to the table
# list to store data
list_models_table_S11 <- as.list(NA)
for(m in 1:nrow(clutches_model2_d6_subset)){
  
  # template to store data
  table_S11_template <- data.frame(coefficient = names(clutches_model2_d6_subset)[1:number_variables]) 
  
  # add model coefficients
  model_coef <- data.frame(coefTable(get.models(clutches_model2_d6_subset, m)[[1]])) %>% 
    mutate(coefficient = rownames(coefTable(get.models(clutches_model2_d6_subset, m)[[1]]))) %>% 
    rename(estimate = Estimate, 
           SE = Std..Error)
  table_S11_00 <- left_join(x = table_S11_template, 
                            y = model_coef %>% select(-df), 
                            by = "coefficient")
  
  ## put table data in right format
  table_S11_01 <- table_S11_00 %>% 
    pivot_wider(names_from = coefficient, values_from = c(estimate, SE)) %>% 
    mutate(k = clutches_model2_d6_subset$df[m],
           AIC = clutches_model2_d6_subset$AIC[m],
           delta = clutches_model2_d6_subset$delta[m]) %>% 
    relocate(`Intercept` = `estimate_(Intercept)`,
             `Intercept SE` = `SE_(Intercept)`,
             `d Number of helping females` = `estimate_cent_female_helpers`,
             `d Number of helping females SE` = `SE_cent_female_helpers`,
             `mean Number of helping females` = `estimate_mean_female_helpers`,
             `mean Number of helping females SE` = `SE_mean_female_helpers`,
             `d Number of helping males` = `estimate_cent_male_helpers`,
             `d Number of helping males SE` = `SE_cent_male_helpers`,
             `mean Number of helping males` = `estimate_mean_male_helpers`,
             `mean Number of helping males SE` = `SE_mean_male_helpers`,
             `Rainfall` = `estimate_rainfall`,            
             `Rainfall SE` = `SE_rainfall`,
             k = k, 
             AIC = AIC, 
             delta = delta)
  
  # save data per model
  list_models_table_S11[[m]] <- table_S11_01           
  
}

# combine data from each model in one table
table_S11_data <- rbindlist(list_models_table_S11)

# remove columns with all NA (remove variables that don't appear in any model in the table)
table_S11_data <- table_S11_data %>%
  select_if(~ !all(is.na(.)))

# form table
table_S11_clean <- table_S11_data %>% 
  gt() %>% 
  fmt_number(columns = c(1:(ncol(table_S11_data)-3), (ncol(table_S11_data)-1):(ncol(table_S11_data))),
             decimals = 3) %>% 
  fmt_number(columns = ncol(table_S11_data)-2,
             decimals = 0) %>% 
  cols_merge_uncert(col_val = `Intercept`, col_uncert = `Intercept SE`) %>% 
  cols_merge_uncert(col_val = `d Number of helping females`, col_uncert =`d Number of helping females SE`) %>% 
  cols_merge_uncert(col_val = `d Number of helping males`, col_uncert =`d Number of helping males SE`) %>% 
  cols_merge_uncert(col_val = `mean Number of helping females`, col_uncert =`mean Number of helping females SE`) %>% 
  cols_merge_uncert(col_val = `mean Number of helping males`, col_uncert =`mean Number of helping males SE`) %>% 
  cols_merge_uncert(col_val = `Rainfall`, col_uncert =`Rainfall SE`) %>% 
  fmt_missing(columns = c(1:(ncol(table_S11_data)-3), (ncol(table_S11_data)-1):(ncol(table_S11_data))), 
              missing_text = " ") %>% 
  cols_label(`d Number of helping females` = html("&Delta; Number of helping females"),
             `d Number of helping males` = html("&Delta; Number of helping males"),
             `mean Number of helping females` = html("&mu; Number of helping females"),
             `mean Number of helping males` = html("&mu; Number of helping males"),
             `delta` = html("&Delta;AIC")) %>% 
  cols_align(align = c("center")) %>% 
  tab_options(table.font.names = "Arial",
              table.font.size = 11,
              table.font.color = "black",
              column_labels.font.weight = "bold",
              table_body.hlines.width = 1,
              table_body.hlines.color = "black",
              table_body.border.bottom.width = 2,
              table_body.border.bottom.color = "black",
              column_labels.border.top.width = 2,
              column_labels.border.top.color = "black",
              column_labels.border.bottom.width = 2,
              column_labels.border.bottom.color = "black")

# TABLE S11
table_S11_clean

# save Table S11 (saved in html, then imported in docx to include in manuscript)
table_S11_clean %>%
  gtsave(filename = "./tables/Table S11.html")





