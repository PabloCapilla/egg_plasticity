###
###
#' 
#' Script for:
#' Mothers front-load their investment to the egg stage when helped in a wild cooperative bird
#' Capilla-Lasheras et al. 
#' Preprint: https://doi.org/10.1101/2021.11.11.468195
#' 
#' Latest update: 2022/08/10
#' 
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list=ls())

##
##
##### Script aim: #####
#' Analysis of total clutches per year.
#' 
##
##

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

##
##
##### total clutches per season and mother dataset #####
##
##
data <- read.csv("./data/total_clutches_per_year.csv") 

## 
##
#### Sixth model - clutch number model - population-level #####
##
##
model6_clutches <- glmer(clutches ~ 
                           female_helpers +
                           male_helpers +
                           rainfall +
                           (1|group_ID) +
                           (1|mother_ID),
                         family = "poisson",
                         na.action="na.fail",
                         data = data)
summary(model6_clutches)
drop1(model6_clutches, test = "Chisq")


## 
##
#### Sixth model - clutch number model - partitioning #####
##
##
model6_clutches_partitioning <- glmer(clutches ~ 
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
summary(model6_clutches_partitioning)
drop1(model6_clutches_partitioning, test = "Chisq")


## 
##
#### Sixth model - clutch number model - testing differences in within- and among-mother slopes #####
##
##
model6_test_slopes <- glmer(clutches ~ 
                              male_helpers +
                              mean_male_helpers +   # this is explicitly testing differences between within and between-mother slopes
                              female_helpers +
                              mean_female_helpers + # this is explicitly testing differences between within and between-mother slopes
                              rainfall +
                              (1|group_ID) +
                              (1|mother_ID),
                            family = "poisson",
                            na.action="na.fail",
                            data = data)
summary(model6_test_slopes)
drop1(model6_test_slopes, test = "Chisq")


## 
##
#### Sixth model - clutch number model - population-level - AIC table #####
##
##

clutches_model6_full_aic_table <- dredge(model6_clutches, rank = "AIC", trace = 3)
clutches_model6_d6_subset <- subset(clutches_model6_full_aic_table, delta < 6) 

##
## create table with AIC results
names(clutches_model6_d6_subset)
number_variables <- 4
col_names <- names(clutches_model6_d6_subset)[1:number_variables]

##
##
##### Code to generate Table S7 #####
##
##

## add each model to the table
# list to store data
list_models_table_S7 <- as.list(NA)
for(m in 1:nrow(clutches_model6_d6_subset)){
  
  # template to store data
  table_S7_template <- data.frame(coefficient = names(clutches_model6_d6_subset)[1:number_variables]) 
  
  # add model coefficients
  model_coef <- data.frame(coefTable(get.models(clutches_model6_d6_subset, m)[[1]])) %>% 
    mutate(coefficient = rownames(coefTable(get.models(clutches_model6_d6_subset, m)[[1]]))) %>% 
    rename(estimate = Estimate, 
           SE = Std..Error)
  table_S7_00 <- left_join(x = table_S7_template, 
                           y = model_coef %>% select(-df), 
                           by = "coefficient")
  
  ## put table data in right format
  table_S7_01 <- table_S7_00 %>% 
    pivot_wider(names_from = coefficient, values_from = c(estimate, SE)) %>% 
    mutate(k = clutches_model6_d6_subset$df[m],
           AIC = clutches_model6_d6_subset$AIC[m],
           delta = clutches_model6_d6_subset$delta[m]) %>% 
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
  list_models_table_S7[[m]] <- table_S7_01           
  
}

# combine data from each model in one table
table_S7_data <- rbindlist(list_models_table_S7)

# remove columns with all NA (remove variables that don't appear in any model in the table)
table_S7_data <- table_S7_data %>%
  select_if(~ !all(is.na(.)))

# form table
table_S7_clean <- table_S7_data %>% 
  gt() %>% 
  fmt_number(columns = c(1:(ncol(table_S7_data)-3), (ncol(table_S7_data)-1):(ncol(table_S7_data))),
             decimals = 3) %>% 
  fmt_number(columns = ncol(table_S7_data)-2,
             decimals = 0) %>% 
  cols_merge_uncert(col_val = `Intercept`, col_uncert = `Intercept SE`) %>% 
  cols_merge_uncert(col_val = `Number of helping females`, col_uncert =`Number of helping females SE`) %>% 
  cols_merge_uncert(col_val = `Number of helping males`, col_uncert =`Number of helping males SE`) %>% 
  cols_merge_uncert(col_val = `Rainfall`, col_uncert =`Rainfall SE`) %>% 
  fmt_missing(columns = c(1:(ncol(table_S7_data)-3), (ncol(table_S7_data)-1):(ncol(table_S7_data))), 
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
table_S7_clean

# save Table S10 (saved in html, then imported in docx to include in manuscript)
table_S7_clean %>%
  gtsave(filename = "./tables/Table S7.html")



## 
##
#### Sixth model - clutch number model - partitioning - AIC table #####
##
##
clutches_model6_part_aic_table <- dredge(model6_clutches_partitioning, rank = "AIC", trace = 3)
clutches_model6_part_d6_subset <- subset(clutches_model6_part_aic_table, delta < 6) 


##
## 
## create table with AIC results
names(clutches_model6_part_d6_subset)
number_variables <- 6
col_names <- names(clutches_model6_part_d6_subset)[1:number_variables]

##
##
##### Code to generate Table S8 #####
##
##

## add each model to the table
# list to store data
list_models_table_S8 <- as.list(NA)
for(m in 1:nrow(clutches_model6_part_d6_subset)){
  
  # template to store data
  table_S8_template <- data.frame(coefficient = names(clutches_model6_part_d6_subset)[1:number_variables]) 
  
  # add model coefficients
  model_coef <- data.frame(coefTable(get.models(clutches_model6_part_d6_subset, m)[[1]])) %>% 
    mutate(coefficient = rownames(coefTable(get.models(clutches_model6_part_d6_subset, m)[[1]]))) %>% 
    rename(estimate = Estimate, 
           SE = Std..Error)
  table_S8_00 <- left_join(x = table_S8_template, 
                            y = model_coef %>% select(-df), 
                            by = "coefficient")
  
  ## put table data in right format
  table_S8_01 <- table_S8_00 %>% 
    pivot_wider(names_from = coefficient, values_from = c(estimate, SE)) %>% 
    mutate(k = clutches_model6_part_d6_subset$df[m],
           AIC = clutches_model6_part_d6_subset$AIC[m],
           delta = clutches_model6_part_d6_subset$delta[m]) %>% 
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
  list_models_table_S8[[m]] <- table_S8_01           
  
}

# combine data from each model in one table
table_S8_data <- rbindlist(list_models_table_S8)

# remove columns with all NA (remove variables that don't appear in any model in the table)
table_S8_data <- table_S8_data %>%
  select_if(~ !all(is.na(.)))

# form table
table_S8_clean <- table_S8_data %>% 
  gt() %>% 
  fmt_number(columns = c(1:(ncol(table_S8_data)-3), (ncol(table_S8_data)-1):(ncol(table_S8_data))),
             decimals = 3) %>% 
  fmt_number(columns = ncol(table_S8_data)-2,
             decimals = 0) %>% 
  cols_merge_uncert(col_val = `Intercept`, col_uncert = `Intercept SE`) %>% 
  cols_merge_uncert(col_val = `d Number of helping females`, col_uncert =`d Number of helping females SE`) %>% 
  cols_merge_uncert(col_val = `d Number of helping males`, col_uncert =`d Number of helping males SE`) %>% 
  cols_merge_uncert(col_val = `mean Number of helping females`, col_uncert =`mean Number of helping females SE`) %>% 
  cols_merge_uncert(col_val = `mean Number of helping males`, col_uncert =`mean Number of helping males SE`) %>% 
  cols_merge_uncert(col_val = `Rainfall`, col_uncert =`Rainfall SE`) %>% 
  fmt_missing(columns = c(1:(ncol(table_S8_data)-3), (ncol(table_S8_data)-1):(ncol(table_S8_data))), 
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

# TABLE S8
table_S8_clean

# save Table S8 (saved in html, then imported in docx to include in manuscript)
table_S8_clean %>%
  gtsave(filename = "./tables/Table S8.html")





