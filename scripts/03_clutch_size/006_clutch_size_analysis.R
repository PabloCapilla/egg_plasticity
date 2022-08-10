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
#' Analysis of clutch size Results presented in Table S8 and S9
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
##### clutch size dataset #####
##
##
data <- read.csv("./data/clutch_size_dataset.csv") 
head(data)

##
##
##### Summary of sample sizes #####
##
##
# sample size in this model
length(unique(data$clutch_ID)) # 344 clutches
length(unique(data$mother_ID)) # 66 mothers
length(unique(data$Group_ID))  # 37 social groups
range(                         # range of number of clutches per mother
  {data %>% 
      group_by(mother_ID, clutch_ID) %>%
      group_by(mother_ID) %>%
      summarise(n_clutch_IDes_female = n())}$n_clutch_IDes_female
)
median(                        # median number of clutches per female
  {data %>% 
      group_by(mother_ID, clutch_ID) %>%
      group_by(mother_ID) %>%
      summarise(n_clutch_IDes_female = n())}$n_clutch_IDes_female
)


## 
##
#### Fifth model - zero-truncated model of clutch size - population level #####
##
##
data$mother_ID <- as.factor(data$mother_ID)
data$Season <- as.factor(data$Season)
data$Group_ID <- as.factor(data$Group_ID)

# model using glmmADMB package
model5_pop_level <- glmmADMB::glmmadmb(clutch_size ~ 
                                         scale(clutch_order) + # improve speed 
                                         female_helpers +
                                         male_helpers + 
                                         (1|mother_ID) + 
                                         (1|Season) + 
                                         (1|Group_ID), 
                                       data=data,
                                       family="truncpoiss")
summary(model5_pop_level)
drop1(model5_pop_level, test = "Chisq")


## 
##
#### Fifth model - zero-truncated model of clutch size - partitioning #####
##
##
model5_partitioning <- glmmADMB::glmmadmb(clutch_size ~ 
                                            scale(clutch_order) + # improve speed 
                                            cent_female_helpers +
                                            cent_male_helpers + 
                                            mean_female_helpers +
                                            mean_male_helpers + 
                                            (1|mother_ID) + 
                                            (1|Season) + 
                                            (1|Group_ID), 
                                          data=data,
                                          family="truncpoiss")
summary(model5_partitioning)
drop1(model5_partitioning, test = "Chisq")

## 
##
#### Fifth model - zero-truncated model of clutch size - testing diffs in within and among-mother effects #####
##
##
model5_test_slopes <- glmmADMB::glmmadmb(clutch_size ~ 
                                           scale(clutch_order) + # improve speed 
                                           male_helpers +
                                           mean_male_helpers +   # this is explicitly testing differences between within and between-mother slopes
                                           female_helpers +
                                           mean_female_helpers + # this is explicitly testing differences between within and between-mother slopes
                                           (1|mother_ID) + 
                                           (1|Season) + 
                                           (1|Group_ID), 
                                         data=data,
                                         family="truncpoiss")
summary(model5_test_slopes)
drop1(model5_test_slopes, test = "Chisq")



## 
##
#### Fifth model - zero-truncated model of clutch size - population level - AIC table #####
##
##

## AIC calculation
clutch_size_model5_full_aic_table <- dredge(model5_pop_level, rank = "AIC", trace = 3)
clutch_size_model5_d6_subset <- subset(clutch_size_model5_full_aic_table, delta < 6) 


##
## 
## create table with AIC results
names(clutch_size_model5_d6_subset)
number_variables <- 4
col_names <- names(clutch_size_model5_d6_subset)[1:number_variables]

##
##
##### Code to generate Table S5 #####
##
##

## add each model to the table
# list to store data
list_models_table_S5 <- as.list(NA)
for(m in 1:nrow(clutch_size_model5_d6_subset)){
  
  # template to store data
  table_S5_template <- data.frame(coefficient = names(clutch_size_model5_d6_subset)[1:number_variables]) 
  
  # add model coefficients
  model_coef <- data.frame(coefTable(get.models(clutch_size_model5_d6_subset, m)[[1]])) %>% 
    mutate(coefficient = rownames(coefTable(get.models(clutch_size_model5_d6_subset, m)[[1]]))) %>% 
    rename(estimate = Estimate, 
           SE = Std..Error)
  table_S5_00 <- left_join(x = table_S5_template, 
                           y = model_coef %>% select(-df), 
                           by = "coefficient")
  
  ## put table data in right format
  table_S5_01 <- table_S5_00 %>% 
    pivot_wider(names_from = coefficient, values_from = c(estimate, SE)) %>% 
    mutate(k = clutch_size_model5_d6_subset$df[m],
           AIC = clutch_size_model5_d6_subset$AIC[m],
           delta = clutch_size_model5_d6_subset$delta[m]) %>% 
    relocate(`Intercept` = `estimate_(Intercept)`,
             `Intercept SE` = `SE_(Intercept)`,
             `Number of helping females` = `estimate_female_helpers`,
             `Number of helping females SE` = `SE_female_helpers`,
             `Number of helping males` = `estimate_male_helpers`,
             `Number of helping males SE` = `SE_male_helpers`,
             `Clutch order` = `estimate_scale(clutch_order)`,            
             `Clutch order SE` = `SE_scale(clutch_order)`,
             k = k, 
             AIC = AIC, 
             delta = delta)
  
  # save data per model
  list_models_table_S5[[m]] <- table_S5_01           
  
}

# combine data from each model in one table
table_S5_data <- rbindlist(list_models_table_S5)

# remove columns with all NA (remove variables that don't appear in any model in the table)
table_S5_data <- table_S5_data %>%
  select_if(~ !all(is.na(.)))

# form table
table_S5_clean <- table_S5_data %>% 
  gt() %>% 
  fmt_number(columns = c(1:(ncol(table_S5_data)-3), (ncol(table_S5_data)-1):(ncol(table_S5_data))),
             decimals = 2) %>% 
  fmt_number(columns = ncol(table_S5_data)-2,
             decimals = 0) %>% 
  cols_merge_uncert(col_val = `Intercept`, col_uncert = `Intercept SE`) %>% 
  cols_merge_uncert(col_val = `Number of helping females`, col_uncert =`Number of helping females SE`) %>% 
  cols_merge_uncert(col_val = `Number of helping males`, col_uncert =`Number of helping males SE`) %>% 
  cols_merge_uncert(col_val = `Clutch order`, col_uncert =`Clutch order SE`) %>% 
  fmt_missing(columns = c(1:(ncol(table_S5_data)-3), (ncol(table_S5_data)-1):(ncol(table_S5_data))), 
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

# TABLE S8
table_S5_clean

# save Table S8 (saved in html, then imported in docx to include in manuscript)
table_S5_clean %>%
  gtsave(filename = "./tables/Table S5.html")



## 
##
#### Fifth model - zero-truncated model of clutch size - partitioning - AIC table #####
##
##

## AIC calculation
clutch_size_model5_part_aic_table <- dredge(model5_partitioning, rank = "AIC", trace = 3)
clutch_size_model5_part_d6_subset <- subset(clutch_size_model5_part_aic_table, delta < 6) 


##
## 
## create table with AIC results
names(clutch_size_model5_part_d6_subset)
number_variables <- 6
col_names <- names(clutch_size_model5_part_d6_subset)[1:number_variables]

##
##
##### Code to generate Table S6 #####
##
##

## add each model to the table
# list to store data
list_models_table_S6 <- as.list(NA)
for(m in 1:nrow(clutch_size_model5_part_d6_subset)){
  
  # template to store data
  table_S6_template <- data.frame(coefficient = names(clutch_size_model5_part_d6_subset)[1:number_variables]) 
  
  # add model coefficients
  model_coef <- data.frame(coefTable(get.models(clutch_size_model5_part_d6_subset, m)[[1]])) %>% 
    mutate(coefficient = rownames(coefTable(get.models(clutch_size_model5_part_d6_subset, m)[[1]]))) %>% 
    rename(estimate = Estimate, 
           SE = Std..Error)
  table_S6_00 <- left_join(x = table_S6_template, 
                           y = model_coef %>% select(-df), 
                           by = "coefficient")
  
  ## put table data in right format
  table_S6_01 <- table_S6_00 %>% 
    pivot_wider(names_from = coefficient, values_from = c(estimate, SE)) %>% 
    mutate(k = clutch_size_model5_part_d6_subset$df[m],
           AIC = clutch_size_model5_part_d6_subset$AIC[m],
           delta = clutch_size_model5_part_d6_subset$delta[m]) %>% 
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
             `Clutch order` = `estimate_scale(clutch_order)`,            
             `Clutch order SE` = `SE_scale(clutch_order)`,
             k = k, 
             AIC = AIC, 
             delta = delta)
  
  # save data per model
  list_models_table_S6[[m]] <- table_S6_01           
  
}

# combine data from each model in one table
table_S6_data <- rbindlist(list_models_table_S6)

# remove columns with all NA (remove variables that don't appear in any model in the table)
table_S6_data <- table_S6_data %>%
  select_if(~ !all(is.na(.)))

# form table
table_S6_clean <- table_S6_data %>% 
  gt() %>% 
  fmt_number(columns = c(1:(ncol(table_S6_data)-3), (ncol(table_S6_data)-1):(ncol(table_S6_data))),
             decimals = 2) %>% 
  fmt_number(columns = ncol(table_S6_data)-2,
             decimals = 0) %>% 
  cols_merge_uncert(col_val = `Intercept`, col_uncert = `Intercept SE`) %>% 
  cols_merge_uncert(col_val = `d Number of helping females`, col_uncert =`d Number of helping females SE`) %>% 
  cols_merge_uncert(col_val = `d Number of helping males`, col_uncert =`d Number of helping males SE`) %>% 
  cols_merge_uncert(col_val = `mean Number of helping females`, col_uncert =`mean Number of helping females SE`) %>% 
  cols_merge_uncert(col_val = `mean Number of helping males`, col_uncert =`mean Number of helping males SE`) %>% 
  cols_merge_uncert(col_val = `Clutch order`, col_uncert =`Clutch order SE`) %>% 
  fmt_missing(columns = c(1:(ncol(table_S6_data)-3), (ncol(table_S6_data)-1):(ncol(table_S6_data))), 
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
table_S6_clean

# save Table S8 (saved in html, then imported in docx to include in manuscript)
table_S6_clean %>%
  gtsave(filename = "./tables/Table S6.html")



