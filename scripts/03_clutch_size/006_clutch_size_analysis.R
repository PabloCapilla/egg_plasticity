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
#' Analysis of clutch size Results presented in Table S8 and S9
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
#### First model - population-level clutch size model #####
##
##
clutch_size_model1 <- lmer(clutch_size ~ 
                             clutch_order +
                             female_helpers +
                             male_helpers + 
                             (1|mother_ID) + 
                             (1|Season) + 
                             (1|Group_ID),
                           REML = FALSE,
                           na.action = "na.fail",
                           data = data)
summary(clutch_size_model1)
hist(residuals(clutch_size_model1)) # check distribution of model residuals

  
# model selection based on AIC
clutch_size_model1_full_aic_table <- dredge(clutch_size_model1, rank = "AIC", trace = 3)
clutch_size_model1_d6_subset <- subset(clutch_size_model1_full_aic_table, delta < 6) 


##
## 
## create table with AIC results
names(model1_d6_subset)
number_variables <- 4
col_names <- names(clutch_size_model1_d6_subset)[1:number_variables]

##
##
##### Code to generate Table S8 #####
##
##

## add each model to the table
# list to store data
list_models_table_S8 <- as.list(NA)
for(m in 1:nrow(clutch_size_model1_d6_subset)){
  
  # template to store data
  table_S8_template <- data.frame(coefficient = names(clutch_size_model1_d6_subset)[1:number_variables]) 
  
  # add model coefficients
  model_coef <- data.frame(coefTable(get.models(clutch_size_model1_d6_subset, m)[[1]])) %>% 
    mutate(coefficient = rownames(coefTable(get.models(clutch_size_model1_d6_subset, m)[[1]]))) %>% 
    rename(estimate = Estimate, 
           SE = Std..Error)
  table_S8_00 <- left_join(x = table_S8_template, 
                         y = model_coef %>% select(-df), 
                         by = "coefficient")
  
  ## put table data in right format
  table_S8_01 <- table_S8_00 %>% 
    pivot_wider(names_from = coefficient, values_from = c(estimate, SE)) %>% 
    mutate(k = clutch_size_model1_d6_subset$df[m],
           AIC = clutch_size_model1_d6_subset$AIC[m],
           delta = clutch_size_model1_d6_subset$delta[m]) %>% 
    relocate(`Intercept` = `estimate_(Intercept)`,
             `Intercept SE` = `SE_(Intercept)`,
             `Number of helping females` = `estimate_female_helpers`,
             `Number of helping females SE` = `SE_female_helpers`,
             `Number of helping males` = `estimate_male_helpers`,
             `Number of helping males SE` = `SE_male_helpers`,
             `Clutch order` = `estimate_clutch_order`,            
             `Clutch order SE` = `SE_clutch_order`,
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
             decimals = 2) %>% 
  fmt_number(columns = ncol(table_S8_data)-2,
             decimals = 0) %>% 
  cols_merge_uncert(col_val = `Intercept`, col_uncert = `Intercept SE`) %>% 
  cols_merge_uncert(col_val = `Number of helping females`, col_uncert =`Number of helping females SE`) %>% 
  cols_merge_uncert(col_val = `Number of helping males`, col_uncert =`Number of helping males SE`) %>% 
  cols_merge_uncert(col_val = `Clutch order`, col_uncert =`Clutch order SE`) %>% 
  fmt_missing(columns = c(1:(ncol(table_S8_data)-3), (ncol(table_S8_data)-1):(ncol(table_S8_data))), 
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
table_S8_clean

# save Table S8 (saved in html, then imported in docx to include in manuscript)
table_S8_clean %>%
  gtsave(filename = "./tables/Table S8.html")



## 
##
#### Second model - partitioning variation in female and male helper number in clutch size model #####
##
##
clutch_size_model2 <- lmer(clutch_size ~ 
                             clutch_order +
                             cent_female_helpers +
                             cent_male_helpers + 
                             mean_female_helpers +
                             mean_male_helpers + 
                             (1|mother_ID) + 
                             (1|Season) + 
                             (1|Group_ID),
                           REML = FALSE,
                           na.action = "na.fail",
                           data = data)

# model selection based on AIC
clutch_size_model2_full_aic_table <- dredge(clutch_size_model2, rank = "AIC", trace = 3)
clutch_size_model2_d6_subset <- subset(clutch_size_model2_full_aic_table, delta < 6) 


##
## 
## create table with AIC results
names(model1_d6_subset)
number_variables <- 6
col_names <- names(clutch_size_model2_d6_subset)[1:number_variables]

##
##
##### Code to generate Table S9 #####
##
##

## add each model to the table
# list to store data
list_models_table_S9 <- as.list(NA)
for(m in 1:nrow(clutch_size_model2_d6_subset)){
  
  # template to store data
  table_S9_template <- data.frame(coefficient = names(clutch_size_model2_d6_subset)[1:number_variables]) 
  
  # add model coefficients
  model_coef <- data.frame(coefTable(get.models(clutch_size_model2_d6_subset, m)[[1]])) %>% 
    mutate(coefficient = rownames(coefTable(get.models(clutch_size_model2_d6_subset, m)[[1]]))) %>% 
    rename(estimate = Estimate, 
           SE = Std..Error)
  table_S9_00 <- left_join(x = table_S9_template, 
                           y = model_coef %>% select(-df), 
                           by = "coefficient")
  
  ## put table data in right format
  table_S9_01 <- table_S9_00 %>% 
    pivot_wider(names_from = coefficient, values_from = c(estimate, SE)) %>% 
    mutate(k = clutch_size_model2_d6_subset$df[m],
           AIC = clutch_size_model2_d6_subset$AIC[m],
           delta = clutch_size_model2_d6_subset$delta[m]) %>% 
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
             `Clutch order` = `estimate_clutch_order`,            
             `Clutch order SE` = `SE_clutch_order`,
             k = k, 
             AIC = AIC, 
             delta = delta)
  
  # save data per model
  list_models_table_S9[[m]] <- table_S9_01           
  
}

# combine data from each model in one table
table_S9_data <- rbindlist(list_models_table_S9)

# remove columns with all NA (remove variables that don't appear in any model in the table)
table_S9_data <- table_S9_data %>%
  select_if(~ !all(is.na(.)))

# form table
table_S9_clean <- table_S9_data %>% 
  gt() %>% 
  fmt_number(columns = c(1:(ncol(table_S9_data)-3), (ncol(table_S9_data)-1):(ncol(table_S9_data))),
             decimals = 2) %>% 
  fmt_number(columns = ncol(table_S9_data)-2,
             decimals = 0) %>% 
  cols_merge_uncert(col_val = `Intercept`, col_uncert = `Intercept SE`) %>% 
  cols_merge_uncert(col_val = `d Number of helping females`, col_uncert =`d Number of helping females SE`) %>% 
  cols_merge_uncert(col_val = `d Number of helping males`, col_uncert =`d Number of helping males SE`) %>% 
  cols_merge_uncert(col_val = `mean Number of helping females`, col_uncert =`mean Number of helping females SE`) %>% 
  cols_merge_uncert(col_val = `mean Number of helping males`, col_uncert =`mean Number of helping males SE`) %>% 
  cols_merge_uncert(col_val = `Clutch order`, col_uncert =`Clutch order SE`) %>% 
  fmt_missing(columns = c(1:(ncol(table_S9_data)-3), (ncol(table_S9_data)-1):(ncol(table_S9_data))), 
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

# TABLE S9
table_S9_clean

# save Table S9 (saved in html, then imported in docx to include in manuscript)
table_S9_clean %>%
  gtsave(filename = "./tables/Table S9.html")


## 
##
#### Third model - zero-trucated model of clutch size #####
##
##
data$mother_ID <- as.factor(data$mother_ID)
data$Season <- as.factor(data$Season)
data$Group_ID <- as.factor(data$Group_ID)

# model using glmmADMB package
truc_model_clutch_size <- glmmADMB::glmmadmb(clutch_size ~ 
                                               scale(clutch_order) + # improve speed 
                                               female_helpers +
                                               male_helpers + 
                                               (1|mother_ID) + 
                                               (1|Season) + 
                                               (1|Group_ID), 
                                             data=data,
                                             family="truncpoiss")
summary(truc_model_clutch_size)

##
## very similar coefficients to Gaussian models, and also no important effects
