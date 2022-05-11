###
###
#' 
#' Script for:
#' Mothers front-load their investment to the egg stage when helped in a wild cooperative bird
#' Capilla-Lasheras et al. 
#' Preprint: https://doi.org/10.1101/2021.11.11.468195
#' 
#' Latest update: 2022/04/01
#' 
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list=ls())

##
##
##### Script aim: #####
#' Analysis of maternal provisioning. Results presented in Table 2 and Table S5 - Table S8
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
##### Maternal provisioning rate dataset #####
##
##
data <- read.csv("./data/maternal_provisioning_dataset.csv")
head(data)

##
##
##### First model of maternal provisioning rates #####
##
##

# sample size in this model
length(unique(data$clutch_ID)) # 118 broods clutches
length(unique(data$mother_ID)) # 50 mothers
length(unique(data$group_ID))  # 34 social groups
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


# full model - including only biologically and justified predictors
model4_full <- lmer (maternal_prov_rate ~ 
                       # rainfall and heat waves
                       poly(rainfall,2)[,1] +
                       poly(rainfall,2)[,2] +
                       scale(temp_above_35) +
                       
                       # interactions
                       female_helpers : scale(brood_size) +
                       male_helpers : scale(brood_size) +
                       
                       # other main effects
                       female_helpers +
                       male_helpers +
                       scale(brood_size) +
                       (1|group_ID) +
                       (1|mother_ID) +
                       (1|season),
                     data = data,
                     na.action = "na.fail",
                     REML = FALSE) # boundary (singular) warning due to group_ID variance = 0. Removing this random effect removes the warning and produces the same estimates
summary(model4_full)

# model selection based on AIC
model4_full_aic_table <- dredge(update(model4_full, 
                                       control = lmerControl( # ignore message for clean trace - the message has been checked in model above
                                         check.conv.singular = .makeCC(action = "ignore", 
                                                                       tol = formals(isSingular)$tol)
                                       )), 
                                subset=dc("poly(rainfall, 2)[, 1]",
                                          "poly(rainfall, 2)[, 2]"),
                                rank = "AIC", 
                                trace = 3)

# AIC results - Table 2 & Table S5
model4_d6_subset <- subset(model4_full_aic_table, # for Table S5
                           delta < 6)
model4_d6_nested_subset <- subset(model4_full_aic_table, # for Table 2
                                  delta < 6 & !nested(.))
summary(get.models(model4_d6_subset,1)[[1]]) # model output from the top model

##
## 
## create table with AIC results
names(model4_d6_nested_subset)
number_variables <- 9
col_names <- names(model4_d6_nested_subset)[1:number_variables]

##
##
##### Code to generate Table 2 #####
##
##

## add each model to the table
# list to store data
list_models_table2 <- as.list(NA)
for(m in 1:nrow(model4_d6_nested_subset)){
  
  # template to store data
  table2_template <- data.frame(coefficient = names(model4_d6_nested_subset)[1:number_variables]) 
  
  # add model coeffiecients
  model_coef <- data.frame(coefTable(get.models(model4_d6_nested_subset, m)[[1]])) %>% 
    mutate(coefficient = rownames(coefTable(get.models(model4_d6_nested_subset, m)[[1]]))) %>% 
    rename(estimate = Estimate, 
           SE = Std..Error)
  table2_00 <- left_join(x = table2_template, 
                         y = model_coef %>% select(-df), 
                         by = "coefficient")
  
  ## put table data in right format
  table2_01 <- table2_00 %>% 
    pivot_wider(names_from = coefficient, values_from = c(estimate, SE)) %>% 
    mutate(k = model4_d6_nested_subset$df[m],
           AIC = model4_d6_nested_subset$AIC[m],
           delta = model4_d6_nested_subset$delta[m]) %>% 
    relocate(`Intercept` = `estimate_(Intercept)`,
             `Intercept SE` = `SE_(Intercept)`,
             `Number of helping females` = `estimate_female_helpers`,
             `Number of helping females SE` = `SE_female_helpers`,
             `Number of helping males` = `estimate_male_helpers`,
             `Number of helping males SE` = `SE_male_helpers`,
             `Brood size` = `estimate_scale(brood_size)`,            # not present in models in Table 1
             `Brood size SE` = `SE_scale(brood_size)`,
             `Rainfall1` = `estimate_poly(rainfall, 2)[, 1]`,
             `Rainfall1 SE` = `SE_poly(rainfall, 2)[, 1]`,
             `Rainfall2` = `estimate_poly(rainfall, 2)[, 2]`,
             `Rainfall2 SE` = `SE_poly(rainfall, 2)[, 2]`,
             `Heat waves` = `estimate_scale(temp_above_35)`,
             `Heat waves SE` = `SE_scale(temp_above_35)`,
             `Number of helping females x Clutch size` = `estimate_female_helpers:scale(brood_size)`, # not present in models in Table 1
             `Number of helping females x Clutch size SE` = `SE_female_helpers:scale(brood_size)`,
             `Number of helping males x Clutch size` = `estimate_male_helpers:scale(brood_size)`,     # not present in models in Table 1
             `Number of helping males x Clutch size SE` = `SE_male_helpers:scale(brood_size)`,
             k = k, 
             AIC = AIC, 
             delta = delta)
  
  # save deta per model
  list_models_table2[[m]] <- table2_01           
  
}

# combine data from each model in one table
table2_data <- rbindlist(list_models_table2)

# remove columns with all NA (remove variables that don't appear in any model in the table)
table2_data <- table2_data %>%
  select_if(~ !all(is.na(.)))

# form table
table2_clean <- table2_data %>% 
  gt() %>% 
  fmt_number(columns = c(1:(ncol(table2_data)-3), (ncol(table2_data)-1):(ncol(table2_data))),
             decimals = 2) %>% 
  fmt_number(columns = ncol(table2_data)-2,
             decimals = 0) %>% 
  cols_merge_uncert(col_val = `Intercept`, col_uncert = `Intercept SE`) %>% 
  cols_merge_uncert(col_val = `Number of helping females`, col_uncert =`Number of helping females SE`) %>% 
  cols_merge_uncert(col_val = `Brood size`, col_uncert =`Brood size SE`) %>% 
  cols_merge_uncert(col_val = Rainfall1, col_uncert =`Rainfall1 SE`) %>% 
  cols_merge_uncert(col_val = Rainfall2, col_uncert =`Rainfall2 SE`) %>% 
  cols_merge_uncert(col_val = `Heat waves`, col_uncert =`Heat waves SE`) %>% 
  fmt_missing(columns = c(1:(ncol(table2_data)-3), (ncol(table2_data)-1):(ncol(table2_data))), 
              missing_text = " ") %>% 
  cols_label(`Rainfall1` = html("Rainfall<sup>1</sup>"),
             `Rainfall2` = html("Rainfall<sup>2</sup>"),
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

# TABLE 1
table2_clean

# save Table 1 (saved in html, then imported in docx to include in manuscript)
table2_clean %>%
  gtsave(filename = "./tables/Table 2.html")


## remove table1 objects
rm(list = c("table2_00", 
            "table2_01", 
            "table2_clean", 
            "list_models_table2", 
            "table2_template", 
            "table2_data"))


##
##
##### Code to generate Table S5 #####
##
##

## add each model to the table
# list to store data
list_models_tableS5 <- as.list(NA)
for(m in 1:nrow(model4_d6_subset)){
  
  # template to store data
  tableS5_template <- data.frame(coefficient = names(model4_d6_subset)[1:number_variables]) 
  
  # add model coeffiecients
  model_coef <- data.frame(coefTable(get.models(model4_d6_subset, m)[[1]])) %>% 
    mutate(coefficient = rownames(coefTable(get.models(model4_d6_subset, m)[[1]]))) %>% 
    rename(estimate = Estimate, 
           SE = Std..Error)
  tableS5_00 <- left_join(x = tableS5_template, 
                          y = model_coef %>% select(-df), 
                          by = "coefficient")
  
  ## put table data in right format
  tableS5_01 <- tableS5_00 %>% 
    pivot_wider(names_from = coefficient, values_from = c(estimate, SE)) %>% 
    mutate(k = model4_d6_subset$df[m],
           AIC = model4_d6_subset$AIC[m],
           delta = model4_d6_subset$delta[m]) %>% 
    relocate(`Intercept` = `estimate_(Intercept)`,
             `Intercept SE` = `SE_(Intercept)`,
             `Number of helping females` = `estimate_female_helpers`,
             `Number of helping females SE` = `SE_female_helpers`,
             `Number of helping males` = `estimate_male_helpers`,
             `Number of helping males SE` = `SE_male_helpers`,
             `Brood size` = `estimate_scale(brood_size)`,            # not present in models in Table 1
             `Brood size SE` = `SE_scale(brood_size)`,
             `Rainfall1` = `estimate_poly(rainfall, 2)[, 1]`,
             `Rainfall1 SE` = `SE_poly(rainfall, 2)[, 1]`,
             `Rainfall2` = `estimate_poly(rainfall, 2)[, 2]`,
             `Rainfall2 SE` = `SE_poly(rainfall, 2)[, 2]`,
             `Heat waves` = `estimate_scale(temp_above_35)`,
             `Heat waves SE` = `SE_scale(temp_above_35)`,
             `Number of helping females x Brood size` = `estimate_female_helpers:scale(brood_size)`, # not present in models in Table 1
             `Number of helping females x Brood size SE` = `SE_female_helpers:scale(brood_size)`,
             `Number of helping males x Brood size` = `estimate_male_helpers:scale(brood_size)`,     # not present in models in Table 1
             `Number of helping males x Brood size SE` = `SE_male_helpers:scale(brood_size)`,
             k = k, 
             AIC = AIC, 
             delta = delta)
  
  # save deta per model
  list_models_tableS5[[m]] <- tableS5_01           
  
}

# combine data from each model in one table
tableS5_data <- rbindlist(list_models_tableS5)

# remove columns with all NA (remove variables that don't appear in any model in the table)
tableS5_data <- tableS5_data %>%
  select_if(~ !all(is.na(.)))

# form table
tableS5_clean <- tableS5_data %>% 
  gt() %>% 
  fmt_number(columns = c(1:(ncol(tableS5_data)-3), (ncol(tableS5_data)-1):(ncol(tableS5_data))),
             decimals = 2) %>% 
  fmt_number(columns = ncol(tableS5_data)-2,
             decimals = 0) %>% 
  cols_merge_uncert(col_val = `Intercept`, col_uncert = `Intercept SE`) %>% 
  cols_merge_uncert(col_val = `Number of helping females`, col_uncert =`Number of helping females SE`) %>% 
  cols_merge_uncert(col_val = `Number of helping males`, col_uncert =`Number of helping males SE`) %>% 
  cols_merge_uncert(col_val = `Brood size`, col_uncert =`Brood size SE`) %>% 
  cols_merge_uncert(col_val = Rainfall1, col_uncert =`Rainfall1 SE`) %>% 
  cols_merge_uncert(col_val = Rainfall2, col_uncert =`Rainfall2 SE`) %>% 
  cols_merge_uncert(col_val = `Heat waves`, col_uncert =`Heat waves SE`) %>% 
  cols_merge_uncert(col_val = `Number of helping females x Brood size`, col_uncert =`Number of helping females x Brood size SE`) %>% 
  cols_merge_uncert(col_val = `Number of helping males x Brood size`, col_uncert =`Number of helping males x Brood size SE`) %>% 
  fmt_missing(columns = 1:(ncol(tableS5_data)-3), 
              missing_text = " ") %>% 
  cols_label(`Rainfall1` = html("Rainfall<sup>1</sup>"),
             `Rainfall2` = html("Rainfall<sup>2</sup>"),
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

# TABLE S5
tableS5_clean

# save Table S5 (saved in html, then imported in docx to include in manuscript)
tableS5_clean %>%
  gtsave(filename = "./tables/Table S5.html")



##
## full model, no model selection - without interactions - results in Table S6
model4_main_effects <- lmer (maternal_prov_rate ~ 
                               # rainfall and heat waves
                               poly(rainfall,2)[,1] +
                               poly(rainfall,2)[,2] +
                               scale(temp_above_35) +
                               
                               # other main effects
                               female_helpers +
                               male_helpers +
                               scale(brood_size) +
                               (1|group_ID) +
                               (1|mother_ID) +
                               (1|season),
                             data = data,
                             na.action = "na.fail",
                             REML = FALSE) # boundary (singular) warning due to group_ID variance = 0. Removing this random effect removes the warning and produces the same estimates
summary(model4_main_effects) # summary of model
drop1(update(object = model4_main_effects, REML = F), 
      test = "Chisq") # Likelihood-ratio test for each predictor

# table with model coefficients - included as Table S7 in manuscript
tab_model(model4_main_effects,
          file="./tables/Table S6 - model coefficients maternal provisioning rate.doc",
          pred.labels = c("Intercept", 
                          "Rainfall^1", 
                          "Rainfall^2",
                          "Heat waves",
                          "Number of female helpers",
                          "Number of male helpers",
                          "Brood size"),
          string.ci = "CI (95%)",
          string.se = "SE",
          show.se = TRUE, 
          show.stat = FALSE,
          show.p = FALSE,
          show.est = TRUE,
          show.intercept = TRUE,
          rm.terms = NULL,
          show.re.var = FALSE,
          show.ngroups = FALSE,
          show.r2 = FALSE,
          show.obs = FALSE,
          ci.hyphen = ",",
          use.viewer = T)
## likelihood-ratio results included in the table above manually


##
##
##### Plotting model predictions from the top performing models above #####
##
##

## 
## Top model (Table 2)
model4_best <- lmer(maternal_prov_rate ~ 
                     # rainfall and heat waves
                     rainfall +
                     I(rainfall^2) +
                     temp_above_35 +
                     
                     # other main effects
                     female_helpers +
                     brood_size +
                     (1|group_ID) +
                     (1|mother_ID) +
                     (1|season),
                   data = data)
summary(model4_best)

plot(fitted(model4_best), residuals(model4_best))
hist( residuals(model4_best))
shapiro.test(residuals(model4_best))

## data table with values to predict
df_predict <- expand.grid(brood_size = mean(as.numeric(data$brood_size)),
                          female_helpers = seq(min(data$female_helpers),max(data$female_helpers),0.25),
                          rainfall = mean(data$rainfall),
                          temp_above_35 = mean(data$temp_above_35))

## model predictions
df_predict$fit <- predict(model4_best, 
                          df_predict, 
                          re.form = NA, 
                          type = "response")

## calculate SE for predictions
mm <- model.matrix(~ # rainfall and heat waves
                     rainfall +
                     I(rainfall^2) +
                     temp_above_35 +
                     
                     # other main effects
                     female_helpers +
                     brood_size,
                   data = df_predict)

pvar1 <- diag(mm %*% tcrossprod(vcov(model4_best),mm))
cmult <- 1 
df_predict <- data.frame(
  df_predict
  , plo = df_predict$fit-cmult*sqrt(pvar1)
  , phi = df_predict$fit+cmult*sqrt(pvar1)
)



##
## model predictions for different number of female helpers, at the mean value of the other variables in the model
model4_best_plot <- ggplot(df_predict,
                      aes(x = female_helpers, 
                          y = fit)) +
  geom_line(size = 1, color = "#a50f15") +
  geom_ribbon(aes(ymin = (plo), ymax = (phi)), 
              alpha = 0.2,
              fill = "#a50f15") +
  geom_point(data = data,
             aes(y = maternal_prov_rate,
                 x = female_helpers),
             size = 2.5,
             position = position_jitter(width = 0.1),
             alpha = 0.20) +
  theme_bw() +
  labs(x = expression(atop("Number of female helpers", " ")),
       y = "Maternal provisioning rate (feeds/hour)") + 
  theme(plot.margin = unit(c(5,5,30,20), "points"),
        panel.grid = element_blank(),
        axis.title.x = element_text(family = "Arial", size = 12, vjust = -2),
        axis.title.y = element_text(family = "Arial", size = 10, vjust = +3),
        axis.text = element_text(family = "Arial", size = 10)) 

ggsave(filename = "./plots/Figure 2c.png", 
       plot = model4_best_plot,
       device = "png", 
       units = "mm",
       width = 73, 
       height = 89)