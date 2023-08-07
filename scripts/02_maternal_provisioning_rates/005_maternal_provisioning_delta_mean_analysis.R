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
               performance,
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

####

##
##
##### partitioning the effects of within and among-mother helper number on maternal provisioning rates #####
##
##

# full model partition within and among-mother variation helper number
model4_main_effects_model <- lmer (maternal_prov_rate ~ 
                                     # rainfall and heat waves
                                     poly(rainfall,2, raw = F)[,1] +
                                     poly(rainfall,2,  raw = F)[,2] +
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
summary(model4_main_effects_model) # boundary (singular) warning due to group and season variances = 0. Removing these two random effects removes the warning and produces the same estimates
confint(model4_main_effects_model)
lmerTest::rand(model4_main_effects_model)# LRT for random effects
drop1(update(object = model4_main_effects_model, REML = F), 
      test = "Chisq") # Likelihood-ratio test for each predictor

# formal test of differences between female helper and male helper effects within-mother from the best model containing both
car::linearHypothesis(model4_main_effects_model, "cent_female_helpers - cent_male_helpers = 0")

# table with model coefficients
tab_model(model4_main_effects_model,
          file="./tables/Table 4 - RAW.doc",
          pred.labels = c("Intercept", 
                          html("Rainfall<sup>1</sup>"),
                          html("Rainfall<sup>2</sup>"),
                          "Heat waves",
                          html("&Delta;Number of female helpers"),
                          html("&mu;Number of female helpers"),
                          html("&Delta;Number of male helpers"),
                          html("&mu;Number of male helpers"),
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
## note that rainfall effects have been back transformed for reporting purposes from the orthogonal scale in which they appear in the table above


##
##
##### Fourth model - partitioning the effects of within and among-mother helper number on maternal provisioning rates - AIC tables #####
##
##


# model selection based on AIC
model4_full_aic_table <- dredge(update(model4_main_effects_model, 
                                       control = lmerControl( # ignore message for clean trace - the message has been checked in model above
                                         check.conv.singular = .makeCC(action = "ignore", 
                                                                       tol = formals(isSingular)$tol)
                                       )), 
                                subset=dc("poly(rainfall, 2)[, 1]",
                                          "poly(rainfall, 2)[, 2]"),
                                rank = "AIC", 
                                trace = 3)

# AIC results
model4_d6_subset <- subset(model4_full_aic_table, 
                           delta < 6)

##
## 
## create table with AIC results
names(model4_d6_subset)
number_variables <- 9
col_names <- names(model4_d6_subset)[1:number_variables]

##
## table

## add each model to the table
# list to store data
list_models_tableS11 <- as.list(NA)
for(m in 1:nrow(model4_d6_subset)){
  
  # template to store data
  tableS11_template <- data.frame(coefficient = names(model4_d6_subset)[1:number_variables]) 
  
  # add model coeffiecients
  model_coef <- data.frame(coefTable(get.models(model4_d6_subset, m)[[1]])) %>% 
    mutate(coefficient = rownames(coefTable(get.models(model4_d6_subset, m)[[1]]))) %>% 
    rename(estimate = Estimate, 
           SE = Std..Error)
  tableS11_00 <- left_join(x = tableS11_template, 
                          y = model_coef %>% select(-df), 
                          by = "coefficient")
  
  ## put table data in right format
  tableS11_01 <- tableS11_00 %>% 
    pivot_wider(names_from = coefficient, values_from = c(estimate, SE)) %>% 
    mutate(k = model4_d6_subset$df[m],
           AIC = model4_d6_subset$AIC[m],
           delta = model4_d6_subset$delta[m]) %>% 
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
             `Brood size` = `estimate_scale(brood_size)`,            # not present in models in Table 1
             `Brood size SE` = `SE_scale(brood_size)`,
             `Rainfall1` = `estimate_poly(rainfall, 2, raw = F)[, 1]`,
             `Rainfall1 SE` = `SE_poly(rainfall, 2, raw = F)[, 1]`,
             `Rainfall2` = `estimate_poly(rainfall, 2, raw = F)[, 2]`,
             `Rainfall2 SE` = `SE_poly(rainfall, 2, raw = F)[, 2]`,
             `Heat waves` = `estimate_scale(temp_above_35)`,
             `Heat waves SE` = `SE_scale(temp_above_35)`,
             k = k, 
             AIC = AIC, 
             delta = delta)
  
  # save deta per model
  list_models_tableS11[[m]] <- tableS11_01           
  
}

# combine data from each model in one table
tableS11_data <- rbindlist(list_models_tableS11)

# remove columns with all NA (remove variables that don't appear in any model in the table)
tableS11_data <- tableS11_data %>%
  select_if(~ !all(is.na(.)))

# form table
tableS11_clean <- tableS11_data %>% 
  gt() %>% 
  fmt_number(columns = c(1:(ncol(tableS11_data)-3), (ncol(tableS11_data)-1):(ncol(tableS11_data))),
             decimals = 2) %>% 
  fmt_number(columns = ncol(tableS11_data)-2,
             decimals = 0) %>% 
  cols_merge_uncert(col_val = `Intercept`, col_uncert = `Intercept SE`) %>% 
  cols_merge_uncert(col_val = `d Number of helping females`, col_uncert =`d Number of helping females SE`) %>% 
  cols_merge_uncert(col_val = `d Number of helping males`, col_uncert =`d Number of helping males SE`) %>% 
  cols_merge_uncert(col_val = `mean Number of helping females`, col_uncert =`mean Number of helping females SE`) %>% 
  cols_merge_uncert(col_val = `mean Number of helping males`, col_uncert =`mean Number of helping males SE`) %>% 
  cols_merge_uncert(col_val = `Brood size`, col_uncert =`Brood size SE`) %>% 
  cols_merge_uncert(col_val = Rainfall1, col_uncert =`Rainfall1 SE`) %>% 
  cols_merge_uncert(col_val = Rainfall2, col_uncert =`Rainfall2 SE`) %>% 
  cols_merge_uncert(col_val = `Heat waves`, col_uncert =`Heat waves SE`) %>% 
  fmt_missing(columns =  c(1:(ncol(tableS11_data)-3), (ncol(tableS11_data)-1):(ncol(tableS11_data))), 
              missing_text = " ") %>% 
  cols_label(`d Number of helping females` = html("&Delta; Number of helping females"),
             `d Number of helping males` = html("&Delta; Number of helping males"),
             `mean Number of helping females` = html("&mu; Number of helping females"),
             `mean Number of helping males` = html("&mu; Number of helping males"),
             `Rainfall1` = html("Rainfall<sup>1</sup>"),
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

tableS11_clean
tableS11_clean %>%
  gtsave(filename = "./tables/Table S11.html")

#####


##
##
##### Plotting model predictions #####
##
##

# re-fit model 4 to facilitate plotting
model4_plot <- lmer (maternal_prov_rate ~ 
                       # rainfall and heat waves
                       rainfall +
                       I(rainfall^2) +
                       temp_above_35 +
                       
                       # other main effects
                       cent_female_helpers +
                       mean_female_helpers +
                       cent_male_helpers +
                       mean_male_helpers +
                       brood_size +
                       (1|group_ID) +
                       (1|mother_ID) +
                       (1|season),
                     data = data,
                     na.action = "na.fail") 
summary(model4_plot)

## data table with values to predict
df_predict <- expand.grid(brood_size = mean(as.numeric(data$brood_size)),
                          mean_male_helpers = mean(data$mean_male_helpers),
                          cent_male_helpers = mean(data$cent_male_helpers),
                          cent_female_helpers = seq(min(data$cent_female_helpers),
                                                    max(data$cent_female_helpers),0.25),
                          mean_female_helpers = mean(data$mean_female_helpers),
                          rainfall = mean(data$rainfall),
                          temp_above_35 = mean(data$temp_above_35))

## model predictions
df_predict$fit <- predict(model4_plot, 
                          df_predict, 
                          re.form = NA, 
                          type = "response")

## calculate SE for predictions
mm <- model.matrix(~ # rainfall and heat waves
                     rainfall +
                     I(rainfall^2) +
                     temp_above_35 +
                     
                     # other main effects
                     cent_female_helpers +
                     mean_female_helpers +
                     cent_male_helpers +
                     mean_male_helpers +
                     brood_size,
                   data = df_predict)

pvar1 <- diag(mm %*% tcrossprod(vcov(model4_plot),mm))
cmult <- 1 
df_predict <- data.frame(
  df_predict
  , plo = df_predict$fit-cmult*sqrt(pvar1)
  , phi = df_predict$fit+cmult*sqrt(pvar1)
)

##
## model predictions for different number of female helpers, at the mean value of the other variables in the model
model4_plot <- ggplot(df_predict,
                           aes(x = cent_female_helpers, 
                               y = fit)) +
  geom_line(size = 1, color = "#a50f15") +
  geom_ribbon(aes(ymin = (plo), ymax = (phi)), 
              alpha = 0.2,
              fill = "#a50f15") +
  geom_point(data = data,
             aes(y = maternal_prov_rate,
                 x = cent_female_helpers),
             size = 2.5,
             position = position_jitter(width = 0.1),
             alpha = 0.20) +
  theme_bw() +
  labs(x = expression(atop(paste(Delta, " female helper number", sep = " "),
                           "(within-mother variation)")), 
       y = "Maternal provisioning rate (feeds/hour)") + 
  theme(plot.margin = unit(c(5,5,30,20), "points"),
        panel.grid = element_blank(),
        axis.title.x = element_text(family = "Arial", size = 12, vjust = -2),
        axis.title.y = element_text(family = "Arial", size = 10, vjust = +3),
        axis.text = element_text(family = "Arial", size = 10)) 

ggsave(filename = "./plots/Figure 2d.png", 
       plot = model4_plot,
       device = "png", 
       units = "mm",
       width = 73, 
       height = 89)


