###
###
#' 
#' Script for:
#' Mothers front-load their investment to the egg stage when helped in a wild cooperative bird
#' Capilla-Lasheras et al. 
#' Preprint: https://doi.org/10.1101/2021.11.11.468195
#' 
#' Latest update: 2022/08/09
#' 
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list=ls())

##
##
##### Script aim: #####
#' Analysis of maternal provisioning. 
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


# sample size in this model
nrow(data)            
length(unique(data$clutch_ID)) # 124 clutches
length(unique(data$mother_ID)) # 50 mothers
length(unique(data$group_ID))  # 34 social groups


# broods per mother
clutches_per_mother <- data %>% 
  group_by(mother_ID, clutch_ID) %>%
  filter(row_number() == 1) %>% 
  group_by(mother_ID) %>%
  summarise(n_clutches_mother = n())

# mother with more than one clutch?
table(clutches_per_mother$n_clutches_mother)
table(clutches_per_mother$n_clutches_mother > 1)

range(clutches_per_mother$n_clutches_mother) # range of clutches per mother
mean(clutches_per_mother$n_clutches_mother) # mean of clutches per mother
median(clutches_per_mother$n_clutches_mother) # median of clutches per mother


## 
## histogram of broods per mother
hist_clutches_mother <- ggplot(data = clutches_per_mother,
                               aes(x = n_clutches_mother)) +
  geom_histogram(fill = "#a50f15") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12)) +
  labs(x = "Number of broods", y = "Count of mothers") 

# save plots
ggsave(filename = "./plots/Figure S6b.png", 
       plot = hist_clutches_mother,
       device = "png", 
       units = "mm",
       width = 90, 
       height = 120)


##
##
##### 1 - Third model - population-level maternal provisioning rates #####
##
##

##
## full model - with interactions - including only biologically and justified predictors
model3_full_model <- lmer (maternal_prov_rate ~ 
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
summary(model3_full_model) # summary of model

# check normality
plot(residuals(model3_full_model))
hist(residuals(model3_full_model), freq=F) # model residuals
lines(density(x = rnorm(n = 10000,         # expected normal distribution
                        mean = 0, 
                        sd = sd(residuals(model3_full_model)))))

##
## LRT to assess importance of interactions
anova(model3_full_model, 
      update(model3_full_model, . ~ . -
               female_helpers : scale(brood_size) -
               male_helpers : scale(brood_size)), 
      test = "Chisq")


##
## full model - without interactions
model3_main_effects_model <- lmer (maternal_prov_rate ~ 
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
summary(model3_main_effects_model) # summary of model
drop1(update(object = model3_main_effects_model, REML = F), 
      test = "Chisq") # Likelihood-ratio test for each predictor

# formal test of differences between female helper and male helper effects within-mother from the best model containing both
car::linearHypothesis(model3_main_effects_model, "female_helpers - male_helpers = 0")




# table with model coefficients 
tab_model(model3_main_effects_model,
          file="./tables/Table 3 - model coefficients maternal provisioning rate.doc",
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
##### 2 - Third model - population-level maternal provisioning rates - AIC tables #####
##
##

# model selection based on AIC
model3_full_aic_table <- dredge(update(model3_full_model, 
                                       control = lmerControl( # ignore message for clean trace - the message has been checked in model above
                                         check.conv.singular = .makeCC(action = "ignore", 
                                                                       tol = formals(isSingular)$tol)
                                       )), 
                                subset=dc("poly(rainfall, 2)[, 1]",
                                          "poly(rainfall, 2)[, 2]"),
                                rank = "AIC", 
                                trace = 3)

# AIC results
model3_d6_subset <- subset(model3_full_aic_table, 
                                  delta < 6)

##
## 
## create table with AIC results
names(model3_full_aic_table)
number_variables <- 9
col_names <- names(model3_full_aic_table)[1:number_variables]

##
##
##### Code to generate Table 2 #####
##
##

## add each model to the table
# list to store data
list_models_tableS3 <- as.list(NA)
for(m in 1:nrow(model3_d6_subset)){
  
  # template to store data
  tableS3_template <- data.frame(coefficient = names(model3_d6_subset)[1:number_variables]) 
  
  # add model coeffiecients
  model_coef <- data.frame(coefTable(get.models(model3_d6_subset, m)[[1]])) %>% 
    mutate(coefficient = rownames(coefTable(get.models(model3_d6_subset, m)[[1]]))) %>% 
    rename(estimate = Estimate, 
           SE = Std..Error)
  tableS3_00 <- left_join(x = tableS3_template, 
                         y = model_coef %>% select(-df), 
                         by = "coefficient")
  
  ## put table data in right format
  tableS3_01 <- tableS3_00 %>% 
    pivot_wider(names_from = coefficient, values_from = c(estimate, SE)) %>% 
    mutate(k = model3_d6_subset$df[m],
           AIC = model3_d6_subset$AIC[m],
           delta = model3_d6_subset$delta[m]) %>% 
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
  list_models_tableS3[[m]] <- tableS3_01           
  
}

# combine data from each model in one table
tableS3_data <- rbindlist(list_models_tableS3)

# remove columns with all NA (remove variables that don't appear in any model in the table)
tableS3_data <- tableS3_data %>%
  select_if(~ !all(is.na(.)))

# form table
tableS3_clean <- tableS3_data %>% 
  gt() %>% 
  fmt_number(columns = c(1:(ncol(tableS3_data)-3), (ncol(tableS3_data)-1):(ncol(tableS3_data))),
             decimals = 2) %>% 
  fmt_number(columns = ncol(tableS3_data)-2,
             decimals = 0) %>% 
  cols_merge_uncert(col_val = `Intercept`, col_uncert = `Intercept SE`) %>% 
  cols_merge_uncert(col_val = `Number of helping females`, col_uncert =`Number of helping females SE`) %>% 
  cols_merge_uncert(col_val = `Number of helping males`, col_uncert =`Number of helping males SE`) %>%
  cols_merge_uncert(col_val = `Brood size`, col_uncert =`Brood size SE`) %>% 
  cols_merge_uncert(col_val = Rainfall1, col_uncert =`Rainfall1 SE`) %>% 
  cols_merge_uncert(col_val = Rainfall2, col_uncert =`Rainfall2 SE`) %>% 
  cols_merge_uncert(col_val = `Heat waves`, col_uncert =`Heat waves SE`) %>% 
  cols_merge_uncert(col_val = `Number of helping females x Brood size`, 
                    col_uncert =`Number of helping females x Brood size SE`) %>% 
  cols_merge_uncert(col_val = `Number of helping males x Brood size`, 
                    col_uncert =`Number of helping males x Brood size SE`) %>% 
  fmt_missing(columns = c(1:(ncol(tableS3_data)-3), (ncol(tableS3_data)-1):(ncol(tableS3_data))), 
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

# TABLE 2
tableS3_clean

# save Table 2 (saved in html, then imported in docx to include in manuscript)
tableS3_clean %>%
  gtsave(filename = "./tables/Table S3.html")


## remove table1 objects
rm(list = c("tableS3_00", 
            "tableS3_01", 
            "tableS3_clean", 
            "list_models_tableS3", 
            "tableS3_template", 
            "tableS3_data"))


##
##
##### Plotting model predictions #####
##
##

# re-fit model 3 to facilitate plotting
model3_plot <- lmer(maternal_prov_rate ~ 
                      # rainfall and heat waves
                      rainfall +
                      I(rainfall^2) +
                      temp_above_35 +
                      
                      # other main effects
                      female_helpers +
                      male_helpers +
                      brood_size +
                      (1|group_ID) +
                      (1|mother_ID) +
                      (1|season),
                    data = data)
summary(model3_plot)


## data table with values to predict
df_predict <- expand.grid(brood_size = mean(as.numeric(data$brood_size)),
                          male_helpers = mean(data$male_helpers),
                          female_helpers = seq(min(data$female_helpers),max(data$female_helpers),0.25),
                          rainfall = mean(data$rainfall),
                          temp_above_35 = mean(data$temp_above_35))

## model predictions
df_predict$fit <- predict(model3_plot, 
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
                     male_helpers +
                     brood_size,
                   data = df_predict)

pvar1 <- diag(mm %*% tcrossprod(vcov(model3_plot),mm))
cmult <- 1 
df_predict <- data.frame(
  df_predict
  , plo = df_predict$fit-cmult*sqrt(pvar1)
  , phi = df_predict$fit+cmult*sqrt(pvar1)
)



##
## model predictions for different number of female helpers, at the mean value of the other variables in the model
model3_plot <- ggplot(df_predict,
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
       plot = model3_plot,
       device = "png", 
       units = "mm",
       width = 73, 
       height = 89)
