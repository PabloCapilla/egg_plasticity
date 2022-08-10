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
#' Analysis of egg volume. Results presented in Table 1, Table S1 and Table S3
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
               MetBrewer,
               ggpubr,
               extrafont,
               MuMIn, 
               lme4)
loadfonts()
# font_import() may be needed first to use font types in ggplot
# check with 'fonts()' that font are available

##
##
##### egg volume dataset #####
##
##
data <- read.csv("./data/egg_volume_dataset.csv") 
head(data)
nrow(data)

##
##
##### Summary of sample sizes and replication of mothers #####
##
##

# sample size in this model
length(unique(data$clutch_ID)) # 271 clutches
length(unique(data$mother_ID)) # 62 mothers
length(unique(data$Group))     # 37 social groups

# eggs per mother
egg_per_mother <- data %>% 
    group_by(mother_ID) %>%
    summarise(n_eggs_mother = n())

range(egg_per_mother$n_eggs_mother) # range of eggs per mother
mean(egg_per_mother$n_eggs_mother) # mean of eggs per mother
median(egg_per_mother$n_eggs_mother) # median of eggs per mother

# clutches per mother
clutches_per_mother <- data %>% 
  group_by(mother_ID, clutch_ID) %>%
  filter(row_number() == 1) %>% 
  group_by(mother_ID) %>%
  summarise(n_clutches_mother = n())

range(clutches_per_mother$n_clutches_mother) # range of clutches per mother
mean(clutches_per_mother$n_clutches_mother) # mean of clutches per mother
median(clutches_per_mother$n_clutches_mother) # median of clutches per mother


## 
## histogram of clutches per mother
hist_egg_mother <- ggplot(data = egg_per_mother,
                          aes(x = n_eggs_mother)) +
  geom_histogram(fill = "#a50f15", binwidth = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12)) +
  labs(x = "Number of eggs", y = "Count of mothers")


# save plot
ggsave(filename = "./plots/Figure S6a.png", 
       plot = hist_egg_mother,
       device = "png", 
       units = "mm",
       width = 90, 
       height = 120)
  
##
##
##### Distribution of number of female and male helpers #####
##
##
df_helpers_summary <- data %>% 
  group_by(clutch_ID) %>% 
  filter(row_number() == 1) %>% 
  select(clutch_ID, 
         female_helpers,
         male_helpers) %>% 
  pivot_longer(cols = 2:3, names_to = "sex", values_to = "number_ind") %>% 
  separate(col = "sex", sep = "_", into = c("sex", "useless"))

distr_helpers <- ggplot(data = df_helpers_summary, 
       aes(x = number_ind, fill = sex)) +
  geom_bar(position = "dodge") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12)) +
  labs(x = "Number of helpers", y = "Count of clutches") +
  scale_fill_manual(values = met.brewer("Java", 2),
                    name = "Helper sex")

##
##
##### 1 - First model - population-level egg volume model #####
##
##

##
## full model - with interactions - including only biologically and justified predictors
model1_full_model <- lmer(egg_volume ~ 
                              # rainfall and heat waves
                              poly(rainfall,2)[,1] +
                              poly(rainfall,2)[,2] +
                              scale(temp_above_35) +
                              
                              # interactions
                              female_helpers : scale(clutch_size) +
                              male_helpers : scale(clutch_size) +
                              female_helpers : scale(egg_position) +
                              male_helpers : scale(egg_position) +
                              
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
                            REML = T)
summary(model1_full_model) # summary of model

# check normality
plot(residuals(model1_full_model))
hist(residuals(model1_full_model), freq=F) # model residuals
lines(density(x = rnorm(n = 10000,          # expected normal distribution
                        mean = 0, 
                        sd = sd(residuals(model1_full_model)))))
##
## LRT to assess importance of interactions
anova(model1_full_model, 
      update(model1_full_model, . ~ . -
               female_helpers : scale(clutch_size) -
               male_helpers : scale(clutch_size) -
               female_helpers : scale(egg_position) -
               male_helpers : scale(egg_position)), 
             test = "Chisq")

##
## full model - without interactions
model1_main_effects_model <- lmer(egg_volume ~ 
                            # rainfall and heat waves
                            poly(rainfall,2)[,1] +
                            poly(rainfall,2)[,2] +
                            scale(temp_above_35) +

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
summary(model1_main_effects_model) # summary of model
drop1(update(object = model1_main_effects_model, REML = F), 
      test = "Chisq") # Likelihood-ratio test for each predictor

##
## repeatability of egg size
as.numeric(summary(model1_full_model)$varcor[[2]]) /      # variance among mothers
  (as.numeric(summary(model1_full_model)$varcor[[1]]) +   # variance among clutches
     as.numeric(summary(model1_full_model)$varcor[[2]]) + # variance among mothers
     as.numeric(summary(model1_full_model)$varcor[[3]]) + # variance among groups
     as.numeric(summary(model1_full_model)$varcor[[4]]) + # variance among seasons
     as.numeric((summary(model1_full_model)$sigma)^2))    # residual variance

# table with model coefficients - included as Table S3 in manuscript
tab_model(model1_main_effects_model,
          file="./tables/Table 1 - model coefficients egg volume - RAW.doc",
          pred.labels = c("Intercept", 
                          "Rainfall^1", 
                          "Rainfall^2",
                          "Heat waves",
                          "Number of female helpers",
                          "Number of male helpers",
                          "Clutch size",
                          "Egg position"),
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


# formal test of differences between female helper and male helper effects
car::linearHypothesis(model1_main_effects_model, "female_helpers - male_helpers = 0")


##
##
##### 2 - First model - population-level egg volume model - AIC tables #####
##
##

# full model 
model1_full <- lmer(egg_volume ~ 
                      # rainfall and heat waves
                      poly(rainfall,2)[,1] +
                      poly(rainfall,2)[,2] +
                      scale(temp_above_35) +
                      
                      # interactions
                      female_helpers : scale(clutch_size) +
                      male_helpers : scale(clutch_size) +
                      female_helpers : scale(egg_position) +
                      male_helpers : scale(egg_position) +
                      
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
summary(model1_full) # boundary (singular) warning due to group and season variances = 0. Removing these two random effects removes the warning and produces the same estimates

# model selection based on AIC
model1_full_aic_table <- dredge(update(model1_full, 
                                       control = lmerControl( # ignore message for clean trace - the message has been checked in model above
                                         check.conv.singular = .makeCC(action = "ignore", 
                                                                       tol = formals(isSingular)$tol)
                                         )), 
                                       trace = 3, 
                                       subset = dc("poly(rainfall, 2)[, 1]", # code to always include main effect when quadratic rainfall effect included
                                                   "poly(rainfall, 2)[, 2]"),
                                       rank = "AIC")
                                
# AIC results - Table S1
model1_d6_subset <- subset(model1_full_aic_table, 
                                  delta < 6)


##
## 
## create table with AIC results
names(model1_d6_subset)
number_variables <- 12
col_names <- names(model1_d6_subset)[1:number_variables]

##
##
##### Code to generate Table 1 #####
##
##

## add each model to the table
# list to store data
list_models_tableS1 <- as.list(NA)
for(m in 1:nrow(model1_d6_subset)){
  
  # template to store data
  tableS1_template <- data.frame(coefficient = names(model1_d6_subset)[1:number_variables]) 
  
  # add model coeffiecients
  model_coef <- data.frame(coefTable(get.models(model1_d6_subset, m)[[1]])) %>% 
    mutate(coefficient = rownames(coefTable(get.models(model1_d6_subset, m)[[1]]))) %>% 
    rename(estimate = Estimate, 
           SE = Std..Error)
  tableS1_00 <- left_join(x = tableS1_template, 
                         y = model_coef %>% select(-df), 
                         by = "coefficient")
  
  ## put table data in right format
  tableS1_01 <- tableS1_00 %>% 
    pivot_wider(names_from = coefficient, values_from = c(estimate, SE)) %>% 
    mutate(k = model1_d6_subset$df[m],
           AIC = model1_d6_subset$AIC[m],
           delta = model1_d6_subset$delta[m]) %>% 
    relocate(`Intercept` = `estimate_(Intercept)`,
             `Intercept SE` = `SE_(Intercept)`,
             `Number of helping females` = `estimate_female_helpers`,
             `Number of helping females SE` = `SE_female_helpers`,
             `Number of helping males` = `estimate_male_helpers`,
             `Number of helping males SE` = `SE_male_helpers`,
             `Clutch size` = `estimate_scale(clutch_size)`,            
             `Clutch size SE` = `SE_scale(clutch_size)`,
             `Egg position` = `estimate_scale(egg_position)`,
             `Egg position SE` = `SE_scale(egg_position)`,
             `Rainfall1` = `estimate_poly(rainfall, 2)[, 1]`,
             `Rainfall1 SE` = `SE_poly(rainfall, 2)[, 1]`,
             `Rainfall2` = `estimate_poly(rainfall, 2)[, 2]`,
             `Rainfall2 SE` = `SE_poly(rainfall, 2)[, 2]`,
             `Heat waves` = `estimate_scale(temp_above_35)`,
             `Heat waves SE` = `SE_scale(temp_above_35)`,
             `Number of helping females x Clutch size` = `estimate_female_helpers:scale(clutch_size)`, 
             `Number of helping females x Clutch size SE` = `SE_female_helpers:scale(clutch_size)`,
             `Number of helping males x Clutch size` = `estimate_male_helpers:scale(clutch_size)`,     
             `Number of helping males x Clutch size SE` = `SE_male_helpers:scale(clutch_size)`,
             `Number of helping females x Egg position` = `estimate_female_helpers:scale(egg_position)`,
             `Number of helping females x Egg position SE` = `SE_female_helpers:scale(egg_position)`,
             `Number of helping males x Egg position` = `estimate_male_helpers:scale(egg_position)`,    
             `Number of helping males x Egg position SE` = `SE_male_helpers:scale(egg_position)`,
             k = k, 
             AIC = AIC, 
             delta = delta)
  
  # save data per model
  list_models_tableS1[[m]] <- tableS1_01           
  
}

# combine data from each model in one table
tableS1_data <- rbindlist(list_models_tableS1)

# remove columns with all NA (remove variables that don't appear in any model in the table)
tableS1_data <- tableS1_data %>%
  select_if(~ !all(is.na(.)))

# form table
tableS1_clean <- tableS1_data %>% 
  gt() %>% 
  fmt_number(columns = c(1:(ncol(tableS1_data)-3), (ncol(tableS1_data)-1):(ncol(tableS1_data))),
             decimals = 2) %>% 
  fmt_number(columns = ncol(tableS1_data)-2,
             decimals = 0) %>% 
  cols_merge_uncert(col_val = `Intercept`, col_uncert = `Intercept SE`) %>% 
  cols_merge_uncert(col_val = `Number of helping females`, col_uncert =`Number of helping females SE`) %>% 
  cols_merge_uncert(col_val = `Number of helping males`, col_uncert =`Number of helping males SE`) %>% 
  cols_merge_uncert(col_val = `Clutch size`, col_uncert =`Clutch size SE`) %>%
  cols_merge_uncert(col_val = `Egg position`, col_uncert =`Egg position SE`) %>% 
  cols_merge_uncert(col_val = Rainfall1, col_uncert =`Rainfall1 SE`) %>% 
  cols_merge_uncert(col_val = Rainfall2, col_uncert =`Rainfall2 SE`) %>% 
  cols_merge_uncert(col_val = `Heat waves`, col_uncert =`Heat waves SE`) %>% 
  cols_merge_uncert(col_val = `Number of helping females x Clutch size`, 
                    col_uncert =`Number of helping females x Clutch size SE`) %>% 
  cols_merge_uncert(col_val = `Number of helping males x Clutch size`, 
                    col_uncert =`Number of helping males x Clutch size SE`) %>% 
  cols_merge_uncert(col_val = `Number of helping females x Egg position`, 
                    col_uncert =`Number of helping females x Egg position SE`) %>% 
  cols_merge_uncert(col_val = `Number of helping males x Egg position`, 
                    col_uncert =`Number of helping males x Egg position SE`) %>% 
  fmt_missing(columns = c(1:(ncol(tableS1_data)-3), (ncol(tableS1_data)-1):(ncol(tableS1_data))), 
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
tableS1_clean

# save Table 1 (saved in html, then imported in docx to include in manuscript)
tableS1_clean %>%
  gtsave(filename = "./tables/Table S1.html")


## remove table1 objects
rm(list = c("table1_00", 
            "table1_01", 
            "table1_clean", 
            "list_models_table1", 
            "table1_template", 
            "tableS1_data"))

##
##
##### Plotting model predictions #####
##
##

# re-fit model 1 to facilitate plotting
model1_plot <- lmer(egg_volume ~ 
                      # rainfall and heat waves
                      rainfall +
                      I(rainfall^2) +
                      temp_above_35 +
                      
                      # other main effects
                      female_helpers +
                      male_helpers +
                      clutch_size +
                      egg_position +
                      
                      # random effects
                      (1|Group) +
                      (1|mother_ID) +
                      (1|Season) + 
                      (1|clutch_ID),
                    data = data, 
                    REML = T) # warning due to large values in Rainfall^2 - note female effect (which is what is going to be plotted) doesn't vary in comparison with the top model above, with rainfall predictor scaled and without this warning

## data table with values to predict
df_predict_ini <- expand.grid(clutch_size = mean(data$clutch_size),
                              egg_position = mean(data$egg_position),
                              male_helpers = mean(data$male_helpers),
                              female_helpers = seq(from = 0, to = 5, 1),
                              rainfall = mean(data$rainfall),
                              temp_above_35 = mean(data$temp_above_35))

## model predictions
df_predict_ini$fit <- predict(model1_plot, 
                              df_predict_ini, 
                              re.form = NA, 
                              type = "response")

## calculate SE for predictions
mm <- model.matrix(~ rainfall +
                     I(rainfall^2) +
                     temp_above_35 +
                     female_helpers +
                     male_helpers +
                     clutch_size +
                     egg_position,
                   data = df_predict_ini)

pvar1 <- diag(mm %*% tcrossprod(vcov(model1_plot),mm))
cmult <- 1 
df_predict_ini <- data.frame(
  df_predict_ini
  , plo = df_predict_ini$fit-cmult*sqrt(pvar1)
  , phi = df_predict_ini$fit+cmult*sqrt(pvar1)
)
df_predict_ini$fit_resp <- df_predict_ini$fit
df_predict_ini$plow_resp <- df_predict_ini$plo
df_predict_ini$phi_resp <- df_predict_ini$phi


##
## model predictions for different number of female helpers, at the mean value of the other variables in the model
model1_plot <- ggplot(df_predict_ini,
                           aes(x = female_helpers, 
                               y = fit_resp/1000)) +
  geom_point(data = data,
             aes(y = egg_volume/1000,
                 x = female_helpers),
             size = 2.5,
             position = position_jitter(width = 0.05),
             alpha = 0.1) +
  geom_line(size = 1, color = "#a50f15") +
  geom_ribbon(aes(ymin = (plow_resp/1000), 
                  ymax = (phi_resp/1000)), 
              alpha = 0.2,
              fill = "#a50f15") +
  
  theme_bw() +
  labs(x = expression(atop('Number of female helpers', " ")),
       y = expression('Egg volume (cm'^3*')')) + 
  theme(plot.margin = unit(c(5,5,30,20), "points"),
        axis.title.x = element_text(family = "Arial", size = 12, vjust = -2),
        axis.title.y = element_text(family = "Arial", size = 12, vjust = +3),
        axis.text.x = element_text(family = "Arial", size = 10),
        axis.text.y = element_text(family = "Arial", size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  scale_x_continuous(limits = c(0,5),
                     labels = seq(0,5,1), 
                     breaks = seq(0,5,1))

ggsave(filename = "./plots/Figure 2a.png", 
       plot = model1_plot,
       device = "png", 
       units = "mm",
       width = 73, 
       height = 89)

