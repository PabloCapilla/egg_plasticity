###
###
#' 
#' Script for:
#' Mothers in a cooperatively breeding bird increase investment per offspring at the pre-natal stage when they will have more help with post-natal care
#' Capilla-Lasheras et al. 
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
               ggExtra,
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
data <- data %>% 
  mutate(egg_volume_cm = egg_volume/1000)
head(data)

##
##
##### Within and among mother distribution of male and female helpers #####
##
##

# among- and within-mother distributions
helper_distr <- ggplot(data = data %>% 
         select(cent_female_helpers, cent_male_helpers, mean_female_helpers, mean_male_helpers) %>% 
         pivot_longer(cols = 1:4, names_to = "variable", values_to = "helper_value") %>% 
         separate(col = "variable", into = c("variable", "sex", "remove"), sep = "_") %>% 
         mutate(variable = recode(variable, cent = "Within mothers", mean = "Among mothers"),
                sex = recode(sex, female = "Female helpers", male = "Male helpers")),
       aes(x = helper_value)) +
  facet_grid(sex~variable, scales = "free") +
  geom_histogram(aes(y = ..density.., fill = sex),
                 color = NA, 
                 alpha = 0.25,) +
  geom_density(aes(color = sex), 
               size = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 12, family = "Arial"),
        axis.text = element_text(size = 10, family = "Arial"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12, family = "Arial")) +
  scale_x_continuous(labels = c(-3:5), breaks = c(-3:5)) +
  scale_color_manual(values = c("#a50f15", "#4575b4")) +
  scale_fill_manual(values = c("#a50f15", "#4575b4")) +
  labs(y = "Density distribition (i.e., normalised histogram)", x = "Number of helpers")
  
ggsave(filename = "./plots/Figure S6.png", 
       plot = helper_distr,
       device = "png", 
       units = "mm",
       width = 125, 
       height = 125)


##
##
##### 1 - Second model - partitioning variation in female and male helper number effects on egg volume model #####
##
##

# full model partition within and among-mother variation helper number
model2_main_effects_model <- lmer(egg_volume_cm ~ 
                                    poly(rainfall,2)[,1] +
                                    poly(rainfall,2)[,2] +
                                    scale(temp_above_35) +
                                    
                                    ## main effects
                                    cent_male_helpers +
                                    mean_male_helpers +
                                    cent_female_helpers +
                                    mean_female_helpers +
                                    scale(egg_position) +
                                    scale(clutch_size) +
                                    
                                    #random effects
                                    (1|Group) +
                                    (1|mother_ID) +
                                    (1|Season) + 
                                    (1|clutch_ID),
                                  data = data,
                                  na.action = "na.fail",
                                  REML = F)
summary(model2_main_effects_model) # boundary (singular) warning due to group and season variances = 0. Removing these two random effects removes the warning and produces the same estimates
lmerTest::rand(model2_main_effects_model)# LRT for random effects
drop1(update(object = model2_main_effects_model, REML = F), 
      test = "Chisq") # Likelihood-ratio test for each predictor

# formal test of differences between female helper and male helper effects within-mother from the best model containing both
car::linearHypothesis(model2_main_effects_model, "cent_female_helpers - cent_male_helpers = 0")

# map poly coefficients to real data scale
source("../prep_scripts/get_poly_orth_map.R")
gamma <- get_poly_orth_map(poly(data$rainfall, 2))
drop(gamma %*% fixef(model2_main_effects_model)[1:3])




# table with model coefficients
tab_model(model2_main_effects_model,
          file="./tables/Table 2 - RAW.doc",
          pred.labels = c("Intercept", 
                          html("Rainfall<sup>1</sup>"),
                          html("Rainfall<sup>2</sup>"),
                          "Heat waves",
                          html("&Delta;Number of male helpers"),
                          html("&mu;Number of male helpers"),
                          html("&Delta;Number of female helpers"),
                          html("&mu;Number of female helpers"),
                          "Egg position",
                          "Clutch size"),
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

#####

##
##
##### 2 - partitioning variation in female and male helper number effects on egg volume model - AIC tables #####
##
##

# model selection based on AIC
model2_full_partitioning_aic_table <- dredge(update(model2_main_effects_model, 
                                                    control = lmerControl( # ignore message for clean trace - the message has been checked in model above
                                                      check.conv.singular = .makeCC(action = "ignore", 
                                                                                    tol = formals(isSingular)$tol)
                                                    )), 
                                             trace = 3, 
                                             subset = dc("poly(rainfall, 2)[, 1]", "poly(rainfall, 2)[, 2]"), # code to always include main effect when quadratic rainfall effect included
                                             rank = "AIC") 

# AIC results 
model2_d6_subset <- subset(model2_full_partitioning_aic_table, # for Table S1
                           delta < 6)


##
## 
## create table with AIC results
names(model2_d6_subset)
number_variables <- 10
col_names <- names(model2_d6_subset)[1:number_variables]

##
##
##### Code to generate Table #####
##
##

## add each model to the table
# list to store data
list_models_tableS5 <- as.list(NA)
for(m in 1:nrow(model2_d6_subset)){
  
  # template to store data
  tableS5_template <- data.frame(coefficient = names(model2_d6_subset)[1:number_variables]) 
  
  # add model coeffiecients
  model_coef <- data.frame(coefTable(get.models(model2_d6_subset, m)[[1]])) %>% 
    mutate(coefficient = rownames(coefTable(get.models(model2_d6_subset, m)[[1]]))) %>% 
    rename(estimate = Estimate, 
           SE = Std..Error)
  tableS5_00 <- left_join(x = tableS5_template, 
                         y = model_coef %>% select(-df), 
                         by = "coefficient")
  
  ## put table data in right format
  tableS5_01 <- tableS5_00 %>% 
    pivot_wider(names_from = coefficient, values_from = c(estimate, SE)) %>% 
    mutate(k = model2_d6_subset$df[m],
           AIC = model2_d6_subset$AIC[m],
           delta = model2_d6_subset$delta[m]) %>% 
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
             `Clutch size` = `estimate_scale(clutch_size)`,            # not present in models in Table 1
             `Clutch size SE` = `SE_scale(clutch_size)`,
             `Egg position` = `estimate_scale(egg_position)`,
             `Egg position SE` = `SE_scale(egg_position)`,
             `Rainfall1` = `estimate_poly(rainfall, 2)[, 1]`,
             `Rainfall1 SE` = `SE_poly(rainfall, 2)[, 1]`,
             `Rainfall2` = `estimate_poly(rainfall, 2)[, 2]`,
             `Rainfall2 SE` = `SE_poly(rainfall, 2)[, 2]`,
             `Heat waves` = `estimate_scale(temp_above_35)`,
             `Heat waves SE` = `SE_scale(temp_above_35)`,
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
  cols_merge_uncert(col_val = `d Number of helping females`, col_uncert =`d Number of helping females SE`) %>% 
  cols_merge_uncert(col_val = `d Number of helping males`, col_uncert =`d Number of helping males SE`) %>% 
  cols_merge_uncert(col_val = `mean Number of helping females`, col_uncert =`mean Number of helping females SE`) %>% 
  cols_merge_uncert(col_val = `mean Number of helping males`, col_uncert =`mean Number of helping males SE`) %>% 
  cols_merge_uncert(col_val = `Clutch size`, col_uncert =`Clutch size SE`) %>% 
  cols_merge_uncert(col_val = `Egg position`, col_uncert =`Egg position SE`) %>% 
  cols_merge_uncert(col_val = Rainfall1, col_uncert =`Rainfall1 SE`) %>% 
  cols_merge_uncert(col_val = Rainfall2, col_uncert =`Rainfall2 SE`) %>% 
  cols_merge_uncert(col_val = `Heat waves`, col_uncert =`Heat waves SE`) %>% 
  fmt_missing(columns =  c(1:(ncol(tableS5_data)-3), (ncol(tableS5_data)-1):(ncol(tableS5_data))), 
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

# TABLE
tableS5_clean

# save Table (saved in html, then imported in docx to include in manuscript)
tableS5_clean %>%
  gtsave(filename = "./tables/Table S5.html")


##
##
##### Plotting model predictions #####
##
##

# re-fit model 2 to facilitate plotting
model2_plot <- lmer(egg_volume ~ 
                      rainfall +
                      I(rainfall^2) +
                      (temp_above_35) +
                      
                      # main effects
                      cent_male_helpers +
                      mean_male_helpers +
                      cent_female_helpers +
                      mean_female_helpers +
                      (egg_position) +
                      (clutch_size) +
                      
                      # random effects
                      (1|Group) +
                      (1|mother_ID) +
                      (1|Season) + 
                      (1|clutch_ID),
                    data = data, REML = F) # warning due to large values in Rainfall^2 - note delta female effect (which is what is going to be plotted) doesn't vary in comparison with the top model above, with rainfall predictor scaled and without this warning
summary(model2_plot)

## data table with values to predict
df_predict_ini <- expand.grid(clutch_size = mean(data$clutch_size),
                              egg_position = mean(data$egg_position),
                              mean_male_helpers = mean(data$mean_male_helpers),
                              cent_male_helpers = mean(data$cent_male_helpers),
                              mean_female_helpers = mean(data$mean_female_helpers),
                              cent_female_helpers = seq(from = min(data$cent_male_helpers),
                                                        to = max(data$cent_male_helpers), 
                                                        by = 0.1),
                              rainfall = mean(data$rainfall),
                              temp_above_35 = mean(data$temp_above_35))

## predictions
df_predict_ini$fit <- predict(model2_plot, 
                              df_predict_ini, 
                              re.form = NA, 
                              type = "response")

## create SE for predictions
mm <- model.matrix(~ rainfall +
                     I(rainfall^2) +
                     (temp_above_35) +
                     
                     # main effects
                     cent_male_helpers +
                     mean_male_helpers +
                     cent_female_helpers +
                     mean_female_helpers +
                     (egg_position) +
                     (clutch_size),
                   data = df_predict_ini)


pvar1 <- diag(mm %*% tcrossprod(vcov(model2_plot),mm))
cmult <- 1
df_predict_ini <- data.frame(
  df_predict_ini
  , plo = df_predict_ini$fit-cmult*sqrt(pvar1)
  , phi = df_predict_ini$fit+cmult*sqrt(pvar1)
)
df_predict_ini$fit_resp <- df_predict_ini$fit/1000  # change egg volume units from mm3 to cm3
df_predict_ini$plow_resp <- df_predict_ini$plo/1000
df_predict_ini$phi_resp <- df_predict_ini$phi/1000


##
## within-female effect on egg volume, at the mean value of the other variables in the model
model2_plot <- ggplot(df_predict_ini,
                           aes(x = cent_female_helpers, 
                               y = fit_resp)) +
  geom_point(data = data,
             aes(y = egg_volume/1000,
                 x = cent_female_helpers),
             size = 2.5,
             alpha = 0.1) +
  geom_line(size = 1, color = "#a50f15") +
  geom_ribbon(aes(ymin = (plow_resp), 
                  ymax = (phi_resp)), 
              alpha = 0.2,
              fill = "#a50f15") +
  theme_bw() +
  labs(x = expression(atop(paste(Delta, " female helper number", sep = " "),
                           "(within-mother variation)")), 
       y = expression('Egg volume (cm'^3*')')) + 
  theme(plot.margin = unit(c(5,5,30,20), "points"),
        axis.title.x = element_text(family = "Arial", size = 12, vjust = -2),
        axis.title.y = element_text(family = "Arial", size = 12, vjust = +3),
        axis.text.x = element_text(family = "Arial", size = 10),
        axis.text.y = element_text(family = "Arial", size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  scale_x_continuous(limits = c(-3,+4),
                     labels = seq(-3,4,1), 
                     breaks = seq(-3,4,1))

ggsave(filename = "./plots/Figure 2b.png", 
       plot = model3_best_plot,
       device = "png", 
       units = "mm",
       width = 73, 
       height = 89)
