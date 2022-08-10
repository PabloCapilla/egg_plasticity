data_old <- read.csv("./data/egg_volume_dataset.csv") 
data_new <- read.csv("./data/zArc/egg_volume_dataset.csv") 


df <- left_join(x = data_old %>% 
                  select(clutch_ID, female_helpers, male_helpers, egg_volume_old = egg_volume) %>% 
                  group_by(clutch_ID) %>% 
                  summarise(egg_volume_old = mean(egg_volume_old),
                            female_old = mean(female_helpers),
                            male_old = mean(male_helpers)),
                y = data_new %>% 
                  select(clutch_ID, female_helpers, male_helpers, egg_volume_new = egg_volume) %>% 
                  group_by(clutch_ID) %>% 
                  summarise(egg_volume_new = mean(egg_volume_new),
                            female_new = mean(female_helpers),
                            male_new = mean(male_helpers)),
                by = "clutch_ID")
head(df)

summary(df$egg_volume_old - df$egg_volume_new)
summary(df$female_old - df$female_new)
summary(df$male_old - df$male_new)


df %>% 
  mutate(diff_female = female_old - female_new) %>% 
  filter(diff_female != 0)

df %>% 
  mutate(diff_male = male_old - male_new) %>% 
  filter(diff_male != 0)
