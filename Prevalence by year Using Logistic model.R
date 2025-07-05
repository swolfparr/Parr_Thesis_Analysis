# Spencer Parr
# 3/30/25
# Prevalence by year Using Logistic model (family Binomial)



if (!require("pacman")) install.packages("pacman")
pacman::p_load(readr, tidyverse, ggplot2, dplyr, knitr, gt, car, gridExtra,
               multcompView, kableExtra, lme4, pwr,mgcv)



raw_clio <- read_csv("Cliona.csv") #load data

#previous data manipulation

#1
cliona <- read_csv("Cliona.csv")%>%
  select(Period,SampleType, Location, SampleYear, Method, Transect, SPP, Cliona, Spo1, Spo1ID, Spo2ID, Spo2)%>%
  mutate(presence=if_else(Cliona>0 | Spo1ID %in% c("CLSP", "BOSP") | Spo2ID %in% c("CLSP", "BOSP"),1,0))%>%
  mutate(presence=replace_na(presence,0))%>%
  filter(Period%in% c("Annual","PostBL")&SampleType%in% c("Permanent","permanent")&Method%in% c("intercept","50 cm belt"))%>%
  filter(!grepl("A",Transect))%>%
  filter(SampleYear>=2005)%>%
  mutate(ID=1:nrow(.),
         SampleYear=as.factor(SampleYear),
         Location=as.factor(Location),
         Transect=as.factor(Transect))


############################################
Prev_by_year<- cliona %>%
  group_by(SampleYear, SPP, ID, Transect) %>%
  summarise(prev = mean(presence), freq = sum(presence))%>%
  ungroup()%>%
  mutate(ID= 1:nrow(.))


Prev_Year_noZeros <- Prev_by_year %>% 
  group_by(SampleYear, ID, SPP, Transect) %>%
  filter(prev>0)


Prev_grouped<- Prev_Year_noZeros %>%
  group_by(SampleYear) %>% 
  summarise(total=sum(prev))


ggplot(Prev_grouped, aes(x = SampleYear, y = total)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Linear Regression of Prevalence by Year",
       x = "Year",
       y = "Prevalence")


ggplot(Prev_grouped, aes(x = as.numeric(SampleYear), y = total)) +
  geom_point(size = 3, color = "steelblue") +
  geom_smooth(method = "loess", se = FALSE, color = "purple")+
  labs(title = "Cliona Observations Over Time",
       x = "Year",
       y = "Total Cliona Observations") +
  theme_minimal()



model_species_year <- glm(prev ~ SampleYear + SPP, data = Prev_Year_noZeros, family = binomial)
summary(model_species_year)











##############################################################
Locationprev<- cliona%>%
  group_by(ID)%>%
  summarise(prev= mean(presence),freq=sum(presence), total=n())


Prevelance_ID<- cliona %>%
  group_by(Location, SampleYear, Transect,SPP, ID) %>%
  summarise(prev = mean(presence), freq = sum(presence), total = n())%>%
  ungroup()%>%
  mutate(ID= 1:nrow(.))



Prev_by_year<- cliona %>%
  group_by(SampleYear, SPP, ID, Transect) %>%
  summarise(prev = mean(presence), freq = sum(presence))%>%
  ungroup()%>%
  mutate(ID= 1:nrow(.))


Prev_Year_noZeros <- Prev_by_year %>% 
 group_by(SampleYear, ID, SPP, Transect) %>%
  filter(prev>0)
  

Prev_grouped<- Prev_Year_noZeros %>%
  group_by(SampleYear) %>% 
  summarise(total=sum(prev))


ggplot(Prev_grouped, aes(x = SampleYear, y = total)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Linear Regression of Prevalence by Year",
       x = "Year",
       y = "Prevalence")

only_year<- cliona %>%
  group_by(SampleYear,Transect)%>%
  summarise(prev = mean(presence), total = n())%>%
  ungroup()%>%
  mutate(ID= 1:nrow(.))

only_year <- only_year %>%
  mutate(prev = ifelse(prev > 0, 1, 0))




# Run a linear regression
Lreg <- lm(prev ~ SampleYear+Transect, data = Prev_by_year)

# View summary of model
summary(Lreg)


ggplot(Lreg, aes(x = SampleYear, y = prev)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Linear Regression of Prevalence by Year",
       x = "Year",
       y = "Prevalence")




# Full model
model_full <- glm(prev ~ SampleYear+Transect, data = only_year, family = binomial)

summary(model_full)

model_reduced <- glm(prev ~ SampleYear, data=only_year, family=binomial)
model_reduced
# Compare using Chi-square (likelihood ratio test)
anova(model_reduced, model_full, test = "Chisq")

model_SPP <- glm(prev ~ SampleYear+SPP, data = Prev_by_year, family = binomial)
model_SPP




ggplot(Prev_by_year, aes(x = as.numeric(as.character(SampleYear)), y = prev)) +
  geom_jitter(height = 0.05, alpha = 0.1) +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"),
              se = FALSE) +
  labs(x = "Sample Year", y = "Prevalence") +
  theme_minimal()




# Plot average prevalence by species
species_summary <- Prev_by_year %>%
  group_by(SPP) %>%
  summarise(
    Cliona_Prevalence = mean(prev),
    n = n()
  ) %>%
  arrange(desc(Cliona_Prevalence))

ggplot(species_summary, aes(x = reorder(SPP, -Cliona_Prevalence), y = Cliona_Prevalence)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Cliona Prevalence by Coral Species (SPP)",
       x = "Coral Species", y = "Prevalence") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))



















