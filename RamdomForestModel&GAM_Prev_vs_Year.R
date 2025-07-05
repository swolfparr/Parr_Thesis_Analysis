#Random forest model for prevalence of cliona by year 

#load library 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readr, tidyverse, ggplot2, dplyr, knitr, gt, car, gridExtra, multcompView, kableExtra, lme4, pwr, mgcv, randomForest, caret, Metrics, gratia, patchwork)


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


#create prevalence column 
cliona_prev <- cliona %>%
  group_by(SampleYear, Transect) %>%
  summarise(
    prevalence = mean(presence, na.rm = TRUE),
    n = n(),  
    .groups = "drop"
  )




#create training test set 
set.seed(123)
trainIndex <- createDataPartition(cliona_prev$prevalence, p = 0.7, list = FALSE)
train_rf <- cliona_prev[trainIndex, ]
test_rf <- cliona_prev[-trainIndex, ]


#train random forest model
rf_model <- randomForest(prevalence ~ SampleYear, data = train_rf, ntree = 500, importance = TRUE)
print(rf_model)

#ntree=500
#mean of squared residuals = 0.0001236339
# % Var explained = 58.74      pretty solid 

#evaluation of model 
pred <- predict(rf_model, newdata = test_rf)
plot(test_rf$SampleYear, pred, main = "Predicted Prevalence by Year")
    #large variation in prevalence per year no exact pattern 


#lets asses the fit then 
rmse <- rmse(test_rf$prevalence, pred)
    #0.015698
    #pretty precise since my prevalence ranges from 0-0.08
rsq <- cor(test_rf$prevalence, pred)^2
    #0.359135
    # explains %36 of the variance meaning moderat generalization 



#ok so lets see what predicted vs observed prevalence is 
pred <- predict(rf_model, newdata = test_rf)

plot(test_rf$prevalence, pred,
     xlab = "Observed Prevalence", ylab = "Predicted Prevalence",
     main = "RF: Observed vs Predicted Prevalence")
abline(0, 1, col = "red")
      #hmmm so its underestimating higher prevalence values 



# ok this is looking good, lets run a GAM and see what it looks like 
    #convert to numeric
train_rf$SampleYear <- as.numeric(as.character(train_rf$SampleYear))
    #run GAM
gam_model <- gam(prevalence ~ s(SampleYear), data = train_rf)
summary(gam_model)
      #edf = 7.24 non-linear trend 
      #F = p<2^-16 highly significant 
      #deviance explained = %53.3 good 
      # adjusted R^2 = 0.483 strong 
      #Cliona prevalence changes significantly over time in a nonlinear fashion, with year explaining a           substantial portion of the trend.

#good now lets plot it 
plot(gam_model, shade = TRUE, seWithMean = TRUE,
     main = "GAM Fit: Cliona Prevalence Over Time",
     xlab = "Year", ylab = "Prevalence")
    #wow so suprised at how this actually is working well in R. love this 


draw(gam_model, residuals = TRUE)  # from {gratia}




#lets look at the derivatives
    #shown by the derivative plot prevalence increased most rapidly 2008, 2014, 2020
smooth_vals <- smooth_estimates(gam_model, select = "s(SampleYear)")

disturbance_years <- c(2005, 2010, 2017, 2020, 2023)

p1 <- ggplot(smooth_vals, aes(x = SampleYear, y = .estimate)) +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = .estimate - .se, ymax = .estimate + .se), alpha = 0.2) +
  geom_vline(xintercept = disturbance_years, linetype = "dashed", color = "black", alpha = 0.5) +
  annotate("text", x = 2005, y = 0.02, label = "2005 Bleaching", angle = 90, vjust = -0.5) +
  annotate("text", x = 2010, y = 0.02, label = "2010 Bleaching", angle = 90, vjust = -0.5) +
  annotate("text", x = 2017, y = 0.02, label = "2017 Hurricanes", angle = 90, vjust = -0.5) +
  annotate("text", x = 2020, y = 0.02, label = "2020 SCTLD", angle = 90, vjust = -0.5) +
  annotate("text", x = 2023, y = 0.02, label = "2023 Bleaching", angle = 90, vjust = -0.5) +
  labs(title = "Smoothed Cliona Prevalence Over Time",
       x = "Year", y = "Estimated Prevalence") +
  theme_minimal()

deriv_vals <- derivatives(gam_model, select = "s(SampleYear)")
names(deriv_vals)
p2 <- ggplot(deriv_vals, aes(x = SampleYear, y = .derivative)) +
  geom_line(color = "darkred") +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = "firebrick", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = disturbance_years, linetype = "dashed", color = "black", alpha = 0.5) +
  annotate("text", x = 2005, y = 0.02, label = "2005 Bleaching", angle = 90, vjust = -0.5) +
  annotate("text", x = 2010, y = 0.02, label = "2010 Bleaching", angle = 90, vjust = -0.5) +
  annotate("text", x = 2017, y = 0.02, label = "2017 Hurricanes", angle = 90, vjust = -0.5) +
  annotate("text", x = 2020, y = 0.02, label = "2020 SCTLD", angle = 90, vjust = -0.5) +
  annotate("text", x = 2023, y = 0.02, label = "2023 Bleaching", angle = 90, vjust = -0.5) +
  labs(title = "Rate of Change in Cliona Prevalence",
       x = "Year", y = "dPrevalence/dYear") +
  theme_minimal()

combined_plot <- p1 / p2  # vertical stack
combined_plot


names(train_rf)
str(train_rf)




###### Adding more predictors ######

cliona_predictors <- read_csv("Cliona.csv")%>%
  select(Period,SampleType, Location, SampleYear, Method, Transect, SPP, Cliona, Spo1, Spo1ID, Spo2ID, Spo2,BL, P, VP, SP, SCTLD)%>%
  mutate(presence=if_else(Cliona>0 | Spo1ID %in% c("CLSP", "BOSP") | Spo2ID %in% c("CLSP", "BOSP"),1,0))%>%
  mutate(presence=replace_na(presence,0))%>%
  filter(Period%in% c("Annual","PostBL")&SampleType%in% c("Permanent","permanent")&Method%in% c("intercept","50 cm belt"))%>%
  filter(!grepl("A",Transect))%>%
  filter(SampleYear>=2005)%>%
  mutate(ID=1:nrow(.),
         SampleYear=as.factor(SampleYear),
         Location=as.factor(Location),
         Transect=as.factor(Transect))

site_impact<- read_csv("SiteMetadata.csv")%>%
  select(Location, ReefComplex, Depth)


# Join site metadata to each coral colony
cliona_full <- left_join(cliona_predictors, site_impact, by = "Location")


PrevWdepth <- cliona_full %>%
  group_by(SampleYear, Transect, Location, Depth, ReefComplex) %>%
  summarise(
    prevalence = mean(presence, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

##Now lets fit the GAM


gam_covariates <- gam(
  prevalence ~ s(as.numeric(as.character(SampleYear))) +
    s(Depth) +
    s(Transect, bs = "re") +
    ReefComplex,
  data = PrevWdepth,
  method = "REML"
)

summary(gam_covariates)
    #Sampl



ggplot(PrevWdepth, aes(x = ReefComplex, y = prevalence)) +
  geom_boxplot() +
  labs(title = "Cliona Prevalence by Reef Complex")










