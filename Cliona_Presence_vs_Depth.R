#Depth Presence vs Depth Logistic Model 
#S. Parr
#5/12/25

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






# Load site depth data (Assuming depth data is in "SiteMetadata.csv")
site_metadata <- read_csv("TCRMP_SiteMetadata.csv") %>%
  select(Location, Depth)  # Ensure it has depth


# Merge depth into cliona data
cliona_depth <- cliona %>%
  left_join(site_metadata, by = "Location") %>%
  filter(!is.na(Depth))  # Drop rows without depth info


# Logistic regression: Cliona presence as a function of depth
glm_depth <- glm(presence ~ Depth,
                 data = cliona_depth,
                 family = binomial)

summary(glm_depth)



# Generate sequence of depth values for prediction
depth_seq <- seq(min(cliona_depth$Depth, na.rm = TRUE),
                 max(cliona_depth$Depth, na.rm = TRUE),
                 length.out = 100)

# Create new data for prediction
pred_data <- data.frame(Depth = depth_seq)
pred_data$predicted <- predict(glm_depth, newdata = pred_data, type = "response")



ggplot(cliona_depth, aes(x = Depth, y = presence)) +
  geom_jitter(height = 0.05, alpha = 0.3, color = "black") +
  geom_line(data = pred_data, aes(x = Depth, y = predicted), color = "blue", size = 1.2) +
  labs(title = "Logistic Regression: Cliona Presence vs. Depth",
       x = "Depth (m)", y = "Predicted Probability of Cliona Presence") +
  theme_minimal()




glm_depth <- glm(presence ~ Depth, data = cliona_depth, family = binomial)



depth_range <- seq(min(cliona_depth$Depth), max(cliona_depth$Depth), length.out = 200)

prediction_df <- data.frame(Depth = depth_range)
prediction_df$Predicted <- predict(glm_depth, newdata = prediction_df, type = "response")


ggplot() +
  geom_point(data = cliona_depth, aes(x = Depth, y = presence),
             color = "red", size = 3, shape = 16, alpha = 0.6) +
  geom_line(data = prediction_df, aes(x = Depth, y = Predicted),
            color = "black", size = 1.2) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(title = "Logistic Regression: Cliona Presence vs. Depth",
       x = "Depth (m)",
       y = "Probability of Presence") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))



# Get depths where at least one presence was detected
depths_with_cliona <- cliona_depth %>%
  group_by(Depth) %>%
  summarise(any_present = any(presence == 1)) %>%
  filter(any_present) %>%
  pull(Depth)

# Filter the dataset: keep only those depths and drop 0s within them
cliona_filtered <- cliona_depth %>%
  filter(Depth %in% depths_with_cliona) %>%
  filter(!(presence == 0))  # Drop 0s from depths where 1s occurred


depths_with_both <- cliona_depth %>%
  group_by(Depth) %>%
  summarise(present = sum(presence == 1),
            absent = sum(presence == 0)) %>%
  filter(present > 0 & absent > 0) %>%
  pull(Depth)

cliona_filtered <- cliona_depth %>%
  filter(Depth %in% depths_with_both)


# Fit logistic regression model on filtered data
glm_filtered <- glm(presence ~ Depth, data = cliona_filtered, family = binomial)

# Generate prediction grid with confidence intervals
depth_seq <- seq(min(cliona_filtered$Depth), max(cliona_filtered$Depth), length.out = 200)
pred_df <- data.frame(Depth = depth_seq)
pred_link <- predict(glm_filtered, newdata = pred_df, type = "link", se.fit = TRUE)

# Convert to probability scale (logit → probability)
pred_df <- pred_df %>%
  mutate(
    fit = plogis(pred_link$fit),
    lower = plogis(pred_link$fit - 1.96 * pred_link$se.fit),
    upper = plogis(pred_link$fit + 1.96 * pred_link$se.fit)
  )

# Plot with jittered binary data and logistic curve + ribbon
ggplot() +
  # Binary data points (jittered)
  geom_jitter(data = cliona_filtered, 
              aes(x = Depth, y = presence), 
              width = 0, height = 0.05, alpha = 0.3, color = "red") +
  
  # Confidence ribbon
  geom_ribbon(data = pred_df, 
              aes(x = Depth, ymin = lower, ymax = upper), 
              fill = "blue", alpha = 0.2) +
  
  # Fitted logistic curve
  geom_line(data = pred_df, 
            aes(x = Depth, y = fit), 
            color = "blue", size = 1.2) +
  
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  
  labs(
    title = "Logistic Regression: Cliona Presence vs. Depth",
    x = "Depth (m)",
    y = "Probability of Presence"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )


















