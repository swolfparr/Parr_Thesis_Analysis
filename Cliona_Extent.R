#CLIONA EXTENT CODE 
#NOV 2024
#SPENCER PARR


library(tidyverse)
library(WRS2)

getwd()
clionaext <- read_csv("Cliona.csv")%>%
  select(Period,SampleType, Location, SampleYear, Method, Transect, SPP, Cliona, Spo1, Spo1ID, Spo2ID, Spo2)%>%
  mutate( CLSP = case_when( Spo1ID == "CLSP" ~ Spo1, Spo2ID == "CLSP" ~ Spo2, TRUE ~ NA_real_),
    BOSP = case_when( Spo1ID == "BOSP" ~ Spo1, Spo2ID == "BOSP" ~ Spo2,
      TRUE ~ NA_real_)) %>%
  select(-Spo1, -Spo1ID, -Spo2, -Spo2ID)%>%
  filter(Period%in% c("Annual","PostBL")&SampleType%in% c("Permanent","permanent")&Method%in% c("intercept","50 cm belt"))%>%
  filter(!grepl("A", "R", "L",Transect))%>%
  filter(SampleYear>=2005)%>%
  mutate(ID=1:nrow(.),
         SampleYear=as.factor(SampleYear),
         Location=as.factor(Location),
         Transect=as.factor(Transect))%>%
  pivot_longer(cols = c(CLSP, BOSP, Cliona), # Combine the columns into a long format
               names_to = "SpongeType", 
               values_to = "Extent")




#Has the extent of Cliona sponge coverage increased or decreased over time? Is there a pattern in prevalence by year?
    #Observing trends over time can provide insights into whether Cliona sponge coverage is growing or declining across the dataset.

ext_time<-clionaext

ext_time <- ext_time %>%
  mutate(SampleYear = as.numeric(as.character(SampleYear)))



#Calculate the total unique corals per year from the original data.


#Calculate Unique Coral Counts
unique_coral_counts <- clionaext %>%
  group_by(SampleYear) %>%
  summarise(total_unique_corals = n_distinct(ID)) 

unique_coral_counts <- unique_coral_counts%>%
  mutate(SampleYear = as.numeric(as.character(SampleYear)))


#Join Unique Counts to ext_time
ext_time_with_counts <- ext_time %>%
  left_join(unique_coral_counts, by = "SampleYear")



#calculate the weighted average using total_unique_corals for each SampleYear.
weighted_summary <- ext_time_with_counts %>%
  group_by(SampleYear) %>%
  summarise(
    total_extent = sum(Extent, na.rm = TRUE),
    count_observations = first(total_unique_corals),  # Use `first` to avoid altering the count
    weighted_avg_extent = total_extent / count_observations
  )


#create a clean time-series plot to visualize trends
ggplot(weighted_summary, aes(x = SampleYear, y = weighted_avg_extent)) +
  geom_line(color = "blue") +
  geom_point(size = 2) +
  labs(
    title = "Weighted Average Cliona Sponge Coverage Over Time",
    x = "Year",
    y = "Weighted Average Extent of Coverage"
  ) +
  theme_minimal()




#How does the extent of coverage differ between coral species each year? Are there species that experience fluctuating or consistent coverage levels over time?
      #Examining how coverage on specific coral species varies over time could reveal interesting temporal patterns and interactions.

spp_summary <- ext_time_with_counts %>%
  group_by(SampleYear, SPP) %>%
  summarise(
    total_extent = sum(Extent, na.rm = TRUE),
    count_observations = first(total_unique_corals),
    weighted_avg_extent = total_extent / count_observations
  )

ggplot(spp_summary, aes(x = SampleYear, y = weighted_avg_extent, color = SPP)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Weighted Average Cliona Sponge Coverage by Coral Species Over Time",
    x = "Year",
    y = "Weighted Average Extent of Coverage"
  ) +
  theme_minimal()



#Determine if Cliona coverage significantly differs between locations.

# Summarize Cliona extent by location
location_summary <- ext_time_with_counts %>%
  group_by(Location) %>%
  summarise(
    mean_extent = mean(Extent, na.rm = TRUE),
    weighted_avg_extent = sum(Extent, na.rm = TRUE) / first(total_unique_corals)
  )

# Plot differences between locations
ggplot(location_summary, aes(x = reorder(Location, -mean_extent), y = mean_extent)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(
    title = "Mean Cliona Sponge Coverage by Location",
    x = "Location",
    y = "Mean Extent of Coverage"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Linear regression for trend analysis
trend_model <- lm(weighted_avg_extent ~ SampleYear, data = weighted_summary)
summary(trend_model)

# Robust trend test (if needed)
library(trend)
mann_kendall <- mk.test(weighted_summary$weighted_avg_extent)
print(mann_kendall)


# Example: Correlation between extent and environmental variable (e.g., temperature)
env_data <- read_csv("EnvironmentalData.csv")  # Load environmental data
ext_with_env <- ext_time_with_counts %>%
  left_join(env_data, by = c("Location", "SampleYear"))

# Correlation
correlation <- cor.test(ext_with_env$Extent, ext_with_env$Temperature, use = "complete.obs")
print(correlation)

# Scatterplot
ggplot(ext_with_env, aes(x = Temperature, y = Extent)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "red") +
  labs(
    title = "Relationship Between Temperature and Cliona Coverage",
    x = "Temperature (°C)",
    y = "Cliona Extent"
  ) +
  theme_minimal()











































#How does the extent of Cliona sponge coverage vary across different coral species?



#What are the average and maximum coverage extents for each coral species? Are some coral species more heavily affected than others?
#This can help identify which coral species might be more susceptible or commonly associated with Cliona sponge coverage.

avg_per_SPP<- clionaext%>%
  group_by(SPP) %>%
  summarise(
    avg_extent = mean(Extent, na.rm = TRUE),
    max_extent = max(Extent, na.rm = TRUE),
    total_extent = sum(Extent, na.rm = TRUE)
  ) %>%
  arrange(desc(avg_extent))








spp_yr<-clionaext
spp_yr %>%
  group_by(SPP, SampleYear) %>%
  summarise(avg_extent = mean(Extent, na.rm = TRUE)) %>%
  ggplot(aes(x = SampleYear, y = avg_extent, color = SPP)) +
  geom_line() +
  labs(title = "Average Cliona Sponge Coverage Over Time by Coral Species",
       x = "Year",
       y = "Average Extent of Coverage",
       color = "Coral Species")


# Calculate the weighted average of Extent by year, with weights based on observation counts
weighted_avg_data <- ext_time %>%
  group_by(SampleYear) %>%
  summarise(
    weighted_avg_extent = sum(Extent, na.rm = TRUE) / n(),  # weighted by count
    total_extent = sum(Extent, na.rm = TRUE),                # total extent for context
    count_observations = n()                                 # number of observations per year
  )



#This will give you the average percent interaction per transect for each year.
# Grouping by Transect, Year, InteractionType, and Location to aggregate data
clionagrouped <- clionalong %>%
  group_by(Location, SampleYear, Transect, InteractionType) %>%
  summarise(avg_interaction = mean(PercentInteraction, na.rm = TRUE),
            total_observations = n())








