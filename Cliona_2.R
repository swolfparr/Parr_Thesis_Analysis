##Cliona Prevalence Analysis 
#Spencer Parr 
#1/12/2025
##Consists of data manipulation, reef complex data analysis, prevalence analysis for locations, GLM, ANOVA, LME



library(tidyverse)
library(WRS2)

getwd()
cliona <- read_csv("Cliona.csv")%>%
  select(Period,SampleType, Location, SampleYear, Method, Transect, SPP, Cliona, Spo1, Spo1ID, Spo2ID, Spo2, SampleDate, SampleMonth)%>%
  mutate(presence=if_else(Cliona>0 | Spo1ID %in% c("CLSP", "BOSP") | Spo2ID %in% c("CLSP", "BOSP"),1,0))%>%
  mutate(presence=replace_na(presence,0))%>%
  filter(Period%in% c("Annual","PostBL")&SampleType%in% c("Permanent","permanent")&Method%in% c("intercept","50 cm belt"))%>%
  filter(!grepl("A",Transect))%>%
  filter(SampleYear>=2005)%>%
  mutate(ID=1:nrow(.),
         SampleYear=as.factor(SampleYear),
         Location=as.factor(Location),
         Transect=as.factor(Transect))
##############################################################







##############################################################
Locationprev<- cliona%>%
  group_by(Location)%>%
  summarise(prev= mean(presence),freq=sum(presence), total=n())


transectprev <- cliona %>%
  group_by(Location, SampleYear, Transect, SampleMonth) %>%
  summarise(prev = mean(presence), freq = sum(presence), total = n())%>%
  ungroup()%>%
  mutate(ID= 1:nrow(.))

#############################################################################
#Using Reef Complex to show prevalence 

reef_complex <- read_csv("SiteMetadata.csv")

#Merge the reef_complex data with your Cliona prevalence data (cliona_data) based on the Location column.
Reef_Types <- transectprev %>%
  left_join(reef_complex, by = "Location")

# Summarize prevalence by location
reef_sum <- Reef_Types %>%
  group_by(Location, ReefComplex) %>%
  summarise(
    mean_prevalence = mean(prev, na.rm = TRUE)  # Calculate mean prevalence
  )%>%
  mutate( ReefComplex = factor(ReefComplex, levels = c("Nearshore", "Offshore", "Mesophotic")))
 
#Visulaize using a Box Plot 
ggplot(reef_sum, aes(x = ReefComplex, y = mean_prevalence, fill = ReefComplex)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +  # Add jitter for data points
  labs(
    title = "Cliona Sponge Prevalence Across Reef Complexes",
    x = "Reef Complex",
    y = "Mean Prevalence"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

####Run ANOVA on this#####
    #The ANOVA shows that the mean prevalence of Cliona sponges differs significantly between reef complexes (𝐹= 6.625 F=6.625,𝑝= 0.00414 p=0.00414). However, it does not specify which groups differ; the Tukey HSD test provides those pairwise comparisons.

# Run ANOVA
anova_result <- aov(mean_prevalence ~ ReefComplex, data = reef_sum)

# Summarize the ANOVA result
summary(anova_result)

# Perform Tukey's HSD test
   #Significant Difference:
 #Mesophotic vs. Nearshore: There is a significant reduction in mean prevalence in Mesophotic sites compared to Nearshore sites (p = 0.0030).
    #No Significant Difference:
 # Offshore vs. Nearshore: No significant difference in mean prevalence (p = 0.1506). 
 # Mesophotic vs. Offshore: No significant difference in mean prevalence (p = 0.2067).
tukey_result <- TukeyHSD(anova_result)

# Print the result
print(tukey_result)

# Visualize the Tukey result
plot(tukey_result)

# Check normality of residuals
    #W = 0.965: This is the Shapiro-Wilk statistic. A value closer to 1 suggests that the data are more consistent with a normal distribution
    #p-value = 0.3555: Since the p-value is greater than the common significance threshold (0.05), you fail to reject the null hypothesis of normality. This means the residuals are consistent with a normal distribution, satisfying the normality assumption for ANOVA.
shapiro.test(residuals(anova_result))




# Generate individual plots
qq_plot <- ggplot(data.frame(residuals = residuals(anova_result)), aes(sample = residuals)) +
  stat_qq() + stat_qq_line(color = "red") + 
  labs(title = "Q-Q Plot of Residuals")

hist_plot <- ggplot(data.frame(residuals = residuals(anova_result)), aes(x = residuals)) +
  geom_histogram(binwidth = 0.005, fill = "lightblue", color = "black") + 
  labs(title = "Histogram of Residuals", x = "Residuals", y = "Frequency") +
  theme_minimal()

density_plot <- ggplot(data.frame(residuals = residuals(anova_result)), aes(x = residuals)) +
  geom_density(color = "blue", fill = "lightblue", alpha = 0.5) +
  stat_function(fun = dnorm, args = list(mean = mean(residuals(anova_result)), 
                                         sd = sd(residuals(anova_result))),
                color = "red", linetype = "dashed") +
  labs(title = "Density Plot with Normal Curve", x = "Residuals", y = "Density")

# Arrange plots
grid.arrange(qq_plot, hist_plot, density_plot, nrow = 1)




# Check homogeneity of variances
    #The Levene's Test result suggests that the p-value is 0.1111, which is greater than the usual significance threshold of 0.05. This means you fail to reject the null hypothesis of homogeneity of variances. Therefore, the assumption of equal variances across reef complexes is satisfied for your ANOVA.
    
library(car)
leveneTest(mean_prevalence ~ ReefComplex, data = reef_sum)

#Normality of Residuals: ✅ The Shapiro-Wilk test confirms that the residuals are normally distributed (p-value > 0.05).
##Homogeneity of Variances: ✅ The Levene's Test shows no significant differences in variances across groups (p-value = 0.1111 > 0.05).

####visualise the resutls#######

# Extract Tukey test results
tukey_result <- TukeyHSD(anova_result, "ReefComplex")

# Convert Tukey results to group letters
library(multcompView)
tukey_labels <- multcompLetters(tukey_result$ReefComplex[, "p adj"])$Letters

# Convert to a data frame for merging
tukey_labels <- data.frame(
  ReefComplex = names(tukey_labels),
  group = tukey_labels
)


# Merge tukey_labels with your reef_sum data
reef_sum <- reef_sum %>%
  left_join(tukey_labels, by = "ReefComplex")



# Plot with Tukey labels
ggplot(reef_sum, aes(x = ReefComplex, y = mean_prevalence, fill = ReefComplex)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  geom_text(
    aes(label = group, x = ReefComplex, y = max(mean_prevalence) + 0.01),
    data = reef_sum, inherit.aes = FALSE, color = "black", size = 4, vjust = -0.5
  ) +
  labs(
    title = "Cliona Sponge Prevalence Across Reef Complexes",
    x = "Reef Complex",
    y = "Mean Prevalence"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )













  
reef_summary <- Reef_Types %>%
  group_by(ReefComplex, Location, SampleYear) %>%
  summarise(
    mean_prevalence = mean(prev, na.rm = TRUE),
    lower_quartile = quantile(prev, 0.25, na.rm = TRUE),
    upper_quartile = quantile(prev, 0.75, na.rm = TRUE),
    iqr = upper_quartile - lower_quartile
  ) %>%
  ungroup()

ggplot(reef_summary, aes(x = Location, y = mean_prevalence, fill = ReefComplex)) +
  # Boxplot-like display using error bars for IQR
  geom_errorbar(aes(ymin = lower_quartile, ymax = upper_quartile), width = 0.2) +
  # Add median points
  geom_point(shape = 21, size = 3, color = "black") +
  # Outliers as jittered points
  geom_jitter(data = Reef_Types, aes(x = Location, y = prev),
              shape = 16, size = 1, alpha = 0.5, color = "black") +
  # Facet for depth categories
  facet_wrap(~ ReefComplex, nrow = 3, scales = "free_x") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Cliona Sponge Prevalence by Depth and Station",
    x = "Reef Sites",
    y = "Median Cliona Prevalence (%)"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "none"
  )


cliona_depth_summary <- reef_summary %>%
  group_by(ReefComplex, Location) %>% 
  summarise(
    median_prevalence = median(mean_prevalence, na.rm = TRUE),
    lower_quartile = quantile(mean_prevalence, 0.25, na.rm = TRUE),
    upper_quartile = quantile(mean_prevalence, 0.75, na.rm = TRUE),
    outliers = list(mean_prevalence[mean_prevalence < quantile(mean_prevalence, 0.25) - 1.5 * IQR(mean_prevalence) |
                                      mean_prevalence > quantile(mean_prevalence, 0.75) + 1.5 * IQR(mean_prevalence)])
  )


ggplot(cliona_depth_summary, aes(x = Location, y = median_prevalence, fill = ReefComplex)) +
  geom_boxplot(aes(lower = lower_quartile, upper = upper_quartile, middle = median_prevalence,
                   ymin = lower_quartile - 1.5 * IQR(median_prevalence), 
                   ymax = upper_quartile + 1.5 * IQR(median_prevalence)),
               stat = "identity", position = position_dodge(width = 0.8)) +
  geom_point(data = cliona_depth_summary %>% unnest(cols = c(outliers)), 
             aes(x = Location, y = outliers), shape = 21, color = "black", fill = "red", size = 2) +
  facet_wrap(~ ReefComplex, scales = "free", nrow = 1) +
  scale_fill_manual(values = c("Nearshore" = "coral", "Offshore" = "skyblue", "Mesophotic" = "lightgreen")) +
  labs(
    title = "Median Cliona Sponge Prevalence by Reef Complex",
    x = "Reef Sites",
    y = "Median Prevalence (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )

###########plot monthly prevalence trends across depth categories (ReefComplex)##########
Reef_month <- Reef_Types %>%
  group_by(ReefComplex, SampleYear, SampleMonth) %>%
  summarise(
    mean_prevalence = mean(prev, na.rm = TRUE),  # Use the prevalence column
    sd_prevalence = sd(prev, na.rm = TRUE),      # Standard deviation
    n = n()                                            # Number of observations
  ) %>%
  mutate(
    se_prevalence = sd_prevalence / sqrt(n)            # Standard error
  )%>%
  mutate(
    SampleMonth = as.numeric(SampleMonth), # Ensure months are numeric
    MonthName = month.abb[SampleMonth],    # Convert numeric months to abbreviations
    YearMonth = factor(paste(SampleYear, MonthName, sep = "-"), 
                       levels = unique(paste(SampleYear, MonthName, sep = "-")))
  )


ggplot(Reef_month, aes(x = YearMonth, y = mean_prevalence, color = ReefComplex, group = ReefComplex)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = mean_prevalence - se_prevalence, ymax = mean_prevalence + se_prevalence),
    width = 0.2, size = 0.5
  ) +
  labs(
    title = "Cliona Sponge Prevalence by Depth and Time",
    x = "Year and Month",
    y = "Cliona Sponge Prevalence (%)",
    color = "Reef Complex"
  ) +
  scale_x_discrete(
    breaks = unique(Reef_month$YearMonth),
    labels = function(x) {
      # Extract the month and year components
      months <- sub("^[0-9]+-", "", x)  # Get the month (e.g., "Jan")
      years <- ifelse(grepl("-Jan", x), sub("-.*$", "", x), "")  # Get the year only for January
      paste0(months, ifelse(years != "", paste0("\n", years), ""))  # Add the year on a new line below January
    }
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.minor.x = element_blank()
  )









ggplot(Reef_month, aes(x = interaction(SampleYear, SampleMonth, sep = " "), y = mean_prevalence, color = ReefComplex, group = ReefComplex)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_prevalence - se_prevalence, ymax = mean_prevalence + se_prevalence), 
                width = 0.2, size = 0.5) +
  labs(
    title = "Cliona Sponge Prevalence by Depth and Time",
    x = "Year and Month",
    y = "Cliona Sponge Prevalence (%)",
    color = "Reef Complex"
  ) +
  theme_minimal() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +  # Ensures clear axis labels
  scale_color_manual(values = c("gold", "green", "blue"))  # Match colors to your example




###############################################################################
# this is the standardized prevalence For all years by Location

#This code calculates, for each location, the average prevalence (proportion) of Cliona sponges observed, the total count of sponge presences, and the total number of observations recorded.
yearsprev <- cliona %>% 
  group_by(Location) %>% 
  summarise(prev = mean(presence), freq = sum(presence), total = n())

# Create a highlight column that marks specific locations
yearsprev <- yearsprev %>%
  mutate(highlight = ifelse(Location %in% c("Brewers Bay", "Magens Bay", "Botany Bay", "Coculus Rock", "Black Point", "Savana"), "Chosen Locations", "Unchosen"))

    #this Highlights promanant sites 
ggplot(yearsprev, aes(x = reorder(Location, -prev), y = prev, fill = highlight)) +
  geom_bar(stat = "identity") +
  labs(title = "Prevalence of Cliona Sponges by Location",
       x = "Location",
       y = "Prevalence") +
  scale_fill_manual(values = c("Chosen Locations" = "red", "Unchosen" = "black")) +  # Assign colors
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

#colored version of Prevaence by Location 
ggplot(yearsprev, aes(x = reorder(Location, -prev), y = prev, fill = Location)) +
  geom_bar(stat = "identity") +
  labs(title = "Prevalence of Cliona Sponges by Location",
       x = "Location",
       y = "Prevalence") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none")  # Hide the legend since it's redundant

#############################################################################

#This code calculates the average prevalence (proportion) of Cliona sponges, the total count of sponge presences, and the total number of observations for each unique combination of Location, SampleYear, and Transect. This provides a detailed summary of sponge presence across different sites and years.

transectprev <- cliona %>%
  group_by(Location, SampleYear, Transect) %>%
  summarise(prev = mean(presence), freq = sum(presence), total = n())%>%
  ungroup()%>%
  mutate(ID= 1:nrow(.))
 
#calculate mean prevalence of cliona for each location
location_summary <- transectprev %>%
  group_by(Location) %>%
  summarise(mean_prev = mean(prev, na.rm = TRUE),
            sd_prev = sd(prev, na.rm = TRUE),
            total_records = n())

###########This ANOVA WORKS!!!######################
# ANOVA (if normality and homogeneity of variance assumptions are met):
  
#   ANOVA can test for significant differences in prev across different Location levels.
# Pairwise comparisons (e.g., Tukey’s test) can follow if ANOVA shows significance.

anova_result <- aov(prev ~ Location, data = transectprev)
summary(anova_result)

par(mfrow=c(2,2))
plot(anova_result)
par(mfrow=c(1,1))


two_way_plot <- ggplot(transectprev, aes(x = Location, y = prev)) +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = 0.2) +      # Error bars for mean and SE
  stat_summary(fun.data = 'mean_se', geom = 'pointrange') +                 # Mean points with error bars
  geom_jitter(data = transectprev, aes(x = Location, y = prev), width = 0.2, alpha = 0.3) +  # Jittered points
  labs(title = "Prevalence of Cliona Sponges by Location",
       x = "Location",
       y = "Prevalence") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 0.2)  # Optional: Adjust y-axis if most values are low

two_way_plot
##Summary of resutls 
    #The ANOVA results indicate a statistically significant effect of location on the prevalence of Cliona sponges. With 32 degrees of freedom associated with Location and 2787 degrees of freedom associated with the residuals, the analysis reveals that the variation in prevalence between locations is substantial. The Sum of Squares for Location is 1.648, and the residual Sum of Squares is 10.196. The Mean Square for Location is calculated as 0.05152, compared to the residual Mean Square of 0.00366. The F-value for Location is 14.08, which is the ratio of the Mean Square for Location to the residual Mean Square. This high F-value suggests that the differences in prevalence across locations are much greater than the within-location variability.

    #The p-value associated with this F-value is less than 2e-16, which is extremely small and highly significant. This result indicates that the differences in prevalence between locations are statistically significant at any conventional significance level (e.g., 0.05, 0.01, 0.001). Therefore, we can confidently reject the null hypothesis, concluding that the prevalence of Cliona sponges is not uniform across locations. These findings suggest that some locations have significantly higher or lower prevalence than others, warranting further investigation to identify specific locations with notably high or low prevalence values. Further pairwise comparisons, such as Tukey’s HSD test, would help determine which specific locations differ in terms of prevalence.

#####Run a Tukey HSD Test##########
# 1. Run Tukey HSD
tukey_result <- TukeyHSD(aov(prev ~ Location, data = transectprev))
tukey_df <- as.data.frame(tukey_result$Location)

# 2. Subset significant comparisons
significant_results <- subset(tukey_df, `p adj` < 0.05)

# 3. Extract and separate rownames
if (nrow(significant_results) > 0) {
  significant_results <- significant_results %>%
    tibble::rownames_to_column("comparison") %>%
    tidyr::separate(comparison, into = c("Location1", "Location2"), sep = "-") %>%
    mutate(Higher_Location = ifelse(diff > 0, Location1, Location2))
  
  # 4. Count which locations had higher prevalence
  summary_high_prevalence <- table(significant_results$Higher_Location)
  summary_high_prevalence <- sort(summary_high_prevalence, decreasing = TRUE)
  
  print(summary_high_prevalence)
} else {
  message("⚠️ No significant differences found in Tukey HSD.")
}




# Run the Tukey test and get CLD letters for Location
tukey_result <- TukeyHSD(anova_result)

# Convert Tukey results into compact letter display (CLD)
tukey_groups <- multcompLetters4(anova_result, tukey_result)

# Extract the letters and create a data frame
cld_df <- data.frame(Location = names(tukey_groups$Location),
                     cld = as.character(tukey_groups$Location))

# Extract the 'Letters' component from the multcompView output
cld_letters <- tukey_groups$Location$Letters

# Convert to a data frame
cld_df <- data.frame(Location = names(cld_letters),
                     cld = as.character(cld_letters))
# Calculate mean prevalence per location for plotting (if needed)
location_means <- transectprev %>%
  group_by(Location) %>%
  summarise(mean_prev = mean(prev, na.rm = TRUE))

# Ensure consistent naming (lowercase and trimmed whitespace) if needed
cld_df$Location <- tolower(trimws(cld_df$Location))
location_means$Location <- tolower(trimws(location_means$Location))

# Join the CLD labels to location_means
location_means <- location_means %>%
  left_join(cld_df, by = "Location")
# Verify the join results
print(location_means)  # Check if 'cld' column has labels without NA

# Ensure both data frames have Location as a factor with the same levels
cliona$Location <- factor(cliona$Location, levels = unique(cliona$Location))
location_means$Location <- factor(location_means$Location, levels = levels(cliona$Location))

# Plot with labels aligned to each location
two_way_plot <- ggplot(cliona_data, aes(x = Location, y = prev)) +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = 0.2) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange') +
  geom_jitter(width = 0.2, alpha = 0.3) +
  geom_text(data = location_means, aes(x = Location, y = mean_prev + 0.05, label = cld), 
            color = "red", fontface = "bold", size = 5) +  # Increased y-offset and size for readability
  labs(title = "Prevalence of Cliona Sponges by Location with Significance Labels",
       x = "Location",
       y = "Prevalence") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

two_way_plot



location_means <- transectprev %>%
  group_by(Location) %>%
  summarise(mean_prev = mean(prev))















##################################
# Using the summary data from Step 1
ggplot(location_summary, aes(x = reorder(Location, -mean_prev), y = mean_prev)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Mean Prevalence of Cliona Sponges by Location",
       x = "Location",
       y = "Mean Prevalence") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))













#This code creates a faceted line plot showing the prevalence of Cliona sponges over years for each transect, with separate panels for each location, allowing for easy comparison of trends across time and locations. The x-axis displays distinct years, and each transect is color-coded for clarity.
      
    # Assuming your dataframe is called transectprev 
 ggplot(transectprev, aes(x = SampleYear, y = prev, color = factor(Transect), group = Transect)) +
  geom_line() +  # Draws lines between points
  geom_point() +  # Adds points for each data point
  facet_wrap(~Location, scales = "free_y") +  # Creates separate plots for each Location
  labs(title = "Prevalence Over Years by Transect",
       x = "Year",
       y = "Prevalence",
       color = "Transect") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate the x-axis labels for better readability

#This code creates a line plot to display the trend of Cliona sponge prevalence over the sample years in the transectprev dataset. The x-axis represents SampleYear, while the y-axis shows the prev (prevalence) values. A blue line (geom_line()) connects these prevalence values across years, representing the overall trend in prevalence, while individual data points are highlighted with geom_point(size = 2). Labels for the title and axes provide context, and theme_minimal() is applied for a clean and simple visual style. This plot helps identify any shifts in sponge prevalence over time.

#plot the prevalence across locations or sample years
ggplot(transectprev, aes(x = SampleYear, y = prev)) +
  geom_line(group = 1, color = "blue") +
  geom_point(size = 2) +
  labs(title = "Prevalence of Cliona Sponges Over Years",
       x = "Sample Year",
       y = "Prevalence") +
  theme_minimal()


#This code creates a grouped bar plot showing the prevalence of Cliona sponges at different locations over sample years. Bars are ordered by prevalence within each location and colored by year, with bars for each year displayed side-by-side. Labels and a minimal theme enhance readability, and x-axis labels are rotated for clarity.
ggplot(transectprev, aes(x = reorder(Location, -prev), y = prev, fill = factor(SampleYear))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Prevalence of Cliona Sponges by Location and Year",
       x = "Location",
       y = "Prevalence",
       fill = "Year") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#This code creates a boxplot to display the distribution of Cliona sponge prevalence across different locations. The x-axis shows locations ordered by decreasing prevalence, and the y-axis represents prevalence values, with each boxplot showing the spread and central tendency of prevalence at each location. The fill color is set to light blue, and x-axis labels are rotated for better readability. This plot highlights the variability in sponge prevalence within each location.
ggplot(transectprev, aes(x = reorder(Location, -prev), y = prev)) +
  geom_boxplot(fill = "lightblue") +
  labs(title = "Prevalence of Cliona Sponges by Location (with Distribution)",
       x = "Location",
       y = "Prevalence") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#This code generates a jitter plot to show the prevalence of Cliona sponges at various locations, with each point representing data from a specific year. Locations are ordered by prevalence, and points are colored by SampleYear to distinguish different years. The geom_jitter() function adds slight horizontal spacing to points, reducing overlap and improving clarity. Labels and a minimal theme enhance readability, with x-axis labels rotated for better visualization.
ggplot(transectprev, aes(x = reorder(Location, -prev), y = prev, color = factor(SampleYear))) +
  geom_jitter(width = 0.2, height = 0, size = 2) +
  labs(title = "Prevalence of Cliona Sponges by Location (Jittered Points)",
       x = "Location", 
       y = "Prevalence", 
       color = "Sample Year") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(transectprev, aes( ))







#This code calculates summary statistics for Cliona sponge prevalence by location and year. For each combination of Location and SampleYear, it computes the average prevalence (avg_prevalence), standard deviation (sd_prevalence), and the total number of transects (total_transects). The resulting data frame, summary_stats, provides a summarized view of prevalence patterns over time and across different locations.

summary_stats <- transectprev %>%
  group_by(Location, SampleYear) %>%
  summarise(
    avg_prevalence = mean(prev, na.rm = TRUE),
    sd_prevalence = sd(prev, na.rm = TRUE),
    total_transects = n()
  )


print(summary_stats)


# Assuming your data frame is called 'transectprev'
transectprev$ID <- 1:nrow(transectprev)



# Run Shapiro-Wilk test for normality on `prev` for each Location and SampleYear
shapiro_results <- transectprev %>%
  group_by(Location, SampleYear) %>%
  summarise(shapiro_p_value = shapiro.test(prev)$p.value)

# View the Shapiro-Wilk test results
print(shapiro_results)


# Fit a mixed-effects model with Transect as a random effect
anova_result <- lmer(prev ~ SampleYear * Location + (1 | Transect), data = transectprev)



# Summarize the model results
summary(anova_result)



#packages
library(rstatix)


anova_test(data = transectprev, formula = prev ~ SampleYear*Location + Error(ID/(SampleYear*Location)))


anova_test(data = transectprev, formula = prev ~SampleYear)



#########Running a Two-Way Repeated Measures ANOVA##############

# #First check the normality of the data using Shapiro Wilks test 
# shapiro_results <- cliona %>%
#   group_by(Location, SampleYear) %>%
#   summarise(shapiro_p_value = if (length(unique(presence)) > 1) {
#     shapiro.test(presence)$p.value
#   } else {
#     NA  # Return NA if all values are identical
#   })
# 
# print(shapiro_results)


#I will flag these zero-only groups separately and exclude them from normality testing while still including them in the broader analysis, such as in an ANOVA or other model.

# Run Shapiro-Wilk test for normality on `prev` for each Location and SampleYear
shapiro_results <- transectprev %>%
  group_by(Location, SampleYear) %>%
  summarise(
    all_zero = all(prev == 0),  # Flag groups with all zero values
    shapiro_p_value = if (!all(prev == 0) && length(unique(prev)) > 1) {
      shapiro.test(prev)$p.value
    } else {
      NA  # Skip the test for all-zero groups or groups with identical values
    }
  )

print(shapiro_results)

#hurddle model turn it into a 1 and do a logistic wquation 

#Given this information, the data shows a mix of normal and non-normal distributions across groups. For groups with p-values ≤ 0.05, consider data transformations (e.g., log transformation) or a non-parametric test if normality cannot be achieved, especially if these groups are critical to your analysis.

####So lets try to Log transform the data####

prev_transform <- transectprev %>%
  mutate(prev_log = log(prev + 1))  # Adding 1 to handle zeros if needed


# Run a robust repeated measures ANOVA with WRS2
# Assuming `prev` is your response variable, `SampleYear` and `Location` are factors,
# and `Transect` is the subject (repeated measure).

result <- rmanova(resp = prev, 
                  factA = SampleYear, 
                  factB = Location, 
                  subj = Transect, 
                  data = transectprev)

# Display the results
print(result)

# Install and load lme4 if not already installed
library(lme4)

# Fit a two-way repeated measures model with Transect as a random effect
anova_model <- lmer(prev ~ SampleYear * Location + (1 | Transect), data = transectprev)

# View the summary of the model
summary(anova_model)


library(lmerTest)
summary(anova_model)

anova(anova_model)
plot(anova_model)


levels(transectprev$SampleYear)


#Need to do a generalized linear model on the raw Presence/absence data 



###running a binomial GLMM

# Load lme4 for GLMM
library(lme4)

# Fit the GLMM with binomial family
glmm_model <- glmer(presence ~ SampleYear * Location + (1 | Transect), 
                    family = binomial, data = short_clio)

# View the model summary
summary(glmm_model)

plot(glmm_model)



