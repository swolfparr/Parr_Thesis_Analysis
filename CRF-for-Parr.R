# 🧪 Cliona Sponge Presence on USVI Coral Reefs
# Focused Analysis: Bleaching (Binary), Depth & Observer Filtering
# ─────────────────────────────────────────────────────────────────────

# 📦 Load Required Libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(glmmTMB)
library(ggeffects)
library(broom.mixed)

# 📥 Load and Prepare Data
cliona_raw <- read_csv("Cliona.csv", show_col_types = FALSE)
site_metadata <- read_csv("TCRMP_SiteMetadata.csv")

# Filter core variables and create Cliona presence flag
cliona <- cliona_raw %>%
  select(Recorder, Period, SampleType, Location, SampleYear, Method, Transect, SPP,
         Cliona, Spo1ID, Spo2ID, SampleDate, BL, Length, Width, Height) %>%
  mutate(
    presence = if_else(
      Cliona > 0 | Spo1ID %in% c("CLSP", "BOSP") | Spo2ID %in% c("CLSP", "BOSP"),
      1, 0
    ),
    presence = replace_na(presence, 0),
    BL = replace_na(BL, 0)
  ) %>%
  filter(
    Period == "Annual",
    SampleType %in% c("Permanent", "permanent"),
    Method %in% c("intercept", "50 cm belt"),
    !grepl("A", Transect),
    SampleYear >= 2005
  ) %>%
  mutate(
    ID = row_number(),
    SampleYear = as.factor(SampleYear),
    Location = as.factor(Location),
    Transect = as.factor(Transect)
  )

# 🐠 Filter to Coral Species Ever Infected by Cliona
infected_species <- cliona %>%
  group_by(SPP) %>%
  summarise(ever_infected = any(presence == 1), .groups = "drop") %>%
  filter(ever_infected) %>%
  pull(SPP)

cliona <- cliona %>%
  filter(SPP %in% infected_species)

# 🌍 Join Site Metadata & Create Environmental Covariates
cliona <- cliona %>%
  left_join(site_metadata, by = "Location") %>%
  mutate(
    bleaching_event = if_else(SampleYear %in% c(2005, 2010), 1, 0),
    high_bleaching = if_else(BL >= 1, 1, 0),  # Binary bleaching threshold
    volume = Length * Width * Height,
    SPP_collapsed = if_else(SPP %in% names(which(table(SPP) >= 50)), SPP, "Other"),
    SPP_collapsed = factor(SPP_collapsed)
  )


# 📈 Fit GLMM: Depth x High Bleaching (Binary)
glmm <- glmmTMB(
  presence ~ Depth + high_bleaching + Depth:high_bleaching + 
    (1 | Code),
  data = cliona,
  family = binomial
)

summary(glmm)

# 📊 Visualize Model: Predicted Cliona Presence by Depth & Bleaching
preds <- ggpredict(glmm, terms = c("Depth", "high_bleaching"))

ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(
    title = "Predicted Probability of Cliona Presence by Depth",
    subtitle = "Stratified by High Bleaching (Binary, BL ≥ 1)",
    x = "Depth (m)",
    y = "Predicted Cliona Presence",
    color = "High Bleaching",
    fill = "High Bleaching"
  ) +
  scale_color_manual(values = c("#1f77b4", "#d62728"), labels = c("No", "Yes")) +
  scale_fill_manual(values = c("#1f77b4", "#d62728"), labels = c("No", "Yes")) +
  theme_minimal()

# 📌 Interpretation:
# - Bleaching increases probability of Cliona presence at most depths.
# - Negative slope of Depth reflects decreasing presence with depth, though moderated by bleaching.

# 🔎 Optional: Prevalence vs Raw Bleaching Percent (continuous)
cliona_tyler_rosmin %>%
  group_by(BL) %>%
  summarise(prevalence = mean(presence, na.rm = TRUE)) %>%
  ggplot(aes(x = BL, y = prevalence)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = TRUE) +
  labs(
    title = "Cliona Prevalence vs Bleaching Percent",
    x = "Bleaching % (BL)",
    y = "Prevalence of Cliona Presence"
  ) +
  theme_minimal()
