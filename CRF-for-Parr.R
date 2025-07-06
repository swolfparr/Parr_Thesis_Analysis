# 🧪 Cliona Sponge Presence on USVI Coral Reefs
# Full Analysis Pipeline — Filtered for Infected Coral Species

# ───────────────────────────────
# 📦 Load Required Libraries
# ───────────────────────────────
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(caret)
library(xgboost)
library(Metrics)
library(glmmTMB)
library(lme4)
library(broom.mixed)
library(ggeffects)

# ───────────────────────────────
# 📥 Load and Preprocess Data
# ───────────────────────────────
cliona_raw <- read_csv("Cliona.csv", show_col_types = FALSE)
site_metadata <- read_csv("TCRMP_SiteMetadata.csv")

# Clean and filter relevant observations
cliona <- cliona_raw %>%
  select(
    Recorder, Period, SampleType, Location, SampleYear, Method,
    Transect, SPP, Cliona, Spo1, Spo1ID, Spo2ID, Spo2,
    SampleDate, SampleMonth, BL, Length, Width, Height
  ) %>%
  mutate(
    presence = if_else(
      Cliona > 0 | Spo1ID %in% c("CLSP", "BOSP") | Spo2ID %in% c("CLSP", "BOSP"),
      1, 0
    ),
    presence = replace_na(presence, 0),
    BL = replace_na(BL, 0)
  ) %>%
  filter(
    Period %in% c("Annual", "PostBL"),
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

# ───────────────────────────────
# 🧑‍🔬 Observer-Level Summaries
# ───────────────────────────────
observer_stats <- cliona %>%
  group_by(Recorder) %>%
  summarise(
    total = n(),
    infected = sum(presence, na.rm = TRUE),
    percent_infected = 100 * infected / total
  ) %>%
  arrange(desc(percent_infected))

# Plot percent infected by observer
ggplot(observer_stats, aes(x = reorder(Recorder, -percent_infected), y = percent_infected)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Percent Infected by Observer",
    x = "Observer",
    y = "Percent Infected"
  ) +
  theme_minimal()

# ───────────────────────────────
# 🧽 Filter Out Uninfected Species
# ───────────────────────────────
infected_species <- cliona %>%
  group_by(SPP) %>%
  summarise(ever_infected = any(presence == 1), .groups = "drop") %>%
  filter(ever_infected) %>%
  pull(SPP)

cliona <- cliona %>%
  filter(SPP %in% infected_species)

# ───────────────────────────────
# 🌍 Add Site Metadata and Disturbance Events
# ───────────────────────────────
cliona_with_site_info <- cliona %>%
  left_join(site_metadata, by = "Location") %>%
  mutate(
    bleaching_event = ifelse(SampleYear %in% c(2005, 2010), 1, 0),
    hurricane_event = ifelse(SampleYear == 2017, 1, 0),
    disease_event   = ifelse(SampleYear == 2020, 1, 0),
    lag1_bleaching  = ifelse(SampleYear %in% c(2006, 2011), 1, 0),
    lag1_hurricane  = ifelse(SampleYear == 2018, 1, 0),
    lag1_disease    = ifelse(SampleYear == 2021, 1, 0),
    volume = Length * Width * Height,
    SPP_collapsed = if_else(SPP %in% names(which(table(SPP) >= 50)), SPP, "Other"),
    SPP_collapsed = factor(SPP_collapsed)
  )

# ───────────────────────────────
# 👤 Filter by Observer
# ───────────────────────────────
cliona_tyler <- cliona_with_site_info %>%
  filter(Recorder == "Tyler Burton Smith")

cliona_tyler_rosmin <- cliona_with_site_info %>%
  filter(Recorder %in% c("Tyler Burton Smith", "Rosmin Ennis"))

# ───────────────────────────────
# 📈 Fit GLMM Models
# ───────────────────────────────

# GLMM 1: Depth, volume, year
glmm_model <- glmmTMB(
  presence ~ Depth + volume + SampleYear + (1 | Island/Code),
  data = cliona_with_site_info,
  family = binomial
)

# GLMM 2: Full model with random effects (including observer)
glmm_model2 <- glmmTMB(
  presence ~ Depth + volume + BL + BL * Depth + 
    (1 | Island/Code) + (1 | SampleYear) + (1 | Recorder),
  data = cliona_with_site_info,
  family = binomial
)

# GLMM 3: Tyler-only data
glmm_model_tyler <- glmmTMB(
  presence ~ Depth + volume + BL + BL * Depth + 
    (1 | Island/Code) + (1 | SampleYear),
  data = cliona_tyler,
  family = binomial
)

# GLMM 4: Tyler + Rosmin only
glmm_model_tyler_and_rosmin <- glmmTMB(
  presence ~ Depth + volume + BL + BL * Depth + 
    (1 | Island/Code) + (1 | SampleYear),
  data = cliona_tyler_rosmin,
  family = binomial
)

# ───────────────────────────────
# 📊 Visualize Fixed Effects
# ───────────────────────────────
broom.mixed::tidy(glmm_model2, effects = "fixed", conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  ggplot(aes(x = term, y = estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(
    title = "Fixed Effects Estimates for Cliona Presence (GLMM)",
    x = "Predictor",
    y = "Estimate (log-odds)"
  ) +
  theme_minimal()

# ───────────────────────────────
# 📉 Predict Presence Across Bleaching and Depth
# ───────────────────────────────
pred <- ggpredict(glmm_model_tyler_and_rosmin, terms = c("BL [all]", "Depth [5,15,25]"), bias_correction = TRUE)

ggplot(pred, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) +
  labs(
    x = "Bleaching Percent (BL)",
    y = "Predicted Cliona Presence Probability",
    color = "Depth (m)",
    fill = "Depth (m)"
  ) +
  theme_minimal()

# ───────────────────────────────
# 📆 Cliona Prevalence Over Time
# ───────────────────────────────
cliona_with_site_info %>%
  mutate(SampleYear = as.numeric(as.character(SampleYear))) %>%
  group_by(SampleYear) %>%
  summarise(
    n = n(),
    prevalence = mean(presence, na.rm = TRUE),
    se = sqrt(prevalence * (1 - prevalence) / n),
    ci_lower = prevalence - 1.96 * se,
    ci_upper = prevalence + 1.96 * se
  ) %>%
  ggplot(aes(x = SampleYear, y = prevalence)) +
  geom_line(color = "darkgreen") +
  geom_point(color = "darkgreen") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "green", alpha = 0.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Cliona Presence Prevalence Over Time (Infected Species Only)",
    x = "Year",
    y = "Prevalence (%)"
  ) +
  theme_minimal()

# ───────────────────────────────
# 🗺️ Site-Level Trends Over Time
# ───────────────────────────────
cliona_with_site_info %>%
  mutate(SampleYear = as.numeric(as.character(SampleYear))) %>%
  group_by(SampleYear, Code) %>%
  summarise(
    n = n(),
    prevalence = mean(presence, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = SampleYear, y = prevalence, color = Code)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Cliona Prevalence Over Time by Site",
    x = "Year",
    y = "Prevalence (%)",
    color = "Site Code"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# ───────────────────────────────
# 🐠 Cliona Prevalence by Coral Species
# ───────────────────────────────
cliona_with_site_info %>%
  group_by(SPP) %>%
  summarise(
    n_total = n(),
    n_cliona = sum(presence, na.rm = TRUE),
    prevalence = n_cliona / n_total
  ) %>%
  filter(n_total > 30) %>%
  ggplot(aes(x = reorder(SPP, prevalence), y = prevalence)) +
  geom_col(fill = "darkorange") +
  coord_flip() +
  labs(
    title = "Cliona Prevalence by Coral Species (Only Infected SPP)",
    x = "Coral Species",
    y = "Prevalence (Proportion Infected)"
  ) +
  theme_minimal()
