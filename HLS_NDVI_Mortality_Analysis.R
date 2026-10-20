###############################################################
# HLS NDVI Mortality Analysis
# Jenny Shepherd
#
# Inputs:
# HLS_monthly_NDVI_final.csv
# master_dataset_cleaned_1.csv
#
# Uses monthly HLS NDVI observations directly
###############################################################

# ============================================================
# 1. PACKAGES
# ============================================================

library(tidyverse)
library(lubridate)
library(lme4)
library(lmerTest)
library(ggeffects)
library(ggplot2)

# ============================================================
# 2. LOAD DATA
# ============================================================
getwd()
setwd("C:/Users/s422478/OneDrive - Cranfield University/Data")
hls <- read_csv(
  "C:/Users/s422478/OneDrive - Cranfield University/Data/Experimental design/New Data/Sentinel2/HLS_monthly_NDVI_final.csv"
)

meta <- read_csv(
  "C:/Users/s422478/OneDrive - Cranfield University/Data/Experimental design/New Data/NDVI_EVI/master_dataset_cleaned_1.csv"
)

# ============================================================
# 3. CREATE CLEAN METADATA TABLE
# ============================================================

meta_clean <- meta %>%
  select(
    plot_id,
    climate,
    forest_type,
    country,
    lat,
    long,
    site_type,
    `%_mort`,
    `mort_year observed`
  ) %>%
  distinct()

meta_clean <- meta_clean %>%
  rename(
    pct_mort = `%_mort`,
    mort_year_observed = `mort_year observed`
  )

# ============================================================
# 4. CLEAN HLS DATA
# ============================================================

hls <- hls %>%
  rename(
    hls_id = site_id,
    ndvi = NDVI,
    years_since = years_from_mortality
  )

# ============================================================
# 4a. LOAD SITE CROSSWALK
# ============================================================

crosswalk <- read_csv(
  "C:/Users/s422478/OneDrive - Cranfield University/Data/Experimental design/Site sampling/site_crosswalk.csv"
)

# ============================================================
# 4b. TRANSLATE HLS SITE NAMES TO PLOT IDs
# ============================================================

hls <- hls %>%
  left_join(
    crosswalk,
    by = c("hls_id" = "site_id")
  )

# check for unmatched sites
unmatched_sites <- hls %>%
  filter(is.na(plot_id)) %>%
  distinct(hls_id)

df <- df %>%
  mutate(
    climate = str_trim(climate),
    climate = str_to_lower(climate)
  )

sort(unique(df$climate))

#standardise upper and lower case
df <- df %>%
  
  mutate(
    climate = str_to_title(climate),
    forest_type = str_to_title(forest_type)
  )

#plot immediate change
#check output

#then you climate summary
climate_change <- plot_changes %>%
  
  group_by(climate) %>%
  
  summarise(
    
    n = n(),
    
    mean_change =
      mean(immediate_change,
           na.rm = TRUE),
    
    sd_change =
      sd(immediate_change,
         na.rm = TRUE),
    
    se =
      sd_change / sqrt(n),
    
    lower =
      mean_change - 1.96 * se,
    
    upper =
      mean_change + 1.96 * se,
    
    .groups = "drop"
  )
#check
print(climate_change)

# test statistical significance
anova_climate <- aov(
  immediate_change ~ climate,
  data = plot_changes
)

summary(anova_climate)

TukeyHSD(anova_climate)

#immediate change by forest type
#summarise by forest type
forest_change <- plot_changes %>%
  
  group_by(forest_type) %>%
  
  summarise(
    
    n = n(),
    
    mean_change =
      mean(
        immediate_change,
        na.rm = TRUE
      ),
    
    sd_change =
      sd(
        immediate_change,
        na.rm = TRUE
      ),
    
    se =
      sd_change / sqrt(n),
    
    lower =
      mean_change - 1.96 * se,
    
    upper =
      mean_change + 1.96 * se,
    
    .groups = "drop"
  )

forest_change

#create a summary table
forest_change %>%
  
  mutate(
    across(
      c(mean_change,
        sd_change,
        lower,
        upper),
      round,
      3
    )
  )

#plot mean change with 95% confidence intervals
library(ggplot2)

ggplot(
  forest_change,
  aes(
    x = forest_type,
    y = mean_change,
    fill = forest_type
  )
) +
  
  geom_col(
    width = 0.7
  ) +
  
  geom_errorbar(
    aes(
      ymin = lower,
      ymax = upper
    ),
    width = 0.2
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  
  labs(
    x = "Forest Type",
    y = "Mean Immediate Change (NDVI)",
    title = "Immediate NDVI Change Following Mortality by Forest Type"
  ) +
  
  theme_minimal(base_size = 8) +
  
  theme(
    legend.position = "none"
  )

# statistical test to formally test whether immediate change differs between forest types
anova_model <- aov(
  immediate_change ~ forest_type,
  data = plot_changes
)

summary(anova_model)
TukeyHSD(anova_model)


# check thresholds

summary(plot_changes$immediate_change)

quantile(plot_changes$immediate_change,
         probs = c(0.05,0.10,0.25,0.5,0.75,0.90,0.95),
         na.rm = TRUE)

sd(plot_changes$immediate_change,
   na.rm = TRUE)

# threshold code

# traffic light colours
colour <- switch(
  
  classification,
  
  "Healthy" = "#2ca25f",
  
  "Watch" = "#fdae61",
  
  "Mortality Alert" = "#d7191c"
)

# risk score from 1 to 100
risk_score <- pmax(
  
  0,
  
  pmin(
    100,
    abs(anomaly / -0.06) * 100
  )
  
)

# leaflet map upgrade (instead of static map)
output$map <- renderLeaflet({
  
  a <- assessment()
  
  colour <- switch(
    a$class,
    "Healthy" = "green",
    "Watch" = "orange",
    "Mortality Alert" = "red"
  )
  
  leaflet() %>%
    addTiles() %>%
    addCircleMarkers(
      lng = input$longitude,
      lat = input$latitude,
      radius = 10,
      color = colour,
      fillOpacity = 0.9,
      popup = paste(
        "Classification:",
        a$class
      )
    )
})



# ============================================================
# 5. JOIN DATASETS
# ============================================================

df <- hls %>%
  left_join(
    meta_clean,
    by = "plot_id"
  )
names(df)
glimpse(df)
# ============================================================
# 6. FILTER TO ANALYSIS WINDOW
# ============================================================

window_df <- df %>%
  filter(
    years_since >= -4,
    years_since <= 5
  )

# ============================================================
# 7. CREATE MONTHS SINCE MORTALITY
# ============================================================

window_df <- window_df %>%
  mutate(
    month = as.numeric(month),
    
    months_since =
      years_since * 12 +
      (month - 1)
  )

# convert month to factor to account for seasonality

window_df$month <- factor(
  window_df$month,
  levels = 1:12,
  labels = c(
    "Jan","Feb","Mar","Apr","May","Jun",
    "Jul","Aug","Sep","Oct","Nov","Dec"
  )
)

# ============================================================
# 8. CHECK DATASET
# ============================================================

cat(
  "Number of sites:",
  n_distinct(window_df$plot_id), "\n"
)

cat(
  "Number of observations:",
  nrow(window_df), "\n"
)

# ======================================
# standardise upper and lower case
# ======================================
window_df <- window_df %>%
  mutate(
    forest_type =
      forest_type %>%
      stringr::str_trim() %>%      # remove extra spaces
      stringr::str_to_lower() %>%  # convert to lower case
      forcats::fct_inorder()
  )

# ============================================================
# 9. OVERALL MONTHLY MODEL
# ============================================================

monthly_model <- lmer(
  ndvi ~
    months_since +
    I(months_since^2) +
    month +
    (1 | plot_id),
  data = window_df
)

summary(monthly_model)

# ============================================================
# 10. FOREST TYPE MODEL
# ============================================================

forest_model <- lmer(
  ndvi ~
    months_since +
    I(months_since^2) +
    month +
    forest_type +
    months_since:forest_type +
    I(months_since^2):forest_type +
    (1 | plot_id),
  data = window_df
)

summary(forest_model)

anova(forest_model)

# ============================================================
# 11. CLIMATE MODEL
# ============================================================

climate_model <- lmer(
  ndvi ~
    months_since +
    I(months_since^2) +
    month +
    climate +
    months_since:climate +
    I(months_since^2):climate +
    (1 | plot_id),
  data = window_df
)

summary(climate_model)

anova(climate_model)

# ============================================================
# 12. MORTALITY SEVERITY
# ============================================================

window_df <- window_df %>%
  mutate(
    mort_group =
      ifelse(
        pct_mort >= 60,
        "High mortality",
        "Low mortality"
      )
  )

mortality_model <- lmer(
  ndvi ~
    months_since +
    I(months_since^2) +
    month +
    mort_group +
    months_since:mort_group +
    I(months_since^2):mort_group +
    (1 | plot_id),
  data = window_df
)

summary(mortality_model)

anova(mortality_model)

# ============================================================
# 13. MONTHLY RECOVERY CURVE
# ============================================================

monthly_summary <- window_df %>%
  group_by(months_since) %>%
  summarise(
    n = n(),
    mean_ndvi = mean(ndvi, na.rm = TRUE),
    sd_ndvi = sd(ndvi, na.rm = TRUE),
    se = sd_ndvi / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    lower = mean_ndvi - 1.96 * se,
    upper = mean_ndvi + 1.96 * se
  )

ggplot(
  monthly_summary,
  aes(months_since, mean_ndvi)
) +
  
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    fill = "lightblue",
    alpha = 0.3
  ) +
  
  geom_line(
    colour = "darkblue",
    linewidth = 0.8
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "red"
  ) +
  
  theme_bw() +
  
  labs(
    x = "Months since mortality",
    y = "Mean NDVI"
  )

# ============================================================
# 14. FOREST TYPE TRAJECTORIES
# ============================================================

forest_summary <- window_df %>%
  group_by(
    forest_type,
    months_since
  ) %>%
  summarise(
    mean_ndvi = mean(ndvi, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(
  forest_summary,
  aes(
    months_since,
    mean_ndvi,
    colour = forest_type
  )
) +
  
  geom_line(
    linewidth = 0.5
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  
  theme_bw(base_size = 9)

# ============================================================
# 15. CLIMATE TRAJECTORIES
# ============================================================

climate_summary <- window_df %>%
  group_by(
    climate,
    months_since
  ) %>%
  summarise(
    mean_ndvi = mean(ndvi, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(
  climate_summary,
  aes(
    months_since,
    mean_ndvi,
    colour = climate
  )
) +
  
  geom_line(
    linewidth = 0.5
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    label = "observed mortality"
  ) +
  
  theme_bw(base_size = 9)

# ============================================================
# 16. IMMEDIATE CHANGE CALCULATION
# ============================================================

plot_changes <- window_df %>%
  group_by(
    plot_id,
    climate,
    forest_type
  ) %>%
  summarise(
    
    pre_mean =
      mean(
        ndvi[
          months_since >= -24 &
            months_since < 0
        ],
        na.rm = TRUE
      ),
    
    post_mean =
      mean(
        ndvi[
          months_since >= 0 &
            months_since <= 24
        ],
        na.rm = TRUE
      ),
    
    immediate_change =
      post_mean - pre_mean,
    
    .groups = "drop"
  )

# ============================================================
# 17. SUMMARY TABLE
# ============================================================

summary_table <- plot_changes %>%
  group_by(
    climate,
    forest_type
  ) %>%
  summarise(
    n_plots = n(),
    mean_change =
      mean(
        immediate_change,
        na.rm = TRUE
      ),
    sd_change =
      sd(
        immediate_change,
        na.rm = TRUE
      ),
    se_change =
      sd_change /
      sqrt(n_plots),
    .groups = "drop"
  )

write_csv(
  summary_table,
  "HLS_Immediate_Change_Summary.csv"
)
print(summary_table)
# ============================================================
# 18. SAVE WORKING DATASET
# ============================================================

write_csv(
  window_df,
  "HLS_window_dataset.csv"
)
print(window_df)

# ============================================
# Visualise immediate change
# ===========================================
#first frequency distribution
ggplot(
  plot_changes,
  aes(x = immediate_change)
) +
  geom_histogram(
    bins = 20,
    fill = "steelblue",
    colour = "black"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "red"
  ) +
  theme_bw() +
  labs(
    x = "Immediate NDVI change",
    y = "Number of plots",
    title = "Distribution of immediate NDVI change"
  )

#climate comparison
ggplot(
  plot_changes,
  aes(
    x = reorder(
      climate,
      immediate_change,
      median,
      na.rm = TRUE
    ),
    y = immediate_change,
    fill = climate
  )
) +
  geom_boxplot(
    alpha = 0.7
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "red"
  ) +
  coord_flip() +
  theme_bw() +
  guides(fill = "none") +
  labs(
    x = "Climate",
    y = "Immediate NDVI change",
    title = "Immediate NDVI change by climate type"
  )
#Forest type comparison
ggplot(
  plot_changes,
  aes(
    x = forest_type,
    y = immediate_change,
    fill = forest_type
  )
) +
  geom_boxplot() +
  geom_hline(
    yintercept = 0,
    colour = "red",
    linetype = "dashed"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  guides(fill = "none") +
  labs(
    x = "Forest type",
    y = "Immediate NDVI change"
  )

#preferred publication graphic mean change +/-95%
#First create summary stats
climate_change <- plot_changes %>%
  group_by(climate) %>%
  summarise(
    n = n(),
    mean_change = mean(
      immediate_change,
      na.rm = TRUE
    ),
    sd_change = sd(
      immediate_change,
      na.rm = TRUE
    ),
    se = sd_change / sqrt(n),
    lower = mean_change - 1.96 * se,
    upper = mean_change + 1.96 * se,
    .groups = "drop"
  )

#then plot
ggplot(
  climate_change,
  aes(
    x = reorder(
      climate,
      mean_change
    ),
    y = mean_change
  )
) +
  
  geom_col(
    fill = "steelblue"
  ) +
  
  geom_errorbar(
    aes(
      ymin = lower,
      ymax = upper
    ),
    width = 0.2
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "red"
  ) +
  
  coord_flip() +
  
  theme_bw() +
  
  labs(
    x = "Climate",
    y = "Mean NDVI change",
    title = "Average immediate NDVI response following mortality"
  )

sort(unique(df$climate))

# ============================================================
# COMPLETE
# ============================================================

cat(
  "\nHLS monthly mortality analysis completed successfully.\n"
)

# =========================================================
# STL decomposition
# ========================================================
install.packages("forecast")
library(forecast)
library(dplyr)

stl_data <- window_df %>%
  arrange(plot_id, date) %>%
  group_by(plot_id) %>%
  group_modify(~{
    
    ts_ndvi <- ts(
      .x$ndvi,
      frequency = 12
    )
    
    fit <- stl(
      ts_ndvi,
      s.window = "periodic"
    )
    
    .x %>%
      mutate(
        trend = as.numeric(fit$time.series[, "trend"]),
        seasonal = as.numeric(fit$time.series[, "seasonal"]),
        residual = as.numeric(fit$time.series[, "remainder"])
      )
    
  })

#LOESS smoothing 
ggplot(
  window_df,
  aes(months_since, ndvi)
) +
  geom_point(
    alpha = 0.1
  ) +
  geom_smooth(
    method = "loess",
    colour = "blue",
    linewidth = 1.2
  ) +
  theme_bw()
#site level smoothing
window_df <- window_df %>%
  group_by(plot_id) %>%
  arrange(date) %>%
  mutate(
    
    ndvi_smooth =
      zoo::rollmean(
        ndvi,
        k = 3,
        fill = NA,
        align = "center"
      )
    
  ) %>%
  ungroup()

#rerun models
monthly_model <- lmer(
  trend ~
    months_since +
    I(months_since^2) +
    (1 | plot_id),
  data = stl_data
)

monthly_model <- lmer(
  ndvi_smooth ~
    months_since +
    I(months_since^2) +
    month +
    (1 | plot_id),
  data = window_df
)
