#run the whole script

# =============================================================================
# NDVI / drought mortality analysis
# -----------------------------------------------------------------------------
# This script does two things, in order:
#   PART 1: cleans the categorical columns (fixes inconsistent capitalisation)
#   PART 2: builds the "years since mortality" variable, restricts to a
#           -2 to +5 year window around each site's mortality event, fits the
#           mixed model, and summarises the immediate change and recovery
#           time to the control baseline, by climate and forest type.
#
# Every step is commented. Run the script top to bottom.
# =============================================================================
install.packages(c("dplyr", "lme4", "purrr","ggplot2"))
install.packages("tidyverse")

# ---- 0. Packages -------------------------------------------------------
library(ggplot2)
library(dplyr)
library(readr)
library(purrr)
library(lme4)        # mixed models
library(magrittr)
library(tidyverse)

# ---- 1. Load the data ---------------------------------------------------
#df <- read_csv(C:\Users\s422478\OneDrive - Cranfield University\Data\Experimental design\New Data\NDVI_EVI\"master_dataset_cleaned_1.csv", show_col_types = FALSE)
df <- read_csv(
  "C:/Users/s422478/OneDrive - Cranfield University/Data/Experimental design/New Data/NDVI_EVI/master_dataset_cleaned_1.csv",
  show_col_types = FALSE
)
# The raw file has some awkward column names (a stray character on the first
# column, a % sign, a space in "mort_year observed"). Rather than fight with
# those names, we just rename all 20 columns in the exact order they appear
# in the file. If your file ever has columns in a different order, check
# this list against your file header first.
colnames(df) <- c(
  "plot_id", "date", "ndvi", "trend", "residual", "year", "doy",
  "country", "event_date", "climate", "forest_type", "quality",
  "pct_mort", "site_type", "start_date", "end_date", "in_drought",
  "lat", "long", "mort_year_observed"
)

ncol(df)


# =============================================================================
# PART 1: fix inconsistent capitalisation
# =============================================================================
# Several categorical columns contain the same category spelled with
# different capitalisation (e.g. "Mortality" and "mortality"), which R would
# otherwise treat as two separate categories. This function capitalises just
# the first letter of a string and leaves the rest untouched, which fixes the
# capitalisation without damaging acronyms like "UK" or "SWUSA".
fix_case <- function(x) {
  x <- trimws(x)                                              # remove stray whitespace
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))     # capitalise first letter only
}

cols_to_fix <- c("site_type", "climate", "forest_type", "country")
df <- df %>%
  mutate(across(all_of(cols_to_fix), fix_case))

# Quick check: print the unique categories for each fixed column
for (col in cols_to_fix) {
  cat(col, ":", paste(sort(unique(df[[col]])), collapse = ", "), "\n")
}

# Note: the "quality" column has separator/spelling inconsistencies as well
# (e.g. "Direct_match" vs "Direct match" vs "Partial"), which are NOT fixed
# here, since that needs a decision on which spelling is canonical rather
# than a simple capitalisation fix.


# =============================================================================
# PART 2: years-since-mortality pipeline
# =============================================================================

# ---- 2.1 Before/after mortality dummy -----------------------------------
# post_mortality = 1 if the observation year is at or after the site's
# recorded mortality year, 0 if before, NA if the site has no recorded
# mortality year at all (17 plots in this dataset).
df <- df %>%
  mutate(post_mortality = if_else(
    is.na(`mort_year_observed`), NA_real_,
    as.numeric(year >= `mort_year_observed`)
  ))


# ---- 2.2 Years relative to the mortality event ---------------------------
# Negative values = years before the event, 0 = the event year itself,
# positive values = years after the event.
df <- df %>%
  mutate(years_since = year - `mort_year_observed`)

# ---- 2.3 Restrict to the (-2) -4 to +5 year window -----------------------------
# This window was chosen because it retains good, balanced plot coverage
# (82% of mortality plots have data across the whole window), unlike wider
# windows where only a handful of long-monitored plots remain at the extremes.
#WE WANT TO EXTEND FROM -2 to -4 years_since to make -4 to +5 window
window_df <- df %>%
  filter(!is.na(years_since), years_since >= -4, years_since <= 5)
#save window_df
write_csv(window_df, "C:/Users/s422478/OneDrive - Cranfield University/Data/Experimental design/New Data/NDVI_EVI/window_df.csv")
names(window_df)

# ---- 2.4 Quadratic mixed model: is there a dip then a rise? --------------
# Fixed effects: years_since and years_since^2 (the squared term lets the
# model curve, rather than assuming a straight line, so it can capture a dip
# followed by a rise).
# Random effect: a random intercept per plot_id, since each plot contributes
# many repeated observations over time.
#WE WANT TO USE THE WHOLEDATASET NOT JUST MORTALITY SITES
#mort_window <- window_df %>% filter(site_type == "Mortality")

quad_model <- lmer(
  trend ~ years_since + I(years_since^2) + (1 | plot_id),
  data = window_df
)
print(summary(quad_model))

# The vertex (turning point) of the fitted curve tells you where, in years
# relative to the event, the model estimates the lowest point of trend NDVI.
b <- fixef(quad_model)
vertex <- -b["years_since"] / (2 * b["I(years_since^2)"])
cat("Estimated turning point (years since event):", round(vertex, 2), "\n")

#turning point = 0.24 (previously 0.33 with only -2 years)

#lets find the reference category for forest type so we can calculte vertex by forest type
fixef(forest_model)

#by forest type
b <- fixef(forest_model)
vertex <- -b["years_since"] / (2 * b["I(years_since^2)"])
cat("Estimated turning point (years since event):", round(vertex, 2), "\n")

#by climate type
b <- fixef(climate_model)
vertex <- -b["years_since"] / (2 * b["I(years_since^2)"])
cat("Estimated turning point (years since event):", round(vertex, 2), "\n")

# ---- 2.5 By climate and forest type: immediate change and recovery time --
# For each climate x forest type combination with mortality plots, this
# computes:
#   - immediate_change: mean trend in years 0-1 minus mean trend in years -2,-1
#   - control_baseline: mean trend of control plots in the same climate x
#     forest type combination, over the same window
#   - recovery_year: the first year (0, 1, 2, ...) at which the mortality
#     plots' mean trend reaches or exceeds the control baseline. Requires at
#     least 2 plots contributing at that year to count, and returns NA if the
#     control baseline is never reached within the available data.

combos <- window_df %>%
  filter(site_type == "Mortality") %>%
  distinct(climate, forest_type)

results <- list()

for (i in seq_len(nrow(combos))) {
  clim  <- combos$climate[i]
  ftype <- combos$forest_type[i]

  mort_sub <- window_df %>%
    filter(site_type == "Mortality", climate == clim, forest_type == ftype)
  ctrl_sub <- window_df %>%
    filter(site_type == "Control", climate == clim, forest_type == ftype)

  n_mort_plots <- n_distinct(mort_sub$plot_id)
  n_ctrl_plots <- n_distinct(ctrl_sub$plot_id)

  pre_mean  <- mort_sub %>% filter(years_since %in% c(-2, -1)) %>% pull(trend) %>% mean(na.rm = TRUE)
  post_mean <- mort_sub %>% filter(years_since %in% c(0, 1))   %>% pull(trend) %>% mean(na.rm = TRUE)
  immediate_change <- post_mean - pre_mean

  if (n_ctrl_plots > 0) {
    control_baseline <- mean(ctrl_sub$trend, na.rm = TRUE)

    recovery_year <- NA_real_
    years_available <- sort(unique(mort_sub$years_since[mort_sub$years_since >= 0]))
    for (yr in years_available) {
      yr_data <- mort_sub %>% filter(years_since == yr)
      if (n_distinct(yr_data$plot_id) < 2) next          # skip years with too few plots
      if (mean(yr_data$trend, na.rm = TRUE) >= control_baseline) {
        recovery_year <- yr
        break
      }
    }
  } else {
    control_baseline <- NA_real_
    recovery_year <- NA_real_   # no control plots exist for this combination
  }

  results[[i]] <- tibble(
    climate = clim,
    forest_type = ftype,
    n_mortality_plots = n_mort_plots,
    n_control_plots = n_ctrl_plots,
    immediate_change = round(immediate_change, 4),
    control_baseline = round(control_baseline, 4),
    recovery_year_since_event = recovery_year
  )
}

summary_table <- bind_rows(results) %>% arrange(immediate_change)
print(summary_table, n = Inf)

# Note: combinations with n_control_plots == 0 have no control baseline, so
# recovery_year_since_event will be NA for those rows, meaning "cannot be
# assessed", not "recovered immediately". Combinations with very few plots
# (check n_mortality_plots and n_control_plots) should be treated as
# indicative only, not statistically robust.

write_csv(summary_table, "C:/Users/s422478/OneDrive - Cranfield University/Data/Experimental design/New Data/NDVI_EVI/climate_forest_type_summary.csv")

#if we want COMBINED SITE TYPE analysis not just mortality sites use this code
combos <- window_df %>%
  distinct(climate, forest_type)

results <- list()

for (i in seq_len(nrow(combos))) {
  
  clim  <- combos$climate[i]
  ftype <- combos$forest_type[i]
  
  sub_df <- window_df %>%
    filter(climate == clim,
           forest_type == ftype)
  
  n_plots <- n_distinct(sub_df$plot_id)
  
  pre_mean <- sub_df %>%
    filter(years_since %in% c(-2, -1)) %>%
    summarise(mean_trend = mean(trend, na.rm = TRUE)) %>%
    pull(mean_trend)
  
  post_mean <- sub_df %>%
    filter(years_since %in% c(0, 1)) %>%
    summarise(mean_trend = mean(trend, na.rm = TRUE)) %>%
    pull(mean_trend)
  
  immediate_change <- post_mean - pre_mean
  
  results[[i]] <- tibble(
    climate = clim,
    forest_type = ftype,
    n_plots = n_plots,
    immediate_change = round(immediate_change, 4)
  )
}

summary_table <- bind_rows(results) %>%
  arrange(immediate_change)

print(summary_table, n = Inf)

write_csv(summary_table, "C:/Users/s422478/OneDrive - Cranfield University/Data/Experimental design/New Data/NDVI_EVI/climate_forest_type_summary_combined.csv")

#OK now we are calculating immediate change for each plot
#then we calculate mean SE and 95% CI across plots within 
#each forest and climate type
#1. Calculate change per plot
library(dplyr)

plot_changes <- window_df %>%
  group_by(plot_id, climate, forest_type) %>%
  summarise(
    pre_mean = mean(trend[years_since %in% c(-2, -1)], na.rm = TRUE),
    post_mean = mean(trend[years_since %in% c(0, 1)], na.rm = TRUE),
    immediate_change = post_mean - pre_mean,
    .groups = "drop"
  )
#create a journal ready table with change by forest and climate type
#including SD and SE change
summary_table <- plot_changes %>%
  group_by(climate, forest_type) %>%
  summarise(
    n_plots = n(),
    mean_change = mean(immediate_change, na.rm = TRUE),
    sd_change = sd(immediate_change, na.rm = TRUE),
    se_change = sd_change / sqrt(n_plots),
    ci_lower = mean_change - 1.96 * se_change,
    ci_upper = mean_change + 1.96 * se_change,
    .groups = "drop"
  ) %>%
  arrange(mean_change)

summary_table

#journal table formatting
journal_table <- summary_table %>%
  mutate(
    estimate =
      sprintf("%.3f ± %.3f",
              mean_change,
              se_change)
  ) %>%
  select(
    climate,
    forest_type,
    n_plots,
    estimate
  )

journal_table
write_csv(summary_table, "C:/Users/s422478/OneDrive - Cranfield University/Data/Experimental design/New Data/NDVI_EVI/journal_table.csv")

#create a word table (better than csv)
install.packages("flextable")
library(flextable)
library(officer)

ft <- flextable(journal_table)

ft <- autofit(ft)

save_as_docx(
  "Table 1" = ft,
  path = "C:/Users/s422478/OneDrive - Cranfield University/Data/Experimental design/New Data/NDVI_EVI/Immediate_Change_Table.docx"
)

#a nicer publication-ready Word table
ft <- flextable(journal_table) %>%
  theme_booktabs() %>%
  autofit()

save_as_docx(
  "Table 1. Immediate change by climate and forest type" = ft,
  path = "C:/Users/s422478/OneDrive - Cranfield University/Data/Experimental design/New Data/NDVI_EVI/Table1_ImmediateChange.docx"
)


#visualise the results with a 
#Figure 1: forest plot
library(ggplot2)

summary_table <- summary_table %>%
  mutate(
    group = paste(climate,
                  forest_type,
                  sep = " | ")
  )

ggplot(summary_table,
       aes(x = mean_change,
           y = reorder(group, mean_change))) +
  
  geom_vline(xintercept = 0,
             linetype = "dashed",
             colour = "grey60") +
  
  geom_errorbarh(
    aes(xmin = ci_lower,
        xmax = ci_upper),
    height = 0.2
  ) +
  
  geom_point(size = 3) +
  
  theme_bw() +
  
  labs(
    x = "Immediate change",
    y = NULL
  )
#Figure 2: Heat map
ggplot(summary_table,
       aes(
         x = forest_type,
         y = climate,
         fill = mean_change
       )) +
  
  geom_tile(colour = "white") +
  
  scale_fill_gradient2(
    low = "firebrick",
    mid = "white",
    high = "darkgreen",
    midpoint = 0
  ) +
  
  theme_bw() +
  
  labs(
    fill = "Immediate\nchange"
  )

#Calculate plot changes
plot_changes <- window_df %>%
  group_by(plot_id, climate, forest_type) %>%
  summarise(
    pre_mean = mean(trend[years_since %in% c(-2, -1)], na.rm = TRUE),
    post_mean = mean(trend[years_since %in% c(0, 1)], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    immediate_change = post_mean - pre_mean
  ) %>%
  filter(is.finite(immediate_change))

#Then plot Figure 3:Climate summary
climate_summary <- plot_changes %>%
  group_by(climate) %>%
  summarise(
    n_plots = n(),
    mean_change = mean(immediate_change),
    se = sd(immediate_change)/sqrt(n_plots),
    .groups = "drop"
  )

#Figure 4: Forest-type summary
forest_summary <- plot_changes %>%
  group_by(forest_type) %>%
  summarise(
    n_plots = n(),
    mean_change = mean(immediate_change),
    se = sd(immediate_change)/sqrt(n_plots),
    .groups = "drop"
  )

#a box plot or violin plot shwoing distrbution of plot level
#immediate changes rather than just means might be helpful
ggplot(plot_changes,
       aes(forest_type,
           immediate_change,
           fill = forest_type)) +
  geom_violin() +
  theme_bw()

#improved figure for publication (including raw value distribution)
ggplot(plot_changes,
       aes(climate,
           immediate_change,
           fill = climate)) +
  
  geom_violin(trim = FALSE,
              alpha = 0.7) +
  
  geom_boxplot(width = 0.12,
               fill = "white",
               outlier.shape = NA) +
  
  geom_jitter(width = 0.08,
              alpha = 0.4,
              size = 1) +
  
  geom_hline(yintercept = 0,
             linetype = "solid",
             colour = "grey40") +
  
  theme_bw() +
theme(
  axis.text.x = element_text(
    angle = 45,
    hjust = 1
  ))

#what about pct_mort?
ggplot(plot_changes,
       aes(pct_mort,
           immediate_change)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  theme_bw()

#we need to join pct_mort from window_df 
library(dplyr)

plot_changes <- plot_changes %>%
  left_join(
    window_df %>%
      select(plot_id, pct_mort) %>%
      distinct(),
    by = "plot_id"
  )
#check it has joined OK
head(plot_changes)
#then create scatter plot
ggplot(plot_changes,
       aes(x = pct_mort,
           y = immediate_change)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm",
              se = TRUE,
              colour = "blue") +
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  theme_bw() +
  labs(
    x = "Mortality (%)",
    y = "Immediate change"
  )

#quantify the relationship
model <- lm(immediate_change ~ pct_mort,
            data = plot_changes)

summary(model)

#create scatterplot with fitted regression line
ggplot(plot_changes,
       aes(x = pct_mort,
           y = immediate_change)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm",
              colour = "blue",
              fill = "lightblue") +
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  theme_bw() +
  labs(
    x = "Mortality (%)",
    y = "Immediate change"
  )
#lets check any relationship between mort and climate
#Does mortality severity have different effects in different climates?
lm(immediate_change ~ pct_mort + climate + forest_type,
   data = plot_changes)
lm(immediate_change ~ pct_mort * climate,
   data = plot_changes)
#run the model summary
summary(
  lm(immediate_change ~ pct_mort * climate,
     data = plot_changes)
)

#plot fitted relationship for each forest type
#===NEW CODE===
library(dplyr)
library(ggplot2)

# Calculate statistics
stats <- plot_changes %>%
  group_by(forest_type) %>%
  summarise(
    R2 = summary(lm(immediate_change ~ pct_mort))$r.squared,
    p = summary(lm(immediate_change ~ pct_mort))$coefficients[2,4],
    r = cor(pct_mort, immediate_change, use = "complete.obs"),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      "R² = ", round(R2, 2),
      "\nr = ", round(r, 2),
      "\np = ", signif(p, 2)
    ),
    x = 65,
    y = c(0.065, 0.033, 0.005)  # adjust as needed
  )

# Plot
ggplot(plot_changes, aes(pct_mort, immediate_change, colour = forest_type)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  geom_text(
    data = stats,
    aes(x = x, y = y, label = label, colour = forest_type),
    hjust = 0,
    show.legend = FALSE
  ) +
  theme_bw(base_size = 8)



#plot fitted relationship for each climate
ggplot(plot_changes,
       aes(x = pct_mort,
           y = immediate_change,
           colour = climate)) +
  
  geom_point(alpha = 0.5) +
  
  geom_smooth(
    method = "lm",
    se = FALSE
  ) +
  
  theme_bw() +
  
  labs(
    x = "Mortality (%)",
    y = "Immediate change"
  )
#simpler additive model
model2 <- lm(
  immediate_change ~ climate + forest_type + pct_mort,
  data = plot_changes
)

summary(model2)
#simple model on forest type-immediate change
summary(
  lm(immediate_change ~ forest_type + pct_mort,
     data = plot_changes)
)
#=========CASE STUDIES====================================
#we want to select some case study sites based on high 
#and low mortality, across forest and climate types.
#first we create one row per plot
plot_changes <- plot_changes %>%
  left_join(
    window_df %>%
      select(
        plot_id,
        pct_mort,
        climate,
        forest_type,
        country,
        lat,
        long
      ) %>%
      distinct(),
    by = "plot_id"
  )
#examine mortality distribution to define high and low
summary(plot_changes$pct_mort)

quantile(
  plot_changes$pct_mort,
  probs = c(0.1, 0.25, 0.5, 0.75, 0.9),
  na.rm = TRUE
)
#rebuild the dataframe
plot_changes_clean <- plot_changes %>%
  transmute(
    plot_id,
    climate = climate,
    forest_type = forest_type,
    pct_mort = pct_mort,
    immediate_change
  )
#check its worked
names(plot_changes_clean)
#select high mort forest plots for case study, one per forest type
high_sites <- plot_changes_clean %>%
  group_by(forest_type) %>%
  slice_max(
    pct_mort,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

high_sites

#select low mortality sites, one per forest type
low_sites <- plot_changes_clean %>%
  group_by(forest_type) %>%
  slice_min(
    pct_mort,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

low_sites
#spread the sites across climate too? high mort sites
#including immediate change values
plot_changes_clean %>%
  arrange(desc(pct_mort)) %>%
  select(
    plot_id,
    climate,
    forest_type,
    pct_mort,
    immediate_change
  ) %>%
  head(15)
#and low mort sites
plot_changes_clean %>%
  arrange(pct_mort) %>%
  select(
    plot_id,
    climate,
    forest_type,
    pct_mort,
    immediate_change
  ) %>%
  head(15)
#lets be objective and select three high and three low 
#mort case study sites automatically
#first we define high and low mort groups
high_threshold <- quantile(
  plot_changes_clean$pct_mort,
  0.75,
  na.rm = TRUE
)

low_threshold <- quantile(
  plot_changes_clean$pct_mort,
  0.25,
  na.rm = TRUE
)

high_candidates <- plot_changes_clean %>%
  filter(pct_mort >= high_threshold)

low_candidates <- plot_changes_clean %>%
  filter(pct_mort <= low_threshold)
#view candidate plots for each category - high
high_candidates %>%
  count(climate, forest_type) %>%
  arrange(desc(n))
#14 temp mixed, 5 med conifer, 4 temp decid

#low candidate plots
low_candidates %>%
  count(climate, forest_type) %>%
  arrange(desc(n))
#13 tropical mixed, 11 med conifer, 4 temp deciduous

#select 3 high mort sites automatically
high_case_studies <- high_candidates %>%
  group_by(forest_type) %>%
  slice_max(
    pct_mort,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

high_case_studies
#   plot_id     climate   forest_type    pct_mort   immediate_change
#  1 ymort       Med       Conifer           80          0.0147 
#  2 kat         Temperate Deciduous         56          0.00903
#  3 bav_control Temperate Mixed             87         -0.0394 
low_case_studies <- low_candidates %>%
  group_by(forest_type) %>%
  slice_min(
    pct_mort,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

low_case_studies
# plot_id   climate  forest_type     pct_mort  immediate_change
#  1 c3      Med      Conifer            0         -0.00263
#  2 c1      Med      Deciduous          0         -0.0319 
#  3 bnt4    Tropical Mixed              1         -0.0172 

#combine in a final case study list

#A more rigorous defendible approach is manually selecting
#a case study from each forest typexclimte combination
high_candidates %>%
  select(
    plot_id,
    climate,
    forest_type,
    pct_mort,
    immediate_change
  ) %>%
  arrange(desc(pct_mort)) 
print(high_candidates, n = 35)
#create a shortlist table
candidate_sites <- bind_rows(
  high_candidates %>% mutate(class = "High"),
  low_candidates %>% mutate(class = "Low")
) %>%
  select(
    class,
    plot_id,
    climate,
    forest_type,
    pct_mort,
    immediate_change
  )

View(candidate_sites)

#make a small candidate table
candidate_sites %>%
  arrange(desc(pct_mort)) %>%
  select(
    plot_id,
    climate,
    forest_type,
    pct_mort,
    immediate_change
  )
#we have selected six plots, 3 high mort 3 low mort
#from a selection of forest type/climate combos
case_ids <- c(
  "harz_mort",
  "pb",
  "ymort",
  "c1",
  "bnt4",
  "craig"
)
#simple publication style panel plot of ndvi profiles for these sites
library(dplyr)
library(ggplot2)

case_data <- window_df %>%
  filter(plot_id %in% case_ids)

ggplot(case_data,
       aes(x = years_since,
           y = trend)) +
  
  geom_line(
    colour = "darkgreen",
    linewidth = 0.8
  ) +
  
  geom_point(
    size = 1.5,
    alpha = 0.8
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "red"
  ) +
  
  facet_wrap(
    ~ plot_id,
    ncol = 3
  ) +
  
  theme_bw() +
  
  labs(
    x = "Years since mortality",
    y = "Detrended NDVI"
  )
#calculate and plot annual means instead (currently multiple)
#DOY for each annual point
case_data_annual <- case_data %>%
  group_by(plot_id, years_since) %>%
  summarise(
    mean_ndvi = mean(trend, na.rm = TRUE),
    .groups = "drop"
  )
#then plot
ggplot(case_data_annual,
       aes(years_since,
           mean_ndvi)) +
  
  annotate(
    "rect",
    xmin = 0,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.1,
    fill = "grey70"
  ) +
  
  geom_line(
    colour = "darkgreen",
    linewidth = 1
  ) +
  
  geom_point(size = 2) +
  
  geom_vline(
    xintercept = 0,
    colour = "red",
    linetype = "dashed"
  ) +
  
  scale_x_continuous(
    breaks = seq(-5, 5, by = 1)
  ) +
  
  facet_wrap(
    ~plot_id,
    ncol = 3
  ) +
  
  theme_bw() +
  
  labs(
    x = "Years since mortality",
    y = "Detrended NDVI"
  )


#create panel labels to make it simple to interpret
panel_info <- case_data %>%
  distinct(
    plot_id,
    climate,
    forest_type,
    pct_mort
  ) %>%
  mutate(
    panel_label = paste0(
      plot_id,
      "\n",
      climate, " | ", forest_type,
      "\nMortality = ", round(pct_mort, 1), "%"
    )
  )
#join labels into plotting data
case_data_annual <- case_data_annual %>%
  left_join(
    panel_info %>%
      select(plot_id, panel_label),
    by = "plot_id"
  )
#use labels in the plot
ggplot(case_data_annual,
       aes(years_since,
           mean_ndvi)) +
  
  annotate(
    "rect",
    xmin = 0,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.1,
    fill = "grey70"
  ) +
  
  geom_line(
    colour = "darkgreen",
    linewidth = 1
  ) +
  
  geom_point(size = 2) +
  
  geom_vline(
    xintercept = 0,
    colour = "red",
    linetype = "dashed"
  ) +
  
  scale_x_continuous(
    breaks = seq(-5, 5, by = 1)
  ) +
  
  facet_wrap(
    ~ panel_label,
    ncol = 3
  ) +
  
  theme_bw() +
  
  labs(
    x = "Years since mortality",
    y = "Detrended NDVI"
  ) +
  
  theme(
    strip.text = element_text(size = 8)
  )
#we want all hi mortality sites on top, low mort on btm row
#first we create a table
#create panel order
panel_order <- c(
  "harz_mort",
  "vr",
  "ymort",
  "c1",
  "bnt4",
  "craig"
)
#extract case study data
case_ids <- panel_order

case_data <- window_df %>%
  filter(plot_id %in% case_ids)
#annual means
case_data_annual <- case_data %>%
  group_by(
    plot_id,
    years_since
  ) %>%
  summarise(
    mean_ndvi = mean(trend, na.rm = TRUE),
    .groups = "drop"
  )
#create labels
panel_info <- tibble(
  plot_id = c(
    "harz_mort","vr","ymort",
    "c1","bnt4","craig"
  ),
  panel_label = c(
    "harz_mort | Temperate | Mixed\n80% mortality",
    "vr | Boreal | Conifer\n80% mortality",
    "ymort | Semi-arid | Conifer\n80% mortality",
    "c1 | Med | Deciduous\n0% mortality",
    "bnt4 | Tropical | Mixed\n1% mortality",
    "craig | Temperate | Deciduous\n0% mortality"
  )
)
#join labels
case_data_annual <- case_data_annual %>%
  left_join(panel_info, by = "plot_id")
#order panels
case_data_annual$panel_label <- factor(
  case_data_annual$panel_label,
  levels = panel_info$panel_label
)
#plot
ggplot(case_data_annual,
       aes(years_since,
           mean_ndvi)) +
  
  annotate(
    "rect",
    xmin = 0,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.1,
    fill = "grey70"
  ) +
  
  geom_line(
    colour = "darkgreen",
    linewidth = 0.8
  ) +
  
  geom_point(
    size = 1.5
  ) +
  
  geom_vline(
    xintercept = 0,
    colour = "red",
    linetype = "solid"
  ) +
  
  scale_x_continuous(
    breaks = seq(-5, 5, 1)
  ) +
  
  facet_wrap(
    ~panel_label,
    ncol = 3
  ) +
  
  theme_bw(base_size = 8) +
  
  labs(
    x = "Years pre and post observed mortality",
    y = "Detrended NDVI"
  ) +
  
  theme(
    strip.text = element_text(
      face = "bold",
      size = 6
    )
  )


#=======Case studies #2======================
# Load packages
library(tidyverse)

# Read data
sites <- read.csv("C:/Users/s422478/OneDrive - Cranfield University/Documents/Discussion Papers/Chapter3_Forest Mortality/case study sites_3.csv")

# View data
print(sites)

# Quick summary
summary(sites)

# Number of sites by mortality category
table(sites$mort_cat)

# Mean mortality by category
sites %>%
  group_by(mort_cat) %>%
  summarise(
    n_sites = n(),
    mean_mortality = mean(pct_mort),
    sd_mortality = sd(pct_mort)
  )
#simple mortality plot
library(ggplot2)

ggplot(sites,
       aes(x = reorder(plot_id, pct_mort),
           y = pct_mort,
           fill = mort_cat)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Site",
    y = "Mortality (%)",
    fill = "Mortality category"
  ) +
  theme_bw()
#world map of case study sites
library(ggplot2)
library(maps)

world <- map_data("world")

ggplot() +
  geom_polygon(
    data = world,
    aes(long, lat, group = group),
    fill = "grey90",
    colour = "grey60",
    linewidth = 0.2
  ) +
  geom_point(
    data = sites,
    aes(long, lat, colour = mort_cat),
    size = 2.5
  ) +
  geom_text(
    data = sites,
    aes(long, lat, label = plot_id),
    nudge_y = 4,
    size = 2.5
  ) +
  coord_fixed(1.3) +
  theme_bw() +
  theme(legend.position = "bottom")
  labs(
    colour = "Mortality category",
    x = "Longitude",
    y = "Latitude"
  )
  
#climate and mortality summary
  sites %>%
    group_by(climate) %>%
    summarise(
      n_sites = n(),
      mean_mortality = mean(pct_mort)
    )
  #forest type and mortality summary
  sites %>%
    group_by(forest_type) %>%
    summarise(
      n_sites = n(),
      mean_mortality = mean(pct_mort)
    )
  
  #produce a publication-ready table
  library(dplyr)
  library(knitr)
  
  site_table <- site_table %>%
    mutate(
      `Climate zone` = tools::toTitleCase(`Climate zone`),
      `Forest type` = tools::toTitleCase(`Forest type`),
      `Mortality class` = tools::toTitleCase(`Mortality class`)
    )
  kable(site_table, align = "l")
  library(flextable)
  library(officer)
  
  ft <- flextable(site_table) |>
    autofit()
  
  save_as_docx(
    "Table 1. Case study sites" = ft,
    path = "site_figures/Table1_CaseStudySites.docx"
  )
   
#AERIAL IMAGERY 
#WE will create 250m buffer, add the ground obs plot size, download ESRI world imagery
#export individual maps and combine in a figure
  install.packages(c(
    "sf",
    "terra",
    "tidyterra",
    "maptiles",
    "ggplot2",
    "patchwork",
    "dplyr"
  ))
  install.packages(c("tidyterra"))
  install.packages("patchwork")
  library(patchwork)
  #################################################
  # CASE STUDY SITE AERIAL IMAGERY FIGURE
  # Creates:
  # 1. 250 m buffers around sites
  # 2. Ground obs plot size
  # 3. Aerial imagery maps for each site
  # 4. Individual PNG files
  # 5. Combined publication figure
  #################################################
  #NEW CODE
  library(sf)
  library(terra)
  library(maptiles)
  
  # Read CSV
  sites <- read.csv("C:/Users/s422478/OneDrive - Cranfield University/Documents/Discussion Papers/Chapter3_Forest Mortality/case study sites_2.csv")
  
  # Create output folder
  dir.create("site_figures", showWarnings = FALSE)
  
  for(i in 1:nrow(sites)) {
    
    site_name <- sites$plot_id[i]
    
    cat("\nProcessing:", site_name, "\n")
    
    tryCatch({
      
      # Create point
      
      site <- st_as_sf(
        sites[i, ],
        coords = c("long", "lat"),
        crs = 4326
      )
      
      # Project to metres
      
      site_m <- st_transform(site, 3857)
      
      # Create 500 m buffer
      buffer250 <- st_buffer(site_m, 250)
      
      # Larger extent for imagery
      extent400 <- st_buffer(site_m, 400)
      
      # Convert extent back to WGS84
      extent4326 <- st_transform(extent400, 4326)
      
      #create ground plot footprint (as a square)
      # Plot area (m2)
      plot_area <- sites$plot_ha[i] * 10000
      
      # Equivalent square side length (m)
      plot_side <- sqrt(plot_area)
      
      # Site coordinates in projected CRS
      coords <- st_coordinates(site_m)
      
      x <- coords[1]
      y <- coords[2]
      
      # Square polygon centred on point
      square <- st_polygon(list(rbind(
        c(x - plot_side/2, y - plot_side/2),
        c(x + plot_side/2, y - plot_side/2),
        c(x + plot_side/2, y + plot_side/2),
        c(x - plot_side/2, y + plot_side/2),
        c(x - plot_side/2, y - plot_side/2)
      )))
      
      square <- st_sfc(square, crs = 3857)
      
      square4326 <- st_transform(square, 4326)
      
      square_v <- terra::vect(square4326)
      
      # Back to WGS84
      plot4326 <- st_transform(plot_buffer, 4326)
      
      # Terra vector
      plot_v <- terra::vect(plot4326)
      # Download imagery
      img <- get_tiles(
        extent4326,
        provider = "Esri.WorldImagery",
        crop = TRUE,
        zoom = 17
      )
      
      # IMPORTANT:
      # Convert buffer and point back to WGS84
      buffer4326 <- st_transform(buffer250, 4326)
      point4326 <- st_transform(site_m, 4326)
      
      # Convert to terra vectors
      buffer_v <- terra::vect(buffer4326)
      point_v <- terra::vect(point4326)
      
      # Save PNG
      png(
        filename = paste0(
          "site_figures/",
          site_name,
          ".png"
        ),
        width = 2000,
        height = 2000,
        res = 300
      )
      
      plotRGB(img)
      #add a legend
      legend(
        "bottomleft",
        legend = c(
          "500 m pixel footprint",
          "Ground plot",
          "Site centre"
        ),
        col = c(circle_col, "yellow", "yellow"),
        lty = c(1, 2, NA),
        pch = c(NA, NA, 16),
        lwd = c(3, 2, NA),
        pt.cex = 1.5,
        bty = "n",
        text.col = "white"
      )
      
      #add buffer circle - red for high, blue for low mort
      circle_col <- ifelse(
        sites$mort_cat[i] == "high",
        "red",
        "deepskyblue"
      )
      plot(
        buffer_v,
        add = TRUE,
        border = circle_col,
        lwd = 4
      )
      # Ground plot outline
      plot(
        square_v,
        add = TRUE,
        border = "yellow",
        lwd = 3,
        lty = 2
      )
      #plot yellow point
      plot(
        point_v,
        add = TRUE,
        pch = 16,
        col = "yellow",
        cex = 2
      )
      #improve labels - add plot size and quality flag info
      #label <- sites$label[i]
      label <- paste0(
        sites$label[i],
        "\n",
        sites$plot_ha[i], " ha",
        " (",
        sites$quality_flag[i],
        ")"
      )
      #increase label spacing
      text(
        x = xmin(ext(img)) +
          0.03 * (xmax(ext(img)) - xmin(ext(img))),
        y = ymax(ext(img)) -
          0.04 * (ymax(ext(img)) - ymin(ext(img))),
        labels = label,
        col = "white",
        font = 2,
        cex = 1.3,
        adj = c(0,1)
      )
      
    
      dev.off()
      
      cat("Completed:", site_name, "\n")
      
    }, error = function(e) {
      
      cat(
        "FAILED:",
        site_name,
        "\n",
        e$message,
        "\n"
      )
      
    })
  }
  
  
#now we create a journal-quality multipanel figure with all PNGs
  #install.packages("magick")
  library(magick)
  
  harz  <- image_read("site_figures/harz_mort.png")
  vr    <- image_read("site_figures/vr.png")
  ymort <- image_read("site_figures/ymort.png")
  
  bnt4  <- image_read("site_figures/bnt4.png")
  c1    <- image_read("site_figures/c1.png")
  craig <- image_read("site_figures/craig.png")
  
  # Make all panels identical size
  harz  <- image_resize(harz,  "1200x1200!")
  vr    <- image_resize(vr,    "1200x1200!")
  ymort <- image_resize(ymort, "1200x1200!")
  bnt4  <- image_resize(bnt4,  "1200x1200!")
  c1    <- image_resize(c1,    "1200x1200!")
  craig <- image_resize(craig, "1200x1200!")
  
  top_row <- image_append(c(harz, vr, ymort))
  bottom_row <- image_append(c(c1, bnt4, craig))
  
  final <- image_append(
    c(top_row, bottom_row),
    stack = TRUE
  )
  
  image_write(
    final,
    path = "site_figures/CaseStudySites_Multipanel3.png",
    format = "png"
  )
  
  cat("Multipanel figure saved.\n")
  
#check whether the file was created   
  list.files(pattern = "CaseStudy")  
  list.files(recursive = TRUE, pattern = "CaseStudy")
 #check all six imgaes exist
  list.files("site_figures")
  library(magick)
  library(magick)
  
  files <- list.files(
    "site_figures",
    pattern = "\\.png$",
    full.names = TRUE
  )
  
  files
  #then
  imgs <- lapply(files, image_read)
  
  final <- image_montage(
    do.call(c, imgs),
    tile = "3x2",
    geometry = "800x800+5+5"
  )
  
  image_write(
    final,
    path = "site_figures/CaseStudySites_Multipanel.png",
    format = "png"
  )
  #lets find where plots are saved (working directory)
  getwd()
  #check what files were created
  list.files("site_figures")
  
    #verify locations
  sites[sites$plot_id == "Ymort", ]
  #view the point on a map
  https://www.google.com/maps?q=31.353,35.044

 #plot interactively before saving
  plotRGB(img)

  
#=========================================
#plot the fitted curve from the mixed model

library(ggplot2)

ggplot(mort_window,
       aes(x = years_since,
           y = trend)) +
  geom_point(alpha = 0.1) +
  stat_smooth(
    method = "lm",
    formula = y ~ x + I(x^2),
    colour = "blue",
    linewidth = 1.2
  ) +
  theme_minimal()
#do the same with the combined site type data window_df
library(ggplot2)

ggplot(window_df,
       aes(x = years_since,
           y = trend)) +
  geom_point(alpha = 0.1) +
  stat_smooth(
    method = "lm",
    formula = y ~ x + I(x^2),
    colour = "blue",
    linewidth = 1.2
  ) +
  theme_minimal()

ls()
#this tells us the model name = quad_model
#this results in a straight line so we will try plotting predictions from the fitted mixed model
pred <- data.frame(
  years_since = seq(
    min(mort_window$years_since),
    max(mort_window$years_since),
    length.out = 100
  )
)

pred$trend <- predict(
  quad_model,
  newdata = pred,
  re.form = NA
)

#now the same for the combined dataset
pred <- data.frame(
  years_since = seq(
    min(window_df$years_since),
    max(window_df$years_since),
    length.out = 100
  )
)

pred$trend <- predict(
  quad_model,
  newdata = pred,
  re.form = NA
)
#Quick check
summary(pred$trend)

#these are very small ranges so the recovery is statistically detectable but ecologically quite small. 
#now we plot it
ggplot(pred, aes(x = years_since, y = trend)) +
  geom_line(colour = "blue", linewidth = 1.2) +
  theme_minimal() +
  labs(
    x = "Years pre and post observed mortality",
    y = "Predicted detrended NDVI",
    #title = "Quadratic mixed model fit"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.8
  ) +
  scale_x_continuous(
    breaks = seq(-4, 5, by = 1)
  ) +
  annotate(
    "text",
    x = 0,
    y = Inf,
    angle = 90,
    label = "Observed Mortality",
    size = 2,
    hjust = 1,
    vjust = 1.5,
    colour = "grey30"
  )
# Another useful plot is to overlay the model prediction on the raw data
ggplot(window_df, aes(years_since, trend)) +
  geom_point(alpha = 0.05) +
  geom_line(
    data = pred,
    aes(years_since, trend),
    colour = "red",
    linewidth = 1.5
  ) +
  theme_minimal()
  

#now we have output lets try the simplest possible plot
plot(pred$years_since,
     pred$trend,
     type = "l",
     col = "red",
     lwd = 2)

#another useful plot would be annual values of observed NDVI
#calculate annual means
library(dplyr)

annual_means <- window_df %>%
  group_by(years_since) %>%
  summarise(
    mean_ndvi = mean(ndvi, na.rm = TRUE)
  )
#calculate confidence intervals
library(dplyr)

annual_summary <- window_df %>%
  group_by(years_since) %>%
  summarise(
    n = n(),
    mean_ndvi = mean(trend, na.rm = TRUE),
    sd_ndvi = sd(trend, na.rm = TRUE),
    se = sd_ndvi / sqrt(n),
    lower = mean_ndvi - 1.96 * se,
    upper = mean_ndvi + 1.96 * se
  )

#plot with error/CI shading
ggplot(annual_summary,
       aes(x = years_since,
           y = mean_ndvi)) +
  
  geom_ribbon(
    aes(ymin = lower,
        ymax = upper),
    fill = "lightblue",
    alpha = 0.4
  ) +
  
  geom_line(
    colour = "blue",
    linewidth = 1
  ) +
  
  geom_point(
    colour = "blue",
    size = 3
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "black"
  ) +
  
  annotate(
    "text",
    x = 0,
    y = max(annual_summary$upper),
    label = "Observed Mortality",
    vjust = -0.5,angle = 90,hjust = 1,
    size = 4
  ) +
  
  scale_x_continuous(
    breaks = seq(-4, 5, by = 1)
  ) +
  
  labs(
    x = "Years since mortality",
    y = "Mean annual detrended NDVI"
  ) +
  
  theme_bw(base_size = 8)

#improved journal ready figure (JAE)
ggplot(
  annual_summary,
  aes(x = years_since,
      y = mean_ndvi)
) +
  
  geom_ribbon(
    aes(ymin = lower,
        ymax = upper),
    fill = "#6BAED6",
    alpha = 0.25
  ) +
  
  geom_line(
    colour = "#08519C",
    linewidth = 0.7
  ) +
  
  geom_point(
    colour = "#08519C",
    size = 1.8
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  
  annotate(
    "text",
    x = 0,
    y = max(annual_summary$upper) * 1.01,
    label = "Observed mortality",
    angle = 90,
    vjust = -0.5,
    hjust = 0.95,
    size = 2.5
  ) +
  
  scale_x_continuous(
    breaks = seq(-4, 5, by = 1)
  ) +
  
  labs(
    x = "Years pre and post observed mortality",
    y = "Detrended NDVI (mean ± 95% CI)"
  ) +
  
  coord_cartesian(clip = "off") +
  
  theme_bw(base_size = 8) +
  
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(
      colour = "black",
      linewidth = 0.5
    ),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    plot.margin = margin(
      t = 10,
      r = 5,
      b = 5,
      l = 5
    )
  )
#export at final journal size vs resizing later
ggsave(
  "C:/Users/s422478/OneDrive - Cranfield University/Data/Experimental design/New Data/NDVI_EVI/Figures/Observed_annual_mean_NDVI_all sites_fontsize8.tiff",
  width = 4,
  height = 3.2,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

#Mean annual NDVI declined in the years preceding tree 
#mortality, reaching its lowest value during the 
#mortality year itself, consistent with progressive 
#canopy stress and loss of photosynthetically active 
#vegetation. NDVI increased slightly during the first 
#years following mortality, although substantial 
#interannual variability remained evident.


#As the effect size of the model is very small (albeit statistically significant) it could be
#that some sites show strong recovery while some show none. i.e. the average response is masking important ecological differences. 
#lets do further exploration of forest type, climate and latitude before we conclude that drought effects are weak.
#"The overall response of NDVI to drought was small, but substantial heterogeneity among forest types and climatic regions suggests that resilience is context-dependent."
#Firstly though a map of sites would be useful


ggplot(df, aes(long, lat)) +
  geom_point(alpha = 0.6) +
  coord_equal() +
  theme_minimal()

#colour points by forest type and climate
ggplot(df,
       aes(long, lat, colour = climate)) +
  geom_point(size = 2,
             alpha = 0.7) +
  scale_colour_viridis_d() +
  coord_equal() +
  theme_minimal() +
  labs(
    x = "Longitude",
    y = "Latitude",
    colour = "Climate"
  )
#now by forest type
ggplot(df,
       aes(long, lat, colour = forest_type)) +
  geom_point(size = 2,
             alpha = 0.7) +
  scale_colour_viridis_d() +
  coord_equal() +
  theme_minimal() +
  labs(
    x = "Longitude",
    y = "Latitude",
    colour = "Forest type"
  )
#a better version with country outlines and colour blind friendly palette.
install.packages("maps")
library(ggplot2)
library(maps)

world <- map_data("world")

ggplot() +
  geom_polygon(
    data = world,
    aes(long, lat, group = group),
    fill = "grey95",
    colour = "grey70",
    linewidth = 0.2
  ) +
  geom_point(
    data = df,
    aes(long, lat, colour = forest_type),
    size = 1.5,
    alpha = 0.8
  ) +
  scale_colour_viridis_d() +
  coord_quickmap() +
  theme(legend.position = "bottom") +
  theme_bw()

#new version black outline around points and legend at the bottom
library(ggplot2)
library(maps)
library(viridis)

world <- map_data("world")

ggplot() +
  geom_polygon(
    data = world,
    aes(long, lat, group = group),
    fill = "grey95",
    colour = "grey70",
    linewidth = 0.2
  ) +
  geom_point(
    data = df,
    aes(long, lat, fill = forest_type),
    shape = 21,
    colour = "black",
    stroke = 0.3,
    size = 2,
    alpha = 0.9
  ) +
  scale_fill_viridis_d(
    option = "D",
    name = "Forest type"
  ) +
  coord_quickmap() +
  theme_bw() +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9)
  )

#Now we explore more analysis by climate, forest type and latitude for context specific variance
#Firstly are some forest types more resilient than others?
#we are using the combined dataset (not just mortality sites)
ggplot(window_df,
       aes(years_since, trend, colour = forest_type)) +
  stat_summary(fun = mean, geom = "line") +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.5
  ) +
  annotate(
    "text",
    x = 0,
    y = Inf,
    label = "Mortality event",
    vjust = -0.75,angle = 90,hjust = 1,
    size = 4,
    colour = "grey30",
  ) +
  labs(
    title = NULL,
    x = "Years since mortality",
    y = "Detrended NDVI",
    colour = NULL,
    fill = NULL
  ) +
  scale_x_continuous(
    breaks = seq(-4, 5, by = 1)
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey85"),
    legend.position = "bottom"
  )
#this reveals very different mean NDVI values and variance in the trend pre and post drought.
#Interestingly they all show a detectable dip pre-drought
#Mixed forests show a slight decrease even before the event and this continues for two years post drought then a slow recovery
#Conifer forests show slight decline pre event and post then a recover after year 2 and decrease afgain by year 4
#Deciduous forests are similar with slight decrease around event but steep increase year 2 then steep decrease 
# The deciduous response could be phenological with early leaf senescence and replacement with new young leave with high chlorophyll levels
#"Although the overall response of NDVI to drought was weak, 
#stratification by forest type revealed marked differences 
#in post-disturbance trajectories. 
#Mixed forests maintained consistently high NDVI values and 
#exhibited only minor fluctuations, indicating strong resistance to 
#drought. Deciduous and conifer forests showed greater temporal variability, 
#with evidence of post-disturbance increases in NDVI several years after drought. 
#These results suggest that drought resilience is forest-type dependent 
#and that aggregated analyses may obscure important ecological differences among forest communities".

#Lets run a statistical test - fit a model with an interaction
library(lme4)

forest_model <- lmer(
  trend ~ years_since +
    I(years_since^2) +
    forest_type +
    years_since:forest_type +
    I(years_since^2):forest_type +
    (1 | plot_id),
  data = window_df
)

summary(forest_model)
#lets plot model predictions by forest type with no title 
#same colours for forest types, vertical line at year zero
#and legend at the bottom
library(ggeffects)
library(ggplot2)

pred <- ggpredict(
  forest_model,
  terms = c("years_since", "forest_type")
)

plot(pred) +
  labs(
    title = NULL,
    x = "Years since mortality",
    y = "Detrended NDVI",
    colour = NULL,
    fill = NULL
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x = 0,
    y = Inf,
    label = "Mortality event",
    vjust = -1.0,angle = 90,hjust = 1.0,
    colour = "grey30"
  ) +
  scale_x_continuous(
    breaks = seq(-4, 5, by = 1)
  ) +
  scale_colour_manual(
    values = c(
      "Mixed" = "blue",
      "Deciduous" = "darkgreen",
      "Conifer" = "darkorange"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Mixed" = "blue",
      "Deciduous" = "darkgreen",
      "Conifer" = "darkorange"
    )
  ) +
  theme(
    legend.position = "bottom"
  )

#lets run anova for forest model 
anova(forest_model)
summary(forest_model)
#lets test with type 3 ANOVA more appropriate for mixed models
library(lmerTest)

anova(forest_model, type = 3)
#type not appropriate for lmer model
#lets try and generate p values from lmer
library(lmerTest)
summary(forest_model)
#still not generating p values
#fit model after loading lmer
library(lmerTest)

forest_model <- lmer(
  trend ~ years_since + I(years_since^2) +
    forest_type +
    years_since:forest_type +
    I(years_since^2):forest_type +
    (1 | plot_id),
  data = window_df
)

summary(forest_model)
anova(forest_model, type = 3)

#this shows significant effects of forest type on trend values
#this code will provide pairwise comparisons
install.packages("emmeans")
library(emmeans)

emmeans(forest_model, pairwise ~ forest_type)


# better emmeans version
ggplot(emm_df,
       aes(x = forest_type,
           y = emmean,
           colour = forest_type)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = asymp.LCL,
        ymax = asymp.UCL),
    width = 0.15
  ) +
  labs(
    title = "Distribution Mean NDVI Trend by Forest Type",
    x = "Forest Type",
    y = "Estimated NDVI Trend"
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "none"
  )
  
# secondly ggmeans
library(emmeans)

emm_forest <- emmeans(forest_model, ~ forest_type)

plot(emm_forest)

#add confidence intervals to the means
ggplot(window_df,
       aes(years_since, trend, colour = forest_type)) +
  stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 1.0
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x = 0,
    y = Inf,
    label = "Observed Mortality",
    vjust = -0.5,angle = 90,hjust = 1.25,
    colour = "grey30"
  ) +
  scale_x_continuous(
    breaks = seq(-4, 5, by = 1)
  ) +
  labs(
    #title = "NDVI Trends Relative to Mortality by Forest Type",
    x = "Years pre and post observed mortality",
    y = "Detrended NDVI",
    colour = "Forest Type",
    legend.position = "bottom"
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey85"),
    legend.position = "bottom"
  )
  
#now we do the same for combined dataset, add confidence interval
#including vertical line at year 0

#best version for journal
library(ggplot2)
library(viridis)

ggplot(window_df,
       aes(years_since,
           trend,
           colour = forest_type)) +
  
  stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 0.8
  ) +
  
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.15,
    linewidth = 0.4
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.4
  ) +
  
  annotate(
    "text",
    x = 0,
    y = Inf,
    label = "Observed mortality",
    vjust = 1.3,
    colour = "grey30",
    size = 2.8
  ) +
  scale_x_continuous(
    breaks = seq(-4, 5, by = 1)
  ) +
  scale_colour_viridis_d(
    option = "D",
    end = 0.9,
    name = "Forest type"
  ) +
  
  labs(
    x = "Years pre and post observed mortality",
    y = "Detrended NDVI"
  ) +
  
  theme_bw(base_size = 9) +
  
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    
    panel.grid.minor = element_blank(),
    
    panel.grid.major = element_line(
      colour = "grey90",
      linewidth = 0.25
    ),
    
    panel.border = element_rect(
      colour = "black",
      linewidth = 0.5
    )
  )

#mean + SE ribbons a style frequently seen in JAE
library(dplyr)

summary_df <- window_df %>%
  group_by(forest_type, years_since) %>%
  summarise(
    mean = mean(trend, na.rm = TRUE),
    se = sd(trend, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

ggplot(summary_df,
       aes(years_since,
           mean,
           colour = forest_type,
           fill = forest_type)) +
  
  geom_ribbon(
    aes(
      ymin = mean - se,
      ymax = mean + se
    ),
    alpha = 0.2,
    colour = NA
  ) +
  
  geom_line(linewidth = 0.8) +
  
  geom_vline(
    xintercept = 0,
    linetype = "solid",
    colour = "grey60",
    linewidth = 0.4
  ) +
  annotate(
    "text",
    x = 0,
    y = Inf,
    label = "Observed mortality",
    vjust = -0.5,
    hjust = 1.50,
    angle = 90,
    colour = "grey50",
    size = 2.5
  ) +
  scale_x_continuous(
    breaks = seq(-4, 5, by = 1)
  ) +
  scale_colour_viridis_d(option = "D") +
  scale_fill_viridis_d(option = "D") +
  
  labs(
    x = "Years pre and post mortality",
    y = "Detrended NDVI",
    colour = "Forest type",
    fill = "Forest type"
  ) +
  
  theme_bw(base_size = 9) +
  
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

#the old version
  ggplot(window_df,
         aes(years_since, trend, colour = forest_type)) +
    stat_summary(
      fun = mean,
      geom = "line",
      linewidth = 1
    ) +
    stat_summary(
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.2
    ) +
    labs(
      x = "Years Since Mortality",
      y = "Detrended NDVI",
      colour = "Forest Type"
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "grey40",
      linewidth = 1
    ) +
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = "Observed Mortality",
      vjust = 1.5,
      colour = "grey30"
    ) +
    scale_colour_manual(
      values = c(
        "Conifer" = "darkorange",
        "Deciduous" = "darkgreen",
        "Mixed" = "blue"
      )
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      
      # White backgrounds
      panel.background = element_rect(fill = "white", colour = "white"),
      plot.background = element_rect(fill = "white", colour = "white"),
      
      # Grey gridlines
      panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
      panel.grid.minor = element_line(colour = "grey92", linewidth = 0.2),
      
      # Keep gridlines visible with theme_classic()
      panel.grid = element_blank()
    )

#Check the sample size for each forest type
table(window_df$forest_type)
#there are even numbers for each forest type (slightly fewer data points for deciduous)
#now we test whether differences are statistically significant as a whole not just individual 
anova(quad_model, forest_model)
#calculate p values to determine statistical significance
anova_out <- anova(forest_model)
anova_out

#comparison with the overall model with the forest type model (Anova) 
#reveals a highly significant difference in post drought trajectory by forest type. 
#Mixed forests have highest NDVI but small changes and minimal drought impact
#This can be interpreted as high resistance
#Deciduous forests have larger temporal fluctuations (phenological response to drought?)
#And strongest evidence of recovery after disturbance - could be canopy replacement with new young growth
#Lower resistance but stronger recovery
#Conifer forests show a response intermediate to the other forest types
#moderate resistance and moderate recovery - due to evergreen phenology?

#resilience appears to be ecosystem specific and different forest types have different trajectories
#this aligns well with Holling's ecological resilience model and Nikinmaa's nested resilience categories
#the general/pooled model reflects the engineering resilience - how quickly does NDVI recover
#the forest type analysis shifts this to ecological resilience - how do different systems respond to diistrubance?

#The narrative becomes:
#1.Overall drought effect is modest.
#2.Forest type significantly improves prediction of NDVI trajectories.
#3.Mixed forests appear most resistant.
#4.Deciduous forests show the strongest recovery response.
#5.Therefore resilience is context-dependent rather than uniform across forests.

#now lets explore variance by climate
ggplot(mort_window,
       aes(years_since, trend, colour = climate)) +
  stat_summary(fun = mean, geom = "line") +
  theme_minimal()

#run statistical test
library(lme4)

climate_model <- lmer(
  trend ~ years_since +
    I(years_since^2) +
    climate +
    years_since:climate +
    I(years_since^2):climate +
    (1 | plot_id),
  data = mort_window
)

summary(climate_model)

#add confidence intervals to the means
ggplot(mort_window,
       aes(years_since, trend, colour = climate)) +
  stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 1
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
  ) +
  labs(
    title = "NDVI Trends Relative to Mortality by Climate Type",
    x = "Years Since Mortality",
    y = "NDVI Trend",
    colour = "Climate Type",
    legend.position = "bottom"
  ) +
  theme(
    plot.title = element_text(size = 10),
    legend.position = "bottom"
  )
theme_minimal()

#check sample size for each climate type
table(mort_window$climate)
#some climate types have lower sample size (Boreal and Continental) but mostly a good spread

#now lets run ANOVA to compare statistical significant differences between all three models
anova(quad_model, forest_model, climate_model)

#plot the fitted trajectories
install.packages("ggeffects")
library(ggeffects)

pred_climate <- ggpredict(
  climate_model,
  terms = c("years_since [all]", "climate")
)

plot(pred_climate)

#lets make tropical climate clearer in the plot (using ggplot) by making it lightgreen
library(ggeffects)
library(ggplot2)

pred_forest_type <- ggpredict(
  climate_model,
  terms = c("years_since [all]", "climate")
)

ggplot(pred_climate,
       aes(x = x, y = predicted,
           colour = group, fill = group)) +
  geom_line(linewidth = 1.0) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, colour = NA) +
  scale_colour_manual(values = c(
    "Continental" = "blue",
    "Med" = "darkgreen",
    "Semi-arid" = "orange",
    "Temperate" = "red",
    "Tropical" = "lightgreen"
  )) +
  scale_fill_manual(values = c(
    "Continental" = "blue",
    "Med" = "darkgreen",
    "Semi-arid" = "orange",
    "Temperate" = "red",
    "Tropical" = "lightgreen"
  )) +
  labs(
    x = "Years since mortality",
    y = "NDVI trend",
    colour = "Climate",
    fill = "Climate",
    title = "Predicted NDVI Recovery"
  ) +
  theme_classic()

#=======================================================
#DIVIDE DATA BY MORTALITY SEVERITY AND REANALYSE
#we want to try rerunning the analysis using two groups of mortality
#this is important because the very wide range of mortality could be
#skewing or smoothing the results? USING COMBINED DATASET
#group 1 = 60% +
#group 2 = 0-59%
window_df$mort_group <- ifelse(
  window_df$pct_mort >= 60,
  "High mortality (≥60%)",
  "Low mortality (<60%)"
)  
#check the number of observations in each group:
table(window_df$mort_group)
#more importantly number of unique sites
table(
  unique(window_df[, c("plot_id", "mort_group")])$mort_group
)

library(dplyr)

window_df %>%
  distinct(plot_id, mort_group) %>%
  count(mort_group)

str(window_df)
window_df %>%
  summarise(
    rows = n(),
    plots = n_distinct(plot_id)
  )

#fit mortality_model60 (to distinguish it from the original model
#based on 40% mortality threshold)

mortality_model60 <- lmer(
  trend ~ years_since +
    I(years_since^2) +
    mort_group +
    years_since:mort_group +
    I(years_since^2):mort_group +
    (1 | plot_id),
  data = window_df,
  REML = FALSE
)

#anova(quad_model, mortality_model60)
#the anova results reveal highly significant results between mortality severity groups
#the stat summary and coefficients will provide more detail

summary(mortality_model60)
#run anova type3 test too
library(lmerTest)

anova(mortality_model60, type = 3)

#summary results show high baseline ndvi trend in mort_groupLow mortality (<40%)
#Estimate = 0.1147
#t = 3.86
#they also have a different temporal rate of change i.e. 
#the trajectory is different not just the level
#I(years_since^2):mort_groupLow mortality (<40>)
#Estimate = -0.0009475
#t = -8.45
#table(analysis_data$mort_group)


#plot predicted trajectories to show whether groups have different trajectories
library(ggeffects)

pred_mort <- ggpredict(
  mortality_model60,
  terms = c("years_since [all]", "mort_group")
)

plot(pred_mort) +
  geom_vline(
    xintercept = 0,
    linetype = "solid",
    colour = "grey60",
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = 0,
    y = Inf,
    label = "Mortality event",
    vjust = -0.5,
    angle = 90,
    hjust = 1.25,
    colour = "grey50"
  ) +
  labs(
    #title = "Predicted NDVI Recovery by Forest Type and Mortality Severity",
    x = "Years pre and post observed mortality",
    y = "Detrended NDVI",
    colour = "Mortality Group",
    fill = "Mortality Group"
  ) +
  scale_x_continuous(
    breaks = seq(-4, 5, by = 1)
  ) +
  theme_bw(base_size = 8) +
  theme(
    axis.title = element_text(size = 8),
    legend.position = "bottom"
  )

#pairwise comparison to test whether the mean NDVI trends differ by mort group
library(emmeans)

emmeans(mortality_model, pairwise ~ mort60)

#check if the 60% threshold produces balanced groups 
window_df%>%
  distinct(plot_id, mort60, pct_mort) %>%
  group_by(mort60) %>%
  summarise(
    n_plots = n(),
    mean_mort = mean(pct_mort, na.rm = TRUE),
    sd_mort = sd(pct_mort, na.rm = TRUE)
  )

#compare mortality model with pooled model and forest type model
AIC(
  quad_model,
  forest_model,
  mortality_model
)
#this show that forest type is the stronger predictor although 
#mortality severity is secondary predictor of NDVI variance
#lets fit a model with both predictors:
combined_model <- lmer(
  trend ~ years_since +
    I(years_since^2) +
    forest_type +
    mort_group +
    years_since:forest_type +
    years_since:mort_group +
    I(years_since^2):forest_type +
    I(years_since^2):mort_group +
    (1 | plot_id),
  data = window_df
)
#and compare AIC
AIC(
  quad_model,
  mortality_model,
  forest_model,
  combined_model
)
#the resulst show that a model including both forest type 
#and mortality severity has the strongest model fit (AIC)
#this indicates that post drought NDVI recovery trajectories 
#are jointly influenced by forest composition and severity of mortality

#lets plot the combined model predictions
#first check the model predictions exist
library(ggeffects)

pred_combined <- ggpredict(
  combined_model,
  terms = c("years_since [all]", "mort_group", "forest_type")
)

head(pred_combined)

#then plot separate panels by forest type
library(ggplot2)

ggplot(pred_combined,
       aes(x = x,
           y = predicted,
           colour = group,
           fill = group)) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(
    aes(ymin = conf.low,
        ymax = conf.high),
    alpha = 0.2,
    colour = NA
  ) +
  facet_wrap(~facet) +
  labs(
    #title = "Predicted NDVI Recovery by Forest Type and Mortality Severity",
    x = "Years pre and post observed mortality",
    y = "Detrended NDVI",
    colour = "Mortality Group",
    fill = "Mortality Group"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "solid",
    colour = "grey60",
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x = 0,
    y = Inf,
    label = "Observed Mortality",
    vjust = -0.5,
    #angle = 90,
    hjust = 1,
    size = 2,
    colour = "grey30"
  ) +
  scale_x_continuous(
    breaks = seq(-4, 5, by = 1)
  ) +
  theme(
    legend.position = "bottom"
  ) +
  theme_classic(base_size = 8) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(
      fill = "white",
      colour = "black"
    ),
    strip.text = element_text(face = "bold")
  )
  


#Fit a model with three way interaction
full_model <- lmer(
  trend ~ years_since +
    I(years_since^2) +
    forest_type * mort_group +
    years_since * forest_type * mort_group +
    I(years_since^2) * forest_type * mort_group +
    (1 | plot_id),
  data = window_df
)
unique(pred_combined$facet)
unique(pred_combined$group)
#quick plot
ggplot(pred_combined,
       aes(x = x,
           y = predicted,
           colour = group)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~facet)
#another one
ggplot(pred_combined,
       aes(x = x,
           y = predicted,
           colour = group,
           fill = group)) +
  geom_line(linewidth = 1.0) +
  geom_ribbon(
    aes(ymin = conf.low,
        ymax = conf.high),
    alpha = 0.20,
    colour = NA
  ) +
  facet_wrap(~facet, nrow = 1) +
  labs(
    #title = "Predicted NDVI Recovery by Forest Type and Mortality Severity",
    x = "Years pre and post observed mortality",
    y = "Detrended NDVI",
    colour = "Mortality Group",
    fill = "Mortality Group"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "solid",
    colour = "grey40",
    linewidth = 0.6
  ) +
  
  annotate(
    "text",
    x = 0,
    y = Inf,
    label = "Observed Mortality",
    angle = 90,
    vjust = -0.5,
    hjust = 1.25,
    size = 4
  ) +
  
  scale_x_continuous(
    breaks = seq(-4, 5, by = 1),
    limits = c(-4, 5)
  ) +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "bottom"
  )

#rerun the chart with three mortality thresholds to test the 
#different response by forest type and how this might help 
#or hinder the detection of forest mortality
#create new mortality classes
window_df <- window_df %>%
  mutate(
    mort_group = case_when(
      `%_mort` >= 75 ~ "≥75%",
      `%_mort` >= 60 ~ "60-74%",
      `%_mort` >= 40 ~ "40-59%",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(mort_group))
#make the order explicit
data$mort_group <- factor(
  data$mort_group,
  levels = c("40-59%", "60-74%", "≥75%")
)
head(window_df)

#run this before fitting combined model
window_df <- window_df %>%
  mutate(
    mort_group = case_when(
      `pct_mort` >= 75 ~ "≥75%",
      `pct_mort` >= 60 & `pct_mort` < 75 ~ "60-74%",
      `pct_mort` >= 40 & `pct_mort` < 60 ~ "40-59%",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(mort_group))


#rerun entire modelling section
combined_model <- lmer(
  trend ~ years_since * mort_group * forest_type +
    (1 | plot_id),
  data = window_df
)
#then
pred_combined <- ggpredict(
  combined_model,
  terms = c(
    "years_since [all]",
    "mort_group",
    "forest_type"
  )
)
#now check the mortality groups
unique(pred_combined$group)


#inspect model dataset
table(data$mort_group)

#order the mortality groups logically
pred_combined$group <- factor(
  pred_combined$group,
  levels = c("40-59%", "60-74%", "≥75%")
)

#plot by forest type and mortality class
ggplot(
  pred_combined,
  aes(
    x = x,
    y = predicted,
    colour = group,
    fill = group
  )
) +
  
  geom_ribbon(
    aes(
      ymin = conf.low,
      ymax = conf.high
    ),
    alpha = 0.20,
    colour = NA
  ) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_vline(
    xintercept = 0,
    colour = "grey",
    linetype = "solid",
    linewidth = 0.6
  ) +
  
  facet_grid(
    group ~ facet
  ) +
  
  scale_x_continuous(
    breaks = -4:5,
    limits = c(-4, 5)
  ) +
  
  labs(
    x = "Years pre and post observed mortality",
    y = "Detrended NDVI"
  ) +
  
  guides(
    colour = "none",
    fill = "none"
  ) +
  
  theme_classic(base_size = 8)  +
  
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold")
  )

#we may need to amend the chart if there are insufficient
#sample sizes in some groups - lets check 
table(window_df$forest_type, window_df$mort_group)
#the output confirms this:
#           ≥75% 40-59% 60-74%
#Conifer   1147    572    876
#Deciduous    0   2069      0
#Mixed      644    876    460

#we amend the code to leave panels blank where there is no data
#with a simple label 'no observations'

#first we need to remove any modelled data in these rows
pred_combined <- pred_combined |>
  dplyr::filter(
    !(facet == "Deciduous" & group %in% c("60-74%", "≥75%"))
  )
#then before plotting set the mortality categories in logical order
pred_combined$group <- factor(
  pred_combined$group,
  levels = c("≥75%", "60-74%", "40-59%")
)
no_data$mortality_class <- factor(
  no_data$group,
  levels = c("≥75%", "60-74%", "40-59%")
)


# Then plot the panels with no observations
# Order rows so that:
# top = 40-59%
# middle = 60-74%
# bottom = ≥75%

# Only annotate the two empty Deciduous panels

# Convert to plain dataframe
pred_df <- as.data.frame(pred_combined)
# Order rows: top=40-59%, middle=60-74%, bottom=≥75%
pred_df$group <- factor(
  pred_df$group,
  levels = c("40-59%", "60-74%", "≥75%")
)
# Labels for empty panels only
no_data <- data.frame(
  group = factor(
    c("60-74%", "≥75%"),
    levels = c("40-59%", "60-74%", "≥75%")
  ),
  facet = c("Deciduous", "Deciduous"),
  x = c(0.5, 0.5),
  y = c(0.60, 0.60),
  label = "No observations\n(n = 0)",
          "No observations\n(n = 0)")

ggplot(
  pred_df,
  aes(
    x = x,
    y = predicted,
    colour = group,
    fill = group
  )
) +
  
  geom_ribbon(
    aes(
      ymin = conf.low,
      ymax = conf.high
    ),
    alpha = 0.15,
    colour = NA
  ) +
  
  geom_line(linewidth = 1) +
  
  geom_vline(
    xintercept = 0,
    colour = "grey",
    linewidth = 0.6
  ) +
  
  geom_text(
    data = no_data,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    colour = "grey40",
    size = 3
  ) +
  
  facet_grid(group ~ facet) +
  
  scale_x_continuous(
    breaks = -4:5,
    limits = c(-4, 5)
  ) +
  
  scale_colour_manual(
    values = c(
      "40-59%" = "#619CFF",
      "60-74%" = "#00BA38",
      "≥75%"   = "#F8766D"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "40-59%" = "#619CFF",
      "60-74%" = "#00BA38",
      "≥75%"   = "#F8766D"
    )
  ) +
  
  labs(
    x = "Years pre and post observed mortality",
    y = "Detrended NDVI"
  ) +
  
  theme_classic(base_size = 8) +
  
  theme(
    legend.position = "none",
    strip.background = element_rect(
      fill = "white",
      colour = "black"
    ),
    strip.text = element_text(face = "bold")
  )


#to create a graph with actual data not modelled/predicted
ggplot(window_df,
       aes(years_since, trend, colour = mort_group)) +
  
  stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 0.8
  ) +
  
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.15,
    linewidth = 0.4
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "solid",
    colour = "grey",
    linewidth = 0.4
  ) +
  
  annotate(
    "text",
    x = 0,
    y = Inf,
    label = "Observed mortality",
    #angle = 90,
    vjust = 1,
    hjust = 0.5,
    size = 2.8
  ) +
  
  scale_x_continuous(
    breaks = seq(-4, 5, by = 1),
    limits = c(-4, 5)
  ) +
  scale_colour_manual(
    values = c(
      "High mortality (≥60%)" = "#D55E00",
      "Low mortality (<60%)" = "#0072B2"
    )
  )+
  labs(
    x = "Years pre and post observed mortality",
    y = "Detrended NDVI",
    colour = "Mortality group"
  ) +
  
  theme_bw(base_size = 9) +
  
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

#create a JAE ready version with ribbons instead of error bars
#create mean and SE values summary then plot from that
library(dplyr)
library(ggplot2)

# Calculate means and SEs
mort_summary <- window_df %>%
  group_by(mort_group, years_since) %>%
  summarise(
    mean_ndvi = mean(trend, na.rm = TRUE),
    se = sd(trend, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    lower = mean_ndvi - se,
    upper = mean_ndvi + se
  )

# Plot
ggplot(
  mort_summary,
  aes(
    x = years_since,
    y = mean_ndvi,
    colour = mort_group,
    fill = mort_group
  )
) +
  
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.20,
    colour = NA
  ) +
  
  geom_line(
    linewidth = 0.9
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "black",
    linewidth = 0.4
  ) +
  
  annotate(
    "text",
    x = 0,
    y = Inf,
    label = "Observed mortality",
    vjust = 1.25,
    hjust = 0.5,
    size = 2.8
  ) +
  
  scale_x_continuous(
    breaks = seq(-4, 5, by = 1),
    limits = c(-4, 5)
  ) +
  
  scale_colour_manual(
    values = c(
      "High mortality (≥60%)" = "#D55E00",
      "Low mortality (<60%)" = "#0072B2"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "High mortality (≥60%)" = "#D55E00",
      "Low mortality (<60%)" = "#0072B2"
    )
  ) +
  
  labs(
    x = "Years pre and post observed mortality",
    y = "Detrended NDVI",
    colour = "Mortality group",
    fill = "Mortality group"
  ) +
  
  theme_bw(base_size = 8) +
  
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )
#try plotting without the manual scales - this reveals three categories
ggplot(
  mort_summary,
  aes(
    x = years_since,
    y = mean_ndvi,
    colour = mort_group,
    fill = mort_group
  )
) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.20,
    colour = NA
  ) +
  
  geom_line(
    linewidth = 0.9
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "solid",
    colour = "grey60",
    linewidth = 0.4
  ) +
  
  annotate(
    "text",
    x = 0,
    y = Inf,
    label = "Observed mortality",
    vjust = 1.25,
    angle= 90,
    hjust = 1.05,
    size = 3
  ) +
  labs(
    x = "Years pre and post observed mortality",
    y = "Detrended NDVI",
    colour = "Mortality group",
    fill = "Mortality group"
  ) +
  scale_x_continuous(
    breaks = seq(-4, 5, by = 1),
    limits = c(-4, 5)) +
  
  #theme_bw(base_size = 8) +
  #theme(
   # legend.position = "inside",
    #legend.position.inside = c(0.5, 0.05),
    #legend.justification = c(0.5, 0)
 # )
  theme_bw(base_size = 12) +
  theme(
    legend.position = c(0.15, 0.05),
    legend.justification = c(0, 0),
    legend.direction = "vertical",
    legend.background = element_rect(
      fill = alpha("white", 0.85),
      colour = "grey80"
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(colour = "black")
  )
 
  ggsave(
    "C:/Users/s422478/OneDrive - Cranfield University/Data/Experimental design/New Data/NDVI_EVI/Figures/3 cats mortality_ndvi.tiff",
    width = 8,
    height = 6.5,
    dpi = 600,
    compression = "lzw"
  )
  

#lets revert to the two mortality categories and run 
#by forest type in a panel plot
ggplot(
  mort_summary,
  aes(
    x = years_since,
    y = mean_ndvi,
    colour = mort_group,
    fill = mort_group
  )
) +
  
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.20,
    colour = NA
  ) +
  
  geom_line(linewidth = 0.9) +
  
  geom_vline(
    xintercept = 0,
    linetype = "solid",
    colour = "grey60",
    linewidth = 0.4
  ) +
  
  facet_wrap(~ forest_type) +
  
  scale_colour_manual(
    values = c(
      "High mortality (≥60%)" = "#D55E00",
      "Low mortality (<60%)" = "#0072B2"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "High mortality (≥60%)" = "#D55E00",
      "Low mortality (<60%)" = "#0072B2"
    )
  ) +
  
  labs(
    x = "Years pre and post observed mortality",
    y = "Detrended NDVI",
    colour = "Mortality group",
    fill = "Mortality group"
  ) +
  
  theme_bw(base_size = 12) +
  
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  "C:/Users/s422478/OneDrive - Cranfield University/Data/Experimental design/New Data/NDVI_EVI/Figures/2 cats mortality_ndvi.tiff",
  width = 8,
  height = 6.5,
  dpi = 600,
  compression = "lzw"
)

#a nice way to show the difference between the groups- this is modelled data
emmip(
  mortality_model60,
  mort_group ~ years_since,
  CIs = TRUE
)
#empirical distribution plot
library(dplyr)
library(ggplot2)

mort_summary <- window_df %>%
  group_by(mort_group) %>%
  summarise(
    mean_ndvi = mean(trend, na.rm = TRUE),
    se = sd(trend, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

ggplot(mort_summary,
       aes(x = mort_group,
           y = mean_ndvi,
           colour = mort_group)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = mean_ndvi - 1.96 * se,
      ymax = mean_ndvi + 1.96 * se
    ),
    width = 0.15
  ) +
  labs(
    #title = "Observed Mean NDVI Trend by Mortality Severity",
    x = "Mortality Group",
    y = "Mean NDVI Trend",
    colour = "Mortality Group"
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom"
  )

#plot by forest type and mortality class
ggplot(pred_combined,
       aes(x = x,
           y = predicted,
           colour = group,
           fill = group)) +
  
  geom_line(linewidth = 1) +
  
  geom_ribbon(
    aes(ymin = conf.low,
        ymax = conf.high),
    alpha = 0.20,
    colour = NA
  ) +
  
  facet_grid(facet ~ group) +
  
  geom_vline(
    xintercept = 0,
    linetype = "solid",
    colour = "grey60"
  ) +
  scale_x_continuous(
    breaks = seq(-4, 5, by = 1),
    limits = c(-4, 5)) +
  labs(
    x = "Years pre and post observed mortality",
    y = "Detrended NDVI",
    colour = "Mortality Class",
    fill = "Mortality Class"
  ) +
  
  theme_classic(base_size = 9)

#statistical test for empirical data
t.test(
  trend ~ mort_group,
  data = mort_window
)
#look at variation amongst sites
window_df %>%
  group_by(plot_id, mort_group) %>%
  summarise(
    mean_trend = mean(trend),
    .groups = "drop"
  ) %>%
  ggplot(aes(mort_group, mean_trend, colour = mort_group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  labs(
    x = "Mortality Group",
    y = "Mean NDVI Trend"
  ) +
  theme_classic()
#tidy up legend and x axis labels
window_df %>%
  group_by(plot_id, mort_group) %>%
  summarise(
    mean_trend = mean(trend),
    .groups = "drop"
  ) %>%
  ggplot(aes(mort_group,
             mean_trend,
             colour = mort_group)) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.7,
    size = 1
  ) +
  labs(
    #title = "Mean NDVI Trend by Mortality Severity",
    x = "Mortality Group",
    y = "Mean NDVI Trend",
    colour = "Mortality Group"
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )
#do similar distribution plot by forest type
window_df %>%
  group_by(plot_id, forest_type) %>%
  summarise(
    mean_trend = mean(trend, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ggplot(aes(
    x = forest_type,
    y = mean_trend,
    colour = forest_type
  )) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.7,
    size = 1
  ) +
  scale_colour_manual(
    values = c(
      "Mixed" = "blue",
      "Deciduous" = "darkgreen",
      "Conifer" = "red"
    )
  ) +
  labs(
    title = "Mean NDVI Trend by Forest Type",
    x = "Forest Type",
    y = "Mean NDVI Trend",
    colour = "Forest Type"
  ) +
  theme_minimal()

#replicate this code for mortality severity
window_df %>%
  group_by(plot_id, mort_group) %>%
  summarise(
    mean_trend = mean(trend, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ggplot(
    aes(
      x = mort_group,
      y = mean_trend,
      fill = mort_group
    )
  ) +
  geom_boxplot(outlier.shape = NA)

#Lets plot some individual sites
#First, calculate years relative to mortality
#If you cant remember the name of your dataframe use
ls()
#that is everything in the directory which doesnt really help!
#instead use:
formula(forest_model)
#this provided variables in the model but not the dataframe name
#try this instead
names(model.frame(forest_model))
#the output from this tells us we can extract data used to 
#fit the model directly from the model itself.
dat <- model.frame(forest_model)

head(dat)
#then create site level plots
library(ggplot2)

ggplot(dat,
       aes(x = years_since,
           y = trend,
           group = plot_id,
           colour = forest_type)) +
  
  geom_line() +
  geom_point() +
  
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  
  facet_wrap(~plot_id) +
  
  labs(
    x = "Years since mortality",
    y = "NDVI trend",
    colour = "Forest type"
  ) +
  
  theme_classic()



library(dplyr)

site_df <- trend %>%
  mutate(
    years_since = year - mort_year_observed
  )

#This doesnt work well because we have too many sites.
#So lets use this code which overlays all trajectories
ggplot(dat,
       aes(x = years_since,
           y = trend,
           group = plot_id,
           colour = forest_type)) +
  
  geom_line(alpha = 0.5) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  
  theme_classic()

#different style
ggplot(dat,
       aes(x = years_since,
           y = trend,
           group = plot_id,
           colour = forest_type)) +
  geom_line(alpha = 0.5) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = 0,
             linetype = "dotted") +
  theme_classic(base_size = 8)

#we want to chart all the trajectories for each forest type on one chart. 
#Mean trajectory as a solid line and individual sites as dotted lines
dat <- as.data.frame(model.frame(forest_model))

ggplot(dat,
       aes(x = years_since,
           y = trend)) +
  
  # Individual sites
  geom_smooth(
    aes(group = plot_id),
    method = "loess",
    se = FALSE,
    linetype = "dotted",
    linewidth = 0.75,
    colour = "grey50"
  ) +
  
  # Mean forest-type trajectory
  geom_smooth(
    aes(colour = forest_type),
    method = "loess",
    se = TRUE,
    linewidth = 1.0
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  
  facet_wrap(~forest_type) +
  
  labs(
    x = "Years since mortality",
    y = "NDVI trend",
    colour = "Forest type"
  ) +
  
  theme_classic()

#this version calculates mean trajectory for each forest type from mixed model 
#and overlays raw site trajectories
ggplot(dat,
       aes(years_since, trend)) +
  
  geom_line(
    aes(group = plot_id),
    colour = "grey70",
    linetype = "dotted",
    alpha = 0.6
  ) +
  
  geom_smooth(
    aes(colour = forest_type),
    method = "lm",
    formula = y ~ poly(x, 2),
    se = TRUE,
    linewidth = 1.0
  ) +
  
  facet_wrap(~forest_type) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  
  theme_classic()

#we are going to rerun the first chart as separate charts with thin solid grey lines
#for indiv sites and thick coloured smooth for forest type mean
#first extract the model data
dat <- as.data.frame(model.frame(forest_model))
#then run separate chart for each forest type, first conifer

library(dplyr)
library(ggplot2)

conifer <- dat %>%
  filter(forest_type == "Conifer")

ggplot(conifer,
       aes(x = years_since,
           y = trend)) +
  
  geom_smooth(
    aes(group = plot_id),
    method = "loess",
    se = FALSE,
    colour = "grey60",
    linewidth = 0.3
  ) +
  
  geom_smooth(
    colour = "darkorange",
    fill = "darkorange",
    method = "loess",
    se = TRUE,
    linewidth = 1.0
  ) +
  
  geom_vline(xintercept = 0,
             linetype = "solid") +
  annotate(
    "text",
    x = 0,
    y = Inf,
    label = "Observed Mortality",
    vjust = -0.5,angle = 90, hjust = 1.25, size = 3,
    colour = "grey30"
  ) +
  scale_x_continuous(
    breaks = seq(-4, 5, by = 1)
  ) +
  labs(
    #title = "Conifer Forests",
    x = "Years pre and post observed mortality",
    y = "Detrended NDVI"
  ) +
  
  theme_classic(base_size = 8)

#next deciduous
deciduous <- dat %>%
  filter(forest_type == "Deciduous")

ggplot(deciduous,
       aes(x = years_since,
           y = trend)) +
  
  geom_smooth(
    aes(group = plot_id),
    method = "loess",
    se = FALSE,
    colour = "grey60",
    linewidth = 0.3
  ) +
  
  geom_smooth(
    colour = "darkgreen",
    fill = "darkgreen",
    method = "loess",
    se = TRUE,
    linewidth = 1.0
  ) +
  
  geom_vline(xintercept = 0,
             linetype = "solid") +
  annotate(
    "text",
    x = 0,
    y = Inf,
    label = "Observed Mortality",
    vjust = -0.5,angle = 90,hjust = 1.25, size = 3,
    colour = "grey30"
  ) +
  scale_x_continuous(
    breaks = seq(-4, 5, by = 1)
  ) +
  
  labs(
    #title = "Deciduous Forests",
    x = "Years pre and post observed mortality",
    y = "Detrended NDVI"
  ) +
  
  theme_classic(base_size = 9)

#finally mixed
mixed <- dat %>%
  filter(forest_type == "Mixed")

ggplot(mixed,
       aes(x = years_since,
           y = trend)) +
  
  geom_smooth(
    aes(group = plot_id),
    method = "loess",
    se = FALSE,
    colour = "grey60",
    linewidth = 0.3
  ) +
  
  geom_smooth(
    colour = "blue",
    method = "loess",
    se = TRUE
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "solid"
  ) +
  
  annotate(
    "text",
    x = 0,
    y = Inf,
    label = "Observed Mortality",
    vjust = -0.5,
    angle = 90,
    hjust = 1.25,
    size = 3,
    colour = "grey30"
  ) +
  
  scale_x_continuous(
    breaks = seq(-4, 5, by = 1)
  ) +
  
  labs(
    x = "Years pre and post observed mortality",
    y = "Detrended NDVI"
  ) +
  
  theme_classic(base_size = 9)



#lets try plotting raw data with smoothed population trend
ggplot(mixed,
       aes(x = years_since,
           y = trend)) +
  
  geom_line(
    aes(group = plot_id),
    colour = "grey70",
    linewidth = 0.3,
    alpha = 0.5
  ) +
  
  geom_point(
    aes(group = plot_id),
    colour = "grey50",
    size = 0.8,
    alpha = 0.5
  ) +
  
  geom_smooth(
    colour = "blue",
    method = "loess",
    se = TRUE
  ) +
  geom_vline(xintercept = 0,
             linetype = "solid") +
  annotate(
    "text",
    x = 0,
    y = Inf,
    label = "Mortality event",
    vjust = -0.5,angle = 90,hjust = 1,size = 3,
    colour = "grey30"
  ) +
  scale_x_continuous(
    breaks = seq(-4, 5, by = 1)
  ) +
  
  labs(
    #title = "Mixed Forests",
    x = "Years since mortality",
    y = "Detrended NDVI trend"
  ) +
  
  theme_classic(base_size = 8)

#looks like the incomplete indiv plot lines are loess smoothing artefacr
#not missing data.
#Lets check number of observations for each plot
mixed %>%
  group_by(plot_id) %>%
  summarise(n_years = n()) %>%
  count(n_years)
#and
mixed %>%
  group_by(plot_id) %>%
  summarise(
    first = min(years_since),
    last = max(years_since),
    n_years = n()
  )

#MISSING DATA? Lets check if data is missing or if its an artefact
ggplot(dat,
       aes(years_since, trend, group = plot_id)) +
  geom_line(alpha = 0.3) +
  geom_point(size = 0.5)

#now lets differentiate between high and low mortality sites
#A figure showing high-mortality (>40%) versus low-mortality (<40%) sites 
#could be very informative and may reveal whether severely 
#affected stands show different recovery trajectories.
dat <- as.data.frame(model.frame(forest_model))
names(dat)
#pct_mort was not included in the model so we need to find original dataset
ls() #list all files in directory
#think it might be df so now we check for pct_mort
names(df)
#now we create mortality class - we are using the 
#statistically established threshold of 60% 
df <- df %>%
  mutate(
    mortality_class = ifelse(
      pct_mort > 60,
      "High mortality",
      "Low mortality"
    ))
    
#plot high (red) vs low (green)
ggplot(df,
           aes(x = years_since,
               y = trend)) +
      
      geom_smooth(
        aes(group = plot_id,
            colour = mortality_class),
        method = "loess",
        se = FALSE,
        linewidth = 0.3
      ) +
      
      geom_smooth(
        aes(colour = mortality_class,
            fill = mortality_class),
        method = "loess",
        se = TRUE,
        linewidth = 1.0
      ) +
      
      scale_colour_manual(values = c(
        "High mortality" = "red",
        "Low mortality" = "darkgreen"
      )) +
      
      scale_fill_manual(values = c(
        "High mortality" = "red",
        "Low mortality" = "darkgreen"
      )) +
      
      geom_vline(xintercept = 0,
                 linetype = "solid") +
      
      theme_classic()
#we have an error message to solve - cant find mortality_class
names(df)
#then we create the mortality class again
library(dplyr)

df <- df %>%
  mutate(
    mortality_class = ifelse(
      pct_mort > 60,
      "High mortality",
      "Low mortality"
    )
  )
#check it worked
table(df$mortality_class)
#then plot = all forest types combined
ggplot(df,
       aes(x = years_since,
           y = trend)) +
  
  geom_smooth(
    aes(group = plot_id,
        colour = mortality_class),
    method = "loess",
    se = FALSE,
    linewidth = 0.3
  ) +
  
  geom_smooth(
    aes(colour = mortality_class,
        fill = mortality_class),
    method = "loess",
    se = TRUE,
    linewidth = 1.2
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "solid"
  ) +
  
  scale_colour_manual(values = c(
    "High mortality" = "red3",
    "Low mortality" = "darkgreen"
  )) +
  
  scale_fill_manual(values = c(
    "High mortality" = "red3",
    "Low mortality" = "darkgreen"
  )) +
  
  labs(
    x = "Years pre and post observed mortality",
    y = "Detrended NDVI",
    colour = "Mortality class",
    fill = "Mortality class"
  ) +
  
  theme_classic(base_size = 9)

#we need to check if this data includes all sites
table(df$site_type)
#this tells us there are contorl and mortality types
#lets check unique sites
unique(df$site_type)
#check how defined
table(df$site_type, df$mortality_class)
#this tells us they are combined mortality and control site plots
table(df$site_type)
#IT CONTAINS ALL COMBINED DATA SO LETS PROCEED


#that is too messy so preferred is separate plots per forest type
ggplot(df,
       aes(x = years_since,
           y = trend)) +
  
  geom_smooth(
    aes(group = plot_id,
        colour = mortality_class),
    method = "loess",
    se = FALSE,
    linewidth = 0.3
  ) +
  
  geom_smooth(
    aes(colour = mortality_class,
        fill = mortality_class),
    method = "loess",
    se = TRUE,
    linewidth = 1.0
  ) +
  
  facet_wrap(~forest_type) +
  
  geom_vline(
    xintercept = 0,
    linetype = "solid"
  ) +
  
  scale_colour_manual(values = c(
    "High mortality" = "red3",
    "Low mortality" = "darkgreen"
  )) +
  
  scale_fill_manual(values = c(
    "High mortality" = "red3",
    "Low mortality" = "darkgreen"
  )) +
  
  labs(
    x = "Years pre and post observed mortality",
    y = "Detrended NDVI",
    colour = "Mortality class",
    fill = "Mortality class"
  ) +
  
  theme_classic(base_size = 9)

#we need to restrict data to the 4 years pre and 5 year post mortality
#we will create a new 'df_short'
df_short <- df %>%
  filter(years_since >= -4,
         years_since <= 5)
#then use df_short in the plotting code 
#lets check data being plotted
mixed <- df_short %>%
  filter(forest_type == "Mixed")

sort(unique(mixed$years_since))
#this shows the data only goes back to -2 years instead of -4
#now lets check where the years disappear
table(df$years_since)
#how many obs for each year and mortality class
mixed %>%
  count(years_since, mortality_class)
table(df_short$years_since)
df_short %>%
  group_by(forest_type) %>%
  summarise(
    min_year = min(years_since),
    max_year = max(years_since),
    n = n()
  )
#first lets check how many plots in df_short
df_short %>%
  group_by(site_type) %>%
  summarise(n_plots = n_distinct(plot_id))
#OK that confirms we have all sites in there:
#control = 59
#mortality = 73

#run separate charts for each forest type - firstly conifer
#try replacing mean loess curve with model consistent quadratic smoother
library(dplyr)
library(ggplot2)

conifer <- df_short %>%
  filter(forest_type == "Conifer")

ggplot(conifer,
       aes(x = years_since,
           y = trend)) +
  
  geom_smooth(
    aes(group = plot_id,
        colour = mortality_class),
    method = "loess",
    se = FALSE,
    linewidth = 0.2
  ) +
  
  geom_smooth(
    aes(colour = mortality_class,
        fill = mortality_class),
    method = "loess",
    se = TRUE,
    linewidth = 1.0,
    alpha = 0.2
  ) +
  
  geom_vline(xintercept = 0,
             linetype = "solid") +
  geom_vline(
    xintercept = 0,
    linetype = "solid",
    colour = "grey60",
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = 0,
    y = Inf,
    label = "Observed Mortality",
    vjust = -0.5,
    angle = 90,
    hjust = 1.25,
    size = 3,
    colour = "grey20"
  ) +
  scale_x_continuous(
    breaks = seq(-4, 5, by = 1)
  ) +
  
  scale_colour_manual(values = c(
    "High mortality" = "red3",
    "Low mortality" = "darkgreen"
  )) +
  
  scale_fill_manual(values = c(
    "High mortality" = "red3",
    "Low mortality" = "darkgreen"
  )) +
  
  labs(
    x = "Years pre and post observed mortality",
    y = "Detrended NDVI",
    colour = "Mortality class",
    fill = "Mortality class",
    subtitle = "Conifer forest "
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(
      fill = "white",
      colour = "black"
    ),
    strip.text = element_text(face = "bold")
  )


#then DECIDUOUS
deciduous <- df_short %>%
  filter(forest_type == "Deciduous")

ggplot(deciduous,
       aes(x = years_since,
           y = trend)) +
  
  geom_smooth(
    aes(group = plot_id,
        colour = mortality_class),
    method = "loess",
    se = FALSE,
    linewidth = 0.3
  ) +
  
  geom_smooth(
    aes(colour = mortality_class,
        fill = mortality_class),
    method = "loess",
    se = TRUE,
    linewidth = 1.0
  ) +
  
  geom_vline(xintercept = 0,
             linetype = "dotted") +
  geom_vline(xintercept = 0,
             linetype = "dotted") +
  geom_vline(
    xintercept = 0,
    linetype = "solid",
    colour = "grey40",
    linewidth = 0.3
  ) +
  annotate(
    "text",
    x = 0,
    y = Inf,
    label = "Observed Mortality",
    vjust = -0.5,
    angle = 90,
    hjust = 1.25,
    size = 3,
    colour = "grey30"
  ) +
  scale_x_continuous(
    breaks = seq(-4, 5, by = 1)
  ) +
  
  scale_colour_manual(values = c(
    "High mortality" = "red3",
    "Low mortality" = "darkgreen"
  )) +
  
  scale_fill_manual(values = c(
    "High mortality" = "red3",
    "Low mortality" = "darkgreen"
  )) +
  
  labs(
    title = "Deciduous Forest",
    x = "Years pre and post observed mortality",
    y = "Detrended NDVI"
  ) +
  
  theme_classic(base_size = 9)+
  theme(
  legend.position = "bottom",
  strip.background = element_rect(
    fill = "white",
    colour = "black"
  ),
  strip.text = element_text(face = "bold")
)

#and MIXED forests
mixed <- df_short %>%
  filter(forest_type == "Mixed")

ggplot(mixed,
       aes(x = years_since,
           y = trend)) +
  
  geom_smooth(
    aes(group = plot_id,
        colour = mortality_class),
    method = "loess",
    se = FALSE,
    linewidth = 0.3
  ) +
  
  geom_smooth(
    aes(colour = mortality_class,
        fill = mortality_class),
    method = "loess",
    se = TRUE,
    linewidth = 1.0
  ) +
  
  geom_vline(xintercept = 0,
             linetype = "dotted") +
  geom_vline(xintercept = 0,
             linetype = "dotted") +
  geom_vline(
    xintercept = 0,
    linetype = "solid",
    colour = "grey40",
    linewidth = 0.4
  ) +
  annotate(
    "text",
    x = 0,
    y = Inf,
    label = "Observed Mortality",
    vjust = -0.5,
    angle = 90,
    hjust = 2.0,
    size = 2.5,
    colour = "grey30"
  ) +
  scale_x_continuous(
    breaks = seq(-4, 5, by = 1)
  ) +
  
  scale_colour_manual(values = c(
    "High mortality" = "red3",
    "Low mortality" = "darkgreen"
  )) +
  
  scale_fill_manual(values = c(
    "High mortality" = "red3",
    "Low mortality" = "darkgreen"
  )) +
  
  labs(
    title = "Mixed Forests",
    x = "Years pre and post observed mortality",
    y = "Detrended NDVI"
  ) +
  
  theme_classic(base_size = 9)+
  theme(
    legend.position = "bottom",
    strip.background = element_rect(
      fill = "white",
      colour = "black"
    ),
    strip.text = element_text(face = "bold")
  )


#Original analysis was only run on mortality category sites
#RERUN ANALYSIS WITH MORTALITY AND CONTROL SITES
quad_model_all <- lmer(
  trend ~ years_since + I(years_since^2) + (1http://127.0.0.1:45802/graphics/24241dd9-699d-4109-bf89-2d0a474dd6a2.png | plot_id),
  data = df_short
)
#to determine whether there is difference in trends between
#mortality and control sites we should fit:
quad_model_type <- lmer(
  trend ~ (years_since + I(years_since^2)) * site_type +
    (1 | plot_id),
  data = df_short
)
#we can then allow different intercepts for control v mortality sites
#different linear trends
#different quadratic trends
#we can then test whether site type matters by comparing models, for example:
anova(quad_model_all, quad_model_type)
#or by examining the coefficients in the quad model type
summary(quad_model_type)
#these highly significant results (p value <0.0001)is 
#telling us that the site type contains important information 
#that was being ignored in the mortality-only model

#next step look at the model summary
summary(quad_model_type)
#Mixed-effects modelling indicated that temporal trends 
#differed between mortality and control sites. 
#Mortality sites had lower overall trend values 
#than control sites and exhibited significantly different 
#temporal responses, with a steeper linear decline and 
#stronger quadratic curvature. 
#Model comparison showed that including site type and 
#its interactions with time significantly improved model fit 
#(χ² = 24.17, df = 3, p < 0.001), indicating that the 
#relationship between trend and years since event 
#was not consistent across site types. 
#The strongest effects were associated with 
#quadratic terms, suggesting that temporal dynamics 
#were predominantly nonlinear and that mortality sites 
#exhibited a more pronounced curved trajectory through time 
#than control sites.

#lets plot the fitted values to see these differences
library(ggplot2)

ggplot(df_short,
       aes(years_since, trend, colour = site_type)) +
  geom_point(alpha = 0.1) +
  geom_smooth(
    method = "lm",
    formula = y ~ x + I(x^2),
    se = TRUE
  )
#this is a very confusing plot with all the raw data too 
#so lets plot only model-fitted curves
library(ggeffects)
library(ggplot2)

pred <- ggpredict(
  quad_model_type,
  terms = c("years_since", "site_type")
)

ggplot(pred,
       aes(x = x, y = predicted, colour = group)) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high,
                  fill = group),
              alpha = 0.2,
              colour = NA) +
  labs(
    x = "Years since",
    y = "Trend",
    colour = "Site type",
    fill = "Site type"
  ) +
  theme_classic()


#lets swap the colours over so the mortality group is red
library(ggeffects)
library(ggplot2)

pred <- ggpredict(
  quad_model_type,
  terms = c("years_since", "site_type")
)

ggplot(pred,
       aes(x = x, y = predicted, colour = group)) +
  geom_line(linewidth = 1.0) +
  geom_ribbon(
    aes(ymin = conf.low,
        ymax = conf.high,
        fill = group),
    alpha = 0.2,
    colour = NA
  ) +
  scale_colour_manual(
    values = c(
      "Control" = "steelblue",
      "Mortality" = "firebrick"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Control" = "steelblue",
      "Mortality" = "firebrick"
    )
  ) +
  labs(
    x = "Years since",
    y = "Trend",
    title = "Control v mortality sites",
    colour = "Site type",
    fill = "Site type"
  ) +
  theme_classic(base_size = 8)
  
  #a few tweaks to the plot e.g. legend at the bottom aaxis label change
  #and larger font size
  library(ggeffects)
  library(ggplot2)
  
  pred <- ggpredict(
    quad_model_type,
    terms = c("years_since", "site_type")
  )
  
  ggplot(pred,
         aes(x = x, y = predicted, colour = group)) +
    geom_line(linewidth = 1.0) +
    geom_ribbon(
      aes(ymin = conf.low,
          ymax = conf.high,
          fill = group),
      alpha = 0.2,
      colour = NA
    ) +
    
    geom_vline(
      xintercept = 0,
      linetype = "solid",
      colour = "grey40",
      linewidth = 0.4
    ) +
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = "Observed Mortality",
      vjust = -0.5,
      angle = 90,
      hjust = 1.25,
      size = 3,
      colour = "grey30"
    ) +
    scale_x_continuous(
      breaks = seq(-4, 5, by = 1)
    ) +
    scale_colour_manual(
      values = c(
        "Control" = "steelblue",
        "Mortality" = "firebrick"
      )
    ) +
    scale_fill_manual(
      values = c(
        "Control" = "steelblue",
        "Mortality" = "firebrick"
      )
    ) +
    labs(
      x = "Years pre and post observed mortality",
      y = "Detrended NDVI",
      title = 'Control v mortality',
      colour = "Site type",
      fill = "Site type"
    ) +
    theme_classic(base_size = 8) +
    theme(
      legend.position = "bottom"
      )
  
#same chart but with empirical data
  library(ggplot2)
  
  ggplot(window_df,
         aes(x = years_since,
             y = trend,
             colour = site_type,
             fill = site_type)) +
    
    stat_summary(
      fun = mean,
      geom = "line",
      linewidth = 1.0
    ) +
    
    stat_summary(
      fun.data = mean_se,
      geom = "ribbon",
      alpha = 0.2,
      colour = NA
    ) +
    
    scale_colour_manual(
      values = c(
        "Control" = "steelblue",
        "Mortality" = "firebrick"
      )
    ) +
    
    scale_fill_manual(
      values = c(
        "Control" = "steelblue",
        "Mortality" = "firebrick"
      )
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "solid",
      colour = "grey40",
      linewidth = 0.4
    ) +
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = "Observed Mortality",
      vjust = -0.5,
      angle = 90,
      hjust = 1,
      size = 2.5,
      colour = "grey30"
    ) +
    scale_x_continuous(
      breaks = seq(-4, 5, by = 1)
    ) +
    
    labs(
      x = "Years since mortality",
      y = "Detrended NDVI",
      colour = "Site type",
      fill = "Site type"
    ) +
    
    theme_classic(base_size = 8) +
    
    theme(
      legend.position = "bottom",
      )  
  
#or how about mean curve by group with individual profiles in grey?
  ggplot() +
    
    geom_line(
      data = window_df,
      aes(years_since, trend, group = plot_id),
      colour = "grey80",
      alpha = 0.3
    ) +
    
    stat_summary(
      data = window_df,
      aes(years_since, trend,
          colour = site_type),
      fun = mean,
      geom = "line",
      linewidth = 1.0
    ) +
    
    scale_colour_manual(
      values = c(
        "Control" = "steelblue",
        "Mortality" = "firebrick"
      )
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "solid",
      colour = "grey40",
      linewidth = 0.4
    ) +
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = "Observed Mortality",
      vjust = -0.5,
      angle = 90,
      hjust = 1,
      size = 2.5,
      colour = "grey30"
    ) +
    scale_x_continuous(
      breaks = seq(-4, 5, by = 1)
    ) +
    
    labs(
      title = "Control v mortality",
      x = "Years since mortality",
      y = "Detrended NDVI",
      colour = "Site type"
    ) +
    
    theme_classic(base_size = 8) +
    
    theme(
      legend.position = "bottom",
          )
  
  #remove vertical lines and smooth the individual curves
  ggplot() +
    
    geom_smooth(
      data = window_df,
      aes(years_since, trend, group = plot_id),
      method = "loess",
      se = FALSE,
      colour = "grey80",
      linewidth = 0.3,
      alpha = 0.3
    ) +
    
    stat_summary(
      data = window_df,
      aes(years_since, trend, colour = site_type),
      fun = mean,
      geom = "line",
      linewidth = 1.0
    ) +
    
    scale_colour_manual(
      values = c(
        "Control" = "steelblue",
        "Mortality" = "firebrick"
      )
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "solid",
      colour = "grey40",
      linewidth = 0.4
    ) +
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = "Observed Mortality",
      vjust = -0.5,
      angle = 90,
      hjust = 1,
      size = 2.5,
      colour = "grey30"
    ) +
    scale_x_continuous(
      breaks = seq(-4, 5, by = 1)
    ) +
    
    labs(
      x = "Years since mortality",
      y = "Detrended NDVI",
      Title = "Control v mortality",
      colour = "Site type"
    ) +
    
    theme_classic(base_size = 8) +
    
    theme(
      legend.position = "bottom",
       )
  
  #smooth the mean lines and the individual site lines
  ggplot(window_df,
         aes(years_since, trend)) +
    
    # Individual site trajectories
    geom_line(
      aes(group = plot_id),
      colour = "grey80",
      alpha = 0.2
    ) +
    
    # Smoothed empirical means
    geom_smooth(
      aes(colour = site_type,
          fill = site_type),
      method = "loess",
      se = TRUE,
      linewidth = 1.0
    ) +
    
    scale_colour_manual(
      values = c(
        "Control" = "steelblue",
        "Mortality" = "firebrick"
      )
    ) +
    
    scale_fill_manual(
      values = c(
        "Control" = "steelblue",
        "Mortality" = "firebrick"
      )
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "solid",
      colour = "grey40",
      linewidth = 0.4
    ) +
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = "Observed Mortality",
      vjust = -0.5,
      angle = 90,
      hjust = 1,
      size = 2.5,
      colour = "grey30"
    ) +
    scale_x_continuous(
      breaks = seq(-4, 5, by = 1)
    ) +
    
    labs(
      x = "Years since mortality",
      y = "Detrended NDVI",
      title = "Control v mortality",
      colour = "Site type",
      fill = "Site type"
    ) +
    
    theme_classic(base_size = 8) +
    
    theme(
      legend.position = "bottom",
      )
  
  #TRY AGAIN
  ggplot(window_df,
         aes(years_since, trend)) +
    
    # Smoothed individual site trajectories
    geom_smooth(
      aes(group = plot_id),
      method = "loess",
      se = FALSE,
      colour = "grey85",
      linewidth = 0.3,
      alpha = 0.15
    ) +
    
    # Smoothed empirical mean curves by site type
    geom_smooth(
      aes(colour = site_type,
          fill = site_type),
      method = "loess",
      se = TRUE,
      linewidth = 1.0
    ) +
    
    scale_colour_manual(
      values = c(
        "Control" = "steelblue",
        "Mortality" = "firebrick"
      )
    ) +
    
    scale_fill_manual(
      values = c(
        "Control" = "steelblue",
        "Mortality" = "firebrick"
      )
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "solid",
      colour = "grey40",
      linewidth = 0.4
    ) +
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = "Observed Mortality",
      vjust = -0.5,
      angle = 90,
      hjust = 1,
      size = 2.5,
      colour = "grey30"
    ) +
    scale_x_continuous(
      breaks = seq(-4, 5, by = 1)
    ) +
    
    labs(
      x = "Years since mortality",
      y = "Detrended NDVI",
      title = "Control v mortality sites",
      colour = "Site type",
      fill = "Site type"
    ) +
    
    theme_classic(base_size = 8) +
    
    theme(
      legend.position = "bottom",
      )
  
  #INCLUDE ALL SITE TYPES _ MORTALITY AND CONTROL
  #now lets rerun the analysis with all sites included
  #we simply need to change database from 'mort_window' to 'window_df'
  
  # ---- 2.4 Quadratic mixed model: is there a dip then a rise? --------------
  # Fixed effects: years_since and years_since^2 (the squared term lets the
  # model curve, rather than assuming a straight line, so it can capture a dip
  # followed by a rise).
  # Random effect: a random intercept per plot_id, since each plot contributes
  # many repeated observations over time.
  
  #mort_window <- window_df %>% filter(site_type == "Mortality")
  
  quad_model <- lmer(
    trend ~ years_since + I(years_since^2) + (1 | plot_id),
    data = window_df
  )
  print(summary(quad_model))
  

  # The vertex (turning point) of the fitted curve tells you where, in years
  # relative to the event, the model estimates the lowest point of trend NDVI.
  b <- fixef(quad_model)
  vertex <- -b["years_since"] / (2 * b["I(years_since^2)"])
  cat("Estimated turning point (years since event):", round(vertex, 2), "\n")
 
  # Including all site types has reduced the turning point (lowest NDVI)
  #from 0.47 to 0.33
  
  
  #lets find the reference category for forest type so we can calculte vertex by forest type
  fixef(forest_model)
  
  #by forest type
  b <- fixef(forest_model)
  vertex <- -b["years_since"] / (2 * b["I(years_since^2)"])
  cat("Estimated turning point (years since event):", round(vertex, 2), "\n")
  
  #by climate type
  b <- fixef(climate_model)
  vertex <- -b["years_since"] / (2 * b["I(years_since^2)"])
  cat("Estimated turning point (years since event):", round(vertex, 2), "\n")
  
  # ---- 2.5 By climate and forest type: immediate change and recovery time --
  # For each climate x forest type combination with mortality plots, this
  # computes:
  #   - immediate_change: mean trend in years 0-1 minus mean trend in years -2,-1
  #   - control_baseline: mean trend of control plots in the same climate x
  #     forest type combination, over the same window
  #   - recovery_year: the first year (0, 1, 2, ...) at which the mortality
  #     plots' mean trend reaches or exceeds the control baseline. Requires at
  #     least 2 plots contributing at that year to count, and returns NA if the
  #     control baseline is never reached within the available data.
 
   
  combos <- window_df %>%
    #filter(site_type == "Mortality") %>%
    distinct(climate, forest_type)
  
  results <- list()
  
  for (i in seq_len(nrow(combos))) {
    clim  <- combos$climate[i]
    ftype <- combos$forest_type[i]
    
    mort_sub <- window_df %>%
      filter(site_type == "Mortality", climate == clim, forest_type == ftype)
    ctrl_sub <- window_df %>%
      filter(site_type == "Control", climate == clim, forest_type == ftype)
    
    n_mort_plots <- n_distinct(mort_sub$plot_id)
    n_ctrl_plots <- n_distinct(ctrl_sub$plot_id)
    
    pre_mean  <- mort_sub %>% filter(years_since %in% c(-2, -1)) %>% pull(trend) %>% mean(na.rm = TRUE)
    post_mean <- mort_sub %>% filter(years_since %in% c(0, 1))   %>% pull(trend) %>% mean(na.rm = TRUE)
    immediate_change <- post_mean - pre_mean
    
    if (n_ctrl_plots > 0) {
      control_baseline <- mean(ctrl_sub$trend, na.rm = TRUE)
      
      recovery_year <- NA_real_
      years_available <- sort(unique(mort_sub$years_since[mort_sub$years_since >= 0]))
      for (yr in years_available) {
        yr_data <- mort_sub %>% filter(years_since == yr)
        if (n_distinct(yr_data$plot_id) < 2) next          # skip years with too few plots
        if (mean(yr_data$trend, na.rm = TRUE) >= control_baseline) {
          recovery_year <- yr
          break
        }
      }
    } else {
      control_baseline <- NA_real_
      recovery_year <- NA_real_   # no control plots exist for this combination
    }
    
    results[[i]] <- tibble(
      climate = clim,
      forest_type = ftype,
      n_mortality_plots = n_mort_plots,
      n_control_plots = n_ctrl_plots,
      immediate_change = round(immediate_change, 4),
      control_baseline = round(control_baseline, 4),
      recovery_year_since_event = recovery_year
    )
  }
  
  summary_table <- bind_rows(results) %>% arrange(immediate_change)
  print(summary_table, n = Inf)
  
  #=====================================================================
  #LETS START AGAIN USING ALL SITE TYPES INSTEAD OF ONLY MORTALITY SITES
  #first we run the overall quadratic model (all sites)
  quad_model_all <- lmer(
    trend ~ years_since + I(years_since^2) +
      (1 | plot_id),
    data = df_short
  )
  
  summary(quad_model_all)
  
  #calculate the turning point (lowest point of ndvi post mortality)
  b <- fixef(quad_model_all)
  
  vertex <- -b["years_since"] /
    (2 * b["I(years_since^2)"])
  
  cat("Estimated turning point:",
      round(vertex, 2), "years\n")
  #this computes as 0.33 years (compred with 0.47 years for mortality sites)
  
  #test whether site type matters
  site_model <- lmer(
    trend ~ (years_since + I(years_since^2)) *
      site_type +
      (1 | plot_id),
    data = window_df
  )
  
  summary(site_model)
  
  anova(quad_model_all, site_model)
  #this reveals statistically significant difference between site model and quad mdel all
  #in other words type of site makes a difference to mort response
  #p = 0.001
  
  #forest type model (all sites)
  forest_model <- lmer(
    trend ~ (years_since + I(years_since^2)) *
      forest_type +
      site_type +
      (1 | plot_id),
    data = window_df
  )
  
  summary(forest_model)
  #test whether adding forest tyep improves fit
  forest_null <- lmer(
    trend ~ years_since + I(years_since^2) +
      site_type +
      (1 | plot_id),
    data = window_df
  )
  
  anova(forest_null, forest_model)
  #significant difference (p=0.001) between null and forest model
  
  #now we run climate model (all sites)
  climate_model <- lmer(
    trend ~ (years_since + I(years_since^2)) *
      climate +
      site_type +
      (1 | plot_id),
    data = window_df
  )
  
  summary(climate_model)
#test if it improves fit
  climate_null <- lmer(
    trend ~ years_since + I(years_since^2) +
      site_type +
      (1 | plot_id),
    data = window_df
  )
  
  anova(climate_null, climate_model)
  #significant difference (p=0.001) between null and climate model
  
  #mortality severity model (all sites)
  mortality_model <- lmer(
    trend ~ (years_since + I(years_since^2)) *
      mort_group +
      site_type +
      (1 | plot_id),
    data = df_short
  )
  
  summary(mortality_model)
  
  #forest type specific turning points
  fixef(forest_model)
  #reference forest type
  b <- fixef(forest_model)
  
  vertex_ref <- -b["years_since"] /
    (2 * b["I(years_since^2)"])
  
  vertex_ref
  
  #to discover which is the reference category 
  levels(window_df$forest_type)
  
  #turning point for deciduous forest tyep
  b <- fixef(forest_model)
  
  lin <- b["years_since"] +
    b["years_since:forest_typeDeciduous"]
  
  quad <- b["I(years_since^2)"] +
    b["I(years_since^2):forest_typeDeciduous"]
  
  vertex_deciduous <- -lin / (2 * quad)
  
  vertex_deciduous
  #the turning point for deciduous forests is -2.194 years
  
  #turning point for conifer forest type
  b <- fixef(forest_model)
  
  lin <- b["years_since"] +
    b["years_since:forest_typeConifer"]
  
  quad <- b["I(years_since^2)"] +
    b["I(years_since^2):forest_typeConifer"]
  
  vertex_conifer <- -lin / (2 * quad)
  
  vertex_conifer
  #the turning point for conifer forests is N/A
  
  #turning point for mixed forests
  b <- fixef(forest_model)
  
  lin <- b["years_since"] +
    b["years_since:forest_typeMixed"]
  
  quad <- b["I(years_since^2)"] +
    b["I(years_since^2):forest_typeMixed"]
  
  vertex_mixed <- -lin / (2 * quad)
  
  vertex_mixed
  #the turning point for mixed forests is 1.58 years
  
  #better approach is automatically calculate for all forest types:
  
  # A tibble: 3 × 4
  
 #forest_type linear       quadratic    turning_point
  #Conifer    -0.000153     0.000376          0.20
  #Deciduous   0.003036     0.000592         -2.57
  #Mixed      -0.003284     0.001035          1.59
  
  #Conifer reaches its minimum almost immediately after the event (~0.2 years).
  #Mixed forests reach their minimum around 1.6 years after the event, 
  #suggesting a more prolonged response.
  #Deciduous forests have a mathematically estimated 
  #minimum at -2.6 years, which is before the event. 
  #In practice this usually means there is no 
  #post-event decline evident in the fitted curve 
  #and the trajectory is generally increasing 
  #across the observed period.
  
  #recommend plotting the fitted curves by forest type 
  #for visual verification and interpretation
  #use ggeffects
  library(ggeffects)
  library(ggplot2)
  
  pred_forest <- ggpredict(
    forest_model,
    terms = c("years_since", "forest_type")
  )
  
  ggplot(pred_forest,
         aes(x = x,
             y = predicted,
             colour = group,
             fill = group)) +
    
    geom_ribbon(
      aes(ymin = conf.low,
          ymax = conf.high),
      alpha = 0.2,
      colour = NA
    ) +
    
    geom_line(linewidth = 1.0) +
    
    scale_colour_manual(
      values = c(
        "Conifer" = "darkorange",
        "Deciduous" = "darkgreen",
        "Mixed" = "blue"
      )
    ) +
    
    scale_fill_manual(
      values = c(
        "Conifer" = "darkorange",
        "Deciduous" = "darkgreen",
        "Mixed" = "blue"
      )
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "solid",
      colour = "grey40",
      linewidth = 0.4
    ) +
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = "Observed Mortality",
      vjust = -0.5,
      angle = 90,
      hjust = 1,
      size = 2.5,
      colour = "grey30"
    ) +
    scale_x_continuous(
      breaks = seq(-4, 5, by = 1)
    ) +
    
    labs(
      x = "Years since mortality",
      y = "Detrended NDVI",
      colour = "Forest type",
      fill = "Forest type"
    ) +
    
    theme_classic(base_size = 8) +
    
    theme(
      legend.position = "bottom",
        )
  
  #add turning points to the plot, first create dataframe
  turning_points <- data.frame(
    forest_type = c("Conifer", "Deciduous", "Mixed"),
    turning_point = c(0.20, -2.57, 1.59)
  )
  #then obtain predicted values at each turning point from fitted curves
  library(dplyr)
  
  tp_plot <- turning_points %>%
    rowwise() %>%
    mutate(
      trend = pred_forest %>%
        filter(group == forest_type) %>%
        slice(which.min(abs(x - turning_point))) %>%
        pull(predicted)
    )
  #now add these layers to the forest type plot:
  ggplot(pred_forest,
         aes(x = x,
             y = predicted,
             colour = group,
             fill = group)) +
    
    geom_ribbon(
      aes(ymin = conf.low,
          ymax = conf.high),
      alpha = 0.2,
      colour = NA
    ) +
    
    geom_line(linewidth = 1.0) +
    
    geom_point(
      data = tp_plot,
      inherit.aes = FALSE,
      aes(
        x = turning_point,
        y = trend,
        colour = forest_type
      ),
      size = 2.0
    ) +
    
    geom_text(
      data = tp_plot,
      inherit.aes = FALSE,
      aes(
        x = turning_point,
        y = trend,
        label = round(turning_point, 2),
        colour = forest_type
      ),
      vjust = -0.8,
      size = 2.0
    ) +
    
    scale_colour_manual(
      values = c(
        "Conifer" = "darkorange",
        "Deciduous" = "darkgreen",
        "Mixed" = "blue"
      ),
      name = "Forest type"
    ) +
    
    scale_fill_manual(
      values = c(
        "Conifer" = "darkorange",
        "Deciduous" = "darkgreen",
        "Mixed" = "blue"
      ),
      name = "Forest type"
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "solid",
      colour = "grey40",
      linewidth = 0.4
    ) +
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = "Observed Mortality",
      vjust = -0.5,
      angle = 90,
      hjust = 1.25,
      size = 2.5,
      colour = "grey30"
    ) +
    scale_x_continuous(
      breaks = seq(-4, 5, by = 1)
    ) +
    
    labs(
      x = "Years pre and post observed mortality",
      y = "Detrended NDVI"
    ) +
    
    theme_bw(base_size = 8) +
    
    theme(
      legend.position = "bottom",
      )
  
  #Now we want to create a BOUNDS CHART for each forest type
  #this constitutes a trajectory plot showing the mean NDVI trend through time, 
  #upper and lower bounds (95% confidence intervals)
  #plotted separately for each forest type
  library(dplyr)
  library(ggplot2)
  
  forest_bounds <- df_short %>%
    group_by(forest_type, years_since) %>%
    summarise(
      mean_trend = mean(trend, na.rm = TRUE),
      se = sd(trend, na.rm = TRUE) / sqrt(n()),
      lower = mean_trend - 1.96 * se,
      upper = mean_trend + 1.96 * se,
      .groups = "drop"
    )
  
  ggplot(forest_bounds,
         aes(years_since, mean_trend)) +
    
    geom_ribbon(
      aes(ymin = lower,
          ymax = upper),
      fill = "lightblue",
      alpha = 0.3
    ) +
    
    geom_line(
      colour = "navy",
      linewidth = 1.2
    ) +
    
    facet_wrap(~ forest_type) +
    
    geom_vline(xintercept = 0,
               linetype = "solid") +
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = "Observed Mortality",
      vjust = -1.0,angle = 90,hjust = 1.25, size = 2.5,
      colour = "grey30"
    ) +
    scale_x_continuous(
      breaks = seq(-4, 5, by = 1)
    ) +
    
    labs(
      x = "Years since mortality",
      y = "Detrended NDVI"
    ) +
    
    theme_bw(base_size = 8) +
    
    theme(
      legend.position = "none",
      )
  
  #check number of sites in windows_df
  df_short %>%
    group_by(site_type) %>%
    summarise(n_plots = n_distinct(plot_id))
  
  #OK good it contains all 132 combined sites
  #control = 59
  #mortality = 73
  
  #we need to rerun model to check all sites are inlcuded
  forest_model <- lmer(
    trend ~ years_since + I(years_since^2) + forest_type +
      (1 | plot_id),
    data = df_short
  )
  
  #then regenerate predictions
  pred_forest <- ggpredict(
    forest_model,
    terms = c("years_since", "forest_type")
  )
  
  #check number of sites again
  n_distinct(df_short$plot_id)
  
  #and after fitting
  n_distinct(model.frame(forest_model)$plot_id)
  
  #now we can rerun the analysis
  
  #model based bounds will match our fitted forest type model
  library(ggeffects)
  
  pred_forest <- ggpredict(
    forest_model,
    terms = c("years_since", "forest_type")
  )
  
  ggplot(pred_forest,
         aes(x = x,
             y = predicted)) +
    
    geom_ribbon(
      aes(ymin = conf.low,
          ymax = conf.high),
      fill = "lightblue",
      alpha = 0.3
    ) +
    
    geom_line(
      colour = "navy",
      linewidth = 1.5
    ) +
    
    facet_wrap(~ group) +
    geom_vline(xintercept = 0,
               linetype = "solid") +
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = "Observed Mortality",
      vjust = -1.0,angle = 90,hjust = 1.25, size = 2.5,
      colour = "grey30"
    ) +
    scale_x_continuous(
      breaks = seq(-4, 5, by = 1)
    ) +
    
    labs(
      x = "Years since mortality",
      y = "Detrended NDVI"
    ) +
    
    theme_bw(base_size = 8) +
    
    
  #we just need to check whether this was run on all sites
  #or just mortality sites
  #we run this to check
  forest_model@frame %>%
    summarise(
      n_sites = n_distinct(plot_id),
      n_rows = n()
    )
  #n = 124 confirming we are using all sites
  
  #now lets colour the sites the usual forest type shade
  library(ggplot2)
  
  ggplot(pred_forest,
         aes(x = x,
             y = predicted,
             colour = group,
             fill = group)) +
    
    geom_ribbon(
      aes(ymin = conf.low,
          ymax = conf.high),
      alpha = 0.25,
      colour = NA
    ) +
    
    geom_line(
      linewidth = 1.5
    ) +
    geom_vline(xintercept = 0,
               linetype = "solid") +
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = "Observed Mortality",
      vjust = -1.0,angle = 90,hjust = 1.15,
      colour = "grey30"
    ) +
    scale_x_continuous(
      breaks = seq(-4, 5, by = 1)
    ) +
    
    facet_wrap(~ group) +
    
    scale_colour_manual(
      values = c(
        "Conifer" = "darkorange",
        "Mixed" = "blue",
        "Deciduous" = "darkgreen"
      )
    ) +
    
    scale_fill_manual(
      values = c(
        "Conifer" = "darkorange",
        "Mixed" = "blue",
        "Deciduous" = "darkgreen"
      )
    ) +
    
    labs(
      x = "Years pre and post mortality",
      y = "Detrended NDVI"
    ) +
    
    theme_bw(base_size = 8) +
    
    theme(
      legend.position = "none",
      )
  
  #Now lets run the bounds charts using the empirical (not modelled) data
  library(dplyr)
  library(ggplot2)
  
  forest_empirical <- window_df %>%
    group_by(forest_type, years_since) %>%
    summarise(
      mean_trend = mean(trend, na.rm = TRUE),
      se = sd(trend, na.rm = TRUE) / sqrt(n()),
      lower = mean_trend - 1.96 * se,
      upper = mean_trend + 1.96 * se,
      .groups = "drop"
    )
  
  ggplot(forest_empirical,
         aes(x = years_since,
             y = mean_trend,
             colour = forest_type,
             fill = forest_type)) +
    
    geom_ribbon(
      aes(ymin = lower,
          ymax = upper),
      alpha = 0.25,
      colour = NA
    ) +
    
    geom_line(
      linewidth = 1.5
    ) +
    geom_vline(xintercept = 0,
               linetype = "solid") +
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = "Observed Mortality",
      vjust = -1.0,angle = 90,hjust = 1.75,
      colour = "grey30"
    ) +
    scale_x_continuous(
      breaks = seq(-4, 5, by = 1)
    ) +
    
    facet_wrap(~ forest_type) +
    
    scale_colour_manual(
      values = c(
        "Conifer" = "darkorange",
        "Mixed" = "blue",
        "Deciduous" = "darkgreen"
      )
    ) +
    
    scale_fill_manual(
      values = c(
        "Conifer" = "darkorange",
        "Mixed" = "blue",
        "Deciduous" = "darkgreen"
      )
    ) +
    
    labs(
      x = "Years pre and post observed mortality",
      y = "Detrended NDVI"
    ) +
    
    theme_bw(base_size = 8) +
    
    theme(
      legend.position = "none",
      )
  
  #smooth the empirical curves
  ggplot(window_df,
         aes(years_since, trend,
             colour = forest_type,
             fill = forest_type)) +
    
    geom_smooth(
      method = "loess",
      linewidth = 1.5,
      alpha = 0.25
    ) +
    
    facet_wrap(~ forest_type) +
    
    scale_colour_manual(
      values = c(
        "Conifer" = "darkorange",
        "Mixed" = "blue",
        "Deciduous" = "darkgreen"
      )
    ) +
    
    scale_fill_manual(
      values = c(
        "Conifer" = "darkorange",
        "Mixed" = "blue",
        "Deciduous" = "darkgreen"
      )
    ) +
    
    labs(
      x = "Years since mortality",
      y = "NDVI trend"
    ) +
    
    theme_bw() +
    
    theme(
      legend.position = "none",
      strip.text = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      panel.grid.major = element_line(
        colour = "grey90",
        linewidth = 0.3
      ),
      panel.grid.minor = element_blank()
    )
  
  #NOW BOUNDS BY FOREST TYPE
  #first check number of sites
  n_distinct(window_df$plot_id)
  #then number of sites by forest tyep
  window_df %>%
    group_by(forest_type) %>%
    summarise(n_sites = n_distinct(plot_id))
  #results = 
  #Conifer          51
  #Deciduous        36
  #Mixed            46
  
  #now create forest-type specific plots
  #one panel per forest type
  ggplot(window_df,
         aes(x = years_since,
             y = trend)) +
    
    # Individual site trajectories
    geom_line(
      aes(group = plot_id),
      colour = "grey70",
      alpha = 0.6,
      linewidth = 0.4
    ) +
    
    # 95% CI ribbon
    stat_summary(
      aes(fill = forest_type),
      fun.data = mean_cl_normal,
      geom = "ribbon",
      alpha = 0.25,
      colour = NA
    ) +
    
    # Mean line
    stat_summary(
      aes(colour = forest_type),
      fun = mean,
      geom = "line",
      linewidth = 1.5
    ) +
    
        geom_vline(xintercept = 0,
               linetype = "dotted") +
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = "Mortality event",
      vjust = -1.0,angle = 90,hjust = 1.25,
      colour = "grey30"
    ) +
    scale_x_continuous(
      breaks = seq(-4, 5, by = 1)
    ) +
    
    facet_wrap(~forest_type, scales = "free_y") +
    
    scale_colour_manual(values = c(
      "Mixed" = "blue",
      "Deciduous" = "darkgreen",
      "Conifer" = "darkorange"
    )) +
    
    scale_fill_manual(values = c(
      "Mixed" = "blue",
      "Deciduous" = "darkgreen",
      "Conifer" = "darkorange"
    )) +
    
    labs(
      x = "Years Since Mortality",
      y = "Detrended NDVI"
    ) +
    
    theme_bw() +
    
    theme(
      legend.position = "none",
      strip.text = element_text(size = 12, face = "bold"),
      panel.grid.major = element_line(colour = "grey90"),
      panel.grid.minor = element_blank()
    )
  #separate chart for each forest type and smoothed indiv site lines
  #Firstly conifer
  ggplot(
    filter(window_df, forest_type == "Conifer"),
    aes(x = years_since, y = trend)
  ) +
    
    # Individual site profiles
    geom_smooth(
      aes(group = plot_id),
      method = "loess",
      se = FALSE,
      colour = "grey70",
      linewidth = 0.4
    ) +
    
    # Mean and 95% CI
    geom_smooth(
      colour = "darkorange",
      fill = "darkorange",
      method = "loess",
      se = TRUE,
      linewidth = 1.5,
      alpha = 0.2
    ) +
    
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "black"
    ) +
    geom_vline(xintercept = 0,
               linetype = "dotted") +
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = "Observed Mortality",
      vjust = -1.0,angle = 90,hjust = 1.25,
      colour = "grey30"
    ) +
    scale_x_continuous(
      breaks = seq(-4, 5, by = 1)
    ) +
    
    labs(
      x = "Years pre and post observed mortality",
      y = "Detrended NDVI"
    ) +
    
    theme_bw(base_size = 9) +
    
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank()
    )
  
  #now deciduous
  ggplot(
    filter(window_df, forest_type == "Deciduous"),
    aes(x = years_since, y = trend)
  ) +
    
    # Individual site profiles
    geom_smooth(
      aes(group = plot_id),
      method = "loess",
      se = FALSE,
      colour = "grey70",
      linewidth = 0.4
    ) +
    
    # Mean and 95% CI
    geom_smooth(
      colour = "darkgreen",
      fill = "darkgreen",
      method = "loess",
      se = TRUE,
      linewidth = 1.5,
      alpha = 0.2
    ) +
    
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "black"
    ) +
    
    labs(
      x = "Years Since Mortality",
      y = "Detrended NDVI"
    ) +
    
    theme_bw() +
    
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank()
    )
  #and mixed
  ggplot(
    filter(window_df, forest_type == "Mixed"),
    aes(x = years_since, y = trend)
  ) +
    
    # Individual site profiles
    geom_smooth(
      aes(group = plot_id),
      method = "loess",
      se = FALSE,
      colour = "grey70",
      linewidth = 0.4
    ) +
    
    # Mean and 95% CI
    geom_smooth(
      colour = "blue",
      fill = "blue",
      method = "loess",
      se = TRUE,
      linewidth = 1.5,
      alpha = 0.2
    ) +
    
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "black"
    ) +
    
    labs(
      x = "Years Since Mortality",
      y = "Detrended NDVI"
    ) +
    
    theme_bw() +
    
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank()
    )
  
  #to calculate and visualise =/- 1SD for mixed forest type
  #which will represent the site variability not just mean
  library(dplyr)
  library(ggplot2)
  
  # Calculate mean ± 1 SD by year - mixed forest type
  mixed_summary <- window_df %>%
    filter(forest_type == "Mixed") %>%
    group_by(years_since) %>%
    summarise(
      mean_ndvi = mean(trend, na.rm = TRUE),
      sd_ndvi = sd(trend, na.rm = TRUE),
      lower = mean_ndvi - sd_ndvi,
      upper = mean_ndvi + sd_ndvi,
      .groups = "drop"
    )
 
  # Calculate mean ± 1 SD by year - conifer forest type
  conifer_summary <- window_df %>%
    filter(forest_type == "Conifer") %>%
    group_by(years_since) %>%
    summarise(
      mean_ndvi = mean(trend, na.rm = TRUE),
      sd_ndvi = sd(trend, na.rm = TRUE),
      lower = mean_ndvi - sd_ndvi,
      upper = mean_ndvi + sd_ndvi,
      .groups = "drop"
    )
  
  #and deciduous
  # Calculate mean ± 1 SD by year - deciduous forest type
  decid_summary <- window_df %>%
    filter(forest_type == "Deciduous") %>%
    group_by(years_since) %>%
    summarise(
      mean_ndvi = mean(trend, na.rm = TRUE),
      sd_ndvi = sd(trend, na.rm = TRUE),
      lower = mean_ndvi - sd_ndvi,
      upper = mean_ndvi + sd_ndvi,
      .groups = "drop"
    )
  
  
  
  
  #calculate average SD and coefficient of variation for this plot
  mixed_summary %>%
    summarise(
      mean_sd = mean(sd_ndvi),
      mean_cv = mean(100 * sd_ndvi / mean_ndvi)
    )
  #mean SD = 0.0816
  #mean coefficient of variation = 11.5
  
  #same for conifer plots
  conifer_summary %>%
    summarise(
      mean_sd = mean(sd_ndvi),
      mean_cv = mean(100 * sd_ndvi / mean_ndvi)
    )
  #mean sd = 0.146
  #mean cv = 26.3
  #n = 
  
  #same for deciduous
  decid_summary %>%
    summarise(
      mean_sd = mean(sd_ndvi),
      mean_cv = mean(100 * sd_ndvi / mean_ndvi)
    )
  #mean sd = 0.160
  #mean cv = 28.5
  #n = 
  
  #add legend and n= plus SD and mean CV - mixed forests
  ggplot() +
    
    # Individual site profiles (grey)
    geom_smooth(
      data = filter(window_df, forest_type == "Mixed"),
      aes(
        x = years_since,
        y = trend,
        group = plot_id,
        colour = "Individual plots"
      ),
      method = "loess",
      se = FALSE,
      linewidth = 0.2,
      alpha = 0.2
    ) +
    
    # SD ribbon (blue)
    geom_ribbon(
      data = mixed_summary,
      aes(
        x = years_since,
        ymin = lower,
        ymax = upper,
        fill = "±1 SD" 
      ),
      alpha = 0.2
    ) +
    
    # Mean line (blue)
    geom_line(
      data = mixed_summary,
      aes(
        x = years_since,
        y = mean_ndvi,
        colour = "Mean NDVI" 
      ),
      linewidth = 0.8
    ) +
    
    # Mortality event
    geom_vline(
      xintercept = 0,
      linetype = "solid",
      colour = "grey30",
      linewidth = 0.4
    ) +
    
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = "Observed Mortality",
      vjust = 1.25,hjust = 0.5,
      size = 3,
      colour = "grey30"
    ) +
    scale_x_continuous(
      breaks = seq(-4, 5, by = 1)
    ) +
    
    # Statistics box
    annotate(
      "label",
      x = 4.95,
      y = 0.95,
      hjust = 1,
      vjust = 1,
      size = 2,
      fill = "white",
      label =
        "46 plots\nMean SD = 0.082\nMean CV = 11.5%"
    ) +
    
    # Force legend colours
    scale_colour_manual(
      name = "",
      values = c(
        "Individual plots" = "grey70",
        "Mean NDVI" = "#08519C"
      )
    ) +
    
    scale_fill_manual(
      name = "",
      values = c(
        "±1 SD" = "#08519C"
      )
    ) +
    
    labs(
      x = "Years pre and post observed mortality",
      y = "Detrended NDVI"
    ) +
    
    theme_bw(base_size = 9) +
    theme(
      #axis.title = element_text(size = 8),
      #axis.text = element_text(size = 7),
      #legend.text = element_text(size = 7),
      legend.position = "bottom"
    )
  #rerun removing incomplete years of individual profiles 
  #first filter out sites without full years data
  complete_mixed <- window_df %>%
    filter(forest_type == "Mixed") %>%
    group_by(plot_id) %>%
    filter(
      min(years_since, na.rm = TRUE) <= -4,
      max(years_since, na.rm = TRUE) >= 5
    ) %>%
    ungroup()
  
  #then replot 
  
  p_mixed <- ggplot() +
    
    # Individual site profiles (grey)
    geom_smooth(
      #data = filter(window_df, forest_type == "Mixed"),
      data = complete_mixed,
      aes(
        x = years_since,
        y = trend,
        group = plot_id,
        colour = "Individual plots"
      ),
      method = "loess",
      se = FALSE,
      linewidth = 0.2,
      alpha = 0.2,
      #show.legend = FALSE
    ) +
    
    # SD ribbon (blue)
    geom_ribbon(
      data = mixed_summary,
      aes(
        x = years_since,
        ymin = lower,
        ymax = upper,
        fill = "±1 SD" 
      ),
      alpha = 0.2
    ) +
    
    # Mean line (blue)
    geom_line(
      data = mixed_summary,
      aes(
        x = years_since,
        y = mean_ndvi,
        colour = "Mean NDVI" 
      ),
      linewidth = 0.8
    ) +
    
    # Mortality event
    geom_vline(
      xintercept = 0,
      linetype = "solid",
      colour = "grey30",
      linewidth = 0.4
    ) +
    
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = "Observed Mortality",
      vjust = 1.50,hjust = 0.53,
      size = 3,
      colour = "grey30"
    ) +
    scale_x_continuous(
      breaks = seq(-4, 5, by = 1)
    ) +
    
    # Statistics box
    annotate(
      "label",
      x = 4.95,
      y = 0.95,
      hjust = 1,
      vjust = 1,
      size = 2,
      fill = "white",
      label =
        "46 plots\nMean SD = 0.082\nMean CV = 11.5%"
    ) +
    
    # Force legend colours
    scale_colour_manual(
      name = "",
      values = c(
        "Individual plots" = "grey70",
        "Mean NDVI" = "#08519C"
      )
    ) +
    
    scale_fill_manual(
      name = "",
      values = c(
        "±1 SD" = "#08519C"
      )
    ) +
    
    labs(
      x = "Years pre and post observed mortality",
      y = "Detrended NDVI"
    ) +
    coord_cartesian(
      ylim = c(0.2, 0.9)
    ) +
    theme_bw(base_size = 9) +
    theme(
      legend.position = "bottom"
    )
  
  #lets do the same for conifer forest type
  
  #first filter out sites without full years data
  complete_conifer <- window_df %>%
    filter(forest_type == "Conifer") %>%
    group_by(plot_id) %>%
    filter(
      min(years_since) <= -4,
      max(years_since) >= 5,
      n_distinct(years_since) == 10
    ) %>%
    ungroup()
  
  #then replot 
  
  p_conifer <- ggplot() +
    
    # Individual site profiles (grey)
    geom_smooth(
      data = complete_conifer,
      aes(
        x = years_since,
        y = trend,
        group = plot_id,
        colour = "Individual plots"
      ),
      method = "loess",
      se = FALSE,
      linewidth = 0.2,
      alpha = 0.2,
      #show.legend = FALSE
    ) +
    
    # SD ribbon (orange)
    geom_ribbon(
      data = conifer_summary,
      aes(
        x = years_since,
        ymin = lower,
        ymax = upper,
        fill = "±1 SD" 
      ),
      alpha = 0.2
    ) +
    
    # Mean line (orange)
    geom_line(
      data = conifer_summary,
      aes(
        x = years_since,
        y = mean_ndvi,
        colour = "Mean NDVI" 
      ),
      linewidth = 0.8
    ) +
    
    # Mortality event
    geom_vline(
      xintercept = 0,
      linetype = "solid",
      colour = "grey30",
      linewidth = 0.4
    ) +
    
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = "Observed Mortality",
      vjust = 1.50,hjust = 0.53,
      size = 3,
      colour = "grey30"
    ) +
    scale_x_continuous(
      breaks = seq(-4, 5, by = 1)
    ) +
    
    # Statistics box
    annotate(
      "label",
      x = 4.95,
      y = 0.95,
      hjust = 1,
      vjust = 1,
      size = 2,
      fill = "white",
      label = "n = 39 plots\nMean SD = 0.146\nMean CV = 26.5%"
    ) +
    
    # Force legend colours
    scale_colour_manual(
      name = "",
      values = c(
        "Individual plots" = "grey70",
        "Mean NDVI" = "darkorange"
      )
    ) +
    
    scale_fill_manual(
      name = "",
      values = c(
        "±1 SD" = "darkorange"
      )
    ) +
    
    labs(
      x = "Years pre and post observed mortality",
      y = "Detrended NDVI"
    ) +
    coord_cartesian(
      ylim = c(0.2, 0.9)
    ) +
    theme_bw(base_size = 9) +
    theme(
      legend.position = "bottom"
    )
  
  #rerun removing incomplete years of individual profiles 
  #first filter out sites without full years data
  complete_deciduous <- window_df %>%
    filter(forest_type == "Deciduous") %>%
    group_by(plot_id) %>%
    filter(
      min(years_since) <= -4,
      max(years_since) >= 5,
      n_distinct(years_since) == 10
    ) %>%
    ungroup()
  
  #then replot 
  
  p_deciduous <- ggplot() +
    
    # Individual site profiles (grey)
    geom_smooth(
      data = complete_deciduous,
      aes(
        x = years_since,
        y = trend,
        group = plot_id,
        colour = "Individual plots"
      ),
      method = "loess",
      se = FALSE,
      linewidth = 0.2,
      alpha = 0.2,
      #show.legend = FALSE
    ) +
    
    # SD ribbon (green)
    geom_ribbon(
      data = decid_summary,
      aes(
        x = years_since,
        ymin = lower,
        ymax = upper,
        fill = "±1 SD" 
      ),
      alpha = 0.2
    ) +
    
    # Mean line (green)
    geom_line(
      data = decid_summary,
      aes(
        x = years_since,
        y = mean_ndvi,
        colour = "Mean NDVI" 
      ),
      linewidth = 0.8
    ) +
    
    # Mortality event
    geom_vline(
      xintercept = 0,
      linetype = "solid",
      colour = "grey30",
      linewidth = 0.4
    ) +
    
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = "Observed Mortality",
      vjust = 1.50,hjust = 0.53,
      size = 3,
      colour = "grey30"
    ) +
    scale_x_continuous(
      breaks = seq(-4, 5, by = 1)
    ) +
    
    # Statistics box
    annotate(
      "label",
      x = 4.95,
      y = 0.95,
      hjust = 1,
      vjust = 1,
      size = 2,
      fill = "white",
      label = "n = 20 plots\nMean SD = 0.16\nMean CV = 28.7%"
    ) +
    
    # Force legend colours
    scale_colour_manual(
      name = "",
      values = c(
        "Individual plots" = "grey70",
        "Mean NDVI" = "darkgreen"
      )
    ) +
    
    scale_fill_manual(
      name = "",
      values = c(
        "±1 SD" = "darkgreen"
      )
    ) +
    
    labs(
      x = "Years Since Mortality",
      y = "Detrended NDVI"
    ) +
    coord_cartesian(
      ylim = c(0.2, 0.9)
    ) +
    theme_bw(base_size = 9) +
    theme(
      legend.position = "bottom"
    )
  
   #check final number of plots after filtering - create database
  complete_plots <- window_df %>%
    group_by(plot_id) %>%
    filter(
      n_distinct(years_since) == 10,
      min(years_since) <= -4,
      max(years_since) >= 5
    ) %>%
    ungroup()
  #then obtain sample sizes
  complete_plots %>%
    group_by(forest_type) %>%
    summarise(
      n_plots = n_distinct(plot_id)
    )
  
  
  #count unique plots per forest type
  library(dplyr)
  
  window_df %>%
    group_by(forest_type) %>%
    summarise(
      n_plots = n_distinct(plot_id)
    )
  
  #and finally deciduous forest type
  #first we remove incomplete plots
  complete_deciduous <- window_df %>%
    filter(forest_type == "Deciduous") %>%
    group_by(plot_id) %>%
    filter(
      min(years_since) <= -4,
      max(years_since) >= 5,
      n_distinct(years_since) == 10
    ) %>%
    ungroup()
  
  #then replot 
  
  ggplot() +
    
    # Individual site profiles (grey)
    geom_smooth(
      #data = filter(window_df, forest_type == "Deciduous"),
      data = complete_deciduous,
      aes(
        x = years_since,
        y = trend,
        group = plot_id,
        colour = "Individual plots"
      ),
      method = "loess",
      se = FALSE,
      linewidth = 0.2,
      alpha = 0.2,
      #show.legend = FALSE
    ) +
    
    # SD ribbon (green)
    geom_ribbon(
      data = deciduous_summary,
      aes(
        x = years_since,
        ymin = lower,
        ymax = upper,
        fill = "±1 SD" 
      ),
      alpha = 0.2
    ) +
    
    # Mean line (green)
    geom_line(
      data = deciduous_summary,
      aes(
        x = years_since,
        y = mean_ndvi,
        colour = "Mean NDVI" 
      ),
      linewidth = 0.8
    ) +
    
    # Mortality event
    geom_vline(
      xintercept = 0,
      linetype = "solid",
      colour = "grey30",
      linewidth = 0.4
    ) +
    
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = "Observed Mortality",
      vjust = 1.50,hjust = 0.53,
      size = 3,
      colour = "grey30"
    ) +
    scale_x_continuous(
      breaks = seq(-4, 5, by = 1)
    ) +
    
    # Statistics box
    annotate(
      "label",
      x = 4.95,
      y = 0.95,
      hjust = 1,
      vjust = 1,
      size = 2,
      fill = "white",
      label = "n = 36 plots\nMean SD = 0.160\nMean CV = 28.5%"
    ) +
    
    # Force legend colours
    scale_colour_manual(
      name = "",
      values = c(
        "Individual plots" = "grey70",
        "Mean NDVI" = "darkgreen"
      )
    ) +
    
    scale_fill_manual(
      name = "",
      values = c(
        "±1 SD" = "darkgreen"
      )
    ) +
    
    labs(
      x = "Years pre and post observed mortality",
      y = "Detrended NDVI"
    ) +
    theme_bw(base_size = 9) +
    theme(
      legend.position = "bottom"
    )
  
  #old plotting code
  ggplot() +
    
    # Individual site profiles (grey)
    geom_smooth(
      data = filter(window_df, forest_type == "Deciduous"),
      aes(
        x = years_since,
        y = trend,
        group = plot_id,
        colour = "Individual plots"
      ),
      method = "loess",
      se = FALSE,
      linewidth = 0.4,
      alpha = 0.5
    ) +
    
    # SD ribbon (blue)
    geom_ribbon(
      data = mixed_summary,
      aes(
        x = years_since,
        ymin = lower,
        ymax = upper,
        fill = "±1 SD"
      ),
      alpha = 0.2
    ) +
    
    # Mean line (blue)
    geom_line(
      data = mixed_summary,
      aes(
        x = years_since,
        y = mean_ndvi,
        colour = "Mean NDVI"
      ),
      linewidth = 1.5
    ) +
    
    # Mortality event
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "black"
    ) +
    geom_vline(xintercept = 0,
               linetype = "dotted") +
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = "Mortality event",
      vjust = -1.0,angle = 90,hjust = 1.25,
      colour = "grey30"
    ) +
    scale_x_continuous(
      breaks = seq(-4, 5, by = 1)
    ) +
    
    # Statistics box
    annotate(
      "label",
      x = 4.95,
      y = 0.95,
      hjust = 1,
      vjust = 1,
      size = 4,
      fill = "white",
      label = "n = 36 plots\nMean SD = 0.16\nMean CV = 28.7%"
    ) +
    
    # Force legend colours
    scale_colour_manual(
      name = "",
      values = c(
        "Individual plots" = "grey70",
        "Mean NDVI" = "darkgreen"
      )
    ) +
    
    scale_fill_manual(
      name = "",
      values = c(
        "±1 SD" = "darkgreen"
      )
    ) +
    
    labs(
      x = "Years Since Mortality",
      y = "Detrended NDVI"
    ) +
    
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal"
    )
  
  #count unique plots per forest type
  library(dplyr)
  
  window_df %>%
    group_by(forest_type) %>%
    summarise(
      n_plots = n_distinct(plot_id)
    )
  #forest_type n_plots
 
 # Conifer          51
 # Deciduous        36
 # Mixed            46
  #===================================================
  #create a panel pplot with all three forest type bounds charts
  #create a plot separately
  p_mixed <- ggplot(...) +
    ... +
    labs(title = "Mixed Forest") +
    theme_bw(base_size = 9)
  
  p_conifer <- ggplot(...) +
    ... +
    labs(title = "Conifer Forest") +
    theme_bw(base_size = 9)
  
  p_deciduous <- ggplot(...) +
    ... +
    labs(title = "Deciduous Forest") +
    theme_bw(base_size = 9)
  #use identical y limits for all three
  coord_cartesian(ylim = c(0.2, 0.9))
  #combine with patchwork
  library(patchwork)
  
  final_plot <-
    (p_mixed /
       p_conifer /
       p_deciduous) +
    plot_annotation(
      tag_levels = "A"
    ) &
    theme(
      legend.position = "bottom"
    )
  
  final_plot
  
  #this version removes legends and has just one x axis 
  #at the bottomleaving more space for the three bounds charts
  library(patchwork)
  
  # Add panel titles and remove repeated x-axis titles
  p_mixed <- p_mixed +
    ggtitle("A. Mixed forest") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      axis.title.x = element_blank()
    )
  
  p_conifer <- p_conifer +
    ggtitle("B. Conifer forest") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      axis.title.x = element_blank()
    )
  
  p_deciduous <- p_deciduous +
    ggtitle("C. Deciduous forest") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0)
    )
  
  # Combine panels
  final_plot <-
    (p_mixed /
       p_conifer /
       p_deciduous) +
    plot_layout(
      guides = "collect"
    ) &
    theme(
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  # Add one shared axis title
  final_plot <- final_plot +
    plot_annotation(
      theme = theme(
        plot.margin = margin(10, 10, 10, 10)
      )
    ) &
    labs(
      x = "Years pre and post observed mortality",
      y = "Detrended NDVI"
    )
  
  final_plot
  
  #this version removes the legend entirely
  final_plot <-
    (p_mixed +
       ggtitle("A. Mixed forest") +
       theme(
         axis.title.x = element_blank(),
         legend.position = "none"
       )) /
    
    (p_conifer +
       ggtitle("B. Conifer forest") +
       theme(
         axis.title.x = element_blank(),
         legend.position = "none"
       )) /
    
    (p_deciduous +
       ggtitle("C. Deciduous forest") +
       theme(
         legend.position = "none"
       ))
  annotate(
    "text",
    x = 0,
    y = 0.88,
    size = 2,
    label = "Observed mortality",
  )
  final_plot &
    theme(
      panel.background = element_rect(fill = "white", colour = "black"),
      plot.background = element_rect(fill = "white", colour = "white"),
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      plot.title = element_text(
        face = "bold",
        size = 10,
        hjust = 0
      )
    )
  
  
 #______________________________________________________ 
  
  # SENSITIVITY TEST MORTALITY THRESHOLD
  #STAT TEST TO FIND OPTIMAL OR BOOTSTRAP?
  
  #1. Likelihood-based threshold search
  thresholds <- seq(10, 80, by = 1)
  
  results <- data.frame()
  
  for(t in thresholds){
    
    temp <- window_df %>%
      mutate(
        mort_class = ifelse(
          pct_mort > t,
          "High",
          "Low"
        )
      )
    
    model <- lmer(
      trend ~ years_since +
        I(years_since^2) +
        mort_class +
        years_since:mort_class +
        I(years_since^2):mort_class +
        (1|plot_id),
      data = temp,
      REML = FALSE
    )
    
    results <- rbind(
      results,
      data.frame(
        threshold = t,
        AIC = AIC(model)
      )
    )
  }
  
  
  #then print results - lowest AIC provides best separation of ndvi trajectories
  results[which.min(results$AIC), ]
  
  #threshold      AIC
  #   75        -108651.1
  #i.e. optimal threshold is 75% with AIC = -108651.1
  #plot the AIC profile
  #first create two objecets before plot
  points_df <- results_all %>%
    filter(threshold %% 5 == 0)
  
  label_df <- best_thresholds %>%
    filter(forest_type != "Overall") %>%
    mutate(
      label = paste0("Optimal = ", threshold, "%")
    )
  #then add to ggplot
  library(ggplot2)
  library(dplyr)
  
  points_df <- results %>%
    filter(threshold %% 5 == 0)
  
  label_df <- data.frame(
    threshold = 75,
    label = "Optimal = 75%"
  )
  
  ggplot(results, aes(threshold, AIC)) +
    
    geom_line(
      linewidth = 0.7
    ) +
    
    geom_point(
      data = points_df,
      size = 2
    ) +
    
    geom_vline(
      xintercept = 75,
      linetype = "dashed",
      colour = "red",
      linewidth = 0.8
    ) +
    
    geom_text(
      data = label_df,
      aes(
        x = threshold,
        y = max(results$AIC),
        label = label
      ),
      angle = 90,
      colour = "red",
      size = 3,
      vjust = -0.5,
      hjust = 1.0
    ) +
    
    scale_x_continuous(
      breaks = seq(10, 80, 10)
    ) +
    
    theme_bw(base_size = 8)
  
  #earlier plot 
  library(ggplot2)
  
  ggplot(results, aes(threshold, AIC)) +
    geom_line() +
    geom_point() +
    geom_vline(xintercept = 75, linetype = "dashed", colour = "red") +
    theme_bw(base_size = 8)
  
  #the profile shows The ability of mortality threshold 
  #to discriminate NDVI trajectories improves substantially 
  #above ~50-55% mortality and reaches a broad optimum 
  #around 75-79%.
  
  #1b - rerun separate threshold graph for each forest type
  #reduce points to every 5% and change axis labels to 
  #every 10%
  #run threshold search by forest type
  library(dplyr)
  library(lme4)
  
  thresholds <- seq(10, 80, by = 1)
  
  results_all <- data.frame()
  
  # Overall analysis
  for(t in thresholds){
    
    temp <- analysis_data %>%
      mutate(
        mort_class = ifelse(pct_mort > t, "High", "Low")
      )
    
    model <- lmer(
      trend ~ years_since +
        I(years_since^2) +
        mort_class +
        years_since:mort_class +
        I(years_since^2):mort_class +
        (1 | plot_id),
      data = temp,
      REML = FALSE
    )
    
    results_all <- rbind(
      results_all,
      data.frame(
        threshold = t,
        AIC = AIC(model),
        forest_type = "Overall"
      )
    )
  }
  
  # Forest-type specific analyses
  for(ft in unique(window_df$forest_type)){
    
    data_ft <- window_df %>%
      filter(forest_type == ft)
    
    for(t in thresholds){
      
      temp <- data_ft %>%
        mutate(
          mort_class = ifelse(pct_mort > t, "High", "Low")
        )
      
      model <- lmer(
        trend ~ years_since +
          I(years_since^2) +
          mort_class +
          years_since:mort_class +
          I(years_since^2):mort_class +
          (1 | plot_id),
        data = temp,
        REML = FALSE
      )
      
      results_all <- rbind(
        results_all,
        data.frame(
          threshold = t,
          AIC = AIC(model),
          forest_type = ft
        )
      )
    }
  }
  
  #plot all four lines together
  library(ggplot2)
  
  ggplot(results_all,
         aes(x = threshold,
             y = AIC,
             colour = forest_type)) +
    
    geom_line(linewidth = 1.2) +
    
    scale_colour_manual(
      values = c(
        "Overall" = "black",
        "Mixed" = "blue",
        "Conifer" = "darkorange",
        "Deciduous" = "darkgreen"
      )
    ) +
    
    scale_x_continuous(
      breaks = seq(10, 80, by = 10),      # labels every 10%
      minor_breaks = seq(10, 80, by = 5)  # gridlines every 5%
    ) +
    
    labs(
      x = "Mortality Threshold (%)",
      y = "AIC",
      colour = "Forest Type"
    ) +
    
    theme_bw()
  
  #add points every 5% to make it easier to see results
  ggplot(results_all,
         aes(threshold, AIC, colour = forest_type)) +
    
    geom_line(linewidth = 1.2) +
    
    geom_point(
      data = subset(results_all, threshold %% 5 == 0),
      size = 2
    ) +
    
    scale_colour_manual(
      values = c(
        "Overall" = "black",
        "Mixed" = "blue",
        "Conifer" = "darkorange",
        "Deciduous" = "darkgreen"
      )
    ) +
    
    scale_x_continuous(
      breaks = seq(10, 80, 10),
      minor_breaks = seq(10, 80, 5)
    ) +
    
    theme_bw()
  
  #now we add optimal thresholds to the chart
  #first compute the optimals thresholds
  library(dplyr)
  
  best_thresholds <- results_all %>%
    group_by(forest_type) %>%
    slice_min(AIC, n = 1) %>%
    ungroup()
  
  best_thresholds
  
  #then plot them with dotted line on each graph
  library(ggplot2)
  
  ggplot(results_all,
         aes(x = threshold,
             y = AIC,
             colour = forest_type)) +
    
    geom_line(linewidth = 1.2) +
    
    geom_point(
      data = subset(results_all, threshold %% 5 == 0),
      size = 2
    ) +
    
    geom_vline(
      data = best_thresholds,
      aes(xintercept = threshold,
          colour = forest_type),
      linetype = "dashed",
      linewidth = 0.8,
      show.legend = FALSE
    ) +
    
    scale_colour_manual(
      values = c(
        "Overall"   = "black",
        "Mixed"     = "blue",
        "Conifer"   = "darkorange",
        "Deciduous" = "darkgreen"
      )
    ) +
    
    scale_x_continuous(
      breaks = seq(10, 80, by = 10),
      minor_breaks = seq(10, 80, by = 5)
    ) +
    
    labs(
      x = "Mortality Threshold (%)",
      y = "AIC",
      colour = "Forest Type"
    ) +
    
    theme_bw() +
    
    theme(
      panel.grid.minor.x = element_line(colour = "grey90")
    )
  
  #because the absolute AIC values vary so much between 
  #forest type (due to sample size) this makes visual 
  #comparison difficult. 
  #Instead we will use deltaAIC to rescale each forest type to zero
  #ΔAIC = AIC - min(AIC)
#Calculate AIC
  library(dplyr)
  
  results_all <- results_all %>%
    group_by(forest_type) %>%
    mutate(
      delta_AIC = AIC - min(AIC)
    ) %>%
    ungroup()
  #find optimal thresholds
  best_thresholds <- results_all %>%
    group_by(forest_type) %>%
    slice_min(AIC, n = 1) %>%
    ungroup()
  
  #plot deltaAIC
  library(ggplot2)
  
  ggplot(results_all,
         aes(x = threshold,
             y = delta_AIC,
             colour = forest_type)) +
    
    geom_line(linewidth = 1.2) +
    
    geom_point(
      data = subset(results_all, threshold %% 5 == 0),
      size = 2
    ) +
    
    geom_vline(
      data = best_thresholds,
      aes(xintercept = threshold,
          colour = forest_type),
      linetype = "dotted",
      linewidth = 0.5,
      show.legend = FALSE
    ) +
    #add deltaAIC = 2 referenceline  
    geom_hline(
      yintercept = 2,
      colour = "red",
      linetype = "dotted"
    )
    scale_colour_manual(
      values = c(
        "Overall" = "grey60",
        "Mixed" = "blue",
        "Conifer" = "darkorange",
        "Deciduous" = "darkgreen"
      )
    ) +
    
    scale_x_continuous(
      breaks = seq(10, 80, by = 10),
      minor_breaks = seq(10, 80, by = 5)
    ) +
    
    labs(
      x = "Mortality Threshold (%)",
      y = expression(Delta*AIC),
      colour = "Forest Type"
    ) +
    
    theme_bw() +
    
    theme(
      panel.grid.minor.x = element_line(colour = "grey90")
    )
  
  #why is there not deciduous forest type in the chart?
    table(analysis_data$forest_type)
    unique(analysis_data$forest_type)
    unique(results_all$forest_type)
  #these diagnostics tell us there are deciduous sites in
  #source data but not in results.
 #lets check they were processed:
    for(ft in unique(analysis_data$forest_type)){
    print(ft)
    }  
  #all three forest types are included
  #lets not check for model failure
    for(ft in unique(analysis_data$forest_type)){
      
      cat("\nProcessing:", ft, "\n")
      
      data_ft <- analysis_data %>%
        filter(forest_type == ft)
      
      for(t in thresholds){
        
        temp <- data_ft %>%
          mutate(
            mort_class = ifelse(pct_mort > t, "High", "Low")
          )
        
        model <- try(
          lmer(
            trend ~ years_since +
              I(years_since^2) +
              mort_class +
              years_since:mort_class +
              I(years_since^2):mort_class +
              (1|plot_id),
            data = temp,
            REML = FALSE
          ),
          silent = TRUE
        )
        
        if(inherits(model, "try-error")){
          cat("Failed:", ft, "threshold =", t, "\n")
          next
        }
        
        results_all <- rbind(
          results_all,
          data.frame(
            threshold = t,
            AIC = AIC(model),
            forest_type = ft
          )
        )
      }
    }
   #check if forest type name is ebing altered
    table(results_all$forest_type)
    #deciduous is not in the table
    #lets check for missing values
    analysis_data %>%
      filter(forest_type == "Deciduous") %>%
      summarise(
        n = n(),
        missing_trend = sum(is.na(trend)),
        missing_years = sum(is.na(years_since)),
        missing_mort = sum(is.na(pct_mort)),
        missing_plot = sum(is.na(plot_id))
      )
    #this confrims there is plenty of deciduous data.
    #lets try fitting a single deciduous model manually:
    data_dec <- analysis_data %>%
      filter(forest_type == "Deciduous") %>%
      mutate(
        mort_class = ifelse(pct_mort > 50, "High", "Low")
      )
    
    model_test <- lmer(
      trend ~ years_since +
        I(years_since^2) +
        mort_class +
        years_since:mort_class +
        I(years_since^2):mort_class +
        (1|plot_id),
      data = data_dec,
      REML = FALSE
    )
    
    summary(model_test)
    
   #this tells us the deciduous model fits well and returns:
    #5302 obs
    #34 plots
    #AIC = -25,232
    
    #so the issue is somewhere in the code that builds results_all
    #quick diagnostic
    for(ft in unique(analysis_data$forest_type)){
      
      cat("\nProcessing:", ft, "\n")
      
      data_ft <- analysis_data %>%
        filter(forest_type == ft)
      
      cat("Rows:", nrow(data_ft), "\n")
      
    }
    #how many results for each forest type?
    results_all %>%
      dplyr::count(forest_type)  
   #try running just deciduous thresholds
    results_dec <- data.frame()
    
    for(t in thresholds){
      
      temp <- analysis_data %>%
        filter(forest_type == "Deciduous") %>%
        mutate(
          mort_class = ifelse(pct_mort > t, "High", "Low")
        )
      
      model <- lmer(
        trend ~ years_since +
          I(years_since^2) +
          mort_class +
          years_since:mort_class +
          I(years_since^2):mort_class +
          (1 | plot_id),
        data = temp,
        REML = FALSE
      )
      
      results_dec <- rbind(
        results_dec,
        data.frame(
          threshold = t,
          AIC = AIC(model)
        )
      )
    }
    
    nrow(results_dec)
    
    count(results_all, forest_type)
    
    #lets rewrite the loop and use a list, its safer and faster
    results_list <- list()
    
    counter <- 1
    
    for(ft in c("Mixed","Conifer","Deciduous")) {
      
      data_ft <- window_df %>%
        filter(forest_type == ft)
      
      for(t in thresholds){
        
        temp <- data_ft %>%
          mutate(
            mort_class = ifelse(
              pct_mort > t,
              "High",
              "Low"
            )
          )
        
        model <- try(
          lmer(
            trend ~ years_since +
              I(years_since^2) +
              mort_class +
              years_since:mort_class +
              I(years_since^2):mort_class +
              (1|plot_id),
            data = temp,
            REML = FALSE
          ),
          silent = TRUE
        )
        
        if(inherits(model,"try-error")) {
          
          cat("Failed:", ft, "threshold =", t, "\n")
          next
          
        }
        
        results_list[[counter]] <- data.frame(
          threshold = t,
          AIC = AIC(model),
          forest_type = ft
        )
        
        counter <- counter + 1
      }
    }
    
    results_all <- bind_rows(results_list)
    
    #final plot code to deal with missing deciduous and 
    #make the chart more journal ready e.g. page width 4" etc
    library(dplyr)
    library(ggplot2)
    
    # Calculate Delta AIC
    results_all <- results_all %>%
      group_by(forest_type) %>%
      mutate(
        delta_AIC = AIC - min(AIC)
      ) %>%
      ungroup()
    
    # One optimum threshold per forest type
    best_thresholds <- results_all %>%
      group_by(forest_type) %>%
      slice_min(
        order_by = AIC,
        n = 1,
        with_ties = FALSE
      ) %>%
      ungroup()
    
    # Plot
    p <- ggplot(
      results_all,
      aes(
        x = threshold,
        y = delta_AIC,
        colour = forest_type
      )
    ) +
      
      geom_line(
        linewidth = 1.0
      ) +
      
      geom_vline(
        data = best_thresholds,
        aes(
          xintercept = threshold,
          colour = forest_type
        ),
        linetype = "dashed",
        linewidth = 1.0,
        show.legend = FALSE
      ) +
      
      geom_hline(
        yintercept = 2,
        linetype = "solid",
        colour = "grey40",
        linewidth = 0.4
      ) +
      
      scale_colour_manual(
        values = c(
          "Overall"   = "grey60",
          "Conifer"   = "darkorange",
          "Mixed"     = "blue",
          "Deciduous" = "darkgreen"
        )
      ) +
      
      scale_x_continuous(
        breaks = seq(10, 80, by = 10),
        minor_breaks = seq(10, 80, by = 5)
      ) +
      
      labs(
        x = "Mortality Threshold",
        y = expression(Delta*AIC),
        colour = "Forest type"
      ) +
      
      theme_bw(base_size = 8) +
      
      theme(
        legend.position = "bottom",
        
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        
        panel.grid.major = element_line(
          colour = "grey85",
          linewidth = 0.3
        ),
        
        panel.grid.minor = element_line(
          colour = "grey93",
          linewidth = 0.2
        ),
        
        panel.border = element_rect(
          colour = "black",
          linewidth = 0.8
        )
      )
    
    p
    
  #the overall results line has disappeared - check
    table(results_all$forest_type)
    unique(results_all$forest_type)
  #we need to append them back into the final dataframe 
  #before calculating deltaAIC
  #recreate overall results
    thresholds <- seq(10, 80, by = 1)
    
    overall_results <- data.frame()
    
    for(t in thresholds){
      
      temp <- window_df %>%
        mutate(
          mort_class = ifelse(
            pct_mort > t,
            "High",
            "Low"
          )
        )
      
      model <- lmer(
        trend ~ years_since +
          I(years_since^2) +
          mort_class +
          years_since:mort_class +
          I(years_since^2):mort_class +
          (1 | plot_id),
        data = temp,
        REML = FALSE
      )
      
      overall_results <- rbind(
        overall_results,
        data.frame(
          threshold = t,
          AIC = AIC(model),
          forest_type = "Overall"
        )
      )
    }   
  #combine with forest type results
    results_all <- bind_rows(
      overall_results,
      results_all
    )
  #sense check
    table(results_all$forest_type)
  #all forest types have results
  #now we calculate deltaAIC
    results_all <- results_all %>%
      group_by(forest_type) %>%
      mutate(
        delta_AIC = AIC - min(AIC)
      ) %>%
      ungroup()
  #find one optimal threshold per forest type
    best_thresholds <- results_all %>%
      group_by(forest_type) %>%
      slice_min(
        AIC,
        n = 1,
        with_ties = FALSE
      ) %>%
      ungroup()
  
    #we want to rerun the chart stopping at 75%
    #and stopping each forest type line at its optimum threshold
    #join optimum thresholds onto the result - c
    library(dplyr)
    
    results_plot <- results_all %>%
      left_join(
        best_thresholds %>%
          select(forest_type, best_thresh = threshold),
        by = "forest_type"
      ) %>%
      filter(
        threshold <= pmin(best_thresh + 10, 75)
      )
    #add optimum points  
    best_points <- results_all %>%
      inner_join(
        best_thresholds %>%
          select(forest_type, threshold),
        by = c("forest_type", "threshold")
      )
    
    #now plot the new chart
    
  #journal ready version
    library(ggplot2)
    
    jae_cols <- c(
      "Overall"   = "grey60",
      "Conifer"   = "#D55E00",
      "Mixed"     = "#0072B2",
      "Deciduous" = "#009E73"
    )
    
    p <- ggplot(
      results_plot,
      aes(
        x = threshold,
        y = delta_AIC,
        colour = forest_type
      )
    ) +
      
      geom_line(
        linewidth = 0.7,
        lineend = "round"
      ) +
      #add a layer for optimum points
      geom_point(
        data = best_points,
        aes(
          x = threshold,
          y = delta_AIC,
          colour = forest_type
        ),
        size = 2.5,
        show.legend = FALSE
      ) +
      
      geom_vline(
        data = best_thresholds %>%
          filter(forest_type != "Overall"),
        aes(
          xintercept = threshold,
          colour = forest_type
        ),
        linetype = "dashed",
        linewidth = 0.8,
        show.legend = FALSE
      ) +
      #add a baseline
      geom_hline(
        yintercept = 0,
        colour = "grey60",
        linewidth = 0.4
      ) +
      
      scale_colour_manual(
        values = jae_cols
      ) +
      
      scale_x_continuous(
        limits = c(10, 75),
        breaks = seq(10, 70, 10),
        minor_breaks = NULL
      ) +
      
      labs(
        x = "Mortality threshold (%)",
        y = expression(Delta*AIC),
        colour = "Forest type"
      ) +
      
      theme_bw(base_size = 8) +
      
      theme(
        legend.position = "top",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        panel.grid.major = element_line(
          colour = "grey88",
          linewidth = 0.25
        ),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(
          linewidth = 0.6,
          colour = "black"
        )
      )
    p
    
    
    #check
    table(results_all$forest_type)
    #this confirms overall is on the results
    #Conifer Deciduous     Mixed   Overall 
    #  70        48        71        71 
    
    #add this to existing code
    geom_vline(
      data = best_thresholds,
      aes(
        xintercept = threshold,
        colour = forest_type
      ),
      linetype = "dashed",
      linewidth = 0.8,
      show.legend = FALSE
    )
    
    
  #and save
    ggsave(
      "C:/Users/s422478/OneDrive - Cranfield University/Data/Experimental design/New Data/NDVI_EVI/Figures/threshold_AIC_forest_type.png",
      p,
      width = 140,
      height = 85,
      units = "mm",
      device = cairo_pdf
    )
  
    #add a supplementary figure showing sample size vs threshold
    sample_sizes <- data.frame(
      threshold = c(40, 50, 60, 70, 75),
      sites = c(28, 21, 13, 8, 5)
    )
    
    ggplot(sample_sizes,
           aes(threshold, sites)) +
      
      geom_line(linewidth = 0.7) +
      
      geom_point(size = 2) +
      
      geom_hline(
        yintercept = 14,
        linetype = "dashed",
        colour = "red"
      ) +
      
      scale_x_continuous(
        breaks = seq(40, 80, 5)
      ) +
      scale_y_continuous(
        breaks = seq(0, 30, 5)
      ) +
      labs(
        x = "Mortality threshold (%)",
        y = "High-mortality sites"
      ) +
      
      theme_bw(base_size = 8)
    
    #we need the same by forest type
    
    #we need to recalculate by counting unique plot ids
    library(dplyr)
    
    thresholds <- seq(10, 80, by = 5)
    
    threshold_counts <- lapply(thresholds, function(th){
      
      df %>%
        
        filter(pct_mort > th) %>%
        
        distinct(plot_id, forest_type) %>%
        
        count(forest_type, name = "high_mortality_sites") %>%
        
        mutate(threshold = th)
      
    }) %>%
      
      bind_rows()
    
    #calculate total sites per forest type
    forest_totals <- df %>%
      
      distinct(plot_id, forest_type) %>%
      
      count(forest_type, name = "total_sites")
    
    forest_totals
    
    #forest_type total_sites
    #<chr>             <int>
   # 1 Conifer              60
   # 2 Deciduous            36
   # 3 Mixed                54
   
    #calculate percentages
    threshold_counts <- threshold_counts %>%
      
      left_join(
        forest_totals,
        by = "forest_type"
      ) %>%
      
      mutate(
        pct_sites = 100 * high_mortality_sites / total_sites,
        sufficient = pct_sites >= 10
      )
    #check 10% rule
    threshold_counts %>%
      
      filter(!sufficient)
    #and
    threshold_counts %>%
      
      filter(sufficient) %>%
      
      group_by(forest_type) %>%
      
      summarise(
        max_threshold = max(threshold)
      )
    
    #generate supplementary figure
    library(ggplot2)
    
    p1 <- ggplot(
      threshold_counts,
      aes(
        x = threshold,
        y = pct_sites,
        colour = forest_type
      )
    ) +
      
      geom_line(linewidth = 0.7) +
      
      geom_point(size = 1) +
      
      geom_hline(
        yintercept = 10,
        linetype = "dashed",
        colour = "red"
      ) +
      
      annotate(
        "text",
        x = 78,
        y = 14,
        label = "10% rule",
        colour = "red",
        size = 2.5
      ) +
      
      scale_x_continuous(
        breaks = seq(10, 80, 10)
      ) +
      
      scale_y_continuous(
        limits = c(0, 100)
      ) +
      
      scale_colour_manual(
        values = c(
          "Conifer" = "#D55E00",
          "Deciduous" = "#009E73",
          "Mixed" = "#0072B2"
        )
      ) +
      
      labs(
        x = "Mortality threshold (%)",
        y = "Sites above threshold (%)",
        colour = "Forest type"
      ) +
      
      theme_bw(base_size = 8) +
      
      theme(
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
      )
    p1
    #and save
    ggsave(
      "C:/Users/s422478/OneDrive - Cranfield University/Data/Experimental design/New Data/NDVI_EVI/Figures/sample_size_forest type.png",
      p1,
      width = 140,
      height = 85,
      units = "mm",
      #device = cairo_png
    )
    
  # caption wording: Figure S1. Sample support across 
    #mortality thresholds for each forest type. 
    #Values show the percentage of sites exceeding 
    #each mortality threshold. The dashed red line 
    #indicates the minimum support criterion of 10% of 
    #available sites within a forest type. 
    #Thresholds below this line were considered poorly 
    #represented and interpreted cautiously.
     
  #and create a journal-ready plot
    library(ggplot2)
    
    p <- ggplot(
      results_all,
      aes(
        x = threshold,
        y = delta_AIC,
        colour = forest_type
      )
    ) +
      
      geom_line(
        linewidth = 1.2
      ) +
      
      geom_vline(
        data = best_thresholds,
        aes(
          xintercept = threshold,
          colour = forest_type
        ),
        linetype = "dashed",
        linewidth = 1.0,
        show.legend = FALSE
      ) +
      
      geom_hline(
        yintercept = 2,
        colour = "grey40",
        linetype = "dotted",
        linewidth = 0.7
      ) +
      
      scale_colour_manual(
        values = c(
          "Overall"   = "black",
          "Conifer"   = "darkorange",
          "Mixed"     = "blue",
          "Deciduous" = "darkgreen"
        )
      ) +
      
      scale_x_continuous(
        breaks = seq(10, 80, 10),
        minor_breaks = seq(10, 80, 5)
      ) +
      
      labs(
        x = "Mortality Threshold",
        y = expression(Delta*AIC),
        colour = "Forest type"
      ) +
      
      theme_bw(base_size = 8) +
      
      theme(
        legend.position = "bottom",
        
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        
        panel.grid.major = element_line(
          colour = "grey85",
          linewidth = 0.3
        ),
        
        panel.grid.minor = element_line(
          colour = "grey93",
          linewidth = 0.2
        ),
        
        panel.border = element_rect(
          colour = "black",
          linewidth = 0.8
        )
      )
    
    p  
      
  #deltaAIC computation looks wrong we need to check it
  results_all <- results_all %>%
    group_by(forest_type) %>%
    mutate(delta_AIC = AIC - min(AIC)) %>%
    ungroup()
  #verify - it is actually all correct
  results_all %>%
    group_by(forest_type) %>%
    summarise(
      min_delta = min(delta_AIC),
      max_delta = max(delta_AIC)
    )
  #results forest_type min_delta max_delta
  
  #1 Conifer             0      501.
  #2 Mixed               0     1205.
  #3 Overall             0     1155.
  
  #2. Bootstrapping will tell us if the apparent optimum is stable
  #resample plots, repeat threshold search, record threshold each time,
  #obtain distribution of optimal thresholds
  #we will bootstrap 121 plot ids (independent units)
  #first create a function 'find_optimal_threshold'
  library(lme4)
  library(dplyr)
  
  find_optimal_threshold <- function(dat,
                                     thresholds = seq(10, 79, by = 1)) {
    
    aic_values <- rep(NA, length(thresholds))
    
    for(i in seq_along(thresholds)) {
      
      th <- thresholds[i]
      
      temp <- dat %>%
        mutate(
          mort_class = ifelse(
            pct_mort > th,
            "High",
            "Low"
          )
        )
      
      # avoid thresholds that create tiny groups
      
      counts <- table(temp$mort_class)
      
      if(any(counts < 10))
        next
      
      fit <- lmer(
        trend ~ years_since +
          I(years_since^2) +
          mort_class +
          years_since:mort_class +
          I(years_since^2):mort_class +
          (1|plot_id),
        data = temp,
        REML = FALSE
      )
      
      aic_values[i] <- AIC(fit)
    }
    
    thresholds[which.min(aic_values)]
  }
  #now the bootstrap code
  set.seed(123)
  
  plot_ids <- unique(analysis_data$plot_id)
  
  B <- 500
  
  boot_thresh <- numeric(B)
  
  for(b in 1:B){
    
    sampled_ids <- sample(
      plot_ids,
      size = length(plot_ids),
      replace = TRUE
    )
    
    boot_data <- bind_rows(
      lapply(sampled_ids,
             function(x)
               analysis_data %>%
               filter(plot_id == x))
    )
    
    boot_thresh[b] <- find_optimal_threshold(boot_data)
    
    if(b %% 20 == 0)
      print(b)
  }
  
  #summarise the bootstrap results
  median(boot_thresh)
  
  quantile(
    boot_thresh,
    probs = c(0.025, 0.25, 0.75, 0.975)
  )
  #and draw a plot
  boot_df <- data.frame(
    threshold = boot_thresh
  )
  
  ggplot(boot_df,
         aes(threshold)) +
    geom_histogram(
      bins = 15,
      fill = "darkgreen",
      colour = "black"
    ) +
    theme_classic() +
    labs(
      x = "Optimal mortality threshold (%)",
      y = "Frequency"
    )
  #check bootstrap results:
  summary(boot_thresh)
  
  table(is.na(boot_thresh))
  #check how many real results we have
  table(boot_thresh)
  #this shows lots of 0
  #remove incomplete runs (from pressing esc during operation)
  boot_thresh_valid <- boot_thresh[boot_thresh > 0]
  
  summary(boot_thresh_valid)
  table(boot_thresh_valid)
  #plot
  ggplot(
    data.frame(threshold = boot_thresh_valid),
    aes(threshold)
  ) +
    geom_histogram(
      bins = 10,
      fill = "darkgreen",
      colour = "black"
    ) +
    theme_classic()
  #there are too few iterations (shoudl be around 500)
  #check how many finished = only 17
  length(boot_thresh[boot_thresh > 0])
  #rerun with this code:
  B <- 50
  #initialise with this code
  boot_thresh <- rep(NA_real_, B)
  #summarise
  boot_thresh_valid <- boot_thresh[!is.na(boot_thresh)]
  
  median(boot_thresh_valid)
  
  quantile(
    boot_thresh_valid,
    probs = c(0.025, 0.5, 0.975)
  )
  #now we calculate thresholds within the AIC plateau i.e.deltaAIC >2
  best_aic <- min(results$AIC)
  
  subset(
    results,
    AIC <= best_aic + 2
  )
  #this results in best thresholds of 75,76,77,78,79
  #so we can conclude optimal threshold is 75-79%
  #for completeness we should check number of sites in each group
  table(
    ifelse(analysis_data$pct_mort > 75,
           "High",
           "Low")
  )
  #number of plots:
  analysis_data %>%
    group_by(plot_id) %>%
    summarise(
      pct_mort = first(pct_mort)
    ) %>%
    mutate(
      mort75 = ifelse(pct_mort > 75, "High", "Low")
    ) %>%
    count(mort75)
  #this reveals only 5 plots in high mortality group (>75%)
  
  #3.an alternative method is to determine sites numbers at different thresholds
  plot_summary <- analysis_data |>
    group_by(plot_id) |>
    summarise(pct_mort = first(pct_mort))
  
  table(plot_summary$pct_mort > 40)
  table(plot_summary$pct_mort > 50)
  table(plot_summary$pct_mort > 60)
  table(plot_summary$pct_mort > 70)
  table(plot_summary$pct_mort > 75)
  #60% could be optimum at 13 sites (v 5 for 75%)
  #fit mortality as a continuous variable
  continuous_model <- lmer(
    trend ~ years_since * pct_mort +
      I(years_since^2) * pct_mort +
      (1 | plot_id),
    data = window_df,
    REML = FALSE
  )
  #examine the results
  summary(continuous_model)
  #compare against threshold model 
  anova(mortality_model, continuous_model)
  AIC(continuous_model, mortality_model)
  #AIC is lower for continuous model suggesting mortality
  #acts as a gradient rather than a sharp threshold
  
  #now lets visualise the model
  #first establish we have all years data in the analysis dataset
  range(window_df$years_since)
  range(newdat$years_since)
  #specify full -4 to 5 yrs data model
  quad_model <- lmer(
    trend ~ years_since + I(years_since^2) + (1|plot_id),
    data = dat,
    REML = FALSE
  )
  #then rerun predictions
  newdat <- expand.grid(
    years_since = seq(-4, 5, length.out = 100),
    pct_mort = c(20, 40, 60, 80),
    plot_id = NA
  )
  #1. Predicted trajectories for several mortality levels
  #newdat <- expand.grid(
   # years_since = seq(
   #   min(analysis_data$years_since),
   #   max(analysis_data$years_since),
   #   length.out = 100
  #  ),
 #   pct_mort = c(20, 40, 60, 80),
 #   plot_id = NA
 # )
 #plot predicted NDVI trend by mortality severity
  newdat$pred <- predict(
    continuous_model,
    newdata = newdat,
    re.form = NA
  )
  #plot
  library(ggplot2)
  
  ggplot(newdat,
         aes(years_since,
             pred,
             colour = factor(pct_mort))) +
    geom_line(linewidth = 1.0) +
    geom_vline(xintercept = 0,
               linetype = "dotted") +
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = "Observed Mortality", size = 2.0,
      vjust = -0.5,hjust = 1.61,angle = 90,
      colour = "grey30"
    ) +
    scale_x_continuous(
      breaks = seq(-4, 5, by = 1)
    ) +
    
    theme_bw(base_size = 8) +
    labs(
      x = "Years Since Mortality",
      y = "Predicted NDVI Trend",
      colour = "% Mortality"
    )
 
  #run panel plot with predicted trends by forest type
  #generate the data
  formula(continuous_model)
  ls()
  formula(pred_forest)
  formula(quad_model)
  formula(forest_model)
  formula(mortality_model60)
  
  #by faceting the existing plot
  library(ggplot2)
  
  ggplot(
    newdat,
    aes(
      x = x,
      y = predicted,
      colour = group,
      fill = group
    )
  ) +
    
    geom_line(linewidth = 0.9) +
    
    geom_ribbon(
      aes(
        ymin = conf.low,
        ymax = conf.high
      ),
      alpha = 0.2,
      colour = NA
    ) +
    
    facet_wrap(
      ~ facet,
      ncol = 2
    ) +
    
    geom_vline(
      xintercept = 0,
      linetype = "dotted",
      colour = "black",
      linewidth = 0.7
    ) +
    
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = "Mortality event",
      angle = 90,
      vjust = 1.5,
      size = 2.5,
      colour = "grey30"
    ) +
    
    labs(
      x = "Years Since Mortality",
      y = "Predicted NDVI Trend",
      colour = "Mortality Group",
      fill = "Mortality Group"
    ) +
    
    scale_x_continuous(
      breaks = seq(-4, 5, 1)
    ) +
    
    theme_classic(base_size = 8) +
    
    theme(
      legend.position = "bottom",
      
      strip.background = element_rect(
        fill = "white",
        colour = "black"
      ),
      
      strip.text = element_text(
        face = "bold"
      ),
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
 
  attr(pred_combined, "terms")
  formula(forest_model)
  formula(combined_model)
  
   #alternative version
  ggplot(pred,
         aes(x = x,
             y = predicted,
             colour = factor(mortality))) +
    
    geom_line(linewidth = 1) +
    
    facet_wrap(
      ~forest_type,
      ncol = 2,
      scales = "fixed"
    ) +
    
    geom_vline(
      xintercept = 0,
      linetype = "dotted",
      colour = "black"
    ) +
    
    scale_colour_brewer(
      palette = "Set1",
      name = "% Mortality"
    ) +
    
    labs(
      x = "Years Since Mortality",
      y = "Predicted NDVI Trend"
    ) +
    
    theme_bw(base_size = 10) +
    
    theme(
      legend.position = "bottom",
      strip.background = element_rect(
        fill = "white",
        colour = "black"
      ),
      strip.text = element_text(
        face = "bold"
      )
    )
  
  #2.raw data and fitted surface (useful exploratory figure)
  ggplot(newdat,
         aes(years_since,
             pred,
             colour = pct_mort)) +
    geom_point(alpha = 0.3) +
    scale_colour_viridis_c() +
    theme_bw()
  
  #3.faceted mortality bands - reveals where trajectories separate
  newdat$mort_band <- cut(
    newdat$pct_mort,
    breaks = c(0, 20, 40, 60, 80, 100),
    include.lowest = TRUE
  )
  
  ggplot(
    newdat,
    aes(
      years_since,
      pred,
      colour = mort_band
    )
  ) +
    geom_line(linewidth = 1.2)
      se = TRUE
     +
    geom_vline(
      xintercept = 0,
      linetype = "dotted"
    ) +
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = "Mortality event",
      size = 4,
      vjust = -1.0,
      angle = 90,
      hjust = 1.75,
      colour = "grey30"
    ) +
    scale_x_continuous(
      breaks = seq(-4, 5, by = 1)
    ) +
  theme_bw()

  