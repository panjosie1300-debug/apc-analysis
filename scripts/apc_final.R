## ==============================================================================
## Comparative age–period–cohort analysis reveals divergent temporal trajectories 
## of breast and colorectal cancer among women in Singapore
## ==============================================================================

## -----------------------------------------------------------------------------
## 0. SETUP: install/load packages, register fonts
## -----------------------------------------------------------------------------

# Install and load packages
required_packages <- list("dplyr", "splines", "tidyr", "ggplot2", "Epi", "stringr", 
                          "colorspace", "ggpubr", "showtext", "epitools", "purrr", "viridis", 
                          "rstudioapi")
packages_to_install <- required_packages[!required_packages %in% installed.packages()]
install.packages(packages_to_install)

## NOTE ON WORKING DIRECTORY:
## This script uses relative paths which assumed that the working directory is the 
## repo's scripts/ folder. Open this .R file in Rstudio and run it from top to bottom.
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load other packages
library(dplyr)
library(splines)
library(tidyr)
library(ggplot2)
library(Epi)
library(stringr)
library(colorspace)
library(ggpubr)
library(showtext)
library(epitools)
library(purrr)


## -----------------------------------------------------------------------------
## 1. LOAD AND CLEAN CANCER INCIDENCE DATA
## -----------------------------------------------------------------------------

# Load and clean cancer incidence data
cancer <- read.csv("../data/cancer.csv", stringsAsFactors = T)

##  Separate year column into start year of period and end year of period and add column for midpoint of period
cancer <- cancer %>% 
  rename(period = year) %>%
  separate(period, into = c("period_start", "period_end"), sep = "-", remove = F, convert = T) %>%
  mutate(period_mid = (period_start + period_end)/2)
cancer$period_mid <- as.numeric(cancer$period_mid)

## Filter data down to include only cases for females aged 25 to 84 until the year 2023
cancer_clean <- cancer %>%
  filter(sex == "F", period_end <= 2023) %>%
  ## Clean age column and split into two numeric variables representing lower and upper bounds of age group
  mutate(
    agedxgrp = str_remove_all(agedxgrp, "Years"),
    agedxgrp = str_trim(agedxgrp),
    age_start = as.numeric(str_extract(agedxgrp, "^\\d+"))
  )
cancer_clean <- cancer_clean %>%
  filter(age_start >= 25 & age_start <= 80) %>%
  mutate(
    age_end   = sapply(str_extract_all(agedxgrp, "\\d+"), function(x) as.numeric(x[2]) ),
    age_mid = (age_start + age_end)/2
  )

##  Rename incidence column to something easy to work with
cancer_clean <- cancer_clean %>% 
  rename(n = n..incidence.count.)
cancer_clean$n <- as.character(cancer_clean$n)

##  Create different incidence count columns where <5 incidences are resolved to different integers (1 through 4)
##  Sensitivity analysis can then be done to observe change in trends when suppression values resolved to different inetgers 
##  Note: incidence count has to be whole number because Poisson model expects integer count data
cancer_clean$n1 <- cancer_clean$n
cancer_clean$n2 <- cancer_clean$n
cancer_clean$n3 <- cancer_clean$n
cancer_clean$n4 <- cancer_clean$n
cancer_clean$n1[cancer_clean$n == "<5"] <- as.character(1)
cancer_clean$n2[cancer_clean$n == "<5"] <- as.character(2)
cancer_clean$n3[cancer_clean$n == "<5"] <- as.character(3)
cancer_clean$n4[cancer_clean$n == "<5"] <- as.character(4)
cancer_clean$n1 <- as.numeric(cancer_clean$n1)
cancer_clean$n2 <- as.numeric(cancer_clean$n2)
cancer_clean$n3 <- as.numeric(cancer_clean$n3)
cancer_clean$n4 <- as.numeric(cancer_clean$n4)


## -----------------------------------------------------------------------------
## 2. DESCRIPTIVE STATISTICS (overall and by race, for breast/colorectal)
## -----------------------------------------------------------------------------

## descriptive stats
breast <- subset(cancer_clean, site == "Breast")
colo <- subset(cancer_clean, site == "Colon & rectum")
descriptive <- tibble(
  total_cases = sum(cancer_clean$n3)
)
descriptive$breast_percent = sum(breast$n3)/descriptive$total_cases*100
descriptive$colo_percent = sum(colo$n3)/descriptive$total_cases*100

# Breast cancer case totals by race group (CN = Chinese, IN = Indian, MY = Malay, XX = other/unknown)
CN_total_breast <- breast %>%
  filter(racegrp == "CN") %>%
  summarise(total = sum(n3)) %>%
  pull(total)
IN_total_breast <- breast %>%
  filter(racegrp == "IN") %>%
  summarise(total = sum(n3)) %>%
  pull(total)
MY_total_breast <- breast %>%
  filter(racegrp == "MY") %>%
  summarise(total = sum(n3)) %>%
  pull(total)
XX_total_breast <- breast %>%
  filter(racegrp == "XX") %>%
  summarise(total = sum(n3)) %>%
  pull(total)
breast_byrace <- tibble(
  total_no = sum(breast$n3),
  CN_percent = CN_total_breast/total_no*100,
  IN_percent = IN_total_breast/total_no*100,
  MY_percent = MY_total_breast/total_no*100,
  XX_percent = XX_total_breast/total_no*100,
) 

# Colorectal cancer case totals by race group
CN_total_colo <- colo %>%
  filter(racegrp == "CN") %>%
  summarise(total = sum(n3)) %>%
  pull(total)
IN_total_colo <- colo %>%
  filter(racegrp == "IN") %>%
  summarise(total = sum(n3)) %>%
  pull(total)
MY_total_colo <- colo %>%
  filter(racegrp == "MY") %>%
  summarise(total = sum(n3)) %>%
  pull(total)
XX_total_colo <- colo %>%
  filter(racegrp == "XX") %>%
  summarise(total = sum(n3)) %>%
  pull(total)
colo_byrace <- tibble(
  total_no = sum(colo$n3),
  CN_percent = CN_total_colo/total_no*100,
  IN_percent = IN_total_colo/total_no*100,
  MY_percent = MY_total_colo/total_no*100,
  XX_percent = XX_total_colo/total_no*100
) 
breast_byrace 
colo_byrace


## -----------------------------------------------------------------------------
## 3. RESTRICT TO 3 MAJOR RACE GROUPS AND DERIVE COHORT MIDPOINTS
## -----------------------------------------------------------------------------

## remove XX racegrp, leaving 3 major race groups in Singapore, Chinese (CN), Malay (MY) and Indian
cancer_clean <- cancer_clean %>%
  filter(racegrp != "XX")

## Add rows for overall data across the 3 ethnicities
overall_data <- cancer_clean %>%
  group_by(age_mid, period_mid, age_start, age_end, period_start, period_end, agedxgrp, period, site, sex) %>%
  summarise(n1 = sum(n1),
            n2 = sum(n2),
            n3 = sum(n3),
            n4 = sum(n4),
            .groups = "drop")
overall_data$racegrp <- "All"
cancer_clean <- bind_rows(cancer_clean, overall_data)

## Derive cohort midpoints
###   cohort = period - age, using midpoints since interval widths match (5-yr both)
cancer_clean <- cancer_clean %>%
  mutate(cohort_mid = period_mid - age_mid)
cancer_clean$cohort_mid <- as.numeric(cancer_clean$cohort_mid)


## -----------------------------------------------------------------------------
## 4. LOAD AND CLEAN POPULATION DATA, THEN MERGE WITH INCIDENCE DATA
## -----------------------------------------------------------------------------

# Load and clean population data
pop <- read.csv("../data/population.csv", stringsAsFactors = T)
pop <- pop %>% 
  rename(agedxgrp = agegrp) %>%
  mutate(
    agedxgrp = str_remove_all(agedxgrp, "Years"),
    agedxgrp = str_trim(agedxgrp),
    age_start = as.numeric(str_extract(agedxgrp, "\\d+"))
  )
pop_clean <- pop %>%
  filter(age_start >= 25 & age_start <= 80, sex == "F", year <= 2023, racegrp != "XX")

## Assign populations to year periods in cancer_clean
# (join each single calendar year of population data to the 5-year period it falls within)
period_lookup <- cancer_clean %>%
  distinct(period, period_start, period_end)
pop_clean <- left_join(pop_clean, period_lookup, by = join_by(year >= period_start, year <= period_end))

## Calculate overall population data across the three ethnicities
all_pop <- pop_clean %>%
  group_by(year, sex, agedxgrp, age_start, period, period_start, period_end) %>%
  summarise (population = sum(population),
             .groups = "drop")

## Add overall data to population dataset
all_pop$racegrp <- "All"
pop_sum <- bind_rows(pop_clean, all_pop) %>%
  group_by(agedxgrp, racegrp, period, sex) %>%
  summarise(pop = sum(population), .groups = "drop")

## Join population data to incidence dataset
data <- cancer_clean %>% left_join(pop_sum, by = c("period", "agedxgrp", "racegrp", "sex"))


## -----------------------------------------------------------------------------
## 5. AGE-STANDARDISED INCIDENCE RATES (ASR) BY PERIOD, STRATIFIED BY ETHNICITY
## -----------------------------------------------------------------------------

# Plot age-standardized incidence rates of colorectal cancer stratified by ethnic group from 1968 to 2023 (as indicated by the midyear of the 5-yearly intervals).
# https://seer.cancer.gov/stdpopulations/stdpop.19ages.html

# Helper: fix the display order of race groups (used consistently across all plots below)
set_race_levels <- function(df) df %>% mutate(racegrp = factor(racegrp, levels = c("All", "IN", "MY", "CN")))

# SEGI world standard population weights, by 5-year age band, used for direct age-standardisation
segi_std <- tibble(
  age_start = c(0, seq(5, 90, by = 5)),
  std_pop = c(24000, 96000, 100000, 90000, 90000, 80000, 80000, 60000, 60000, 
              60000, 60000, 50000, 40000, 40000, 30000, 20000, 10000, 5000, 5000)
)

# Compute a directly age-standardised rate with gamma-distribution confidence intervals
compute_asr <- function(count, pop, stdpop, conf.level = 0.95) {
  rate <- count / pop
  weight <- stdpop / 1000000   
  y <- sum(weight * rate)
  v <- sum((weight^2) * (count / pop^2))
  wM <- max(weight / pop)
  alpha <- 1 - conf.level
  lci <- qgamma(alpha/2, shape = y^2/v, rate = y/v)
  uci <- qgamma(1 - alpha/2, 
                shape = (y + wM)^2 / (v + wM^2), 
                rate  = (y + wM)  / (v + wM^2) )
  c(crude.rate = sum(count) / sum(pop), adj.rate = y, lci = lci, uci = uci)
}

# Apply the SEGI standard weights to a single site/racegrp/period group and scale to per-100,000
apply_asr <- function(df, df_call_n) {
  df <- df %>% 
    left_join(segi_std, by = "age_start")
  
  
  result <- compute_asr(
    count  = df_call_n,
    pop    = df$pop,
    stdpop = df$std_pop
  )
  
  tibble(
    total_count = sum(df_call_n, na.rm = TRUE),
    total_pop   = sum(df$pop, na.rm = TRUE),
    value       = result["adj.rate"] * 100000,
    lowercl     = result["lci"] * 100000,
    uppercl     = result["uci"] * 100000
  )
}

# Compute ASR summaries separately for each of the 4 possible imputation values (n1-n4)
asr_summary_1 <- data %>%
  group_by(site, racegrp, period_mid) %>%
  group_split() %>%
  map_dfr(function(df_group) {
    grp_name <- df_group %>% select(site, racegrp, period_mid) %>% slice(1)
    result <- apply_asr(df_group, df_group$n1)
    bind_cols(grp_name, result)
  })
asr_summary_1 <- set_race_levels(asr_summary_1)

asr_summary_2 <- data %>%
  group_by(site, racegrp, period_mid) %>%
  group_split() %>%
  map_dfr(function(df_group) {
    grp_name <- df_group %>% select(site, racegrp, period_mid) %>% slice(1)
    result <- apply_asr(df_group, df_group$n2)
    bind_cols(grp_name, result)
  })
asr_summary_2 <- set_race_levels(asr_summary_2)

asr_summary_3 <- data %>%
  group_by(site, racegrp, period_mid) %>%
  group_split() %>%
  map_dfr(function(df_group) {
    grp_name <- df_group %>% select(site, racegrp, period_mid) %>% slice(1)
    result <- apply_asr(df_group, df_group$n3)
    bind_cols(grp_name, result)
  })
asr_summary_3 <- set_race_levels(asr_summary_3)

asr_summary_4 <- data %>%
  group_by(site, racegrp, period_mid) %>%
  group_split() %>%
  map_dfr(function(df_group) {
    grp_name <- df_group %>% select(site, racegrp, period_mid) %>% slice(1)
    result <- apply_asr(df_group, df_group$n4)
    bind_cols(grp_name, result)
  })
asr_summary_4 <- set_race_levels(asr_summary_4)

# Smooth the ASR trend over time (per site/racegrp) with Loess smoothing
smooth_asr <- function(df, span = 0.75, n_out = 200) {
  df %>%
    group_by(site, racegrp) %>%
    group_modify(~ {
      grid <- seq(min(.x$period_mid), max(.x$period_mid), length.out = n_out)
      
      fit_val <- loess(value   ~ period_mid, data = .x, weights = total_count, span = span)
      fit_lo  <- loess(lowercl ~ period_mid, data = .x, weights = total_count, span = span)
      fit_hi  <- loess(uppercl ~ period_mid, data = .x, weights = total_count, span = span)
      
      tibble(
        period_mid = grid,
        value      = predict(fit_val, newdata = tibble(period_mid = grid)),
        lowercl    = predict(fit_lo,  newdata = tibble(period_mid = grid)),
        uppercl    = predict(fit_hi,  newdata = tibble(period_mid = grid))
      )
    }) %>%
    ungroup() %>%
    mutate(lowercl = pmax(0, lowercl))
}

asr_smoothed_1 <- smooth_asr(asr_summary_1)
asr_smoothed_2 <- smooth_asr(asr_summary_2)
asr_smoothed_3 <- smooth_asr(asr_summary_3)
asr_smoothed_4 <- smooth_asr(asr_summary_4)

# Shared colour/linetype scheme for race groups, used across all plots below
mycolours <- c("All" = "#002D62", 
               "IN" = "#785EF0", 
               "MY" = "#DC267F",
               "CN" = "#FE6100")
mycolours_light <- lighten(mycolours, amount = 0.6)  
mylines <- c("All" = "solid",     
             "IN"  = "42",
             "MY"  = "4212",   
             "CN"  = "11") 

# Plot smoothed breast cancer ASR trend over calendar year, by race
breast_stand_smooth <- function(df) {
  df %>%
    filter(site == "Breast") %>%
    ggplot(aes(x = period_mid, y = value, colour = racegrp, linetype = racegrp)) +
    geom_ribbon(aes(ymin = lowercl, ymax = uppercl, fill = racegrp),
                colour = NA, alpha = 0.25) +
    geom_line(linewidth = 1) +
    scale_colour_manual(values = mycolours) +
    scale_linetype_manual(values = mylines) +
    scale_fill_manual(values = mycolours_light) +
    ylim(0,110) +
    labs(x = "Year", 
         y = "Breast cancer ASR",
         colour = "Race",
         linetype = "Race",
         fill = "Race") +
    theme_classic() +
    theme(text = element_text(face = "bold"),
          axis.title.y = element_text(margin = margin(r=5), size = 14),
          axis.title.x = element_text(margin = margin(t=10), size = 14),
          axis.text = element_text(size = 12, colour = "black"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.line = element_line(linewidth = 1),
          axis.ticks = element_line(linewidth = 1),
          legend.text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5),
          plot.margin = margin(t = 25))
}

# Plot smoothed colorectal cancer ASR trend over calendar year, by race
colo_stand_smooth <- function(df) {
  df %>%
    filter(site == "Colon & rectum") %>%
    ggplot(aes(x = period_mid, y = value, colour = racegrp, linetype = racegrp)) +
    geom_ribbon(aes(ymin = lowercl, ymax = uppercl, fill = racegrp),
                colour = NA, alpha = 0.25) +
    geom_line(linewidth =1) +
    scale_colour_manual(values = mycolours) +
    scale_linetype_manual(values = mylines) +
    scale_fill_manual(values = mycolours_light) +
    ylim(0,110) +
    labs(x = "Year", 
         y = "Colorectal cancer ASR",
         colour = "Race",
         linetype = "Race", 
         fill = "Race") +
    theme_classic() +
    theme(text = element_text(face = "bold"),
          axis.title.y = element_text(margin = margin(r=5), size = 14),
          axis.title.x = element_text(margin = margin(t=10), size = 14),
          axis.text = element_text(size = 12, colour = "black"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.line = element_line(linewidth = 1),
          axis.ticks = element_line(linewidth = 1),
          legend.text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5),
          plot.margin = margin(t = 25))
}

# Generate each combination of site x suppression-resolution scheme (1-4)
breast_stand_smooth_1 <- breast_stand_smooth(asr_smoothed_1)
colo_stand_smooth_1 <- colo_stand_smooth(asr_smoothed_1)
breast_stand_smooth_2 <- breast_stand_smooth(asr_smoothed_2)
colo_stand_smooth_2 <- colo_stand_smooth(asr_smoothed_2)
breast_stand_smooth_3 <- breast_stand_smooth(asr_smoothed_3)
colo_stand_smooth_3 <- colo_stand_smooth(asr_smoothed_3)
breast_stand_smooth_4 <- breast_stand_smooth(asr_smoothed_4)
colo_stand_smooth_4 <- colo_stand_smooth(asr_smoothed_4)

## Plot ASR graphs for sensitivity analysis
# (grid comparing all 4 suppression-resolution schemes side by side, for both cancer sites)

asr_labels <- c( "a) Breast, <5 suppressed to 1", "", "b) Colorectal, <5 suppressed to 1", 
                 "", "", "",
                 "c) Breast, <5 suppressed to 2", "", "d) Colorectal, <5 suppressed to 2",
                 "", "", "",
                 "e) Breast, <5 suppressed to 3", "", "f) Colorectal, <5 suppressed to 3", 
                 "", "", "",
                 "g) Breast, <5 suppressed to 4", "", "h) Colorectal, <5 suppressed to 4")
ggarrange(breast_stand_smooth_1, NULL, colo_stand_smooth_1, 
          NULL, NULL, NULL,
          breast_stand_smooth_2, NULL, colo_stand_smooth_2, 
          NULL, NULL, NULL,
          breast_stand_smooth_3, NULL, colo_stand_smooth_3, 
          NULL, NULL, NULL,
          breast_stand_smooth_4, NULL, colo_stand_smooth_4, 
          labels = asr_labels, 
          font.label = list(size = 14, color = "black", face = "bold"),
          hjust = 0, widths = c(1,0.1,1), heights = c(1, 0.05, 1, 0.05, 1, 0.05, 1),
          ncol = 3, nrow = 7, common.legend = T)

ggsave("../results/diff_suppression_counts_smoothed_asr.pdf", width = 8, height = 11)
dev.off()

## Plot Figure 1. Age-standardised breast (left) and colorectal (right) cancer incidence rates in females 
##                of different races in Singapore aged 25 to 84, 1968 - 2023. Year refers to midyears of 
##                5-year periods; CN = Chinese, IN = Indian, MY = Malay, All = CN+IN+MY. Shading indicates
##                95 % confidence intervals; curves and confidence intervals smoothed with Loess smoothing, span = 0.75.

# Use imputation value 3 (i.e. n3, <5 counts -> 3) as the main analysis
breast_stand <- breast_stand_smooth_3 + ylim(0,90)
colo_stand <- colo_stand_smooth_3 + ylim(0,90)
ggarrange(breast_stand, NULL, colo_stand,
          labels = c("a)", "", "b)"),
          widths = c(1, 0.1, 1),
          font.label = list(size = 14, color = "black", face = "bold"),
          ncol = 3, nrow = 1,
          common.legend = T)
ggsave("../results/ASR_smoothed.pdf", width = 8, height = 4.5)
dev.off()


## -----------------------------------------------------------------------------
## 6. AGE-SPECIFIC RATES BY ETHNICITY (POOLED ACROSS ALL PERIODS) - FIGURE 2
## -----------------------------------------------------------------------------

# Plot age-specific cancer rates stratified by ethnicity and 5-year period i.e. Figure 2
## Figure 2. Age-specific rates of breast (left) and colorectal (right) cancer in females of different 
##           ethnicities aged 25 to 84 in Singapore. Age refers to midyears of 5-year age groups; 
##           CN = Chinese, IN = Indian, MY = Malay, All = CN+IN+MY. Shading indicates 95 % confidence 
##           intervals; curves and confidence intervals smoothed with Loess smoothing, span = 0.75

# Smooth an age-specific rate curve (with CIs) per race group, using a weighted Loess fit
smooth_age <- function(df, span = 0.75, n_out = 200) {
  df %>%
    group_by(racegrp) %>%
    group_modify(~{
      grid <- seq(min(.x$age_mid), max(.x$age_mid), length.out = n_out)
      fit_val <- loess(value ~ age_mid, span = span, weights = total_cases, data = .x)
      fit_lo <- loess(lowercl ~ age_mid, span = span, weights = total_cases, data = .x)
      fit_hi <- loess(uppercl ~ age_mid, span = span, weights = total_cases, data = .x)
      
      tibble(
        age_mid = grid,
        value = predict(fit_val, newdata = tibble(age_mid = grid)),
        lowercl = predict(fit_lo, newdata = tibble(age_mid = grid)),
        uppercl = predict(fit_hi, newdata = tibble(age_mid = grid))
      )}
    ) %>%
    ungroup() %>%
    mutate(lowercl = pmax(0, lowercl),
           value = pmax(0, value),
           uppercl = pmax(0, uppercl))
}

# Compute a raw incidence rate per 100,000 with exact Poisson (gamma-based) confidence intervals
poisson_ci <- function(count, pop, multiplier = 100000, conf.level = 0.95) {
  alpha <- 1 - conf.level
  rate <- count / pop * multiplier
  lci <- ifelse(count == 0, 0, qgamma(alpha/2, shape = count, rate = 1) / pop * multiplier)
  uci <- qgamma(1 - alpha/2, shape = count + 1, rate = 1) / pop * multiplier
  tibble(value = rate, lowercl = lci, uppercl = uci)
}

# Colorectal: combine cases/population across all periods, by age band and race, then compute rates
age_spec_colo_3 <- data %>%
  filter(site == "Colon & rectum") %>%
  group_by(age_mid, racegrp) %>%
  summarise(
    total_cases = sum(n3),
    total_pop   = sum(pop),
    .groups = "drop"
  ) %>%
  bind_cols(poisson_ci(count = .$total_cases, pop = .$total_pop))

age_spec_colo_3 <- set_race_levels(age_spec_colo_3) %>%
  smooth_age()

# Breast: same as above, combined across all periods
age_spec_breast_3 <- data %>%
  filter(site == "Breast") %>%
  group_by(age_mid, racegrp) %>%
  summarise(
    total_cases = sum(n3),
    total_pop   = sum(pop),
    .groups = "drop"
  ) %>%
  bind_cols(poisson_ci(count = .$total_cases, pop = .$total_pop))

age_spec_breast_3 <- set_race_levels(age_spec_breast_3) %>%
  smooth_age()


# Plot smoothed age-specific colorectal cancer incidence curve, by race
age_spec_smooth_colo <- function(df) {
  ggplot(df, aes(x = age_mid, y = value, colour = racegrp, linetype = racegrp)) +
    geom_ribbon(aes(ymin = lowercl, ymax = uppercl, fill = racegrp), colour = NA, alpha = 0.3) +
    geom_line(linewidth = 1) +
    scale_colour_manual(values = mycolours) +
    scale_linetype_manual(values = mylines) +
    scale_fill_manual(values = mycolours_light) +
    labs(x = "Age (years)", y = "Colorectal cancer incidence per 100,000", 
         colour = "Race", linetype = "Race", fill = "Race") +
    theme_classic() +
    scale_x_continuous(breaks = seq(20, 80, by = 10))+
    theme(text = element_text(face = "bold"),
          axis.title.y = element_text(margin = margin(r=10), size = 14),
          axis.title.x = element_text(margin = margin(t=10), size = 14),
          axis.text = element_text(size = 12, colour = "black"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.line = element_line(linewidth = 1),
          axis.ticks = element_line(linewidth = 1),
          legend.text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5))
}

# Plot smoothed age-specific breast cancer incidence curve, by race
age_spec_smooth_breast <- function(df) {
  ggplot(df, aes(x = age_mid, y = value, colour = racegrp, linetype = racegrp)) +
    geom_ribbon(aes(ymin = lowercl, ymax = uppercl, fill = racegrp), colour = NA, alpha = 0.3) +
    geom_line(linewidth = 1) +
    scale_colour_manual(values = mycolours) +
    coord_cartesian(xlim = c(30, 84)) +
    scale_linetype_manual(values = mylines) +
    scale_fill_manual(values = mycolours_light) +
    labs(x = "Age (years)", y = "Breast cancer incidence per 100,000", 
         colour = "Race", linetype = "Race", fill = "Race") +
    theme_classic() +
    scale_x_continuous(breaks = seq(20, 80, by = 10))+
    theme(text = element_text(face = "bold"),
          axis.title.y = element_text(margin = margin(r=10), size = 14),
          axis.title.x = element_text(margin = margin(t=10), size = 14),
          axis.text = element_text(size = 12, colour = "black"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.line = element_line(linewidth = 1),
          axis.ticks = element_line(linewidth = 1),
          legend.text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5))
}

age_spec_smooth_colo_3 <- age_spec_smooth_colo(age_spec_colo_3)
age_spec_smooth_breast_3 <- age_spec_smooth_breast(age_spec_breast_3)


ggarrange(age_spec_smooth_breast_3, age_spec_smooth_colo_3, ncol = 2, nrow = 1, common.legend = T)
ggsave("../results/smooth_age_specific_rate_with_age.pdf", width = 8, height = 4.5)
dev.off()


## ---------------------------------------------------------------------------------
## 7. AGE-SPECIFIC RATES OVER SUCCESSIVE PERIODS (racegrp == "All" only) - FIGURE 3
## ---------------------------------------------------------------------------------

## To look at how age-incidence changes over periods for breast and colorectal cancer

# Overall age-specific colorectal cancer rates stratified by ethnicity and 5-yearly diagnosis period from 1968 – 2023 (as indicated by the first year of the 5-yearly intervals).
poisson_ci <- function(count, pop, multiplier = 100000, conf.level = 0.95) {
  alpha <- 1 - conf.level
  rate <- count / pop * multiplier
  lci <- ifelse(count == 0, 0, qgamma(alpha/2, shape = count, rate = 1) / pop * multiplier)
  uci <- qgamma(1 - alpha/2, shape = count + 1, rate = 1) / pop * multiplier
  tibble(value = rate, lowercl = lci, uppercl = uci)
}

# Colorectal: age-specific rates computed separately for each period (not combined), race group "All" only
age_spec_colo_3_overyears <- data %>%
  filter(site == "Colon & rectum") %>%
  group_by(age_mid, racegrp, period_mid) %>%
  summarise(
    total_cases = sum(n3),
    total_pop   = sum(pop),
    .groups = "drop"
  ) %>%
  bind_cols(poisson_ci(count = .$total_cases, pop = .$total_pop))

age_spec_colo_overyears <- set_race_levels(age_spec_colo_3_overyears) %>%
  filter(racegrp == "All")

# Breast: same as above, by period, race group "All" only
age_spec_breast_3_overyears <- data %>%
  filter(site == "Breast") %>%
  group_by(age_mid, racegrp, period_mid) %>%
  summarise(
    total_cases = sum(n3),
    total_pop   = sum(pop),
    .groups = "drop"
  ) %>%
  bind_cols(poisson_ci(count = .$total_cases, pop = .$total_pop))

age_spec_breast_overyears <- set_race_levels(age_spec_breast_3_overyears) %>%
  filter(racegrp == "All")

# Smooth age-specific rate curve (with CIs) separately within each period, using Loess smoothing
smooth_age2 <- function(df, span = 0.75, n_out = 200) {
  df %>%
    group_by(period_mid) %>%
    group_modify(~{
      grid <- seq(min(.x$age_mid), max(.x$age_mid), length.out = n_out)
      fit_val <- loess(value ~ age_mid, span = span, weights = total_cases, data = .x)
      fit_lo <- loess(lowercl ~ age_mid, span = span, weights = total_cases, data = .x)
      fit_hi <- loess(uppercl ~ age_mid, span = span, weights = total_cases, data = .x)
      
      tibble(
        age_mid = grid,
        value = predict(fit_val, newdata = tibble(age_mid = grid)),
        lowercl = predict(fit_lo, newdata = tibble(age_mid = grid)),
        uppercl = predict(fit_hi, newdata = tibble(age_mid = grid))
      )}
    ) %>%
    ungroup() %>%
    mutate(lowercl = pmax(0, lowercl),
           value = pmax(0, value),
           uppercl = pmax(0, uppercl))
}
age_spec_breast_overyears <- smooth_age2(age_spec_breast_overyears)
age_spec_colo_overyears <- smooth_age2(age_spec_colo_overyears)


# Plot age-specific colorectal cancer incidence, one curve per period, coloured by period
library(viridis)
plot_age_spec_colo <- function(df) {
  ggplot(df, aes(x = age_mid, y = value, colour = period_mid, group = period_mid)) +
    geom_ribbon(aes(ymin = lowercl, ymax = uppercl, fill = period_mid, group = period_mid), colour = NA, alpha = 0.3) +
    geom_line(linewidth = 1) +
    scale_colour_viridis() +
    scale_fill_viridis() +
    labs(x = "Age (years)", y = "Colorectal cancer incidence per 100,000", 
         colour = "Period midyear", fill = "Period midyear") +
    theme_classic() +
    scale_x_continuous(breaks = seq(20, 80, by = 10))+
    theme(text = element_text(face = "bold"),
          axis.title.y = element_text(margin = margin(r=10), size = 14),
          axis.title.x = element_text(margin = margin(t=10), size = 14),
          axis.text = element_text(size = 12, colour = "black"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.line = element_line(linewidth = 1),
          axis.ticks = element_line(linewidth = 1),
          legend.text = element_text(size = 12),
          legend.key.width = unit(2.5, "cm"),
          plot.title = element_text(hjust = 0.5))
}

# Plot age-specific breast cancer incidence, one curve per period, coloured by period
plot_age_spec_breast <- function(df) {
  ggplot(df, aes(x = age_mid, y = value, colour = period_mid, group = period_mid)) +
    geom_ribbon(aes(ymin = lowercl, ymax = uppercl, fill = period_mid, group = period_mid), colour = NA, alpha = 0.3) +
    geom_line(linewidth = 1) +
    scale_colour_viridis() +
    scale_fill_viridis() +
    coord_cartesian(xlim = c(30, 84)) +
    labs(x = "Age (years)", y = "Breast cancer incidence per 100,000", 
         colour = "Period midyear", fill = "Period midyear") +
    theme_classic() +
    scale_x_continuous(breaks = seq(20, 80, by = 10))+
    theme(text = element_text(face = "bold"),
          axis.title.y = element_text(margin = margin(r=10), size = 14),
          axis.title.x = element_text(margin = margin(t=10), size = 14),
          axis.text = element_text(size = 12, colour = "black"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.line = element_line(linewidth = 1),
          axis.ticks = element_line(linewidth = 1),
          legend.text = element_text(size = 12),
          legend.key.width = unit(2.5, "cm"),
          plot.title = element_text(hjust = 0.5))
}

age_spec_colo_overyears_plot <- plot_age_spec_colo(age_spec_colo_overyears)
age_spec_breast_overyears_plot <- plot_age_spec_breast(age_spec_breast_overyears)

## Figure 3. Age-incidence relationship for breast (left) and colorectal (right) cancer in females 
##           aged 25 to 84 in Singapore overall in different periods. Age refers to midyears of 
##           5-year age groups; shading indicates 95 % confidence intervals; curves and confidence
##           intervals smoothed with Loess smoothing, span = 0.75

ggarrange(age_spec_breast_overyears_plot, age_spec_colo_overyears_plot, ncol = 2, nrow = 1, common.legend = T)
ggsave("../results/smooth_age_specific_rate_overtheyears.pdf", width = 8, height = 4.5)
dev.off()


## -----------------------------------------------------------------------------
## 8. FIT AGE-PERIOD-COHORT (APC) MODELS, PER RACE GROUP, FOR EACH CANCER SITE
## -----------------------------------------------------------------------------

# Fit APC model for breast and colorectal cancer in females for each race 

## set references to be lower bounds of ranges
ref_period <- 1970
ref_cohort <- 1890

# Fit a natural-spline APC model (Epi::apc.fit) to a single site/racegrp subset of the data
fit_apc <- function(data) {
  apc.fit(
    ref.c = ref_cohort,
    ref.p = ref_period,
    model = "ns",    
    A = data$age_mid,
    P = data$period_mid,
    D = data$n3,
    Y= data$pop
  )
}

# Fit one APC model per race group, for breast cancer
breast_grouped <- data %>% filter(site == "Breast") %>% group_by(racegrp)
apc_results_breast <- breast_grouped %>%
  group_split() %>%
  setNames(group_keys(breast_grouped)$racegrp) %>%
  lapply(fit_apc)

# Fit one APC model per race group, for colorectal cancer
colo_grouped <- data %>% filter(site == "Colon & rectum") %>% group_by(racegrp)
apc_results_colorectal <- colo_grouped %>%
  group_split() %>%
  setNames(group_keys(colo_grouped)$racegrp) %>%
  lapply(fit_apc)


## -----------------------------------------------------------------------------
## 9. EXTRACT AND PLOT AGE / PERIOD / COHORT COMPONENT EFFECTS - FIGURE 4
## -----------------------------------------------------------------------------

## Component effects for a) age, b) period and c) cohort from full age-period-cohort model.
menopause_age <- 49

###  Age effects -- colorectal cancer

# Pull a named effect component ("Age"/"Per"/"Coh") out of a fitted apc.fit object for a given site
extract_effect <- as.data.frame(function(racegrp, effect, site) {
  if (site == "Breast"){
    apc_results_breast[[racegrp]][[effect]]
  } else if (site == "Colon & rectum"){
    apc_results_colorectal[[racegrp]][[effect]]
  }
})

# Colorectal cancer: age effect, per race group, combined into one data frame
CN_age_effect_c <- extract_effect("CN", "Age", "Colon & rectum") %>%
  mutate(racegrp = "CN")
IN_age_effect_c <- extract_effect("IN", "Age", "Colon & rectum")  %>%
  mutate(racegrp = "IN")
MY_age_effect_c <- extract_effect("MY", "Age", "Colon & rectum")  %>%
  mutate(racegrp = "MY")
All_age_effect_c <- extract_effect("All", "Age", "Colon & rectum")  %>%
  mutate(racegrp = "All")

colorectal_age_effects <- bind_rows(CN_age_effect_c, 
                                    IN_age_effect_c, 
                                    MY_age_effect_c,
                                    All_age_effect_c)


###  Period effects -- colorectal cancer
CN_per_effect_c <- extract_effect("CN", "Per", "Colon & rectum") %>%
  mutate(racegrp = "CN")
IN_per_effect_c <- extract_effect("IN", "Per", "Colon & rectum")  %>%
  mutate(racegrp = "IN")
MY_per_effect_c <- extract_effect("MY", "Per", "Colon & rectum")  %>%
  mutate(racegrp = "MY")
All_per_effect_c <- extract_effect("All", "Per", "Colon & rectum")  %>%
  mutate(racegrp = "All")

colorectal_period_effects <- bind_rows(CN_per_effect_c, 
                                       IN_per_effect_c, 
                                       MY_per_effect_c,
                                       All_per_effect_c)

###  Cohort effects -- colorectal cancer
CN_coh_effect_c <- extract_effect("CN", "Coh", "Colon & rectum") %>%
  mutate(racegrp = "CN")
IN_coh_effect_c <- extract_effect("IN", "Coh", "Colon & rectum")  %>%
  mutate(racegrp = "IN")
MY_coh_effect_c <- extract_effect("MY", "Coh", "Colon & rectum")  %>%
  mutate(racegrp = "MY")
All_coh_effect_c <- extract_effect("All", "Coh", "Colon & rectum")  %>%
  mutate(racegrp = "All")

colorectal_cohort_effects <- bind_rows(CN_coh_effect_c, 
                                       IN_coh_effect_c, 
                                       MY_coh_effect_c,
                                       All_coh_effect_c)

### Generation cutoffs by birth year 
generation_bounds <- data.frame(
  x_min = c(-Inf, 1946, 1965, 1981, 1997),
  x_max = c(1945, 1964, 1980, 1996, Inf),
  label = c("Silent", "Boomers", "Gen X", "Millennials", "Gen Z"),
  label_pos = c(1940, 1955, 1973, 1989, 2003)
)

### plot colorectal cancer age effects
colorectal_age_effects    <- set_race_levels(colorectal_age_effects)
colorectal_period_effects <- set_race_levels(colorectal_period_effects)
colorectal_cohort_effects <- set_race_levels(colorectal_cohort_effects)

# Colorectal age effect: incidence per 100,000 as a function of age, faceted by race, marking menopause age
colorectal_age <-  ggplot(colorectal_age_effects, aes(x = value.Age, y = value.Rate, colour = racegrp, fill = racegrp, linetype = racegrp)) +
  geom_ribbon(aes(ymin = value.2.5., ymax = value.97.5.), colour = NA, alpha = 0.3) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~racegrp) +
  geom_vline(xintercept = menopause_age, linetype = "dashed", colour = "grey30") +
  annotate("text", x = menopause_age - 12, y = 0.011, 
           label = "< 49 years", size = 4) +
  annotate("text", x = menopause_age + 15, y = 0.011, 
           label = "> 49 years", size = 4) +
  labs(title = "Age effect, colorectal", 
       x = "Age (years)", 
       y = "Incidence per 100,000",
       colour = "Race",
       linetype = "Race",
       fill = "Race") +
  theme_classic() + 
  scale_colour_manual(values = mycolours) + 
  scale_fill_manual(values = mycolours_light) +
  scale_linetype_manual(values = mylines) +
  theme(text = element_text(face = "bold"),
        plot.title = element_text(margin = margin(b=10), size = 12),
        axis.title.y = element_text(margin = margin(r=10), size = 12),
        axis.title.x = element_text(margin = margin(t=5), size = 12),
        axis.text = element_text(size = 11, colour = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1),
        legend.text = element_text(size = 11))

### plot colorectal cancer period effects
# Colorectal period effect: rate ratio relative to ref_period, faceted by race
colorectal_period <- ggplot(colorectal_period_effects, aes(x = value.Per, y = value.P.RR, colour = racegrp, fill = racegrp, linetype = racegrp)) +
  geom_ribbon(aes(ymin = value.2.5., ymax = value.97.5.), colour = NA, alpha = 0.3) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  labs(title = "Period effect, colorectal", x = "Year", y = "Rate ratio") +
  theme_classic() + 
  facet_wrap(~racegrp) +
  scale_colour_manual(values = mycolours) + 
  scale_fill_manual(values = mycolours_light) +
  scale_linetype_manual(values = mylines) +
  theme(text = element_text(face = "bold"),
        plot.title = element_text(margin = margin(b=10), size = 12),
        axis.title.y = element_text(margin = margin(r=10), size = 12),
        axis.title.x = element_text(margin = margin(t=5), size = 12),
        axis.text = element_text(size = 11, colour = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1),
        legend.text = element_text(size = 11),
        axis.text.x = element_text(hjust = 1, angle = 45)) 

### plot colorectal cancer cohort effects
# Colorectal cohort effect: rate ratio relative to ref_cohort, faceted by race, with generation bands
colorectal_cohort <- ggplot(colorectal_cohort_effects, aes(x = value.Coh, y = value.C.RR, colour = racegrp, fill = racegrp, linetype = racegrp)) +
  geom_vline(data = generation_bounds, aes(xintercept = x_min), 
             linetype = "dotted", colour = "grey50", inherit.aes = FALSE) +
  coord_cartesian(xlim = c(1888,2000), clip = "off") +
  scale_x_continuous(breaks = c(1890, 1910, 1930, 1950, 1970, 1990)) +
  geom_ribbon(aes(ymin = value.2.5., ymax = value.97.5.), colour = NA, alpha = 0.3) +
  geom_line(linewidth = 0.8) +
  ylim(0,15) +
  facet_wrap(~racegrp) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_text(
    data = generation_bounds,
    aes(x = label_pos, y = 10, label = label),
    inherit.aes = FALSE,
    angle = 90,
    size = 3.5,
    hjust = 0) +
  labs(title = "Cohort effect, colorectal", x = "Birth cohort", y = "Rate ratio") +
  theme_classic() + 
  scale_colour_manual(values = mycolours) + 
  scale_fill_manual(values = mycolours_light) +
  scale_linetype_manual(values = mylines) +
  theme(text = element_text(face = "bold"),
        plot.title = element_text(margin = margin(b=10), size = 12),
        axis.title.y = element_text(margin = margin(r=10), size = 12),
        axis.title.x = element_text(margin = margin(t=5), size = 12),
        axis.text = element_text(size = 11, colour = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1),
        legend.text = element_text(size = 11),
        plot.margin = margin(r = 10),
        axis.text.x = element_text(hjust = 1, angle = 45))

### Age effects -- breast cancer
CN_age_effect_b <- extract_effect("CN", "Age", "Breast") %>%
  mutate(racegrp = "CN")
IN_age_effect_b <- extract_effect("IN", "Age", "Breast")  %>%
  mutate(racegrp = "IN")
MY_age_effect_b <- extract_effect("MY", "Age", "Breast")  %>%
  mutate(racegrp = "MY")
All_age_effect_b <- extract_effect("All", "Age", "Breast")  %>%
  mutate(racegrp = "All")

breast_age_effects <- bind_rows(CN_age_effect_b, 
                                IN_age_effect_b, 
                                MY_age_effect_b,
                                All_age_effect_b)


### Period effects -- breast cancer
CN_per_effect_b <- extract_effect("CN", "Per", "Breast") %>%
  mutate(racegrp = "CN")
IN_per_effect_b <- extract_effect("IN", "Per", "Breast")  %>%
  mutate(racegrp = "IN")
MY_per_effect_b <- extract_effect("MY", "Per", "Breast")  %>%
  mutate(racegrp = "MY")
All_per_effect_b <- extract_effect("All", "Per", "Breast")  %>%
  mutate(racegrp = "All")

breast_period_effects <- bind_rows(CN_per_effect_b, 
                                   IN_per_effect_b, 
                                   MY_per_effect_b,
                                   All_per_effect_b)

### Cohort effects -- breast cancer
CN_coh_effect_b <- extract_effect("CN", "Coh", "Breast") %>%
  mutate(racegrp = "CN")
IN_coh_effect_b <- extract_effect("IN", "Coh", "Breast")  %>%
  mutate(racegrp = "IN")
MY_coh_effect_b <- extract_effect("MY", "Coh", "Breast")  %>%
  mutate(racegrp = "MY")
All_coh_effect_b <- extract_effect("All", "Coh", "Breast")  %>%
  mutate(racegrp = "All")

breast_cohort_effects <- bind_rows(CN_coh_effect_b, 
                                   IN_coh_effect_b, 
                                   MY_coh_effect_b,
                                   All_coh_effect_b)
breast_age_effects        <- set_race_levels(breast_age_effects)
breast_period_effects     <- set_race_levels(breast_period_effects)
breast_cohort_effects     <- set_race_levels(breast_cohort_effects)

### plot breast cancer age effects
# Breast age effect: incidence per 100,000 as a function of age, faceted by race, marking menopause age
breast_age <-  ggplot(breast_age_effects, aes(x = value.Age, y = value.Rate, colour = racegrp, fill = racegrp, linetype = racegrp)) +
  geom_ribbon(aes(ymin = value.2.5., ymax = value.97.5.), colour = NA, alpha = 0.3) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~racegrp) +
  geom_vline(xintercept = menopause_age, linetype = "dashed", colour = "grey30") +
  annotate("text", x = menopause_age - 12, y = 0.0025, 
           label = "< 49 years", size = 4) +
  annotate("text", x = menopause_age + 15, y = 0.0025, 
           label = "> 49 years", size = 4) +
  labs(title = "Age effect, breast", 
       x = "Age (years)", 
       y = "Incidence per 100,000",
       colour = "Race",
       linetype = "Race",
       fill = "Race") +
  theme_classic() + 
  scale_colour_manual(values = mycolours) + 
  scale_fill_manual(values = mycolours_light) +
  scale_linetype_manual(values = mylines) +
  theme(text = element_text(face = "bold"),
        plot.title = element_text(margin = margin(b=10), size = 12),
        axis.title.y = element_text(margin = margin(r=10), size = 12),
        axis.title.x = element_text(margin = margin(t=5), size = 12),
        axis.text = element_text(size = 11, colour = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1),
        legend.text = element_text(size = 11))

### plot breast cancer period effects
# Breast period effect: rate ratio relative to ref_period, faceted by race
breast_period <- ggplot(breast_period_effects, aes(x = value.Per, y = value.P.RR, colour = racegrp, fill = racegrp, linetype = racegrp)) +
  geom_ribbon(aes(ymin = value.2.5., ymax = value.97.5.), colour = NA, alpha = 0.3) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~racegrp) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  labs(title = "Period effect, breast", x = "Year", y = "Rate ratio") +
  theme_classic() + 
  scale_colour_manual(values = mycolours) + 
  scale_fill_manual(values = mycolours_light) +
  scale_linetype_manual(values = mylines) +
  theme(text = element_text(face = "bold"),
        plot.title = element_text(margin = margin(b=10), size = 12),
        axis.title.y = element_text(margin = margin(r=10), size = 12),
        axis.title.x = element_text(margin = margin(t=5), size = 12),
        axis.text = element_text(size = 11, colour = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1),
        legend.text = element_text(size = 11),
        axis.text.x = element_text(hjust = 1, angle = 45)) 

### plot breast cancer cohort effects
# Breast cohort effect: rate ratio relative to ref_cohort, faceted by race, with generation bands
breast_cohort <- ggplot(breast_cohort_effects, aes(x = value.Coh, y = value.C.RR, colour = racegrp, fill = racegrp, linetype = racegrp)) +
  geom_vline(data = generation_bounds, aes(xintercept = x_min), 
             linetype = "dotted", colour = "grey50", inherit.aes = FALSE) +
  coord_cartesian(xlim = c(1888,2000), clip = "off") +
  scale_x_continuous(breaks = c(1890, 1910, 1930, 1950, 1970, 1990)) +
  geom_ribbon(aes(ymin = value.2.5., ymax = value.97.5.), colour = NA, alpha = 0.3) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~racegrp) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_text(
    data = generation_bounds,
    aes(x = label_pos, y = 16.5, label = label),
    inherit.aes = FALSE,
    angle = 90,
    size = 3.5,
    hjust = 0)+
  labs(title = "Cohort effect, breast", x = "Birth cohort", y = "Rate ratio") +
  theme_classic() + 
  scale_colour_manual(values = mycolours) + 
  scale_fill_manual(values = mycolours_light) +
  scale_linetype_manual(values = mylines) +
  theme(text = element_text(face = "bold"),
        plot.title = element_text(margin = margin(b=10), size = 12),
        axis.title.y = element_text(margin = margin(r=10), size = 12),
        axis.title.x = element_text(margin = margin(t=5), size = 12),
        axis.text = element_text(size = 11, colour = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1),
        legend.text = element_text(size = 11),
        plot.margin = margin(r = 10),
        axis.text.x = element_text(hjust = 1, angle = 45))

## arrange all component effects together in 2 panels
## Figure 4. Age, period and cohort effects on breast (left) and colorectal (right) cancer incidence 
##           in females aged 25 to 84 of different ethnicities, 1968-2023. Age, year, and birth cohort 
##          refer to 5-year midyears; CN = Chinese, IN = Indian, MY = Malay, All = CN+IN+MY

plot_labels <- c(
  "a)", "b)", "c)", "d)", "e)", "f)"
)
ggarrange(breast_age, colorectal_age,
          breast_period, colorectal_period,
          breast_cohort,   colorectal_cohort, 
          ncol = 2, nrow = 3,
          common.legend = TRUE, 
          legend = "top",
          heights = c(0.8, 0.8, 1),
          labels = plot_labels,
          font.label = list(size = 14, color = "black", face = "bold"))
ggsave("../results/component_effects.pdf",
       width = 10, height = 15)
dev.off()


## -----------------------------------------------------------------------------
## 10. MODEL DIAGNOSTICS: MULTIPLE-TESTING CORRECTION AND OVERDISPERSION CHECK
## -----------------------------------------------------------------------------

# Apply BH correction across ALL p-values collected (excluding the NA baseline Age rows)

## Collect every Anova table across sites and race groups
all_anova <- bind_rows(
  lapply(names(apc_results_breast), function(r) {
    apc_results_breast[[r]]$Anova %>% mutate(racegrp = r, site = "Breast")
  }),
  lapply(names(apc_results_colorectal), function(r) {
    apc_results_colorectal[[r]]$Anova %>% mutate(racegrp = r, site = "Colon & rectum")
  })
)

## Apply correction
all_anova <- all_anova %>%
  mutate(p_adj_BH = p.adjust(`Pr(>Chi)`, method = "BH"))

## Table 1. Goodness-of-fit measures for APC analyses of colorectal and breast cancer incidence 
##          in females of different ethnicities aged 25 - 84 in Singapore, 1968–2023. AIC refers 
##          to Akaike Information Criterion, df refers to degrees of freedom used in the model, 
##          p refers to the p-value from the likelihood ratio test between models in adjacent 
##          lines corrected with Benjamini-Hochberg method, Dev refers to model deviance, and χ²/df 
##          refers to Pearson χ²/df dispersion ratio.
all_anova 

# Check for overdispersion in the model
check_dispersion <- function(model_fit) {
  mod <- model_fit$Model
  pearson_chisq <- sum(residuals(mod, type = "pearson")^2)
  df <- mod$df.residual
  ratio <- pearson_chisq / df
  tibble(pearson_chisq = pearson_chisq, df = df, dispersion_ratio = ratio)
}

dispersion_breast <- lapply(names(apc_results_breast), function(x) {
  check_dispersion(apc_results_breast[[x]]) %>% mutate(racegrp = x, site = "Breast")
}) %>% bind_rows()

dispersion_colorectal <- lapply(names(apc_results_colorectal), function(x) {
  check_dispersion(apc_results_colorectal[[x]]) %>% mutate(racegrp = x, site = "Colon & rectum")
}) %>% bind_rows()

dispersion_all <- bind_rows(dispersion_breast, dispersion_colorectal)

dispersion_all %>%
  select(racegrp, site, dispersion_ratio)


## -----------------------------------------------------------------------------
## 11. PRE- VS POST-MENOPAUSAL COHORT EFFECT IN BREAST CANCER - FIGURE 5
## -----------------------------------------------------------------------------

# See how cohort effects differ by women under average menopause age of 49 vs 49+ women

# For each race group, test (via nested quasi-Poisson GLMs) whether the cohort effect on breast
# cancer incidence differs significantly between women under 49 and women 49+
cohort_by_age_group_test <- lapply(unique(data$racegrp), function(race) {
  filtered_data <- data %>% filter(site == "Breast", racegrp == race) %>%
    mutate(menopause_grp = if_else(age_mid < 49, "pre", "post"))
  
  no_meno <- glm(n3 ~ ns(cohort_mid, df = 5) + ns(age_mid, df = 5),
                 offset = log(pop), family = quasipoisson, data = filtered_data)
  with_meno <- glm(n3 ~ ns(cohort_mid, df = 5) * menopause_grp + ns(age_mid, df = 5),
                   offset = log(pop), family = quasipoisson, data = filtered_data)
  
  test <- anova(no_meno, with_meno, test = "F")
  tibble(racegrp = race, F_stat = test$F[2], p_value = test$`Pr(>F)`[2])
}) %>% bind_rows()

cohort_by_age_group_test

# For each race group, fit both the null and menopause-interaction models and generate predicted
# incidence rates (with 95% CIs) across birth cohort, separately for the pre- and post-49 age groups
all_pred <- lapply(unique(data$racegrp), function(race) {
  filtered_data <- data %>%
    filter(site == "Breast", racegrp == race) %>%
    mutate(
      menopause_grp = if_else(age_mid < 49, "pre", "post"),
      log_pop = log(pop)
    )
  
  no_meno <- glm(n3 ~ ns(cohort_mid, df = 5) + ns(age_mid, df = 5),
                 offset = log_pop, family = quasipoisson, data = filtered_data)
  with_meno <- glm(n3 ~ ns(cohort_mid, df = 5) * menopause_grp + ns(age_mid, df = 5),
                   offset = log_pop, family = quasipoisson, data = filtered_data)
  
  age_pre  <- median(filtered_data$age_mid[filtered_data$menopause_grp == "pre"])
  age_post <- median(filtered_data$age_mid[filtered_data$menopause_grp == "post"])
  
  grid_pre <- expand.grid(
    cohort_mid = sort(unique(filtered_data %>% filter(menopause_grp == "pre") %>% pull(cohort_mid))),
    menopause_grp = "pre",
    age_mid = age_pre
  ) %>% mutate(log_pop = log(100000))
  
  grid_post <- expand.grid(
    cohort_mid = sort(unique(filtered_data %>% filter(menopause_grp == "post") %>% pull(cohort_mid))),
    menopause_grp = "post",
    age_mid = age_post
  ) %>% mutate(log_pop = log(100000))
  
  grid <- bind_rows(grid_pre, grid_post)
  
  
  predict_from <- function(model, model_label) {
    p <- predict(model, newdata = grid, type = "link", se.fit = TRUE)
    grid %>%
      mutate(
        fit = as.vector(p$fit),
        se = as.vector(p$se.fit),
        rate = exp(fit),
        lowercl = exp(fit - 1.96 * se),
        uppercl = exp(fit + 1.96 * se),
        racegrp = race,
        model = model_label
      )
  }
  
  bind_rows(
    predict_from(no_meno, "Null (no menopause interaction)"),
    predict_from(with_meno, "Stratified (menopause interaction)")
  )
  
}) %>%
  bind_rows() %>%
  mutate(model = factor(model, levels = c("Null (no menopause interaction)",
                                          "Stratified (menopause interaction)")))

## Plot Figure 5. Overall predicted breast cancer incidence among women under the average menopause age 
##                of 49, and women aged 49 and above. Predicted breast cancer incidence for each age 
##                group was derived from a generalised linear model with quasi-Poisson distribution, 
##                including age group (under 49 or 49 +) as an interaction variable, and using median age 
##                for each group as the reference age. Birth cohort refers to 5-year midyears; 
##                CN = Chinese, IN = Indian, MY = Malay, All = CN+IN+MY; Shading indicates 95 % confidence 
##                intervals; curves and confidence intervals smoothed with Loess smoothing, span = 0.75

# Restrict to race group "All" and plot Loess-smoothed predicted incidence by cohort, split by age group
ggplot(all_pred %>% filter(racegrp == "All"), aes(x = cohort_mid, y = rate, fill = menopause_grp, colour = menopause_grp, group = menopause_grp)) +
  geom_smooth(method = "loess", linewidth = 1, alpha = 0.2, span = 0.75) +
  theme_classic() +
  scale_colour_manual(values = c("pre" = "#F8766D", "post" = "#00BFC4"),
                      labels = c("pre" = "Under 49", "post" = "49 +")) +
  scale_fill_manual(values = c("pre" = "#F8766D", "post" = "#00BFC4"),
                    labels = c("pre" = "Under 49", "post" = "49 +")) +
  theme(text = element_text(face = "bold"),
        axis.title.y = element_text(margin = margin(r=5), size = 14),
        axis.title.x = element_text(margin = margin(t=10), size = 14),
        axis.text = element_text(size = 12, colour = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  labs(
    x = "Birth cohort",
    y = "Predicted breast cancer incidence rate per 100,000",
    colour = "Age group",
    fill = "Age group"
  )
ggsave("../results/pre_vs_post_cohort_effects.pdf", width = 7, height = 5)
dev.off()