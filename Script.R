# =====================================================

# =====================================================
# FDI and Carbon Emissions in South Asia
# Independent Econometric Analysis
# Author: Your Name
# =====================================================


# 1. Libraries
library(tidyverse)
library(plm)
library(lmtest)
library(sandwich)
library(ggplot2)

# 2. Load data
co2 <- read.csv("data/co2_emissions.csv")
fdi <- read.csv("data/FDI.csv")
gdp <- read.csv("data/GDP per capita.csv")
energy <- read.csv("data/Energy Use.csv")
urban <- read.csv("data/Urban Population.csv")
reg <- read.csv("data/Regulatory quality.csv")
trade <- read.csv("data/Trade openness.csv")
dominv <- read.csv("data/Domestic Investment.csv")

# 3. Cleaning function
clean_wdi <- function(df, value_name){
  year_cols <- names(df)[5:ncol(df)]
  names(df)[5:ncol(df)] <- substr(gsub("[^0-9]", "", year_cols), 1, 4)
  
  df %>%
    rename(
      Country = Country.Name,
      CountryCode = Country.Code,
      SeriesName = Series.Name,
      SeriesCode = Series.Code
    ) %>%
    pivot_longer(
      cols = 5:ncol(.),
      names_to = "Year",
      values_to = value_name
    ) %>%
    mutate(Year = as.integer(Year))
}

# 4. Convert all datasets to long
co2_long <- clean_wdi(co2, "CO2")
fdi_long <- clean_wdi(fdi, "FDI")
gdp_long <- clean_wdi(gdp, "GDP")
energy_long <- clean_wdi(energy, "Energy")
urban_long <- clean_wdi(urban, "UrbanPop")
reg_long <- clean_wdi(reg, "RegQuality")
trade_long <- clean_wdi(trade, "Trade_openness")
dominv_long <- clean_wdi(dominv, "DomesticInv")

# 5. Merge panel
panel_full <- reduce(
  list(co2_long, fdi_long, gdp_long, energy_long,
       urban_long, reg_long, trade_long, dominv_long),
  left_join,
  by = c("Country", "CountryCode", "Year")
)

# 6. Panel setup
pdata <- pdata.frame(panel_full, index = c("Country", "Year"))

# 7. Models
ols <- lm(CO2 ~ FDI, data = panel_full)
fe <- plm(CO2 ~ FDI, data = pdata, model = "within")
re <- plm(CO2 ~ FDI, data = pdata, model = "random")

# 8. Robust FE
coeftest(fe, vcov = vcovHC(fe, type = "HC1"))

# 9. Plot (Pakistan)
ggplot(
  subset(panel_full, Country == "Pakistan"),
  aes(x = Year)
) +
  geom_line(aes(y = CO2, color = "CO2")) +
  geom_line(aes(y = scales::rescale(FDI, to = range(CO2)), color = "FDI")) +
  labs(title = "FDI vs CO2 Emissions in Pakistan")
# =========================
# B/C: PANEL MODELS (plm)
# =========================

# Packages used for panel regressions + robust SEs
library(plm)
library(lmtest)
library(sandwich)

# 1) Make sure you use the dataset that has NO duplicates
# (if you already fixed duplicates, use that object here)
pdata <- pdata.frame(panel_full, index = c("Country", "Year"))

# --- Model 1: Fixed Effects (FDI only) ---
fe <- plm(CO2 ~ FDI, data = pdata, model = "within")

# Robust SEs (cluster by country)
fe_robust <- coeftest(fe, vcov = vcovHC(fe, type = "HC1", cluster = "group"))
fe_robust

# --- Model 2: Fixed Effects + controls + RegQuality ---
fe_controls <- plm(
  CO2 ~ FDI + GDP + Energy + Trade_openness + UrbanPop + DomesticInv + RegQuality,
  data = pdata, model = "within"
)

fe_controls_robust <- coeftest(
  fe_controls,
  vcov = vcovHC(fe_controls, type = "HC1", cluster = "group")
)
fe_controls_robust

# --- Model 3: Fixed Effects + controls + Sustainability ---
fe_sustain <- plm(
  CO2 ~ FDI + GDP + Energy + Trade_openness + UrbanPop + DomesticInv + Sustainability,
  data = pdata, model = "within"
)

fe_sustain_robust <- coeftest(
  fe_sustain,
  vcov = vcovHC(fe_sustain, type = "HC1", cluster = "group")
)
fe_sustain_robust
#Removing duplicates
library(dplyr)
panel_full_clean <- panel_full %>%
  +     group_by(Country, Year) %>%
  +     summarise(
    +         across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
    +         across(where(~!is.numeric(.x)), ~ dplyr::first(.x)),
    +         .groups = "drop")
 panel_full_clean <- panel_full[!duplicated(panel_full[c("Country","Year")]), ]
 table(duplicated(panel_full_clean[c("Country","Year")]))
#Pdata frame
library(plm)
panel_full <- panel_full %>% mutate(Country = as.factor(Country),     Year = as.integer(Year))
sum(duplicated(panel_full[c("Country", "Year")]))
pdata<- pdata.frame(panel_full, index=c ("Country","Year"), drop.index=TRUE, row.names=TRUE)
pdim(pdata)
head(index(pdata))
ols <- lm(CO2 ~ FDI, data = panel_full)
summary(ols)
fe <- plm(CO2 ~ FDI, data = pdata, model = "within")
lmtest::coeftest(fe, vcov = sandwich::vcovHC(fe, type="HC1", cluster="group"))
fe_controls <- plm(CO2 ~ FDI + GDP + Energy + Trade_openness + UrbanPop + DomesticInv + RegQuality,
                   data = pdata, model = "within")
lmtest::coeftest(fe_controls, vcov = sandwich::vcovHC(fe_controls, type="HC1", cluster="group"))
fe_sustain <- plm(CO2 ~ FDI + GDP + Energy + Trade_openness + UrbanPop + DomesticInv + Sustainability,
                  data = pdata, model = "within")
lmtest::coeftest(fe_sustain, vcov = sandwich::vcovHC(fe_sustain, type="HC1", cluster="group"))
re <- plm(CO2 ~ FDI, data = pdata, model = "random")
phtest(fe, re)
sink("outputs/reg_results.txt")
print(summary(ols))
print(lmtest::coeftest(fe, vcov = sandwich::vcovHC(fe, type="HC1", cluster="group")))
print(lmtest::coeftest(fe_controls, vcov = sandwich::vcovHC(fe_controls, type="HC1", cluster="group")))
print(lmtest::coeftest(fe_sustain, vcov = sandwich::vcovHC(fe_sustain, type="HC1", cluster="group")))
sink()
file.show("outputs/reg_results.txt")
pakistan_data <- subset(panel_full, Country == "Pakistan")

lm_pak <- lm(
  CO2 ~ FDI + GDP + Energy + Trade_openness + DomesticInv,
  data = pakistan_data)
summary(lm_pak)
sink("outputs/reg_results_pakistan.txt")
summary(lm_pak)
sink()
library(sandwich)
library(lmtest)
coeftest(lm_pak, vcov = vcovHC(lm_pak, type = "HC1"))
library(modelsummary)
vc <- function(m) sandwich::vcovHC(m, type = "HC1", cluster = "group")

modelsummary(list("OLS" = ols,
      +         "FE (FDI only)" = fe,
      +         "FE + controls" = fe_controls,
      +         "FE + sustainability" = fe_sustain),
    +     vcov = vc,
    +     stars = TRUE,
    +     output = "outputs/reg_table.tex"
    library(plm)
  fe_controls_twfe <- plm(  CO2 ~ FDI + GDP + Energy + Trade_openness + UrbanPop + DomesticInv + RegQuality,  data = pdata,  model = "within", effect = "twoways")        
  lmtest::coeftest(fe_controls_twfe,vcov = sandwich::vcovHC(fe_controls_twfe, type = "HC1", cluster = "group"))
  library(modelsummary)
  modelsummary(list("OLS" = ols,
      +         "FE (FDI only)" = fe,
      +         "FE + controls" = fe_controls,
      +         "TWFE (main)" = fe_controls_twfe ),    vcov = vc,
    +     stars = TRUE,
    +     output = "outputs/reg_table.tex" )
  modelsummary( list("OLS" = ols,
      +         "FE (country only)" = fe_controls,
      +         "FE (country + year)" = fe_controls_twfe ),
    +     vcov = vc,
    +     stars = TRUE,
    +     gof_omit = "AIC|BIC|Log.Lik.|RMSE|Std.Errors",
    +     output = "outputs/reg_table.tex")
  # =========================
  # DATA SECTION: SUMMARY TABLES
  # Plug-and-play script
  # =========================
  
  # Packages (install if missing)
  pkgs <- c("tidyverse", "skimr", "knitr", "kableExtra", "modelsummary")
  to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if (length(to_install) > 0) install.packages(to_install)
  
  library(tidyverse)
  library(skimr)
  library(knitr)
  library(kableExtra)
  library(modelsummary)
  
  # -------------------------
  # 0) Quick sanity checks
  # -------------------------
  stopifnot(exists("panel_full"))
  
  # Make sure types are right
  panel_full <- panel_full %>%
    mutate(
      Year = as.integer(Year),
      Country = as.factor(Country),
      across(c(CO2, FDI, GDP, Energy, Trade_openness, UrbanPop, DomesticInv,
               RegQuality, Sustainability),
             ~ suppressWarnings(as.numeric(.)))
    )
  
  # Optional: drop rows missing core vars for the main model (edit if needed)
  core_vars <- c("CO2","FDI","GDP","Energy","Trade_openness","UrbanPop","DomesticInv")
  df_core <- panel_full %>% filter(if_all(all_of(core_vars), ~ !is.na(.)))
  
  # -------------------------
  # 1) Table A: Sample coverage (N by country, years covered)
  # -------------------------
  sample_coverage <- panel_full %>%
    group_by(Country) %>%
    summarise(
      N_total = n(),
      Year_min = min(Year, na.rm = TRUE),
      Year_max = max(Year, na.rm = TRUE),
      N_years = n_distinct(Year),
      N_CO2 = sum(!is.na(CO2)),
      N_FDI = sum(!is.na(FDI)),
      .groups = "drop"
    ) %>%
    arrange(Country)
  
  kable(sample_coverage, caption = "Sample coverage by country") %>%
    kable_styling(full_width = FALSE)
  
  # -------------------------
  # 2) Table B: Missingness overview (share missing per variable)
  # -------------------------
  missingness <- panel_full %>%
    summarise(across(where(is.numeric),
                     ~ mean(is.na(.)))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "share_missing") %>%
    arrange(desc(share_missing))
  
  kable(missingness, digits = 3, caption = "Missing data (share missing by variable)") %>%
    kable_styling(full_width = FALSE)
  
  # -------------------------
  # 3) Table C: Descriptive statistics (overall)
  # -------------------------
  desc_vars <- c("CO2","FDI","GDP","Energy","Trade_openness","UrbanPop","DomesticInv",
                 "RegQuality","Sustainability")
  
  desc_overall <- panel_full %>%
    select(any_of(desc_vars)) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    group_by(variable) %>%
    summarise(
      N = sum(!is.na(value)),
      Mean = mean(value, na.rm = TRUE),
      SD = sd(value, na.rm = TRUE),
      Min = min(value, na.rm = TRUE),
      P25 = quantile(value, 0.25, na.rm = TRUE),
      Median = median(value, na.rm = TRUE),
      P75 = quantile(value, 0.75, na.rm = TRUE),
      Max = max(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  kable(desc_overall, digits = 3, caption = "Descriptive statistics (overall)") %>%
    kable_styling(full_width = FALSE)
  
  # -------------------------
  # 4) Table D: Descriptive statistics by country (mean/sd)
  # -------------------------
  desc_by_country <- panel_full %>%
    group_by(Country) %>%
    summarise(across(all_of(desc_vars),
                     list(mean = ~ mean(.x, na.rm = TRUE),
                          sd   = ~ sd(.x, na.rm = TRUE)),
                     .names = "{.col}_{.fn}"),
              .groups = "drop") %>%
    arrange(Country)
  
  kable(desc_by_country, digits = 3, caption = "Descriptive statistics by country (mean, sd)") %>%
    kable_styling(full_width = FALSE) %>%
    scroll_box(height = "450px")
  
  # -------------------------
  # 5) Export tables (optional)
  # -------------------------
  # Write CSVs (handy for appendix or checking)
  write_csv(sample_coverage, "table_sample_coverage.csv")
  write_csv(missingness, "table_missingness.csv")
  write_csv(desc_overall, "table_desc_overall.csv")
  write_csv(desc_by_country, "table_desc_by_country.csv")
  
  message("Done. Tables printed + CSVs exported to your working directory.")