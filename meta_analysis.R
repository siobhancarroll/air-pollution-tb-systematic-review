library(readODS)
library(tidyverse)
library(meta)
library(here)

# Load the extraction sheets ##################################################
HAP_data <- read_ods(here("HAP_final_dataset.ods"))
AAP_data <- read_ods(here("AAP_final_dataset.ods"))

# Set output directory for all tables
if (!dir.exists(here("figures"))) {
  dir.create(here("figures"), recursive = TRUE)
}

fp <- here("figures")

###############################################################################
################################# AAP #########################################
###############################################################################

# PM2.5 #######################################################################
PM2.5_incidence <- AAP_data %>%
  filter(Pollutant == "PM2.5",
         primary_outcomes == "TB")

PM2.5_incidence_overall <- AAP_data %>%
  filter(Pollutant == "PM2.5" & 
           primary_outcomes == "TB" &
           (!is.na(per) & per != "unclear") &
           (flag == 1 | is.na(flag)) &
           (type_of_effect_estimate == "RR" | type_of_effect_estimate == "OR") &
           (lag_specific_or_cumulative == "cumulative" | is.na(lag_specific_or_cumulative))) %>%
  mutate(effect_estimate_standardized = 
           case_when(
             per == "1 ug/m3" ~ effect_estimate^10,
             per == "34.7 ug/m3" ~ effect_estimate^(10/34.7),
             per == "50 ug/m3" ~ effect_estimate^(10/50),
             per == "7.78 ug/m3" ~ effect_estimate^(10/7.78),
             per == "10 ug/m3" ~ effect_estimate)) %>%
  mutate(lower_ci_standardized = 
           case_when(
             per == "1 ug/m3" ~ lower_ci^10,
             per == "34.7 ug/m3" ~ lower_ci^(10/34.7),
             per == "50 ug/m3" ~ lower_ci^(10/50),
             per == "7.78 ug/m3" ~ lower_ci^(10/7.78),
             per == "10 ug/m3" ~ lower_ci)) %>%
  mutate(upper_ci_standardized = 
           case_when(
             per == "1 ug/m3" ~ upper_ci^10,
             per == "34.7 ug/m3" ~ upper_ci^(10/34.7),
             per == "50 ug/m3" ~ upper_ci^(10/50),
             per == "7.78 ug/m3" ~ upper_ci^(10/7.78),
             per == "10 ug/m3" ~ upper_ci))

PM2.5_incidence_shortterm <- PM2.5_incidence_overall %>%
  filter(short_or_long_term == "short")
PM2.5_incidence_longterm <- PM2.5_incidence_overall %>%
  filter(short_or_long_term == "long")

PM2.5_incidence_shortterm_meta <- 
  metagen(data = PM2.5_incidence_shortterm,
          TE = log(effect_estimate_standardized),
          lower = log(lower_ci_standardized),
          upper = log(upper_ci_standardized),
          studlab = study_ID,
          sm = "RR",
          common = FALSE,
          random = TRUE,
          method.tau = "DL",
          method.bias = "Egger")
summary(PM2.5_incidence_shortterm_meta)

PM2.5_incidence_longterm_meta <- 
  metagen(data = PM2.5_incidence_longterm,
          TE = log(effect_estimate_standardized),
          lower = log(lower_ci_standardized),
          upper = log(upper_ci_standardized),
          studlab = study_ID,
          sm = "RR",
          common = FALSE,
          random = TRUE,
          method.tau = "DL",
          method.bias = "Egger")
summary(PM2.5_incidence_longterm_meta)

path <- paste0(fp, "/forest-plot-PM2.5-incidence-shortterm.png")
png(file = path, height = "3.5", width = "11", units = "in", res = 300)
forest(PM2.5_incidence_shortterm_meta,
            sortvar = effect_estimate_standardized,
            col.square = "turquoise",
            col.square.lines = "black",
       leftcols = c("studlab", "country", "age_group", "Sex", "study_design"),
       leftlabs = c("Study", "Location", "Age Group", "Sex", "Study Design"),
            overall = TRUE,
            hetstat = TRUE,
       print.tau2 = FALSE,
       print.pval.Q = FALSE)
dev.off()

path <- paste0(fp, "/forest-plot-PM2.5-incidence-shortterm-simple.png")
png(file = path, height = "3.5", width = "7", units = "in", res = 300)
forest(PM2.5_incidence_shortterm_meta,
       sortvar = effect_estimate_standardized,
       col.square = "turquoise",
       col.square.lines = "black",
       leftcols = c("studlab"),
       leftlabs = c("Study"),
       overall = TRUE,
       hetstat = TRUE,
       print.tau2 = FALSE,
       print.pval.Q = FALSE)
dev.off()

path <- paste0(fp, "/forest-plot-PM2.5-incidence-longterm.png")
png(file = path, height = "3.5", width = "11", units = "in", res = 300)
forest(PM2.5_incidence_longterm_meta,
       sortvar = effect_estimate_standardized,
       col.square = "turquoise",
       col.square.lines = "black",
       leftcols = c("studlab", "country", "age_group", "Sex", "study_design"),
       leftlabs = c("Study", "Location", "Age Group", "Sex", "Study Design"),
       overall = TRUE,
       hetstat = TRUE,
       print.tau2 = FALSE,
       print.pval.Q = FALSE)
dev.off()

path <- paste0(fp, "/forest-plot-PM2.5-incidence-longterm-simple.png")
png(file = path, height = "3.5", width = "7", units = "in", res = 300)
forest(PM2.5_incidence_longterm_meta,
       sortvar = effect_estimate_standardized,
       col.square = "turquoise",
       col.square.lines = "black",
       leftcols = c("studlab"),
       leftlabs = c("Study"),
       overall = TRUE,
       hetstat = TRUE,
       print.tau2 = FALSE,
       print.pval.Q = FALSE)
dev.off()

# Add a grouping variable to each object
PM2.5_incidence_shortterm_meta$group <- "Short-Term"
PM2.5_incidence_longterm_meta$group <- "Long-Term"


# PM10 ########################################################################
PM10_incidence <- AAP_data %>%
  filter(Pollutant == "PM10",
         primary_outcomes == "TB")

PM10_incidence_overall <- AAP_data %>%
  filter(Pollutant == "PM10" & 
           primary_outcomes == "TB" &
           (!is.na(per) & per != "unclear") &
           (flag == 1 | is.na(flag)) &
           (type_of_effect_estimate == "RR" | type_of_effect_estimate == "OR") &
           (lag_specific_or_cumulative == "cumulative" | is.na(lag_specific_or_cumulative))) %>%
  mutate(effect_estimate_standardized = 
           case_when(
             per == "1 ug/m3" ~ effect_estimate^10,
             per == "interquartile increase (4.1 ug/m3)" ~ effect_estimate^(10/4.1),
             per == "5.63 ug/m3" ~ effect_estimate^(10/5.63),
             per == "17.35 ug/m3" ~ effect_estimate^(10/17.35),
             per == "97 ug/m3" ~ effect_estimate^(10/97),
             per == "10 ug/m3" ~ effect_estimate)) %>%
  mutate(lower_ci_standardized = 
           case_when(
             per == "1 ug/m3" ~ lower_ci^10,
             per == "interquartile increase (4.1 ug/m3)" ~ lower_ci^(10/4.1),
             per == "5.63 ug/m3" ~ lower_ci^(10/5.63),
             per == "17.35 ug/m3" ~ lower_ci^(10/17.35),
             per == "97 ug/m3" ~ lower_ci^(10/97),
             per == "10 ug/m3" ~ lower_ci)) %>%
  mutate(upper_ci_standardized = 
           case_when(
             per == "1 ug/m3" ~ upper_ci^10,
             per == "interquartile increase (4.1 ug/m3)" ~ upper_ci^(10/4.1),
             per == "5.63 ug/m3" ~ upper_ci^(10/5.63),
             per == "17.35 ug/m3" ~ upper_ci^(10/17.35),
             per == "97 ug/m3" ~ upper_ci^(10/97),
             per == "10 ug/m3" ~ upper_ci))

PM10_incidence_shortterm <- PM10_incidence_overall %>%
  filter(short_or_long_term == "short")
PM10_incidence_longterm <- PM10_incidence_overall %>%
  filter(short_or_long_term == "long")

PM10_incidence_longterm_meta <- 
  metagen(data = PM10_incidence_longterm,
          TE = log(effect_estimate_standardized),
          lower = log(lower_ci_standardized),
          upper = log(upper_ci_standardized),
          studlab = study_ID,
          sm = "RR",
          common = FALSE,
          random = TRUE,
          method.tau = "DL",
          method.bias = "Egger")
summary(PM10_incidence_longterm_meta)

path <- paste0(fp, "/forest-plot-PM10-incidence-longterm.png")
png(file = path, height = "3.5", width = "12", units = "in", res = 300)
forest(PM10_incidence_longterm_meta,
       sortvar = effect_estimate_standardized,
       col.square = "turquoise",
       col.square.lines = "black",
       leftcols = c("studlab", "country", "age_group", "Sex", "study_design"),
       leftlabs = c("Study", "Location", "Age Group", "Sex", "Study Design"),
       overall = TRUE,
       hetstat = TRUE,
       prediction = TRUE,
       print.tau2 = FALSE,
       print.pval.Q = FALSE)
dev.off()

# SO2 #########################################################################
SO2_incidence <- AAP_data %>%
  filter(Pollutant == "SO2",
         primary_outcomes == "TB")

SO2_incidence_overall <- AAP_data %>%
  filter(Pollutant == "SO2" & 
           primary_outcomes == "TB" &
           (!is.na(per) & per != "unclear") &
           (flag == 1 | is.na(flag)) &
           (type_of_effect_estimate == "RR" | type_of_effect_estimate == "OR")&
           (lag_specific_or_cumulative == "cumulative" | is.na(lag_specific_or_cumulative))) %>%
  mutate(effect_estimate_standardized = 
           case_when(
             per == "1 ug/m3" ~ effect_estimate^10,
             per == "5 ug/m3" ~ effect_estimate^2,
             per == "8.59 ug/m3" ~ effect_estimate^(10/8.59),
             per == "10 ug/m3" ~ effect_estimate)) %>%
  mutate(lower_ci_standardized = 
           case_when(
             per == "1 ug/m3" ~ lower_ci^10,
             per == "5 ug/m3" ~ lower_ci^2,
             per == "8.59 ug/m3" ~ lower_ci^(10/8.59),
             per == "10 ug/m3" ~ lower_ci)) %>%
  mutate(upper_ci_standardized = 
           case_when(
             per == "1 ug/m3" ~ upper_ci^10,
             per == "5 ug/m3" ~ upper_ci^2,
             per == "8.59 ug/m3" ~ lower_ci^(10/8.59),
             per == "10 ug/m3" ~ upper_ci))

SO2_incidence_shortterm <- SO2_incidence_overall %>%
  filter(short_or_long_term == "short")
SO2_incidence_longterm <- SO2_incidence_overall %>%
  filter(short_or_long_term == "long") %>%
  filter(study_ID != "Li JX 2023")

SO2_incidence_shortterm_meta <- 
  metagen(data = SO2_incidence_shortterm,
          TE = log(effect_estimate_standardized),
          lower = log(lower_ci_standardized),
          upper = log(upper_ci_standardized),
          studlab = study_ID,
          sm = "RR",
          common = FALSE,
          random = TRUE,
          method.tau = "DL",
          method.bias = "Egger")
summary(SO2_incidence_shortterm_meta)

path <- paste0(fp, "/forest-plot-SO2-incidence-shortterm.png")
png(file = path, height = "3.5", width = "11", units = "in", res = 300)
forest(SO2_incidence_shortterm_meta,
       sortvar = effect_estimate_standardized,
       col.square = "turquoise",
       col.square.lines = "black",
       leftcols = c("studlab", "country", "age_group", "Sex", "study_design"),
       leftlabs = c("Study", "Location", "Age Group", "Sex", "Study Design"),
       overall = TRUE,
       hetstat = TRUE,
       prediction = TRUE,
       print.tau2 = FALSE,
       print.pval.Q = FALSE)
dev.off()

SO2_incidence_longterm_meta <- 
  metagen(data = SO2_incidence_longterm,
          TE = log(effect_estimate_standardized),
          lower = log(lower_ci_standardized),
          upper = log(upper_ci_standardized),
          studlab = study_ID,
          sm = "RR",
          common = FALSE,
          random = TRUE,
          method.tau = "DL",
          method.bias = "Egger")
summary(SO2_incidence_longterm_meta)

path <- paste0(fp, "/forest-plot-SO2-incidence-longterm.png")
png(file = path, height = "3.5", width = "12", units = "in", res = 300)
forest(SO2_incidence_longterm_meta,
       sortvar = effect_estimate_standardized,
       col.square = "turquoise",
       col.square.lines = "black",
       leftcols = c("studlab", "country", "age_group", "Sex", "study_design"),
       leftlabs = c("Study", "Location", "Age Group", "Sex", "Study Design"),
       overall = TRUE,
       hetstat = TRUE,
       prediction = TRUE,
       print.tau2 = FALSE,
       print.pval.Q = FALSE)
dev.off()

# NO2 #########################################################################
NO2_incidence <- AAP_data %>%
  filter(Pollutant == "NO2",
         primary_outcomes == "TB")

NO2_incidence_overall <- AAP_data %>%
  filter(Pollutant == "NO2" & 
           primary_outcomes == "TB" &
           (!is.na(per) & per != "unclear") &
           (flag == 1 | is.na(flag)) &
           (type_of_effect_estimate == "RR" | type_of_effect_estimate == "OR") &
           (lag_specific_or_cumulative == "cumulative" | is.na(lag_specific_or_cumulative))) %>%
  mutate(effect_estimate_standardized = 
           case_when(
             per == "1 ug/m3" ~ effect_estimate^10,
             per == "5 ug/m3" ~ effect_estimate^2,
             per == "1.79 ug/m3" ~ effect_estimate^(10/1.79),
             per == "10 ug/m3" ~ effect_estimate)) %>%
  mutate(lower_ci_standardized = 
           case_when(
             per == "1 ug/m3" ~ lower_ci^10,
             per == "5 ug/m3" ~ lower_ci^2,
             per == "1.79 ug/m3" ~ lower_ci^(10/1.79),
             per == "10 ug/m3" ~ lower_ci)) %>%
  mutate(upper_ci_standardized = 
           case_when(
             per == "1 ug/m3" ~ upper_ci^10,
             per == "5 ug/m3" ~ upper_ci^2,
             per == "1.79 ug/m3" ~ upper_ci^(10/1.79),
             per == "10 ug/m3" ~ upper_ci))

NO2_incidence_shortterm <- NO2_incidence_overall %>%
  filter(short_or_long_term == "short")
NO2_incidence_longterm <- NO2_incidence_overall %>%
  filter(short_or_long_term == "long") %>%
  filter(study_ID != "Chen KY 2016" & study_ID != "Hwang 2014")

NO2_incidence_shortterm_meta <- 
  metagen(data = NO2_incidence_shortterm,
          TE = log(effect_estimate_standardized),
          lower = log(lower_ci_standardized),
          upper = log(upper_ci_standardized),
          studlab = study_ID,
          sm = "RR",
          common = FALSE,
          random = TRUE,
          method.tau = "DL",
          method.bias = "Egger")
summary(NO2_incidence_shortterm_meta)

path <- paste0(fp, "/forest-plot-NO2-incidence-shortterm.png")
png(file = path, height = "4", width = "11", units = "in", res = 300)
forest(NO2_incidence_shortterm_meta,
       sortvar = effect_estimate_standardized,
       col.square = "turquoise",
       col.square.lines = "black",
       leftcols = c("studlab", "country", "age_group", "Sex", "study_design"),
       leftlabs = c("Study", "Location", "Age Group", "Sex", "Study Design"),
       overall = TRUE,
       hetstat = TRUE,
       prediction = TRUE,
       print.tau2 = FALSE,
       print.pval.Q = FALSE)
dev.off()

NO2_incidence_longterm_meta <- 
  metagen(data = NO2_incidence_longterm,
          TE = log(effect_estimate_standardized),
          lower = log(lower_ci_standardized),
          upper = log(upper_ci_standardized),
          studlab = study_ID,
          sm = "RR",
          common = FALSE,
          random = TRUE,
          method.tau = "DL",
          method.bias = "Egger")
summary(NO2_incidence_longterm_meta)

path <- paste0(fp, "/forest-plot-NO2-incidence-longterm.png")
png(file = path, height = "3", width = "11", units = "in", res = 300)
forest(NO2_incidence_longterm_meta,
       sortvar = effect_estimate_standardized,
       col.square = "turquoise",
       col.square.lines = "black",
       leftcols = c("studlab", "country", "age_group", "Sex", "study_design"),
       leftlabs = c("Study", "Location", "Age Group", "Sex", "Study Design"),
       overall = TRUE,
       hetstat = TRUE,
       prediction = TRUE,
       print.tau2 = FALSE,
       print.pval.Q = FALSE)
dev.off()

# CO ##########################################################################
CO_incidence <- AAP_data %>%
  filter(Pollutant == "CO",
         primary_outcomes == "TB")

CO_incidence_overall <- AAP_data %>%
  filter(Pollutant == "CO" & 
           primary_outcomes == "TB" &
           (!is.na(per) & per != "unclear") &
           (flag == 1 | is.na(flag)) &
           (type_of_effect_estimate == "RR" | type_of_effect_estimate == "OR") &
           (lag_specific_or_cumulative == "cumulative" | is.na(lag_specific_or_cumulative))) %>%
  mutate(effect_estimate_standardized = 
           case_when(
             per == "1 ug/m3" ~ effect_estimate^10,
             per == "100 ug/m3" ~ effect_estimate^(10/100),
             per == "10 ug/m3" ~ effect_estimate)) %>%
  mutate(lower_ci_standardized = 
           case_when(
             per == "1 ug/m3" ~ lower_ci^10,
             per == "100 ug/m3" ~ lower_ci^(10/100),
             per == "10 ug/m3" ~ lower_ci)) %>%
  mutate(upper_ci_standardized = 
           case_when(
             per == "1 ug/m3" ~ upper_ci^10,
             per == "100 ug/m3" ~ upper_ci^(10/100),
             per == "10 ug/m3" ~ upper_ci))

CO_incidence_shortterm <- CO_incidence_overall %>%
  filter(short_or_long_term == "short")
CO_incidence_longterm <- CO_incidence_overall %>%
  filter(short_or_long_term == "long")

# O3 ##########################################################################
O3_incidence <- AAP_data %>%
  filter(Pollutant == "O3",
         primary_outcomes == "TB")

O3_incidence_overall <- AAP_data %>%
  filter(Pollutant == "O3" & 
           primary_outcomes == "TB" &
           (!is.na(per) & per != "unclear") &
           (flag == 1 | is.na(flag)) &
           (type_of_effect_estimate == "RR" | type_of_effect_estimate == "OR") &
           (lag_specific_or_cumulative == "cumulative" | is.na(lag_specific_or_cumulative))) %>%
  mutate(effect_estimate_standardized = 
           case_when(
             per == "1 ug/m3" ~ effect_estimate^10,
             per == "5 ug/m3" ~ effect_estimate^2,
             per == "10 mg/m3" ~ effect_estimate^(1/1000),
             per == "10 ug/m3" ~ effect_estimate)) %>%
  mutate(lower_ci_standardized = 
           case_when(
             per == "1 ug/m3" ~ lower_ci^10,
             per == "5 ug/m3" ~ lower_ci^2,
             per == "10 mg/m3" ~ lower_ci^(1/1000),
             per == "10 ug/m3" ~ lower_ci)) %>%
  mutate(upper_ci_standardized = 
           case_when(
             per == "1 ug/m3" ~ upper_ci^10,
             per == "5 ug/m3" ~ upper_ci^2,
             per == "10 mg/m3" ~ upper_ci^(1/1000),
             per == "10 ug/m3" ~ upper_ci))

O3_incidence_shortterm <- O3_incidence_overall %>%
  filter(short_or_long_term == "short")
O3_incidence_longterm <- O3_incidence_overall %>%
  filter(short_or_long_term == "long")

O3_incidence_shortterm_meta <- 
  metagen(data = O3_incidence_shortterm,
          TE = log(effect_estimate_standardized),
          lower = log(lower_ci_standardized),
          upper = log(upper_ci_standardized),
          studlab = study_ID,
          sm = "RR",
          common = FALSE,
          random = TRUE,
          method.tau = "DL",
          method.bias = "Egger")
summary(O3_incidence_shortterm_meta)

path <- paste0(fp, "/forest-plot-O3-incidence-shortterm.png")
png(file = path, height = "3", width = "11", units = "in", res = 300)
forest(O3_incidence_shortterm_meta,
       sortvar = effect_estimate_standardized,
       col.square = "turquoise",
       col.square.lines = "black",
       leftcols = c("studlab", "country", "age_group", "Sex", "study_design"),
       leftlabs = c("Study", "Location", "Age Group", "Sex", "Study Design"),
       overall = TRUE,
       hetstat = TRUE,
       prediction = TRUE,
       print.tau2 = FALSE,
       print.pval.Q = FALSE)
dev.off()


###############################################################################
# MORTALITY ###################################################################
###############################################################################

PM2.5_mortality <- AAP_data %>%
  filter(Pollutant == "PM2.5", primary_outcomes == "Death")

CO_mortality <- AAP_data %>%
  filter(Pollutant == "CO", primary_outcomes == "Death")

O3_mortality <- AAP_data %>%
  filter(Pollutant == "O3", primary_outcomes == "Death")

SO2_mortality <- AAP_data %>%
  filter(Pollutant == "SO2", primary_outcomes == "Death")

NO2_mortality <- AAP_data %>%
  filter(Pollutant == "NO2", primary_outcomes == "Death")

dust_mortality <- AAP_data %>%
  filter((Pollutant == "Dust" | 
            Pollutant == "SPM" | 
            Pollutant == "PM2.5-10") & primary_outcomes == "Death")

traffic_mortality <- AAP_data %>%
  filter(Pollutant == "Traffic", primary_outcomes == "Death")

###############################################################################
# Drug Resistance #############################################################
###############################################################################

PM2.5_resistance <- AAP_data %>%
  filter(Pollutant == "PM2.5",
         primary_outcomes == "Drug resistance")

PM10_resistance <- AAP_data %>%
  filter(Pollutant == "PM10",
         primary_outcomes == "Drug resistance")

SO2_resistance <- AAP_data %>%
  filter(Pollutant == "SO2",
         primary_outcomes == "Drug resistance")

NO2_resistance <- AAP_data %>%
  filter(Pollutant == "NO2",
         primary_outcomes == "Drug resistance")

CO_resistance <- AAP_data %>%
  filter(Pollutant == "CO",
         primary_outcomes == "Drug resistance")

O3_resistance <- AAP_data %>%
  filter(Pollutant == "O3",
         primary_outcomes == "Drug resistance")

dust_resistance <- AAP_data %>%
  filter(Pollutant == "Dust",
         primary_outcomes == "Drug resistance")

###############################################################################
################################## HAP ########################################
###############################################################################

HAP_incidence <- HAP_data %>%
  filter(primary_outcomes == "TB") %>%
  filter(!is.na(effect_est) & !is.na(lower.ci) & !is.na(upper.ci)) %>%
  filter(exposure_category != "direct measurement") %>%
  filter(type_of_effect_est == "OR" | type_of_effect_est == "RR")

HAP_incidence_overall <- HAP_incidence %>%
  filter(flag == 1 | is.na(flag))

HAP_incidence_overall_meta <- 
  metagen(data = HAP_incidence_overall,
          TE = log(effect_est),
          lower = log(lower.ci),
          upper = log(upper.ci),
          studlab = study_ID,
          sm = "RR",
          common = FALSE,
          random = TRUE,
          method.tau = "DL",
          method.bias = "Egger")
summary(HAP_incidence_overall_meta)

path <- paste0(fp, "/forest-plot-HAP-incidence-overall.png")
png(file = path, height = "10", width = "13", units = "in", res = 300)
forest(HAP_incidence_overall_meta,
       sortvar = effect_est,
       col.square = "turquoise",
       col.square.lines = "black",
       leftcols = c("studlab", "country", "age_group", "sex", "study_design", "simplified_n"),
       leftlabs = c("Study", "Location", "Age Group", "Sex", "Study Design", "Sample Size"),
       overall = TRUE,
       hetstat = TRUE,
       print.tau2 = FALSE,
       print.pval.Q = FALSE,
       prediction = TRUE)
dev.off()

path <- paste0(fp, "/forest-plot-HAP-incidence-overall-simple.png")
png(file = path, height = "10", width = "9", units = "in", res = 300)
forest(HAP_incidence_overall_meta,
       sortvar = effect_est,
       col.square = "turquoise",
       col.square.lines = "black",
       leftcols = c("studlab", "simplified_n"),
       leftlabs = c("Study", "Sample Size"),
       overall = TRUE,
       hetstat = TRUE,
       print.tau2 = FALSE,
       print.pval.Q = FALSE,
       prediction = TRUE)
dev.off()

HAP_incidence_female <- HAP_incidence %>%
  filter(sex == "Female")

HAP_incidence_female <- HAP_incidence_female[-c(2, 3, 9, 10, 11, 13, 14), ]

HAP_incidence_female_meta <- 
  metagen(data = HAP_incidence_female,
          TE = log(effect_est),
          lower = log(lower.ci),
          upper = log(upper.ci),
          studlab = study_ID,
          sm = "RR",
          common = FALSE,
          random = TRUE,
          method.tau = "DL",
          method.bias = "Egger")
summary(HAP_incidence_female_meta)

path <- paste0(fp, "/forest-plot-HAP-incidence-female.png")
png(file = path, height = "3", width = "13", units = "in", res = 300)
forest(HAP_incidence_female_meta,
       sortvar = effect_est,
       col.square = "turquoise",
       col.square.lines = "black",
       leftcols = c("studlab", "country", "age_group", "sex", "study_design", "simplified_n"),
       leftlabs = c("Study", "Location", "Age Group", "Sex", "Study Design", "Sample Size"),
       overall = TRUE,
       hetstat = TRUE,
       print.tau2 = FALSE,
       print.pval.Q = FALSE,
       prediction = TRUE)
dev.off()

HAP_incidence_children <- HAP_incidence %>%
  filter(age_group == "Children")

HAP_incidence_children_meta <- 
  metagen(data = HAP_incidence_children,
          TE = log(effect_est),
          lower = log(lower.ci),
          upper = log(upper.ci),
          studlab = study_ID,
          sm = "RR",
          common = FALSE,
          random = TRUE,
          method.tau = "DL",
          method.bias = "Egger")
summary(HAP_incidence_children_meta)

path <- paste0(fp, "/forest-plot-HAP-incidence-children.png")
png(file = path, height = "3.5", width = "13", units = "in", res = 300)
forest(HAP_incidence_children_meta,
       sortvar = effect_est,
       col.square = "turquoise",
       col.square.lines = "black",
       leftcols = c("studlab", "country", "age_group", "sex", "study_design", "simplified_n"),
       leftlabs = c("Study", "Location", "Age Group", "Sex", "Study Design", "Sample Size"),
       overall = TRUE,
       hetstat = TRUE,
       print.tau2 = FALSE,
       print.pval.Q = FALSE,
       prediction = TRUE)
dev.off()

HAP_incidence_adults <- HAP_incidence %>%
  filter(age_group == "Adults") %>%
  filter(flag == 1 | is.na(flag))

HAP_incidence_adults_meta <- 
  metagen(data = HAP_incidence_adults,
          TE = log(effect_est),
          lower = log(lower.ci),
          upper = log(upper.ci),
          studlab = study_ID,
          sm = "RR",
          common = FALSE,
          random = TRUE,
          method.tau = "DL",
          method.bias = "Egger")
summary(HAP_incidence_adults_meta)

path <- paste0(fp, "/forest-plot-HAP-incidence-adults.png")
png(file = path, height = "7", width = "13", units = "in", res = 300)
forest(HAP_incidence_adults_meta,
       sortvar = effect_est,
       col.square = "turquoise",
       col.square.lines = "black",
       leftcols = c("studlab", "country", "age_group", "sex", "study_design", "simplified_n"),
       leftlabs = c("Study", "Location", "Age Group", "Sex", "Study Design", "Sample Size"),
       overall = TRUE,
       hetstat = TRUE,
       print.tau2 = FALSE,
       print.pval.Q = FALSE,
       prediction = TRUE)
dev.off()

HAP_incidence_allages <- HAP_incidence %>%
  filter(age_group == "All") %>%
  filter(flag == 1 | is.na(flag))

HAP_incidence_allages_meta <- 
  metagen(data = HAP_incidence_allages,
          TE = log(effect_est),
          lower = log(lower.ci),
          upper = log(upper.ci),
          studlab = study_ID,
          sm = "RR",
          common = FALSE,
          random = TRUE,
          method.tau = "DL",
          method.bias = "Egger")
summary(HAP_incidence_allages_meta)

path <- paste0(fp, "/forest-plot-HAP-incidence-allages.png")
png(file = path, height = "3.5", width = "13", units = "in", res = 300)
forest(HAP_incidence_allages_meta,
       sortvar = effect_est,
       col.square = "turquoise",
       col.square.lines = "black",
       leftcols = c("studlab", "country", "age_group", "sex", "study_design", "simplified_n"),
       leftlabs = c("Study", "Location", "Age Group", "Sex", "Study Design", "Sample Size"),
       overall = TRUE,
       hetstat = TRUE,
       print.tau2 = FALSE,
       print.pval.Q = FALSE,
       prediction = TRUE)
dev.off()

HAP_incidence_minconf <- HAP_incidence %>%
  filter(min_confounder_adjustment == "Yes") %>%
  filter(flag == 1 | is.na(flag))

HAP_incidence_minconf_meta <- 
  metagen(data = HAP_incidence_minconf,
          TE = log(effect_est),
          lower = log(lower.ci),
          upper = log(upper.ci),
          studlab = study_ID,
          sm = "RR",
          common = FALSE,
          random = TRUE,
          method.tau = "DL",
          method.bias = "Egger")
summary(HAP_incidence_minconf_meta)

path <- paste0(fp, "/forest-plot-HAP-incidence-minconf.png")
png(file = path, height = "4", width = "13", units = "in", res = 300)
forest(HAP_incidence_minconf_meta,
       sortvar = effect_est,
       col.square = "turquoise",
       col.square.lines = "black",
       leftcols = c("studlab", "country", "age_group", "sex", "study_design", "simplified_n"),
       leftlabs = c("Study", "Location", "Age Group", "Sex", "Study Design", "Sample Size"),
       overall = TRUE,
       hetstat = TRUE,
       print.tau2 = FALSE,
       print.pval.Q = FALSE,
       predicton = TRUE)
dev.off()

# Level 1
HAP_incidence_level1 <- HAP_incidence %>%
  filter(diagnosis_level == 1) %>%
  filter(flag == 1 | is.na(flag))

HAP_incidence_level1_meta <- 
  metagen(data = HAP_incidence_level1,
          TE = log(effect_est),
          lower = log(lower.ci),
          upper = log(upper.ci),
          studlab = study_ID,
          sm = "RR",
          common = FALSE,
          random = TRUE,
          method.tau = "DL",
          method.bias = "Egger")
summary(HAP_incidence_level1_meta)

path <- paste0(fp, "/forest-plot-HAP-incidence-level1.png")
png(file = path, height = "3", width = "13", units = "in", res = 300)
forest(HAP_incidence_level1_meta,
       sortvar = effect_est,
       col.square = "turquoise",
       col.square.lines = "black",
       leftcols = c("studlab", "country", "age_group", "sex", "study_design", "simplified_n"),
       leftlabs = c("Study", "Location", "Age Group", "Sex", "Study Design", "Sample Size"),
       overall = TRUE,
       hetstat = TRUE,
       print.tau2 = FALSE,
       print.pval.Q = FALSE,
       predicton = TRUE)
dev.off()

# Level 2
HAP_incidence_level2 <- HAP_incidence %>%
  filter(diagnosis_level == 2) %>%
  filter(flag == 1 | is.na(flag))

HAP_incidence_level2_meta <- 
  metagen(data = HAP_incidence_level2,
          TE = log(effect_est),
          lower = log(lower.ci),
          upper = log(upper.ci),
          studlab = study_ID,
          sm = "RR",
          common = FALSE,
          random = TRUE,
          method.tau = "DL",
          method.bias = "Egger")
summary(HAP_incidence_level2_meta)

path <- paste0(fp, "/forest-plot-HAP-incidence-level2.png")
png(file = path, height = "5", width = "13", units = "in", res = 300)
forest(HAP_incidence_level2_meta,
       sortvar = effect_est,
       col.square = "turquoise",
       col.square.lines = "black",
       leftcols = c("studlab", "country", "age_group", "sex", "study_design", "simplified_n"),
       leftlabs = c("Study", "Location", "Age Group", "Sex", "Study Design", "Sample Size"),
       overall = TRUE,
       hetstat = TRUE,
       print.tau2 = FALSE,
       print.pval.Q = FALSE,
       predicton = TRUE)
dev.off()

# Level 3
HAP_incidence_level3 <- HAP_incidence %>%
  filter(diagnosis_level == 3) %>%
  filter(flag == 1 | is.na(flag))

HAP_incidence_level3_meta <- 
  metagen(data = HAP_incidence_level3,
          TE = log(effect_est),
          lower = log(lower.ci),
          upper = log(upper.ci),
          studlab = study_ID,
          sm = "RR",
          common = FALSE,
          random = TRUE,
          method.tau = "DL",
          method.bias = "Egger")
summary(HAP_incidence_level3_meta)

path <- paste0(fp, "/forest-plot-HAP-incidence-level3.png")
png(file = path, height = "5", width = "13", units = "in", res = 300)
forest(HAP_incidence_level3_meta,
       sortvar = effect_est,
       col.square = "turquoise",
       col.square.lines = "black",
       leftcols = c("studlab", "country", "age_group", "sex", "study_design", "simplified_n"),
       leftlabs = c("Study", "Location", "Age Group", "Sex", "Study Design", "Sample Size"),
       overall = TRUE,
       hetstat = TRUE,
       print.tau2 = FALSE,
       print.pval.Q = FALSE,
       predicton = TRUE)
dev.off()

# LMIC
HAP_incidence_lmic <- HAP_incidence %>%
  filter(lmic == "yes") %>%
  filter(flag == 1 | is.na(flag))

HAP_incidence_lmic_meta <- 
  metagen(data = HAP_incidence_lmic,
          TE = log(effect_est),
          lower = log(lower.ci),
          upper = log(upper.ci),
          studlab = study_ID,
          sm = "RR",
          common = FALSE,
          random = TRUE,
          method.tau = "DL",
          method.bias = "Egger")
summary(HAP_incidence_lmic_meta)

path <- paste0(fp, "/forest-plot-HAP-incidence-lmic.png")
png(file = path, height = "7", width = "13", units = "in", res = 300)
forest(HAP_incidence_lmic_meta,
       sortvar = effect_est,
       col.square = "turquoise",
       col.square.lines = "black",
       leftcols = c("studlab", "country", "age_group", "sex", "study_design", "simplified_n"),
       leftlabs = c("Study", "Location", "Age Group", "Sex", "Study Design", "Sample Size"),
       overall = TRUE,
       hetstat = TRUE,
       print.tau2 = FALSE,
       print.pval.Q = FALSE,
       predicton = TRUE)
dev.off()

# Case-control studies
HAP_incidence_casecontrol <- HAP_incidence %>%
  filter(study_design == "Case-control") %>%
  filter(flag == 1 | is.na(flag))

HAP_incidence_casecontrol_meta <- 
  metagen(data = HAP_incidence_casecontrol,
          TE = log(effect_est),
          lower = log(lower.ci),
          upper = log(upper.ci),
          studlab = study_ID,
          sm = "RR",
          common = FALSE,
          random = TRUE,
          method.tau = "DL",
          method.bias = "Egger")
summary(HAP_incidence_casecontrol_meta)

path <- paste0(fp, "/forest-plot-HAP-incidence-casecontrol.png")
png(file = path, height = "7", width = "13", units = "in", res = 300)
forest(HAP_incidence_casecontrol_meta,
       sortvar = effect_est,
       col.square = "turquoise",
       col.square.lines = "black",
       leftcols = c("studlab", "country", "age_group", "sex", "study_design", "simplified_n"),
       leftlabs = c("Study", "Location", "Age Group", "Sex", "Study Design", "Sample Size"),
       overall = TRUE,
       hetstat = TRUE,
       print.tau2 = FALSE,
       print.pval.Q = FALSE,
       predicton = TRUE)
dev.off()

# Cross-sectional studies
HAP_incidence_crosssec <- HAP_incidence %>%
  filter(study_design == "Cross-sectional") %>%
  filter(flag == 1 | is.na(flag))

HAP_incidence_crosssec_meta <- 
  metagen(data = HAP_incidence_crosssec,
          TE = log(effect_est),
          lower = log(lower.ci),
          upper = log(upper.ci),
          studlab = study_ID,
          sm = "RR",
          common = FALSE,
          random = TRUE,
          method.tau = "DL",
          method.bias = "Egger")
summary(HAP_incidence_crosssec_meta)

path <- paste0(fp, "/forest-plot-HAP-incidence-crosssec.png")
png(file = path, height = "5", width = "13", units = "in", res = 300)
forest(HAP_incidence_crosssec_meta,
       sortvar = effect_est,
       col.square = "turquoise",
       col.square.lines = "black",
       leftcols = c("studlab", "country", "age_group", "sex", "study_design", "simplified_n"),
       leftlabs = c("Study", "Location", "Age Group", "Sex", "Study Design", "Sample Size"),
       overall = TRUE,
       hetstat = TRUE,
       print.tau2 = FALSE,
       print.pval.Q = FALSE,
       predicton = TRUE)
dev.off()

# Current use
HAP_incidence_current <- HAP_incidence %>%
  filter(current_past_both == "current") %>%
  filter(flag == 1 | is.na(flag))

HAP_incidence_current_meta <- 
  metagen(data = HAP_incidence_current,
          TE = log(effect_est),
          lower = log(lower.ci),
          upper = log(upper.ci),
          studlab = study_ID,
          sm = "RR",
          common = FALSE,
          random = TRUE,
          method.tau = "DL",
          method.bias = "Egger")
summary(HAP_incidence_current_meta)

path <- paste0(fp, "/forest-plot-HAP-incidence-current.png")
png(file = path, height = "8", width = "13", units = "in", res = 300)
forest(HAP_incidence_current_meta,
       sortvar = effect_est,
       col.square = "turquoise",
       col.square.lines = "black",
       leftcols = c("studlab", "country", "age_group", "sex", "study_design", "simplified_n"),
       leftlabs = c("Study", "Location", "Age Group", "Sex", "Study Design", "Sample Size"),
       overall = TRUE,
       hetstat = TRUE,
       print.tau2 = FALSE,
       print.pval.Q = FALSE,
       predicton = TRUE)
dev.off()

# Past use
HAP_incidence_both <- HAP_incidence %>%
  filter(current_past_both == "both") %>%
  filter(flag == 1 | is.na(flag))

HAP_incidence_current_meta <- 
  metagen(data = HAP_incidence_current,
          TE = log(effect_est),
          lower = log(lower.ci),
          upper = log(upper.ci),
          studlab = study_ID,
          sm = "RR",
          common = FALSE,
          random = TRUE,
          method.tau = "DL",
          method.bias = "Egger")
summary(HAP_incidence_current_meta)

path <- paste0(fp, "/forest-plot-HAP-incidence-current.png")
png(file = path, height = "8", width = "13", units = "in", res = 300)
forest(HAP_incidence_current_meta,
       sortvar = effect_est,
       col.square = "turquoise",
       col.square.lines = "black",
       leftcols = c("studlab", "country", "age_group", "sex", "study_design", "simplified_n"),
       leftlabs = c("Study", "Location", "Age Group", "Sex", "Study Design", "Sample Size"),
       overall = TRUE,
       hetstat = TRUE,
       print.tau2 = FALSE,
       print.pval.Q = FALSE,
       predicton = TRUE)
dev.off()

# FUNNEL PLOTS ################################################################
# Define fill colors for contour
col.contour = c("gray75", "gray85", "gray95")

funnel(HAP_incidence_overall_meta, 
       studlab = TRUE,
       contour = c(0.9, 0.95, 0.99),
       col.contour = col.contour)
