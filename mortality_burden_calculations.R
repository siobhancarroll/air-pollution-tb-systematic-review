# Load packages and data
library(tidyverse)
library(readxl)
library(data.table)
library(officer)
library(flextable)
library(here)
library(matrixStats)

set.seed(123)

# Set output directory for all tables
if (!dir.exists(here("tables"))) {
  dir.create(here("tables"), recursive = TRUE)
}
output_dir <- here("tables")
doc <- read_docx()

# Load and Clean PAFs
pafs <- read_csv("pafs_to_share_v3.csv") %>%
  # Remove exact row duplicates
  distinct() %>%
  # Ensure only one estimate per unique strata
  group_by(location_id, year_id, grouping, variable, type) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  rename(paf_estimate = mean,
         paf_lower = lower,
         paf_upper = upper) %>%
  filter(type == "yll")

# Load and Clean TB Counts
tb_counts <- read_csv("tb_counts_v2.csv") %>%
  # Remove exact row duplicates
  distinct() %>%
  # Ensure only one estimate per unique strata
  group_by(location_id, year, age_id, sex_id, measure_name) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  filter(measure_name == "Deaths") %>%
  rename(count_estimate = val,
         count_lower = lower,
         count_upper = upper)

# Create course grouping in TB counts
tb_counts <- tb_counts %>%
  mutate(grouping = case_when(
    age_id == 39 ~ "child",
    age_id %in% c(24, 41, 234) & sex_id == 1 ~ "male",
    age_id %in% c(24, 41, 234) & sex_id == 2 ~ "female",
    .default = "other"
  )) %>%
  mutate(year_id = year)

# Perform a merge with PAFs and TB counts using coarse grouping
burden <- inner_join(
  tb_counts,
  pafs %>% select(location_id, grouping, year_id, paf_estimate, paf_lower, 
                  paf_upper, variable),
  by = c("location_id", "grouping", "year_id"),
  relationship = "many-to-many"
)

# Load the location mapping file
region_map <- read_csv("locmeta_GBD2023.csv")

region_map <- region_map %>% 
  select(location_id, parent_id, level, location_name) %>%
  rename(name = location_name) %>%
  mutate(across(c(location_id, parent_id, level), as.integer)) %>%
  as.data.table()

# Copy only necessary columns for the walk
ancestor_dt_level2 <- copy(region_map[, .(location_id, parent_id, level, name)]) 

# Initialize ancestor_id with location_id
ancestor_dt_level2[, ancestor_id := location_id]

# Create a named vector for parent lookup
parent_vec <- setNames(region_map$parent_id, region_map$location_id)
level_vec  <- setNames(region_map$level, region_map$location_id)
name_vec   <- setNames(region_map$name, region_map$location_id)

# --- FIND LEVEL 2 ANCESTOR ---
ancestor_vec_level2 <- ancestor_dt_level2$ancestor_id
ancestor_level_level2 <- ancestor_dt_level2$level

# Keep climbing until all rows are <= level 2 or NA
repeat {
  to_update <- !is.na(ancestor_vec_level2) & ancestor_level_level2 > 2
  if (!any(to_update)) break
  
  # Update ancestor IDs in one step
  ancestor_vec_level2[to_update] <- 
    parent_vec[as.character(ancestor_vec_level2[to_update])]
  ancestor_level_level2[to_update] <- 
    level_vec[as.character(ancestor_vec_level2[to_update])]
}

# Assign final level 2 info
ancestor_dt_level2[, level2_id := ifelse(ancestor_level_level2 == 2, 
                                         ancestor_vec_level2, 
                                         NA_integer_)]
ancestor_dt_level2[, level2_name := ifelse(is.na(level2_id), 
                                           NA_character_,
                                           name_vec[as.character(level2_id)])]

# Final location map
level2_map <- 
  ancestor_dt_level2[, .(location_id, 
                         level2_id, 
                         level2_name)]

# Define burden_prep globally for reuse with SEs and map ---
burden_prep <- burden %>%
  # Calculate SEs on the combined/coarse data
  mutate(se_paf = (paf_upper - paf_lower) / (2 * 1.96),
         se_count = (count_upper - count_lower) / (2 * 1.96)) %>%
  # Join with the location map
  left_join(level2_map, 
            by = "location_id")

# --------------------------------------------------------------------

# GLOBAL FUNCTION
calculate_grouped_summary_tidy <- function(data, paf_variable, n_draws, years, grouping_col = NULL, filter_col = NULL, filter_na_group = TRUE) {
  
  burden_filtered <- data %>%
    filter(variable == paf_variable)
  
  if (!is.null(filter_col) && filter_na_group) {
    burden_filtered <- burden_filtered %>% filter(!is.na(.data[[filter_col]]))
  }
  
  final_results <- list()
  group_symbols <- syms(grouping_col)
  
  for (yr in years) {
    burden_year <- burden_filtered %>% filter(year_id == yr)
    n_strata <- nrow(burden_year)
    
    if (n_strata == 0) next
    
    # Tidy Draw Generation
    draw_idx <- rep(1:n_draws, each = n_strata)
    strata_idx <- rep(1:n_strata, times = n_draws)
    
    # Create base tibble with grouping column(s) and draw index
    tidy_draws <- tibble(draw_id = draw_idx)
    if (!is.null(grouping_col)) {
      for (col in grouping_col) {
        tidy_draws <- tidy_draws %>% mutate(!!sym(col) := burden_year[[col]][strata_idx])
      }
    }
    
    tidy_draws <- tidy_draws %>%
      mutate(
        # Generate random draws for both counts and PAFs
        count_draw = rnorm(n_draws * n_strata, 
                           mean = burden_year$count_estimate[strata_idx], 
                           sd = burden_year$se_count[strata_idx]),
        paf_draw = rnorm(n_draws * n_strata, 
                         mean = burden_year$paf_estimate[strata_idx], 
                         sd = burden_year$se_paf[strata_idx])
      ) %>%
      # Enforce bounds
      mutate(
        count_draw = pmax(0, count_draw),
        paf_draw = pmax(0, pmin(1, paf_draw))
      ) %>%
      # Calculate burden for each draw
      mutate(burden_draw = count_draw * paf_draw)
    
    # Tidy Aggregation (Sum across strata for each draw and group)
    grouped_draw_sums <- tidy_draws %>%
      group_by(!!!group_symbols, draw_id) %>%
      summarise(
        sum_count = sum(count_draw, na.rm = TRUE),
        sum_burden = sum(burden_draw, na.rm = TRUE),
        .groups = "drop_last"
      )
    
    # Summarize the draws to get the final mean and 95% UI for both metrics
    summary_results <- grouped_draw_sums %>%
      reframe(
        year = yr,
        paf_variable = paf_variable,
        
        # Count Summary
        count_estimate = mean(sum_count),
        count_lower = quantile(sum_count, probs = 0.025),
        count_upper = quantile(sum_count, probs = 0.975),
        
        # Burden Summary
        burden_estimate = mean(sum_burden),
        burden_lower = quantile(sum_burden, probs = 0.025),
        burden_upper = quantile(sum_burden, probs = 0.975),
        
        # Fraction Summary
        fraction_mean = mean(sum_burden / sum_count),
        fraction_lower = quantile(sum_burden / sum_count, 0.025, na.rm = TRUE),
        fraction_upper = quantile(sum_burden / sum_count, 0.975, na.rm = TRUE)
      )
    
    final_results[[as.character(yr)]] <- summary_results
    
    rm(tidy_draws, grouped_draw_sums)
    gc()
  }
  
  return(bind_rows(final_results))
}

# --- GLOBAL PARAMETERS ---
paf_variables_to_run <- c("paf_pm", "paf_hap", "paf_ambient")
years <- c(seq(1990, 2020, by = 5), 2022, 2023)
num_draws <- 1000

# -----------------------------------------------------------------------------
# 1. TOTAL BY YEAR (No Grouping Column)
# -----------------------------------------------------------------------------

burden_prep <- burden_prep %>%
  filter(sex_name != "Both",
         age_name != "All ages")

# Calculate burden for all three PAF variables
all_burden_summaries_list <- lapply(paf_variables_to_run, function(v) {
  calculate_grouped_summary_tidy(
    data = burden_prep,
    paf_variable = v,
    n_draws = num_draws,
    years = years,
    grouping_col = NULL, 
    filter_col = NULL,
    filter_na_group = FALSE 
  )
})
all_burden_summaries <- bind_rows(all_burden_summaries_list)

# Extract total counts from the paf_pm run
total_count_by_year <- all_burden_summaries %>%
  filter(paf_variable == "paf_pm") %>%
  select(year, count_estimate, count_lower, count_upper) %>%
  mutate(paf_variable = "count_total") # Dummy variable name for joining

variable_order <- c(
  "Total TB Deaths",
  "PM2.5-Attributable Deaths",
  "Household PM2.5-Attributable Deaths",
  "Outdoor PM2.5-Attributable Deaths"
)

burden_long <- all_burden_summaries %>%
  filter(paf_variable != "count_total") %>%
  mutate(year = as.character(year)) %>%
  mutate(
    Variable = case_when(
      paf_variable == "paf_pm" ~ "PM2.5-Attributable Deaths",
      paf_variable == "paf_hap" ~ "Household PM2.5-Attributable Deaths",
      paf_variable == "paf_ambient" ~ "Outdoor PM2.5-Attributable Deaths",
      TRUE ~ paf_variable
    ),
    Estimate = format(round(burden_estimate), big.mark = ",", trim = TRUE),
    Lower_95UI = format(round(burden_lower), big.mark = ",", trim = TRUE),
    Upper_95UI = format(round(burden_upper), big.mark = ",", trim = TRUE),
    Frac_Est = round(fraction_mean * 100, 1),
    `Estimate (95% UI)` = paste0(Estimate, " (", Lower_95UI, "-", Upper_95UI, 
                                 ")\n", "[", Frac_Est, "% (",
                                 round(fraction_lower * 100, 1), "-",
                                 round(fraction_upper * 100, 1), "%)]"
    )
  ) %>%
  select(Year = year, Variable, `Estimate (95% UI)`)

count_long <- total_count_by_year %>%
  mutate(year = as.character(year)) %>%
  mutate(
    Variable = "Total TB Deaths",
    `Estimate (95% UI)` = paste0(
      format(round(count_estimate), big.mark = ",", trim = TRUE), " (",
      format(round(count_lower), big.mark = ",", trim = TRUE), "-",
      format(round(count_upper), big.mark = ",", trim = TRUE), ")"
    )
  ) %>%
  select(Year = year, Variable, `Estimate (95% UI)`)

final_long_table <- bind_rows(burden_long, count_long) %>%
  mutate(Variable = factor(Variable, levels = variable_order)) %>%
  arrange(Year, Variable) %>%
  mutate(Variable = case_when(
    Variable == "Household PM2.5-Attributable Deaths" ~ "   Household PM2.5-Attributable",
    Variable == "Outdoor PM2.5-Attributable Deaths"   ~ "   Outdoor PM2.5-Attributable",
    TRUE ~ as.character(Variable)
  ))

ft <- flextable(final_long_table) %>%
  merge_v(j = "Year") %>%
  valign(j = "Year", valign = "top", part = "body") %>%
  
  bold(i = ~ Variable == "Total TB Deaths", part = "body") %>%
  fontsize(i = ~ Variable == "Total TB Deaths", part = "body") %>%
  
  bold(i = ~ Variable == "PM2.5-Attributable Deaths", part = "body") %>%
  padding(i = ~ Variable == "PM2.5-Attributable Deaths", j = "Variable", padding.left = 15, part = "body") %>%
  
  padding(i = ~ Variable %in% c("   Household PM2.5-Attributable", "   Outdoor PM2.5-Attributable"),
          j = "Variable", padding.left = 30, part = "body") %>%
  bg(i = ~ Variable %in% c("   Household PM2.5-Attributable", "   Outdoor PM2.5-Attributable"),
     bg = "#F5F5F5", part = "body") %>%
  
  border(i = ~ Variable == "Total TB Deaths",
         border.top = officer::fp_border(color = "#999999", width = 1.5),
         part = "body") %>%
  
  autofit() %>%
  align(align = "left", part = "all")

# doc output for Total By Year 
doc <- doc %>%
  body_add_flextable(ft) %>%
  body_add_par("", style = "Normal")

# -----------------------------------------------------------------------------
# 2. ESTIMATES BY REGION
# -----------------------------------------------------------------------------

calculate_combined_regional_summary_tidy <- function(data, paf_variable, n_draws, years) {
  return(calculate_grouped_summary_tidy(
    data = data,
    paf_variable = paf_variable,
    n_draws = n_draws,
    years = years,
    grouping_col = c("level2_id", "level2_name"),
    filter_col = "level2_id",
    filter_na_group = TRUE
  ) %>% rename(level2_name = `level2_name`)) 
}

years_to_analyze <- c(1990, 1995, 2000, 2005, 2010, 2015, 2020, 2022, 2023) 
paf_variable_name <- "paf_pm"

regional_totals_wide <- calculate_combined_regional_summary_tidy(
  data = burden_prep,
  paf_variable = paf_variable_name,
  n_draws = num_draws,
  years = years_to_analyze
)

# Prepare the Count Table
count_formatted_regional <- regional_totals_wide %>%
  mutate(
    count_estimate = round(count_estimate),
    count_lower = round(count_lower),
    count_upper = round(count_upper)
  ) %>%
  mutate(
    'Total Deaths (95% UI)' = paste0(
      format(count_estimate, big.mark = ",", trim = TRUE),
      " (", format(count_lower, big.mark = ",", trim = TRUE), ", ",
      format(count_upper, big.mark = ",", trim = TRUE), ")"
    )
  ) %>%
  select(Region = level2_name, Year = year, 'Total Deaths (95% UI)')


# Prepare the Burden Table
burden_formatted_regional <- regional_totals_wide %>%
  mutate(Attributable_Estimate_for_Sort = fraction_mean) %>% 
  mutate(
    burden_estimate = round(burden_estimate),
    burden_lower = round(burden_lower),
    burden_upper = round(burden_upper)
  ) %>%
  mutate(
    'PM2.5-Attributable Deaths (95% UI)' = paste0(
      format(burden_estimate, big.mark = ",", trim = TRUE),
      " (", format(burden_lower, big.mark = ",", trim = TRUE), ", ",
      format(burden_upper, big.mark = ",", trim = TRUE), ")"
    )
  ) %>%
  select(Region = level2_name, Year = year, 
         'PM2.5-Attributable Deaths (95% UI)', Attributable_Estimate_for_Sort,
         fraction_mean, fraction_lower, fraction_upper)


# Join, arrange, and finalize the table
final_table_with_mortality <- burden_formatted_regional %>%
  left_join(count_formatted_regional, by = c("Region", "Year")) %>%
  mutate(Fraction = 
           paste0(
             round(fraction_mean * 100, 1), "% (",
             round(fraction_lower * 100, 1), "-",
             round(fraction_upper * 100, 1), "%)"
           )) %>%
  arrange(Year, desc(Attributable_Estimate_for_Sort)) %>%
  select(Region, Year, 'PM2.5-Attributable Deaths (95% UI)', 
         'Total Deaths (95% UI)', Fraction)

# Loop through each year to create and save a separate flextable
years_to_table <- unique(final_table_with_mortality$Year)

for (yr in years_to_table) {
  yearly_data <- final_table_with_mortality %>%
    filter(Year == yr) %>%
    select(-Year)
  
  doc <- doc %>%
    body_add_par(paste0("Year: ", yr), style = "heading 3") %>%
    body_add_flextable(
      flextable(yearly_data) %>%
        set_table_properties(layout = "fixed") %>%
        width(j = "Region", width = 1.5) %>%
        width(j = "PM2.5-Attributable Deaths (95% UI)", width = 1.7) %>%
        width(j = "Total Deaths (95% UI)", width = 1.6) %>%
        width(j = "Fraction", width = 1.7) %>% 
        align(align = "left", part = "all") %>%
        fontsize(size = 9, part = "all")
    ) %>%
    body_add_par("", style = "Normal")
}


# -----------------------------------------------------------------------------
# 3. ESTIMATES BY SEX 
# -----------------------------------------------------------------------------

global_sex_summary_df <- calculate_grouped_summary_tidy(
  data = burden_prep, 
  paf_variable = "paf_pm",
  n_draws = num_draws,
  years = years,
  grouping_col = "sex_name", 
  filter_col = "sex_name",
  filter_na_group = TRUE
) %>% rename(Sex = sex_name)

# Apply formatting and sorting logic
final_sex_data <- global_sex_summary_df %>%
  filter(Sex %in% c("Male", "Female")) %>%
  mutate(Year = as.character(year)) %>%
  mutate(
    burden_estimate_round = round(burden_estimate),
    count_estimate_round = round(count_estimate),
    'PM2.5-Attributable Deaths (95% UI)' = paste0(
      format(burden_estimate_round, big.mark = ",", trim = TRUE),
      " (", format(round(burden_lower), big.mark = ",", trim = TRUE), "-",
      format(round(burden_upper), big.mark = ",", trim = TRUE), ")"
    ),
    'Total Deaths (95% UI)' = paste0(
      format(count_estimate_round, big.mark = ",", trim = TRUE),
      " (", format(round(count_lower), big.mark = ",", trim = TRUE), "-",
      format(round(count_upper), big.mark = ",", trim = TRUE), ")"),
    Fraction = 
             paste0(
               round(fraction_mean * 100, 1), "% (",
               round(fraction_lower * 100, 1), "-",
               round(fraction_upper * 100, 1), "%)"
             )) %>%
  # 3. Final selection and sorting
  select(
    Year,
    Sex,
    'PM2.5-Attributable Deaths (95% UI)',
    'Total Deaths (95% UI)',
    Fraction
  ) %>%
  arrange(Year, Sex) 

# Loop through male/female to create and save a separate flextable
sex_groups <- unique(final_sex_data$Sex)

for (sx in sex_groups) {
  yearly_data_sex <- final_sex_data %>%
    filter(Sex == sx) %>%
    select(-Sex) 
  
  doc <- doc %>%
    body_add_par(paste0("Sex: ", sx), style = "heading 3") %>%
    body_add_flextable(
      flextable(yearly_data_sex) %>%
        set_table_properties(layout = "fixed") %>%
        # Proportional widths adding up to 6.5
        width(j = 1, width = 1.2) %>% 
        width(j = 2, width = 1.8) %>%
        width(j = 3, width = 1.8) %>% 
        width(j = 4, width = 1.7) %>% 
        fontsize(size = 9, part = "all")
    ) %>%
    body_add_par("", style = "Normal")
}

# -----------------------------------------------------------------------------
# 4. ESTIMATES BY AGE
# -----------------------------------------------------------------------------

raw_age_summary_df <- calculate_grouped_summary_tidy(
  data = burden_prep,
  paf_variable = "paf_pm",
  n_draws = num_draws,
  years = years,
  grouping_col = "age_name",
  filter_col = "age_name",
  filter_na_group = TRUE
) %>% rename(`Age Group` = age_name)

# Define the Age Order for the final table
age_order <- c("0-14 years", "15-49 years", "50-74 years", "75+ years")

final_age_data <- raw_age_summary_df %>%
  mutate(Year = as.character(year)) %>%
  mutate(
    burden_estimate_round = round(burden_estimate),
    count_estimate_round = round(count_estimate),
    'PM2.5-Attributable Deaths (95% UI)' = paste0(
      format(burden_estimate_round, big.mark = ",", trim = TRUE),
      " (", format(round(burden_lower), big.mark = ",", trim = TRUE), "-",
      format(round(burden_upper), big.mark = ",", trim = TRUE), ")"
    ),
    'Total Deaths (95% UI)' = paste0(
      format(count_estimate_round, big.mark = ",", trim = TRUE),
      " (", format(round(count_lower), big.mark = ",", trim = TRUE), "-",
      format(round(count_upper), big.mark = ",", trim = TRUE), ")"
    ),
    Fraction = 
      paste0(
        round(fraction_mean * 100, 1), "% (",
        round(fraction_lower * 100, 1), "-",
        round(fraction_upper * 100, 1), "%)"
      )
  ) %>%
  mutate(`Age Group` = factor(`Age Group`, levels = age_order, ordered = TRUE)) %>%
  select(
    Year,
    `Age Group`,
    'PM2.5-Attributable Deaths (95% UI)',
    'Total Deaths (95% UI)',
    Fraction
  ) %>%
  arrange(Year, `Age Group`) 

# Loop through each year to create and save a separate flextable
for (yr in years_to_analyze) {
  yearly_data_age <- final_age_data %>%
    filter(Year == yr) %>%
    select(-Year)
  
  doc <- doc %>%
    body_add_par(paste0("Year: ", yr), style = "heading 3") %>%
    body_add_flextable(
      flextable(yearly_data_age) %>%
        set_table_properties(layout = "fixed") %>%
        # Proportional widths adding up to 6.5
        width(j = 1, width = 1.2) %>% 
        width(j = 2, width = 1.8) %>% 
        width(j = 3, width = 1.8) %>% 
        width(j = 4, width = 1.7) %>% 
        fontsize(size = 9, part = "all") 
    ) %>%
    body_add_par("", style = "Normal")
}

# -------------------------------
# Save Word document with tables
# -------------------------------
output_file <- file.path(output_dir, "mortality_burden_tables.docx")
print(doc, target = output_file)
