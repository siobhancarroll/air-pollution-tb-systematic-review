# Load packages and data
library(tidyverse)
library(readxl)
library(data.table)
library(officer)
library(flextable)
library(here)

# Enforce working directory to the project root
setwd(here())

pafs <- read_csv(here("pafs_to_share.csv"))
tb_counts <- read_csv(here("tb_counts.csv"))

# Set output directory for all tables
if (!dir.exists(here("tables"))) {
  dir.create(here("tables"), recursive = TRUE)
}

output_dir <- here("tables")

# Rename PAF and TB count estimates and UIs to avoid confusion
pafs <- pafs %>%
  rename(paf_estimate = mean,
         paf_lower = lower,
         paf_upper = upper) %>%
  # Filter TB counts to the incidence figures
  filter(type == "yll")

tb_counts <- tb_counts %>%
  rename(count_estimate = val,
         count_lower = lower,
         count_upper = upper)

# Create course grouping in TB counts
tb_counts <- tb_counts %>%
  mutate(grouping = case_when(
    age_group_id %in% c(388, 238, 34, 389) ~ "child",
    sex_id == 1 ~ "male",
    sex_id == 2 ~ "female",
    .default = "other"
  ))

# Perform a merge with PAFs and TB counts using coarse grouping
burden <- inner_join(
  tb_counts,
  pafs %>% select(location_id, grouping, year_id, paf_estimate, paf_lower, 
                  paf_upper, variable),
  by = c("location_id", "grouping", "year_id"),
  relationship = "many-to-many"
)

# Load the location mapping file
region_map <- read_csv(here("locmeta_GBD2023.csv"))
region_map <- region_map %>% 
  select(location_id, parent_id, level, location_name) %>%
  rename(name = location_name) %>%
  mutate(across(c(location_id, parent_id, level), as.integer)) %>%
  as.data.table()

# Copy only necessary columns for the walk
ancestor_dt <- copy(region_map[, .(location_id, parent_id, level, name)]) 

# Initialize ancestor_id with location_id
ancestor_dt[, ancestor_id := location_id]

# Create a named vector for parent lookup
parent_vec <- setNames(region_map$parent_id, region_map$location_id)
level_vec  <- setNames(region_map$level, region_map$location_id)
name_vec   <- setNames(region_map$name, region_map$location_id)

# Walk the hierarchy vectorized
ancestor_vec <- ancestor_dt$ancestor_id
ancestor_level <- ancestor_dt$level

# Keep climbing until all rows are <= level 2 or NA
repeat {
  to_update <- !is.na(ancestor_vec) & ancestor_level > 2
  if (!any(to_update)) break
  
  # Update ancestor IDs in one step
  ancestor_vec[to_update] <- parent_vec[as.character(ancestor_vec[to_update])]
  ancestor_level[to_update] <- level_vec[as.character(ancestor_vec[to_update])]
}

# Assign final level 2 info
ancestor_dt[, level2_id := ifelse(ancestor_level == 2, ancestor_vec, NA_integer_)]
ancestor_dt[, level2_name := ifelse(is.na(level2_id), NA_character_,
                                    name_vec[as.character(level2_id)])]

# Result
level2_map <- ancestor_dt[, .(location_id, level2_id, level2_name)]

# Write function that will calculate attributable burden by year for total PM, HAP, and OAP
calculate_burden_summary <- function(data, paf_variable, n_draws, years) {
  
  # 1. Filter the burden data for the specific PAF variable
  burden_filtered <- data %>%
    filter(variable == paf_variable)
  
  # 2. Estimate the SE for each PAF and TB count estimate from the UIs
  burden_filtered <- burden_filtered %>%
    mutate(se_paf = (paf_upper - paf_lower) / (2 * 1.96),
           se_count = (count_upper - count_lower) / (2 * 1.96))
  
  # 3. Create a list to store summarized burden results
  total_burden_by_year <- list()
  
  # 4. Loop through years to generate draws, calculate burden, and summarize
  for (yr in years) {
    set.seed(123) # Set seed for reproducibility
    burden_year <- burden_filtered %>% filter(year_id == yr)
    n_strata <- nrow(burden_year)
    
    # Generate Monte Carlo draws for counts and PAFs for the year
    count_draws <- matrix(
      rnorm(n_draws * n_strata, mean = burden_year$count_estimate, sd = burden_year$se_count),
      nrow = n_draws, ncol = n_strata
    )
    count_draws[count_draws < 0] <- 0
    
    paf_draws <- matrix(
      rnorm(n_draws * n_strata, mean = burden_year$paf_estimate, sd = burden_year$se_paf),
      nrow = n_draws, ncol = n_strata
    )
    paf_draws[paf_draws < 0] <- 0
    paf_draws[paf_draws > 1] <- 1
    
    # Compute burden draws and sum them
    burden_draws <- count_draws * paf_draws
    total_draws <- rowSums(burden_draws) # Sum across all strata for each draw
    
    # Store the summary for the year
    total_burden_by_year[[as.character(yr)]] <- data.frame(
      year = as.character(yr),
      paf_variable = paf_variable, # Store the variable name
      burden_estimate = mean(total_draws),
      burden_lower = quantile(total_draws, probs = 0.025),
      burden_upper = quantile(total_draws, probs = 0.975)
    )
    
    # Clear large matrices to free memory
    rm(count_draws, paf_draws, burden_draws, total_draws, burden_year)
    gc() 
  }
  
  return(bind_rows(total_burden_by_year))
}

paf_variables_to_run <- c("paf_pm", "paf_hap", "paf_ambient")
years <- c(seq(1990, 2020, by = 5), 2022)

# Run the calculation for each paf variable and store the results in a list
all_burden_summaries_list <- lapply(paf_variables_to_run, function(v) {
  burden_prep <- burden %>%
    mutate(se_paf = (paf_upper - paf_lower) / (2 * 1.96),
           se_count = (count_upper - count_lower) / (2 * 1.96)) %>%
    left_join(level2_map, 
              by = "location_id")
  
  calculate_burden_summary(
    data = burden_prep,
    paf_variable = v,
    n_draws = 1000,
    years = years
  )
})

all_burden_summaries <- bind_rows(all_burden_summaries_list)

# Get the total mortality count for each year
count_draw_matrices_by_year_all <- list()
n_draws = 1000

burden_for_count <- burden %>% filter(variable == "paf_pm")

for (yr in years) {
  set.seed(123)
  burden_year <- burden_for_count %>% filter(year_id == yr)
  n_strata <- nrow(burden_year)
  
  burden_year <- burden_year %>%
    mutate(se_count = (count_upper - count_lower) / (2 * 1.96))
  
  count_draws <- matrix(
    rnorm(n_draws * n_strata, mean = burden_year$count_estimate, sd = burden_year$se_count),
    nrow = n_draws, ncol = n_strata
  )
  count_draws[count_draws < 0] <- 0
  
  count_draw_matrices_by_year_all[[as.character(yr)]] <- count_draws
}

total_count_by_year <- lapply(names(count_draw_matrices_by_year_all), function(yr) {
  draws <- count_draw_matrices_by_year_all[[yr]]
  
  # Sum across all strata for each draw
  total_draws <- rowSums(draws)
  
  # Summarize
  data.frame(
    year = as.character(yr),
    count_estimate = mean(total_draws),
    count_lower = quantile(total_draws, probs = 0.025),
    count_upper = quantile (total_draws, probs = 0.975)
  )
}) %>% bind_rows()

# Clean up columns for presentation
burden_formatted <- all_burden_summaries %>%
  mutate(
    burden_presentation = paste0(
      round(burden_estimate),
      " (",
      round(burden_lower),
      ", ",
      round(burden_upper),
      ")"
    )
  ) %>%
  select(
    Year = year,
    paf_variable,
    burden_presentation
  )

burden_wide <- burden_formatted %>%
  pivot_wider(names_from = "paf_variable",
              values_from = "burden_presentation")

count_formatted <- total_count_by_year %>%
  mutate(
    'Total Mortality Estimate (95% UI)' = paste0(
      round(count_estimate),
      " (",
      round(count_lower),
      ", ",
      round(count_upper),
      ")"
    )
  ) %>%
  select(
    Year = year,
    'Total Mortality Estimate (95% UI)'
  )

# Combine the results
final_table_data <- burden_wide %>%
  left_join(count_formatted, by = "Year")

# Rename columns for final presentation
final_table_data <- final_table_data %>%
  rename(
    'Total Attributable Burden (95% UI)' = paf_pm,
    'Household Attributable Burden (95% UI)' = paf_hap,
    'Outdoor Attributable Burden (95% UI)' = paf_ambient
  )

# REGIONAL ESTIMATES ##########################################################
# Initialize a list to store aggregated results
aggregated_burden_by_year_region <- list()

for (yr in years) {
  # Get the data for the current year
  burden_year <- burden %>% filter(year_id == yr)
  
  # Get the burden draws matrix for the current year
  burden_draws <- draw_matrices_by_year[[as.character(yr)]]
  
  # Extract the region grouping 
  region_grouping <- burden_year %>%
    select(level2_id, level2_name) %>%
    distinct() 
  
  # Create a vector to map each column index (stratum) to its region ID
  region_map_vector <- burden_year$level2_id
  
  # Aggregate draws by summing up the burden draws for all strata within each region
  
  # 1. Convert the burden draws matrix to a data.table for efficient grouping
  draws_dt <- as.data.table(t(burden_draws)) # Transpose so draws are columns
  
  # 2. Add the region ID as the first column
  draws_dt[, level2_id := region_map_vector]
  
  # 3. Group by region and sum the draws (columns V1 to V1000)
  # The columns V1...Vn represent the n_draws (1000)
  agg_draws_dt <- draws_dt[, lapply(.SD, sum), by = level2_id, .SDcols = paste0("V", 1:n_draws)]
  
  # 4. Convert back to a matrix, removing the region ID column
  aggregated_draws_matrix <- as.matrix(agg_draws_dt[, -c("level2_id")])
  
  # 5. Summarize the draws for each region using data.table
  agg_draws_dt[, c("estimate", "lower", "upper") := list(
    rowMeans(.SD),
    apply(.SD, 1, quantile, probs = 0.025),
    apply(.SD, 1, quantile, probs = 0.975)
  ), .SDcols = paste0("V", 1:n_draws)]
  
  region_summary <- agg_draws_dt %>%
    select(level2_id, estimate, lower, upper) %>%
    # Use 'yr' from the loop directly
    mutate(year = yr) %>% 
    # Merge back the region name
    left_join(distinct(burden_year %>% select(level2_id, level2_name)), 
              by = "level2_id")
  
  # Store results in the list
  aggregated_burden_by_year_region[[as.character(yr)]] <- region_summary
}

# Combine all year-region summaries into one data frame
aggregated_burden_df <- bind_rows(aggregated_burden_by_year_region)

# Final cleanup and presentation
aggregated_burden_clean <- aggregated_burden_df %>%
  mutate(
    estimate = round(estimate),
    lower    = round(lower),
    upper    = round(upper)
  ) %>%
  rename(
    Region = level2_name,
    Year = year,
    Estimate = estimate,
    `Lower 95% UI` = lower,
    `Upper 95% UI` = upper
  ) %>%
  arrange(Year, desc(Estimate)) %>%
  # Select and order final columns
  select(Region, Year, Estimate, `Lower 95% UI`, `Upper 95% UI`)

# ESTIMATES BY SEX ############################################################
sex_totals_by_year <- list()

for (yr in years) {
  burden_year <- burden %>% filter(year_id == yr)
  draws <- draw_matrices_by_year[[as.character(yr)]]
  
  # Map each column to its sex
  sex_vec <- burden_year$sex  # length = ncol(draws)
  unique_sexes <- unique(sex_vec)
  
  # Initialize storage for this year
  year_summary <- data.frame(category = character(), estimate = numeric(),
                             lower = numeric(), upper = numeric(), stringsAsFactors = FALSE)
  
  # Total burden
  total_draws <- rowSums(draws)
  year_summary <- rbind(year_summary, data.frame(
    category = "Total",
    estimate = mean(total_draws),
    lower = quantile(total_draws, 0.025),
    upper = quantile(total_draws, 0.975)
  ))
  
  # Sex-specific totals
  for (sx in unique_sexes) {
    cols <- which(sex_vec == sx)
    sex_draws <- rowSums(draws[, cols, drop = FALSE])
    year_summary <- rbind(year_summary, data.frame(
      category = sx,
      estimate = mean(sex_draws),
      lower = quantile(sex_draws, 0.025),
      upper = quantile(sex_draws, 0.975)
    ))
  }
  
  year_summary$year <- yr
  sex_totals_by_year[[as.character(yr)]] <- year_summary
}

# Combine all years
final_sex_table <- bind_rows(sex_totals_by_year) %>%
  select(year, category, estimate, lower, upper) %>%
  mutate(across(c(estimate, lower, upper), round))

# Identify categories to ensure correct column names later
# Filter out "other" if it shouldn't be included in the final table
sex_categories <- c("Female", "Male")
ordered_categories <- c(sex_categories, "Total")

# Prepare the data
sex_table <- final_sex_table %>%
  mutate(
    category = factor(category, levels = ordered_categories),
    burden_with_ui = paste0(estimate, " (", lower, "–", upper, ")")
  ) %>%
  select(year, category, burden_with_ui) %>%
  pivot_wider(
    names_from = category,
    values_from = burden_with_ui
  ) %>%
  select(year, all_of(ordered_categories))

# Quick checks
print(colnames(sex_table))
print(head(sex_table))

# ESTIMATES BY AGE ############################################################

# Initialize a list to store aggregated results by age group
aggregated_burden_by_year_age <- list()

for (yr in years) {
  # Filter burden data for this year
  burden_year <- burden %>% filter(year_id == yr)
  
  # Get the burden draws matrix
  burden_draws <- draw_matrices_by_year[[as.character(yr)]]
  
  # Extract age group info
  age_group_vector <- burden_year$age_group_name  # length = ncol(draws)
  
  # Convert the draws matrix to a data.table (transpose so draws are columns)
  draws_dt <- as.data.table(t(burden_draws))
  
  # Add the age group as the first column
  draws_dt[, age_group_name := age_group_vector]
  
  # Aggregate draws by summing over each age group
  agg_draws_dt <- draws_dt[, lapply(.SD, sum), by = age_group_name, .SDcols = paste0("V", 1:n_draws)]
  
  # Summarize draws for each age group
  age_summary <- agg_draws_dt %>%
    rowwise() %>%
    mutate(
      year = yr,
      estimate = mean(c_across(starts_with("V"))),
      lower    = quantile(c_across(starts_with("V")), probs = 0.025),
      upper    = quantile(c_across(starts_with("V")), probs = 0.975)
    ) %>%
    select(age_group_name, year, estimate, lower, upper) %>%
    ungroup()
  
  # Store in the list
  aggregated_burden_by_year_age[[as.character(yr)]] <- age_summary
}

# Combine all years
aggregated_burden_age_df <- bind_rows(aggregated_burden_by_year_age) %>%
  mutate(
    estimate = round(estimate),
    lower    = round(lower),
    upper    = round(upper)
  ) %>%
  rename(
    `Age Group` = age_group_name,
    Year = year,
    Estimate = estimate,
    `Lower 95% UI` = lower,
    `Upper 95% UI` = upper
  ) %>%
  mutate(`Age Group` = factor(`Age Group`,
                              levels = c("1-5 months",
                                         "6-11 months",
                                         "12 to 23 months",
                                         "2 to 4",
                                         "5 to 9",
                                         "10 to 14",
                                         "15 to 19",
                                         "20 to 24", 
                                         "25 to 29",
                                         "30 to 34",
                                         "35 to 39",
                                         "40 to 44", 
                                         "45 to 49",
                                         "50 to 54",
                                         "55 to 59",
                                         "60 to 64",
                                         "65 to 69", 
                                         "70 to 74",
                                         "75 to 79",
                                         "80 to 84",
                                         "85 to 89", 
                                         "90 to 94",
                                         "95 plus"),
                              ordered = TRUE))

# CREATE WORD TABLES FOR PUBLICATION ##########################################

# Create a new Word document
doc <- read_docx()

# -------------------------------
# Total TB burden table by year
# -------------------------------
doc <- doc %>%
  body_add_par("Table E1. Total TB deaths attributable to PM2.5 by year", style = "heading 2") %>%
  body_add_flextable(
    flextable(final_table_data) %>%
      autofit() %>%
      align(align = "left", part = "all")
  ) %>%
  body_add_par("", style = "Normal")

# -------------------------------
# Regional burden tables
# -------------------------------
doc <- doc %>%
  body_add_par("Table E2. Regional TB deaths attributable to PM2.5 by year", style = "heading 2")

for (yr in sort(unique(aggregated_burden_clean$Year))) {
  yearly_data <- aggregated_burden_clean %>% filter(Year == yr) %>% select(-Year)
  doc <- doc %>%
    body_add_par(paste0("Year: ", yr), style = "heading 3") %>%
    body_add_flextable(
      flextable(yearly_data) %>%
        autofit() %>%
        align(align = "left", part = "all")
    ) %>%
    body_add_par("", style = "Normal")
}

# -------------------------------
# Sex-specific burden table
# -------------------------------
doc <- doc %>%
  body_add_par("Table E3. Sex-specific TB deaths attributable to PM2.5 by year", style = "heading 2") %>%
  body_add_flextable(
    flextable(sex_table) %>%
      autofit() %>%
      align(align = "left", part = "all")
  ) %>%
  body_add_par("", style = "Normal")

# -------------------------------
# Age-specific burden tables
# -------------------------------
doc <- doc %>%
  body_add_par("Table E4. Age-specific TB deaths attributable to PM2.5 by year", style = "heading 2")

for (yr in sort(unique(aggregated_burden_age_df$Year))) {
  yearly_data <- aggregated_burden_age_df %>% filter(Year == yr) %>% select(-Year)
  
  # Append total row
  correct_total_row <- total_by_year_clean %>%
    filter(Year == yr) %>%
    select(Year, Estimate, `Lower 95% UI`, `Upper 95% UI`) %>%
    mutate(`Age Group` = "Total") %>%
    select(`Age Group`, Estimate, `Lower 95% UI`, `Upper 95% UI`)
  
  yearly_data <- bind_rows(yearly_data, correct_total_row)
  
  doc <- doc %>%
    body_add_par(paste0("Year: ", yr), style = "heading 3") %>%
    body_add_flextable(
      flextable(yearly_data) %>%
        autofit() %>%
        align(align = "left", part = "all")
    ) %>%
    body_add_par("", style = "Normal")
}

# -------------------------------
# Save Word document
# -------------------------------
output_file <- file.path(output_dir, "tb_burden_tables.docx")
print(doc, target = output_file)