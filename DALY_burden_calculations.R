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
tb_counts <- read_csv(here("tb_morbidity.csv"))

# Set output directory for all tables
if (!dir.exists(here("tables"))) {
  dir.create(here("tables"), recursive = TRUE)
}

output_dir <- here("tables")
doc <- read_docx()

# Rename PAF and TB count estimates and UIs to avoid confusion
pafs <- pafs %>%
  rename(paf_estimate = mean,
         paf_lower = lower,
         paf_upper = upper) %>%
  # Filter TB counts to the incidence figures
  filter(type == "yld")

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

# TOTAL BY YEAR ###############################################################
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

burden_long <- all_burden_summaries %>%
  mutate(year = as.character(year)) %>%
  mutate(
    Variable = case_when(
      paf_variable == "paf_pm" ~ "PM2.5-Attributable DALYs",
      paf_variable == "paf_hap" ~ "Household PM2.5-Attributable DALYs",
      paf_variable == "paf_ambient" ~ "Outdoor PM2.5-Attributable DALYs",
      .default = paf_variable
    ),
    Estimate = round(burden_estimate),
    Lower_95UI = round(burden_lower),
    Upper_95UI = round(burden_upper)
  ) %>%
  select(Year = year, Variable, Estimate, Lower_95UI, Upper_95UI)

count_long <- total_count_by_year %>%
  mutate(year = as.character(year)) %>%
  mutate(
    Variable = "Total DALYs",
    Estimate = round(count_estimate),
    Lower_95UI = round(count_lower),
    Upper_95UI = round(count_upper)
  ) %>%
  select(Year = year, Variable, Estimate, Lower_95UI, Upper_95UI)

final_long_table <- bind_rows(burden_long, count_long) %>%
  arrange(Year, Variable) %>%
  mutate(
    'Estimate (95% CI)' = paste0(
      format(Estimate, big.mark = ",", trim = TRUE),
      " (",
      format(Lower_95UI, big.mark = ",", trim = TRUE),
      "–", 
      format(Upper_95UI, big.mark = ",", trim = TRUE),
      ")"
    )
  ) %>%
  select(
    Year,
    Variable,
    'Estimate (95% CI)'
  )

variable_order <- c(
  "Total DALYs",
  "PM2.5-Attributable DALYs",
  "Household PM2.5-Attributable DALYs",
  "Outdoor PM2.5-Attributable DALYs"
)

final_long_table$Variable <- factor(final_long_table$Variable, levels = variable_order)
final_long_table <- arrange(final_long_table, Year, Variable)

doc <- doc %>%
  body_add_par("Table E1. Total TB DALYs attributable to PM2.5 and Total DALYs by year", style = "heading 2") %>%
  body_add_flextable(
    flextable(final_long_table) %>%
      merge_v(j = 1) %>%
      valign(j = 1, valign = "top", part = "body") %>%
      border_remove() %>%
      
      # Add thick border to top/bottom of the whole table
      hline_top(part = "all", border = fp_border(color = "black", width = 1.5)) %>%
      hline_bottom(part = "all", border = fp_border(color = "black", width = 1.5)) %>%
      
      hline(i = ~ Variable == "Outdoor PM2.5-Attributable DALYs", 
            part = "body") %>%
      
      # Optional: Clean up vertical borders between merged cells
      vline(j = 1, border = fp_border(width = 0), part = "body") %>%
      
      autofit() %>%
      align(align = "left", part = "all")
  ) %>%
  body_add_par("", style = "Normal")

# ESTIMATES BY REGION #########################################################
n_draws <- 1000
years <- c(seq(1990, 2020, by = 5), 2022)

burden_prep <- burden %>%
  # Estimate the SE for each PAF and TB count estimate from the UIs
  mutate(se_paf = (paf_upper - paf_lower) / (2 * 1.96),
         se_count = (count_upper - count_lower) / (2 * 1.96)) %>%
  # Merge the map into the burden data
  left_join(level2_map,
            by = "location_id")

calculate_burden_summary_with_region <- function(data, paf_variable, n_draws, years) {
  
  burden_filtered <- data %>%
    filter(variable == paf_variable)
  
  global_burden_by_year <- list()
  regional_burden_by_year <- list()
  
  for (yr in years) {
    set.seed(123)
    burden_year <- burden_filtered %>% filter(year_id == yr)
    n_strata <- nrow(burden_year)
    
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
    
    burden_draws <- count_draws * paf_draws
    
    # --- REGIONAL AGGREGATION ---
    region_map_vector <- burden_year$level2_id
    
    draws_dt <- as.data.table(t(burden_draws))
    draws_dt[, level2_id := region_map_vector]
    
    agg_draws_dt <- draws_dt[, lapply(.SD, sum), by = level2_id, .SDcols = paste0("V", 1:n_draws)]
    
    region_summary <- agg_draws_dt[, c("estimate", "lower", "upper") := list(
      rowMeans(.SD),
      apply(.SD, 1, quantile, probs = 0.025),
      apply(.SD, 1, quantile, probs = 0.975)
    ), .SDcols = paste0("V", 1:n_draws)] %>%
      select(level2_id, estimate, lower, upper) %>%
      mutate(year = yr, paf_variable = paf_variable) %>%
      # Merge back the region name
      left_join(distinct(burden_year %>% select(level2_id, level2_name)),
                by = "level2_id")
    
    regional_burden_by_year[[as.character(yr)]] <- region_summary
    
    # --- GLOBAL AGGREGATION ---
    total_draws <- rowSums(burden_draws)
    
    global_burden_by_year[[as.character(yr)]] <- data.frame(
      year = as.character(yr),
      paf_variable = paf_variable,
      burden_estimate = mean(total_draws),
      burden_lower = quantile(total_draws, probs = 0.025),
      burden_upper = quantile(total_draws, probs = 0.975)
    )
    
    # Memory clean up
    rm(count_draws, paf_draws, burden_draws, total_draws, burden_year, draws_dt, agg_draws_dt)
    gc()
  }
  
  return(list(
    global_summary = bind_rows(global_burden_by_year),
    regional_summary = bind_rows(regional_burden_by_year)
  ))
}

calculate_regional_count_summary <- function(data, n_draws, years) {
  burden_for_count <- data %>%
    filter(variable == "paf_pm")
  
  regional_count_by_year <- list()
  
  for (yr in years) {
    set.seed(123)
    burden_year <- burden_for_count %>% filter(year_id == yr)
    n_strata <- nrow(burden_year)
    
    count_draws <- matrix(
      rnorm(n_draws * n_strata, mean = burden_year$count_estimate, sd = burden_year$se_count),
      nrow = n_draws, ncol = n_strata
    )
    count_draws[count_draws < 0] <- 0
    
    # --- REGIONAL AGGREGATION ---
    region_map_vector <- burden_year$level2_id
    
    draws_dt <- as.data.table(t(count_draws))
    draws_dt[, level2_id := region_map_vector]
    
    agg_draws_dt <- draws_dt[, lapply(.SD, sum), by = level2_id, .SDcols = paste0("V", 1:n_draws)]
    
    region_summary <- agg_draws_dt[, c("count_estimate", "count_lower", "count_upper") := list(
      rowMeans(.SD),
      apply(.SD, 1, quantile, probs = 0.025),
      apply(.SD, 1, quantile, probs = 0.975)
    ), .SDcols = paste0("V", 1:n_draws)] %>%
      select(level2_id, count_estimate, count_lower, count_upper) %>%
      mutate(year = yr) %>%
      left_join(distinct(burden_year %>% select(level2_id, level2_name)),
                by = "level2_id")
    
    regional_count_by_year[[as.character(yr)]] <- region_summary
    
    # Memory Clean-up
    rm(count_draws, burden_year, draws_dt, agg_draws_dt)
    gc()
  }
  
  return(bind_rows(regional_count_by_year))
}

paf_variables_to_run <- c("paf_pm") # Only PM

pm_results_list <- lapply(paf_variables_to_run, function(v) {
  calculate_burden_summary_with_region(
    data = burden_prep, 
    paf_variable = v,
    n_draws = 1000,
    years = years
  )
})

pm_regional_summary <- bind_rows(lapply(pm_results_list, `[[`, "regional_summary"))

regional_count_summary_df <- calculate_regional_count_summary(
  data = burden_prep, 
  n_draws = 1000,
  years = years
)

burden_formatted_regional <- pm_regional_summary %>%
  mutate(year = as.character(year)) %>%
  mutate(estimate = round(estimate),
         lower = round(lower),
         upper = round(upper)) %>%
  mutate(
    'PM2.5-Attributable DALYs (95% UI)' = paste0(
      format(estimate, big.mark = ",", trim = TRUE),
      " (",
      format(lower, big.mark = ",", trim = TRUE),
      ", ",
      format(upper, big.mark = ",", trim = TRUE),
      ")"
    )
  ) %>%
  select(
    Year = year,
    Region = level2_name,
    'PM2.5-Attributable DALYs (95% UI)',
    Attributable_Estimate = estimate
  )

# Prepare Count Table
count_formatted_regional <- regional_count_summary_df %>%
  mutate(year = as.character(year)) %>%
  mutate(count_estimate = round(count_estimate),
         count_lower = round(count_lower),
         count_upper = round(count_upper)) %>%
  mutate(
    'Total DALYs (95% UI)' = paste0(
      format(count_estimate, big.mark = ",", trim = TRUE),
      " (",
      format(count_lower, big.mark = ",", trim = TRUE),
      ", ",
      format(count_upper, big.mark = ",", trim = TRUE),
      ")"
    )
  ) %>%
  select(
    Year = year,
    Region = level2_name,
    'Total DALYs (95% UI)'
  )

final_table_with_mortality <- burden_formatted_regional %>%
  left_join(
    count_formatted_regional,
    by = c("Year", "Region")
  ) %>%
  arrange(Year, desc(Attributable_Estimate)) %>%
  select(
    Region,
    Year,
    'PM2.5-Attributable DALYs (95% UI)',
    'Total DALYs (95% UI)'
  ) 

# Get the list of unique years
years_to_table <- unique(final_table_with_mortality$Year)

# Loop through each year to create and save a separate flextable
for (yr in years_to_table) {
  
  yearly_data <- final_table_with_mortality %>%
    filter(Year == yr) %>%
    select(-Year)
  
  doc <- doc %>%
    body_add_par(paste0("Year: ", yr), style = "heading 3") %>%
    body_add_flextable(
      flextable(yearly_data) %>%
        autofit() %>%
        align(align = "left", part = "all")
    ) %>%
    body_add_par("", style = "Normal")
}

# ESTIMATES BY SEX ############################################################
calculate_global_summary_by_sex <- function(data, paf_variable, n_draws, years) {
  
  burden_filtered <- data %>%
    filter(variable == paf_variable)
  
  summary_by_year_sex <- list()
  
  for (yr in years) {
    set.seed(123)
    burden_year <- burden_filtered %>% filter(year_id == yr)
    n_strata <- nrow(burden_year)
    
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
    
    burden_draws <- count_draws * paf_draws
    
    grouping_map_dt <- burden_year %>%
      select(grouping) %>%
      as.data.table()
    
    burden_draws_dt <- as.data.table(t(burden_draws))
    burden_draws_dt[, grouping := grouping_map_dt$grouping]
    
    agg_burden_dt <- burden_draws_dt[, lapply(.SD, sum), by = "grouping", .SDcols = paste0("V", 1:n_draws)]
    
    burden_summary <- agg_burden_dt[, c("burden_estimate", "burden_lower", "burden_upper") := list(
      rowMeans(.SD),
      apply(.SD, 1, quantile, probs = 0.025),
      apply(.SD, 1, quantile, probs = 0.975)
    ), .SDcols = paste0("V", 1:n_draws)] %>%
      select(grouping, burden_estimate, burden_lower, burden_upper)
    
    count_draws_dt <- as.data.table(t(count_draws))
    count_draws_dt[, grouping := grouping_map_dt$grouping]
    
    agg_count_dt <- count_draws_dt[, lapply(.SD, sum), by = "grouping", .SDcols = paste0("V", 1:n_draws)]
    
    count_summary <- agg_count_dt[, c("count_estimate", "count_lower", "count_upper") := list(
      rowMeans(.SD),
      apply(.SD, 1, quantile, probs = 0.025),
      apply(.SD, 1, quantile, probs = 0.975)
    ), .SDcols = paste0("V", 1:n_draws)] %>%
      select(grouping, count_estimate, count_lower, count_upper)
    
    yearly_summary <- burden_summary %>%
      left_join(count_summary, by = "grouping") %>%
      mutate(Year = yr, paf_variable = paf_variable)
    
    summary_by_year_sex[[as.character(yr)]] <- yearly_summary
    
    # Memory clean up
    rm(count_draws, paf_draws, burden_draws, burden_year, grouping_map_dt,
       burden_draws_dt, count_draws_dt, agg_burden_dt, agg_count_dt)
    gc()
  }
  
  return(bind_rows(summary_by_year_sex))
}

global_sex_summary_df <- calculate_global_summary_by_sex(
  data = burden_prep, 
  paf_variable = "paf_pm",
  n_draws = n_draws,
  years = years
)

final_sex_data <- global_sex_summary_df %>%
  mutate(Year = as.character(Year)) %>%
  mutate(burden_estimate = round(burden_estimate),
         burden_lower = round(burden_lower),
         burden_upper = round(burden_upper)) %>%
  mutate(count_estimate = round(count_estimate),
         count_lower = round(count_lower),
         count_upper = round(count_upper)) %>%
  filter(grouping %in% c("male", "female")) %>% 
  mutate(
    'PM2.5-Attributable DALYs (95% UI)' = paste0(
      format(burden_estimate, big.mark = ",", trim = TRUE),
      " (",
      format(burden_lower, big.mark = ",", trim = TRUE),
      ", ",
      format(burden_upper, big.mark = ",", trim = TRUE),
      ")"
    ),
    'Total DALYs (95% UI)' = paste0(
      format(count_estimate, big.mark = ",", trim = TRUE),
      " (",
      format(count_lower, big.mark = ",", trim = TRUE),
      ", ",
      format(count_upper, big.mark = ",", trim = TRUE),
      ")"
    )
  ) %>%
  select(
    Year,
    Group = grouping,
    'PM2.5-Attributable DALYs (95% UI)',
    'Total DALYs (95% UI)'
  ) %>%
  arrange(Year)

# Define the groups to iterate over
sex_groups <- c("male", "female")

for (sx in sex_groups) {
  yearly_data_sex <- final_sex_data %>%
    filter(Group == sx) %>%
    select(-Group) 
  
  doc <- doc %>%
    body_add_par(paste0(sx), style = "heading 3") %>%
    body_add_flextable(
      flextable(yearly_data_sex) %>%
        autofit() %>%
        align(align = "left", part = "all")
    ) %>%
    body_add_par("", style = "Normal")
}

# ESTIMATES BY AGE ############################################################

calculate_age_burden_summary <- function(data, paf_variable, n_draws, years) {
  
  burden_filtered <- data %>%
    filter(variable == paf_variable)
  
  aggregated_burden_by_year_age <- list()
  
  for (yr in years) {
    set.seed(123)
    burden_year <- burden_filtered %>% filter(year_id == yr)
    n_strata <- nrow(burden_year)
    
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
    burden_draws <- count_draws * paf_draws
    
    # --- Age Group Aggregation ---
    age_group_vector <- burden_year$age_group_name
    
    # Convert burden draws to data.table for efficient aggregation
    draws_dt <- as.data.table(t(burden_draws)) 
    
    # Add the age group as the first column
    draws_dt[, age_group_name := age_group_vector]
    
    # Aggregate draws by summing over each age group
    agg_draws_dt <- draws_dt[, lapply(.SD, sum), by = age_group_name, .SDcols = paste0("V", 1:n_draws)]
    
    # Summarize draws for each age group (using data.table for efficiency)
    age_summary <- agg_draws_dt[, c("estimate", "lower", "upper") := list(
      rowMeans(.SD),
      apply(.SD, 1, quantile, probs = 0.025),
      apply(.SD, 1, quantile, probs = 0.975)
    ), .SDcols = paste0("V", 1:n_draws)] %>%
      select(age_group_name, estimate, lower, upper) %>%
      mutate(year = yr)
    
    # Store results in the list
    aggregated_burden_by_year_age[[as.character(yr)]] <- age_summary
    
    # Memory Clean-up (Crucial)
    rm(count_draws, paf_draws, burden_draws, burden_year, draws_dt, agg_draws_dt)
    gc()
  }
  
  return(bind_rows(aggregated_burden_by_year_age))
}

calculate_age_count_summary <- function(data, n_draws, years) {
  burden_for_count <- data %>%
    filter(variable == "paf_pm") 
  
  age_count_by_year <- list()
  
  for (yr in years) {
    set.seed(123)
    burden_year <- burden_for_count %>% filter(year_id == yr)
    n_strata <- nrow(burden_year)
    
    count_draws <- matrix(
      rnorm(n_draws * n_strata, mean = burden_year$count_estimate, sd = burden_year$se_count),
      nrow = n_draws, ncol = n_strata
    )
    count_draws[count_draws < 0] <- 0
    
    age_group_vector <- burden_year$age_group_name
    
    draws_dt <- as.data.table(t(count_draws))
    draws_dt[, age_group_name := age_group_vector]
    
    agg_draws_dt <- draws_dt[, lapply(.SD, sum), by = age_group_name, .SDcols = paste0("V", 1:n_draws)]
    
    age_summary <- agg_draws_dt[, c("count_estimate", "count_lower", "count_upper") := list(
      rowMeans(.SD),
      apply(.SD, 1, quantile, probs = 0.025),
      apply(.SD, 1, quantile, probs = 0.975)
    ), .SDcols = paste0("V", 1:n_draws)] %>%
      select(age_group_name, count_estimate, count_lower, count_upper) %>%
      mutate(year = yr)
    
    age_count_by_year[[as.character(yr)]] <- age_summary
    
    # Memory Clean-up
    rm(count_draws, burden_year, draws_dt, agg_draws_dt)
    gc()
  }
  
  return(bind_rows(age_count_by_year))
}

# Run the calculation for Total PM burden
raw_age_burden_df <- calculate_age_burden_summary(
  data = burden_prep,
  paf_variable = "paf_pm",
  n_draws = n_draws,
  years = years
)

# Calculate the total number of deaths (TB count) by age group
raw_age_count_df <- calculate_age_count_summary(
  data = burden_prep,
  n_draws = n_draws,
  years = years
)

# Prepare the count table for joining and formatting
count_formatted_age <- raw_age_count_df %>%
  mutate(year = as.character(year)) %>%
  mutate(count_estimate = round(count_estimate),
         count_lower = round(count_lower),
         count_upper = round(count_upper)) %>%
  mutate(
    'Total DALYs (95% UI)' = paste0(
      format(count_estimate, big.mark = ",", trim = TRUE),
      " (",
      format(count_lower, big.mark = ",", trim = TRUE),
      ", ",
      format(count_upper, big.mark = ",", trim = TRUE),
      ")"
    )
  ) %>%
  select(
    Year = year,
    `Age Group` = age_group_name,
    'Total DALYs (95% UI)'
  )

# Final formatting and cleanup
final_age_data <- raw_age_burden_df %>%
  mutate(year = as.character(year)) %>%
  mutate(estimate = round(estimate),
         lower = round(lower),
         upper = round(upper)) %>%
  mutate(
    'PM2.5-Attributable DALYs (95% UI)' = paste0(
      format(estimate, big.mark = ",", trim = TRUE),
      " (",
      format(lower, big.mark = ",", trim = TRUE),
      ", ",
      format(upper, big.mark = ",", trim = TRUE),
      ")"
    )
  ) %>%
  rename(
    `Age Group` = age_group_name,
    Year = year
  ) %>%
  left_join(count_formatted_age,
            by = c("Year", "Age Group")) %>%
  select(`Age Group`, Year, `PM2.5-Attributable DALYs (95% UI)`, `Total DALYs (95% UI)`) %>%
  mutate(`Age Group` = factor(`Age Group`,
                              levels = c("1-5 months", "6-11 months", "12 to 23 months", "2 to 4", "5 to 9", 
                                         "10 to 14", "15 to 19", "20 to 24", "25 to 29", "30 to 34", 
                                         "35 to 39", "40 to 44", "45 to 49", "50 to 54", "55 to 59", 
                                         "60 to 64", "65 to 69", "70 to 74", "75 to 79", "80 to 84", 
                                         "85 to 89", "90 to 94", "95 plus"),
                              ordered = TRUE)) %>%
  arrange(`Age Group`)

# Loop through each year to create and save a separate flextable
for (yr in years_to_table) {
  
  yearly_data_age <- final_age_data %>%
    filter(Year == yr) %>%
    select(-Year)
  
  doc <- doc %>%
    body_add_par(paste0("Year: ", yr), style = "heading 3") %>%
    body_add_flextable(
      flextable(yearly_data_age) %>%
        autofit() %>%
        align(align = "left", part = "all")
    ) %>%
    body_add_par("", style = "Normal")
}

# -------------------------------
# Save Word document with tables
# -------------------------------
output_file <- file.path(output_dir, "DALY_burden_tables.docx")
print(doc, target = output_file)

