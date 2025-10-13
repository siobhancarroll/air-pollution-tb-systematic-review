# Load packages and data
library(tidyverse)
library(readxl)
library(data.table)
library(kableExtra)
library(webshot2)
library(chromote)
library(officer)
library(flextable)

pafs <- read_csv("/home/siobhancarroll/Documents/PhD/systematic_review/PAF_data/pafs_to_share.csv")
tb_counts <- read_csv("/home/siobhancarroll/Documents/PhD/systematic_review/PAF_data/tb_counts.csv")

# Set output directory for all tables
output_dir <- "/home/siobhancarroll/Documents/PhD/systematic_review/final_figures_tables/supplemental/"

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
burden_data <- inner_join(
  tb_counts,
  pafs %>% select(location_id, grouping, year_id, paf_estimate, paf_lower, 
                  paf_upper, variable),
  by = c("location_id", "grouping", "year_id")
)

# Final filtering
burden <- burden_data %>%
  filter(variable == "paf_pm")

# UNCERTAINTY CALCULATIONS ####################################################
# Estimate the SE for each PAF and TB count estimate from the UIs
# Assuming a Normal distribution since true draws aren't available
burden <- burden %>%
  mutate(se_paf = (paf_upper - paf_lower) / (2 * 1.96),
         se_count = (count_upper - count_lower) / (2 * 1.96))

# Load the location mapping file
region_map <- read_csv("/home/siobhancarroll/Documents/PhD/systematic_review/locmeta_GBD2023.csv") 
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

# Merge the map into the burden data
burden <- burden %>%
  left_join(level2_map, 
             by = "location_id")

# Generate Monte Carlo draws
n_draws <- 1000
n_strata <- nrow(burden)

# Years of interest
years <- c(seq(1990, 2020, by = 5), 2022)

draw_matrices_by_year <- list()

for (yr in years) {
  set.seed(123)
  burden_year <- burden %>% filter(year_id == yr)
  n_strata <- nrow(burden_year)
  
  # Generate all draws for this year in one go
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
  
  # Compute burden draws
  burden_draws <- count_draws * paf_draws
  
  # Store in the list
  draw_matrices_by_year[[as.character(yr)]] <- burden_draws
}

# Summarize total burden per year
total_by_year <- lapply(names(draw_matrices_by_year), function(yr) {
  draws <- draw_matrices_by_year[[yr]]
  
  # Sum across all strata for each draw
  total_draws <- rowSums(draws)
  
  # Summarize
  data.frame(
    year = as.integer(yr),
    estimate = mean(total_draws),
    lower    = quantile(total_draws, probs = 0.025),
    upper    = quantile(total_draws, probs = 0.975)
  )
}) %>% bind_rows()

total_by_year <- total_by_year %>%
  mutate(
    estimate = round(estimate),
    lower    = round(lower),
    upper    = round(upper)
  )

# Clean up column names for presentation
total_by_year_clean <- total_by_year %>%
  rename(
    Year = year,
    Estimate = estimate,
    `Lower 95% UI` = lower,
    `Upper 95% UI` = upper
  )

# Save total table as PNG
temp_html_file <- tempfile(fileext = ".html")
total_by_year_clean %>%
  kbl(
    caption = "Estimated TB Deaths Attributable to PM2.5 per Year",
    row.names = FALSE,        # suppress row names
    align = "c"               # center-align all columns
  ) %>%
  kable_styling(
    full_width = FALSE,
    position = "center",
    bootstrap_options = c("striped", "hover", "condensed")
  ) %>%
  save_kable(file = temp_html_file)

webshot2::webshot(
  temp_html_file,
  vwidth = 500,
  delay = 0.5,
  zoom = 3,
  file = file.path(output_dir, "total_tb_burden.png")
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

# Get the unique years
unique_years <- sort(unique(aggregated_burden_clean$Year))

# Iterate through each year and save a separate table as PNG
for (yr in unique_years) {
  
  # 1. Filter the data for the current year
  yearly_data <- aggregated_burden_clean %>%
    filter(Year == yr) %>%
    select(-Year)
  
  # 2. Create the kableExtra object and render it to HTML
  table_html_object <- yearly_data %>%
    kbl(
      caption = paste0("Estimated TB Deaths Attributable to PM2.5 by GBD Region (Year: ", yr, ")"),
      row.names = FALSE,
      align = "c"
    ) %>%
    kable_styling(
      full_width = FALSE,
      position = "center",
      bootstrap_options = c("striped", "hover", "condensed")
    )
  
  # 3. Save the HTML object to a temporary file
  # Note: The output from kableExtra is an HTML object, which must be saved before webshot2 can render it
  temp_html_file <- tempfile(fileext = ".html")
  save_kable(table_html_object, file = temp_html_file)
  
  # 4. Define the output file path for the PNG
  file_name <- paste0("regional_burden_tb_", yr, ".png")
  output_path <- file.path(output_dir, file_name)
  
  # 5. Save the table as a PNG file by providing the temporary HTML file path
  webshot2::webshot(
      temp_html_file, 
      vwidth = 600,        
      delay = 0.5,           
      zoom = 3,
      file = output_path 
    )
}

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

# Create the Kable using the correct final column names
col_names_header <- c("Year", ordered_categories)

# Save sex-specific table as PNG
temp_html_file <- tempfile(fileext = ".html")
sex_table %>%
  kbl(
    caption = "Estimated TB Deaths Attributable to PM2.5 by Sex",
    col.names = col_names_header,
    align = "c"
  ) %>%
  add_header_above(c(" " = 1, "Burden Estimate (95% UI)" = 3)) %>%
  kable_styling(full_width = FALSE, position = "center") %>%
  save_kable(file = temp_html_file)

webshot2::webshot(
  temp_html_file,
  vwidth = 800,
  delay = 0.5,
  zoom = 3,
  file = file.path(output_dir, "sex_tb_burden.png")
)

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

# Save separate tables per year
unique_years <- sort(unique(aggregated_burden_age_df$Year))

for (yr in unique_years) {
  yearly_data <- aggregated_burden_age_df %>%
    filter(Year == yr) %>%
    select(-Year)
  
  # Retrieve the pre-calculated total from the total_by_year data frame
  correct_total_row <- total_by_year %>%
    filter(year == yr) %>%
    select(year, Estimate = estimate, `Lower 95% UI` = lower, `Upper 95% UI` = upper) %>%
    mutate(`Age Group` = "Total") %>%
    select(`Age Group`, everything())

  yearly_data <- bind_rows(yearly_data, correct_total_row %>% select(-year))

  # Make sure Total is last
  yearly_data <- yearly_data %>%
    mutate(`Age Group` = factor(`Age Group`, levels = c(levels(aggregated_burden_age_df$`Age Group`), "Total"), ordered = TRUE)) %>%
    arrange(`Age Group`)
  
  # ---- Create and save table ----
  table_html_object <- yearly_data %>%
    kbl(
      caption = paste0("Estimated TB Deaths Attributable to PM2.5 by Age Group (Year: ", yr, ")"),
      row.names = FALSE,
      align = "c"
    ) %>%
    kable_styling(
      full_width = FALSE,
      position = "center",
      bootstrap_options = c("striped", "hover", "condensed")
    )
  
  temp_html_file <- tempfile(fileext = ".html")
  save_kable(table_html_object, file = temp_html_file)
  
  file_name <- paste0("age_group_burden_tb_", yr, ".png")
  output_path <- file.path(output_dir, file_name)
  
  webshot2::webshot(
    temp_html_file,
    vwidth = 600,
    delay = 0.5,
    zoom = 3,
    file = output_path
  )
}

# CREATE WORD TABLES FOR PUBLICATION ##########################################

# Create a new Word document
doc <- read_docx()

# -------------------------------
# Total TB burden table by year
# -------------------------------
doc <- doc %>%
  body_add_par("Table E1. Total TB deaths attributable to PM2.5 by year", style = "heading 2") %>%
  body_add_flextable(
    flextable(total_by_year_clean) %>%
      autofit() %>%
      align(align = "center", part = "all")
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
        align(align = "center", part = "all")
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
      align(align = "center", part = "all")
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
        align(align = "center", part = "all")
    ) %>%
    body_add_par("", style = "Normal")
}

# -------------------------------
# Save Word document
# -------------------------------
output_file <- file.path(output_dir, "tb_burden_tables.docx")
print(doc, target = output_file)

