# Packages
library(tidyverse)
library(readODS)
library(officer)
library(flextable)
library(here)

# Load data
HAP_studies <- read_ods(here("HAP_final_dataset.ods"))
AAP_studies <- read_ods(here("AAP_final_dataset.ods"))

# Set output directory for all tables
if (!dir.exists(here("tables"))) {
  dir.create(here("tables"), recursive = TRUE)
}

output_dir <- here("tables")

# Make data sets that contain the information for the tables
HAP_collapsed <- HAP_studies %>% 
  group_by(study_ID) %>%
  summarise(
    Exposure = str_c(unique(exposure), collapse = ", "),
    Outcome = str_c(unique(primary_outcomes), collapse = ", "),
    `Study Design` = first(study_design),
    Country = first(name_long),
    `Sample Size` = first(simplified_n),
    `Age Group` = first(age_group),
    Sex = first(sex),
    `Variables Adjusted For` = first(adjusted_variables),
    .groups = "drop"
  ) %>%
  rename(Study = study_ID)

AAP_collapsed <- AAP_studies %>% 
  group_by(study_ID) %>%
  summarise(
    Exposure = str_c(unique(Pollutant), collapse = ", "),
    Outcome = str_c(unique(primary_outcomes), collapse = ", "),
    `Study Design` = first(study_design),
    Country = first(name_long),
    `Age Group` = first(age_group),
    Sex = first(Sex),
    `Variables Adjusted For` = first(adjusted_variables),
    .groups = "drop"
  ) %>%
  rename(Study = study_ID)

# Turn these into flextables
ft_HAP <- flextable(HAP_collapsed) %>%
  autofit() %>%
  set_header_labels(
    Study = "Study",
    Exposure = "Exposure(s)",
    Outcome = "Outcome(s)",
    `Study Design` = "Study Design",
    Country = "Country",
    `Sample Size` = "Sample Size",
    `Age Group` = "Age Group",
    Sex = "Sex",
    `Variables Adjusted For` = "Variables Adjusted For"
  ) %>%
  theme_vanilla() %>%
  align(j = c("Study", "Exposure", "Outcome", "Study Design", "Country", 
              "Age Group", "Sex", "Variables Adjusted For"),
        align = "left", part = "all") %>%
  align(j = "Sample Size", align = "right", part = "all")

ft_AAP <- flextable(AAP_collapsed) %>%
  autofit() %>%
  set_header_labels(
    Study = "Study",
    Exposure = "Exposure(s)",
    Outcome = "Outcome(s)",
    `Study Design` = "Study Design",
    Country = "Country",
    `Age Group` = "Age Group",
    Sex = "Sex",
    `Variables Adjusted For` = "Variables Adjusted For"
  ) %>%
  theme_vanilla() %>%
  align(j = c("Study", "Exposure", "Outcome", "Study Design", "Country", 
              "Age Group", "Sex", "Variables Adjusted For"),
        align = "left", part = "all")

# Define landscape orientation
landscape_section <- prop_section(
  page_size = page_size(orient = "landscape")
)

# Make a Word document with the tables
doc <- read_docx() %>%
  body_set_default_section(landscape_section) %>%
  body_add_flextable(ft_HAP) %>%
  body_add_par("", style = "Normal") %>%
  body_add_flextable(ft_AAP)

# Save the document
print(doc, target = paste0(output_dir, "descriptive_studies_tables.docx"))
