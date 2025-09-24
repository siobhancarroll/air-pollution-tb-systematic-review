# Packages
library(rnaturalearth)
library(readODS)
library(tidyverse)
library(tmap)
library(xtable)
library(here)
library(PrettyCols)

# Overall map #################################################################

world <- ne_countries(scale = "medium", returnclass = "sf")

HAP_studies <- read_ods(here("HAP_final_dataset.ods"))
countries_HAP <- HAP_studies %>%
  filter(!duplicated(study_ID)) %>%
  separate_rows(name_long, sep = ", ") %>%
  select(name_long) %>%
  group_by(name_long) %>%
  summarize(freq = n())

AAP_studies <- read_ods(here("AAP_final_dataset.ods"))
countries_AAP <- AAP_studies %>%
  filter(!duplicated(study_ID)) %>%
  separate_rows(name_long, sep = ", ") %>%
  select(name_long) %>%
  group_by(name_long) %>%
  summarize(freq = n())

HAP_countries <- HAP_studies %>% select(study_ID, name_long)
AAP_countries <- AAP_studies %>% select(study_ID, name_long)

studies <- rbind(HAP_countries, AAP_countries)
countries_all <- studies %>%
  filter(!duplicated(study_ID)) %>%
  separate_rows(name_long, sep = ", ") %>%
  select(name_long) %>%
  group_by(name_long) %>%
  summarize(freq = n())
countries_all_countries <- as.data.frame(countries_all$name_long)

#countries <- merge(countries_HAP, countries_AAP, by = "name_long", all = TRUE) 
#countries$freq.x[is.na(countries$freq.x)] <- 0
#countries$freq.y[is.na(countries$freq.y)] <- 0
#countries$freq <- countries$freq.x + countries$freq.y
#countries_countries <- as.data.frame(countries$name_long)

all_studies <- merge(world, countries_all, by = "name_long", all.x = TRUE)
all_studies[50, 1] <- "Cote d'Ivoire" # Change because the French special characters cause problems when merging with the IHME TB data
all_studies$freq[is.na(all_studies$freq)] <- 0

map_labels <- all_studies[(all_studies$freq != 0 & !is.na(all_studies$freq)),]
map_labels_countries <- as.data.frame(map_labels$name_long)

IHME_GBD_2021 <- read.csv(here("IHME-GBD_2021_DATA-fbef2080-1.csv"),
                          encoding = "latin1")
all_studies <- merge(all_studies, IHME_GBD_2021, by = "name_long")
tuberculosis_data_check <- all_studies %>%
  filter(is.na(val))

# Add 1 to the freq column to avoid issues with the log transformation in the map scale

map_labels$freq <- as.numeric(map_labels$freq)
all_studies$freq <- as.numeric(all_studies$freq)
all_studies$val <- as.numeric(all_studies$val)

all_studies$val[is.na(all_studies$val)] <- 0

tmap_mode("plot")
overall_map <- tm_shape(all_studies) + 
  tm_polygons(fill = "val", 
              fill.scale = tm_scale_continuous(values = "brewer.reds"),
              fill.legend = tm_legend(title = "TB incidence per 100,000",
                                 orientation = "landscape")) +
  tm_shape(map_labels) +
  tm_bubbles(size = "freq",
             size.scale = tm_scale_continuous(trans = "sqrt",
                                              ticks = c(1, 5, 25, 50),
                                              values = c(2, 5, 12, 20)),
             fill = "grey25",
             col = "black",
             size.legend = tm_legend(title = "Number of studies",
                                orientation = "landscape")) +
  tm_layout(
    legend.outside = TRUE,
    legend.outside.position = "bottom",
    legend.stack = "horizontal",
    legend.outside.size = 0.1,  
    legend.hist.height = 0.5,   
    legend.title.size = 1,
    legend.text.size = 1
  ) +
  tm_crs("ESRI:54030")
overall_map

tmap_save(overall_map, filename = "overall_map.png", height = 4.5, width = 7.5,
          units = "in")

# Table 2

table2data <- read_excel("Table Data/table2data.xlsx")

table2data <- table2data %>%
  select(Study, 
         'Study design', 
         Region, Population, 
         'Exposure(s) of interest',
         'Outcome(s) of interest')

xtable(table2data, auto = TRUE, longtable = TRUE)
