# Packages
library(rnaturalearth)
library(readODS)
library(tidyverse)
library(tmap)
library(xtable)
library(PrettyCols)

# Set output directory for all tables
if (!dir.exists(figures)) {
  dir.create(figures, recursive = TRUE)
}

output_dir <- here("figures")

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

tmap_save(overall_map, filename = paste0(output_dir, "overall_map.png"), height = 4.5, width = 7.5,
          units = "in")

# HAP map #####################################################################

countries_hap <- HAP_countries %>%
  filter(!duplicated(study_ID)) %>%
  separate_rows(name_long, sep = ", ") %>%
  select(name_long) %>%
  group_by(name_long) %>%
  summarize(freq = n())
countries_hap_countries <- as.data.frame(countries_hap$name_long)

hap_studies <- merge(world, countries_hap, by = "name_long", all.x = TRUE)
hap_studies[50, 1] <- "Cote d'Ivoire" # Change because the French special characters cause problems when merging with the IHME TB data
hap_studies$freq[is.na(hap_studies$freq)] <- 0

map_labels <- hap_studies[(hap_studies$freq != 0 & !is.na(hap_studies$freq)),]
map_labels_countries <- as.data.frame(map_labels$name_long)

hap_studies <- merge(hap_studies, IHME_GBD_2021, by = "name_long")

map_labels$freq <- as.numeric(map_labels$freq)
hap_studies$freq <- as.numeric(hap_studies$freq)
hap_studies$val <- as.numeric(hap_studies$val)

hap_studies$val[is.na(hap_studies$val)] <- 0

tmap_mode("plot")
hap_map <- tm_shape(hap_studies) + 
  tm_polygons(fill = "val", 
              fill.scale = tm_scale_continuous(values = "brewer.reds"),
              fill.legend = tm_legend(title = "TB incidence per 100,000",
                                      orientation = "landscape")) +
  tm_shape(map_labels) +
  tm_bubbles(size = "freq",
             size.scale = tm_scale_continuous(trans = "sqrt",
                                              ticks = c(1, 3, 5, 20),
                                              values = c(1, 5, 15, 20)),
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
hap_map

tmap_save(hap_map, filename = paste0(output_dir, "hap_map.png"), height = 4.5, width = 7.5,
          units = "in")

# AAP map #####################################################################

countries_aap <- AAP_countries %>%
  filter(!duplicated(study_ID)) %>%
  separate_rows(name_long, sep = ", ") %>%
  select(name_long) %>%
  group_by(name_long) %>%
  summarize(freq = n())
countries_aap_countries <- as.data.frame(countries_aap$name_long)

aap_studies <- merge(world, countries_aap, by = "name_long", all.x = TRUE)
aap_studies[50, 1] <- "Cote d'Ivoire" # Change because the French special characters cause problems when merging with the IHME TB data
aap_studies$freq[is.na(aap_studies$freq)] <- 0

map_labels <- aap_studies[(aap_studies$freq != 0 & !is.na(aap_studies$freq)),]
map_labels_countries <- as.data.frame(map_labels$name_long)

aap_studies <- merge(aap_studies, IHME_GBD_2021, by = "name_long")

map_labels$freq <- as.numeric(map_labels$freq)
aap_studies$freq <- as.numeric(aap_studies$freq)
aap_studies$val <- as.numeric(aap_studies$val)

aap_studies$val[is.na(aap_studies$val)] <- 0

tmap_mode("plot")
aap_map <- tm_shape(aap_studies) + 
  tm_polygons(fill = "val", 
              fill.scale = tm_scale_continuous(values = "brewer.reds"),
              fill.legend = tm_legend(title = "TB incidence per 100,000",
                                      orientation = "landscape")) +
  tm_shape(map_labels) +
  tm_bubbles(size = "freq",
             size.scale = tm_scale_continuous(trans = "sqrt",
                                              ticks = c(1, 3, 5, 50),
                                              values = c(1, 5, 15, 20)),
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
aap_map

tmap_save(aap_map, filename = paste0(output_dir, "aap_map.png"), height = 4.5, width = 7.5,
          units = "in")
