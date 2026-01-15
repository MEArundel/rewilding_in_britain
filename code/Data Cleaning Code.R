#################################################
## Preamble ##
#################################################

# Clear Memory
rm(list=ls())

# Load Packages
library(dplyr)
library(tidyr)
library(xlsx)

# Setting seed
set.seed(13) #ensures that any randomness is controlled between MAC and Windows and between sessions

# Paths
data_path <- "./data/"
out_dir <- "./data/"

# Read in raw data
home_range_data <- read.csv("./data/HomeRangeData_2025_04_11.csv")

target_species <- c(
  "Castor fiber", "Bison bonasus", "Ursus arctos", "Lynx lynx", "Alces alces",
  "Martes martes", "Sus scrofa", "Felis silvestris", "Canis lupus", "Gulo gulo"
)

# Filter data
home_range_data <- home_range_data %>%
  filter(
    Species %in% target_species,
    Context %in% c("Wild", "Reintroduced"),
    Longitude >= -25, Longitude <= 45,
    Latitude >= 34, Latitude <= 72,
    !grepl("familiaris|italicus", subspecies, ignore.case = TRUE)
  )
  
# Tidy data
home_range_data <- home_range_data %>%
  # Split Species into Genus and Species
  separate(Species, into = c("Genus", "Species"), sep = " ") %>%
  # Add species column right after Species
  mutate(species = case_when(
    Genus == "Castor" & Species == "fiber" ~ "Beaver",
    Genus == "Bison" & Species == "bonasus" ~ "Bison",
    Genus == "Ursus" & Species == "arctos" ~ "Brown Bear",
    Genus == "Lynx" & Species == "lynx" ~ "Lynx",
    Genus == "Alces" & Species == "alces" ~ "Elk",
    Genus == "Martes" & Species == "martes" ~ "Pine Marten",
    Genus == "Sus" & Species == "scrofa" ~ "Wild Boar",
    Genus == "Felis" & Species == "silvestris" ~ "Wildcat",
    Genus == "Canis" & Species == "lupus" ~ "Wolf",
    Genus == "Gulo" & Species == "gulo" ~ "Wolverine",
    TRUE ~ NA_character_
  )) %>%
  relocate(species, .after = Species) %>%
  # Remove unwanted columns
  select(-Study_ID, -subspecies, -Ind_ID, -Sex, -Life_Stage, -Reproductive_Status, 
         -Body_mass_kg, -Locomotion, -dayStart, -monthStart, -dayEnd, -monthEnd)

# Define the desired order for common names
species_order <- c(
  "Beaver", "Bison", "Brown Bear", "Lynx", "Elk",
  "Pine Marten", "Wild Boar", "Wildcat", "Wolf", "Wolverine"
)

# Apply ordering
home_range_data <- home_range_data %>%
  mutate(species = factor(species, levels = species_order)) %>%
  arrange(species)

# Define output path
output_file <- file.path(out_dir, "Supplementary Table 5 - Home Range Values.xlsx")
# Write to Excel
write.xlsx(home_range_data, output_file)

# Calculate mean home range for each species
mean_home_ranges <- home_range_data %>%
  group_by(species) %>%
  summarise(mean_home_range_km2 = mean(Home_Range_km2, na.rm = TRUE)) %>%
  arrange(factor(species, levels = species_order))

# Add mean home range column (in m²)
home_range_data_export <- mean_home_ranges %>%
  group_by(species) %>%
  mutate(mean_home_range_m2 = mean(mean_home_range_km2, na.rm = TRUE) * 1e6) %>%
  ungroup()

# Define output path
output_file <- file.path(out_dir, "Supplementary Table 6 - Mean Home Range Values.xlsx")
# Write to Excel
write.xlsx(home_range_data_export, output_file)

# Calculate mean No_Individuals for species with HR_Level = "Group" or "Group (mean)"
home_range_data_group <- home_range_data %>%
  filter(HR_Level %in% c("Group", "Group (mean)"))%>%
  mutate(No_Individuals = as.integer(No_Individuals)) %>%
  filter(!is.na(No_Individuals))

# Calculate mean No_Individuals for each species 
mean_no_individuals_by_species <- home_range_data_group %>%
  group_by(species) %>%
  summarise(mean_no_individuals = mean(No_Individuals)) %>%
  arrange(factor(species, levels = species_order))

# Define output path
output_file <- file.path(out_dir, "Supplementary Table 7 - Mean Group Size Values.xlsx")
# Write to Excel
write.xlsx(mean_no_individuals_by_species, output_file)
