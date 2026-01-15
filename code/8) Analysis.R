
###################################################
# Summary Table Combined Stats FOR SPLIT WITHOUT BARRIERS#
###################################################

# Load packages
library(readxl)
library(dplyr)
library(writexl)

# Set folder path
folder_path <- "./data/Suitable Habitats Buffered"

# List all Excel files in the folder
files <- list.files(folder_path, pattern = "\\.xlsx$", full.names = TRUE)

# Helper function to summarise each species file
summarise_species <- function(file) {
  df <- read_excel(file)
  
  # Skip empty or invalid files
  if (nrow(df) == 0 || all(is.na(df$area_m2))) {
    message("⚠️ Skipping empty file: ", basename(file))
    return(NULL)
  }
  
  summary_df <- df %>%
    group_by(species) %>%
    summarise(
      n_patches = n(),
      total_area_km2 = sum(area_m2, na.rm = TRUE) / 1e6,
      mean_patch_km2 = mean(area_m2, na.rm = TRUE) / 1e6,
      median_patch_km2 = median(area_m2, na.rm = TRUE) / 1e6,
      sd_patch_km2 = sd(area_m2, na.rm = TRUE) / 1e6,
      largest_patch_km2 = max(area_m2, na.rm = TRUE) / 1e6,
      perc_largest_patch = (max(area_m2, na.rm = TRUE) / sum(area_m2, na.rm = TRUE)) * 100,
      fragmentation_index = 1 - (max(area_m2, na.rm = TRUE) / sum(area_m2, na.rm = TRUE))
    ) %>%
    ungroup()
  
  # Round all numeric values except n_patches to 2 decimal places
  summary_df <- summary_df %>%
    mutate(across(where(is.numeric) & !matches("n_patches"), ~ round(., 2)))
  
  return(summary_df)
}

# Read and summarise all valid files
summary_list <- lapply(files, summarise_species)
summary_list <- summary_list[!sapply(summary_list, is.null)]  # remove skipped species

# Combine summaries and sort alphabetically
patch_summary <- bind_rows(summary_list) %>%
  arrange(species)

# Save final overview
output_path <- file.path(folder_path, "Patch_Summary_Overview.xlsx")
write_xlsx(patch_summary, output_path)

message("✅ Patch Summary Overview saved to: ", output_path)

###################################################
# Summary Table Combined Stats FOR SPLIT WITH BARRIERS#
###################################################

# Load packages
library(readxl)
library(dplyr)
library(writexl)

# Set folder path
folder_path <- "./data/Suitable Habitats with Barriers"

# List all Excel files in the folder
files <- list.files(folder_path, pattern = "\\.xlsx$", full.names = TRUE)

# Helper function to summarise each species file
summarise_species <- function(file) {
  df <- read_excel(file)
  
  # Skip empty or invalid files
  if (nrow(df) == 0 || all(is.na(df$area_m2))) {
    message("⚠️ Skipping empty file: ", basename(file))
    return(NULL)
  }
  
  summary_df <- df %>%
    group_by(species) %>%
    summarise(
      n_patches = n(),
      total_area_km2 = sum(area_m2, na.rm = TRUE) / 1e6,
      mean_patch_km2 = mean(area_m2, na.rm = TRUE) / 1e6,
      median_patch_km2 = median(area_m2, na.rm = TRUE) / 1e6,
      sd_patch_km2 = sd(area_m2, na.rm = TRUE) / 1e6,
      largest_patch_km2 = max(area_m2, na.rm = TRUE) / 1e6,
      perc_largest_patch = (max(area_m2, na.rm = TRUE) / sum(area_m2, na.rm = TRUE)) * 100,
      fragmentation_index = 1 - (max(area_m2, na.rm = TRUE) / sum(area_m2, na.rm = TRUE))
    ) %>%
    ungroup()
  
  # Round all numeric values except n_patches to 2 decimal places
  summary_df <- summary_df %>%
    mutate(across(where(is.numeric) & !matches("n_patches"), ~ round(., 2)))
  
  return(summary_df)
}

# Read and summarise all valid files
summary_list <- lapply(files, summarise_species)
summary_list <- summary_list[!sapply(summary_list, is.null)]  # remove skipped species

# Combine summaries and sort alphabetically
patch_summary <- bind_rows(summary_list) %>%
  arrange(species)

# Save final overview
output_path <- file.path(folder_path, "Patch_Summary_Overview.xlsx")
write_xlsx(patch_summary, output_path)

message("✅ Patch Summary Overview saved to: ", output_path)

###################################################
# Hotspot Patch Summary: Primary & Secondary
###################################################

# Load packages
library(terra)
library(sf)
library(dplyr)
library(writexl)

# -----------------------
# Paths
# -----------------------
hotspot_dir <- "C:/Users/jx23973/OneDrive - University of Bristol/PhD/Coding/rewilding_in_britain/data/Rewilding Hotspots"

hotspot_files <- c(
  file.path(hotspot_dir, "8_8_Rewilding_Hotspots_Renumbered.tif"),
  file.path(hotspot_dir, "8_7_Rewilding_Hotspots_Renumbered.tif")
)

out_path <- file.path(hotspot_dir, "Hotspot_Patch_Summary.xlsx")

# -----------------------
# Helper function
# -----------------------
summarise_hotspot <- function(raster_file) {
  
  if (!file.exists(raster_file)) {
    warning("⚠️ File not found: ", raster_file)
    return(NULL)
  }
  
  message("📊 Processing: ", basename(raster_file))
  
  r <- rast(raster_file)
  
  # Convert patches to polygons (NA cells ignored automatically)
  patches_vect <- as.polygons(r, dissolve = TRUE)
  
  # Convert to sf for area calculation
  patches_sf <- st_as_sf(patches_vect)
  patches_sf$area_m2 <- as.numeric(st_area(patches_sf))
  
  # Identify hotspot type from filename
  hotspot_type <- if (grepl("8_8", basename(raster_file))) {
    "Primary"
  } else if (grepl("8_7", basename(raster_file))) {
    "Secondary"
  } else {
    "Unknown"
  }
  
  # Patch summary
  summary_df <- patches_sf |>
    st_drop_geometry() |>
    summarise(
      hotspot = hotspot_type,
      n_patches = n(),
      total_area_km2 = sum(area_m2) / 1e6,
      mean_patch_km2 = mean(area_m2) / 1e6,
      median_patch_km2 = median(area_m2) / 1e6,
      sd_patch_km2 = sd(area_m2) / 1e6,
      largest_patch_km2 = max(area_m2) / 1e6,
      perc_largest_patch = (max(area_m2) / sum(area_m2)) * 100,
      fragmentation_index = 1 - (max(area_m2) / sum(area_m2))
    ) |>
    mutate(across(where(is.numeric) & !matches("n_patches"), ~ round(.x, 2)))
  
  return(summary_df)
}

# -----------------------
# Run summaries
# -----------------------
hotspot_summary_list <- lapply(hotspot_files, summarise_hotspot)
hotspot_summary_list <- hotspot_summary_list[!sapply(hotspot_summary_list, is.null)]

hotspot_summary <- bind_rows(hotspot_summary_list)

# -----------------------
# Save Excel
# -----------------------
write_xlsx(hotspot_summary, out_path)

message("✅ Hotspot patch summary saved to: ", out_path)

###################################################
# Identify & Plot Largest Secondary Hotspot Patch
###################################################

library(terra)
library(sf)
library(dplyr)
library(ggplot2)

# -----------------------
# Load secondary hotspot raster
# -----------------------
hotspot_dir <- "C:/Users/jx23973/OneDrive - University of Bristol/PhD/Coding/rewilding_in_britain/data/Rewilding Hotspots"
secondary_file <- file.path(hotspot_dir, "8_7_Rewilding_Hotspots_Renumbered.tif")

if (!file.exists(secondary_file)) stop("⚠️ File not found: ", secondary_file)

r <- rast(secondary_file)

# Convert patches to polygons
patches_vect <- as.polygons(r, dissolve = TRUE)
patches_sf <- st_as_sf(patches_vect)
patches_sf$area_m2 <- as.numeric(st_area(patches_sf))

# Identify the largest patch
largest_patch <- patches_sf %>%
  filter(area_m2 == max(area_m2))

# -----------------------
# Plot largest patch using base R and save as JPEG
# -----------------------
jpeg_filename <- file.path(hotspot_dir, "largest_secondary_patch_base.jpg")

# Open JPEG device
jpeg(filename = jpeg_filename, width = 2000, height = 2000, res = 300)

# Plot all patches in light grey
plot(patches_sf$geometry, col = "grey90", border = "grey50", main = "Largest Secondary Hotspot Patch")

# Overlay the largest patch in green
plot(largest_patch$geometry, col = "forestgreen", border = "black", add = TRUE)

# Close device
dev.off()

cat("✅ JPEG saved to:", jpeg_filename, "\n")
