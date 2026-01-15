#########################################
# Formation of habitat suitability maps #
#########################################
# Description - This script generates binary habitat suitability maps for focal
# species in Britain by reclassifying a national land cover raster using species-
# specific suitability rules. For each species, land cover classes are converted to
# suitable (1) or unsuitable (0) habitat based on a reference table, masked to a
# dissolved England–Scotland–Wales boundary, and exported as both a GeoTIFF and a
# static JPEG map. The script additionally compiles all species-specific suitability
# rasters into a single multi-panel figure for comparative visualisation.

# Clean session (optional)
# rm(list = ls())

# Load packages
library(terra)

# -----------------------
# Settings / I-O paths
# -----------------------
set.seed(13) # included for consistency across scripts

in_dir <- "./data/"
out_dir_data <- "./data/Suitable Habitats/"
out_dir_res <- "./results/Suitable Habitats/"

# Ensure output directories exist
dir.create(out_dir_data, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_res, recursive = TRUE, showWarnings = FALSE)

# -------------------------------
# Inputs and checks
# -------------------------------
landcover_path <- file.path(in_dir, "Land Coverage Map Reclassified Raster.tif")
hsref_path     <- file.path(in_dir, "Species HS Reference Data.csv")
boundary_path  <- file.path(in_dir, "Britain Boundary Polygon.geojson")

stopifnot(file.exists(landcover_path))
stopifnot(file.exists(hsref_path))
stopifnot(file.exists(boundary_path))

# Load raster layer
land_cover_reclassified <- rast(landcover_path)

# Load species HS reference data
species_hs_ref <- read.csv(hsref_path)

# ------------
# Boundary
# ------------
boundary_vect <- vect(boundary_path)

# Project boundary to the raster CRS
boundary_vect <- project(boundary_vect, crs(land_cover_reclassified))

# Dissolve internal borders (creates one unified polygon)
boundary_vect <- aggregate(boundary_vect)

# -------------
# Species list
# -------------
species_list <- unique(species_hs_ref$Species)

# =========================================================
# Main loop: process each species raster independently
# =========================================================
for (species in species_list) {
  
  # Subset suitability rules for species
  spp_data <- subset(species_hs_ref, Species == species)
  
  # Create a 2-column reclassification matrix: landcover_class -> suitability (0/1)
  spp_matrix <- as.matrix(spp_data[, 2:3])
  
  # Reclassify raster to binary suitability
  species_habitat <- classify(land_cover_reclassified, spp_matrix)
  
  # Mask to Britain boundary (ensures outside = NA)
  species_habitat <- mask(species_habitat, boundary_vect)
  
  # File names (Fix 6: consistent file.path usage)
  jpeg_filename <- file.path(out_dir_res, paste0(species, " Hs Britain Map.jpg"))
  tif_filename  <- file.path(out_dir_data, paste0(species, " Hs Britain Map.tif"))
  
  # ---------------------
  # Save single-species GeoTIFF
  # ---------------------
  writeRaster(species_habitat, filename = tif_filename, overwrite = TRUE)
}

# ---------------------------------------------
# Plot all suitable habitat maps in one figure
# ---------------------------------------------
species_files <- list.files(path = out_dir_data,
                            pattern = "Hs Britain Map\\.tif$",
                            full.names = TRUE)

species_rasters <- lapply(species_files, rast)

# Clean species names for titles
species_names <- basename(species_files)
species_names <- gsub(" Hs Britain Map\\.tif$", "", species_names)

# Layout
n_cols <- 5
n_rows <- ceiling(length(species_rasters) / n_cols)

jpeg_filename_all <- file.path(out_dir_res, "All_Species_Suitable_Habitats.jpg")
jpeg(jpeg_filename_all, width = 5200, height = 3100, res = 200)

par(mfrow = c(n_rows, n_cols),
    mar = c(0.5, 0.5, 0.5, 0.5),
    oma = c(0, 0, 6, 0))

for (i in seq_along(species_rasters)) {
  r <- mask(species_rasters[[i]], boundary_vect)
  
  plot(r,
       col = c("white", "forestgreen"),
       legend = FALSE,
       axes = FALSE,
       box = FALSE)
  
  plot(boundary_vect, add = TRUE, border = "black", lwd = 1.8)
  
  text(x = (xmin(r) + xmax(r)) / 2,
       y = ymax(r) - 0.02 * (ymax(r) - ymin(r)),
       labels = species_names[i],
       cex = 2.5,
       font = 2,
       col = "black",
       adj = c(0.5, 1))
}

mtext("Species Suitable Habitats", outer = TRUE, cex = 4, line = 2, font = 2)

dev.off()
message("✅ Combined habitat map saved to: ", jpeg_filename_all)
