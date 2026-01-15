#############################################################
# Splitting buffered suitable habitat patches by barriers   #
#############################################################
# Description - This script refines species-specific buffered habitat maps by
# explicitly accounting for dispersal barriers and minimum patch size constraints.
# For each species, buffered suitable habitat rasters are converted to polygons,
# convex hulls are generated for contiguous patches, and these patches are split
# using multiple barrier types (roads, rivers, canals, railways, and coastline).
# Resulting habitat fragments are rasterised, reclassified, and filtered to retain
# only patches exceeding species-specific minimum area requirements. The workflow
# outputs final filtered habitat rasters, spatial patch layers, patch-level and
# species-level summary tables, and publication-ready maps with dispersal barriers
# overlaid.

# Clean session (optional)
# rm(list = ls())

# Load Packages
library(terra)
library(dplyr)
library(readxl)
library(sf)
library(writexl)

# -----------------------
# Settings / I/O paths
# -----------------------
set.seed(13)  # kept for consistency across scripts

# Use more RAM to reduce disk thrash; put tempdir on a fast SSD
dir.create("./_terra_tmp", showWarnings = FALSE, recursive = TRUE)

terraOptions(
  memfrac = 0.7,                     # adjust for your machine
  tempdir = "./_terra_tmp"
)

data_path <- "./data/Suitable Habitats Buffered/"
out_dir   <- "./data/Suitable Habitats with Barriers/"
out_dir_2 <- "./results/Suitable Habitats Buffered and with Boundaries/"

dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_2, showWarnings = FALSE, recursive = TRUE)

# -----------------------
# Species parameters
# -----------------------
min_area <- read_excel("./data/Species Minimum Area Required.xlsx")
min_area_m2 <- setNames(min_area$`min_area (m2)`, min_area$species)
rm(min_area)

# -----------------------
# Boundary
# -----------------------
boundary_vect <- vect("./data/Britain Boundary Polygon.geojson")
boundary_vect <- project(boundary_vect, "EPSG:27700")

# Dissolve Britain boundary to remove internal borders
boundary_union <- terra::aggregate(boundary_vect)

# -----------------------
# Find suitability rasters
# -----------------------
hs_files <- list.files(
  data_path,
  pattern = " Suitable Habitat Buffered\\.tif$",
  full.names = TRUE
)

if (length(hs_files) == 0) stop("No 'Suitable Habitat Buffered.tif' files found in ", data_path)

# -----------------------
# Dispersal Barriers
# -----------------------

barriers_all <- terra::vect("./data/Combined Dispersal Barriers.geojson") %>%
  terra::makeValid()

# Only add coastline here if it is NOT already present
if (!("type" %in% names(barriers_all))) {
  stop("Combined Dispersal Barriers is missing the expected 'type' attribute.")
}

if (!("coast" %in% unique(barriers_all$type))) {
  boundary_line <- terra::as.lines(boundary_union)
  boundary_line$type <- "coast"
  barriers_all <- rbind(barriers_all, boundary_line)
}

# Create merged version for plotting (optional)
barriers_merged <- terra::aggregate(barriers_all)

# =====================================================
# Main loop: process each species raster independently #
# =====================================================

# Process all species
for (hs_path in hs_files) {
  
  message("\n==============================")
  message("Processing: ", basename(hs_path))
  message("==============================")
  
  # -----------------------
  # Normalize species name from filename
  # -----------------------
  species_name_raw <- gsub(" Suitable Habitat Buffered\\.tif$", "", basename(hs_path))
  
  # Convert underscores to spaces and trim
  species_name_norm <- gsub("_", " ", species_name_raw)
  species_name_norm <- trimws(species_name_norm)
  
  # For robust matching: lowercase
  species_name_lower <- tolower(species_name_norm)
  
  # Normalize min_area table names for flexible matching
  min_area_names_lower <- tolower(gsub("_", " ", names(min_area_m2)))
  
  # -----------------------
  # Match species name in the minimum area table
  # -----------------------
  if (!(species_name_lower %in% min_area_names_lower)) {
    stop("Species '", species_name_raw, "' not found in minimum area table.")
  }
  
  min_area_val <- min_area_m2[which(min_area_names_lower == species_name_lower)]
  
  # Use underscored version for file naming
  species_filename <- gsub(" ", "_", species_name_norm)
  
  # -----------------------
  # Load raster
  # -----------------------
  valid_cluster_raster <- rast(hs_path)
  
  # Ensure raster CRS matches EPSG:27700
  if (is.na(crs(valid_cluster_raster, proj = TRUE)) || crs(valid_cluster_raster) != "EPSG:27700") {
    valid_cluster_raster <- project(valid_cluster_raster, "EPSG:27700", method = "near")
  }
  
  plot(valid_cluster_raster)
  
  # -----------------------
  # Convert raster to polygons and create convex hulls
  # -----------------------
  patch_polys <- terra::as.polygons(valid_cluster_raster)
  hulls <- terra::hull(patch_polys, "convex", by = "patches")
  plot(patch_polys)
  plot(hulls)
  
  # Clip hulls to Britain outline
  hulls_clipped <- terra::intersect(hulls, boundary_union)
  hulls_clipped <- terra::makeValid(hulls_clipped)
  hulls_clipped <- terra::buffer(hulls_clipped, 0)
  
  plot(hulls_clipped)
  
  # -----------------------
  # Split by barrier types
  # -----------------------
  non_rail_types <- setdiff(unique(barriers_all$type), "rail")
  split_patches <- hulls_clipped
  
  for (bt in non_rail_types) {
    message("Splitting by barrier type: ", bt)
    b <- barriers_all[barriers_all$type == bt, ] %>% terra::makeValid()
    split_patches <- terra::split(split_patches, b)
  }
  
  # Handle rail barriers incrementally
  rail_barriers <- barriers_all[barriers_all$type == "rail", ] %>%
    terra::makeValid()
  
  # Add a small buffer to ensure intersections register
  rail_barriers <- terra::buffer(rail_barriers, width = 25)
  
  new_split <- list()
  for (i in 1:length(split_patches)) {
    patch <- split_patches[i]
    pieces <- tryCatch({
      terra::split(patch, rail_barriers)
    }, error = function(e) patch)
    new_split <- c(new_split, pieces)
  }
  
  # Ensure we produced something to rasterize
  if (length(new_split) == 0) {
    warning("No split patches produced for ", species_name_norm, "; skipping outputs.")
    next
  }
  
  split_patches <- do.call(c, new_split)
  
  # -----------------------
  # Rasterize split patches
  # -----------------------
  all_patches <- do.call(rbind, new_split)
  all_patches <- terra::makeValid(all_patches)
  all_patches$patch_id <- seq_len(nrow(all_patches))
  
  split_raster <- terra::rasterize(all_patches, valid_cluster_raster,
                                   field = "patch_id", touches = TRUE)
  
  # -----------------------
  # Reclassify and mask
  # -----------------------
  vals_orig <- values(valid_cluster_raster)
  vals_split <- values(split_raster)
  vals_orig[!is.na(vals_orig)] <- vals_split[!is.na(vals_orig)]
  reclassified_raster <- valid_cluster_raster
  values(reclassified_raster) <- vals_orig
  masked_raster <- terra::mask(reclassified_raster, boundary_vect)
  
  # -----------------------
  # Filter by minimum area
  # -----------------------
  cell_area <- terra::cellSize(masked_raster, unit = "m")
  cluster_areas <- terra::zonal(cell_area, masked_raster, fun = "sum", na.rm = TRUE)
  colnames(cluster_areas) <- c("cluster", "area")
  valid_clusters <- cluster_areas %>% filter(area >= min_area_val)
  valid_cluster_ids <- as.numeric(valid_clusters$cluster)
  
  masked_raster <- terra::as.int(masked_raster)
  
  # === FAIL-SAFE CHECK ===
  if (length(valid_cluster_ids) == 0 || all(is.na(valid_cluster_ids))) {
    message("⚠️ No valid clusters remain after area filtering for ", species_name_norm)
    filtered_final_raster <- masked_raster
    filtered_final_raster[] <- NA
  } else {
    filtered_final_raster <- terra::ifel(masked_raster %in% valid_cluster_ids,
                                         masked_raster,
                                         NA)
  }
  
  # -----------------------
  # Handle case where no patches remain
  # -----------------------
  if (length(valid_cluster_ids) == 0) {
    message("⚠️ No suitable patches remaining for ", species_name_norm,
            " after applying minimum area filter. Creating empty outputs.")
    
    # Empty raster
    filtered_final_raster[] <- NA
    
    # Empty polygon layer (Britain outline)
    filtered_polys <- boundary_union
    filtered_polys$patch_id <- NA
    filtered_polys$area_m2 <- NA
    filtered_polys$species <- species_name_norm
    
    # Empty patch/summary tables
    patch_df <- data.frame(species = character(0), patch_id = numeric(0), area_m2 = numeric(0))
    summary_df <- data.frame(
      species = species_name_norm,
      n_patches = 0,
      mean_patch_size_m2 = NA_real_
    )
    
  } else {
    # Normal case
    filtered_polys <- terra::as.polygons(filtered_final_raster, dissolve = TRUE, na.rm = TRUE)
    names(filtered_polys) <- "patch_id"
    filtered_polys$area_m2 <- terra::expanse(filtered_polys, unit = "m")
    filtered_polys$species <- species_name_norm
    
    patch_df <- as.data.frame(filtered_polys)[, c("species", "patch_id", "area_m2")]
    summary_df <- data.frame(
      species = species_name_norm,
      n_patches = nrow(patch_df),
      mean_patch_size_m2 = mean(patch_df$area_m2, na.rm = TRUE)
    )
  }
  
  # -----------------------
  # Save raster outputs
  # -----------------------
  tif_filename  <- paste0(out_dir, species_filename, " Suitable Habitat Split and Filtered.tif")
  jpeg_filename <- paste0(out_dir_2, species_filename, " Suitable Habitat Split and Filtered.jpg")
  
  writeRaster(filtered_final_raster, tif_filename, overwrite = TRUE)
  
  jpeg(jpeg_filename, width = 1200, height = 1200, res = 150)
  plot(filtered_final_raster, col = rep("forestgreen", 100), legend = FALSE,
       main = paste0(species_name_norm, " Split and Filtered Suitable Habitat"),
       axes = FALSE, box = FALSE)
  plot(boundary_vect, add = TRUE, border = "black", lwd = 1.5)
  dev.off()
  
  # -----------------------
  # Patch summary outputs
  # -----------------------
  gpkg_filename  <- paste0(out_dir, species_filename, " Patch Summary.gpkg")
  excel_filename <- paste0(out_dir, species_filename, " Patch Summary.xlsx")
  
  terra::writeVector(filtered_polys, gpkg_filename, filetype = "GPKG", overwrite = TRUE)
  
  writexl::write_xlsx(list(
    Patches = patch_df,
    Summary = summary_df
  ), excel_filename)
  
  message("✅ Done: ", species_name_norm)
  
}

# =====================================================
# Plot all final habitat maps with barriers and legend
# =====================================================

message("\n==============================")
message("Plotting all final habitat maps with barriers")
message("==============================")

# Find all final rasters
final_tifs <- list.files(out_dir,
                         pattern = "Suitable Habitat Split and Filtered\\.tif$",
                         full.names = TRUE)

if (length(final_tifs) == 0) stop("No final rasters found in ", out_dir)

# Define consistent colors for each barrier type
barrier_colors <- c(
  river = "#0571b0",
  road  = "#ca0020",
  rail  = "#7b3294",
  canal = "#e66101",
  coast = "gray30"
)

# Check which types exist in your data and match colors
barrier_types <- intersect(names(barrier_colors), unique(barriers_all$type))

for (tif_path in final_tifs) {
  
  species_name <- gsub(" Suitable Habitat Split and Filtered\\.tif$", "", basename(tif_path))
  species_name <- gsub("_", " ", species_name)
  
  message("Plotting: ", species_name)
  
  # Load raster
  r <- terra::rast(tif_path)
  
  # Define output JPEG path
  jpeg_filename <- paste0(out_dir_2, gsub(" ", "_", species_name),
                          " Suitable Habitat Split and Filtered with Barriers.jpg")
  
  # Open graphics device
  jpeg(jpeg_filename, width = 1200, height = 1200, res = 150)
  
  # Expand margins to make space for the legend
  par(mar = c(4, 4, 4, 4))  # bottom, left, top, right
  
  # Plot the suitable habitat
  plot(r,
       col = rep("forestgreen", 100),
       legend = FALSE,
       main = paste0(species_name, " Suitable Habitat with Barriers"),
       axes = FALSE, box = FALSE)
  
  # Add Britain boundary
  plot(boundary_vect, add = TRUE, border = "black", lwd = 1.5)
  
  # Plot barriers
  for (bt in barrier_types) {
    if (bt %in% barriers_all$type) {
      plot(barriers_all[barriers_all$type == bt, ],
           add = TRUE, col = barrier_colors[bt], lwd = 1)
    }
  }
  
  # Add legend (inset slightly above bottom-left corner)
  legend("topright",
         inset = c(5, 0),  # left as-is per request
         legend = c("Suitable Habitat", names(barrier_colors[barrier_types])),
         col = c("forestgreen", barrier_colors[barrier_types]),
         lwd = c(10, rep(2, length(barrier_types))),
         bty = "n", cex = 0.9,
         title = "Map Features")
  
  dev.off()
}

message("✅ All habitat maps plotted with barriers overlay and legend.")
