############################################################
# Species habitat maps, patch maps, and overlap heatmap (GB)
############################################################
# Description - This script visualises final species habitat outputs by (1) plotting
# binary suitable habitat maps for all species in a multi-panel figure, (2) plotting
# patch-ID maps (one panel per species) with dispersal barriers overlaid, and (3)
# computing and mapping the number of species overlapping per pixel across Great
# Britain as a heatmap. Outputs are saved as publication-ready JPEGs.

# Clean session (optional)
# rm(list = ls())

# ----------------
# Load packages
# ----------------
library(terra)
library(dplyr)
library(RColorBrewer)

# -----------------------
# Settings / I/O paths
# -----------------------
set.seed(13)  # kept for consistency across scripts

dir.create("./_terra_tmp", showWarnings = FALSE, recursive = TRUE)
terraOptions(memfrac = 0.5, tempdir = "./_terra_tmp")

data_path <- "./data/Suitable Habitats with Barriers/"
map_dir   <- "./results/Suitable Habitats Buffered and with Boundaries/"
out_dir_2 <- "./results/Rewilding Hotspots/"

dir.create(map_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_2, showWarnings = FALSE, recursive = TRUE)

# -----------------------
# Boundary
# -----------------------
boundary_vect <- vect("./data/Britain Boundary Polygon.geojson")
boundary_vect <- project(boundary_vect, "EPSG:27700")
boundary_vect <- aggregate(boundary_vect)  # dissolved GB polygon

# -----------------------
# Barriers
# -----------------------
barriers_vect <- vect("./data/Combined Dispersal Barriers.geojson")
barriers_vect <- project(barriers_vect, "EPSG:27700")
barriers_vect <- intersect(barriers_vect, boundary_vect)

# ---------------------
# Load species rasters (binary presence)
# ---------------------
hs_files <- list.files(
  data_path,
  pattern = "Suitable Habitat Split and Filtered\\.tif$",
  full.names = TRUE
)
if (length(hs_files) == 0) stop("No rasters found in ", data_path)

species_rasters <- list()
for (f in hs_files) {
  species_name <- gsub(" Suitable Habitat Split and Filtered\\.tif$", "", basename(f))
  species_name <- gsub("_", " ", species_name)
  message("Loading species: ", species_name)
  
  r <- rast(f)
  r <- ifel(r > 0, 1, NA)   # binary presence for the multi-panel map + overlap heatmap
  species_rasters[[species_name]] <- r
}

# Optional exclusion
all_species     <- names(species_rasters)
species_to_plot <- setdiff(all_species, c("Wolf", "Wolverine"))

# ==========================================================
# 1) Multi-panel map: binary suitable habitat for each species
# ==========================================================
message("📊 Plotting all species suitable habitats...")

n_species <- length(species_rasters)
if (n_species == 0) stop("No suitable habitat rasters loaded!")

species_names <- names(species_rasters)

n_cols <- 5
n_rows <- ceiling(n_species / n_cols)

jpeg_filename_all <- file.path(map_dir, "All_Species_Suitable_Habitats.jpg")
jpeg(jpeg_filename_all, width = 5200, height = 3100, res = 200)

par(mfrow = c(n_rows, n_cols),
    mar = c(0, 0, 0, 0),
    oma = c(0, 0, 6, 0))

for (sp in species_names) {
  r <- mask(species_rasters[[sp]], boundary_vect)
  
  plot(r, col = c("white", "forestgreen"), legend = FALSE, axes = FALSE, box = FALSE)
  plot(boundary_vect, add = TRUE, border = "black", lwd = 1.8)
  
  text(x = (xmin(r) + xmax(r)) / 2,
       y = ymax(r) - 0.02 * (ymax(r) - ymin(r)),
       labels = sp,
       cex = 2.5, font = 2, col = "black", adj = c(0.5, 1))
}

mtext("Species Specific Rewilding Sites", outer = TRUE, cex = 4, line = 2, font = 2)
dev.off()
message("✅ Saved: ", jpeg_filename_all)

# ==========================================================
# 2) Multi-panel map: patch-ID rasters with barriers (Paired)
# ==========================================================
message("🎨 Plotting patch-ID maps with barriers...")

# Read patch-ID rasters from files (not binary)
species_rasters_filtered <- list()
for (sp in species_to_plot) {
  raster_file <- file.path(
    data_path,
    paste0(gsub(" ", "_", sp), " Suitable Habitat Split and Filtered.tif")
  )
  
  if (file.exists(raster_file)) {
    species_rasters_filtered[[sp]] <- rast(raster_file)
  } else {
    message("Raster file not found for ", sp, ", skipping.")
  }
}

species_to_plot2 <- names(species_rasters_filtered)
if (length(species_to_plot2) == 0) stop("No patch-ID rasters found to plot.")

n_cols <- 4
n_rows <- ceiling(length(species_to_plot2) / n_cols)

jpeg_filename_patches <- file.path(map_dir, "All_Species_Suitable_Habitats_Patch_IDs_Paired_Barriers.jpg")
jpeg(jpeg_filename_patches, width = 5200, height = 3100, res = 200)

par(mfrow = c(n_rows, n_cols),
    mar = c(0, 0, 0, 0),
    oma = c(0, 0, 6, 0))

paired12 <- brewer.pal(12, "Paired")

for (sp in species_to_plot2) {
  
  message("Processing: ", sp)
  r <- mask(species_rasters_filtered[[sp]], boundary_vect)
  
  ids <- terra::freq(r)$value
  ids <- ids[!is.na(ids)]
  n_id <- length(ids)
  if (n_id == 0) next
  
  cols <- rep(paired12, length.out = n_id)
  cols <- sample(cols)
  patch_cols <- setNames(cols, as.character(ids))
  
  plot(r, col = patch_cols, legend = FALSE, axes = FALSE, box = FALSE)
  plot(boundary_vect, add = TRUE, border = "black", lwd = 1.3)
  
  # Add barriers (CRS-safe)
  if (!is.null(barriers_vect)) {
    barriers_to_plot <- barriers_vect
    if (!terra::same.crs(r, barriers_to_plot)) {
      barriers_to_plot <- terra::project(barriers_to_plot, crs(r))
    }
    if (terra::relate(ext(r), ext(barriers_to_plot), "intersects")) {
      plot(barriers_to_plot, add = TRUE, col = "#404040", border = "#404040", lwd = 0.6)
    }
  }
  
  text(x = (xmin(r) + xmax(r)) / 2,
       y = ymax(r) - 0.02 * (ymax(r) - ymin(r)),
       labels = paste0(sp, " - ", n_id, " patch", ifelse(n_id > 1, "es", "")),
       cex = 2.5, font = 2, col = "black", adj = c(0.5, 1))
}

mtext("Species-Specific Rewilding Sites", outer = TRUE, cex = 4, line = 1.5, font = 2)
dev.off()
message("✅ Saved: ", jpeg_filename_patches)

# ==========================================================
# 3) Heatmap: number of overlapping species per pixel
# ==========================================================
message("🔥 Building and plotting species overlap heatmap...")

template <- rast(ext(species_rasters[[1]]), resolution = 50, crs = "EPSG:27700")
species_stack <- rast(lapply(species_rasters, function(r) resample(r, template, method = "near")))

overlap_raster <- sum(species_stack, na.rm = TRUE)
overlap_raster_masked <- mask(overlap_raster, overlap_raster, maskvalues = 0)

max_overlap <- global(overlap_raster_masked, "max", na.rm = TRUE)[[1]]
max_overlap <- min(max_overlap, 8)

cols <- brewer.pal(n = max_overlap, name = "YlOrRd")

overlap_factor <- as.factor(round(overlap_raster_masked))
levels(overlap_factor) <- data.frame(
  value = 1:max_overlap,
  label = as.character(1:max_overlap)
)

jpeg(file.path(out_dir_2, "Species_Overlap_Heatmap_All_Species.jpeg"),
     width = 2000, height = 1800, res = 300)

plot(overlap_factor,
     col = cols,
     main = "Species Overlap",
     axes = FALSE,
     box = FALSE,
     legend = TRUE,
     colNA = "white",
     plg = list(
       title = "Number of species",
       cex = 1.5,
       title.cex = 0.9,
       title.font = 2
     ))

plot(boundary_vect, add = TRUE, border = "black", lwd = 1.5)
dev.off()

message("✅ Saved overlap heatmap to: ", file.path(out_dir_2, "Species_Overlap_Heatmap_All_Species.jpeg"))
