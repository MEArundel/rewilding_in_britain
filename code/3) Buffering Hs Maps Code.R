####################################################
# Assigning species dispersal distance as a buffer #
####################################################
# Description - This script expands species-specific suitable habitat maps by
# applying dispersal-distance buffers to account for landscape connectivity and
# potential movement. For each species, binary habitat suitability rasters are
# buffered by species-specific dispersal distances, processed separately within
# England, Scotland, and Wales, and mosaicked into a single Great Britain–wide
# raster. Resulting habitat patches are then clustered and filtered to retain
# only those exceeding species-specific minimum area requirements. The script
# exports buffered and filtered habitat rasters as GeoTIFFs, generates individual
# JPEG maps for each species, and produces a combined multi-panel figure for
# comparative visualisation.

# Clean session (optional)
# rm(list = ls())

# Load packages
library(terra)
library(readxl)

# -----------------------
# Settings / I-O paths
# -----------------------
# set.seed(13) # no randomness used; keep only if you want consistency across scripts

in_dir            <- "./data/"
in_hs_dir         <- file.path(in_dir, "Suitable Habitats")
out_hs_buf_dir    <- file.path(in_dir, "Suitable Habitats Buffered")
out_hs_buf_mapdir <- "./results/Suitable Habitats Buffered/"
terra_tmp_dir     <- "./_terra_tmp"

dir.create(out_hs_buf_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(out_hs_buf_mapdir, showWarnings = FALSE, recursive = TRUE)
dir.create(terra_tmp_dir,     showWarnings = FALSE, recursive = TRUE)

# Use more RAM to reduce disk thrash; put tempdir on a fast SSD
terraOptions(
  memfrac = 0.7,               # adjust for your machine
  tempdir = terra_tmp_dir
)

# -----------------------
# Species parameters
# -----------------------
min_area_path   <- file.path(in_dir, "Species Minimum Area Required.xlsx")
dispersal_path  <- file.path(in_dir, "Species Dispersal Distances.xlsx")
boundary_path   <- file.path(in_dir, "Britain Boundary Polygon.geojson")

stopifnot(file.exists(min_area_path))
stopifnot(file.exists(dispersal_path))
stopifnot(file.exists(boundary_path))

min_area <- read_excel(min_area_path)
min_area_m2 <- setNames(min_area$`min_area (m2)`, min_area$species)
rm(min_area)

dispersal_df <- read_excel(dispersal_path)
dispersal_m  <- setNames(dispersal_df$dispersal_dist_m, dispersal_df$species)
rm(dispersal_df)

# -----------------------
# Boundary + regions
# -----------------------
boundary_vect <- vect(boundary_path)
boundary_vect <- project(boundary_vect, "EPSG:27700")

# IMPORTANT: subset regions BEFORE any operations that might drop attributes
if (!("shape1" %in% names(boundary_vect))) {
  stop("Expected attribute 'shape1' not found in boundary file: ", boundary_path)
}

regions <- list(
  England  = boundary_vect[boundary_vect$shape1 == "England", ],
  Scotland = boundary_vect[boundary_vect$shape1 == "Scotland", ],
  Wales    = boundary_vect[boundary_vect$shape1 == "Wales", ]
)

# Dissolved GB outline for plotting and masking where needed
boundary_dissolved <- aggregate(boundary_vect)

# -----------------------
# Find suitability rasters
# -----------------------
hs_files <- list.files(
  in_hs_dir,
  pattern = " Hs Britain Map\\.tif$",
  full.names = TRUE
)
if (length(hs_files) == 0) stop("No 'Hs Britain Map.tif' files found in ", in_hs_dir)

# -----------------------
# Helper: write compressed GTiff
# -----------------------
write_gtiff <- function(x, filename) {
  wopt <- list(
    datatype = "INT4S",
    gdal     = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER")
  )
  writeRaster(x, filename = filename, overwrite = TRUE, wopt = wopt)
}

# =========================================================
# Main loop: process each species raster independently
# =========================================================
for (hs_path in hs_files) {
  
  message("\n==============================")
  message("Processing: ", basename(hs_path))
  message("==============================")
  
  # -----------------------
  # Normalize species name from filename
  # -----------------------
  species_name_raw <- gsub(" Hs Britain Map\\.tif$", "", basename(hs_path))
  species_name_norm <- trimws(gsub("_", " ", species_name_raw))
  species_name_lower <- tolower(species_name_norm)
  
  # Normalize table names for robust matching
  dispersal_names_lower <- tolower(gsub("_", " ", names(dispersal_m)))
  min_area_names_lower  <- tolower(gsub("_", " ", names(min_area_m2)))
  
  # Lookup species parameters
  if (!(species_name_lower %in% dispersal_names_lower)) {
    stop("Species '", species_name_raw, "' not found in dispersal table.")
  }
  if (!(species_name_lower %in% min_area_names_lower)) {
    stop("Species '", species_name_raw, "' not found in minimum area table.")
  }
  
  buffer_width <- dispersal_m[which(dispersal_names_lower == species_name_lower)]
  min_area_val <- min_area_m2[which(min_area_names_lower == species_name_lower)]
  
  # Use underscored name for file outputs
  species_filename <- gsub(" ", "_", species_name_norm)
  
  # -----------------------
  # Load raster, ensure EPSG:27700 and binary (1/NA)
  # -----------------------
  r <- rast(hs_path)
  
  if (is.na(crs(r, proj = TRUE)) || crs(r) != "EPSG:27700") {
    r <- project(r, "EPSG:27700", method = "near")
  }
  
  r <- ifel(r == 1, 1, NA)
  r <- mask(r, boundary_dissolved)
  
  if (global(!is.na(r), "sum", na.rm = TRUE)[1, 1] == 0) {
    warning("No suitable habitat (value==1) in ", basename(hs_path), "; skipping.")
    next
  }
  
  # -----------------------
  # Regional dilation -> mosaic
  # -----------------------
  buf_list <- list()
  
  for (region_name in names(regions)) {
    message("  Region: ", region_name)
    
    region <- regions[[region_name]]
    region_ext <- buffer(region, width = buffer_width)
    
    rr <- crop(r, region_ext)
    rr <- mask(rr, region_ext)
    
    if (is.null(rr) || ncell(rr) == 0 ||
        global(!is.na(rr), "sum", na.rm = TRUE)[1, 1] == 0) {
      message("    (no suitable cells here)")
      next
    }
    
    message("    computing distance()...")
    d <- distance(rr)
    
    message("    thresholding to buffer halo...")
    rr_buf <- ifel(d <= buffer_width, 1, NA)
    names(rr_buf) <- paste0("buf_", region_name)
    
    buf_list[[region_name]] <- rr_buf
  }
  
  if (length(buf_list) == 0) {
    warning("No buffered regions produced for ", species_name_raw, "; skipping.")
    next
  }
  
  message("  Mosaicing regional buffers...")
  r_buf_gb <- do.call(mosaic, c(unname(buf_list), list(fun = "max")))
  
  # -----------------------
  # GB-wide clumps
  # -----------------------
  message("  Clumping GB-wide (patches on dilated raster)...")
  r_clump <- patches(r_buf_gb, directions = 8)
  r_clump <- extend(r_clump, r)
  cluster_ids <- mask(r_clump, r)
  
  # -----------------------
  # Filter by minimum area
  # -----------------------
  message("  Computing cluster areas...")
  cell_area <- cellSize(cluster_ids, unit = "m")
  cl_area <- zonal(cell_area, cluster_ids, fun = "sum", na.rm = TRUE)
  colnames(cl_area) <- c("cluster", "area_m2")
  
  keep_ids <- cl_area[, "cluster"][cl_area[, "area_m2"] >= min_area_val]
  
  if (length(keep_ids) == 0) {
    warning("No clusters meet area threshold for ", species_name_raw, "; skipping outputs.")
    next
  }
  
  message("  Filtering valid clusters (", length(keep_ids), " kept)...")
  keep_mask <- app(cluster_ids, fun = function(x) {
    x[!(x %in% keep_ids)] <- NA
    x[!is.na(x)] <- 1
    x
  })
  
  valid_cluster_raster <- mask(cluster_ids, keep_mask)
  
  # -----------------------
  # Save outputs
  # -----------------------
  tif_filename  <- file.path(out_hs_buf_dir,    paste0(species_filename, " Suitable Habitat Buffered.tif"))
  jpeg_filename <- file.path(out_hs_buf_mapdir, paste0(species_filename, " Suitable Habitat Buffered.jpg"))
  
  message("  Writing GeoTIFF: ", tif_filename)
  write_gtiff(valid_cluster_raster, tif_filename)
  
  message("  Writing JPEG: ", jpeg_filename)
  jpeg(jpeg_filename, width = 1200, height = 1200, res = 150)
  plot(valid_cluster_raster,
       main = paste0(species_name_norm, " Buffered & Filtered Suitable Habitat"),
       axes = FALSE, box = FALSE)
  plot(boundary_dissolved, add = TRUE, border = "black", lwd = 1.2)
  dev.off()
  
  message("Done: ", species_name_norm)
}

message("\nAll species processed.")

# ------------------------------------------------------
# Plot all buffered suitable habitat maps in one figure
# ------------------------------------------------------
# Note: this summary plot is binary and mainly for quick comparison

species_files <- list.files(
  path = out_hs_buf_dir,
  pattern = " Suitable Habitat Buffered\\.tif$",
  full.names = TRUE
)

sd_maps <- lapply(species_files, rast)

species_names <- basename(species_files)
species_names <- gsub(" Suitable Habitat Buffered\\.tif$", "", species_names)
species_names <- gsub("_", " ", species_names)

n_cols <- 5
n_rows <- ceiling(length(sd_maps) / n_cols)

jpeg_filename_all <- file.path(out_hs_buf_mapdir, "All_Species_Buffered_Suitable_Habitats.jpg")
jpeg(jpeg_filename_all, width = 5200, height = 3100, res = 200)

par(mfrow = c(n_rows, n_cols),
    mar = c(0.5, 0.5, 0.5, 0.5),
    oma = c(0, 0, 6, 0))

for (i in seq_along(sd_maps)) {
  rr <- mask(sd_maps[[i]], boundary_dissolved)
  
  plot(rr,
       col = c("white", "forestgreen"),
       legend = FALSE,
       axes = FALSE,
       box = FALSE)
  
  plot(boundary_dissolved, add = TRUE, border = "black", lwd = 1.8)
  
  text(x = (xmin(rr) + xmax(rr)) / 2,
       y = ymax(rr) - 0.02 * (ymax(rr) - ymin(rr)),
       labels = species_names[i],
       cex = 2.5,
       font = 2,
       col = "black",
       adj = c(0.5, 1))
}

mtext("Species Buffered Suitable Habitats", outer = TRUE, cex = 4, line = 2, font = 2)

dev.off()
message("✅ Combined habitat map saved to: ", jpeg_filename_all)
