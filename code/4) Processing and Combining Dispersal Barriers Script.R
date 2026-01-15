###############################################
# Dispersal barriers: build + combine (GB)    #
###############################################
# Description - This script (1) identifies and exports four dispersal barrier
# layers (major roads, canals, railways, wide rivers) from source datasets,
# produces a JPEG map for each barrier type, and (2) combines all barrier layers
# (plus coastline) into a single GeoJSON and a combined JPEG map.

# Clean session (optional)
# rm(list = ls())

# ----------------
# Load packages
# ----------------
library(terra)
library(sf)
library(tidyterra)
library(dplyr)

# -----------------------
# Settings / I/O paths
# -----------------------
set.seed(13)  # kept for consistency across scripts

# Base paths
boundary_geojson <- "./data/Britain Boundary Polygon.geojson"
data_dir         <- "./data/"
plots_path       <- "./results/Dispersal Barriers/"

# -----------------------
# Inputs and checks
# -----------------------
# Barrier-specific input dirs
roads_dir   <- "./data/Road Data/"
rivers_dir  <- "./data/River Data/"
rails_dir   <- "./data/Railway Network shp Data/"

# Ensure output dirs exist
dir.create(plots_path, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(boundary_geojson))
stopifnot(dir.exists(data_dir))
stopifnot(dir.exists(roads_dir))
stopifnot(dir.exists(rivers_dir))
stopifnot(dir.exists(rails_dir))

# -----------------------
# Boundary (GB dissolved)
# -----------------------
boundary_vect <- vect(boundary_geojson)
boundary_vect <- project(boundary_vect, "EPSG:27700")
boundary_vect <- aggregate(boundary_vect)  # dissolved GB polygon

# sf version for sf ops
boundary_sf <- st_as_sf(boundary_vect)

# ============================================================
# 1) Identification of major road networks (A roads + motorways)
# ============================================================

roadlink_files <- list.files(
  path = roads_dir,
  pattern = "RoadLink.*\\.shp$",
  full.names = TRUE
)
if (length(roadlink_files) == 0) stop("No RoadLink shapefiles found in: ", roads_dir)

roads_list <- lapply(roadlink_files, vect)
roads <- do.call(rbind, roads_list)

# Ensure CRS is EPSG:27700
if (is.na(crs(roads, proj = TRUE)) || crs(roads) != "EPSG:27700") {
  roads <- project(roads, "EPSG:27700")
}

# Filter: Keep only Motorways and specified A Roads
filtered_roads <- roads[
  roads$`function` == "Motorway" |
    (roads$`function` == "A Road" & roads$formOfWay %in% c("Collapsed Dual Carriageway", "Dual Carriageway")),
]

# Filter by length (grouped by roadNumber)
filtered_roads <- filtered_roads %>%
  group_by(roadNumber) %>%
  mutate(length_group = dplyr::if_else(sum(length) >= 500, "long", "short")) %>%
  ungroup()

# Remove 'short' roads
filtered_roads <- filtered_roads[filtered_roads$length_group == "long", ]

# Save roads
writeVector(filtered_roads, file.path(data_dir, "A Roads and Motorways.geojson"), overwrite = TRUE)

# Plot roads
jpeg(filename = file.path(plots_path, "Roads Map.jpg"),
     width = 20, height = 15, units = "cm", res = 300)

par(mar = c(3, 1, 3, 1) + 0.1)
plot(boundary_vect, col = NA, border = "black", lwd = 0.75, main = "", axes = FALSE)
plot(filtered_roads, col = "#ca0020", add = TRUE, lwd = 0.75)
mtext("Dual Carriageways and Motorways", side = 3, outer = FALSE, cex = 1.5, line = 1)
dev.off()

# =====================================
# 2) Identification of canal networks
# =====================================

watercourse_files <- list.files(
  path = rivers_dir,
  pattern = "WatercourseLink\\.shp$",
  full.names = TRUE
)
if (length(watercourse_files) == 0) stop("No WatercourseLink.shp found in: ", rivers_dir)

rivers_list <- lapply(watercourse_files, vect)
rivers <- do.call(rbind, rivers_list)

# Ensure CRS is EPSG:27700
if (is.na(crs(rivers, proj = TRUE)) || crs(rivers) != "EPSG:27700") {
  rivers <- project(rivers, "EPSG:27700")
}

# Filter by form + length
canals <- rivers %>%
  tidyterra::filter(form == "canal" & length >= 500)

# Save canals
writeVector(canals, file.path(data_dir, "Canals.geojson"), overwrite = TRUE)

# Plot canals
jpeg(filename = file.path(plots_path, "Canals Map.jpg"),
     width = 20, height = 15, units = "cm", res = 300)

par(mar = c(3, 1, 3, 1) + 0.1)
plot(boundary_vect, col = NA, border = "black", lwd = 0.75, main = "", axes = FALSE)
plot(canals, col = "#e66101", add = TRUE, lwd = 0.75)
mtext("Canals", side = 3, outer = FALSE, cex = 1.5, line = 1)
dev.off()

# ======================================
# 3) Identification of railway networks
# ======================================

rail_shp <- file.path(rails_dir, "NetworkLinks.shp")
if (!file.exists(rail_shp)) stop("Railway shapefile not found: ", rail_shp)

railway_data_vect <- vect(rail_shp)

# Ensure CRS is EPSG:27700
if (is.na(crs(railway_data_vect, proj = TRUE)) || crs(railway_data_vect) != "EPSG:27700") {
  railway_data_vect <- project(railway_data_vect, "EPSG:27700")
}

# Save railways
writeVector(railway_data_vect, file.path(data_dir, "Railways.geojson"), overwrite = TRUE)

# Plot railways
jpeg(filename = file.path(plots_path, "Railways Map.jpg"),
     width = 20, height = 15, units = "cm", res = 300)

par(mar = c(3, 1, 3, 1) + 0.1)
plot(boundary_vect, col = NA, border = "black", lwd = 0.75, main = "", axes = FALSE)
plot(railway_data_vect, col = "#7b3294", add = TRUE, lwd = 0.75)
mtext("Railways", side = 3, outer = FALSE, cex = 1.5, line = 1)
dev.off()

# ==========================================
# 4) Identification of major river networks
# ==========================================

river_shp_files <- list.files(rivers_dir, pattern = "\\.shp$", full.names = TRUE)
if (length(river_shp_files) == 0) stop("No shapefiles found in: ", rivers_dir)

river_shp_list <- lapply(river_shp_files, st_read, quiet = TRUE)

# Reproject to EPSG:27700
river_shp_list <- lapply(river_shp_list, function(x) {
  if (st_crs(x) != 27700) st_transform(x, 27700) else x
})

# Keep only linework (ignore points)
river_lines <- river_shp_list[sapply(river_shp_list, function(x) {
  any(st_geometry_type(x) %in% c("LINESTRING", "MULTILINESTRING"))
})]
if (length(river_lines) == 0) stop("No LINESTRING/MULTILINESTRING layers found in: ", rivers_dir)

# Combine
river_sf <- bind_rows(river_lines)

# Clip to boundary
british_rivers <- st_intersection(river_sf, boundary_sf)

# Filter by width + length
wide_rivers <- british_rivers %>%
  filter(width_mean >= 100,
         Shape_Leng >= 0.5)

wide_rivers$form <- "river"

# Convert back to terra
wide_rivers_vect <- vect(wide_rivers)

# Ensure CRS is EPSG:27700
if (is.na(crs(wide_rivers_vect, proj = TRUE)) || crs(wide_rivers_vect) != "EPSG:27700") {
  wide_rivers_vect <- project(wide_rivers_vect, "EPSG:27700")
}

# Save wide rivers
writeVector(wide_rivers_vect, file.path(data_dir, "Wide Rivers.geojson"), overwrite = TRUE)

# Plot wide rivers
jpeg(filename = file.path(plots_path, "Wide Rivers Map.jpg"),
     width = 20, height = 15, units = "cm", res = 300)

par(mar = c(3, 1, 3, 1) + 0.1)
plot(boundary_vect, col = NA, border = "black", lwd = 0.75, axes = FALSE)
plot(wide_rivers_vect, col = "#0571b0", add = TRUE, lwd = 0.75)
mtext("Wide Rivers", side = 3, outer = FALSE, cex = 1.5, line = 1)
dev.off()

# =================================
# 5) Combining dispersal barriers
# =================================

geojson_files <- file.path(
  data_dir,
  c("Canals.geojson", "A Roads and Motorways.geojson", "Railways.geojson", "Wide Rivers.geojson")
)
if (!all(file.exists(geojson_files))) {
  missing_files <- geojson_files[!file.exists(geojson_files)]
  stop("Missing expected GeoJSON outputs:\n", paste(missing_files, collapse = "\n"))
}

barrier_list <- lapply(geojson_files, vect)
combined_barriers <- do.call(rbind, barrier_list)

# Ensure CRS is EPSG:27700
if (is.na(crs(combined_barriers, proj = TRUE)) || crs(combined_barriers) != "EPSG:27700") {
  combined_barriers <- project(combined_barriers, "EPSG:27700")
}

# Initial filtering
combined_barriers <- combined_barriers %>%
  select(form, `function`, L_SYSTEM)

# Add "type" column based on conditions
combined_barriers <- combined_barriers %>%
  mutate(type = case_when(
    form == "canal" ~ "canal",
    form == "river" ~ "river",
    `function` %in% c("Motorway", "A Road") ~ "road",
    L_SYSTEM %in% c("M", "K") ~ "rail",
    TRUE ~ NA_character_
  ))

# -----------------------
# Add coastline as barrier
# -----------------------
boundary_line <- as.lines(boundary_vect)
boundary_line$type <- "coast"

# Combine coastline with other barriers
barriers_all <- rbind(combined_barriers, boundary_line)

# Save combined barriers
writeVector(barriers_all, file.path(data_dir, "Combined Dispersal Barriers.geojson"), overwrite = TRUE)

# -----------------------
# Define colors
# -----------------------
barrier_colors <- c(
  river = "#0571b0",
  road  = "#ca0020",
  rail  = "#7b3294",
  canal = "#e66101",
  coast = "gray30"
)

# -----------------------
# Plot combined barriers
# -----------------------
jpeg(filename = file.path(plots_path, "Combined Dispersal Barriers.jpg"),
     width = 20, height = 15, units = "cm", res = 300)

par(mar = c(3, 1, 3, 1) + 0.1)

# Plot coastline/boundary outline first
plot(boundary_line, col = "gray30", lwd = 0.75, axes = FALSE, main = "")

# Plot each barrier type
types_present <- intersect(names(barrier_colors), unique(barriers_all$type))
types_present <- types_present[!is.na(types_present)]

if (length(types_present) > 0) {
  first_type <- types_present[1]
  
  plot(barriers_all[barriers_all$type == first_type, ],
       col = barrier_colors[first_type],
       lwd = 0.75,
       add = TRUE)
  
  if (length(types_present) > 1) {
    for (type_name in types_present[-1]) {
      subset_geom <- barriers_all[barriers_all$type == type_name, ]
      if (nrow(subset_geom) > 0) {
        plot(subset_geom,
             col = barrier_colors[type_name],
             add = TRUE,
             lwd = 0.75)
      }
    }
  }
}

mtext("Dispersal Barriers Map", side = 3, outer = FALSE, cex = 1.5, line = 1)

legend("topright",
       legend = c("Coastline", "Rivers", "Roads", "Railways", "Canals"),
       fill   = c("gray30", "#0571b0", "#ca0020", "#7b3294", "#e66101"),
       title  = "Barrier Types",
       cex    = 0.8,
       bty    = "n")

dev.off()
