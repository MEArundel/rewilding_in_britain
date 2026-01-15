####################################################################################
# Reclassification of the UKCEH 21 land coverage types into to 13 simplified types #
####################################################################################
# Description - This script processes a national land cover raster for Britain by
# masking it to a dissolved England–Scotland–Wales boundary, aggregating the raster
# to a lower spatial resolution using a custom modal function, and reclassifying
# 21 original land cover classes into 13 broader categories. The workflow outputs
# a reclassified raster and a static JPEG map.

# Clean session (optional)
# rm(list = ls())

# Load packages
library(terra)

# -------------------
# Settings / I-O paths
# -------------------
set.seed(13) # included for consistency across scripts

in_dir <- "./data/"
out_dir_data <- "./data/"
out_dir_res <- "./results/Land Coverage Map/"

# Ensure output directory exists
dir.create(out_dir_res, recursive = TRUE, showWarnings = FALSE)

# -----------------------
# Inputs and checks
# -----------------------
lc_path <- file.path(in_dir, "Land Coverage Type.tif")
bd_path <- file.path(in_dir, "Britain Boundary Polygon.geojson")

stopifnot(file.exists(lc_path))
stopifnot(file.exists(bd_path))

# Load land coverage data
land_cover <- rast(file.path(in_dir, "Land Coverage Type.tif"))
land_cover <- land_cover[[1]]

# ------------
# Boundary
# ------------
boundary_vect <- vect(file.path(in_dir, "Britain Boundary Polygon.geojson"))

# Project boundary to the raster CRS
boundary_vect <- project(boundary_vect, crs(land_cover))

# Dissolve internal borders (creates one unified polygon)
boundary_vect <- aggregate(boundary_vect)

# -----------------------
# Masking and aggregation
# -----------------------
# Mask land_cover
land_cover_masked <- mask(land_cover, boundary_vect)

# Written function to replace buggy terra::modal
fmod <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  uv <- unique(x)
  uv[which.max(tabulate(match(x, uv)))]
}

# Aggregate to lower resolution
land_cover_agg <- aggregate(land_cover_masked, fact = 5, fun = fmod)
# fact = number of pixels (in each direction) combined into one output cell

# ------------------
# Reclassification
# ------------------
reclass_matrix <- matrix(
  c(0, 0,
    1, 1,
    2, 1,
    3, 2,
    4, 3,
    5, 3,
    6, 4,
    7, 5,
    8, 6,
    9, 7,
    10, 7,
    11, 6,
    12, 8,
    13, 9,
    14, 10,
    15, 11,
    16, 11,
    17, 12,
    18, 12,
    19, 12,
    20, 13,
    21, 13),
  ncol = 2, byrow = TRUE
)

land_cover_reclassified <- classify(land_cover_agg, reclass_matrix)

# Optional sanity check: confirm output values are within expected range
# print(sort(unique(values(land_cover_reclassified))))

# ---------------------
# Plotting and saving
# ---------------------
palette <- c("white", "limegreen", "yellow2", "darkorange",
             "goldenrod1", "darkgoldenrod3", "maroon1", "red",
             "azure2", "navyblue", "blue", "plum",
             "brown", "black")

labels <- c("NA", "Woodland", "Arable and Horticulture",
            "Improved and Neutral Grassland", "Calcareous Grassland", "Acid Grassland",
            "Fen, Marsh, Swamp and Bog", "Heather Grassland", "Inland Rock", "Saltwater",
            "Freshwater", "Supralittoral Rock and Sediment",
            "Littoral Rock and Sediment and Saltmarsh", "Urban and Suburban")

# Save JPEG map
jpeg(
  filename = file.path(out_dir_res, "Land Coverage Map Reclassified.jpg"),
  width = 20, height = 15, units = "cm", res = 300
)

layout(matrix(c(1, 2), nrow = 1), widths = c(1, 1))

par(mar = c(5, 4, 4, 1) + 0.1)
plot(land_cover_reclassified, col = palette, legend = FALSE, axes = FALSE,
     main = "Land Cover in Britain")

par(mar = c(5, 0, 4, 2) + 0.1)
plot.new()
legend("center", legend = labels, fill = palette, title = "Land Cover Types", cex = 0.8)

dev.off()

# Save reclassified raster
writeRaster(
  land_cover_reclassified,
  filename = file.path(out_dir_data, "Land Coverage Map Reclassified Raster.tif"),
  overwrite = TRUE
)
