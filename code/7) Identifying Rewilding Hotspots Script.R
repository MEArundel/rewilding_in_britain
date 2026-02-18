###################################################
# Identify rewilding hotspots (barrier polygons)  #
###################################################
# Description - This script identifies multi-species rewilding hotspots by
# partitioning Great Britain into barrier-defined polygons and quantifying how
# many rewilding species can be supported within each polygon. First, a dissolved
# Britain boundary is split progressively using dispersal barrier line layers
# (roads, rivers, canals, railways, plus coastline) to generate a set of barrier
# polygons with unique IDs and areas. Each barrier polygon is then assigned to
# England, Scotland, or Wales based on maximum spatial overlap. Next, for each
# species’ final patch-ID raster, patches are converted to polygons and intersected
# with the barrier polygons to calculate the proportional overlap of each patch
# with each polygon; a polygon is considered to support a species if at least one
# patch overlaps the polygon by a specified threshold (default = 0.95). The script
# outputs long and wide summary tables of species support per barrier polygon,
# ranks polygons supporting ≥7 species (prioritising more species, then smaller
# polygon area), and produces a publication-ready ranked hotspot map highlighting
# primary (8-species) and secondary (7-species) hotspot polygons.

# Clean session (optional)
# rm(list = ls())

# Packages
library(terra)
library(dplyr)
library(tidyr)
library(writexl)

# -----------------------
# Helpers
# -----------------------
flatten_to_vect <- function(x) {
  if (inherits(x, "SpatVector")) return(x)
  
  if (inherits(x, "SpatVectorCollection")) {
    xs <- as.list(x)
    xs <- xs[!vapply(xs, is.null, logical(1))]
    if (length(xs) == 0) return(NULL)
    return(do.call(rbind, xs))
  }
  
  stop("Unknown object type: ", paste(class(x), collapse = ", "))
}

split_one_poly <- function(poly, splitter) {
  out <- tryCatch(terra::split(poly, splitter), error = function(e) poly)
  flatten_to_vect(out)
}

assert_has_field <- function(v, field) {
  if (!(field %in% names(v))) {
    stop("Missing field '", field, "' in object. Fields are: ", paste(names(v), collapse = ", "))
  }
}

# -----------------------
# Settings / I/O paths
# -----------------------
set.seed(13)

# Fix 1: gate debug/interactive plots
do_checks <- FALSE

tmp_dir <- "./_terra_tmp"
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
terraOptions(memfrac = 0.5, tempdir = tmp_dir, progress = 1)

path_barriers <- "./data/Combined Dispersal Barriers.geojson"        # barrier LINES with a 'type' column
path_boundary <- "./data/Britain Boundary Polygon.geojson"           # England/Scotland/Wales polygons with 'shape1'
species_dir   <- "./data/Suitable Habitats with Barriers/"           # contains *patch-ID* rasters
out_dir       <- "./results/Rewilding Hotspots/Barrier Polygon Summary/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

crs_target <- "EPSG:27700"

# Threshold: patch considered "within" a barrier polygon
overlap_threshold <- 0.95
message("Using overlap_threshold = ", overlap_threshold)

# Fix 2: input existence checks early
stopifnot(file.exists(path_barriers))
stopifnot(file.exists(path_boundary))
stopifnot(dir.exists(species_dir))

# -----------------------
# Load boundary + regions
# -----------------------
boundary_vect <- terra::vect(path_boundary)
boundary_vect <- terra::project(boundary_vect, crs_target)
boundary_vect <- terra::makeValid(boundary_vect)

assert_has_field(boundary_vect, "shape1")

regions <- list(
  England  = boundary_vect[boundary_vect$shape1 == "England", ],
  Scotland = boundary_vect[boundary_vect$shape1 == "Scotland", ],
  Wales    = boundary_vect[boundary_vect$shape1 == "Wales", ]
)

# One Britain polygon (dissolve internal borders)
boundary_union <- terra::aggregate(boundary_vect)
boundary_union <- terra::makeValid(boundary_union)

# -----------------------
# Load barrier LINES and build barrier POLYGONS by splitting Britain
# -----------------------
barrier_lines_raw <- terra::vect(path_barriers)
barrier_lines_raw <- terra::project(barrier_lines_raw, crs_target)
barrier_lines_raw <- terra::makeValid(barrier_lines_raw)

assert_has_field(barrier_lines_raw, "type")

# Add coastline line so faces close nicely
coast_line <- terra::as.lines(boundary_union)
coast_line$type <- "coast"

barrier_lines <- rbind(barrier_lines_raw, coast_line)
barrier_lines <- terra::makeValid(barrier_lines)

# Start with the Britain polygon, then split progressively
barrier_polys <- boundary_union

# Split by each non-rail barrier type
non_rail_types <- setdiff(unique(barrier_lines$type), "rail")

for (bt in non_rail_types) {
  message("Splitting by: ", bt)
  b <- barrier_lines[barrier_lines$type == bt, ]
  if (nrow(b) == 0) next
  b <- terra::makeValid(b)
  
  barrier_polys <- flatten_to_vect(barrier_polys)
  new_polys <- vector("list", nrow(barrier_polys))
  
  for (i in seq_len(nrow(barrier_polys))) {
    new_polys[[i]] <- split_one_poly(barrier_polys[i], b)
  }
  barrier_polys <- do.call(rbind, new_polys)
}

# Rail split using a small buffer
rail_lines <- barrier_lines[barrier_lines$type == "rail", ]
if (nrow(rail_lines) > 0) {
  message("Splitting by: rail (buffered)")
  rail_lines <- terra::makeValid(rail_lines)
  rail_buf   <- terra::buffer(rail_lines, width = 25)
  
  barrier_polys <- flatten_to_vect(barrier_polys)
  new_polys <- vector("list", nrow(barrier_polys))
  
  for (i in seq_len(nrow(barrier_polys))) {
    new_polys[[i]] <- split_one_poly(barrier_polys[i], rail_buf)
  }
  barrier_polys <- do.call(rbind, new_polys)
}

# Finalize barrier polygons
barrier_polys <- terra::makeValid(barrier_polys)

message("Barrier polys geomtype: ", terra::geomtype(barrier_polys))
if (terra::geomtype(barrier_polys) != "polygons") {
  stop("barrier_polys is not polygons (geomtype = ", terra::geomtype(barrier_polys), "). Check your barrier input.")
}

# IDs + areas
barrier_polys$barrier_poly_id <- seq_len(nrow(barrier_polys))
barrier_polys$area_m2  <- terra::expanse(barrier_polys, unit = "m")
barrier_polys$area_km2 <- barrier_polys$area_m2 / 1e6

message("Barrier polygons created: ", nrow(barrier_polys))
summary(barrier_polys$area_km2)

# Fix 1: optional check plotting (gated)
if (do_checks) {
  plot(barrier_polys)
  plot(boundary_union, add = TRUE, border = "black", lwd = 1.5)
}

# -----------------------
# Assign each barrier polygon to a single country (largest overlap)
# -----------------------
country_assignments <- lapply(names(regions), function(ctry) {
  
  inter <- terra::intersect(barrier_polys, regions[[ctry]])
  if (is.null(inter) || nrow(inter) == 0) return(NULL)
  
  if (!("barrier_poly_id" %in% names(inter))) {
    stop("barrier_poly_id not present after intersect for ", ctry, ".")
  }
  
  inter$country <- ctry
  inter$overlap_km2 <- terra::expanse(inter, unit = "km")
  
  inter[, c("barrier_poly_id", "country", "overlap_km2")]
})

country_assignments <- do.call(rbind, country_assignments)

country_df <- as.data.frame(country_assignments) %>%
  group_by(barrier_poly_id) %>%
  slice_max(overlap_km2, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(barrier_poly_id, country)

# -----------------------
# Find final species rasters
# -----------------------
# Fix 5: deterministic ordering
species_tifs <- sort(list.files(
  species_dir,
  pattern = "Suitable Habitat Split and Filtered\\.tif$",
  full.names = TRUE
))

if (length(species_tifs) == 0) {
  stop("No final species rasters found in ", species_dir)
}

# -----------------------
# Patch-level species support inside each barrier polygon
# -----------------------
support_tables <- list()

for (tif_path in species_tifs) {
  
  sp <- gsub(" Suitable Habitat Split and Filtered\\.tif$", "", basename(tif_path))
  sp <- gsub("_", " ", sp)
  message("\n--- Species: ", sp, " ---")
  
  r <- terra::rast(tif_path)
  
  # CRS ensure (robust)
  if (is.na(terra::crs(r))) stop("Raster has NA CRS: ", tif_path)
  if (!terra::same.crs(r, crs_target)) {
    r <- terra::project(r, crs_target, method = "near")
  }
  
  # Skip if empty
  if (terra::global(!is.na(r), "sum", na.rm = TRUE)[1, 1] == 0) {
    message("  (no patches; skipping)")
    next
  }
  
  # Raster -> polygons, dissolve by patch id
  patch_polys <- terra::as.polygons(r, dissolve = TRUE, na.rm = TRUE)
  names(patch_polys) <- "patch_id"
  patch_polys <- terra::makeValid(patch_polys)
  
  # Patch total area
  patch_polys$patch_area_m2 <- terra::expanse(patch_polys, unit = "m")
  
  # Intersect patches with barrier polygons
  pieces <- terra::intersect(patch_polys, barrier_polys)
  
  if (is.null(pieces) || nrow(pieces) == 0) {
    message("  (no intersections with barrier polygons; skipping)")
    rm(r, patch_polys, pieces); gc()
    next
  }
  
  if (!("barrier_poly_id" %in% names(pieces))) stop("barrier_poly_id missing in pieces for species: ", sp)
  if (!("patch_id" %in% names(pieces))) stop("patch_id missing in pieces for species: ", sp)
  if (!("patch_area_m2" %in% names(pieces))) stop("patch_area_m2 missing in pieces for species: ", sp)
  
  pieces <- terra::makeValid(pieces)
  pieces$piece_area_m2 <- terra::expanse(pieces, unit = "m")
  
  pieces_df <- as.data.frame(pieces)[, c("barrier_poly_id", "patch_id", "patch_area_m2", "piece_area_m2")]
  
  patch_in_barrier <- pieces_df %>%
    group_by(barrier_poly_id, patch_id) %>%
    summarise(
      patch_area_m2   = first(patch_area_m2),
      overlap_area_m2 = sum(piece_area_m2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Fix 3: safe overlap proportion (avoid divide-by-zero)
    mutate(overlap_prop = ifelse(patch_area_m2 > 0, overlap_area_m2 / patch_area_m2, NA_real_))
  
  supported_by_poly <- patch_in_barrier %>%
    group_by(barrier_poly_id) %>%
    summarise(
      species = sp,
      supports_species  = any(overlap_prop >= overlap_threshold, na.rm = TRUE),
      best_overlap_prop = max(overlap_prop, na.rm = TRUE),
      best_patch_id     = patch_id[which.max(overlap_prop)][1],
      .groups = "drop"
    )
  
  support_tables[[sp]] <- supported_by_poly
  
  rm(r, patch_polys, pieces, pieces_df, patch_in_barrier, supported_by_poly)
  gc()
}

# Fix 4: guard against empty outputs
if (length(support_tables) == 0) {
  stop("No species produced support results; check your inputs (rasters/barriers) and overlap_threshold.")
}

# Combine all species results into one long table
support_long <- bind_rows(support_tables)

# Wide table: one row per barrier polygon, columns per species TRUE/FALSE
support_wide <- support_long %>%
  select(barrier_poly_id, species, supports_species) %>%
  tidyr::pivot_wider(
    names_from  = species,
    values_from = supports_species,
    values_fill = FALSE
  )

# Ensure all barrier polygons represented
support_wide <- as.data.frame(barrier_polys)[, "barrier_poly_id", drop = FALSE] %>%
  left_join(support_wide, by = "barrier_poly_id") %>%
  mutate(across(where(is.logical), ~ tidyr::replace_na(.x, FALSE)))

# Add totals per barrier polygon
hotspot_summary <- support_wide %>%
  mutate(n_species_supported = rowSums(across(where(is.logical)))) %>%
  left_join(
    as.data.frame(barrier_polys)[, c("barrier_poly_id", "area_km2")],
    by = "barrier_poly_id"
  ) %>%
  left_join(country_df, by = "barrier_poly_id") %>%
  relocate(barrier_poly_id, country, area_km2, n_species_supported)

# -----------------------
# Save outputs
# -----------------------
write.csv(support_long,
          file.path(out_dir, "BarrierPolygon_SpeciesSupport_LONG.csv"),
          row.names = FALSE)

write.csv(hotspot_summary,
          file.path(out_dir, "BarrierPolygon_HotspotSummary_WIDE.csv"),
          row.names = FALSE)

hotspot_summary_nonzero <- hotspot_summary[hotspot_summary$n_species_supported >= 1, ]
write.csv(hotspot_summary_nonzero,
          file.path(out_dir, "BarrierPolygon_HotspotSummary_WIDE_nonzero.csv"),
          row.names = FALSE)

# -----------------------
# Filter: keep only patches supporting >= 7 species
# -----------------------
hotspot_7plus <- hotspot_summary[hotspot_summary$n_species_supported >= 7, , drop = FALSE]

# -----------------------
# Rank: more species is better; for ties, smaller area is better
# -----------------------
hotspot_7plus <- hotspot_7plus[order(-hotspot_7plus$n_species_supported, hotspot_7plus$area_km2), , drop = FALSE]
hotspot_7plus$rank_best <- seq_len(nrow(hotspot_7plus))

write.csv(hotspot_7plus,
          file.path(out_dir, "Rewilding_Hotspots.csv"),
          row.names = FALSE)

message("\n✅ Done. Output written to: ", out_dir)

# --- Subset polygons to only ranked hotspots (use barrier_poly_id) ---
hotspot_polys <- barrier_polys[barrier_polys$barrier_poly_id %in% hotspot_7plus$barrier_poly_id, ]

if (is.null(hotspot_polys) || nrow(hotspot_polys) == 0) {
  stop("No hotspot polygons matched barrier_poly_id. Check IDs / filtering.")
}

# Match order to hotspot_7plus and attach attributes
ord <- match(hotspot_polys$barrier_poly_id, hotspot_7plus$barrier_poly_id)
hotspot_polys$rank_best           <- hotspot_7plus$rank_best[ord]
hotspot_polys$n_species_supported <- hotspot_7plus$n_species_supported[ord]
hotspot_polys$rank_best <- as.integer(hotspot_polys$rank_best)

# Split by species class
is_primary   <- hotspot_polys$n_species_supported == 8
is_secondary <- hotspot_polys$n_species_supported == 7

# Palettes (dark = best rank, light = worse rank)
cols_primary   <- colorRampPalette(c("#3182bd", "#9ecae1"))(sum(is_primary))   # blue
cols_secondary <- colorRampPalette(c("#de2d26", "#fc9272"))(sum(is_secondary)) # orange

hotspot_polys$col <- NA_character_

if (any(is_primary)) {
  idx <- which(is_primary)
  hotspot_polys$col[idx] <- cols_primary[
    rank(hotspot_polys$rank_best[idx], ties.method = "first")
  ]
}

if (any(is_secondary)) {
  idx <- which(is_secondary)
  hotspot_polys$col[idx] <- cols_secondary[
    rank(hotspot_polys$rank_best[idx], ties.method = "first")
  ]
}

# ---------- SAVE AS JPEG ----------
jpeg(
  filename = file.path(out_dir, "Rewilding_Hotspots_Ranked.jpg"),
  width    = 3000,
  height   = 2400,
  res      = 300
)

par(mar = c(6, 1.5, 1.5, 1.5))

ext  <- terra::ext(boundary_union)
dx   <- ext$xmax - ext$xmin
dy   <- ext$ymax - ext$ymin
xpad <- 0.02 * dx
ypad <- 0.02 * dy

plot(
  boundary_union,
  col        = NA,
  border     = "black",
  lwd        = 1.5,
  xlim       = c(ext$xmin - xpad, ext$xmax + xpad),
  ylim       = c(ext$ymin - ypad, ext$ymax + ypad),
  axes       = FALSE,
  frame.plot = FALSE
)

plot(
  hotspot_polys,
  col    = hotspot_polys$col,
  border = "grey40",
  add    = TRUE
)

cc <- terra::centroids(hotspot_polys)
text(cc, labels = hotspot_polys$rank_best, cex = 0.85, col = "black")

par(xpd = NA)

legend(
  x      = ext$xmin + 0.65 * dx,
  y      = ext$ymax - 0.12 * dy,
  legend = c("Primary hotspots", "Secondary hotspots"),
  fill   = c("#3182bd", "#de2d26"),
  border = NA,
  bty    = "n",
  cex    = 0.95
)

mtext(
  "Darker colour = higher rank (more species, then smaller area)",
  side = 1,
  line = 3,
  cex  = 0.95
)

dev.off()
