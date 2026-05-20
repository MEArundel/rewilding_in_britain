**Rewilding Hotspots as Spatial Priorities for Trophic Restoration**

This repository contains the full spatial modelling and analysis workflow used to identify potential rewilding hotspots across Great Britain. Using species distribution modelling, habitat suitability mapping, dispersal buffering, and barrier-constrained patch analysis, the project evaluates where 
large-scale trophic rewilding and species reintroductions may be ecologically viable. The repository includes scripts for habitat processing, connectivity modelling, hotspot identification, visualisation, and statistical analysis for multiple focal rewilding species.

Code Descriptions:

Code 1) Land Coverage Reclassification Script
Reclassification of the UKCEH 21 land coverage types into 13 simplified types
Description - This script processes a national land cover raster for Britain by masking it to a dissolved England–Scotland–Wales boundary, aggregating the raster to a lower spatial resolution using a custom modal function, and reclassifying 21 original land cover classes into 13 broader categories. 
The workflow outputs a reclassified raster and a static JPEG map.

Code 2) Hs Code 
Formation of habitat suitability maps (to create Habitat Fragments (HFs))
Description - This script generates binary habitat suitability maps (referred to as HFs) for focal species in Britain by reclassifying a national land cover raster using species-specific suitability rules. For each species, land cover classes are converted to suitable (1) or unsuitable (0) habitat 
based on a reference table, masked to a dissolved England–Scotland–Wales boundary, and exported as both a GeoTIFF and a static JPEG map. The script additionally compiles all species-specific suitability rasters into a single multi-panel figure for comparative visualisation.

Code 3) Buffering Hs Maps Code
Assigning species dispersal distance as a buffer (to create Dispersal Constrained Patches (DCPs))
Description - This script expands species-specific suitable habitat maps by applying dispersal-distance buffers to account for landscape connectivity and potential movement. For each species, binary habitat suitability rasters are buffered by species-specific dispersal distances, processed separately 
within England, Scotland, and Wales, and mosaicked into a single Great Britain–wide raster. Resulting habitat patches are then clustered and filtered to retain only those exceeding species-specific minimum area requirements. The script exports buffered and filtered habitat rasters as GeoTIFFs, 
generates individual JPEG maps for each species, and produces a combined multi-panel figure for comparative visualisation.

Code 4) Processing and Combining Dispersal Barriers Script
Dispersal barriers: build + combine (GB)
Description - This script (1) identifies and exports four dispersal barrier layers (major roads, canals, railways, wide rivers) from source datasets, produces a JPEG map for each barrier type, and (2) combines all barrier layers (plus coastline) into a single GeoJSON and a combined JPEG map.

Code 5) Splitting Centroid Patches with Barriers
Splitting buffered suitable habitat patches by barriers (to create Barrier Constrained Patches (BCPs))
Description - This script refines species-specific buffered habitat maps by explicitly accounting for dispersal barriers and minimum patch size constraints. For each species, buffered suitable habitat rasters are converted to polygons, convex hulls are generated for contiguous patches, and these 
patches are split using multiple barrier types (roads, rivers, canals, railways, and coastline). Resulting habitat fragments are rasterised, reclassified, and filtered to retain only patches exceeding species-specific minimum area requirements. The workflow outputs final filtered habitat rasters, 
spatial patch layers, patch-level and species-level summary tables, and publication-ready maps with dispersal barriers overlaid.

Code 6) Heat Map and Visulisation Script
Species habitat maps, patch maps, and overlap heatmap (GB) (to create Viable Patches (VPs))
Description - This script visualises final species habitat outputs by (1) plotting binary suitable habitat maps for all species in a multi-panel figure, (2) plotting patch-ID maps (one panel per species) with dispersal barriers overlaid, and (3) computing and mapping the number of species overlapping 
per pixel across Great Britain as a heatmap. Outputs are saved as publication-ready JPEGs.

Code 7) Identifying Rewilding Hotspots Script
Identify rewilding hotspots (barrier polygons)
Description - This script identifies multi-species rewilding hotspots by partitioning Great Britain into barrier-defined polygons and quantifying how many rewilding species can be supported within each polygon. First, a dissolved Britain boundary is split progressively using dispersal barrier line 
layers (roads, rivers, canals, railways, plus coastline) to generate a set of barrier polygons with unique IDs and areas. Each barrier polygon is then assigned to England, Scotland, or Wales based on maximum spatial overlap. Next, for each species’ final patch-ID raster, patches are converted to 
polygons and intersected with the barrier polygons to calculate the proportional overlap of each patch with each polygon; a polygon is considered to support a species if at least one patch overlaps the polygon by a specified threshold (default = 0.95). The script outputs long and wide summary tables 
of species support per barrier polygon, ranks polygons supporting ≥7 species (prioritising more species, then smaller polygon area), and produces a publication-ready ranked hotspot map highlighting primary (8-species) and secondary (7-species) hotspot polygons.

Code 8) Analysis
Summary table combined stats (split with barriers)
Description - This script summarises species-level habitat patch statistics from Excel outputs generated after suitable habitat has been split by dispersal barriers. For each species file, the script calculates the number of remaining habitat patches, total suitable habitat area, mean, median and 
standard deviation of patch area, largest patch size, the percentage of total habitat contained within the largest patch, and a domination index describing how evenly habitat area is distributed among patches. Empty or invalid files are skipped, and all valid species summaries are combined, sorted 
alphabetically, and exported as a single overview Excel table.
AND
Dispersal distance and minimum area relationship
Description - This script analyses the scaling relationship between species dispersal distance and the minimum habitat area required to support a viable population. Species-level minimum area and dispersal distance tables are joined by species name, filtered to remove missing or invalid values, and
converted into km² and km for interpretation. A Pearson correlation test and linear regression model are then fitted in log-log space to assess whether species with greater dispersal distances also require larger habitat areas. The script produces a labelled scatterplot with a fitted regression line 
and saves the figure as a high-resolution PNG.


Download Locations:

Spatial Downloads:
Land Coverage (UKCEH) Map -  Morton, R.D., Marston, C.G., O’Neil, A.W., Rowland, C.S., 2024. Land Cover Map 2023 (land parcels, GB). https://doi.org/10.5285/50B344EB-8343-423B-8B2F-0E9800E34BBD
The cleaning and processing of these data can be found in code 1) Land Coverage Reclassification Script
European Rivers - Pope, A., 2017. European River Data. https://doi.org/10.7488/DS/1886
British Railways - Pope, A., 2017a. GB Railways and stations. https://doi.org/10.7488/DS/1773
British Roads - Ordnance Survey, 2024. OS Open Roads
The cleaning and processing of these data can be found in code 4) Processing and Combining Dispersal Barriers Script

Species Specific Downloads:
Minimum Viable Populations - Traill, L.W., Bradshaw, C.J.A., Brook, B.W., 2007. Minimum viable population size: A meta-analysis of 30 years of published estimates. Biological Conservation 139, 159–166. https://doi.org/10.1016/j.biocon.2007.06.011
Population Densities - Santini, L., Benítez‐López, A., Dormann, C.F., Huijbregts, M.A.J., 2022. Population density estimates for terrestrial mammal species. Global Ecol Biogeogr 31, 978–994. https://doi.org/10.1111/geb.13476
Home Range - Broekman, M.J.E., Hilbers, J.P., Hoeks, S., Huijbregts, M.A.J., Schipper, A.M., Tucker, M.A., 2024. Environmental drivers of global variation in home range size of terrestrial and marine mammals. Journal of Animal Ecology 93, 488–500. https://doi.org/10.1111/1365-2656.14073


