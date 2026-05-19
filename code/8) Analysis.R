# Clean session (optional)
# rm(list = ls())

# Load packages
library(terra)
library(readxl)
library(ggplot2)
library(dplyr)
library(writexl)

###################################################
# Summary Table Combined Stats FOR SPLIT WITH BARRIERS#
###################################################

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
      domination_index = 1 - (max(area_m2, na.rm = TRUE) / sum(area_m2, na.rm = TRUE))
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

# -----------------------
# Settings / I-O paths
# -----------------------
in_dir        <- "./data/"
out_dir       <- "./results/analysis/"
terra_tmp_dir <- "./_terra_tmp"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(terra_tmp_dir, showWarnings = FALSE, recursive = TRUE)

terraOptions(
  memfrac = 0.7,
  tempdir = terra_tmp_dir
)

# -----------------------
# Species parameters
# -----------------------
min_area_path  <- file.path(in_dir, "Species Minimum Area Required.xlsx")
dispersal_path <- file.path(in_dir, "Species Dispersal Distances.xlsx")

min_area_df <- read_excel(min_area_path)
min_area_m2 <- setNames(min_area_df$`min_area (m2)`, min_area_df$species)
rm(min_area_df)

dispersal_df <- read_excel(dispersal_path)
dispersal_m  <- setNames(dispersal_df$dispersal_dist_m, dispersal_df$species)
rm(dispersal_df)

# -----------------------
# Build analysis dataframe (JOIN BY SPECIES NAME)
# -----------------------
common_species <- intersect(names(min_area_m2), names(dispersal_m))

df <- data.frame(
  species     = common_species,
  min_area_m2 = as.numeric(min_area_m2[common_species]),
  dispersal_m = as.numeric(dispersal_m[common_species])
)

# Drop missing / non-positive values (log10 requires > 0)
df <- df[is.finite(df$min_area_m2) & is.finite(df$dispersal_m) &
           df$min_area_m2 > 0 & df$dispersal_m > 0, ]

# Unit conversions for plotting and interpretability
df$min_area_km2 <- df$min_area_m2 / 1e6
df$dispersal_km <- df$dispersal_m / 1000

# -----------------------
# Statistics: correlation + regression in log-log space
# -----------------------
cor_test <- cor.test(
  log10(df$dispersal_km),
  log10(df$min_area_km2),
  method = "pearson"
)

cor_text <- paste0(
  "log-log r = ", round(unname(cor_test$estimate), 2),
  ", p = ", signif(cor_test$p.value, 2)
)

fit <- lm(log10(min_area_km2) ~ log10(dispersal_km), data = df)
r2_text <- paste0("R² = ", round(summary(fit)$r.squared, 2))

# Optional: print full model output to console
print(cor_test)
print(summary(fit))

# -----------------------
# Label nudges (edit species names to match exactly)
# -----------------------
vjust_vec <- ifelse(df$species %in% c("Pine_Marten", "Lynx", "Elk"), 1.5, -0.8)

# -----------------------
# Plot
# -----------------------
p <- ggplot(df, aes(x = dispersal_km, y = min_area_km2, label = species)) +
  geom_point(size = 3, color = "forestgreen") +
  geom_text(vjust = vjust_vec, size = 2.75) +
  scale_y_log10() +
  scale_x_log10() +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  labs(
    x = "Dispersal Distance (km, log scale)",
    y = "Minimum Area Required (km², log scale)",
    title = "Scaling Relationship between Dispersal Distance and Minimum Area",
    subtitle = paste(cor_text, "|", r2_text)
  ) +
  theme_minimal()

print(p)

# -----------------------
# Save plot (optional)
# -----------------------
ggsave(
  filename = file.path(out_dir, "dispersal_vs_min_area_loglog.png"),
  plot = p,
  width = 7, height = 5, dpi = 300
)
