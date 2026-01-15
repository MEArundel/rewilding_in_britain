#################################################
## Preamble ##
#################################################

# Clear Memory
rm(list = ls())

# Clear Memory EXCEPT TEST AREAS
#rm(list = setdiff(ls(), c("test_region", "regions")))

# Load Packages
library(ggplot2)

# Setting seed
set.seed(13) #ensures that any randomness is controlled between MAC and Windows and between sessions

# Data
df <- data.frame(
  species = c("Beaver", "Bison", "Brown Bear", "Elk", "Lynx",
              "Pine Marten", "Wild Boar", "Wildcat", "Wolf", "Wolverine"),
  min_area = c(1.26, 212.096774, 727.272727, 7102.65, 847.457627,
               61.6, 70.92, 518, 37825.498, 11545.54),
  dispersal_dist_km = c(5, 16, 73.5, 44.5, 32.5, 11.35, 10.55, 10, 65.8, 55.5)
)

# Compute log-log correlation
cor_val <- cor(log10(df$dispersal_dist_km), log10(df$min_area))
cor_text <- paste0("log-log r = ", round(cor_val, 2))

# Fit log-log linear model
fit <- lm(log10(min_area) ~ log10(dispersal_dist_km), data = df)
r2_val <- summary(fit)$r.squared
r2_text <- paste0("R² = ", round(r2_val, 2))

# Create vjust vector for labels
vjust_vec <- ifelse(df$species %in% c("Pine Marten", "Lynx"), 1.5, -0.8)

# Plot
ggplot(df, aes(x = dispersal_dist_km, y = min_area, label = species)) +
  geom_point(size = 3, color = "forestgreen") +
  geom_text(vjust = vjust_vec, size = 3) +
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
