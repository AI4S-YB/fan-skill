# theme_plant_scientific.R
# A clean, publication-ready ggplot2 theme for plant biology
# Based on Tufte principles: maximize data-ink ratio

library(ggplot2)

theme_plant_scientific <- function(base_size = 8, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      # Remove top and right borders (Tufte data-ink ratio)
      panel.border = element_rect(fill = NA, color = "black", size = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(size = 0.2, color = "grey90"),

      # Clean axis lines
      axis.line = element_blank(),
      axis.ticks = element_line(size = 0.3),
      axis.ticks.length = unit(0.05, "cm"),

      # Typography hierarchy
      axis.title = element_text(size = base_size + 2, face = "plain"),
      axis.text = element_text(size = base_size),
      strip.text = element_text(size = base_size + 1, face = "bold"),
      strip.background = element_rect(fill = "grey95", color = NA),
      legend.text = element_text(size = base_size),
      legend.title = element_text(size = base_size + 1),
      legend.key = element_rect(fill = NA, color = NA),
      legend.key.size = unit(0.4, "cm"),

      # Faceting
      panel.spacing = unit(0.8, "lines"),

      # Default aspect ratio (can be overridden per-chart)
      aspect.ratio = 0.75
    )
}

# Colorblind-safe palette (Okabe-Ito, 8 colors)
# Suitable for up to 8 groups; for more, use viridis
okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000"
)

# Viridis discrete palette (reduced to 6 distinct colors for clarity)
viridis_d <- c("#440154", "#414487", "#2A788E", "#22A884", "#7AD151", "#FDE725")

# Plant-specific color palettes

# Green-brown palette for field/soil/plant contexts
plant_field <- c(
  "#228B22",  # ForestGreen (healthy plant)
  "#8B4513",  # SaddleBrown (soil)
  "#6B8E23",  # OliveDrab (mature/stressed plant)
  "#DAA520",  # Goldenrod (harvest/yield)
  "#556B2F",  # DarkOliveGreen (canopy)
  "#CD853F"   # Peru (dry soil)
)

# Photosynthesis / light-response palette (green-to-yellow gradient)
plant_photo <- c(
  "#00441B", "#006D2C", "#238B45", "#41AB5D",
  "#74C476", "#A1D99B", "#C7E9C0", "#F7FCF5"
)

# Drought/heat stress palette (green to brown to red)
plant_stress <- c(
  "#1B9E77", "#66A61E", "#E6AB02", "#D95F02", "#A6761D", "#7570B3"
)
