# Script to create Figure 1 Australia sample map

library(tidyverse)
library(ggrepel)
library(cowplot)
library(rnaturalearth)
library(rnaturalearthdata)
library(grid)

setwd("analyses/data/")
locations <- read.csv("sampling_geo_data.csv")

pie_data <- locations %>%
  group_by(site, lat, long, animal) %>%
  mutate(total_samples = sum(num_samples)) %>%
  spread(key = sample_type, value = num_samples, fill = 0) %>%
  ungroup()

# Define colors for sample_type
sample_colors <- c(
  "faecal" = "#6C73B1", 
  "rectal_swab" = "#A5C2CA", 
  "nasal_swab" = "#B2BA72", 
  "tissue" = "#C7A1BC"
)

# Dummy dataset for color legend
dummy_data <- data.frame(
  sample_type = names(sample_colors),
  x = 1:length(sample_colors),
  y = 1
)

# Create a dataset for the size legend
size_legend_data <- data.frame(
  x = 1:3,
  y = 1,
  total_samples = c(min(pie_data$total_samples), median(pie_data$total_samples), max(pie_data$total_samples))
)

# Create a pie chart grob for each site
create_pie_grob <- function(data_row) {
  pie_data <- data.frame(
    sample_type = c("faecal", "rectal_swab", "nasal_swab", "tissue"),
    proportion = c(data_row$faecal, data_row$rectal_swab, data_row$nasal_swab, 
                   data_row$tissue)
  )
  pie_data <- pie_data[pie_data$proportion > 0, ]  # Remove zeros
  pie <- ggplot(pie_data, aes(x = "", y = proportion, fill = sample_type)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = sample_colors) +
    theme_void() +
    theme(legend.position = "none")
  ggplotGrob(pie)
}

# Get map data for Australia
australia_map_data <- ne_countries(scale = "medium", country = "Australia", returnclass = "sf")
australia_states <- ne_states(country = "Australia", returnclass = "sf")

# Base plot
main_plot <- ggplot() +
  geom_sf(data = australia_map_data, fill = "grey95", color = "black") +
  geom_sf(data = australia_states, fill = NA, color = "black", linetype = "dotted") +
  # Add labels for each site
  geom_text_repel(data = pie_data,
                  aes(x = long, y = lat, label = site),
                  size = 3, box.padding = 0.5, max.overlaps = Inf) +
  coord_sf(xlim = c(110, 155), ylim = c(-45, -10), expand = FALSE) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_blank()
  )

# Add pie charts to the base plot
for (i in 1:nrow(pie_data)) {
  pie_grob <- create_pie_grob(pie_data[i, ])
  radius <- sqrt(pie_data$total_samples[i]) / 8  # Scale pies
  main_plot <- main_plot +
    annotation_custom(
      pie_grob,
      xmin = pie_data$long[i] - radius, xmax = pie_data$long[i] + radius,
      ymin = pie_data$lat[i] - radius, ymax = pie_data$lat[i] + radius
    )
}

# Add color legend
main_plot <- main_plot +
  geom_point(data = dummy_data, aes(x = x, y = y, fill = sample_type), shape = 21, size = 5) +
  scale_fill_manual(values = sample_colors, name = "Sample Type")

# Add size legend
main_plot <- main_plot +
  geom_point(data = size_legend_data, aes(x = x, y = y, size = sqrt(total_samples) / 8), fill = "grey50", shape = 21) +
  scale_size_continuous(
    name = "Number of Samples",
    breaks = sqrt(size_legend_data$total_samples) / 8, 
    labels = size_legend_data$total_samples
  )

print(main_plot)

ggsave(plot = main_plot, filename = "map.svg", height = 12, width = 12, limitsize = FALSE, dpi = 1000)




