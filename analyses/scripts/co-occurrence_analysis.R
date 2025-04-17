# Script to generate Figure 7 virus co-occurrence plot
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(readxl)

setwd("analyses/data/")
co_occurence_data <- read_xlsx("co_occurence_data.xlsx")

# Data cleaning and assembly
co_occurence_data$virus_species <- gsub("_", " ", co_occurence_data$virus_species)
sorted_levels <- sort(unique(co_occurence_data$virus_species))
co_occurence_data$virus_species <- factor(co_occurence_data$virus_species, levels = sorted_levels)
library_matrix <- table(co_occurence_data$library_summarised, co_occurence_data$virus_species)
total_occurrences <- colSums(library_matrix)

# Filter out viruses found in only one library
virus_filter <- total_occurrences >= 2
filtered_library_matrix <- library_matrix[, virus_filter]
filtered_sorted_levels <- names(virus_filter)[virus_filter]

# Recalculate the matrix for combined occurrences of virus pairs
library_occurrence_matrix <- t(filtered_library_matrix) %*% filtered_library_matrix

# Calculate percentages of co-occurrence relative to individual occurrences
percentage_matrix <- sweep(library_occurrence_matrix, 2, colSums(filtered_library_matrix), FUN="/") * 100
percentage_matrix <- t(percentage_matrix)  # Transpose for correct orientation for Figure 7 

combined_df <- expand.grid(Virus1 = filtered_sorted_levels, Virus2 = filtered_sorted_levels)
combined_df$Count <- as.vector(library_occurrence_matrix)
combined_df$Percentage <- as.vector(percentage_matrix)

my_palette <- colorRampPalette(brewer.pal(9, "Blues"))(100)

# Legend size breaks
size_breaks <- c(0, 8, 17)
size_actual <- c(0, 8, 17) 

# Figure 7
circle_plot <- ggplot(combined_df, aes(x = Virus2, y = Virus1, size = Count, fill = Percentage)) +
  geom_point(shape = 21, color = "black", stroke = 0.5) + 
  scale_size_continuous(
    breaks = size_breaks,
    range = c(1, 17),
    guide = guide_legend(title = "Co-occurrence Count", override.aes = list(size = size_actual), labels = size_breaks)) +  
  scale_fill_gradientn(
    colors = my_palette, 
    limits = c(0, 100), 
    oob = scales::squish,
    guide = guide_colourbar(ticks = FALSE, title.position = "top", title.hjust = 0.5)) +
  scale_x_discrete(limits = filtered_sorted_levels) + 
  scale_y_discrete(limits = rev(filtered_sorted_levels)) +
  labs(x = "Virus Species", y = "Virus Species", fill = "Co-occurrence %", size = "Co-occurrence Count") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.key = element_rect(color = "white")
  )

#ggsave(circle_plot, "circle_plot.svg", dpi = 1000, width = 8, height = 6)
