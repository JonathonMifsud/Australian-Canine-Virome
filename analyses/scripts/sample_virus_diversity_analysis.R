# Stats and plotting underlying Figure 2 and Supplmenetary Table 6 and 7
library(tidyverse)
library(readxl)
library(gridExtra)
library(cowplot)
library(FSA)
library(scales)

setwd("analyses/data/")

# Load the virus table and the sample table
virus_data <- readxl::read_xlsx("Supplementary_Table_5.xlsx")
sample_data_original <- readxl::read_xlsx("Supplementary_Table_1.xlsx", sheet = "Pool")

sample_data_cleaned <- sample_data_original %>%
  rename(library = `Sample Name`,
         location = Location,
         sample_type = `Sample Type`,
         host_sex = `Host Sex`,
         host_age = `Host Age`,
         diseases_and_conditions = `Summary of Indivdual Diseases and Conditions of Samples`,
         individuals_per_pool = `Individuals Per Pool`)

# Merge virus data with sample data
combined_data_merged <- sample_data_cleaned %>%
  left_join(virus_data, by = "library")

# clean the merged data
combined_data_cleaned <- combined_data_merged %>%
  select(library, location.x, sample_type.x, host_sex.x, host_age.x, individuals_per_pool.x, 
         virus_species, percentage_of_library_reads) %>%
  rename(location = location.x, sample_type = sample_type.x, host_sex = host_sex.x, host_age = host_age.x,
         individuals_per_pool = individuals_per_pool.x)

# replace missing values with 0 and calculate virus number
combined_data_final <- combined_data_cleaned %>%
  mutate(virus_number = ifelse(is.na(virus_species), 0, 1),
         virus_species = ifelse(is.na(virus_species), "None", virus_species),
         percentage_of_library_reads = ifelse(is.na(percentage_of_library_reads), 0, percentage_of_library_reads))

# filter out those with missing metadata and dingo libs
data_filtered_final <- combined_data_final %>%
  filter(!is.na(sample_type) & sample_type != "Unknown",
         !is.na(host_sex) & host_sex != "Unknown",
         !is.na(host_age) & host_age != "Unknown")
  
# Summarize the data at the library level to get virus count, diversity, and abundance per library
# Group by the library and summarize virus data without being grouped by virus species
data_filtered_abundance <- data_filtered_final %>%
  group_by(library, sample_type, host_sex, host_age) %>%
  summarise(percentage_of_library_reads = sum(percentage_of_library_reads),
            individuals_per_pool = first(individuals_per_pool), 
            .groups = 'drop')

data_filtered_count <-  data_filtered_final %>%
  group_by(library, sample_type, host_sex, host_age) %>%
  distinct(virus_species, .keep_all = T) %>%
  ungroup() %>% # as we are distinct by library and virus_species we need to increase the counter
  # for those libs in which there were multiple of the "same" virus species within a single lib
  # we can do that manually
  mutate(virus_number = case_when((library == "hopevale_4" & virus_species == "canine_astrovirus") ~ 3,
                                  (library == "perth_5" & virus_species == "canine_calicivirus") ~ 3,
                                  (library == "huntervalley_6" & virus_species == "canine_calicivirus") ~ 3,
                                  (library == "perth_2" & virus_species == "canine_picodicistrovirus") ~ 2,
                                  TRUE ~ as.numeric(.$virus_number))) %>%
  group_by(library, sample_type, host_sex, host_age) %>%
  summarise(virus_count = sum(virus_number),
            virus_diversity = n_distinct(virus_species[virus_species != "None"]), .groups = 'drop')

# Virus Number
kruskal_sample_type_count <- kruskal.test(virus_count ~ sample_type, data = data_filtered_count)
dunn_sample_type_count <- dunnTest(virus_count ~ sample_type, data = data_filtered_count)
kruskal_host_sex_count <- kruskal.test(virus_count ~ host_sex, data = data_filtered_count)
kruskal_host_age_count <- kruskal.test(virus_count ~ host_age, data = data_filtered_count)

# Virus Diversity
kruskal_sample_type_diversity <- kruskal.test(virus_diversity ~ sample_type, data = data_filtered_count)
dunn_sample_type_diversity <- dunnTest(virus_diversity ~ sample_type, data = data_filtered_count)
kruskal_host_sex_diversity <- kruskal.test(virus_diversity ~ host_sex, data = data_filtered_count)
kruskal_host_age_diversity <- kruskal.test(virus_diversity ~ host_age, data = data_filtered_count)

# Virus Abundance
# Calculate virus abundance based on percentage_of_library_reads
kruskal_sample_type_abundance <- kruskal.test(percentage_of_library_reads ~ sample_type, data = data_filtered_abundance)
dunn_sample_type_abundance <- dunnTest(percentage_of_library_reads ~ sample_type, data = data_filtered_abundance)
kruskal_host_sex_abundance <- kruskal.test(percentage_of_library_reads ~ host_sex, data = data_filtered_abundance)
kruskal_host_age_abundance <- kruskal.test(percentage_of_library_reads ~ host_age, data = data_filtered_abundance)

# Stat summary table Supplementary Table 6
kruskal_test_results <- data.frame(
  Test = c("Sample Type count", 
           "Host Sex count", 
           "Host Age count", 
           "Sample Type abundance", 
           "Host Sex abundance", 
           "Host Age abundance",
           "Sample Type diversity", 
           "Host Sex diversity", 
           "Host Age diversity"),
  
  Chi_Squared = c(kruskal_sample_type_count$statistic,
                  kruskal_host_sex_count$statistic,
                  kruskal_host_age_count$statistic,
                  kruskal_sample_type_abundance$statistic,
                  kruskal_host_sex_abundance$statistic,
                  kruskal_host_age_abundance$statistic,
                  kruskal_sample_type_diversity$statistic,
                  kruskal_host_sex_diversity$statistic,
                  kruskal_host_age_diversity$statistic),
  
  p_value = c(kruskal_sample_type_count$p.value,
              kruskal_host_sex_count$p.value,
              kruskal_host_age_count$p.value,
              kruskal_sample_type_abundance$p.value,
              kruskal_host_sex_abundance$p.value,
              kruskal_host_age_abundance$p.value,
              kruskal_sample_type_diversity$p.value,
              kruskal_host_sex_diversity$p.value,
              kruskal_host_age_diversity$p.value)
)

# Function to tidy dunnTest results
tidy_dunn <- function(dunn_result, variable, metric) {
  dunn_result$res %>%
    select(variable, metric, Comparison, Z, P.unadj, P.adj)
}

# Tidy each Dunn test result
dunn_sample_type_count_df <- tidy_dunn(dunn_sample_type_count, "sample_type", "virus_count")
dunn_sample_type_diversity_df <- tidy_dunn(dunn_sample_type_diversity, "sample_type", "virus_diversity")
dunn_sample_type_abundance_df <- tidy_dunn(dunn_sample_type_abundance, "sample_type", "virus_abundance")

# Combine all Dunn test results into Supplementary Table 7
dunn_summary_table <- bind_rows(
  dunn_sample_type_count_df,
  dunn_sample_type_diversity_df,
  dunn_sample_type_abundance_df
)

# Plotting
color_palette <- c("#666EB3", "#D7C770", "#A3BA72")
# Plot theme
themeGeneral1 <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(
    colour = "black",
    fill = NA,
    size = 1
  ),
  axis.ticks.length = unit(-1.5, "mm"),
  axis.ticks = element_line(size = 1),
  axis.line.x = element_line(),
  axis.line.y = element_line(),
  axis.title.x = element_text(
    face = "bold",
    size = 13,
    margin = margin(10, 0, 0, 0),
    family="Helvetica"
  ),
  strip.text.x = element_text(size = 12),
  axis.title.y = element_text(
    face = "bold",
    size = 13,
    margin = margin(0, 10, 0, 0),
    family="Helvetica"
  ),
  axis.text.x = element_text(
    face = "bold",
    size = 11,
    colour = "black",
    margin = margin(15, 0, 0, 0, "pt"),
    family="Helvetica"
  ),
  axis.text.y = element_text(
    face = "bold",
    size = 11,
    colour = "black",
    margin = margin(0, 15, 0, 0, "pt"),
    family="Helvetica"
  ),
  panel.spacing.y = unit(1, "points")
)

# give a very small minimun value to those with percentages so small that they have been rounded by excel for visualisation
# this is so that the log transformation can work
data_filtered_abundance_for_plot <- data_filtered_abundance %>%
  mutate(percentage_of_library_reads = ifelse(percentage_of_library_reads == 0, 1e-7, percentage_of_library_reads))

# Virus number by Sample Type
boxplot_number_sample_type <- ggplot(data_filtered_count, aes(x = sample_type, y = virus_count, fill = sample_type)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitter(width = 0.1, height = 0), alpha = 0.5, size = 1) +
  labs(y = "Virus number") +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  themeGeneral1 +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

# Virus number by Host Sex
boxplot_number_host_sex <- ggplot(data_filtered_count, aes(x = host_sex, y = virus_count, fill = host_sex)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitter(width = 0.1, height = 0), alpha = 0.5, size = 1) +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  themeGeneral1 +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())

# Virus number by Host Age
boxplot_number_host_age <- ggplot(data_filtered_count, aes(x = host_age, y = virus_count, fill = host_age)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitter(width = 0.1, height = 0), alpha = 0.5, size = 1) +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  themeGeneral1 +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())

# Virus Abundance by Sample Type 
boxplot_abundance_sample_type <- ggplot(data_filtered_abundance_for_plot, aes(x = sample_type, y = percentage_of_library_reads, fill = sample_type)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitter(width = 0.1, height = 0), alpha = 0.5, size = 1) +
  scale_y_log10(labels = label_log()) +
  annotation_logticks(sides = "l") +
  labs(y = "Virus Abundance (Log Percentage of Library Reads)") +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  themeGeneral1 +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

# Virus Abundance by Host Sex
boxplot_abundance_host_sex <- ggplot(data_filtered_abundance_for_plot, aes(x = host_sex, y = percentage_of_library_reads, fill = host_sex)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitter(width = 0.1, height = 0), alpha = 0.5, size = 1) +
  scale_y_log10(labels = label_log()) +
  annotation_logticks(sides = "l") + 
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  themeGeneral1 +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())

# Virus Abundance by Host Age
boxplot_abundance_host_age <- ggplot(data_filtered_abundance_for_plot, aes(x = host_age, y = percentage_of_library_reads, fill = host_age)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitter(width = 0.1, height = 0), alpha = 0.5, size = 1) +
  scale_y_log10(labels = label_log()) +
  annotation_logticks(sides = "l") +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  themeGeneral1 +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())

# Virus Diversity by Sample Type
boxplot_diversity_sample_type <- ggplot(data_filtered_count, aes(x = sample_type, y = virus_diversity, fill = sample_type)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitter(width = 0.1, height = 0), alpha = 0.5, size = 1) +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  themeGeneral1

# Virus Diversity by Host Sex
boxplot_diversity_host_sex <- ggplot(data_filtered_count, aes(x = host_sex, y = virus_diversity, fill = host_sex)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitter(width = 0.1, height = 0), alpha = 0.5, size = 1) +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  themeGeneral1 +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank())

# Virus Diversity by Host Age
boxplot_diversity_host_age <- ggplot(data_filtered_count, aes(x = host_age, y = virus_diversity, fill = host_age)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitter(width = 0.1, height = 0), alpha = 0.5, size = 1) +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  themeGeneral1 +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank())


final_plot <- plot_grid(
  boxplot_number_sample_type + theme(legend.position = "none"), 
  boxplot_number_host_sex + theme(legend.position = "none"), 
  boxplot_number_host_age + theme(legend.position = "none"),
  boxplot_abundance_sample_type + theme(legend.position = "none"), 
  boxplot_abundance_host_sex + theme(legend.position = "none"), 
  boxplot_abundance_host_age + theme(legend.position = "none"),
  boxplot_diversity_sample_type + theme(legend.position = "none"),
  boxplot_diversity_host_sex + theme(legend.position = "none"),
  boxplot_diversity_host_age + theme(legend.position = "none"),
  ncol = 3, nrow = 3
)

ggsave(plot = final_plot, filename = "boxplot.svg", height = 12, width = 12, limitsize = FALSE, dpi = 1200)

