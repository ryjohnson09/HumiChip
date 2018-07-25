####################################################################
# Name: Humichip_NumSpecies_Visit.R
# Author: Ryan Johnson
# Date Created: 25 July 2018
# Purpose: Plot the number of detected species by visit number
####################################################################

library(tidyverse)

# Read in merged Humichip Data
Humichip <- read_tsv("data/processed/Merged_Humichip_Tidy.tsv", 
                     col_types = "cciccccccccd",
                     progress = TRUE)


# Extract the strain/species probes
Humichip_str_spe <- Humichip %>%
  filter(gene =="STR_SPE") %>%
  
  # Remove any rows with no signal
  filter(!is.na(Signal))


# Set threshold to determine if strain/species is present
Humichip_strain_signal <- Humichip_str_spe %>%
  group_by(glomics_ID, study_id, visit_number, species) %>%
  summarise(n = n(), mean_signal = mean(Signal)) %>%
  ungroup() %>%
  filter(n >= 8) # must have at least # probes present

# Calculate number of species present by studyID and visit number
num_species <- Humichip_strain_signal %>%
  ungroup() %>%
  group_by(study_id, visit_number) %>%
  summarize(n = n())


# Plot
num_species_plot <- ggplot(data = num_species, aes(x = factor(visit_number), y = n)) +
  geom_jitter(width = 0.2, height = 0, size = 2) +
  geom_boxplot(aes(fill = factor(visit_number)), alpha = 0.7, outlier.shape = NA) +
  xlab("Visit Number") +
  ylab("Number of Species Detected") +
  ggtitle("Number of Species Detected by Visit") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "none"
  )

ggsave(plot = num_species_plot, filename = "results/figures/NumSpecies_Visit.png")