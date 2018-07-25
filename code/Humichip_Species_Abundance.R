######################################################################
# Name: Humichip_Species_Abundance.R
# Author: Ryan Johnson
# Date Created: 25 July 2018
# Purpose: Create a heatmap showing the most prevalent species among
#          all samples in the HumiChip dataset
######################################################################

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
  filter(n >= 7) # must have at least 5 probes present


# Extract only strains/species shared by multiple samples
highly_prev_species <- Humichip_strain_signal %>%
  select(glomics_ID, species) %>%
  distinct() %>%
  group_by(species) %>%
  summarise(n = n()) %>%
  filter(n >= 10) %>% # threshold for how many samples must contain species
  pull(species)


# Subset the strain signal tibble to only contain
# highly abundant strains/species
Humichip_strain_signal_highly_prev <- Humichip_strain_signal %>%
  filter(species %in% highly_prev_species) %>%
  select(glomics_ID, species, mean_signal) %>%
  complete(glomics_ID, species) # needed so geom_tile will have blanks

# Read in decoder
ID_list <- read.csv("data/processed/ID_Decoder.csv")

# Join to Humichip data
Humichip_strain_signal_highly_prev <- Humichip_strain_signal_highly_prev %>%
  left_join(., ID_list, by = c("glomics_ID")) %>%
  select(glomics_ID, study_id, visit_number, everything())


# Plot as geom_tile
Humichip_species_plot <- 
  ggplot(data = Humichip_strain_signal_highly_prev, aes(x = glomics_ID, y = species)) +
  geom_tile(aes(fill = mean_signal)) + 
  scale_fill_viridis_c(name = "Normalized\nSignal\nIntensity", na.value = "black") +
  facet_grid(. ~ visit_number, scales = "free", space = "free") +
  xlab("Glomics ID") +
  ylab("Species/Strain") +
  ggtitle("Species Normalized Signal") +
  theme(
    axis.text.x = element_text(size = 12, angle = 90),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    strip.text.x = element_text(size = 18)
  )


Humichip_species_plot

ggsave(plot = Humichip_species_plot, filename = "results/figures/Species_Abundance_All_Samples.png", width = 25, height = 14)
