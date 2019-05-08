# EAEC_vs_Virus_ord.R
# Author: Ryan Johnson
# Date: 7 May 2019
# Purpose: Use ordination methods to detect changes in EAEC vs 
#  viral TD

library(tidyverse)
library(vegan)
library(plotly)

probe_type <- "Functional"
visits <- c(1)

## Read in data -----------------------------------------------------------------------------------------
humichip <- suppressWarnings(suppressMessages(read_tsv("data/processed/Merged_humichip.tsv")))
treat <- suppressWarnings(suppressMessages(read_csv("data/processed/TrEAT_Clinical_Metadata_tidy.csv")))
ID_Decoder <- suppressWarnings(suppressMessages(read_csv("data/processed/ID_Decoder_Humichip.csv")))


## Select for matched isolates --------------------------------------------------
# Only include matched isolates from `visits`
ID_Decoder_matched <- ID_Decoder %>%
  filter(visit_number %in% visits) %>%
  group_by(study_id) %>% 
  filter(n() == length(visits)) 

rm(ID_Decoder)

# Merge treat to ID_Decoder
treat <- ID_Decoder_matched %>% 
  left_join(., treat, by = c("study_id" = "STUDY_ID"))

## Filter treat to include samples of interest -----------------------------------------------
# Remove LOP and PLA
treat_filter <- treat %>% 
  filter(!Treatment %in% c("LOP", "PLA"))

# Extract EAEC positive samples
EAEC_patients_either <- treat_filter %>% 
  filter(EAEC_both == "yes") %>% 
  pull(glomics_ID)

# Extract Norovirus positive samples
NoroV_patients_either <- treat_filter %>% 
  filter(Norovirus_either == "yes") %>% 
  pull(glomics_ID)

# Final filter of treat
treat_filter <- treat_filter %>% 
  filter(glomics_ID %in% c(EAEC_patients_either, NoroV_patients_either))


rm(treat)

## Filter the humichip data to include only samples in treat_filter --------------------------------
humichip_filtered <- humichip %>% 
  select_at(c("Genbank_UniqueID", "gene", "species", "lineage",
              "annotation", "geneCategory", "subcategory1",
              "subcategory2", treat_filter$glomics_ID))

rm(humichip)

## Probe Filtering ---------------------------------------------------------
if(probe_type == "Functional"){
  humi_probes_filtered <- humichip_filtered %>%
    filter(gene != "STR_SPE")
} else if (probe_type == "Strain/Species"){
  humi_probes_filtered <- humichip_filtered %>%
    filter(gene == "STR_SPE")
} else if (probe_type == "All"){
  humi_probes_filtered <- humichip_filtered
} else {
  stop("Error filter probes")
}

## Prepare Data for Ordination Analysis ----------------------------------------
humi_matrix <- humi_probes_filtered %>% 
  
  # Set NA's to 0 and values not NA to original value
  select(starts_with("X")) %>%
  mutate_all(list(~ifelse(is.na(.), 0, .)))

# Remove rows that equal 0
humi_matrix <- humi_matrix[rowSums(humi_matrix) != 0,]

# Return matrix
humi_matrix <- as.matrix(humi_matrix)


## Ordination Analysis ----------------------------------------------------------
humi_PCA <- vegan::rda(t(humi_matrix))

# Extract coordinates as tibble
humi_coordinates <- as.data.frame(scores(humi_PCA, display = "sites")) %>%
  rownames_to_column(var = "glomics_ID")



# Get Proportion explained
ord_prop_expl <- summary(eigenvals(humi_PCA))[2,] * 100



## Merge Ordination Analysis with Metadata --------------------------------------
humi_PCA_metadata <- humi_coordinates %>%
  # Add study ID's
  full_join(., treat_filter, by = "glomics_ID")

## Add in Custom Category column -----------------------------------------------
# Viral or EAEC
humi_PCA_metadata <- humi_PCA_metadata %>% 
  mutate(viral_bacteria= ifelse(glomics_ID %in% EAEC_patients_either & 
                                glomics_ID %in% NoroV_patients_either, "Both", 
                         ifelse(glomics_ID %in% EAEC_patients_either, "EAEC", 
                         ifelse(glomics_ID %in% NoroV_patients_either, "NoroV", "Unknown"))))


## Plot ---------------------------------------------------------------------------------
# Factor Columns
humi_PCA_metadata <- humi_PCA_metadata %>% 
  mutate(visit_number = factor(visit_number)) %>%
  mutate(Impact_of_illness_on_activity_level = factor(Impact_of_illness_on_activity_level))


# Aesthetic sizes
axis_title_size <- 18
axis_text_size <- 16
title_size <- 20
legend_text_size <- 13
point_size <- 4


# Set up base plot
humi_PCA_base <- 
  ggplot(humi_PCA_metadata) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = axis_title_size),
    axis.text.x = element_text(size = axis_text_size, hjust = 1),
    axis.text.y = element_text(size = axis_text_size),
    axis.title.y = element_text(size = axis_title_size),
    plot.title = element_text(size = title_size, face = "bold"),
    legend.text = element_text(size = legend_text_size),
    legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size=7)))

# PCA Plot
humi_PCA_plot <- humi_PCA_base +
  xlab(paste0("PC1(", round(ord_prop_expl[[1]], 2), "%)")) +
  ylab(paste0("PC2(", round(ord_prop_expl[[2]], 2), "%)")) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 color = LLS_severity,
                 text = study_id),
                 alpha = 1, size = point_size) +
  ggtitle("PCA Analysis - Humichip")

ggplotly(humi_PCA_plot)
