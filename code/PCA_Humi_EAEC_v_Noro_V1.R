##################################################################
# Name: 
# Author: Ryan Johnson
# Date Created: 
##################################################################

library(tidyverse)
library(vegan)

probe_type <- "Functional"
visits <- c(1)
detection_method <- "culture"
pathogens <- c("EAEC", "Norovirus")


## Read in data -----------------------------------------------------------------------------------------
humichip <- suppressWarnings(suppressMessages(read_tsv("data/processed/Merged_humichip_Renormalized.tsv")))
treat <- suppressWarnings(suppressMessages(read_csv("data/processed/TrEAT_Clinical_Metadata_tidy.csv")))
ID_Decoder <- suppressWarnings(suppressMessages(read_csv("data/processed/ID_Decoder_Humichip.csv")))

# Only include matched isolate for visit 1 and visit 4 or 5
ID_Decoder_select <- ID_Decoder %>%
  filter(visit_number %in% visits) %>%
  group_by(study_id) %>% 
  filter(n() == length(visits)) 


# Merge treat to ID_Decoder
treat <- ID_Decoder_select %>% 
  left_join(., treat, by = c("study_id" = "STUDY_ID"))


## Filter treat to include treat samples of interest -----------------------------------------------

treat_filter <- treat %>% 
  # Remove LOP and PLA
  filter(!Treatment %in% c("LOP", "PLA")) %>% 
  ungroup()


## Only extract those that had JUST EAEC and JUST Norovirus (no coinfections)
treat_eaec_noro <- treat_filter %>% 
  select(study_id, glomics_ID, ends_with(paste0("_", detection_method))) %>% 
  mutate_at(vars(ends_with(paste0("_", detection_method))), list(~ifelse(. == "yes", 1, 0))) %>% 
  mutate(orig_rowSum = rowSums(.[-1:-2], na.rm = TRUE)) %>% 
  select_at(c("study_id", "glomics_ID", "orig_rowSum", paste(pathogens, detection_method, sep = "_"))) %>% 
  mutate(new_rowSum = rowSums(.[-1:-3])) %>% 
  filter(new_rowSum == orig_rowSum) %>% 
  filter(new_rowSum != 0) %>%
  filter(new_rowSum != length(pathogens)) %>% 
  select(-new_rowSum, -orig_rowSum) %>% 
  map2_df(., names(.), ~  replace(.x, .x == 1, gsub(x = .y, 
                                                    pattern = paste0("_", detection_method), 
                                                    replacement = ""))) %>% 
  mutate_at(vars(ends_with(detection_method)), list(~ifelse(. == 0, "", .))) %>% 
  unite(col = "detection_results", sep = "", paste(pathogens, detection_method, sep = "_"))



## Filter the humichip data to include only samples in treat_filter --------------------------------
humichip_filtered <- humichip %>% 
  select_at(c("Genbank.ID", "Gene", "Organism", "Lineage",
              "Gene_category", "Subcategory1", "Subcategory2", treat_eaec_noro$glomics_ID))


## Probe Filtering ---------------------------------------------------------
if(probe_type == "Functional"){
  humi_probes_filtered <- humichip_filtered %>%
    filter(Gene != "STR_SPE")
} else if (probe_type == "Strain/Species"){
  humi_probes_filtered <- humichip_filtered %>%
    filter(Gene == "STR_SPE")
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
  full_join(., treat_eaec_noro, by = "glomics_ID")



## Plot ---------------------------------------------------------------------------------

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
    legend.position = "bottom",
    legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size=7)))

# PCA Plot
humi_PCA_plot <- humi_PCA_base +
  xlab(paste0("PC1(", round(ord_prop_expl[[1]], 2), "%)")) +
  ylab(paste0("PC2(", round(ord_prop_expl[[2]], 2), "%)")) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 color = detection_results),
             alpha = 1, size = point_size, fill = "black") +
  stat_ellipse(aes(x = PC1, y = PC2, color = detection_results), 
               level = 0.75, size = 1.5, linetype = "dashed")

humi_PCA_plot

# ggsave(plot = humi_PCA_plot, filename = "results/figures/Manuscript_Figures/Humi_V1_impact_PCA.png", 
#        height = 7, width = 8)
