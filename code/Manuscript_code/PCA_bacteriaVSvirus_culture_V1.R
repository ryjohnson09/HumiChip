##################################################################
# Name: PCA_bacteriaVSvirus_culture_V1.R
# Author: Ryan Johnson
# Date Created: 20 Aug 2019
# Purpose: Functional probe PCA to determine functional
#  differences between bacterial vs viral infections
##################################################################

library(tidyverse)
library(vegan)

probe_type <- "Functional"
visits <- c(1)
detection_method <- "culture"


## Read in data -----------------------------------------------------------------------------------------
humichip <- suppressWarnings(suppressMessages(read_tsv("data/processed/Merged_humichip_Renormalized.tsv")))
treat <- suppressWarnings(suppressMessages(read_csv("data/processed/TrEAT_Clinical_Metadata_tidy.csv")))
ID_Decoder <- suppressWarnings(suppressMessages(read_csv("data/processed/ID_Decoder_Humichip.csv")))

# Extract visit 1 samples (don't need to be matched with 4 or 5)
ID_Decoder_select <- ID_Decoder %>%
  filter(visit_number %in% visits)


# Merge treat to ID_Decoder
treat <- ID_Decoder_select %>% 
  left_join(., treat, by = c("study_id" = "STUDY_ID"))


## Filter treat to include treat samples of interest -----------------------------------------------

treat_filter <- treat %>% 
  # Remove LOP and PLA
  filter(!Treatment %in% c("LOP", "PLA"))



## Select patient ID's that contain pathogens of interest -----------------------------------------
# Bacteria vs Viral
bacteria_list <- c("ETEC|EAEC|EPEC|EHEC|EIEC|C_diff|STEC|Shigella|Salmonella|Vibrio|EcoliO157H7|Aeromonas|Yersinia|Campylobacter")


treat_filter_bacVSvirus <- treat_filter %>% 
  select(glomics_ID, study_id, ends_with(paste0("_", detection_method))) %>% 
  mutate_at(vars(ends_with(paste0("_", detection_method))), list(~ifelse(. == "yes", 1, 0))) %>% 
  # Remove parasites
  select_at(vars(-matches("Giardia|E_histo|Crypto|Cyclospora"))) %>% 
  # Identify bacteria versus virus patients
  mutate(bacteria = rowSums(select_at(., vars(matches(bacteria_list))))) %>% 
  mutate(virus = rowSums(select_at(., vars(contains("virus"))))) %>% 
  # Select bactria only and virus only
  filter(bacteria > 0 & virus == 0 | bacteria == 0 & virus > 0) %>% 
  mutate(bact_or_virus = ifelse(bacteria > 0, "Bacteria", "Virus")) %>% 
  select(glomics_ID, study_id, bact_or_virus)






## Filter the humichip data to include only samples in treat_filter --------------------------------
humichip_filtered <- humichip %>% 
  select_at(c("Genbank.ID", "Gene", "Organism", "Lineage",
              "Gene_category", "Subcategory1", "Subcategory2", treat_filter_bacVSvirus$glomics_ID))




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
  full_join(., treat_filter_bacVSvirus, by = "glomics_ID")



## Plot ---------------------------------------------------------------------------------
# Factor Columns


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
                 color = bact_or_virus),
             alpha = 1, size = point_size, fill = "black") +
  stat_ellipse(aes(x = PC1, y = PC2, color = bact_or_virus), 
               level = 0.75, size = 1.5, linetype = "dashed")

humi_PCA_plot

ggsave(plot = humi_PCA_plot, filename = "results/figures/Manuscript_Figures/PCA_bacteriaVSvirus_culture_V1.png",
       height = 7,
       width = 8)


## Perform Statistics ----------------------------------------------

# Remove samples for which the groupings are NA values
new_ord_meta <- humi_PCA_metadata %>%
    filter(!is.na(bact_or_virus))

new_matrix <- humi_matrix[ ,colnames(humi_matrix) %in% new_ord_meta$glomics_ID]


# Display the groupings
stats_groups <- paste("Groupings:", paste(unique(new_ord_meta$bact_or_virus), collapse = ", "))

adonis_results <-
  # Ensure that new_ord_meta$glomics_ID is in the
  #  same order as colnames(new_matrix
  if (!all(new_ord_meta$glomics_ID == colnames(new_matrix))){
    stopApp("Error calculating adonis P-value")
  } else {
    
    # Calculate adonis results
    
    adonis_temp <-
      vegan::adonis(t(new_matrix) ~ new_ord_meta$bact_or_virus,
                    method = "bray",
                    perm = 99)
    
    paste("adonis p-value: ", adonis_temp$aov.tab$`Pr(>F)`[1])
  }



mrpp_results <-
  # Ensure that new_ord_meta$glomics_ID is in the
  #  same order as colnames(new_matrix())
  if (!all(new_ord_meta$glomics_ID == colnames(new_matrix))){
    stopApp("Error calculating mrpp P-value")
  } else {

    # Calculate mrpp results
      mrpp_temp <-  vegan::mrpp(t(new_matrix),
                                new_ord_meta$bact_or_virus)
    
    paste("mrpp p-value: ", mrpp_temp$Pvalue)
  }


## Print stats --------------------------------------------------
stats_groups

adonis_results

mrpp_results
