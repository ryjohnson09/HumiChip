################################################
# Name: ordination_plot_functions.R
# Date Created: 27 August 2018
# Author: Ryan Johnson
# Purpose: This script will contain functions
#   and commands to generate various ordination
#   plots for the humichip data.
################################################

## Ordination Function
humichip_ord <- function(humichip_path,
                         ID_decoder_path,
                         treat_DB_path,
                         keep_unmatched_patients = TRUE,
                         remove_LOP_PLA = FALSE,
                         visit_numbers = c(1,4,5),
                         treatment_groups = "All",
                         pathogen_groups = "All",
                         probe_type = "All",
                         gene_category = "All",
                         ord_type = "PCoA"
                         ){
  
  ## Libraries ---------------------------------------------------------------
  library(tidyverse)
  library(vegan)
  library(ape)
  library(ggExtra)
 
  
   ## Dependencies Read In --------------------------------------------------
  ID_decoder <- suppressMessages(suppressWarnings(read_csv(ID_decoder_path)))
  treat <- suppressMessages(suppressWarnings(read_csv(treat_DB_path)))
  humichip <- suppressMessages(suppressWarnings(read_tsv(humichip_path)))
  
  
  ## Modify Input Data -------------------------------------------------------
  
  # Filter ID_decoder to only include samples in humichip
  ID_decoder <- ID_decoder %>%
    filter(glomics_ID %in% colnames(humichip))
  
  
  # Add glomics IDs to the treat DB
  treat <- treat %>%
    left_join(., ID_decoder, by = c("STUDY_ID" = "study_id"))
  
  
  # Only keep patients with matched visits 
  if(keep_unmatched_patients == FALSE){
    # Get vector of patients that have matched visits
    matched_samples <- ID_decoder %>%
      filter(visit_number %in% visit_numbers) %>%
      count(study_id) %>%
      filter(n == length(visit_numbers)) %>%
      pull(study_id)
    
    # Get glomics ID's that correspond to the patients
    # with matched samples
    matched_glomics <- ID_decoder %>%
      filter(study_id %in% matched_samples) %>%
      pull(glomics_ID)
    
    # Subset the humichip data
    humichip <- humichip %>%
      select_if(colnames(.) %in% c("Genbank.ID", "gene", "species", "lineage",
                                   "annotation", "geneCategory", "subcategory1",
                                   "subcategory2", matched_glomics))
    
    rm(matched_glomics, matched_samples)
      
    } else if (keep_unmatched_patients == TRUE){
      humichip <- humichip
    } else {
      stop("Error when subsetting matched/unmatched patients")
    }
  
  
  # Remove samples from humichip that are in the LOP & PLA tx groups
  if (isTRUE(remove_LOP_PLA)) {
    LOP_PLA_samples <- treat %>%
      filter(Treatment %in% c("LOP", "PLA")) %>%
      pull(glomics_ID)
    
    humichip <- humichip %>%
      select(-one_of(LOP_PLA_samples))
    
    rm(LOP_PLA_samples)
  }
  
  
  # Visit Number Filter
  visit_samples <- treat %>%
    filter(visit_number %in% visit_numbers) %>%
    pull(glomics_ID)
  
  humichip <- humichip %>%
    select(Genbank.ID, gene, species, lineage, annotation, 
           geneCategory, subcategory1, subcategory2, one_of(visit_samples))
  
  rm(visit_samples)
  
  
  # Treatment Group Filter
  if (all(treatment_groups %in% c("RIF", "AZI", "LEV"))){
    treatment_samples <- treat %>%
      filter(Treatment %in% treatment_groups) %>%
      pull(glomics_ID)
    
    humichip <- humichip %>%
      select(Genbank.ID, gene, species, lineage, annotation, 
             geneCategory, subcategory1, subcategory2, one_of(treatment_samples))
    
    rm(treatment_samples)
    
  } else if (treatment_groups == "All"){
    humichip <- humichip
  } else {
    stop("Invalid treatment group(s)")
  }
  
  
  # Pathogen Filter
  if (pathogen_groups == "All"){
    humichip <- humichip
  } else if (pathogen_groups %in% colnames(treat)){
    pathogen <- sym(pathogen_groups)
    pathogen_samples <- treat %>%
      filter(!!pathogen == "yes") %>%
      pull(glomics_ID)
    
    humichip <- humichip %>%
      select(Genbank.ID, gene, species, lineage, annotation, 
             geneCategory, subcategory1, subcategory2, one_of(pathogen_samples))
    
    rm(pathogen_samples, pathogen)
  } else {
    stop("Invalid pathogen group")
  }
  
  
  # Probe Filter
  if(probe_type == "Functional" & gene_category == "All"){
    humichip <- humichip %>%
      filter(gene != "STR_SPE")
  } else if (probe_type == "Functional" & probe_type != "All"){
    humichip <- humichip %>%
      filter(gene != "STR_SPE") %>%
      filter(geneCategory == input$geneCategory)
  } else if (probe_type == "Strain/Species"){
    humichip <- humichip %>%
      filter(gene == "STR_SPE")
  } else if (probe_type == "All"){
    humichip <- humichip
  } else {
    stop("Invalid probe type")
  }
  
  
  ## Ordination analysis ----------------------------------------------------------
  
  # Remove descriptive columns
  humichip <- humichip %>%
    select(-Genbank.ID, -gene, -species, -lineage, -annotation, 
           -geneCategory, -subcategory1, -subcategory2)
  
  # Set NA's to 0 and values not NA to original value
  humichip <- humichip %>%
    mutate_all(funs(ifelse(is.na(.), 0, .)))
  
  # Remove rows that equal 0
  humichip <- humichip[rowSums(humichip) != 0,]
  
  # Convert to matrix
  humichip <- as.matrix(humichip)
  
  
  # Ordination Analsys
  if(ord_type == "PCA"){
    ord_results <- vegan::rda(humichip)
    
    # Perform DCA analysis using vegan
  } else if (ord_type == "DCA") {
    ord_results <- vegan::decorana(humichip)
    
    # Perform PCoA using ape package
  } else if (ord_type == "PCoA"){
    humichip_dist <- vegan::vegdist(as.matrix(t(humichip)))
    ord_results <- ape::pcoa(humichip_dist, correction = "none")
    rm(humichip_dist)
  } else {
    stop("Invalid Ordination Type")
  }
  
  
  ## Merge ordination results into treat
  # PCA and DCA coordinates
  if(ord_type == "PCA" | ord_type == "DCA"){
    
    # Get the coordinated for the samples
    ord_coords <- scores(ord_results, display = "species")
    
    # Make tibble
    ord_coords <- as.data.frame(ord_coords) %>%
      rownames_to_column(var = "glomics_ID")
    
    ord_coords
    
    # PCoA Coords
  } else if (ord_type == "PCoA") {
    
    ord_coords <- ord_results$vectors[,1:2]
    
    # Make tibble
    ord_coords <- as.data.frame(ord_coords) %>%
      rownames_to_column(var = "glomics_ID")
  } else {
    stop("Error when extracting ordination coordinates")
  }
  
  
  # Extract proportion explained by first couple PC's for PCA and PCoA
  if (ord_type == "PCA"){
    prop_exp <- summary(eigenvals(ord_results))[2,] * 100
  } else if (ord_type == "PCoA"){
    prop_exp <- ord_results$values$Relative_eig * 100
  } else {
    prop_exp <- NA
  }
  
  
  # Merge into treat
    
  treat <- ord_coords %>%
    left_join(., treat, by = c("glomics_ID")) %>%
    
    # Factor columns
    mutate(visit_number = factor(visit_number)) %>%
    mutate(Impact_of_illness_on_activity_level = factor(Impact_of_illness_on_activity_level))

    
  # Return object
  
  humichip_return <- list(humichip_data = treat)
  humichip_return$prop_exp <- prop_exp
  humichip_return
}




#########################################
ord_plot <- function(ord_object,
                     ord_type, 
                     color_points_by = visit_number,
                     plot_title = "",
                     y_lab = "", 
                     x_lab = ""){
  
  if (ord_type == "PCoA"){
    x_axis = "Axis.1"
    y_axis = "Axis.2"
  } else if (ord_type == "PCA"){
    x_axis = "PC1"
    y_axis = "PC2"
  } else if (ord_type == "DCA"){
    x_axis = "DCA1"
    y_axis = "DCA2"
  } else {
    stop("Invalid x/y axes")
  }
  
  ord_plot <- ggplot(ord_object$humichip_data, 
                      aes_string(x = x_axis, y = y_axis)) +
    xlab(x_lab) +
    ylab(y_lab) +
    geom_point(aes_string(color = color_points_by), alpha = 0.8, size = 6) +
    geom_point(pch = 21, size = 6, alpha = 0.6) +
    ggtitle(plot_title) + 
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 14), 
      axis.text.x = element_text(size = 12, hjust = 1),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 14), 
      plot.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 12),
      legend.title = element_blank())
  ggMarginal(ord_plot, groupFill = TRUE, groupColour = TRUE)
}



###### Generate Figures #########################

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#

pcoa_LLSseverity <- humichip_ord(humichip_path = "data/processed/Merged_humichip.tsv",
                          ID_decoder_path = "data/processed/ID_Decoder.csv",
                          treat_DB_path = "data/processed/TrEAT_Clinical_Metadata_tidy.csv",
                          remove_LOP_PLA = TRUE,
                          visit_numbers = 1,
                          treatment_groups = "All",
                          pathogen_groups = "All",
                          probe_type = "All",
                          gene_category = "All",
                          ord_type = "PCoA")

png("results/figures/ordination_plots/pcoa_LLS_severity_visit1.png", height = 8, width = 9, units = "in", res = 300)
ord_plot(ord_object = pcoa_LLSseverity, 
         ord_type = "PCoA", 
         color_points_by = "LLS_severity", 
         plot_title = "PCoA Analysis", 
         x_lab = paste0("PCoA1 (", round(pcoa_LLSseverity$prop_exp[1], 2), "%)"),
         y_lab = paste0("PCoA2 (", round(pcoa_LLSseverity$prop_exp[2], 2), "%)"))
dev.off()

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#

pcoa_impact <- humichip_ord(humichip_path = "data/processed/Merged_humichip.tsv",
                                 ID_decoder_path = "data/processed/ID_Decoder.csv",
                                 treat_DB_path = "data/processed/TrEAT_Clinical_Metadata_tidy.csv",
                                 remove_LOP_PLA = TRUE,
                                 visit_numbers = 1,
                                 treatment_groups = "All",
                                 pathogen_groups = "All",
                                 probe_type = "All",
                                 gene_category = "All",
                                 ord_type = "PCoA")

png("results/figures/ordination_plots/pcoa_impact_visit1.png", height = 8, width = 9, units = "in", res = 300)
ord_plot(ord_object = pcoa_impact, 
         ord_type = "PCoA", 
         color_points_by = "Impact_of_illness_on_activity_level", 
         plot_title = "PCoA Analysis", 
         x_lab = paste0("PCoA1 (", round(pcoa_impact$prop_exp[1], 2), "%)"),
         y_lab = paste0("PCoA2 (", round(pcoa_impact$prop_exp[2], 2), "%)"))
dev.off()

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#

RIF <- humichip_ord(humichip_path = "data/processed/Merged_humichip.tsv",
                    ID_decoder_path = "data/processed/ID_Decoder.csv",
                    treat_DB_path = "data/processed/TrEAT_Clinical_Metadata_tidy.csv",
                    keep_unmatched_patients = FALSE,
                    remove_LOP_PLA = TRUE,
                    visit_numbers = c(1,5),
                    treatment_groups = "RIF",
                    pathogen_groups = "All",
                    probe_type = "All",
                    gene_category = "All",
                    ord_type = "PCoA")

png("results/figures/ordination_plots/pcoa_rif.png", height = 8, width = 9, units = "in", res = 300)
ord_plot(ord_object = RIF, 
         ord_type = "PCoA", 
         color_points_by = "visit_number", 
         plot_title = "PCoA Analysis", 
         x_lab = paste0("PCoA1 (", round(RIF$prop_exp[1], 2), "%)"),
         y_lab = paste0("PCoA2 (", round(RIF$prop_exp[2], 2), "%)"))
dev.off()
