##############################################################
# Name: pcoa_rif_visit15.R
# Date Created: 28 August 2018
# Author: Ryan Johnson
# Purpose: Plot script for humichip ordination.
##############################################################

source("code/ordination_plots/ordination_function.R")
source("code/ordination_plots/ordination_plot.R")


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

png("results/figures/ordination_plots/pcoa_rif_visit15.png", height = 8, width = 9, units = "in", res = 300)
ord_plot(ord_object = RIF, 
         ord_type = "PCoA", 
         color_points_by = "visit_number", 
         plot_title = "PCoA Analysis", 
         x_lab = paste0("PCoA1 (", round(RIF$prop_exp[1], 2), "%)"),
         y_lab = paste0("PCoA2 (", round(RIF$prop_exp[2], 2), "%)"))
dev.off()