##############################################################
# Name: pcoa_impact_visit1.R
# Date Created: 28 August 2018
# Author: Ryan Johnson
# Purpose: Plot script for humichip ordination.
##############################################################

source("code/ordination_plots/ordination_function.R")
source("code/ordination_plots/ordination_plot.R")

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