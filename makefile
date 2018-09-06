
# List all targes by typing make list
.PHONY: list
list:
	@$(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$' | xargs -n 1




# Create decoder that links Glomics ID's to TrEAT Subject IDs
# Depends on:	data/raw/IDCRP_Glomics_Subject_ID_List_11-21-17.xlsx
#		code/ID_Decoder.R
# Produces:	data/processed/ID_Decoder.csv
data/processed/ID_Decoder.csv : data/raw/IDCRP_Glomics_Subject_ID_List_11-21-17.xlsx\
				code/ID_Decoder.R
	R -e "source('code/ID_Decoder.R')"




# Merge all of the Humichip datasets into one
# Depends on:	data/raw/HumiChip_New/HumiChip-1-LTO.txt
#		data/raw/HumiChip_New/HumiChip-2-LTO.txt
#		data/raw/HumiChip_New/HumiChip-3-LTO.txt
#		data/raw/HumiChip_New/HumiChip-4-LTO.txt
#		data/raw/HumiChip_New/HumiChip-5-LTO.txt
#		code/Merge_Humichip.R
# Produces:	data/processed/Merged_humichip.tsv
data/processed/Merged_humichip.tsv : data/raw/HumiChip_New/HumiChip-1-LTO.txt\
		                     data/raw/HumiChip_New/HumiChip-2-LTO.txt\
		                     data/raw/HumiChip_New/HumiChip-3-LTO.txt\
		                     data/raw/HumiChip_New/HumiChip-4-LTO.txt\
		                     data/raw/HumiChip_New/HumiChip-5-LTO.txt\
				     code/Merge_Humichip.R
	R -e "source('code/Merge_Humichip.R', echo=T)"



# Make the Merged Humichip data long/tidy
# Depends on:	data/processed/Merged_Humichip.tsv
#		data/processed/ID_Decoder.csv
#		code/Merged_Humichip_Tidy.R
# Produces:	data/processed/Merged_Humichip_Tidy.tsv
data/processed/Merged_Humichip_Tidy.tsv : data/processed/Merged_humichip.tsv\
			                  data/processed/ID_Decoder.csv\
			                  code/Merged_Humichip_Tidy.R
	R -e "source('code/Merged_Humichip_Tidy.R', echo=T)"



# Create Clinical Metadata Table Extracted from TrEAT DB
# Depends on:   data/processed/Merged_humichip.tsv
#               data/raw/TrEAT_Merge_2018.06.27.XLSX
#               data/raw/TrEAT_Merge_DataDictionary_2018.06.27.XLSX
#               data/raw/IDCRP_Glomics_Subject_ID_List_11-21-17.xlsx
#               code/Create_Clin_Metadata.R
# Produces:     data/processed/TrEAT_Clinical_Metadata_tidy.csv
data/processed/TrEAT_Clinical_Metadata_tidy.csv : data/processed/Merged_humichip.tsv\
                                                  data/raw/TrEAT_Merge_2018.06.27.XLSX\
                                                  data/raw/TrEAT_Merge_DataDictionary_2018.06.27.XLSX\
                                                  data/raw/IDCRP_Glomics_Subject_ID_List_11-21-17.xlsx\
                                                  code/Create_Clin_Metadata.R
	R -e "source('code/Create_Clin_Metadata.R', echo=T)"






# Visualize the Humichip Signal Distribution
# Depends on:	data/processed/Merged_Humichip_Tidy.tsv
#		code/Humichip_Signal_Distribution.R
# Produces:	results/figures/Humichip_Signal_Distribution.png
results/figures/Humichip_Signal_Distribution.png : data/processed/Merged_Humichip_Tidy.tsv\
				                   code/Humichip_Signal_Distribution.R
	R -e "source('code/Humichip_Signal_Distribution.R', echo=T)"




###########################
### Category Abundances ###
###########################


# Plot the Relative Abundance of all Gene Categories for all samples
# Depends on:	data/processed/Merged_Humichip_Tidy.tsv
#		code/Gene_Category_Abundance.R
# Produces:	results/figures/Humichipl_Gene_Abundance.png
results/figures/Humichip_Gene_Abundance.png : data/processed/Merged_Humichip_Tidy.tsv\
					      code/Gene_Category_Abundance.R
	R -e "source('code/Gene_Category_Abundance.R', echo=T)"



# Plot the Relative Abundance of Gene Categories by visit number
# Depends on:	data/processed/Merged_Humichip_Tidy.tsv
#               code/Gene_Category_Visit_Number.R
# Produces:	results/figures/Humichip_Gene_Abundance_Visit_Number.png
results/figures/Humichip_Gene_Abundance_Visit_Number.png : data/processed/Merged_Humichip_Tidy.tsv\
			       			            code/Gene_Category_Visit_Number.R
	R -e "source('code/Gene_Category_Visit_Number.R', echo=T)"






###################################
#### Humichip Response Ratios #####
###################################



# Humichip Response Ratios for Gene Categories based on Visit Number
# Depends on:	data/processed/Merged_Humichip.tsv
#		data/processed/ID_Decoder.csv
#               data/processed/TrEAT_Clinical_Metadata_tidy.csv
# 	     	code/Humichip_ResRatio_Visit.R
# Produces:	results/figures/Humichip_RespRatio_Visit.png
results/figures/Humichip_RespRatio_Visit.png : data/processed/Merged_Humichip.tsv\
					       data/processed/ID_Decoder.csv\
					       data/processed/TrEAT_Clinical_Metadata_tidy.csv\
					       code/Humichip_ResRatio_Visit.R
	R -e "source('code/Humichip_ResRatio_Visit.R', echo=T)"


# Humichip Response Ratios for Genes based on Visit Number
# Depends on:   data/processed/Merged_Humichip.tsv
#               data/processed/ID_Decoder.csv
#               data/processed/TrEAT_Clinical_Metadata_tidy.csv
#               code/Humichip_ResRatio_Visit.R
# Produces:     results/figures/Humichip_RespRatio_Visit_gene.png
results/figures/Humichip_RespRatio_Visit_gene.png : data/processed/Merged_Humichip.tsv\
                                             	    data/processed/ID_Decoder.csv\
                                              	    data/processed/TrEAT_Clinical_Metadata_tidy.csv\
                                              	     code/Humichip_ResRatio_Visit_gene.R
	R -e "source('code/Humichip_ResRatio_Visit_gene.R', echo=T)"



# Humichip Response Ratios for subcategory1 based on Visit Number
# Depends on:   data/processed/Merged_Humichip.tsv
#               data/processed/ID_Decoder.csv
#               data/processed/TrEAT_Clinical_Metadata_tidy.csv
#               code/Humichip_ResRatio_Visit_subcategory1.R
# Produces:     results/figures/Humichip_RespRatio_Visit_subcateogry1.png
results/figures/Humichip_RespRatio_Visit_subcategory1.png : data/processed/Merged_Humichip.tsv\
	                                                    data/processed/ID_Decoder.csv\
             	 	                                    data/processed/TrEAT_Clinical_Metadata_tidy.csv\
                                               		    code/Humichip_ResRatio_Visit_subcategory1.R
	R -e "source('code/Humichip_ResRatio_Visit_subcategory1.R', echo=T)"




# Humichip Response Ratios for subcategory2 based on Visit Number
# Depends on:   data/processed/Merged_Humichip.tsv
#               data/processed/ID_Decoder.csv
#               data/processed/TrEAT_Clinical_Metadata_tidy.csv
#               code/Humichip_ResRatio_Visit_subcategory2.R
# Produces:     results/figures/Humichip_RespRatio_Visit_subcateogry2.png
results/figures/Humichip_RespRatio_Visit_subcategory2.png : data/processed/Merged_Humichip.tsv\
                                                            data/processed/ID_Decoder.csv\
                                                            data/processed/TrEAT_Clinical_Metadata_tidy.csv\
                                                            code/Humichip_ResRatio_Visit_subcategory2.R
	R -e "source('code/Humichip_ResRatio_Visit_subcategory2.R', echo=T)"


# Humichip Response Ratios for species based on Visit Number
# Depends on:   data/processed/Merged_Humichip.tsv
#               data/processed/ID_Decoder.csv
#               data/processed/TrEAT_Clinical_Metadata_tidy.csv
#               code/Humichip_ResRatio_Visit_species.R
# Produces:     results/figures/Humichip_RespRatio_Visit_species.png
results/figures/Humichip_RespRatio_Visit_species.png : data/processed/Merged_Humichip.tsv\
                                                            data/processed/ID_Decoder.csv\
                                                            data/processed/TrEAT_Clinical_Metadata_tidy.csv\
                                                            code/Humichip_ResRatio_Visit_species.R
	R -e "source('code/Humichip_ResRatio_Visit_species.R', echo=T)"



###########################
#### Species Abundance ####
###########################

# Heatmap of species present in a specific number of samples
# Color indicates Normalized Signal Intensity
# Depends on:	data/processed/Merged_Humichip.tsv
#		code/Humichip_Species_Abundance.R
# Produces:	results/figures/Species_Abundance_All_Samples.png
results/figures/Species_Abundance_All_Samples.png : data/processed/Merged_Humichip.tsv\
						    code/Humichip_Species_Abundance.R
	R -e "source('code/Humichip_Species_Abundance.R', echo=T)"


# Scatter/Boxplot showing number of species per Visit
# Depends on:	data/processed/Merged_Humichip.tsv
#		code/Humichip_NumSpecies_Visit.R
# Produces:	results/figures/NumSpecies_Visit.png
results/figures/NumSpecies_Visit.png : data/processed/Merged_Humichip.tsv\
				       code/Humichip_NumSpecies_Visit.R
	R -e "source('code/Humichip_NumSpecies_Visit.R', echo=T)"





####################################
#### Ordination Plots/Functions ####
####################################

# Humichip ordination plots and functions. Will place all figures
# in the results/figures/ordination_plots directory
# Depends on:	data/processed/Merged_Humichip.tsv
#		data/processed/ID_Decoder.csv
#		data/processed/TrEAT_Clinical_Metadata_tidy.csv
#		code/ordination_plots/ordination_function.R
#		code/ordination_plots/ordination_plot.R
#		code/ordination_plots/*.R
# Produces:	results/figures/ordination_plots/*.png
results/figures/ordination_plots/pcoa_LLS_severity_visit1.png : data/processed/Merged_Humichip.tsv\
                				                data/processed/ID_Decoder.csv\
                               					data/processed/TrEAT_Clinical_Metadata_tidy.csv\
                               					code/ordination_plots/ordination_function.R\
                               					code/ordination_plots/ordination_plot.R\
        	                       				code/ordination_plots/pcoa_LLS_severity_visit1.R
	R -e "source('code/ordination_plots/pcoa_LLS_severity_visit1.R', echo=T)"


results/figures/ordination_plots/pcoa_impact_visit1.png : data/processed/Merged_Humichip.tsv\
                                                          data/processed/ID_Decoder.csv\
                                                          data/processed/TrEAT_Clinical_Metadata_tidy.csv\
                                                          code/ordination_plots/ordination_function.R\
                                                          code/ordination_plots/ordination_plot.R\
                                                          code/ordination_plots/pcoa_impact_visit1.R
	R -e "source('code/ordination_plots/pcoa_impact_visit1.R', echo=T)"


results/figures/ordination_plots/pcoa_rif_visit15.png : data/processed/Merged_Humichip.tsv\
                                                        data/processed/ID_Decoder.csv\
                                                        data/processed/TrEAT_Clinical_Metadata_tidy.csv\
                                                        code/ordination_plots/ordination_function.R\
                                                        code/ordination_plots/ordination_plot.R\
                                                        code/ordination_plots/pcoa_rif_visit15.R
	R -e "source('code/ordination_plots/pcoa_rif_visit15.R', echo=T)"


results/figures/ordination_plots/pcoa_lev_visit15.png : data/processed/Merged_Humichip.tsv\
                                                        data/processed/ID_Decoder.csv\
                                                        data/processed/TrEAT_Clinical_Metadata_tidy.csv\
                                                        code/ordination_plots/ordination_function.R\
                                                        code/ordination_plots/ordination_plot.R\
                                                        code/ordination_plots/pcoa_lev_visit15.R
	R -e "source('code/ordination_plots/pcoa_lev_visit15.R', echo=T)"


results/figures/ordination_plots/pcoa_azi_visit15.png : data/processed/Merged_Humichip.tsv\
                                                        data/processed/ID_Decoder.csv\
                                                        data/processed/TrEAT_Clinical_Metadata_tidy.csv\
                                                        code/ordination_plots/ordination_function.R\
                                                        code/ordination_plots/ordination_plot.R\
                                                        code/ordination_plots/pcoa_azi_visit15.R
	R -e "source('code/ordination_plots/pcoa_azi_visit15.R', echo=T)"
