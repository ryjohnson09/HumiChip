
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
