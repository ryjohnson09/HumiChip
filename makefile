
# List all targes by typing make list
.PHONY: list
list:
	@$(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$' | xargs -n 1




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


# Create decoder that links Glomics ID's to TrEAT Subject IDs
# Depends on:   data/raw/IDCRP_Glomics_Subject_ID_List_11-21-17.xlsx
#               code/ID_Decoder_Humichip.R
#               data/processed/Merged_humichip.tsv
# Produces:     data/processed/ID_Decoder.csv
data/processed/ID_Decoder_Humichip.csv : data/raw/IDCRP_Glomics_Subject_ID_List_11-21-17.xlsx\
                                	 code/ID_Decoder_Humichip.R\
					 data/processed/Merged_humichip.tsv
	R -e "source('code/ID_Decoder_Humichip.R', echo=T)"
	
	
# Process TaqMan data for visit 1 and visit 5
# Depends on:	data/raw/Taqman_results.xlsx
#		data/processed/ID_Decoder_Humichip.csv
#		code/process_taq_humi.R
# Produces:	data/processed/Taq_tidy.csv
data/processed/Taq_tidy.csv : data/raw/Taqman_results.xlsx\
		              data/processed/ID_Decoder_Humichip.csv\
			      code/process_taq_humi.R
	R -e "source('code/process_taq_humi.R', echo=T)"



# Create Clinical Metadata Table Extracted from TrEAT DB
# Depends on:   data/processed/Merged_humichip.tsv
#               data/raw/TrEAT_Merge_ESBL_2018.09.13_v2.XLSX
#               data/raw/TrEAT_Merge_DataDictionary_2018.06.27.XLSX
#               data/raw/IDCRP_Glomics_Subject_ID_List_11-21-17.xlsx
#               code/Create_Clin_Metadata.R
#		data/processed/Taq_tidy.csv
# Produces:     data/processed/TrEAT_Clinical_Metadata_tidy.csv
data/processed/TrEAT_Clinical_Metadata_tidy.csv : data/processed/Merged_humichip.tsv\
                                                  data/raw/TrEAT_Merge_ESBL_2018.09.13_v2.XLSX\
                                                  data/raw/TrEAT_Merge_DataDictionary_2018.06.27.XLSX\
                                                  data/raw/IDCRP_Glomics_Subject_ID_List_11-21-17.xlsx\
						  data/processed/Taq_tidy.csv\
                                                  code/Create_Clin_Metadata.R
	R -e "source('code/Create_Clin_Metadata.R', echo=T)"


