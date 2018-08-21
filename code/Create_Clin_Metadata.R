##################################################################
# Name: Create_Clin_Metadata.R
# Author: Ryan Johnson
# Date Created: 3 Aug, 2018
# Purpose: Extract data from the Merged TrEAT DB and compile into 
#          clean tidy format that can be appended to data files.
#          Tailored to the HumiChip data.
##################################################################

library(tidyverse)
library(readxl)

## Read in TrEAT and Clean ----------------------------------------------

# Read in TrEAT DB
treat <- read_excel("data/raw/TrEAT_Merge_2018.06.27.XLSX")

# Import data that explains treat column names
treat_explain <- read_excel("data/raw/TrEAT_Merge_DataDictionary_2018.06.27.XLSX")

# Remove spaces from truncated study ID
treat_explain <- treat_explain %>% 
  mutate(`Truncated Study ID` = gsub(" ", "_", `Truncated Study ID`)) %>% 
  mutate(`Truncated Study ID` = gsub("\\/", "_", `Truncated Study ID`)) %>%
  mutate(`Truncated Study ID` = gsub("\\(", "", `Truncated Study ID`)) %>%
  mutate(`Truncated Study ID` = gsub("\\)", "", `Truncated Study ID`))

# Ensure all colnames in treat match the values in treat_explain
foo <- tibble(STUDY_ID_TRUNC = colnames(treat[-1])) %>% 
  full_join(., treat_explain, by = "STUDY_ID_TRUNC")

# Check if labels match
ifelse(nrow(foo) == ncol(treat[,-1]), "Clear", stop("labels don't match"))

# Assign explanations to columns
colnames(treat) <- c("STUDY_ID", foo$`Truncated Study ID`)

# Clean
rm(treat_explain, foo)








## Import Merged humichip data -----------------------------

# Read in merged humichip data (just first instance, only need column names)
humichip_samples <- read_tsv("data/processed/Merged_humichip.tsv", n_max = 1)







## Import Glomics ID decoder list --------------------------------------

# Read in decoder
ID_list <- read.csv("data/processed/ID_Decoder.csv")

# Filter ID_list to only include samples from Glomics Data
ID_list <- ID_list %>%
  filter(glomics_ID %in% colnames(humichip_samples))

# Then, filter treat to only include samples in newly filtered ID_list
treat <- treat %>%
  filter(STUDY_ID %in% ID_list$study_id)

rm(ID_list, humichip_samples) # clean up





## Create Clinical Metadata Table ------------------------------

################################################
# As more clinical data becomes relevant, must 
# add it to the treat_clin
# tibble, then recode it below (if necessary)
################################################

# Extract relevant columns from treat DB
treat_clin <- treat %>% 
  select(STUDY_ID, 
         
         # Diarrhea and Fever
         Diarrhea_classification,
         Fever_present_at_presentation,
         
         # Stool
         Maximum_number_of_loose_liquid_stools_in_any_24_hours_prior_to_presentation, 
         Number_of_loose_liquid_stools_in_last_8_hours_prior_to_presentation, 
         Number_of_loose_liquid_stools_since_the_start_of_symptoms_prior_to_presentation,
         
         # Impact on activity
         Impact_of_illness_on_activity_level, 
         
         # Vomit
         Vomiting_present, 
         Number_of_vomiting_episodes,
         
         # Ab Cramps and Gas
         Abdominal_cramps_present_at_presentation, 
         Excesssive_gas_flatulence_present_at_presentation,
         
         # Nausea
         Nausea_present_at_presentation, 
         
         # Stool passing pain
         Ineffective_and_or_paiful_straining_to_pass_a_stool_at_presentation,
         
         # Tenesmus
         Tenesmus_present_at_presentation,
         
         # Malaise
         Malaise_present_at_presentation,
         
         # Incontenent / Constipation
         Fecal_incontinence_present_at_presentation,
         Constipation_present_at_presentation,
         
         # Other Symptoms
         Other_symptom_present_at_presentation,
         
         # Time from admittance to last unformed stool
         Time_to_last_unformed_stool,
         
         # Treatment
         Treatment,
         
         # Time to cure
         Time_to_cure,
         
         # Gross blood in stool
         Gross_blood_in_stool,
         
         # Occult blood test
         Occult_blood_result)


# Recode to make more legible
treat_clin <- treat_clin %>% 
  
  mutate(Diarrhea_classification = 
           ifelse(Diarrhea_classification == 1, "AWD", "Febrile")) %>%
  
  mutate(Fever_present_at_presentation = 
           ifelse(Fever_present_at_presentation == 0, "No", 
                  ifelse(Fever_present_at_presentation == 1, "Yes", NA))) %>%
  
  mutate(Vomiting_present =
           ifelse(Vomiting_present == 0, "No", 
                  ifelse(Vomiting_present == 1, "Yes", NA))) %>%
  
  mutate(Abdominal_cramps_present_at_presentation = 
           ifelse(Abdominal_cramps_present_at_presentation == 0, "No",
                  ifelse(Abdominal_cramps_present_at_presentation == 1, "Yes", NA))) %>%
  
  mutate(Excesssive_gas_flatulence_present_at_presentation = 
           ifelse(Excesssive_gas_flatulence_present_at_presentation == 0, "No",
                  ifelse(Excesssive_gas_flatulence_present_at_presentation == 1, "Yes", NA))) %>%
  
  mutate(Nausea_present_at_presentation =
           ifelse(Nausea_present_at_presentation == 0, "No", 
                  ifelse(Nausea_present_at_presentation == 1, "Yes", NA))) %>%
  
  mutate(Ineffective_and_or_paiful_straining_to_pass_a_stool_at_presentation = 
           ifelse(Ineffective_and_or_paiful_straining_to_pass_a_stool_at_presentation == 0, "No",
                  ifelse(Ineffective_and_or_paiful_straining_to_pass_a_stool_at_presentation == 1, "Yes", NA))) %>%
  
  mutate(Tenesmus_present_at_presentation =
           ifelse(Tenesmus_present_at_presentation == 0, "No",
                  ifelse(Tenesmus_present_at_presentation == 1, "Yes", NA))) %>%
  
  mutate(Malaise_present_at_presentation = 
           ifelse(Malaise_present_at_presentation == 0, "No", 
                  ifelse(Malaise_present_at_presentation == 1, "Yes", NA))) %>%
  
  mutate(Fecal_incontinence_present_at_presentation = 
           ifelse(Fecal_incontinence_present_at_presentation == 0, "No",
                  ifelse(Fecal_incontinence_present_at_presentation == 1, "Yes", NA))) %>%
  
  mutate(Constipation_present_at_presentation =
           ifelse(Constipation_present_at_presentation == 0, "No", 
                  ifelse(Constipation_present_at_presentation == 1, "Yes", NA))) %>%
  
  mutate(Other_symptom_present_at_presentation = 
           ifelse(Other_symptom_present_at_presentation == 0, "No", 
                  ifelse(Other_symptom_present_at_presentation == 1, "Yes", NA))) %>%
  
  mutate(Occult_blood_result =
           ifelse(Occult_blood_result == "N/A", NA, Occult_blood_result)) %>%
  
  
  #############################
  ### Add in custom columns ###
  #############################
  mutate(LLS_severity = 
    ifelse(Number_of_loose_liquid_stools_in_last_8_hours_prior_to_presentation <= 1, "mild",
    ifelse(Number_of_loose_liquid_stools_in_last_8_hours_prior_to_presentation %in% c(2,3,4), "moderate",
    ifelse(Number_of_loose_liquid_stools_in_last_8_hours_prior_to_presentation >= 5, "severe", NA))))
         

##########################
### Pathogen Detection ###
##########################
treat_path <- treat %>% 
  select(STUDY_ID,
        
         # C. Diff
         C._difficile_tcdA_result_from_taq_stool_sample,
         C._difficile_tcdA_result_from_taq_stool_card_sample,
         C._difficile_tcdB_result_from_taq_stool_sample,
         C._difficile_tcdB_result_from_taq_stool_card_sample,
         
         # EAEC
         `E._coli-like_isolate_1_is_EAEC_routine_microbiology`,
         `E._coli-like_isolate_2_is_EAEC_routine_microbiology`,
         `E._coli-like_isolate_3_is_EAEC_routine_microbiology`,
         `E._coli-like_isolate_4_is_EAEC_routine_microbiology`,
         `E._coli-like_isolate_5_is_EAEC_routine_microbiology`,
         EAEC_aaiC_result_from_taq_stool_sample,
         EAEC_aaiC_result_from_taq_stool_card_sample,
         EAEC_aatA_result_from_taq_stool_sample,
         EAEC_aatA_result_from_taq_stool_card_sample,
         EAEC_aggR_result_from_taq_stool_sample,
         EAEC_aggR_result_from_taq_stool_card_sample,
         
         # ETEC
         `E._coli-like_isolate_1_is_ETEC_routine_microbiology`,
         `E._coli-like_isolate_2_is_ETEC_routine_microbiology`,
         `E._coli-like_isolate_3_is_ETEC_routine_microbiology`,
         `E._coli-like_isolate_4_is_ETEC_routine_microbiology`,
         `E._coli-like_isolate_5_is_ETEC_routine_microbiology`,
         ETEC_LT_result_from_taq_stool_sample,
         ETEC_LT_result_from_taq_stool_card_sample,
         ETEC_STh_result_from_taq_stool_sample,
         ETEC_STh_result_from_taq_stool_card_sample,
         ETEC_STp_result_from_taq_stool_sample,
         ETEC_STp_result_from_taq_stool_card_sample,
         
         # EPEC
         `E._coli-like_isolate_1_is_EPEC_routine_microbiology`,
         `E._coli-like_isolate_2_is_EPEC_routine_microbiology`,
         `E._coli-like_isolate_3_is_EPEC_routine_microbiology`,
         `E._coli-like_isolate_4_is_EPEC_routine_microbiology`,
         `E._coli-like_isolate_5_is_EPEC_routine_microbiology`,
         EPEC_eae_result_from_taq_stool_sample,
         EPEC_eae_result_from_taq_stool_card_sample,
         
         # EHEC
         `E._coli-like_isolate_1_is_EHEC_routine_microbiology`,
         `E._coli-like_isolate_2_is_EHEC_routine_microbiology`,
         `E._coli-like_isolate_3_is_EHEC_routine_microbiology`,
         `E._coli-like_isolate_4_is_EHEC_routine_microbiology`,
         `E._coli-like_isolate_5_is_EHEC_routine_microbiology`,
         
         # EIEC
         `E._coli-like_isolate_1_is_EIEC_routine_microbiology`,
         `E._coli-like_isolate_2_is_EIEC_routine_microbiology`,
         `E._coli-like_isolate_3_is_EIEC_routine_microbiology`,
         `E._coli-like_isolate_4_is_EIEC_routine_microbiology`,
         `E._coli-like_isolate_5_is_EIEC_routine_microbiology`,
         
         # STEC
         STEC_stx1_result_from_taq_stool_sample,
         STEC_stx1_result_from_taq_stool_card_sample,
         STEC_stx2_result_from_taq_stool_sample,
         STEC_stx2_result_from_taq_stool_card_sample,
         
         # Shigella
         Shigella_routine_microbiology,
         S._dysenteriae_I_result_from_taq_stool_sample,
         S._dysenteriae_I_result_from_taq_stool_card_sample,
         S._flexneri_6_result_from_taq_stool_sample,
         S._flexneri_6_result_from_taq_stool_card_sample,
         `S._flexneri_non-6_result_from_taq_stool_sample`,
         `S._flexneri_non-6_result_from_taq_stool_card_sample`,
         S._sonnei_result_from_taq_stool_sample,
         S._sonnei_result_from_taq_stool_card_sample,
         
         # Other Bacteria
         Aremonas_result_from_taq_stool_sample,
         Aremonas_result_from_taq_stool_card_sample,
         Campylobacter_routine_microbiology,
         Campylobacter_result_from_taq_stool_sample,
         Campylobacter_result_from_taq_stool_card_sample,
         C._coli_result_from_taq_stool_sample,
         C._coli_result_from_taq_stool_card_sample,
         C._jejuni_result_from_taq_stool_sample,
         C._jejuni_result_from_taq_stool_card_sample,
         C._jejuni_coli_result_from_taq_stool_sample,
         C._jejuni_coli_result_from_taq_stool_card_sample,
         E._coli_O157_H7_result_from_taq_stool_sample,
         E._coli_O157_H7_result_from_taq_stool_card_sample,
         Salmonella_routine_microbiology,
         Salmonella_result_from_taq_stool_sample,
         Salmonella_result_from_taq_stool_card_sample,
         Vibrio_routine_microbiology,
         V._parahaemolyticus_result_from_taq_stool_sample,
         V._parahaemolyticus_result_from_taq_stool_card_sample,
         Yersinia_routine_microbiology,
         
         # Viral
         Adenovirus_result_from_taq_stool_sample,
         Adenovirus_result_from_taq_stool_card_sample,
         Astrovirus_result_from_taq_stool_sample,
         Astrovirus_result_from_taq_stool_card_sample,
         Rotavirus_routine_microbiology,
         Rotavirus_result_from_taq_stool_sample,
         Rotavirus_result_from_taq_stool_card_sample,
         Sapovirus_result_from_taq_stool_sample,
         Sapovirus_result_from_taq_stool_card_sample,
         Norovirus_GI_routine_microbiology,
         Norovirus_GII_routine_microbiology,
         Norovirus_other_routine_microbiology,
         Norovirus_GI_result_from_taq_stool_sample,
         Norovirus_GI_result_from_taq_stool_card_sample,
         Norovirus_GII_result_from_taq_stool_sample,
         Norovirus_GII_result_from_taq_stool_card_sample,
         
         # Parasite
         C._parvum_routine_microbiology,
         Cryptosporidium_result_from_taq_stool_sample,
         Cryptosporidium_result_from_taq_stool_card_sample,
         Cyclospora_result_from_taq_stool_sample,
         Cyclospora_result_from_taq_stool_card_sample,
         E._histolytica_routine_microbiology,
         E._histolytica_result_from_taq_stool_sample,
         E._histolytica_result_from_taq_stool_card_sample,
         G._lamblia_routine_microbiology,
         Giardia_result_from_taq_stool_sample,
         Giardia_result_from_taq_stool_card_sample
         )

# Convert all columns to numeric (will coerce 'undetermined' to NA)
treat_path <- treat_path %>%
  mutate_at(vars(-STUDY_ID, -ends_with("_routine_microbiology")), as.numeric) %>%
  # Convert all values under 30 to 1 and over 35 to NA
  mutate_at(vars(-STUDY_ID, -ends_with("_routine_microbiology")), funs(
    ifelse(. < 30, 1, 0)))



# Set Rules for Pathogen detection
treat_path <- treat_path %>%
  
  #EAEC
  mutate(EAEC_culture = ifelse(
    `E._coli-like_isolate_1_is_EAEC_routine_microbiology` %in% c("Yes", "yes", "1") |
    `E._coli-like_isolate_2_is_EAEC_routine_microbiology` %in% c("Yes", "yes", "1") |
    `E._coli-like_isolate_3_is_EAEC_routine_microbiology` %in% c("Yes", "yes", "1") |
    `E._coli-like_isolate_4_is_EAEC_routine_microbiology` %in% c("Yes", "yes", "1") |
    `E._coli-like_isolate_5_is_EAEC_routine_microbiology` %in% c("Yes", "yes", "1"), "yes", "no")) %>% 
    
  mutate(EAEC_taq = ifelse( 
    EAEC_aaiC_result_from_taq_stool_sample == 1 | 
    EAEC_aaiC_result_from_taq_stool_card_sample == 1 |
    EAEC_aatA_result_from_taq_stool_sample == 1 |
    EAEC_aatA_result_from_taq_stool_card_sample == 1 |
    EAEC_aggR_result_from_taq_stool_sample == 1 | 
    EAEC_aggR_result_from_taq_stool_card_sample == 1, "yes", "no")) %>%
  
  mutate(EAEC_either = ifelse(
    EAEC_culture == "yes" | EAEC_taq == "yes", "yes", "no")) %>%
  
  mutate(EAEC_both = ifelse(
    EAEC_culture == "yes" & EAEC_taq == "yes", "yes", "no")) %>%
  
  
  # ETEC
  mutate(ETEC_culture = ifelse(
    `E._coli-like_isolate_1_is_ETEC_routine_microbiology` %in% c("Yes", "yes", "1") |
    `E._coli-like_isolate_2_is_ETEC_routine_microbiology` %in% c("Yes", "yes", "1") |
    `E._coli-like_isolate_3_is_ETEC_routine_microbiology` %in% c("Yes", "yes", "1") |
    `E._coli-like_isolate_4_is_ETEC_routine_microbiology` %in% c("Yes", "yes", "1") |
    `E._coli-like_isolate_5_is_ETEC_routine_microbiology` %in% c("Yes", "yes", "1"), "yes", "no")) %>%
  
  mutate(ETEC_taq = ifelse(
    ETEC_LT_result_from_taq_stool_sample == 1 | 
    ETEC_LT_result_from_taq_stool_card_sample == 1 |
    ETEC_STh_result_from_taq_stool_sample == 1 | 
    ETEC_STh_result_from_taq_stool_card_sample == 1 |
    ETEC_STp_result_from_taq_stool_sample == 1 | 
    ETEC_STp_result_from_taq_stool_card_sample == 1, "yes", "no")) %>%
  
  mutate(ETEC_either = ifelse(
    ETEC_culture == "yes" | ETEC_taq == "yes", "yes", "no")) %>%
  
  mutate(ETEC_both = ifelse(
    ETEC_culture == "yes" & ETEC_taq == "yes", "yes", "no")) %>%
  
  
  # EPEC
  mutate(EPEC_culture = ifelse(
    `E._coli-like_isolate_1_is_EPEC_routine_microbiology` %in% c("Yes", "yes", "1") |
    `E._coli-like_isolate_2_is_EPEC_routine_microbiology` %in% c("Yes", "yes", "1") |
    `E._coli-like_isolate_3_is_EPEC_routine_microbiology` %in% c("Yes", "yes", "1") |
    `E._coli-like_isolate_4_is_EPEC_routine_microbiology` %in% c("Yes", "yes", "1") |
    `E._coli-like_isolate_5_is_EPEC_routine_microbiology` %in% c("Yes", "yes", "1"), "yes", "no")) %>%
  
  mutate(EPEC_taq = ifelse(
    EPEC_eae_result_from_taq_stool_sample == 1 |
    EPEC_eae_result_from_taq_stool_card_sample == 1, "yes", "no")) %>%
  
  mutate(EPEC_either = ifelse(
    EPEC_culture == "yes" | EPEC_taq == "yes", "yes", "no")) %>%
  
  mutate(EPEC_both = ifelse(
    EPEC_culture == "yes" & EPEC_taq == "yes", "yes", "no")) %>%
  
  
  # EHEC
  mutate(EHEC_culture = ifelse(
    `E._coli-like_isolate_1_is_EHEC_routine_microbiology` %in% c("Yes", "yes", "1") |
    `E._coli-like_isolate_2_is_EHEC_routine_microbiology` %in% c("Yes", "yes", "1") |
    `E._coli-like_isolate_3_is_EHEC_routine_microbiology` %in% c("Yes", "yes", "1") |
    `E._coli-like_isolate_4_is_EHEC_routine_microbiology` %in% c("Yes", "yes", "1") |
    `E._coli-like_isolate_5_is_EHEC_routine_microbiology` %in% c("Yes", "yes", "1"), "yes", "no")) %>%
  
  mutate(EHEC_either = ifelse(
    EHEC_culture == "yes", "yes", "no")) %>%
  
  mutate(EHEC_both = NA) %>%
  

  # EIEC
  mutate(EIEC_culture = ifelse(
    `E._coli-like_isolate_1_is_EIEC_routine_microbiology` %in% c("Yes", "yes", "1") |
    `E._coli-like_isolate_2_is_EIEC_routine_microbiology` %in% c("Yes", "yes", "1") |
    `E._coli-like_isolate_3_is_EIEC_routine_microbiology` %in% c("Yes", "yes", "1") |
    `E._coli-like_isolate_4_is_EIEC_routine_microbiology` %in% c("Yes", "yes", "1") |
    `E._coli-like_isolate_5_is_EIEC_routine_microbiology` %in% c("Yes", "yes", "1"), "yes", "no")) %>%
  
  mutate(EIEC_either = ifelse(
    EIEC_culture == "yes", "yes", "no")) %>%
  
  mutate(EIEC_both = NA) %>%
  
  # C.diff
  mutate(C_diff_taq = ifelse(
    C._difficile_tcdA_result_from_taq_stool_sample == 1 |
    C._difficile_tcdA_result_from_taq_stool_card_sample == 1 |
    C._difficile_tcdB_result_from_taq_stool_sample == 1 |
    C._difficile_tcdB_result_from_taq_stool_card_sample == 1, "yes", "no")) %>%
  
  mutate(C_diff_either = ifelse(
    C_diff_taq == "yes", "yes", "no")) %>%
  
  mutate(C_diff_both = NA) %>%
  
  # STEC
  mutate(STEC_taq = ifelse(
    STEC_stx1_result_from_taq_stool_sample == 1 | 
    STEC_stx1_result_from_taq_stool_card_sample == 1 | 
    STEC_stx2_result_from_taq_stool_sample == 1 | 
    STEC_stx2_result_from_taq_stool_card_sample == 1, "yes", "no")) %>%
  
  mutate(STEC_either = ifelse(
    STEC_taq == "yes", "yes", "no")) %>%
  
  mutate(STEC_both = NA) %>%
  
  # Shigella
  mutate(Shigella_culture = ifelse(
    Shigella_routine_microbiology %in% c("Yes", "yes", "1"), "yes", "no")) %>%
  
  mutate(Shigella_taq = ifelse(
    S._dysenteriae_I_result_from_taq_stool_sample == 1 | 
    S._dysenteriae_I_result_from_taq_stool_card_sample == 1 | 
    S._flexneri_6_result_from_taq_stool_sample == 1 | 
    S._flexneri_6_result_from_taq_stool_card_sample == 1 | 
    `S._flexneri_non-6_result_from_taq_stool_sample` == 1 | 
    `S._flexneri_non-6_result_from_taq_stool_card_sample` == 1 | 
    S._sonnei_result_from_taq_stool_sample == 1 | 
    S._sonnei_result_from_taq_stool_card_sample == 1, "yes", "no")) %>%
  
  mutate(Shigella_either = ifelse(
    Shigella_culture == "yes" | Shigella_taq == "yes", "yes", "no")) %>%
  
  mutate(Shigella_both = ifelse(
    Shigella_culture == "yes" & Shigella_taq == "yes", "yes", "no")) %>%
  
  # Aeromonas
  mutate(Aeromonas_taq = ifelse(
    Aremonas_result_from_taq_stool_sample == 1 |
    Aremonas_result_from_taq_stool_card_sample == 1, "yes", "no")) %>%
  
  mutate(Aeromonas_either = ifelse(
    Aeromonas_taq == "yes", "yes", "no")) %>%
  
  mutate(Aeromonas_both = NA) %>%
  
  # Campylobacter
  mutate(Campylobacter_culture = ifelse(
    Campylobacter_routine_microbiology %in% c("Yes", "yes", "1"), "yes", "no")) %>%
  
  mutate(Campylobacter_taq = ifelse(
    Campylobacter_result_from_taq_stool_sample == 1 |
    Campylobacter_result_from_taq_stool_card_sample == 1 |
    C._coli_result_from_taq_stool_sample == 1 |
    C._coli_result_from_taq_stool_card_sample == 1 |
    C._jejuni_result_from_taq_stool_sample == 1 |
    C._jejuni_result_from_taq_stool_card_sample == 1 |
    C._jejuni_coli_result_from_taq_stool_sample == 1 |
    C._jejuni_coli_result_from_taq_stool_card_sample == 1, "yes", "no")) %>%
  
  mutate(Campylobacter_either = ifelse(
    Campylobacter_culture == "yes" | Campylobacter_taq == "yes", "yes", "no")) %>%
  
  mutate(Shigella_both = ifelse(
    Campylobacter_culture == "yes" & Campylobacter_taq == "yes", "yes", "no")) %>%
  
  
  # EcoliO157H7
  mutate(EcoliO157H7_taq = ifelse(
    E._coli_O157_H7_result_from_taq_stool_sample == 1 |
    E._coli_O157_H7_result_from_taq_stool_card_sample == 1, "yes", "no")) %>%
  
  mutate(EcoliO157H7_either = ifelse(
    EcoliO157H7_taq == "yes", "yes", "no")) %>%
  
  mutate(EcoliO157H7_both = NA) %>%
  
 
   # Salmonella
  mutate(Salmonella_culture = ifelse(
    Salmonella_routine_microbiology %in% c("Yes", "yes", "1"), "yes", "no")) %>%
  
  mutate(Salmonella_taq = ifelse(
    Salmonella_result_from_taq_stool_sample == 1 |
    Salmonella_result_from_taq_stool_card_sample == 1, "yes", "no")) %>%
  
  mutate(Salmonella_either = ifelse(
    Salmonella_culture == "yes" | Salmonella_taq == "yes", "yes", "no")) %>%
  
  mutate(Shigella_both = ifelse(
    Salmonella_culture == "yes" & Salmonella_taq == "yes", "yes", "no")) %>%
  
  
  # Vibrio
  mutate(Vibrio_culture = ifelse(
    Vibrio_routine_microbiology  %in% c("Yes", "yes", "1"), "yes", "no")) %>%
  
  mutate(Vibrio_taq = ifelse(
    V._parahaemolyticus_result_from_taq_stool_sample == 1 |
    V._parahaemolyticus_result_from_taq_stool_card_sample == 1, "yes", "no")) %>%
  
  mutate(Vibrio_either = ifelse(
    Vibrio_culture == "yes" | Vibrio_taq == "yes", "yes", "no")) %>%
  
  mutate(Vibrio_both = ifelse(
    Vibrio_culture == "yes" & Vibrio_taq == "yes", "yes", "no")) %>%
  
  # Adenovirus
  mutate(Adenovirus_taq = ifelse(
    Adenovirus_result_from_taq_stool_sample == 1 |
    Adenovirus_result_from_taq_stool_card_sample == 1, "yes", "no")) %>%
  
  mutate(Adenovirus_either = ifelse(
    Adenovirus_taq == "yes", "yes", "no")) %>%
  
  mutate(Adenovirus_both = NA) %>%
  
  # Astrovirus
  mutate(Astrovirus_taq = ifelse(
    Astrovirus_result_from_taq_stool_sample == 1 |
    Astrovirus_result_from_taq_stool_card_sample == 1, "yes", "no")) %>%
  
  mutate(Astrovirus_either = ifelse(
    Astrovirus_taq == "yes", "yes", "no")) %>%
  
  mutate(Astrovirus_both = NA) %>%
  
  # Rotavirus
  mutate(Rotavirus_culture = ifelse(
    Rotavirus_routine_microbiology %in% c("Yes", "yes", "1"), "yes", "no")) %>%
  
  mutate(Rotavirus_taq = ifelse(
    Rotavirus_result_from_taq_stool_sample == 1 |
    Rotavirus_result_from_taq_stool_card_sample == 1, "yes", "no")) %>%
  
  mutate(Rotavirus_either = ifelse(
    Rotavirus_culture == "yes" | Rotavirus_taq == "yes", "yes", "no")) %>%
  
  mutate(Rotavirus_both = ifelse(
    Rotavirus_culture == "yes" & Rotavirus_taq == "yes", "yes", "no")) %>%
  
  
  # Sapovirus
  mutate(Sapovirus_taq = ifelse(
    Sapovirus_result_from_taq_stool_sample == 1 |
    Sapovirus_result_from_taq_stool_card_sample == 1, "yes", "no")) %>%
  
  mutate(Sapovirus_either = ifelse(
    Sapovirus_taq == "yes", "yes", "no")) %>%
  
  mutate(Sapovirus_both = NA) %>%
  
  # Norovirus
  mutate(Norovirus_culture = ifelse(
    Norovirus_GI_routine_microbiology %in% c("Yes", "yes", "1") |
    Norovirus_GII_routine_microbiology %in% c("Yes", "yes", "1") |
    Norovirus_other_routine_microbiology %in% c("Yes", "yes", "1"), "yes", "no")) %>%
  
  mutate(Norovirus_taq = ifelse(
    Norovirus_GI_result_from_taq_stool_sample == 1 |
    Norovirus_GI_result_from_taq_stool_card_sample == 1 |
    Norovirus_GII_result_from_taq_stool_sample == 1 |
    Norovirus_GII_result_from_taq_stool_card_sample == 1, "yes", "no")) %>%
  
  mutate(Norovirus_either = ifelse(
    Norovirus_culture == "yes" | Norovirus_taq == "yes", "yes", "no")) %>%
  
  mutate(Norovirus_both = ifelse(
    Norovirus_culture == "yes" & Norovirus_taq == "yes", "yes", "no")) %>%

  # Cryptosporidium
  mutate(Cryptosporidium_culture = ifelse(
    C._parvum_routine_microbiology  %in% c("Yes", "yes", "1"), "yes", "no")) %>%
  
  mutate(Cryptosporidium_taq = ifelse(
    Cryptosporidium_result_from_taq_stool_sample == 1 |
    Cryptosporidium_result_from_taq_stool_card_sample == 1, "yes", "no")) %>%
  
  mutate(Cryptosporidium_either = ifelse(
    Cryptosporidium_culture == "yes" | Cryptosporidium_taq == "yes", "yes", "no")) %>%
  
  mutate(Cryptosporidium_both = ifelse(
    Cryptosporidium_culture == "yes" & Cryptosporidium_taq == "yes", "yes", "no")) %>%
  
  
  # Cyclospora
  mutate(Cyclospora_taq = ifelse(
    Cyclospora_result_from_taq_stool_sample == 1 |
    Cyclospora_result_from_taq_stool_card_sample == 1, "yes", "no")) %>%
  
  mutate(Cyclospora_either = ifelse(
    Cyclospora_taq == "yes", "yes", "no")) %>%
  
  mutate(Cyclospora_both = NA) %>%
  
  # E_histolytica
  mutate(E_histolytica_culture = ifelse(
    E._histolytica_routine_microbiology %in% c("Yes", "yes", "1"), "yes", "no")) %>%
  
  mutate(E_histolytica_taq = ifelse(
    E._histolytica_result_from_taq_stool_sample == 1 |
    E._histolytica_result_from_taq_stool_card_sample == 1, "yes", "no")) %>%
  
  mutate(E_histolytica_either = ifelse(
    E_histolytica_culture == "yes" | E_histolytica_taq == "yes", "yes", "no")) %>%
  
  mutate(E_histolytica_both = ifelse(
    E_histolytica_culture == "yes" & E_histolytica_taq == "yes", "yes", "no")) %>%
  
  
  # Giardia
  mutate(Giardia_culture = ifelse(
    G._lamblia_routine_microbiology %in% c("Yes", "yes", "1"), "yes", "no")) %>%
  
  mutate(Giardia_taq = ifelse(
    Giardia_result_from_taq_stool_sample == 1 |
    Giardia_result_from_taq_stool_card_sample == 1, "yes", "no")) %>%
  
  mutate(Giardia_either = ifelse(
    Giardia_culture == "yes" | Giardia_taq == "yes", "yes", "no")) %>%
  
  mutate(Giardia_both = ifelse(
    Giardia_culture == "yes" & Giardia_taq == "yes", "yes", "no")) %>%
  
  # Yersinia
  mutate(Yersinia_culture = ifelse(
    Yersinia_routine_microbiology  %in% c("Yes", "yes", "1"), "yes", "no")) %>%
  
  mutate(Yersinia_either = ifelse(
    Yersinia_culture == "yes", "yes", "no")) %>%
  
  mutate(Yersinia_both = NA)








# Select columns of interest
treat_path <- treat_path %>%
  select(STUDY_ID, EAEC_culture, EAEC_taq, EAEC_either, EAEC_both, ETEC_culture, ETEC_taq, 
         ETEC_either, ETEC_both, EPEC_culture, EPEC_taq, EPEC_either, EPEC_both, 
         EHEC_culture, EHEC_either, EHEC_both, EIEC_culture, EIEC_either, EIEC_both, 
         C_diff_taq, C_diff_either, C_diff_both, STEC_taq, STEC_either, STEC_both, 
         Shigella_culture, Shigella_taq, Shigella_either, Shigella_both, Aeromonas_taq, 
         Aeromonas_either, Aeromonas_both, Campylobacter_culture, Campylobacter_taq, 
         Campylobacter_either, Shigella_both, EcoliO157H7_taq, EcoliO157H7_either, 
         EcoliO157H7_both, Salmonella_culture, Salmonella_taq, Salmonella_either, 
         Shigella_both, Vibrio_culture, Vibrio_taq, Vibrio_either, Vibrio_both, 
         Adenovirus_taq, Adenovirus_either, Adenovirus_both, Astrovirus_taq, 
         Astrovirus_either, Astrovirus_both, Rotavirus_culture, Rotavirus_taq, 
         Rotavirus_either, Rotavirus_both, Sapovirus_taq, Sapovirus_either, Sapovirus_both, 
         Norovirus_culture, Norovirus_taq, Norovirus_either, Norovirus_both, 
         Cryptosporidium_culture, Cryptosporidium_taq, Cryptosporidium_either, 
         Cryptosporidium_both, Cyclospora_taq, Cyclospora_either, Cyclospora_both, 
         E_histolytica_culture, E_histolytica_taq, E_histolytica_either, E_histolytica_both, 
         Giardia_culture, Giardia_taq, Giardia_either, Giardia_both, Yersinia_culture, 
         Yersinia_either, Yersinia_both)




## Mege Treat Clin and Treat Path -----------------------------------------------
treat_full <- treat_clin %>%
  full_join(., treat_path, by = c("STUDY_ID"))


## Write to processed data ---------------------------------
write_csv(treat_full, "data/processed/TrEAT_Clinical_Metadata_tidy.csv")


## Clean ----------------
rm(treat, treat_clin, treat_path)
