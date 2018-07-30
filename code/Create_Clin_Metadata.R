##################################################################
# Name: Create_Clin_Metadata.R
# Author: Ryan Johnson
# Date Created: 26 July, 2018
# Purpose: Extract data from the Merged TrEAT DB and compile into 
#          clean tidy format that can be appended to data files
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

# Read in the ID list
ID_list <- read_excel("data/raw/IDCRP_Glomics_Subject_ID_List_11-21-17.xlsx") %>%
  select_if(~ !all(is.na(.)))

glomics_ID <- unlist(select(ID_list, starts_with("X")), # stack all glomics ID values
                     use.names = FALSE)

study_ID <- unlist(select(ID_list, starts_with("study")), # stack all study ID values
                   use.names = FALSE)

ID_list <- tibble(glomics_ID, study_ID) %>% # merge into tibble, drop NA rows
  drop_na()

# Split study_ID to "study_ID" and "visit_numer
ID_list <- ID_list %>%
  separate(study_ID, into = c("study_id", "visit_number"), "(?<=\\d{4})-") %>%
  mutate(glomics_ID = paste0("X", glomics_ID)) # Add "X" in front of glomics_ID

rm(glomics_ID, study_ID) # clean up

# Filter for observations in ID_list that are found in the humichip_samples
ID_list_humichip <- ID_list %>% 
  filter(glomics_ID %in% colnames(humichip_samples))

# Extract unique patient ID's
humichip_patient_IDs <- unique(ID_list_humichip$study_id)

# Filter treat to include only those in humichip_patient_id
treat_humichip <- treat %>% 
  filter(STUDY_ID %in% humichip_patient_IDs)

# Clean
rm(humichip_patient_IDs)



## Create Clinical Metadata Table ------------------------------

# Extract relevant columns from treat DB
treat_humichip_clin <- treat_humichip %>% 
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
treat_humichip_clin <- treat_humichip_clin %>% 
  
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
           ifelse(Occult_blood_result == "N/A", NA, Occult_blood_result))




## Write to processed data ---------------------------------
write_csv(treat_humichip_clin, "data/processed/TrEAT_Clinical_Metadata_tidy.csv")

