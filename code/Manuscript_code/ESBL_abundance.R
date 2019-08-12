######################################################################
# Name: ESBL_abundance
# Author: Ryan Johnson
# Date Created: 12 Aug 2019
# Purpose: Determine the ESBL abundance at visit 1 to visit 5
#   using TAC.
######################################################################

library(tidyverse)
library(readxl)
library(gt)

## Load data --------------------------------------------------------------------------------------------
taq_data_visit1 <- read_xlsx("data/raw/Taqman_results.xlsx", sheet = 1)
taq_data_visit5 <- read_xlsx("data/raw/Taqman_results.xlsx", sheet = 2)
ID_Decoder <- read_csv("data/processed/ID_Decoder_Humichip.csv")
metadata <- read_csv("data/processed/TrEAT_Clinical_Metadata_tidy.csv")


# Set Ct threshold ---------------------------------------------------------------------------------
taq_threshold <- 35


## Process Taq data ---------------------------------------------------------
taq_data_visit1_tidy <- taq_data_visit1 %>%
  # Remove if SITE == "N/A
  filter(SITE != "N/A") %>%
  # Remove non-needed columns
  select(-OBS, -SITE, -DAYS_BETWEEN:-FTACARD_STOOLSAMPLE_PAIRED, -SAMPLE_COLLECTED, -SAMPLE_RUN) %>%
  # Remove Blanks
  filter(SPECIMENS_1 != "Blank") %>%
  # Make Tidy
  gather(key = Target, value = taq_value, -SPECIMENS_1, -STUDYID) %>%
  # Add Visit Column
  mutate(Visit = "Visit_1")

taq_data_visit5_tidy <- taq_data_visit5 %>%
  # Remove if SITE == "N/A
  filter(SITE != "N/A") %>%
  # Remove non-needed columns
  select(-OBS, -SITE, -STOOL_ALONE) %>%
  # Remove Blanks
  filter(SPECIMENS_1 != "Blank") %>%
  # Strip ending from study_id
  mutate(STUDYID = str_replace(STUDYID, "-5$", "")) %>%
  # Make Tidy
  gather(key = Target, value = taq_value, -SPECIMENS_1, -STUDYID) %>%
  # Add Visit Column
  mutate(Visit = "Visit_5")



# Merge the two visits 
taq_data_full <- taq_data_visit1_tidy %>%
  full_join(., taq_data_visit5_tidy) %>%
  # Fix taq values so indeterminates are NA, Undetermined is 0, and any
  #  other values are as is
  mutate(taq_value = as.numeric(
    ifelse(taq_value == "Undetermined", 0, 
           ifelse(taq_value %in% c("Indeterminate", "Indeterminat", "."), NA, taq_value)))) %>%
  select(STUDYID, SPECIMENS_1, Target, Visit, taq_value) %>% 
  # Filter for HuMiChip study IDs
  filter(STUDYID %in% ID_Decoder$study_id)


## Add in Clin Metadata --------------------------------------------------------------
taq_data_full_treat <- taq_data_full %>%
  left_join(., metadata, by = c("STUDYID" = "STUDY_ID")) %>%
  # Not interested in samples not in a treatment groups #!#!#!#!#!
  filter(!is.na(Treatment))


## Prep data -------------------------------------------
taq_data_select <- taq_data_full_treat %>%
  select(STUDYID, SPECIMENS_1, Target, Visit, taq_value, Treatment, ESBL_V1, ESBL_V5) %>% 
  
  # Filter just for ESBL targets (and 16S if needed)
  filter(Target %in% c("CMY", "CTX", "KPC", "NDM", "SHV", "TEM")) %>%
  
  # Filter for just the three tx groups of interest
  filter(Treatment %in% c("RIF", "LEV", "AZI")) %>% 
  
  # Don't have card data for visit 5
  filter(SPECIMENS_1 == "Stool") %>%
  
  # Remove invalide (NA) taq values
  filter(!is.na(taq_value)) %>% 
  
  # Remove non-matched isolates
  group_by(STUDYID, Target) %>%
  filter(n() == 2) %>%
  ungroup()


# Set rules for detection
taq_ESBL_stool_matched <- taq_data_select %>%
  mutate(detection = ifelse(taq_value == 0, 0, 
                            ifelse(taq_value > 0 & taq_value < taq_threshold, 1, 0)))


## Calculate Percentage of ESBL positives ----------------------------------
ESBL_percentages <- taq_ESBL_stool_matched %>%
  group_by(Treatment, Target, Visit) %>%
  summarise(total_ESBLs = sum(detection), total_in_tx = n()) %>% 
  mutate(percent_ESBL = round((total_ESBLs / total_in_tx) * 100, 2)) %>% 
  ungroup()


# Create table
ESBL_percentages %>% 
  filter(!Target %in% c("NDM", "KPC")) %>% 
  # Rename columns
  mutate(Treatment = ifelse(Treatment == "AZI", "Azithromycin",
                            ifelse(Treatment == "LEV", "Levofloxacin", "Rifaximin"))) %>% 
  mutate(Visit = ifelse(Visit == "Visit_1", "Visit 1", "Visit 5")) %>% 
  select(-total_ESBLs, -total_in_tx) %>% 
  spread(Visit, percent_ESBL) %>% 
  group_by(Treatment) %>% 
  gt() %>% 
  # Add spanner
  tab_spanner(
    label = "Percent Detection",
    columns = vars(
      "Visit 1", "Visit 5"
    )
  ) %>% 
  # Remove target column header
  cols_label(
    Target = ""
  )
