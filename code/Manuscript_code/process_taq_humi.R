# Name: process_taq_humi.R
# Author: Ryan Johnson
# Date Created: 12 Aug 2019
# Purpose: Process the raw visit 1 and visit 5 taq data

# Load Libraries ---------------------------------------
library(tidyverse)
library(readxl)

# Read in data -------------------------------------------
taq_data_visit1 <- read_xlsx("data/raw/Taqman_results.xlsx", sheet = 1)
taq_data_visit5 <- read_xlsx("data/raw/Taqman_results.xlsx", sheet = 2)
ID_Decoder <- read_csv("data/processed/ID_Decoder_Humichip.csv")
metadata <- read_csv("data/processed/TrEAT_Clinical_Metadata_tidy.csv")


# Process raw data ----------------------------------------
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



# Merge the two visits -------------------------------------
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


# Add in Clin Metadata --------------------------------------
taq_data_full_treat <- taq_data_full %>%
  left_join(., metadata, by = c("STUDYID" = "STUDY_ID")) %>%
  # Not interested in samples not in a treatment groups #!#!#!#!#!
  filter(!is.na(Treatment))


## Write to csv ----------------------------------------------
write_csv(x = taq_data_full_treat, path = "data/processed/Taq_tidy_humi.csv")

