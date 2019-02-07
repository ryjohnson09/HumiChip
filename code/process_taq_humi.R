# Name: process_taq_humi.R
# Author: Ryan Johnson
# Date Created: 6 February 2019
# Purpose: Process the raw visit 1 and visit 5 taq data

# Load Libraries ---------------------------------------
library(tidyverse)
library(readxl)

# Read in data -------------------------------------------
taq_data_visit1 <- read_xlsx("data/raw/Taqman_results.xlsx", sheet = 1)
taq_data_visit5 <- read_xlsx("data/raw/Taqman_results.xlsx", sheet = 2)
ID_Decoder <- read_csv("data/processed/ID_Decoder_Humichip.csv")


# Process raw data ----------------------------------------
taq_data_visit1_format <- taq_data_visit1 %>%
  # Remove if SITE == "N/A
  filter(SITE != "N/A") %>%
  # Remove non-needed columns
  select(STUDYID, SPECIMENS_1, CMY, CTX, KPC, NDM, SHV, TEM) %>%
  # Remove Blanks
  filter(SPECIMENS_1 != "Blank") %>%
  # Make Tidy
  gather(key = Target, value = taq_value, -SPECIMENS_1, -STUDYID) %>%
  # Add Visit Column
  mutate(visit_number = 1) %>%
  # Merge the Specimen and target
  unite(col = Target, c("Target", "SPECIMENS_1")) %>%
  # filter for samples in humichip data
  filter(STUDYID %in% ID_Decoder$study_id)


taq_data_visit5_format <- taq_data_visit5 %>%
  # Remove if SITE == "N/A
  filter(SITE != "N/A") %>%
  # Remove non-needed columns
  select(STUDYID, SPECIMENS_1, CMY, CTX, KPC, NDM, SHV, TEM) %>%
  # Remove Blanks
  filter(SPECIMENS_1 != "Blank") %>%
  # Strip ending from study_id
  mutate(STUDYID = str_replace(STUDYID, "-5$", "")) %>%
  # Make Tidy
  gather(key = Target, value = taq_value, -SPECIMENS_1, -STUDYID) %>%
  # Add Visit Column
  mutate(visit_number = 5) %>% 
  # Merge the Specimen and target
  unite(col = Target, c("Target", "SPECIMENS_1")) %>%
  # filter for samples in humichip data
  filter(STUDYID %in% ID_Decoder$study_id)
  


# Merge the two visits, then merge intoID Decoder ------------------------------
taq_data_full <- taq_data_visit1_format %>%
  full_join(., taq_data_visit5_format) %>% 
  # Fix taq values so indeterminates are NA, Undetermined is 0, and any
  #  other values are as is
  mutate(taq_value = as.numeric(
    ifelse(taq_value == "Undetermined", 0, 
    ifelse(taq_value %in% c("Indeterminate", "Indeterminat", "."), NA, taq_value))))


# Merge ont ID_Decoder so that only humichip samples are included
taq_data_humi <- ID_Decoder %>% 
  # remove visit 4 from ID_Decoder
  filter(visit_number != 4) %>% 
  # join to taq full dat
  left_join(., taq_data_full, by = c("study_id" = "STUDYID", "visit_number")) %>% 
  # Set threshold for what is present vs not-detected
  mutate(taq_value = ifelse(is.na(taq_value), "Not Available",
                     ifelse(taq_value == 0, "No", 
                     ifelse(taq_value > 0 & taq_value < 35, "Yes", 
                     ifelse(taq_value >= 35, "No", "Not Available")))))


# Spread into format compatable with clinical metadata
taq_data_humi <- taq_data_humi %>% 
  # Remove NA Target due to missing visit 5 taq data
  filter(!is.na(Target)) %>% 
  # Spread
  spread(key = Target, value = taq_value) %>%
  # For all NA's that are added, change to "Not Available"
  mutate_at(vars(ends_with("Card")), funs(ifelse(is.na(.), "Not Available", .))) %>%
  mutate_at(vars(ends_with("Stool")), funs(ifelse(is.na(.), "Not Available", .)))


# Add Either column if detected by either stool or card
taq_data_humi_either <- taq_data_humi %>% 
  # CMY
  mutate(CMY_either = ifelse(CMY_Card == "Yes" | CMY_Stool == "Yes", "Yes",
                      ifelse(CMY_Card == "Not Available" & CMY_Stool == "Not Available", "Not Available", "No"))) %>%
                      
  # CTX
  mutate(CTX_either = ifelse(CTX_Card == "Yes" | CTX_Stool == "Yes", "Yes",
                      ifelse(CTX_Card == "Not Available" & CTX_Stool == "Not Available", "Not Available", "No"))) %>%
                      
  # KPC
  mutate(KPC_either = ifelse(KPC_Card == "Yes" | KPC_Stool == "Yes", "Yes",
                      ifelse(KPC_Card == "Not Available" & KPC_Stool == "Not Available", "Not Available", "No"))) %>%
                      
  # NDM
  mutate(NDM_either = ifelse(NDM_Card == "Yes" | NDM_Stool == "Yes", "Yes",
                      ifelse(NDM_Card == "Not Available" & NDM_Stool == "Not Available", "Not Available", "No"))) %>%
                      
  # SHV
  mutate(SHV_either = ifelse(SHV_Card == "Yes" | SHV_Stool == "Yes", "Yes",
                      ifelse(SHV_Card == "Not Available" & SHV_Stool == "Not Available", "Not Available", "No"))) %>%
                      
  # TEM
  mutate(TEM_either = ifelse(TEM_Card == "Yes" | TEM_Stool == "Yes", "Yes",
                      ifelse(TEM_Card == "Not Available" & TEM_Stool == "Not Available", "Not Available", "No"))) %>%
                      
  # Arrange
  select(study_id, visit_number, CMY_Card, CMY_Stool, CMY_either,
         CTX_Card, CTX_Stool, CTX_either, KPC_Card, KPC_Stool, KPC_either,
         NDM_Card, NDM_Stool, NDM_either, SHV_Card, SHV_Stool, SHV_either,
         TEM_Card, TEM_Stool, TEM_either)

# Make so that visit is incorporated into columns
# Reminder that there was no visit 5 card detection, just stool
taq_data_humi_either <- taq_data_humi_either %>% 
  # Change visit to V1 and V5
  mutate(visit_number = ifelse(visit_number == 1, "V1",
                        ifelse(visit_number == 5, "V5", NA))) %>% 
  gather(key = Target, value = taq_value, -study_id, -visit_number) %>% 
  unite(Target, c("Target", "visit_number")) %>%
  spread(key = Target, value = taq_value)


# Write to data/processed -----------------------------------
write_csv(x = taq_data_humi_either, path = "data/processed/Taq_tidy.csv")



