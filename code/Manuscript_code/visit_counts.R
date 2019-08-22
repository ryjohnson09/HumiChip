##########################################
# Title: visit_counts.R
# Author: Ryan Johnson
# Date Created: 21 August 2019
# Purpose: Calculate patient visits
##########################################

library(tidyverse)

# Read in data
treat <- read_csv("data/processed/TrEAT_Clinical_Metadata_tidy.csv")
ID_Decoder <- read_csv("data/processed/ID_Decoder_Humichip.csv")

# Merge treat to ID_Decoder
treat <- ID_Decoder %>% 
  left_join(., treat, by = c("study_id" = "STUDY_ID"))

# Remove LOP and PLA
treat_filter <- treat %>% 
  filter(!Treatment %in% c("LOP", "PLA"))

# Get number of matched vs non-matched
treat_number <- treat_filter %>% 
  select(study_id, visit_number)


# Generate overlap/non-overlap numbers
visit_1 <- treat_number %>% 
  group_by(study_id) %>% 
  filter(n() == 1) %>% 
  filter(visit_number == 1) %>% 
  nrow()

visit_4 <- treat_number %>% 
  group_by(study_id) %>% 
  filter(n() == 1) %>% 
  filter(visit_number == 4) %>% 
  nrow()

visit_5 <- treat_number %>% 
  group_by(study_id) %>% 
  filter(n() == 1) %>% 
  filter(visit_number == 5) %>% 
  nrow()

visit_14 <- treat_number %>% 
  group_by(study_id) %>% 
  filter(n() == 2) %>% 
  filter(visit_number %in% c(1,4)) %>% 
  filter(n() == 2) %>% 
  pull(study_id) %>% 
  unique() %>% 
  length()

visit_15 <- treat_number %>% 
  group_by(study_id) %>% 
  filter(n() == 2) %>% 
  filter(visit_number %in% c(1,5)) %>% 
  filter(n() == 2) %>% 
  pull(study_id) %>% 
  unique() %>% 
  length()

visit_45 <- treat_number %>% 
  group_by(study_id) %>% 
  filter(n() == 2) %>% 
  filter(visit_number %in% c(4,5)) %>% 
  filter(n() == 2) %>% 
  pull(study_id) %>% 
  unique() %>% 
  length()

visit_145 <- treat_number %>% 
  group_by(study_id) %>% 
  filter(n() == 3) %>% 
  pull(study_id) %>% 
  unique() %>% 
  length()

# Create final table
visit_groups <- tibble(visit_1, visit_4, visit_5, 
                       visit_14, visit_15, visit_45,
                       visit_145)
