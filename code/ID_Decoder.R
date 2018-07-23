################################################
# Name: ID_Decoder.R
# Author: Ryan Johnson
# Date Created: 8 July 2018
# Purpose: Generate a clean table that links
#          Glomics ID's to TrEAT Subject IDs
################################################

library(tidyverse)
library(readxl)

## Read in the GeoChip Decoder -----------------------------
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


## Write to processed data -----------------------
write_csv(ID_list, "data/processed/ID_Decoder.csv")