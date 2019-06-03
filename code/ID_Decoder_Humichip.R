################################################
# Name: ID_Decoder.R
# Author: Ryan Johnson
# Date Created: 8 July 2018
# Purpose: Generate a clean table that links
#          Glomics ID's to TrEAT Subject IDs
################################################

library(tidyverse)
library(readxl)

## Read in the Decoder ------------------------------------
ID_list <- suppressMessages(suppressWarnings(read_excel("data/raw/IDCRP_Glomics_Subject_ID_List_11-21-17.xlsx")))

## Modify Decoder ---------------------------------------

# Replace column names
colnames(ID_list) <- 1:ncol(ID_list)

# Remove columns that are all NA
ID_list <- ID_list %>%
  select_if(~!all(is.na(.)))

# Stack all glomics ID values
glomics_ID <- na.omit(unlist(select_if(ID_list, is.double), 
                     use.names = FALSE))
# Stack all study ID values
study_ID <- na.omit(unlist(select_if(ID_list, is.character), 
                           use.names = FALSE))

# merge into tibble
ID_list <- tibble(glomics_ID, study_ID)

## Split study_ID to "study_ID" and "visit_numer" ---------------------------
ID_list <- ID_list %>%
  separate(study_ID, into = c("study_id", "visit_number"), "(?<=\\d{4})-") %>%
  mutate(glomics_ID = paste0("X", glomics_ID)) # Add "X" in front of glomics_ID

rm(glomics_ID, study_ID) # clean up

# Filter for samples in Humichip data
humichip_samples <- suppressMessages(suppressWarnings(read_tsv("data/processed/Merged_humichip_Renormalized.tsv", n_max = 1))) %>%
  select(starts_with("X")) %>%
  colnames()

ID_list <- ID_list %>%
  filter(glomics_ID %in% humichip_samples)


## Write to processed data -----------------------
write_csv(ID_list, "data/processed/ID_Decoder_Humichip.csv")
