####################################################
# Name: Merged_Humichip_Tidy.R
# Author: Ryan Johnson
# Date Created: 8 July 2018
# Purpose: Take the merged Humichip data and convert
#          it to long/tidy format
####################################################

library(tidyverse)

# Read in merged Humichip Data
Humichip <- read_tsv("data/processed/Merged_Humichip.tsv", guess_max = 100000, progress = TRUE)

# Make long/tidy
Humichip_long <- Humichip %>%
  gather(key = glomics_ID, 
         value = Signal, 
         colnames(.[8:ncol(.)]), -Genbank.ID, -gene,
         -species, -lineage, -annotation, -geneCategory, -subcategory1, 
         -subcategory2)

rm(Humichip) # Clean

# Read in decoder
ID_list <- read.csv("data/processed/ID_Decoder.csv")

# Join to Humichip data
Humichip_long <- Humichip_long %>%
  left_join(., ID_list, by = c("glomics_ID")) %>%
  select(glomics_ID, study_id, visit_number, everything())

# Write to `data/processed`
write_tsv(Humichip_long, "data/processed/Merged_Humichip_Tidy.tsv")
