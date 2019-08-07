
#!#!#!#!#!# INCOMPLETE CODE !#!#!#!#!#!#!#

##########################################
# Title: V1_diversity.R
# Author: Ryan Johnson
# Date Created: 5 August 2019
# Purpose: Calculate diversity statistics
#  for V1 samples
##########################################

library(tidyverse)
library(gt)
library(vegan)

## Read in and filter data ----------------------------------------------
humichip <- read_tsv("data/processed/Merged_humichip_Renormalized.tsv")
treat <- read_csv("data/processed/TrEAT_Clinical_Metadata_tidy.csv")
ID_Decoder <- read_csv("data/processed/ID_Decoder_Humichip.csv")

# Filter Humichip
humichip_filt <- humichip %>% 
  # filter for just functional probes
  filter(Gene != "STR_SPE") %>% 
  # filter out any "NA" in Gene column
  filter(!is.na(Gene)) %>% 
  # Select just glomics samples
  select(starts_with("X")) %>%
  # Replace NA's with 0
  replace(., is.na(.), 0)


## Calculate diversity metrics -------------------------------------------
# Shannon
humi_shannon <- diversity(t(humichip_filt), "shannon")

# Count probes detected
humichip_count_probe <- humichip_filt %>% 
  mutate_all(list(~ifelse(. == 0, 0, 1)))

humichip_count_probe <- colSums(humichip_count_probe)

# Convert to tibble
humi_diversity <- tibble(glomics_ID = names(humi_shannon), 
                         shannon = humi_shannon, 
                         probe_count = humichip_count_probe)


## Merge with Metadata ----------------------------------------------------

# V1 isolates
ID_Decoder_V1 <- ID_Decoder %>%
  filter(visit_number %in% c(1, 4, 5)) %>% 
  group_by(study_id) %>% 
  filter(n() == 3)

# Add treat to humi_diversity
humi_diversity_treat <- ID_Decoder_V1 %>% 
  left_join(., humi_diversity, by = "glomics_ID") %>%
  left_join(., treat, by = c("study_id" = "STUDY_ID")) %>%
  select(glomics_ID, study_id, visit_number, country, Impact, probe_count, shannon, everything())

# Filter out the LOP and PLA treatment groups
humi_diversity_treat <- humi_diversity_treat %>%
  filter(!Treatment %in% c("LOP", "PLA")) %>% 
  ungroup()


## Metrics for Table ---------------------------------------
# By visit
probe_count_visit <- humi_diversity_treat %>% 
  group_by(visit_number) %>% 
  summarise(mean_probe_count = round(mean(probe_count), 1),
            sd_probe_count = round(sd(probe_count), 1),
            mean_shannon = round(mean(shannon), 2),
            sd_shannon = round(sd(shannon), 2)) %>% 
gt() %>% 
  # Merge Columns
  cols_merge(
    col_1 = vars(mean_probe_count), 
    col_2 = vars(sd_probe_count), 
    pattern = "{1} ± {2}") %>% 
  cols_merge(
    col_1 = vars(mean_shannon), 
    col_2 = vars(sd_shannon), 
    pattern = "{1} ± {2}")


# Get probe counts by visit by Tx
probe_count_visit_tx <- humi_diversity_treat %>% 
  group_by(visit_number, Treatment) %>% 
  summarise(mean_probe_count = round(mean(probe_count), 1),
            sd_probe_count = round(sd(probe_count), 1),
            mean_shannon = round(mean(shannon), 2),
            sd_shannon = round(sd(shannon), 2)) %>% 
  gt() %>% 
  # Merge Columns
  cols_merge(
    col_1 = vars(mean_probe_count), 
    col_2 = vars(sd_probe_count), 
    pattern = "{1} ± {2}") %>% 
  cols_merge(
    col_1 = vars(mean_shannon), 
    col_2 = vars(sd_shannon), 
    pattern = "{1} ± {2}")




probe_count_visit 
probe_count_visit_tx


  

