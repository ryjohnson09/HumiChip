##########################################
# Title: Patient_demographics.R
# Author: Ryan Johnson
# Date Created: 20 August 2019
# Purpose: Calculate patient demographics
##########################################

library(tidyverse)
library(gt)

## Read in and filter data ----------------------------------------------
treat <- read_csv("data/processed/TrEAT_Clinical_Metadata_tidy.csv")
ID_Decoder <- read_csv("data/processed/ID_Decoder_Humichip.csv")

# Merge treat to ID_Decoder
treat <- ID_Decoder %>% 
  left_join(., treat, by = c("study_id" = "STUDY_ID"))

# Remove LOP and PLA
treat_filter <- treat %>% 
  filter(!Treatment %in% c("LOP", "PLA"))


## By Treatment Groups -------------------------------------------------
all_samples <- treat_filter %>% 
  # Ensure that all patients IDs included
  # even if they only gave a single sample
  select(-glomics_ID, -visit_number) %>% 
  unique() %>% 
  summarise(
    
    # n
    N = n(),
    
    # Age
    `Age, median (IQR)` = paste0(median(Age), " (", 
                                 quantile(Age, 0.25), "-", 
                                 quantile(Age, 0.75), ")"),
    
    # Sex
    `Male, No. (%)` = paste0(sum(Sex == 1), " (",
                             round((sum(Sex == 1) / N) * 100, 1), ")"),
    
    # Race
    `White, No. (%)` = paste0(sum(Race == "White"), " (",
                              round((sum(Race == "White") / N) * 100, 1), ")"),
    `Black, No. (%)` = paste0(sum(Race == "Black"), " (",
                              round((sum(Race == "Black") / N) * 100, 1), ")"),
    `Other, No. (%)` = paste0(sum(Race == "Other"), " (",
                              round((sum(Race == "Other") / N) * 100, 1), ")"),
    
    
    # Military Affilation
    `US Military, No. (%)` = paste0(sum(Military_affiliation == 1), " (",
                                    round((sum(Military_affiliation == 1) / N) * 100, 1), ")"),
    
    
    # Duty Station (Country)
    `Kenya, No. (%)` = paste0(sum(country == "Kenya"), " (",
                              round((sum(country == "Kenya") / N) * 100, 1), ")"),
    `Djibouti, No. (%)` = paste0(sum(country == "Djibouti"), " (",
                                 round((sum(country == "Djibouti") / N) * 100, 1), ")"),
    `Honduras, No. (%)` = paste0(sum(country == "Honduras"), " (",
                                 round((sum(country == "Honduras") / N) * 100, 1), ")"),
    
    # Disease Severity
    `Maximum LLS in 24 hours (all cases), mean (SD), No. of stools` = 
      paste0(round(mean(Maximum_number_of_loose_liquid_stools_in_any_24_hours_prior_to_presentation), 1), 
             " (", round(sd(Maximum_number_of_loose_liquid_stools_in_any_24_hours_prior_to_presentation), 1),
             ")"),
    `Total LLS prior to presentation, mean (SD), No. of stools` =
      paste0(round(mean(Number_of_loose_liquid_stools_since_the_start_of_symptoms_prior_to_presentation), 1),
             " (", round(sd(Number_of_loose_liquid_stools_since_the_start_of_symptoms_prior_to_presentation), 1),
             ")"),
    `Total LLS in 8 hours prior to presentation, mean (SD), No. of stools` = 
      paste0(round(mean(Number_of_loose_liquid_stools_in_last_8_hours_prior_to_presentation), 1),
             " (", round(sd(Number_of_loose_liquid_stools_in_last_8_hours_prior_to_presentation), 1),
             ")"),
    `No impact on activity level, No. (%)` = 
      paste0(sum(Impact_of_illness_on_activity_level == 0), " (",
             round((sum(Impact_of_illness_on_activity_level == 0) / N) * 100, 1), ")"),
    `Decreased (≤50%) impact on activity level, No. (%)` = 
      paste0(sum(Impact_of_illness_on_activity_level == 1), " (",
             round((sum(Impact_of_illness_on_activity_level == 1) / N) * 100, 1), ")"),
    `Decreased (>50%) impact on activity level, No. (%)` = 
      paste0(sum(Impact_of_illness_on_activity_level == 2), " (",
             round((sum(Impact_of_illness_on_activity_level == 2) / N) * 100, 1), ")"),
    `Complete inability to function, No. (%)` = 
      paste0(sum(Impact_of_illness_on_activity_level == 3), " (",
             round((sum(Impact_of_illness_on_activity_level == 3) / N) * 100, 1), ")"),
  )

gt_table <- all_samples %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("Characteristic") %>% 
  gt()

gt_table






## By Treatment Groups -----------------------------------------------------------------
by_tx_group <- treat_filter %>% 
  # Ensure that all patients IDs included
  # even if they only gave a single sample
  select(-glomics_ID, -visit_number) %>% 
  unique() %>% 
  group_by(Treatment) %>% 
  summarise(
    
    # n
    N = n(),
  
    # Age
    `Age, median (IQR)` = paste0(median(Age), " (", 
                            quantile(Age, 0.25), "-", 
                            quantile(Age, 0.75), ")"),
  
    # Sex
    `Male, No. (%)` = paste0(sum(Sex == 1), " (",
                            round((sum(Sex == 1) / N) * 100, 1), ")"),
    
    # Race
    `White, No. (%)` = paste0(sum(Race == "White"), " (",
                              round((sum(Race == "White") / N) * 100, 1), ")"),
    `Black, No. (%)` = paste0(sum(Race == "Black"), " (",
                              round((sum(Race == "Black") / N) * 100, 1), ")"),
    `Other, No. (%)` = paste0(sum(Race == "Other"), " (",
                              round((sum(Race == "Other") / N) * 100, 1), ")"),
    
    
    # Military Affilation
    `US Military, No. (%)` = paste0(sum(Military_affiliation == 1), " (",
                          round((sum(Military_affiliation == 1) / N) * 100, 1), ")"),
    
    
    # Duty Station (Country)
    `Kenya, No. (%)` = paste0(sum(country == "Kenya"), " (",
                          round((sum(country == "Kenya") / N) * 100, 1), ")"),
    `Djibouti, No. (%)` = paste0(sum(country == "Djibouti"), " (",
                             round((sum(country == "Djibouti") / N) * 100, 1), ")"),
    `Honduras, No. (%)` = paste0(sum(country == "Honduras"), " (",
                             round((sum(country == "Honduras") / N) * 100, 1), ")"),
    
    # Disease Severity
    `Maximum LLS in 24 hours (all cases), mean (SD), No. of stools` = 
      paste0(round(mean(Maximum_number_of_loose_liquid_stools_in_any_24_hours_prior_to_presentation), 1), 
             " (", round(sd(Maximum_number_of_loose_liquid_stools_in_any_24_hours_prior_to_presentation), 1),
             ")"),
    `Total LLS prior to presentation, mean (SD), No. of stools` =
      paste0(round(mean(Number_of_loose_liquid_stools_since_the_start_of_symptoms_prior_to_presentation), 1),
             " (", round(sd(Number_of_loose_liquid_stools_since_the_start_of_symptoms_prior_to_presentation), 1),
             ")"),
    `Total LLS in 8 hours prior to presentation, mean (SD), No. of stools` = 
      paste0(round(mean(Number_of_loose_liquid_stools_in_last_8_hours_prior_to_presentation), 1),
             " (", round(sd(Number_of_loose_liquid_stools_in_last_8_hours_prior_to_presentation), 1),
             ")"),
    `No impact on activity level, No. (%)` = 
      paste0(sum(Impact_of_illness_on_activity_level == 0), " (",
        round((sum(Impact_of_illness_on_activity_level == 0) / N) * 100, 1), ")"),
    `Decreased (≤50%) impact on activity level, No. (%)` = 
      paste0(sum(Impact_of_illness_on_activity_level == 1), " (",
             round((sum(Impact_of_illness_on_activity_level == 1) / N) * 100, 1), ")"),
    `Decreased (>50%) impact on activity level, No. (%)` = 
      paste0(sum(Impact_of_illness_on_activity_level == 2), " (",
             round((sum(Impact_of_illness_on_activity_level == 2) / N) * 100, 1), ")"),
    `Complete inability to function, No. (%)` = 
      paste0(sum(Impact_of_illness_on_activity_level == 3), " (",
             round((sum(Impact_of_illness_on_activity_level == 3) / N) * 100, 1), ")"),
  )

gt_table_tx_groups <- by_tx_group %>% 
  column_to_rownames("Treatment") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("Characteristic") %>% 
  gt()

gt_table_tx_groups


