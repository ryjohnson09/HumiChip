# Name: array_run.R
# Author: Ryan Johnson
# Date: 8 May 2019
# Purpose: Assign each glomics ID to whichever array "run" it came from.
#  Can be used to assess changes in arrays.

## Load Libraries ---------------------------------
library(tidyverse)

## Load Data -------------------------------------
humichip_files <- c(
  "data/raw/HumiChip_New/HumiChip-1-LTO.txt",
  "data/raw/HumiChip_New/HumiChip-2-LTO.txt",
  "data/raw/HumiChip_New/HumiChip-3-LTO.txt", 
  "data/raw/HumiChip_New/HumiChip-4-LTO.txt", 
  "data/raw/HumiChip_New/HumiChip-5-LTO.txt"
)


# Start loop that pulls the colnames from each humichip file
for (hchip in humichip_files){
  # Read in hchip
  humichip_data <- read_tsv(hchip, n_max = 1) %>%
    select(starts_with("X")) %>% 
    gather() %>% 
    pull(key)
  
  # Remove any extraneous text at end of sample headers
  humichip_data <- str_replace(string = humichip_data, 
                                        pattern = "(^[Xx]\\d{1,3}).*", 
                                        replacement = "\\1")
  
  # Capitalize "x" if lowercase
  humichip_data <- str_replace(string = humichip_data,
                                        pattern = "^x",
                                        replacement = "X")
  
  assign(basename(hchip), humichip_data)
}

## Remove duplicates and excluded samples ----------------------------------------

# This comes from the Merge_Humichip.R script
# Remove first X75
`HumiChip-2-LTO.txt` <- `HumiChip-2-LTO.txt`[!`HumiChip-2-LTO.txt` == "X75"]
# Remove second X193
`HumiChip-2-LTO.txt` <- `HumiChip-2-LTO.txt`[!`HumiChip-2-LTO.txt` == "X193"]
# Remove first X278
`HumiChip-1-LTO.txt` <- `HumiChip-1-LTO.txt`[!`HumiChip-1-LTO.txt` == "X248"]
# Remove X10
`HumiChip-1-LTO.txt` <- `HumiChip-1-LTO.txt`[!`HumiChip-1-LTO.txt` == "X10"]

## Create Data Frame --------------------------------------------
run_data <- data.frame("glomics_ID" = c(`HumiChip-1-LTO.txt`, `HumiChip-2-LTO.txt`, `HumiChip-3-LTO.txt`,
                                        `HumiChip-4-LTO.txt`, `HumiChip-5-LTO.txt`),
                       "Array_Run" = c(rep("Run1", length(`HumiChip-1-LTO.txt`)),
                                       rep("Run2", length(`HumiChip-2-LTO.txt`)),
                                       rep("Run3", length(`HumiChip-3-LTO.txt`)),
                                       rep("Run4", length(`HumiChip-4-LTO.txt`)),
                                       rep("Run5", length(`HumiChip-5-LTO.txt`))))

## Write to data/processed
write_csv(run_data, path = "data/processed/array_runs.csv")
