####################################################
# Name: Merge_Humichip.R
# Author: Ryan Johnson
# Date Created: 23 July, 2018
# Purpose: Take all Humichip data sets and combine 
#          them into a single data frame
####################################################

library(tidyverse)

# list vector of humichip files
humichip_files <- c(
  "data/raw/HumiChip_New/HumiChip-1-LTO.txt",
  "data/raw/HumiChip_New/HumiChip-2-LTO.txt",
  "data/raw/HumiChip_New/HumiChip-3-LTO.txt", 
  "data/raw/HumiChip_New/HumiChip-4-LTO.txt", 
  "data/raw/HumiChip_New/HumiChip-5-LTO.txt"
)

# Create empty data frame
humichip_data <- tibble()

# Start loop that will read in each humichip separately
for (hchip in humichip_files){
  
  # If tibble is empty (first occurence)
  if(is_empty(humichip_data)){
    
    # Read in hchip
    humichip_data <- read_tsv(hchip, guess_max = 100000) %>%
      select(`Genbank ID`, gene, species, lineage, annotation, geneCategory, 
             subcategory1, subcategory2, starts_with("X"))
  } else {
    
    # Read in hchip and merge into humichip_data
    humichip_temp <-  read_tsv(hchip, guess_max = 100000) %>%
      select(`Genbank ID`, gene, species, lineage, annotation, geneCategory, 
             subcategory1, subcategory2, starts_with("X"))
    
    humichip_data <- full_join(humichip_data, humichip_temp, 
                              by = c("Genbank ID", "gene", "species", "lineage", "annotation", 
                                     "geneCategory", "subcategory1", "subcategory2"))
  }
}

# Clean
rm(humichip_temp, hchip, humichip_files)


# Remove any extraneous text at end of sample headers
colnames(humichip_data) <- str_replace(string = colnames(humichip_data), 
                                        pattern = "(^[Xx]\\d{1,3}).*", 
                                        replacement = "\\1")

# Capitalize "x" if lowercase
colnames(humichip_data) <- str_replace(string = colnames(humichip_data),
                                      pattern = "^x",
                                      replacement = "X")



## Make sample names Unique ---------------------------------------------------

# Check for duplicates and make unique
humichip_data <- humichip_data %>%
  setNames(make.names(names(.), unique = TRUE))


#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
# Remove Duplicates ****MAY NEED TO BE CHANGED******* 
humichip_data <- humichip_data %>%
  select(-X75, -X193, -X248)

colnames(humichip_data) <- str_replace(string = colnames(humichip_data), 
                                      pattern = "\\.1", 
                                      replacement = "")
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#




## Remove samples that are excluded from the study ----------------------------
humichip_data <- humichip_data %>%
  select(-X10) # 48-2406




## Fixing Category names ------------------------------------------------------

# geneCategory
humichip_data$geneCategory <- toupper(humichip_data$geneCategory)
humichip_data$geneCategory <- gsub(pattern = " ", replacement = "_", humichip_data$geneCategory)

# annotation
humichip_data$annotation <- toupper(humichip_data$annotation)
humichip_data$annotation <- gsub(pattern = " ", replacement = "_", humichip_data$annotation)

# subcategory1
humichip_data$subcategory1 <- toupper(humichip_data$subcategory1)
humichip_data$subcategory1 <- gsub(pattern = " ", replacement = "_", humichip_data$subcategory1)

# subcategory2
humichip_data$subcategory2 <- toupper(humichip_data$subcategory2)
humichip_data$subcategory2 <- gsub(pattern = " ", replacement = "_", humichip_data$subcategory2)




# Return compiled data frame
write_tsv(humichip_data, "data/processed/Merged_humichip.tsv")

