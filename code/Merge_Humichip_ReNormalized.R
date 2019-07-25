####################################################
# Name: Merge_Humichip_ReNormalized.R
# Author: Ryan Johnson
# Date Created: 13 May, 2019
# Purpose: Take all Humichip data sets and combine 
#  them into a single data frame. This script
#  tidy's up the new renormalized data (13 May 2019)
####################################################

library(tidyverse)

# vector of humichip files
humichip_files <- c(
  "data/raw/HumiChip_ReNormalized_Updated/Norm.Riddle.HumiChip.All.Nooutlier.Log.txt"
)

# Read in hchip
humichip_data <- read_tsv(humichip_files, col_types = cols(uniqueID = col_character(),
                                                  gene = col_character(),
                                                  proteinGI = col_character(),
                                                  annotation = col_character(),
                                                  species = col_character(),
                                                  lineage = col_character(),
                                                  geneCategory = col_character(),
                                                  subcategory1 = col_character(),
                                                  subcategory2 = col_character(),
                                                  .default = col_double())) %>%
  select(uniqueID, gene, species, lineage, geneCategory, subcategory1, subcategory2, starts_with("X"))

# Rename to make compatable with code
colnames(humichip_data)[1:7] <- c("Genbank ID", "Gene", "Organism", "Lineage", "Gene_category", 
                                  "Subcategory1", "Subcategory2")
  

# Clean
rm(humichip_files)



# Check that the Genbank ID column has no duplicates
ensure_no_dups <- sum(duplicated(humichip_data$`Genbank ID`) + 
                        duplicated(humichip_data$`Genbank ID`, fromLast = TRUE))

# Stop execution if duplicates found
if(ensure_no_dups != 0){
  stop("Non unique probe identifiers!")
}







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
# Remove Duplicates ****Remove sample with fewer values******* 
humichip_data <- humichip_data %>%
  select(-X75, -X193.1, -X248)

colnames(humichip_data) <- str_replace(string = colnames(humichip_data), 
                                       pattern = "\\.1", 
                                       replacement = "")
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#




## Remove samples that are excluded from the study ----------------------------
humichip_data <- humichip_data %>%
  select(-X10) %>%  # 48-2406 (Removed due to pill vomitting)
  select(-X66) %>% # Removed due to low biomass
  select(-X319) # Removed due to low biomass




## Fixing Category names ------------------------------------------------------

# Gene_category
humichip_data$Gene_category <- toupper(humichip_data$Gene_category)
humichip_data$Gene_category <- gsub(pattern = " ", replacement = "_", humichip_data$Gene_category)

# annotation
# humichip_data$annotation <- toupper(humichip_data$annotation)
# humichip_data$annotation <- gsub(pattern = " ", replacement = "_", humichip_data$annotation)

# Subcategory1
humichip_data$Subcategory1 <- toupper(humichip_data$Subcategory1)
humichip_data$Subcategory1 <- gsub(pattern = " ", replacement = "_", humichip_data$Subcategory1)

# Subcategory2
humichip_data$Subcategory2 <- toupper(humichip_data$Subcategory2)
humichip_data$Subcategory2 <- gsub(pattern = " ", replacement = "_", humichip_data$Subcategory2)




# Return compiled data frame
write_tsv(humichip_data, "data/processed/Merged_humichip_Renormalized.tsv")

