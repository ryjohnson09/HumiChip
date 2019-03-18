#######################################################
# Name: Functional_Category_Country.R
# Author: Ryan Johnson
# Date Created: 11 March 2019
# Purpose: Plot the relative abundance of functional 
#           gene categories for all Humichip samples
#           by country
#######################################################

library(tidyverse)

probe_type <- "Functional"


## Read in data -----------------------------------------------------------------------------------------
humichip <- suppressWarnings(suppressMessages(read_tsv("data/processed/Merged_humichip.tsv")))
treat <- suppressWarnings(suppressMessages(read_csv("data/processed/TrEAT_Clinical_Metadata_tidy.csv")))
ID_Decoder <- suppressWarnings(suppressMessages(read_csv("data/processed/ID_Decoder_Humichip.csv")))

# Only include matched isolate for visit 1 and visit 4 or 5
ID_Decoder_select <- ID_Decoder %>%
  filter(visit_number %in% c(1)) %>%
  group_by(study_id) %>% 
  filter(n() == 1) 

rm(ID_Decoder)

# Merge treat to ID_Decoder
treat <- ID_Decoder_select %>% 
  left_join(., treat, by = c("study_id" = "STUDY_ID"))


## Filter treat to include treat samples of interest -----------------------------------------------

treat_filter <- treat %>% 
  # Remove LOP and PLA
  filter(!Treatment %in% c("LOP", "PLA"))

rm(treat)

## Filter the humichip data to include only samples in treat_filter --------------------------------
humichip_filtered <- humichip %>% 
  select_if(colnames(.) %in% c("Genbank.ID", "gene", "species", "lineage",
                               "annotation", "geneCategory", "subcategory1",
                               "subcategory2", treat_filter$glomics_ID))

rm(humichip)

## Probe Filtering ---------------------------------------------------------
if(probe_type == "Functional"){
  humi_probes_filtered <- humichip_filtered %>%
    filter(gene != "STR_SPE")
} else if (probe_type == "Strain/Species"){
  humi_probes_filtered <- humichip_filtered %>%
    filter(gene == "STR_SPE")
} else if (probe_type == "All"){
  humi_probes_filtered <- humichip_filtered
} else {
  stop("Error filter probes")
}






#!#!#!#! Below needs fixing #!#!#!!#!

## Convert to Relative Abundance and group by Gene Category --------------------------------
humichip_RA <- humi_probes_filtered %>%
  
  # Select columns of interest
  select(geneCategory, starts_with("X")) %>% 
  
  # if na convert to 0
  mutate_at(vars(starts_with("X")), funs(ifelse(is.na(.), 0, .))) %>% 
  
  
  group_by(geneCategory) %>%
  mutate_at(vars(starts_with("X")), funs(sum(.))) %>% # convert to relative abundance
  distinct() %>% 
  ungroup() %>% 
  
  # Make tidy
  gather(glomics_ID, Signal, -geneCategory) %>% 
  
  # Convert to Relative Abundance
  group_by(glomics_ID) %>% 
  mutate(probe_total = sum(Signal)) %>% 
  mutate(percent_abundance = (Signal / probe_total) * 100) %>% 
  ungroup() %>% 
  
  # Add in treat metadata
  left_join(., treat_filter, by = "glomics_ID") %>% 
  select(glomics_ID, visit_number, country, geneCategory, percent_abundance) %>% 
  
  # Summary stats
  group_by(country, geneCategory) %>% 
  summarise(mean_relabun = mean(percent_abundance),
              sd_relabun = sd(percent_abundance)) %>% 
  
  # Make labels pretty
  mutate(pretty_cat = str_to_title(geneCategory)) %>%
  mutate(pretty_cat = str_replace_all(pretty_cat, "_"," "))



## Creat Plot ---------------------------------------------------------
humichip_plot <- ggplot(humichip_RA, aes(x = pretty_cat, 
                                         y = mean_relabun,
                                         fill = country)) +
  
  # Bar Plot
  geom_bar(stat = "identity", color = "black", position = position_dodge(0.9)) +
  
  # Error Bar
  geom_errorbar(aes(ymin = mean_relabun - sd_relabun, 
                    ymax = mean_relabun + sd_relabun), 
                width=.2, 
                position=position_dodge(0.9)) +
  
  # Aesthetics
  scale_fill_discrete(name = "Country") +
  ylab("Relative Abundance") +
  xlab("Gene Category") +
  ggtitle("Relative Abundance of HumiChip Gene Categories") +
  theme_minimal() +
  theme(
    #legend.position = "none",
    axis.title.x = element_text(size = 14), 
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14), 
    plot.title = element_text(size = 16, face = "bold"),
    plot.margin = unit(c(0.5,0.5,0.5,5), "cm")
  )

humichip_plot

## Save Plot ----------------------------------------------------------------------------------------------
ggsave("results/figures/Humichip_Gene_Abundance.png", width = 12, height = 7)
