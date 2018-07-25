#############################################################
# Name: Functional_Category_Visit_Number.R
# Author: Ryan Johnson
# Date Created: 11 July 2018
# Purpose: Visualize differences in gene categories
#           abundance between the various visit numbers
#############################################################

library(tidyverse)

## Read in Merged Humichip Data --------------------------------------------------------------------------
Humichip <- read_tsv("data/processed/Merged_Humichip_Tidy.tsv", 
                    col_types = "cciccccccccd", 
                    progress = TRUE)



## Convert to Relative abundance and summarize by geneCategory and visit_number ------------------------
Humichip_RA_visitnumber <- Humichip %>%
  
  # Remove the STR_SPE probes
  filter(gene != "STR_SPE") %>%
  
  group_by(glomics_ID) %>%
  mutate(Signal_Relative_Abundance = Signal / sum(Signal, na.rm = TRUE)) %>% # Relative Abundance Calc
  group_by(geneCategory, glomics_ID, visit_number) %>%
  summarise(geneCategory_relative_abundance = sum(Signal_Relative_Abundance, na.rm = TRUE)) %>%
  group_by(geneCategory, visit_number) %>%
  summarise(mean_sig = mean(geneCategory_relative_abundance),
            sd_sig = sd(geneCategory_relative_abundance))


## Create Plot -------------------------------------------------------------------------------------------

Humichip_plot <- ggplot(Humichip_RA_visitnumber, aes(x = geneCategory, 
                                                   y = mean_sig, 
                                                   fill = factor(visit_number))) +
  geom_errorbar(aes(ymax = mean_sig + sd_sig,
                    ymin = mean_sig - (mean_sig/2)), 
                width=.2, 
                position=position_dodge(.9)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  scale_fill_manual(name = "Visit Number", values = c('#d7191c','#fdae61','#abd9e9')) +
  
  ylab("Relative Abundance") +
  xlab("Gene Category") +
  ggtitle("Relative Abundance of Gene Categories by Visit") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14), 
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 14), 
    plot.title = element_text(size = 16, face = "bold"),
    plot.margin=unit(c(0.5,0.5,0.5,5),"cm")
  )

## Save Plot ----------------------------------------------------------------------------------------------
ggsave("results/figures/Humichip_Gene_Abundance_Visit_Number.png", width = 12, height = 7)