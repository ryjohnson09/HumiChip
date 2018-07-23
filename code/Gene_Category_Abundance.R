#######################################################
# Name: Functional_Category_Abundance.R
# Author: Ryan Johnson
# Date Created: 10 July 2018
# Purpose: Plot the relative abundance of functional 
#           gene categories for all Humichip samples
#######################################################

library(tidyverse)

## Read in Merged Humichip Data -----------------------------
Humichip <- read_tsv("data/processed/Merged_Humichip_Tidy.tsv", 
                     col_types = "cciccccccccd", 
                     progress = TRUE)




## Convert to Relative Abundance and group by Gene Category --------------------------------
Humichip_RA <- Humichip %>%
  group_by(glomics_ID) %>%
  mutate(Signal_Relative_Abundance = Signal / sum(Signal, na.rm = TRUE)) %>% # convert to relative abundance
  group_by(geneCategory, glomics_ID) %>%
  
  # Add all relative abundance values together per gene category
  summarise(geneCategory_relative_abundance = sum(Signal_Relative_Abundance, na.rm = TRUE)) %>%
  group_by(geneCategory) %>%
  
  # Report mean gene category relative abundance with standard deviation
  summarise(mean_sig = mean(geneCategory_relative_abundance),
            sd_sig = sd(geneCategory_relative_abundance))




## Creat Plot ---------------------------------------------------------
Humichip_plot <- ggplot(Humichip_RA, aes(x = geneCategory, 
                                             y = mean_sig,
                                             fill = geneCategory)) +
  
  # Bar Plot
  geom_bar(stat = "identity", color = "black") +
  
  # Error Bar
  geom_errorbar(aes(ymin = mean_sig - sd_sig, 
                    ymax = mean_sig + sd_sig), 
                width=.2, 
                position=position_dodge(.9)) +
  
  # Aesthetics
  scale_fill_discrete(name = "Glomics ID") +
  ylab("Relative Abundance") +
  xlab("Gene Category") +
  ggtitle("Relative Abundance of HumiChip Gene Categories") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 14), 
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14), 
    plot.title = element_text(size = 16, face = "bold"),
    plot.margin = unit(c(0.5,0.5,0.5,5), "cm")
  )

Humichip_plot

## Save Plot ----------------------------------------------------------------------------------------------
ggsave("results/figures/Humichip_Gene_Abundance.png", width = 12, height = 7)