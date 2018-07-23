#############################################################
# Name: Humichip_Signal_Distribution.R
# Author: Ryan Johnson
# Date Created: 10 July 2018
# Purpose: Visualize the distribution of signal intensities
#          for the merged Humichip Data
#############################################################

library(tidyverse)

# Read in only the Signal Column from Humichip Data
Humichip <- read_tsv("data/processed/Merged_Humichip_Tidy.tsv", progress = TRUE,
                    col_types = cols_only(Signal = "d")) %>% na.omit()

# Create Plot
Humichip_plot <- ggplot(Humichip, aes(x = Signal)) +
  geom_histogram(binwidth = 0.1) +
  ggtitle("Signal Distribution") +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold")
  )

ggsave("results/figures/Humichip_Signal_Distribution.png")

