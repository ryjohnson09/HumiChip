######################################################################
# Name: Humichip_ResRatio_Visit.R
# Author: Ryan Johnson
# Date Created: 24 July 2018
# Purpose: Determine which genes categories are significantly altered
#          from visit 1 to visit 4 and 5 by calculating
#          response ratios with 95% CI using relative
#          abundance.
######################################################################

library(tidyverse)

# Read in merged Humichip Data
Humichip <- read_tsv("data/processed/Merged_Humichip_Tidy.tsv", 
                     col_types = "cciccccccccd",
                     progress = TRUE)

Humichip_RR <- Humichip %>%
  
  # Remove the STR_SPE probes
  filter(gene != "STR_SPE") %>%
  
  group_by(glomics_ID) %>%
  mutate(Signal_Relative_Abundance = (Signal / sum(Signal, na.rm = TRUE)* 100)) %>%
  
  # Remove columns not needed
  select(-gene, -Genbank.ID, -subcategory1,  
         -Signal, -species, -subcategory2, -annotation, -lineage) %>%
  
  # Remove any rows with NA in the signal category or geneCategory
  filter(!is.na(Signal_Relative_Abundance)) %>%
  filter(!is.na(geneCategory)) %>%

  
  # Calculate mead, sd, and counts (n)
  group_by(geneCategory, glomics_ID, visit_number) %>%
  summarise(geneCategory_relative_abundance = sum(Signal_Relative_Abundance, na.rm = TRUE)) %>%
  group_by(geneCategory, visit_number) %>%
  summarise(mean_signal = mean(geneCategory_relative_abundance),
            sd_signal = sd(geneCategory_relative_abundance),
            n = sum(!is.na(geneCategory_relative_abundance))) %>%


  # Spread the signal mean by visit number
  group_by(geneCategory) %>%
  spread(visit_number, mean_signal) %>%
  
  # Rename visit mean columns
  rename(Visit1_mean = `1`,
         Visit4_mean = `4`,
         Visit5_mean = `5`) %>% 
  
  # Spread the sd and n columns by visit
  mutate(sd_visit1 = ifelse(!is.na(Visit1_mean), sd_signal, NA)) %>%
  mutate(sd_visit4 = ifelse(!is.na(Visit4_mean), sd_signal, NA)) %>%
  mutate(sd_visit5 = ifelse(!is.na(Visit5_mean), sd_signal, NA)) %>%
  mutate(n_visit1 = ifelse(!is.na(Visit1_mean), n, NA)) %>%
  mutate(n_visit4 = ifelse(!is.na(Visit4_mean), n, NA)) %>%
  mutate(n_visit5 = ifelse(!is.na(Visit5_mean), n, NA)) %>%
  select(-sd_signal, -n) %>%
  

  # Compress NAs
  group_by(geneCategory) %>%
  summarise_all(funs(sum(., na.rm = T))) %>%
  
  # Must have at least 10 observations in each subcategory
  filter(n_visit1 >= 10) %>%
  filter(n_visit4 >= 10) %>%
  filter(n_visit5 >= 10) %>%
  
  # Calculate SEM for each mean
  mutate(SEM_visit1 = sd_visit1 / sqrt(n_visit1)) %>%
  mutate(SEM_visit4 = sd_visit4 / sqrt(n_visit4)) %>%
  mutate(SEM_visit5 = sd_visit5 / sqrt(n_visit5)) %>%
  
  # Calculate the Response Ratio (RR) for each gene (V4-V1, V5-V1)
  mutate(RR_41 = log(Visit4_mean / Visit1_mean)) %>%
  mutate(RR_51 = log(Visit5_mean / Visit1_mean)) %>%
  
  # Calculate the Standard error for the RR
  mutate(SE_RR_41 = sqrt((SEM_visit1**2 / Visit1_mean**2) + (SEM_visit4**2 / Visit4_mean**2))) %>%
  mutate(SE_RR_51 = sqrt((SEM_visit1**2 / Visit1_mean**2) + (SEM_visit5**2 / Visit5_mean**2))) %>%
  
  # Calcualte the 95% confidence interval for each RR
  mutate(CI95_41 = abs(1.96 * SE_RR_41)) %>%
  mutate(CI95_51 = abs(1.96 * SE_RR_51))


# Plot
Humichip_RR_plot <- ggplot(data = arrange(Humichip_RR, desc(RR_41))) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  # 1 v 4
  geom_point(aes(x = RR_41, y = geneCategory), size = 4, color = "black") +
  geom_errorbarh(aes(xmin = RR_41 - CI95_41, xmax = RR_41 + CI95_41, y = geneCategory), 
                 color = "black") +
  # 1 v 5
  geom_point(aes(x = RR_51, y = geneCategory), size = 4, color = "red") +
  geom_errorbarh(aes(xmin = RR_51 - CI95_51, xmax = RR_51 + CI95_51, y = geneCategory), 
                 color = "red") +
  ggtitle("Response Ratio: Gene Categories by Visit\nBlack = Visit 4 vs Visit 1\nRed = Visit 5 vs Visit 1") +
  ylab("Gene Category") +
  xlab("Response Ratio") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    plot.margin = unit(c(0.5, 2, 0.5, 0.5), "cm")
  )


ggsave("results/figures/Humichip_RespRatio_Visit.png", height = 10, width = 10)
