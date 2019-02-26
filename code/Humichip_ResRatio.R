######################################################################
# Name: Humichip_ResRatio.R
# Author: Ryan Johnson
# Date Created: 7 February 2019
# Purpose: Determine which genes categories are significantly altered
#          from visit 1 to visit 4/5 by calculating
#          response ratios with 95% CI using relative
#          abundance.
######################################################################

library(tidyverse)

## Load data --------------------------------------------------------------------------------------------
humichip <- suppressMessages(suppressWarnings(read_tsv("data/processed/Merged_humichip.tsv")))
treat <- suppressMessages(suppressWarnings(read_csv("data/processed/TrEAT_Clinical_Metadata_tidy.csv")))
ID_Decoder <- suppressMessages(suppressWarnings(read_csv("data/processed/ID_Decoder_Humichip.csv")))

# Only include matched isolate for visit 1 and visit 4 or 5
ID_Decoder_matched <- ID_Decoder %>%
  filter(visit_number %in% c(1,5)) %>%
  group_by(study_id) %>% 
  filter(n() == 2) 

rm(ID_Decoder)

# Merge treat to ID_Decoder
treat <- ID_Decoder_matched %>% 
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

## Convert the Humichip data to relative abundance -------------------------------------------------
humi_relabun <- humichip_filtered %>%
  # only consider probes with a gene category designation
  filter(!is.na(geneCategory)) %>%
  # Convert all values to 1 and 0
  mutate_at(vars(starts_with("X")), funs(ifelse(is.na(.), 0, 1))) %>% 
  # Convert all values to relative abundance
  mutate_at(vars(starts_with("X")), funs((./sum(.)) * 100))

rm(humichip_filtered)

## Summarize Data ----------------------------------------------------------------------------------
humi_grouped <- humi_relabun %>%
  # Select grouping column of interest and remove rest
  select(geneCategory, starts_with("X")) %>%
  # Make long
  gather(key = glomics_ID, value = rel_abun_value, -geneCategory) %>%
  # Group by category of interest
  group_by(glomics_ID, geneCategory) %>%
  # Calculate total relative abundance for each category
  summarise(category_abundance = sum(rel_abun_value)) %>%
  # Add in metadata
  left_join(., treat_filter, by = "glomics_ID") %>% 
  ungroup()
  
rm(humi_relabun)

## Calculate Response Ratio ------------------------------------------------------------------------

humi_RR <- humi_grouped %>% 
  # Group by category and visit and Treatment
  group_by(geneCategory, visit_number, Treatment) %>% 
  summarise(mean_signal = mean(category_abundance),
            sd_signal = sd(category_abundance),
            n = sum(!is.na(category_abundance))) %>% 
  
  # Spread the signal mean by visit number and Treatment
  ungroup() %>%
  group_by(geneCategory, Treatment) %>%
  spread(visit_number, mean_signal) %>% 
  
  # Rename mean columns
  rename(group1_mean = colnames(.)[length(colnames(.)) - 1],
         group2_mean = colnames(.)[length(colnames(.))]) %>% 
  
  # Spread the sd and n columns by visit
  mutate(sd_group1 = ifelse(!is.na(group1_mean), sd_signal, NA)) %>%
  mutate(sd_group2 = ifelse(!is.na(group2_mean), sd_signal, NA)) %>%
  mutate(n_group1 = ifelse(!is.na(group1_mean), n, NA)) %>%
  mutate(n_group2 = ifelse(!is.na(group2_mean), n, NA)) %>%
  select(-sd_signal, -n) %>% 
  
  # Compress NAs
  ungroup() %>%
  group_by(geneCategory, Treatment) %>%
  summarise_all(funs(sum(., na.rm = T))) %>% 
  
  # Must have at least __ observations in each subcategory
  filter(n_group1 >= 10) %>%
  filter(n_group2 >= 10) %>%
  
  # Calculate SEM for each mean
  mutate(SEM_group1 = sd_group1 / sqrt(n_group1)) %>%
  mutate(SEM_group2 = sd_group2 / sqrt(n_group2)) %>%
  
  # Calculate the Response Ratio (RR)
  mutate(RR = log(group2_mean / group1_mean)) %>%
  
  # Calculate the Standard error for the RR
  mutate(SE_RR = sqrt((SEM_group1**2 / group1_mean**2) + (SEM_group2**2 / group2_mean**2))) %>%
  
  # Calcualte the 95% confidence interval for each RR
  mutate(CI95 = abs(1.96 * SE_RR)) %>%
  
  # Add in keeper column if does not overlap 0
  mutate(keeper = ifelse(0 > (RR - CI95) & 0 < (RR + CI95), "No", "Yes")) %>%
  
  # Make labels pretty
  mutate(pretty_cat = str_to_title(geneCategory)) %>%
  mutate(pretty_cat = str_replace_all(pretty_cat, "_"," ")) %>%
  
  # Factor columns
  #mutate(pretty_cat = fct_reorder(pretty_cat, RR)) %>%
  ungroup()

## Remove any cateogories where the mean was below threshold ----------------
humi_RR_filtered <- humi_RR %>% 
  filter(group1_mean > 0.01) %>% 
  filter(group2_mean > 0.01)


## Plot ---------------------------------------------------------------------

RR_plot <- ggplot(data = humi_RR_filtered) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  
  # points and error bar
  geom_point(aes(x = pretty_cat, y = RR, color = Treatment), 
             size = 4,
             position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = RR - CI95, 
                    ymax = RR + CI95, 
                    x = pretty_cat,
                    color = Treatment),
                width = 0.25,
                position = position_dodge(width = 0.4)) +
  
  # Group labels
  annotate(geom = "text", label = "Visit 1", x = Inf, y = -Inf, hjust = 0, vjust = 1, 
           size = 5, color = "red", fontface = 2) +
  annotate(geom = "text", label = "Visit 5", x = Inf, y = Inf, hjust = 1, vjust = 1, 
           size = 5, color = "red", fontface = 2) +
  
  # plot labels
  labs(title = "Response Ratio",
       x = "Category",
       y = "Response Ratio") +
  
  theme_minimal() +
  coord_flip() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0.5)
  )

RR_plot

ggsave(plot = RR_plot, filename = "results/figures/Humi_geneCat_Treatment.png", height = 12, width = 7)

