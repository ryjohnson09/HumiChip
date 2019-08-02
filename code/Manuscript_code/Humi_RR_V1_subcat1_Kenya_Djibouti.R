######################################################################
# Name: Humi_RR_V1_subcat1_Kenya_Djibouti.R
# Author: Ryan Johnson
# Date Created: June 11 2019
# Purpose: Determine which genes categories are significantly altered.
#   This script is highly ammendable depending on what you are looking
#   at (i.e. country, visit, category, etc).
######################################################################

library(tidyverse)

countries <- c("Kenya", "Djibouti")
visits <- c(1)

## Load data --------------------------------------------------------------------------------------------
humichip <- suppressMessages(suppressWarnings(read_tsv("data/processed/Merged_humichip_Renormalized.tsv")))
treat <- suppressMessages(suppressWarnings(read_csv("data/processed/TrEAT_Clinical_Metadata_tidy.csv")))
ID_Decoder <- suppressMessages(suppressWarnings(read_csv("data/processed/ID_Decoder_Humichip.csv")))

# Only include matched isolates
ID_Decoder_matched <- ID_Decoder %>%
  filter(visit_number %in% visits) %>%
  group_by(study_id) %>% 
  filter(n() == length(visits)) 


# Merge treat to ID_Decoder
treat <- ID_Decoder_matched %>% 
  left_join(., treat, by = c("study_id" = "STUDY_ID"))


## Filter treat to include treat samples of interest -----------------------------------------------

treat_filter <- treat %>% 
  # Remove LOP and PLA
  filter(!Treatment %in% c("LOP", "PLA"))


## Filter the humichip data to include only samples in treat_filter --------------------------------
humichip_filtered <- humichip %>% 
  select_at(c("Genbank.ID", "Gene", "Organism", "Lineage",
              "Gene_category", "Subcategory1",
              "Subcategory2", treat_filter$glomics_ID))


## Convert the Humichip data to relative abundance -------------------------------------------------
humi_relabun <- humichip_filtered %>%
  # only consider probes with a gene category designation (removes STR_SPE probes)
  filter(!is.na(Gene_category)) %>%
  # Convert all NA values to 0
  mutate_at(vars(starts_with("X")), list(~ifelse(is.na(.), 0, .))) %>% 
  # Convert all values to relative abundance
  mutate_at(vars(starts_with("X")), list(~(./sum(.)) * 100))


## Summarize Data ----------------------------------------------------------------------------------
humi_grouped <- humi_relabun %>%
  # Select specific category?
  #filter(Subcategory1 == "ANTIBIOTIC_RESISTANCE") %>% 
  # Select grouping column of interest and remove rest
  select(Gene_category, Subcategory1, starts_with("X")) %>%
  # Make long
  gather(key = glomics_ID, value = rel_abun_value, -Subcategory1, -Gene_category) %>%
  # Group by category of interest
  group_by(glomics_ID, Gene_category, Subcategory1) %>%
  # Calculate total relative abundance for each category
  summarise(category_abundance = sum(rel_abun_value)) %>%
  # Add in metadata
  left_join(., treat_filter, by = "glomics_ID") %>% 
  ungroup()


## Select only certain countries at a time ------------------------------------------------------------
humi_grouped <- humi_grouped %>%
  filter(country %in% countries)

## Calculate Response Ratio ------------------------------------------------------------------------

humi_RR <- humi_grouped %>% 
  # Group by Subcategory1 and country
  group_by(Gene_category, Subcategory1, country) %>% 
  summarise(mean_signal = mean(category_abundance),
            sd_signal = sd(category_abundance),
            n = sum(!is.na(category_abundance))) %>% 
  
  # Spread the signal mean by country
  ungroup() %>%
  group_by(Gene_category, Subcategory1, country) %>%
  spread(country, mean_signal) %>% 
  
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
  group_by(Gene_category, Subcategory1) %>%
  summarise_all(funs(sum(., na.rm = T))) %>% 
  
  # Must have at least __ observations in each subcategory
  filter(n_group1 >= 7) %>%
  filter(n_group2 >= 7) %>%
  
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
  mutate(pretty_cat = str_to_title(Subcategory1)) %>%
  mutate(pretty_genecat = str_to_title(Gene_category)) %>% 
  mutate(pretty_cat = str_replace_all(pretty_cat, "_"," ")) %>%
  mutate(pretty_genecat = str_replace_all(pretty_genecat, "_"," ")) %>%
  
  # Factor columns
  #mutate(pretty_cat = fct_reorder(pretty_cat, RR)) %>%
  ungroup()

## Remove any cateogories where the mean was below threshold ----------------
humi_RR_filtered <- humi_RR %>% 
  filter(group1_mean > 0.01) %>% 
  filter(group2_mean > 0.01) %>% 
  # Remove any NA categories from subcats
  filter(!is.na(Subcategory1)) %>% 
  # Remove an non-significant categories?
  filter(keeper == "Yes")


## Plot ---------------------------------------------------------------------

RR_plot <- ggplot(data = humi_RR_filtered) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  
  # points and error bar
  geom_point(aes(x = pretty_cat, y = RR), 
             size = 3,
             position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = RR - CI95, 
                    ymax = RR + CI95, 
                    x = pretty_cat),
                width = 0.25, 
                position = position_dodge(width = 0.4)) +
  
  # Group labels
  # annotate(geom = "text", label = "country1", x = Inf, y = -Inf, hjust = 0, vjust = 1,
  #          size = 5, color = "red", fontface = 2) +
  # annotate(geom = "text", label = "country2", x = Inf, y = Inf, hjust = 1, vjust = 1,
  #          size = 5, color = "red", fontface = 2) +
  
  # plot labels
  labs(y = "Response Ratio",
       x = "") +
  scale_x_discrete(position = "left") +
  
  theme_linedraw() +
  facet_wrap(~pretty_genecat, 
             scales = "free_y", ncol = 1) +
  coord_flip() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0.5),
    strip.text = element_text(size = 12)
  )

RR_plot

ggsave(plot = RR_plot, filename = "results/figures/Manuscript_Figures/RR_V1_subcat1_Kenya_Djibouti.png", height = 12, width = 7)

