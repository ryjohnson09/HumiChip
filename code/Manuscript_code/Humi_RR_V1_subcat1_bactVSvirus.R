######################################################################
# Name: Humi_RR_V1_subcat1_bactVSvirus.R
# Author: Ryan Johnson
# Date Created: 20 Aug 2019
# Purpose: Functional probe RR to determine functional
#  differences between bacterial vs viral infections @ visit 1
######################################################################

library(tidyverse)

visits <- c(1)
detection_method <- "culture"

## Load data --------------------------------------------------------------------------------------------
humichip <- suppressMessages(suppressWarnings(read_tsv("data/processed/Merged_humichip_Renormalized.tsv")))
treat <- suppressMessages(suppressWarnings(read_csv("data/processed/TrEAT_Clinical_Metadata_tidy.csv")))
ID_Decoder <- suppressMessages(suppressWarnings(read_csv("data/processed/ID_Decoder_Humichip.csv")))

# Extract visit 1 samples (don't need to be matched with 4 or 5)
ID_Decoder_select <- ID_Decoder %>%
  filter(visit_number %in% visits)


# Merge treat to ID_Decoder
treat <- ID_Decoder_select %>% 
  left_join(., treat, by = c("study_id" = "STUDY_ID"))


## Filter treat to include treat samples of interest -----------------------------------------------

treat_filter <- treat %>% 
  # Remove LOP and PLA
  filter(!Treatment %in% c("LOP", "PLA"))



## Select patient ID's that contain pathogens of interest -----------------------------------------
# Bacteria vs Viral
bacteria_list <- c("ETEC|EAEC|EPEC|EHEC|EIEC|C_diff|STEC|Shigella|Salmonella|Vibrio|EcoliO157H7|Aeromonas|Yersinia|Campylobacter")


treat_filter_bacVSvirus <- treat_filter %>% 
  select(glomics_ID, study_id, ends_with(paste0("_", detection_method))) %>% 
  mutate_at(vars(ends_with(paste0("_", detection_method))), list(~ifelse(. == "yes", 1, 0))) %>% 
  # Remove parasites
  select_at(vars(-matches("Giardia|E_histo|Crypto|Cyclospora"))) %>% 
  # Identify bacteria versus virus patients
  mutate(bacteria = rowSums(select_at(., vars(matches(bacteria_list))))) %>% 
  mutate(virus = rowSums(select_at(., vars(contains("virus"))))) %>% 
  # Select bactria only and virus only
  filter(bacteria > 0 & virus == 0 | bacteria == 0 & virus > 0) %>% 
  mutate(bact_or_virus = ifelse(bacteria > 0, "Bacteria", "Virus")) %>% 
  select(glomics_ID, study_id, bact_or_virus)






## Filter the humichip data to include only samples in treat_filter --------------------------------
humichip_filtered <- humichip %>% 
  select_at(c("Genbank.ID", "Gene", "Organism", "Lineage",
              "Gene_category", "Subcategory1", "Subcategory2", treat_filter_bacVSvirus$glomics_ID))


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
  left_join(., treat_filter_bacVSvirus, by = "glomics_ID") %>% 
  ungroup()


## Calculate Response Ratio ------------------------------------------------------------------------

humi_RR <- humi_grouped %>% 
  # Group by Subcategory1
  group_by(Gene_category, Subcategory1, bact_or_virus) %>% 
  summarise(mean_signal = mean(category_abundance),
            sd_signal = sd(category_abundance),
            n = sum(!is.na(category_abundance))) %>% 
  
  # Spread the signal mean by Impact
  ungroup() %>%
  group_by(Gene_category, Subcategory1, bact_or_virus) %>%
  spread(bact_or_virus, mean_signal) %>% 
  
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
  summarise_all(list(~sum(., na.rm = T))) %>% 
  
  # Must have at least __ observations in each subcategory
  filter(n_group1 >= 3) %>%
  filter(n_group2 >= 3) %>%
  
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
  
  
  # plot labels
  labs(y = "Response Ratio",
       x = "") +
  scale_x_discrete(position = "left") +
  
  theme_linedraw() +
  coord_flip() +
  facet_wrap(~pretty_genecat, ncol = 1, scales = "free_y") +
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

ggsave(plot = RR_plot, filename = "results/figures/Manuscript_Figures/RR_V1_subcat1_bactVSvirus.png",
       height = 6,
       width = 5)

