###################################
#### Humichip Shiny App Server ####
###################################

## Load Libraries -----------------------------------
library(shiny)
library(tidyverse)
library(vegan)
library(ggExtra)
library(ape)


# Define server -------------------------------------
shinyServer(function(input, output){
  
  
  
  ## Read in Raw Data with Progress Bar --------------------------------
  withProgress(message = "Reading in Data:", {
               
               # Humichip
               incProgress(amount = 1/3, detail = "Reading Humichip Data")
               humichip <- suppressWarnings(suppressMessages(read_tsv("Merged_humichip_Renormalized.tsv")))
               
               
               # Treat
               incProgress(amount = 1/3, detail = "Reading Metadata Data")
               treat <- suppressWarnings(suppressMessages(read_csv("TrEAT_Clinical_Metadata_tidy.csv")))
               
               
               # ID Decoder
               incProgress(amount = 1/3, detail = "Reading Decoder")
               ID_decoder <- suppressWarnings(suppressMessages(read_csv("ID_Decoder_Humichip.csv")))
      })
  
  
  
  
  ## Filter Patients from Data -------------------------------------------------------
  
  ##############################
  ### Visit Number / Matched ###
  ##############################
  ID_v <- reactive({
    
    # Ensure that at least visit is selected
    validate(need(input$visit_number, 'Please select at least one visit number'))
    
    # If non-matched
    if (!input$matched){
      filter(ID_decoder, visit_number %in% as.numeric(input$visit_number))
    
    # if matched 
    } else if (input$matched){
      matched_samples <- ID_decoder %>%
        filter(visit_number %in% as.numeric(input$visit_number)) %>%
        count(study_id) %>%
        filter(n == length(as.numeric(input$visit_number))) %>%
        pull(study_id)
      
      # Get glomics ID's that correspond to the patients
      # with matched samples
      matched_glomics <- ID_decoder %>%
        filter(visit_number %in% as.numeric(input$visit_number)) %>%
        filter(study_id %in% matched_samples) %>%
        pull(glomics_ID)
      
      # Filter ID Decoder for only matched
      filter(ID_decoder, glomics_ID %in% matched_glomics)
    }
  })
  
  
  #########################
  ### Filter by Country ###
  #########################
  ID_v_c <- reactive({
    
    # Ensure that at least country is selected
    validate(need(input$country, 'Please select at least one country'))
    
    # Select patients by country from treat
    country_patients <- treat %>%
      filter(country %in% input$country) %>%
      pull(STUDY_ID)
    
    # Filter ID_Decoder
    ID_v() %>%
      filter(study_id %in% country_patients)
  })
  
  
  ##########################################
  ### Treatment Group / Remove LOP & PLA ###
  ##########################################
  ID_v_c_t <- reactive({
    
    # Ensure that at least treatment group is selected
    validate(need(input$treatment_groups, 'Please select at least one treatment group'))
    
    treat_studyIDs <- treat %>%
      # Remove LOP and PLA samples
      filter(!Treatment %in% c("LOP", "PLA")) %>%
      # Samples in select treatment groups
      filter(Treatment %in% input$treatment_groups) %>%
      pull(STUDY_ID)
    
    # Filter ID_decoder
    ID_v_c() %>%
      filter(study_id %in% treat_studyIDs)
  })
  
  
  
  ########################
  ### Disease Severity ###
  ########################
  ID_v_c_t_d <- reactive({
    
    # Ensure that at least disease severity is selected
    validate(need(input$disease_severity, 'Please select at least one disease severity group'))
    
    treat_disease <- treat %>%
      filter(LLS_severity %in% input$disease_severity) %>%
      pull(STUDY_ID)
    
    # Filter ID_decoder
    ID_v_c_t() %>%
      filter(study_id %in% treat_disease)
  })
  
  
  
  ##########################
  ### Pathogen Detection ###
  ##########################
  treat_pathogens <- reactive({
    treat_select <- treat %>%
      # Select columns with appropriate detection method
      select(STUDY_ID, ends_with(input$detection_method)) %>%
      # Determine if co-infection or not
      mutate(num_pathogens = rowSums(. == "yes", na.rm = TRUE))
    
    # Allow Coinfections?
    if (!input$allow_coinfections){
      treat_select_coinfections <- treat_select %>%
        filter(num_pathogens == 1)
    } else if (input$allow_coinfections){
      treat_select_coinfections <- treat_select
    }
    
    treat_new <- treat_select_coinfections %>%
      # Select for pathogens of interest
      select(STUDY_ID, paste(input$pathogens, input$detection_method, sep = "_")) %>%
      # Replace yes no with Pathogen
      mutate_at(vars(matches(input$detection_method)),
                list(~ifelse(. == "yes", gsub(x = deparse(substitute(.)), pattern = "_.*$", replacement = ""), ""))) %>%
      # If sample has an NA, remove it (wasn't tested)
      mutate(has_na = rowSums(is.na(.))) %>%
      filter(has_na == 0) %>%
      select(-has_na) %>%
      # If no pathogens of intereste detected, remove it
      mutate(path_present = rowSums(. == "")) %>%
      filter(path_present != length(input$pathogens)) %>%
      select(-path_present) %>%
      # Merge pathogens into one column
      unite(pathogens, paste(input$pathogens, input$detection_method, sep = "_"), sep = "_") %>%
      mutate(pathogens = str_replace_all(pathogens, "^_*", "")) %>%
      mutate(pathogens = str_replace_all(pathogens, "_*$", "")) %>%
      mutate(pathogens = str_replace_all(pathogens, "_+", "-"))
    
    treat_new
  })
  
  # Filter ID if selecting pathogens
  ID_v_c_t_d_p <- reactive({
    if (input$pathogen_select){
      
      # Ensure that at least one pathogen is selected
      validate(need(input$pathogens, 'Please select at least one Pathogen'))
      
      ID_v_c_t_d() %>%
        filter(study_id %in% treat_pathogens()$STUDY_ID)
    } else if (!input$pathogen_select){
      ID_v_c_t_d()
    }
  })
  
  
  #####################################
  ### Filter for patients with ESBL ###
  #####################################
  treat_esbls <- reactive({
    
    # Create tibble of just ESBL columns and modify
    treat_esbl_select <- treat %>% 
      # Select columns of interest
      select(STUDY_ID, ESBL_V1, ESBL_V5, ends_with("either_V1"), ends_with("either_V5")) %>%
      # Change ESBL columns to Yes/No
      mutate_at(vars(starts_with("ESBL")), list(~ifelse(. == "Positive", "Yes",
                                                 ifelse(. == "Negative", "No", NA)))) %>% 
      # Replace yes no with ESBL
      mutate_at(vars(-STUDY_ID),
                list(~ifelse(. == "Yes", deparse(substitute(.)), ""))) %>% 
      # If no ESBLs of interest detected, remove it
      mutate(esbl_present = rowSums(. != "", na.rm = TRUE)) %>%
      filter(esbl_present > 1) %>%
      select(-esbl_present) %>% 
      # Merge pathogens into one column
      unite(all_esbls, -STUDY_ID, sep = ":") %>%
      mutate(all_esbls = str_replace_all(all_esbls, "NA", "")) %>% 
      mutate(all_esbls = str_replace_all(all_esbls, "^:*", "")) %>%
      mutate(all_esbls = str_replace_all(all_esbls, ":*$", "")) %>%
      mutate(all_esbls = str_replace_all(all_esbls, ":+", "; "))
    
    # Select for ESBL of interest
    treat_esbl <- treat_esbl_select %>%
      rowwise() %>% 
      # change if you want to select if contain ALL ESBL selected
      # or select if contain ANY ESBL selected
      #filter(all(str_detect(all_esbls, input$esbls))) %>% # all
      filter(str_detect(all_esbls, paste(input$esbls, collapse = "|"))) %>% 
      # remove either from column
      mutate(all_esbls = str_replace_all(string = all_esbls, pattern = "_either", replacement = ""))
    
    treat_esbl
  })
  
  # Filter ID if selecting ESBLs
  ID_v_c_t_d_p_e <- reactive({
    if (input$esbl_select){
      
      # Ensure that at least one esbl is selected
      validate(need(input$esbls, 'Please select at least one ESBL'))
      
      ID_v_c_t_d_p() %>%
        filter(study_id %in% treat_esbls()$STUDY_ID)
    } else if (!input$esbl_select){
      ID_v_c_t_d_p()
    }
  })
    
  
  
  
  
  ## Probe Filtering ---------------------------------------------------------
  humi_probes_filtered <- reactive({
    
    if(input$probe_type == "Functional"){
      # Ensure that at least one functional group is selected
      validate(need(input$Gene_category, 'Please select at least one Functional Group'))
      
      humichip %>%
        filter(Gene_category %in% input$Gene_category)
    } else if (input$probe_type == "Strain/Species"){
      humichip %>%
        filter(Gene == "STR_SPE")
    } else if (input$probe_type == "All"){
      humichip
    }
  })
  
  
  
  ## Filter patients from Humichip Data --------------------------------------
  humi_probes_patient_filtered <- eventReactive(input$action, {
    humi_probes_filtered() %>%
      select_at(c("Genbank.ID", "Gene", "Organism", "Lineage",
                  "Gene_category", "Subcategory1",
                  "Subcategory2", ID_v_c_t_d_p_e()$glomics_ID))
  })
  
  
  ## Prepare Data for Ordination Analysis ----------------------------------------
  humi_matrix <- eventReactive(input$action, {

    #validate(need(ncol(humi_probes_patient_filtered()) > 7, "No Patients meet criteria"))

    # Set NA's to 0 and values not NA to original value
    humi1 <- humi_probes_patient_filtered() %>%
      select_at(ID_v_c_t_d_p_e()$glomics_ID) %>%
      mutate_all(list(~ifelse(is.na(.), 0, .)))

    # Remove rows that equal 0
    humi1 <- humi1[rowSums(humi1) != 0,]

    # Return matrix
    as.matrix(humi1)
  })


  ## Ordination Analysis ----------------------------------------------------------
  humi_ordination_results <- eventReactive(input$action, {
    withProgress(message = "Performing Ordination: ", {

      # Perform PCA analysis using vegan
      if(input$ordination == "PCA"){
        incProgress(amount = 1/2, detail = "PCA analysis")
        vegan::rda(t(humi_matrix()))
      # Perform PCoA analysis using ape
      } else if (input$ordination == "PCoA"){
        incProgress(amount = 1/3, detail = "Calculating distance matrix")
        humichip_dist <- vegan::vegdist(as.matrix(t(humi_matrix())))
        incProgress(amount = 1/3, detail = "PCoA analysis")
        ape::pcoa(humichip_dist, correction = "none")
      # Perform DCA analysis using vegan
      } else if (input$ordination == "DCA"){
        incProgress(amount = 1/2, detail = "DCA analysis")
        vegan::decorana(t(humi_matrix()))
      } else {
        stopApp("Error in running ordination")
      }
    })
  })



  # Extract coordinates as tibble
  humi_coordinates <- eventReactive(input$action, {

    # PCA or DCA
    if (input$ordination == "PCA" | input$ordination == "DCA"){
      as.data.frame(scores(humi_ordination_results(), display = "sites")) %>%
        rownames_to_column(var = "glomics_ID")
    # If PCoA
    } else if (input$ordination == "PCoA"){
      as.data.frame(humi_ordination_results()$vectors[,1:2]) %>%
        rownames_to_column(var = "glomics_ID")
    } else {
      stopApp("Error in extracting ordination coordinates")
    }
  })



  # Get Proportion explained
  ord_prop_expl <- eventReactive(input$action, {
    if (input$ordination == "PCA"){
      summary(eigenvals(humi_ordination_results()))[2,] * 100
    } else if (input$ordination == "PCoA"){
      humi_ordination_results()$values$Relative_eig * 100
    } else {
      stopApp("Error in extraction percent explaind by axes")
    }
  })




  ## Merge Ordination Analysis with Metadata --------------------------------------
  humi_ordination_metadata <- eventReactive(input$action, {
    humi_coordinates() %>%
      # Add study ID's
      full_join(., ID_v_c_t_d_p_e(), by = "glomics_ID") %>%
      # Add in pathogen list from treat_pathogen()
      left_join(., treat_pathogens(), by = c("study_id" = "STUDY_ID")) %>%
      # Add in ESBL list from treat_esbl()
      left_join(., treat_esbls(), by = c("study_id" = "STUDY_ID")) %>%
      # Add in remaining treat metadata
      left_join(., treat, by = c("study_id" = "STUDY_ID")) %>%

      # Factor Columns
      mutate(visit_number = factor(visit_number)) %>%
      mutate(Impact_of_illness_on_activity_level = factor(Impact_of_illness_on_activity_level))
  })





  ## Plot --------------------------------------------------------------
  # Aesthetic sizes
  axis_title_size <- 18
  axis_text_size <- 16
  title_size <- 20
  legend_text_size <- 13
  point_size <- 4

  # Fill
  my_fill <- reactive({
    ifelse(input$point_color == "None", "NULL", input$point_color)
  })


  # Set up Base Plot
  humi_ord_base_plot <- reactive({

    ggplot(humi_ordination_metadata()) +
      theme_minimal() +
      theme(
        axis.title.x = element_text(size = axis_title_size),
        axis.text.x = element_text(size = axis_text_size, hjust = 1),
        axis.text.y = element_text(size = axis_text_size),
        axis.title.y = element_text(size = axis_title_size),
        plot.title = element_text(size = title_size, face = "bold"),
        legend.text = element_text(size = legend_text_size),
        legend.title = element_blank()) +
      guides(fill = guide_legend(override.aes = list(size=7)))
  })

  # eventReactive to chose plot type
  plot_type <- eventReactive(input$action, {
    ifelse(input$ordination == "PCA", "PCA",
           ifelse(input$ordination == "DCA", "DCA", "PCoA"))
  })



  plotInput <- reactive({


    ################
    ### PCA PLOT ###
    ################
    if (plot_type() == "PCA"){
      humi_ord_plot <- humi_ord_base_plot() +
        xlab(paste0("PC1(", round(ord_prop_expl()[[1]], 2), "%)")) +
        ylab(paste0("PC2(", round(ord_prop_expl()[[2]], 2), "%)")) +
        geom_point(aes_string(x = "PC1", y = "PC2", color = my_fill()),
                   pch = 1, alpha = 1, size = point_size) +
        geom_point(aes_string(x = "PC1", y = "PC2", color = my_fill()),
                   pch = 19, alpha = 0.8, size = point_size) +
        ggtitle("PCA Analysis")

      # Connect lines?
      if (input$lines == "Yes"){
        humi_ord_plot <- humi_ord_plot +
          geom_line(aes_string(x = "PC1", y = "PC2", group = "study_id"), linetype = 1, color = "black", size = 0.8)
      } else {
        humi_ord_plot <- humi_ord_plot
      }
      # Plot
      ggMarginal(humi_ord_plot, groupColour = TRUE, groupFill = TRUE)

    #################
    ### PCoA PLOT ###
    #################
    } else if (plot_type() == "PCoA"){
      humi_ord_plot <- humi_ord_base_plot() +
        xlab(paste0("PCoA1(", round(ord_prop_expl()[[1]], 2), "%)")) +
        ylab(paste0("PCoA2(", round(ord_prop_expl()[[2]], 2), "%)")) +
        geom_point(aes_string(x = "Axis.1", y = "Axis.2", color = my_fill()),
                   pch = 1, alpha = 1, size = point_size) +
        geom_point(aes_string(x = "Axis.1", y = "Axis.2", color = my_fill()),
                   pch = 19, alpha = 0.8, size = point_size) +
        ggtitle("PCoA Analysis")

      # Connect lines?
      if (input$lines == "Yes"){
        humi_ord_plot <- humi_ord_plot +
          geom_line(aes_string(x = "Axis.1", y = "Axis.2", group = "study_id"), linetype = 1, color = "black", size = 0.8)
      } else {
        humi_ord_plot <- humi_ord_plot
      }
      # Plot
      ggMarginal(humi_ord_plot, groupColour = TRUE, groupFill = TRUE)

    ################
    ### DCA PLOT ###
    ################
    } else if (plot_type() == "DCA"){
      humi_ord_plot <- humi_ord_base_plot() +
        xlab("DCA1") +
        ylab("DCA2") +
        geom_point(aes_string(x = "DCA1", y = "DCA2", color = my_fill()),
                   pch = 1, alpha = 1, size = point_size) +
        geom_point(aes_string(x = "DCA1", y = "DCA2", color = my_fill()),
                   pch = 19, alpha = 0.8, size = point_size) +
        ggtitle("DCA Analysis")

      # Connect lines?
      if (input$lines == "Yes"){
        humi_ord_plot <- humi_ord_plot +
          geom_line(aes_string(x = "DCA1", y = "DCA2", group = "study_id"), linetype = 1, color = "black", size = 0.8)
      } else {
        humi_ord_plot <- humi_ord_plot
      }
      # Plot
      ggMarginal(humi_ord_plot, groupColour = TRUE, groupFill = TRUE)

    } else {stopApp("Error generating plot")}
  })



  ## Display Plot ---------------------------------------------------
  output$humi_plot <- renderPlot({
    print(plotInput())
  })

  ## Download plot --------------------------------------------------
  output$downloadPlot <- downloadHandler(
    filename = function(){paste("shiny_plot",'.png',sep='')},
    content = function(file){
      ggsave(file, plot=plotInput(), width = 10, height = 9)
    })




  ## Perform Statistics ----------------------------------------------

  # Remove samples for which the groupings are NA values
  new_ord_meta <- reactive({
    humi_ordination_metadata() %>%
      filter(!is.na(!!sym(input$point_color)))
  })

  new_matrix <- reactive({
    humi_matrix()[ ,colnames(humi_matrix()) %in% new_ord_meta()$glomics_ID]
  })


  # Display the groupings
  stats_groups <- reactive({
    if (input$stat_calc){
      paste("Groupings:", paste(unique(new_ord_meta()[[input$point_color]]), collapse = ", "))
    }
  })

  adonis_results <- reactive({
    # Ensure that new_ord_meta()$glomics_ID is in the
    #  same order as colnames(new_matrix())
    if (!all(new_ord_meta()$glomics_ID == colnames(new_matrix()))){
      stopApp("Error calculating adonis P-value")
    }

    if (input$stat_calc){

      # Calculate adonis results
      withProgress(message = "Performing adonis test: ", value = 0.33, {
        adonis_temp <-  vegan::adonis(t(new_matrix()) ~ new_ord_meta()[[input$point_color]],
                                      method = "bray",
                                      perm = 99)
      })

      paste("adonis p-value: ", adonis_temp$aov.tab$`Pr(>F)`[1])
    }
  })


  mrpp_results <- reactive({
    # Ensure that new_ord_meta()$glomics_ID is in the
    #  same order as colnames(new_matrix())
    if (!all(new_ord_meta()$glomics_ID == colnames(new_matrix()))){
      stopApp("Error calculating mrpp P-value")
    }

    if (input$stat_calc){

      # Calculate mrpp results
      withProgress(message = "Performing mrpp test: ", value = 0.66, {
        mrpp_temp <-  vegan::mrpp(t(new_matrix()),
                                  new_ord_meta()[[input$point_color]])
      })

      paste("mrpp p-value: ", mrpp_temp$Pvalue)
    }
  })


  ## Print stats --------------------------------------------------

  output$stats_groupings <- renderText({
    stats_groups()
  })

  output$adonis_pvalue <- renderText({
    adonis_results()
  })

  output$mrpp_pvalue <- renderText({
    mrpp_results()
  })
  
    

  ## Show Table ----------------------------------------------------
  output$humi_table <- renderTable({humi_ordination_metadata()})
  
})
