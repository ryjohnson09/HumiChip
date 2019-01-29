###################################
#### Humichip Shiny App Server ####
###################################

## Load Libraries -----------------------------------
library(shiny)
library(tidyverse)
library(vegan)
library(ggExtra)


# Define server -------------------------------------
shinyServer(function(input, output){
  
  
  
  ## Read in Raw Data with Progress Bar --------------------------------
  withProgress(message = "Reading in Data:", {
               
               # Humichip
               incProgress(amount = 1/3, detail = "Reading Humichip Data")
               humichip <- suppressWarnings(suppressMessages(read_tsv("Merged_humichip.tsv")))
               
               
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
                funs(ifelse(. == "yes", gsub(x = deparse(substitute(.)), pattern = "_.*$", replacement = ""), ""))) %>%
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
  
  ## Probe Filtering ---------------------------------------------------------
  humi_probes_filtered <- reactive({
    
    if(input$probe_type == "Functional"){
      # Ensure that at least one functional group is selected
      validate(need(input$geneCategory, 'Please select at least one Functional Group'))
      
      humichip %>%
        filter(geneCategory %in% input$geneCategory)
    } else if (input$probe_type == "Strain/Species"){
      humichip %>%
        filter(gene == "STR_SPE")
    } else if (input$probe_type == "All"){
      humichip
    }
  })
  
  
  
  ## Filter patients from Humichip Data --------------------------------------
  humi_probes_patient_filtered <- eventReactive(input$action, {
    humi_probes_filtered() %>%
      select_if(colnames(.) %in% c("Genbank.ID", "gene", "species", "lineage",
                                   "annotation", "geneCategory", "subcategory1",
                                   "subcategory2", ID_v_c_t_d_p()$glomics_ID))
  })
  
  
  ## Prepare Data for Ordination Analysis ----------------------------------------
  humi_matrix <- eventReactive(input$action, {
    
    validate(need(ncol(humi_probes_patient_filtered()) > 8, "No Patients meet criteria"))
   
    # Set NA's to 0 and values not NA to original value
    humi1 <- humi_probes_patient_filtered() %>%
      select_if(colnames(.) %in% ID_v_c_t_d_p()$glomics_ID) %>%
      mutate_all(funs(ifelse(is.na(.), 0, .)))
    
    # Remove rows that equal 0
    humi1 <- humi1[rowSums(humi1) != 0,]
    
    # Return matrix
    as.matrix(humi1)
  })
  
  
  ## Ordination Analysis ----------------------------------------------------------
  humi_ordination_results <- eventReactive(input$action, {
    withProgress(message = "Performing Ordination:", {
    
    ###########
    ### PCA ###
    ###########
    incProgress(amount = 1/2)
    humi_PCA <- vegan::rda(t(humi_matrix()))
    })
  })
  
  # Extract coordinates as tibble
  humi_PCA_coordinates <- eventReactive(input$action, {
    as.data.frame(scores(humi_ordination_results(), display = "sites")) %>%
      rownames_to_column(var = "glomics_ID")
  })
    
  # Get Proportion explained
  PCA_prop_expl <- eventReactive(input$action, {
    summary(eigenvals(humi_ordination_results()))[2,] * 100
  })

  
  
  
  ## Merge Ordination Analysis with Metadata --------------------------------------
  humi_ordination_metadata <- eventReactive(input$action, {
    humi_PCA_coordinates() %>%
      # Add study ID's
      full_join(., ID_v_c_t_d_p(), by = "glomics_ID") %>%
      # Add in pathogen list from treat_pathogen()
      left_join(., treat_pathogens(), by = c("study_id" = "STUDY_ID")) %>%
      # Add in remaining treat metadata
      left_join(., treat, by = c("study_id" = "STUDY_ID")) %>%
      
      # Factor Columns
      mutate(visit_number = factor(visit_number)) %>%
      mutate(Impact_of_illness_on_activity_level = factor(Impact_of_illness_on_activity_level))
  })
    
  
  
  
    
  ## Plot --------------------------------------------------------------
  
  ################
  ### PCA PLOT ###
  ################
  plotInput <- reactive({
    
    # Color of points
    my_fill <- ifelse(input$point_color == "None", "NULL", input$point_color)
    
    # Aesthetic sizes
    axis_title_size <- 18
    axis_text_size <- 16
    title_size <- 20
    legend_text_size <- 13
    point_size <- 5
    
    # PCA plot
    pca_plot <- ggplot(humi_ordination_metadata(),
                       aes_string(x = "PC1", 
                                  y = "PC2", 
                                  color = my_fill)) +
      
      # Set up proportion explained
      xlab(paste0("PC1(", round(PCA_prop_expl()[[1]], 2), "%)")) +
      ylab(paste0("PC2(", round(PCA_prop_expl()[[2]], 2), "%)")) +
      
      
      geom_point(pch = 1, alpha = 1, size = point_size) +
      geom_point(pch = 19, alpha = 0.8, size = point_size) +
      ggtitle("PCA Analysis") +
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
    
    ggMarginal(pca_plot, groupColour = TRUE, groupFill = TRUE)
  })
  
  
  ####################
  ### Display Plot ###
  ####################
  output$humi_plot <- renderPlot({
    print(plotInput())
  })
  
  #####################
  ### Download plot ###
  #####################
  
  output$downloadPlot <- downloadHandler(
    filename = function(){paste("shiny_plot",'.png',sep='')},
    content = function(file){
      ggsave(file, plot=plotInput())
    })
    
  
  
  
  
  
  
  
  # Show Humi Data Table
  output$humi_table <- renderTable({humi_ordination_metadata()})
  
  output$random_text <- renderText({ 
    paste(ncol(humi_probes_patient_filtered()), " columns")
  })

})
