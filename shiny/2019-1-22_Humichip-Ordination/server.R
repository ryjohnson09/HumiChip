###################################
#### Humichip Shiny App Server ####
###################################

## Load Libraries -----------------------------------
library(shiny)
library(tidyverse)
library(vegan)

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
  
  
  ##########################################
  ### Treatment Group / Remove LOP & PLA ###
  ##########################################
  ID_v_t <- reactive({
    treat_studyIDs <- treat %>%
      # Remove LOP and PLA samples
      filter(!Treatment %in% c("LOP", "PLA")) %>%
      # Samples in select treatment groups
      filter(Treatment %in% input$treatment_groups) %>%
      pull(STUDY_ID)
    
    # Filter ID_decoder
    ID_v() %>%
      filter(study_id %in% treat_studyIDs)
  })
  
  
  
  ########################
  ### Disease Severity ###
  ########################
  ID_v_t_d <- reactive({
    treat_disease <- treat %>%
      filter(LLS_severity %in% input$disease_severity) %>%
      pull(STUDY_ID)
    
    # Filter ID_decoder
    ID_v_t() %>%
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
      mutate(pathogens = str_replace(pathogens, "^_*", "")) %>%
      mutate(pathogens = str_replace(pathogens, "_*$", "")) %>%
      mutate(pathogens = str_replace(pathogens, "_+", "-"))
    
    treat_new
  })
  
  # Filter ID if selecting pathogens
  ID_v_t_d_p <- reactive({
    if (input$pathogen_select){
    ID_v_t_d() %>%
      filter(study_id %in% treat_pathogens()$STUDY_ID)
    } else if (!input$pathogen_select){
      ID_v_t_d()
    }
  })
  
  ## Probe Filtering ---------------------------------------------------------
  humi_probes_filtered <- reactive({
    
    if(input$probe_type == "Functional"){
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
                                   "subcategory2", ID_v_t_d_p()$glomics_ID))
  })
  
  
  
  
  
  # Show Humi Data Table
  output$humi_table <- renderTable({head(humi_probes_patient_filtered())})

})
