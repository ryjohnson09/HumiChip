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
  
  
  
  
  ## Filter Patients from Data -----------------------------------------
  
  # Filter ID_decoder by visit number (matched or non-matched)
  ID_visit <- reactive({
    
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
  
  
  # Filter Treat based on treatment group, then modify ID_visit
  ID_visit_treatment <- reactive({
    treat_studyIDs <- treat %>%
      # Remove LOP and PLA samples
      filter(!Treatment %in% c("LOP", "PLA")) %>%
      # Samples in select treatment groups
      filter(Treatment %in% input$treatment_groups) %>%
      pull(STUDY_ID)
    
    # Filter ID_decoder
    ID_visit() %>%
      filter(study_id %in% treat_studyIDs)
  })
  
  
  
  # Filter from Humichip Data
  humi_filtered <- eventReactive(input$action, {
    humichip %>%
      select_if(colnames(.) %in% c("Genbank.ID", "gene", "species", "lineage",
                                   "annotation", "geneCategory", "subcategory1",
                                   "subcategory2", ID_visit_treatment()$glomics_ID))
  })
  
  
  
  
  # Show Humi Data Table
  output$humi_table <- renderTable({head(humi_filtered())})

})
