# Humichip Response Ratio App

library(shiny)
library(tidyverse)

#######################
#### Read in Files ####
#######################

ID_decoder <- suppressWarnings(suppressMessages(read_csv("ID_Decoder.csv")))
treat <- suppressWarnings(suppressMessages(read_csv("TrEAT_Clinical_Metadata_tidy.csv")))
humichip <- suppressWarnings(suppressMessages(read_tsv("Merged_humichip.tsv")))


######################
### Set up choices ###
######################

matched_choices <- list("All Samples" = "all_samples",
                        "Matched Samples" = "matched_samples")

tx_choices <- c("RIF", "LEV", "AZI", "All")

visit_choices <- list()

phylum_choices <- humichip %>%
  filter(str_detect(lineage, "Bacteria")) %>%
  filter(str_detect(lineage, ";phylum")) %>%
  mutate(Phylum = gsub(x = lineage, # Phylum column
                       pattern = ".*;phylum:(\\w*\\s*[-]*\\w*);.*", 
                       replacement = "\\1")) %>%
  select(Phylum) %>%
  distinct() %>%
  pull(Phylum)

gene_cat_choices <- c()

detection_choices <- c("Culture", "Taq", "Both", "Either")

culture_choices <- c("All", colnames(treat)[grep(pattern = "_culture$", colnames(treat))])

taq_choices <- c("All", colnames(treat)[grep(pattern = "_taq$", colnames(treat))])

both_choices <- c("All", colnames(treat)[grep(pattern = "_both$", colnames(treat))])

either_choices <- c("All", colnames(treat)[grep(pattern = "_either$", colnames(treat))])





## UI ------------------------------------------------------------------------------------

ui <- fluidPage(
  
  titlePanel("Humichip Ordination"),
  sidebarLayout(
    sidebarPanel(
      
      # Matched or All samples
      radioButtons("matched", label = "Matched or All samples",
                   choices = matched_choices, inline = TRUE, selected = "matched_samples"),
      helpText("Matched = Only patients that provided samples for all selected visits"),
      
    ##########################################
    ### Inputs that Select Certain Samples ###
    ##########################################
    fluidRow(
      h3("Sample Selection:"),
      
      column(12, 
             wellPanel(
               # Treatment Choices
               selectInput("tx_group", "Treatment Groups", choices = tx_choices, selected = "All"),
               helpText("Select patients samples in specified treatment groups"),
               
               # Visit
               checkboxGroupInput('Visit_Number', 'Visit:', choices = visit_choices, selected = c(1, 4, 5), inline = TRUE),
               helpText("Select patient samples from specified visit number"),
               
               # Pathogen Detection
               selectInput("path_detection", "Pathogen Detection Method", choices = detection_choices, selected = "Both"),
               helpText("How are samples determined to be positive for pathogen"),
               
               # Pathogen Choices
               uiOutput("secondSelection"),
               helpText("If", code("All"), ", then all samples included.")))),
    
    ########################
    ### Y axis categories ##
    ########################
    fluidRow(
      h3("Functional Group Selection"),
      
      column(12, 
             wellPanel(
               # gene category vs subcateogry1 vs subcategory2
               selectInput("gene_cat", "Gene Categories", choices = gene_cat_choices, selected = "All"),
               helpText("Select patients samples in specified treatment groups"),
               ))),
    
    
    fluidRow(
      h3("Ordination and Aesthetics"),
      
      column(12, 
             wellPanel(
               ))),
    
      
      downloadButton('downloadPlot','Download Plot')),
    
    
    #sidebar width
    width = 4),
  
  # Plot
  mainPanel(
    plotOutput("plot", width = "800px", height = "800px"),
    
    # Table to see patients (not needed, but useful for troubleshooting)
    #fluidRow(column(12,tableOutput('table'))),
    
    ################################
    ### Notes Regarding Analysis ###
    ################################
    helpText("")
  )))



## Server --------------------------------------------------------------------------------

server <- function(input, output){
  
  #####################################################
  ### Add glomics ID's and visit number to treat DB ###
  #####################################################
  
  # Filter ID_decoder to only include samples in humichip
  ID_decoder <- ID_decoder %>%
    filter(glomics_ID %in% colnames(humichip))
  
  treat <- treat %>%
    left_join(., ID_decoder, by = c("STUDY_ID" = "study_id"))
  
  #########################################
  ### Remove samples from humichip that ###
  ### are in the LOP & PLA tx groups ######
  #########################################
  
  LOP_PLA_samples <- treat %>%
    filter(Treatment %in% c("LOP", "PLA")) %>%
    pull(glomics_ID)
  
  humichip <- humichip %>%
    select(-one_of(LOP_PLA_samples))
  
  rm(LOP_PLA_samples)
  
  ##############################################
  ### Only keep patients with matched visits ###
  ##############################################
  
  humichip_select <- reactive({
    if(input$matched == "matched_samples"){
      # Get vector of patients that have matched visits
      matched_samples <- ID_decoder %>%
        filter(visit_number %in% as.numeric(input$Visit_Number)) %>%
        count(study_id) %>%
        filter(n == length(as.numeric(input$Visit_Number))) %>%
        pull(study_id)
      
      # Get glomics ID's that correspond to the patients
      # with matched samples
      matched_glomics <- ID_decoder %>%
        filter(study_id %in% matched_samples) %>%
        pull(glomics_ID)
      
      # Subset the humichip data
      humichip %>%
        select_if(colnames(.) %in% c("Genbank.ID", "gene", "species", "lineage",
                                     "annotation", "geneCategory", "subcategory1",
                                     "subcategory2", matched_glomics))
      
      #rm(matched_glomics, matched_samples)
      
    } else {
      humichip
    }
  })
  
  
  
  #######################################
  ### Render pathogen detection list ###
  #######################################
  output$secondSelection <- renderUI({
    if(input$path_detection == "Culture"){
      selectInput("pathogens", "Pathogens", choices = culture_choices)
    } else if (input$path_detection == "Taq"){
      selectInput("pathogens", "Pathogens", choices = taq_choices)
    } else if (input$path_detection == "Either"){
      selectInput("pathogens", "Pathogens", choices = either_choices)
    } else {
      selectInput("pathogens", "Pathogens", choices = both_choices)
    }
  })
  
  
  
}







shinyApp(server = server, ui = ui)