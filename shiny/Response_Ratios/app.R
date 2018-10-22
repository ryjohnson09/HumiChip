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

visit_choices <- c("Visit 1 vs 4",
                   "Visit 1 vs 5",
                   "Acute vs Conv")

phylum_choices <- humichip %>%
  filter(str_detect(lineage, "Bacteria")) %>%
  filter(str_detect(lineage, ";phylum")) %>%
  mutate(Phylum = gsub(x = lineage, # Phylum column
                       pattern = ".*;phylum:(\\w*\\s*[-]*\\w*);.*", 
                       replacement = "\\1")) %>%
  select(Phylum) %>%
  distinct() %>%
  pull(Phylum)

cat_choices <- list("Gene Category" = "geneCategory",
                    "Subcategory 1" = "subcategory1",
                    "Subcategory 2" = "subcategory2",
                    "Phylum" = "Phylum")

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
               radioButtons('Visit_Number', 'Visit:', choices = visit_choices, selected = "Visit 1 vs 5", inline = TRUE),
               helpText("Select patient groups to compare by visit"),
               
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
               selectInput("cat_choice", "Y-axis Categories", choices = cat_choices, selected = "geneCategory"),
               helpText("Select Categories to Compare")))),
    
    
    fluidRow(
      h3("Ordination and Aesthetics"),
      
      column(12, 
             wellPanel(
               )),
    
      
      downloadButton('downloadPlot','Download Plot')),
    
    
    #sidebar width
    width = 4),
  
  # Plot
  mainPanel(
    #plotOutput("plot", width = "800px", height = "800px"),
    
    # Table to see patients (not needed, but useful for troubleshooting)
    fluidRow(column(12,tableOutput('table')))
    
    ################################
    ### Notes Regarding Analysis ###
    ################################
    #helpText("")
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
      matched_samples <- 
        
        #Visit 1 and 4
        if (input$Visit_Number == "Visit 1 vs 4"){
          filter(ID_decoder, visit_number %in% c(1,4))
          
        # Visit 1 and 5
        } else if (input$Visit_Number == "Visit 1 vs 5"){
          filter(ID_decoder, visit_number %in% c(1,5))
          
        # Acute vs Conv
        } else if (input$Visit_Number == "Acute vs Conv"){
          filter(ID_decoder, visit_number %in% c(1,4,5))
            
        } else {
          stopApp("Error with matched samples")
        }
    
      matched_samples <- matched_samples %>% 
        count(study_id) %>%
        filter(n == length(unique(matched_samples$visit_number))) %>%
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
    } else {
      humichip
    }
  })
  
  
  ###############################################
  ### Select for Functional or STR_SPE probes ###
  ###############################################
  
  humichip_probe <- reactive({
    
    # If cat_choice == Phylum, filter for STR_SPE
    # and add Phylum category
    if(input$cat_choice == "Phylum"){
      humichip_select() %>%
        filter(gene == "STR_SPE") %>%
        filter(str_detect(lineage, "Bacteria")) %>%
        filter(str_detect(lineage, ";phylum")) %>%
        mutate(Phylum = gsub(x = lineage, # Phylum column
                             pattern = ".*;phylum:(\\w*\\s*[-]*\\w*);.*", 
                             replacement = "\\1"))
      
    # If cat_choice == Functional group
    } else if (input$cat_choice != "Phylum"){
      humichip_select() %>%
        filter(gene != "STR_SPE")
    } else {
      stopApp("Problem filtering STR_SPE vs Functional Probes")
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
  
  
  #######################
  ### Treatment Group ###
  #######################
  
  # Filter for treatment Group
  treat_tx <- reactive({
    if(input$tx_group == "AZI"){
      filter(treat, Treatment == "AZI")
      
    } else if (input$tx_group == "RIF"){
      filter(treat, Treatment == "RIF")
      
    } else if (input$tx_group == "LEV"){
      filter(treat, Treatment == "LEV")
      
    } else {
      return(treat)
    }
  })
  
  
  #################
  ##### Visit #####
  #################
  
  # Filter ID_decoder based on Visit Number
  treat_visit_filtered <- reactive({
    if(input$Visit_Number == "Visit 1 vs 4"){
      filter(treat_tx(), visit_number %in% c(1,4))
    } else if (input$Visit_Number == "Visit 1 vs 5"){
      filter(treat_tx(), visit_number %in% c(1,5))
    } else if (input$Visit_Number == "Acute vs Conv"){
      filter(treat_tx(), visit_number %in% c(1,4,5))
    } else {
      stopApp("Problem Filtering by Visit")
    }
  })
  

  #####################################
  #### TrEAT DB filter by pathogen ####
  #####################################
  
  # Filter for pathogen of interest
  treat_pathogen_filtered <- reactive({
    
    # Include all patients
    if (input$pathogens == "All"){
      treat_visit_filtered()
    } else {
      
      # Select for patients with specified pathogen
      pathogen <- sym(input$pathogens)
      treat_visit_filtered() %>%
        filter(!!pathogen == "yes")
    }
  })
  
  
  ############################################
  ### Filter Humichip for certain patients ###
  ############################################
  
  
  humichip_final <- reactive({
    
    humichip_probe()
      
  
  })
  
  output$table <- renderTable({head(treat_visit_filtered(), 25)})
  
}







shinyApp(server = server, ui = ui)