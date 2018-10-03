# SPECIES HEATMAP

library(shiny)
library(tidyverse)
library(readxl)

#######################
#### Read in Files ####
#######################

ID_decoder <- suppressWarnings(suppressMessages(read_csv("ID_Decoder.csv")))
treat <- suppressWarnings(suppressMessages(read_csv("TrEAT_Clinical_Metadata_tidy.csv")))
humichip <- suppressWarnings(suppressMessages(read_tsv("Merged_humichip.tsv")))
total_probes <- suppressWarnings(suppressMessages(read_xlsx("Species_probe_numbers.xlsx")))


######################
### Set up choices ###
######################

matched_choices <- list("All Samples" = "all_samples",
                        "Matched Samples" = "matched_samples")

visit_choices <- list("Visit 1" = 1,
                      "Visit 4" = 4,
                      "Visit 5" = 5)

tx_choices <- c("RIF", "LEV", "AZI", "All")

detection_choices <- c("Culture", "Taq", "Both", "Either")

culture_choices <- c("All", colnames(treat)[grep(pattern = "_culture$", colnames(treat))])

taq_choices <- c("All", colnames(treat)[grep(pattern = "_taq$", colnames(treat))])

both_choices <- c("All", colnames(treat)[grep(pattern = "_both$", colnames(treat))])

either_choices <- c("All", colnames(treat)[grep(pattern = "_either$", colnames(treat))])

facet_choices <- list(
  "Visit" = "visit_number",
  "Treatment" = "Treatment",
  "Impact on Activity" = "Impact_of_illness_on_activity_level",
  "Past 8 hours diarrhea severity" = "LLS_severity",
  "Alternative Categories" = "Alternative_Categories",
  "ESBL Visit 1" = "ESBL_V1",
  "ESBL Visit 5" = "ESBL_V5")



## UI ------------------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Species Heatmap"),
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
      fluidRow(
        h3("Visualization Parameters"),
        
        column(12, 
               wellPanel(
                 # Ordination Type
                 radioButtons("facet", label = h3("Facet by:"),
                              choices = facet_choices, selected = "visit_number"),
                 
                 # Probe per species
                 sliderInput("probe_number", "Minimum Probe # per species:",
                             min = 0, max = 17,
                             value = 10),
                 
                 # Probe Threshold
                 sliderInput("probe_threshold", "% probes needed to be consider species present:",
                             min = 0, max = 100,
                             value = 60),
                 
                 # Species in samples
                 sliderInput("species_threshold", "# samples positive to include species in heatmap:",
                             min = 0, max = 10,
                             value = 3))),
        
        
      
      downloadButton('downloadPlot','Download Plot')),
    
    
    #sidebar width
    width = 4),
  
  # Plot
  mainPanel(
    plotOutput("plot", width = "1000px", height = "800px"),
    
    # Table to see patients (not needed, but useful for troubleshooting)
    fluidRow(column(12,tableOutput('table'))),
    
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
  ### Remove functional probes, keep STR_SPE ###
  ##############################################
  humichip <- humichip %>%
    filter(gene == "STR_SPE")
  
  
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
    if(all(input$Visit_Number %in% "1")){
      filter(treat_tx(), visit_number == 1)
      
    } else if (all(input$Visit_Number %in% "4")){
      filter(treat_tx(), visit_number == 4)
      
    } else if (all(input$Visit_Number %in% "5")){
      filter(treat_tx(), visit_number == 5)
      
    } else if (all(input$Visit_Number %in% c("4","5"))){
      filter(treat_tx(), visit_number %in% c(4,5))
      
    } else if (all(input$Visit_Number %in% c("1","4"))){
      filter(treat_tx(), visit_number %in% c(1,4))
      
    } else if (all(input$Visit_Number %in% c("1","5"))){
      filter(treat_tx(), visit_number %in% c(1,5))
      
    } else if (is.null(input$Visit_Number)){
      return(treat_tx())
      
    } else {
      return(treat_tx())
    }
  })
  
  
  ##################
  #### TrEAT DB ####
  ##################
  
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
  
  
  #######################################
  ### Convert humichip to tidy format ###
  #######################################
  
  humichip_tidy <- reactive({
    humichip_select() %>%
      # Remove unwanted columns
      select(-Genbank.ID, -gene, -annotation, -geneCategory, -subcategory1, -subcategory2) %>%
      
      # Gather into long format
      gather(key = glomics_ID, value = signal, -species, -lineage) %>%
      
      # remove row if NA value present in signal
      filter(!is.na(signal)) %>%
      
      # Group by species and glomics ID, count probes and report average signal
      group_by(species, lineage, glomics_ID) %>%
      summarise(detected_probes = n(), avg_signal = mean(signal)) %>%
      
      # Merge in study_id's and visit number
      left_join(., ID_decoder, by = "glomics_ID") %>%
      
      # Merge in treat DB
      left_join(., treat, by = c("study_id" = "STUDY_ID", "visit_number", "glomics_ID")) %>%
      
      # Filter out samples based on Treat filtering
      filter(glomics_ID %in% treat_pathogen_filtered()$glomics_ID)
  })
  
  
  ########################
  ### Probe Thresholds ###
  ########################
  
  species_percent_probe <- reactive({
    
    humichip_tidy() %>%
    left_join(., total_probes, by = c("species" = "SPECIES")) %>%
    rename(total_probes = `Number of probes`) %>%
    filter(!is.na(total_probes)) %>%
    mutate(percent_probe = (detected_probes / total_probes) * 100) %>%
    select(species, lineage, glomics_ID, 
           study_id, visit_number, 
           detected_probes, total_probes, percent_probe, avg_signal, everything()) %>%
    filter(total_probes >= input$probe_number) %>% # How many probes each species must have
    filter(percent_probe >= input$probe_threshold) # what percent of probes must be present to consider positive
  })
  
  
  ##############################
  ### Heatmap data/factoring ###
  ##############################
  
  heatmap_data <- reactive({
    
    heatmap_filter_patients <- species_percent_probe() %>%
      select(species, lineage, study_id, visit_number, avg_signal, Treatment,
           Diarrhea_classification, LLS_severity, ESBL_V1, ESBL_V5, 
           Impact_of_illness_on_activity_level) %>%
      group_by(species) %>%
      filter(n() > input$species_threshold) %>% # how many samples must have species to be considered present
      mutate(Phylum = gsub(x = lineage, pattern = ".*;phylum:(\\w*\\s*\\w*);.*", replacement = "\\1")) # add phylum column
    
    # Align species by most to least abundant
    species_order <- heatmap_filter_patients %>%
      ungroup() %>%
      group_by(species) %>%
      summarise(n = n(), mean_sig = mean(avg_signal)) %>%
      arrange(desc(n)) %>%
      pull(species)
    
    # factor species order
    heatmap_filter_patients$species <- factor(heatmap_filter_patients$species, levels = rev(species_order))
    
    heatmap_filter_patients
  })
    
    
  
  # Get list of phylum colors
  heatmap_phylum_colors <- reactive({
    heatmap_data() %>%
      ungroup() %>%
      mutate(phylum_colors = ifelse(Phylum == "Actinobacteria", "red", 
                             ifelse(Phylum == "Bacteroidetes", "blue",
                             ifelse(Phylum == "Firmicutes", "green",
                             ifelse(Phylum == "Proteobacteria", "purple", "black"))))) %>%
      select(species, phylum_colors) %>%
      unique() %>%
      arrange(species) %>%
      pull(phylum_colors)
  })
  
  
  
  #####################
  ### Final Heatmap ###
  #####################
  
  plotInput <- reactive({
    
    # Set variables
    
    
    ggplot(heatmap_data(), aes(x = study_id, y = species)) +
      geom_tile(aes(fill = avg_signal)) +
      scale_y_discrete(expand=c(0,0)) +
      scale_x_discrete(expand=c(0,0)) +
      facet_grid(as.formula(paste("~", input$facet)), scales = "free") +
      scale_fill_viridis_c(name = "Normalized\nSignal\nIntensity", na.value = "black") +
      labs(caption = "Purple = Proteobacteria\nBlue = Bacteroidetes\nRed = Actinobacteria\nGreen = Firmicutes") +
      theme(
        plot.background = element_blank(), 
        panel.background = element_rect(fill = "black"),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, size = 7, vjust = 0.5),
        axis.text.y = element_text(color = heatmap_phylum_colors(), size = 12),
        plot.caption = element_text(hjust = 0),
        strip.text.x = element_text(size = 12)
      )
  })
  
  
  ####################
  ### Display Plot ###
  ####################
  output$plot <- renderPlot({
    print(plotInput())
  })
  
  #output$table <- renderTable({heatmap_data()})
  
  
  #####################
  ### Download plot ###
  #####################
  
  output$downloadPlot <- downloadHandler(
    filename = function(){paste("shiny_plot",'.png',sep='')},
    content = function(file){
      ggsave(file, plot=plotInput(), width = 17)
    })
  
}


shinyApp(server = server, ui = ui)