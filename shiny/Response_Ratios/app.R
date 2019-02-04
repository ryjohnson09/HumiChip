# Humichip Response Ratio App

library(shiny)
library(tidyverse)

#######################
#### Read in Files ####
#######################

ID_decoder <- suppressWarnings(suppressMessages(read_csv("ID_Decoder_Humichip.csv")))
treat <- suppressWarnings(suppressMessages(read_csv("TrEAT_Clinical_Metadata_tidy.csv")))
humichip <- suppressWarnings(suppressMessages(read_tsv("Merged_humichip.tsv")))


######################
### Set up choices ###
######################

matched_choices <- list("All Samples" = "all_samples",
                        "Matched Samples" = "matched_samples")

tx_choices <- c("RIF", "LEV", "AZI", "All")

visit_choices <- list("Visit 1" = 1,
                      "Visit 4" = 4,
                      "Visit 5" = 5)

compare_choices <- c("visit_number", 
                     "ESBL_V1",
                     "ESBL_V5")

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
  
  titlePanel("Humichip Response Ratios"),
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
               checkboxGroupInput('Visit_Number', 'Visit:', choices = visit_choices, selected = c(1, 5), inline = TRUE),
               helpText("Select patient samples from specified visit number"),
               helpText("Only select two visit numbers"),
               
               # Comparing Groups
               radioButtons('compare_groups', 'Groups to Compare:', choices = compare_choices, selected = "visit_number", inline = TRUE),
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
               helpText("Select Categories to Compare"),
               
               # Number of minimum samples that must have Category
               sliderInput("cat_min_number", "Minimum Observation in Category:",
                           min = 1, max = 20,
                           value = 10)))),
    
    #######################
    ### Probe Selection ###
    #######################
    fluidRow(
      h3("Probe Selection"),
      
      column(12, 
             wellPanel(
               
               # Select Phylum probes?
               h5("Phylum Specific Probes"),
               checkboxInput("select_phylum", label = "Select Phylum Probes?", value = FALSE),
               helpText("If selected, only probes from selected bacterial phlya will be included in analysis"),
               
               # Phyla output
               uiOutput("phyla"),
               helpText("Select probes based on bacterial phyla"),
               
               # Should signal intensity be converted to 0 and 1
               checkboxInput("signal", "Convert probe signal to 0 or 1?", value = FALSE)
               ))),
    
    
    fluidRow(
      h3("Aesthetics"),
      
      column(12, 
             wellPanel(
               # Only keep significant RR's
               checkboxInput("keeper", label = "Remove non-significant RR's", value = FALSE)
               )),
    
      
      downloadButton('downloadPlot','Download Plot')),
    
    
    #sidebar width
    width = 4),
  
  # Plot
  mainPanel(
    plotOutput("plot", width = "800px", height = "800px"),
    br()
    
    # Table to see patients (not needed, but useful for troubleshooting)
    #fluidRow(column(12,tableOutput('table')))
    
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
  
  
  ########################
  ### Filter by Phylum ###
  ########################
  
  humichip_phylum <- reactive({
    
    # If filtering by phylum, and phylum column (only bacteria with phylum designation)
    if(input$select_phylum){
      humichip_select() %>%
        filter(str_detect(lineage, "Bacteria")) %>%
        filter(str_detect(lineage, ";phylum")) %>%
        mutate(Phylum = gsub(x = lineage, # Phylum column
                             pattern = ".*;phylum:(\\w*\\s*[-]*\\w*);.*", 
                             replacement = "\\1")) %>%
        filter(Phylum %in% input$phylum) # Filter step
        
    } else if (!input$select_phylum) {
      humichip_select() %>%
        filter(str_detect(lineage, "Bacteria")) %>%
        filter(str_detect(lineage, ";phylum")) %>%
        mutate(Phylum = gsub(x = lineage, # Phylum column
                             pattern = ".*;phylum:(\\w*\\s*[-]*\\w*);.*", 
                             replacement = "\\1"))
    } else {
      stopApp()
    }
  })
  
  #########################################
  ### Render list of Phyla if applicable ##
  #########################################
  output$phyla <- renderUI({
    
    if(input$select_phylum){
      checkboxGroupInput("phylum", "Select Phyla:",
                         choices = phylum_choices, inline = TRUE)
    }
  })
  
  
  ###############################################
  ### Select for Functional or STR_SPE probes ###
  ###############################################
  
  humichip_probe <- reactive({
    
    # If cat_choice == Phylum, filter for STR_SPE
    # and add Phylum category
    if(input$cat_choice == "Phylum"){
      humichip_phylum() %>%
        filter(gene == "STR_SPE")
      
    # If cat_choice == Functional group
    } else if (input$cat_choice != "Phylum"){
      humichip_phylum() %>%
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
  
  # Filter treat based on Visit Number
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
  ### And Make Long ##########################
  ############################################
  
  
  humichip_long <- reactive({
    # Subset Samples based on treat_pathogen_filtered
    humichip_probe() %>%
      select_if(colnames(.) %in% c("Genbank.ID", "gene", "species", "lineage",
                                   "annotation", "geneCategory", "subcategory1",
                                   "subcategory2", "Phylum",
                                   treat_pathogen_filtered()$glomics_ID)) %>%
      # Make Long
      gather(key = glomics_ID,
             value = Signal,
             -Genbank.ID, -gene,
             -species, -lineage, -annotation, -geneCategory, -subcategory1,
             -subcategory2, -Phylum) %>%
      
      # Merge in ID decoder
      left_join(., ID_decoder, by = c("glomics_ID")) %>%
      select(glomics_ID, study_id, visit_number, everything())
      
  })
  
  
  #######################################
  ### Rel Abundance by probe presence ###
  ### or probe signal ############### ###
  #######################################
  
  humichip_final <- reactive({
    if (input$signal){
      humichip_long() %>%
        mutate(Signal = ifelse(is.na(Signal), 0, 1))
    } else {
      humichip_long()
    }
  })
  
  #################################
  ### Calculate Rel Abundance #####
  ### Merge Humichip with TrEAT ###
  #################################
  humichip_treat <- reactive({
    
    y_axis_cat <- sym(input$cat_choice)
    grouping_cat <- sym(input$compare_groups)
    
    humichip_final() %>%
      
      group_by(glomics_ID) %>%
      mutate(Signal_Relative_Abundance = (Signal / sum(Signal, na.rm = TRUE)* 100)) %>%
      
      # Remove columns not needed
      select(glomics_ID, visit_number, !!y_axis_cat, Signal_Relative_Abundance) %>%
      
      
      # Remove any rows with NA in the signal category or y_axis_cat
      filter(!is.na(Signal_Relative_Abundance)) %>%
      filter(!is.na(!!y_axis_cat)) %>%
      
      
      # Calculate mead, sd, and counts (n)
      ungroup() %>%
      group_by(!!y_axis_cat, glomics_ID, visit_number) %>%
      summarise(cat_relative_abundance = sum(Signal_Relative_Abundance, na.rm = TRUE)) %>%
      
      # Add in TrEAT metadata
      left_join(., treat_pathogen_filtered(), by = c("glomics_ID", "visit_number")) %>%
      select(!!y_axis_cat, glomics_ID, visit_number, cat_relative_abundance, !!grouping_cat) %>%
      filter(!is.na(!!grouping_cat)) %>%
      ungroup()
  })
  
  
  
  ################################
  ### Calculate Response Ratio ###
  ################################
  

  
  humichip_RR <- reactive({
    
    y_axis_cat <- sym(input$cat_choice)
    grouping_cat <- sym(input$compare_groups)
    
     humichip_treat() %>%
      
      # Calculate mean, sd, and n
      group_by(!!y_axis_cat, !!grouping_cat) %>%
      summarise(mean_signal = mean(cat_relative_abundance),
                sd_signal = sd(cat_relative_abundance),
                n = sum(!is.na(cat_relative_abundance))) %>%
        
        
      # Spread the signal mean by visit number
      ungroup() %>%
      group_by(!!y_axis_cat) %>%
      spread(!!grouping_cat, mean_signal) %>%
       
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
      group_by(!!y_axis_cat) %>%
      summarise_all(funs(sum(., na.rm = T))) %>%

      # Must have at least __ observations in each subcategory
      filter(n_group1 >= input$cat_min_number) %>%
      filter(n_group2 >= input$cat_min_number) %>%

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
      mutate(pretty_cat = str_to_title(!!y_axis_cat)) %>%
      mutate(pretty_cat = str_replace_all(pretty_cat, "_"," ")) %>%
      
      # Factor columns
      mutate(pretty_cat = fct_reorder(pretty_cat, RR)) %>%
      ungroup()
  })
  
  
  ##########################################
  ### Filter for RR that don't overlap 0 ###
  ##########################################
  
  humichip_significant <- reactive({
    
    if (input$keeper){
      humichip_RR() %>%
        filter(keeper == "Yes")
    } else if (!input$keeper){
      humichip_RR()
    } else {
      stopApp("Failure at 'keeper' column creation")
    }
  })
      
    
    
  ##############################
  #### Create plots/tables #####
  ##############################
  
  plotInput <- reactive({
    
    y_axis_cat <- sym(input$cat_choice)
    grouping_cat <- sym(input$compare_groups)
    
    group_labels <- humichip_treat() %>%
      select(!!grouping_cat) %>%
      distinct() %>%
      pull(!!grouping_cat)
      
    
    ggplot(data = humichip_significant()) +
      geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
      
      # points and error bar
      geom_point(aes(x = pretty_cat, y = RR), size = 4) +
      geom_errorbar(aes(ymin = RR - CI95, 
                        ymax = RR + CI95, 
                        x = pretty_cat),
                    width = 0.25) +
      
      # Group labels
      annotate(geom = "text", label = group_labels[1], x = Inf, y = -Inf, hjust = 0, vjust = 1, 
               size = 5, color = "red", fontface = 2) +
      annotate(geom = "text", label = group_labels[2], x = Inf, y = Inf, hjust = 1, vjust = 1, 
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
    
  })
  
  
  
  ####################
  ### Display Plot ###
  ####################
  output$plot <- renderPlot({
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

  
  
  
  
  
  #output$table <- renderTable({head(humichip_significant(), 25)})
  
}







shinyApp(server = server, ui = ui)