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
                   "Visit 1 vs 5")

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
               helpText("Select Categories to Compare"),
               
               # Number of minimum samples that must have Category
               sliderInput("cat_min_number", "Minimum Observation in Category:",
                           min = 1, max = 20,
                           value = 10)))),
    
    
    fluidRow(
      h3("Ordination and Aesthetics"),
      
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
        filter(gene != "STR_SPE") %>%
        mutate(Phylum = NA)
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
  ### And Make Long ##########################
  ############################################
  
  
  humichip_final <- reactive({
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
      
      # Merge in ID_Decoder
      left_join(., ID_decoder, by = c("glomics_ID")) %>%
      select(glomics_ID, study_id, visit_number, everything())
  })
  
  
  
  
  
  
  ################################
  ### Calculate Response Ratio ###
  ################################
  

  
  humichip_RR <- reactive({
    
    y_axis_cat <- sym(input$cat_choice)
    
    humichip_final() %>%
    
      group_by(glomics_ID) %>%
      mutate(Signal_Relative_Abundance = (Signal / sum(Signal, na.rm = TRUE)* 100)) %>%
      
      # Remove columns not needed
      select(glomics_ID, visit_number, !!y_axis_cat, Signal_Relative_Abundance) %>%
    
      
      # Remove any rows with NA in the signal category or y_axis_cat
      filter(!is.na(Signal_Relative_Abundance)) %>%
      filter(!is.na(!!y_axis_cat)) %>%
       
       
      # Calculate mead, sd, and counts (n)
      group_by(!!y_axis_cat, glomics_ID, visit_number) %>%
      summarise(cat_relative_abundance = sum(Signal_Relative_Abundance, na.rm = TRUE)) %>%
      group_by(!!y_axis_cat, visit_number) %>%
      summarise(mean_signal = mean(cat_relative_abundance),
                sd_signal = sd(cat_relative_abundance),
                n = sum(!is.na(cat_relative_abundance))) %>%
       
       
      # Spread the signal mean by visit number
      group_by(!!y_axis_cat) %>%
      spread(visit_number, mean_signal) %>%
      
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
      group_by(!!y_axis_cat) %>%
      summarise_all(funs(sum(., na.rm = T))) %>%
      
      # Must have at least 10 observations in each subcategory
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
      mutate(pretty_cat = str_replace_all(pretty_cat, "_"," "))
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
    
    ggplot(data = humichip_significant()) +
      geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
      
      # points and error bar
      geom_point(aes(x = pretty_cat, y = RR), size = 4) +
      geom_errorbar(aes(ymin = RR - CI95, 
                        ymax = RR + CI95, 
                        x = pretty_cat),
                    width = 0.25) +
      
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

  
  
  
  
  
  output$table <- renderTable({head(humichip_significant(), 25)})
  
}







shinyApp(server = server, ui = ui)