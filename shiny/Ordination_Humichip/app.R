# Humichip TrEAT App

library(shiny)
library(tidyverse)
library(vegan)
library(ape)
library(ggExtra)


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

ordination_choices <- list("PCA" = "PCA",
                           "DCA" = "DCA", 
                           "PCoA" = "PCoA")

probe_choices <- c("Functional", "Strain/Species", "All")

geneCategory_choices <- c("All", unique(na.omit(humichip$geneCategory)))

visit_choices <- list("Visit 1" = 1,
                      "Visit 4" = 4,
                      "Visit 5" = 5)

tx_choices <- c("RIF", "LEV", "AZI", "All")

phylum_choices <- humichip %>%
  filter(str_detect(lineage, "Bacteria")) %>%
  filter(str_detect(lineage, ";phylum")) %>%
  mutate(Phylum = gsub(x = lineage, # Phylum column
                       pattern = ".*;phylum:(\\w*\\s*[-]*\\w*);.*", 
                       replacement = "\\1")) %>%
  select(Phylum) %>%
  distinct() %>%
  pull(Phylum)

detection_choices <- c("Culture", "Taq", "Both", "Either")

culture_choices <- c("All", colnames(treat)[grep(pattern = "_culture$", colnames(treat))])

taq_choices <- c("All", colnames(treat)[grep(pattern = "_taq$", colnames(treat))])

both_choices <- c("All", colnames(treat)[grep(pattern = "_both$", colnames(treat))])

either_choices <- c("All", colnames(treat)[grep(pattern = "_either$", colnames(treat))])

color_choices <- list(
  "None" = "None",
  "Visit" = "visit_number",
  "Treatment" = "Treatment",
  "Impact on Activity" = "Impact_of_illness_on_activity_level",
  "Past 8 hours diarrhea severity" = "LLS_severity",
  "Alternative Categories" = "Alternative_Categories",
  "ESBL Visit 1" = "ESBL_V1",
  "ESBL Visit 5" = "ESBL_V5")
  
color_choices1 <- c("Fever_present_at_presentation", "Vomiting_present", 
    "Abdominal_cramps_present_at_presentation", 
    "Excesssive_gas_flatulence_present_at_presentation", 
    "Nausea_present_at_presentation", 
    "Ineffective_and_or_paiful_straining_to_pass_a_stool_at_presentation", 
    "Tenesmus_present_at_presentation", "Malaise_present_at_presentation", 
    "Fecal_incontinence_present_at_presentation", 
    "Constipation_present_at_presentation", "Gross_blood_in_stool", 
    "Occult_blood_result")

point_size_choices <- list(
  "None" = "None",
  "TTLUS" = "Time_to_last_unformed_stool",
  "Loose Stools 8 hrs Prior to Presentation" = "Number_of_loose_liquid_stools_in_last_8_hours_prior_to_presentation",
  "Time to Cure" = "Time_to_cure",
  "Inverse Simpson" = "invsimpson",
  "Shannon Index" = "shannon"
)








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
  
  fluidRow(
    h3("Probe Selection"),
    
    column(12, 
           wellPanel(
             # Probe Type
             selectInput('probe', 'Probe Type:', choices = probe_choices, selected = "All"),
             helpText(code("Functional:"), " probes specific for functional genes", br(),
                      code("Strain/Species:"), " probes specific for microbial species and strains", br(),
                      code("All:"), " all probes included in analysis"),
             br(),
             # Probe Functional Category
             selectInput("geneCategory", "Probe Functional Category", choices = geneCategory_choices, selected = "All"),
             helpText("If ", code("Probe Type "), "= ", strong("Functional"), 
                      ", can select by functional group"),
             
             # Select Phylum probes?
             h5("Phylum Specific Probes"),
             checkboxInput("select_phylum", label = "Select Phylum Probes?", value = FALSE),
             helpText("If selected, only probes from selected bacterial phlya will be included in analysis"),
             
             # Phyla output
             uiOutput("phyla"),
             helpText("Select probes based on bacterial phyla")))),
  
  
  fluidRow(
    h3("Ordination and Aesthetics"),
    
    column(12, 
           wellPanel(
             # Ordination Type
             radioButtons("ordination", label = h3("Ordination Type:"),
                          choices = ordination_choices, selected = "PCA"),
    
             # Color Choices (main choices)
             radioButtons("color", label = h3("Color Points By:"),
                          choices = color_choices, inline = TRUE, selected = "visit_number"),
             # Color Choices (alternative)
             uiOutput("alt_cats"),
    
             # Point Size
             radioButtons("point_size", label = h3("Point Size"),
                          choices = point_size_choices, inline = TRUE, selected = "None"),
             
             # Connect points?
             checkboxInput("lines", label = "Connect Points by study id", value = FALSE))),
    
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
  
  
  #############################
  ### Functional or STR_SPE ###
  #############################
  
  humichip_probe <- reactive({
    if(input$probe == "Functional" & input$geneCategory == "All"){
      humichip_select() %>%
        filter(gene != "STR_SPE")
    } else if (input$probe == "Functional" & input$geneCategory != "All"){
      humichip_select() %>%
        filter(gene != "STR_SPE") %>%
        filter(geneCategory == input$geneCategory)
    } else if (input$probe == "Strain/Species"){
      humichip_select() %>%
        filter(gene == "STR_SPE")
    } else {
      humichip_select()
    }
  })
  
  
  ########################
  ### Filter by Phylum ###
  ########################
  
  humichip_phylum <- reactive({
    
    if(input$select_phylum){
      humichip_probe() %>%
        filter(!is.na(lineage)) %>%
        filter(str_detect(lineage, "Bacteria")) %>%
        filter(str_detect(lineage, ";phylum")) %>%
        mutate(Phylum = gsub(x = lineage, # Phylum column
                             pattern = ".*;phylum:(\\w*\\s*\\w*);.*", 
                             replacement = "\\1")) %>%
        filter(Phylum %in% input$phylum)
      
    } else if (!input$select_phylum) {
      humichip_probe()
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
  
  
  #######################
  ### Merged Humichip ###
  #######################
  
  # Filer Humichip to only include specified Visit_Number samples
   humichip_filtered <- reactive({
     humichip1 <- humichip_phylum()[,which(colnames(humichip_phylum()) %in% treat_pathogen_filtered()$glomics_ID)]
     
     # Set NA's to 0 and values not NA to original value
     humichip1 <- humichip1 %>%
       mutate_all(funs(ifelse(is.na(.), 0, .)))
     
     # Remove rows that equal 0
     humichip1 <- humichip1[rowSums(humichip1) != 0,]
     
     as.matrix(humichip1)
    })
  
  
  
  ################################
  ### Alpha Diveristy Analysis ###
  ################################
  
  # Remove all columns that are not samples
  humichip_alpha <- reactive({
    
    # Inverse Simpson
    humichip_invs <- as.tibble(diversity(t(humichip_filtered()), "inv")) %>%
      rownames_to_column("glomics_ID") %>%
      rename(invsimpson = value)
    
    # Shannon
    humichip_shann <- as.tibble(diversity(t(humichip_filtered()), "shannon")) %>%
      rownames_to_column("glomics_ID") %>%
      rename(shannon = value)
    
    # Merge into One
    humichip_alpha_diversity <- humichip_invs %>%
      full_join(., humichip_shann, by = c("glomics_ID"))
    
    humichip_alpha_diversity
  })
  
  
  #############################
  #### Ordination Analysis ####
  #############################
  
  # Ordination choices
  humichip_ord <- reactive({
    # Perform PCA analysis using vegan
    if(input$ordination == "PCA"){
      vegan::rda(t(humichip_filtered()))
    
    # Perform DCA analysis using vegan
    } else if (input$ordination == "DCA") {
      vegan::decorana(t(humichip_filtered()))
    
    # Perform PCoA using ape package
    } else if (input$ordination == "PCoA"){
      
      humichip_dist <- vegan::vegdist(as.matrix(t(humichip_filtered())))
      ape::pcoa(humichip_dist, correction = "none")
    } else {
      stopApp()
    }
  })
  
  
  humichip_coords <- reactive({
    
    # PCA and DCA coordinates
    if(input$ordination == "PCA" | input$ordination == "DCA"){
    
    # Get the coordinated for the samples
    ord_coords <- scores(humichip_ord(), display = "sites")

    # Make tibble
    ord_coords <- as.data.frame(ord_coords) %>%
      rownames_to_column(var = "glomics_ID")
    
    ord_coords
    
    # PCoA Coords
    } else if (input$ordination == "PCoA") {
      
    ord_coords <- humichip_ord()$vectors[,1:2]
    
    # Make tibble
    ord_coords <- as.data.frame(ord_coords) %>%
      rownames_to_column(var = "glomics_ID")
    
    ord_coords
    } else {
      stopApp()
    }
  })
  
  
  humichip_prop_exp <- reactive({
    
    # Extract proportion explained by first couple PC's for PCA
    if (input$ordination == "PCA"){
      summary(eigenvals(humichip_ord()))[2,] * 100
    } else if (input$ordination == "PCoA"){
      humichip_ord()$values$Relative_eig * 100
    } else {
      stopApp()
    }
  })
  
  
  ########################################################
  #### Merge coords and alpha diversity into to treat ####
  ########################################################
  
  humichip_coords_treat <- reactive({
    
    humichip_coords() %>%
      left_join(., treat_pathogen_filtered(), by = c("glomics_ID")) %>%
      right_join(., humichip_alpha(), by = c("glomics_ID")) %>%
      
    # Factor columns
      mutate(visit_number = factor(visit_number)) %>%
      mutate(Impact_of_illness_on_activity_level = factor(Impact_of_illness_on_activity_level))
  })

  
  
  
  ###################################
  ### Set up alt color categories ###
  ###################################
  
  output$alt_cats <- renderUI({
    if(input$color == "Alternative_Categories"){
      selectInput("color_alt", "Alternative Categories:", choices = color_choices1, width = "500px")
    }
  })
  
  
  ##############################
  #### Create plots/tables #####
  ##############################
  
  plotInput <- reactive({
    
    # Get color input for points
    if(input$color != "Alternative_Categories"){
      my_fill <- ifelse(input$color == "None", "NULL", input$color)
    } else {
      my_fill <- input$color_alt
    }
    
    # Get size input for points
    my_size <- ifelse(input$point_size == "None", "4", input$point_size)
    
    
    # Aesthetic sizes
    axis_title_size <- 16
    axis_text_size <- 14
    title_size <- 18
    legend_text_size <- 12
    
      
    # PCA Plot
    if(input$ordination == "PCA"){
    pca_plot <- ggplot(humichip_coords_treat(),
            aes_string(x = "PC1", 
                       y = "PC2", 
                       color = my_fill,
                       size = my_size)) +
        
      # Set up proportion explained
      xlab(paste0("PC1(", round(humichip_prop_exp()[[1]], 2), "%)")) +
      ylab(paste0("PC2(", round(humichip_prop_exp()[[2]], 2), "%)")) +
      geom_point(pch = 1, alpha = 1) +
      geom_point(pch = 19, alpha = 0.8) +
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
      guides(fill = guide_legend(override.aes = list(size=7))) # make legend points larger
    
    # Connect lines?
    if (input$lines == TRUE){
      pca_plot <- pca_plot + geom_line(aes(group = STUDY_ID), linetype = 1, color = "black", size = 0.8)
    }
    
    ggMarginal(pca_plot, groupColour = TRUE, groupFill = TRUE)
    
    
  } else if (input$ordination == "DCA"){
  
      
    # DCA plot
    dca_plot <- ggplot(humichip_coords_treat(),
           aes_string(x = "DCA1",
                      y = "DCA2", 
                      color = my_fill,
                      size = my_size)) +
      xlab("DCA1") +
      ylab("DCA2") +
      geom_point(pch = 1, alpha = 1) +
      geom_point(pch = 19, alpha = 0.8) +
      ggtitle("DCA Analysis") +
      theme_minimal() +
      theme(
        axis.title.x = element_text(size = axis_title_size),
        axis.text.x = element_text(size = axis_text_size, hjust = 1),
        axis.text.y = element_text(size = axis_text_size),
        axis.title.y = element_text(size = axis_title_size),
        plot.title = element_text(size = title_size, face = "bold"),
        legend.text = element_text(size = legend_text_size),
        legend.title = element_blank()) +
      guides(fill = guide_legend(override.aes = list(size=7))) # make legend points larger
    
    # Connect lines?
    if (input$lines == TRUE){
      dca_plot <- dca_plot + geom_line(aes(group = STUDY_ID), linetype = 1, color = "black", size = 0.8)
    }
    
    ggMarginal(dca_plot, groupColour = TRUE, groupFill = TRUE)
  
    
  } else if (input$ordination == "PCoA"){
    
    # PCoA Plot
    pcoa_plot <- ggplot(humichip_coords_treat(), 
                aes_string(x = "Axis.1", 
                           y = "Axis.2", 
                           color = my_fill, 
                           size = my_size)) +
      scale_fill_discrete(name = "Visit Number") +
      xlab(paste0("PCoA1(", round(humichip_prop_exp()[1], 2), "%)")) +
      ylab(paste0("PCoA2(", round(humichip_prop_exp()[2], 2), "%)")) +
      geom_point(pch = 1, alpha = 1) +
      geom_point(pch = 19, alpha = 0.8) + 
      ggtitle("PCoA Analysis") + 
      theme_minimal() +
      theme(
        axis.title.x = element_text(size = axis_title_size), 
        axis.text.x = element_text(size = axis_text_size, hjust = 1),
        axis.text.y = element_text(size = axis_text_size),
        axis.title.y = element_text(size = axis_title_size), 
        plot.title = element_text(size = title_size, face = "bold"),
        legend.text = element_text(size = legend_text_size),
        legend.title = element_blank()) +
      guides(fill = guide_legend(override.aes = list(size=5))) # make legend points larger
    
    # Connect lines?
    if (input$lines == TRUE){
      pcoa_plot <- pcoa_plot + geom_line(aes(group = STUDY_ID), linetype = 1, color = "black", size = 0.8)
    }
    
    ggMarginal(pcoa_plot, groupColour = TRUE, groupFill = TRUE)
    }
  })
  
  
  ####################
  ### Display Plot ###
  ####################
  output$plot <- renderPlot({
    print(plotInput())
  })
  
  #output$table <- renderTable({humichip_coords_treat()})
  
  #####################
  ### Download plot ###
  #####################
  
  output$downloadPlot <- downloadHandler(
    filename = function(){paste("shiny_plot",'.png',sep='')},
    content = function(file){
      ggsave(file, plot=plotInput())
    })
  
  
}



shinyApp(server = server, ui = ui)