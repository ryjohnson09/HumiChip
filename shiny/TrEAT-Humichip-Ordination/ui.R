###############################
#### Humichip Shiny App UI ####
###############################

## Load Libraries -----------------------------------
library(shiny)


## Define UI ----------------------------------------
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Humichip Ordination"),
  
  # Sidebar 
  sidebarLayout(
  sidebarPanel(
      
  # Sample Selection
  fluidRow(column(12,
                  h3("Sample Selection"),
                  #############
                  ### Visit ###
                  #############
                  checkboxGroupInput('visit_number', 'Visit:', choices = list("Visit 1" = 1,"Visit 4" = 4,"Visit 5" = 5), 
                  selected = c(1, 4, 5), inline = TRUE),
                  
                  ###############
                  ### Matched ###
                  ###############
                  checkboxInput("matched", label = "Matched Samples?", value = TRUE),
                  
                  ########################
                  ### Country of Origin ##
                  ########################
                  checkboxGroupInput('country', 'Country:', 
                                     choices = c("Kenya", "Honduras", "Djibouti"), 
                                     selected = c("Kenya", "Honduras", "Djibouti"), inline = TRUE),
                  
                  ########################
                  ### Treatment Groups ###
                  ########################
                  checkboxGroupInput('treatment_groups', 'Treatment Groups:', 
                                     choices = c("LEV", "RIF", "AZI"), 
                                     selected = c("LEV", "RIF", "AZI"), inline = TRUE),
                  
                  ########################
                  ### Disease Severity ###
                  ########################
                  checkboxGroupInput("disease_severity", "Disease Severity:",
                                     choices = c("mild", "moderate", "severe"),
                                     selected = c("mild", "moderate", "severe")),
                  
                  ########################
                  ### Detection Method ###
                  ########################
                  br(),
                  h4("Select by Pathogen"),
                  checkboxInput("pathogen_select", label = "Select by Pathogen?", value = FALSE),
                  radioButtons("detection_method", "Detection Method:", 
                               choices = c("taq", "culture", "either", "both"),
                               selected = "both", inline = TRUE),
                  
                  ##########################
                  ### Pathogen Selection ###
                  ##########################
                  checkboxGroupInput("pathogens", "Select Pathogens:",
                                     choices = c("EAEC", "ETEC", "EPEC", "Shigella", "Norovirus"),
                                     selected = c("EAEC", "ETEC", "EPEC", "Shigella", "Norovirus"),
                                     inline = TRUE),
                  checkboxInput("allow_coinfections", label = "Allow Coinfections?", value = FALSE),
                  
                  ######################
                  ### ESBL Selection ###
                  ######################
                  br(),
                  h4("ESBL Selection"),
                  checkboxInput("esbl_select", label = "Select by ESBL?", value = FALSE),
                  checkboxGroupInput("esbls", label = "Which ESBLs?", choices = list(
                    "ESBL+ by phenotype, Visit 1" = "ESBL_V1",
                    "ESBL+ by phenotype, Visit 5" = "ESBL_V5",
                    "CMY, Visit 1" = "CMY_either_V1",
                    "CMY, Visit 5" = "CMY_either_V5",
                    "CTX, Visit 1" = "CTX_either_V1",
                    "CTX, Visit 5" = "CTX_either_V5",
                    "SHV, Visit 1" = "SHV_either_V1",
                    "SHV, Visit 5" = "SHV_either_V5",
                    "TEM, Visit 1" = "TEM_either_V1",
                    "TEM, Visit 5" = "TEM_either_V5")
                    , selected = "ESBL_V1") 
  )),
  
  # Probe Selection
  fluidRow(column(12,
                  h3("Probe Selection"),
                  
                  ##################
                  ### Probe Type ###
                  ##################
                  selectInput('probe_type', 'Probe Type:', choices = c("Functional", "Strain/Species", "All"), 
                              selected = "All"),
                  checkboxGroupInput("geneCategory", "Gene Category:", 
                                     choices = c("VIRULENCE", "PURINE_METABOLISM", 
                                                 "GLYCAN_BIOSYNTHESIS_AND_METABOLISM", 
                                                 "AMINO_ACID_SYNTHESIS", "EXOTIC_METABOLISMS", 
                                                 "COFACTOR_BIOSYNTHESIS", "GLYCOSAMINOGLYCAN_DEGRADATION", 
                                                 "FEEDER_PATHWAYS_TO_GLYCOLYSIS", 
                                                 "CENTRAL_CARBON_METABOLISM_PATHWAYS", "NITROGEN_METABOLISM", 
                                                 "GLYCAN_STRUCTURES_-_DEGRADATION", 
                                                 "AMINO_ACID_TRANSPORT_AND_METABOLISM", "ORGANIC_ACIDS", 
                                                 "PYRIMIDINE_METABOLISM", "ISOPRENOID_BIOSYNTHESIS", 
                                                 "COMPLEX_CARBOHYDRATES", "FATTY_ACID_METABOLISM", 
                                                 "GLYCEROLIPID_METABOLISM", "FATTY_ACID_BIOSYNTHESIS", 
                                                 "RESPIRATION", "N-GLYCAN_DEGRADATION"))
  )),
  
  # Ordination Type
  fluidRow(column(12,
                  h3("Ordination Type"),
                  
                  ##################
                  ### Ordination ###
                  ##################
                  radioButtons("ordination", "Ordination:", choices = c("PCA", "PCoA", "DCA"))
  )),
  
  
  # Plot Aesthetics
  fluidRow(column(12,
                  h3("Plot Aesthetics"),
                  
                  ####################
                  ### Point Colors ###
                  ####################
                  radioButtons("point_color", "Color Points By:", choices =  list(
                    "None" = "None",
                    "Visit" = "visit_number",
                    "Country" = "country",
                    "Treatment" = "Treatment",
                    "Pathogens" = "pathogens",
                    "Impact on Activity" = "Impact_of_illness_on_activity_level",
                    "Past 8 hours diarrhea severity" = "LLS_severity",
                    "ESBL Visit 1" = "ESBL_V1",
                    "ESBL Visit 5" = "ESBL_V5",
                    "ESBL Either 1 or 5" = "ESBL_either",
                    "TEM_V1" = "TEM_either_V1",
                    "TEM_V5" = "TEM_either_V5",
                    "CTX_V1" = "CTX_either_V1",
                    "CTX_V5" = "CTX_either_V5",
                    "CMY_V1" = "CMY_either_V1",
                    "CMY_V5" = "CMY_either_V5",
                    "SHV_V1" = "SHV_either_V1",
                    "SHV_V5" = "SHV_either_V5"),
                    selected = "Treatment"),
                  
                  ######################
                  ### Connect Points ###
                  ######################
                  radioButtons("lines", "Connect points by study ID?", 
                               choices = c("Yes", "No"), selected = "No")
  )),
                  
                  
                  
  
  # Perform Ordination Button
  actionButton("action", label = "Analyze"),
  
  # Download plot button
  downloadButton('downloadPlot','Download Plot'),

  # Side Bar Width
  width = 4),
    
  # Visualize Ordination Plot/data table
  mainPanel(
    plotOutput("humi_plot", width = "800px", height = "800px"),
    br(),
    verbatimTextOutput("stats_groupings"),
    br(),
    verbatimTextOutput("adonis_pvalue"),
    br(),
    verbatimTextOutput("mrpp_pvalue"),
    br(),
    fluidRow(column(12,tableOutput("humi_table")))
))))